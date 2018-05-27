# *****************************************************************************
# steady_frame.jl
#
# Author: Jonathan Martinez 
# Email: jon.martinez@colostate.edu
# Julia version: 0.6.0
#
# This script contains functions required to compute the most steady frame of 
# reference from Doppler radar data as described by Matejka (2002).
#
# Function list:
# steadyframe
# *****************************************************************************

#==============================================================================
steadyframe

This function computes the most steady frame of reference given Doppler radar
data following Matejka (2002).
Three Cartesian gridded analyses are required which can be obtained from three 
sequential volumes using Radx2Grid.

** Note **
The function assumes that each of the analyses are on the same grid (i.e., use
the same Cartesian coordinates).

Currenlty supporting stationary radars with support for airborne platforms in 
progress.

Required input vars

x - Cartesian x-coordinate array (meters)
y - Cartesian y-coordinate array (meters)
z - Cartesian z-coordinate array (meters)
delta_t1 - Time of first analysis minus time of second analysis (seconds)
delta_t3 - Time of third analysis minus time of second analysis
vr1 - Gridded Doppler velocities for the first analysis (m/s)
vr2 - Gridded Doppler velocities for the second analysis (m/s)
vr3 - Gridded Doppler velocities for the third analysis (m/s)
U - 1d U-motion array, ideally with best guess in the specified range
V - 1d V-motion array, " " 
xr - Cartesian x-location of the radar (meters)
yr - Cartesian y-location of the radar (meters)
zr - Radar altitude which can be obtained from metadata in soloii, solo3, or
     online (meters)  
rho - 3-D density for the region covered by the grid. A simple density profile
      such as assuming an isothermal atmosphere should suffice given that the 
      sole purpose of this variable is to give less weight to observations at 
      high altitudes

Output vars:

Q = Matrix of values computed following eq. (42) for either mobile (eq. 36) or 
    stationary (eq. 39) platforms 
numer_out = Numerator of eq. (42) computed following eq. (42) for either mobile 
            (eq. 36) or stationary (eq. 39) platforms 
==============================================================================#

# Stationary radars

function steadyframe(x::AbstractVector{<:Real},y::AbstractVector{<:Real},z::AbstractVector{<:Real},
                     delta_t1::Real,delta_t3::Real,vr1::AbstractArray{<:Real,3},
                     vr2::AbstractArray{<:Real,3},vr3::AbstractArray{<:Real,3},
                     U::AbstractVector{<:Real},V::AbstractVector{<:Real},xr::Real,yr::Real,zr::Real,
                     rho::AbstractArray{<:Real,3})
                                                        
    # Get zr from metadata in solo and import into function at start
    # Compute the expected squared error terms
    # See last paragragph on pg 1039 for reference 
    sigsq_vr = 1.0
    sigsq_udotdot = 1e-12
    # Allocate necessary arrays -- this can't be the most efficient way to do this..
    Q = Array{Float64}(length(U),length(V))
    numer_out = Array{Float64}(length(U),length(V))
    chisq = Array{Float64}(length(x),length(y),length(z))
    sigsq = Array{Float64}(length(x),length(y),length(z))
    for v in eachindex(V)
        for u in eachindex(U) 
            # Compute x_i and y_i for each analysis time (z stays the same)
            # See equation (13)
            # ** Can use same x and y assuming all grids are identical
            x1 = x + U[u] * delta_t1
            x2 = x + U[u] * 0.
            x3 = x + U[u] * delta_t3
            y1 = y + V[v] * delta_t1
            y2 = y + V[v] * 0.
            y3 = y + V[v] * delta_t3
            # Compute R at each time (see beneath equation 5)
            for k in eachindex(z)
                for j in eachindex(y) 
                    for i in eachindex(x)
                        R1 = sqrt((x1[i] - xr)^2 + (y1[j] - yr)^2 + (z[k] - zr)^2)
                        R2 = sqrt((x2[i] - xr)^2 + (y2[j] - yr)^2 + (z[k] - zr)^2)
                        R3 = sqrt((x3[i] - xr)^2 + (y3[j] - yr)^2 + (z[k] - zr)^2)
                        # Compute the p_hats (equation 40)
                        p_hat12 = -delta_t1 / (R1 * R2)
                        p_hat13 = (delta_t3 - delta_t1) / (R1 * R3)
                        p_hat23 = delta_t3 / (R2 * R3)
                        # Compute the alpha_hats (equation 41)
                        alpha_hatx = (delta_t1 * delta_t3 * (delta_t1 - delta_t3) * 
                                     ( (x2[i] - xr) + (U[u] * (delta_t1 + delta_t3)) )) / 
                                     (2.0 * R1 * R2 * R3)
                        alpha_haty = (delta_t1 * delta_t3 * (delta_t1 - delta_t3) * 
                                     ( (y2[j] - yr) + (V[v] * (delta_t1 + delta_t3)) )) / 
                                     (2.0 * R1 * R2 * R3)
                        alpha_hatz = (delta_t1 * delta_t3 * (delta_t1 - delta_t3) * z[k]) / 
                                     (2.0 * R1 * R2 * R3)
                        # Compute the chi squared and sigma squared terms in (42) -- refer to (39)
                        chisq[i,j,k] = (vr1[i,j,k] * p_hat23 - vr2[i,j,k] * p_hat13 + vr3[i,j,k] * p_hat12)^2 
                        sigsq[i,j,k] = sigsq_vr * (p_hat23^2 + p_hat13^2 + p_hat12^2) + 
                                       sigsq_udotdot * (alpha_hatx^2 + alpha_haty^2 + alpha_hatz^2)
                    end
                end
            end
            # Compute the numerator in (42) 
            numer = rho .* chisq ./ sigsq
            # Use numer to mask values in rho where NaNs are present
            rho[findin(numer,NaN)] = NaN
            # Integrate the numerator using the trapezoidal method 
            numer_out[u,v] = trapz3d(x,y,z,numer)
            # Integrate the denominator 
            rho_norm = trapz3d(x,y,z,rho)
            # Compute Q for each U and V (see equations 39 and 42) 
            Q[u,v] = numer_out[u,v] / rho_norm
        end
    end 
    # Create 2D arrays for U and V since Q is 2D
    U2d,V2d = grid2d(U,V)
    # Find the index of minimum Q
    ind_Qmin = findmin(Q)[2]
    # Find the U,V corresponding to the minimum Q
    # This is the most steady frame of reference!
    minU = U2d[ind_Qmin]
    minV = V2d[ind_Qmin]
    println("Q is minimized for U = ", minU, " and V = ", minV)
    return(minU,minV,Q,numer_out)
end


