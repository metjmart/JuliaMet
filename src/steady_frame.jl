# *****************************************************************************
# steady_frame.jl
#
# Author: Jonathan Martinez 
# Email: jon.martinez@colostate.edu
# Julia version: 0.4.5
#
# This script contains functions required to compute the most steady frame of 
# reference from Doppler radar data as described by Matejka (2002).
#
# Function list:
# stationary_frame
# More to come ...
# *****************************************************************************

#==============================================================================
stationary_frame

This function computes the most steady frame of reference for one stationary 
Doppler radar following Matejka (2002) equations (39) and (42). The function 
assumes that the radar is always located at the center of the grid (x,y = 0,0).
Three Cartesian gridded analyses are required which can be obtained from three 
sequential volumes using Radx2Grid.

Required input vars:
t1file = Gridded Doppler radar analysis for first chronological volume
t2file = Gridded Doppler radar analysis for second " "
t3file = Gridded Doppler radar analysis for third " "
U = 1d U-motion array, ideally with best guess in the specified range
    ** Set the resolution of the array to 1 m/s to iterate over large ranges
V = 1d V-motion array, " " 
zr = Radar altitude which can be obtained from metadata in soloii or solo3
   
Output vars:
Q = Matrix of values computed from (39) and (42) for U,V 
U = From input: U-velocity 
V = From input: V-velocity 
numer_out = Numerator of (42) evaluated
rho_norm = Normalizing factor (denominator) of (42) evaluated
==============================================================================#

function stationary_frame{Ta<:Real,Tb<:Real,Tc<:Real}(t1file::AbstractString,
                         t2file::AbstractString,t3file::AbstractString,
                         U::AbstractVector{Ta},V::AbstractVector{Tb},
                         zr::Tc)
                                                        
    # Define the required input vars for each time 
    # *** Only importing x,y,z once assuming they're identical 
    # for each time (same grid used, radar centered)

    t1vars = ["x0","y0","z0","start_time","stop_time","VG"]
    t2vars = ["start_time","stop_time","VG"]
    t3vars = ["start_time","stop_time","VG"]

    # Read in the vars 

    x,y,z,t1_start,t1_stop,vr1 = read_ncvars(t1file,t1vars)
    t2_start,t2_stop,vr2 = read_ncvars(t2file,t2vars)
    t3_start,t3_stop,vr3 = read_ncvars(t3file,t3vars)

    # Convert x,y,z to meters

    xm = x .* 1e3
    ym = y .* 1e3
    zm = z .* 1e3

    # Create 3d x,y,z arrays

    x3d,y3d,z3d = grid_3d(x,y,z)    

    # Convert the 3d array units to meters

    x3dm = x3d .* 1e3
    y3dm = y3d .* 1e3
    z3dm = z3d .* 1e3

    # Compute the time of each analysis as the mid point between start/end

    t1 = (t1_start + t1_stop) / 2
    t2 = (t2_start + t2_stop) / 2
    t3 = (t3_start + t3_stop) / 2

    # Compute the delta_t's

    delta_t1 = t1 - t2
    delta_t3 = t3 - t2

    # Define the location of the radar 
    # ** Assumes that radar is at center of grid

    xr = 0
    yr = 0

    # Get zr from metadata in solo and import into function at start

    # z2 only needs to be defined once and is not dependent on U,V 
    # see equation (13)

    z2 = z3d

    # Compute the expected squared error terms
    # See last paragragph on pg 1039 for reference 

    sigsquared_vr = 1.0
    sigsquared_udotdot = 10.0 ^ -12

    # Define the normalizing factor as the volume in equation (39)

    vol_norm = length(x) * length(y) * length(z)

    # Compute the density weighing factor
    # Assume density only varies w/ z, use isothermal atmosphere
    # Should be sufficient as it adds less weight at upper levels

    rho_o = 1.225
    H = 8000 # Scale height in meters
    rho = rho_o .* exp(-z3dm/H)

    # Now calculate the normalizing factor in equation (42)
    # Using trapezoidal integration, make sure x,y,z are in meters!
 
    rho_norm = trapz_3d(xm,ym,zm,rho)

    # Create the array for Q and compute it for a subset of U and V

    Q = Array(Float64,length(U),length(V))

    # Output array for numerator in equation (42)

    numer_out = Array(Float64,length(U),length(V))

    for i in eachindex(U)
        for j in eachindex(V) 

            # Compute x, y, and z for each analysis time
            # See equation (13)
            # ** Can use same x and y assuming all grids are identical
            
            x1 = x3dm .+ U[i]*(t1-t2)
            x2 = x3dm .+ U[i]*(t2-t2)
            x3 = x3dm .+ U[i]*(t3-t2)
    
            y1 = y3dm .+ V[j]*(t1-t2)
            y2 = y3dm .+ V[j]*(t2-t2)
            y3 = y3dm .+ V[j]*(t3-t2)
            # z2 was defined above

            # Compute R at each time (see beneath equation 5)
    
            R1 = ((x1 - xr).^2) + ((y1 - yr).^2) + ((z2 - zr).^2)
            R2 = ((x2 - xr).^2) + ((y2 - yr).^2) + ((z2 - zr).^2)
            R3 = ((x3 - xr).^2) + ((y3 - yr).^2) + ((z2 - zr).^2)
            
            # Compute the p_hats (equation 40)
    
            p_hat12 = -delta_t1 ./ (R1 .* R2)
            p_hat13 = (delta_t3 - delta_t1) ./ (R1 .* R3)
            p_hat23 = delta_t3 ./ (R2 .* R3)
            
            # Compute the alpha_hats (equation 41)
    
            alpha_hatx = (delta_t1 .* delta_t3 .* (delta_t1 - delta_t3) .* 
                         ( (x2 - xr) .+ (U[i] .* (delta_t1 + delta_t3)) )) ./ 
                         (2.0 .* R1 .* R2 .* R3)

            alpha_haty = (delta_t1 .* delta_t3 .* (delta_t1 - delta_t3) .* 
                         ( (y2 - yr) .+ (V[j] .* (delta_t1 + delta_t3)) )) ./ 
                         (2.0 .* R1 .* R2 .* R3)
    
            alpha_hatz = (delta_t1 .* delta_t3 .* (delta_t1 - delta_t3) .* z2) ./ 
                         (2.0 .* R1 .* R2 .* R3)
            
            # Evaluate the numerator in (39) and include density as in (42)
            # ** Not integrating yet!!!
 
            numer_intg = ( rho .* ( ((vr1 .* p_hat23) - (vr2 .* p_hat13) + (vr3 .* p_hat12)) .^2 ) ) ./ 
                         ( (sigsquared_vr .* (p_hat23 .^2 + p_hat13 .^2 + p_hat12 .^2) ) + 
                         (sigsquared_udotdot .* (alpha_hatx .^ 2 + alpha_haty .^2 + alpha_hatz .^2)) )

            # Integrate the numerator using the trapezoidal method
            # Make sure x,y,z are in meters!

            numer = trapz_3d(xm,ym,zm,numer_intg)

            # Store the numerator for each U,V

            numer_out[i,j] = numer

            # Compute Q for each U and V
            # See equations (39) and (42)
    
            Q[i,j] = numer / rho_norm
            
        end
    end 

    # Create 2D arrays for U and V since Q is 2D
    
    U2d,V2d = grid_2d(U,V)
    
    # Find the index of minimum Q
    
    ind_Qmin = findin(Q,minimum(Q))[1]

    # Find the U,V corresponding to the minimum Q
    # This is the most steady frame of reference!
    
    minU = U2d[ind_Qmin][1]
    minV = V2d[ind_Qmin][1]
    
    println(" ")
    println("Q is minimized for U = ", minU, " and V = ", minV)
    
    return(Q,U,V,numer_out,rho_norm)

end



