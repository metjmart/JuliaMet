# *****************************************************************************
# cyclone.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.4.5
#
# This script contains functions that are more or less specific to tropical 
# cyclone related applications.
# 
# calc_rmw
# calc_azmean
# *****************************************************************************

#==============================================================================
calc_rmw

Determine the azimuthal mean radius of maximum tangential winds at each vertical 
level. Assumes that the input vt is the azimuthally averaged vt and has 
dimensions of r,z.
==============================================================================#

function calc_rmw{Ta<:Real,Tb<:Real,Tc<:Real}(r::AbstractVector{Ta},
                                              z::AbstractVector{Tb}, 
                                              azmean_vt::AbstractArray{Tc,2})
    # Create an array with the radial values at each height
    r_new = Array(Float64,length(r),length(z))
    for i in eachindex(z)
        r_new[:,i] = r
    end
    # Compute the azimuthal mean RMW at each vertical level
    rmw = similar(z,Float64)
    fill!(rmw, NaN)
    vtmax = findmax(azmean_vt,1)
    for k in eachindex(z)
        # Only store the RMW values if they're present
        if vtmax[2][k] != 0
            rmw[k] = r_new[vtmax[2][k]]
        end
    end

    return rmw

end

#==============================================================================
calc_azmean

Compute the azimuthal mean of the input variable. 
** Be sure to use the nanmean function since NaNs may be present
** var needs to be a three dimensional field that's indexed r,theta,z
==============================================================================#

function calc_azmean{Ta<:Real,Tb<:Real,Tc<:Real}(r::AbstractVector{Ta},
                                                 z::AbstractVector{Tb},
                                                 var::AbstractArray{Tc,3})
    # Define an array for the azimuthal mean
    azmean_var = Array(Float64,length(r),length(z)) 
    fill!(azmean_var, NaN)
    # Calculate the azimuthal mean 
    for i in eachindex(r)
        for j in eachindex(z)
            azmean_var[i,j] = nanmean(var[i,:,j])
        end
    end

    return azmean_var

end

#==============================================================================
uv2urvt

Convert u (east-west) and v (north-south) winds to radial and tangential winds
** Assumes u and v input winds are on an RTZ grid (radius,theta,height).
** Input theta vector assumed to be in units of degrees
==============================================================================#

function uv2urvt{Ta<:Real,Tb<:Real,Tc<:Real}(theta::AbstractVector{Ta},
                                             u::AbstractArray{Tb,3},
                                             v::AbstractArray{Tc,3})
    # Create the ur and vt arrays
    ur = Array(Float64,size(u))
    vt = Array(Float64,size(v))
    # Convert u,v to ur,vt
    for i in eachindex(theta)
        ur[:,i,:] =  u[:,i,:] .* cos(deg2rad(theta[i])) + 
                     v[:,i,:] .* sin(deg2rad(theta[i]))
        vt[:,i,:] = -u[:,i,:] .* sin(deg2rad(theta[i])) + 
                     v[:,i,:] .* cos(deg2rad(theta[i]))
    end
 
    return ur,vt

end
