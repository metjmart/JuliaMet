# *****************************************************************************
# regrid.jl
#
# Author: Jonathan Martinez 
# Email: jon.martinez@colostate.edu
# Julia version: 0.6.0
#
# This script contains functions to regrid data using various methods.
#
# Function list:
# closest_ind
# grid2d
# grid3d
# xy2rp
# regrid_xy2rp 
# regrid_xyz2rpz
# regrid_pol2cart
# *****************************************************************************

#==============================================================================
# closest_ind

Search for the index within a 1-D array which most closely corresponds to the 
specified value. 
==============================================================================#

function closest_ind(arr::AbstractVector{<:Real},val::Real)
    idiff = abs.(arr-val)
    return findin(idiff,minimum(idiff))[1]
end

#==============================================================================
# grid2d

This function is the equivalent of Matlab's ndgrid. It creates a two
dimensional grid from the one-dimensional input coordinate arrays.
** x and y do not have to be the same size.
==============================================================================#

function grid2d(x::AbstractVector{<:Real},y::AbstractVector{<:Real})
    return @fastmath x .* ones(y)', ones(x) .* y'
end

#==============================================================================
# grid3d

3D formulation of grid2d
** x, y, and z do not have to be the same size.
==============================================================================#

function grid3d(x::AbstractVector{<:Real},y::AbstractVector{<:Real},
                z::AbstractVector{<:Real})

    # Create new arrays
    x_3d = Array{Float64}(length(x),length(y),length(z))
    y_3d = Array{Float64}(length(x),length(y),length(z))
    z_3d = Array{Float64}(length(x),length(y),length(z))
    # 3d x-array
    for i in eachindex(x)
        x_3d[i,:,:] = x[i]
    end
    # 3d y-array
    for j in eachindex(y)
        y_3d[:,j,:] = y[j]
    end
    # 3d z-array
    for k in eachindex(z)
        z_3d[:,:,k] = z[k]
    end
    return x_3d,y_3d,z_3d
end

#==============================================================================
xy2rp

This function will convert x and y arrays to radius and phi arrays given the 
center of a TC. See center_finding.jl for methods to determine the center.
The purpose of requiring the center is if there is stretching in the domain,
this will use the x and y grid spacing nearest to the center (i.e., without 
stretching). 
==============================================================================#

function xy2rp(cx::Real,cy::Real,x::AbstractVector{<:Real},
               y::AbstractVector{<:Real})

    # Pinpoint the center of grid that has no stretching
    # Has no impact on grids with constant spacing
    x0 = closest_ind(x - cx,0.0)
    y0 = closest_ind(y - cy,0.0)
    # Determine the max radius of polar coordinate grid and the phi increment
    rmax = ceil(sqrt((maximum(abs.(x)))^2 + (maximum(abs.(y))^2)))
    phi_inc = floor(atan2(y[y0]-y[y0-1],maximum(abs.(x))))
    phi_inc < pi/180 ? phi_inc = pi/180 : nothing
    # Return radius and azimuth arrays
    return collect(0:(x[x0]-x[x0-1]):rmax), collect(0:phi_inc:2*pi - phi_inc)
end

#==============================================================================
# regrid_xy2rp

This function takes two-dimensional Cartesian data and interpolates it to a 
polar coordinate grid. It is specifically designed to work with NetCDF data
by working around issues with fill values when interpolating the data.
** Follows http://mathworld.wolfram.com/PolarCoordinates.html for converting 
   Cartesian to polar coordinates
 
rp_out - Option to output r and phi arrays. Default is false.
==============================================================================#

function regrid_xy2rp(cx::Real,cy::Real,x::AbstractVector{<:Real},
                      y::AbstractVector{<:Real},field::AbstractArray{<:Real,2})

    # Define r and phi
    r,phi = xy2rp(cx,cy,x,y)
    # Interpolate the data from the Cartesian grid to the polar grid 
    field_rp = Array{Float64}(length(r),length(phi))
    field_itp = extrapolate(interpolate((x-cx, y-cy),field,Gridded(Linear())),NaN)
    for j in eachindex(phi)
        for i in eachindex(r)
             @inbounds field_rp[i,j] = field_itp[r[i] * cos(phi[j]), r[i] * sin(phi[j])]
        end
    end
    return field_rp
end

#==============================================================================
regrid_xyz2rpz

This function interpolates three-dimensional Cartesian data to a cylindrical 
coordinate grid by executing regrid_xy2rp at each vertical level.

rp_out - Option to output r and phi arrays, default is false.
==============================================================================#

function regrid_xyz2rpz(cx::Real,cy::Real,x::AbstractVector{<:Real},
                        y::AbstractVector{<:Real},z::AbstractVector{<:Real},
                        field::AbstractArray{<:Real,3})

    # Define r and phi
    r,phi = xy2rp(cx,cy,x,y)
    # Define the dimensions of the new var
    field_rpz = Array{Float64}(length(r),length(phi),length(z))
    # Call regrid_xy2rp at each vertical level 
    for k in eachindex(z)
        @inbounds field_rpz[:,:,k] = regrid_xy2rp(cx,cy,x,y,field[:,:,k])
    end
    # Return field_rpz
    return field_rpz
end

#==============================================================================
regrid_pol2cart

This function will regrid data on a polar grid to Cartesian coordinates
*** Assume data is centered on polar grid 
*** Currently not handling presence of NaNs -- use caution!
==============================================================================#

function regrid_pol2cart(r::AbstractVector{<:Real},phi::AbstractVector{<:Real},
                         field::AbstractArray{<:Real,2})

    x = collect(-r[end]:r[2]-r[1]:r[end])
    y = collect(-r[end]:r[2]-r[1]:r[end])
    # Create 2-D grids
    xx,yy = grid2d(x,y)
    # Transform Cartesian to polar reference points for azimuth
    # due to atan2 returning negative values
    pp_cart = atan2.(yy,xx)
    pp_cart[pp_cart .< 0] = pp_cart[pp_cart .< 0] + 2 * pi
    # Create the interpolation object
    field_xy = Array{Float64}(length(x),length(y))
    field_itp = extrapolate(interpolate((r,phi),field,Gridded(Linear())),NaN)
    # Interpolate from polar to Cartesian
    for j in eachindex(y)
        for i in eachindex(x)
            @inbounds field_xy[i,j] = field_itp[sqrt(x[i]^2 + y[j]^2), pp_cart[i,j]]
        end
    end
    return field_xy
end


