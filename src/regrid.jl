# *****************************************************************************
# regrid.jl
#
# Author: Jonathan Martinez 
# Email: jon.martinez@colostate.edu
# Julia version: 0.4.5
#
# This script contains functions to regrid data using various methods.
#
# Function list:
# grid_2d
# grid_3d
# regrid_xy2rt 
# regrid_xyz2rtz
# *****************************************************************************

#==============================================================================
# grid_2d

This function is the equivalent of Matlab's ndgrid. It creates a two
dimensional grid from the one-dimensional input coordinate arrays.
** x and y do not have to be the same size.
==============================================================================#

function grid_2d{Ta<:Real,Tb<:Real}(x::Vector{Ta},y::Vector{Tb})
    # Create new arrays
    x_2d = Array(Float64,length(x),length(y))
    y_2d = Array(Float64,length(x),length(y))
    # 2d x-array
    for i in eachindex(x)
        x_2d[i,:] = x[i]
    end
    # 2d y-array
    for j in eachindex(y)
        y_2d[:,j] = y[j]
    end
    
    return x_2d,y_2d
end

#==============================================================================
# grid_3d

3D formulation of grid_2d
** x, y, and z do not have to be the same size.
==============================================================================#

function grid_3d{Ta<:Real,Tb<:Real,Tc<:Real}(x::Vector{Ta},y::Vector{Tb},
                                             z::Vector{Tc})
    # Create new arrays
    x_3d = Array(Float64,length(x),length(y),length(z))
    y_3d = Array(Float64,length(x),length(y),length(z))
    z_3d = Array(Float64,length(x),length(y),length(z))
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
# regrid_xy2rt

This function takes two-dimensional Cartesian data and interpolates it to a 
polar coordinate grid. It is specifically designed to work with NetCDF data
by working around issues with fill values when interpolating the data.
** Follows http://mathworld.wolfram.com/PolarCoordinates.html for converting 
   Cartesian to polar coordinates
 
** The CoordInterpGrid function can handle the presence of NaNs

rt_out - Option to output r and theta arrays, default is false.
==============================================================================#

function regrid_xy2rt{Ta<:Real,Tb<:Real,Tc<:Real}(x::Vector{Ta},y::Vector{Tb},
                                                  vardata::Array{Tc,2},
                                                  rt_out::Bool=false)
    # Determine the max radius of polar coordinate grid and the theta increment
    rmax = ceil(sqrt(((x[end]-x[1])/2.0)^2 + ((y[end]-y[1])/2.0)^2))
    theta_inc = floor(atan2(y[2]-y[1],(x[end]-x[1])/2.0)/pi*180.0)
    if theta_inc<1.0
        theta_inc=1.0
    end
    # Define r and theta
    r = collect(0:(x[2]-x[1]):rmax)
    theta = collect(0:theta_inc:360-theta_inc)
    # Create two-dimensional arrays for r and theta components of the polar grid
    # Define the Cartesian points in terms of the 2-d polar coordinates
    r_2d,theta_2d = grid_2d(r,theta)
    x_polar = r_2d .* cos(deg2rad(theta_2d))
    y_polar = r_2d .* sin(deg2rad(theta_2d))
    # Specify the x and y range of the data
    x_range = x[1]:(x[2]-x[1]):x[end]
    y_range = y[1]:(y[2]-y[1]):y[end]
    # Interpolate the data from the Cartesian grid to the polar grid 
    field_rt = Array(Float64,size(x_polar))
    field_interp = CoordInterpGrid((x_range,y_range), vardata, BCnan, InterpLinear);
    for i in eachindex(r)
        for j in eachindex(theta)
             field_rt[i,j] = field_interp[x_polar[i,j],y_polar[i,j]]
        end
    end
    # Return field_rt, r, and theta if rt_out is true
    if rt_out == true
        return r,theta,field_rt
    else
        return field_rt
    end
end

#==============================================================================
regrid_xyz2rtz

This function interpolates two-dimensional Cartesian data to a cylindrical 
coordinate grid by executing regrid_xy2rt at each vertical level.

rt_out - Option to output r and theta arrays, default is false.
==============================================================================#

function regrid_xyz2rtz{Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}(
                        x::Vector{Ta},
                        y::Vector{Tb},
                        z::Vector{Tc},
                        vardata::Array{Td,3},
                        rt_out::Bool=false)
    # Determine the max radius of polar coordinate grid and the theta increment
    rmax = ceil(sqrt(((x[end]-x[1])/2.0)^2 + ((y[end]-y[1])/2.0)^2))
    theta_inc = floor(atan2(y[2]-y[1],(x[end]-x[1])/2.0)/pi*180.0)
    if theta_inc<1.0
        theta_inc=1.0
    end
    # Define r and theta
    r = collect(0:(x[2]-x[1]):rmax)
    theta = collect(0:theta_inc:360-theta_inc)
   # Define the dimensions of the new var
    field_rtz = Array(Float64,length(r),length(theta),length(z))
    # Call regridxy2rt at each vertical level 
    for k in eachindex(z)
        field_rtz[:,:,k] = regrid_xy2rt(x,y,vardata[:,:,k])
    end
    # Return field_rtz, r, and theta if rt_out is true
    if rt_out == true
        return r,theta,field_rtz
    else
        return field_rtz
    end
end


