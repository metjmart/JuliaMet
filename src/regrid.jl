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
# xy2rt
# regrid_xy2rt 
# regrid_xyz2rtz
# regrid_gfrelxz
# regrid_gfrelxyz
# *****************************************************************************

#==============================================================================
# closest_ind

Search for the index within a 1-D array which most closely corresponds to the 
specified value. 
==============================================================================#

function closest_ind(arr::AbstractVector{<:Real},val::Real)

        idiff = abs.(arr-val)
        ibest = findin(idiff,minimum(idiff))[1]

        return ibest
end

#==============================================================================
# grid2d

This function is the equivalent of Matlab's ndgrid. It creates a two
dimensional grid from the one-dimensional input coordinate arrays.
** x and y do not have to be the same size.
==============================================================================#

function grid2d(x::AbstractVector{<:Real},y::AbstractVector{<:Real})

    # Create new arrays
    x_2d = Array{Float64}(length(x),length(y))
    y_2d = Array{Float64}(length(x),length(y))
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
xy2rt

This function will convert x and y arrays to radius and theta arrays given the 
center of a TC. See center_finding.jl for methods to determine the center.
The purpose of requiring the center is if there is stretching in the domain,
this will use the x and y grid spacing nearest to the center (i.e., without 
stretching). 
==============================================================================#

function xy2rt(cx::Real,cy::Real,x::AbstractVector{<:Real},
               y::AbstractVector{<:Real})

    # Re-define x and y arrays based on specified center
    xn = x - cx
    yn = y - cy
    # Pinpoint the center of grid that has no stretching
    # Has no impact on grids with constant spacing
    x0 = closest_ind(x,0.0)
    y0 = closest_ind(y,0.0)
    # Determine the max radius of polar coordinate grid and the theta increment
    rmax = ceil(sqrt((maximum(abs.(x)))^2 + (maximum(abs.(y))^2)))
    theta_inc = floor(atan2(y[y0]-y[y0-1],maximum(abs.(x)))/pi*180.0)
    if theta_inc<1.0
        theta_inc=1.0
    end
    # Define r and theta
    r = collect(0:(x[x0]-x[x0-1]):rmax)
    theta = collect(0:theta_inc:360-theta_inc)

    return r,theta

end

#==============================================================================
# regrid_xy2rt

This function takes two-dimensional Cartesian data and interpolates it to a 
polar coordinate grid. It is specifically designed to work with NetCDF data
by working around issues with fill values when interpolating the data.
** Follows http://mathworld.wolfram.com/PolarCoordinates.html for converting 
   Cartesian to polar coordinates
 
rt_out - Option to output r and theta arrays. Default is false.
==============================================================================#

function regrid_xy2rt(cx::Real,cy::Real,x::AbstractVector{<:Real},
                      y::AbstractVector{<:Real},var::AbstractArray{<:Real,2},
                      rt_out::Bool=false)

    # Re-define x and y arrays based on specified center
    xn = x - cx
    yn = y - cy
    # Pinpoint the center of grid that has no stretching
    # Has no impact on grids with constant spacing
    x0 = closest_ind(x,0.0)
    y0 = closest_ind(y,0.0)
    # Determine the max radius of polar coordinate grid and the theta increment
    rmax = ceil(sqrt((maximum(abs.(x)))^2 + (maximum(abs.(y))^2)))
    theta_inc = floor(atan2(y[y0]-y[y0-1],maximum(abs.(x)))/pi*180.0)
    if theta_inc<1.0
        theta_inc=1.0
    end
    # Define r and theta
    r = collect(0:(x[x0]-x[x0-1]):rmax)
    theta = collect(0:theta_inc:360-theta_inc)
    # Create two-dimensional arrays for r and theta components of the polar grid
    # Define the Cartesian points in terms of the 2-d polar coordinates
    r_2d,theta_2d = grid2d(r,theta)
    x_polar = r_2d .* cos.(deg2rad.(theta_2d))
    y_polar = r_2d .* sin.(deg2rad.(theta_2d))
    # Interpolate the data from the Cartesian grid to the polar grid 
    field_rt = Array{Float64}(size(x_polar))
    field_interp = extrapolate(interpolate((xn,yn), var, Gridded(Linear())), NaN)
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

function regrid_xyz2rtz(cx::Real,cy::Real,x::AbstractVector{<:Real},
                        y::AbstractVector{<:Real},z::AbstractVector{<:Real},
                        vardata::AbstractArray{<:Real,3},rt_out::Bool=false)

    # Re-define x and y arrays based on specified center
    xn = x - cx
    yn = y - cy
    # Pinpoint the center of grid that has no stretching
    # Has no impact on grids with constant spacing
    x0 = closest_ind(x,0.0)
    y0 = closest_ind(y,0.0)
    # Determine the max radius of polar coordinate grid and the theta increment
    rmax = ceil(sqrt((maximum(abs.(x)))^2 + (maximum(abs.(y))^2)))
    theta_inc = floor(atan2(y[y0]-y[y0-1],maximum(abs.(x)))/pi*180.0)
    if theta_inc<1.0
        theta_inc=1.0
    end
    # Define r and theta
    r = collect(0:(x[x0]-x[x0-1]):rmax)
    theta = collect(0:theta_inc:360-theta_inc)
    # Define the dimensions of the new var
    field_rtz = Array{Float64}(length(r),length(theta),length(z))
    # Call regridxy2rt at each vertical level 
    for k in eachindex(z)
        field_rtz[:,:,k] = regrid_xy2rt(cx,cy,x,y,vardata[:,:,k])
    end
    # Return field_rtz, r, and theta if rt_out is true
    if rt_out == true
        return r,theta,field_rtz
    else
        return field_rtz
    end
end

#==============================================================================
regrid_gfrelxz

This "algorithm" interpolates a two-dimensional Cartesian squall line to gust 
front relative coordinates. It searches for the minimum in y-vorticity at the 
surface and then searches for the gridpoint nearest to -1 K potential 
temperature perturbation ahead of the minimum. The data at each vertical level 
are then interpolated to a new grid with respect to the surface gust front 
location.

** The algorithm was created for cm1 output of a quasi-linear squall line. If
   the system begins to bow, additional complexity may need to be added (e.g.,
   reconsider the x-buffer).

** This algorithm won't work if there's no thermodynamic data. Consider using 
   only y-vort in those cases.

** Input must be two dimensional [x,z] (i.e., specify a time stamp).

xsec_out - Option to output the x-dimension for the cross section
==============================================================================#

function regrid_gfrelxz(x::AbstractVector{<:Real},z::AbstractVector{<:Real},
                        yvort::AbstractArray{<:Real,2},thpert::AbstractArray{<:Real,2},
                        var::AbstractArray{<:Real,2},xsec_out::Bool=false) 

    # Start by finding the minimum in y-vorticity(η) at the surface
    minyvort = minimum(yvort[:,1])
    minvort_ind = findin(yvort[:,1],minyvort)[1]
    # Find the x-index where minyvort is located
    x_min = x[minvort_ind]
    # Search for the -1K thpert within the x buffer range
    indxgf = closest_ind(thpert[minvort_ind:end,1],-1.0)
    xgf = x[minvort_ind:end][indxgf]
    xgfnorm = x - xgf
    # Interpolate the var to the location of the gust front 
    # Define the universal x-grid for the interpolation
    gfrel_var = Array{Float64}(length(xgfnorm),length(z))
    for k in eachindex(z)
        fieldinterp = interpolate((xgfnorm,),var[:,k],Gridded(Linear()))
        for i in eachindex(x)
            gfrel_var[i,k] = fieldinterp[xgfnorm[i]]
        end
    end
    # Set the min and max range of the gfrel x-values (-50.0,20.0 for now)
    xgfmin = closest_ind(xgfnorm,-50.0)
    xgfmax = closest_ind(xgfnorm,20.0)
    # Define the x-section and var within the x-section
    xsec = xgfnorm[xgfmin:xgfmax]
    gfrel_var_xsec = gfrel_var[xgfmin:xgfmax,:]
    # Return the xsec only if xsec=true
    if xsec_out == true
        return xsec,gfrel_var_xsec
    else
        return gfrel_var_xsec
    end
end

#==============================================================================
regrid_gfrelxyz

This "algorithm" interpolates a three-dimensional Cartesian squall line to gust 
front relative coordinates. It searches for the minimum in y-vorticity at the 
surface and then searches for the gridpoint nearest to -1 K potential 
temperature perturbation ahead of the minimum. The data at each vertical level 
are then interpolated to a new grid with the location of the gust front at the 
center.

** The algorithm was created for cm1 output of a quasi-linear squall line. If
   the system begins to bow, additional complexity may need to be added (e.g.,
   reconsider the x-buffer).

** Input must be thee-dimensional [x,y,z] (i.e., specify a time stamp)

xsec_out - Option to output the x-dimension for the cross section
==============================================================================#

function regrid_gfrelxyz(x::AbstractVector{<:Real},y::AbstractVector{<:Real},
                         z::AbstractVector{<:Real},yvort::AbstractArray{<:Real,3},
                         thpert::AbstractArray{<:Real,3},var::AbstractArray{<:Real,3},
                         xsec_out::Bool=false)

    # Start by finding the minimum in y-vorticity(η) at the surface
    # Then find the x-index where minyvort is located at each y
    minyvort = similar(y,Float64)
    minvort_ind = Array{Int64}(length(y))
    x_min = similar(y)
    xbuff = Array{Float64}(21,length(y))
    fill!(x_min, NaN)
    fill!(minyvort, NaN)
    for j in eachindex(y)
        minyvort[j] = minimum(yvort[:,j,1])
        minvort_ind[j] = findin(yvort[:,j,1],minyvort[j])[1]
        x_min[j] = x[minvort_ind[j]]
        xbuff[:,j] = x[minvort_ind[j]-10:minvort_ind[j]+10]
    end
    # Search for the -1K thpert within the x buffer range
    xgf = similar(y,Float64)
    xgfnorm = Array{Float64}(length(x),length(y))
    fill!(xgf, NaN)
    fill!(xgfnorm, NaN)
    for j in eachindex(y)
        indxgf = closest_ind(thpert[minvort_ind[j]:end,j,1],-1.0)
        xgf[j] = x[minvort_ind[j]:end][indxgf]
        xgfnorm[:,j] = x-xgf[j]
    end
    # Interpolate the var to the location of the gust front 
    # Define the universal x-grid for the interpolation
    xunorm = collect(-170.0:129.0)
    gfrel_var = Array{Float64}(length(xunorm),length(y),length(z))
    for k in eachindex(z)
        for j in eachindex(y)
            fieldinterp = interpolate((xgfnorm[:,j],),var[:,j,k],Gridded(Linear()))
            for i in eachindex(x)
                gfrel_var[i,j,k] = fieldinterp[xunorm[i]]
            end
        end
    end
    # Set the min and max range of the gfrel x-values (-50.0,20.0 for now)
    xgfmin = closest_ind(xunorm,-50.0)
    xgfmax = closest_ind(xunorm,20.0)
    # Define the x-section and var within the x-section
    xsec = xunorm[xgfmin:xgfmax]
    gfrel_var_xsec = gfrel_var[xgfmin:xgfmax,:,:]
    # Return the xsec only if xsec=true
    if xsec_out == true
        return xsec,gfrel_var_xsec
    else
        return gfrel_var_xsec
    end
end

