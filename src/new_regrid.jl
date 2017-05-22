# *****************************************************************************
# new_regrid.jl
#
# Author: Jonathan Martinez 
# Email: jon.martinez@colostate.edu
# Julia version: 0.4.5
#
# ** Functions seem to be handling NaNs correctly with added 
# extrapolate functionality. Preliminary tests showed no differences using 
# new_regrid_xyz2rtz and regrid_xyz2rtz from regrid.jl. Thresholds no longer 
# required but are left commented out for reference.
#
# This script contains functions to interpolate data between Cartesian and 
# polar coordinates. See each individual function below for required data 
# input/output and formats.
#
# Function list:
# new_grid_2d
# new_grid_3d
# new_regrid_xy2rt 
# new_regrid_xyz2rtz
# regrid_gfrelxz
# regrid_gfrelxyz
#
# ** Note that the "new" prefixes are there just to avoid any conflicts when 
# loading this script and the original regrid.jl script in the same code
# *****************************************************************************

#==============================================================================
# new_grid_2d

This function is the equivalent of Matlab's ndgrid. It creates a two
dimensional grid from the one-dimensional input coordinate arrays.
** x and y do not have to be the same size.
==============================================================================#

function new_grid_2d{T<:Real}(x::Vector{T},y::Vector{T})
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
# new_grid_3d

3D formulation of new_grid_2d
** x, y, and z do not have to be the same size.
==============================================================================#

function new_grid_3d{T<:Real}(x::Vector{T},y::Vector{T},z::Vector{T})
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
# new_regrid_xy2rt

This function takes two-dimensional Cartesian data and interpolates it to a 
polar coordinate grid. It is specifically designed to work with NetCDF data
by working around issues with fill values when interpolating the data.
** Follows http://mathworld.wolfram.com/PolarCoordinates.html for converting 
   Cartesian to polar coordinates
 
rt_out - Option to output r and theta arrays. Default is false.
==============================================================================#

function new_regrid_xy2rt{T<:Real}(x::Vector{T},y::Vector{T},
                                   vardata::Array{T,2},rt_out::Bool=false)
    # Define lower/upper thresholds based on min/max values in orginal array
    #thresh_min = minimum(vardata[vardata .!= -999.0])
    #thresh_max = maximum(vardata[vardata .!= -999.0])
    # Determine the radius of the polar coordinate grid and the theta increment
    rmax = ceil(sqrt(((x[end]-x[1])/2.0)^2 + ((y[end]-y[1])/2.0)^2))
    r = collect(0:(x[2]-x[1]):rmax)
    theta_inc = floor(atan2(y[2]-y[1],(x[end]-x[1])/2.0)/pi*180.0)
    if theta_inc<1.0
        theta_inc=1.0
    end
    # Define theta
    theta = collect(0:theta_inc:360-theta_inc)
    # Create two-dimensional arrays for r and theta components of the polar grid
    r_2d,theta_2d = new_grid_2d(r,theta)
    # Define the Cartesian points in terms of the 2-d polar coordinates
    x_polar = r_2d .* cos(deg2rad(theta_2d))
    y_polar = r_2d .* sin(deg2rad(theta_2d))
    # Interpolate the data from the Cartesian grid to the polar grid 
    field_rt = Array(Float64,size(x_polar))
    field_interp = extrapolate(interpolate((x,y), vardata, Gridded(Linear())), NaN)
    for i in eachindex(r)
        for j in eachindex(theta)
             field_rt[i,j] = field_interp[x_polar[i,j],y_polar[i,j]]
        end
    end 
    # Now replace all interpolated values less/greater than 
    # thresh_min/thresh_max with NaN. This will effectively replace all
    # -999.0 values with NaN.
    #field_rt[field_rt .< thresh_min] = NaN
    #field_rt[field_rt .> thresh_max] = NaN
    # Return field_rt, r, and theta if rt_out is true
    if rt_out == true 
        return r,theta,field_rt
    else 
        return field_rt
    end
end

#==============================================================================
new_regrid_xyz2rtz

This function interpolates two-dimensional Cartesian data to a cylindrical 
coordinate grid by executing new_regrid_xy2rt at each vertical level.

rt_out - Option to output r and theta arrays. Default is false.
==============================================================================#

function new_regrid_xyz2rtz{T<:Real}(x::Vector{T},y::Vector{T},z::Vector{T},
                                     vardata::Array{T,3},rt_out::Bool=false)
    # Determine the radius of the polar coordinate grid and the theta increment
    rmax = ceil(sqrt(((x[end]-x[1])/2.0)^2 + ((y[end]-y[1])/2.0)^2))
    r = collect(0:(x[2]-x[1]):rmax)
    theta_inc = floor(atan2(y[2]-y[1],(x[end]-x[1])/2.0)/pi*180.0)
    if theta_inc<1
        theta_inc=1
    end
    # Define theta
    theta = collect(0:theta_inc:360-theta_inc)
   # Define the dimensions of the new var
    field_rtz = Array(Float64,length(r),length(theta),length(z))
    # Call regridxy2rt at each vertical level 
    for k in eachindex(z)
        field_rtz[:,:,k] = new_regrid_xy2rt(x,y,vardata[:,:,k])
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

** The algorithm uses the Interpolations pkg since there are no missing data.
   May have to adjust to CoordInterpGrid if missing data are present or 
   adjust the script to handle NaNs.

** Input must be two dimensional [x,z] (i.e., specify a time stamp).

xsec_out - Option to output the x-dimension for the cross section
==============================================================================#

function regrid_gfrelxz{T<:Real}(x::Vector{T},z::Vector{T},yvort::Array{T,2},
                                 thpert::Array{T,2},var::Array{T,2},
                                 xsec_out::Bool=false)
    # Start by finding the minimum in y-vorticity(η) at the surface
    minyvort = minimum(yvort[:,1])
    minvort_ind = findin(yvort[:,1],minyvort)[1]
    # Find the x-index where minyvort is located
    x_min = x[minvort_ind]
    # Define a function that will find the closest grid point to -1K
    function closest_index(arr,val)
        idiff = abs(arr-val);
        mindiff = minimum(idiff);
        ibest = findin(idiff,mindiff)[1]
       
        return ibest
    end
    # Search for the -1K thpert within the x buffer range
    indxgf = closest_index(thpert[minvort_ind:end,1],-1.0)
    xgf = x[minvort_ind:end][indxgf]
    xgfnorm = x - xgf
    # Interpolate the var to the location of the gust front 
    # Define the universal x-grid for the interpolation
    gfrel_var = Array(Float64,(length(xgfnorm),length(z)));
    for k in eachindex(z)
        fieldinterp = interpolate((xgfnorm,),var[:,k],Gridded(Linear()));
        for i in eachindex(x)
            gfrel_var[i,k] = fieldinterp[xgfnorm[i]]
        end    
    end   
    # Set the min and max range of the gfrel x-values (-50.0,20.0 for now)
    xgfmin = closest_index(xgfnorm,-50.0)
    xgfmax = closest_index(xgfnorm,20.0)
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

** The algorithm uses the Interpolations pkg since there are no missing data.
   May have to adjust to CoordInterpGrid if missing data are present or 
   adjust the script to handle NaNs.

** Input must be thee-dimensional [x,y,z] (i.e., specify a time stamp)

xsec_out - Option to output the x-dimension for the cross section
==============================================================================#

function regrid_gfrelxyz{T<:Real}(x::Vector{T},y::Vector{T},z::Vector{T},
                                  yvort::Array{T,3},thpert::Array{T,3},
                                  var::Array{T,3},xsec_out::Bool=false)
    # Start by finding the minimum in y-vorticity(η) at the surface
    # Then find the x-index where minyvort is located at each y
    minyvort = Array(Float64,length(y)) .* NaN
    minvort_ind = Array(Int64,length(y))
    x_min = Array(Float64,length(y)) .* NaN
    xbuff = Array(Float64,21,length(y)) 
    for j in eachindex(y)
        minyvort[j] = minimum(yvort[:,j,1])
        minvort_ind[j] = findin(yvort[:,j,1],minyvort[j])[1]
        x_min[j] = x[minvort_ind[j]]
        xbuff[:,j] = x[minvort_ind[j]-10:minvort_ind[j]+10]
    end 
    # Define a function that will find the closest grid point to -1K
    function closest_index(arr,val)
        idiff = abs(arr-val);
        mindiff = minimum(idiff);
        ibest = findin(idiff,mindiff)[1]
       
        return ibest
    end
    # Search for the -1K thpert within the x buffer range
    xgf = Array(Float64,length(y)) .* NaN
    xgfnorm = Array(Float64,length(x),length(y)) .* NaN    
    for j in eachindex(y)
        indxgf = closest_index(thpert[minvort_ind[j]:end,j,1],-1.0)
        xgf[j] = x[minvort_ind[j]:end][indxgf]
        xgfnorm[:,j] = x-xgf[j]
    end
    # Interpolate the var to the location of the gust front 
    # Define the universal x-grid for the interpolation
    xunorm = collect(-170.0:129.0);
    gfrel_var = Array(Float64,(length(xunorm),length(y),length(z)));
    for k in eachindex(z)
        for j in eachindex(y)
            fieldinterp = interpolate((xgfnorm[:,j],),var[:,j,k],Gridded(Linear()));
            for i in eachindex(x)
                gfrel_var[i,j,k] = fieldinterp[xunorm[i]]
            end    
        end   
    end 
    # Set the min and max range of the gfrel x-values (-50.0,20.0 for now)
    xgfmin = closest_index(xunorm,-50.0)
    xgfmax = closest_index(xunorm,20.0)
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


