# *****************************************************************************
# regrid.jl
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# This script contains functions to regrid data using various methods.
#
# Function list
# closest_ind
# grid2d
# grid3d
# xy2rp
# regrid_xy2rp
# regrid_xyz2rpz
# regrid_pol2cart
# unstagger
# curv2rect 
# curv2rect_wgts
# *****************************************************************************

#==============================================================================
# closest_ind

Search for the index within a 1-D array which most closely corresponds to the
specified value.
==============================================================================#

function closest_ind(arr::AbstractVector{T},val::Real) where T<:Real
    return argmin(abs.(arr .- val))
end

#==============================================================================
# grid2d

This function is the equivalent of Matlab's ndgrid. It creates a two
dimensional grid from the one-dimensional input coordinate arrays.
** x and y do not have to be the same size.
==============================================================================#

function grid2d(x::AbstractVector{Ta},y::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}
    return @fastmath x .* ones(length(y))', ones(length(x)) .* y'
end

#==============================================================================
# grid3d

3D formulation of grid2d
** x, y, and z do not have to be the same size.
==============================================================================#

function grid3d(x::AbstractVector{Ta},y::AbstractVector{Tb},
                z::AbstractVector{Tc}) where {Ta<:Real,Tb<:Real,Tc<:Real}

    # Create new arrays
    x_3d = Array{Float64}(undef,length(x),length(y),length(z))
    y_3d = Array{Float64}(undef,length(x),length(y),length(z))
    z_3d = Array{Float64}(undef,length(x),length(y),length(z))
    # 3d x-array
    for i in eachindex(x)
        x_3d[i,:,:] .= x[i]
    end
    # 3d y-array
    for j in eachindex(y)
        y_3d[:,j,:] .= y[j]
    end
    # 3d z-array
    for k in eachindex(z)
        z_3d[:,:,k] .= z[k]
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

function xy2rp(cx::Real,cy::Real,x::AbstractVector{Ta},
               y::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real,Tc<:Real}

    # Pinpoint the center of grid that has no stretching
    # Has no impact on grids with constant spacing
    x0 = closest_ind(x .- cx,0.0)
    y0 = closest_ind(y .- cy,0.0)
    # Determine the max radius of polar coordinate grid and the phi increment
    rmax = floor(sqrt((maximum(abs.(x)))^2 + (maximum(abs.(y))^2)))
    phi_inc = floor(atan(y[y0]-y[y0-1],maximum(abs.(x))))
    phi_inc < pi/180 ? phi_inc = pi/180 : nothing
    # Return radius and azimuth arrays
    r = collect(0:(round(x[x0],digits=2) - round(x[x0-1],digits=2)):rmax)
    phi = collect(0:phi_inc:2*pi - phi_inc)
    return r,phi
end

#==============================================================================
# regrid_xy2rp

This function takes two-dimensional Cartesian data and interpolates it to a
polar coordinate grid. It is specifically designed to work with NetCDF data
by working around issues with fill values when interpolating the data.
** Follows http://mathworld.wolfram.com/PolarCoordinates.html for converting
   Cartesian to polar coordinates
==============================================================================#

function regrid_xy2rp(cx::Real,cy::Real,x::AbstractVector{Ta},y::AbstractVector{Tb},
                      field::AbstractArray{Tc,2}) where {Ta<:Real,Tb<:Real,Tc<:Real}

    # Define r and phi
    r,phi = xy2rp(cx,cy,x,y)
    # Interpolate the data from the Cartesian grid to the polar grid
    field_rp = Array{Float64}(undef,length(r),length(phi))
    field_itp = extrapolate(interpolate((x.-cx, y.-cy),field,Gridded(Linear())),NaN)
    for j in eachindex(phi)
        for i in eachindex(r)
             @inbounds field_rp[i,j] = field_itp(r[i] * cos(phi[j]), r[i] * sin(phi[j]))
        end
    end
    return field_rp
end

# Allow r and phi as input arguments to the function

function regrid_xy2rp(cx::Real,cy::Real,x::AbstractVector{Ta},y::AbstractVector{Tb},
                      r::AbstractVector{Tc},phi::AbstractVector{Td},
                      field::AbstractArray{Te,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real}

    # Interpolate the data from the Cartesian grid to the polar grid
    field_rp = Array{Float64}(undef,length(r),length(phi))
    field_itp = extrapolate(interpolate((x.-cx, y.-cy),field,Gridded(Linear())),NaN)
    for j in eachindex(phi)
        for i in eachindex(r)
             @inbounds field_rp[i,j] = field_itp(r[i] * cos(phi[j]), r[i] * sin(phi[j]))
        end
    end
    return field_rp
end

# If Cartesian domain is already centered at (0,0), no need for cx and cy input

function regrid_xy2rp(x::AbstractVector{Ta},y::AbstractVector{Tb},
                      r::AbstractVector{Tc},phi::AbstractVector{Td},
                      field::AbstractArray{Te,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real}

    # Interpolate the data from the Cartesian grid to the polar grid
    field_rp = Array{Float64}(undef,length(r),length(phi))
    field_itp = extrapolate(interpolate((x,y),field,Gridded(Linear())),NaN)
    for j in eachindex(phi)
        for i in eachindex(r)
             @inbounds field_rp[i,j] = field_itp(r[i] * cos(phi[j]), r[i] * sin(phi[j]))
        end
    end
    return field_rp
end

#==============================================================================
regrid_xyz2rpz

This function interpolates three-dimensional Cartesian data to a cylindrical
coordinate grid by executing regrid_xy2rp at each vertical level.
==============================================================================#

function regrid_xyz2rpz(cx::Real,cy::Real,x::AbstractVector{Ta},
                        y::AbstractVector{Tb},z::AbstractVector{Tc},
                        field::AbstractArray{Td,3}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}

    # Define r and phi
    r,phi = xy2rp(cx,cy,x,y)
    # Define the dimensions of the new var
    field_rpz = Array{Float64}(undef,length(r),length(phi),length(z))
    # Call regrid_xy2rp at each vertical level
    for k in eachindex(z)
        @inbounds field_rpz[:,:,k] = regrid_xy2rp(cx,cy,x,y,field[:,:,k])
    end
    # Return field_rpz
    return field_rpz
end

# Allow r and phi as input arguments to the function

function regrid_xyz2rpz(cx::Real,cy::Real,x::AbstractVector{Ta},
                        y::AbstractVector{Tb},z::AbstractVector{Tc},
                        r::AbstractVector{Td},phi::AbstractVector{Te},
                        field::AbstractArray{Tf,3}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real,Tf<:Real}

    # Define the dimensions of the new var
    field_rpz = Array{Float64}(undef,length(r),length(phi),length(z))
    # Call regrid_xy2rp at each vertical level
    for k in eachindex(z)
        @inbounds field_rpz[:,:,k] = regrid_xy2rp(cx,cy,x,y,r,phi,field[:,:,k])
    end
    # Return field_rpz
    return field_rpz
end

# If Cartesian domain is already centered at (0,0), no need for cx and cy input

function regrid_xyz2rpz(x::AbstractVector{Ta},y::AbstractVector{Tb},z::AbstractVector{Tc},
                        r::AbstractVector{Td},phi::AbstractVector{Te},
                        field::AbstractArray{Tf,3}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real,Tf<:Real}

    # Define the dimensions of the new var
    field_rpz = Array{Float64}(undef,length(r),length(phi),length(z))
    # Call regrid_xy2rp at each vertical level
    for k in eachindex(z)
        @inbounds field_rpz[:,:,k] = regrid_xy2rp(x,y,r,phi,field[:,:,k])
    end
    # Return field_rpz
    return field_rpz
end

#==============================================================================
regrid_pol2cart

This function will regrid data on a polar grid to Cartesian coordinates
==============================================================================#

# Interpolate from radius to (x,y) -- axisymmetric application

function regrid_pol2cart(cx::Real,cy::Real,r::AbstractVector{Ta},
                         field::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}


    # Specify the x and y arrays based on r
    x = collect(-r[end]:r[2]-r[1]:r[end])
    y = collect(-r[end]:r[2]-r[1]:r[end])
    # Create the interpolation object
    field_xy = Array{Float64}(undef,length(x),length(y))
    field_itp = extrapolate(interpolate((r,),field,Gridded(Linear())),NaN)
    # Interpolate from axisymmetric polar to Cartesian
    for j in eachindex(y)
        for i in eachindex(x)
            @inbounds field_xy[i,j] = field_itp(sqrt((x[i]-cx)^2 + (y[j]-cy)^2))
        end
    end
    return field_xy
end

# Allow x and y as input arguments to the axisymmetric function

function regrid_pol2cart(cx::Real,cy::Real,r::AbstractVector{Ta},
                         x::AbstractVector{Tb},y::AbstractVector{Tc},
                         field::AbstractVector{Td}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}

    # Create the interpolation object
    field_xy = Array{Float64}(undef,length(x),length(y))
    field_itp = extrapolate(interpolate((r,),field,Gridded(Linear())),NaN)
    # Interpolate from axisymmetric polar to Cartesian
    for j in eachindex(y)
        for i in eachindex(x)
            @inbounds field_xy[i,j] = field_itp(sqrt((x[i]-cx)^2 + (y[j]-cy)^2))
        end
    end
    return field_xy
end

# Interpolate from (r,phi) to (x,y)
# ** Assumes data on polar grid are centered at (0,0)

function regrid_pol2cart(r::AbstractVector{Ta},phi::AbstractVector{Tb},
                         field::AbstractArray{Tc,2}) where {Ta<:Real,Tb<:Real,Tc<:Real}

    x = collect(-r[end]:r[2]-r[1]:r[end])
    y = collect(-r[end]:r[2]-r[1]:r[end])
    # Create 2-D grids
    xx,yy = grid2d(x,y)
    # Transform Cartesian to polar reference points for azimuth
    # due to atan2 returning negative values
    pp_cart = atan.(yy,xx)
    pp_cart[pp_cart .< 0] = pp_cart[pp_cart .< 0] .+ 2*pi
    # Create the interpolation object
    field_xy = Array{Float64}(undef,length(x),length(y))
    field_itp = extrapolate(interpolate((r,phi),field,Gridded(Linear())),NaN)
    # Interpolate from polar to Cartesian
    for j in eachindex(y)
        for i in eachindex(x)
            @inbounds field_xy[i,j] = field_itp(sqrt(x[i]^2 + y[j]^2), pp_cart[i,j])
        end
    end
    return field_xy
end

# Allow x and y as input (still assumes data on polar grid are centered at (0,0))

function regrid_pol2cart(r::AbstractVector{Ta},phi::AbstractVector{Tb},
                         x::AbstractVector{Tc},y::AbstractVector{Td},
                         field::AbstractArray{Te,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real}

    # Create 2-D grids
    xx,yy = grid2d(x,y)
    # Transform Cartesian to polar reference points for azimuth
    # due to atan2 returning negative values
    pp_cart = atan.(yy,xx)
    pp_cart[pp_cart .< 0] = pp_cart[pp_cart .< 0] .+ 2*pi
    # Create the interpolation object
    field_xy = Array{Float64}(undef,length(x),length(y))
    field_itp = extrapolate(interpolate((r,phi),field,Gridded(Linear())),NaN)
    # Interpolate from polar to Cartesian
    for j in eachindex(y)
        for i in eachindex(x)
            @inbounds field_xy[i,j] = field_itp(sqrt(x[i]^2 + y[j]^2), pp_cart[i,j])
        end
    end
    return field_xy
end

# Interpolate from (r,phi,z) to (x,y,z)
# ** Assumes data at each vertical level are centered at (0,0)

function regrid_pol2cart(r::AbstractVector{Ta},phi::AbstractVector{Tb},
                         x::AbstractVector{Tc},y::AbstractVector{Td},
                         z::AbstractVector{Te},field::AbstractArray{Tf,3}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real,Tf<:Real}

    # Define the dimensions of the new var
    field_xyz = Array{Float64}(undef,length(x),length(y),length(z))
    # Interpolate from polar to Cartesian
    for k in eachindex(z)
        @inbounds field_xyz[:,:,k] = regrid_pol2cart(r,phi,x,y,field[:,:,k])
    end
    return field_xyz
end

#==============================================================================
unstagger

This function will take u, v, or w (vectors) on a 3-D grid and unstagger them,
placing them on the same grid as scalar variables
==============================================================================#

function unstagger(grid::AbstractArray{<:Real},method::Symbol=:u)
    if method == :u
        gridout = 0.5 * (grid[1:end-1,:,:] + grid[2:end,:,:])
        return gridout
    elseif method == :v
        gridout = 0.5 * (grid[:,1:end-1,:] + grid[:,2:end,:])
        return gridout
    elseif method == :w
        gridout = 0.5 * (grid[:,:,1:end-1] + grid[:,:,2:end])
        return gridout
    end
end

#==============================================================================
curv2rect

Interpolate a 2-D curvilinear (Cartesian) grid to a rectilinear (Cartesian)
grid using inverse-distance-squared weighting
- Note - Currently don't have any constraints to prevent function from going 
         out of bounds --> Make the (x,y) coordinate arrays at least one index
         smaller than corresponding extrema for (x2d,y2d)

(x2d, y2d) = 2-D curvilinear coordinate arrays for x and y dimensions
(x, y) = User-specified 1-D rectilinear coordinate arrays for x and y dimensions
field = 2-D variable mapped on the curvilinear grid

Output
field(x,y) = 2-D variable interpolated to the specified rectilinear grid
==============================================================================#

function curv2rect(x2d::AbstractArray{Ta,2},y2d::AbstractArray{Tb,2},
                   x::AbstractVector{Tc},y::AbstractVector{Td},
                   field::AbstractArray{Te,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real}

    field_out = zeros(length(x),length(y))
    for j in eachindex(y)
        for i in eachindex(x)
            # Find closest point on the curvilinear grid
            r = sqrt.((x[i] .- x2d).^2 + (y[j] .- y2d).^2)
            rmin_ind = Tuple(argmin(r))
            # If r = 0, then no weighting is necessary
            if r[rmin_ind...] == 0.
                field_out[i,j] = field[i,j]
            else 
                xmin_ind = rmin_ind[1]
                ymin_ind = rmin_ind[2]
                # Define indices for the neighboring points in each cardinal direction
                east_ind = (xmin_ind+1, ymin_ind)
                west_ind = (xmin_ind-1, ymin_ind)
                north_ind = (xmin_ind, ymin_ind+1)
                south_ind = (xmin_ind, ymin_ind-1)
                # Pull center point and neighboring points from each cardinal direction
                field_cen = field[rmin_ind...]
                field_east = field[east_ind...]
                field_west = field[west_ind...]
                field_north = field[north_ind...]
                field_south = field[south_ind...]
                # Calculate the radius weights using inverse distance squared
                inv_rsq_cen = inv(sqrt( (x[i] - x2d[rmin_ind...])^2 + (y[j] - y2d[rmin_ind...])^2)^2 )
                inv_rsq_east = inv(sqrt( (x[i] - x2d[east_ind...])^2 + (y[j] - y2d[east_ind...])^2)^2 )
                inv_rsq_west = inv(sqrt( (x[i] - x2d[west_ind...])^2 + (y[j] - y2d[west_ind...])^2)^2 )
                inv_rsq_north = inv(sqrt( (x[i] - x2d[north_ind...])^2 + (y[j] - y2d[north_ind...])^2)^2 )
                inv_rsq_south = inv(sqrt( (x[i] - x2d[south_ind...])^2 + (y[j] - y2d[south_ind...])^2)^2 )
                # Calculate the weighted average for the given rectilinear grid point
                field_num = (field_cen * inv_rsq_cen) + (field_east * inv_rsq_east) + (field_west * inv_rsq_west) +
                            (field_north * inv_rsq_north) + (field_south * inv_rsq_south)
                field_denom = inv_rsq_cen + inv_rsq_east + inv_rsq_west + inv_rsq_north + inv_rsq_south
                # Interpolated variable
                field_out[i,j] = field_num/field_denom
            end # r = 0 condition
        end # x[i]
    end # y[j]
    return field_out
end

#==============================================================================
curv2rect_wgts

Create the stencil and weights for an n x n inverse-distance-squared weighted
average
- Requires 2-D curvilinear grid points and desired rectilinear grid points
  to be in Cartesian coordinates

- Note - Currently don't have any constraints to prevent function from going 
         out of bounds --> Make the desired (x,y) coordinate arrays at least 
         one index smaller than corresponding extrema for (x2d,y2d)

Input
(x2d, y2d) = 2-D curvilinear coordinate arrays for x and y dimensions
(x, y) = User-specified 1-D rectilinear coordinate arrays for x and y dimensions
stencil = Vector containing indices to construct the square weighting stencil
- Note - Currently limited to square stencils (n x n)

Output
wgt_inds = Array containing indices for each square stencil created for each
           (x,y) grid point
wgt_vals = Array containing the inverse-distance-squared weight values for each
           index on the square stencils
==============================================================================#

function curv2rect_wgts(x2d::AbstractArray{Ta,2},y2d::AbstractArray{Tb,2},
                        x::AbstractVector{Tc},y::AbstractVector{Td};
                        stencil::AbstractVector{Te}=[-1,0,1]) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real}

    # Construct n x n wgt arrays for n=length(stencil)
    n = length(stencil)
    wgt_inds = zeros(CartesianIndex{2},n,n,length(x),length(y))
    wgt_vals = zeros(n,n,length(x),length(y))
    for j in eachindex(y)
        for i in eachindex(x)
            # Find closest point on the curvilinear grid
            r = sqrt.((x[i] .- x2d).^2 + (y[j] .- y2d).^2)
            rmin_ind = argmin(r)
            for jj in 1:n
                for ii in 1:n
                    wgt_inds[ii,jj,i,j] = CartesianIndex(rmin_ind[1] + stencil[ii], rmin_ind[2] + stencil[jj])
                end 
            end 
            # r = 0 condition
            if r[rmin_ind] == 0.
                mid_ind = div(n,2)+1
                wgt_vals[mid_ind,mid_ind,i,j] += 1. 
            else 
                wgt_vals[:,:,i,j] .+= inv.(sqrt.( (x[i] .- x2d[wgt_inds[:,:,i,j]]).^2 .+ (y[j] .- y2d[wgt_inds[:,:,i,j]]).^2).^2 )
            end
        end # x[i]
    end # y[j]
    return wgt_inds, wgt_vals
end
