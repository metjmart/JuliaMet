# *****************************************************************************
# derivative4.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia Version: 1.0.0
#
# This script contains several functions used to compute derivatives with
# finite differencing using fourth-order accurate (5-point) stencils.
# The units of the input coordinate arrays will determine the units of the
# differentiated field. E.g., for a coordinate array in meters (m), the
# differentiated field of temperature (K) will have units of K/m
#
# Function list
# fd_weights4
# finite_ds4
# finite_dsx4
# finite_dsy4
# finite_dsz4
# *****************************************************************************

#==============================================================================
fd_weights4

Determine the fourth-order finite difference weights for either uniform or
non-uniform grids by solving a linear system of equations
Note: This is really sloppy. Need to consider using composite types to pass
around the stencils and matrices
See Eq. (2.53) in Dave Randall's textbook: An Introduction to Numerical 
Modeling of the Atmosphere
==============================================================================#

# Fourth-order accurate first derivative stencil

function fd_weights4(x::AbstractVector{Ta};order::Int=4) where Ta<:Real

    # Create the stencil given the order of accuracy
    stencil = zeros(order+1)
    # Array b constructed from Kronecker delta given 0 <= m <= n + l - 1
    # m = l ? delta = 1 : delta = 0 (where n = order; l = lth derivative)
    b = [0.,1.,0.,0.,0.]
    wgts = Array{Float64}(undef,length(stencil),length(x))
    for j in eachindex(x)
        if j == 1
            stencil = [0,1,2,3,4]
        elseif j == 2
            stencil = [-1,0,1,2,3]
        elseif j == length(x)
            stencil = [-4,-3,-2,-1,0]
        elseif j == length(x)-1
            stencil = [-3,-2,-1,0,1]
        else
            stencil = [-2,-1,0,1,2]
        end
        # Matrix must be constructed for each grid point
        a = ones(order+1,order+1)
        for k in eachindex(stencil)
            for m in 0:order
                a[m+1,k] = (x[j+stencil[k]] - x[j])^m
            end
        end
        wgts[:,j] = a\b
    end
    return wgts
end

#==============================================================================
finite_ds4

General one-dimensional application
Fourth-order accurate first derivative finite difference scheme for uniform
or non-uniform grids
==============================================================================#

function finite_ds4(wgts::AbstractArray{Ta,2},field::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}

    # Create the stencil given the order of accuracy
    stencil = zeros(size(wgts)[1])
    dvards = similar(field,Float64)
    for i in eachindex(field)
        if i == 1
            stencil = [i,i+1,i+2,i+3,i+4]
        elseif i == 2
            stencil = [i-1,i,i+1,i+2,i+3]
        elseif i == length(field)
            stencil = [i-4,i-3,i-2,i-1,i]
        elseif i == length(field)-1
            stencil = [i-3,i-2,i-1,i,i+1]
        else
            stencil = [i-2,i-1,i,i+1,i+2]
        end
        dvards[i] = sum(wgts[:,i] .* field[stencil])
    end
    return dvards
end

#==============================================================================
finite_dsx4

Fourth-order accurate first derivative finite difference scheme for uniform
or non-uniform grids
Apply to the x-dimension (first dimension) of an xyz field
==============================================================================#

# x-dimension

function finite_dsx4(wgts::AbstractArray{Ta,2},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

    # Two-dimensional variable
    if ndims(field) == 2
        if size(wgts)[2] != size(field)[1]
            error("There should be one set of weights for each grid point along
                   the dimension to be differentiated")
        end
        # Create the stencil given the order of accuracy
        stencil = zeros(size(wgts)[1])
        dvards = similar(field,Float64)
        for i in 1:size(field)[1]
            if i == 1
                stencil = [i,i+1,i+2,i+3,i+4]
            elseif i == 2
                stencil = [i-1,i,i+1,i+2,i+3]
            elseif i == size(field)[1]
                stencil = [i-4,i-3,i-2,i-1,i]
            elseif i == size(field)[1]-1
                stencil = [i-3,i-2,i-1,i,i+1]
            else
                stencil = [i-2,i-1,i,i+1,i+2]
            end
            for j in 1:size(field)[2]
                dvards[i,j] = sum(wgts[:,i] .* field[stencil,j])
            end
        end
        return dvards
    # Three-dimensional variable
    elseif ndims(field) == 3
        if size(wgts)[2] != size(field)[1]
            error("There should be one set of weights for each grid point along
                   the dimension to be differentiated")
        end
        # Create the stencil given the order of accuracy
        stencil = zeros(size(wgts)[1])
        dvards = similar(field,Float64)
        for i in 1:size(field)[1]
            if i == 1
                stencil = [i,i+1,i+2,i+3,i+4]
            elseif i == 2
                stencil = [i-1,i,i+1,i+2,i+3]
            elseif i == size(field)[1]
                stencil = [i-4,i-3,i-2,i-1,i]
            elseif i == size(field)[1]-1
                stencil = [i-3,i-2,i-1,i,i+1]
            else
                stencil = [i-2,i-1,i,i+1,i+2]
            end
            for k in 1:size(field)[3]
                for j in 1:size(field)[2]
                    dvards[i,j,k] = sum(wgts[:,i] .* field[stencil,j,k])
                end
            end
        end
        return dvards
    else
        error("Input variable to be differentiated cannot exceed 3 dimensions")
    end
end

#==============================================================================
finite_dsy4

Fourth-order accurate first derivative finite difference scheme for uniform
or non-uniform grids
Apply to the y-dimension (second dimension) of an xyz field
==============================================================================#

# y-dimension

function finite_dsy4(wgts::AbstractArray{Ta,2},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

    # Two-dimensional variable
    if ndims(field) == 2
        if size(wgts)[2] != size(field)[2]
            error("There should be one set of weights for each grid point along
                   the dimension to be differentiated")
        end
        # Create the stencil given the order of accuracy
        stencil = zeros(size(wgts)[1])
        dvards = similar(field,Float64)
        for j in 1:size(field)[2]
            if j == 1
                stencil = [j,j+1,j+2,j+3,j+4]
            elseif j == 2
                stencil = [j-1,j,j+1,j+2,j+3]
            elseif j == size(field)[2]
                stencil = [j-4,j-3,j-2,j-1,j]
            elseif j == size(field)[2]-1
                stencil = [j-3,j-2,j-1,j,j+1]
            else
                stencil = [j-2,j-1,j,j+1,j+2]
            end
            for i in 1:size(field)[1]
                dvards[i,j] = sum(wgts[:,j] .* field[i,stencil])
            end
        end
        return dvards
    # Three-dimensional variable
    elseif ndims(field) == 3
        if size(wgts)[2] != size(field)[2]
            error("There should be one set of weights for each grid point along
                   the dimension to be differentiated")
        end
        # Create the stencil given the order of accuracy
        stencil = zeros(size(wgts)[1])
        dvards = similar(field,Float64)
        for j in 1:size(field)[2]
            if j == 1
                stencil = [j,j+1,j+2,j+3,j+4]
            elseif j == 2
                stencil = [j-1,j,j+1,j+2,j+3]
            elseif j == size(field)[2]
                stencil = [j-4,j-3,j-2,j-1,j]
            elseif j == size(field)[2]-1
                stencil = [j-3,j-2,j-1,j,j+1]
            else
                stencil = [j-2,j-1,j,j+1,j+2]
            end
            for k in 1:size(field)[3]
                for i in 1:size(field)[1]
                    dvards[i,j,k] = sum(wgts[:,j] .* field[i,stencil,k])
                end
            end
        end
        return dvards
    else
        error("Input variable to be differentiated cannot exceed 3 dimensions")
    end
end

#==============================================================================
finite_dsz4

Fourth-order accurate first derivative finite difference scheme for uniform
or non-uniform grids
Apply to the z-dimension (third dimension) of an xyz field
==============================================================================#

# z-dimension

function finite_dsz4(wgts::AbstractArray{Ta,2},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

    # Two-dimensional variable
    if ndims(field) == 2
        if size(wgts)[2] != size(field)[2]
            error("There should be one set of weights for each grid point along
                   the dimension to be differentiated")
        end
        # Create the stencil given the order of accuracy
        stencil = zeros(size(wgts)[1])
        dvards = similar(field,Float64)
        for k in 1:size(field)[2]
            if k == 1
                stencil = [k,k+1,k+2,k+3,k+4]
            elseif k == 2
                stencil = [k-1,k,k+1,k+2,k+3]
            elseif k == size(field)[2]
                stencil = [k-4,k-3,k-2,k-1,k]
            elseif k == size(field)[2]-1
                stencil = [k-3,k-2,k-1,k,k+1]
            else
                stencil = [k-2,k-1,k,k+1,k+2]
            end
            for j in 1:size(field)[1]
                dvards[j,k] = sum(wgts[:,k] .* field[j,stencil])
            end
        end
        return dvards
    # Three-dimensional variable
    elseif ndims(field) == 3
        if size(wgts)[2] != size(field)[3]
            error("There should be one set of weights for each grid point along
                   the dimension to be differentiated")
        end
        # Create the stencil given the order of accuracy
        stencil = zeros(size(wgts)[1])
        dvards = similar(field,Float64)
        for k in 1:size(field)[3]
            if k == 1
                stencil = [k,k+1,k+2,k+3,k+4]
            elseif k == 2
                stencil = [k-1,k,k+1,k+2,k+3]
            elseif k == size(field)[3]
                stencil = [k-4,k-3,k-2,k-1,k]
            elseif k == size(field)[3]-1
                stencil = [k-3,k-2,k-1,k,k+1]
            else
                stencil = [k-2,k-1,k,k+1,k+2]
            end
            for j in 1:size(field)[2]
                for i in 1:size(field)[1]
                    dvards[i,j,k] = sum(wgts[:,k] .* field[i,j,stencil])
                end
            end
        end
        return dvards
    else
        error("Input variable to be differentiated cannot exceed 3 dimensions")
    end
end




