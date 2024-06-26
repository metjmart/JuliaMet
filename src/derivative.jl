# *****************************************************************************
# derivative.jl
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# This script contains several functions used to compute finite difference 
# approximations to derivatives
#
# The units of the input coordinate arrays will determine the units of the
# differentiated field. E.g., for a coordinate array in meters (m), the
# differentiated field of temperature (K) will have units of K/m
#
# Function list
# finite_dr
# finite_dx
# finite_dy
# finite_dz
# finite_laplacian
# kron_delta
# fd_weights
# finite_ds
# finite_dsx
# finite_dsy
# finite_dsz
# diff1
# spdiff1
# Laplacian
# *****************************************************************************

#==============================================================================
finite_dr

Calculate the radial derivative of a specified rpz field.
*** Assumes the radial dimension is the first dimension in the input array!
==============================================================================#

function finite_dr(r::AbstractVector{Ta},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

    if ndims(field) == 1
        if length(r) != length(field)
            error("Vector r and input variable to be differentiated must be of
                   same length")
        end
        dvardr = similar(field,Float64)
        # Loop to compute radial derivative
        for i in eachindex(r)
            # Forward difference at innermost boundary
            if i==1
                dvardr[i] = (-3.0*field[i] + 4.0*field[i+1] - field[i+2]) /
                            (r[i+2]-r[i])
            # Reverse difference at outermost boundary
            elseif i==length(r)
                dvardr[i] = (3.0*field[i] - 4.0*field[i-1] + field[i-2]) /
                            (r[i]-r[i-2])
            # Centered difference at all other grid points
            else
                dvardr[i] = (field[i+1]-field[i-1]) / (r[i+1]-r[i-1])
            end
        end

        return dvardr
    elseif ndims(field) == 2
        if length(r) != size(field)[1]
            error("Vector r and first dimension of variable to be differentiated
                   must be of same length")
        end
        dvardr = similar(field,Float64)
        # Loop over all dimensions to compute the radial derivative
        for j in 1:size(field)[2]
            for i in eachindex(r)
                # Forward difference at innermost boundary
                if i==1
                    dvardr[i,j] = (-3.0*field[i,j] + 4.0*field[i+1,j] -
                                   field[i+2,j]) / (r[i+2]-r[i])
                # Reverse difference at outermost boundary
                elseif i==length(r)
                    dvardr[i,j] = (3.0*field[i,j] - 4.0*field[i-1,j] +
                                   field[i-2,j]) / (r[i]-r[i-2])
                # Centered difference at all other grid points
                else
                    dvardr[i,j] = (field[i+1,j]-field[i-1,j]) / (r[i+1]-r[i-1])
                end
            end
        end

        return dvardr
    elseif ndims(field) == 3
        if length(r) != size(field)[1]
            error("Vector r and first dimension of variable to be differentiated
                   must be of same length")
        end
        dvardr = similar(field,Float64)
        # Loop over all dimensions to compute the radial derivative
        for k in 1:size(field)[3]
            for j in 1:size(field)[2]
                for i in eachindex(r)
                    # Forward difference at innermost boundary
                    if i==1
                        dvardr[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i+1,j,k] -
                                         field[i+2,j,k]) / (r[i+2]-r[i])
                    # Reverse difference at outermost boundary
                    elseif i==length(r)
                        dvardr[i,j,k] = (3.0*field[i,j,k] - 4.0*field[i-1,j,k] +
                                         field[i-2,j,k]) / (r[i]-r[i-2])
                    # Centered difference at all other grid points
                    else
                        dvardr[i,j,k] = (field[i+1,j,k]-field[i-1,j,k]) /
                                        (r[i+1]-r[i-1])
                    end
                end
            end
        end

        return dvardr
    else
        error("Input variable to be differentiated cannot exceed 3 dimensions")
    end
end

#==============================================================================
finite_dx

Calculate the x-horizontal derivative of a specified xyz field.
*** Assumes the x-dimension is the first dimension in the input array!
==============================================================================#

function finite_dx(x::AbstractVector{Ta},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

    if ndims(field) == 1
        if length(x) != length(field)
            error("Vector x and input variable to be differentiated must be of
                   same length")
        end
        dvardx = similar(field,Float64)
        # Loop to compute x-derivative
        for i in eachindex(x)
            # Forward difference at innermost boundary
            if i==1
                dvardx[i] = (-3.0*field[i] + 4.0*field[i+1] - field[i+2]) /
                            (x[i+2]-x[i])
            # Reverse difference at outermost boundary
            elseif i==length(x)
                dvardx[i] = (3.0*field[i] - 4.0*field[i-1] + field[i-2]) /
                            (x[i]-x[i-2])
            # Centered difference at all other grid points
            else
                dvardx[i] = (field[i+1]-field[i-1]) / (x[i+1]-x[i-1])
            end
        end

        return dvardx
    elseif ndims(field) == 2
        if length(x) != size(field)[1]
            error("Vector x and first dimension of variable to be differentiated
                   must be of same length")
        end
        dvardx = similar(field,Float64)
        # Loop over all dimensions to compute the x-derivative
        for j in 1:size(field)[2]
            for i in eachindex(x)
                # Forward difference at innermost boundary
                if i==1
                    dvardx[i,j] = (-3.0*field[i,j] + 4.0*field[i+1,j] -
                                   field[i+2,j]) / (x[i+2]-x[i])
                # Reverse difference at outermost boundary
                elseif i==length(x)
                    dvardx[i,j] = (3.0*field[i,j] - 4.0*field[i-1,j] +
                                   field[i-2,j]) / (x[i]-x[i-2])
                # Centered difference at all other grid points
                else
                    dvardx[i,j] = (field[i+1,j]-field[i-1,j]) / (x[i+1]-x[i-1])
                end
            end
        end

        return dvardx
    elseif ndims(field) == 3
        if length(x) != size(field)[1]
            error("Vector x and first dimension of variable to be differentiated
                   must be of same length")
        end
        dvardx = similar(field,Float64)
        # Loop over all dimensions to compute the x-derivative
        for k in 1:size(field)[3]
            for j in 1:size(field)[2]
                for i in eachindex(x)
                    # Forward difference at innermost boundary
                    if i==1
                        dvardx[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i+1,j,k] -
                                         field[i+2,j,k]) / (x[i+2]-x[i])
                    # Reverse difference at outermost boundary
                    elseif i==length(x)
                        dvardx[i,j,k] = (3.0*field[i,j,k] - 4.0*field[i-1,j,k] +
                                         field[i-2,j,k]) / (x[i]-x[i-2])
                    # Centered difference at all other grid points
                    else
                        dvardx[i,j,k] = (field[i+1,j,k]-field[i-1,j,k]) /
                                        (x[i+1]-x[i-1])
                    end
                end
            end
        end

        return dvardx
    else
        error("Input variable to be differentiated cannot exceed 3 dimensions")
    end
end

#==============================================================================
finite_dy

Calculate the y-horizontal derivative of a specified xyz field.
*** Assumes the y-dimension is the first dimension in the input array!
==============================================================================#

function finite_dy(y::AbstractVector{Ta},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

    if ndims(field) == 1
        if length(y) != length(field)
            error("Vector y and input variable to be differentiated must be of
                   same length")
        end
        dvardy = similar(field,Float64)
        # Loop over all dimensions to compute the y-derivative
        for i in eachindex(y)
            # Forward difference at innermost boundary
            if j==1
                dvardy[j] = (-3.0*field[j] + 4.0*field[j+1] - field[j+2]) /
                                (y[j+2]-y[j])
            # Reverse difference at outermost boundary
            elseif j==length(y)
                dvardy[j] = (3.0*field[j] - 4.0*field[j-1] + field[j-2]) /
                                (y[j]-y[j-2])
            # Centered difference at all other grid points
            else
                dvardy[j] = (field[j+1]-field[j-1]) / (y[j+1]-y[j-1])
            end
        end

        return dvardy
    elseif ndims(field) == 2
        if length(y) != size(field)[2]
            error("Vector y and second dimension of the variable to differentiated
                   must be of same length")
        end
        dvardy = similar(field,Float64)
        # Loop over all dimensions to compute the y-derivative
        for j in eachindex(y)
            for i in 1:size(field)[1]
                # Forward difference at innermost boundary
                if j==1
                    dvardy[i,j] = (-3.0*field[i,j] + 4.0*field[i,j+1] -
                                   field[i,j+2]) / (y[j+2]-y[j])
                # Reverse difference at outermost boundary
                elseif j==length(y)
                    dvardy[i,j] = (3.0*field[i,j] - 4.0*field[i,j-1] +
                                   field[i,j-2]) / (y[j]-y[j-2])
                # Centered difference at all other grid points
                else
                    dvardy[i,j] = (field[i,j+1] - field[i,j-1]) /
                                  (y[j+1]-y[j-1])
                end
            end
        end

        return dvardy
    elseif ndims(field) == 3
        if length(y) != size(field)[2]
            error("Vector y and second dimension of the variable to differentiated
                   must be of same length")
        end
        dvardy = similar(field,Float64)
        # Loop over all dimensions to compute the y-derivative
        for k in 1:size(field)[3]
            for j in eachindex(y)
                for i in 1:size(field)[1]
                    # Forward difference at innermost boundary
                    if j==1
                        dvardy[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i,j+1,k] -
                                         field[i,j+2,k]) / (y[j+2]-y[j])
                    # Reverse difference at outermost boundary
                    elseif j==length(y)
                        dvardy[i,j,k] = (3.0*field[i,j,k] - 4.0*field[i,j-1,k] +
                                         field[i,j-2,k]) / (y[j]-y[j-2])
                    # Centered difference at all other grid points
                    else
                        dvardy[i,j,k] = (field[i,j+1,k] - field[i,j-1,k]) /
                                        (y[j+1]-y[j-1])
                    end
                end
            end
        end

        return dvardy
    else
        error("Input variable to be differentiated cannot exceed 3 dimensions")
    end
end

#==============================================================================
finite_dz

Calculate the vertical derivative of a specified rpz or xyz field.
*** Assumes the vertical dimension is the third dimension in the input array!
==============================================================================#

function finite_dz(z::AbstractVector{Ta},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

    if ndims(field) == 1
        if length(z) != length(field)
            error("Vector z and input variable to be differentiated must be of
                   same length")
        end
        dvardz = similar(field,Float64)
        # Loop over all dimensions to compute the vertical derivative
        for k in eachindex(z)
            # Forward difference at lower boundary
            if k==1
                dvardz[k] = (-3.0*field[k] + 4.0*field[k+1] - field[k+2]) /
                            (z[k+2]-z[k])
            # Reverse difference at upper boundary
            elseif k==length(z)
                dvardz[k] = (3.0*field[k] - 4.0*field[k-1] + field[k-2]) /
                            (z[k]-z[k-2])
            # Centered difference at all other grid points
            else
                dvardz[k] = (field[k+1]-field[k-1]) / (z[k+1]-z[k-1])
            end
        end

        return dvardz
    elseif ndims(field) == 2
        if length(z) != size(field)[2]
            error("Vector z and second dimension of variable to be differentiated
                   must be of same length")
        end
        dvardz = similar(field,Float64)
        # Loop over all dimensions to compute the vertical derivative
        for k in eachindex(z)
            for j in 1:size(field)[1]
                # Forward difference at lower boundary
                if k==1
                    dvardz[j,k] = (-3.0*field[j,k] + 4.0*field[j,k+1] -
                                   field[j,k+2]) / (z[k+2]-z[k])
                # Reverse difference at upper boundary
                elseif k==length(z)
                    dvardz[j,k] = (3.0*field[j,k] - 4.0*field[j,k-1] +
                                   field[j,k-2]) / (z[k]-z[k-2])
                # Centered difference at all other grid points
                else
                    dvardz[j,k] = (field[j,k+1]-field[j,k-1]) /
                                  (z[k+1]-z[k-1])
                end
            end
        end

        return dvardz
    elseif ndims(field) == 3
        if length(z) != size(field)[3]
            error("Vector z and third dimension of variable to be differentiated
                   must be of same length")
        end
        dvardz = similar(field,Float64)
        # Loop over all dimensions to compute the vertical derivative
        for k in eachindex(z)
            for j in 1:size(field)[2]
                for i in 1:size(field)[1]
                    # Forward difference at lower boundary
                    if k==1
                        dvardz[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i,j,k+1] -
                                         field[i,j,k+2]) / (z[k+2]-z[k])
                    # Reverse difference at upper boundary
                    elseif k==length(z)
                        dvardz[i,j,k] = (3.0*field[i,j,k] - 4.0*field[i,j,k-1] +
                                         field[i,j,k-2]) / (z[k]-z[k-2])
                    # Centered difference at all other grid points
                    else
                        dvardz[i,j,k] = (field[i,j,k+1]-field[i,j,k-1]) /
                                        (z[k+1]-z[k-1])
                    end
                end
            end
        end

        return dvardz
    else
        error("Input variable to be differentiated cannot exceed 3 dimensions")
    end
end

#==============================================================================
finite_dt

Calculate the time derivative of a 1-D variable with the following criteria:

- Forward difference = First-order accurate
- Centered diffeference = Second-order accurate
- Backward difference = First-order accurate
==============================================================================#

function finite_dt(ti::AbstractVector{Ta},field::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}

    dvardt = zeros(length(ti))
    for n in eachindex(ti)
        if n == 1
            dvardt[n] = (field[n+1] - field[n]) / (ti[n+1] - ti[n]) 
        elseif n == length(ti)
            dvardt[n] = (field[n] - field[n-1]) / (ti[n] - ti[n-1])
        else 
            dvardt[n] = (field[n+1] - field[n-1]) / (ti[n+1] - ti[n-1])
        end
    end 
    return dvardt
end

#==============================================================================
finite_laplacian

Calculate the laplacian of a specified field on a Cartesian grid.
- Assumes x and y grid spacing are uniform (not necessarily equal)
- Assumes the x-dimension is the first dimension and that the y-dimension is
  the second dimension of the input array
==============================================================================#

function finite_laplacian(x::AbstractVector{Ta},y::AbstractVector{Tb},
                          field::AbstractArray{Tc}) where {Ta<:Real,Tb<:Real,Tc<:Real}

    if ndims(field) == 2
        if length(x) != size(field)[1] || length(y) != size(field)[2]
            error("Input variable to be differentiated must have dimensions of
                  (x,y)")
        end
        dx = (x[2]-x[1])^2
        dy = (y[2]-y[1])^2
        out = zeros(length(x),length(y))
        # Loop over all dimensions to compute the Laplacian
        for j in 2:length(y)-1
            for i in 2:length(x)-1
                out[i,j] = (field[i+1,j] - 2. * field[i,j] + field[i-1,j])/dx +
                           (field[i,j+1] -2. * field[i,j] + field[i,j-1])/dy
            end
        end
        return out
    elseif ndims(field) == 3
        if length(x) != size(field)[1] || length(y) != size(field)[2]
            error("Input variable to be differentiated must have dimensions of
                   (x,y,z)")
        end
        dx = (x[2]-x[1])^2
        dy = (y[2]-y[1])^2
        out = zeros(length(x),length(y))
        # Loop over all dimensions to compute the Laplacian
        for k in 1:size(field)[3]
            for j in 2:length(y)-1
                for i in 2:length(x)-1
                   out[i,j,k] = (field[i+1,j,k] -2. * field[i,j] + field[i-1,j,k])/dx +
                                (field[i,j+1,k] -2. * field[i,j] + field[i,j-1,k])/dy
                end
            end
        end
        return out
    end
end

#==============================================================================
kron_delta

Define the Kronecker Delta for constructing an lth-order derivative with
nth order of accuracy $\delta_{m,l}$ 

m = l ? delta = 1 : delta = 0 (where n = order; l = lth derivative)

- l = order of derivative (e.g., 1 for first derivative)
- n = order of accuracy (e.g., 2 for second-order accurate)
==============================================================================#

kron_delta(l::Int,n::Int) = [m == l ? 1. : 0. for m in 0:n+l-1]

#==============================================================================
fd_weights

Determine the second-order finite difference weights for either uniform or
non-uniform grids by solving a linear system of three equations
==============================================================================#

function fd_weights(x::AbstractVector{Ta};l::Int=1,n::Int=2) where Ta<:Real

    n == 2 ? nothing : error("Currently only supporting second-order accurate finite differences.")
    # Create the stencil given the order of accuracy
    stencil = zeros(n+1)
    # Array b constructed from Kronecker delta for lth-order derivative with nth-order accuracy
    b = kron_delta(l,n)
    wgts = zeros(n+1,length(x))
    for j in eachindex(x)
        if j == 1
            stencil = [0,1,2]
        elseif j == length(x)
            stencil = [0,-1,-2]
        else
            stencil = [-1,0,1]
        end
        # Matrix must be constructed for each grid point
        a = ones(n+1,n+1)
        for k in eachindex(stencil)
            for m in 0:n
                a[m+1,k] = (x[j+stencil[k]] - x[j])^m
            end
        end
        wgts[:,j] = a\b
    end
    return wgts
end

#==============================================================================
finite_ds

Calculate the second-order finite difference of a specified xyz field for
grids that are either uniform or non-uniform
Apply the second-order finite difference weights to a given field to determine
its derivative
General one-dimensional application
==============================================================================#

function finite_ds(wgts::AbstractArray{Ta,2},field::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}

    if size(wgts)[2] != length(field)
        error("There should be one set of weights for each grid point along
                   the dimension to be differentiated")
    end
    # Create the stencil given the order of accuracy
    stencil = zeros(size(wgts)[1])
    dvards = similar(field,Float64)
    for i in eachindex(field)
        if i == 1
            stencil = [i,i+1,i+2]
        elseif i == length(field)
            stencil = [i,i-1,i-2]
        else
            stencil = [i-1,i,i+1]
        end
        dvards[i] = sum(wgts[:,i] .* field[stencil])
    end
    return dvards
end

# x-dimension

function finite_dsx(wgts::AbstractArray{Ta,2},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

    # Two-dimensional variable
    if ndims(field) == 2
        if size(wgts)[2] != size(field)[1]
            error("There should be one set of weights for each grid point along
                   the dimension to be differentiated")
        end
        # Create the stencil given the order of accuracy
        stencil = zeros(size(wgts)[1])
        dvards = similar(field,Float64)
        for j in 1:size(field)[2]
            for i in 1:size(field)[1]
                if i == 1
                    stencil = [i,i+1,i+2]
                elseif i == size(field)[1]
                    stencil = [i,i-1,i-2]
                else
                    stencil = [i-1,i,i+1]
                end
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
        for k in 1:size(field)[3]
            for j in 1:size(field)[2]
                for i in 1:size(field)[1]
                    if i == 1
                        stencil = [i,i+1,i+2]
                    elseif i == size(field)[1]
                        stencil = [i,i-1,i-2]
                    else
                        stencil = [i-1,i,i+1]
                    end
                    dvards[i,j,k] = sum(wgts[:,i] .* field[stencil,j,k])
                end
            end
        end
        return dvards
    else
        error("Input variable to be differentiated cannot exceed 3 dimensions")
    end
end

# y-dimension

function finite_dsy(wgts::AbstractArray{Ta,2},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

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
                stencil = [j,j+1,j+2]
            elseif j == size(field)[2]
                stencil = [j,j-1,j-2]
            else
                stencil = [j-1,j,j+1]
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
        for k in 1:size(field)[3]
            for j in 1:size(field)[2]
                if j == 1
                    stencil = [j,j+1,j+2]
                elseif j == size(field)[2]
                    stencil = [j,j-1,j-2]
                else
                    stencil = [j-1,j,j+1]
                end
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

# z-dimension

function finite_dsz(wgts::AbstractArray{Ta,2},field::AbstractArray{Tb}) where {Ta<:Real,Tb<:Real}

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
                stencil = [k,k+1,k+2]
            elseif k == size(field)[2]
                stencil = [k,k-1,k-2]
            else
                stencil = [k-1,k,k+1]
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
                stencil = [k,k+1,k+2]
            elseif k == size(field)[3]
                stencil = [k,k-1,k-2]
            else
                stencil = [k-1,k,k+1]
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

#==============================================================================
diff1

Construct the (m+1) x m first derivative matrix D, not including the 1/ds factor
==============================================================================#

diff1(m::Real) = [[1.0 zeros(1,m-1)]; diagm(1=>ones(m-1)) - I]

#==============================================================================
spdiff1

Construct the sparse first derivative matrix D
==============================================================================#

spdiff1(m::Real) = sparse(diff1(m))

#==============================================================================
Laplacian

Construct the sparse, discrete 2-D Laplacian operator 
==============================================================================#

function Laplacian(x::AbstractVector{Ta},y::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    nx = length(x)
    ny = length(y)
    Dx = spdiff1(nx) / dx
    Dy = spdiff1(ny) / dy
    Ax = -transpose(Dx) * Dx
    Ay = -transpose(Dy) * Dy
    return kron(sparse(1.0I,ny,ny), Ax) + kron(Ay, sparse(1.0I,nx,nx))
end
