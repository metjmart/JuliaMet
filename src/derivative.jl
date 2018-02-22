# *****************************************************************************
# derivative.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia Version: 0.6.0
#
# This script contains several functions used to compute derivatives with 
# finite differencing using second-order accurate (3-point) stencils. 
# The grid points are assumed to be uniform.
# 
# Function list:
# finite_dr 
# finite_dx
# finite_dy
# finite_dz 
# finite_laplacian
# *****************************************************************************

#==============================================================================
finite_dr

Calculate the radial derivative of a specified rtz field.
The function requires the input (r) to have units of meters.
*** Assumes the radial dimension is the first dimension in the input array!
==============================================================================#

function finite_dr(r::AbstractVector{<:Real},field::AbstractArray{<:Real})

    if ndims(field) == 1 
        if length(r) != length(field)
            error("Vector r and input variable to be differentiated must be of 
                   same length")
        end
        dvardr = similar(field,Float64)
        fill!(dvardr,NaN)
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
        # Define the dimensions of the input variable
        d1,d2 = size(field)
        if length(r) != d1
            error("Vector r and first dimension of variable to be differentiated 
                   must be of same length")
        end 
        dvardr = similar(field,Float64)
        fill!(dvardr, NaN)
        # Loop over all dimensions to compute the radial derivative
        for j in 1:d2
            for i in 1:d1
                # Forward difference at innermost boundary
                if i==1
                    dvardr[i,j] = (-3.0*field[i,j] + 4.0*field[i+1,j] - 
                                   field[i+2,j]) / (r[i+2]-r[i])
                # Reverse difference at outermost boundary
                elseif i==d1
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
        # Define the dimensions of the input variable
        d1,d2,d3 = size(field)
        if length(r) != d1
            error("Vector r and first dimension of variable to be differentiated 
                   must be of same length")
        end 
        dvardr = similar(field,Float64)
        fill!(dvardr, NaN)
        # Loop over all dimensions to compute the radial derivative
        for k in 1:d3
            for j in 1:d2
                for i in 1:d1
                    # Forward difference at innermost boundary
                    if i==1
                        dvardr[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i+1,j,k] - 
                                         field[i+2,j,k]) / (r[i+2]-r[i])
                    # Reverse difference at outermost boundary
                    elseif i==d1
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
The function requires the input (x) to have units of meters.
*** Assumes the x-dimension is the first dimension in the input array!
==============================================================================#

function finite_dx(x::AbstractVector{<:Real},field::AbstractArray{<:Real})

    if ndims(field) == 1 
        if length(x) != length(field)
            error("Vector x and input variable to be differentiated must be of 
                   same length")
        end
        dvardx = similar(field,Float64)
        fill!(dvardx,NaN)
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
        # Define the dimensions of the input variable
        d1,d2 = size(field)
        if length(x) != d1
            error("Vector x and first dimension of variable to be differentiated 
                   must be of same length")
        end 
        dvardx = similar(field,Float64)
        fill!(dvardx, NaN)
        # Loop over all dimensions to compute the x-derivative
        for j in 1:d2
            for i in 1:d1
                # Forward difference at innermost boundary
                if i==1
                    dvardx[i,j] = (-3.0*field[i,j] + 4.0*field[i+1,j] - 
                                   field[i+2,j]) / (x[i+2]-x[i])
                # Reverse difference at outermost boundary
                elseif i==d1
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
        # Define the dimensions of the input variable
        d1,d2,d3 = size(field)
        if length(x) != d1
            error("Vector x and first dimension of variable to be differentiated 
                   must be of same length")
        end 
        dvardx = similar(field,Float64)
        fill!(dvardx, NaN)
        # Loop over all dimensions to compute the x-derivative
        for k in 1:d3
            for j in 1:d2
                for i in 1:d1
                    # Forward difference at innermost boundary
                    if i==1
                        dvardx[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i+1,j,k] - 
                                         field[i+2,j,k]) / (x[i+2]-x[i])
                    # Reverse difference at outermost boundary
                    elseif i==d1
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
The function requires the input (y) to have units of meters.
*** Assumes the y-dimension is the first dimension in the input array!
==============================================================================#

function finite_dy(y::AbstractVector{<:Real},field::AbstractArray{<:Real})

    if ndims(field) == 1
        if length(y) != length(field)
            error("Vector y and input variable to be differentiated must be of 
                   same length")
        end 
        dvardy = similar(field,Float64)
        fill!(dvardy, NaN)
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
        # Define the dimensions of the input variable
        d1,d2 = size(field)
        if length(y) != d2 
            error("Vector y and second dimension of the variable to differentiated
                   must be of same length")
        end 
        dvardy = similar(field,Float64)
        fill!(dvardy, NaN)
        # Loop over all dimensions to compute the y-derivative
        for j in 1:d2
            for i in 1:d1
                # Forward difference at innermost boundary
                if j==1
                    dvardy[i,j] = (-3.0*field[i,j] + 4.0*field[i,j+1] - 
                                   field[i,j+2]) / (y[j+2]-y[j])
                # Reverse difference at outermost boundary
                elseif j==d2
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
        # Define the dimensions of the input variable
        d1,d2,d3 = size(field)
        if length(y) != d2 
            error("Vector y and second dimension of the variable to differentiated
                   must be of same length")
        end 
        dvardy = similar(field,Float64)
        fill!(dvardy, NaN)
        # Loop over all dimensions to compute the y-derivative
        for k in 1:d3
            for j in 1:d2
                for i in 1:d1
                    # Forward difference at innermost boundary
                    if j==1
                        dvardy[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i,j+1,k] - 
                                         field[i,j+2,k]) / (y[j+2]-y[j])
                    # Reverse difference at outermost boundary
                    elseif j==d2
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

Calculate the vertical derivative of a specified rtz or xyz field. 
The function requires the input (z) to have units of meters.
*** Assumes the vertical dimension is the third dimension in the input array!
==============================================================================#

function finite_dz(z::AbstractVector{<:Real},field::AbstractArray{<:Real})

    if ndims(field) == 1 
        if length(z) != length(field)
            error("Vector z and input variable to be differentiated must be of 
                   same length")
        end 
        dvardz = similar(field,Float64)
        fill!(dvardz, NaN)
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
        # Define the dimensions of the input variable
        d1,d2 = size(field)
        if length(z) != d2
            error("Vector z and second dimension of variable to be differentiated 
                   must be of same length")
        end 
        dvardz = similar(field,Float64)
        fill!(dvardz, NaN)
        # Loop over all dimensions to compute the vertical derivative
        for k in 1:d2
            for j in 1:d1
                # Forward difference at lower boundary
                if k==1
                    dvardz[j,k] = (-3.0*field[j,k] + 4.0*field[j,k+1] - 
                                   field[j,k+2]) / (z[k+2]-z[k])
                # Reverse difference at upper boundary
                elseif k==d2
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
        # Define the dimensions of the input variable
        d1,d2,d3 = size(field)
        if length(z) != d3 
            error("Vector z and third dimension of variable to be differentiated 
                   must be of same length")
        end 
        dvardz = similar(field,Float64)
        fill!(dvardz, NaN)
        # Loop over all dimensions to compute the vertical derivative
        for k in 1:d3
            for j in 1:d2
                for i in 1:d1
                    # Forward difference at lower boundary
                    if k==1
                        dvardz[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i,j,k+1] - 
                                         field[i,j,k+2]) / (z[k+2]-z[k])
                    # Reverse difference at upper boundary
                    elseif k==d3
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
finite_laplacian

Calculate the laplacian of a specified field on a Cartesian grid.
The function requires the input (x) and (y) to have units of meters.
*** Assumes x and y grid spacing are identical
*** Assumes the x-dimension is the first dimension and that the y-dimension is 
    the second dimension of the input array!
==============================================================================#

function finite_laplacian(x::AbstractVector{<:Real},y::AbstractVector{<:Real},
                          field::AbstractArray{<:Real})

    if ndims(field) == 2
        # Define the dimensions of the input variable
        d1,d2 = size(field) 
        if length(x) != d1 || length(y) != d2
            error("Input variable to be differentiated must have dimensions of  
                  [x,y]")
        end
        hx = (x[2]-x[1])^2
        hy = (y[2]-y[1])^2
        out = similar(field,Float64)
        fill!(out, NaN)
        # Loop over all dimensions to compute the Laplacian
        for j in 2:d2-1
            for i in 2:d1-1
                out[i,j] = (field[i+1,j] + field[i-1,j])/hx + 
                           (field[i,j+1] + field[i,j-1])/hy 
            end
        end

        return out
    elseif ndims(field) == 3 
        # Define the dimensions of the input variable
        d1,d2,d3 = size(field) 
        if length(x) != d1 || length(y) != d2
            error("Input variable to be differentiated must have dimensions of  
                   [x,y,z]")
        end 
        hx = (x[2]-x[1])^2
        hy = (y[2]-y[1])^2
        out = similar(field,Float64)
        fill!(out, NaN)
        # Loop over all dimensions to compute the Laplacian
        for k in 1:d3
            for j in 2:d2-1
                for i in 2:d1-1
                   out[i,j,k] = (field[i+1,j,k] + field[i-1,j,k])/hx + 
                                (field[i,j+1,k] + field[i,j-1,k])/hy 
                end
            end
        end

        return out
    else 
        error("Input variable to be differentiated must be 2-D or 3-D")
    end 
end

