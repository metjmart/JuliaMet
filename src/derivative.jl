# *****************************************************************************
# derivative.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia Version: 0.4.5
#
# This script contains several functions used to calculate derivatives with 
# second-order accurate (3-point) stencils. 
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

function finite_dr{Ta<:Real,Tb<:Real}(r::AbstractVector{Ta},
                                      field::AbstractArray{Tb,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field)
    dvardr = similar(field,Float64)
    fill!(dvardr, NaN)
    # Loop over all dimensions to compute the radial derivative
    for k in 1:d3
        for j in 1:d2
            for i in 1:d1
                # Forward difference at innermost boundary
                if i==1
                    dvardr[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i+1,j,k] - field[i+2,j,k])/(r[i+2]-r[i])
                # Reverse difference at outermost boundary
                elseif i==d1
                    dvardr[i,j,k] = (3.0*field[i,j,k] - 4.0*field[i-1,j,k] + field[i-2,j,k])/(r[i]-r[i-2])
                # Centered difference at all other grid points
                else
                    dvardr[i,j,k] = (field[i+1,j,k]-field[i-1,j,k])/(r[i+1]-r[i-1])
                end
            end
        end
    end

   return dvardr
end

#==============================================================================
finite_dx

Calculate the x-horizontal derivative of a specified xyz field. 
The function requires the input (x) to have units of meters.
*** Assumes the x-dimension is the first dimension in the input array!
==============================================================================#

function finite_dx{Ta<:Real,Tb<:Real}(x::AbstractVector{Ta},
                                      field::AbstractArray{Tb,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field)
    dvardx = similar(field,Float64)
    fill!(dvardx, NaN)
    # Loop over all dimensions to compute the x-derivative
    for k in 1:d3
        for j in 1:d2
            for i in 1:d1
                # Forward difference at innermost boundary
                if i==1
                    dvardx[i,j,k] = (-3.0*field[i+1,j,k] + 4.0*field[i+1,j,k] - field[i+2,j,k])/(x[i+2]-x[i])
                # Reverse difference at outermost boundary
                elseif i==d1
                    dvardx[i,j,k] = (3.0*field[i,j,k] - 4.0*field[i-1,j,k] + field[i-2,j,k])/(x[i]-x[i-2])
                # Centered difference at all other grid points
                else
                    dvardx[i,j,k] = (field[i+1,j,k]-field[i-1,j,k])/(x[i+1]-x[i-1])
                end
            end
        end
    end

   return dvardx
end

#==============================================================================
finite_dy

Calculate the y-horizontal derivative of a specified xyz field.
The function requires the input (y) to have units of meters.
*** Assumes the y-dimension is the first dimension in the input array!
==============================================================================#

function finite_dy{Ta<:Real,Tb<:Real}(y::AbstractVector{Ta},
                                      field::AbstractArray{Tb,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field)
    dvardy = similar(field,Float64)
    fill!(dvardy, NaN)
    # Loop over all dimensions to compute the y-derivative
    for k in 1:d3
        for j in 1:d2
            for i in 1:d1
                # Forward difference at innermost boundary
                if j==1
                    dvardy[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i,j+1,k] - field[i,j+2,k])/(y[j+2]-y[j])
                # Reverse difference at outermost boundary
                elseif j==d2
                    dvardy[i,j,k] = (3.0*field[i,j,k] - 4.0*field[i,j-1,k] + field[i,j-2,k])/(y[j]-y[j-2])
                # Centered difference at all other grid points
                else
                    dvardy[i,j,k] = (field[i,j+1,k]-field[i,j-1,k])/(y[j+1]-y[j-1])
                end
            end
        end
    end

   return dvardy
end

#==============================================================================
finite_dz

Calculate the vertical derivative of a specified rtz or xyz field. 
The function requires the input (z) to have units of meters.
*** Assumes the vertical dimension is the third dimension in the input array!
==============================================================================#

function finite_dz{Ta<:Real,Tb<:Real}(z::AbstractVector{Ta},
                                      field::AbstractArray{Tb,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field)
    dvardz = similar(field,Float64)
    fill!(dvardz, NaN)
    # Loop over all dimensions to compute the vertical derivative
    for k in 1:d3
        for j in 1:d2
            for i in 1:d1
                # Forward difference at lower boundary
                if k==1
                    dvardz[i,j,k] = (-3.0*field[i,j,k] + 4.0*field[i,j,k+1] - field[i,j,k+2])/(z[k+2]-z[k])
                # Reverse difference at upper boundary
                elseif k==d3
                    dvardz[i,j,k] = (3.0*field[i,j,k] - 4.0*field[i,j,k-1] + field[i,j,k-2])/(z[k]-z[k-2])
                # Centered difference at all other grid points
                else
                    dvardz[i,j,k] = (field[i,j,k+1]-field[i,j,k-1])/(z[k+1]-z[k-1])
                end
            end
        end
    end

    return dvardz
end

#==============================================================================
finite_laplacian

Calculate the laplacian of a specified xyz field.
The function requires the input (x) and (y) to have units of meters.
*** Assumes the x-dimension is the first dimension and that the y-dimension is 
    the second dimension of the input array!
*** Assumes x and y grid spacing are identical
==============================================================================#

function finite_laplacian{Ta<:Real,Tb<:Real,Tc<:Real}(
                          x::AbstractVector{Ta},
                          y::AbstractVector{Tb},
                          field::AbstractArray{Tc,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field) 
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
end

