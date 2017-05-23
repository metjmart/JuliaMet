# *****************************************************************************
# derivative.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia Version: 0.4.5
#
# This script contains several functions used to calculate derived quantities. 
# These quantities may be derivatives of a given field or simply a quantity
# derived from available vars. Read each function description for specific 
# details
# 
# Function list:
# finite_dz 
# finite_dr 
# finite_dx
# finite_dy
# *****************************************************************************

#==============================================================================
finite_dz

Calculate the vertical derivative of a specified rtz or xyz field using finite 
differencing. The function requires the input (z) to have units of meters.
*** Assumes the vertical dimension is the third dimension in the input array!
==============================================================================#

function finite_dz{T<:Real}(z::Vector{T},field::Array{T,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field)
    dvardz = Array(Float64,size(field)) .* NaN
    # Loop over all dimensions to compute the vertical derivative
    for k in collect(1:d3)
        for j in collect(1:d2)
            for i in collect(1:d1)
                # Forward difference at lower boundary
                if k==1
                    dvardz[i,j,k] = (field[i,j,k+1]-field[i,j,k])/(z[k+1]-z[k])
                # Reverse difference at upper boundary
                elseif k==d3
                    dvardz[i,j,k] = (field[i,j,k]-field[i,j,k-1])/(z[k]-z[k-1])
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
finite_dr

Calculate the radial derivative of a specified rtz field using finite 
differencing. The function requires the input (r) to have units of meters.
*** Assumes the radial dimension is the first dimension in the input array!
==============================================================================#

function finite_dr{T<:Real}(r::Vector{T},field::Array{T,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field)
    dvardr = Array(Float64,size(field)) .* NaN;
    # Loop over all dimensions to compute the radial derivative
    for i in collect(1:d1)
        for j in collect(1:d2)
            for k in collect(1:d3)
                # Forward difference at innermost boundary
                if i==1
                    dvardr[i,j,k] = (field[i+1,j,k]-field[i,j,k])/(r[i+1]-r[i])
                # Reverse difference at outermost boundary
                elseif i==d1
                    dvardr[i,j,k] = (field[i,j,k]-field[i-1,j,k])/(r[i]-r[i-1])
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

Calculate the x-horizontal derivative of a specified xyz field using finite 
differencing. The function requires the input (x) to have unites of meters.
*** Assumes the x-dimension is the first dimension in the input array!
==============================================================================#

function finite_dx{T<:Real}(x::Vector{T},field::Array{T,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field)
    dvardx = Array(Float64,size(field)) .* NaN;
    # Loop over all dimensions to compute the radial derivative
    for i in collect(1:d1)
        for j in collect(1:d2)
            for k in collect(1:d3)
                # Forward difference at innexmost boundary
                if i==1
                    dvardx[i,j,k] = (field[i+1,j,k]-field[i,j,k])/(x[i+1]-x[i])
                # Reverse difference at outermost boundary
                elseif i==d1
                    dvardx[i,j,k] = (field[i,j,k]-field[i-1,j,k])/(x[i]-x[i-1])
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

Calculate the y-horizontal derivative of a specified xyz field using finite 
differencing. The function requires the input (y) to have units of meters.
*** Assumes the y-dimension is the first dimension in the input array!
==============================================================================#

function finite_dy{T<:Real}(y::Vector{T},field::Array{T,3})
    # Define the dimensions of the input variable
    d1,d2,d3 = size(field)
    dvardy = Array(Float64,size(field)) .* NaN;
    # Loop over all dimensions to compute the radial derivative
    for i in collect(1:d1)
        for j in collect(1:d2)
            for k in collect(1:d3)
                # Forward difference at innermost boundary
                if j==1
                    dvardy[i,j,k] = (field[i,j+1,k]-field[i,j,k])/(y[j+1]-y[j])
                # Reverse difference at outermost boundary
                elseif j==d2
                    dvardy[i,j,k] = (field[i,j,k]-field[i,j-1,k])/(y[j]-y[j-1])
                # Centered difference at all other grid points
                else
                    dvardy[i,j,k] = (field[i,j+1,k]-field[i,j-1,k])/(y[j+1]-y[j-1])
                end
            end
        end
    end

   return dvardy
end




