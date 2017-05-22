# *****************************************************************************
# integrate.jl
#
# Author: Jonathan Martinez 
# Email: jon.martinez@colostate.edu
# Julia version: 0.4.5
#
# This script contains a set of functions to integrate in different dimensions 
# using the trapezoidal method.
#
# Function list:
# trapz_1d
# trapz_2d
# trapz_3d
# *****************************************************************************

#==============================================================================
trapz_1d

Trapezoidal integration of a 1D radial or Cartesian grid that accounts for 
uneven grid spacing
See https://en.wikipedia.org/wiki/Trapezoidal_rule under the "Non-uniform grid"
** Can't be used for azimuthal integration (requires r*dphi)
** Assumes that x is in units of meters
Required input:
x = vector for x-dimension
y = dependent variable (vector with same size as x)
==============================================================================#

function trapz_1d{T<:Real}(x::Vector{T},y::Vector{T})
    local len::Int = length(y)
    if (len != length(x))
        error("Vectors must be of same length")
    end    
    sum = 0.0
    for i in 2:length(x)
        if !isnan(y[i])
            sum += (x[i] - x[i-1]) * (y[i] + y[i-1])
        end
    end
    
    return sum/2.0
end

#==============================================================================
trapz_2d

Trapezoidal integration of a 2D Cartesian grid that accounts for uneven grid 
spacing.
** Can't be used for polar coordinate integration (requires r*dr*dphi)
** Assumes x and y are in units of meters
Required input:
x   = vector for x-dimension
y   = vector for y-dimension
var = dependent variable with size [x,y]
==============================================================================#

function trapz_2d{T<:Real}(x::Vector{T},y::Vector{T},var::Array{T,2})
    int1 = Array(Float64,length(x)) .* NaN
    # Integrate along y-axis 
    for i in eachindex(x)
        int1[i] = trapz_1d(y,var[i,:])
    end
    # Integrate along x-axis
    int2 = trapz_1d(x,int1)

    return int2            
end

#==============================================================================
trapz_3d

Trapezoidal integration of a 3D Cartesian grid that accounts for uneven grid 
spacing.
** Can't be used for polar coordinate integration (requires r*dr*dphi*dz)
** Assumes x, y, and z are in units of meters
Required input:
x   = vector for x-dimension
y   = vector for y-dimension
z   = vector for z-dimension
var = dependent variable with size [x,y,z]
==============================================================================#

function trapz_3d{T<:Real}(x::Vector{T},y::Vector{T},z::Vector{T},
                           var::Array{T,3})
    int1 = Array(Float64,length(x),length(y)) .* NaN
    int2 = Array(Float64,length(x)) .* NaN
    # Integrate along z-axis 
    for i in eachindex(x)
        for j in eachindex(y)
            int1[i,j] = trapz_1d(z,var[i,j,:])
        end
        # Integrate along y-axis
        int2[i] = trapz_1d(y,int1[i,:])
    end
    # Integrate along x-axis
    int3 = trapz_1d(x,int2)
    return int3            
end



