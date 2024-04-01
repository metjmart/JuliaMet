# *****************************************************************************
# integrate.jl

# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# This script contains a set of functions for integrating discretized data
#
# Function list
# newtoncotes
# trapz2d
# trapz3d
# *****************************************************************************

#==============================================================================
newtoncotes

Family of integration methods from the closed Newton-Cotes formulae -- see
http://mathworld.wolfram.com/Newton-CotesFormulas.html for reference.
** These integration methods assume that the data are evenly spaced.

** Need to proceed with caution if attempting to integrate observational
   data where NaNs are present. The higher-degree integration methods require
   a certain number of sub-intervals for the integration to be carried out
   accurately. If the data are not continuously available, these requirements
   will likely be invalidated and result in an incorrect solution. Therefore,
   it is recommended that the trapezoidal method is used for integrating
   data interspersed with NaNs as the trapezoidal method does not place a
   restriction on the number of sub-intervals required. For this reason, only
   the trapezoidal method works in the presence of NaNs. If NaNs are only
   present at the start and end of an analysis with continuous data in between,
   the user can splice the data such that the NaNs are removed and the number
   of sub-intervals present from start to end match the required number of
   sub-intervals for a given integration method.

Methods

:trapz   - trapezoidal (2-point) rule
:simps13 - Simpson's 1/3 (3-point) rule (requires even sub-intervals)
:simps38 - Simpson's 3/8 (4-point) rule (requires multiple of 3 sub-intervals)
:boole   - Boole's (5-point) rule (requires multiple of 4 sub-intervals)

==============================================================================#

function newtoncotes(x::AbstractVector{Ta},y::AbstractVector{Tb},
                     method::Symbol=:trapz) where {Ta<:Real,Tb<:Real}

    length(y) == length(x) ? nothing :
        throw(DimensionMismatch("Input vectors must have same length"))
    # Define the number of sub-intervals
    n = length(y)-1
    numer = 0.0
    # Trapezoidal rule
    if method == :trapz
        n < 1 ? error("Must have at least one sub-interval for integration") : nothing
        for i in 2:length(x)
            if !isnan(y[i]) && !isnan(y[i-1])
                @inbounds numer += (x[i] - x[i-1]) * (y[i] + y[i-1])
            end
        end
        return numer/2.0
    # Simpson's 1/3 rule
    elseif method == :simps13
        if !iseven(n)
            error("Number of sub-intervals must be even")
        elseif any(isnan.(y))
            error("NaNs are currently not handled in Simpson's 1/3 rule,
                   consider using the trapezoidal rule")
        end
        for i in 3:2:length(x)
            @fastmath @inbounds numer += (x[i] - x[i-2]) * (y[i] + 4.0*y[i-1] + y[i-2])
        end
        return numer/6.0
    # Simpson's 3/8 rule
    elseif method == :simps38
        if n % 3 != 0
            error("Number of sub-intervals must be a multiple of 3")
        elseif any(isnan.(y))
            error("NaNs are currently not handled in Simpson's 3/8 rule,
                   consider using the trapezoidal rule")
        end
        for i in 4:3:length(x)
            @fastmath @inbounds numer += (x[i] - x[i-3]) * (y[i] + 3.0*y[i-1] + 3.0*y[i-2] + y[i-3])
        end
        return numer/8.0
    # Boole's rule
    elseif method == :boole
        if n % 4 != 0
            error("Number of sub-intervals must be a multiple of 4")
        elseif any(isnan.(y))
            error("NaNs are currently not handled in Boole's rule,
                   consider using the trapezoidal rule")
        end
        for i in 5:4:length(x)
            @fastmath @inbounds numer += (x[i] - x[i-4]) * (7.0*y[i] + 32.0*y[i-1] + 12.0*y[i-2] +
                                        32.0*y[i-3] + 7.0*y[i-4])
        end
        return numer/90.0
    end
end

#==============================================================================
trapz2d

Trapezoidal integration of a 2D Cartesian grid
** Assumes x and y are in units of meters
Required input:
x   = vector for x-dimension
y   = vector for y-dimension
field = dependent variable with size (x,y)
==============================================================================#

function trapz2d(x::AbstractVector{Ta},y::AbstractVector{Tb},
                 field::AbstractArray{Tc,2}) where {Ta<:Real,Tb<:Real,Tc<:Real}

    size(field)[1] == length(x) && size(field)[2] == length(y) ? nothing :
        throw(DimensionMismatch("Input variable to be integrated must have
                                 dimensions of (x,y)"))
    int1 = similar(x,Float64)
    fill!(int1, NaN)
    # Integrate along y-axis
    for i in eachindex(x)
        @inbounds int1[i] = newtoncotes(y,field[i,:])
    end
    # Integrate along x-axis
    int2 = newtoncotes(x,int1)
    return int2
end

#==============================================================================
trapz3d

Trapezoidal integration of a 3D Cartesian grid
** Assumes x, y, and z are in units of meters
Required input:
x   = vector for x-dimension
y   = vector for y-dimension
z   = vector for z-dimension
field = dependent variable with size (x,y,z)
==============================================================================#

function trapz3d(x::AbstractVector{Ta},y::AbstractVector{Tb},
                 z::AbstractVector{Tc},field::AbstractArray{Td,3}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}

    size(field)[1] == length(x) && size(field)[2] == length(y) && size(field)[3] == length(z) ? nothing :
        throw(DimensionMismatch("Input variable to be integrated must have dimensions
                                 of (x,y,z)"))
    int1 = Array{Float64}(undef,length(x),length(y))
    fill!(int1, NaN)
    int2 = similar(x)
    fill!(int2, NaN)
    # Integrate along z-axis
    for i in eachindex(x)
        for j in eachindex(y)
            @inbounds int1[i,j] = newtoncotes(z,field[i,j,:])
        end
        # Integrate along y-axis
        @inbounds int2[i] = newtoncotes(y,int1[i,:])
    end
    # Integrate along x-axis
    int3 = newtoncotes(x,int2)
    return int3
end
