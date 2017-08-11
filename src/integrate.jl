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
# newton_cotes
# trapz_1d
# trapz_2d
# trapz_3d
# *****************************************************************************

#==============================================================================
newton_cotes

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

function newton_cotes{Ta<:Real,Tb<:Real}(x::AbstractVector{Ta},
                                         y::AbstractVector{Tb},
                                         method::Symbol=:trapz)
    len::Int = length(y)
    if (len != length(x))
       error("Vectors x and y must be of same length")
    end
    n::Int = length(y)-1
    sum = 0.0

    if method == :trapz
        if n < 2
            error("Must have at least one sub-interval for integration")
        end
        for i in 2:length(x)
            if !isnan(y[i]) && !isnan(y[i-1])
                sum += (x[i] - x[i-1]) * (y[i] + y[i-1])
            end
        end

    return sum/2.0

    elseif method == :simps13
        if !iseven(n)
            error("Number of sub-intervals must be even")
        elseif any(isnan(y))
            error("NaNs are currently not handled in Simpson's 1/3 rule, 
                   consider using the trapezoidal rule")
        end
        for i in 3:2:length(x)
            sum += (x[i] - x[i-2]) * (y[i] + 4.0*y[i-1] + y[i-2])
        end
    
        return sum/6.0

    elseif method == :simps38
        if n % 3 != 0
            error("Number of sub-intervals must be a multiple of 3")
        elseif any(isnan(y))
            error("NaNs are currently not handled in Simpson's 3/8 rule, 
                   consider using the trapezoidal rule")
        end 
        for i in 4:3:length(x)
            sum += (x[i] - x[i-3]) * (y[i] + 3.0*y[i-1] + 3.0*y[i-2] + y[i-3])
        end  

        return sum/8.0 

    elseif method == :boole
        if n % 4 != 0 
            error("Number of sub-intervals must be a multiple of 4")
        elseif any(isnan(y))
            error("NaNs are currently not handled in Boole's rule, 
                   consider using the trapezoidal rule")
        end  
        for i in 5:4:length(x)
            sum += (x[i] - x[i-4]) * (7.0*y[i] + 32.0*y[i-1] + 12.0*y[i-2] + 
                                      32.0*y[i-3] + 7.0*y[i-4])
        end
    
        return sum/90.0

    end

end

#==============================================================================
trapz_1d

Trapezoidal integration of a 1D radial or Cartesian grid 
See https://en.wikipedia.org/wiki/Trapezoidal_rule 
** Assumes that x is in units of meters
Required input:
x = vector for x-dimension
y = dependent variable to be integrated (vector with same size as x)
==============================================================================#

function trapz_1d{Ta<:Real,Tb<:Real}(x::AbstractVector{Ta},
                                     y::AbstractVector{Tb})
    len::Int = length(y)
    if (len != length(x))
        error("Vectors must be of same length")
    end    
    sum = 0.0
    for i in 2:length(x)
        if !isnan(y[i]) && !isnan(y[i-1])
            sum += (x[i] - x[i-1]) * (y[i] + y[i-1])
        end
    end
    
    return sum/2.0
end

#==============================================================================
trapz_2d

Trapezoidal integration of a 2D Cartesian grid 
** Assumes x and y are in units of meters
Required input:
x   = vector for x-dimension
y   = vector for y-dimension
var = dependent variable with size [x,y]
==============================================================================#

function trapz_2d{Ta<:Real,Tb<:Real,Tc<:Real}(x::AbstractVector{Ta},
                                              y::AbstractVector{Tb},
                                              var::AbstractArray{Tc,2})
    int1 = similar(x,Float64)
    fill!(int1, NaN)
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

Trapezoidal integration of a 3D Cartesian grid 
** Assumes x, y, and z are in units of meters
Required input:
x   = vector for x-dimension
y   = vector for y-dimension
z   = vector for z-dimension
var = dependent variable with size [x,y,z]
==============================================================================#

function trapz_3d{Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}(x::AbstractVector{Ta},
                                                       y::AbstractVector{Tb},
                                                       z::AbstractVector{Tc},
                                                       var::AbstractArray{Td,3})
    int1 = Array(Float64,length(x),length(y)) 
    fill!(int1, NaN)
    int2 = similar(x)
    fill!(int2, NaN)
    # Integrate along z-axis 
    for i in eachindex(x)
        for j in eachindex(y)
            int1[i,j] = trapz_1d(z,squeeze(squeeze(var[i,j,:],2),1))
        end
        # Integrate along y-axis
        int2[i] = trapz_1d(y,squeeze(int1[i,:],1))
    end
    # Integrate along x-axis
    int3 = trapz_1d(x,int2)
    return int3            
end



