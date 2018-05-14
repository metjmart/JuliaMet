# *****************************************************************************
# vortexprofs.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.6.0
#
# This script contains functions to create radial profiles of tangential wind 
# for vortices widely used to initiate model simulations of tropical cyclones.
# References are included where necessary.
# 
# rankine
# modrankine
# hermite
# re87
# cw87
# wc04
# *****************************************************************************

#==============================================================================
rankine

This function will generate a 1-D Rankine vortex for a given rmax and vmax
*** Assumes r and rmax are in units of meters and vmax is m/s!!!
===============================================================================#

function rankine(radii::AbstractVector{<:Real},rmax::Real,vmax::Real,method::Symbol=:vt)
    if method == :vt
        return @fastmath [r < rmax ? vmax * r / rmax : rmax * vmax / r for r in radii]
    elseif method == :vort
        return @fastmath [r < rmax ? 2. * vmax / rmax : 0. for r in radii]
    end
end

#==============================================================================
modrankine

This function will generate a 1-D modified Rankine vortex for a given rmax, 
vmax, and decay parameter (alpha)
*** Assumes r and rmax are in units of meters, vmax is m/s, and alpha is postive
===============================================================================#

function modrankine(radii::AbstractVector{<:Real},rmax::Real,vmax::Real,alpha::Real)
    return @fastmath [r < rmax ? vmax * r / rmax : vmax * (r/rmax)^-alpha for r in radii]
end

#==============================================================================
hermite

Basic cubic Hermite shape function satisfying S(0) = 1, S(1) = 0, and 
S'(0) = S'(1) = 0.
Can be used to provide smooth transition regions for, e.g., vorticity rings.
See Schubert et al. (1999; JAS) for examples
===============================================================================#

# Single data points

function hermite(radius::Real)
    return @fastmath 1. - 3. * radius^2 + 2. * radius^3
end

# Vectors

function hermite(radius::AbstractVector{<:Real})
    return @fastmath 1. - 3. * radius.^2 + 2. * radius.^3
end

#==============================================================================
re87

This function will generate a 1-D vortex profile 
following Rotunno and Emanuel (1987, JAS - equation 37)
Need to specify radius array, rmax, r0, vmax, coriolis force
*** Assumes r, rmax, and r0 are in units of meters, vmax is m/s, and coriolis 
    force is s^-1
===============================================================================#

function re87(radii::AbstractVector{<:Real},rmax::Real,r0::Real,vmax::Real,fcor::Real)
    @fastmath [r < r0 ? sqrt( vmax^2 * (r / rmax)^2 * ( ( (2. * rmax) ./ (r + rmax) )^3 - ( (2. * rmax) / (r0 + rmax) )^3 ) +
                              (fcor^2 * r^2 / 4.) ) - ( (fcor * r) / 2. ) : 0. for r in radii]
end

#==============================================================================
cw87

This function will generate a 1-D vortex profile following 
Chan and Williams (1987, JAS - equation 2.10)
Need to specify radius array, rmax, vmax, and b
*** Assumes r and rmax are in units of meters, vmax is m/s
    b is a unitless parameter to control the shape of the outer wind profile
===============================================================================#

function cw87(r::AbstractVector{<:Real},rmax::Real,vmax::Real,b::Real)
    return @fastmath vmax * (r/rmax) .* exp.( (1/b) * (1 - (r/rmax).^b) )
end

#==============================================================================
wc04

This function will generate a 1-D vortex profile following 
Wong and Chan (2004, JAS - equation 2.2)
Need to specify radius array, rmax, vmax, b, and rtc
*** Assumes r, rmax, and rtc  are in units of meters, vmax is m/s
    b is a unitless parameter to control the shape of the outer wind profile
===============================================================================#

function wc04(radii::AbstractVector{<:Real},rmax::Real,vmax::Real,b::Real,rtc::Real)
    return [r < rtc ? vmax * (r/rmax) .* ( exp.( (1 - ((r/rmax)).^b)/b ) - 
                      abs.(r - rmax) ./ (rtc-rmax) .* exp.( (1 - ((rtc/rmax)).^b)/b ) ) : 0. for r in radii]
end




