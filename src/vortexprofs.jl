# *****************************************************************************
# vortexprofs.jl
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# This script contains functions to create radial profiles of tangential wind
# for vortices widely used to initiate model simulations of tropical cyclones.
# References are included where necessary.
#
# Function list
#
# rankine
# smrankine
# modrankine
# hermite
# re87
# cw87
# wc04
# *****************************************************************************

#==============================================================================
rankine

This function will generate the tangential velocity or relative vorticity
profile for a 1-D Rankine vortex

Input (and units)
radii - vector containing radius grid points (m)
rmax - radius of maximum tangential winds (m)
vmax - maximum tangential velocity (m/s)
method - output 1-D tangential velocity profile (:vt) or relative vorticity
         profile (:vort)
===============================================================================#

function rankine(radii::AbstractVector{T},rmax::Real,vmax::Real,method::Symbol=:vt) where T<:Real
    if method == :vt
        return @fastmath [r < rmax ? vmax * r / rmax : rmax * vmax / r for r in radii]
    elseif method == :vort
        return @fastmath [r < rmax ? 2. * vmax / rmax : 0. for r in radii]
    end
end

#==============================================================================
modrankine

This function will generate the tangential velocity or relative vorticity
profile for a 1-D modified Rankine vortex

Input (and units)
radii - vector containing radius grid points (m)
rmax - radius of maximum tangential winds (m)
vmax - maximum tangential velocity (m/s)
alpha - decay parameter (unitless; positive value)
method - output 1-D tangential velocity profile (:vt) or relative vorticity
         profile (:vort)
===============================================================================#

function modrankine(radii::AbstractVector{T},rmax::Real,vmax::Real,alpha::Real,method::Symbol=:vt) where T<:Real
    if method == :vt
        return @fastmath [r < rmax ? vmax * r / rmax : vmax * (r/rmax)^-alpha for r in radii]
    elseif method == :vort
        return @fastmath [r < rmax ? 2 .* vmax / rmax : (1 - alpha) * (vmax * rmax^alpha * r^-(1+alpha)) for r in radii]
    end
end

#==============================================================================
smrankine

This function will generate the relative vorticity profile for a 1-D smoothed
Rankine vortex following Schecter and Montgomery (2004; Phys. Fluids, eq. 37)

Input (and units)
radii - vector containing radius grid points (m)
rknot - radius where central vorticity decreases (m; exact if delta = 0,
        yielding a pure Rankine vortex)
zknot - central vorticity (s^-1)
delta - smoothing parameter (unitless; larger values = smoother profile)
===============================================================================#

function smrankine(radii::AbstractVector{T},rknot::Real,zknot::Real,delta::Real) where T<:Real
    return @fastmath zknot/2.0*[1.0 - tanh((r - rknot)/(rknot * delta))  for r in radii]
end

#==============================================================================
hermite

Basic cubic Hermite shape function satisfying S(0) = 1, S(1) = 0, and
S'(0) = S'(1) = 0.
Can be used to create smooth transition regions (e.g., vorticity rings)
See Schubert et al. (1999; JAS) for examples
===============================================================================#

# Single data points--use dot broadcasting for vectors and arrays

function hermite(radii::Real)
    return @fastmath 1. - 3. * radii^2 + 2. * radii^3
end

#==============================================================================
re87

This function will generate the tangential velocity profile for
a 1-D vortex following Rotunno and Emanuel (1987; JAS, equation 37)

Input (and units)
radii - vector containing radius grid points (m)
rmax - radius of maximum tangential winds (m)
r0 - radius where tangential velocity decreases to zero (m)
vmax - maximum tangential velocity (m/s)
fcor - Coriolis force (s^-1)
===============================================================================#

function re87(radii::AbstractVector{T},rmax::Real,r0::Real,vmax::Real,fcor::Real) where T<:Real
    return @fastmath [r < r0 ? sqrt( vmax^2 * (r / rmax)^2 * ( ( (2. * rmax) ./ (r + rmax) )^3 - ( (2. * rmax) / (r0 + rmax) )^3 ) +
                                     (fcor^2 * r^2 / 4.) ) - ( (fcor * r) / 2. ) : 0. for r in radii]
end

#==============================================================================
cw87

This function will generate the tangential velocity or relative vorticity
profile for a 1-D vortex following Chan and Williams (1987; JAS, equation 2.10)

Input (and units)
radii - vector containing radius grid points (m)
rmax - radius of maximum tangential winds (m)
vmax - maximum tangential velocity (m/s)
b - decay parameter  (unitless; positive value)
method - output 1-D tangential velocity profile (:vt) or relative vorticity
         profile (:vort)
===============================================================================#

function cw87(radii::AbstractVector{T},rmax::Real,vmax::Real,b::Real,method::Symbol=:vt) where T<:Real
    if method == :vt
        return @fastmath vmax * (radii/rmax) .* exp.( (1/b) * (1 .- (radii/rmax).^b) )
    elseif method == :vort
        return (2 * vmax) / rmax * (1 .- 0.5 * (radii/rmax).^b) .* exp.( (1/b) * (1 .- (radii/rmax).^b) )
    end
end

#==============================================================================
wc04

This function will generate the tangential velocity profile for a 1-D vortex
following Wong and Chan (2004; JAS, equation 2.2)

Input (and units)
radii - vector containing radius grid points (m)
rmax - radius of maximum tangential winds (m)
vmax - maximum tangential velocity (m/s)
b - decay parameter  (unitless; positive value)
rtc - radius of the tropical cyclone circulation (m)
===============================================================================#

function wc04(radii::AbstractVector{T},rmax::Real,vmax::Real,b::Real,rtc::Real) where T<:Real
    return [r < rtc ? vmax * (r/rmax) .* ( exp.( (1 - ((r/rmax)).^b)/b ) -
                      abs.(r - rmax) ./ (rtc-rmax) .* exp.( (1 - ((rtc/rmax)).^b)/b ) ) : 0. for r in radii]
end
