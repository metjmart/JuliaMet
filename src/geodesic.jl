# *****************************************************************************
# geodesic.jl
#
# Functions to calculate the distance between two points on a sphere
# *****************************************************************************

#==============================================================================
haversine(radius::Real,λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)

Define the haversine function to calculate the great circle distance between two
points on a sphere
See Snyder (1987), pg. 30

(λ_c, ϕ_c) = lon, lat of center point (in radians)
(λ, ϕ) = lon, lat of second point (in radians)
radius = radius of sphere (units in = units out)
==============================================================================#

function haversine(radius::Real,λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)
    return 2. * radius * asin(sqrt(sin((ϕ - ϕ_c) / 2) ^ 2 +
           cos(ϕ_c) * cos(ϕ) * sin((λ - λ_c) / 2) ^ 2))
end

#==============================================================================
haversined(radius::Real,λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)

Same as haversine(), but input coordinates in degrees

Define the haversine function to calculate the great circle distance between two
points on a sphere
See Snyder (1987), pg. 30

(λ_c, ϕ_c) = lon, lat of center point (in degrees)
(λ, ϕ) = lon, lat of second point (in degrees)
radius = radius of sphere (units in = units out)
==============================================================================#

function haversined(radius::Real,λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)
    return 2. * radius * asin(sqrt(sind((ϕ - ϕ_c) / 2) ^ 2 +
           cosd(ϕ_c) * cosd(ϕ) * sind((λ - λ_c) / 2) ^ 2))
end

#==============================================================================
angular_dist(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)

Define a function to calculate the angular distance between two points using
the Law of Cosines
See Snyder (1987), pg. 30

Note - This is just the haversine formula absent the radius factor

(λ_c, ϕ_c) = lon, lat of center point (in radians)
(λ, ϕ) = lon, lat of second point (in radians)
Output = angular distance (in radians)
==============================================================================#

function angular_dist(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)
    return 2. * asin(sqrt(sin((ϕ - ϕ_c) / 2) ^ 2 +
           cos(ϕ_c) * cos(ϕ) * sin((λ - λ_c) / 2) ^ 2))
end

#==============================================================================
atan_full(ang::Real)

Define a function to return atan values on the full unit circle (0, 2π)

ang = Input azimuth angle to be converted (in radians)
==============================================================================#

function atan_full(ang::Real)
    return ang < 0. ? ang += 2. * pi : ang
end

#==============================================================================
azimuth_equidist(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)

Define a function to calculate the azimuth angle between two points using an
azimuthal equidistant map projection
Return azimuth angles on the full unit circle (0, 2π)

(λ_c, ϕ_c) = lon, lat of center point (in radians)
(λ, ϕ) = lon, lat of second point (in radians)
==============================================================================#

function azimuth_equidist(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)
    c = angular_dist(λ_c,ϕ_c,λ,ϕ)
    x = c * inv(sin(c)) * cos(ϕ) * sin(λ-λ_c)
    y = c * inv(sin(c)) * (cos(ϕ_c)*sin(ϕ) - sin(ϕ_c)*cos(ϕ)*cos(λ-λ_c))
    return atan_full(atan(y,x))
end

#==============================================================================
gc_azimuth(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)

Define a function to calculate the azimuth angle between two points using the
inaccurate great circle approach
Return azimuth angles on the full unit circle (0, 2π)

(λ_c, ϕ_c) = lon, lat of center point (in radians)
(λ, ϕ) = lon, lat of second point (in radians)
==============================================================================#

function gc_azimuth(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)
    x = (λ-λ_c)*cos(ϕ)
    y = ϕ-ϕ_c
    return atan_full(atan(y,x))
end
