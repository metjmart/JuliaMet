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
Output = distance between the two points (given by input units for radius)
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
Output = distance between the two points (given by input units for radius)
==============================================================================#

function haversined(radius::Real,λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)
    return 2. * radius * asin(sqrt(sind((ϕ - ϕ_c) / 2) ^ 2 +
           cosd(ϕ_c) * cosd(ϕ) * sind((λ - λ_c) / 2) ^ 2))
end

#==============================================================================
invhaversine(sphere_radius::Real,radius::Real,λ_c::Real,ϕ_c::Real,θ::Real)

Invert the haversine equation to calculate a (lon, lat) destination point on the
sphere given a (lon, lat) origin point, the distance traveled, and the heading

sphere_radius = radius of sphere (units in = units out)
radius = distance traveled on sphere (same units as sphere_radius)
(λ_c, ϕ_c) = lon, lat of origin point (in radians)
theta = polar angle for the heading (note: this follows the mathematical
convention for polar angles where east = 0; north = π/2; west = π/2,
and south = 3π/2. The function will convert the polar angle to the meteorological
azimuth angle)
Output: (λ, ϕ) = Tuple with lon, lat of destination point (in radians)
==============================================================================#

function invhaversine(sphere_radius::Real,radius::Real,λ_c::Real,ϕ_c::Real,θ::Real)
    az = pi/2. - θ
    c = radius/sphere_radius
    ϕ = asin(sin(ϕ_c) * cos(c) + cos(ϕ_c) * sin(c) * cos(az))
    λ = λ_c + atan(sin(az) * sin(c) * cos(ϕ_c), cos(c) - sin(ϕ_c) * sin(ϕ))
    return λ, ϕ
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
polar_azimuth_equidist(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)

Define a function to calculate the polar angle between two points using an
azimuthal equidistant map projection

(λ_c, ϕ_c) = lon, lat of center point (in radians)
(λ, ϕ) = lon, lat of second point (in radians)
Output = polar angle between the two points (in radians)
==============================================================================#

function polar_azimuth_equidist(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)
    c = angular_dist(λ_c,ϕ_c,λ,ϕ)
    x = c * inv(sin(c)) * cos(ϕ) * sin(λ-λ_c)
    y = c * inv(sin(c)) * (cos(ϕ_c)*sin(ϕ) - sin(ϕ_c)*cos(ϕ)*cos(λ-λ_c))
    return atan(y,x)
end

#==============================================================================
polar_gc(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)

Define a function to calculate the azimuth angle between two points using the
inaccurate great circle approach

(λ_c, ϕ_c) = lon, lat of center point (in radians)
(λ, ϕ) = lon, lat of second point (in radians)
Output = polar angle between the two points
==============================================================================#

function polar_gc(λ_c::Real,ϕ_c::Real,λ::Real,ϕ::Real)
    x = (λ-λ_c)*cos(ϕ)
    y = ϕ-ϕ_c
    return atan(y,x)
end
