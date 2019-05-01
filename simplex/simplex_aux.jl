#*******************************************************************************
# simplex.jl
#
# This script will house the necessary functions to run the objective center
# finding algorithm described by Bell and Lee (2012)
#*******************************************************************************

#===============================================================================
init_config

Initialize configuration for center finding algorithm given an initial center
guess position, an rmw guess, and the number of center guesses to accompany
the primary guess (3x3, 4x4, or 5x5 in Bell and Lee 2012; can technically
search as many as you want but more grid points = more expensive).
E.g., nguesses = 3 for 3x3, 4 for 4x4, 5 for 5x5
Range of radii is set to +/- 5 km of RMW guess in increments of 1 km
===============================================================================#

function init_config(xguess::Real,yguess::Real,rmwguess::Real,nguesses::Real)
    xinit = collect(range(xguess - 10., length=nguesses, stop=xguess + 10.))
    yinit = collect(range(yguess - 10., length=nguesses, stop=yguess + 10.))
    radii = collect(rmwguess - 5.:1.:rmwguess + 5.)
    return xinit, yinit, radii
end

#===============================================================================
nanmeanvt_annulus

Compute the mean tangential wind within an annulus of +/- 2 km from a specified
radius
===============================================================================#

function nanmeanvt_annulus(loc::AbstractVector{Ta},r::Real,u_glob::AbstractArray{Tb,2},
                           v_glob::AbstractArray{Tc,2}) where {Ta<:Real,Tb<:Real,Tc<:Real}
    # Define center from input
    xc = loc[1]
    yc = loc[2]
    # Compute the tangential wind
    vt = -u_glob .* sin.(atan.(ones(length(x)) .* (y .- yc)',(x .- xc) .* ones(length(y))')) +
          v_glob .* cos.(atan.(ones(length(x)) .* (y .- yc)',(x .- xc) .* ones(length(y))'))
    # Create an accessible interpolation object via bi-linear interpolation
    vt_itp = extrapolate(interpolate((x .- xc,y .- yc), vt, Gridded(Linear())), NaN)
    phi = collect(0:pi/180.:2*pi - pi/180.)
    # Define the range of radii for the annulus as +/- 2 km in 1 km intervals
    rrange = collect(r - 2.:1.:r + 2.)
    vt_rings = Array{Float64}(undef,length(rrange),length(phi))
    # Save vt for rings and also compute mean vt for rings
    for j in eachindex(phi)
        for i in eachindex(rrange)
            @inbounds vt_rings[i,j] = vt_itp(rrange[i] * cos(phi[j]), rrange[i] * sin(phi[j]))
        end
    end
    # Return the negative of the mean vt within the annulus
    vt_bar_annulus = nanmean(vt_rings)
    return -vt_bar_annulus
end

#===============================================================================
rm_outliers

Remove outliers from the simplex centers
===============================================================================#

function rm_outliers(xinit::AbstractVector{Ta},yinit::AbstractVector{Tb},
                     simplex_xcenters::AbstractArray{Tc,2},
                     simplex_ycenters::AbstractArray{Td,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}
    conv_xcenters = Float64[]
    conv_ycenters = Float64[]
    xbar_centers = mean(simplex_xcenters)
    ybar_centers = mean(simplex_ycenters)
    stdv_centers = sqrt(var(simplex_xcenters) + var(simplex_ycenters))
        for j in eachindex(yinit)
            for i in eachindex(xinit)
                @inbounds dist = sqrt( (simplex_xcenters[i,j] - xbar_centers)^2 +
                                       (simplex_ycenters[i,j] - ybar_centers)^2 )
                if dist < stdv_centers
                    push!(conv_xcenters,simplex_xcenters[i,j])
                    push!(conv_ycenters,simplex_ycenters[i,j])
                end
            end
        end
    prelim_xbar_centers = mean(conv_xcenters)
    prelim_ybar_centers = mean(conv_ycenters)
    prelim_stdv_centers = sqrt(var(conv_xcenters) + var(conv_ycenters))
    return prelim_xbar_centers, prelim_ybar_centers, prelim_stdv_centers
end

#===============================================================================
nanmeanvt_ring

Compute the mean tangential wind at a specified radius
===============================================================================#

# Define a function to determine the center which maximizes vt

function nanmeanvt_ring(x::AbstractVector{Ta},y::AbstractVector{Tb},prelim_xc::Real,
                        prelim_yc::Real,radius::Real,u::AbstractArray{Tc,2},
                        v::AbstractArray{Td,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}
    # Compute the tangential wind
    vt = -u .* sin.(atan.(ones(length(x)) .* (y .- prelim_yc)',(x .- prelim_xc) .* ones(length(y))')) +
          v .* cos.(atan.(ones(length(x)) .* (y .- prelim_yc)',(x .- prelim_xc) .* ones(length(y))'))
    # Create an accessible interpolation object via bi-linear interpolation
    vt_itp = extrapolate(interpolate((x .- prelim_xc,y .- prelim_yc), vt, Gridded(Linear())), NaN)
    phi = collect(0:pi/180.:2*pi - pi/180.)
    vt_ring = Array{Float64}(undef,length(phi))
    # Loop over all phi and store vt at given radius
    for j in eachindex(phi)
        @inbounds vt_ring[j] = vt_itp(radius * cos(phi[j]), radius * sin(phi[j]))
    end
    return nanmean(vt_ring)
end
