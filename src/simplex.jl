#*******************************************************************************
# simplex_aux.jl
#
# Author:
#     Jonathan Martinez
#
#*******************************************************************************

"""
    init_config(xc_guess, yc_guess, rmw_guess, nguesses, dx)

Construct the initial config for the simplex algorithm.

The configuration defines nguesses^2 simplexes evenly distributed between\\
+/- 10dx of the initial center guess and 10 RMW guesses centered on the\\
initial RMW guess ranging from +/- 5dx in increments of dx.

# Arguments
- `xc_guess::Real`: x-location of initial center guess
- `yc_guess::Real`: y-location of initial center guess
- `rmw_guess::Real`: Initial RMW guess 
- `nguesses::Real`: sqrt(number of simplexes)--e.g. nguesses=5 gives 25 simplexes
- `dx::Real`: Grid spacing 
# Output
- `xinit::AbstractVector{<:Real}`: x-locations to initiate simplex searches
- `yinit::AbstractVector{<:Real}`: y-locations to initiate simplex searches
- `radii::AbstractVector{<:Real}`: Radii centered on rmwguess to initiate simplex searches
"""
function init_config(
    xc_guess::Real,
    yc_guess::Real,
    rmw_guess::Real,
    nguesses::Real,
    dx::Real)

    xinit = collect(range(xc_guess - 10. *dx, length=nguesses, stop=xc_guess + 10. *dx))
    yinit = collect(range(yc_guess - 10. *dx, length=nguesses, stop=yc_guess + 10. *dx))
    # Remove radii < 2dx:
    # Required to compute the tangential circulation within a +/- 2dx annulus
    radii = [r for r in rmw_guess-5*dx:dx:rmw_guess+5*dx if r > 2*dx]
    return xinit, yinit, radii
end

"""
    meanvta(loc, r, x, y, u, v) 

Compute the mean tangential wind from (u,v) wind components within an annulus.

# Arguments
- `loc::AbstractVector{<:Real}`: 2-element vector containing the (x,y) location of the center point
- `r::Real`: radius to center the annulus 
- `x::AbstractVector{<:Real}`: x-coordinate 
- `y::AbstractVector{<:Real}`: y-coordinate
- `u::AbstractArray{<:Real,2}`: east-west wind component (size = len(x),len(y))
- `v::AbstractArray{<:Real,2}`: north-south wind component (size = len(x),len(y))
# Output 
- `-vt_mean::Real` - minus the mean tangential wind within the annulus
"""
function meanvta(
    loc::AbstractVector{Ta},
    r::Real,
    x::AbstractVector{Tb},
    y::AbstractVector{Tc},
    u::AbstractArray{Td,2},
    v::AbstractArray{Te,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real}

    # Define center from input
    xc,yc = loc
    # Define the annulus for +/- 2*dx
    dx = x[2]-x[1]
    annulus_rinner = r - 2*dx
    annulus_router = r + 2*dx
    # Calculate mean vt within the annulus
    vt_mean = 0.
    n = 0.
    for j in eachindex(y)
        for i in eachindex(x)
            radius = sqrt((x[i] - xc)^2 + (y[j] - yc)^2)
            if (radius >= annulus_rinner) & (radius <= annulus_router)
                theta = atan(y[j] - yc, x[i] - xc)
                vt_mean += -u[i,j] * sin(theta) + v[i,j] * cos(theta)
                n += 1
            end
        end
    end
    vt_mean /= n
    # Return minus the mean vt within the annulus
    return -vt_mean
end

"""
    nanmeanvta(
        loc::AbstractVector{<:Real},
        r::Real,
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        u::AbstractArray{<:Real,2},
        v::AbstractArray{<:Real,2})

Compute the mean tangential wind from (u,v) wind components within an annulus.

Ignore NaNs in the data when calculating the annulus-mean tangential winds.\\
This approach is useful for observational (e.g., radar) data.

# Arguments
- `loc::AbstractVector{<:Real}`: 2-element vector containing the (x,y) location of the center point
- `r::Real`: radius to center the annulus 
- `x::AbstractVector{<:Real}`: x-coordinate 
- `y::AbstractVector{<:Real}`: y-coordinate
- `u::AbstractArray{<:Real,2}`: east-west wind component (size = len(x),len(y))
- `v::AbstractArray{<:Real,2}`: north-south wind component (size = len(x),len(y))
# Output 
- `-vt_mean::Real` - minus the mean tangential wind within the annulus
"""
function nanmeanvta(
    loc::AbstractVector{Ta},
    r::Real,
    x::AbstractVector{Tb},
    y::AbstractVector{Tc},
    u::AbstractArray{Td,2},
    v::AbstractArray{Te,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real}

    # Define center from input
    xc,yc = loc
    # Define the annulus for +/- 2*dx
    dx = x[2]-x[1]
    annulus_rinner = r - 2*dx
    annulus_router = r + 2*dx
    # Calculate mean vt within the annulus
    vt_mean = 0.
    n = 0.
    for j in eachindex(y)
        for i in eachindex(x)
            if (!isnan(u[i,j])) & (!isnan(v[i,j]))
                radius = sqrt((x[i] - xc)^2 + (y[j] - yc)^2)
                if (radius >= annulus_rinner) & (radius <= annulus_router)
                    theta = atan(y[j] - yc, x[i] - xc)
                    vt_mean += -u[i,j] * sin(theta) + v[i,j] * cos(theta)
                    n += 1
                end
            end
        end
    end
    vt_mean /= n
    # Return minus the mean vt within the annulus
    return -vt_mean
end

"""
    get_prelim_center(sxcen, sycen)

Compute the preliminary simplex center.

Average the subset of input simplex center points within 1σ of of the mean\\
input simplex centers.

# Arguments
- `sxcen::AbstractArray{<:Real,2}`: Array containing the x-locations of the simplex center points
- `sycen::AbstractArray{<:Real,2}`: Array containing the y-locations of the simplex center points
# Output
- `prelim_xbar::Real`: Preliminary simplex center x-location
- `prelim_ybar::Real`: Preliminary simplex center y-location
- `prelim_stdv::Real`: Standard deviation of the preliminary simplex center locations
"""
function get_prelim_center(
    sxcen::AbstractArray{Ta,Na},
    sycen::AbstractArray{Tb,Nb}) where {Ta<:Real,Tb<:Real,Na,Nb}
    #sxcen::AbstractArray{Ta,2},
    #sycen::AbstractArray{Tb,2}) where {Ta<:Real,Tb<:Real}

    xbar = mean(sxcen)
    ybar = mean(sycen)
    stdv = sqrt(var(sxcen) + var(sycen))
    dist = sqrt.((sxcen .- xbar).^2 + (sycen .- ybar).^2)
    cond = dist .< stdv
    prelim_xbar = mean(sxcen[cond])
    prelim_ybar = mean(sycen[cond])
    prelim_stdv = sqrt(var(sxcen[cond]) + var(sycen[cond]))
    return prelim_xbar, prelim_ybar, prelim_stdv
end

"""
    meanvtr(xc, yc, r, x, y, u, v)

Compute the mean tangential wind from (u,v) wind components for a given radius.

# Arguments 
- `xc::Real`: x-location of the center point
- `yc::Real`: y-location of the center point
- `r::Real`: Radius at which to calculate mean tangential wind
- `x::Abstract{<:Real}`: x-coordinate 
- `y::AbstractVector{<:Real}`: y-coordinate
- `u::AbstractArray{<:Real,2}`: east-west wind component (size = len(x),len(y))
- `v::AbstractArray{<:Real,2}`: north-south wind component (size = len(x),len(y))
# Output 
- `vt_mean::Real`: mean tangential wind at the specified radius
"""
function meanvtr(
    xc::Real,
    yc::Real,
    r::Real,
    x::AbstractVector{Ta},
    y::AbstractVector{Tb},
    u::AbstractArray{Tc,2},
    v::AbstractArray{Td,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}

    # Calculate mean vt within r +/- dx
    dx = x[2]-x[1]
    rinner = r - 0.5*dx
    router = r + 0.5*dx
    vt_mean = 0.
    n = 0.
    for j in eachindex(y)
        for i in eachindex(x)
            radius = sqrt((x[i] - xc)^2 + (y[j] - yc)^2)
            if (radius >= rinner) & (radius <= router)
                theta = atan(y[j] - yc, x[i] - xc)
                vt_mean += -u[i,j] * sin(theta) + v[i,j] * cos(theta)
                n += 1
            end
        end
    end
    vt_mean /= n
    return vt_mean
end

"""
    nanmeanvtr(xc, yc, r, x, y, u, v)

Compute the mean tangential wind from (u,v) wind components for a given radius.

Ignore NaNs in the data when calculating the mean tangential winds.\\
This approach is useful for observational (e.g., radar) data.

# Arguments 
- `xc::Real`: x-location of the center point
- `yc::Real`: y-location of the center point
- `r::Real`: Radius at which to calculate mean tangential wind
- `x::Abstract{<:Real}`: x-coordinate 
- `y::AbstractVector{<:Real}`: y-coordinate
- `u::AbstractArray{<:Real,2}`: east-west wind component (size = len(x),len(y))
- `v::AbstractArray{<:Real,2}`: north-south wind component (size = len(x),len(y))
# Output 
- `vt_mean::Real`: mean tangential wind at the specified radius
"""
function nanmeanvtr(
    xc::Real,
    yc::Real,
    r::Real,
    x::AbstractVector{Ta},
    y::AbstractVector{Tb},
    u::AbstractArray{Tc,2},
    v::AbstractArray{Td,2}) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real}

    # Calculate mean vt within r +/- dx
    dx = x[2]-x[1]
    rinner = r - 0.5*dx
    router = r + 0.5*dx
    vt_mean = 0.
    n = 0.
    for j in eachindex(y)
        for i in eachindex(x)
            if (!isnan(u[i,j])) & (!isnan(v[i,j]))
                radius = sqrt((x[i] - xc)^2 + (y[j] - yc)^2)
                if (radius >= rinner) & (radius <= router)
                    theta = atan(y[j] - yc, x[i] - xc)
                    vt_mean += -u[i,j] * sin(theta) + v[i,j] * cos(theta)
                    n += 1
                end
            end
        end
    end
    vt_mean /= n
    return vt_mean
end

"""
    objsimplex(xc_guess, yc_guess, rmw_guess, nguesses, coord, u, v)

Find the center of a tropical cyclone by maximizing the tangential circulation.

Following Lee and Marks (2000; MWR) and Bell and Lee (2012; JAMC) section 2,
determine the "optimal" center location of a tropical cyclone. The objective 
simplex center-finding algorithm will launch several simplexes surrounding an
initial center guess and maximize the tangential circulation within annuli 
centered on a set of initial guesses for the radius of maximum winds.

# Arguments 
- `xc_guess::Real`: x-location of initial center guess
- `yc_guess::Real`: y-location of initial center guess
- `rmw_guess::Real`: Initial guess for the radius of maximum winds`
- `nguesses::Real`: sqrt(number of simplexes)--e.g. nguesses=5 gives 25 simplexes
- `coord::Any`: Composite data type (struct) containing x, y, and dx
# Output
- `sxc::Real`: x-location of the objective simplex center
- `syc::Real`: y-location of the objective simplex center
- `srmw::Real`: RMW accompanying the objective simplex center
- `sstd::Real`: Standard deviation of the preliminary simplex center locations 
"""
function objsimplex(
    xc_guess::Real,
    yc_guess::Real,
    rmw_guess::Real,
    nguesses::Real,
    coord,
    u::AbstractArray{Ta,2},
    v::AbstractArray{Tb,2}) where {Ta<:Real,Tb<:Real}
    
    xinit,yinit,radii = init_config(xc_guess,yc_guess,rmw_guess,nguesses,coord.dx)
    sxc = 0.
    syc = 0.
    srmw = 0.
    sstd = 0.
    max_mean_vt = -Inf
    for ir in eachindex(radii)
        sxcen = Float64[]
        sycen = Float64[]
        for j in eachindex(yinit)
            for i in eachindex(xinit)
                soln = optimize(loc -> meanvta(loc,radii[ir],coord.x,coord.y,u,v),
                                [xinit[i], yinit[j]], NelderMead(), Optim.Options(g_tol=1e-4, iterations=30))
                push!(sxcen, soln.minimizer[1])
                push!(sycen, soln.minimizer[2])
            end
        end
        # Remove the outliers for each radius and return the prelim center
        pxc,pyc,pstd = get_prelim_center(sxcen,sycen)
        # Compute the azimuthal mean tangential wind for the prelim center
        # If mean_vt could be NaN, add a joint conditional or outer conditional statement
        mean_vt = meanvtr(pxc,pyc,radii[ir],coord.x,coord.y,u,v)
        if mean_vt > max_mean_vt 
            sxc = pxc
            syc = pyc 
            srmw = radii[ir]
            sstd = pstd
            max_mean_vt = mean_vt
        end
    end
    return sxc, syc, srmw, sstd
end

"""
    nanobjsimplex(xc_guess, yc_guess, rmw_guess, nguesses, coord, u, v)

Find the center of a tropical cyclone by maximizing the tangential circulation.

Following Lee and Marks (2000; MWR) and Bell and Lee (2012; JAMC) section 2,
determine the "optimal" center location of a tropical cyclone. The objective 
simplex center-finding algorithm will launch several simplexes surrounding an
initial center guess and maximize the tangential circulation within annuli 
centered on a set of initial guesses for the radius of maximum winds.

Ignore NaNs in the data when calculating the mean tangential winds.\\
This approach is useful for observational (e.g., radar) data.

# Arguments 
- `xc_guess::Real`: x-location of initial center guess
- `yc_guess::Real`: y-location of initial center guess
- `rmw_guess::Real`: Initial guess for the radius of maximum winds`
- `nguesses::Real`: sqrt(number of simplexes)--e.g. nguesses=5 gives 25 simplexes
- `coord::Any`: Composite data type (struct) containing x, y, and dx
# Output
- `sxc::Real`: x-location of the objective simplex center
- `syc::Real`: y-location of the objective simplex center
- `srmw::Real`: RMW accompanying the objective simplex center
- `sstd::Real`: Standard deviation of the preliminary simplex center locations 
"""
function nanobjsimplex(
    xc_guess::Real,
    yc_guess::Real,
    rmw_guess::Real,
    nguesses::Real,
    coord,
    u::AbstractArray{Ta,2},
    v::AbstractArray{Tb,2}) where {Ta<:Real,Tb<:Real}
    
    xinit,yinit,radii = init_config(xc_guess,yc_guess,rmw_guess,nguesses,coord.dx)
    sxc = 0.
    syc = 0.
    srmw = 0.
    sstd = 0.
    max_mean_vt = -Inf
    for ir in eachindex(radii)
        sxcen = Float64[]
        sycen = Float64[]
        for j in eachindex(yinit)
            for i in eachindex(xinit)
                soln = optimize(loc -> nanmeanvta(loc,radii[ir],coord.x,coord.y,u,v),
                                [xinit[i], yinit[j]], NelderMead(), Optim.Options(g_tol=1e-4, iterations=30))
                push!(sxcen, soln.minimizer[1])
                push!(sycen, soln.minimizer[2])
            end
        end
        # Remove the outliers for each radius and return the prelim center
        pxc,pyc,pstd = get_prelim_center(sxcen,sycen)
        # Compute the azimuthal mean tangential wind for the prelim center
        mean_vt = nanmeanvtr(pxc,pyc,radii[ir],coord.x,coord.y,u,v)
        if mean_vt > max_mean_vt 
            sxc = pxc
            syc = pyc 
            srmw = radii[ir]
            sstd = pstd
            max_mean_vt = mean_vt
        end
    end
    return sxc, syc, srmw, sstd
end
