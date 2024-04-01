# *****************************************************************************
# center_finding.jl

# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# This script contains functions for various methods of determining the center
# of a tropical cyclone.
#
# Function list
# filtloop
# psi_centroid_xy
# p_centroid (Adapted from Ellie Delap)
# p_centroid_
# *****************************************************************************

#==============================================================================
# filtloop

Smooth a given field n times with a 1-2-1 filter (by default)

filtloop(field::AbstractVector{T},n::Int;g=[0.25,0.5,0.25]) where T<:Real

Input
field - One-dimensional field along the axis to be smoothed 
n - Number of times the filter is applied
g - Filter weights - default is a 1-2-1 filter, but can be modified
Output
field - Smoothed field 
==============================================================================#

function filtloop(field::AbstractVector{T},n::Int;g=[0.25,0.5,0.25]) where T<:Real
    for i in 1:n
        field = filtfilt(g,field)
    end
    return field
end

#==============================================================================
psi_centroid_xy

Generic function to calculate the centroid of a specified variable in Cartesian
(x, y) coordinates

psi_centroid_xy(xc::Real,yc::Real,r0::Real,x::AbstractVector{Ta},y::AbstractVector{Tb},
                         psi::AbstractArray{Tc,2}) 

Input
Note - xc, yc, R0, x, and y should all have the same units!
xc - x-coordinate for initial center guess 
yc - y-coordinate for initial center guess 
r0 - radius of centroid 
x - x-coordinate vector 
y - y-coordinate vector 
Output
psi - 2-D (x, y) variable for computing centroid
==============================================================================#

function psi_centroid_xy(xc::Real,yc::Real,r0::Real,x::AbstractVector{Ta},y::AbstractVector{Tb},
                         psi::AbstractArray{Tc,2}) where {Ta<:Real,Tb<:Real,Tc<:Real}

    psi_x = 0.
    psi_y = 0.
    denom = 0.
    for j in eachindex(y)
        for i in eachindex(x)
            R = sqrt((x[i]-xc)^2 + (y[j]-yc)^2)
            if R <= r0
                psi_x += psi[i,j] * x[i]
                psi_y += psi[i,j] * y[j]
                denom += psi[i,j]
            end
        end
    end
    x_bar = psi_x/denom
    y_bar = psi_y/denom
    i_bar = 0.
    j_bar = 0.
    itr = 0
    while !isapprox(x_bar,i_bar,atol=1e-4) && !isapprox(y_bar,j_bar,atol=1e-4)
        i_bar = x_bar
        j_bar = y_bar
        itr += 1
        psi_x = 0.
        psi_y = 0.
        denom = 0. 
        #nr = 0
        for j in eachindex(y)
            for i in eachindex(x)
                R = sqrt((x[i]-i_bar)^2 + (y[j]-j_bar)^2)
                if R <= r0
                    psi_x += psi[i,j] * x[i]
                    psi_y += psi[i,j] * y[j]
                    denom += psi[i,j]
                    #nr += 1
                end
            end
        end
        x_bar = psi_x/denom
        y_bar = psi_y/denom
        itr > 30 && break
    end
    return x_bar, y_bar
end

#==============================================================================
# p_centroid

This function determines the center of a tropical cyclone from model output.
The method follows Nguyen et al. (2014), specifically equations 2 & 3 and the
discussion at the end of section 7 (pg 4337, last paragraph).
The script computes the pressure centroid at each vertical level in the domain
and then uses the 2-8 km layer centroid average as the center.
The function can handle one time-step so the user must create a loop if
the pressure centroid for each time step is desired.
==============================================================================#

function p_centroid(x::AbstractVector{Ta},y::AbstractVector{Tb},
                    z::AbstractVector{Tc},prs::AbstractArray{Td,3},
                    u::AbstractArray{Te,3},v::AbstractArray{Tf,3},
                    prnt_cents::Bool=false) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real,Tf<:Real}

    # First guess is minimum pressure centroid using R = 100 km
    # at each vertical level
    xbar0 = similar(z,Float64)
    ybar0 = similar(z,Float64)
    pp    = similar(prs)
    z2km  = closest_ind(z,2.0) # Only find centroids for
    z8km  = closest_ind(z,8.0) # the 2-8 km layer
    r,phi = xy2rp(0.,0.,x,y)
    for k in eachindex(z)
        if k >= z2km && k <= z8km
            min_ind = findmin(prs[:,:,k])[2]
            pminx = x[min_ind[1]]
            pminy = y[min_ind[2]]
            # Compute pp as p - p_env
            prs_rpz = regrid_xy2rp(pminx,pminy,x,y,r,phi,prs[:,:,k])
            azmean_p = nanmean(prs_rpz,2,true)
            # p_env is azimuthally averaged p at r = 500 km
            r500 = closest_ind(r,500)
            p_env = azmean_p[r500]
            pp[:,:,k] = prs[:,:,k] .- p_env
            r0 = 100.
            circle_pp = Float64[]
            circle_x  = Float64[]
            circle_y  = Float64[]
            for j in eachindex(y)
                for i in eachindex(x)
                    R = sqrt( (x[i] - pminx).^2 + (y[j] - pminy).^2)
                    if R <= r0
                        push!(circle_pp,pp[i,j,k])
                        push!(circle_x,x[i])
                        push!(circle_y,y[j])
                    end
                end
            end
            pc_den = sum(circle_pp)
            pc_numx = sum(circle_pp .* circle_x)
            pc_numy = sum(circle_pp .* circle_y)
            xbar0[k] = pc_numx / pc_den;
            ybar0[k] = pc_numy / pc_den;
        end
    end
    # Now use the first guess at 2 km to interpolate data
    # and find azimuthal mean RMW
    u_rpz = regrid_xy2rp(xbar0[z2km],ybar0[z2km],x,y,r,phi,u[:,:,z2km])
    v_rpz = regrid_xy2rp(xbar0[z2km],ybar0[z2km],x,y,r,phi,v[:,:,z2km])
    # Convert u,v at 2 km to ur,vt
    ur,vt = uv2urvt(phi,u_rpz,v_rpz)
    # Compute the azimuthal mean tangential wind
    azmean_vt = nanmean(vt,2,true)
    # Compute the RMW and 2R80
    vtmax = nanfindmax(azmean_vt)
    rmw = r[vtmax[2]]
    r80 = closest_ind(azmean_vt[1:vtmax[2]],0.8*vtmax[1])
    two_r80 = 2.0 * r[r80]
    # Now use 2r80 as the radius threshold
    for k in eachindex(z)
        if k >= z2km && k <= z8km
            ixbar = 0.
            iybar = 0.
            itr = 0
            while xbar0[k] != ixbar && ybar0[k] != iybar
                ixbar = xbar0[k]
                iybar = ybar0[k]
                itr += 1
                circle_pp = Float64[]
                circle_x  = Float64[]
                circle_y  = Float64[]
                for j in eachindex(y)
                    for i in eachindex(x)
                        R = sqrt( (x[i]-ixbar).^2 + (y[j]-iybar).^2)
                        if R <= two_r80
                            push!(circle_pp,pp[i,j,k])
                            push!(circle_x,x[i])
                            push!(circle_y,y[j])
                        end
                    end
                end
                pc_den = sum(circle_pp)
                pc_numx = sum(circle_pp .* circle_x)
                pc_numy = sum(circle_pp .* circle_y)
                xbar0[k] = pc_numx / pc_den
                ybar0[k] = pc_numy / pc_den
           end
        end
    end
   if prnt_cents == true
        for k in eachindex(z)
            if k >= z2km && k <= z8km
                println("z = ", z[k], ", xbar = ", xbar0[k], ", ybar = ", ybar0[k])
            end
        end
    end
    # Now take the 2-8 km layer average for the center
    pcx = mean(xbar0[z2km:z8km])
    pcy = mean(ybar0[z2km:z8km])
    return pcx, pcy
end

#==============================================================================
# p_centroid_

As in p_centroid, but compute pressure centroids for range of specified
vertical levels
Use a range of indices such as 2:8 for zlevs, where 2:8 should correspond
to indices in the altitude array
==============================================================================#

function p_centroid_(x::AbstractVector{Ta},y::AbstractVector{Tb},
                      zlevs::UnitRange{Int},z::AbstractVector{Tc},prs::AbstractArray{Td},
                      u::AbstractArray{Te},v::AbstractArray{Tf},
                      prnt_cents::Bool=false) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real,Tf<:Real}

    # First guess is minimum pressure centroid using R = 100 km
    # at each vertical level
    z2km = closest_ind(z,2.0)
    z2km in zlevs ? nothing : error("2-km altitude must be in the range of zlevs")
    size(prs)[3] == size(u)[3] == size(v)[3] == length(zlevs) ? nothing :
        error("Only input prs, u, and v arrays for the range given by zlevs")
    pp    = similar(prs)
    r,phi = xy2rp(0.,0.,x,y)
    xbar0 = similar(zlevs,Float64)
    ybar0 = similar(zlevs,Float64)
    for k in eachindex(zlevs)
        min_ind = findmin(prs[:,:,k])[2]
        pminx = x[min_ind[1]]
        pminy = y[min_ind[2]]
        # Compute pp as p - p_env
        prs_rpz = regrid_xy2rp(pminx,pminy,x,y,r,phi,prs[:,:,k])
        azmean_p = nanmean(prs_rpz,2,true)
        # p_env is azimuthally averaged p at r = 500 km
        r500 = closest_ind(r,500)
        p_env = azmean_p[r500]
        pp[:,:,k] = prs[:,:,k] .- p_env
        # Use R = 100 km as the first guess inner-core region
        r0 = 100.
        circle_pp = Float64[]
        circle_x  = Float64[]
        circle_y  = Float64[]
        for j in eachindex(y)
            for i in eachindex(x)
                R = sqrt( (x[i] - pminx).^2 + (y[j] - pminy).^2)
                if R <= r0
                    push!(circle_pp,pp[i,j,k])
                    push!(circle_x,x[i])
                    push!(circle_y,y[j])
                end
            end
        end
        pc_den = sum(circle_pp)
        pc_numx = sum(circle_pp .* circle_x)
        pc_numy = sum(circle_pp .* circle_y)
        xbar0[k] = pc_numx / pc_den;
        ybar0[k] = pc_numy / pc_den;
    end
    # Now use the first guess at 2 km to interpolate data
    # and find azimuthal mean RMW
    z2km_zlev = findall(x->x==z2km,zlevs)[1]
    u_rpz = regrid_xy2rp(xbar0[z2km_zlev],ybar0[z2km_zlev],x,y,r,phi,u[:,:,z2km_zlev])
    v_rpz = regrid_xy2rp(xbar0[z2km_zlev],ybar0[z2km_zlev],x,y,r,phi,v[:,:,z2km_zlev])
    # Convert u,v at 2 km to ur,vt
    ur,vt = uv2urvt(phi,u_rpz,v_rpz)
    # Compute the azimuthal mean tangential wind
    azmean_vt = nanmean(vt,2,true)
    # Compute the RMW and 2R80
    vtmax = nanfindmax(azmean_vt)
    rmw = r[vtmax[2]]
    r80 = closest_ind(azmean_vt[1:vtmax[2]],0.8*vtmax[1])
    two_r80 = 2.0 * r[r80]
    # Now use 2r80 as the radius threshold
    for k in eachindex(zlevs)
        ixbar = 0.
        iybar = 0.
        itr = 0
        while xbar0[k] != ixbar && ybar0[k] != iybar
            ixbar = xbar0[k]
            iybar = ybar0[k]
            itr += 1
            circle_pp = Float64[]
            circle_x  = Float64[]
            circle_y  = Float64[]
            for j in eachindex(y)
                for i in eachindex(x)
                    R = sqrt( (x[i]-ixbar).^2 + (y[j]-iybar).^2)
                    if R <= two_r80
                        push!(circle_pp,pp[i,j,k])
                        push!(circle_x,x[i])
                        push!(circle_y,y[j])
                    end
                end
            end
            pc_den = sum(circle_pp)
            pc_numx = sum(circle_pp .* circle_x)
            pc_numy = sum(circle_pp .* circle_y)
            xbar0[k] = pc_numx / pc_den
            ybar0[k] = pc_numy / pc_den
        end
    end
    if prnt_cents == true
        for k in zlevs
            println("z = ", z[k], ", xbar = ", xbar0[k], ", ybar = ", ybar0[k])
        end
    end
    return xbar0, ybar0
end
