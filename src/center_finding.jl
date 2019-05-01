# *****************************************************************************
# center_finding.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.6.0
#
# This script contains functions for various methods of determining the center
# of a tropical cyclone.
#
# Function list:
#
# p_centroid (Adapted from Ellie Delap)
# More to come ...
# *****************************************************************************

#==============================================================================
# filtloop

This function smooths a given field (e.g., pressure or vorticity) n times with
a 1-2-1 filter. Can be used to determine the center via the location of the
smoothed mininum (in the case of pressure) or maximum (in the case of vorticity)

field
Two-dimensional field (x-y plane) - can be pressure, vorticity, geopotential, etc.

n
Number of times the filter is applied

g
Filter weights - default is a 1-2-1 filter, but can be modified
==============================================================================#

function filtloop(field::AbstractArray{T,2},n::Int;g=[0.25,0.5,0.25]) where T<:Real
    for i in 1:n
        #println("1", " ", extrema(field))
        field = filtfilt(g,field)
        #println("2", " ", extrema(field))
    end
    return field
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
            #println("==============================")
            #println("z = ",z[k]," km, k = ",k)
            #println("==============================")
            while xbar0[k] != ixbar && ybar0[k] != iybar
                ixbar = xbar0[k]
                iybar = ybar0[k]
                itr += 1
                #println("itr = ", itr)
                #println("ixbar = ", ixbar)
                #println("iybar = ", iybar)
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
                #println("new xbar0 = ", xbar0[k])
                #println("new ybar0 = ", ybar0[k])
            end
            #println("z = ", z[k], "Total itrs = ",itr)
        end
    end
    #println("---------------")
    #println(" ")
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
                      zlevs::UnitRange{Int};z=AbstractVector{zlevs},prs=AbstractArray{Tc,[:,:,zlevs]},
                      u=AbstractArray{Td,[:,:,zlevs]},v=AbstractArray{Te,[:,:,zlevs]},
                      prnt_cents::Bool=false) where {Ta<:Real,Tb<:Real,Tc<:Real,Td<:Real,Te<:Real}

    # First guess is minimum pressure centroid using R = 100 km
    # at each vertical level
    z2km = closest_ind(z,2.0)
    z2km in zlevs ? nothing : error("2-km altitude must be in the subset of altitudes")
    pp    = similar(prs)
    r,phi = xy2rp(0.,0.,x,y)
    xbar0 = similar(zlevs,Float64)
    ybar0 = similar(zlevs,Float64)
    zind = 0
    for k in zlevs
        zind += 1
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
        xbar0[zind] = pc_numx / pc_den;
        ybar0[zind] = pc_numy / pc_den;
    end
    # Now use the first guess at 2 km to interpolate data
    # and find azimuthal mean RMW
    u_rpz = regrid_xy2rp(xbar0[z2km-zlevs[1]+1],ybar0[z2km-zlevs[1]+1],x,y,r,phi,u[:,:,z2km])
    v_rpz = regrid_xy2rp(xbar0[z2km-zlevs[1]+1],ybar0[z2km-zlevs[1]+1],x,y,r,phi,v[:,:,z2km])
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
    zind =0
    for k in zlevs
        zind += 1
        ixbar = 0.
        iybar = 0.
        itr = 0
        #println("==============================")
        #println("z = ",z[k]," km, k = ",k)
        #println("==============================")
        while xbar0[zind] != ixbar && ybar0[zind] != iybar
            ixbar = xbar0[zind]
            iybar = ybar0[zind]
            itr += 1
            #println("itr = ", itr)
            #println("ixbar = ", ixbar)
            #println("iybar = ", iybar)
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
            xbar0[zind] = pc_numx / pc_den
            ybar0[zind] = pc_numy / pc_den
            #println("new xbar0 = ", xbar0[k])
            #println("new ybar0 = ", ybar0[k])
        end
        #println("z = ", z[k], "Total itrs = ",itr)
    end
    #println("---------------")
    #println(" ")
    if prnt_cents == true
        for k in zlevs
            println("z = ", z[k], ", xbar = ", xbar0[k], ", ybar = ", ybar0[k])
        end
    end
    return xbar0, ybar0
end
