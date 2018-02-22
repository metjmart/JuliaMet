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
# p_centroid

This function determines the center of a tropical cyclone from model output. 
The method follows Nguyen et al. (2014), specifically equations 2 & 3 and the 
discussion at the end of section 7 (pg 4337, last paragraph).
The script computes the pressure centroid at each vertical level in the domain
and then uses the 2-8 km layer centroid average as the center. 
The function can handle one time-step so the user must create a loop if 
the pressure centroid for each time step is desired.
==============================================================================#

function p_centroid(x::AbstractVector{<:Real},y::AbstractVector{<:Real},
                    z::AbstractVector{<:Real},prs::AbstractArray{<:Real,3},
                    u::AbstractArray{<:Real,3},v::AbstractArray{<:Real,3},
                    prnt_cents::Bool=false)
    
    # First guess is minimum pressure centroid using R = 100 km 
    # at each vertical level
    xx,yy = grid2d(x,y)
    xbar0 = Array{Float64}(length(z))
    ybar0 = Array{Float64}(length(z))
    pp    = similar(prs)
    z2km  = closest_ind(z,2.0) # Only find centroids for 
    z8km  = closest_ind(z,8.0) # the 2-8 km layer
    for k in eachindex(z)
        if k >= z2km && k <= z8km
            min_ind = findmin(prs[:,:,k])[2]
            pminx = x[findin(x,xx[min_ind])[1]]
            pminy = y[findin(y,yy[min_ind])[1]]
            # Compute pp as p - p_env
            r,phi,prs_rpz = regrid_xy2rt(pminx,pminy,x,y,prs[:,:,k],true)
            azmean_p = azmean(r,z,prs_rpz)
            # p_env is azimuthally averaged p at r = 500 km
            r500 = closest_ind(r,500)
            p_env = azmean_p[r500]
            pp[:,:,k] = prs[:,:,k] - p_env
            r0 = 100
            circle_pp = []
            circle_x  = []
            circle_y  = []
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
    r,phi,u_rpz = regrid_xy2rt(xbar0[z2km],ybar0[z2km],x,y,u[:,:,z2km],true)
    v_rpz = regrid_xy2rt(xbar0[z2km],ybar0[z2km],x,y,v[:,:,z2km])
    # Convert u,v at 2 km to ur,vt
    ur,vt = uv2urvt(phi,u_rpz,v_rpz)
    # Compute the azimuthal mean tangential wind
    azmean_vt = azmean(r,z,vt)
    # Compute the RMW and 2R80
    vtmax = findmax(azmean_vt)
    rmw = r[vtmax[2]]
    r80 = closest_ind(azmean_vt[1:vtmax[2]],0.8*vtmax[1])
    two_r80 = 2.0 * r[r80]
    # Now use 2r80 as the radius threshold
    for k in eachindex(z)
        if k >= z2km && k <= z8km 
            ixbar = 0
            iybar = 0
            itr = 0
            #println("==============================")
            #println("z = ",z[k]," km, k = ",k)
            #println("==============================")
            while xbar0[k] != ixbar && ybar0[k] != iybar
                ixbar = xbar0[k]
                iybar = ybar0[k]
                itr = itr + 1
                #println("itr = ", itr)
                #println("ixbar = ", ixbar)
                #println("iybar = ", iybar)
                circle_pp = []
                circle_x  = []
                circle_y  = []
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
       
