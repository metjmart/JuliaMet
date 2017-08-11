# *****************************************************************************
# cyclone.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.4.5
#
# This script contains functions that are more or less specific to tropical 
# cyclone related applications.
# 
# calc_rmw
# calc_azmean
# uv2urvt
# p3swploc
# *****************************************************************************

#==============================================================================
calc_rmw

Determine the azimuthal mean radius of maximum tangential winds at each vertical 
level. Assumes that the input vt is the azimuthally averaged vt and has 
dimensions of r,z.
==============================================================================#

function calc_rmw{Ta<:Real,Tb<:Real,Tc<:Real}(r::AbstractVector{Ta},
                                              z::AbstractVector{Tb}, 
                                              azmean_vt::AbstractArray{Tc,2})
    # Create an array with the radial values at each height
    r_new = Array(Float64,length(r),length(z))
    for i in eachindex(z)
        r_new[:,i] = r
    end
    # Compute the azimuthal mean RMW at each vertical level
    rmw = similar(z,Float64)
    fill!(rmw, NaN)
    vtmax = findmax(azmean_vt,1)
    for k in eachindex(z)
        # Only store the RMW values if they're present
        if vtmax[2][k] != 0
            rmw[k] = r_new[vtmax[2][k]]
        end
    end

    return rmw

end

#==============================================================================
calc_azmean

Compute the azimuthal mean of the input variable. 
** Be sure to use the nanmean function since NaNs may be present
** var needs to be a three dimensional field that's indexed r,theta,z
==============================================================================#

function calc_azmean{Ta<:Real,Tb<:Real,Tc<:Real}(r::AbstractVector{Ta},
                                                 z::AbstractVector{Tb},
                                                 var::AbstractArray{Tc,3})
    # Define an array for the azimuthal mean
    azmean_var = Array(Float64,length(r),length(z)) 
    fill!(azmean_var, NaN)
    # Calculate the azimuthal mean 
    for i in eachindex(r)
        for j in eachindex(z)
            azmean_var[i,j] = nanmean(var[i,:,j])
        end
    end

    return azmean_var

end

#==============================================================================
uv2urvt

Convert u (east-west) and v (north-south) winds to radial and tangential winds
** Assumes u and v input winds are on an RTZ grid (radius,theta,height).
** Input theta vector assumed to be in units of degrees
==============================================================================#

function uv2urvt{Ta<:Real,Tb<:Real,Tc<:Real}(theta::AbstractVector{Ta},
                                             u::AbstractArray{Tb,3},
                                             v::AbstractArray{Tc,3})
    # Create the ur and vt arrays
    ur = Array(Float64,size(u))
    vt = Array(Float64,size(v))
    # Convert u,v to ur,vt
    for i in eachindex(theta)
        ur[:,i,:] =  u[:,i,:] .* cos(deg2rad(theta[i])) + 
                     v[:,i,:] .* sin(deg2rad(theta[i]))
        vt[:,i,:] = -u[:,i,:] .* sin(deg2rad(theta[i])) + 
                     v[:,i,:] .* cos(deg2rad(theta[i]))
    end
 
    return ur,vt

end

#==============================================================================
p3swploc

Determine the spatial location of the fore and aft radar beams based on the 
location of the P3. The function requires the 1-sec flight level data from 
the HRD site (for example, 
ftp://ftp.aoml.noaa.gov/hrd/pub/data/flightlevel/2015/patricia/20151022I1.1sec.txt)
** It is assumed that the header of the 1-sec flight level data text file has 
already been removed.
** The extent of the fore/aft beams is assumed to be 67 km (the range of data
retained when using the medium-threshold editing script)

==================================
Example of how to use the function
==================================

    filein = "/Users/jmart23/TCI/Patricia/noaa_flight_level/20151023I1/20151023I1.sec";
    
    fl_time,fl_lon,fl_lat,foreR_lon,foreR_lat,foreL_lon,foreL_lat,
    aftR_lon,aftR_lat,aftL_lon,aftL_lat = p3swploc(filein);

    # Specify a time for the P3 location

    p3time = 173839.0
    p3ind = findin(fl_time,p3time)[1]

    p3_lon = fl_lon[p3ind]
    p3_lat = fl_lat[p3ind]

    # Plot the analysis and overlay the swp locations for the specified time
    
    plt.contourf(lon,lat,v[:,:,9]',levels=collect(-80:10:80),
                 cmap="radar_carbone2",extend="both")
    plt.colorbar()
    
    plt.plot(fl_lon[p3ind],fl_lat[p3ind],color="k",markersize=4)
    plt.plot([fl_lon[p3ind],foreR_lon[p3ind]],[fl_lat[p3ind],foreR_lat[p3ind]],
             color="C0",linewidth=2)
    plt.plot([fl_lon[p3ind],foreL_lon[p3ind]],[fl_lat[p3ind],foreL_lat[p3ind]],
             color="C0",linewidth=2)
    plt.plot([fl_lon[p3ind],aftR_lon[p3ind]],[fl_lat[p3ind],aftR_lat[p3ind]],
             color="C1",linewidth=2)
    plt.plot([fl_lon[p3ind],aftL_lon[p3ind]],[fl_lat[p3ind],aftL_lat[p3ind]]
             ,color="C1",linewidth=2)
    plt.scatter(fl_lon[p3ind],fl_lat[p3ind],color="k",linewidth=4,zorder=10)

==============================================================================#

function p3swploc(filein::AbstractString)

    # Read in the flight level data 
    fl_data = readdlm(filein)
    # Read in the required variables
    fl_time = Array{Float64}(fl_data[1:end-1,1])
    fl_lat  = Array{Float64}(fl_data[1:end-1,2])
    fl_lon  = Array{Float64}(fl_data[1:end-1,3])
    heading = Array{Float64}(fl_data[1:end-1,4]) 
    # Specify the angles of the fore/aft beams given the heading of the plane
    foreR = heading + 70.0
    foreL = heading - 70.0
    aftR  = heading + 110.0
    aftL  = heading - 110.0
    # Constrain the range of azimuth values to 0-360 degrees
    for i in eachindex(foreR)
        if foreR[i] > 360.0
            foreR[i] -= 360.0
        end
        if foreL[i] < 0.0
            foreL[i] += 360.0
        end
        if aftR[i] > 360.0
            aftR[i] -= 360.0
        end
        if aftL[i] < 0.0
            aftL[i] += 360.0
        end
    end 
    # Convert the azimuth angles to math degrees
    foreRdeg = (foreR - 90.0) .* -1.0
    foreLdeg = (foreL - 90.0) .* -1.0
    aftRdeg  = (aftR - 90.0) .* -1.0
    aftLdeg  = (aftL - 90.0) .* -1.0
    # Determine the x- and y-distances from the location of aircraft
    # out to 67 km (range of data in sweep files)
    foreR_xdist = 67.0 .* cos(deg2rad(foreRdeg))
    foreR_ydist = 67.0 .* sin(deg2rad(foreRdeg)) 
    foreL_xdist = 67.0 .* cos(deg2rad(foreLdeg))
    foreL_ydist = 67.0 .* sin(deg2rad(foreLdeg)) 
    aftR_xdist = 67.0 .* cos(deg2rad(aftRdeg))
    aftR_ydist = 67.0 .* sin(deg2rad(aftRdeg)) 
    aftL_xdist = 67.0 .* cos(deg2rad(aftLdeg))
    aftL_ydist = 67.0 .* sin(deg2rad(aftLdeg))
    # Determine the distance of 1 deg lat and lon at the aircraft location
    # ** In units of km
    fl_latrad = deg2rad(fl_lat);
    fac_lat = 111.13209 - 0.56605 .* cos(2.0 .* fl_latrad)
        + 0.00012 .* cos(4.0 .* fl_latrad) - 0.000002 .* cos(6.0 .* fl_latrad)
    fac_lon = 111.41513 .* cos(fl_latrad)
        - 0.09455 .* cos(3.0 .* fl_latrad) + 0.00012 .* cos(5.0 .* fl_latrad)
    # Take the x- and y-distances and convert them to a lat and lon 
    # position relative to the aircraft location
    foreR_lon = (foreR_xdist ./ fac_lon) + fl_lon
    foreR_lat = (foreR_ydist ./ fac_lat) + fl_lat
    foreL_lon = (foreL_xdist ./ fac_lon) + fl_lon
    foreL_lat = (foreL_ydist ./ fac_lat) + fl_lat
    aftR_lon = (aftR_xdist ./ fac_lon) + fl_lon
    aftR_lat = (aftR_ydist ./ fac_lat) + fl_lat
    aftL_lon = (aftL_xdist ./ fac_lon) + fl_lon
    aftL_lat = (aftL_ydist ./ fac_lat) + fl_lat

    return fl_time,fl_lon,fl_lat,foreR_lon,foreR_lat,foreL_lon,foreL_lat,
           aftR_lon,aftR_lat,aftL_lon,aftL_lat

end

