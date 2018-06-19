# *****************************************************************************
# balanced_vortex.jl
# 
# This script will compute the balanced state of a vortex following 
# Smith (2006; Tellus). 

# Description of variables to be modified by user
# 
# fin, fout 
# Input file name (fin), output filename (fout) if user specifes
# writeout = true (see below)
# 
# mode 
# "hydrostatic" - create a hydrostatic (horizontally homogeneous) reference 
#                 state based on ambient profile
# "vortex_cart" - compute the balanced state (gradient wind and hydrostatic) 
#                 of a vortex initially on a Cartesian coordinate grid
# "vortex_cyln" - compute the balanced state (gradient wind and hydrostatic) 
#                 of a vortex initially on a cylindrical coordinate grid
# 
# cenfile, cx, cy 
# Specify the center of the storm for a gridded Cartesian analysis
# Only applies to mode = "vortex_cart" 
# If center is at (0,0) manually override and make cx = 0, cy = 0
#
# paramBL
# Parameterize the boundary layer when computing the balanced vortex? 
# Only applies to mode = "vortex_cart" and mode = "vortex_cylnd"
# Note: Current parameterization simply sets winds below 2 km altitude to
#       the value at 2 km altitude
#
# f 
# Specify the Coriolis parameter for the analysis
#
# writeout
# Option to write output to NetCDF file in the format required for running 
# the thermodynamic retrieval
# *****************************************************************************

using JuliaMet, Interpolations, NetCDF, DataStructures
include("./balanced_aux.jl")

#=========================== BEGIN MODIFY ====================================#
# Specify the input file (currently, script formulated specific to SAMURAI)
    fin = "./samurai_XYZ_analysis.nc"
# Choose which mode to run
    mode = "vortex_cart"
    paramBL = true
# Specify the center location
    cx = 0.
    cy = 0.
# Specify NetCDF output file name in format required by thermodynamic retrieval
    fout = "./20151022_test_1800.nc"
# Specify the Coriolis parameter for case at hand
    f = 0.0000451 # Coriolis
# Write output to NetCDF file?
    writeout = false
#=========================== END MODIFY ======================================#

#==============================================================================
Read in all required variables for computing the balanced state
==============================================================================#

    balvarnames = ["longitude","latitude","x","y","altitude","DBZ","U","V","T","THETA","QV","RHOA"]
    gradvarnames = ["W","DUDX","DVDX","DWDX","DUDY","DVDY","DWDY","DUDZ","DVDZ","DWDZ"]

    lon,lat,x,y,z,dbz,u,v,temp,theta,qv,rhoa = read_ncvars(fin,balvarnames)
    w,dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz = read_ncvars(fin,gradvarnames,false)

    # Specify an offset based on length of x-array (avoid boundaries)
    offset = convert(Int64,floor(0.1*length(x))) 

    xs = x[offset+1]
    xe = x[end-offset]
    ys = y[offset+1]
    ye = y[end-offset]

    xr = xe - xs + 1
    yr = ye - ys + 1

    lon = lon[offset+1:end-offset]
    lat = lat[offset+1:end-offset]
    x = x[offset+1:end-offset]
    y = y[offset+1:end-offset]
    dbz = dbz[offset+1:end-offset,offset+1:end-offset,:]
    u = u[offset+1:end-offset,offset+1:end-offset,:]
    v = v[offset+1:end-offset,offset+1:end-offset,:]
    w = w[offset+1:end-offset,offset+1:end-offset,:]
    temp = temp[offset+1:end-offset,offset+1:end-offset,:]
    theta = theta[offset+1:end-offset,offset+1:end-offset,:]
    qv = qv[offset+1:end-offset,offset+1:end-offset,:]
    rhoa = rhoa[offset+1:end-offset,offset+1:end-offset,:]
    dudx = dudx[offset+1:end-offset,offset+1:end-offset,:]
    dvdx = dvdx[offset+1:end-offset,offset+1:end-offset,:]
    dwdx = dwdx[offset+1:end-offset,offset+1:end-offset,:]
    dudy = dudy[offset+1:end-offset,offset+1:end-offset,:]
    dvdy = dvdy[offset+1:end-offset,offset+1:end-offset,:]
    dwdy = dwdy[offset+1:end-offset,offset+1:end-offset,:]
    dudz = dudz[offset+1:end-offset,offset+1:end-offset,:]
    dvdz = dvdz[offset+1:end-offset,offset+1:end-offset,:]
    dwdz = dwdz[offset+1:end-offset,offset+1:end-offset,:]

    # Define constants
    g = 9.81 # Gravity constant 
    Cp = 1005. # Specific heat of dry air at constant pressure

    # Now regrid the data to cylindrical coordinates
    r,phi = xy2rp(cx,cy,x,y)
    dbz_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,dbz)
    dbzlin_rpz = 10.0.^(dbz_rpz ./ 10)
    u_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,u)
    v_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,v)
    temp_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,temp)
    theta_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,theta)
    qv_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,qv)
    rhoa_rpz = regrid_xyz2rpz(cx,cy,x,y,z,r,phi,rhoa)

    # Convert u,v to ur,vt 
    ur_rpz,vt_rpz = uv2urvt(phi,u_rpz,v_rpz)

    # Compute the azimuthal means    
    azmean_dbzlin = nanmean(dbzlin_rpz,2,true)
    azmean_dbz = 10.0 .* log10.(azmean_dbzlin)
    azmean_temp = nanmean(temp_rpz,2,true)
    azmean_theta = nanmean(theta_rpz,2,true)
    azmean_qv = nanmean(qv_rpz,2,true)
    azmean_rhoa = nanmean(rhoa_rpz,2,true)
    azmean_vt = nanmean(vt_rpz,2,true)

    # Specify the radius to compute the ambient profile
    # Following Annette's convention
    ind_prof = findfirst(r,abs(x[1]))

    # Compute the ambient profile 
    amb_exner = Array{Float64}(length(z))
    amb_thetarho = Array{Float64}(length(z))     

    for k in eachindex(z)
        amb_exner[k],amb_thetarho[k] = ambientprof(z[k],azmean_dbz[ind_prof,k],azmean_rhoa[ind_prof,k],
                                                   azmean_temp[ind_prof,k],azmean_theta[ind_prof,k],
                                                   azmean_qv[ind_prof,k])
    end

    # Convert to SI units prior to integration 
    rmet = r .* 1e3
    zmet = z .* 1e3

    # Integrate the hydrostatic equation for the ambient profile
    # This will give the hydrostatic pressure at each altitude
    amb_hexner = Array{Float64}(length(z))

    for k in eachindex(z)
        if k == 1
            amb_hexner[k] = amb_exner[1]
        else
            if any(isnan.(amb_thetarho[1:k]))
                amb_hexner[k] = NaN
            else
                amb_hexner[k] = amb_exner[1] - (g/Cp) * newtoncotes(zmet[1:k], 1 ./ amb_thetarho[1:k])
            end
        end
    end

#==============================================================================
Branch off to the specified mode
==============================================================================#

    # Define arrays for pib, trb, dpibdx, dpibdy

    pib = Array{Float64}(size(azmean_vt))
    trb = Array{Float64}(size(azmean_vt))
    pib_xy = similar(u)
    trb_xy = similar(u)

    fill!(pib,NaN)
    fill!(trb,NaN)

#============================== HYDROSTATIC ==================================#

    if mode == "hydrostatic"
        # Use the ambient profiles to create hydrostatic base state
        # Can skip ahead and broadcast to pib_xy, trb_xy
        for j in eachindex(x)
            for i in eachindex(y)
                pib_xy[i,j,:] = amb_hexner
                trb_xy[i,j,:] = amb_thetarho
            end
        end
        # Compute dpibdx and dpibdy
        # ** NOTE **
        # dpibdx and dpibdy are required in units of 1/km and conversion to 
        # SI units takes place in NetCDF_XYZ.cpp of thermodynamic retrieval
        dpibdx = finite_dx(x,pib_xy)
        dpibdy = finite_dy(y,pib_xy)
        # Replace all NaNs with -999.0 required by thermodynamic retrieval
        pib_xy[isnan.(pib_xy)] = -999.0
        trb_xy[isnan.(trb_xy)] = -999.0
        dpibdx[isnan.(dpibdx)] = -999.0
        dpibdy[isnan.(dpibdy)] = -999.0
    end

#============================== VORTEX_CART ==================================#

    if mode == "vortex_cart"
        if paramBL
            # Set vt below 2 km to the value at 2 km 
            # This forces the boundary layer parameterization
            BLtop = closest_ind(z,2.0)
            azmean_vt[:,1:BLtop] .= azmean_vt[:,BLtop]
        end

        # Create interpolation objects for hydrostatic exner function
        # and thetarho at ambient profile
        amb_hexner_itp = interpolate((zmet,),amb_hexner,Gridded(Linear()))
        amb_thetarho_itp = interpolate((zmet,),amb_thetarho,Gridded(Linear()))
    
        # Compute required variables
        C = azmean_vt.^2 ./ rmet
        C[1,:] = NaN # Undefined at r = 0
        dCdz = finite_dz(zmet,C)

        # **Note: The integration only goes out to ind_prof (radius of ambient profile!)

        for k in eachindex(z)
            for i in 1:length(r[1:ind_prof])-1
                # Don't allow NaNs along the integration path
                if any(isnan.(C[i:ind_prof,k])) || any(isnan.(dCdz[i:ind_prof,k])) 
                    pib[i,k] = NaN
                    trb[i,k] = NaN
                else
                    z_R = zmet[k] + newtoncotes(rmet[i:ind_prof],C[i:ind_prof,k]) / g
                    pib[i,k] = amb_hexner_itp[z_R]
                    trb[i,k] = exp(log(amb_thetarho_itp[z_R]) - (1. / g) * newtoncotes(rmet[i:ind_prof],dCdz[i:ind_prof,k]))
                end
            end
        end 

        # Interpolate pib and trb to Cartesian grid
        for k in eachindex(z)
            # Interpolate to the x-y grid
            pib_xy[:,:,k] = regrid_pol2cart(cx,cy,r[1:ind_prof],x,y,pib[1:ind_prof,k])
            trb_xy[:,:,k] = regrid_pol2cart(cx,cy,r[1:ind_prof],x,y,trb[1:ind_prof,k])
        end
        
        # Compute dpibdx and dpibdy
        # ** NOTE **
        # dpibdx and dpibdy are required in units of 1/km and conversion to 
        # SI units takes place in NetCDF_XYZ.cpp of thermodynamic retrieval
        dpibdx = finite_dx(x,pib_xy)
        dpibdy = finite_dy(y,pib_xy)

        # Replace all NaNs with -999.0 required by thermodynamic retrieval
        pib_xy[isnan.(pib_xy)] = -999.0
        trb_xy[isnan.(trb_xy)] = -999.0
        dpibdx[isnan.(dpibdx)] = -999.0
        dpibdy[isnan.(dpibdy)] = -999.0
    end 
    
#==============================================================================
# Write all the variables to the NetCDF file
# Rename dictionary entries according to what's required by thermodynamic 
# retrieval and then push additional variables
==============================================================================#

if writeout

# Set u,v missing values to -999.0

    u[isnan.(u)] = -999.0
    v[isnan.(v)] = -999.0

# Place the required coordinate arrays in a dict and variables in another dict
# Need to add the singleton (time) dimension back to variables

    coord_dict = OrderedDict()
    thermo_dict = OrderedDict()

    coord_dict["lon"] = lon
    coord_dict["lat"] = lat
    coord_dict["x"] = x
    coord_dict["y"] = y
    coord_dict["z"] = z

    thermo_dict["u"] = reshape(u,(length(x),length(y),length(z),1))
    thermo_dict["v"] = reshape(v,(length(x),length(y),length(z),1))
    thermo_dict["w"] = reshape(w,(length(x),length(y),length(z),1))
    thermo_dict["dudx"] = reshape(dudx,(length(x),length(y),length(z),1))
    thermo_dict["dvdx"] = reshape(dvdx,(length(x),length(y),length(z),1))
    thermo_dict["dwdx"] = reshape(dwdx,(length(x),length(y),length(z),1))
    thermo_dict["dudy"] = reshape(dudy,(length(x),length(y),length(z),1))
    thermo_dict["dvdy"] = reshape(dvdy,(length(x),length(y),length(z),1))
    thermo_dict["dwdy"] = reshape(dwdy,(length(x),length(y),length(z),1))
    thermo_dict["dudz"] = reshape(dudz,(length(x),length(y),length(z),1))
    thermo_dict["dvdz"] = reshape(dvdz,(length(x),length(y),length(z),1))
    thermo_dict["dwdz"] = reshape(dwdz,(length(x),length(y),length(z),1))
    thermo_dict["trb"] = reshape(trb_xy,(length(x),length(y),length(z),1))
    thermo_dict["pib"] = reshape(pib_xy,(length(x),length(y),length(z),1))
    thermo_dict["dpibdx"] = reshape(dpibdx,(length(x),length(y),length(z),1))
    thermo_dict["dpibdy"] = reshape(dpibdy,(length(x),length(y),length(z),1))

# Define an array of strings for the vars needed by thermodynamic retrieval
# (excluding lon,lat,x,y,z coordinate arrays)

    varnames = ["u","v","w","dudx","dvdx","dwdx","dudy","dvdy","dwdy",
                "dudz","dvdz","dwdz","trb","pib","dpibdx","dpibdy"]

# Write to output NetCDF file

    writeout_ncvars(fout,coord_dict,thermo_dict)

end



