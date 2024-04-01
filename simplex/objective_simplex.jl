# *****************************************************************************
# objective_simplex
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# Top-level script
# Following Lee and Marks (2000; MWR) and section 2 of Bell and Lee (2012; JAMC)
# Determine the optimal center location of a tropical cyclone circulation
# The objective simplex center-finding algorithm will launch several simplexes
# at a specified number of grid points for an initial center and RMW guess.
# Note - This script was developed to handle observational analyses that 
#        contain missing values (NaNs)
# *****************************************************************************

using JuliaMet, Interpolations, Optim, NetCDF, DelimitedFiles, Statistics
include("./simplex_aux.jl")

#==============================================================================
Interactive portion - The algorithm requires the following variables

x - east-west coordinate array (units = km)
y - north-south coordinate array (units = km)
z - vertical coordinate array (units = km)
u_glob - u-wind (east-west); size(u_glob) = (length(x),length(y),length(z));
         units = m/s
v_glob - v-wind (north-south); size(v_glob) = (length(x),length(y),length(z));
         units = m/s
zs - index of the lowest vertical level to run the algorithm
ze - index of highest vertical level to run the algorithm

init_config is a function that requires the following arguments

x-location of first center guess (km)
y-location of first center guess (km)
First-guess radius of maximum tangential winds (km)
sqrt(# of grid points to initialize simplexes) - e.g., 3 for 3x3, 5 for 5x5

The following variables are produced by the init_config function
*Note that the length of xinit and yinit correspond to the example provided
below where there are nine (3x3) initial simplex locations

xinit - x-location for each initial center guess (km); length(xinit) = 3
yinit - y-location for each initial center guess (km); length(yinit) = 3
radii - radius rings centered on the first-guess RMW (km); length(radii) = 11
        (default is +/- 5 km in 1 km increments; can be changed in simplex_aux.jl)

fout - Output NetCDF file name
==============================================================================#

    fin = "./test/samurai_XYZ_analysis.nc"
    vars = ["x","y","altitude","U","V"]

    const x,y,z,u_glob,v_glob = read_ncvars(fin,vars)
    const zs = closest_ind(z,2.0)
    const ze = closest_ind(z,8.0)
    const xinit,yinit,radii = init_config(0.,0.,15.,3)
    const fout = "./20151022_pass1_test3_center.nc"
    zlevs = zs:ze

#==============================================================================
# Initialize all required variables
==============================================================================#

    simplex_xcenters = Array{Float64}(undef,length(xinit),length(yinit),length(radii),length(zlevs))
    simplex_ycenters = Array{Float64}(undef,length(xinit),length(yinit),length(radii),length(zlevs))
    simplex_numiters = Array{Float64}(undef,length(xinit),length(yinit),length(radii),length(zlevs))

    prelim_xbar_centers = Array{Float64}(undef,length(radii),length(zlevs))
    prelim_ybar_centers = Array{Float64}(undef,length(radii),length(zlevs))
    prelim_stdv_centers = Array{Float64}(undef,length(radii),length(zlevs))

    final_xc = Array{Float64}(undef,length(zlevs))
    final_yc = Array{Float64}(undef,length(zlevs))
    final_rmw = Array{Float64}(undef,length(zlevs))
    final_std = Array{Float64}(undef,length(zlevs))
    mean_vts = Array{Float64}(undef,length(radii),length(zlevs))

# Compute dzmet

    zmet = z*1e3

    dzmet = zeros(length(z))

    for k in 1:length(z)
        if k == 1
            dzmet[k] = zmet[k] - 0.
        else
            dzmet[k] = zmet[k] - zmet[k-1]
        end
    end

# Confine dzmet to zlevs

    dzmet_zlayer = dzmet[zlevs]

#==============================================================================
# Run the algorithm
!!! Do not modify below !!!
==============================================================================#

    ind = 0
    for k in zlevs
        global ind += 1
        for ir in eachindex(radii)
            for j in eachindex(yinit)
                for i in eachindex(xinit)
                    solution = optimize(loc -> nanmeanvt_annulus(loc,radii[ir],u_glob[:,:,k],v_glob[:,:,k]),
                                        [xinit[i], yinit[j]], NelderMead(), Optim.Options(g_tol=1e-4, iterations=30))
                    @inbounds simplex_xcenters[i,j,ir,ind] = solution.minimizer[1]
                    @inbounds simplex_ycenters[i,j,ir,ind] = solution.minimizer[2]
                    @inbounds simplex_numiters[i,j,ir,ind] = solution.iterations
                end
            end
            # Remove the outliers for each vertical level and radius ring
            # and return the preliminary center information
            @inbounds prelim_xbar_centers[ir,ind],
                      prelim_ybar_centers[ir,ind],
                      prelim_stdv_centers[ir,ind] =
                      rm_outliers(xinit,yinit,simplex_xcenters[:,:,ir,ind],simplex_ycenters[:,:,ir,ind])
            # Compute the azimuthal mean tangential wind for the preliminary centers
            @inbounds mean_vts[ir,ind] = nanmeanvt_ring(x,y,
                                                        prelim_xbar_centers[ir,ind],
                                                        prelim_ybar_centers[ir,ind],radii[ir],
                                                        u_glob[:,:,k],v_glob[:,:,k])
        end

        # Find the radial index that maximizes the mean vt at each vertical level
        ind_rmax = nanargmax(mean_vts[:,ind])
        # Use index to store final center, rmw, and std
        final_xc[ind] = prelim_xbar_centers[ind_rmax,ind]
        final_yc[ind] = prelim_ybar_centers[ind_rmax,ind]
        final_rmw[ind] = radii[ind_rmax]
        final_std[ind] = prelim_stdv_centers[ind_rmax,ind]
        println("Final Center @ z = ", z[k], " km: ", final_xc[ind], ", ", final_yc[ind],
                ", std = ", final_std[ind], ", rmw = ", final_rmw[ind])
    end

# Compute the layer-avearged simplex center

    final_xc_mean = dropdims(sum(final_xc .* dzmet_zlayer,dims=1),dims=1) / sum(dzmet_zlayer)
    final_yc_mean = dropdims(sum(final_yc .* dzmet_zlayer,dims=1),dims=1) / sum(dzmet_zlayer)

# Print the final result

    println(" ")
    println(z[zs], "-", z[ze], " km layer-average center: (", final_xc_mean, ", ", final_yc_mean, ")")

# Store desired output in a NetCDF file

# Define dimensions

    xinit_dim = NcDim("nxinit", length(xinit))
    yinit_dim = NcDim("nyinit", length(yinit))
    rdim = NcDim("nr", length(radii))
    zdim = NcDim("nz",length(z[zlevs]))

# Create vars

    xinit_var = NcVar("xinit",[xinit_dim],t=Float64)
    yinit_var = NcVar("yinit",[yinit_dim],t=Float64)
    rvar = NcVar("radii",[rdim],t=Float64)
    zlayervar = NcVar("zlayer",[zdim],t=Float64)
    numitersvar = NcVar("numiters",[xinit_dim,yinit_dim,rdim,zdim],t=Float64)
    final_xcvar = NcVar("final_xc",[zdim],t=Float64)
    final_ycvar = NcVar("final_yc",[zdim],t=Float64)
    final_stdvar = NcVar("final_std",[zdim],t=Float64)
    final_rmwvar = NcVar("final_rmw",[zdim],t=Float64)
    dzmet_zlayervar = NcVar("dzmet_zlayer",[zdim],t=Float64)

# Remove any preexisting NetCDF file with desired name

    isfile(fout) ? rm(fout) : nothing
    nc = NetCDF.create(fout,xinit_var,yinit_var,rvar,zlayervar,numitersvar,
                       final_xcvar,final_ycvar,final_stdvar,final_rmwvar,
                       dzmet_zlayervar)

# Store the data in the NetCDF file

    println("Storing data in NetCDF file ...")

    NetCDF.putvar(nc,"xinit",xinit)
    NetCDF.putvar(nc,"yinit",yinit)
    NetCDF.putvar(nc,"radii",radii)
    NetCDF.putvar(nc,"zlayer",z[zlevs])
    NetCDF.putvar(nc,"numiters",simplex_numiters)
    NetCDF.putvar(nc,"final_xc",final_xc)
    NetCDF.putvar(nc,"final_yc",final_yc)
    NetCDF.putvar(nc,"final_std",final_std)
    NetCDF.putvar(nc,"final_rmw",final_rmw)
    NetCDF.putvar(nc,"dzmet_zlayer",dzmet_zlayer)

    NetCDF.ncclose(fout)

    println("Complete!")
