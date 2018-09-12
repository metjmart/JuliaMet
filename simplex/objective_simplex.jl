# *****************************************************************************
# objective_simplex
#
# Top-level script
# Following Bell and Lee (2012), the objective simplex center-finding algorithm
# will launch several simplexes at a specified number of grid points for an
# initial center and RMW guess.
# *****************************************************************************

using JuliaMet, Interpolations, Optim, NetCDF, DelimitedFiles, Statistics
include("./simplex_aux.jl")

#==============================================================================
Interactive portion: Read in Cartesian coordinate arrays, (u,v), and specify
the initial center guess/RMW.
zs,ze - specify the lowest and highest vertical levels to run the algorithm 
        These will be the closest indices to the specified values
Final init_config arg is sqrt(# of grid points to initialize simplices)
E.g., 3 for 3x3, 4 for 4x4, 5 for 5x5
fout - Output NetCDF file name
==============================================================================#

    fin = "./tests/samurai_XYZ_analysis.nc"
    vars = ["x","y","altitude","U","V"]
    
    const x,y,z,u_glob,v_glob = read_ncvars(fin,vars)
    const zs = closest_ind(z,2.0)
    const ze = closest_ind(z,8.0)
    const xinit,yinit,radii = init_config(0.,0.,15.,3)
    const fout = "./20151022_pass1_test3_center.nc"

#==============================================================================
# Initialize all required variables
==============================================================================#

    simplex_xcenters = Array{Float64}(undef,length(xinit),length(yinit),length(radii),length(z[zs:ze]))
    simplex_ycenters = Array{Float64}(undef,length(xinit),length(yinit),length(radii),length(z[zs:ze]))
    simplex_numiters = Array{Float64}(undef,length(xinit),length(yinit),length(radii),length(z[zs:ze]))

    prelim_xbar_centers = Array{Float64}(undef,length(radii),length(z[zs:ze]))
    prelim_ybar_centers = Array{Float64}(undef,length(radii),length(z[zs:ze]))
    prelim_stdv_centers = Array{Float64}(undef,length(radii),length(z[zs:ze]))

    final_xc = Array{Float64}(undef,length(z[zs:ze]))
    final_yc = Array{Float64}(undef,length(z[zs:ze]))
    final_rmw = Array{Float64}(undef,length(z[zs:ze]))
    final_std = Array{Float64}(undef,length(z[zs:ze]))
    mean_vts = Array{Float64}(undef,length(radii),length(z[zs:ze]))

#==============================================================================
# Run the algorithm
!!! Do not modify below !!!
==============================================================================#

    ind = 0
    for k in zs:ze # 2-8 km
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
        ind_rmax = nanargmax(mean_vts[:,ind])[2]
        # Use index to store final center, rmw, and std
        final_xc[ind] = prelim_xbar_centers[ind_rmax,ind]
        final_yc[ind] = prelim_ybar_centers[ind_rmax,ind]
        final_rmw[ind] = radii[ind_rmax]
        final_std[ind] = prelim_stdv_centers[ind_rmax,ind]
        println("Final Center @ z = ", z[k], " km: ", final_xc[ind], ", ", final_yc[ind],
                ", std = ", final_std[ind], ", rmw = ", final_rmw[ind])
    end

# Print the final result

    println(" ")
    println("2-8 km layer average center: (", mean(final_xc), ", ", mean(final_yc), ")")

# Store desired output in a NetCDF file

# Define dimensions

    xinit_dim = NcDim("nxinit", length(xinit))
    yinit_dim = NcDim("nyinit", length(yinit))
    rdim = NcDim("nr", length(radii))
    zdim = NcDim("nz",length(z[zs:ze]))

# Create vars

    xinit_var = NcVar("xinit",[xinit_dim],t=Float64)
    yinit_var = NcVar("yinit",[yinit_dim],t=Float64)
    rvar = NcVar("radii",[rdim],t=Float64)
    zvar = NcVar("zlayer",[zdim],t=Float64)
    numitersvar = NcVar("numiters",[xinit_dim,yinit_dim,rdim,zdim],t=Float64)
    final_xcvar = NcVar("final_xc",[zdim],t=Float64)
    final_ycvar = NcVar("final_yc",[zdim],t=Float64)
    final_stdvar = NcVar("final_std",[zdim],t=Float64)
    final_rmwvar = NcVar("final_rmw",[zdim],t=Float64)

# Remove any preexisting NetCDF file with desired name

    isfile(fout) ? rm(fout) : nothing
    nc = NetCDF.create(fout,xinit_var,yinit_var,rvar,zvar,numitersvar,
                            final_xcvar,final_ycvar,final_stdvar,final_rmwvar)

# Store the data in the NetCDF file

    println("Storing data in NetCDF file ...")

    NetCDF.putvar(nc,"xinit",xinit)
    NetCDF.putvar(nc,"yinit",yinit)
    NetCDF.putvar(nc,"radii",radii)
    NetCDF.putvar(nc,"zlayer",z[zs:ze])
    NetCDF.putvar(nc,"numiters",simplex_numiters)
    NetCDF.putvar(nc,"final_xc",final_xc)
    NetCDF.putvar(nc,"final_yc",final_yc)
    NetCDF.putvar(nc,"final_std",final_std)
    NetCDF.putvar(nc,"final_rmw",final_rmw)

    NetCDF.ncclose(fout)

    println("Complete!")
