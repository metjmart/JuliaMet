__precompile__()

module JuliaMet

    using DataStructures, Interpolations, LinearAlgebra, SparseArrays, NetCDF, Statistics, DelimitedFiles, DSP

    # Define a function to search a given path for all files containing the given key
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

    export searchdir
    export filtloop, psi_centroid_xy, p_centroid, p_centroid_
    export calc_rmw, uv2urvt, p3swploc
    export finite_dr, finite_dx, finite_dy, finite_dz, finite_dt, finite_laplacian, kron_delta,
           fd_weights, finite_ds, finite_dsx, finite_dsy, finite_dsz, diff1, spdiff1, Laplacian
    export fd_weights4, finite_ds4, finite_dsx4, finite_dsy4, finite_dsz4
    export haversine, haversined, invhaversine, angular_dist, polar_azimuth_equidist, polar_gc
    export fourierols, wavecoeffs
    export newtoncotes, trapz2d, trapz3d
    export nansum, nanmean, nanvar, nanstd, nanmin, nanmax, nanextrema, nanfindmax,
           nanargmax, nanfindmin, nanargmin
    export read_ncvars
    export closest_ind, grid2d, grid3d, xy2rp, regrid_xy2rp, regrid_xyz2rpz,
           regrid_pol2cart, unstagger, curv2rect, curv2rect_wgts
    export obslocs_rtz
    export get_bolton_thetae, get_sat_vap_prs, rslf, rsif
    export rankine, smrankine, modrankine, hermite, re87, cw87, wc04

    include("center_finding.jl")
    include("cyclone.jl")
    include("derivative.jl")
    include("derivative4.jl")
    include("geodesic.jl")
    include("harmonics.jl")
    include("integrate.jl")
    include("nanstats.jl")
    include("nc_tools.jl")
    include("regrid.jl")
    include("samurai_obslocs.jl")
    include("thermo.jl")
    include("vortexprofs.jl")

end # module
