__precompile__()

module JuliaMet

    using DataStructures, Interpolations, LinearAlgebra, NetCDF, Statistics, DelimitedFiles, DSP

    export filtloop, p_centroid, p_centroid_
    export calc_rmw, uv2urvt, p3swploc
    export finite_dr, finite_dx, finite_dy, finite_dz, finite_laplacian,
           fd_weights, finite_ds, finite_dsx, finite_dsy, finite_dsz
    export fd_weights4, finite_ds4, finite_dsx4, finite_dsy4, finite_dsz4
    export haversine, haversined, angular_dist, atan_full, azimuth_equidist, gc_azimuth
    export rmnan, fourierols, wavecoeffs
    export newtoncotes, trapz2d, trapz3d
    export nansum, nanmean, nanvar, nanstd, nanmin, nanmax, nanextrema, nanfindmax,
           nanargmax, nanfindmin, nanargmin
    export read_ncvars
    export closest_ind, grid2d, grid3d, xy2rp, regrid_xy2rp, regrid_xyz2rpz,
           regrid_pol2cart, unstagger
    export obslocs_rtz
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
    include("vortexprofs.jl")

end # module
