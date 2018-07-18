__precompile__()

module JuliaMet

    using NetCDF, DataStructures, Interpolations

    export p_centroid, p_centroid2
    export calc_rmw, uv2urvt, p3swploc
    export finite_dr, finite_dx, finite_dy, finite_dz, finite_laplacian
    export rmnan, fourierols, wavecoeffs
    export newtoncotes, trapz2d, trapz3d
    export nansum, nanmean, nanvar, nanstd, nanmin, nanmax, nanextrema
    export read_ncvars
    export closest_ind, grid2d, grid3d, xy2rp, regrid_xy2rp, regrid_xyz2rpz, 
           regrid_pol2cart 
    export obslocs_rtz
    export steadyframe
    export rankine, modrankine, hermite, re87, cw87, wc04

    include("center_finding.jl")
    include("cyclone.jl")
    include("derivative.jl")    
    include("harmonics.jl")
    include("integrate.jl")
    include("nanstats.jl")
    include("nc_tools.jl")
    include("regrid.jl")
    include("samurai_obslocs.jl")
    include("steady_frame.jl")
    include("vortexprofs.jl")

end # module
