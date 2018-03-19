module JuliaMet

    using NetCDF, DataStructures, Interpolations, Colors, PyPlot

    export p_centroid
    export extract_rgbs, create_cmap
    export calc_rmw, uv2urvt, rankine, hermite, p3swploc
    export finite_dr, finite_dx, finite_dy, finite_dz, finite_laplacian
    export rmnan, fourierols, wavecoeffs
    export newtoncotes, trapz2d, trapz3d
    export nansum, nanmean, nanvar, nanstd, nanmin, nanmax, nanextrema
    export read_ncvars
    export closest_ind, grid2d, grid3d, xy2rp, regrid_xy2rp, regrid_xyz2rpz, 
           regrid_pol2cart 
    export obslocs_rtz
    export stationary_frame


    include("center_finding.jl")
    include("cmaps.jl")
    include("cyclone.jl")
    include("derivative.jl")    
    include("harmonics.jl")
    include("integrate.jl")
    include("nanstats.jl")
    include("nc_tools.jl")
    include("regrid.jl")
    include("samurai_obslocs.jl")
    include("steady_frame.jl")

end # module
