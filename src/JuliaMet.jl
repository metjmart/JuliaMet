module JuliaMet

    using NetCDF, DataStructures, Grid, Interpolations, Colors, PyPlot
    plt = PyPlot

    export read_ncvars
    export finite_dz, finite_dr, finite_dx, finite_dy
    export newton_cotes, trapz_1d, trapz_2d, trapz_3d
    export grid_2d, grid_3d, regrid_xy2rt, regrid_xyz2rtz
    export new_grid2d, new_grid3d, new_regrid_xy2rt, new_regrid_xyz2rtz, 
           regrid_gfrelxz, regrid_gfrel_xyz
    export calc_rmw, calc_azmean, uv2urvt
    export nanmax, nanmin, nansum, nanmean, nanmedian, nanvar, nanstd,
           nanskewness, nankurtosis
    export obslocs_rtz
    export create_colormap, extractrgbs, get_rgbs, make_line_colors, kwr_cmap, 
           kwk_cmap, blue_cmap, blue_r_cmap, red_cmap, red_r_cmap, radar_camp, 
           radar2_cmap, radar3_cmap, radar_carbone_cmap, radar_carbone2_cmap, 
           pvwcolors_cmap, pvcolors_cmap, ssec_cmap, bluegie_cmap
           

    include("nc_tools.jl")
    include("derivative.jl")    
    include("integrate.jl")
    include("regrid.jl")
    include("new_regrid.jl")
    include("cyclone.jl")
    include("nanstats.jl")
    include("samurai_obslocs.jl")
    include("colormaps.jl")

end # module
