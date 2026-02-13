# *****************************************************************************
# wrfcoordvars.jl
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.10.4
#
# This script contains functions and data types to extract coordinate arrays
# from WRF NetCDF file output
#
# *****************************************************************************

"""
    WRFCoordVars(filepath::String)

Store unstaggered coordinate arrays and their respective dimensions from WRF simulations

This struct is specifically designed to work with WRF idealized TC simulations.\\
Note that the vertical grid is excluded, but its unstaggered size is included.

# Fields
- `dxm`: x-coordinate grid spacing (meters)
- `dym`: y-coordinate grid spacing (meters)
- `xm`: x-coordinate (meters)
- `ym`: y-coordinate (meters)
- `dx`: x-coordinate grid spacing (km)
- `dy`: y-coordinate grid spacing (km)
- `x`: x-coordinate (km)
- `y`: y-coordinate (km)
- `nx`: number of unstaggered x-coordinate grid cells 
- `ny`: number of unstaggered y-coordinate grid cells 
- `nz`: number of unstaggered z-coordinate grid cells 

# Construction 

The struct is created by providing a string that points to the WRF NetCDF 
filepath from which to extract variables:

    WRFCoordVars(filepath::String)

# Example

```julia
filepath = "path/to/wrfout_dxx_yyyy-mm-dd-hh:mm"
coordvars = WRFCoordVars(filepath)
```
"""
struct WRFCoordVars
    dxm::Float64
    dym::Float64
    xm::Vector{Float64}
    ym::Vector{Float64}
    dx::Float64
    dy::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    nx::Int64
    ny::Int64
    nz::Int64
end

"""
    WRFCoordVars(filepath::String)

Outer constructor for the WRFCoordVars struct

# Arguments
- `filepath::String`: path to WRF output file 
"""
function WRFCoordVars(filepath::String)
    WRFCoordVars(_WRFCoordVars(filepath)...)
end

"""
Extract unstaggered coordinate arrays and respective dimensions for WRFCoordVars
"""
function _WRFCoordVars(filepath::String)
    # get dx, dy from metadata
    dxm = ncgetatt(filepath, "Global", "DX")
    dym = ncgetatt(filepath, "Global", "DY")
    # get number of unstaggered u,v,w grid points from metadata
    nx = ncgetatt(filepath, "Global", "WEST-EAST_GRID_DIMENSION")
    ny = ncgetatt(filepath, "Global", "SOUTH-NORTH_GRID_DIMENSION")
    nz = ncgetatt(filepath, "Global", "BOTTOM-TOP_GRID_DIMENSION")
    # define x,y arrays 
    _xm = zeros(nx)
    _ym = zeros(ny)
    for i in 1:nx
        _xm[i] = (i-1) * dxm
    end
    for j in 1:ny
        _ym[j] = (j-1) * dym
    end
    # center
    _xm .-= _xm[div(nx,2)+1]
    _ym .-= _ym[div(ny,2)+1]
    # unstagger 
    xm = unstagger(_xm,dims=1)
    ym = unstagger(_ym,dims=1)
    return dxm,dym,xm,ym,dxm/1e3,dym/1e3,xm./1e3,ym./1e3,nx-1,ny-1,nz-1
end

"""
    WRFPostProcCoordVars(filepath::String)

Store post-processed coordinate arrays and their respective dimensions from WRF simulations

This struct is designed to work with post-processed WRF idealized TC simulations.  
Post-processed variables have been unstaggered and interpolated to altitude (vertical) coordinates.

# Fields
- `dxm`: x-coordinate grid spacing (meters)
- `dym`: y-coordinate grid spacing (meters)
- `dzm`: z-coordinate grid spacing (meters)
- `xm`: x-coordinate (meters)
- `ym`: y-coordinate (meters)
- `zm`: z-coordinate (meters)
- `dx`: x-coordinate grid spacing (km)
- `dy`: y-coordinate grid spacing (km)
- `dz`: z-coordinate grid spacing (km)
- `x`: x-coordinate (km)
- `y`: y-coordinate (km)
- `z`: z-coordinate (km)
- `nx`: number of unstaggered x-coordinate grid cells
- `ny`: number of unstaggered y-coordinate grid cells
- `nz`: number of unstaggered z-coordinate grid cells

# Construction 

The struct is created by providing a string that points to the post-processed  
WRF NetCDF filepath from which to extract variables:

    WRFPostProcCoordVars(filepath::String)

# Example

```julia
filepath = "path/to/wrfout_dxx_yyyy-mm-dd-hh:mm"
coordvars = WRFPostProcCoordVars(filepath)
```
"""
struct WRFPostProcCoordVars
    dxm::Float64
    dym::Float64
    dzm::Float64
    xm::Vector{Float64}
    ym::Vector{Float64}
    zm::Vector{Float64}
    dx::Float64
    dy::Float64
    dz::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    nx::Int64
    ny::Int64
    nz::Int64
end

"""
    WRFPostProcCoordVars(filepath::String)

Outer constructor for the WRFPostProcCoordVars struct

# Arguments
- `filepath::String`: path to post-processed WRF output file 
"""
function WRFPostProcCoordVars(filepath::String)
    WRFPostProcCoordVars(_WRFPostProcCoordVars(filepath)...)
end

"""
Extract unstaggered coordinate arrays and respective dimensions for WRFPostProcCoordVars
"""
function _WRFPostProcCoordVars(filepath::String)
    # get nx, ny, nz
    ds = NCDataset(filepath)
    nx, ny, nz = [ds.dim[key] for key in keys(ds.dim)]
    x = ds["x"][:]
    y = ds["y"][:]
    z = ds["z"][:]
    xm = x * 1e3
    ym = y * 1e3
    zm = z * 1e3
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    dz = z[2]-z[1]
    dxm = dx * 1e3
    dym = dy * 1e3
    dzm = dz * 1e3
    return dxm,dym,dzm,xm,ym,zm,dx,dy,dz,x,y,z,nx,ny,nz
end

"""
    get_wrf_files(filepath, domain; prefix)

Get the paths to WRF output files

# Arguments
- `filepath::String`: Path to WRF output files
- `domain::Int`: WRF domain (e.g., 1,2,3)
- `prefix::String`: prefix to WRF output file pattern  
    - Use for files other than raw WRF output but still contain the "wrfout_d0x" naming convention)
# Returns
- `fins::Vector{String}`: Vector of strings with the paths to WRF output files
"""
function get_wrf_files(filepath::String, domain::Int; prefix::String="^")
    # Convert domain to search for regex matching WRF output files
    pattern = "$(prefix)wrfout_d$(@sprintf("%02d",domain)).*"
    pattern = Regex(pattern)
    fins = joinpath.(filepath, [fin for fin in readdir(filepath) if !isnothing(match(pattern,fin))])
    return fins
end

"""
    get_wrf_files()

Parse command line arguments to get the paths to WRF output files 

# Returns
- `fins`: Array of strings with the paths to WRF output files
"""
function get_wrf_files()
    function parseargs()
        s = ArgParseSettings()
        @add_arg_table s begin
            "--filepath"
                help = "Path to WRF output files"
                required = true
                arg_type = String
            "--domain"
                help = "WRF domain (e.g., 1,2,3)"
                required = true
                arg_type = Int
        end
        return parse_args(s)
    end
    args = parseargs()
    filepath = args["filepath"]
    domain = args["domain"]
    # Convert domain to search for regex matching WRF output files
    pattern = "^wrfout_d$(@sprintf("%02d",domain)).*"
    pattern = Regex(pattern)
    fins = joinpath.(filepath, [fin for fin in readdir(filepath) if !isnothing(match(pattern,fin))])
    return fins
end

"""
    get_wrf_simulation(cfgpath::String)

Parse command line arguments to configure analyses for a WRF simulation

# Returns
- `sim`: Name of simulation parsed from CLI
- `dom`: Domain parsed from CLI
"""
function get_wrf_simulation(cfgpath::String)
    # get parent path from cfg
    parpath = TOML.parsefile(cfgpath)["config"]["parpath"]
    function parseargs()
        s = ArgParseSettings()
        @add_arg_table s begin
            "--simulation"
                help = "WRF simulation directory name"
                required = true
                arg_type = String
            "--domain"
                help = "WRF domain (e.g., 1,2,3)"
                required = true
                arg_type = Int
        end
        return parse_args(s)
    end
    args = parseargs()
    sim = args["simulation"]
    dom = args["domain"]
    return parpath, sim, dom
end

"""
    match_wrf_files(wrfoutfiles, xfiles)

Verify that the date/time of a list of WRF files match a corresponding list of files

For each WRF file in a list, get the date/time regex and match it to the
corresponding index of a second input list of strings. If the date/time regex
does not match, return an error

# Arguments
- `wrfout_files::AbstractVector{<:String}`: Vector of strings with WRF file names
- `xfiles::AbstractVector{<:String}`: Vector of strings to match to WRF files
# Output
- nothing (raise an error if a mismatch occurs)
"""
function match_wrf_files(wfiles::Vector{String}, xfiles::Vector{String})
    # YYYY-MM-DD_HH:MM:SS
    pattern = r"\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}"
    for (i,wfile) in enumerate(wfiles)
        datetime_pattern = match(pattern, wfile)
        if !occursin(datetime_pattern.match, xfiles[i])
            error("datetime $(datetime_pattern.match) does not match file $(xfiles[i])")
        end
    end
end

"""
    get_wrf_time_specs()

Parse command line arguments to get user-defined time specifications from a WRF simulation

# Returns
- `ts::Int`: start time index
- `te::Int`: end time index
- `ri::String`
    - input: flag to read/use rapid intensification index from a `rapid_init.toml` file as `te`
    - output: if `ri == true`: string with prefix for figure names "ri_"; if `ri == false`: empty string ""
- `ri_path::String`: path to `rapid_init.toml` file containing index `ri_ind`
"""
function get_wrf_time_specs()
    function parseargs()
        s = ArgParseSettings()
        @add_arg_table s begin
            "--ts"
                help = "Start time index"
                arg_type = Int
            "--te"
                help = "End time index"
                arg_type = Int
            "--ri"
                help = "Read rapid_init.toml ri_ind to ts"
                arg_type = String
            "--ri_path"
                help = "Path to rapid_init.toml file (default is directory where parent script is run)"
                arg_type = String
                default = ""
        end
        return parse_args(s)
    end
    args = parseargs()
    ts = args["ts"]
    te = args["te"]
    ri = lowercase(args["ri"])
    ri_path = args["ri_path"]
    if ri == "true"
        config = TOML.parsefile(joinpath(ri_path, "rapid_init.toml"))
        ri_ind = config["ri_ind"]
        ts = ri_ind
        desc = "ri"
    else
        desc = ""
    end
    if ts > te
        error("start time index ts cannot be larger than end time index te")
    end
    return ts, te, desc
end

"""
    WRFTimeSpecs

Store user-defined time specifications from a WRF simulation

This struct currently operates on indexes, not time units (e.g., hours or mins)

# Fields
- `ts::Int`: start time index
- `te::Int`: end time index
- `desc::String`: string with prefix for figure names

# Construction 

Call the outer-constructor with no args to parse command-line args  

```julia
    WRFTimeSpecs()
````

Or, pass args directly (note: this gives the option to define `desc` differently)

```julia
    WRFTimeSpecs(1, 168, "ri")
```

"""
struct WRFTimeSpecs
    ts::Int
    te::Int
    desc::String
end

"""
    WRFTimeSpecs

Outer constructor to generate the WRFTimeSpecs struct via command-line args from `get_wrf_time_specs()`
"""
function WRFTimeSpecs()
    WRFTimeSpecs(get_wrf_time_specs()...)
end
