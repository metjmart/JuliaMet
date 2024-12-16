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
    WRFCoordVars

Store unstaggered coordinate arrays and their respective dimensions from WRF simulations.\\
This struct is specifically designed to work with WRF idealized TC simulations.\\
Note that the vertical grid is excluded, but its unstaggered size is included.

# Fields
- `dxm`: x-coordinate grid spacing (meters)
- `dym`: y-coordinate grid spacing (meters)
- `xm`: x-coordinate (meters)
- `ym`: y-coordinate (meters)
- `nx`: number of unstaggered x-coordinate grid cells 
- `ny`: number of unstaggered y-coordinate grid cells 
- `nz`: number of unstaggered z-coordinate grid cells 

# Construction 

The struct is created by providing a string that points to the path/to/file
for the WRF NetCDF file from which to extract variables:

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

# Outer constructor for the WRFCoordVars struct

function WRFCoordVars(filepath::String)
    dxm,dym,xm,ym,dx,dy,x,y,nx,ny,nz = _WRFCoordVars(filepath)
    WRFCoordVars(dxm,dym,xm,ym,dx,dy,x,y,nx,ny,nz)
end

# Extract the unstaggered coordinate arrays and their respective dimensions
# for the WRFCoordVars struct 

function _WRFCoordVars(filepath::String)
    # Get dx, dy from metadata
    dxm = ncgetatt(filepath, "Global", "DX")
    dym = ncgetatt(filepath, "Global", "DY")
    # Get number of unstaggered u,v,w grid points from metadata
    nx = ncgetatt(filepath, "Global", "WEST-EAST_GRID_DIMENSION")
    ny = ncgetatt(filepath, "Global", "SOUTH-NORTH_GRID_DIMENSION")
    nz = ncgetatt(filepath, "Global", "BOTTOM-TOP_GRID_DIMENSION")
    # Define x,y arrays 
    _xm = zeros(nx)
    _ym = zeros(ny)
    for i in 1:nx
        _xm[i] = (i-1) * dxm
    end
    for j in 1:ny
        _ym[j] = (j-1) * dym
    end
    # Center on (_nx-1,_ny-1)
    _xm .-= _xm[div(nx,2)+1]
    _ym .-= _ym[div(ny,2)+1]
    # Unstagger 
    xm = unstagger(_xm,dims=1)
    ym = unstagger(_ym,dims=1)
    return dxm,dym,xm,ym,dxm/1e3,dym/1e3,xm./1e3,ym./1e3,nx-1,ny-1,nz-1
end
