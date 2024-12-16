# *****************************************************************************
# constants.jl
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.10.4
#
# This script contains data types to store constants for different applications
#
# *****************************************************************************

"""
    WRFConstants

Store constants from WRF model environment.

Constants are defined following WRF definitions under shared/module_constants.F.\\
An instance `wrfcs` is generated to precompile the WRF constants.

# Defined fields
- `r_a`: Gas constant for dry air (J K^-1 kg^-1)
- `r_v`: Gas constant for vapor (J K^-1 kg^-1)
- `g`: Gravitational constant
- `p_00`: Reference pressure level (Pa)
- `t_00`: Reference potential temperature (K)

# Derived fields

- `c_pa`: Specific heat capacity for dry air at constant pressure (J K^-1 kg^-1)
- `epsilon`: r_a/r_v
- `kappa`: r_a/c_pa
"""
struct WRFConstants{T}
    # Defined constants
    r_a::T
    r_v::T
    g::T
    p_00::T
    t_00::T
    # Derived constants
    c_pa::T
    epsilon::T
    kappa::T
    function WRFConstants(r_a::T, r_v::T, g::T, p_00::T, t_00::T) where T<:Real
        c_pa = 7. * r_a / 2
        epsilon = r_a / r_v
        kappa = r_a / c_pa
        new{T}(r_a, r_v, g, p_00, t_00, c_pa, epsilon, kappa)
    end
  end
  
  wrfcs = WRFConstants(
    287.,
    461.6,
    9.81,
    100000.0,
    300.
  )
