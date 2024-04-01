# *****************************************************************************
# thermo.jl
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# This script contains several functions to compute thermodynamic quantities
# and carry out conversions
#
# Function list
# get_bolton_thetae
# get_sat_vap_prs 
# *****************************************************************************

#=============================================================================
get_bolton_thetae

Calculate the pseudo-adiabatic equivalent potential temperature following the
empirical method outlined by Bolton (1980)

Steps
Compute the vapor pressure
Compute the dewpoint temperature
Compute temperature at LCL
Compute potential temperature (Eq. 7 with rv in units of g/kg, or equivalently,
Eq. 6 from Bryan (2008) with rv in units of kg/kg)

Variable units within function

# T_K = temperature (K)
# T_D = dewpoint temperature (deg C)
# T_DK = dewpoint temperature (K)
# prs = pressure (hPa)
# rv = water vapor mixing ratio (kg/kg)
# R_d = gas constant for dry air (J K^-1 kg^-1)
# R_v = gas constant for water vapor (J K^-1 kg^-1)
# epsilon (Rd/Rv) = 0.622
=============================================================================#

function get_bolton_thetae(T_K::Real,T_D::Real,prs::Real)

    # Define constants 
    R_d = 287.04
    R_v = 461.5
    epsilon = R_d/R_v
    # Calculate the vapor pressure (vap_prs) from Eq. 16 in Bolton (1980)
    # rv units = kg/kg
    # prs units determines units of vap_prs
    vap_prs = prs*rv/(epsilon + rv)
    # Eq. 7 in Bolton (1980) for potential temperature
    # rv units = kg/kg
    # prs units = hPa
    theta = T_K * (1000/prs)^(0.2854 * (1 - 0.28 * rv))
    # Invert Eq. 10 in Bolton (1980) to get dewpoint temperature (deg C)
    # vap_prs units = hPa
    T_D = 243.5 / ( (17.67 / log(vap_prs/6.112)) - 1 )
    T_DK = T_D + 273.15
    # Eq. 15 in Bolton (1980) for the temperature at the LCL
    # All temperatures in units of K
    T_L = 1 / ( (1 / (T_DK - 56)) + log(T_K/T_DK)/800 ) + 56
    # Eq. 38 from Bolton (1980) or equivalently Eq. 6 from Bryan (2008)
    # theta units = K
    # TL units = K
    # rv units = kg/kg
    theta_e = theta * exp( (3376/T_L - 2.54) * rv * (1 + 0.81*rv))
    return theta_e
end

#=============================================================================
get_sat_vap_prs 

Calculate the saturation vapor pressure following Bolton (1980) - Eq. (10)

Input 
T_K - temperature (Kelvin)

Output 
Saturation vapor pressure (Pa)
=============================================================================#

function get_sat_vap_prs(T_K::Real)
    return 611.2 * exp( (17.67 * (T_K - 273.15)) / (T_K - 29.65) )
end 

#=============================================================================
Calculate the saturation mixing ratio with respect to liquid water

- Eigth order polynomial following Flatau et al. (1992; Table 4, relative 
  error norm)
- Cross-referenced with MPAS, Thompson Microphysics, and Morrison Microphysics

Input
prs - Pressure (Pa)
temp - Temperature (K)

Output
Saturation mixing ratio r_s (kg/kg)
=============================================================================#

function rslf(prs::Real,temp::Real)

    a0 = 0.611583699e03
    a1 = 0.444606896e02
    a2 = 0.143177157e01
    a3 = 0.264224321e-1
    a4 = 0.299291081e-3
    a5 = 0.203154182e-5
    a6 = 0.702620698e-8
    a7 = 0.379534310e-11
    a8 =-0.321582393e-13
    x = max(-80.,temp-273.16)
    e_sl = a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*(a6+x*(a7+x*a8)))))))
    r_sl = 0.622*e_sl/(prs-e_sl)
    return r_sl
end 

#=============================================================================
Calculate the saturation mixing ratio with respect to ice

- Eigth order polynomial following Flatau et al. (1992; Table 4, relative 
  error norm)
- Cross-verified with MPAS, Thompson Microphysics, and Morrison Microphysics

Input
prs - Pressure (Pa)
temp - Temperature (K)

Output
Saturation mixing ratio r_s (kg/kg)
=============================================================================#

function rsif(prs::Real,temp::Real)

    a0 = 0.609868993e03
    a1 = 0.499320233e02
    a2 = 0.184672631e01
    a3 = 0.402737184e-1
    a4 = 0.565392987e-3
    a5 = 0.521693933e-5
    a6 = 0.307839583e-7
    a7 = 0.105785160e-9
    a8 = 0.161444444e-12
    x = max(-80.,temp-273.16)
    e_si = a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*(a6+x*(a7+x*a8)))))))
    r_si = 0.622*e_si/(prs-e_si)
    return r_si
end 
