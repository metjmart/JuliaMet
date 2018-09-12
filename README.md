# JuliaMet

## Introduction

Functions developed to analyze gridded meteorological data and model output. 

Operating Julia version: 1.0.0

## Package dependencies 

The following packages are not included in either Julia's base or stdlib, and 
should be installed prior to using JuliaMet

* [DataStructures](https://github.com/JuliaCollections/DataStructures.jl)
* [Interpolations](https://github.com/JuliaMath/Interpolations.jl)
* [NetCDF](https://github.com/JuliaGeo/NetCDF.jl)

This can be accomplished by first opening the Julia REPL and entering the Pkg 
REPL-mode by hitting the `]` key.
Then, packages can be added via the following command
```julia
add <pkgname>
```
where `<pkgname>` is replaced by any of the package names listed above.

## Installation

Begin by using `git clone` to clone the repository to a local directory
```
git clone https://github.com/metjmart/JuliaMet.git
```

Since JuliaMet is not (yet) an official Julia package, it must be loaded as 
a module. This can be done by manually extending the global variable 
`LOAD_PATH`. This variable contains the directories which Julia will search 
for modules when calling `require`. 

To extend the `LOAD_PATH` variable, open the `~/.julia/config/startup.jl` file 
(or create it in the above specified directory if it doesn't exist) and add the following to the file
(see the [Modules documentation](https://docs.julialang.org/en/stable/manual/modules/))
```julia
push!(LOAD_PATH, "/path/to/JuliaMet/src/")
```
where `"/path/to/"` is replaced by the path to the user's local copy of JuliaMet.
This will extend the `LOAD_PATH` variable on every Julia initialization. 
Alternatively, one could define an environment variable `JULIA_LOAD_PATH` in 
the `~/.bashrc` file (or which ever shell is used) and then source the file
```
export JULIA_LOAD_PATH="/path/to/JuliaMet/src/"
source ~/.bashrc
```
Either method will suffice! We should now be able to load JuliaMet as a custom
module by adding the following to the beginning of our code
```julia
using JuliaMet
```

## Utility

### center_finding.jl
* Functions to determine the center of a tropical cyclone 

### cyclone.jl 
* Functions to compute relevant quantities for gridded tropical cyclone data

### derivative.jl 
* Numerical methods for computing derivatives using finite differencing

### harmonics.jl
* Functions to compute discrete Fourier transforms via ordinary least squares regression

### integrate.jl
* Numerical methods for integrating discretized data using the closed Newtown-Cotes formulae

### nanstats.jl 
* Redefined statistical functions in Julia that handle the presence of NaNs

### nc_tools.jl 
* Extensive NetCDF data reading function (NetCDF writing function in progress)

### regrid.jl 
* Re-grid data between Cartesian and cylindrical coorindates using the Interpolations.jl package

### samurai_obslocs.jl 
* Determine the spatial location of data input into [SAMURAI](https://github.com/mmbell/samurai)

### steady_frame.jl
* Functions used to compute the [most steady frame of reference](https://journals.ametsoc.org/doi/abs/10.1175/1520-0426(2002)019%3C1035%3AETMSFO%3E2.0.CO%3B2) given Doppler velocity observations

### vortexprofs.jl
* Functions to create radial profiles of tangential velocity widely used to initiate tropical cyclone simulations

**Disclaimer**: The functions included in JuliaMet have been primarily tested by myself and therefore may contain bugs.
Users are encouraged to examine the source code prior to using any function in JuliaMet.
Please contribute if you discover bugs, find a more efficient way of doing something, or if you have functions you'd like to include.
Documentation and officialization coming soon!
