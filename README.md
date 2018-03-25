# JuliaMet

## Introduction

Functions developed to analyze gridded meteorological data and model output. 

Operating Julia version: 0.6.0

## Package dependencies 

The following packages should be installed prior to using JuliaMet:

* [NetCDF](https://github.com/JuliaGeo/NetCDF.jl)
* [DataStructures](https://github.com/JuliaCollections/DataStructures.jl)
* [Interpolations](https://github.com/JuliaMath/Interpolations.jl)

This can be accomplished by opening the Julia REPL and running
```julia
Pkg.add("Package")
```
where "Package" is replaced by the package names listed above.

## Installation

Begin by using `git clone` to clone the repository to a local directory:
```
git clone https://github.com/metjmart/JuliaMet.git
```

Since JuliaMet is not an official Julia package, it won't be located in 
Julia's `LOAD_PATH`. Therefore, we have to manually extend `LOAD_PATH` to 
include the location of JuliaMet. There are a couple ways to do this. 

First, it can be extended by opening the `~/.juliarc.jl` file (or creating it 
if it doesn't exist) and adding the following to the file
(see the [Modules documentation](https://docs.julialang.org/en/stable/manual/modules/)):
```julia
push!(LOAD_PATH, "/path/to/JuliaMet/src/")
```
This will extend the `LOAD_PATH` on every Julia initialization. Alternatively, 
we could append additional directories to our `LOAD_PATH` evironment variable 
by adding the following to our `~/.bashrc` (or which ever shell is used) and 
sourcing `~/.bashrc`:
```
export JULIA_LOAD_PATH="/path/to/JuliaMet/src/"
source ~/.bashrc
```
Either method will suffice! We should now be able to load JuliaMet as a custom
module by adding the following to the beginning of our code:
```julia
using JuliaMet
```

## Utility

### center_finding.jl
* Functions to determine the center of a tropical cyclone 

### cyclone.jl 
* Functions to compute relevant quantities from tropical cyclone data

### derivative.jl 
* Numerical methods for computing derivatives using finite differencing

### integrate.jl
* Numerical methods for integrating discretized data using the closed Newtown-Cotes formulae

### nanstats.jl 
* Redefined statistical functions in Julia that handle the presence of NaNs

### nc_tools.jl 
* Extensive NetCDF data reading function (NetCDF writing function in progress)

### regrid.jl 
* Regridding functions using the Interpolations.jl package

### samurai_obslocs.jl 
* Determine the spatial location of data input into SAMURAI

### steady_frame.jl
* Functions used to compute the most steady frame of reference given Doppler velocity observations


Disclaimer: The functions included in JuliaMet have largely been tested by myself and therefore may contain bugs.
Users are encouraged to examine the source code prior to using any function in JuliaMet.
Please contribute if you find bugs, find a more efficient way of doing something, or if you have your own functions!
Documentation and officialization coming soon!
