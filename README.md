# JuliaMet

## Introduction

Simple module developed to analyze gridded meteorological data and model output. 

Operating Julia version: 1.0.0

## Installation

Begin by using `git clone` to clone the repository to a local directory
```
git clone https://github.com/metjmart/JuliaMet.git
```

Since JuliaMet is not an official Julia package, it must be loaded as 
a module. This can be done by manually extending the global variable 
`LOAD_PATH`. This variable contains the directories which Julia will search 
for modules when calling `require`. 

To extend the `LOAD_PATH` variable, open the `~/.julia/config/startup.jl` file 
(or create it in the above specified directory if it doesn't exist) and add the following to the file
(see the [Modules documentation](https://docs.julialang.org/en/v1/manual/modules/index.html))
```julia
push!(LOAD_PATH, "/path/to/JuliaMet/src/")
```
where `"/path/to/"` is replaced by the path to the user's local copy of JuliaMet.
This will extend the `LOAD_PATH` variable on every Julia initialization. 

## Package dependencies 

`Project.toml` and `Manifest.toml` files have been included to instantiate an
environment with all the necessary package dependencies for running JuliaMet.
Note that JuliaMet has not been updated since Julia v1.0. 

To create the JuliaMet environment, initiate the Julia REPL from the `JuliaMet/` 
directory and hit the `]` key to enter the Pkg REPL-mode. Then, run the following
sequence of commands
```julia
activate .
instantiate
```
Once the JuliaMet environment has been resolved, subsequent executions of Julia code
dependent on JuliaMet can be executed from the command line via
```bash
julia --project=/path/to/JuliaMet
```
where `"/path/to/"` is replaced by the path to the user's local copy of JuliaMet.

## Run

We should now be able to load JuliaMet as a custom
module by adding the following to the beginning of our code
```julia
using JuliaMet
```

## Utility

### center_finding.jl
* Functions to determine the center of a fluid circulation

### cyclone.jl 
* Functions to compute relevant quantities for gridded tropical cyclone data

### derivative.jl 
* Numerical methods for computing the finite difference approximation to first-order derivatives with second-order accuracy 

### geodesic.jl
* Functions to compute distances and angles between points on a sphere

### harmonics.jl
* Functions to compute discrete Fourier transforms via ordinary least squares regression

### integrate.jl
* Numerical methods for approximating integrals using the closed Newtown-Cotes formulae

### nanstats.jl 
* Redefined statistical functions in Julia that handle the presence of NaNs

### nc_tools.jl 
* Extensive NetCDF data reading function (NetCDF writing function in progress)

### regrid.jl 
* Re-grid data between Cartesian and cylindrical coorindates using the Interpolations.jl package

### samurai_obslocs.jl 
* Determine the spatial location of data input into [SAMURAI](https://github.com/mmbell/samurai)

### thermo.jl
* Functions to calculate and convert various thermodynamic quantities 

### vortexprofs.jl
* Functions to create analytical vortex profiles of tangential (azimuthal) velocity 
