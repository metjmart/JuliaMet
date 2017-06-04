# JuliaMet

## Introduction

Functions developed to analyze gridded meteorological data and model output. 

Julia versions tested: 0.4.5, 0.5.2 

## Package dependencies 

The following packages should be installed prior to using JuliaMet:

* [NetCDF](https://github.com/JuliaGeo/NetCDF.jl)
* [DataStructures](https://github.com/JuliaCollections/DataStructures.jl)
* [Grid](https://github.com/timholy/Grid.jl)
* [Interpolations](https://github.com/JuliaMath/Interpolations.jl)
* [Colors](https://github.com/JuliaGraphics/Colors.jl)
* [PyPlot](https://github.com/JuliaPy/PyPlot.jl)

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
```
push!(LOAD_PATH, "/path/to/JuliaMet.jl/")
```
This will extend the `LOAD_PATH` on every Julia initialization. Alternatively, 
we could append additional directories to our `LOAD_PATH` evironment variable 
by adding the following to our `~/.bashrc` (or which ever shell is used) and 
sourcing `~/.bashrc`:
```
export JULIA_LOAD_PATH="/path/to/JuliaMet.jl/"
source ~/.bashrc
```
Either method will suffice! We should now be able to load JuliaMet as a custom
module by adding the following to the beginning of our code:
```julia
using JuliaMet
```

## Functions

### nc_tools.jl 
* Extensive NetCDF data reading function (NetCDF writing function in progress)

### derivative.jl 
* Functions to aid in numerical analyses (e.g., derivatives, means, etc.)

### integrate.jl
* Simple integration methods

### regrid.jl 
* Interpolate functions using the deprecated CoordInterpGrid.jl package 

### new_regrid.jl 
* Interpolate functions using the new Interpolations.jl package

### cyclone.jl 
* Simple functions to calculate relevant quantities with tropical cyclone data

### colormaps.jl 
* Access and create colorbars in Julia

### samurai_obslocs.jl 
* Determine the spatial location of data input into SAMURAI



