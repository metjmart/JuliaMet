# JuliaMet

## Introduction

Functions developed to analyze gridded meteorological data and model output. 

Julia version: 0.4.5

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
