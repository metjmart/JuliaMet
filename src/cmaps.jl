# *****************************************************************************
# cmaps.jl
#
# Author: Jonathan Martinez 
# Email: jon.martinez@colostate.edu
# Julia version: 0.6.0
#
# Functions designed to work with colormaps in Julia. 
#
# Function list:
# extract_rgbs
# create_cmap
# *****************************************************************************

#==============================================================================
extract_rgbs

This function will extract a specified number of RGB triplets from a colorbar 
that was already created and stored in a matplotlib colormap.

cbar_name - string defining the colorbar name defined by the colormap
inc - number of increments (i.e., number of RGB triplets)
==============================================================================#

function extract_rgbs(cbar_name::AbstractString,inc::Real)

    julia_cmap = matplotlib[:cm][:get_cmap](cbar_name)
    rgbs = []
    rgb_vals = []
    for i in 1:inc
        itr = (i-1)/(inc-1) 
        push!(rgbs,julia_cmap(itr))
        push!(rgb_vals,collect(rgbs[i][1:3]))
        println(round.(Int,rgb_vals[i].*255))
    end

    return rgb_vals
end

#==============================================================================
create_cmap

This function is used locally to create the colormaps from the pre-defined 
RGB triplets in the colormaps directory. Users can also utilize the function
to create their own colormaps by adding a delimited file to the colormaps 
directory containing RGB triplets (see colormaps directory for examples).

Input args:
fin = file with delimited list of RGB triplets 
ncolors = number of discrete colors to be generated from input RGBs
==============================================================================#

function create_cmap(name::AbstractString,ncolors::Any)
   
    # Define path to JuliaMet and extend to colormaps directory
    path = []
    for ipath in LOAD_PATH
        if ismatch(r"JuliaMet",ipath)
            path = ipath
        end
    end
    fin = path*"colormaps/"*name
    rgb = readdlm(fin)
    n   = size(rgb)[1]
    rs  = Array{Tuple{Float64,Float64,Float64}}(0)
    gs  = Array{Tuple{Float64,Float64,Float64}}(0)
    bs  = Array{Tuple{Float64,Float64,Float64}}(0)
    for i in 1:n
        push!(rs,((i-1)/(n-1),rgb[i,1]/255.,rgb[i,1]/255.))
        push!(gs,((i-1)/(n-1),rgb[i,2]/255.,rgb[i,2]/255.))
        push!(bs,((i-1)/(n-1),rgb[i,3]/255.,rgb[i,3]/255.))
    end
    cmap = PyPlot.ColorMap(name,rs,gs,bs,ncolors,1.0) 
    PyPlot.register_cmap(name, cmap)

end

#==============================================================================
Create colormaps from current RGB triplets in colormaps directory
==============================================================================#

create_cmap("bluegie",28)
create_cmap("carbone",56)
create_cmap("carbone2",56)
create_cmap("pyro",40)
create_cmap("pyro_div",40)
create_cmap("radar",56)
create_cmap("radar2",56)
create_cmap("ssec",32)
create_cmap("ssec2",28)

#==============================================================================
Manually create colormaps 
==============================================================================#

radar3 = PyPlot.ColorMap("radar3",
[(0.0,0.0000,0.0000),(0.1,0.0000,0.0000),(0.2,0.03529,0.03529),(0.3,0.03137,0.03137),(0.5,1.00000,1.00000),(0.6,1.00000,1.00000),(0.7,1.0,1.0),(0.8,0.73725,0.73725),(0.9,0.47450,0.47450),(1.0,0.76470,0.76470)],
[(0.0,0.6157,0.6157),(0.1,0.0000,0.0000),(0.2,0.50980,0.50980),(0.3,0.68627,0.68627),(0.5,0.83921,0.83921),(0.6,0.59607,0.59607),(0.7,0.0,0.0),(0.8,0.00000,0.00000),(0.9,0.00000,0.00000),(1.0,0.63921,0.63921)],
[(0.0,1.0000,1.0000),(0.1,1.0000,1.0000),(0.2,0.68627,0.68627),(0.3,0.07843,0.07843),(0.5,0.00000,0.00000),(0.6,0.00000,0.00000),(0.7,0.0,0.0),(0.8,0.21176,0.21176),(0.9,0.42745,0.42745),(1.0,0.83137,0.83137)],56,1.0)
PyPlot.register_cmap("radar3", radar3)

