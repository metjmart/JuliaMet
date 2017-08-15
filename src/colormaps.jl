# *****************************************************************************
# colormaps.jl
#
# Author: Jonathan Martinez 
# Email: jon.martinez@colostate.edu
# Julia version: 0.5.2
#
# -- Adapted/updated from a previous version by Annette Foerster --
#
# Several functions designed to work with color bars in Julia. To create a 
# colorbar manually, see the method introduced at the end.
# *****************************************************************************

function create_colormap(colorCodes::Array,mapname)
    dim1 = length(colorCodes)
    rs  = (Float64,Float64,Float64)[]
    gs  = (Float64,Float64,Float64)[]
    bs  = (Float64,Float64,Float64)[]

    for i in [1:1:dim1]
        push!(rs,((i-1)/(dim1-1),Colors.color(colorCodes[i]).r,Colors.color(colorCodes[i]).r))
        push!(gs,((i-1)/(dim1-1),Colors.color(colorCodes[i]).g,Colors.color(colorCodes[i]).g))
        push!(bs,((i-1)/(dim1-1),Colors.color(colorCodes[i]).b,Colors.color(colorCodes[i]).b))
    end

    cmap = plt.ColorMap(mapname,rs,gs,bs,dim1,1.0)

    plt.register_cmap(mapname, cmap)
    return true
end

function create_colormap(reds::Array,greens::Array,blues::Array,mapname)
    dim1 = length(reds)
    rs  = (Float64,Float64,Float64)[]
    gs  = (Float64,Float64,Float64)[]
    bs  = (Float64,Float64,Float64)[]

    for i in [1:1:dim1]
        push!(rs,((i-1)/(dim1-1),reds[i]/255,reds[i]/255))
        push!(gs,((i-1)/(dim1-1),greens[i]/255,greens[i]/255))
        push!(bs,((i-1)/(dim1-1),blues[i]/255,blues[i]/255))
    end

    cmap = plt.ColorMap(mapname,rs,gs,bs,dim1,1.0)

    plt.register_cmap(mapname, cmap)
    return true
end

#==============================================================================
extract_rgbs

This function will extract a specified number of RGB triplets from a colorbar 
that was already created and stored in a matplotlib colormap.

cbar_name - string defining the colorbar name defined by the colormap
inc - number of increments (i.e., number of RGB triplets)
==============================================================================#

function extract_rgbs(cbar_name::AbstractString,inc::Real)

    julia_cmap = matplotlib[:cm][:get_cmap](cbar_name)
    num_colors = Array(Float64,length(1:inc))

    for i in eachindex(num_colors)
        num_colors[i] = i/inc
    end

    rgbs = []
    rgb_vals = []

    for i in eachindex(num_colors)
        push!(rgbs,julia_cmap(num_colors[i]))
        push!(rgb_vals,collect(rgbs[i][1:3]))
        println(round(rgb_vals[i].*255))
    end

    return rgb_vals

end

#==============================================================================
# Extract the rgb values from the format that colormaps are stored in
==============================================================================#

function get_rgbs(r_in,g_in,b_in)
    dim1 = length(r_in)
    rgb_out = Array[]
    for i in [1:1:dim1]
        push!(rgb_out,[r_in[i][2],g_in[i][2],b_in[i][2]])
    end
    return rgb_out
end

#==============================================================================
# Use the colors of the colormap to color the lines of your lineplot
==============================================================================#

function make_line_colors(colors_in,nrOfColors::Number)
    colors_out = typeof(colors_in[1])[]
    repeat = ceil(nrOfColors/size(colors_in)[1])
    for i in [1:1:size(colors_in)[1]]
        for m in [1:1:repeat]
            push!(colors_out,colors_in[i])
        end
    end
    return colors_out
end

#==============================================================================
# Specify different colormaps to be created
#
# This function will create a Julia colorbar based on specified colors and 
# gradients. The required structure is as follows:
#
# plt.ColorMap("cbarname,
# [(num,grad,rs)],
# [(num,grad,gs)],
# [(num,grad,bs)], interval, 1.0)
#
# num - specifies the fractional number of the color based on the total 
#       number of colors (e.g., 0.0, 0.1, .... 1.0 for an 11-color colorbar)
# grad - specifies the gradient to be used for the color. Typically set 
#        this to the decimal RGB value described below
# rs, gs, bs - specify the decimal RGB values for each respective color
#              For example, for a red hue of 125, rs = 125/255 = 0.4902
# interval - specifies the number of colors to be placed on the cbar
# 1.0 - (not sure what this is for? always set as 1.0...)
#
# plt.register_cmap("cbarname", cbarobject)
#
# "cbarname" - specify the name of your colorbar
# "cbarobject - specify the cbar object you created from plt.ColorMap
#
==============================================================================#

kwr_cmap = plt.ColorMap("kwr",
[(0.0,0.0,0.0),(0.42,1.0,1.0),(0.58,1.0,1.0),(1.0,0.7,0.7)],
[(0.0,0.0,0.0),(0.42,1.0,1.0),(0.58,1.0,1.0),(1.0,0.0,1.0)],
[(0.0,0.0,0.0),(0.42,1.0,1.0),(0.58,1.0,1.0),(1.0,0.0,0.0)],33,1.0)
plt.register_cmap("kwr", kwr_cmap)

kwk_cmap = plt.ColorMap("kwk",
[(0.0,0.0,0.0),(0.42,1.0,1.0),(0.58,1.0,1.0),(1.0,0.0,0.0)],
[(0.0,0.0,0.0),(0.42,1.0,1.0),(0.58,1.0,1.0),(1.0,0.0,0.0)],
[(0.0,0.0,0.0),(0.42,1.0,1.0),(0.58,1.0,1.0),(1.0,0.0,0.0)],33,1.0)
plt.register_cmap("kwk", kwk_cmap)

blue_cmap = plt.ColorMap("blue",
[(0.0,0.0,0.0),(0.9,1.0,1.0),(1.0,1.0,1.0)],
[(0.0,0.0,0.0),(0.9,1.0,1.0),(1.0,1.0,1.0)],
[(0.0,0.7,0.7),(0.9,1.0,1.0),(1.0,1.0,1.0)],33,1.0)
plt.register_cmap("blue", blue_cmap)

blue_r_cmap = plt.ColorMap("blue_r",
[(0.0,1.0,1.0),(0.1,1.0,1.0),(1.0,0.0,0.0)],
[(0.0,1.0,1.0),(0.1,1.0,1.0),(1.0,0.0,0.0)],
[(0.0,1.0,1.0),(0.1,1.0,1.0),(1.0,0.7,0.7)],33,1.0)
plt.register_cmap("blue_r", blue_r_cmap)

red_cmap = plt.ColorMap("red",
[(0.0,0.7,0.7),(0.9,1.0,1.0),(1.0,1.0,1.0)],
[(0.0,0.0,0.0),(0.9,1.0,1.0),(1.0,1.0,1.0)],
[(0.0,0.0,0.0),(0.9,1.0,1.0),(1.0,1.0,1.0)],33,1.0)
plt.register_cmap("red", red_cmap)

red_r_cmap = plt.ColorMap("red_r",
[(0.0,1.0,1.0),(0.1,1.0,1.0),(1.0,0.7,0.7)],
[(0.0,1.0,1.0),(0.1,1.0,1.0),(1.0,0.0,0.0)],
[(0.0,1.0,1.0),(0.1,1.0,1.0),(1.0,0.0,0.0)],33,1.0)
plt.register_cmap("red_r", red_r_cmap)

radar_cmap = plt.ColorMap("radar",
[(0.0,0.0,0.0),(0.09,0.0000,0.0000),(0.17,0.0000,0.0000),(0.25,0.03529,0.03529),(0.33,0.03137,0.03137),(0.41,1.00000,1.00000),(0.49,1.00000,1.00000),(0.57,1.0,1.0),(0.65,0.86666,0.86666),(0.73,0.73725,0.73725),(0.81,0.47450,0.47450),(0.89,0.47450,0.47450),(1.0,0.76470,0.76470)],
[(0.0,0.9,0.9),(0.09,0.6157,0.6157),(0.17,0.0000,0.0000),(0.25,0.50980,0.50980),(0.33,0.68627,0.68627),(0.41,0.83921,0.83921),(0.49,0.59607,0.59607),(0.57,0.0,0.0),(0.65,0.00000,0.00000),(0.73,0.00000,0.00000),(0.81,0.00000,0.00000),(0.89,0.20000,0.20000),(1.0,0.63921,0.63921)],
[(0.0,1.0,1.0),(0.09,1.0000,1.0000),(0.17,1.0000,1.0000),(0.25,0.68627,0.68627),(0.33,0.07843,0.07843),(0.41,0.00000,0.00000),(0.49,0.00000,0.00000),(0.57,0.0,0.0),(0.65,0.10588,0.10588),(0.73,0.21176,0.21176),(0.81,0.42745,0.42745),(0.89,0.62745,0.62745),(1.0,0.83137,0.83137)],56,1.0)
plt.register_cmap("radar", radar_cmap)

radar2_cmap = plt.ColorMap("radar2",
[(0.0,0.0,0.0),(0.1,0.0000,0.0000),(0.2,0.0000,0.0000),(0.3,0.03529,0.03529),(0.4,0.03137,0.03137),(0.5,1.00000,1.00000),(0.6,1.00000,1.00000),(0.7,1.0,1.0),(0.8,0.73725,0.73725),(0.9,0.47450,0.47450),(1.0,0.76470,0.76470)],
[(0.0,0.9,0.9),(0.1,0.6157,0.6157),(0.2,0.0000,0.0000),(0.3,0.50980,0.50980),(0.4,0.68627,0.68627),(0.5,0.83921,0.83921),(0.6,0.59607,0.59607),(0.7,0.0,0.0),(0.8,0.00000,0.00000),(0.9,0.00000,0.00000),(1.0,0.63921,0.63921)],
[(0.0,1.0,1.0),(0.1,1.0000,1.0000),(0.2,1.0000,1.0000),(0.3,0.68627,0.68627),(0.4,0.07843,0.07843),(0.5,0.00000,0.00000),(0.6,0.00000,0.00000),(0.7,0.0,0.0),(0.8,0.21176,0.21176),(0.9,0.42745,0.42745),(1.0,0.83137,0.83137)],56,1.0)
plt.register_cmap("radar2", radar2_cmap)

radar3_cmap = plt.ColorMap("radar3",
[(0.0,0.0000,0.0000),(0.1,0.0000,0.0000),(0.2,0.03529,0.03529),(0.3,0.03137,0.03137),(0.5,1.00000,1.00000),(0.6,1.00000,1.00000),(0.7,1.0,1.0),(0.8,0.73725,0.73725),(0.9,0.47450,0.47450),(1.0,0.76470,0.76470)],
[(0.0,0.6157,0.6157),(0.1,0.0000,0.0000),(0.2,0.50980,0.50980),(0.3,0.68627,0.68627),(0.5,0.83921,0.83921),(0.6,0.59607,0.59607),(0.7,0.0,0.0),(0.8,0.00000,0.00000),(0.9,0.00000,0.00000),(1.0,0.63921,0.63921)],
[(0.0,1.0000,1.0000),(0.1,1.0000,1.0000),(0.2,0.68627,0.68627),(0.3,0.07843,0.07843),(0.5,0.00000,0.00000),(0.6,0.00000,0.00000),(0.7,0.0,0.0),(0.8,0.21176,0.21176),(0.9,0.42745,0.42745),(1.0,0.83137,0.83137)],56,1.0)
plt.register_cmap("radar3", radar3_cmap)

radar_carbone_cmap = plt.ColorMap("radar_carbone",
[(0.0,0.46667,0.46667),(0.1,0.29412,0.29412),(0.2,0.25882,0.25882),(0.3,0.01961,0.01961),(0.4,0.44706,0.44706),(0.5,0.87059,0.87059),(0.6,0.93725,0.93725),(0.7,0.77647,0.77647),(0.8,0.545098,0.545098),(0.9,0.70980,0.70980),(1.0,1.00000,1.00000)],
[(0.0,0.01961,0.01961),(0.1,0.06275,0.06275),(0.2,0.57647,0.57647),(0.3,0.55294,0.55294),(0.4,0.74118,0.74118),(0.5,0.89019,0.89019),(0.6,0.73725,0.73725),(0.7,0.54509,0.54905),(0.8,0.34901,0.34901),(0.9,0.07451,0.07451),(1.0,0.011765,0.011765)],
[(0.0,0.63921,0.63921),(0.1,0.77647,0.77647),(0.2,0.55294,0.55294),(0.3,0.01961,0.01961),(0.4,0.44706,0.44706),(0.5,0.75294,0.75294),(0.6,0.03921,0.03921),(0.7,0.16863,0.16863),(0.8,0.24705,0.24705),(0.9,0.16471,0.16471),(1.0,0.00000,0.00000)],56,1.0)
plt.register_cmap("radar_carbone", radar_carbone_cmap)

radar_carbone2_cmap = plt.ColorMap("radar_carbone2",
[(0.0,0.46667,0.46667),(0.1,0.29412,0.29412),(0.2,0.25882,0.25882),(0.3,0.01961,0.01961),(0.4,0.44706,0.44706),(0.5,0.87059,0.87059),(0.6,0.93725,0.93725),(0.7,0.77647,0.77647),(0.8,0.545098,0.545098),(0.9,0.72549,0.72549),(1.0,0.88235,0.88235)],
[(0.0,0.01961,0.01961),(0.1,0.06275,0.06275),(0.2,0.57647,0.57647),(0.3,0.55294,0.55294),(0.4,0.74118,0.74118),(0.5,0.89019,0.89019),(0.6,0.73725,0.73725),(0.7,0.54509,0.54905),(0.8,0.34901,0.34901),(0.9,0.09019,0.09019),(1.0,0.211765,0.211765)],
[(0.0,0.63921,0.63921),(0.1,0.77647,0.77647),(0.2,0.55294,0.55294),(0.3,0.01961,0.01961),(0.4,0.44706,0.44706),(0.5,0.75294,0.75294),(0.6,0.03921,0.03921),(0.7,0.16863,0.16863),(0.8,0.24705,0.24705),(0.9,0.19216,0.19216),(1.0,0.35294,0.35294)],56,1.0)
plt.register_cmap("radar_carbone2", radar_carbone2_cmap)

pvwcolors_cmap = plt.ColorMap("pvwcolors",
[(0.0,0.6314,0.6314),(0.0769,0.6784,0.6784),(0.1538,0.1764,0.1764),(0.2307,0.3529,0.3529),(0.3076,0.5294,0.5294),(0.3845,0.7059,0.7059),(0.4614,1.000,1.000),(0.5383,1.000,1.000),(0.6152,1.0000,1.0000),(0.6921,1.0000,1.0000),(0.7690,1.0000,1.0000),(0.8459,1.0000,1.0000),(0.9228,1.0000,1.0000),(1.0,0.9647,0.9647)],
[(0.0,0.1490,0.1490),(0.0769,0.2588,0.2588),(0.1538,0.2745,0.2745),(0.2307,0.4274,0.4274),(0.3076,0.5451,0.5451),(0.3845,0.7412,0.7412),(0.4614,1.000,1.000),(0.5383,1.000,1.000),(0.6152,0.9333,0.9333),(0.6921,0.8118,0.8118),(0.7690,0.6890,0.6890),(0.8459,0.4549,0.4549),(0.9228,0.1725,0.1725),(1.0,0.0000,0.0000)],
[(0.0,1.0000,1.0000),(0.0769,1.0000,1.0000),(0.1538,1.0000,1.0000),(0.2307,1.0000,1.0000),(0.3076,1.0000,1.0000),(0.3845,1.0000,1.0000),(0.4614,1.000,1.000),(0.5383,1.000,1.000),(0.6152,0.3882,0.3882),(0.6921,0.0549,0.0549),(0.7690,0.0549,0.0549),(0.8459,0.0549,0.0549),(0.9228,0.1529,0.1529),(1.0,0.2706,0.2706)],28,1.0)
plt.register_cmap("pvwcolors",pvwcolors_cmap)

pvcolors_cmap = plt.ColorMap("pvcolors",
[(0.0,0.6314,0.6314),(0.0909,0.6784,0.6784),(0.1818,0.1764,0.1764),(0.2727,0.3529,0.3529),(0.3636,0.5294,0.5294),(0.4545,0.7059,0.7059),(0.5454,1.0000,1.0000),(0.6363,1.0000,1.0000),(0.7272,1.0000,1.0000),(0.8181,1.0000,1.0000),(0.9090,1.0000,1.0000),(1.0,0.9647,0.9647)],
[(0.0,0.1490,0.1490),(0.0909,0.2588,0.2588),(0.1818,0.2745,0.2745),(0.2727,0.4274,0.4274),(0.3636,0.5451,0.5451),(0.4545,0.7412,0.7412),(0.5454,0.9333,0.9333),(0.6363,0.8118,0.8118),(0.7272,0.6890,0.6890),(0.8181,0.4549,0.4549),(0.9090,0.1725,0.1725),(1.0,0.0000,0.0000)],
[(0.0,1.0000,1.0000),(0.0909,1.0000,1.0000),(0.1818,1.0000,1.0000),(0.2727,1.0000,1.0000),(0.3636,1.0000,1.0000),(0.4545,1.0000,1.0000),(0.5454,0.3882,0.3882),(0.6363,0.0549,0.0549),(0.7272,0.0549,0.0549),(0.8181,0.0549,0.0549),(0.9090,0.1529,0.1529),(1.0,0.2706,0.2706)],28,1.0)
plt.register_cmap("pvcolors",pvcolors_cmap)

ssec_cmap = plt.ColorMap("ssec_colors",
[(0.0000,0.0000,0.0000),(0.0667,0.0000,0.0000),(0.1333,0.0000,0.0000),(0.2000,0.0000,0.0000),(0.2667,0.0000,0.0000),(0.3333,0.0000,0.0000),(0.4000,0.3059,0.3059),(0.4667,0.5843,0.5843),(0.5333,0.8588,0.8588),(0.6000,1.0000,1.0000),(0.6667,1.0000,1.0000),(0.7333,1.0000,1.0000),(0.8000,1.0000,1.0000),(0.8667,0.8980,0.8980),(0.9333,0.7921,0.7921),(1.0000,0.6980,0.6980)],
[(0.0000,0.0353,0.0353),(0.0667,0.2000,0.2000),(0.1333,0.3647,0.3647),(0.2000,0.5333,0.5333),(0.2667,0.6980,0.6980),(0.3333,0.8627,0.8627),(0.4000,0.9294,0.9294),(0.4667,0.9569,0.9569),(0.5333,0.9843,0.9843),(0.6000,0.9255,0.9255),(0.6667,0.7765,0.7765),(0.7333,0.6275,0.6275),(0.8000,0.4784,0.4784),(0.8667,0.3490,0.3490),(0.9333,0.2235,0.2235),(1.0000,0.1137,0.1137)],
[(0.0000,0.2039,0.2039),(0.0667,0.3373,0.3373),(0.1333,0.4706,0.4706),(0.2000,0.6039,0.6039),(0.2667,0.7373,0.7373),(0.3333,0.8706,0.8706),(0.4000,0.6196,0.6196),(0.4667,0.3725,0.3725),(0.5333,0.1216,0.1216),(0.6000,0.0000,0.0000),(0.6667,0.0000,0.0000),(0.7333,0.0000,0.0000),(0.8000,0.0000,0.0000),(0.8667,0.0000,0.0000),(0.9333,0.0000,0.0000),(1.0000,0.0000,0.0000)],32.0,1.0)
plt.register_cmap("ssec_colors",ssec_cmap)

bluegie_cmap = plt.ColorMap("bluegie",
[(0.0000,0.0000,0.0000),(0.0769,0.0000,0.0000),(0.1538,0.0000,0.0000),(0.2308,0.0000,0.0000),(0.3077,0.0000,0.0000),(0.3846,0.0000,0.0000),(0.4615,0.0000,0.0000),(0.5385,0.0000,0.0000),(0.6154,0.1686,0.1686),(0.6923,0.3059,0.3059),(0.7692,0.4471,0.4471),(0.8462,0.5843,0.5843),(0.9231,0.7569,0.7569),(1.0000,1.0000,1.0000)],
[(0.0000,0.3530,0.0353),(0.0769,0.1176,0.1176),(0.1538,0.2000,0.2000),(0.2308,0.2824,0.2824),(0.3077,0.3647,0.3647),(0.3846,0.4471,0.4471),(0.4615,0.6157,0.6157),(0.5385,0.8627,0.8627),(0.6154,0.9176,0.9176),(0.6923,0.9294,0.9294),(0.7692,0.9059,0.9059),(0.8462,0.9569,0.9569),(0.9231,0.9725,0.9725),(1.0000,0.9608,0.9608)],
[(0.0000,0.2039,0.2039),(0.0769,0.2706,0.2706),(0.1538,0.3373,0.3373),(0.2308,0.4039,0.4039),(0.3077,0.4706,0.4706),(0.3846,0.5373,0.5373),(0.4615,0.6706,0.6706),(0.5385,0.8706,0.8706),(0.6154,0.7451,0.7451),(0.6923,0.6196,0.6196),(0.7692,0.4980,0.4980),(0.8462,0.3725,0.3725),(0.9231,0.2157,0.2157),(1.0000,0.0000,0.0000)],28.0,1.0)
plt.register_cmap("bluegie",bluegie_cmap)

#=
Another possible radar colormap
colors = [[ .469, .020, .640 ], [ .164, .055, .582 ],
         [ .352, .141, .898 ], [ .445, .559, .996 ], [ .004, .445, .000 ],
         [ .059, .586, .059 ], [ .289, .680, .289 ], [ .633, .816, .633 ],
         [ .938, .906, .703 ], [ .938, .812, .000 ], [ .926, .672, .086 ],
         [ .816, .578, .148 ], [ .648, .438, .242 ], [ .485, .328, .297 ],
         [ .625, .003, .000 ], [ .879, .211, .355 ]]
=#
