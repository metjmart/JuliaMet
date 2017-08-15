# *****************************************************************************
# samurai_obslocs.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.5.2
#
# -- Adapted from Ellie Delap -- 
#
# This script contains functions that will be utilized in determining and 
# plotting the location of SAMURAI observations. 
# 
# Function List:
# obslocs_rtz
# More to come...
# *****************************************************************************

#==============================================================================
obslocs_rtz

This function reads in the SAMURAI_QC_analysis.out file in cylindrical 
coordinates and stores the rtz location of each data type in a dictionary. The
dictionary keys are the data types (integers) and the subkeys are the names 
of the data types (e.g. obs, inverse_error, etc.). The function will output 
the dictionary and automatically create specified figures.
==============================================================================#
 
function obslocs_rtz(filein,name::AbstractString,year::AbstractString,
                     month::AbstractString,day::AbstractString,rmax,zmax,
                     data_out=false,plots_out=false)
    # Read in the data 
    dlmdata = readdlm(filein, '\t'); 
    # Define the dimension sizes
    nrows,ncols = size(dlmdata)
    # Import the data from each individual column and convert to dtype{Float64}
    obs = Array{Float64}(dlmdata[2:nrows,1])
    inverse_error = Array{Float64}(dlmdata[2:nrows,2])
    radius = Array{Float64}(dlmdata[2:nrows,3])
    theta = Array{Float64}(dlmdata[2:nrows,4])
    z = Array{Float64}(dlmdata[2:nrows,5])
    dtype = Array{Float64}(dlmdata[2:nrows,6])
    # Need to import time as a string and convert to type Int
    time_str = Array{AbstractString}(dlmdata[2:nrows,7])
    time = Array(Int64,size(time_str))
    for i in collect(1:length(time))
        time[i] = parse(Int,time_str[i][1:3]*time_str[i][5:6]*time_str[i][8:9])
    end  
    rhou = Array{Float64}(dlmdata[2:nrows,8])
    rhov = Array{Float64}(dlmdata[2:nrows,9])
    rhow = Array{Float64}(dlmdata[2:nrows,10])
    tprime = Array{Float64}(dlmdata[2:nrows,11])
    qvprime = Array{Float64}(dlmdata[2:nrows,12])
    rhoaprime = Array{Float64}(dlmdata[2:nrows,13])
    qr = Array{Float64}(dlmdata[2:nrows,14])
    analysis = Array{Float64}(dlmdata[2:nrows,15])
    background = Array{Float64}(dlmdata[2:nrows,16])
    # There's a blank 17th column not being used here
    # Create a dictionary that contains keys for the data dtypes 
    # values and subkeys for the names of the data dtypes
    data_dict=Dict()
    for i in collect(1:length(dtype))
    # Only create keys for each unique dtype value
        if !haskey(data_dict,dtype[i]) 
            # Create subkeys for the names of the data dtypes
            # and set them as open arrays 
            data_dict[dtype[i]] = Dict()
            data_dict[dtype[i]]["obs"] = Float64[]
            data_dict[dtype[i]]["inverse_error"] = Float64[]
            data_dict[dtype[i]]["radius"] = Float64[]
            data_dict[dtype[i]]["theta"] = Float64[]
            data_dict[dtype[i]]["z"] = Float64[]
            data_dict[dtype[i]]["time"] = Int64[]
            data_dict[dtype[i]]["rhou"] = Float64[]
            data_dict[dtype[i]]["rhov"] = Float64[]
            data_dict[dtype[i]]["rhow"] = Float64[]
            data_dict[dtype[i]]["tprime"] = Float64[]
            data_dict[dtype[i]]["qvprime"] = Float64[]
            data_dict[dtype[i]]["rhoaprime"] = Float64[]
            data_dict[dtype[i]]["qr"] = Float64[]
            data_dict[dtype[i]]["analysis"] = Float64[]
            data_dict[dtype[i]]["background"] = Float64[]
        end
        # Place the data values in their corresponding dictionaries by dtype and name
        append!(data_dict[dtype[i]]["obs"],[obs[i]]) 
        append!(data_dict[dtype[i]]["inverse_error"],[inverse_error[i]])
        append!(data_dict[dtype[i]]["radius"],[radius[i]])
        append!(data_dict[dtype[i]]["theta"],[theta[i]])
        append!(data_dict[dtype[i]]["z"],[z[i]])
        append!(data_dict[dtype[i]]["time"],[time[i]])
        append!(data_dict[dtype[i]]["rhou"],[rhou[i]])
        append!(data_dict[dtype[i]]["rhov"],[rhov[i]])
        append!(data_dict[dtype[i]]["rhow"],[rhow[i]])
        append!(data_dict[dtype[i]]["tprime"],[tprime[i]])
        append!(data_dict[dtype[i]]["qvprime"],[qvprime[i]])
        append!(data_dict[dtype[i]]["rhoaprime"],[rhoaprime[i]])
        append!(data_dict[dtype[i]]["qr"],[qr[i]])
        append!(data_dict[dtype[i]]["analysis"],[analysis[i]])
        append!(data_dict[dtype[i]]["background"],[background[i]])
    end
    # Return the data_dict if data_out is true 
    if data_out == true
        return data_dict  
    end
    # Create png figures if plots_out is true
    #if plots_out == true

        # Figure 1 - Dropsonde, flight level, sfmr in r-z
        # Plot the locations of dropsonde and flight level if both
        # are present. If both are not present, plot the data that is
        # ** Will eventually include SFMR 
        
        # Plot data in the r-z plane
        if haskey(data_dict,0) && haskey(data_dict,1)
            # Plot the dropsonde and flight level data 
            plt.scatter(data_dict[0]["radius"], data_dict[0]["z"], marker="o", color="red")
            plt.scatter(data_dict[1]["radius"], data_dict[1]["z"], marker="o", color="blue")
            # Set plot configuration
            ax = plt.axes()
            ax[:axis]([0,rmax,0,zmax])
            ax[:legend](["Dropsonde", "Flight Level"])
            ax[:set_title]("Non-radar observations")
            ax[:set_xlabel]("Radius (km)")
            ax[:set_ylabel]("Altitude (km)")
            # Save the figure
            plt.savefig("./" * name * year * month * day * "_dropsonde_fl_rz.png")
            println("Succesfully plotted dropsonde and flight level data in r-z plane!")
        elseif haskey(data_dict,0) && !haskey(data_dict,1)
            # Plot just the dropsonde data if flight level is missing
            plt.scatter(data_dict[0]["radius"], data_dict[0]["z"], marker="o", color="red")
            # Set plot configuration
            ax = plt.axes()
            ax[:axis]([0,rmax,0,zmax])
            ax[:legend](["Dropsonde"])
            ax[:set_title]("Non-radar observations")
            ax[:set_xlabel]("Radius (km)")
            ax[:set_ylabel]("Altitude (km)")
            plt.savefig("./" * name * year * month * day * "_dropsonde_rz.png")
            println("Succesfully plotted dropsonde data in r-z plane!")
        elseif !haskey(data_dict,0) && haskey(data_dict,1)
            # Plot just the flight level data if dropsondes are missing
            plt.scatter(data_dict[1]["radius"], data_dict[1]["z"], marker="o", color="blue")
            # Set plot configuration
            ax = plt.axes()
            ax[:axis]([0,rmax,0,zmax])
            ax[:legend](["Flight Level"])
            ax[:set_title]("Non-radar observations")
            ax[:set_xlabel]("Radius (km)")
            ax[:set_ylabel]("Altitude (km)")
            # Save the figure
            plt.savefig("./" * name * year * month * day * "_fl_rz.png")
            println("Succesfully plotted flight level data in r-z plane!")
        end
    
        # Figure 2 - Dropsonde, flight level, sfmr in r-theta
        # Plot the locations of dropsonde and flight level if both
        # are present. If both are not present, plot the data that is
        # ** Will eventually include SFMR 
        
        # Plot data in the r-theta plane
        if haskey(data_dict,0) && haskey(data_dict,1)
            # Plot the flight level and dropsonde data
            ax = plt.axes(polar="true")
            plt.scatter((deg2rad(data_dict[0]["theta"])-pi/2).*-1, data_dict[0]["radius"], marker="o", color="red")
            plt.scatter((deg2rad(data_dict[1]["theta"])-pi/2).*-1, data_dict[1]["radius"], marker="o", color="blue")
            # Set plot configuration
            ax[:axis]([0,360,0,rmax])
            ax[:set_theta_zero_location]("N")
            ax[:set_theta_direction](-1)
            ax[:set_rlabel_position](75)
            ax[:legend](["Dropsonde", "Flight Level"])
            ax[:set_title]("Non-radar observations",y=1.06)
            # Save the figure
            plt.savefig("./" * name * year * month * day * "_dropsonde_fl_rt.png")
            println("Succesfully plotted dropsonde and flight level data in r-theta plane!")
        elseif haskey(data_dict,0) && !haskey(data_dict,1)
            # Plot just the dropsonde data if flight level is missing
            ax = plt.axes(polar="true")
            plt.scatter((deg2rad(data_dict[0]["theta"])-pi/2).*-1, data_dict[0]["radius"], marker="o", color="red")
            ax[:axis]([0,360,0,rmax])
            # Set plot configuration
            ax[:axis]([0,360,0,rmax])
            ax[:set_theta_zero_location]("N")
            ax[:set_theta_direction](-1)
            ax[:set_rlabel_position](75)
            ax[:legend](["Dropsonde"])
            ax[:set_title]("Non-radar observations",y=1.06)
            plt.savefig("./" * name * year * month * day * "_dropsonde_rt.png")
            println("Succesfully plotted dropsonde data in r-theta plane!")
        elseif !haskey(data_dict,0) && haskey(data_dict,1)
            # Plot just the flight level data if dropsondes are missing
            ax = plt.axes(polar="true")
            plt.scatter((deg2rad(data_dict[1]["theta"])-pi/2).*-1, data_dict[1]["radius"], marker="o", color="blue")
            ax[:axis]([0,360,0,rmax])
            ax[:set_theta_zero_location]("N")
            ax[:set_theta_direction](-1)
            ax[:set_rlabel_position](75)
            ax[:legend](["Flight Level"])
            ax[:set_title]("Non-radar observations",y=1.06)
            plt.savefig("./" * name * year * month * day * "_fl_rt.png")
            println("Succesfully plotted flight level data in r-theta plane!")
        end
    
        # Figure 3 - Radar observations in the r-z plane
    
        # Plot the radar data in the r-z plane if it's present
        if haskey(data_dict,2) 
            plt.scatter(data_dict[2]["radius"], data_dict[2]["z"], marker=".", color="black")
            # Set plot configuration
            ax = plt.axes()
            ax[:axis]([0,rmax,0,zmax])
            ax[:set_title]("Radar observations")
            ax[:set_xlabel]("Radius (km)")
            ax[:set_ylabel]("Altitude (km)")
            # Save the figure
            plt.savefig("./" * name * year * month * day * "_radar_rz.png")
            println("Succesfully plotted radar data in r-z plane!")
        end
    
        # Figure 4 - Radar observations in r-theta plane
    
        # Plot the radar data in the r-theta plane if it's present
        if haskey(data_dict,2)
            ax = plt.axes(polar="true")
            plt.scatter((deg2rad(data_dict[2]["theta"])-pi/2).*-1, data_dict[2]["radius"], marker=".", color="black")
            ax[:axis]([0,360,0,rmax])
            # Set plot configuration
            ax[:axis]([0,360,0,rmax])
            ax[:set_theta_zero_location]("N")
            ax[:set_theta_direction](-1)
            ax[:set_rlabel_position](75)
            ax[:set_title]("Radar observations",y=1.06)
            plt.savefig("./" * name * year * month * day * "_radar_rt.png")
            println("Succesfully plotted radar data in r-theta plane!")
        end 

    #end
 
end

#==============================================================================
obslocs_xyz

This function reads in the SAMURAI_QC_analysis.out file in cartesian
coordinates and stores the xyz location of each data type in a dictionary. The
dictionary keys are the data types (integers) and the subkeys are the names 
of the data types (e.g. obs, inverse_error, etc.). The function will output 
the dictionary and automatically create specified figures.
==============================================================================#




