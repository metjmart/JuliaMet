# *****************************************************************************
# nc_tools.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.4.5
#
# This script contains functions that will handle reading and writing NetCDF 
# files in Julia.
#
# Function List:
# read_ncvars
# *****************************************************************************

#==============================================================================
read_ncvars

This function is designed to read in specified vars from NetCDF files.
It will remove single-dimensions when necessary and replace fill values 
with NaNs if necessary. 

Format for function use:

Define varnames as a string or an array of strings.
Ex: varnames = "DBZ" -- Single var
Ex: varnames = ["x","y","altitude","DBZ"] -- Array of vars
==============================================================================#

function read_ncvars(ncfile::AbstractString,varnames::Array,mask_opt::Bool=true,
                     dict_opt::Bool=false)
    # Determine if one or more vars needs to be read in
    # If one var, type is string
    if typeof(varnames) <: AbstractString
        println("Reading in " * varnames * " ...")
        vardata = ncread(ncfile, varnames)
        if ndims(vardata) == 1
            # Determine if values need to be masked
            if mask_opt == true
                # If they do, determine the fill/miss values for the variable
                fillval = ncgetatt(ncfile,var,"_FillValue")
                missval = ncgetatt(ncfile,var,"missing_value")
                # Only replace fill/mask values with NaN if fill/mask values 
                # exist for the variable
                if typeof(fillval) != Void
                    vardata[findin(vardata,fillval)] = NaN
                elseif typeof(fillval) == Void && typeof(missval) != Void
                    vardata[findin(vardata,missval)] = NaN
                end
                println("Succesfully read in " * varnames * "!")
                return vardata
            else
                # If no masking, just return vardata
                println("Succesfully read in " * varnames * "!")
                return vardata
            end
        # If not 1-d, determine if vardata has a single-dimension
        elseif ndims(vardata) > 1        
            # Remove single-dimension from multi-dimensional array
            for i in collect(1:ndims(vardata))
                if size(vardata)[i] == 1
                    vardata = squeeze(vardata,i)
                end
            end
            # Determine if fill values need to be masked
            if mask_opt == true
                # If they do, determine the fill/mask values for the variable
                fillval = ncgetatt(ncfile,var,"_FillValue")
                missval = ncgetatt(ncfile,var,"missing_value")
                # Only replace fill/miss values with NaN if fill/miss values 
                # exist for the variable
                if typeof(fillval) != Void
                    vardata[findin(vardata,fillval)] = NaN
                elseif typeof(fillval) == Void && typeof(missval) != Void
                    vardata[findin(vardata,missval)] = NaN
                end
                println("Successfully read in " * varnames * "!")
                return vardata
            else
                # If no masking, just return svardata
                println("Successfully read in " * varnames * "!")
                return vardata
            end
        end 
    # If more than one var, type is list
    elseif typeof(varnames) <: Array
        # Create an ordered dict for the vars
        varsdata = OrderedDict()
        svarsdata = OrderedDict()
        for var in varnames
            println("Reading in " * var * " ...")
            varsdata[var] = ncread(ncfile,var)
            # Create a second variable to overwrite data that is squeezed
            svarsdata[var] = ncread(ncfile,var)
            # Determine if any var in varsdata has a single-dimension
            if ndims(varsdata[var]) > 1
                # Remove single-dimension from multi-dimensional array
                for i in collect(1:ndims(varsdata[var]))
                    if size(varsdata[var])[i] == 1
                        svarsdata[var] = squeeze(varsdata[var],i)
                    #else
                        #svarsdata[var] = varsdata[var] 
                    end
                end
            #else
                #svarsdata[var] = varsdata[var]
            end
            # Determine if values need to be masked 
            if mask_opt == true
                # If they do, determine the fill/miss values for the variable
                fillval = ncgetatt(ncfile,var,"_FillValue")
                missval = ncgetatt(ncfile,var,"missing_value")
                # Only replace fill/miss values with NaN if fill/miss values 
                # exist for the variable
                if typeof(fillval) != Void
                    svarsdata[var][findin(svarsdata[var],fillval)] = NaN
                elseif typeof(fillval) == Void && typeof(missval) != Void
                    svarsdata[var][findin(svarsdata[var],missval)] = NaN
                end
            else
                # If no masking, just re-store the data
                svarsdata[var] = svarsdata[var]
            end
        end
        # Determine varsdata output type: Values or OrderedDict 
        if dict_opt == true
            for var in varnames
                println("Successfully read in " * var * "!")
            end
            return svarsdata
        else
            for var in varnames
                println("Successfully read in " * var * "!")
            end 
            return collect(values(svarsdata))
        end
    else  
       println("Failed to read in var(s), be sure to define them as a string or an array of strings!")
    end 
end    




