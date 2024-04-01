# *****************************************************************************
# nc_tools.jl
#
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# Functions for NetCDF I/O in Julia
#
# Function list
# read_ncvars
# *****************************************************************************

#==============================================================================
read_ncvars

Read in variables from a NetCDF file and squeeze singleton dimensions

ncfile = string pointing to NetCDF file location/name
varnames = array of strings with variable names (1-element array for 1 variable)

Keyword arguments

mask_opt = if true--replace fill/miss values with NaN if fill/miss values
           exist for the variable (opt for fillval first, then missval)

dict_opt = if true--return variable(s) in an ordered dictionary
==============================================================================#

function read_ncvars(ncfile::AbstractString,varnames::AbstractArray;
                     mask_opt::Bool,dict_opt::Bool)

    # Create an ordered dict for the vars
    vardata = OrderedDict()
    for var in varnames
        vararr = ncread(ncfile,var)
        # Determine if the variable has a single-dimension
        if ndims(vararr) > 1
            # If any, remove single-dimension from multi-dimensional array
            if any(size(vararr) .== 1)
                ddims = Tuple(findall(size(vararr) .== 1))
                vardata[var] = dropdims(vararr,dims=ddims) 
            else 
                vardata[var] = vararr
            end
        else 
            vardata[var] = vararr
        end
        # Determine if values need to be masked
        if mask_opt 
            # If yes, determine the fill/miss values for the variable
            fillval = ncgetatt(ncfile,var,"_FillValue")
            missval = ncgetatt(ncfile,var,"missing_value")
            # Only replace fill/miss values with NaN if fill/miss values
            # exist for the variable
            # Opt for fillval first, then missval
            if fillval !== nothing
                vardata[var][findall(fillval,vardata[var])] .= NaN
            elseif fillval === nothing && missval !== nothing
                vardata[var][findall(missval,vardata[var])] .= NaN
            end
        end
    end
    # Determine varsdata output type: OrderedDict or arrays
    if dict_opt 
        return vardata
    else
        return collect(values(vardata))
    end
end
