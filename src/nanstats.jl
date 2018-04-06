# *****************************************************************************
# nanstats.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.6.0
#
# Statistical functions that properly handle the presence of NaNs.
# Re-defining nanstats functions to more generalized applications.
# ** Note: NaNMath package was not the solution given that you can't specify 
#         a region
#
# Adapted from Tamas_papp @: https://discourse.julialang.org/t/nanmean-options/4994
# 
# Function list:
# nansum
# nanmean
# nanvar
# nanstd
# nanmin
# nanmax
# nanextrema
# Working on matrix multiplication and covariance in presence of NaNs
# *****************************************************************************

#==============================================================================
filtnan

Generic function to remove NaNs from an array using filter
=============================================================================#

filtnan(x::AbstractArray{<:Real}) = filter(!isnan,x)

#==============================================================================
nansum

Compute the sum of an array excluding NaNs with the option of specifying the
region of a multi-dimensional array
==============================================================================#

# General nansum function over entire array x

function nansum(x::AbstractArray{<:Real})
    return length(filtnan(x)) == 0 ? NaN : sum(filtnan(x))
end

# Apply nansum over a specific region in array x and option to squeeze the 
# output (sout) -- default behavior is the same as Base.sum

function nansum(x::AbstractArray{<:Real},region::Int64,sout::Bool=false)
 
    # Compute the sum over the desired region  
    xsum = mapslices(nansum,x,region)
    # Squeeze singleton dim if sout = true
    if sout
        return squeeze(xsum,findin(size(xsum),1)[1])  
    else
        return xsum
    end
end

#==============================================================================
nanmean

Compute the mean of an array excluding NaNs with the option of specifying the
region of a multi-dimensional array
==============================================================================#

# General nanmean function over entire array x

function nanmean(x::AbstractArray{<:Real})
    return mean(filtnan(x))
end

# Apply nanmean over a specific region in array x and option to squeeze the 
# output (sout) -- default behavior is the same as Base.mean

function nanmean(x::AbstractArray{<:Real},region::Int64,sout::Bool=false)
 
    # Compute the mean over the desired region  
    xmean = mapslices(nanmean,x,region)
    # Squeeze singleton dim if sout = true
    if sout
        return squeeze(xmean,findin(size(xmean),1)[1])  
    else
        return xmean
    end
end

#==============================================================================
nanvar

Compute the variance of an array excluding NaNs with the option of 
specifying the region of a multi-dimensional array
** Note: Using the default call to Base.var which uses a bias corrected 
         estimator (1/N-1)
==============================================================================#

# General nanvar function over entire array x

function nanvar(x::AbstractArray{<:Real})
    return var(filtnan(x))
end

# Apply nanvar over a specific region in array x and option to squeeze the 
# output (sout) -- default behavior is the same as Base.var

function nanvar(x::AbstractArray{<:Real},region::Int64,sout::Bool=false)
 
    # Compute the var over the desired region  
    xvar = mapslices(nanvar,x,region)
    # Squeeze singleton dim if sout = true
    if sout
        return squeeze(xvar,findin(size(xvar),1)[1])  
    else
        return xvar
    end
end

#==============================================================================
nanstd

Compute the standard deviation of an array excluding NaNs with the option of 
specifying the region of a multi-dimensional array
** Note: Using the default call to Base.std which uses a bias corrected 
         estimator (1/N-1)
==============================================================================#

# General nanstd function over entire array x

function nanstd(x::AbstractArray{<:Real})
    return std(filtnan(x))
end

# Apply nanstd over a specific region in array x and option to squeeze the 
# output (sout) -- default behavior is the same as Base.std

function nanstd(x::AbstractArray{<:Real},region::Int64,sout::Bool=false)
 
    # Compute the std over the desired region  
    xstd = mapslices(nanstd,x,region)
    # Squeeze singleton dim if sout = true
    if sout
        return squeeze(xstd,findin(size(xstd),1)[1])  
    else
        return xstd
    end
end

#==============================================================================
nanmin

Compute the minimum value of an array excluding NaNs with the option of 
specifying the region of a multi-dimensional array
==============================================================================#

# General nanmin function over entire array x

function nanmin(x::AbstractArray{<:Real})
    return minimum(filtnan(x))
end

# Apply nanmin over a specific region in array x 

function nanmin(x::AbstractArray{<:Real},region::Int64)
 
    # Compute the min over the desired region  
    xmin = mapslices(nanmin,x,region)
    return xmin
end

#==============================================================================
nanmax

Compute the maximum value of an array excluding NaNs with the option of 
specifying the region of a multi-dimensional array
==============================================================================#

# General nanmax function over entire array x

function nanmax(x::AbstractArray{<:Real})
    return maximum(filtnan(x))
end

# Apply nanmax over a specific region in array x 

function nanmax(x::AbstractArray{<:Real},region::Int64)
 
    # Compute the min over the desired region  
    xmax = mapslices(nanmax,x,region)
    return xmax
end

#==============================================================================
nanextrema

Compute the extrema of an array excluding NaNs with the option of 
specifying the region of a multi-dimensional array
==============================================================================#

# General nanextrema function over entire array x

function nanextrema(x::AbstractArray{<:Real})
    return extrema(filtnan(x))
end

# Apply nanextrema over a specific region in array x 

function nanextrema(x::AbstractArray{<:Real},region::Int64)
 
    # Compute the extrema over the desired region  
    xextrema = mapslices(nanextrema,x,region)
    return xextrema
end




