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

filtnan(x::AbstractArray{T}) where T<:Real = filter(!isnan,x)

#==============================================================================
nansum

Compute the sum of an array excluding NaNs with the option of specifying the
region of a multi-dimensional array
==============================================================================#

# General nansum function over entire array x

function nansum(x::AbstractArray{T}) where T<:Real
    return length(filtnan(x)) == 0 ? NaN : sum(filtnan(x))
end

# Apply nansum over a specific region in array x and option to squeeze the
# output (sout) -- default behavior is the same as Base.sum

function nansum(x::AbstractArray{T},region::Int,sout::Bool=false) where T<:Real

    # Compute the sum over the desired region
    xsum = mapslices(nansum,x,dims=region)
    # Squeeze singleton dim if sout = true
    if sout
        return dropdims(xsum,dims=findmin(size(xsum))[2])
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

function nanmean(x::AbstractArray{T}) where T<:Real
    return length(filtnan(x)) == 0 ? NaN : mean(filtnan(x))
end

# Apply nanmean over a specific region in array x and option to squeeze the
# output (sout) -- default behavior is the same as Base.mean

function nanmean(x::AbstractArray{T},region::Int,sout::Bool=false) where T<:Real

    # Compute the mean over the desired region
    xmean = mapslices(nanmean,x,dims=region)
    # Squeeze singleton dim if sout = true
    if sout
        return dropdims(xmean,dims=findmin(size(xmean))[2])
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

function nanvar(x::AbstractArray{T}) where T<:Real
    return length(filtnan(x)) == 0 ? NaN : var(filtnan(x))
end

# Apply nanvar over a specific region in array x and option to squeeze the
# output (sout) -- default behavior is the same as Base.var

function nanvar(x::AbstractArray{T},region::Int,sout::Bool=false) where T<:Real

    # Compute the var over the desired region
    xvar = mapslices(nanvar,x,dims=region)
    # Squeeze singleton dim if sout = true
    if sout
        return dropdims(xvar,dims=findmin(size(xvar))[2])
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

function nanstd(x::AbstractArray{T}) where T<:Real
    return length(filtnan(x)) == 0 ? NaN : std(filtnan(x))
end

# Apply nanstd over a specific region in array x and option to squeeze the
# output (sout) -- default behavior is the same as Base.std

function nanstd(x::AbstractArray{T},region::Int,sout::Bool=false) where T<:Real

    # Compute the std over the desired region
    xstd = mapslices(nanstd,x,dims=region)
    # Squeeze singleton dim if sout = true
    if sout
        return squeeze(xstd,dims=findmin(size(xstd))[2])
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

function nanmin(x::AbstractArray{T}) where T<:Real
    return minimum(filtnan(x))
end

# Apply nanmin over a specific region in array x

function nanmin(x::AbstractArray{T},region::Int) where T<:Real

    # Compute the min over the desired region
    xmin = mapslices(nanmin,x,dims=region)
    return xmin
end

#==============================================================================
nanmax

Compute the maximum value of an array excluding NaNs with the option of
specifying the region of a multi-dimensional array
==============================================================================#

# General nanmax function over entire array x

function nanmax(x::AbstractArray{T}) where T<:Real
    return maximum(filtnan(x))
end

# Apply nanmax over a specific region in array x

function nanmax(x::AbstractArray{T},region::Int) where T<:Real

    # Compute the min over the desired region
    xmax = mapslices(nanmax,x,dims=region)
    return xmax
end

#==============================================================================
nanextrema

Compute the extrema of an array excluding NaNs with the option of
specifying the region of a multi-dimensional array
==============================================================================#

# General nanextrema function over entire array x

function nanextrema(x::AbstractArray{T}) where T<:Real
    return extrema(filtnan(x))
end

# Apply nanextrema over a specific region in array x

function nanextrema(x::AbstractArray{T},region::Int) where T<:Real

    # Compute the extrema over the desired region
    xextrema = mapslices(nanextrema,x,dims=region)
    return xextrema
end

#==============================================================================
nanfindmax

Find the maximum value in an array excluding NaNs and return both the value and
corresponding index as a tuple. Can handle multi-dimensional arrays.
This code was taken from Base and modified to allow the presence of NaNs.
* Currently does not support dimension specification.
==============================================================================#

function _nanfindmax(a::AbstractArray{T}, ::Colon) where T<:Real
    p = pairs(a)
    y = iterate(p)
    if y === nothing
        throw(ArgumentError("collection must be non-empty"))
    end
    (mi, m), s = y
    isnan(m) ? m = -Inf : nothing
    i = mi
    while true
        y = iterate(p, s)
        y === nothing && break
        (i, ai), s = y
        isnan(ai) ? ai = -Inf : nothing
        if isless(m, ai)
            m = ai
            mi = i
        end
    end
    return (m, mi)
end

nanfindmax(a) = _nanfindmax(a,:)

#==============================================================================
nanargmax

Find the index corresponding to the maximum value in an array, excluding NaNs.
Can handle multi-dimensional arrays.
This code was taken from Base and modified to allow the presence of NaNs.
* Currently does not support dimension specification.
==============================================================================#

nanargmax(a) = nanfindmax(a)[2]

#==============================================================================
nanfindmin

Find the minimum value in an array excluding NaNs and return both the value and
corresponding index as a tuple. Can handle multi-dimensional arrays.
This code was taken from Base and modified to allow the presence of NaNs.
* Currently does not support dimension specification.
==============================================================================#

function _nanfindmin(a::AbstractArray{T}, ::Colon) where T<:Real
    p = pairs(a)
    y = iterate(p)
    if y === nothing
        throw(ArgumentError("collection must be non-empty"))
    end
    (mi, m), s = y
    i = mi
    while true
        y = iterate(p, s)
        y === nothing && break
        (i, ai), s = y
        if isless(ai,m)
            m = ai
            mi = i
        end
    end
    return (m, mi)
end

nanfindmin(a) = _nanfindmin(a,:)

#==============================================================================
nanargmin

Find the index corresponding to the minimum value in an array, excluding NaNs.
Can handle multi-dimensional arrays.
This code was taken from Base and modified to allow the presence of NaNs.
* Currently does not support dimension specification.
==============================================================================#

nanargmin(a) = nanfindmin(a)[2]
