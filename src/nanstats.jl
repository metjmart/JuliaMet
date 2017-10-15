# *****************************************************************************
# nanstats.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.6.0
#
# Adapted from milktrader: https://github.com/JuliaStats/StatsBase.jl/issues/3
#
# Statistical functions that properly handle the presence of NaNs.
# 
# *****************************************************************************

#==============================================================================
nanstats

Define a list of statistical functions that can handle the presence of NaNs.
==============================================================================#

for(nam, func) = ((:nanmax, :maximum), (:nanmin, :minimum), (:nanextrema, :extrema),
                  (:nansum, :sum),(:nanmean, :mean), (:nanmedian, :median), 
                  (:nanvar, :var), (:nanstd, :std))

    @eval begin
        function ($nam)(arr::Array)
            newarr = typeof(arr[1])[]
            for i in collect(1:length(arr))
                if ~isnan(arr[i])
                    push!(newarr, arr[i])
                end
            end

            ($func)(newarr)
        end
    end
end

