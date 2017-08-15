# *****************************************************************************
# nanstats.jl
#
# Author: Jonathan Martinez
# Email: jon.martinez@colostate.edu
# Julia version: 0.5.2
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

for(nam, func) = ((:nanmax, :max), (:nanmin, :min), (:nansum, :sum),
                  (:nanmean, :mean), (:nanmedian, :median), (:nanvar, :var),
                  (:nanstd, :std), (:nanskewness, :skewness),
                  (:nankurtosis, :kurtosis))
    @eval begin
        function ($nam)(var::Array)
            newvar = typeof(var[1])[]
            for i in collect(1:length(var))
                if ~isnan(var[i])
                    push!(newvar, var[i])
                end
            end

            ($func)(newvar)
        end
    end
end

