# *****************************************************************************
# harmonics
#
# This script contains functions that are utilized to decompose the azimuth 
# of a TC into a harmonic series. This is done via ordinary least squares 
# regression where the basis functions are a series of sin and cos waves.
# The functions can also be generalized for other applications where Fourier
# decomposition is desired.
# *****************************************************************************

#==============================================================================
rmnan

For a given independent variable (x) and dependent variable (y), remove
all indices where y == NaN
Compresses the size of original arrays to only include indices without NaN
Ex: x = [1.,2.,3.,4.], y = [1.,2.,3.,NaN]
c,d = rmnan(x,y) returns c = [1.,2.,3.], d = [1.,2.,3.]
==============================================================================#

# Only remove NaNs from one input vector

function rmnan(y::AbstractVector{<:Real})

    # Construct new array absent NaN indices    
    newy = Float64[]
    for i in eachindex(y)
        if !isnan(y[i])
            append!(newy, y[i])
        end
    end
    return newy
end

# Remove NaNs from both x and y based on their presence in y

function rmnan(x::AbstractVector{<:Real},y::AbstractVector{<:Real})

    length(x) == length(y) ? nothing :
        throw(DimensionMismatch("Input vectors must have same length"))
    # Construct new arrays absent NaN indices    
    newx = Float64[]
    newy = Float64[]
    for i in eachindex(x)
        if !isnan(y[i])
            append!(newx, x[i])
            append!(newy, y[i])
        end
    end
    return newx, newy
end

#==============================================================================
fourierols

Fourier decomposition of azimuth with Ordinary Least Squares (OLS) regression.
The function begins by constructing the design matrix A that contains a linear
combination of M basis functions and N data points.
In Fourier analysis, the basis functions are given by the harmonic series:
cos(k*phi) + sin(k*phi) where k is the integer wavenumber and phi is the 
azimuth angle from 0 to 2pi.
The function returns the vector "a" which contains the amplitudes of the 
coefficients.
** Note: If input wavenumber is nk = 2, the function will return coefficients
   for wavenumbers 0-2 (i.e., 3 coefficients)
** Note that if the independent variable in your data is not azimuth but say 
   time, you can still use this function. Just input phi with a normalization 
   factor.
   Ex: You have a time series (t) that spans (0 <= t <= T), then
   phi = 2 * pi * t/T
==============================================================================#

function fourierols(nk::Int,phi::AbstractVector{<:Real},b::AbstractVector{<:Real})

    length(phi) == length(b) ? nothing :
        throw(DimensionMismatch("Input vectors must have same length"))
    # Ensure that no wavenumber higher than the Nyquist is permitted
    nk <= (length(b) - 1)/2 ? nothing : error("Wavenumbers exceed the Nyquist frequency")
    # Construct the design matrix with M basis functions
    M = nk * 2 + 1
    A = ones(length(phi),M)
    for j in 2:M
        for i in eachindex(phi)
            @fastmath @inbounds iseven(j) ? A[i,j] = cos.(j/2*phi[i]) : A[i,j] = sin.((j-1)/2*phi[i])
        end 
    end
    # If det(A^T * A) = 0, no unique solution -- throw error
    det(At_mul_B(A,A)) == 0 ? error("No unique solution exists") : nothing 
    # Solve the normal equations (A^T * A) * a = A^T * b 
    a = (At_mul_B(A,A)) \ (At_mul_B(A,b))
    return a
end

#==============================================================================
wavecoeffs

Given an input set of Fourier coefficients, reconstruct the various wavenubmers
of a decomposed field
==============================================================================#

function wavecoeffs(a::AbstractVector{<:Real},phi::AbstractVector{<:Real})

    isodd(length(a)) ? nothing : error("Number of coefficients should be odd 
        (e.g., 3 coefficients required for wavenumbers 0 and 1)")
    nk = div(length(a) + 1, 2)
    nkwaves = zeros(nk,length(phi))
    # Wavenumber-0 is a[1] (mean)
    nkwaves[1,:] = a[1]
    for j in 1:nk-1
        ind = j + 1
        @fastmath @inbounds nkwaves[ind,:] = a[j*2] .* cos.(j*phi) + a[j*2+1] .* sin.(j*phi)
    end
    return nkwaves
end

