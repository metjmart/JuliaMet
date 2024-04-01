# *****************************************************************************
# harmonics
# 
# Author:
#       Jonathan Martinez
#
# Julia version: 
#       1.0.0
#
# Functions to construct a Fourier series via ordinary least squares regression
# and retrieve the harmonics corresponding to each integer wavenumber 
#
# Function list
# fourierols
# regress_harmonics
# *****************************************************************************

#==============================================================================
fourierols(nk::Int,b::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}

Construct a Fourier series via Ordinary Least Squares (OLS) regression.
- Construct the design matrix A that contains a linear combination of M basis 
  functions and N data points
- The basis functions are given by the harmonic series cos(k*ω) + sin(k*ω), 
  where k is the integer wavenumber and ω is the angular frequency given by
  ω = 2*pi*i/N, for i = 0:N-1

Input
nk = largest integer wavenumber to retrieve
b = data vector

Output
a = vector containing the amplitudes of the Fourier coefficients for each 
    wavenumber in k = 0:nk (A_k = cos coefficient; B_k = sin coefficient)
    k_0 = first index (wavenumber 0)
    A_k = even-valued indices
    B_k = odd-valued indices
==============================================================================#

function fourierols(nk::Int,b::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}

    N = length(b)
    # Ensure that no wavenumber higher than the Nyquist is permitted
    nk <= (N-1)/2 ? nothing : error("Wavenumbers exceed the Nyquist frequency")
    # Construct the design matrix with M basis functions
    M = nk * 2 + 1
    A = ones(N,M)
    for j in 2:M
        for (ind,i) in enumerate(0:N-1)
            # Below is equivalent to 2*pi*(j/2)*i/N, where j/2 = wavenumber k 
            @fastmath @inbounds iseven(j) ? A[ind,j] = cos(j*pi*i/N) : A[ind,j] = sin((j-1)*pi*i/N) 
        end
    end
    # If det(A^T * A) = 0, no unique solution -- throw error
    det(transpose(A)*A) == 0 ? error("No unique solution exists") : nothing
    # Solve the normal equations (A^T * A) * a = A^T * b
    a = (transpose(A)*A) \ (transpose(A)*b)
    return a
end

#==============================================================================
regress_harmonics(a::AbstractVector{Ta},b::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}

Given an input set of Fourier coefficients from fourierols, retrieve the 
regression onto each harmonic 

Input 
a = vector containing the amplitudes of the Fourier coefficients for each 
    wavenumber in k = 0:nk (A_k = cos coefficient; B_k = sin coefficient)
    k_0 = first index (wavenumber 0)
    A_k = even-valued indices
    B_k = odd-valued indices
b = data vector

Output 
h = regression of b onto each harmonic 
==============================================================================#

function regress_harmonics(a::AbstractVector{Ta},b::AbstractVector{Tb}) where {Ta<:Real,Tb<:Real}

    N = length(b)
    nk = div(length(a)-1, 2)
    h = zeros(nk+1,N)
    # Wavenumber 0 is a[1] 
    h[1,:] .+= a[1]
    for k in 1:nk
        for (ind,i) in enumerate(0:N-1)
            @fastmath @inbounds h[k+1,ind] = a[k*2] * cos(k*2*pi*i/N) + a[k*2+1] * sin(k*2*pi*i/N) 
        end
    end
    return h
end
