
include( "complex.jl" )

# Wavelet functions.
function morlet(ω::defFloat; ω0::defFloat=5.0)::defFloat
    return (1/π)^(0.25)*exp( -0.5*(ω - ω0)^2 )
end

# Scaling transform.
function scale(s::defFloat, ω::defFloat; ω0::defFloat=0.6, wavelet::Function=morlet)::defFloat
    return wavelet( s*ω; ω0=ω0 )*√s
end

# Class functions.
# mutable struct Wavelet
#     # Fourier modes.
#     M::FourierModes
# end
