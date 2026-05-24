
include( "complex.jl" )

# Wavelet functions.
function morlet(ω::defFloat; ω0::defFloat=5.0)::defFloat
    return (1/π)^(0.25)*exp( -0.5*(ω - ω0)^2 )
end

# Scaling function.
function scale(F::ComplexFourier, s::defFloat; wavelet::Function=morlet)::Vector{defComp}
    # Compute the scaled wavelet in the frequency domain.
    φ = √s.*conj.( wavelet.( s.*F.ω ) )
    return F.C.*φ
end

mutable struct Wavelet
    # Fourier variable.
    F::ComplexFourier

    # Wavelet "mother" function.
    φ::Function
end

Base.broadcastable(W::Wavelet) = Ref( W )

function Wavelet(F::ComplexFourier; φ::Function=morlet)::Wavelet
    return Wavelet( F, φ )
end

function solve(W::Wavelet, s::defFloat; wavelet::Function=morlet)::Vector{defFloat}
    # Compute the scaled wavelet coefficients.
    Ĉ = scale( W.F, s; wavelet=wavelet )

    # Solve the wavelet coefficients.
    return solve( ComplexFourier( W.F.M, Ĉ, W.F.X, W.F.Y, W.F.ω ); X=W.F.X )
end

println( "Loaded wavelet.jl class file." )
