
include( "power.jl" )

# Wavelet functions.
function morlet(ω::defFloat; ω0::defFloat=5.0)::defFloat
    σ = (1 + exp( -ω0^2 ) - 2exp( -0.75*ω0^2 ) )^(-0.5)
    κ = exp( -0.5*ω0^2 )
    return σ*(1/π)^(0.25)*(exp( -0.5*(ω - ω0)^2 ) - κ*exp( -0.5*ω^2 ))
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
    # Return the wavelet variable.
    return Wavelet( F, φ )
end

function solve(W::Wavelet, s::defFloat;
    wavelet::Function=morlet, X::Vector{defFloat}=W.F.X)::Vector{defFloat}
    # Compute the scaled wavelet coefficients.
    Ĉ = scale( W.F, s; wavelet=wavelet )

    # Solve the wavelet coefficients.
    return solve( ComplexFourier( W.F.M, Ĉ, W.F.X, W.F.Y, W.F.ω ); X=X )
end

println( "Loaded wavelet.jl class file." )
