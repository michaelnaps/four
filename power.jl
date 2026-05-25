
include( "complex.jl" )

using StatsBase

# Structure and helper variables for the power spectrum.
mutable struct PowerSpectrum
    # Mate the power spectrum to its Fourier transform.
    F::Union{Fourier,ComplexFourier}

    # Power spectrum variables.
    P::Matrix{defFloat}
    R::Vector{defFloat}
    R̂::defFloat

    # Variable that sorts the spectrum.
    i::Vector{defInt}
end

function PowerSpectrum(F::Union{Fourier,ComplexFourier})::PowerSpectrum
    # Compute the power spectrum for sin/cos individually.
    N = length( F )

    if typeof( F ) == Fourier
        P = (1/2).*[F.A.^2 F.B.^2]  # Parseval normalization.
    elseif typeof( F ) == ComplexFourier
        P = hcat( (1/2).*real.( F.C.*conj.( F.C ) ) )
    end

    # Sum the power spectrum and sort large → small.
    R = vec( sum( P; dims=2 ) )
    i = sortperm( R; rev=true )

    # Return structure.
    return PowerSpectrum( F, P, R, sum( R ), i )
end

function Base.max(P::PowerSpectrum)::defFloat
    return P.F.ω[P.i[1]]
end

function StatsBase.mean(P::PowerSpectrum)::defFloat
    return sum( P.R.*P.F.ω )./P.R̂
end

function Plots.plot!(plt::Plots.Plot, P::PowerSpectrum; norm::Bool=true, args...)::Plots.Plot
    R̂ = norm ? P.R̂ : 1.0
    return plot!( plt, P.F.ω, P.R./R̂; args... )
end

function Plots.plot!(plt::Plots.Plot, P::PowerSpectrum; norm::Bool=true, args...)::Plots.Plot
    R̂ = norm ? P.R̂ : 1.0
    return plot!( plt, P.F.ω, P.R./R̂; args... )
end

println( "Loaded power.jl class file." )
