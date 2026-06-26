
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

function PowerSpectrum(F::Fourier)::PowerSpectrum
    # Parseval-consistent power for real Fourier (cos/sin components).
    P = 0.5.*[F.A.^2 F.B.^2]

    # Compute the power and rank.
    R = vec( sum( P; dims=2 ) )
    i = sortperm( R; rev=true )

    # Return the power spectrum.
    return PowerSpectrum( F, P, R, sum( R ), i )
end


function PowerSpectrum(F::ComplexFourier; Δω::Real=π)::PowerSpectrum
    # Proper complex power spectrum
    P = (Δω/2).*abs2.( F.C )

    # Compute the power and rank.
    R = vec( sum( P; dims=2 ) )
    i = sortperm( R; rev=true )

    # Return the power spectrum.
    return PowerSpectrum( F, P, R, sum( R ), i )
end

function Base.max(P::PowerSpectrum; ι::Union{defInt,Vector{defInt}}=1)::Union{defFloat,Vector{defFloat}}
    return P.F.ω[P.i[ι]]
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
