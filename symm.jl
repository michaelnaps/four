
include( "four.jl" )

function perturb(F::Fourier, P::PowerSpectrum, φ::defFloat; K::defInt=length( F ))::Fourier
    # Copy the given Fourier series.
    F̂ = copy( F )

    # Perturb that K-top frequencies.
    for i ∈ P.i[1:K]
        # Unpack coefficients and frequency.
        a, b, ω = F.A[i], F.B[i], F.ω[i]

        # Translate frequencies.
        F̂.A[i] = a*cos( ω*φ ) - b*sin( ω*φ )
        F̂.B[i] = a*sin( ω*φ ) + b*cos( ω*φ )
        # NOTE: Removed manual check for ||a b|| = 0 (see Python code).
    end

    # Return the perturbed Fourier transform.
    return F̂
end

# Compute the evenness score of a given Fourier transform.
function evenness(F::Fourier)::defFloat
    B = F.B.^2
    AB = F.A.^2 .+ F.B.^2
    return sum( B[2:end] )/sum( AB[2:end] )
end

# Compute the best evenness score for the Fourier transform.
function evenness(F::Fourier, P::PowerSpectrum; δφ::defFloat=F.M.δx, φ̂::defFloat=F.X[end]
    )::Tuple{Vector{defFloat},Vector{defFloat}}
    # Initialize list of evenness scores.
    Φ = 0:δφ:φ̂;  N = length( Φ )
    Σ = Vector{defFloat}( undef, N )

    # Iterate through list, computing evenness score for each.
    for i ∈ 1:N
        F̂ = perturb( F, P, Φ[i] )
        Σ[i] = evenness( F̂ )
    end

    # Return the evenness list.
    return Φ, Σ
end

println( "Loaded symm.jl class file." )
