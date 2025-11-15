
include( "args.jl" )


# Data set structure and helpers.
struct Data
    # Data-related variables.
    X::Vector{defFloat}
    Y::Vector{defFloat}
end

function Base.length(D::Data)
    return length( D.X )
end

function Plots.plot!(plt::Plots.Plot, D::Data; args...)::Plots.Plot
    return plot!( plt, D.X::AbstractVector, D.Y::AbstractVector; args... )
end

# Transform structure and helpers.
struct FourierModes
    # Data set.
    D::Data

    # Hyper-parameters.
    N::defInt       # Number of modes.
    δx::defFloat    # Difference between x-points.
    δω::defFloat    # "Step" of the frequency spectrum.

    # Frequency and power spectrum variables.
    ω::Vector{defFloat}
end

function FourierModes(X::Vector{defFloat}, Y::Vector{defFloat})::FourierModes
    # Truncate data set if given an odd number of points.
    ε = length( X ) % 2 == 0
    X̂ = ε ? X : X[1:end-1]
    Ŷ = ε ? Y : Y[1:end-1]
    D = Data( X̂, Ŷ )

    # Frequency spectrum hyper-parametrs.
    N = defInt( length( D )/2 )
    δx = defFloat( X̂[2] - X̂[1] )  # Assume evenly spaced data points.
    δω = defFloat( 2*N*δx )

    # Construct frequency spectrum.
    ω = defFloat.( 2π/δω.*(0:N) )

    return FourierModes( D, N, δx, δω, ω )
end

function Base.length(M::FourierModes)::defInt
    return M.N + 1
end

function Plots.plot!(plt::Plots.Plot, M::FourierModes; args...)::Plots.Plot
    return plot!( plt, M.D; args... )
end

# Real, discrete Fourier transform structure and helpers.
mutable struct Fourier
    # Fourier modes.
    M::FourierModes

    # Real transform coefficient list.
    A::Vector{defFloat}
    B::Vector{defFloat}

    # Power spectrum variable.
    P::Vector{defFloat}

    # Top-level data references.
    X::Vector{defFloat}
    Y::Vector{defFloat}
    ω::Vector{defFloat}
end

function Fourier(X::Vector{defFloat}, Y::Vector{defFloat})::Fourier
    # Construct the discrete Fourier modes.
    M = FourierModes( copy( X ), copy( Y ) )

    # Initialize Fourier amplitudes
    T = length( M )
    A = Vector{defFloat}( undef, T )
    B = Vector{defFloat}( undef, T )
    P = Vector{defFloat}( undef, T )

    return Fourier( M, A, B, P, M.D.X, M.D.Y, M.ω )
end

function Base.length(F::Fourier)::defInt
    return length( F.M )
end

function Plots.plot!(plt::Plots.Plot, F::Fourier; args...)::Plots.Plot
    return plot!( plt, F.M; args... )
end

function serialize(F::Fourier; X::Vector{defFloat}=F.X)::Tuple{Matrix{defFloat},Matrix{defFloat}}
    # Create serialized set of waves from frequency list.
    xSin = sin.( F.ω.*X' )
    xCos = cos.( F.ω.*X' )
    return xSin, xCos
end

function dft!(F::Fourier)::Fourier
    # Extract seriealized time-series.
    N = length( F ) - 1
    xSin, xCos = serialize( F )

    # Solve for 0-frequency.
    F.A[1] = 0.0
    F.B[1] = 1/(2*N)*sum( F.Y )  # Mean value of the data points.

    # Compute the Fourier coefficients.
    for i ∈ 2:length( F )-1
        F.A[i] = 1/N*sum( F.Y.*xSin[i,:] )
        F.B[i] = 1/N*sum( F.Y.*xCos[i,:] )
    end

    # Solve for the Nyquists frequency.
    F.A[end] = 0.0
    F.B[end] = 1/(2*N)*sum( F.Y.*xCos[end,:] )

    # Return the discrete Fourier transform.
    return F
end

function solve(F::Fourier; X::Vector{defFloat}=F.X)::Vector{defFloat}
    # Extract serialized wave lists.
    xSin, xCos = serialize( F; X )

    # Return Fourier series estimates.
    return vec( F.A'*xSin + F.B'*xCos )
end

function error(F::Fourier; X::Vector{defFloat}=F.X, Y::Vector{defFloat}=F.Y)::defFloat
    return √sum( (Y .- solve( F; X=X )).^2 )
end

println( "Loaded four.jl class file." )
