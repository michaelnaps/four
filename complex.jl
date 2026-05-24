
include( "args.jl" )
include( "four.jl" )

mutable struct ComplexFourier
    # Fourier modes.
    M::FourierModes

    # Imaginary-valued coefficients.
    C::Vector{defComp}

    # Top-level data references.
    X::Vector{defFloat}
    Y::Vector{defFloat}
    ω::Vector{defFloat}
end

Base.broadcastable(F::ComplexFourier) = Ref( F )

function converter(A::Vector{defFloat}, B::Vector{defFloat})::Vector{defComp}
    # Zero-frequency component.
    C0 = complex( B[1], 0.0 )

    # Mid-frequency components.
    Ck = 0.5.*(B[2:end-1] .- im.*A[2:end-1])

    # Nyquist frequency component.
    CN = complex( B[end], 0.0 )

    # Return coefficients.
    return [C0;Ck;CN]
end

function ComplexFourier(F::Fourier)::ComplexFourier
    # Compute the complex Fourier coefficients.
    C = converter( F.A, F.B )

    # Return the new complex variable.
    return ComplexFourier( F.M, C, F.X, F.Y, F.ω )
end

function solve(F::ComplexFourier; X::Vector{defFloat}=F.X)::Vector{defComp}
    # Cross the points with the coefficients and their conjugates.
    ε = im.*(F.ω.*X')

    # Zero-frequency component.
    Y0 = F.C[1].*exp.( ε[1,:] )

    # Mid-frequency components.
    Yk = vec( sum( F.C[2:end-1].*exp.( ε[2:end-1,:] ); dims=1 ) )
    Ŷk = vec( sum( conj.( F.C[2:end-1] ).*exp.( -ε[2:end-1,:] ); dims=1 ) )

    # Nyquist frequency component.
    YN = F.C[end].*exp.( ε[end,:] )
    return Y0 .+ Yk .+ Ŷk .+ YN
end

# ---- HELPER FOR CONVERTING BACK TO REAL DOMAIN ---- #
function converter(C::Vector{defComp})::Tuple{Vector{defFloat},Vector{defFloat}}
    # Zero-frequency components.
    A0 = 0.0
    B0 = real( C[1] )

    # Mid-frequency components.
    Ak = -2.0.*imag.( C[2:end-1] )
    Bk =  2.0.*real.( C[2:end-1] )

    # Nyquist frequency.
    AN = 0.0
    BN = real.( C[end] )

    # Return coefficients.
    return ([A0;Ak;AN], [B0;Bk;BN])
end

function Fourier(F::ComplexFourier)::Fourier
    # Compute the complex Fourier coefficients.
    A, B = converter( F.C )

    # Return the new complex variable.
    return Fourier( F.M, A, B, F.X, F.Y, F.ω )
end

println( "Loaded complex.jl class file." )
