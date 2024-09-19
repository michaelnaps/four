
from args import *

if __name__ == "__main__":
    # Generate random data set with equally spaced points.
    Nx = 20   # Chosen arbitrarily.
    L = 10    # Chosen arbitrarily.
    dx = np.random.rand()
    Xlearn = np.array( [[i*dx for i in range( Nx )]] )
    Ylearn = 2*L*np.random.rand( 1,Nx ) - L

    # Results from class solver.
    fvar = RealFourier( Xlearn, Ylearn ).dft()
    print( fvar )
    print( realcentroid( fvar ) )

    # Solve over range using real Fourier.
    X = np.linspace( 0, dx*(Nx - 1), 1000 )[None]
    Yr = fvar.solve( X )

    # Test real -> complex function.
    cvar = ComplexFourier( Xlearn, Ylearn ).dft().powerspec()
    print( cvar )
    print( complexcentroid( cvar ) )

    # Solve over range using complex Fourier.
    Yc = cvar.solve( X )

    # Plot results.
    fig1, axs1 = plt.subplots()
    axs1.plot( Xlearn.T, Ylearn.T, marker='o', label='true' )
    axs1.plot( X.T, Yr.T, linestyle='--', label='real' )
    axs1.plot( X.T, np.real( Yc.T ), linestyle=':', label='complex' )
    axs1.legend()
    axs1.grid( 1 )

    fig2, axs2 = plt.subplots()
    axs2.plot( fvar.w, fvar.R.T, linestyle='--', label='real' )
    axs2.plot( cvar.w[cvar.N-1:], cvar.R.T[cvar.N-1:], linestyle=':', label='complex' )
    axs2.legend()

    plt.show()