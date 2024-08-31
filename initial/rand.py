
from args import *

if __name__ == "__main__":
    # Generate random data set with equally spaced points.
    Nx = 20  # Chosen arbitrarily.
    L = 10   # Chosen arbitrarily.
    dx = np.random.rand()
    Xlearn = np.array( [[i*dx for i in range( Nx )]] )
    Ylearn = 2*L*np.random.rand( 1,Nx ) - L

    # Results from class solver.
    fvar = RealFourier( Xlearn, Ylearn )
    fvar.dft()
    print( fvar )

    # Solve over range using real Fourier.
    X = np.linspace( 0, dx*(Nx - 1), 1000 )[None]
    Yr = fvar.solve( X )

    # Test real -> complex function.
    cvar = ComplexFourier( Xlearn, Ylearn ).dft( verbose=1 )
    print( cvar )

    # Solve over range using complex Fourier.
    Yc = cvar.solve( X )

    # Plot results.
    fig, axs = plt.subplots()
    plt.plot( Xlearn.T, Ylearn.T, marker='o', label='true' )
    plt.plot( X.T, Yr.T, linestyle='--', label='real' )
    plt.plot( X.T, np.real( Yc.T ), linestyle=':', label='complex' )
    plt.legend()
    plt.grid( 1 )
    plt.show()
