
from args import *

if __name__ == "__main__":
    # Generate random data set with equally spaced points.
    Nx = 20  # Chosen arbitrarily.
    L = 10   # Chosen arbitrarily.
    dx = np.random.rand()
    X = np.array( [[i*dx for i in range( Nx )]] )
    Y = 2*L*np.random.rand( 1,Nx ) - L

    # Results from class solver.
    fvar = RealFourier( X, Y )
    fvar.dft()
    print( fvar )
    print( fvar.A )
    print( fvar.B )

    # Solve over range.
    Xdft = np.linspace( 0, dx*(Nx - 1), 1000 )[None]
    Ydft = fvar.solve( Xdft )

    # Plot results.
    fig, axs = plt.subplots()
    plt.plot( X.T, Y.T, marker='o', label='true' )
    plt.plot( Xdft.T, Ydft.T, linestyle=':', label='solver' )
    plt.legend()
    plt.grid( 1 )
    plt.show()