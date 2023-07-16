import numpy as np
import matplotlib.pyplot as plt

from FOUR.Fouriers import *

if __name__ == "__main__":
    # Generate random data set with equally spaced points.
    Nx = 20;  # chosen arbitrarily.
    dx = np.random.rand();
    X = np.array( [[i*dx for i in range( Nx )]] );
    Y = 2*np.random.rand( 1,Nx ) - 1;
    Xdft = np.linspace(0,dx*(Nx-1),1000).reshape( 1,1000 );

    # Results from class solver.
    fvar = RealFourier( X, Y );
    fvar.dft();
    print( 'Error:', np.linalg.norm( Y - fvar.solve( X ) ) );
    print( 'A:', fvar.A );
    print( 'B:', fvar.B );
    Ysolver = fvar.solve( Xdft );

    # Plot results.
    fig, axs = plt.subplots();
    plt.plot( X.T, Y.T, marker='o', label='true' );
    plt.plot( Xdft.T, Ysolver.T, linestyle=':', label='solver' );
    plt.legend();
    plt.grid( 1 );
    plt.show();