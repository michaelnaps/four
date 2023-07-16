import numpy as np
import matplotlib.pyplot as plt

from FOUR.Fouriers import *

if __name__ == "__main__":
    X = np.array( [[i for i in range( 10 )]] );
    Y = np.random.rand( 1, 10 );
    Xdft = np.linspace(0,9,1000).reshape( 1,1000 );

    # # Results from handwritten notes.
    # f = lambda x: 1/2 + 1/2*np.sin( np.pi*x/2 ) - 1/2*np.cos( np.pi*x/2 );
    # Ydft = f( Xdft );

    # Results from class solver.
    fvar = RealFourier( X, Y );
    fvar.dft();
    print( 'A:', fvar.A );
    print( 'B:', fvar.B );
    Ysolver = fvar.solve( Xdft );

    fig, axs = plt.subplots()
    plt.plot( X.T, Y.T, marker='o', label='true' );
    # plt.plot( Xdft.T, Ydft.T, linestyle='--', label='dft' );
    plt.plot( Xdft.T, Ysolver.T, linestyle=':', label='solver' );

    plt.legend();
    plt.grid( 1 );
    plt.show();