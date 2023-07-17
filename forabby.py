import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from FOUR.Fouriers import *

datafile = 'data/abbydrawing.csv';

if __name__ == "__main__":
    # Import data set and create X/Y lists.
    csvcontents = pd.read_csv( datafile );
    points = csvcontents[ ['x','y'] ].to_numpy();

    # Extract and form data vectors.
    Nx = len( points );
    T = np.array( [ [i for i in range( Nx )] ] );
    X = points[:,0].reshape( 1,Nx );
    Y = points[:,1].reshape( 1,Nx );
    XY = np.vstack( (X, Y) );

    # Initialize Transform variables.
    fvar = RealFourier( T, XY );
    fvar.dft();

    # print( 'X:' );
    # print( '\tA:', fxvar.A );
    # print( '\tB:', fxvar.B );
    # print( 'Y:' );
    # print( '\tA:', fyvar.A );
    # print( '\tB:', fyvar.B );

    # Compute comparison vectors.
    dt = 0.01;
    Nf = round( Nx/dt );
    Tf = np.array( [[ i*dt for i in range( Nf ) ]] );
    XYf = fvar.solve( Tf );

    # Plot results.
    fig, axs = plt.subplots();
    axs.plot( X.T, Y.T, label='Drawing' );
    axs.plot( XYf[0], XYf[1], linestyle='--', label='Fourier' );
    plt.legend();
    plt.grid( 1 );
    plt.show();