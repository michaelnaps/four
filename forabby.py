import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from FOUR.Fouriers import *

datafile = 'draw/xytest01.csv';

if __name__ == "__main__":
    # Import data set and create X/Y lists.
    csvcontents = pd.read_csv( datafile );
    points = csvcontents[ ['x','y'] ].to_numpy();

    # Extract and form data vectors.
    Nx = len( points );
    T = np.array( [ [i for i in range( Nx )] ] );
    X = points[:,0].reshape( 1,Nx );
    Y = points[:,1].reshape( 1,Nx );

    # Initialize Transform variables.
    fxvar = RealFourier( T, X );
    fyvar = RealFourier( T, Y );

    # Solve the DFT Problem.
    fxvar.dft();
    fyvar.dft();

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
    Xf = fxvar.solve( Tf );
    Yf = fyvar.solve( Tf );

    # Plot results.
    fig, axs = plt.subplots();
    axs.plot( X.T, Y.T, label='Drawing' );
    axs.plot( Xf.T, Yf.T, linestyle='--', label='Fourier' );
    plt.legend();
    plt.grid( 1 );
    plt.show();