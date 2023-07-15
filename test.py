import sys
from os.path import expanduser
sys.path.insert(0, expanduser('~')+'/prog/fourier');

import numpy as np
import matplotlib.pyplot as plt

# Personal classes.
from FOUR.Fouriers import *

# Hyper parameter(s).
beta = 0.01;

# Square wave initialization.
def wave(x):
    return np.sign( x );

if __name__ == '__main__':
    # Generate x-data for square wave.
    T = 0.03;  Nt = round( T/beta ) + 1;
    X = np.array( [[beta*(i-Nt+1) for i in range( 2*Nt-1 )]] );
    Y = wave( X );

    print( 'X:\n', X );
    print( 'Y:\n', Y );

    # Test Fourier method class.
    fvar = RealFourier( h=beta );
    fvar.dft( X, Y );

    print( fvar.generateSeries( X )[1] );

    # # Results data
    # print( fvar.A.shape );
    # print( fvar.B.shape );
    # xSinList, xCosList = fvar.liftData( X );
    # Yf = fvar.A@xSinList+ fvar.B@xCosList;

    # print( X );
    # print( Y );
    # print( Yf );

    print ( 'A:\n', fvar.A );
    print ( 'B:\n', fvar.B );

    Yf = fvar.propagate( X );

    # Plot results
    fig, axs = plt.subplots();
    axs.plot( X.T, Y.T );
    axs.plot( X.T, Yf.T );
    axs.grid( 1 );
    plt.show();