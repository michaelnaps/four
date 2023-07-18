import sys
from os.path import expanduser
sys.path.insert(0, expanduser('~')+'/prog/kman')

import numpy as np
import matplotlib.pyplot as plt

# Personal classes.
from FOUR.Transforms import *

# Hyper parameter(s).
Nmax = 20
dN = 4
beta = 0.01

# Square wave initialization.
def wave(x):
    return np.sign( x )

if __name__ == '__main__':
    # Generate x-data for square wave.
    T = 1;  Nt = round( T/beta ) + 1
    X = np.array( [[beta*(i-Nt+1) for i in range( 2*Nt )]] )
    Y = wave( X )

    # Initialize solver variable.
    fvar = RealFourier( X, Y )

    # Plot results.
    fig, axs = plt.subplots()
    axs.plot( X.T, Y.T, color='r', label='Model' )

    for n in range( 0, Nmax+1, dN ):
        if n == 0:  n = n + 1;
        fvar.ls( N=n )
        Yf = fvar.solve( X )
        axs.plot( X.T, Yf.T, linestyle=None, label=('N=%i' % n) )

    plt.grid( 1 )
    plt.legend()
    plt.show()
