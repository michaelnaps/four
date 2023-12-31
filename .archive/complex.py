import sys
from os.path import expanduser
sys.path.insert(0, expanduser('~')+'/prog/kman')

import numpy as np
import matplotlib.pyplot as plt

# Personal classes.
from KMAN.Operators import *

# Hyper parameter(s).
Nmax = 25
dN = 5
beta = 0.01

# Square wave initialization.
def wave(x):
    return x**2

# Observables
def thetaN(X, N=1):
    Theta = np.empty( (N+1, X.shape[1]) )
    for i in range( -N, N+1 ):
        Theta[i,:] = np.exp( i*X )  # not complex
    print( Theta )
    return Theta

if __name__ == '__main__':
    # Generate x-data for square wave.
    T = 1.0;  Nt = round( T/beta ) + 1
    X = np.array( [[beta*(i-Nt+1) for i in range( 2*Nt-1 )]] )
    Y = wave( X )

    # Plot results.
    fig, axs = plt.subplots()
    axs.plot( X.T, Y.T, color='r', label='Model' )

    for n in range( 0, Nmax+1, dN ):
        theta = lambda x=None: thetaN( x, N=n )

        solver = Regressor( theta( X ), Y )
        C, _ = solver.dmd()

        print( C )
        print( '---------' )

        Yf = C@theta( X )
        axs.plot( X.T, Yf.T, linestyle='-.', label=('N=%i' % n) )

    plt.grid( 1 )
    plt.legend()
    plt.show()
