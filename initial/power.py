
from args import *

# Hyper parameter(s).
Nmax = 100
dN = 10
beta = 0.001

# Square wave initialization.
def model(x):
    return x**3

if __name__ == '__main__':
    # Generate x-data for square wave.
    T = 2.5;  Nt = round( T/beta ) + 1
    X = np.array( [[beta*(i-Nt+1) for i in range( 2*Nt )]] )
    Y = model( X )

    # Initialize solver variable.
    fvar = RealFourier( X, Y )

    # Plot results.
    fig, axs = plt.subplots()
    axs.plot( X.T, Y.T, color='r', label='Model' )

    Xf = np.linspace( -T-1, T+1, 2*Nt )[None]
    for n in range( 1, Nmax+1, dN ):
        N = n-1 if n > 1 else n
        fvar.dmd( N=N )
        Yf = fvar.solve( Xf )
        axs.plot( Xf.T, Yf.T, linestyle=None, label=('N=%i' % N) )
        print( fvar )

    plt.grid( 1 )
    plt.legend()
    plt.show()
