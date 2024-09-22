
from args import *

# Hyper parameter(s).
Nmax = 20
dN = 4
beta = 1e-3

# Square wave initialization.
def wave(x):
    return np.sign( x )

if __name__ == '__main__':
    # Generate x-data for square wave.
    T = 1;  Nt = round( T/beta ) + 1
    X = np.array( [[beta*(i-Nt+1) for i in range( 2*Nt )]] )
    Y = wave( X )

    # Initialize solver variable.
    fdft = RealFourier( X, Y ).dft().powerspec()
    # fptb = perturbseries( fdft, eps=5e-1 )

    # Plot results.
    fig, axs = plt.subplots()
    axs.plot( X.T, Y.T, color='r', label='Model' )

    Xf = np.linspace( -T, T, 1000 )[None]
    axs.plot( Xf.T, fdft.solve( Xf ).T, color='k', label='DFT' )
    # axs.plot( Xf.T, fptb.solve( Xf ).T, label='DFT perturbed' )

    fvar = RealFourier( X, Y )
    for n in range( 0, Nmax+1, dN ):
        if n == 0:  n = n + 1;
        fvar.dmd( N=n )
        Yf = fvar.solve( X )
        axs.plot( X.T, Yf.T, linestyle=None, label=('N=%i' % n) )

    axs.grid( 1 )
    axs.legend()

    fig2, axs2 = plt.subplots()
    axs2.plot( fdft.F, fdft.R.T )

    plt.show()
