
from args import *

filename = 'model_fit_exp_ct.txt'
data = np.loadtxt( filename, dtype=float )

if __name__ == '__main__':
    # Perform Fourier transform.
    fvar = RealFourier( data[:,0][None], data[:,1][None] ).dft()

    # Generate fit data.
    xlist = np.linspace( np.min( data[:,0] ), np.max( data[:,0] ), 1000 )[None]
    ylist = fvar.solve( xlist )

    # Plot results.
    fig, axs1 = plt.subplots()
    axs1.plot( data[:,0], data[:,1] )
    axs1.plot( xlist.T, ylist.T, linestyle='--' )

    fig, axs2 = plt.subplots()
    axs2.plot( fvar.freqlist().F.T, fvar.powerspec().P.T, marker='x' )

    plt.show()
