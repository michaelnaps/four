import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from FOUR.Transforms import *

datafile = 'data/abbydrawing.csv';

if __name__ == "__main__":
    # Import data set and create X/Y lists.
    data = pd.read_csv( datafile );

    # Data sets fror each object.
    Xlist = (
        data[ data['z']=='a' ].to_numpy()[:,:2].T,
        data[ data['z']=='b' ].to_numpy()[:,:2].T,
        data[ data['z']=='c' ].to_numpy()[:,:2].T,
        data[ data['z']=='d' ].to_numpy()[:,:2].T
    );

    # Size of data sets.
    Nlist = [ X.shape[1] for X in Xlist ];
    print( Nlist );

    # Generate time series.
    Tlist = [ np.array( [[i for i in range( N )]] ) for N in Nlist ];

    # Fourier series lists.
    flist = [ RealFourier( X, T ) for X, T in zip( Xlist, Tlist ) ];
    for fvar in flist:
        fvar.dft();
        print( 1 );


    # # Plot data.
    # fig, axs = plt.subplots();
    # axs.plot( hat1[:,0], hat1[:,1], label='hat1' );
    # axs.plot( mike[:,0], mike[:,1], label='mike' );
    # axs.plot( abby[:,0], abby[:,1], label='abby' );
    # plt.legend();
    # plt.show();