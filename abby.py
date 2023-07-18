# pyinstaller --onefile mult.py
import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/mpc' );

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from MPC.Vehicle2D import *
from FOUR.Transforms import *

datafile = 'data/abbydrawing03.csv';

if __name__ == "__main__":
    # Import data set and create X/Y lists.
    data = pd.read_csv( datafile );

    # Data sets fror each object.
    Xlist = [
        data[ data['z']=='a' ].to_numpy()[:,:2].T,
        data[ data['z']=='b' ].to_numpy()[:,:2].T,
        data[ data['z']=='c' ].to_numpy()[:,:2].T,
        data[ data['z']=='d' ].to_numpy()[:,:2].T
    ];

    # Size of data sets - make sets even.
    Nlist = [ X.shape[1] for X in Xlist ];
    for i, N in enumerate( Nlist ):
        if N % 2 != 0:
            Nlist[i] = N - 1;
            Xlist[i] = Xlist[i][:,:-1];

    # Generate time series.
    Tlist = [ np.array( [[i for i in range( N )]] ) for N in Nlist ];

    # Fourier series lists.
    flist = [ RealFourier( T, X ) for X, T in zip( Xlist, Tlist ) ];
    for fvar in flist:
        # print( fvar.N );
        # fvar.dft();
        fvar.ls( N=100 );

    # Initial conditions.
    t = np.array( [[0]] );
    xa = flist[0].solve( t );
    xm = flist[2].solve( t );
    xh = flist[3].solve( t );

    # Create vehicles.
    fig, axs = plt.subplots();
    abby = Vehicle2D( None, xa, fig=fig, axs=axs, vhc_color='k', tail_length=2500 );
    mike = Vehicle2D( None, xm, fig=fig, axs=axs, vhc_color='k', tail_length=2500 );
    hat1 = Vehicle2D( None, xh, fig=fig, axs=axs, vhc_color='k', tail_length=2500 );
    abby.setLimits( xlim=(0,500), ylim=(0,400) );

    axs.grid( 0 );
    abby.draw();

    # Simulate.
    dt = 1;  t = t + dt;
    ans = input( "Press ENTER to start simulation..." );
    while t < 10000 and ans != 'n':
        xa = flist[0].solve( 4.0*t );
        xm = flist[2].solve( 1.5*t );
        xh = flist[3].solve( 1.0*t );

        abby.update( xa, pause=0 );
        mike.update( xm, pause=0 );
        hat1.update( xh );

        t = t + dt;
    if ans != 'n':
        input( "Press ENTER to end program..." );