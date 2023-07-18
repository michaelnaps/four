# pyinstaller --onefile mult.py
import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/mpc' );

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from MPC.Vehicle2D import *
from FOUR.Transforms import *

datafile = 'data/abbydrawing01.csv';

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
        fvar.ls( N=150 );

    # Initial conditions.
    t = np.array( [[0]] );
    xh = flist[0].solve( t );
    xm = flist[1].solve( t );
    xa = flist[3].solve( t );

    # Create vehicles.
    fig, axs = plt.subplots();
    hat1 = Vehicle2D( None, xh, fig=fig, axs=axs, vhc_color='k', tail_length=2500 );
    mike = Vehicle2D( None, xm, fig=fig, axs=axs, vhc_color='k', tail_length=2500 );
    abby = Vehicle2D( None, xa, fig=fig, axs=axs, vhc_color='k', tail_length=2500 );
    abby.setLimits( xlim=(100,550), ylim=(40,450) );

    axs.grid( 0 );
    abby.draw();

    # Simulate.
    dt = 0.5;
    t = t + dt;
    ans = input( "Press ENTER to start simulation..." );
    while t < 10000 and ans != 'n':
        xh = flist[0].solve( t );
        xm = flist[1].solve( t );
        xa = flist[3].solve( t );

        hat1.update( xh, pause=0 );
        mike.update( xm, pause=0 );
        abby.update( xa );

        t = t + dt;
    if ans != 'n':
        input( "Press ENTER to end program..." );