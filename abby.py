# pyinstaller --onefile mult.py
# File: abby.py
# Created by: Michael Napoli
# Created for: Abby Feldmann
# Purpose: To draw a line-based sketch using Fourier series of
#   a boy and a girl.

import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/mpc' )

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from MPC.Vehicle2D import *
from FOUR.Transforms import *

datafile = 'data/abbydrawing03.csv'

if __name__ == "__main__":
    # Import data set and create X/Y lists.
    data = pd.read_csv( datafile )

    # Data sets fror each object.
    Xlist = [
        data[ data['z']=='a' ].to_numpy()[:,:2].T,
        data[ data['z']=='b' ].to_numpy()[:,:2].T,
        data[ data['z']=='c' ].to_numpy()[:,:2].T,
        data[ data['z']=='d' ].to_numpy()[:,:2].T
    ]

    # Size of data sets - make sets even.
    Nlist = [ X.shape[1] for X in Xlist ]
    for i, N in enumerate( Nlist ):
        if N % 2 != 0:
            Nlist[i] = N - 1
            Xlist[i] = Xlist[i][:,:-1]

    # Generate time series.
    Tlist = [ np.array( [[i for i in range( N )]] ) for N in Nlist ]

    # Fourier series lists.
    flist = [ RealFourier( T, X ) for X, T in zip( Xlist, Tlist ) ]
    for fvar in flist:
        # print( fvar.N )
        # fvar.dft()
        fvar.ls( N=75 )

    # Initial conditions.
    glist = [ 4.0, 0.0, 1.5, 1.0 ]
    t = np.array( [[0]] )
    xa = flist[0].solve( t )
    xm = flist[2].solve( t )
    xh = flist[3].solve( t )

    # Create vehicles.
    fig, axs = plt.subplots()
    axs.set_title( 'Abby and Michael' )

    abby = Vehicle2D( xa, fig=fig, axs=axs,
        vhc_color='plum', tail_length=round( Nlist[0]/glist[0] )-5 )
    mike = Vehicle2D( xm, fig=fig, axs=axs,
        vhc_color='cornflowerblue', tail_length=round( Nlist[2]/glist[2] )-5 )
    hat1 = Vehicle2D( xh, fig=fig, axs=axs,
        vhc_color='sandybrown', tail_length=round( Nlist[3]/glist[3] )-5 )
    abby.setFigureDimensions( w=5, h=4 )
    abby.setLimits( xlim=(0,500), ylim=(0,400) )

    fig.tight_layout()
    axs.axes.xaxis.set_ticklabels( [] )
    axs.axes.yaxis.set_ticklabels( [] )
    axs.grid( 0 )
    abby.draw()

    # Simulate.
    dt = 1;  t = t + dt
    ans = input( "Press ENTER to start simulation..." )
    while t < 10000 and ans != 'n':
        xa = flist[0].solve( glist[0]*t )
        xm = flist[2].solve( glist[2]*t )
        xh = flist[3].solve( glist[3]*t )

        abby.update( xa, pause=0 )
        mike.update( xm, pause=0 )
        hat1.update( xh )

        t = t + dt
    if ans != 'n':
        input( "Press ENTER to end program..." )