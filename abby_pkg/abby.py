# pyinstaller --onefile mult.py
# File: abby.py
# Created by: Michael Napoli
# Created for: Abby Feldmann
# Purpose: To draw a line-based sketch using Fourier series of
#   a boy and a girl.

import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/four' )
sys.path.insert( 0, expanduser('~')+'/prog/mpc' )

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from MPC.Vehicle2D import *
from FOUR.Transforms import *

datafile = 'sketchdata.csv'
plt.style.use( 'dark_background' )

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
    sim_glasses = 0
    glist = 1.0*np.array( [3.5, 0.75, 2.0, 3.0] )
    t = np.array( [[0]] )
    xa = flist[0].solve( t )
    if sim_glasses:
        xg = flist[1].solve( t )
    xm = flist[2].solve( t )
    xh = flist[3].solve( t )

    # Initialize plots.
    fig, axs = plt.subplots()
    axs.set_title( 'Abby and Michael' )

    # Vehicle variables.
    abby = Vehicle2D( xa, fig=fig, axs=axs, zorder=200,
        vhc_color='plum', tail_length=round( Nlist[0]/glist[0] )-5 )
    if sim_glasses:
        glss = Vehicle2D( xg, fig=fig, axs=axs, zorder=150,
        vhc_color='yellowgreen', tail_length=round( Nlist[1]/glist[1] )-5 )
    mike = Vehicle2D( xm, fig=fig, axs=axs, zorder=50,
        vhc_color='cornflowerblue', tail_length=round( Nlist[2]/glist[2] )-5 )
    hat1 = Vehicle2D( xh, fig=fig, axs=axs, zorder=100,
        vhc_color='sandybrown', tail_length=round( Nlist[3]/glist[3] )-5 )
    abby.setFigureDimensions( w=4.75, h=4.40 )
    abby.setLimits( xlim=(-10,425), ylim=(-20,375) )

    # Axis edits and draw.
    fig.tight_layout()
    axs.axes.xaxis.set_ticklabels( [] )
    axs.axes.yaxis.set_ticklabels( [] )
    axs.grid( 0 )
    abby.draw()

    # Simulate.
    dt = 1;  t = t + dt
    ans = input( "Press ENTER to start simulation..." )
    while t < 2500 and ans != 'n':
        xa = flist[0].solve( glist[0]*t )
        if sim_glasses and t > 250:
            xg = flist[1].solve( glist[1]*(t - 250) )
            glss.update( xg, pause=0 )
        xm = flist[2].solve( glist[2]*t )
        xh = flist[3].solve( glist[3]*t )

        abby.update( xa, pause=0 )
        mike.update( xm, pause=0 )
        hat1.update( xh )

        t = t + dt
    if ans != 'n':
        input( "Press ENTER to end program..." )