# pyinstaller --onefile mult.py
# File: abby.py
# Created by: Michael Napoli
# Created for: Abby Feldmann
# Purpose: To draw a line-based sketch using Fourier series of
#   a boy and a girl.

import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/four' )
sys.path.insert( 0, expanduser('~')+'/prog/geom' )

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from FOUR.Transforms import *
from GEOM.Vectors import *
from GEOM.Vehicle2D import *

datafile = 'sketchdata.csv'
# plt.style.use( 'dark_background' )

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
    fList = [ RealFourier( T, X ) for X, T in zip( Xlist, Tlist ) ]
    for fvar in fList:
        fvar.dmd( N=50 )

    # Initial conditions and plotting parameters.
    t = np.array( [[0]] )
    gList = 2.0*np.array( [3.5, 0.75, 2.5, 3.5] )
    xList = [ f.solve( t ) for f in fList ]
    zorderList = [ 200, 150, 50, 100 ]
    colorList = [ 'plum', 'yellowgreen', 'cornflowerblue', 'sandybrown' ]
    tailList = [ round( N/g ) - 5 for N, g in zip( Nlist, gList ) ]

    # Initialize plots.
    fig, axs = plt.subplots()
    axs.set_title( 'Abby and Michael' )
    plt.show( block=0 )

    # Vehicle variables.
    sim_glasses = 1
    vhcList = [ Vehicle2D( x, fig=fig, axs=axs,
        zorder=z, color=c, tail_length=l )
        for x, z, c, l in zip( xList, zorderList, colorList, tailList ) ]

    # Creating vector entities.
    offsetList = [ [600, 600, -175, -175], [-150, -150, -150, 500] ]
    vecList = [ f.vectors( t ) for f in fList ]
    vxList = [ Vectors( v[0]+[[0],[oy]], fig=fig, axs=axs, color=c )
        for v, c, oy in zip( vecList, colorList, offsetList[0] ) ]
    vyList = [ Vectors( np.flipud( v[1] )+[[ox],[0]], fig=fig, axs=axs, color=c )
        for v, c, ox in zip( vecList, colorList, offsetList[1] ) ]
    dxList = [ Vectors( np.hstack( (v[0][:,-1,None]+[[0],[oy]], x) ), fig=fig, axs=axs, color='grey' )
        for v, x, oy in zip( vecList, xList, offsetList[1] ) ]
    dyList = [ Vectors( np.hstack( (np.flipud( v[1][:,-1,None] )+[[ox],[0]], x) ), fig=fig, axs=axs, color='grey' )
        for v, x, ox in zip( vecList, xList, offsetList[0] ) ]

    iList = [ 0, 1, 2, 3 ]
    for i, vhc, vx, vy, dx, dy in \
        zip( iList, vhcList, vxList, vyList, dxList, dyList ):
        if (i != 1) or (i == 1 and sim_glasses):
            vhc.draw()
            vx.draw()
            vy.draw()
            dx.setLineStyle( ':' );  dx.setLineWidth( 1.25 )
            dy.setLineStyle( ':' );  dy.setLineWidth( 1.25 )
            dx.draw()
            dy.draw()

    # Axis edits and draw.
    # axs.axes.xaxis.set_ticklabels( [] )
    # axs.axes.yaxis.set_ticklabels( [] )
    axs.axis( 'equal' )
    axs.grid( 0 )

    # Simulate.
    dt = 1;  t = t + dt
    ans = input( "Press ENTER to start simulation..." )
    while t < 2500 and ans != 'n':
        for i, f, g, ox, oy in zip( iList, fList, gList, offsetList[0], offsetList[1] ):
            if (i != 1) or (i == 1 and sim_glasses):
                xList[i] = f.solve( g*t )
                vecList[i] = f.vectors( g*t )
                vxList[i].setVertices( vecList[i][0] )
                vxList[i].transform( dx=[[0],[oy]] )
                vyList[i].setVertices( np.flipud( vecList[i][1] ) )
                vyList[i].transform( dx=[[ox],[0]] )
                dxList[i].setVertices( np.hstack( (vecList[i][0,:,-1,None]+[[0],[oy]], xList[i]) ) )
                dyList[i].setVertices( np.hstack( (np.flipud( vecList[i][1,:,-1,None] )+[[ox],[0]], xList[i]) ) )
        for x, vhc, vx, vy, dx, dy in zip( xList, vhcList, vxList, vyList, dxList, dyList ):
            if (i != 1) or (i == 1 and sim_glasses):
                vhc.update( x )
                vx.update()
                vy.update()
                dx.update()
                dy.update()
        plt.pause( 1e-3 )
        t = t + dt
    if ans != 'n':
        input( "Press ENTER to end program..." )