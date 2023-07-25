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
    fList = [ RealFourier( T, X ) for X, T in zip( Xlist, Tlist ) ]
    for fvar in fList:
        # print( fvar.N )
        # fvar.dft()
        fvar.ls( N=75 )

    # Initial conditions and plotting parameters.
    t = np.array( [[0]] )
    gList = 1.5*np.array( [3.5, 0.75, 2.0, 3.0] )
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
    vhcList[0].setLimits( xlim=(-50,800), ylim=(-300,400) )

    # Creating vector entities.
    offsetList = [ [600, 0, 0, 0], [-150, 0, 0, 0] ]
    vecList = [ f.vectors( t ) for f in fList ]
    vxList = [ Vectors( v[0], fig=fig, axs=axs, color=c )
        for v, c in zip( vecList, colorList ) ]
    vyList = [ Vectors( np.flipud( v[1] ), fig=fig, axs=axs, color=c )
        for v, c in zip( vecList, colorList ) ]
    # cabby = [
    #     Vectors( np.hstack( (avectors[0,:,-1,None] , xa) ), fig=fig, axs=axs, color='grey' ),
    #     Vectors( np.hstack( (np.flipud( avectors[1,:,-1,None] ) , xa) ), fig=fig, axs=axs, color='grey' )
    # ]

    # Axis edits and draw.
    fig.tight_layout()
    # axs.axes.xaxis.set_ticklabels( [] )
    # axs.axes.yaxis.set_ticklabels( [] )
    axs.grid( 0 )

    for vhc, vx, vy in zip( vhcList, vxList, vyList ):
        vhc.draw()
        vx.draw( new=0 )
        vy.draw( new=0 )

    # Simulate.
    iList = [ 0, 1, 2, 3 ]
    dt = 1;  t = t + dt
    ans = input( "Press ENTER to start simulation..." )
    while t < 2500 and ans != 'n':
        for i, f, g in zip( iList, fList, gList ):
            xList[i] = f.solve( g*t )
            vecList[i] = f.vectors( g*t )
            vxList[i].setVertices( vecList[i][0] )
            vyList[i].setVertices( np.flipud( vecList[i][1] ) )


        avectors = fList[0].vectors( gList[0]*t )
        for i, v, c, vList in zip( [0,1], vabby, cabby, avectors ):
            if i == 0:
                vList[1] = vList[1] - 150
            if i == 1:
                vList = np.flipud( vList )
                vList[0] = vList[0] + 600
            v.setVertices( vList )
            c.setVertices( np.hstack( (vList[:,-1,None], xa) ) )
            v.update()
            c.update()

        abby.update( xa, pause=0 )
        mike.update( xm, pause=0 )
        hat1.update( xh )

        t = t + dt
    if ans != 'n':
        input( "Press ENTER to end program..." )