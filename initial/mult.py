# pyinstaller --onefile mult.py
import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/geom' )

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from GEOM.Vectors import *
from GEOM.Vehicle2D import *

from args import *

datafile = 'data/xytest01.csv'

if __name__ == "__main__":
    # Import data set and create X/Y lists.
    csvcontents = pd.read_csv( datafile )
    points = csvcontents[ ['x','y'] ].to_numpy()

    # Extract and form data vectors.
    Nx = len( points )
    T = np.array( [ [i for i in range( Nx )] ] )
    X = points[:,0].reshape( 1,Nx )
    Y = points[:,1].reshape( 1,Nx )
    XY = np.vstack( (X, Y) )

    # Initialize Transform variables.
    fvar = RealFourier( T, XY ).dmd( N=25 )

    # Plot results.
    fig, axs = plt.subplots()
    axs.plot( X.T, Y.T, label='Drawing' )

    # Initialize vehicle.
    dt = 0.25
    t = np.array( [[0]] )
    x = fvar.solve( t )

    vhc = Vehicle2D( x, fig=fig, axs=axs,
        radius=10, tail_length=1000 ).draw()

    vxList, vyList = fvar.vectors( t )
    xconnector = np.hstack( (vxList[:,-1,None], x) )
    yconnector = np.hstack( (np.flipud( vyList[:,-1,None] ), x) )
    vxvar = Vectors( vxList, fig=fig, axs=axs ).draw()
    vyvar = Vectors( np.flipud( vyList ), fig=fig, axs=axs ).draw()
    cxvar = Vectors( xconnector, fig=fig, axs=axs, color='grey' ).draw()
    cyvar = Vectors( yconnector, fig=fig, axs=axs, color='grey' ).draw()
    axs.axis( 'equal' )

    while t < 500:
        t = t + dt
        x = fvar.solve( t )

        vxList, vyList = fvar.vectors( t )
        xconnector = np.hstack( (vxList[:,-1,None], x) )
        yconnector = np.hstack( (np.flipud( vyList[:,-1,None] ), x) )

        vxvar.setVertices( vxList )
        vyvar.setVertices( np.flipud( vyList ) )
        cxvar.setVertices( xconnector )
        cyvar.setVertices( yconnector )

        vxvar.update()
        vyvar.update()
        cxvar.update()
        cyvar.update()
        vhc.update( x )
        plt.pause( 1e-3 )
    input( "Press ENTER to end program..." )
