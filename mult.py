# pyinstaller --onefile mult.py
import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/geom' )

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from FOUR.Transforms import *
from GEOM.Vectors import *
from GEOM.Vehicle2D import *

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
    fvar = RealFourier( T, XY )
    fvar.dft()

    # Plot results.
    fig, axs = plt.subplots()
    axs.plot( X.T, Y.T, label='Drawing' )

    # Initialize vehicle.
    dt = 0.1
    t = np.array( [[0]] )
    x = fvar.solve( t )

    vhc = Vehicle2D( x, fig=fig, axs=axs,
        radius=10, tail_length=1000 )
    # vvar = Vectors( fvar.vectors( t ), fig=fig, axs=axs )
    vhc.setLimits( xlim=(-50, 700), ylim=(-50, 500) )

    # vvar.draw()
    vhc.draw()
    while t < 500:
        t = t + dt
        x = fvar.solve( t )

        # vvar.setVertices( fvar.vectors( t ) )
        # vvar.update()
        vhc.update( x )
