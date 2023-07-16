import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from FOUR.Fouriers import *

datafile = 'data/xytest01.csv';

if __name__ == "__main__":
    # Import data set and create X/Y lists.
    csvcontents = pd.read_csv( datafile );
    points = csvcontents[ ['x','y'] ].to_numpy();

    # Extract and form data vectors.
    Nx = len( points  );
    T = np.array( [[i for i in range( Nx )]] );
    X = points[:,0].reshape( 1,Nx );
    Y = points[:,1].reshape( 1,Nx );
    XY = np.vstack( (X, Y) );