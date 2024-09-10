# File: Preprocess.py
# Created by: Michael Napoli
# Created on: September 10, 2024
# Purpose: For methods of preprocessing before use with
#   Transform() class.

import numpy as np

def simplemovingaverage(T, X, W, dT=0.1):
    # Extract bounds of the new time-series.
    MT = T.shape[0];  MX = X.shape[0]
    Tmin = np.min( T );  Tmax = np.max( T )
    NT = round( (Tmax - Tmin)/dT ) + 1

    # Initialize new time-series.
    Tm = np.empty( (MT, NT) )
    Xm = np.empty( (MX, NT) )

    # Calculate SMA over the data set at selected points.
    for k in range( NT ):
        Tm[:,k] = k*dT + Tmin
        Xm[:,k] = np.mean( X[((Tm[:,k] - W) < T) & (T < (Tm[:,k] + W))] )

    # Return the equally spaced, averaged data.
    return Tm, Xm
