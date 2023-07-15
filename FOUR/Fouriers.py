# File: LearningStrategies.py
# Created by: Michael Napoli
# Created on: Jul 11, 2023
# Purpose: To develop and study various methods of expanding
#   and solving real/complex Fourier series.

import numpy as np

class Fourier:
    def __init__(self, X=None, Y=None, N=-1, h=1e-3):
        self.X = X;
        self.Y = Y;
        self.h = h;
        self.N = N;
        self.A = None;
        self.B = None;

    def setDataSets(self, X, Y):
        self.X = X;
        self.Y = Y;
        # Return instance of self.
        return self;

    def setCoefficientNumber( self, N ):
        self.N = N;
        # Return instance of self.
        return self;


class RealFourier( Fourier ):
    def __init__(self, X=None, Y=None, N=-1, h=1e-3):
        Fourier.__init__( self, X=X, Y=Y, N=N, h=h );

    def generateSeries(self, X):
        assert self.N != -1, "ERROR: Coefficient number, N, is not set."
        M = X.shape[1];
        xSinList = np.ones( (M, 2*self.N) );
        xCosList = np.empty( (M, 2*self.N) );
        for k in range( 2*self.N ):
            xSinList[:,k] = np.sin( 2*np.pi*k*X/self.h );
            xCosList[:,k] = np.cos( 2*np.pi*k*X/self.h );
        return xSinList, xCosList;

    def propagate(self, X):
        assert self.A is not None, "ERROR: A coefficient vector is empty."
        assert self.B is not None, "ERROR: B coefficient vector is empty."
        xSinList, xCosList = self.generateSeries( X );
        print( xSinList.shape, xCosList.shape )
        Y = self.A@xSinList + selfB@xCostList;

    def dft(self, X=None, Y=None, h=None):
        if X is not None:
            self.setDataSets( X, Y );
        if h is not None:
            self.h = h;

        # Set expansion number (strict in DFT).
        self.N = int( self.X.shape[1]/2 );

        self.A = np.empty( (1, self.N) );
        self.B = np.empty( (1, self.N) );

        self.A[0][0] = 0;
        self.B[0][0] = 1/(2*self.N)*np.sum( self.Y );

        xSinList, xCosList = self.generateSeries( X );
        for k in range( 1, self.N ):
            for y, xSin, xCos in zip( self.Y[0], xSinList[0], xCosList[0] ):
                self.A[0][k] += y*xSin;
                self.B[0][k] += y*xCos;

        self.A[0][-1] = 0;
        self.B[0][-1] = 0; # 1/(2*self.N)*np.sum( [ self.Y[0][j]*np.cos( 2*np.pi*self.N*self.X[0][j]/self.h ) for j in range( self.N ) ] );

        # Return instance of self.
        return self;
