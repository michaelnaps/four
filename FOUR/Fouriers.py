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

    def setDataSets(self, X, Y):
        self.X = X;
        self.Y = Y;
        # Return instance of self.
        return self;

    def setCoefficientNumber(self, N):
        self.N = N;
        # Return instance of self.
        return self;


class RealFourier( Fourier ):
    def __init__(self, X=None, Y=None, N=-1, h=1e-3):
        self.A = None;
        self.B = None;
        Fourier.__init__( self, X=X, Y=Y, N=N, h=h );

    def generateSeries(self, X):
        assert self.N != -1, "ERROR: Coefficient number, N, is not set."
        M = X.shape[1];
        xSinList = np.ones( (M, self.N) );
        xCosList = np.empty( (M, self.N) );
        for k in range( self.N ):
            xSinList[:,k] = np.sin( 2*np.pi*k*X/self.h );
            xCosList[:,k] = np.cos( 2*np.pi*k*X/self.h );
        # Return sin(X) and cox(X) lists separately.
        return xSinList, xCosList;

    def propagate(self, X):
        assert self.A is not None, "ERROR: A coefficient vector is empty."
        assert self.B is not None, "ERROR: B coefficient vector is empty."
        xSinList, xCosList = self.generateSeries( X );
        print( xSinList.shape, xCosList.shape )
        # Return approximate solution using coefficient matrices.
        return self.A@xSinList.T + self.B@xCosList.T;

    def dft(self, X=None, Y=None, h=None):
        if X is not None:
            self.setDataSets( X, Y );
        if h is not None:
            self.h = h;

        assert self.X is not None, "ERROR: X data set is empty."
        assert self.Y is not None, "ERROR: Y data set is empty."

        # Set expansion number (strict in DFT).
        self.N = int( self.X.shape[1]/2 );

        self.A = np.zeros( (1, self.N) );
        self.B = np.zeros( (1, self.N) );

        self.A[0][0] = 0;
        self.B[0][0] = 1/(2*self.N)*np.sum( self.Y );

        xSinList, xCosList = self.generateSeries( X );
        for k in range( 1, self.N ):
            for y, xSin, xCos in zip( self.Y[0], xSinList[:,k], xCosList[:,k] ):
                self.A[0][k] += 1/self.N*y*xSin;
                self.B[0][k] += 1/self.N*y*xCos;

        self.A[0][-1] = 0;
        self.B[0][-1] = 1/(2*self.N)*np.sum( Y[:,-1]*xCosList[:,-1] );

        # Return instance of self.
        return self;


class ComplexFourier( Fourier ):
    def __init__(self, X=None, Y=None, N=-1, h=1e-3):
        self.C = None;
        Fourier.__init__( self, X=X, Y=Y, N=N, h=h );