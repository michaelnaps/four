# File: LearningStrategies.py
# Created by: Michael Napoli
# Created on: Jul 11, 2023
# Purpose: To develop and study various methods of expanding
#   and solving real/complex Fourier series.

import numpy as np

# Class: Transform()
class Transform:
    def __init__(self, T, X, N=None, dt=None):
        self.T = T;
        self.X = X;

        if N is None:
            self.N = round( T.shape[1]/2 );
        else:
            self.N = N;

        if dt is None:
            self.dt = self.T[0,1] - self.T[0,0];
        else:
            dt = dt;

        self.tau = 2*self.N*self.dt;

    def setDataSets(self, T, X):
        self.__init__( T, X );
        # Return instance of self.
        return self;

    def setCoefficientNumber(self, N):
        self.N = N;
        # Return instance of self.
        return self;


# Class: RealFourier()
class RealFourier( Transform ):
    def __init__(self, T, X, N=None, dt=None):
        # Intialize coefficient matrices to None.
        self.A = None;
        self.B = None;
        self.Nx = X.shape[0];

        # Initialize Transform() as parent class.
        Transform.__init__( self, T, X, N=N, dt=dt );

    def serialize(self, T=None):
        # If data set is given use instead of 'default'.
        if T is None:
            T = self.T
            M = 2*self.N;
        else:
            M = T.shape[1];

        # Initialize serialized set matrices.
        tSin = np.zeros( (self.N+1, M) );
        tCos = np.ones(  (self.N+1, M) );

        # Iterate through for varying frequencies.
        for k in range( 1, self.N+1 ):
            tSin[k,:] = np.sin( 2*np.pi*k*T[0,:]/self.tau );
            tCos[k,:] = np.cos( 2*np.pi*k*T[0,:]/self.tau );

        # Return the serialized sets.
        return tSin, tCos;

    def dft(self):
        # Serialize the given data set.
        tSin, tCos = self.serialize();

        # Initialize coefficient vectors.
        self.A = np.empty( (self.Nx, self.N+1) );
        self.B = np.empty( (self.Nx, self.N+1) );

        for i in range( self.Nx ):
            # Solve for when k=0.
            self.A[i,0] = 0;
            self.B[i,0] = 1/(2*self.N)*sum( self.X[i,:] );

            # Solve for when 0 < k < N.
            for k in range( 1,self.N ):
                self.A[i,k] = 1/self.N*sum( self.X[i,:]*tSin[k,:] );
                self.B[i,k] = 1/self.N*sum( self.X[i,:]*tCos[k,:] );

            # Solve for when k = N.
            self.A[i,self.N] = 0;
            self.B[i,self.N] = 1/(2*self.N)*sum( self.X[i,:]*tCos[self.N,:] );

        # Return instance of self.
        return self;

    def solve(self, T=None):
        # Is given set is none, use default.
        if T is None:
            T = self.T;
            Nt = 2*self.N;
        else:
            Nt = T.shape[1];

        # Get serialized form of data set.
        tSin, tCos = self.serialize( T );

        # Return approtimation using coefficient matrices.
        Y = np.empty( (self.Nx, Nt) );
        for i in range( self.Nx ):
            Y[i,:] = self.A[i,:]@tSin + self.B[i,:]@tCos;
        return Y;


# Class: ComplexFourier()
class ComplexFourier( Transform ):
    def __init__(self, T, X, N=None, dt=None):
        # Intialize coefficient matrices to None.
        self.C = None;

        # Initialize Transform() as parent class.
        Transform.__init__( self, T, X, N=N, dt=dt );