# File: LearningStrategies.py
# Created by: Michael Napoli
# Created on: Jul 11, 2023
# Purpose: To develop and study various methods of expanding
#   and solving real/complex Fourier series.

import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/kman' )

import numpy as np
from KMAN.Regressors import *

# Class: Transform()
class Transform:
    def __init__(self, T, X, N=None, dt=None):
        # Initialize transform variables and dimensions.
        self.T = T      # Input data.
        self.X = X      # Output data.
        self.F = None   # List of frequencies.
        self.K = T.shape[0]
        self.N = round( T.shape[1]/2 ) if N is None else N
        self.dt = self.T[0,1] - self.T[0,0] if dt is None else dt
        self.tau = 2*self.N*self.dt
        self.err = None

    def setDataSets(self, T, X):
        self.__init__( T, X )
        # Return instance of self.
        return self

    def setCoefficientNumber(self, N):
        self.N = N
        # Return instance of self.
        return self

    def frequencies(self):
        self.F = np.empty( (self.K*(self.N+1),1) )
        for k in range( self.K*(self.N+1) ):
            self.F[k,:] = 2*np.pi*k/self.tau

        # Return instance of self.
        return self

# Class: RealFourier()
class RealFourier( Transform ):
    def __init__(self, T, X, N=None, dt=None):
        # Intialize coefficient matrices to None.
        self.A = None   # sin(t) coefficients.
        self.B = None   # cos(t) coefficients.
        self.P = None   # Individual power spectrum values.
        self.R = None   # Power spectrum.
        self.Nt = X.shape[0]

        # Initialize Transform() as parent class.
        Transform.__init__( self, T, X, N=N, dt=dt )

    # Default print function.
    def __str__(self):
        assert self.A is not None or self.B is not None, \
            "\nERROR: RealFourier.A, or RealFourier.B has not been set...\n"
        line1 = 'Error: %.5e\n' % (-1 if self.err is None else self.err)
        line2 = '\tA.shape: (' + str(self.A.shape[0]) + ', ' + str(self.A.shape[1]) + ')\n'
        line3 = '\tB.shape: (' + str(self.B.shape[0]) + ', ' + str(self.B.shape[1]) + ')\n'
        return line1 + line2 + line3

    def serialize(self, T=None):
        # If data set is given use instead of 'default'.
        T = self.T if T is None else T

        # Create serialized set from frequency list.
        self.frequencies()
        tSin = np.sin( self.F*T )
        tCos = np.cos( self.F*T )

        # Return the serialized sets.
        return tSin, tCos

    def dft(self, verbose=0):
        # Check that system is single input.
        assert self.K == 1, \
            "\nERROR: RealFourier.dft() requires that system be single input.\n"

        # Serialize the given data set.
        if verbose:
            print( 'Creating sin and cos series over data set...' )
        tSin, tCos = self.serialize()

        # Initialize coefficient vectors.
        self.A = np.empty( (self.Nt, self.N+1) )
        self.B = np.empty( (self.Nt, self.N+1) )

        for i in range( self.Nt ):
            if verbose:
                print( 'Calculating coefficients for state space %i/%i.' % (i, self.Nt) )

            # Solve for when k=0.
            self.A[i,0] = 0
            self.B[i,0] = 1/(2*self.N)*sum( self.X[i,:] )

            # Solve for when 0 < k < N.
            for k in range( 1,self.N ):
                self.A[i,k] = 1/self.N*sum( self.X[i,:]*tSin[k,:] )
                self.B[i,k] = 1/self.N*sum( self.X[i,:]*tCos[k,:] )
                if verbose:
                    print( '\tCoefficients %i/%i: (A,B) = (%.3e, %.3e).'
                        % (k, self.N, self.A[i,k], self.B[i,k]) )

            # Solve for when k = N.
            self.A[i,self.N] = 0
            self.B[i,self.N] = 1/(2*self.N)*sum( self.X[i,:]*tCos[self.N,:] )

        # Return instance of self.
        self.powerspec()
        self.resError( self.T, self.X, save=1 )
        return self

    def dmd(self, N=None):
        # Set number of cos/sin terms.
        if N is not None:
            self.N = N

        # Serialize and stack the data set.
        TSinCos = np.vstack( self.serialize() )

        # Initialize the regressor variable and solve.
        regr = Regressor( TSinCos, self.X )
        C, _ = regr.dmd()

        # Set coefficient vectors.
        self.A = C[:,:self.K*(self.N+1)]
        self.B = C[:,self.K*(self.N+1):]

        # Return instance of self.
        self.powerspec()
        self.resError( self.T, self.X, save=1 )
        return self

    def vectors(self, t):
        # Check that system is single input.
        assert self.K == 1, \
            "\nERROR: RealFourier.vectors() requires that system be single input.\n"

        # Expand sin/cos functions around point.
        tSin, tCos = self.serialize( t )

        # Initialize vector matrices.
        Vx = np.empty( (1, 2*(self.N+1)) )
        Vy = np.empty( (1, 2*(self.N+1)) )
        V = np.zeros( (self.Nt, 2, 2*(self.N+1)) )
        for i in range( self.Nt ):
            # Calculate x and y components of vectors.
            k = 0
            for j in range( self.N+1 ):
                Vx[0,k:k+2] = np.array( [ self.B[i,j]*tCos.T[0,j], self.A[i,j]*tSin.T[0,j] ] )
                Vy[0,k:k+2] = np.array( [ self.B[i,j]*tSin.T[0,j], self.A[i,j]*tCos.T[0,j] ] )
                k += 2
            V[i,:,:] = np.vstack( (Vx, Vy) )

            # Sum vectors together for plotting.
            for j in range( 1, 2*(self.N+1) ):
                V[i,:,j] = V[i,:,j] + V[i,:,j-1]

        # Return vector list.
        return V

    def powerspec(self):
        # Quit if coefficients are not set.
        assert self.A is not None or self.B is not None, \
            "\nERROR: RealFourier.A, or RealFourier.B has not been set...\n"

        # Calculate power series from sin and cos functions.
        self.P = 1/(4*(self.N + 1))*np.array( [self.A**2, self.B**2] )
        self.R = self.P[0] + self.P[1]

        # Return instance of self.
        return self

    def solve(self, T=None):
        # Is given set is none, use default.
        T = self.T if T is None else T

        # Get serialized form of data set.
        tSin, tCos = self.serialize( T )

        # Return approximation from coefficient matrices.
        Y = self.A@tSin + self.B@tCos
        return Y

    def resError(self, T=None, X=None, save=0):
        # Quit if coefficients are not set.
        assert self.A is not None or self.B is not None, \
            "\nERROR: RealFourier.A, or RealFourier.B has not been set...\n"

        # Initialize data matrix (default: training data).
        T = self.T if T is None else T
        X = self.X if X is None else X

        # Solve for approximation of set.
        Y = self.solve( T )

        # Calculate residual error.
        err = np.linalg.norm( X - Y )**2

        # Save if requested and return.
        if save:
            self.err = err
        return err

# Class: ComplexFourier()
class ComplexFourier( Transform ):
    def __init__(self, T, X, N=None, dt=None):
        # Intialize coefficient matrices to None.
        self.Cn = None   # Complex coefficients (negative).
        self.Cp = None   # Complex coefficients (positive).
        self.P = None   # Power spectrum.

        # Initialize Transform() as parent class.
        Transform.__init__( self, T, X, N=N, dt=dt )

    # Default print function.
    def __str__(self):
        assert self.Cn is not None or self.Cp is not None, \
            "\nERROR: ComplexFourier.C has not been set...\n"
        line1 = 'Error: %.5e\n' % (-1 if self.err is None else self.err)
        line2 = '\tC.shape: (' + str(self.C.shape[0]) + ', ' + str(self.C.shape[1]) + ')\n'
        return line1 + line2

    @property
    def C(self):
        # Return cumulative coefficient list.
        return np.hstack( (self.Cn, self.Cp) )

    def serialize(self, T=None):
        # If data set is given use instead of 'default'.
        T = self.T if T is None else T

        # Create serialized set from frequency list.
        self.frequencies()
        tExpN = np.exp( -1j*self.F*T )
        tExpP = np.exp(  1j*self.F*T )

        # Return the serialized sets.
        return tExpN, tExpP

    def dft(self, verbose=0):
        # Check that system is single input.
        assert self.K == 1, \
            "\nERROR: ComplexFourier.dft() requires that system be single input.\n"

        # Serialize the given data set.
        if verbose:
            print( 'Creating sin and cos series over data set...' )
        tExpN, tExpP = self.serialize()

        # Initialize coefficient vectors.
        self.Cn = np.empty( (self.Nt, self.N+1) )
        self.Cp = np.empty( (self.Nt, self.N+1) )

        # for i in range( self.Nt ):
        #     if verbose:
        #         print( 'Calculating coefficients for state space %i/%i.' % (i, self.Nt) )

        #     # Solve for when k=0.
        #     self.A[i,0] = 0
        #     self.B[i,0] = 1/(2*self.N)*sum( self.X[i,:] )

        #     # Solve for when 0 < k < N.
        #     for k in range( 1,self.N ):
        #         self.A[i,k] = 1/self.N*sum( self.X[i,:]*tSin[k,:] )
        #         self.B[i,k] = 1/self.N*sum( self.X[i,:]*tCos[k,:] )
        #         if verbose:
        #             print( '\tCoefficients %i/%i: (A,B) = (%.3e, %.3e).'
        #                 % (k, self.N, self.A[i,k], self.B[i,k]) )

        #     # Solve for when k = N.
        #     self.A[i,self.N] = 0
        #     self.B[i,self.N] = 1/(2*self.N)*sum( self.X[i,:]*tCos[self.N,:] )

        # # Return instance of self.
        # self.powerspec()
        # self.resError( self.T, self.X, save=1 )
        return self

    def solve(self, T=None):
        # Is given set is none, use default.
        T = self.T if T is None else T

        # Get serialized form of data set.
        tExpN, tExpP = self.serialize( T )

        # Return approximation from coefficient matrices.
        Y = self.Cn@tExpN + self.Cp@tExpP
        return Y

    def RtoC(self, fvar):
        # Convert sin/cos series to pos/neg coefficient groups.
        self.Cn = 1/2*fvar.B + 1j/2*fvar.A
        self.Cp = 1/2*fvar.B - 1j/2*fvar.A

        # Return instance of self.
        return self
