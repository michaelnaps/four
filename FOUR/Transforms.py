# File: Transforms.py
# Created by: Michael Napoli
# Created on: Jul 11, 2023
# Purpose: To develop and study various methods of expanding
#   and solving real/complex Fourier series.

import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/kman' )

import numpy as np
from KMAN.Regressors import *

def conjugate(X):
    return np.real( X ) - np.imag( X )*1j

def autocorrelate(cvar=None, fvar=None, llist=None):
    assert cvar is not None or fvar is not None, \
        'ERROR: Autocorrelate requires ComplexFourier or RealFourier variable.'

    # Convert to complex Fourier to match notes.
    if cvar is not None:
        tvar = cvar
        fac = lambda f, fD: f@conjugate( fD ).T
    elif fvar is not None:
        tvar = fvar
        fac = lambda f, fD: f@fD.T

    # Initialize sets.
    llist = tvar.T if llist is None else llist
    flist = np.empty( llist.shape, dtype=complex )

    # Iterate through lag list and calculate correlate.
    for i, l in enumerate( llist.T ):
        f = tvar.solve( fvar.T )
        fD = tvar.solve( fvar.T - l )
        flist[:,i] = fac( f, fD )
    flist = flist/flist[:,0]

    return llist, flist

def reverse_autocorrelate(cvar=None, fvar=None, llist=None):
    assert cvar is not None or fvar is not None, \
        'ERROR: Autocorrelate requires ComplexFourier or RealFourier variable.'

    # Convert to complex Fourier to match notes.
    if cvar is None:
        cvar = ComplexFourier( fvar.T, fvar.X ).RtoC( fvar )

    # Initialize sets.
    llist = cvar.T if llist is None else llist
    flist = np.empty( llist.shape, dtype=complex )

    # Iterate through lag list and calculate correlate.
    for i, l in enumerate( llist.T ):
        f = cvar.solve( fvar.T )
        fC = cvar.solve( -(fvar.T - l) ).T
        flist[:,i] = f@conjugate( fC )
    flist = flist/np.max( flist, axis=1 )

    return llist, flist

# Class: Transform()
class Transform:
    def __init__(self, T, X, N=None, dt=None):
        # Initialize transform variables and dimensions.
        self.T = T      # Input data.
        self.X = X      # Output data.
        self.F = None   # List of frequencies.
        self.R = None   # Power spectrum values.
        self.K = T.shape[0]
        self.N = round( T.shape[1]/2 ) if N is None else N
        self.Nt = X.shape[0]
        self.dt = self.T[0,1] - self.T[0,0] if dt is None else dt
        self.tau = 2*self.N*self.dt
        self.err = None

        # Values from spectral analysis.
        self.Fmax = None
        self.Fmean = None
        self.Tmax = None
        self.Tmean = None

    def setDataSets(self, T, X):
        self.__init__( T, X )
        # Return instance of self.
        return self

    def setCoefficientNumber(self, N):
        self.N = N
        # Return instance of self.
        return self

    def frequencies(self, scale=None):
        # Scale may be used to change frequency spacing.
        scale = 1 if scale is None else scale

        # Generate frequency list.
        self.F = scale*2*np.pi/self.tau*np.arange( self.K*(self.N + 1) )[:,None]

        # Return instance of self.
        return self

    def centroidfreq(self, P=None):
        assert not (self.F is None and self.R is None), \
            "\nERROR: Either Transform.F or Transform.R have not been set.\n"

        # Frequency from weighted mean and spectral max.
        self.Fmax = self.F[self.sort[:,-1]]
        self.Fmean = self.R@self.F/np.sum( self.R, axis=1 )

        # Period from weighted mean and spectral max.
        self.Tmax = 2*np.pi/self.Fmax
        self.Tmean = 2*np.pi/self.Fmean

        # Return instance of self.
        return self

# Class: RealFourier()
class RealFourier( Transform ):
    def __init__(self, T, X, N=None, dt=None):
        # Intialize coefficient matrices to None.
        self.A = None   # sin(t) coefficients.
        self.B = None   # cos(t) coefficients.
        self.P = None   # Power spectrum split by A/B coefficients.
        self.Rmax = None    # Maximum spectral coefficient.

        # Initialize Transform() as parent class.
        Transform.__init__( self, T, X, N=N, dt=dt )

    # Default print function.
    def __str__(self):
        assert not (self.A is None and self.B is None), \
            "\nERROR: RealFourier.A, or RealFourier.B has not been set...\n"
        line1 = 'Error: %.5e\n' % (-1 if self.err is None else self.err)
        line2 = 'Centroid frequency: ' + str( self.Fmean.T ) + '\n'
        line3 = 'Average period:     ' + str( self.Tmean.T ) + '\n'
        line4 = '\tA.shape: (' + str(self.A.shape[0]) + ', ' + str(self.A.shape[1]) + ')\n'
        line5 = '\tB.shape: (' + str(self.B.shape[0]) + ', ' + str(self.B.shape[1]) + ')'
        return line1 + line2 + line3 + line4 + line5

    def serialize(self, T=None, scale=None):
        # If data set is given use instead of 'default'.
        T = self.T if T is None else T

        # Create serialized set from frequency list.
        self.frequencies( scale=scale )
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
            self.B[i,0] = 1/(2*self.N)*np.sum( self.X[i,:] )

            # Solve for when 0 < k < N.
            for k in range( 1,self.N ):
                self.A[i,k] = 1/self.N*np.sum( self.X[i,:]*tSin[k,:] )
                self.B[i,k] = 1/self.N*np.sum( self.X[i,:]*tCos[k,:] )
                if verbose:
                    print( '\tCoefficients %i/%i: (A,B) = (%.3e, %.3e).'
                        % (k, self.N, self.A[i,k], self.B[i,k]) )

            # Solve for when k = N.
            self.A[i,self.N] = 0
            self.B[i,self.N] = 1/(2*self.N)*np.sum( self.X[i,:]*tCos[self.N,:] )

        # Return instance of self.
        self.powerspec()
        self.centroidfreq()
        self.resError( self.T, self.X, save=1 )
        return self

    def dmd(self, N=None):
        # Set number of cos/sin terms.
        self.N
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
        self.centroidfreq()
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

        # Divide all coefficients by maximum and sub for power spectrum.
        self.R = np.sum( self.P, axis=0 )
        self.Rmax = np.max( self.R )

        # Normalize components to power spectrum.
        self.R = self.R/self.Rmax
        self.P = self.P/self.Rmax

        # Create sorted list of most significant coefficient terms.
        self.sort = np.argsort( self.R, kind='quicksort' )

        # Return instance of self.
        return self

    def solve(self, T=None, N=None):
        # Truncate signal list if requested
        N = self.N if N is None else N
        assert N <= self.N, \
            "\nERROR: Requested number of coefficients exceeds length of current list.\n"

        # Return solution over training data unless given.
        T = self.T if T is None else T

        # Get serialized form of data set.
        tSin, tCos = self.serialize( T )

        # Return approximation from coefficient matrices.
        Y = np.empty( (self.Nt, T.shape[1]) )
        for i, isort in enumerate( self.sort ):
            Y[i] = self.A[i,isort[-N:]]@tSin[isort[-N:]] + self.B[i,isort[-N:]]@tCos[isort[-N:]]
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
        self.Cn = np.empty( (self.Nt, self.N+1), dtype=complex )
        self.Cp = np.empty( (self.Nt, self.N+1), dtype=complex )

        for i in range( self.Nt ):
            if verbose:
                print( 'Calculating coefficients for state space %i/%i.' % (i, self.Nt) )

            # Solve for when k=0.
            self.Cn[i,0] = 1/(2*self.N)*sum( self.X[i,:] )
            self.Cp[i,0] = 1/(2*self.N)*sum( self.X[i,:] )

            # Solve for when 0 < k < N.
            for k in range( 1,self.N ):
                self.Cn[i,k] = 1/(2*self.N)*sum( self.X[i,:]*tExpN[k,:] )
                self.Cp[i,k] = 1/(2*self.N)*sum( self.X[i,:]*tExpP[k,:] )
                if verbose:
                    print( '\tCoefficients %i/%i: (Cn,Cp) = (' % (k, self.N),
                        self.Cn[i,k], ',', self.Cp[i,k], ')' )

            # Solve for when k = N.
            self.Cn[i,self.N] = 0
            self.Cp[i,self.N] = 1/(2*self.N)*sum( self.X[i,:]*tCos[self.N,:] )

        # # Return instance of self.
        # self.powerspec()
        # self.resError( self.T, self.X, save=1 )
        return self

        # Return instance of self.
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
