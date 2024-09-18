
from FOUR.Transforms import *

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
        line5 = '\tB.shape: (' + str(self.B.shape[0]) + ', ' + str(self.B.shape[1]) + ')\n'
        line6 = self.cwave.__str__()
        return line1 + line2 + line3 + line4 + line5 + line6

    def serialize(self, T=None):
        # If data set is given use instead of 'default'.
        T = self.T if T is None else T

        # Create serialized set from frequency list.
        self.frequencies()
        tSin = np.sin( self.F*T )
        tCos = np.cos( self.F*T )

        # Return the serialized sets.
        return tSin, tCos

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

    def autocorrelate(self, llist=None, reverse=0):
        # Select either reverse/forward AC function.
        if reverse:
            fauto = lambda tlist, l: (self.solve( tlist ), self.solve( l - self.T ))
        else:
            fauto = lambda tlist, l: (self.solve( tlist ), self.solve( self.T - l ))

        # Initialize sets.
        llist = self.T if llist is None else llist
        flist = np.empty( llist.shape )

        # Iterate through lag list and calculate correlate.
        for i, l in enumerate( llist.T ):
            f, fD = fauto( self.T, l )
            flist[:,i] = f@fD.T/(np.sqrt( f@f.T )*np.sqrt( fD@fD.T ))

        return llist, flist

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
        # Y = np.empty( (self.Nt, T.shape[1]) )
        # for i, isort in enumerate( self.sort ):
        #     Y[i] = self.A[i,isort[-N:]]@tSin[isort[-N:]] + self.B[i,isort[-N:]]@tCos[isort[-N:]]
        # return Y
        return self.A@tSin + self.B@tCos

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

    def CtoR(self, cvar):
        self.A = np.imag( cvar.Cn - cvar.Cp )
        self.B = np.real( cvar.Cn + cvar.Cp )

        # Return instance of self.
        return self

class Characterize:
    def __init__(self, fvar=None, cvar=None):
        # Check and save Fourier variable and type.
        assert not (fvar is None and cvar is None), \
            'ERROR: Characterize class requires RealFourier() or ComplexFourier() variables.'
        self.fvar = cvar if fvar is None else fvar
        self.type = 'complex' if fvar is None else 'real'

    def centroidalwave(self):
        fvar = self.fvar

        if self.type == 'real':
            print( '--- real' )
            freq = fvar.R@fvar.F/np.sum( fvar.R, axis=1 )
            period = 2*np.pi/freq

            A = fvar.R@fvar.A.T/np.sum( fvar.R, axis=1 )
            B = fvar.R@fvar.B.T/np.sum( fvar.R, axis=1 )
            C = np.sqrt( A**2 + B**2 )
            ampl = 2*np.pi/np.arccos( A/C )

            wave = CharacteristicWave( ampl, freq, period )

        elif self.type == 'complex':
            print( '--- complex' )
            F = np.hstack( (-fvar.F, fvar.F) )

            print( F.shape, fvar.R.shape )
            freq = fvar.R@F/np.sum( fvar.R, axis=1 )
            period = 2*np.pi/freq

            wave = CharacteristicWave( 1, freq, period )

        return wave
