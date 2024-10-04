
from FOUR.Transforms import *

def conjugate(X):
    return np.real( X ) - np.imag( X )*1j

def complexiwave(cvar, i=0):
    # Perform power spectrum calculations if necessary.
    if cvar.R is None:
        cvar.powerspec()

    # Get parameters from series.
    cn = cvar.Cn[:,cvar.sort[:,i]]
    cp = cvar.Cp[:,cvar.sort[:,i]]

def complexmaximum(cvar):
    # Perform power spectrum calculations if necessary.
    if cvar.R is None:
        cvar.powerspec()

    # Get parameters from series variable.
    Cmax = cvar.C[:,cvar.sort[:,-1]]
    ampl = 2*np.real( np.sqrt( Cmax*conjugate( Cmax ) ) )
    freq = np.abs( cvar.w[cvar.sort[:,-1]] )
    phase = 0

    return CharacteristicWave( ampl, freq, phase )

def complexcentroid(cvar):
    assert isinstance( cvar, ComplexFourier ), \
        'ERROR: Incorrect variable type given to complexcentroid().'

    # Perform power spectrum calculations if neccesary.
    if cvar.R is None:
        cvar.powerspec()

    # Get parameters from series variable.
    N = cvar.N
    C = cvar.C
    w = cvar.w[N:]
    R = cvar.R[:,N:]

    freq = R@w/np.sum( R, axis=1 )
    phase = 0

    return CharacteristicWave( 1, freq, phase )

# Class: ComplexFourier()
class ComplexFourier(Transform):
    def __init__(self, T, X, N=None, dt=None):
        # Intialize coefficient matrices to None.
        self.Cn = None   # Complex coefficients (negative).
        self.Cp = None   # Complex coefficients (positive).

        # Initialize Transform() as parent class.
        Transform.__init__( self, T, X, N=N, dt=dt )

    @property
    def check(self):
        assert self.Cn is not None and self.Cp is not None, \
            'ComplexFourier.C has not been set.'
        return True

    @property
    def C(self):
        # Return cumulative coefficient list.
        return np.hstack( (self.Cn[:,::-1], self.Cp) )

    @property
    def w(self):
        return np.vstack( (-self.F[::-1], self.F) )

    # Default print function.
    def __str__(self):
        self.check
        line1 = 'Error: %.5e\n' % (-1 if self.err is None else self.err)
        line2 = '\tC.shape: (' + str(self.C.shape[0]) + ', ' + str(self.C.shape[1]) + ')'
        return line1 + line2

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
            # Solve for when k=0.
            self.Cn[i,0] = 1/(4*self.N)*np.sum( self.X[i,:] )
            self.Cp[i,0] = self.Cn[i,0]

            # Solve for when 0 < k < N.
            for k in range( 1,self.N ):
                self.Cn[i,k] = 1/(2*self.N)*np.sum( self.X[i,:]*tExpP[k,:] )

            # Solve for when k = N.
            self.Cn[i,self.N] = 1/(4*self.N)*np.sum( self.X[i,:]*np.cos( self.F[-1]*self.T ) )
        self.Cp[:,1:] = conjugate( self.Cn )[:,1:]

        # # Return instance of self.
        self.resError( self.T, self.X, save=1 )
        return self

    def solve(self, T=None):
        # Is given set is none, use default.
        T = self.T if T is None else T

        # Get serialized form of data set.
        tExpN, tExpP = self.serialize( T )

        # Return approximation from coefficient matrices.
        Y = self.Cn@tExpN + self.Cp@tExpP
        return Y

    def powerspec(self):
        self.check

        # Calculate power spectrum.
        self.R = np.real( 1/self.N*self.C*conjugate( self.C ) )

        # Normalize to maximum spectrum value.
        self.Rmax = np.max( self.R )
        self.R = self.R/self.Rmax

        # Create sorted list of most significant terms.
        self.sort = np.argsort( self.R, kind='quicksort' )

        # Return instance of self.
        return self

    def RtoC(self, fvar):
        # Convert sin/cos series to pos/neg coefficient groups.
        self.Cn = 1/2*fvar.B + 1j/2*fvar.A
        self.Cp = 1/2*fvar.B - 1j/2*fvar.A

        # Return instance of self.
        return self
