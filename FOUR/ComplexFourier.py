
from FOUR.Transforms import *

# Class: ComplexFourier()
class ComplexFourier(Transform):
    def __init__(self, T, X, N=None, dt=None):
        # Intialize coefficient matrices to None.
        self.Cn = None   # Complex coefficients (negative).
        self.Cp = None   # Complex coefficients (positive).
        self.P = None    # Power spectrum split by p/n coefficients.

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

    def powerspec(self):
        self.P = 1/self.N*np.array( [self.Cn, self.Cp] )**2
        self.R = np.hstack( self.P )

        # Return instance of self.
        return self

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
        self.powerspec()
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
