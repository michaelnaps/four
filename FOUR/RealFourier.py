
from FOUR.Transforms import *

def evenness(fvar):
    asum = np.abs( np.sum( fvar.A ) )
    bsum = np.abs( np.sum( fvar.B ) )
    return bsum/(asum + bsum)

def realcentroid(fvar):
    assert isinstance( fvar, RealFourier ), \
        'ERROR: Incorrect variable type given to realcentroid().'

    # Perform power spectrum calculations if neccesary.
    if fvar.R is None:
        fvar.powerspec()

    # Get parameters from series variable.
    A = fvar.A
    B = fvar.B
    w = fvar.w
    R = fvar.R

    # Calculate centroid frequency from power spectrum.
    freq = R@w/np.sum( R, axis=1 )

    a = R@A.T/np.sum( R, axis=1 )
    b = R@B.T/np.sum( R, axis=1 )
    ampl = np.sqrt( a**2 + b**2 )

    # Get weighted phase.
    phase = R@np.arccos( A.T/np.sqrt( A.T**2 + B.T**2 ) )/np.sum( R, axis=1 )
    # plist = R@np.arctan( A.T/B.T )/np.sum( R, axis=1 )

    wave = CharacteristicWave( ampl, freq, phase )

    return wave

def realiwave(fvar, i=0, wave_type='cos'):
    # Perform power spectrum calculations if necessary.
    if fvar.R is None:
        fvar.powerspec()

    # Get parameters from series.
    a = fvar.A[:,fvar.sort[:,i]]
    b = fvar.B[:,fvar.sort[:,i]]
    freq = fvar.w[fvar.sort[:,i]]
    ampl = np.sqrt( a**2 + b**2 )

    # Calculate phase based on selection of form.
    if wave_type == 'sin':
        phase = np.arctan( b/a ) if np.abs( a ) > 0 else np.sign( b )*np.pi/2
    elif wave_type == 'cos':
        phase = np.arctan( -a/b ) if np.abs( b ) > 0 else np.sign( -a )*np.pi/2

    return CharacteristicWave( ampl, freq, phase, wave_type=wave_type )

def realmaximum(fvar, wave_type='cos'):
    return realiwave( fvar, i=-1, wave_type=wave_type )

def phasedistr(fvar, wave_type='cos'):
    plist = np.empty( (fvar.K, fvar.N + 1) )
    for i in range( fvar.N + 1 ):
        plist[:,i] = realiwave( fvar, i=i, wave_type=wave_type ).phase
    return plist

def perturbseries(fvar, eps, imin=0, imax=1, unit_time=1):
    imax = min( imax, fvar.N + 1 ) if imin < imax else imin + 1

    # Perform power spectrum calculations if necessary.
    if fvar.R is None:
        fvar.powerspec()

    # Create copy of series.
    fptb = deepcopy( fvar )

    # Create index list.
    ilist = [-(i + 1) for i in range( imin, imax )]

    # Use appropriate unit.
    if unit_time:
        wlist  = fvar.F
    else:
        wlist = np.ones( fvar.F.shape )

    for i in ilist:
        j = fvar.sort[:,i]
        a, b, w = fvar.A[:,j], fvar.B[:,j], wlist[j]
        if np.linalg.norm( [a, b] ) == 0:
            fptb.A[:,j] = fptb.B[:,j] = 0
        else:
            fptb.A[:,j] = a*np.cos( w*eps ) - b*np.sin( w*eps )
            fptb.B[:,j] = a*np.sin( w*eps ) + b*np.cos( w*eps )

    fptb.powerspec()
    fptb.resError( save=1 )
    return fptb

def offsetseries(fvar, phi, unit_time=1):
    return perturbseries( fvar, eps=phi, imax=np.inf, unit_time=unit_time )

# Class: RealFourier()
class RealFourier( Transform ):
    def __init__(self, T, X, N=None, dt=None):
        # Intialize coefficient matrices to None.
        self.A = None   # sin(t) coefficients.
        self.B = None   # cos(t) coefficients.
        self.P = None   # Power spectrum split by A/B coefficients.

        # Initialize Transform() as parent class.
        Transform.__init__( self, T, X, N=N, dt=dt )

    @property
    def check(self):
        assert self.A is not None or self.B is not None, \
            'RealFourier.A or RealFourier.B has not been set.'
        return True

    @property
    def w(self):
        return self.F

    # Default print function.
    def __str__(self):
        self.check
        line1 = 'Error: %.5e\n' % (-1 if self.err is None else self.err)
        # line2 = 'Centroid frequency: ' + str( self.Fmean.T ) + '\n'
        # line3 = 'Average period:     ' + str( self.Tmean.T ) + '\n'
        line4 = '\tA.shape: (' + str(self.A.shape[0]) + ', ' + str(self.A.shape[1]) + ')\n'
        line5 = '\tB.shape: (' + str(self.B.shape[0]) + ', ' + str(self.B.shape[1]) + ')'
        # line6 = self.cwave.__str__()
        return line1 + line4 + line5

    def __deepcopy__(self, memo={}):
        # Create new copies of variable parameters.
        T = deepcopy( self.T )
        X = deepcopy( self.X )
        N = deepcopy( self.N )
        dt = deepcopy( self.dt )

        # Initialize new variable using copied parameters.
        fnew = RealFourier( T, X, N=N, dt=dt )

        # Set series coefficients and power spectrum.
        fnew.A = deepcopy( self.A )
        fnew.B = deepcopy( self.B )
        fnew.P = deepcopy( self.P )
        fnew.F = deepcopy( self.F )
        fnew.sort = deepcopy( self.sort )

        # Return new copy...
        return fnew

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
        # Check that system is solved and single input.
        self.check
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

        # # Return approximation from coefficient matrices.
        # Y = np.empty( (self.Nt, T.shape[1]) )
        # for i, isort in enumerate( self.sort ):
        #     Y[i] = self.A[i,isort[-N:]]@tSin[isort[-N:]] + self.B[i,isort[-N:]]@tCos[isort[-N:]]
        # return Y
        return self.A@tSin + self.B@tCos

    def powerspec(self):
        self.check

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

    def CtoR(self, cvar):
        self.A = np.imag( cvar.Cn - cvar.Cp )
        self.B = np.real( cvar.Cn + cvar.Cp )

        # Return instance of self.
        return self
