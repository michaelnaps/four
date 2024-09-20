# File: Transforms.py
# Created by: Michael Napoli
# Created on: Jul 11, 2023
# Purpose: To develop and study various methods of expanding
#   and solving real/complex Fourier series.

import sys
from os.path import expanduser
sys.path.insert( 0, expanduser('~')+'/prog/kman' )

import numpy as np
from copy import deepcopy
from FOUR.Preprocess import *
from KMAN.Regressors import *

# Class: CharacteristicWave()
# Purpose: To be used to save the characteristic wave form found through the Transform() class.
class CharacteristicWave:
    def __init__(self, ampl=None, freq=None, phase=None):
        # Initialize class variables.
        self.ampl = ampl
        self.freq = freq
        self.phase = phase

    @property
    def period(self):
        return 2*np.pi/self.freq

    def __str__(self):
        line1 = 'Characteristic wave: '
        line2 = 'x(t) = %.3f sin( %.3f t + %.3f )' % (self.ampl, self.freq, self.phase)
        return line1 + line2

    def setCharacteristics(self, ampl=None, freq=None, phase=None):
        # Update components if given values are not None.
        self.ampl = self.ampl if ampl is None else ampl
        self.freq = self.freq if freq is None else freq
        self.phase = self.phase if phase is None else phase

        # Return instance of self.
        return self

    def solve(self, T):
        assert None not in (self.freq, self.phase, self.ampl), \
            'ERROR: One (or all) of CharacteristicWave.freq/phase/ampl is not set.'
        # Return solution of wave form.
        return self.ampl*np.sin( self.freq*T + self.phase )

# Class: Transform()
class Transform:
    def __init__(self, T, X, N=None, dt=None):
        # Initialize transform variables and dimensions.
        self.T = T      # Input data.
        self.X = X      # Output data.

        self.F = None   # List of frequencies.
        self.P = None   # Power spectrum list split by coefficient members.
        self.R = None   # Power spectrum.
        self.Rmax = None    # Max value from power spectrum.

        # Dimensions and frequency parameter setup for the transform.
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

        # Initialize parent class.
        self.cwave = CharacteristicWave()
        self.cmean = CharacteristicWave()

    def setDataSets(self, T, X):
        self.__init__( T, X )
        # Return instance of self.
        return self

    def setCoefficientNumber(self, N):
        self.N = N
        # Return instance of self.
        return self

    def frequencies(self):
        # Generate frequency list.
        self.F = 2*np.pi/self.tau*np.arange( self.K*(self.N + 1) )[:,None]

        # Return instance of self.
        return self

    def centroidfreq(self, P=None):
        assert not (self.F is None and self.R is None), \
            "\nERROR: Either Transform.F or Transform.R have not been set.\n"

        self.Fmax = self.F[self.sort[:,-1]]
        self.Tmax = 2*np.pi/self.Fmax
        Cmax = np.sqrt( self.A[:,self.sort[0,-1]]**2 + self.B[:,self.sort[0,-1]]**2 )
        pmax = 2*np.pi/np.arccos( self.A[:,self.sort[0,-1]]/Cmax )
        self.cwave.setCharacteristics( Cmax, self.Fmax, pmax )

        # Mean characteristic wave.
        self.Fmean = self.R@self.F/np.sum( self.R, axis=1 )
        self.Tmean = 2*np.pi/self.Fmean
        Amean = self.R@self.A.T/np.sum( self.R, axis=1 )
        Bmean = self.R@self.B.T/np.sum( self.R, axis=1 )
        Cmean = np.sqrt( Amean**2 + Bmean**2 )
        pmean = 2*np.pi/np.arccos( Amean/Cmean )
        self.cmean.setCharacteristics( Cmean, self.Fmean, pmean )

        # Return instance of self.
        return self

    def resError(self, T=None, X=None, save=0):
        self.check

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
