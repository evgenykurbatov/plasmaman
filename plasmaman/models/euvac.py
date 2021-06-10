# -*- coding: utf-8 -*-
"""
EUVAC. See ``[RFT94]_``.

.. [RFT94] Richards P G, Fenelly J A, Torr D G. EUVAC: A solar EUV flux model
   for aeronomic calculations // Journal of Geophysical Research,
   99(A5):8981-8992 (1994).
   https://ui.adsabs.harvard.edu/abs/1994JGR....99.8981R
"""

import numpy as np

from ..const import h, c, Angstrom



class EUVAC(object):
    """Richards et al. model for solar EUV.

    See ``[RFT94]_``.

    Flux scaling:
    .. math:: F = F74113 [1 + A (P - 80)]

    Proxy:
    .. math:: P = (F_{10.7} + F_{10.7A})/2


    Attributes
    ----------
    lambda_min : ndarray
        Lower bin edges in wave length [Angstrom].
    lambda_max : ndarray
        Upper bin edges in wave length [Angstrom].
    F74113 : ndarray
        Measured reference solar spectrum at 1 AU.
        Flux per bin [10^9 photons cm^{-2} sec^{-1}].
    A : ndarray
        Flux scaling factor.
    dE : ndarray
        Energy measure in the bins [erg].


    Methods
    -------
    Phi(P)
        EUV flux [photons cm^{-2} sec^{-1}] per bin.
    F(P)
        EUV flux [erg cm^{-2} sec^{-1}] per bin.
    """

    ## Wave length bins [Angstrom]
    lamda_min = None
    lamda_max = None
    ## Flux per bin [10^9 photons cm^{-2} sec^{-1}]
    F74113 = None
    ## Multiplier for solar activity rescaling
    A = None


    def __init__(self):
        ## Load data
        import os
        location = os.path.dirname(os.path.realpath(__file__))
        self.lambda_min, self.lambda_max, self.F74113, self.A \
            = np.loadtxt(os.path.join(location, 'euvac.csv'), delimiter=';', comments='#', unpack=True)

        ## Calc energy bins
        dlambda = self.lambda_max - self.lambda_min
        self.dE = h*c/(dlambda*Angstrom) * np.log(self.lambda_max/self.lambda_min)


    def Phi(self, P):
        """Photon flux surface density

        Parameters
        ----------
        P : float
            Proxy for solar activity level.

        Returns
        -------
        array of binned fluxes [photons cm^{-2} sec^{-1}].
        """

        return self.F74113 * (1 + self.A*(P-80)) * 1e9


    def F(self, P):
        """Energy flux surface density

        Parameters
        ----------
        P : float
            Proxy for solar activity level.

        Returns
        -------
        array of binned fluxes [erg cm^{-2} sec^{-1}].
        """

        return self.Phi(P) * self.dE
