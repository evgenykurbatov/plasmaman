# -*- coding: utf-8 -*-
"""
Parker model for solar wind. See ``[P58]_``.

.. [P58] Parker E N. Dynamics of the Interplanetary Gas and Magnetic Fields //
   ApJ 128:664 (1958).
   https://ui.adsabs.harvard.edu/abs/1958ApJ...128..664P
"""

import numpy as np



class Parker1958(object):
    """Parker model for solar wind.

    See ``[P58]_``.

    Attributes
    ----------
    T : list, (7,)
        Temperature grid [K].
    r : ndarray
        Radial coordinate [cm].
    N : ndarray
        Number density of protons [cm^{-3}].
    v : ndarray
        Radial velocity of the wind [cm s^{-1}].
    """

    ## Temperature grid [K]
    T = [0.5e6, 0.75e6, 1e6, 1.5e6, 2e6, 3e6, 4e6]

    ## Radial coordinate [cm]
    r = None
    ## Number density of protons [cm^{-3}]
    N = None
    ## Radial velocity of the wind [cm s^{-1}]
    v = None


    def __init__(self):
        import os
        location = os.path.dirname(os.path.realpath(__file__))
        data = np.loadtxt(os.path.join(location, 'solwind_parker1958.csv'), delimiter=';', comments='#')

        self.r = data[:,0::3].T
        self.N = data[:,1::3].T
        self.v = data[:,2::3].T
