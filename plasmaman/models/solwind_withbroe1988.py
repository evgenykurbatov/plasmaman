# -*- coding: utf-8 -*-
"""
Withbroe model for corona and inner solar wind. See ``[W88]_``.

.. [W88] Withbroe G L. The Temperature Structure, Mass, and Energy Flow in the Corona
   and Inner Solar Wind // ApJ 325:442 (1988)
   ADS URL: https://ui.adsabs.harvard.edu/abs/1988ApJ...325..442W
"""

import numpy as np

from ..const import R_sol



class Withbroe1988(object):
    """Withbroe model for corona and inner solar wind.

    See ``[W88]_``.

    Attributes
    ----------
    QS : class SolWind
        Quiet Sun.
    PCHmax : class SolWind
        Polar Coronal Hole in the maximum.
    ECHmin : class SolWind
        Equatorial Coronal Hole in the minimum.
    PCHmin : class SolWind
        Polar Coronal Hole in the minimum.
    """

    class SolWind:
        """
        Attributes
        ----------
        r : ndarray
            Radial coordinate [cm].
        T : ndarray
            Temperature [K].
        N_e : ndarray
            Number density of electrons [cm^{-3}].
        v : ndarray
            Radial velocity of the wind [cm s^{-1}].
        """

        ## Radial coordinate [cm]
        r = None
        ## Temperature [K]
        T = None
        ## Number density of protons [cm^{-3}]
        N = None
        ## Radial velocity of the wind [cm s^{-1}]
        v = None


    QS     = SolWind()
    PCHmax = SolWind()
    ECHmin = SolWind()
    PCHmin = SolWind()


    def __init__(self):
        import os
        location = os.path.dirname(os.path.realpath(__file__))
        data = np.loadtxt(os.path.join(location, 'solwind_withbroe1988.csv'), delimiter=';', comments='#')

        r = data[:,0] * R_sol

        ## Quiet Sun
        self.QS.r   = r
        self.QS.T   = data[:,1]
        self.QS.N_e = data[:,2]
        self.QS.v   = data[:,3]

        ## Polar Coronal Hole (max)
        self.PCHmax.r   = r
        self.PCHmax.T   = data[:,4]
        self.PCHmax.N_e = data[:,5]
        self.PCHmax.v   = data[:,6]

        ## Equatorial Coronal Hole (min)
        self.ECHmin.r   = r
        self.ECHmin.T   = data[:,7]
        self.ECHmin.N_e = data[:,8]
        self.ECHmin.v   = data[:,9]

        ## Polar Coronal Hole (min)
        self.PCHmin.r   = r
        self.PCHmin.T   = data[:,10]
        self.PCHmin.N_e = data[:,11]
        self.PCHmin.v   = data[:,12]
