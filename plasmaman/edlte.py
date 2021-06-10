# -*- coding: utf-8 -*-
"""
EM field + LTE plasma

Z is atomic number
A is mass number, A >= Z
"""

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.special

from .const import *
from .edmicro import *



## Multiplier in free-free absorption coefficient
absorption_ff_mul = 4*pi*e**6/(3*m_e**2*c*h) * sqrt(2*m_e/(3*pi))
## Multiplier in free-free emission coefficient
emission_ff_mul = absorption_ff_mul * 2*h/c**2



##
## Photon-LTE plasma interactions
##


def planck_nu(T, nu):
    """Planck function [erg/cm^2/Hz/s]"""
    chi = h*nu/(k_B*T)
    return 2*(k_B*T)**3/(h*c)**2 * chi**2 * exp(-chi) / sp.special.exprel(-chi)


def dplanck_nu_dkT(T, nu):
    """Derivative of Planck function by 'k_B T' [1/cm^2/Hz/s]"""
    chi = h*nu/(k_B*T)
    return 2*(k_B*T/(h*c))**2 * chi**2 * exp(-chi) / sp.special.exprel(-chi)**2


def planck_mean(f, T):
    """Calculate Planck mean for a function 'f(nu)' and a given temperature"""
    integrand = lambda x, T : exp(x) * f(exp(x)) * 4*pi/c*planck_nu(T, exp(x))
    tmp = sp.integrate.quad(integrand, log(1e-6*k_B*T/h), log(5*k_B*T/h), args=(T), epsabs=0)[0]
    return tmp / (a_rad*T**4)


def rosseland_mean(f, T):
    """Calculate Rosseland mean for a function 'f(nu)' and a given temperature"""
    integrand = lambda x, T : exp(x) / f(exp(x)) * dplanck_nu_dkT(T, exp(x))
    tmp_1 = sp.integrate.quad(integrand, log(1e-6*k_B*T/h), log(5*k_B*T/h), args=(T), epsabs=0)[0]
    integrand = lambda x, T : exp(x) * dplanck_nu_dkT(T, exp(x))
    tmp_2 = sp.integrate.quad(integrand, log(1e-6*k_B*T/h), log(5*k_B*T/h), args=(T), epsabs=0)[0]
    return tmp_2 / tmp_1


def absorption_ff_nu(T, nu):
    """Free-free absorption coefficient [cm^5]"""
    chi = h*nu/(k_B*T)
    gaunt_ff = sqrt(3.0)/pi * sp.special.k0e(0.5*chi)
    return h**3/(k_B*T)**3.5 * absorption_ff_mul * gaunt_ff/chi**2 * sp.special.exprel(-chi)


def emission_ff_nu(T, nu):
    """Free-free emission coefficient [erg cm^3/Hz/s]"""
    ## See Ch.15 of Shu. Physics of Astrophysics I. Radiation (1991)
    chi = h*nu/(k_B*T)
    gaunt_ff = sqrt(3.0)/pi * sp.special.k0e(0.5*chi)
    return emission_ff_mul/sqrt(k_B*T) * gaunt_ff * exp(-chi)


def transition_fb_n(T, n, Z=1):
    """Free-bound transition rate to the 'n'-th state, e.g. recombination coefficient [cm^3/s]"""
    kT = k_B*T
    eta = ionization_energy(n, Z) / kT
    ## Some workaround to avoid overflow in exponent
    exp_E1_exact  = lambda eta : exp(eta) * sp.special.exp1(eta)
    exp_E1_approx = lambda eta : 0.5*log(1+2/eta)
    exp_E1 = np.array([ exp_E1_exact(eta_) if eta_ < 100 else exp_E1_approx(eta_)  for eta_ in eta ])

    ## g_n, g_e, g_plus = 2*n**2, 2, 1
    ## mul = 2*g_n/(g_e*g_plus)
    mul = 2*n**2
    mul *= sqrt(2/pi) * (kT / (m_e*c**2))**1.5 * c*sigma_bf_1
    return mul * n*eta**3 * exp_E1


def transition_fb(T, n_min=1, n_max=6, Z=1):
    """Free-bound transition rate, e.g. recombination coefficient [cm^3/s]"""
    return np.array([ transition_fb_n(T, n_, Z)  for n_ in range(n_min, n_max) ]).sum(axis=0)


def cooling_fb_n(T, n, Z=1):
    """Cooling rate in free-bound transitions to the 'n'-th state [erg cm^3/s]"""
    kT = k_B*T
    eta = ionization_energy(n, Z) / kT
    ## Some workaround to avoid overflow in exponent
    exp_E1_exact  = lambda eta : exp(eta) * sp.special.exp1(eta)
    exp_E1_approx = lambda eta : 0.5*log(1+2/eta)
    exp_E1 = np.array([ exp_E1_exact(eta_) if eta_ < 100 else exp_E1_approx(eta_)  for eta_ in eta ])

    ## g_n, g_e, g_plus = 2*n**2, 2, 1
    ## mul = 2*g_n/(g_e*g_plus)
    mul = 2*n**2
    mul *= sqrt(2/pi) * (kT / (m_e*c**2))**1.5 * c*sigma_bf_1
    return mul * kT * n*eta**3 * (1 - eta*exp_E1)


def cooling_fb(T, n, n_max=6, Z=1):
    """"""
    return np.array([ cooling_fb_n(T, n_, Z)  for n_ in range(1, n_max) ]).sum(axis=0)


def emission_fb_n(T, n, Z=1):
    """Emission coefficient for free-bound transitions to the 'n'-th state [erg cm^3/s]"""
    kT = k_B*T
    eta = ionization_energy(n, Z) / kT

    ## g_n, g_e, g_plus = 2*n**2, 2, 1
    ## mul = 2*g_n/(g_e*g_plus)
    mul = 2*n**2
    mul *= sqrt(2/pi) * (kT / (m_e*c**2))**1.5 * c*sigma_bf_1
    return mul/(4*pi) * kT * n*eta**3
