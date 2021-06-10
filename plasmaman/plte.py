# -*- coding: utf-8 -*-
"""
Particle + LTE plasma

Z is atomic number
A is mass number, A >= Z
"""

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.special

from .const import *
from .pmicro import *



##
## LTE things
##


def maxwellian_avel_pdf(v, T, m):
    """Maxwellian PDF for absolute vellocity of a particle [(cm/s)^{-1}]"""
    tmp = 0.5*m/(k_B*T)
    v2 = v**2
    return 4*pi*v2 *(tmp/pi)**1.5 * exp(-tmp*v2)


def maxwellian_energy_pdf(E, T):
    """Maxwellian PDF for kinetic energy of a particle [erg^{-1}]"""
    x = E/(k_B*T)
    return 2/(k_B*T)*sqrt(x/pi) * exp(-x)


def maxwellian_energy_cdf(E, T):
    """Maxwellian CDF for kinetic energy of a particle"""
    x = E/(k_B*T)
    return sp.special.erf(x) - 2*sqrt(x/pi) * exp(-x)


def prb_vel(T, m):
    """Most probably velocity for Maxwellian distribution (maximum of it) [cm/s]"""
    return sqrt(2*k_B*T/m)


def rms_vel(T, m):
    """Root of mean square 3-D velocity for Maxwellian distribution, sqrt(<v_x^2> + <v_y^2> + <v_z^2>) [cm/s]"""
    return sqrt(3*k_B*T/m)


def plasma_gas_pressure(N_i, T, Z=1):
    """Neutral fully ionized plasma pressure [erg/cm^3]"""
    return (Z+1) * N_i * k_B*T


def sound_vel(T, m, gamma=5/3):
    """Sound velocity in a gas of uniform composition and given 'gamma' [cm/s]"""
    return sqrt(gamma/3) * rms_vel(T, m)



##
## Particle-plasma interactions
##


def plasma_fundamental_freq(N_e):
    """Plasma (or fundamental, or Langmuir) frequency for electrons [rad/s]"""
    return sqrt(4*pi*e**2*N_e/m_e)


def debye_length(N_e, T):
    """Debye length [cm]"""
    return sqrt(k_B*T/m_e) / plasma_fundamental_freq(N_e)


def coulomb_section_ei(N_e, T, Z=1):
    """Coulomb scatter section for electron-ion collisions [cm^2]"""
    kin_energy = 3*k_B*T
    pot_energy = Z*e**2 / debye_length(N_e, T)
    Lambda = log(kin_energy/pot_energy)
    return 4*pi * (Z*e**2/(3*k_B*T))**2 * Lambda


def collision_coef_ei(N_e, T, Z=1):
    """Collision coefficient for electron-ion collisions [cm^3/s]"""
    return coulomb_section_ei(N_e, T, Z) * rms_vel(T, m_e)


def free_path_time_ei(N_i, T, Z=1):
    """Free path time for electron-ion collisions in neutral fully ionized plasma [s]"""
    return 1 / (N_i * collision_coef_ei(Z*N_i, T, Z))


def ohm_conductivity_ei(N_i, T, Z=1):
    """Ohm conductivity in electron-ion collisions"""
    return e**2/m_e * Z*N_i * free_path_time_ei(N_i, T, Z)


def collision_coef_in():
    """Collision coefficient for ion-neutral collisions [cm^3/s]
    Blaes & Balbus (1994) https://ui.adsabs.harvard.edu/abs/1994ApJ...421..163B
    """
    return 1.9e-9


def collision_coef_en(T):
    """Collision coefficient for electron-neutral collisions [cm^3/s]
    Blaes & Balbus (1994) https://ui.adsabs.harvard.edu/abs/1994ApJ...421..163B
    """
    return 1e-15 * sqrt(128*k_B*T/(9*pi*m_e))


def thermal_conductivity(T):
    """Coefficient of thermal conductivity for an ionized gas [g cm/s^3/K]
    """
    return 3e9 * (T/2e6)**2.5
