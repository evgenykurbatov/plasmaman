# -*- coding: utf-8 -*-
"""
EM field + particle

Z is atomic number
A is mass number, A >= Z
"""

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

from .const import *



## Ionization energy for hydrogen atom from the 1S energy level [erg]
ionization_energy_1 = (e**2/(hbar*c))**2 * m_e*c**2/2

## Electron-photon scatter section [cm^2]
section_es = 8*pi/3 * (e**2/(m_e*c**2))**2



##
## Photon-particle interactions
##


def ionization_energy(n, Z=1):
    """Ionization energy for hydrogen-like atom from the 'n'-th energy level [erg]"""
    return ionization_energy_1 * (Z/n)**2


def gaunt_bf(n, nu, Z=1):
    """Gaunt factor for bound-free transition from the 'n'-th energy level
    Walter A, Cox A N. Allen's Astrophysical Quantities (eq. 6.22 at p.133)
    """
    u = h*nu / ionization_energy(n, Z)
    tmp = 1 / (u+1)**0.6667
    return 1 + 0.1728*(u-1)*tmp/n - 0.0496*(u**2 + 1.3333*u + 1)*tmp**2/n
    #return (1 - (0.1728 + 0.0496)/n) * np.ones_like(nu)


def section_bf_nu(n, nu, Z=1):
    """Bound-free scatter section from the 'n'-th level [cm^2]"""
    E = h*nu
    chi_E = ionization_energy(n, Z) / E
    return n*gaunt_bf(n, nu, Z) * sigma_bf_1/Z**2 * chi_E**3 * np.heaviside(1 - chi_E, np.ones_like(E))


def section_fb_nu(n, nu, v, Z=1):
    """Free-bound scatter section to the 'n'-th level, by an electron [cm^2]"""
    ## g_n, g_e, g_plus = 2*n**2, 2, 1
    ## mul = 2*g_n/(g_e*g_plus)
    mul = 2*n**2
    return mul * (h*nu/(m_e*c*v))**2 * section_bf_nu(n, nu, Z)


def gaunt_ff(nu, v, Z=1):
    """Gaunt factor for free-free transitions
    Rybicki G B, Lightman A P.Radiative Processes in Astrophysics (eq. 5.12 at p.159)
    """
    b_max = v/(2*pi*nu)
    b_min_1 = 4*Z*e**2/(pi*m_e*v**2)
    b_min_2 = h/(m_e*v)
    b_min = np.where(b_min_1 > b_min_2, b_min_1, b_min_2)
    return sqrt(3)/pi * log(b_max/b_min)



##
## Particle in magnetic field
##


def cyclotron_freq(B, m=m_e):
    """Cyclotronic frequency [rad/s]"""
    return e*B/(m*c)


def cyclotron_radius(v, B, m=m_e):
    """Cyclotronic radius [cm]"""
    return v / cyclotron_freq(B, m)



##
## If the source is executed as a main program
##

if __name__ == "__main__":

    ## From mm to Ly-limit photons
    nu = np.logspace(log10(c/0.1), 2+log10(ionization_energy_1/h), 200)

    from plte import *
    T = 1e4  ## [K]
    v = prb_vel(T, m_H)



    ##
    ## Plot
    ##

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=1, ncols=2)

    ax_ = ax[0]
    ax_.set_title("Gaunt factors for 'bf' transitions")
    ax_.semilogx(nu, gaunt_ff(nu, 1e6), 'k', label=r"ff")
    for n in range(1, 6):
        ax_.semilogx(nu, gaunt_bf(n, nu), label=(r"n = %d" % n))
    ax_.set_xlabel(r"$\nu$ [Hz]")
    ax_.legend()

    ax_ = ax[1]
    ax_.set_title("Scatter sections\nfor 'bf' (solid) and 'fb' (dashed) transitions")
    ax_.loglog(nu, section_es * np.ones_like(nu), ':k', label=r"es")
    for n in range(1, 6):
        l = ax_.loglog(nu, section_bf_nu(n, nu), ls='-', label=(r"n = %d" % n))
        ax_.loglog(nu, section_fb_nu(n, nu, v), ls='--', c=l[0].get_color())
    ax_.set_xlabel(r"$\nu$ [Hz]")
    ax_.set_ylabel(r"[cm$^2$]")
    ax_.legend()

    plt.tight_layout()
    plt.show()
