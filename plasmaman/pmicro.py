# -*- coding: utf-8 -*-
"""
Particle + particle

Z is atomic number
A is mass number, A >= Z
"""

from .const import *



def ion_mass(Z=1, A=1):
    """Ion mass"""
    return Z*m_p + (A-Z)*m_n


def atom_mass(Z=1, A=1):
    """Atom mass, rho/N_i"""
    return Z*m_e + ion_mass(Z, A)


def plasma_mean_part_mass(Z=1, A=1):
    """Mean particle mass in neutral fully ionized plasma"""
    return atom_mass(Z, A) / (Z+1)


def plasma_mass_density(N_i, Z=1, A=1):
    """Mass density in neutral fully ionized plasma"""
    return atom_mass(Z, A) * N_i
