# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.integrate

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from plasmaman.const import *
from plasmaman.edlte import *



## From mm to X-ray photons
nu = np.logspace(log10(c/0.1), 2+log10(ionization_energy_1/h), 100)

##
N_i = 1e12
## From PPD to BH disks
T = np.logspace(2, 6, 100)

## Planck-averaged free-free scatter section
section_ff_P = np.array([ planck_mean(lambda nu : N_i*absorption_ff_nu(T_, nu), T_)
                          for T_ in T ])

## Rosseland-averaged free-free and electron scatter section
section_ff_R = np.array([ rosseland_mean(lambda nu : N_i*absorption_ff_nu(T_, nu) + section_es, T_)
                          for T_ in T ])

## Planck-averaged bound-free scatter section
section_bf_P_1 = np.array([ planck_mean(lambda nu : section_bf_nu(1, nu), T_)
                            for T_ in T ])

## Planck-averaged bound-free scatter section
section_bf_P_2 = np.array([ planck_mean(lambda nu : section_bf_nu(2, nu), T_)
                            for T_ in T ])

##
emission_ff_XUV = np.array([ sp.integrate.quad(lambda x, T_ : emission_ff_nu(T_, exp(x)),
                                               log(ionization_energy_1/h), log(100*ionization_energy_1/h),
                                               args=(T_), epsabs=0)[0]
                             for T_ in T ])

##
emission_fb_XUV = np.array([ [ emission_fb_n(T_, n)
                               for n in range(1, 6) ]
                             for T_ in T ]).sum(axis=1)



##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=2, ncols=2)

ax_ = ax[0,0]
for T_ in [1e2, 1e3, 1e4, 1e5, 1e6]:
    ax_.loglog(nu, nu*planck_nu(T_, nu), label=(r"$T = %g$ [K]" % T_))
ax_.set_xlabel(r"$\nu$ [Hz]")
ax_.set_ylabel(r"$\nu B_\nu(T)$ [erg$/$cm$^2/$s]")
ax_.set_ylim(nu[0]*planck_nu(T[0], nu[0]), max(nu*planck_nu(T[-1], nu)))
ax_.legend()

ax_ = ax[0,1]
l = ax_.loglog(T, ionization_energy_1*transition_fb(T, 1), label=r"$h \nu_1 \alpha_\mathrm{A}$")
ax_.loglog(T, ionization_energy_1*transition_fb(T, 2), label=r"$h \nu_1 \alpha_\mathrm{B}$",
           c=l[0].get_color(), ls='--')
l = ax_.loglog(T, cooling_fb(T, 1), label=r"$k_\mathrm{B}T \beta_\mathrm{A}$")
ax_.loglog(T, cooling_fb(T, 2), label=r"$k_\mathrm{B}T \beta_\mathrm{B}$",
           c=l[0].get_color(), ls='--')
ax_.set_xlabel(r"$T$ [K]")
ax_.set_ylabel(r"[erg cm$^3/$s]")
ax_.legend()

ax_ = ax[1,0]
ax_.loglog(T, section_es * np.ones_like(T), ls=':', label=r"$\sigma^\mathrm{es}$")
ax_.loglog(T, section_ff_P, label=(r"$\sigma_\mathrm{P}^\mathrm{ff}$, $N_\mathrm{i} = %g$ cm$^{-3}$" % N_i))
ax_.loglog(T, section_ff_R, label=(r"$\sigma_\mathrm{R}^\mathrm{ff}$, $N_\mathrm{i} = %g$ cm$^{-3}$" % N_i))
ax_.loglog(T, section_bf_P_1, label=r"$\sigma_\mathrm{P,1}^\mathrm{bf}$")
ax_.loglog(T, section_bf_P_2, label=r"$\sigma_\mathrm{P,2}^\mathrm{bf}$")
ax_.set_xlabel(r"$T$ [K]")
ax_.set_ylabel(r"[cm$^2$]")
ax_.legend()

ax_ = ax[1,1]
alpha_Spitzer = lambda T : 2.07e-11/sqrt(T) * np.where(log10(T) < 3, 3.0, 1.5)
ax_.loglog(T, alpha_Spitzer(T), ':k', label=r"Spitzer (1978)")
ax_.loglog(T, transition_fb(T, 1), label=r"$\alpha_\mathrm{A}$")
ax_.loglog(T, transition_fb(T, 2), label=r"$\alpha_\mathrm{B}$")
ax_.set_xlabel(r"$T$ [K]")
ax_.set_ylabel(r"[cm$^3/$s]")
ax_.legend()

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
