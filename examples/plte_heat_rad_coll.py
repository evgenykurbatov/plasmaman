# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.special

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from plasmaman.const import *
from plasmaman.plte import *



## Black-body radiative flux per unit surface:
## q_rad  = sigma_SB T_gas^4
## Collisional kinetic energy flux per unit surface:
## q_coll = v_th N_gas (3/2) k_B T_gas
## Gas number density in equillibrium  q_rad = q_coll:
## N_gas = sigma_SB T_gas^4 / [ v_th (3/2) k_B T_gas ]

N_gas = lambda T : sigma_SB * T**4 / ( rms_vel(T, m_p) * (3/2)*k_B*T )

T_gas = np.logspace(1, 5, 100)



##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=1, ncols=1)

ax_ = ax
ax_.set_title(r"Radiation v.s. collision heating rate")
ax_.loglog(T_gas, N_gas(T_gas))
ax_.set_xlabel(r"$T_\mathrm{gas}$ [K]")
ax_.set_ylabel(r"$N_\mathrm{gas}$ [cm$^{-3}$]")
ax_.text(1e2, 1e17, r"$\Gamma_\mathrm{rad} < \Gamma_\mathrm{coll}$", ha='center', va='center')
ax_.text(1e4, 1e13, r"$\Gamma_\mathrm{rad} > \Gamma_\mathrm{coll}$", ha='center', va='center')

plt.tight_layout()
plt.show()
