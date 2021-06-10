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



N_e = np.logspace(3, 16, 400)
T = np.logspace(2, 9, 200)

N_mesh, T_mesh = np.meshgrid(N_e, T)

kin_energy = 3*k_B*T_mesh
pot_energy = e**2 / debye_length(N_mesh, T_mesh)
Lambda = log(kin_energy/pot_energy)
print("min(Lambda) = %g,  max(Lambda) = %g" % (np.min(Lambda), np.max(Lambda)))



##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=1, ncols=1)

ax_ = ax
ax_.set_title(r"Coulomb's $\Lambda$")
cs = ax_.contour(log10(N_mesh), log10(T_mesh), Lambda)
ax_.clabel(cs, fmt="%g", inline=1)
ax_.set_xlabel(r"$\lg N_\mathrm{e}$ [cm$^{-3}$]")
ax_.set_ylabel(r"$\lg T$ [K]")

plt.tight_layout()
plt.show()
