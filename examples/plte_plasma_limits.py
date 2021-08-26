# -*- coding: utf-8 -*-

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from plasmaman.const import *
from plasmaman.plte import *



## Log number density [cm-3]
N = np.logspace(3, 16, 400)
## Log ionization fraction
x = np.logspace(-3, 0, 200)

N_mesh, x_mesh = np.meshgrid(N, x)

## Fundamental (Langmuir) wave length
lam_pe = 2*pi*c / plasma_fundamental_freq(x_mesh*N_mesh)
## Distance between ions
lam_elem = (x_mesh*N_mesh)**(-1/3)



##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=1, ncols=1)

ax_ = ax
ax_.loglog(N, 2*pi*c / plasma_fundamental_freq(N), label="Langmuir limit")
ax_.loglog(N, N**(-1/3), label="Elementary limit")
ax_.set_xlabel(r"$N$ [cm$^{-3}$]")
ax_.set_ylabel(r"[cm]")
ax_.legend()

ax__ = ax_.twinx()
ax__.set_yscale('log')
lam_min, lam_max = ax_.get_ylim()
ax__.set_ylim(c/lam_min, c/lam_max)
ax__.set_ylabel(r"[Hz]")

plt.tight_layout()
plt.show()
