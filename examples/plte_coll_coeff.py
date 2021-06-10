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



N_e = [1e0, 1e8, 1e16]
T = np.logspace(2, 8, 200)



##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=1, ncols=1)

ax_ = ax
ax_.set_title(r"Collision coefficient")
ax_.loglog(T, collision_coef_in() * np.ones_like(T), '--', label=r"'in'")
ax_.loglog(T, collision_coef_en(T), '--', label=r"'en'")
for N_e_ in N_e:
    ax_.loglog(T, collision_coef_ei(N_e_, T), label=(r"'ei', $N_\mathrm{e} = 10^{%g}$ [cm$^{-3}$]" % log10(N_e_)))
ax_.legend()
ax_.set_xlabel(r"$T$ [K]")
ax_.set_ylabel(r"$\langle \sigma_{\!\dots}\,v \rangle$ [cm$^3/$s]")

plt.tight_layout()
plt.show()
