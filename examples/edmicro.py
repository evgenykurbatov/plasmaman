# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from plasmaman.const import *
from plasmaman.edmicro import *
from plasmaman.plte import *



## From mm to Ly-limit photons
nu = np.logspace(log10(c/0.1), 2+log10(ionization_energy_1/h), 200)

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
