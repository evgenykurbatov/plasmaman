# -*- coding: utf-8 -*-
"""
EUVAC. See ``[RFT94]_``.

.. [RFT94] Richards P G, Fenelly J A, Torr D G. EUVAC: A solar EUV flux model
   for aeronomic calculations // Journal of Geophysical Research,
   99(A5):8981-8992 (1994).
   https://ui.adsabs.harvard.edu/abs/1994JGR....99.8981R
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np

from plasmaman.const import m_p, AU, mbar
from plasmaman.plte import sound_vel
from plasmaman.models.solwind_parker1958 import Parker1958



wind = Parker1958()



##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=1, ncols=3)


ax_ = ax[0]
for i, T in enumerate(wind.T):
    r = wind.r[i]
    N = wind.N[i]
    ax_.loglog(r/AU, N, label=(r"$T = %.1e [K]$" % T))
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_ylabel(r"$N$ [cm$^{-3}$]")
ax_.legend()


ax_ = ax[1]
for i, T in enumerate(wind.T):
    r = wind.r[i]
    v = wind.v[i]
    l = ax_.loglog(r/AU, v, label=(r"$T = %.1e [K]$" % T))
    ax_.axhline(sound_vel(T, m_p, gamma=5/3), c=l[0].get_color(), ls='--')
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_ylabel(r"$v$ [cm s$^{-1}$]")
ax_.legend()


ax_ = ax[2]
for i, T in enumerate(wind.T):
    r = wind.r[i]
    N = wind.N[i]
    v = wind.v[i]
    ax_.loglog(r/AU, N*v**2, label=(r"$T = %.1e [K]$" % T))
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_ylabel(r"$N v^2$ [cm$^{-1}$ s$^{-2}$]")
ax_.legend()
ax2_ = ax_.twinx()
Nv2_min, Nv2_max = ax_.get_ylim()
ax2_.set_ylim((m_p*Nv2_min/(1e-3*mbar), m_p*Nv2_max/(1e-3*mbar)))
ax2_.set_yscale('log')
ax2_.set_ylabel(r"$m_\mathrm{p} N v^2$ [$\mu$bar]")


plt.tight_layout()
plt.show()
