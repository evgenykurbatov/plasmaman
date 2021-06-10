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

from plasmaman.const import h, c, Angstrom, eV, AU, L_sol
from plasmaman.models.euvac import EUVAC



sol = EUVAC()
print("%d wave length bins from %g to %g Angstrom" % \
      (len(sol.lambda_min), min(sol.lambda_min), max(sol.lambda_max)))

## Scale fluxes to certain levels of solar activity
Phi_80  = sol.Phi(80)
Phi_140 = sol.Phi(140)
Phi_200 = sol.Phi(200)


##
## Integral EUV flux at the solar surface
##

L_EUV_80  = 4*np.pi*AU**2 * np.sum(Phi_80 * sol.dE)
L_EUV_140 = 4*np.pi*AU**2 * np.sum(Phi_140 * sol.dE)
L_EUV_200 = 4*np.pi*AU**2 * np.sum(Phi_200 * sol.dE)
print("L_sol = %.2e [erg/s]" % L_sol)
print("L_EUV_80  = %.2e [erg/s] = %.2e L_sol" % (L_EUV_80,  L_EUV_80/L_sol))
print("L_EUV_140 = %.2e [erg/s] = %.2e L_sol" % (L_EUV_140, L_EUV_140/L_sol))
print("L_EUV_200 = %.2e [erg/s] = %.2e L_sol" % (L_EUV_200, L_EUV_200/L_sol))


##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=1, ncols=1)

dlambda = sol.lambda_max - sol.lambda_min
ax.bar(sol.lambda_min, Phi_200/1e9, width=0.9*dlambda, align='edge', color='m')
ax.bar(sol.lambda_min, Phi_140/1e9, width=0.9*dlambda, align='edge', color='c')
ax.bar(sol.lambda_min, Phi_80/1e9,  width=0.9*dlambda, align='edge', color='b')
ax.set_xlabel(r"[Angstrom]")
ax.set_ylabel(r"[$10^9$ photons cm$^{-2}$ sec$^{-1}$]")
ax.legend(labels=[r"$P = 200$", r"$P = 140$", r"$P = 80$"])
ax2 = ax.twiny()
lambda_min_, lambda_max_ = ax.get_xlim()
ax2.set_xlim((h*c/(lambda_min_*Angstrom)/(1e3*eV), h*c/(lambda_max_*Angstrom)/(1e3*eV)))
ax2.set_xlabel(r"[keV]")

plt.tight_layout()
plt.show()
