# -*- coding: utf-8 -*-
"""
Withbroe model for corona and inner solar wind. See ``[W88]_``.

.. [W88] Withbroe G L. The Temperature Structure, Mass, and Energy Flow in the Corona
   and Inner Solar Wind // ApJ 325:442 (1988)
   ADS URL: https://ui.adsabs.harvard.edu/abs/1988ApJ...325..442W
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from numpy import sqrt

from plasmaman.const import G, m_p, k_B, M_sol, R_sol, AU, mbar
from plasmaman.plte import sound_vel
from plasmaman.models.solwind_withbroe1988 import Withbroe1988



wind = Withbroe1988()



##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

## rc settings (see http://matplotlib.sourceforge.net/users/customizing.html#customizing-matplotlib)
mpl.rc('font', family='serif')
mpl.rc('font', size='8.0')
mpl.rc('text', usetex=True)
mpl.rc('lines', linewidth=0.5)
mpl.rc('axes', linewidth=0.5)
mpl.rc('legend', frameon=False)
mpl.rc('legend', handlelength=2.5)

figwidth = 24.0 / 2.54           ## convert cm to in
figheight = 16.0 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=2, ncols=3)


ax_ = ax[0,0]
ax_.loglog(wind.QS.r/R_sol,     wind.QS.T,     label="QS")
ax_.loglog(wind.PCHmax.r/R_sol, wind.PCHmax.T, label="PCH (max)")
ax_.loglog(wind.ECHmin.r/R_sol, wind.ECHmin.T, label="ECH (min)")
ax_.loglog(wind.PCHmin.r/R_sol, wind.PCHmin.T, label="PCH (min)")
ax_.set_xlabel(r"$r/R_\odot$")
ax_.set_ylabel(r"$T$ [K]")
ax_.set_ylim(4e5, 2e6)
ax_.legend()


ax_ = ax[0,1]
ax_.loglog(wind.QS.r/R_sol,     wind.QS.N_e,     label="QS")
ax_.loglog(wind.PCHmax.r/R_sol, wind.PCHmax.N_e, label="PCH (max)")
ax_.loglog(wind.ECHmin.r/R_sol, wind.ECHmin.N_e, label="ECH (min)")
ax_.loglog(wind.PCHmin.r/R_sol, wind.PCHmin.N_e, label="PCH (min)")
ax_.set_xlabel(r"$r/R_\odot$")
ax_.set_ylabel(r"$N_\mathrm{e}$ [cm$^{-3}$]")
ax_.legend()


ax_ = ax[0,2]
ax_.loglog(wind.QS.r/R_sol, sqrt(2*G*M_sol/wind.QS.r), ':k', label="Escape vel.")
l = ax_.loglog(wind.QS.r/R_sol,     wind.QS.v,     label="QS")
ax_.plot(wind.QS.r/R_sol, sound_vel(wind.QS.T,     m_p, gamma=5/3), c=l[0].get_color(), ls='--')
l = ax_.loglog(wind.PCHmax.r/R_sol, wind.PCHmax.v, label="PCH (max)")
ax_.plot(wind.QS.r/R_sol, sound_vel(wind.PCHmax.T, m_p, gamma=5/3), c=l[0].get_color(), ls='--')
l = ax_.loglog(wind.ECHmin.r/R_sol, wind.ECHmin.v, label="ECH (min)")
ax_.plot(wind.QS.r/R_sol, sound_vel(wind.ECHmin.T, m_p, gamma=5/3), c=l[0].get_color(), ls='--')
l = ax_.loglog(wind.PCHmin.r/R_sol, wind.PCHmin.v, label="PCH (min)")
ax_.plot(wind.QS.r/R_sol, sound_vel(wind.PCHmin.T, m_p, gamma=5/3), c=l[0].get_color(), ls='--')
ax_.set_xlabel(r"$r/R_\odot$")
ax_.set_ylabel(r"$v$ [cm$/$s]")
ax_.legend()


ax_ = ax[1,0]
ax_.loglog(wind.QS.r/R_sol,     wind.QS.r**2 * wind.QS.N_e * wind.QS.v,             label="QS")
ax_.loglog(wind.PCHmax.r/R_sol, wind.PCHmax.r**2 * wind.PCHmax.N_e * wind.PCHmax.v, label="PCH (max)")
ax_.loglog(wind.ECHmin.r/R_sol, wind.ECHmin.r**2 * wind.ECHmin.N_e * wind.ECHmin.v, label="ECH (min)")
ax_.loglog(wind.PCHmin.r/R_sol, wind.PCHmin.r**2 * wind.PCHmin.N_e * wind.PCHmin.v, label="PCH (min)")
ax_.set_xlabel(r"$r/R_\odot$")
ax_.set_ylabel(r"$r^2 N_\mathrm{e} v$ [s$^{-1}$]")
ax_.legend()


ax_ = ax[1,1]
J = wind.QS.r**2 * wind.QS.N_e * wind.QS.v
ax_.semilogx(wind.QS.r/R_sol,     J[-1]/J,             label="QS")
J = wind.PCHmax.r**2 * wind.PCHmax.N_e * wind.PCHmax.v
ax_.semilogx(wind.PCHmax.r/R_sol, J[-1]/J, label="PCH (max)")
J = wind.ECHmin.r**2 * wind.ECHmin.N_e * wind.ECHmin.v
ax_.semilogx(wind.ECHmin.r/R_sol, J[-1]/J, label="ECH (min)")
J = wind.PCHmin.r**2 * wind.PCHmin.N_e * wind.PCHmin.v
ax_.semilogx(wind.PCHmin.r/R_sol, J[-1]/J, label="PCH (min)")
ax_.set_xlabel(r"$r/R_\odot$")
ax_.set_ylabel(r"solid angle (relative units)")
ax_.legend()


ax_ = ax[1,2]
l = ax_.loglog(wind.QS.r/R_sol,     wind.QS.N_e * wind.QS.v**2,         label="QS")
ax_.loglog(wind.QS.r/R_sol,     wind.QS.N_e * (k_B/m_p)*wind.QS.T,         c=l[0].get_color(), ls='--')
l = ax_.loglog(wind.PCHmax.r/R_sol, wind.PCHmax.N_e * wind.PCHmax.v**2, label="PCH (max)")
ax_.loglog(wind.PCHmax.r/R_sol, wind.PCHmax.N_e * (k_B/m_p)*wind.PCHmax.T, c=l[0].get_color(), ls='--')
l = ax_.loglog(wind.ECHmin.r/R_sol, wind.ECHmin.N_e * wind.ECHmin.v**2, label="ECH (min)")
ax_.loglog(wind.ECHmin.r/R_sol, wind.ECHmin.N_e * (k_B/m_p)*wind.ECHmin.T, c=l[0].get_color(), ls='--')
l = ax_.loglog(wind.PCHmin.r/R_sol, wind.PCHmin.N_e * wind.PCHmin.v**2, label="PCH (min)")
ax_.loglog(wind.PCHmin.r/R_sol, wind.PCHmin.N_e * (k_B/m_p)*wind.PCHmin.T, c=l[0].get_color(), ls='--')
ax_.set_xlabel(r"$r/R_\odot$")
ax_.set_ylabel(r"$N_\mathrm{e} v^2$ [cm$^{-1}$ s$^{-2}$]")
ax_.legend()
ax2_ = ax_.twinx()
Nv2_min, Nv2_max = ax_.get_ylim()
ax2_.set_ylim((m_p*Nv2_min/(1e-3*mbar), m_p*Nv2_max/(1e-3*mbar)))
ax2_.set_yscale('log')
ax2_.set_ylabel(r"$m_\mathrm{p} N_\mathrm{e} v^2$ [$\mu$bar]")


plt.tight_layout()
plt.show()
#plt.savefig('solwind_withbroe1988.pdf')
