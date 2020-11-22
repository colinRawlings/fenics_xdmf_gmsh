# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap


"""

from __future__ import absolute_import
from __future__ import print_function
import fenics_helpers as fh
import particle_helpers as ph
import numpy as np
import matplotlib.pyplot as plt
from six.moves import range

Tdeg = 21
lD_nm = 17.5  # Debye length
pH_nom = 6.7

# --- surface models
surf = dict()
surf['Silanol'] = {'pK': 7.5, 'Gamma': 9e18, 'C_stern': 2.5}  # Silica
surf['Carboxyl'] = {'pK': 4.9, 'Gamma': 0.25e18, 'C_stern': 100}  # Carboxyl

SDs = {'vol': {'sph': 3, 'poly': 2, 'Si': 1, 'fluid': 0},
       'surf': {'edge': 4, 'sph': 3, 'poly': 1, 'Si': 2}}  # domain numbering

# --- calc ND surf pars
surf_ND = dict()
for mat in ['Silanol', 'Carboxyl']:
    surf_ND[mat] = ph.calc_surf_ND_pars(surf[mat], Tdeg)

# ---- fig2 data
pH_beh = np.arange(4, 11, 1)
scaleFac = -50e-3 / 55.5
sigma_beh = {'Silanol': np.asarray([0, 2, 5, 11, 24, 47, float('nan')]) * scaleFac,
             'Carboxyl': np.asarray([2, 5.5, 12, 22, 36, 43, 44.5]) * scaleFac}
const = ph.load_const()
l0_m = ph.calc_l0_m(Tdeg)

# --- calc fig 2
pH_r = np.linspace(3, 11, 50)

sigma_ND_r = dict()
phi_d_mV_r = dict()
N_it = 0  # count total evaluations required

for mat in ['Silanol', 'Carboxyl']:
    sigma_ND_r[mat] = np.zeros(np.shape(pH_r))
    phi_d_mV_r[mat] = np.zeros(np.shape(pH_r))
    s0 = -1

    for p in range(len(pH_r)):
        sigma_ND_r[mat][p], phi_d, phi_d_mV_r[mat][p] = ph.solve_isol_pl_Qreg(pH_r[p],
                                                                              surf_ND[
                                                                                  mat],
                                                                              lD_nm, Tdeg,
                                                                              s0_ND=s0)
        s0 = sigma_ND_r[mat][p]

# --- calculate case for pH_nom
sigma_ND, phi_d, phi_d_mV = ph.solve_isol_pl_Qreg(pH_nom, surf_ND['Silanol'], lD_nm, Tdeg,
                                                  s0_ND=-1)
sigma_SI = sigma_ND * const['e'] / l0_m ** 2
print(('pH: {}, sigma: {:.2f} mC /m^2, phi_d: {:.2f} mV'.format(
    pH_nom, sigma_SI * 1e3, phi_d_mV)))

##--- display
fig, ax = fh.cfigw(ncols=3)

q = 0
for mat in ['Silanol', 'Carboxyl']:
    col = fh.chot(q, 2)
    ax[0].plot(pH_r, sigma_ND_r[mat], '-', label='fenics: ' + mat, color=col)
    ax[1].plot(pH_r, sigma_ND_r[mat] * const['e'] / l0_m ** 2, '-', \
               label='fenics: ' + mat, color=col)
    ax[1].plot(pH_beh, sigma_beh[mat], 'o', color=col,
               label='Behrens: $\lambda_D=$ 9.6nm')

    ax[2].plot(pH_r, phi_d_mV_r[mat], '-', label='fenics: ' + mat, color=col)
    q += 1

ax[0].set_xlabel('pH')
ax[0].set_ylabel('$\sigma l_0^2/e$ ')
ax[0].legend()

ax[1].set_xlabel('pH')
ax[1].set_ylabel('$\sigma$ (C/m$^2$)')
ax[1].legend()
ax[1].set_ylim((-55e-3, 0))
ax[1].set_title(
    '$n_0$: %.1e M ($\lambda_D$: %.1f nm)' % (ph.lD_nm_to_n0_M(lD_nm, Tdeg), lD_nm))

ax[2].set_xlabel('pH')
ax[2].set_ylabel('$\phi_d$ (mV)')

plt.pause(0.1)
