# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap

cf. directly with Stankovich96 case for figure 6a 

psi_s = e/kT
psi_p = 3e/kT

a = lD

"""

from dolfin import *
import particle_helpers as ph
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
import fenics_helpers as fh
import pickle

plt.ion()

parameters.form_compiler.quadrature_degree = 4
parameters.allow_extrapolation = True

fName_tail = r'/home/fenics/legacy/twoD_cyl/mesh/cyl_sph_plane_tail_01.geo'
fName_log = r'/home/fenics/legacy/twoD_cyl/logFile.txt'

#--- Stank96 variables (6a)
psi_s_ND = 1
psi_p_ND = 3
a_ND = 1 

#-- arb scaling vars
lD_nm = 30 # Debye length
Tdeg = 21

#-- ideally insensitive parameters
Rs_mesh = 0.3 # cell size as a function of 1/k0 at the charge sheets
Nh = 10 # number of different particle positions to try
w = 10*lD_nm # field size

const={'e': 1.602e-19,'e0': 8.85e-12,'kB': 1.38e-23,\
           'Na': 6.022e23}

psi_mV={'Si':1e3*psi_p_ND/const['e']*const['kB']*(Tdeg+273),\
    'sph':1e3*psi_s_ND/const['e']*const['kB']*(Tdeg+273),'poly':1e-3}
geom_nm={'w':w,'L':w,'h':2*lD_nm,'a':lD_nm,'sb':2*lD_nm,'st':lD_nm,'dsr':Rs_mesh*lD_nm} # cylindrical symmetry
er_dry = {'sph': 1000, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars

#--- calc hr
hr_nm = np.linspace(0.5*lD_nm,2*lD_nm,Nh)
#hr_nm = np.asarray([5*lD_nm])
        
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)

#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}

Fr_int = np.zeros(len(hr_nm))

for p in range(len(hr_nm)):
    geom_ND['h'] = hr_nm[p]*1e-9/geom_ND['l0_m']    
    Fr_int[p],Fch,Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,SDs,fName_tail,doPlots=(p==0),cylSim=True)
    print('Finished: %d/%d'%(p+1,len(hr_nm)))


#======-- Stankovich1996 soln
Hr = np.arange(0.5,2.1,0.1)
f = 100./94*np.asarray([24,21,18,15.5,13.5,12,11,10,8.5,7.5,7.0,6,6,5.5,5.0,4])
N = 100

Hrc = np.linspace(Hr[0],Hr[np.alen(Hr)-1],N)
P = np.polyfit(Hr,f,4);
fc = np.polyval(P,Hrc)

I = np.arange(N-1,-1,-1)
fclr = fc[I]

Fci = np.cumsum(fclr)*(Hrc[1]-Hrc[0])
Fci=Fci[I]

#--- convert to units of kT
Fi_kT = Fci*params_ND['er_fluid']**1.5/np.sqrt(2*params_ND['n0'])
#====================================

#--- do XX
f0,axs = plt.subplots(nrows=1,ncols=3)
# fh.make_cfigw()
# fh.place_fig(f0)

fh.plotXX(u_,(0,-geom_ND['sb']),(0,geom_ND['L']+geom_ND['st']),IV='x1',axs=axs[0],style='-k')
axs[0].set_xlabel('$z/l_0$')
axs[0].set_ylabel('$e\phi/k_B T$')

fh.plotXX(u_,(0,geom_ND['L']-geom_ND['h']-geom_ND['a']),\
        (geom_ND['w'],geom_ND['L']-geom_ND['h']-geom_ND['a']),IV='x0',axs=axs[1],style='-k')
axs[1].set_xlabel('$r/l_0$')
axs[1].set_ylabel('$e\phi/k_B T$')

#-- add in F energy curve
# f,axs2=fh.pfig()
f0,axs2 = plt.subplots(nrows=1,ncols=1)

axs2.plot(Hrc,Fi_kT+(np.interp(2*lD_nm,hr_nm,Fr_int)-np.min(Fr_int)),'-r',\
    label = 'Stankovich96')
axs2.plot(hr_nm/lD_nm,Fr_int-np.min(Fr_int),'.b',label='fenics')
axs2.legend(loc='best')
fh.cxlbl(axs2,'x/\lambda_D','F/(k_BT)')
axs2.set_title('$a=%.1f\lambda_D$, $e\phi/k_BT$: sph=%.0f, pl=%.0f'\
    %(a_ND,psi_s_ND,psi_p_ND))

plt.show(block=True)
