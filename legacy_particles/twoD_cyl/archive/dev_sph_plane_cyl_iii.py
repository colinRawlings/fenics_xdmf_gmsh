# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap


"""

from dolfin import *
import particle_helpers as ph
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
import fenics_helpers as fh

reload(ph)
reload(fh)

plt.close('all')

parameters.form_compiler.quadrature_degree = 4
parameters.allow_extrapolation = True

fName_tail = r'/home/cra/dolfin/particles/twoD_cyl/mesh/cyl_sph_plane_tail_01.geo'
fName_log = r'/home/cra/dolfin/particles/twoD_cyl/logFile.txt'


Rs_mesh = 0.3 # cell size as a function of 1/k0 at the charge sheets

Nh = 5 # number of different particle positions to try

lD_nm = 30 # Debye length
Tdeg = 21
w = 200 # field size

psi_mV={'Si':67,'sph':58,'poly':67}
geom_nm={'w':w,'L':160,'h':50,'a':30,'sb':60,'st':30,'dsr':Rs_mesh*lD_nm} # cylindrical symmetry
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars

#--- calc hr
h_min = lD_nm # don't let the sphere actually touch the bounding wall
hr_nm = np.linspace(h_min,geom_nm['L']-geom_nm['a']*2-h_min,Nh)
        
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)

#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}

Fr_int = np.zeros(len(hr_nm))

def x_gap_nm(h,geom_nm):
    return geom_nm['L']*0.5-h-geom_nm['a']

for p in range(len(hr_nm)):
    geom_ND['h'] = hr_nm[p]*1e-9/geom_ND['l0_m']    
    Fr_int[p],Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,SDs,fName_tail,doPlots=False,cylSim=True)

#--- do XX
f0,axs = plt.subplots(nrows=1,ncols=3)
fh.make_cfigw()
fh.place_fig(f0)

#-- add in 3d solution
R = pickle.load(open(r'/home/cra/dolfin/particles/thD_ratchet/benchmark_cyl_results20160606-162211/XX.pickle','rb'))
axs[0].plot(R['IVz']['x1'],R['DVz'][0],'or')
axs[1].plot((R['IVr']['x0']-R['geom_ND']['x0']),R['DVr'][0],'or')

fh.plotXX(u_,(0,-geom_ND['sb']),(0,geom_ND['L']+geom_ND['st']),IV='x1',axs=axs[0],style='-k')
axs[0].set_xlabel('$z$ (nm)')
axs[0].set_ylabel('$\phi$')

fh.plotXX(u_,(0,geom_ND['L']-geom_ND['h']-geom_ND['a']),\
        (geom_ND['w'],geom_ND['L']-geom_ND['h']-geom_ND['a']),IV='x0',axs=axs[1],style='-k')
axs[1].set_xlabel('$r$ (nm)')
axs[1].set_ylabel('$\phi$')
axs[0].legend(('3D','axi'))

#-- add in F energy curve
axs[2].plot(x_gap_nm(hr_nm,geom_nm),Fr_int-np.min(Fr_int),'.k')
axs[2].set_xlabel('$x$ (nm)')
axs[2].set_ylabel('$F/(k_BT)$')
axs[2].plot(x_gap_nm(R['hr_nm'],geom_nm),R['Fr']-np.min(R['Fr']),'or')
axs[2].legend(('axis','3D'))

plt.draw()
