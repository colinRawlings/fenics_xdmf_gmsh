# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

stripped down version of sph_plane_divide_opt_vary_gap which just calculates 
free energy for some different particle positions

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

Nh = 10 # number of different particle positions to try
lD_nm = 30 # Debye length
Tdeg = 21
w = 200 # field size

psi_mV={'Si':67,'sph':58,'poly':67}
geom_nm={'w':w,'L':160,'h':50,'a':30,'sb':60,'st':30,'dsr':Rs_mesh*lD_nm} # cylindrical symmetry
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars
        
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)

#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':2,'Si':1}}


F_int,Fc_int,u_ = ph.solve_ratchet_problem(geom_ND,params_ND,SDs,fName_tail,doPlots=True,cylSim=True)

#--- do XX
f0,axs = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()
fh.place_fig(f0)

#-- add in 3d solution
R = pickle.load(open(r'/home/cra/dolfin/particles/thD_ratchet/benchmark_cyl_results20160603-100222/XX.pickle','rb'))
axs[0].plot(R['zXX'],R['uXXz'],'or')
axs[1].plot((R['rXX']-R['geom_ND']['x0']),R['uXXr'],'or')

fh.plotXX(u_,(0,-geom_ND['sb']),(0,geom_ND['L']+geom_ND['st']),IV='x1',axs=axs[0],style='-k')
axs[0].set_xlabel('$z$ (nm)')
axs[0].set_ylabel('$\phi$')

fh.plotXX(u_,(0,geom_ND['L']-geom_ND['h']-geom_ND['a']),\
        (geom_ND['w'],geom_ND['L']-geom_ND['h']-geom_ND['a']),IV='x0',axs=axs[1],style='-k')
axs[1].set_xlabel('$r$ (nm)')
axs[1].set_ylabel('$\phi$')
axs[0].legend(('3D','axi'))

plt.draw()
