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
geom_nm={'w':w,'L':160,'h':50,'a':30,'sb':60,'st':30,'dsr':Rs_mesh*lD_nm}
         
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,1,Tdeg,psi_mV)



#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':2,'Si':1}}

#-- construct basis
mesh,sd_MF,bnd_MF=fh.params_to_mesh(geom_ND,fName_tail,'out.geo')
#fh.paraview_mesh(mesh,sd_MF=sd_MF)

V = FunctionSpace(mesh,'CG',2)
u_ = Function(V)
u = TrialFunction(V)
v = TestFunction(V)

#-- define equation
ds = Measure('dS')[bnd_MF]
dx = Measure('dx')[sd_MF]
dx_dry = dx(SDs['vol']['sph'])+dx(SDs['vol']['Si'])+dx(SDs['vol']['poly'])
dx_wet = dx(SDs['vol']['fluid'])  

r,z = SpatialCoordinate(mesh) 

F = params_ND['sigma_sph']*v('+')*ds(SDs['surf']['sph'])\
    + params_ND['sigma_Si']*v('+')*ds(SDs['surf']['Si'])\
    + params_ND['sigma_poly']*v('+')*ds(SDs['surf']['poly'])\
    + v/(r+1e-6)*Dx(u,0)*params_ND['er_fluid']*dx_wet\
    + v/(r+1e-6)*Dx(u,0)*params_ND['er_dry']*dx_dry\
    - 2*params_ND['n0']*sinh(u)*v*dx_wet\
    - dot(params_ND['er_fluid']*grad(u),grad(v))*dx_wet\
    - dot(params_ND['er_dry']*grad(u),grad(v))*dx_dry
    
F = action(F,u_)
J = derivative(F,u_,u)

#- n.b. with |J| added 
u_el = Constant(params_ND['er_fluid']*0.5)*dot(grad(u_),grad(u_))*2*pi*r*dx_wet+\
       Constant(params_ND['er_dry']*0.5)*dot(grad(u_),grad(u_))*2*pi*r*dx_dry
Tds = 2*Constant(params_ND['n0'])*(cosh(u_)-1.0-u_*sinh(u_))*2*pi*r*dx_wet

f = u_el - Tds

problem = NonlinearVariationalProblem(F,u_,[],J)

solver = NonlinearVariationalSolver(problem)
solver.parameters['newton_solver']['linear_solver'] ='gmres'
solver.parameters['newton_solver']['error_on_nonconvergence'] = False
solver.solve()

F_int=assemble(f) 

#fh.paraview_fld(u_) 


#--- make plots
Np = 100
f0,axs = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()
fh.place_fig(f0)


zXX = np.linspace(-geom_ND['sb'],geom_ND['L']+geom_ND['st'],Np)
rXX = np.linspace(0,geom_ND['w'],Np)

uXXz = np.zeros(Np)
uXXr = np.zeros(Np)

for q in range(Np):
    uXXz[q] = u_(0,zXX[q])
    uXXr[q] = u_(rXX[q],geom_ND['L']-geom_ND['h']-geom_ND['a'])

#-- plots
axs[0].plot(zXX*geom_ND['l0_m']*1e9,uXXz,'-k')
axs[0].set_xlabel('$z$ (nm)')
axs[0].set_ylabel('$\phi$')

axs[1].plot(rXX*geom_ND['l0_m']*1e9,uXXr,'-k')
axs[1].set_xlabel('$r$ (nm)')
axs[1].set_ylabel('$\phi$')

plt.draw()

#-- add in 3d solution
R = pickle.load(open(r'/home/cra/dolfin/particles/thD_ratchet/benchmark_cyl_results20160603-100222/XX.pickle','rb'))
axs[0].plot(R['zXX']*R['geom_ND']['l0_m']*1e9,R['uXXz'],'or')
axs[1].plot((R['rXX']-R['geom_ND']['x0'])*R['geom_ND']['l0_m']*1e9,R['uXXr'],'or')

