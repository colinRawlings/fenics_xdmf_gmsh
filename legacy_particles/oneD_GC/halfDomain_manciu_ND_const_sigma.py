# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:30:19 2015

@author: cra
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from dolfin import *

#-- constants  
const={'e': 1.602e-19,'e0': 8.85e-12,'kB': 1.38e-23,\
       'Na': 6.022e23}

#-- dim variables
er_fluid = 80
Tkelv = 300
ce = 0.01*1000*const['Na'] # electrolyte concentration
KD = 1
N = 6.25E15

Nl = 10
lr_m = np.linspace(1e-9,10e-9,Nl)

#----- derived
kD_m = np.sqrt(const['e']**2*ce/(er_fluid*const['e0']*const['kB']*Tkelv))
lD_m = 1./kD_m

sigma_dim = -const['e']*N #-- assuming large K


#--- pars
Rmesh = 10 # elements/Debye length

#-- calc ND variables
l0_m = const['e']**2/(const['e0']*const['kB']*Tkelv)
n0 = er_fluid/(2*lD_m**2*l0_m) 
n0_ND = n0*l0_m**3

k0 = np.sqrt(n0*const['e']**2/(er_fluid*const['e0']*const['kB']*Tkelv))
k0_ND = k0*l0_m

sigma_ND = sigma_dim*l0_m**2/const['e']

lr_ND = lr_m/l0_m
lD_ND = lD_m/l0_m

#---- mesh
mesh = IntervalMesh(int(Rmesh*np.max(lr_ND/2./lD_ND)),0,lr_ND[8]/2)
s0 = Constant(sigma_ND)

V = FunctionSpace(mesh,'CG',2)
u_ = Function(V)
u = TrialFunction(V)
v = TestFunction(V)

#--- setup BC
bound_parts = MeshFunction('size_t',mesh,mesh.topology().dim()-1)
bound_parts.set_all(0) # init  

class LeftBound(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x<1e3*DOLFIN_EPS

lb = LeftBound()
lb.mark(bound_parts,1)
ds = Measure('ds')[bound_parts]

i = indices(1)

F = s0/er_fluid*v*ds(1)-dot(grad(u),grad(v))*dx-2*n0/er_fluid*sinh(u)*v*dx

F = action(F,u_)
J = derivative(F,u_,u)

problem = NonlinearVariationalProblem(F,u_,[],J)
solver = NonlinearVariationalSolver(problem)
solver.solve()


##-- anl solution
#kv = k*sqrt(2)
#C = np.log(np.tanh(u_(0)/4))
#
#u_anl = linspace(0,u_(0),100)
#x_num = linspace(0,R_domain/k-DOLFIN_EPS,100)
#x_anl = u_anl.copy()
#u_num = u_anl.copy()
#for p in range(np.alen(u_anl)):
#    x_anl[p] = (C-np.log(np.tanh(u_anl[p]/4)))/kv
#    u_num[p] = u_(x_num[p])
#
#plt.figure()
#plt.plot(x_num,u_num,'or',x_anl,u_anl,'-k')
#plt.legend(('fem','anl'))
#
#
##-- integrals
##U_el_int = assemble(-Constant(n0)*u_*sinh(u_)*dx)
#U_el_int = assemble(Constant(er*0.5)*Dx(u_,0)**2*dx)
#T_dS_int = assemble(2*Constant(n0)*(cosh(u_)-1-u_*sinh(u_))*dx)
#
#U_el_anl = sqrt(er*n0*8)*(np.cosh(u_(0)/2)-1)
#min_T_dS_anl = np.sqrt(er*n0*8)*(3-3*np.cosh(u_(0)/2)+u_(0)*np.sinh(u_(0)/2))
#
#print '-TdS: num:', -T_dS_int,' anl: ', min_T_dS_anl
#print 'U_el: num:', U_el_int,' anl: ', U_el_anl