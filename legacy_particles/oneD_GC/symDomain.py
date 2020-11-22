# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:30:19 2015

@author: cra
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from dolfin import *

#-- variables  (all ND)
n0 = 3.
er = 2.
sigma = 40.

k = sqrt(n0/er)

L = 1.

mesh = IntervalMesh(200,0,0.5*L)
s0 = Constant(sigma)

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

F = sigma/er*v*ds(1)-dot(grad(u),grad(v))*dx-2*n0/er*sinh(u)*v*dx

F = action(F,u_)
J = derivative(F,u_,u)

problem = NonlinearVariationalProblem(F,u_,[],J)
solver = NonlinearVariationalSolver(problem)
solver.solve()


#-- anl solution
kv = k*sqrt(2)
C = np.log(np.tanh(u_(0)/4))

u_anl = linspace(0,u_(0),100)
x_num = linspace(0,0.5*L,100)
x_anl = u_anl.copy()
u_num = u_anl.copy()
for p in range(np.alen(u_anl)):
    x_anl[p] = (C-np.log(np.tanh(u_anl[p]/4)))/kv
    u_num[p] = u_(x_num[p])

#plt.figure()
plt.plot(x_num,u_num,'or',x_anl,u_anl,'-k')
plt.legend(('fem','anl'))


#-- integrals
#U_el_int = assemble(-Constant(n0)*u_*sinh(u_)*dx)
U_el_int = assemble(Constant(er*0.5)*Dx(u_,0)**2*dx)
T_dS_int = assemble(2*Constant(n0)*(cosh(u_)-1-u_*sinh(u_))*dx)

U_el_anl = sqrt(er*n0*8)*(np.cosh(u_(0)/2)-1)
min_T_dS_anl = np.sqrt(er*n0*8)*(3-3*np.cosh(u_(0)/2)+u_(0)*np.sinh(u_(0)/2))

print '-TdS: num:', -T_dS_int,' anl: ', min_T_dS_anl
print 'U_el: num:', U_el_int,' anl: ', U_el_anl