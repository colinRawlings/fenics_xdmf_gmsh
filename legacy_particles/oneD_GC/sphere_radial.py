# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 13:22:33 2015

@author: cra
"""

from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

# solve sphere problem in 1D to relate surface potential 
# to surface charge

n0 = 3.
er = 2.

A1 = 5.
I1 = 7.8747 # scaled surface charge

R_balloon = 3. # 
R_mesh = 0.05 # mesh size as a function of debye length

k0 = sqrt(n0/er)

a = A1/(sqrt(2)*k0)
sigma = er*k0*sqrt(2)*I1

mesh = IntervalMesh(int(np.ceil(R_balloon/R_mesh)),a,a+R_balloon/k0)
V = FunctionSpace(mesh,'CG',2)

#-- setup for BC
bnd_mf = MeshFunction('size_t',mesh,0)
class initMF(SubDomain):
    def inside(self,x,on_boundary):
        return True
        
class LeftB(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and np.abs(x[0]-a) < 1e2*DOLFIN_EPS

init_mf = initMF()
lb = LeftB()
init_mf.mark(bnd_mf,0)
lb.mark(bnd_mf,1)
ds = Measure('ds')[bnd_mf]

#-- specify problem
u_ = Function(V)
u = TrialFunction(V)
v = TestFunction(V)

r = SpatialCoordinate(mesh)

F = r**2/er*Constant(sigma)*v*ds(1)-\
    r**2*Dx(u,0)*Dx(v,0)*dx-\
    2*n0*r**2/er*sinh(u)*v*dx

F = action(F,u_)
J = derivative(F,u_,u)

problem = NonlinearVariationalProblem(F,u_,[],J)
solver = NonlinearVariationalSolver(problem)
solver.solve()

plot(u_)
print 'surface u: ',u_(a)