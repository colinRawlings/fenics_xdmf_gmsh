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
sigma = 4.

k = sqrt(n0/er)

R_domain = 10. 

mesh = IntervalMesh(200,0,R_domain/k)
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

plot(u_)