# !python

###########################################################
# Imports
###########################################################

import matplotlib.pyplot as plt

import dolfin as fn

import fenics_utils as fu

###########################################################
# Main
###########################################################

mesh = fn.UnitSquareMesh(10, 10)

V = fn.FunctionSpace(mesh, "CG", 2)

u = fn.Function(V)
v = fn.TrialFunction(V)

bc = fn.DirichletBC(V, fn.Constant(0), fu.IsBoundary())

F = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
F += u**3 * v * fn.dx + fn.Constant(100) * v * fn.dx  # type: ignore

fn.solve(F == 0,
         u,
         bc,
         solver_parameters={
             "nonlinear_solver": "newton",
             "newton_solver": {
                 "linear_solver": "gmres",
                 "preconditioner": "ilu"
             }
         })

plt.figure()
p = fn.plot(u, title="u")
plt.colorbar(p)
