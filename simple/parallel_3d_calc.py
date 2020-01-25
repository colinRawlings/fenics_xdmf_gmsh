# !python3
"""
Here:

- use the stable container since parallel execution
seems broken in the MeshView dev container.
- use fenicstools to evaluate solutions for parallel executions

run as:

mpirun -np 4 /usr/bin/python3 /workspaces/fenics_xdmf_gmsh/simple/parallel_3d_calc.py

"""

###########################################################
# imports
###########################################################

import dolfin as fn

import fenicstools as ft
import numpy as np

###########################################################
# Classes
###########################################################


class IsBoundary(fn.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


###########################################################
# Main
###########################################################

N = 5

mesh = fn.UnitCubeMesh(N, N, N)

V = fn.FunctionSpace(mesh, 'CG', 2)

u = fn.TrialFunction(V)
v = fn.TestFunction(V)

a = fn.Constant(-1) * fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
L = fn.Constant(1) * v * fn.dx

# construct for anl solution
x = fn.SpatialCoordinate(mesh)
anl_soln = fn.Constant(1 / 6) * (x[0]**2 + x[1]**2 + x[2]**2)  # type: ignore
bc = fn.DirichletBC(V, anl_soln, IsBoundary())

#

print(
    f"mpi process: rank: {fn.MPI.rank(fn.MPI.comm_world)}, size: {fn.MPI.size(fn.MPI.comm_world)}")

u = fn.Function(V)
fn.solve(a == L, u, bc)

# post pro

u_anl = fn.project(anl_soln, V)

x = np.array([[0.5, 0.5, 0.5]])

probes_u = ft.Probes(x.flatten(), V)
probes_u(u)
u_fem = probes_u.array()

probes_u_expected = ft.Probes(x.flatten(), V)
probes_u_expected(u_anl)
u_expected = probes_u_expected.array()

if fn.MPI.rank(fn.MPI.comm_world) == 0:  # is root process
    assert abs(u_fem - u_expected) < 1e-5 * abs(u_expected)  # type: ignore
    print("Success!")
