# !python3
"""
Print the system matrices for a one dimensional problem:

    dO: u'' - 1 = 0
    lhs: u = 0.5
    rhs: u' = 0

treated using CG1 elements such that it may be easily compared 
with  the expected system matrices.

n.b. the use of the legacy container since one d meshes seem
broken in the the MeshView dev container.
"""

###########################################################
# imports
###########################################################

import os

import dolfin as fn
import matplotlib.pyplot as plt
import numpy as np

import typing as ty

###########################################################
# Definitions
###########################################################

NUM_CELLS = 4
F_VAL = 1
SECTION_DISPLAY_WIDTH = 20

###########################################################
# Functions
###########################################################


def np_array_from_fenics(Af: ty.Any) -> np.array:

    if isinstance(Af, fn.Matrix):
        A = np.zeros((Af.size(0), Af.size(1)))

        for p in range(Af.size(0)):
            sparsity_pattern, values = Af.getrow(p)
            A[p, sparsity_pattern] = values

    elif isinstance(Af, fn.Vector):
        A = Af.get_local()
    else:
        raise ValueError(f"Unexpected type: {type(Af)}")

    return A


###########################################################
# Main
###########################################################

# mesh

mesh = fn.UnitIntervalMesh(NUM_CELLS)

V = fn.FunctionSpace(mesh, "CG", 1)
u = fn.TrialFunction(V)
v = fn.TestFunction(V)

plt.figure()
ax = plt.gca()
fn.plot(mesh, title="mesh")

for dof_n, dof_x in enumerate(V.tabulate_dof_coordinates()):
    ax.annotate(f"dof: {dof_n}", (dof_x, 0))

def is_lhs(x, on_boundary):
    condition = abs(x[0] - 0.0) < 10 * fn.DOLFIN_EPS
    return condition

bcs = [fn.DirichletBC(V, fn.Constant(0.5), is_lhs, method='pointwise')]

a = fn.Constant(-1) * fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
L = fn.Constant(F_VAL) * v * fn.dx

u = fn.Function(V)

fn.solve(a == L, u, bcs)

plt.figure()
ax = plt.gca()
line_u = fn.plot(u, title="solution", label="fenics")
x = np.linspace(0, 1)
line_anl = ax.plot(x, 0.5 * x**2 - x + 0.5, ":k", label="anl")
plt.legend()
ax.set_xlabel("x")
ax.set_ylabel("u")

if "PYTEST_CURRENT_TEST" not in os.environ:
    plt.show()

#

print(f"{'='*SECTION_DISPLAY_WIDTH}\nMatrices\n{'='*SECTION_DISPLAY_WIDTH}")

dx_anl = 1 / NUM_CELLS

print(f"\nBefore BC applied\n{'-'*SECTION_DISPLAY_WIDTH}")

Af = fn.assemble(a)
Lf = fn.assemble(L)

print(f"A:\n{np_array_from_fenics(Af)}")
print(f"anl coeffs: (-1 *) (2 *) 1/dx * 1/dx * dx: (-1 *) (2 *) {1/dx_anl}")

print(f"L:\n{np_array_from_fenics(Lf).reshape(-1, 1)}")
print(f"anl coeffs: (2 *) 0.5*dx*f(x): (2 *) {0.5 * dx_anl * F_VAL}")

print(f"\nAfter BC applied\n{'-'*SECTION_DISPLAY_WIDTH}")

bcs[0].apply(Af)
bcs[0].apply(Lf)

print(f"Af: \n{np_array_from_fenics(Af)}")
print(f"Lf: \n{np_array_from_fenics(Lf).reshape(-1, 1)}")
