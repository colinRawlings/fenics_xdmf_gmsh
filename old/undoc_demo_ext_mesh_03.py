# !python3

###########################################################
# Imports
###########################################################

import os

import matplotlib.pyplot as plt

import dolfin as fn

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.abspath("geo")

###########################################################
# Functions/Classes
###########################################################


# Sub domains for Dirichlet boundary conditions
class DirichletBoundary1(fn.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


class DirichletBoundary2(fn.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


def newton_solver_parameters():
    return {"nonlinear_solver": "newton", "newton_solver": {"linear_solver": "gmres"}}


###########################################################
# Main
###########################################################

# Create meshes

msh_filepath = os.path.join(GEO_DIR, "circle_square_center_pt.geo")
mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(msh_filepath, do_plots=False)
mesh = mesh_data["mesh"]

submesh1 = fn.MeshView.create(mesh_data["subdomain_mesh_func"], 1)
submesh2 = fn.MeshView.create(mesh_data["subdomain_mesh_func"], 2)

mesh1 = mesh
mesh2 = submesh2

plt.figure()
fn.plot(mesh1, title="mesh2")

plt.figure()
fn.plot(mesh2, title="mesh2")
plt.show()

# Create function spaces
W1 = fn.FunctionSpace(mesh, "Lagrange", 1)
W2 = fn.FunctionSpace(submesh2, "Lagrange", 1)
# Define the mixed function space W = W1 x W2
W = fn.MixedFunctionSpace(W1, W2)

# Define boundary conditions
g = fn.Constant(1.0)
bc1 = fn.DirichletBC(W1, g, DirichletBoundary1())
bc2 = fn.DirichletBC(W2, g, DirichletBoundary2())
bcs = [bc1, bc2]

f = fn.Expression("100*x[0]*sin(x[1])", degree=2)

# Define mixed-domains variational problem
(v1, v2) = fn.TestFunctions(W)
u = fn.Function(W)
u1 = u.sub(0)
u2 = u.sub(1)

dx1 = fn.Measure("dx",
                 domain=W.sub_space(0).mesh(),
                 subdomain_data=mesh_data["subdomain_mesh_func"])
dx2 = fn.Measure("dx", domain=W.sub_space(1).mesh())

F1 = fn.inner((1 + u1**2) * fn.grad(u1), fn.grad(v1)) * dx1 - f * v1 * dx1(2)
F2 = fn.inner((1 + u2**2) * fn.grad(u2), fn.grad(v2)) * dx2 - 100 * u1 * v2 * dx2
F = F1 + F2

# Compute solution - mixed-domains problem
fn.solve(F == 0, u, bcs, solver_parameters=newton_solver_parameters())

plt.figure()
p = fn.plot(u1, title="u1")
plt.colorbar(p)

plt.figure()
p = fn.plot(u2, title="u2")
plt.colorbar(p)