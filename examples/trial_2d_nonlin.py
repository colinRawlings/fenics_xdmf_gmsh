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
# Main
###########################################################

    
# Create meshes

msh_filepath = os.path.join(GEO_DIR, "circle_square_center_pt.geo")
mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(msh_filepath, do_plots=False)
mesh = mesh_data["mesh"]

mesh_data1 = fu.create_mesh_view(mesh_data["mesh"], mesh_data["subdomain_mesh_func"])
mesh_data2 = fu.create_mesh_view(mesh_data["mesh"], mesh_data["subdomain_mesh_func"], 2)

plt.figure()
fn.plot(mesh_data1["mesh"], title="mesh1")

plt.figure()
fn.plot(mesh_data2["mesh"], title="mesh2")
plt.show()

# Create function spaces
W1 = fn.FunctionSpace(mesh_data1["mesh"], "Lagrange", 1)
W2 = fn.FunctionSpace(mesh_data2["mesh"], "Lagrange", 1)

W = fn.MixedFunctionSpace(W1, W2)

# Define boundary conditions
bc1 = fn.DirichletBC(W1, fn.Constant(0), fu.IsBoundary())
bc2 = fn.DirichletBC(W2, fn.Constant(0), fu.IsBoundary())
bcs = [bc1, bc2]

f = fn.Expression("x[0]*sin(x[1])", degree=2)

# Define mixed-domains variational problem
(v1, v2) = fn.TestFunctions(W)  # type: ignore
u = fn.Function(W)
u1 = u.sub(0)
u2 = u.sub(1)

dx1 = fn.Measure("dx", domain=W.sub_space(0).mesh(), subdomain_data=mesh_data1["subdomain_mesh_func"])
dx2 = fn.Measure("dx", domain=W.sub_space(1).mesh(), subdomain_data=mesh_data2["subdomain_mesh_func"])

F1 = fn.inner((fn.Constant(1) + u1 * u1) * fn.grad(u1), fn.grad(v1)) * dx1  # type: ignore
F1 -= fn.Constant(400) * (fn.Constant(1) + u2) * f * v1 * dx1(2)
F2 = fn.inner((fn.Constant(1) + u2 * u2) * fn.grad(u2), fn.grad(v2)) * dx2  # type: ignore
F2 -= fn.Constant(400) * (fn.Constant(1) + u1) * f * v2 * dx2
F = F1 + F2

# Compute solution - mixed-domains problem

fn.solve(F == 0, u, bcs, solver_parameters=fu.newton_solver_parameters())

# plot solution

plt.figure()
p = fn.plot(u1, title="u1")
plt.colorbar(p)

plt.figure()
p = fn.plot(u2, title="u2")
plt.colorbar(p)