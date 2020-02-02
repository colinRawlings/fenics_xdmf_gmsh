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

GEO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "shared_geo"))

###########################################################
# Main
###########################################################

# Load Mesh

geo_filepath = os.path.join(GEO_DIR, "circle_square_center_pt.geo")
labelled_mesh = fu.convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath, do_plots=False, geo_params={"radius": 1})

labelled_mesh1 = fu.create_mesh_view(labelled_mesh)
labelled_mesh2 = fu.create_mesh_view(labelled_mesh, 2)

x1 = fn.SpatialCoordinate(labelled_mesh1.mesh)
x2 = fn.SpatialCoordinate(labelled_mesh2.mesh)

#

plt.figure()
fn.plot(labelled_mesh1.mesh, title="mesh1")

plt.figure()
fn.plot(labelled_mesh2.mesh, title="mesh2")
plt.show()

# Spaces

W1 = fn.FunctionSpace(labelled_mesh1.mesh, "CG", 2)
W2 = fn.FunctionSpace(labelled_mesh2.mesh, "CG", 2)

W = fn.MixedFunctionSpace(W1, W2)

# Define BC

bc1 = fn.DirichletBC(W1, fn.Constant(0), fu.IsBoundary())
bc2 = fn.DirichletBC(W2, fn.Constant(0), fu.IsBoundary())
bcs = [bc1, bc2]

# Define PDE

(v1, v2) = fn.TestFunctions(W)  # type: ignore
u = fn.Function(W)
u1 = u.sub(0)
u2 = u.sub(1)

dx1 = fn.Measure("dx",
                 domain=W.sub_space(0).mesh(),
                 subdomain_data=labelled_mesh1.subdomain_mesh_func)
dx2 = fn.Measure("dx",
                 domain=W.sub_space(1).mesh(),
                 subdomain_data=labelled_mesh2.subdomain_mesh_func)

f1 = x1[0] * fn.sin(x1[1])  # type: ignore
f2 = x2[0] * fn.sin(x2[1])  # type: ignore

F1 = fn.inner((fn.Constant(1) + u1 * u1) * fn.grad(u1), fn.grad(v1)) * dx1  # type: ignore
F1 -= fn.Constant(400) * (fn.Constant(1) + u2) * f1 * v1 * dx1(2)
F2 = fn.inner((fn.Constant(1) + u2 * u2) * fn.grad(u2), fn.grad(v2)) * dx2  # type: ignore
F2 -= fn.Constant(400) * (fn.Constant(1) + u1) * f2 * v2 * dx2
F = F1 + F2

# Compute Solution

fn.solve(F == 0, u, bcs, solver_parameters=fu.newton_solver_parameters())

# Plot Solution

plt.figure()
p = fn.plot(u1, title="u1")
plt.colorbar(p)

plt.figure()
p = fn.plot(u2, title="u2")
plt.colorbar(p)
