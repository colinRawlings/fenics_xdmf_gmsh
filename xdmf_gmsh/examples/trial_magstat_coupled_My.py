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

GEO_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, "geo"))

###########################################################
# Functions
###########################################################


def is_center_pt(x, on_boundary):
    condition = (abs(x[0] - 0.0) < 10 * fn.DOLFIN_EPS
                 and abs(x[1] - 0.0) < 10 * fn.DOLFIN_EPS)  # type: ignore
    return condition


###########################################################
# Main
###########################################################

# Load Mesh

geo_filepath = os.path.join(GEO_DIR, "circle_square_center_pt.geo")
labelled_mesh = fu.convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath, {
    "radius": 2,
    "dx_inner_mesh": 0.1,
    "dx_outer_mesh": 0.4
},
                                                      do_plots=False)

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

bcs = [bc1]

# Define PDE

(v_phi, v_My) = fn.TestFunctions(W)  # type: ignore

u = fn.Function(W)
u_phi = u.sub(0)
u_My = u.sub(1)

dx1 = fn.Measure("dx",
                 domain=W.sub_space(0).mesh(),
                 subdomain_data=labelled_mesh1.subdomain_mesh_func)
dx2 = fn.Measure("dx",
                 domain=W.sub_space(1).mesh(),
                 subdomain_data=labelled_mesh2.subdomain_mesh_func)

M = fn.as_vector((fn.Constant(1), u_My))


F1 = fn.inner(fn.grad(u_phi), fn.grad(v_phi)) * dx1
F1 += fn.inner(M, fn.grad(v_phi)) * dx1(2)

F2 = fn.Constant(0.1)*fn.inner(fn.grad(u_My), fn.grad(v_My)) * dx2
F2 += fn.Constant(-1) * fn.Dx(u_phi, 1) * v_My * dx2

F = F1 + F2

# Compute Solution

fn.solve(F == 0, u, bcs, solver_parameters=fu.newton_solver_parameters())

# Plot Solution

plt.figure()
p = fn.plot(u_phi, title="u_phi")
plt.colorbar(p)

plt.figure()
p = fn.plot(u_My, title="u_My")
plt.colorbar(p)

fn.plot(M)