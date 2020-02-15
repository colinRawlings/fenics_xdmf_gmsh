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
    os.path.join(os.path.dirname(__file__), os.pardir, "shared", "geo"))

invA = 0.01
Hx_app = 0.0

###########################################################
# Functions
###########################################################


def is_center_pt(x, on_boundary):
    condition = (abs(x[0] - 0.0) < 10 * fn.DOLFIN_EPS and abs(x[1] - 0.0) < 10 * fn.DOLFIN_EPS)  # type: ignore
    return condition

Mx_init = fn.Expression('-x[1]/sqrt(pow(x[0],2)+pow(x[1],2)+0.001)*(1-exp(-(x[0]*x[0]+x[1]*x[1])/0.04))', degree=1)
My_init = fn.Expression('x[0]/sqrt(pow(x[0],2)+pow(x[1],2)+0.001)*(1-exp(-(x[0]*x[0]+x[1]*x[1])/0.04))', degree=1)
Mz_init = fn.Expression('(exp(-(x[0]*x[0]+x[1]*x[1])/0.1))', degree=1)

###########################################################
# Main
###########################################################

# Load Mesh

geo_filepath = os.path.join(GEO_DIR, "cylinder_sphere.geo")
labelled_mesh = fu.convert_3d_gmsh_geo_to_fenics_mesh(geo_filepath, {
    "sphere_radius": 2,
    "dx_mesh_inner": 0.05,
    "dx_mesh_outer": 0.4
})

labelled_mesh1 = fu.create_mesh_view(labelled_mesh)
labelled_mesh2 = fu.create_mesh_view(labelled_mesh, 2)

#

plt.figure()
fn.plot(labelled_mesh1.mesh, title="mesh1")

plt.figure()
fn.plot(labelled_mesh2.mesh, title="mesh2")
plt.show()

# Spaces

W_phi = fn.FunctionSpace(labelled_mesh1.mesh, "CG", 1)
W_Mx = fn.FunctionSpace(labelled_mesh2.mesh, "CG", 1)
W_My = fn.FunctionSpace(labelled_mesh2.mesh, "CG", 1)
W_Mz = fn.FunctionSpace(labelled_mesh2.mesh, "CG", 1)
W_l = fn.FunctionSpace(labelled_mesh2.mesh, "CG", 1)

W = fn.MixedFunctionSpace(W_phi, W_Mx, W_My, W_Mz, W_l)

# Define BC

bc1 = fn.DirichletBC(W_phi, fn.Constant(0), fu.IsBoundary())

# Define PDE

(v_phi, v_Mx, v_My, v_Mz, v_l) = fn.TestFunctions(W)  # type: ignore

u = fn.Function(W)
u_phi = u.sub(0)
u_Mx = u.sub(1)
u_My = u.sub(2)
u_Mz = u.sub(3)
u_l = u.sub(4)

x = fn.SpatialCoordinate(labelled_mesh2.mesh)

M = fn.as_vector((u_Mx, u_My, u_Mz))

u_Mx.interpolate(Mx_init)
u_My.interpolate(My_init)
u_Mz.interpolate(Mz_init)
u_l.interpolate(fn.Constant(0.01))

dx1 = fn.Measure("dx", domain=W.sub_space(0).mesh(), subdomain_data=labelled_mesh1.subdomain_mesh_func)
dx2 = fn.Measure("dx", domain=W.sub_space(1).mesh(), subdomain_data=labelled_mesh2.subdomain_mesh_func)

F0 = fn.inner(fn.grad(u_phi), fn.grad(v_phi)) * dx1
F0 += fn.inner(M, fn.grad(v_phi)) * dx1(2)

F1 = fn.Constant(invA) * fn.inner(fn.grad(u_Mx), fn.grad(v_Mx)) * dx2
F1 += fn.Constant(-1) * (fn.Dx(u_phi, 0)+fn.Constant(Hx_app)) * v_Mx * dx2 - u_l * u_Mx * v_Mx * dx2  # type: ignore

F2 = fn.Constant(invA) * fn.inner(fn.grad(u_My), fn.grad(v_My)) * dx2
F2 += fn.Constant(-1) * fn.Dx(u_phi, 1) * v_My * dx2 - u_l * u_My * v_My * dx2  # type: ignore

F3 = fn.Constant(invA) * fn.inner(fn.grad(u_Mz), fn.grad(v_Mz)) * dx2
F3 += fn.Constant(-1) * fn.Dx(u_phi, 2) * v_Mz * dx2 - u_l * u_Mz * v_Mz * dx2  # type: ignore

F4 = +v_l * (u_Mx * u_Mx + u_My * u_My + u_Mz * u_Mz - 1.0) * dx2  # type: ignore

F = F0 + F1 + F2 + F3 + F4

# Compute Solution

fn.solve(F == 0,
            u,
            bc1,
            solver_parameters={
                "nonlinear_solver": "newton",
                "newton_solver": {
                    "linear_solver": "superlu",
                    "relaxation_parameter": 1.0
                }
            })


# Plot Solution

comm = fn.MPI.comm_world

VM = fn.VectorFunctionSpace(labelled_mesh2.mesh, 'CG', 1)
u_Mm = fn.project(M, VM)
file = fn.File(comm, "ms_3d_genlin_M.pvd")
file << u_Mm

file = fn.File(comm, "ms_3d_genlin_Mz.pvd")
file << u_Mz

file = fn.File(comm, "ms_3d_genlin_u_phi.pvd")
file << u_phi