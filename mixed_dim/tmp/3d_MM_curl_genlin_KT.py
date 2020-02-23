# !python3

###########################################################
# Imports
###########################################################

import os
import time

import matplotlib.pyplot as plt

import dolfin as fn
import ufl

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, "shared", "geo"))

RESULTS_DIR = fu.get_clean_results_dir(__file__)

A = 0.02
Hx_app = 0.0
core_radius_guess = 0.2
element_order = 1
cylinder_height = 0.5

def is_cylinder_pt(x, on_boundary):
    condition = (abs(x[0] - 0.5) < 10 * fn.DOLFIN_EPS
                 and abs(x[1] - 0.0) < 10 * fn.DOLFIN_EPS  # type: ignore
                 and abs(x[2] - 0.0 * cylinder_height) < 10 * fn.DOLFIN_EPS)  # type: ignore

    if condition==True:
        print("!!!!!!!!!!!!!! POINT FOUND !!!!!!!!!!!")

    return condition

###########################################################
# Main
###########################################################

# Load Mesh

geo_filepath = os.path.join(GEO_DIR, "cylinder_sphere.geo")
labelled_mesh = fu.convert_3d_gmsh_geo_to_fenics_mesh(geo_filepath, {
    "sphere_radius": 1.0,
    "cylinder_height": 0.5,
    "dx_mesh_inner": 0.1,
    "dx_mesh_outer": 0.4
})

mesh_phi = fu.create_mesh_view(labelled_mesh)
mesh_M = fu.create_mesh_view(labelled_mesh, 2)
mesh_kt = fn.MeshView.create(labelled_mesh.boundary_mesh_func, 1)


# Spaces

W_phi = fn.FunctionSpace(mesh_phi.mesh, "CG", 1)
W_M = fn.VectorFunctionSpace(mesh_M.mesh, "CG", element_order)
W_l = fn.FunctionSpace(mesh_M.mesh, "CG", element_order)
W_kt = fn.FunctionSpace(mesh_kt, "CG", 1)
W_phiE = fn.FunctionSpace(mesh_phi.mesh, "CG", 1)

W = fn.MixedFunctionSpace(W_phi, W_M, W_l, W_kt, W_phiE)

# Define BC

# bc1 = fn.DirichletBC(W_phi, fn.Constant(0), fu.IsBoundary())
# bc1 = []

bc1 = fn.DirichletBC(W.sub_space(4),
                   fn.Constant(0),
                   is_cylinder_pt,
                   method='pointwise'),

# Define terms

(v_phi, v_M, v_l, v_Lkt, v_phiE) = fn.TestFunctions(W)

u = fn.Function(W)
u_phi = u.sub(0)
u_M = u.sub(1)
u_l = u.sub(2)
u_Lkt = u.sub(3)  # the lagrange multiplier for pinning phi internal to phi external
u_phiE = u.sub(4)

i, j = ufl.indices(2)
x = fn.SpatialCoordinate(mesh_M.mesh)
dx_phi = fn.Measure("dx",
                    domain=mesh_phi.mesh,
                    subdomain_data=mesh_phi.subdomain_mesh_func)
dx_M = fn.Measure("dx",
                  domain=mesh_M.mesh,
                  subdomain_data=mesh_M.subdomain_mesh_func)
dL = fn.Measure("dx", domain=mesh_kt)


H_app = fn.as_vector((fn.Constant(Hx_app), fn.Constant(0), fn.Constant(0)))

# Initial condition (curling state)

u_M_init = fn.TrialFunction(W_M)
v_M_init = fn.TestFunction(W_M)

rho = fn.sqrt(x[0] * x[0] + x[1] * x[1])
Mz_guess = fn.exp(-(rho * rho) / fn.Constant(core_radius_guess**2))
xy_rescaling = fn.sqrt(1 - Mz_guess * Mz_guess)
Mx_guess = -x[1] / (rho + 0.001) * xy_rescaling
My_guess = x[0] / (rho + 0.001) * xy_rescaling

M_init = fn.as_vector((Mx_guess, My_guess, Mz_guess))

a_init = u_M_init[i] * v_M_init[i] * fn.dx
L_init = M_init[i] * v_M_init[i] * fn.dx

fn.solve(a_init == L_init, u_M, [])

# Define PDE

F__phi_M = fn.Dx(u_phi, i) * fn.Dx(v_phi, i) * dx_phi
F__phi_M += u_M[i] * fn.Dx(v_phi, i) * dx_phi(2)

F__M_phi = fn.Constant(A) * fn.Dx(u_M[i], j) * fn.Dx(v_M[i], j) * dx_M
F__M_phi += fn.Constant(-1) * (fn.Dx(u_phi, i) + H_app[i]) * v_M[i] * dx_M

F__M_l = fn.Constant(-1) * u_l * u_M[i] * v_M[i] * dx_M
F__M_l += (u_M[i] * u_M[i] - 1) * v_l * dx_M

F__phiE_Lkt = fn.inner(fn.grad(u_phiE), fn.grad(v_phiE))*fn.dx
F__phiE_Lkt += (u_phi - u_phiE) * v_Lkt * dL + u_Lkt * (v_phi - v_phiE) * dL


F = F__phi_M + F__M_phi + F__M_l + F__phiE_Lkt

# Compute Solution

start_time = time.time()
try:
    fn.solve(F == 0,
            u,
            bc1,
            solver_parameters={
                "nonlinear_solver": "newton",
                "newton_solver": {
                    "linear_solver": "lu",
                    "relaxation_parameter": 1.0,
                    "maximum_iterations": 10,
                }
            })
except Exception as e:
    print(repr(e))


print(f"Solving took: {time.time() - start_time}s")

# Export Solution

fu.save_function(u_phi, os.path.join(RESULTS_DIR, "u_phi.pvd"))
fu.save_function(u_phiE, os.path.join(RESULTS_DIR, "u_phiE.pvd"))
fu.save_function(u_M.sub(2), os.path.join(RESULTS_DIR, "Mz.pvd"))
fu.save_function(u_M, os.path.join(RESULTS_DIR, "M.pvd"))
