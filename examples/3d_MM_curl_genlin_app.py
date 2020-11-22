# !python3

###########################################################
# Imports
###########################################################

import os

import numpy as np
import dolfin as fn
import ufl
import tqdm

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "geo"))

RESULTS_DIR = fu.get_clean_results_dir(__file__)

A = 0.02
Hx_app = 0.01
core_radius_guess = 0.2
element_order = 2

###########################################################
# Main
###########################################################

# Create Mesh

geo_filepath = os.path.join(GEO_DIR, "cylinder_sphere.geo")
labelled_mesh = fu.convert_3d_gmsh_geo_to_fenics_mesh(
    geo_filepath,
    {
        "sphere_radius": 2,
        "dx_mesh_inner": 0.10,
        "dx_mesh_outer": 0.4,
        "cylinder_height": 0.5,
    },
)

mesh_phi = fu.create_mesh_view(labelled_mesh)
mesh_M = fu.create_mesh_view(labelled_mesh, 2)

# Create Function Spaces

W_phi = fn.FunctionSpace(mesh_phi.mesh, "CG", 1)
W_M = fn.VectorFunctionSpace(mesh_M.mesh, "CG", element_order)
W_l = fn.FunctionSpace(mesh_M.mesh, "CG", element_order)

W = fn.MixedFunctionSpace(W_phi, W_M, W_l)

# Define BC

bc1 = fn.DirichletBC(W_phi, fn.Constant(0), fu.IsBoundary())

# Define terms

u = fn.Function(W)
(v_phi, v_M, v_l) = fn.TestFunctions(W)  # type: ignore

u_phi, u_M, u_l = u.sub(0), u.sub(1), u.sub(2)

i, j = ufl.indices(2)
x = fn.SpatialCoordinate(mesh_M.mesh)
dx_phi = fn.Measure(
    "dx", domain=mesh_phi.mesh, subdomain_data=mesh_phi.subdomain_mesh_func
)
dx_M = fn.Measure("dx", domain=mesh_M.mesh, subdomain_data=mesh_M.subdomain_mesh_func)


# Setup Initial condition (curling state)

u_M_init = fn.TrialFunction(W_M)
v_M_init = fn.TestFunction(W_M)

rho = fn.sqrt(x[0] * x[0] + x[1] * x[1])  # type: ignore
Mz_guess = fn.exp(-(rho * rho) / fn.Constant(core_radius_guess ** 2))
xy_rescaling = fn.sqrt(1 - Mz_guess * Mz_guess)
Mx_guess = -x[1] / (rho + 0.001) * xy_rescaling  # type: ignore
My_guess = x[0] / (rho + 0.001) * xy_rescaling  # type: ignore

M_init = fn.as_vector((Mx_guess, My_guess, Mz_guess))

a_init = u_M_init[i] * v_M_init[i] * fn.dx  # type: ignore
L_init = M_init[i] * v_M_init[i] * fn.dx  # type: ignore

fn.solve(a_init == L_init, u_M, [])

# Compute Solutions

Hx_app_all = np.linspace(0, 0.15, 100)

for index, Hx_app in enumerate(tqdm.tqdm(Hx_app_all)):

    # Define PDE

    H_app = fn.as_vector((fn.Constant(Hx_app), fn.Constant(0), fn.Constant(0)))

    F__phi_M = fn.Dx(u_phi, i) * fn.Dx(v_phi, i) * dx_phi
    F__phi_M += u_M[i] * fn.Dx(v_phi, i) * dx_phi(2)  # type: ignore

    F__M_phi = (
        fn.Constant(A) * fn.Dx(u_M[i], j) * fn.Dx(v_M[i], j) * dx_M
    )  # type: ignore
    F__M_phi += (
        fn.Constant(-1) * (fn.Dx(u_phi, i) + H_app[i]) * v_M[i] * dx_M
    )  # type: ignore

    F__M_l = fn.Constant(-1) * u_l * u_M[i] * v_M[i] * dx_M  # type: ignore
    F__M_l += (u_M[i] * u_M[i] - 1) * v_l * dx_M  # type: ignore

    F = F__phi_M + F__M_phi + F__M_l

    # Solve

    fn.solve(
        F == 0,
        u,
        bc1,
        solver_parameters={
            "nonlinear_solver": "newton",
            "newton_solver": {
                "linear_solver": "superlu",
                "maximum_iterations": 50,
                "relaxation_parameter": 1.0,
            },
        },
    )

    # Export Solution

    H = fu.calculate_gradient(u_phi)

    fu.save_function(H, os.path.join(RESULTS_DIR, f"H_{index}.pvd"))
    fu.save_function(u_phi, os.path.join(RESULTS_DIR, f"phi_{index}.pvd"))
    fu.save_function(u_M.sub(2), os.path.join(RESULTS_DIR, f"Mz_{index}.pvd"))
    fu.save_function(u_M, os.path.join(RESULTS_DIR, f"M_{index}.pvd"))
