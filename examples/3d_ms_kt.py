# !python3

###########################################################
# imports
###########################################################

import dolfin as fn
import os

import fenics_utils as fu
import numpy as np

import matplotlib.pyplot as plt

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "geo"))
RESULTS_DIR = fu.get_clean_results_dir(__file__)

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

# import mesh

msh_filepath = os.path.join(GEO_DIR, "sph_sph_center_pt.geo")
labelled_mesh = fu.convert_3d_gmsh_geo_to_fenics_mesh(msh_filepath)

# magnetic equation

ob_mesh = fn.MeshView.create(labelled_mesh.boundary_mesh_func, 1)

V = fn.FunctionSpace(labelled_mesh.mesh, "CG", 2)
VE = fn.FunctionSpace(labelled_mesh.mesh, "CG", 2)
VL = fn.FunctionSpace(ob_mesh, "CG", 2)
W = fn.MixedFunctionSpace(V, VE, VL)

(u, ue, ul) = fn.TrialFunctions(W)  # type: ignore
(v, ve, vl) = fn.TestFunctions(W)  # type: ignore

bcs = [
    fn.DirichletBC(W.sub_space(1),
                   fn.Constant(0),
                   is_center_pt,
                   method='pointwise'),
]

M0 = fn.Expression(("0", "1", "0"), degree=2)

dx_mf = fn.dx(subdomain_data=labelled_mesh.subdomain_mesh_func)
ds_mf = fn.dx(subdomain_data=labelled_mesh.boundary_mesh_func)
dL = fn.Measure("dx", domain=ob_mesh)

a = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
a += fn.inner(fn.grad(ue), fn.grad(ve)) * fn.dx
L = fn.inner(M0, fn.grad(v)) * dx_mf(2)
a += (u - ue) * vl * dL + ul * (v - ve) * dL
L += fn.Constant(0) * vl * dL

#

w = fn.Function(W)
fn.solve(a == L, w, [], solver_parameters={"linear_solver": "lu"})

#

file_u = fn.File(os.path.join(RESULTS_DIR, "u.pvd"))
file_u << w.sub(0)

file_ue = fn.File(os.path.join(RESULTS_DIR, "ue.pvd"))
file_ue << w.sub(1)

#

yv = np.linspace(-1, 1)
uv = [w.sub(0)(0, yvv, 0) for yvv in yv]
uev = [w.sub(1)(0, yvv, 0) for yvv in yv]

plt.figure()
plt.plot(yv, uv)

plt.figure()
plt.plot(yv, uev)

print("done")
