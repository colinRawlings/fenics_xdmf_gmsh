# !python3

###########################################################
# imports
###########################################################

import os

import numpy as np
import matplotlib.pyplot as plt

import dolfin as fn

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "shared_geo"))

###########################################################
# Functions
###########################################################


def is_center_pt(x, on_boundary):
    condition = (abs(x[0] - 0.0) < 10 * fn.DOLFIN_EPS and abs(x[1] - 0.0) < 10 * fn.DOLFIN_EPS)  # type: ignore
    return condition


###########################################################
# Main
###########################################################

# import mesh

geo_filepath = os.path.join(GEO_DIR, "circle_square_center_pt.geo")
labelled_mesh = fu.convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath, do_plots=False, geo_params={"radius": 1})

# magnetic equation

ob_mesh = fn.MeshView.create(labelled_mesh.boundary_mesh_func, 1)

V = fn.FunctionSpace(labelled_mesh.mesh, "CG", 2)
VE = fn.FunctionSpace(labelled_mesh.mesh, "CG", 2)
VL = fn.FunctionSpace(ob_mesh, "CG", 2)
W = fn.MixedFunctionSpace(V, VE, VL)

(u, ue, ul) = fn.TrialFunctions(W)  # type: ignore
(v, ve, vl) = fn.TestFunctions(W)  # type: ignore

bcs = [
    fn.DirichletBC(W.sub_space(1), fn.Constant(0), is_center_pt, method='pointwise'),
]

M0 = fn.Expression(("0", "1"), degree=2)

dx_mf = fn.dx(subdomain_data=labelled_mesh.subdomain_mesh_func)
ds_mf = fn.dx(subdomain_data=labelled_mesh.boundary_mesh_func)
dL = fn.Measure("dx", domain=ob_mesh)

a_u = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx + ul * ve * dL
L_u = fn.inner(M0, fn.grad(v)) * dx_mf(2)
a_ue = fn.inner(fn.grad(ue), fn.grad(ve)) * fn.dx - ul * v * dL
a_vl = (u - ue) * vl * dL
L_vl = fn.Constant(0) * vl * dL

#

w = fn.Function(W)

fn.solve(a_u + a_ue + a_vl == L_u + L_vl, w, bcs, solver_parameters={'linear_solver': 'lu'})

u, ue = w.sub(0), w.sub(1)

#

plt.figure()
c = fn.plot(w.sub(0), title="u")
plt.colorbar(c)

plt.figure()
c = fn.plot(w.sub(1), title="ue")
plt.colorbar(c)

yv = np.linspace(-1, 1)
uv = [u(0, yvv) for yvv in yv]
uev = [ue(0, yvv) for yvv in yv]

plt.figure()
plt.plot(yv, uv)

plt.figure()
plt.plot(yv, uev)

plt.show()

#

# TODO: cf with analytical solution for uniformly magnetised rod!!!
# expected_ue_pt = 0.040048882089143374
# assert abs(ue(0, 0.25) - expected_ue_pt) / expected_ue_pt < 1e-4, "Solution seems strange"

#

print("Done")