# !python3

###########################################################
# imports
###########################################################

import dolfin as fn
import os

import numpy as np
import matplotlib.pyplot as plt

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "geo"))

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

msh_filepath = os.path.join(GEO_DIR, "circle_square_center_pt.geo")
mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(msh_filepath, do_plots=False)

# magnetic equation

ob_mesh = fn.MeshView.create(mesh_data["boundary_mesh_func"], 1)

V = fn.FunctionSpace(mesh_data["mesh"], "CG", 2)
VE = fn.FunctionSpace(mesh_data["mesh"], "CG", 2)
VL = fn.FunctionSpace(ob_mesh, "CG", 2)
W = fn.MixedFunctionSpace(V, VE, VL)

(u, ue, ul) = fn.TrialFunctions(W)  # type: ignore
(v, ve, vl) = fn.TestFunctions(W)  # type: ignore

bcs = [
    fn.DirichletBC(W.sub_space(1), fn.Constant(0), is_center_pt, method='pointwise'),
]

M0 = fn.Expression(("0", "1"), degree=2)

dx_mf = fn.dx(subdomain_data=mesh_data["subdomain_mesh_func"])
ds_mf = fn.dx(subdomain_data=mesh_data["boundary_mesh_func"])
dL = fn.Measure("dx", domain=ob_mesh)

a = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
a += fn.inner(fn.grad(ue), fn.grad(ve)) * fn.dx
L = fn.inner(M0, fn.grad(v)) * dx_mf(2)
a += (u - ue) * vl * dL + ul * (ve - v) * dL
L += fn.Constant(0) * vl * dL

#

w = fn.Function(W)
fn.solve(a == L, w, bcs)

#

plt.figure()
c = fn.plot(w.sub(0), title="u")
plt.colorbar(c)

plt.figure()
c = fn.plot(w.sub(1), title="ue")
plt.colorbar(c)

yv = np.linspace(-1, 1)
uv = [w.sub(0)(0, yvv) for yvv in yv]
uev = [w.sub(1)(0, yvv) for yvv in yv]

plt.figure()
plt.plot(yv, uv)

plt.figure()
plt.plot(yv, uev)