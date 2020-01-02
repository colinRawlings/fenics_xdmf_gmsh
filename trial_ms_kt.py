# !python3

###########################################################
# imports
###########################################################

import dolfin as fn
import os

import fenics_utils as fu

import matplotlib.pyplot as plt

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.abspath("geo")

###########################################################
# Main
###########################################################

if __name__ == "__main__":

    # import mesh

    msh_filepath = os.path.join(GEO_DIR, "circle_square_center_pt.geo")
    mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(msh_filepath, do_plots=False)

    def edge_pt(x, on_boundary):
        condition = (abs(x[0]-0.0) < 10*fn.DOLFIN_EPS and abs(x[1]-0.0) < 10*fn.DOLFIN_EPS)
        return condition

    # magnetic equation

    ob_mesh = fn.MeshView.create(mesh_data["boundary_mesh_func"], 1)

    V = fn.FunctionSpace(mesh_data["mesh"], "CG", 2)
    VE = fn.FunctionSpace(mesh_data["mesh"], "CG", 2)
    VL = fn.FunctionSpace(ob_mesh, "CG", 2)
    W = fn.MixedFunctionSpace(V, VE, VL)

    (u, ue, ul) = fn.TrialFunctions(W)
    (v, ve, vl) = fn.TestFunctions(W)

    bcs = [
        fn.DirichletBC(W.sub_space(1), fn.Constant(0), edge_pt, method='pointwise'),
    ]

    M0 = fn.Expression(("0", "1"), degree=2)

    dx_mf = fn.dx(subdomain_data=mesh_data["subdomain_mesh_func"])
    ds_mf = fn.dx(subdomain_data=mesh_data["boundary_mesh_func"])
    dL = fn.Measure("dx", domain=ob_mesh)

    a = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
    a += fn.inner(fn.grad(ue), fn.grad(ve)) * fn.dx
    L = fn.inner(M0, fn.grad(v)) * dx_mf(2)
    a += (u - ue) * vl * dL + ul * ve * dL
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