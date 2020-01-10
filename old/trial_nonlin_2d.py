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
    mesh = mesh_data["mesh"]

    # Solve a dummy problem

    el1 = fn.FiniteElement("CG", mesh.ufl_cell(), 2)  # probably one needs to set dim=3?
    el2 = fn.FiniteElement("CG", mesh.ufl_cell(), 2)
    W = fn.FunctionSpace(mesh, fn.MixedElement(el1, el2))

    w = fn.Function(W)
    u, ui = fn.split(w)
    v, vi = fn.TestFunctions(W)

    bcs = [
        fn.DirichletBC(W.sub(0), fn.Constant(0), mesh_data["boundary_mesh_func"], 1),
        fn.DirichletBC(W.sub(1), fn.Constant(0), mesh_data["boundary_mesh_func"], 1)
    ]

    dx_mf = fn.dx(subdomain_data=mesh_data["subdomain_mesh_func"])

    F = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
    F -= fn.Constant(1) * v * dx_mf(2)

    F += fn.inner(fn.grad(ui), fn.grad(vi)) * fn.dx
    F -= 100*u * vi * dx_mf(2)
    F += (ui + 500*ui**3) * vi * dx_mf(2)

    #

    fn.solve(F == 0, w, bcs)

    # Plot solution
    plt.figure()
    c = fn.plot(w.sub(0), title="u")
    plt.colorbar(c)

    plt.figure()
    c = fn.plot(w.sub(1), title="ui")
    plt.colorbar(c)
