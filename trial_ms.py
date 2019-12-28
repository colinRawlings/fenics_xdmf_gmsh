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

    msh_filepath = os.path.join(GEO_DIR, "square_in_circle.geo")
    mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(msh_filepath, do_plots=False)

    # magnetic equation

    V = fn.FunctionSpace(mesh_data["mesh"], "CG", 2)
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    bcs = [fn.DirichletBC(V, fn.Constant(0), mesh_data["boundary_mesh_func"], 1)]

    M0 = fn.Expression(("0", "1"), degree=2)

    dx_mf = fn.dx(subdomain_data=mesh_data["subdomain_mesh_func"])
    ds_mf = fn.dx(subdomain_data=mesh_data["boundary_mesh_func"])

    n = fn.FacetNormal(mesh_data["mesh"])

    F = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
    F -= fn.inner(M0, fn.grad(v)) * dx_mf(2)

    #

    u = fn.Function(V)
    u = fu.solve_weak_form(u, F, bcs)

    plt.figure()
    c = fn.plot(u)
    plt.colorbar(c)
