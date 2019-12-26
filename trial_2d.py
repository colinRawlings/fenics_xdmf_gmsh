# !python3

###########################################################
# imports
###########################################################

import dolfin as fn
import os
import pkg_resources as rc

import fenics_utils as fu

import matplotlib.pyplot as plt

###########################################################
# Definitions
###########################################################

MSH_DIR = rc.resource_filename("fenics_utils", "assets")

###########################################################
# Main
###########################################################

if __name__ == "__main__":

    # import mesh

    msh_filepath = os.path.join(MSH_DIR, "simple_poisson_2d.msh")
    mesh_data = fu.convert_2d_gmsh_msh_to_fenics_msh(msh_filepath, do_plots=True)

    # Solve a dummy problem

    V = fn.FunctionSpace(mesh_data["mesh"], "CG", 2)
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    bcs = [fn.DirichletBC(V, fn.Constant(0), mesh_data["boundary_mesh_func"], 2)]

    dx_mf = fn.dx(subdomain_data=mesh_data["subdomain_mesh_func"])

    F = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
    F += fn.Constant(10) * u * v * fn.dx
    F -= fn.Constant(1) * v * dx_mf(2)

    #

    u = fn.Function(V)
    u = fu.solve_weak_form(u, F, bcs)

    # Plot solution(s)
    plt.figure()
    c = fn.plot(u, title="Solution")
    plt.colorbar(c)
