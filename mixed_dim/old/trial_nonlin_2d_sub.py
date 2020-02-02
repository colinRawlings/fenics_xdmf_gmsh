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

    sq_mesh = fn.MeshView.create(mesh_data["subdomain_mesh_func"], 2)

    plt.figure()
    fn.plot(sq_mesh)
    plt.show()

    # Solve a dummy problem

    el = fn.FiniteElement("CG", mesh.ufl_cell(), 2)  # probably one needs to set dim=3?
    eli = fn.FiniteElement("CG", sq_mesh.ufl_cell(), 2)

    V = fn.FunctionSpace(mesh, el)
    Vi = fn.FunctionSpace(sq_mesh, eli)

    W = fn.MixedFunctionSpace(V, Vi)

    w = fn.Function(W)
    u = w.sub(1)
    ui = w.sub(0)
    v, vi = fn.TestFunctions(W)

    class DirichletBoundary2(fn.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    bcs = [
        fn.DirichletBC(V, fn.Constant(0), mesh_data["boundary_mesh_func"], 1),
        fn.DirichletBC(Vi, fn.Constant(0), DirichletBoundary2())
    ]

    dx_mf = fn.dx(subdomain_data=mesh_data["subdomain_mesh_func"])

    dx1 = fn.Measure("dx", domain=W.sub_space(0).mesh())
    dxi = fn.Measure("dx", domain=W.sub_space(1).mesh())

    F = fn.inner(fn.grad(u), fn.grad(v)) * dx1
    F -= fn.Constant(1) * v * dx1

    F += fn.inner(fn.grad(ui), fn.grad(vi)) * dxi
    F -= 100 * u * vi * dxi

    # F += (ui + 500*ui**3) * vi * dx_mf(2)

    #


    def newton_solver_parameters():
        return {"nonlinear_solver": "newton", "newton_solver": {"linear_solver": "gmres"}}

    fn.solve(F == 0, w, bcs, solver_parameters=newton_solver_parameters())

    # Plot solution
    plt.figure()
    c = fn.plot(w.sub(0), title="u")
    plt.colorbar(c)

    plt.figure()
    c = fn.plot(w.sub(1), title="ui")
    plt.colorbar(c)
