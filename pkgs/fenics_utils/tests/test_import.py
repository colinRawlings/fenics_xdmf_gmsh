# !python3

###########################################################
# imports
###########################################################

import os

import pkg_resources as rc

import dolfin as fn

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = rc.resource_filename("fenics_utils", "assets")  # type: ignore

###########################################################
# Fixures
###########################################################


###########################################################
# Test
###########################################################

def test_2d_geo_import_and_solve():
    geo_filepath = os.path.join(GEO_DIR, "simple_poisson_2d.geo")
    labelled_mesh = fu.convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath)

    # Solve a dummy problem

    V = fn.FunctionSpace(labelled_mesh.mesh, "CG", 2)
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    bcs = [fn.DirichletBC(V, fn.Constant(0), labelled_mesh.boundary_mesh_func, 2)]

    dx_mf = fn.dx(subdomain_data=labelled_mesh.subdomain_mesh_func)

    F = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
    F += fn.Constant(10) * u * v * fn.dx
    F -= fn.Constant(1) * v * dx_mf(2)

    #

    u = fn.Function(V)
    u = fu.solve_weak_form(u, F, bcs)

    # 

    assert round(0.0175472, 5) == round(u(0.8, 0.8), 5)

def test_3d_geo_import_and_solve():
    geo_filepath = os.path.join(GEO_DIR, "simple_poisson_3d.geo")
    labelled_mesh = fu.convert_3d_gmsh_geo_to_fenics_mesh(geo_filepath)

    # Solve a dummy problem

    V = fn.FunctionSpace(labelled_mesh.mesh, "CG", 2)
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    bcs = [fn.DirichletBC(V, fn.Constant(0), labelled_mesh.boundary_mesh_func, 2)]

    dx_mf = fn.dx(subdomain_data=labelled_mesh.subdomain_mesh_func)

    F = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
    F += fn.Constant(10) * u * v * fn.dx
    F -= fn.Constant(1) * v * dx_mf(2)

    #

    u = fn.Function(V)
    u = fu.solve_weak_form(u, F, bcs)

    assert round(0.0146543, 5) == round(u(0.8, 0.8, 0.5), 5)


