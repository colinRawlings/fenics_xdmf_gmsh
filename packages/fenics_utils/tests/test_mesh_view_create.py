# !python3

###########################################################
# imports
###########################################################

import os
import pkg_resources as rc

import numpy as np

import dolfin as fn

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = rc.resource_filename("fenics_utils", "assets")  # type: ignore

###########################################################
# Test
###########################################################


def test_view_on_full_mesh_2d():

    msh_filepath = os.path.join(GEO_DIR, "simple_poisson_2d.geo")
    mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(msh_filepath, do_plots=False)

    mesh_data1 = fu.create_mesh_view(mesh_data["mesh"], mesh_data["subdomain_mesh_func"])

    smf_view = mesh_data1["subdomain_mesh_func"]
    smf_vals = np.asarray(smf_view.array(), dtype=np.int32)  # type: ignore

    vals = np.sort(np.unique(smf_vals))
    assert (vals == np.asarray([1, 2])).all()

def test_view_on_sub_mesh_2d():

    msh_filepath = os.path.join(GEO_DIR, "simple_poisson_2d.geo")
    mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(msh_filepath, do_plots=False)

    mesh_data1 = fu.create_mesh_view(mesh_data["mesh"], mesh_data["subdomain_mesh_func"], 2)

    smf_view = mesh_data1["subdomain_mesh_func"]
    smf_vals = np.asarray(smf_view.array(), dtype=np.int32)  # type: ignore

    vals = np.sort(np.unique(smf_vals))
    assert (vals == np.asarray([2])).all()

def test_view_on_full_mesh_3d():

    msh_filepath = os.path.join(GEO_DIR, "simple_poisson_3d.geo")
    mesh_data = fu.convert_3d_gmsh_geo_to_fenics_mesh(msh_filepath)

    mesh_data1 = fu.create_mesh_view(mesh_data["mesh"], mesh_data["subdomain_mesh_func"])

    smf_view = mesh_data1["subdomain_mesh_func"]
    smf_vals = np.asarray(smf_view.array(), dtype=np.int32)  # type: ignore

    vals = np.sort(np.unique(smf_vals))
    assert (vals == np.asarray([1, 2])).all()

def test_view_on_sub_mesh_3d():

    msh_filepath = os.path.join(GEO_DIR, "simple_poisson_3d.geo")
    mesh_data = fu.convert_3d_gmsh_geo_to_fenics_mesh(msh_filepath)

    mesh_data1 = fu.create_mesh_view(mesh_data["mesh"], mesh_data["subdomain_mesh_func"], 2)

    smf_view = mesh_data1["subdomain_mesh_func"]
    smf_vals = np.asarray(smf_view.array(), dtype=np.int32)  # type: ignore

    vals = np.sort(np.unique(smf_vals))
    assert (vals == np.asarray([2])).all()
    