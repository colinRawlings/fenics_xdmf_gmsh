# !python3

"""
Interestingly plots get scrambled if plt.show() is not called
and the same name is given to each mesh_data instance even with
the copy construction ...
"""

###########################################################
# Imports
###########################################################

import os
import tempfile

import matplotlib.pyplot as plt
import pyvista as pv

import dolfin as fn

import fenics_utils as fu

import pyvista as pv

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "geo")

###########################################################
# Main
###########################################################

geo_filepath = os.path.join(GEO_DIR, "param_geo.geo")
assert os.path.isfile(geo_filepath)

for geom_params in [
    {"dx_inner_mesh": 0.5},
    {"radius": 3.5, "x0": -2},
    {"radius": 2.5, "x0": 0, "y0": -1},
]:

    mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath, geom_params)

    if "PYTEST_CURRENT_TEST" not in os.environ:

        with tempfile.TemporaryDirectory() as tmp_dir:

            tmp_xdmf_filepath = os.path.join(tmp_dir, "tmp.xdmf")

            with fn.XDMFFile(tmp_xdmf_filepath) as tmp_xd:
                tmp_xd.write(mesh_data.subdomain_mesh_func)

            mesh = pv.read(tmp_xdmf_filepath)

            plotter = pv.Plotter()
            plotter.add_mesh_clip_box(mesh, show_edges=True)

            plotter.show()
