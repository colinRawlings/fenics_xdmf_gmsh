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
# import subprocess as sp
import matplotlib.pyplot as plt

import dolfin as fn

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "geo")

###########################################################
# Main
###########################################################

geo_filepath = os.path.join(GEO_DIR, "param_geo.geo")
assert os.path.isfile(geo_filepath)

for geom_params in [{"dx_inner_mesh": 0.5}, {"radius": 3.5, "x0": -2}, {"radius": 2.5, "x0": 0, "y0": -1}]:

    mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath, geom_params)
    p = fn.plot(mesh_data.subdomain_mesh_func)
    fn.plot(mesh_data.mesh)
    plt.title(f"geom params: {geom_params}")
    plt.colorbar(p)
    plt.show()

