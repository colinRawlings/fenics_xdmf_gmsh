# !python3

###########################################################
# Imports
###########################################################

import os
import matplotlib.pyplot as plt

import dolfin as fn

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = os.path.join(os.path.dirname(__file__), "mesh")

###########################################################
# Main
###########################################################

geo_filepath = os.path.join(GEO_DIR, "cyl_sph_plane_tail_parametric.geo")
assert os.path.isfile(geo_filepath), f"Failed to find: {geo_filepath}"

mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath)
p = fn.plot(mesh_data.subdomain_mesh_func)
fn.plot(mesh_data.mesh)
plt.colorbar(p)
plt.show()