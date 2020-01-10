# !python3
"""
Create mesh view of the global mesh with a preserved mesh function

create expression from the global mesh function.

use this expression as a means of transferring the global mesh
function to the view's cells

"""

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

GEO_DIR = os.path.abspath("geo")

###########################################################
# Main
###########################################################

# Create meshes

msh_filepath = os.path.join(GEO_DIR, "circle_square_center_pt.geo")
mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(msh_filepath, do_plots=False)
global_mesh = mesh_data["mesh"]
global_smf = mesh_data["subdomain_mesh_func"]

mv_data = fu.create_mesh_view(global_mesh, global_smf, 1)

plt.figure()
fn.plot(mv_data["mesh"], title="mesh view")

plt.figure()
p = fn.plot(mv_data["subdomain_mesh_func"], title="smf view")
plt.colorbar(p)
