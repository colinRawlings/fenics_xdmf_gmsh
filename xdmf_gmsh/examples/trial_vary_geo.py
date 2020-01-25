# !python3

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
# tmp_msh_filepath = "out.msh"

assert os.path.isfile(geo_filepath)

# if os.path.isfile(tmp_msh_filepath):
#     print(f"removing: {tmp_msh_filepath}")
#     os.remove(tmp_msh_filepath)

# result = sp.run(["gmsh", "-2", "-o", tmp_msh_filepath, "-setnumber", "lc", "3.5", geo_filepath])
# assert result.returncode == 0, f"{result.stdout.decode()}, \n{result.stderr.decode()}"

# mesh_data = fu.convert_2d_gmsh_msh_to_fenics_mesh(tmp_msh_filepath)

mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath, {"lc": 3.5, "x0": -0.5})

plt.figure()
fn.plot(mesh_data.mesh)
