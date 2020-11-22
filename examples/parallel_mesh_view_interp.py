#!python3

###############################################################
# Imports
###############################################################

import os
from time import sleep
from datetime import datetime
import tempfile as tmp
import subprocess as sp

import dolfin as fn
import matplotlib.pyplot as plt


import fenics_utils as fu

###############################################################
# Definitions
###############################################################

GEO_FILEPATH = os.path.join(os.path.dirname(__file__), os.pardir, "geo",
                            "circle_square_center_pt.geo")
assert os.path.isfile(GEO_FILEPATH)

RESULTS_DIR = fu.get_clean_results_dir(__file__)

###############################################################
# Functions
###############################################################


def log(comm: fn.MPI.comm_world, msg: str) -> None:
    print(f"{datetime.now()}: log: rank {fn.MPI.rank(comm)}: {msg}")


###############################################################
# Main
###############################################################

print(
    f"mpi process: rank: {fn.MPI.rank(fn.MPI.comm_world)}, size: {fn.MPI.size(fn.MPI.comm_world)}")


comm = fn.MPI.comm_world
rank = fn.MPI.rank(comm)

geo_params = {}

lbl_mesh = fu.convert_2d_gmsh_geo_to_fenics_mesh(GEO_FILEPATH)

plt.figure()
fn.plot(lbl_mesh.subdomain_mesh_func)

if "PYTEST_CURRENT_TEST" not in os.environ:
    plt.show()

lbl_meshview = fu.create_mesh_view(lbl_mesh)

file = fn.File(comm, os.path.join(RESULTS_DIR, "mpi_meshview_out_2d.pvd"))
file << lbl_meshview.subdomain_mesh_func