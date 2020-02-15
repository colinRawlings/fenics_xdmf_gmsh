#!python3

###############################################################
# Imports
###############################################################

import os
from datetime import datetime

import dolfin as fn

import fenics_utils as fu

###############################################################
# Definitions
###############################################################

GEO_FILEPATH = os.path.join(os.path.dirname(__file__), os.pardir, "geo",
                            "circle_square_center_pt.geo")
assert os.path.isfile(GEO_FILEPATH)

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

#
# Solve
#

V = fn.FunctionSpace(lbl_mesh.mesh, "CG", 2)
u = fn.Function(V)
v = fn.TrialFunction(V)

dx_mf = fn.dx(subdomain_data=lbl_mesh.subdomain_mesh_func)
bc = fn.DirichletBC(V, fn.Constant(0), lbl_mesh.boundary_mesh_func, 1)

F = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
F += fn.Constant(1) * v * dx_mf(2)

fn.solve(F == 0, u, bc)

file = fn.File(comm, "mpi_out_2d.pvd")
file << u