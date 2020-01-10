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
import numpy

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

marker = fn.MeshFunction("size_t", global_mesh,
                         global_mesh.topology().dim(), 0)  # mark entirety of the full mesh
mesh1 = fn.MeshView.create(marker, 0)

plt.figure()
p = fn.plot(global_smf, title="smf")
plt.colorbar(p)

# make expression

V = fn.FunctionSpace(mesh1, "DG", 0)
u_smf = fn.Function(V)

helper = numpy.asarray(global_smf.array(), dtype=numpy.int32)

dm = V.dofmap()
for cell in fn.cells(mesh1):
    helper[dm.cell_dofs(cell.index())] = global_smf[cell]

u_smf.vector()[:] = helper

plt.figure()
p = fn.plot(u_smf, title="smf func")
plt.colorbar(p)

view_smf = fn.MeshFunction("size_t", mesh1, mesh1.topology().dim(), 0)
for c in fn.cells(mesh1):

    if mesh1.topology().dim() == 2:
        cell_midpoint = (c.midpoint().x(), c.midpoint().y())
    elif mesh1.topology().dim() == 3:
        cell_midpoint = (c.midpoint().x(), c.midpoint().y(), c.midpoint.z())
    else:
        assert False, "Unexpected condition"

    view_smf[c] = int(u_smf(*cell_midpoint))

plt.figure()
p = fn.plot(view_smf, title="view smf func")
plt.colorbar(p)
