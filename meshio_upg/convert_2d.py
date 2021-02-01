import meshio
import numpy as np

msh = meshio.read("square_in_circle.msh")

triangle_cells = []
line_cells = []
for cell in msh.cells:
    if cell.type == "triangle":
        if len(triangle_cells) == 0:
            triangle_cells = cell.data
        else:
            triangle_cells = np.vstack([triangle_cells, cell.data])
    elif cell.type == "line":
        if len(line_cells) == 0:
            line_cells = cell.data
        else:
            line_cells = np.vstack([line_cells, cell.data])

line_data = []
triangle_data = []
for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "line":
        if len(line_data) == 0:
            line_data = msh.cell_data_dict["gmsh:physical"][key]
        else:
            line_data = np.vstack([line_data, msh.cell_data_dict["gmsh:physical"][key]])
    elif key == "triangle":
        if len(triangle_data) == 0:
            triangle_data = msh.cell_data_dict["gmsh:physical"][key]
        else:
            triangle_data = np.vstack(
                [triangle_data, msh.cell_data_dict["gmsh:physical"][key]]
            )

if len(triangle_data) > 0:
    dom_cell_data = {"name_to_read": [triangle_data]}
else:
    dom_cell_data = {}

triangle_mesh = meshio.Mesh(
    points=msh.points[:, :2],
    cells={"triangle": triangle_cells},
    cell_data=dom_cell_data,
)

if len(line_data) > 0:
    bnd_cell_data = {"name_to_read": [line_data]}
else:
    bnd_cell_data = {}

line_mesh = meshio.Mesh(
    points=msh.points[:, :2], cells=[("line", line_cells)], cell_data=bnd_cell_data
)

meshio.write("mesh.xdmf", triangle_mesh)

meshio.xdmf.write("bmf.xdmf", line_mesh)

# load

import dolfin as fn

mesh = fn.Mesh()
with fn.XDMFFile(fn.MPI.comm_world, "mesh.xdmf") as xdmf_infile:
    xdmf_infile.read(mesh)

#

mvc_dom = fn.MeshValueCollection("size_t", mesh, 2)
with fn.XDMFFile(fn.MPI.comm_world, "mesh.xdmf") as xdmf_infile:
    xdmf_infile.read(mvc_dom, "name_to_read")

mf_dom = fn.MeshFunction("size_t", mesh, mvc_dom)

#

mvc_bnd = fn.MeshValueCollection("size_t", mesh, 1)
with fn.XDMFFile(fn.MPI.comm_world, "bmf.xdmf") as xdmf_infile:
    xdmf_infile.read(mvc_bnd, "name_to_read")

mf_bnd = fn.MeshFunction("size_t", mesh, mvc_bnd)

# plot

# import matplotlib.pyplot as plt

# plt.figure()
# c = fn.plot(mf_dom, title="domain mesh function")
# plt.colorbar(c)

# plt.figure()
# c = fn.plot(mf_bnd, title="domain mesh function")


# plt.show()
