import pyvista as pv

# mesh = pv.read("mesh.xdmf")

# plotter = pv.Plotter()
# plotter.add_mesh(mesh)

# plotter.show()

mesh = pv.read("/home/fenics/work/tmp3.xdmf")

plotter = pv.Plotter()
plotter.add_mesh(mesh)

plotter.show()
