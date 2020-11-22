#!python3
# %%
###############################################################
# Imports
###############################################################

import os
import pathlib as pt
from time import perf_counter
from configparser import ConfigParser

import numpy as np
import dolfin as fn
import matplotlib.pyplot as plt

import fenics_utils as fu

fu.reset()

# %%
###########################################################
# Definitions: Geometry
###########################################################

assert pt.Path("config.ini").is_file()
parser = ConfigParser()
parser.read(pt.Path("config.ini"))
PROJECT_DIR = pt.Path(parser["Config"]["project_dir"])

inner_square_boundary_index = 1
inner_circle_boundary_index = 2
outer_circle_boundary_index = 3

inner_square_domain_index = 1
inner_circle_domain_index = 2
outer_circle_domain_index = 3

square_side = 0.5
R_inner_circle = 1.0
R_outer_circle = 20

geo_filepath = PROJECT_DIR / "geo" / "square_in_circles.geo"
assert geo_filepath.is_file(), f"{geo_filepath}"

mesh_data = fu.convert_2d_gmsh_geo_to_fenics_mesh(
    str(geo_filepath),
    geo_params={
        "radius_1": R_inner_circle,
        "radius_2": R_outer_circle,
        "square_side": square_side,
    },
)

mesh_full = fu.create_mesh_view(mesh_data)
mesh_outer = fu.create_mesh_view(mesh_data, inner_circle_domain_index)

plt.figure()

plt.subplot(1, 2, 1)
fn.plot(mesh_full.subdomain_mesh_func)
fn.plot(mesh_full.mesh)

plt.subplot(1, 2, 2)
fn.plot(mesh_outer.subdomain_mesh_func)
fn.plot(mesh_outer.mesh)

# %%
###########################################################
# Definitions: Function Spaces
###########################################################

V = fn.FunctionSpace(mesh_full.mesh, "CG", 2)
axt = fn.TrialFunction(V)
vxt = fn.TestFunction(V)

V_outer = fn.FunctionSpace(mesh_outer.mesh, "CG", 2)
axt_outer = fn.TrialFunction(V_outer)
vxt_outer = fn.TestFunction(V_outer)

# %%
###########################################################
# Definitions: Equation ax (ref)
###########################################################

bc = [fn.DirichletBC(V, fn.Constant(0), fu.IsBoundary())]

dx_mf = fn.dx(subdomain_data=mesh_full.subdomain_mesh_func)

x, y = fn.SpatialCoordinate(mesh_full.mesh)

jx = x
lhs = -1 * fn.inner(fn.grad(axt), fn.grad(vxt)) * fn.dx  # type: ignore
rhs = vxt * jx * dx_mf(inner_square_domain_index)  # type: ignore

ax = fn.Function(V)
fn.solve(
    lhs == rhs,
    ax,
    bc,
    solver_parameters={"linear_solver": "gmres", "preconditioner": "ilu"},
)

# %%
###########################################################
# Definitions: Equation ax_outer: Dirichlet
###########################################################

bc_outer_dirichlet = [
    fn.DirichletBC(
        V_outer,
        fn.Constant(0),
        mesh_outer.boundary_mesh_func,
        inner_circle_boundary_index,
    ),
    fn.DirichletBC(
        V_outer, ax, mesh_outer.boundary_mesh_func, inner_square_boundary_index
    ),
]

lhs = -1 * fn.inner(fn.grad(axt_outer), fn.grad(vxt_outer)) * fn.dx  # type: ignore
rhs = vxt_outer * fn.Constant(0) * fn.dx  # type: ignore

ax_outer_dirichlet = fn.Function(V_outer)
fn.solve(
    lhs == rhs,
    ax_outer_dirichlet,
    bc_outer_dirichlet,
    solver_parameters={"linear_solver": "gmres", "preconditioner": "ilu"},
)

# %%
###########################################################
# Definitions: Equation ax_outer: Neumann
###########################################################

bc_outer_neumann = [
    fn.DirichletBC(
        V_outer, ax, mesh_outer.boundary_mesh_func, inner_square_boundary_index
    )
]

lhs = -1 * fn.inner(fn.grad(axt_outer), fn.grad(vxt_outer)) * fn.dx  # type: ignore
rhs = vxt_outer * fn.Constant(0) * fn.dx  # type: ignore

ax_outer_neumann = fn.Function(V_outer)
fn.solve(
    lhs == rhs,
    ax_outer_neumann,
    bc_outer_neumann,
    solver_parameters={"linear_solver": "gmres", "preconditioner": "ilu"},
)

# %%
###########################################################
# PostProcessing:
###########################################################

ax_outer_ref = fn.project(ax, V_outer)

plt.figure()

plt.subplot(4, 2, 1)
p = fn.plot(ax, title="ref")
b = plt.colorbar(p)

plt.subplot(4, 2, 2)
p = fn.plot(ax_outer_ref, title="ref")
plt.colorbar(p)

plt.subplot(4, 2, 3)
p = fn.plot(ax_outer_dirichlet, title="outer-dirichlet")
plt.colorbar(p)

plt.subplot(4, 2, 4)
p = fn.plot(ax_outer_neumann, title="outer-neumann")
plt.colorbar(p)

plt.subplot(4, 2, 5)
p = fn.plot(
    0.5 * ax_outer_dirichlet + 0.5 * ax_outer_neumann, title="outer-average"
)  # type: ignore
plt.colorbar(p)

plt.subplot(4, 2, 6)
p = fn.plot(
    ax_outer_ref - 0.5 * ax_outer_dirichlet - 0.5 * ax_outer_neumann,  # type: ignore
    title="error",
)
plt.colorbar(p)

# Cross sections

# y = 0

x_xXX_outer = np.linspace(0.5 * square_side + 0.01, R_inner_circle - 0.01)

ax_outer_dirichlet_xXX = [ax_outer_dirichlet(x_val, 0) for x_val in x_xXX_outer]
ax_outer_neumann_xXX = [ax_outer_neumann(x_val, 0) for x_val in x_xXX_outer]

ax_avg_xXX = [
    0.5 * ax_outer_dirichlet(x_val, 0) + 0.5 * ax_outer_neumann(x_val, 0)
    for x_val in x_xXX_outer
]

ax_xXX = [ax(x_val, 0) for x_val in x_xXX_outer]

plt.subplot(4, 2, 7)

plt.plot(x_xXX_outer, ax_xXX, "-b", label="full")
plt.plot(x_xXX_outer, ax_outer_dirichlet_xXX, ":g", label="dirichlet")
plt.plot(x_xXX_outer, ax_outer_neumann_xXX, "--g", label="neumann")
plt.plot(x_xXX_outer, ax_avg_xXX, ":r", label="avg")
plt.xlabel("$x$")
plt.ylabel("$a_x(x,0)$")
plt.legend()

# x = y

x_xyXX_outer = np.linspace(
    0.5 * square_side + 0.01, R_inner_circle / np.sqrt(2) - 0.01
)  # type: ignore

ax_outer_dirichlet_xyXX = [ax_outer_dirichlet(x_val, x_val) for x_val in x_xyXX_outer]
ax_outer_neumann_xyXX = [ax_outer_neumann(x_val, x_val) for x_val in x_xyXX_outer]

ax_avg_xyXX = [
    0.5 * ax_outer_dirichlet(x_val, x_val) + 0.5 * ax_outer_neumann(x_val, x_val)
    for x_val in x_xyXX_outer
]

ax_xyXX = [ax(x_val, x_val) for x_val in x_xyXX_outer]

plt.subplot(4, 2, 8)

plt.plot(x_xyXX_outer, ax_xyXX, "-b", label="full")
plt.plot(x_xyXX_outer, ax_outer_dirichlet_xyXX, ":g", label="dirichlet")
plt.plot(x_xyXX_outer, ax_outer_neumann_xyXX, "--g", label="neumann")
plt.plot(x_xyXX_outer, ax_avg_xyXX, ":r", label="avg")
plt.xlabel("$x$")
plt.ylabel("$a_x(x, x)$")
plt.legend()
plt.tight_layout()

if "PYTEST_CURRENT_TEST" not in os.environ:
    plt.show()
