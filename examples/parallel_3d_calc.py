# !python3
"""
Here:

- trial using fenicstools to access solution values for a
  parallellised calculation

run as:

mpirun -np 4 /usr/bin/python3 <this_file>.py

"""

###########################################################
# imports
###########################################################

from time import time
import dolfin as fn
import numpy as np

import fenics_utils as fu
import argparse as ap
import json

###########################################################
# Main
###########################################################

parser = ap.ArgumentParser('Simulate trivial poisson problem')
parser.add_argument('-N',
                    '--num-cells',
                    default=20,
                    type=int,
                    help="number of cells (per cube side) in the mesh.")
parser.add_argument('--linear-solver', default='gmres', help="Linear solver (e.g. lu).")
parser.add_argument('--preconditioner', default='none', help="Preconditioner to use (e.g. none)")
parser.add_argument('--check-solution', action='store_true')
parser.add_argument('--results-filepath',
                    default='',
                    help='json file to save the result data to (if empty no data is saved).')

args, _ = parser.parse_known_args()

# solve

comm = fn.MPI.comm_world
rank = fn.MPI.rank(comm)

mesh = fn.UnitCubeMesh(args.num_cells, args.num_cells, args.num_cells)

V = fn.FunctionSpace(mesh, 'CG', 2)

u = fn.TrialFunction(V)
v = fn.TestFunction(V)

a = fn.Constant(-1) * fn.inner(fn.grad(u), fn.grad(v)) * fn.dx
L = fn.Constant(1) * v * fn.dx

# construct for anl solution
x = fn.SpatialCoordinate(mesh)
anl_soln = fn.Constant(1 / 6) * (x[0]**2 + x[1]**2 + x[2]**2)  # type: ignore
bc = fn.DirichletBC(V, anl_soln, fu.IsBoundary())

#

print(
    f"mpi process: rank: {fn.MPI.rank(fn.MPI.comm_world)}, size: {fn.MPI.size(fn.MPI.comm_world)}")

t_start_solve = time()
u = fn.Function(V)
fn.solve(a == L,
         u,
         bc,
         solver_parameters={
             'linear_solver': args.linear_solver,
             'preconditioner': args.preconditioner
         })
solve_time_s = time() - t_start_solve

# report

if rank == 0:
    print(f"Solution took: {solve_time_s: .3f} s for {V.dim()} DOF")

if rank == 0 and args.results_filepath:
    with open(args.results_filepath, 'w') as f:
        json.dump({'dof': V.dim(), 'solve_time_rank0_s': solve_time_s}, f, sort_keys=True, indent=4)

# post pro:

if rank == 0 and args.check_solution:
    import fenicstools as ft
    u_anl = fn.project(anl_soln, V)

    x = np.array([[0.5, 0.5, 0.5]])

    probes_u = ft.Probes(x.flatten(), V)
    probes_u(u)
    u_fem = probes_u.array()

    probes_u_expected = ft.Probes(x.flatten(), V)
    probes_u_expected(u_anl)
    u_expected = probes_u_expected.array()

    assert abs(u_fem - u_expected) < 1e-3 * abs(u_expected)  # type: ignore

if rank == 0:  # is root process
    print("Success!")
