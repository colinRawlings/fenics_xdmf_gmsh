import matplotlib.pyplot as plt

import dolfin as fn


def newton_solver_parameters():
    return {
        "nonlinear_solver": "newton",
        "newton_solver": {
            "linear_solver": "gmres"
        }
    }

class IsBoundary(fn.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


mesh = fn.UnitSquareMesh(10, 10)

V = fn.FunctionSpace(mesh, "CG", 2)

u = fn.Function(V)
v = fn.TrialFunction(V)

bc = fn.DirichletBC(V, fn.Constant(0), IsBoundary())

F = fn.inner(fn.grad(u), fn.grad(v)) * fn.dx + u**3 * v * fn.dx + fn.Constant(100) * v * fn.dx

fn.solve(F == 0, u, bc, solver_parameters=newton_solver_parameters())

plt.figure()
p = fn.plot(u, title="u")
plt.colorbar(p)
