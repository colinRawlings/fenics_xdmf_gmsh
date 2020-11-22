# !python

###########################################################
# Imports
###########################################################

import os

import dolfin as fn
import ufl

###########################################################
# Definitions
###########################################################


def detect_function_degree(u: fn.Function) -> int:
    """Look up the degree of the element providing support
    for the supplied function
    
    Arguments:
        u {fn.Function} -- [description]
    
    Returns:
        int -- [description]
    """

    return u.ufl_element().degree()


#############################################################################


def refine_and_project_to_cg1(u: fn.Function, max_refinements: int = 2) -> fn.Function:
    """Project a scalar or vector function supported
    on higher order elements onto first order elements 
    defined on a suitably refined mesh.
    
    Arguments:
        u {fn.Function} -- [description]
    
    Returns:
        fn.Function -- [description]
    """

    num_refinements = detect_function_degree(u) - 1
    if num_refinements < 1:
        return u

    num_refinements = min(num_refinements, max_refinements) + 1

    refined_mesh = u.function_space().mesh()

    while num_refinements > 0:
        refined_mesh = fn.refine(refined_mesh)
        num_refinements -= 1

    if u.num_sub_spaces() == 0:  # scalar
        Vr = fn.FunctionSpace(refined_mesh, "CG", 1)
    else:
        assert (
            u.num_sub_spaces() == u.geometric_dimension()
        ), "Only scalar or vector functions supported"
        Vr = fn.VectorFunctionSpace(refined_mesh, "CG", 1)

    ur = fn.Function(Vr)
    ur.interpolate(u)

    return ur


#############################################################################


def save_function(u: fn.Function, filepath: str):

    assert os.path.splitext(filepath)[-1] == ".pvd", "Expected a paraview file"

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    ur = refine_and_project_to_cg1(u)

    file = fn.File(fn.MPI.comm_world, filepath)
    file << ur


#############################################################################


def calculate_gradient(u: fn.Function) -> fn.Function:
    """Compute the gradient for use during post processing
    
    Arguments:
        u {fn.Function} -- [description]
    
    Returns:
        fn.Function -- [description]
    """

    assert u.num_sub_spaces() == 0, "Function does not appear to be a scalar"

    V = fn.VectorFunctionSpace(u.function_space().mesh(), "CG", 1)
    i = ufl.Index()

    u_H = fn.TrialFunction(V)
    v_H = fn.TestFunction(V)

    a = u_H[i] * v_H[i] * fn.dx  # type: ignore
    L = fn.Dx(u, i) * v_H[i] * fn.dx  # type: ignore

    u_H = fn.Function(V)
    fn.solve(a == L, u_H, [])

    return u_H
