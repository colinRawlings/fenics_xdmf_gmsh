# !python

###########################################################
# Imports
###########################################################

import os

import dolfin as fn

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


def refine_and_project_to_cg1(u: fn.Function,
                              max_refinements: int = 2) -> fn.Function:
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
        Vr = fn.FunctionSpace(refined_mesh, 'CG', 1)
    else:
        assert u.num_sub_spaces() == u.geometric_dimension(),\
        "Only scalar or vector functions supported"
        Vr = fn.VectorFunctionSpace(refined_mesh, 'CG', 1)

    ur = fn.Function(Vr)
    ur.interpolate(u)

    return ur


def save_function(u: fn.Function, filepath: str):

    assert os.path.splitext(filepath)[-1] == '.pvd', "Expected a paraview file"

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    ur = refine_and_project_to_cg1(u)

    file = fn.File(fn.MPI.comm_world, filepath)
    file << ur