"""
Some helpful code for simulations with fenics
"""

#############################################################################
# IMPORTS
#############################################################################

import os
import shutil
import subprocess as sp
import typing as ty
import tempfile

import matplotlib.pyplot as plt

import meshio
import dolfin as fn

#############################################################################
# CONSTANTS
#############################################################################

PARAVIEW_TMP_FOLDERPATH = os.path.abspath("./paraview_tmp")

#############################################################################
# FUNCTIONS
#############################################################################


def reset(matplotlib_mode='inline'):
    """
    If a notebook session call %reset and set matplotlib magic mode.
    This must be called _first_ as it will "erase" the history of any
    previous imports (if called in a notebook - otherwise it will do
    nothing).

    :param matplotlib_mode: {'inline', 'notebook', 'ipympl', ... }
    :return:
    """

    try:
        import IPython

        if IPython.get_ipython():
            IPython.get_ipython().magic('reset -f')

            import IPython

            IPython.get_ipython().magic(f'matplotlib {matplotlib_mode}')

            IPython.get_ipython().magic('reload_ext autoreload')
            IPython.get_ipython().magic('autoreload 2')
            IPython.get_ipython().magic('aimport -dolfin -fenics')

            print("Configured magics")
        else:
            print("Could not configure magics")

    except ImportError:
        print("Could not import IPython")


#############################################################################


def clean_output_folder(mesh_folderpath: str) -> None:
    """

    :param mesh_folderpath:
    :return:
    """

    if os.path.isdir(mesh_folderpath):
        shutil.rmtree(mesh_folderpath)

    os.mkdir(mesh_folderpath)


#############################################################################


def close_all_paraview() -> None:
    sp.call(['killall', 'paraview'])


#############################################################################


def convert_2d_gmsh_msh_to_fenics_msh(msh_filepath: str, do_plots=False) -> ty.Dict:
    """[convert_2d_gmsh_msh_to_fenics_msh]
    
    Arguments:
        msh_filepath {str} --
        do_plots {bool} 
    
    Returns:
        ty.Dict -- dict(mesh=mesh, subdomain_mesh_func=mf_dom, boundary_mesh_func=mf_bnd)
    """

    msh = meshio.read(msh_filepath)

    #

    with tempfile.TemporaryDirectory() as temp_dir:

        tmp_domain_filepath = os.path.join(os.path.join(temp_dir, "domains.xdmf"))
        tmp_boundary_filepath = os.path.join(os.path.join(temp_dir, "boundaries.xdmf"))

        meshio_dom = meshio.Mesh(
            points=msh.points[:, :2],  # Converting to 2D
            cells={"triangle": msh.cells["triangle"]},
            cell_data={"triangle": {
                "subdomain": msh.cell_data["triangle"]["gmsh:physical"]
            }},
            field_data=msh.field_data)

        meshio.write(tmp_domain_filepath, meshio_dom)

        #

        meshio_bnd = meshio.Mesh(
            points=msh.points[:, :2],  # Converting to 2D
            cells={"line": msh.cells["line"]},
            cell_data={"line": {
                "boundaries": msh.cell_data["line"]["gmsh:physical"]
            }})

        meshio.write(tmp_boundary_filepath, meshio_bnd)

        # Load into fenics

        mesh = fn.cpp.mesh.Mesh()
        with fn.XDMFFile(fn.MPI.comm_world, tmp_domain_filepath) as xdmf_infile:
            xdmf_infile.read(mesh)

        #

        mvc_dom = fn.MeshValueCollection("size_t", mesh, 2)
        with fn.XDMFFile(fn.MPI.comm_world, tmp_domain_filepath) as xdmf_infile:
            xdmf_infile.read(mvc_dom, "subdomain")

        mf_dom = fn.MeshFunction("size_t", mesh, mvc_dom)

        #

        mvc_bnd = fn.MeshValueCollection("size_t", mesh, 1)
        with fn.XDMFFile(fn.MPI.comm_world, tmp_boundary_filepath) as xdmf_infile:
            xdmf_infile.read(mvc_bnd, "boundaries")

        mf_bnd = fn.MeshFunction("size_t", mesh, mvc_bnd)

    if do_plots:
        plt.figure()
        c = fn.plot(mf_dom, title="domain mesh function")
        plt.colorbar(c)

        plt.figure()
        c = fn.plot(mesh, title="mesh")


    return dict(mesh=mesh, subdomain_mesh_func=mf_dom, boundary_mesh_func=mf_bnd)


#############################################################################


def convert_3d_gmsh_msh_to_fenics_msh(msh_filepath: str) -> ty.Dict:
    """[convert_3d_gmsh_msh_to_fenics_msh]
    
    Arguments:
        msh_filepath {str} -- 
    
    Returns:
        ty.Dict -- dict(mesh=mesh, subdomain_mesh_func=mf_dom, boundary_mesh_func=mf_bnd)
    """

    msh = meshio.read(msh_filepath)

    #

    with tempfile.TemporaryDirectory() as temp_dir:

        tmp_domain_filepath = os.path.join(os.path.join(temp_dir, "domains.xdmf"))
        tmp_boundary_filepath = os.path.join(os.path.join(temp_dir, "boundaries.xdmf"))

        meshio_dom = meshio.Mesh(
            points=msh.points,  # Converting to 2D
            cells={"tetra": msh.cells["tetra"]},
            cell_data={"tetra": {
                "subdomain": msh.cell_data["tetra"]["gmsh:physical"]
            }},
            field_data=msh.field_data)

        meshio.write(tmp_domain_filepath, meshio_dom)

        #

        meshio_bnd = meshio.Mesh(
            points=msh.points,  # Converting to 2D
            cells={"triangle": msh.cells["triangle"]},
            cell_data={"triangle": {
                "boundaries": msh.cell_data["triangle"]["gmsh:physical"]
            }})

        meshio.write(tmp_boundary_filepath, meshio_bnd)

        # Load into fenics

        mesh = fn.cpp.mesh.Mesh()
        with fn.XDMFFile(fn.MPI.comm_world, tmp_domain_filepath) as xdmf_infile:
            xdmf_infile.read(mesh)

        #

        mvc_dom = fn.MeshValueCollection("size_t", mesh, 2)
        with fn.XDMFFile(fn.MPI.comm_world, tmp_domain_filepath) as xdmf_infile:
            xdmf_infile.read(mvc_dom, "subdomain")

        mf_dom = fn.MeshFunction("size_t", mesh, mvc_dom)

        #

        mvc_bnd = fn.MeshValueCollection("size_t", mesh, 1)
        with fn.XDMFFile(fn.MPI.comm_world, tmp_boundary_filepath) as xdmf_infile:
            xdmf_infile.read(mvc_bnd, "boundaries")

        mf_bnd = fn.MeshFunction("size_t", mesh, mvc_bnd)

    return dict(mesh=mesh, subdomain_mesh_func=mf_dom, boundary_mesh_func=mf_bnd)


#############################################################################


def solve_weak_form(u, F, boundary_conditions=[]):
    """

    A flexible solution strategy (not optimal in case of a linear problem)

    :param u: field
    :param F: weak form F(u, w) = 0
    :param boundary_conditions: List of boundary conditions
    :return: u (n.b. not a copy just a reference to the field supplied as input)
    """

    import dolfin as fn

    R = fn.action(F, u)
    dR = fn.derivative(R, u)

    problem = fn.NonlinearVariationalProblem(R, u, boundary_conditions, dR)
    solver = fn.NonlinearVariationalSolver(problem)
    solver.parameters['newton_solver']['relative_tolerance'] = 1e-6
    solver.parameters['newton_solver']['relaxation_parameter'] = 1.
    solver.parameters['newton_solver']['error_on_nonconvergence'] = False
    solver.solve()

    return u


#############################################################################


def save_mesh(mesh, mesh_filepath=None, open_file=False):
    """

    :param mesh:
    :param mesh_filepath:
    :param open_file: try to open the saved mesh in paraview
    :return:
    """
    import dolfin as fn

    if not isinstance(mesh, fn.Mesh):
        raise ValueError("Expected a dolfin.Mesh")

    if mesh_filepath is None:
        mesh_filepath = os.path.join(PARAVIEW_TMP_FOLDERPATH, "mesh.pvd")

    mesh_filepath = os.path.abspath(mesh_filepath)

    if os.path.isfile(mesh_filepath):
        os.remove(mesh_filepath)

    mesh_file = fn.File(mesh_filepath)

    mesh_file << mesh

    print(f"Wrote mesh to: {mesh_filepath}")

    if open_file:
        sp.Popen(['paraview', mesh_filepath])

    return mesh_filepath


#############################################################################


def save_function(u, function_filepath=None, open_file=False):
    """

    :param u: (Function or Boundary Mesh Function or Subdomain Mesh Function)
    :param function_filepath: location to store the temporary pvd file
    :param open_file: try to open the saved function in paraview

    :return: function_filepath (absolute)
    """

    import dolfin as fn

    if function_filepath is None:
        function_filepath = os.path.join(PARAVIEW_TMP_FOLDERPATH, "function.pvd")

    function_filepath = os.path.abspath(function_filepath)

    if os.path.isfile(function_filepath):
        os.remove(function_filepath)

    result_file = fn.File(function_filepath)

    result_file << u

    print(f"Wrote function to: {function_filepath}")

    if open_file:
        sp.Popen(['paraview', function_filepath])

    return function_filepath
