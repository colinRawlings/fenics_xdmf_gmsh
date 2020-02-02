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

import numpy as np
import matplotlib.pyplot as plt

import meshio
import dolfin as fn

from typing import Optional

#############################################################################
# CONSTANTS
#############################################################################

PARAVIEW_TMP_FOLDERPATH = os.path.abspath("./paraview_tmp")

#############################################################################
# DEFINITIONS
#############################################################################


class LabelledMesh():
    def __init__(self,
                 mesh: fn.Mesh,
                 subdomain_mesh_func: fn.MeshFunction,
                 boundary_mesh_func: Optional[fn.MeshFunction] = None) -> None:

        # copy inputs

        self._mesh = fn.Mesh(mesh)
        assert id(self._mesh) != id(mesh)

        #

        with tempfile.TemporaryDirectory() as temp_dir:

            tmp_mf_filepath = os.path.join(str(temp_dir), "domain_mf.xdmf")

            with fn.XDMFFile(fn.MPI.comm_world,
                             tmp_mf_filepath) as tmp_file_out:
                tmp_file_out.write(subdomain_mesh_func)

            self._subdomain_mesh_func = fn.MeshFunction(
                "size_t", self._mesh, self._mesh.geometric_dimension())
            with fn.XDMFFile(fn.MPI.comm_world, tmp_mf_filepath) as tmp_file_in:
                tmp_file_in.read(self._subdomain_mesh_func)

        assert id(self._subdomain_mesh_func) != id(subdomain_mesh_func)

        #

        self._boundary_mesh_func = None

        if boundary_mesh_func is None:
            return

        with tempfile.TemporaryDirectory() as temp_dir:

            tmp_mf_filepath = os.path.join(str(temp_dir), "boundary_mf.xdmf")

            with fn.XDMFFile(fn.MPI.comm_world,
                             tmp_mf_filepath) as tmp_file_out:
                tmp_file_out.write(boundary_mesh_func)

            self._boundary_mesh_func = fn.MeshFunction(
                "size_t", self._mesh,
                self._mesh.geometric_dimension() - 1)
            with fn.XDMFFile(fn.MPI.comm_world, tmp_mf_filepath) as tmp_file_in:
                tmp_file_in.read(self._boundary_mesh_func)

        assert id(self._boundary_mesh_func) != id(boundary_mesh_func)

    @property
    def mesh(self) -> fn.Mesh:
        return self._mesh

    @property
    def subdomain_mesh_func(self) -> fn.MeshFunction:
        return self._subdomain_mesh_func

    @property
    def boundary_mesh_func(self) -> Optional[fn.MeshFunction]:
        return self._boundary_mesh_func


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


def convert_2d_gmsh_msh_to_fenics_mesh(msh_filepath: str,
                                       do_plots=False) -> LabelledMesh:
    """[convert_2d_gmsh_msh_to_fenics_mesh]
    
    Arguments:
        msh_filepath {str} --
        do_plots {bool} 
    
    Returns:
        LabelledMesh
    """

    msh = meshio.read(msh_filepath)

    #

    with tempfile.TemporaryDirectory() as temp_dir:

        tmp_domain_filepath = os.path.join(
            os.path.join(str(temp_dir), "domains.xdmf"))
        tmp_boundary_filepath = os.path.join(
            os.path.join(str(temp_dir), "boundaries.xdmf"))

        meshio_dom = meshio.Mesh(
            points=msh.points[:, :2],  # Converting to 2D
            cells={"triangle": msh.cells["triangle"]},
            cell_data={
                "triangle": {
                    "subdomain": msh.cell_data["triangle"]["gmsh:physical"]
                }
            },
            field_data=msh.field_data)

        meshio.write(tmp_domain_filepath, meshio_dom)

        #

        meshio_bnd = meshio.Mesh(
            points=msh.points[:, :2],  # Converting to 2D
            cells={"line": msh.cells["line"]},
            cell_data={
                "line": {
                    "boundaries": msh.cell_data["line"]["gmsh:physical"]
                }
            })

        meshio.write(tmp_boundary_filepath, meshio_bnd)

        # Load into fenics

        mesh = fn.Mesh()
        with fn.XDMFFile(fn.MPI.comm_world, tmp_domain_filepath) as xdmf_infile:
            xdmf_infile.read(mesh)

        #

        mvc_dom = fn.MeshValueCollection("size_t", mesh, 2)
        with fn.XDMFFile(fn.MPI.comm_world, tmp_domain_filepath) as xdmf_infile:
            xdmf_infile.read(mvc_dom, "subdomain")

        mf_dom = fn.MeshFunction("size_t", mesh, mvc_dom)

        #

        mvc_bnd = fn.MeshValueCollection("size_t", mesh, 1)
        with fn.XDMFFile(fn.MPI.comm_world,
                         tmp_boundary_filepath) as xdmf_infile:
            xdmf_infile.read(mvc_bnd, "boundaries")

        mf_bnd = fn.MeshFunction("size_t", mesh, mvc_bnd)

    if do_plots:
        plt.figure()
        c = fn.plot(mf_dom, title="domain mesh function")
        plt.colorbar(c)

        plt.figure()
        c = fn.plot(mesh, title="mesh")

    return LabelledMesh(mesh=mesh,
                        subdomain_mesh_func=mf_dom,
                        boundary_mesh_func=mf_bnd)


#############################################################################


def convert_3d_gmsh_msh_to_fenics_mesh(msh_filepath: str) -> LabelledMesh:
    """[convert_3d_gmsh_msh_to_fenics_mesh]
    
    Arguments:
        msh_filepath {str} -- 
    
    Returns:
        LabelledMesh
    """

    msh = meshio.read(msh_filepath)

    #

    with tempfile.TemporaryDirectory() as temp_dir:

        tmp_domain_filepath = os.path.join(
            os.path.join(str(temp_dir), "domains.xdmf"))
        tmp_boundary_filepath = os.path.join(
            os.path.join(str(temp_dir), "boundaries.xdmf"))

        meshio_dom = meshio.Mesh(
            points=msh.points,  # Converting to 2D
            cells={"tetra": msh.cells["tetra"]},
            cell_data={
                "tetra": {
                    "subdomain": msh.cell_data["tetra"]["gmsh:physical"]
                }
            },
            field_data=msh.field_data)

        meshio.write(tmp_domain_filepath, meshio_dom)

        #

        meshio_bnd = meshio.Mesh(
            points=msh.points,  # Converting to 2D
            cells={"triangle": msh.cells["triangle"]},
            cell_data={
                "triangle": {
                    "boundaries": msh.cell_data["triangle"]["gmsh:physical"]
                }
            })

        meshio.write(tmp_boundary_filepath, meshio_bnd)

        # Load into fenics

        mesh = fn.Mesh()
        with fn.XDMFFile(fn.MPI.comm_world, tmp_domain_filepath) as xdmf_infile:
            xdmf_infile.read(mesh)

        #

        mvc_dom = fn.MeshValueCollection("size_t", mesh, 2)
        with fn.XDMFFile(fn.MPI.comm_world, tmp_domain_filepath) as xdmf_infile:
            xdmf_infile.read(mvc_dom, "subdomain")

        mf_dom = fn.MeshFunction("size_t", mesh, mvc_dom)

        #

        mvc_bnd = fn.MeshValueCollection("size_t", mesh, 1)
        with fn.XDMFFile(fn.MPI.comm_world,
                         tmp_boundary_filepath) as xdmf_infile:
            xdmf_infile.read(mvc_bnd, "boundaries")

        mf_bnd = fn.MeshFunction("size_t", mesh, mvc_bnd)

    return LabelledMesh(mesh=mesh,
                        subdomain_mesh_func=mf_dom,
                        boundary_mesh_func=mf_bnd)


#############################################################################


def _check_mesh_conversion_result(result: sp.CompletedProcess,
                                  tmp_msh_filepath: str) -> None:
    """Check the result of running gmsh to create a msh file from a geo file

        n.b. gmsh does not give an non-zero exit code on meshing failure)

    
    Arguments:
        
        result {sp.CompletedProcess} -- [description]
        tmp_msh_filepath {str} -- [description]

    """

    msg = "gmsh failed to create msh\n"
    msg += f"stdout:\n{result.stdout.decode()}\n"
    msg += f"stderr:\n{result.stderr.decode()}"
    assert result.returncode == 0, msg

    msh = meshio.read(tmp_msh_filepath)
    assert msh.points.shape[0] != 0, msg


#############################################################################


def _construct_param_args(
        geo_params: ty.Optional[ty.Dict["str", float]]) -> ty.List[str]:
    """convert the params in the supplied dict into the corresponding
    command line args for gmsh
    
    Arguments:
        geo_params {ty.Optional[ty.Dict[} -- [description]
    
    Returns:
        ty.List[str] -- [description]
    """

    param_args = []

    if not geo_params:
        return param_args

    for geo_param in geo_params:
        param_args += ["-setnumber", geo_param, str(geo_params[geo_param])]

    return param_args


#############################################################################


def convert_2d_gmsh_geo_to_fenics_mesh(geo_filepath: str,
                                       geo_params: ty.Optional[ty.Dict[
                                           "str", float]] = None,
                                       do_plots: bool = False) -> LabelledMesh:
    """[convert_2d_gmsh_geo_to_fenics_mesh]
    
    Arguments:
        geo_filepath {str} -- 
    
    Returns:
        LabelledMesh
    """

    assert os.path.isfile(geo_filepath), f"{geo_filepath} not found"

    with tempfile.TemporaryDirectory() as temp_dir:
        tmp_msh_filepath = os.path.join(str(temp_dir), "tmp_msh.msh")

        param_args = _construct_param_args(geo_params)
        cmd_list = ["gmsh", "-2"
                    ] + param_args + ["-o", tmp_msh_filepath, geo_filepath]

        result = sp.run(cmd_list, stdout=sp.PIPE, stderr=sp.PIPE)

        _check_mesh_conversion_result(result, tmp_msh_filepath)

        mesh_data = convert_2d_gmsh_msh_to_fenics_mesh(tmp_msh_filepath,
                                                       do_plots=do_plots)

    return mesh_data


#############################################################################


def convert_3d_gmsh_geo_to_fenics_mesh(
        geo_filepath: str,
        geo_params: ty.Optional[ty.Dict["str", float]] = None) -> LabelledMesh:
    """[convert_3d_gmsh_geo_to_fenics_mesh]
    
    Arguments:
        geo_filepath {str} -- 
    
    Returns:
        LabelledMesh
    """

    assert os.path.isfile(geo_filepath), f"{geo_filepath} not found"

    with tempfile.TemporaryDirectory() as temp_dir:
        tmp_msh_filepath = os.path.join(str(temp_dir), "tmp_msh.msh")

        param_args = _construct_param_args(geo_params)
        cmd_list = ["gmsh", "-3"
                    ] + param_args + ["-o", tmp_msh_filepath, geo_filepath]

        result = sp.run(cmd_list, stdout=sp.PIPE, stderr=sp.PIPE)

        _check_mesh_conversion_result(result, tmp_msh_filepath)

        mesh_data = convert_3d_gmsh_msh_to_fenics_mesh(tmp_msh_filepath)

    return mesh_data


#############################################################################


def create_mesh_view(labelled_mesh: LabelledMesh,
                     domain_index: ty.Optional[int] = None) -> LabelledMesh:
    """[summary]
    
    Arguments:
        labelled_mesh {LabelledMesh} -- [description]
        domain_index {int} -- [description] If none return entire mesh as view
    
    Returns:
        LabelledMesh
    """

    if domain_index is None:
        marker = fn.MeshFunction("size_t", labelled_mesh.mesh,
                                 labelled_mesh.mesh.topology().dim(),
                                 0)  # mark entirety of the full mesh
        mesh_view = fn.MeshView.create(marker, 0)
    else:
        mesh_view = fn.MeshView.create(labelled_mesh.subdomain_mesh_func,
                                       domain_index)

    # make expression for original mesh

    V = fn.FunctionSpace(labelled_mesh.mesh, "DG", 0)
    u_smf = fn.Function(V)

    helper = np.asarray(labelled_mesh.subdomain_mesh_func.array(),
                        dtype=np.int32)  # type: ignore

    dm = V.dofmap()
    for cell in fn.cells(labelled_mesh.mesh):
        helper[dm.cell_dofs(cell.index())] = labelled_mesh.subdomain_mesh_func[
            cell]  # type: ignore

    u_smf.vector()[:] = helper

    # evaluate on mesh view

    view_smf = fn.MeshFunction("size_t", mesh_view,
                               mesh_view.topology().dim(), 0)
    for c in fn.cells(mesh_view):

        if mesh_view.topology().dim() == 2:
            cell_midpoint = (c.midpoint().x(), c.midpoint().y())
        elif mesh_view.topology().dim() == 3:
            cell_midpoint = (c.midpoint().x(), c.midpoint().y(),
                             c.midpoint().z())
        else:
            assert False, "Unexpected condition"

        view_smf[c] = int(u_smf(*cell_midpoint))  # type: ignore

    return LabelledMesh(mesh=mesh_view, subdomain_mesh_func=view_smf)


#############################################################################


def newton_solver_parameters():
    return {
        "nonlinear_solver": "newton",
        "newton_solver": {
            "linear_solver": "petsc",
            "relative_tolerance": 1e-2
        }
    }


#############################################################################
# DEPRECATED!!!
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
        function_filepath = os.path.join(PARAVIEW_TMP_FOLDERPATH,
                                         "function.pvd")

    function_filepath = os.path.abspath(function_filepath)

    if os.path.isfile(function_filepath):
        os.remove(function_filepath)

    result_file = fn.File(function_filepath)

    result_file << u

    print(f"Wrote function to: {function_filepath}")

    if open_file:
        sp.Popen(['paraview', function_filepath])

    return function_filepath


###########################################################
# Classes
###########################################################


class IsBoundary(fn.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary
