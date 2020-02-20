# !python

###########################################################
# Imports
###########################################################

import os

import pytest
import dolfin as fn

import fenics_utils as fu

###########################################################
# Definitions
###########################################################


###########################################################
# Fixtures
###########################################################

@pytest.fixture
def mesh_2d():
    return fn.UnitSquareMesh(10, 10)

@pytest.fixture
def mesh_3d():
    return fn.UnitCubeMesh(10, 10, 10)

@pytest.fixture
def cg2_2d_scalar_func(mesh_2d):
    V = fn.FunctionSpace(mesh_2d, 'CG', 2)
    u = fn.Function(V)
    return u

@pytest.fixture
def cg2_2d_vector_func(mesh_2d):
    V = fn.VectorFunctionSpace(mesh_2d, 'CG', 2)
    u = fn.Function(V)
    return u


###########################################################
# Tests
###########################################################

def test_detect_function_degree1(mesh_2d):

    V1 = fn.FunctionSpace(mesh_2d, 'CG', 1)
    u1 = fn.Function(V1)

    assert fu.detect_function_degree(u1) == 1

def test_detect_function_degree2(mesh_2d):

    V2 = fn.FunctionSpace(mesh_2d, 'CG', 2)
    u2 = fn.Function(V2)

    assert fu.detect_function_degree(u2) == 2

def test_detect_vector_function_degree2(mesh_3d):

    V2 = fn. VectorFunctionSpace(mesh_3d, 'CG', 2)
    u2 = fn.Function(V2)

    assert fu.detect_function_degree(u2) == 2

def test_refinement_degree2(cg2_2d_scalar_func):
    ur = fu.refine_and_project_to_cg1(cg2_2d_scalar_func)

def test_refinement_degree2(cg2_2d_vector_func):
    ur = fu.refine_and_project_to_cg1(cg2_2d_vector_func)