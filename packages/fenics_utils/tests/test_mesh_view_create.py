# !python3

###########################################################
# imports
###########################################################

import os

import pkg_resources as rc

import dolfin as fn

import fenics_utils as fu

###########################################################
# Definitions
###########################################################

GEO_DIR = rc.resource_filename("fenics_utils", "assets")  # type: ignore

###########################################################
# Test
###########################################################
