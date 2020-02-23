# !python3
"""
Shamelessly stolen from: SO56807698
"""

###########################################################
# Imports
###########################################################

import os
import glob
import runpy
import pytest
import sys

###########################################################
# Defintions
###########################################################

EXAMPLE_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "examples")
SHARED_EXAMPLE_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "examples")

###########################################################
# Tests
###########################################################


@pytest.mark.parametrize('script', glob.glob(os.path.join(EXAMPLE_SCRIPTS_DIR, "*.py")))
def test_script_execution(script):

    assert(os.path.isdir(EXAMPLE_SCRIPTS_DIR))

    print(f"testing: {script}")

    if script.endswith('app.py'):
        return

    runpy.run_path(script)

@pytest.mark.parametrize('script', glob.glob(os.path.join(SHARED_EXAMPLE_SCRIPTS_DIR, "*.py")))
def test_shared_script_execution(script):

    assert(os.path.isdir(SHARED_EXAMPLE_SCRIPTS_DIR))

    print(f"testing: {script}")
    runpy.run_path(script)