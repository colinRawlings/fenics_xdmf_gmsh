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

EXAMPLE_SCRIPTS_BLACKLIST = ["3d_MM_curl_genlin_app.py",]
SHARED_EXAMPLE_SCRIPTS_BLACKLIST = ["parallel_3d_calc.py"]

###########################################################
# Tests
###########################################################


@pytest.mark.parametrize('script', glob.glob(os.path.join(EXAMPLE_SCRIPTS_DIR, "*.py")))
def test_script_execution(script):

    assert (os.path.isdir(EXAMPLE_SCRIPTS_DIR))

    if os.path.split(script)[-1] in EXAMPLE_SCRIPTS_BLACKLIST:
        return

    print(f"testing: {script}")
    runpy.run_path(script)


@pytest.mark.parametrize('script', glob.glob(os.path.join(SHARED_EXAMPLE_SCRIPTS_DIR, "*.py")))
def test_shared_script_execution(script):

    assert (os.path.isdir(SHARED_EXAMPLE_SCRIPTS_DIR))

    if os.path.split(script)[-1] in SHARED_EXAMPLE_SCRIPTS_BLACKLIST:
        return

    print(f"testing: {script}")
    runpy.run_path(script)