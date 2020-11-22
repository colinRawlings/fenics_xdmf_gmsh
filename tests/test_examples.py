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

EXAMPLE_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "examples")

EXAMPLE_SCRIPTS_BLACKLIST = ["3d_MM_curl_genlin_app.py"]
EXAMPLE_SCRIPTS_WHITELIST = ["calc_system_matrices.py"]

###########################################################
# Tests
###########################################################


@pytest.mark.parametrize("script", glob.glob(os.path.join(EXAMPLE_SCRIPTS_DIR, "*.py")))
def test_script_execution(script):

    assert os.path.isdir(EXAMPLE_SCRIPTS_DIR)

    script_name = os.path.split(script)[-1]
    if script_name in EXAMPLE_SCRIPTS_BLACKLIST:
        return

    print(f"testing: {script}")
    runpy.run_path(script)
