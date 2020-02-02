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

EXAMPLE_SCRIPTS = glob.glob(os.path.join(os.path.dirname(__file__), os.pardir, "examples", "trial*.py"))


###########################################################
# Tests
###########################################################


@pytest.mark.parametrize('script', EXAMPLE_SCRIPTS)
def test_script_execution(script):
    runpy.run_path(script)

