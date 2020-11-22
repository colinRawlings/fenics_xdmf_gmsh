# !python3

###########################################################
# imports
###########################################################

import sys
import os
import subprocess as sp
import re
import glob

import pytest

###########################################################
# Defintions
###########################################################

EXAMPLE_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "examples")


###########################################################
# Tests
###########################################################

@pytest.mark.skipif(True, reason="parallel running not working with mixed dimensional?")
@pytest.mark.parametrize('script', glob.glob(os.path.join(EXAMPLE_SCRIPTS_DIR, "parallel*.py")))
def test_parallel_execution_with_mpi(script):

    NUM_THREADS = 2

    parallel_calc_script_filepath = os.path.join(os.path.dirname(__file__),
                                                 os.pardir, "examples",
                                                 script)

    assert os.path.isfile(parallel_calc_script_filepath)

    result = sp.run([
        "mpirun", "-np",
        str(NUM_THREADS), sys.executable, parallel_calc_script_filepath
    ],
                    stdout=sp.PIPE,
                    stderr=sp.PIPE)

    assert result.returncode == 0

    stdout_lines = result.stdout.decode().splitlines()

    reported_pool_size = None
    for line in stdout_lines:
        if line.startswith("mpi process"):
            match = re.search(r"size: (\d+)", line)
            assert match
            reported_pool_size=int(match[1])

        if reported_pool_size is not None:
            break
    
    assert reported_pool_size is not None
    assert reported_pool_size == NUM_THREADS