# !python3

###########################################################
# Imports
###########################################################

import sys
import os
import subprocess as sp
import re

###########################################################
# Tests
###########################################################


def test_parallel_execution_with_mpi():

    NUM_THREADS = 2

    parallel_calc_script_filepath = os.path.join(os.path.dirname(__file__),
                                                 os.pardir, "examples",
                                                 "parallel_3d_calc.py")

    assert os.path.isfile(parallel_calc_script_filepath)

    result = sp.run([
        "mpirun", "-np",
        str(NUM_THREADS), sys.executable, parallel_calc_script_filepath
    ],
                    stdout=sp.PIPE,
                    stderr=sp.PIPE)

    assert result.returncode == 0

    stdout_lines = result.stdout.decode().splitlines()

    size_found = False
    for line in stdout_lines:
        if line.startswith("mpi process"):
            match = re.search(r"size: (\d+)", line)
            assert match
            assert match[1] == str(NUM_THREADS)
            size_found = True

        if size_found:
            break
