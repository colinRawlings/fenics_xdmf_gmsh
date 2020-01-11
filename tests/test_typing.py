# !python3

###########################################################
# imports
###########################################################

import os
import subprocess as sp

###########################################################
# Tests
###########################################################


def test_typing_with_pyright():

    config_filepath = os.path.join(os.path.dirname(__file__), os.pardir, "pyrightconfig.json")
    assert os.path.isfile(config_filepath)

    result = sp.run(["which", "pyright"])
    assert result.returncode == 0, "No pyright found"

    result = sp.run(["pyright", "-p", config_filepath], stdout=sp.PIPE, stderr=sp.PIPE)

    msg = f"typing check failed with:\n{result.stdout.decode()}\n{result.stderr.decode()}"

    assert result.returncode == 0, msg
