#!/usr/bin/env python3
"""
conftest.py for testing the PopART-IBM

W. Probert, 2020
"""

import subprocess, shutil, os
from os.path import join
import numpy as np, pandas as pd, sys
import pytest 

sys.path.append("./src/popart_ibm")
from parameters import ParameterSet, SamplingSchedule

from utils import adjust_macros
import constants as c

###################################
# Add --runslow command; by default
# suppresses running slow tests, 
# those marked with the decorator
# @pytest.mark.slow.  --runslow will
# run all tests.  
# ---------------------------------

def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )

def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


##########################################
# Setup/teardown functions defined
# at the test-class and test-method levels
# ----------------------------------------

# class-level setup/teardown
@pytest.fixture(scope = "class", autouse = True)
def compile_popart_ibm(request):
    """
    Set up function: call the IBM so as to generate test output of Annual_outputs files
    """
    print("Creating temporary code directory")
    # Make a temporary copy of the code (remove this temporary directory if it already exists)
    shutil.rmtree(c.IBM_DIR_TEST, ignore_errors = True)
    shutil.copytree(c.IBM_DIR, c.IBM_DIR_TEST)

    # Turn on the macros for the output files above
    macros_one = [
        "PRINT_EACH_RUN_OUTPUT", 
        "PRINT_ALL_RUNS", 
        "WRITE_EVERYTIMESTEP"]

    macro_file = join(c.IBM_DIR_TEST, "constants.h")
    [adjust_macros(macro_file, macro, "1", macro_file) for macro in macros_one]

    # Construct the compilation command and compile
    compile_command = "make clean; make all"
    completed_compilation = subprocess.run([compile_command], 
        shell = True, cwd = c.IBM_DIR_TEST, capture_output = True)
    
    def fin():
        """Remove the temporary code directory"""
        shutil.rmtree(c.IBM_DIR_TEST, ignore_errors = True)
    request.addfinalizer(fin)


# Method ("function") setup/teardown
@pytest.fixture(scope = "function", autouse = True)
def setup_popart_ibm_methods(request):
    """
    Before each method, create a temporary data directory
    """
    print("Creating temporary data directory")
    # Make a temporary data directory, copy parameters to this dir
    #shutil.copytree(c.INPUTDIR, c.DATA_DIR_TEST)
    p = ParameterSet()
    p.read_test_parameters(
        parameter_transpose_file = c.parameter_transpose_file,
        fertility_file = c.fertility_file, 
        mortality_file = c.mortality_file, 
        patchinfo_file = c.patchinfo_file, 
        python_seed_file = c.python_seed_file)
    
    p.write_param_set(output_dir = c.DATA_DIR_TEST)
    
    # Remove the Output folder and recreate it
    shutil.rmtree(join(c.DATA_DIR_TEST, "Output"), ignore_errors = True)
    os.mkdir(join(c.DATA_DIR_TEST, "Output"))
    
    # Call the model
    completed_run = subprocess.run([c.command, c.DATA_DIR_TEST, str(c.NRUNS)], 
        capture_output = True)
    
    def fin():
        """
        After each method, remove the temporary data directory
        """
        shutil.rmtree(c.DATA_DIR_TEST, ignore_errors = True)
    request.addfinalizer(fin)


@pytest.fixture(scope = "function")
def basic_model_output():
    df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))
    df_timestep = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_timestep))
    return df_annual, df_timestep
