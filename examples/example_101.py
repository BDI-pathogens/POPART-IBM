#!/usr/bin/env python3
"""
Basic example for organising input parameters from the test parameter set and running PopART-IBM.  

Steps 1 and 2 set up the parameters from the test parameters.  
Steps 3 and 4 run the model (these can be done from the command-line too).  

Author: p-robot, 2021
"""

import sys, os, subprocess, shutil
from os.path import join

sys.path.append(join("..", "src", "popart_ibm"))
from parameters import ParameterSet

MODEL_DIR = ".."
IBM_DIR_TEST = join(MODEL_DIR, "src")
TEST_PARAM_DIR = join(MODEL_DIR, "tests", "data")
EXAMPLE_PARAM_DIR = join(MODEL_DIR, "examples", "PARAMS_COMMUNITY5")

NRUNS = 1
EXE = "popart-simul.exe"
COMMAND = f"{IBM_DIR_TEST}/{EXE} {EXAMPLE_PARAM_DIR} {NRUNS}"

parameter_transpose_file = join(TEST_PARAM_DIR, "parameters_transpose.csv")
fertility_file = join(TEST_PARAM_DIR, "parameters_fertility.csv")
mortality_file = join(TEST_PARAM_DIR, "parameters_mortality.csv")
patchinfo_file = join(TEST_PARAM_DIR, "patchinfo.txt")
python_seed_file = join(TEST_PARAM_DIR, "python_seed.txt")

# 1. Create a ParameterSet object and read the test data
p = ParameterSet()
p.read_test_parameters(
    parameter_transpose_file = parameter_transpose_file,
    fertility_file = fertility_file, 
    mortality_file = mortality_file, 
    patchinfo_file = patchinfo_file,
    python_seed_file = python_seed_file)

# 2. Write parameters in the form required for POPART-IBM
p.write_param_set(output_dir = EXAMPLE_PARAM_DIR)

# Create an "Output" directory within the parameter directory 
# (remove it if it already exists)
# (this is where model output will be stored)
shutil.rmtree(join(EXAMPLE_PARAM_DIR, "Output"), ignore_errors = True)
os.mkdir(join(EXAMPLE_PARAM_DIR, "Output"))

# 3. Construct the compilation command and compile
# (Note: this can be done from the terminal/console using
# the command 'cd src; make clean; make all)
compile_command = "make clean; make all"
completed_compilation = subprocess.run([compile_command], 
    shell = True, cwd = IBM_DIR_TEST, capture_output = True)

# 4. Call/run the model on the example parameter set
completed_run = subprocess.run([COMMAND], shell = True, capture_output = True)
