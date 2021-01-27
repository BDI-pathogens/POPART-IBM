#!/usr/bin/env python3
"""
Tests of PopART-IBM associated with timing using pytest. 

Several checks are performed by adjusting input parameters and comparing to expected outputs in
the `Annual_outputs*.csv` files.  A single simulation is checked (and only the intervention patch).

W. Probert, 2019
"""
import subprocess, shutil, os, pytest, sys, numpy as np, pandas as pd
from os.path import join

sys.path.append("./src/popart_ibm")
from parameters import ParameterSet, SamplingSchedule

from utils import adjust_macros
import constants as c


######################################################
# Override setup_popart_ibm_methods() in conftest.py 
# (since the setup/teardown protocol in setup_popart_ibm_methods() creates a new data directory for
# each method; here, we compile and run the model once at the level of the TestClass object)
# ----------------------------------------------------

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
    
    def fin():
        """
        After each method, remove the temporary data directory
        """
        shutil.rmtree(c.DATA_DIR_TEST, ignore_errors = True)
    request.addfinalizer(fin)


class TestClass(object):
    def test_late_hiv_introduction(self):
        """
        Set HIV seeding date to be 1990 and make sure there are no HIV positive individuals before
        """
        
        # Adjust HIV seeding time
        for patch in [0, 1]:

            param_filename = "param_processed_patch" + str(patch) + "_times.csv"

            df_times_params = pd.read_csv(join(c.DATA_DIR_TEST, param_filename), \
                delim_whitespace = True)

            df_times_params['start_time_hiv'] = 1990

            df_times_params.to_csv(join(c.DATA_DIR_TEST, param_filename), sep = " ", index = False)
        
        # Call the model
        completed_run = subprocess.run([c.command, c.DATA_DIR_TEST, str(c.NRUNS)], 
            capture_output = True)
        
        # Read the Annual_output* dataframes
        df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))

        # Check if there were any HIV cases
        np.testing.assert_equal(np.sum(df_annual.NumberPositive[df_annual.Year < 1990].values), 0)


    def test_start_end_simul(self):
        """
        Change start/end time of the simulation (simulation should end 1st Jan of end_time_simul)
        """
        
        start = 1910; end = 2022
        
        # Adjust start/end time of the simulation
        for patch in [0, 1]:

            param_filename = "param_processed_patch" + str(patch) + "_times.csv"

            df_times_params = pd.read_csv(join(c.DATA_DIR_TEST, param_filename), \
                delim_whitespace = True)

            df_times_params['start_time_simul'] = start
            df_times_params['end_time_simul'] = end

            df_times_params.to_csv(join(c.DATA_DIR_TEST, param_filename), sep = " ", index = False)
        
        # Call the model
        completed_run = subprocess.run([c.command, c.DATA_DIR_TEST, str(c.NRUNS)], 
            capture_output = True)
        
        # Read the Annual_output* dataframes
        df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))

        # Check if there were any HIV cases
        np.testing.assert_array_equal([df_annual.Year.min(), df_annual.Year.max()], [start, end-1])
