#!/usr/bin/env python3
"""
Tests of PopART-IBM associated with male circumcision using pytest. 

Several checks are performed by adjusting input parameters and comparing to expected outputs in
the `Annual_outputs*.csv` files.  A single simulation is checked (and only the intervention patch).

W. Probert, 2019
"""

import subprocess, shutil, os, pytest, sys
from os.path import join
import numpy as np, pandas as pd

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
    def test_tmc_prevalence(self):
        """
        TMC prevalence should be roughly (within 0.05) the same as the overall prevalence if
        there's no VMMC.  
        """
        
        TMC_PREVALENCE = 0.7
        TOLERANCE = 0.05 
        
        # Adjust TMC prevalence to be TMC_PREVALENCE
        for patch in [0, 1]:

            param_filename = "param_processed_patch" + str(patch) + "_HIV.csv"

            df_hiv_params = pd.read_csv(join(c.DATA_DIR_TEST, param_filename), \
                delim_whitespace = True, engine = 'python')

            df_hiv_params['p_child_circ'] = TMC_PREVALENCE

            df_hiv_params.to_csv(join(c.DATA_DIR_TEST, param_filename), sep = " ", index = False)
        
        
        # Adjust VMMC acceptance probability to be zero
        for patch in [0, 1]:

            param_filename = "param_processed_patch" + str(patch) + "_cascade.csv"

            df_cascade_params = pd.read_csv(join(c.DATA_DIR_TEST, param_filename), \
                delim_whitespace = True, engine = 'python')

            df_cascade_params['p_circ_popart'] = 0.0
            df_cascade_params['p_circ_nopopart'] = 0.0

            df_cascade_params.to_csv(join(c.DATA_DIR_TEST, param_filename), 
                sep = " ", index = False)

        # Call the model
        completed_run = subprocess.run([c.command, c.DATA_DIR_TEST, str(c.NRUNS)], 
            capture_output = True)

        # Read the Annual_output* dataframes
        df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))

        # Check if any of these are less than 0.1 from TMC_PREVALENCE
        assertion = any((df_annual.PropMenCirc - TMC_PREVALENCE) > TOLERANCE)
        np.testing.assert_equal(assertion, False)

    def test_no_vmmc_no_tmc(self):
        """
        VMMC prevalence to zero, TMC prevalence to zero (no proportion of circumcision)
        should mean there is no circumcised men in the population (PropMenCirc in Annual outputs
        should sums zero)
        """

        # Adjust TMC prevalence to be zero
        for patch in [0, 1]:

            param_filename = "param_processed_patch" + str(patch) + "_HIV.csv"

            df_hiv_params = pd.read_csv(join(c.DATA_DIR_TEST, param_filename), \
                delim_whitespace = True, engine = 'python')

            df_hiv_params['p_child_circ'] = 0.0

            df_hiv_params.to_csv(join(c.DATA_DIR_TEST, param_filename), sep = " ", index = False)


        # Adjust VMMC acceptance probability to be zero
        for patch in [0, 1]:

            param_filename = "param_processed_patch" + str(patch) + "_cascade.csv"

            df_cascade_params = pd.read_csv(join(c.DATA_DIR_TEST, param_filename), \
                delim_whitespace = True, engine = 'python')

            df_cascade_params['p_circ_popart'] = 0.0
            df_cascade_params['p_circ_nopopart'] = 0.0

            df_cascade_params.to_csv(join(c.DATA_DIR_TEST, param_filename), 
                sep = " ", index = False)

        # Call the model
        completed_run = subprocess.run([c.command, c.DATA_DIR_TEST, str(c.NRUNS)], 
            capture_output = True)

        # Read the Annual_output* dataframes
        df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))

        np.testing.assert_equal(np.sum(df_annual.PropMenCirc.values), 0)

    def test_tmc_eff_1(self):
        """
        TMC effectiveness of 1 and TMC prevalence of 100% should mean there is no HIV in the 
        male population after the seeding years (so incidence of zero).  
        """
        
        # Adjust TMC effectiveness to be 1 and TMC prevalence to be 100%
        for patch in [0, 1]:
            
            param_filename = "param_processed_patch" + str(patch) + "_HIV.csv"
            
            df_hiv_params = pd.read_csv(join(c.DATA_DIR_TEST, param_filename), \
                delim_whitespace = True, engine = 'python')
            
            df_hiv_params['p_child_circ'] = 1.0
            df_hiv_params['eff_circ_tmc'] = 1.0
            df_hiv_params['eff_circ_vmmc'] = 1.0 # need to change for internal model checks
            
            df_hiv_params.to_csv(join(c.DATA_DIR_TEST, param_filename), sep = " ", index = False)
        
        # Call the model
        completed_run = subprocess.run([c.command, c.DATA_DIR_TEST, str(c.NRUNS)], 
            capture_output = True)
        print(completed_run)
        # Read the Annual_output* dataframes
        df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))
        
        cols = [c for c in df_annual.columns if "IncM" in c]
        incidence_by_age = df_annual[df_annual.Year > 2000][cols]
        
        np.testing.assert_equal(incidence_by_age.to_numpy().sum(), 0)
