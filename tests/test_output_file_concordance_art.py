#!/usr/bin/env python3
"""
Test concordance of timestep outputs file and ART_status file

W. Probert, 2018
"""
import subprocess, shutil, os, pytest, sys
from os.path import join
import numpy as np, pandas as pd

sys.path.append("./src/popart_ibm")
from parameters import ParameterSet, SamplingSchedule

from utils import adjust_macros
import constants as c

@pytest.fixture
def extended_model_output():
    df_timestep = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_timestep))
    df_art_status = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_art_status))
    
    return df_timestep, df_art_status

######################################################
# Override setup_popart_ibm_methods() in conftest.py 
# (since the setup/teardown protocol in setup_popart_ibm_methods() creates a new data directory for
# each method; here, we compile and run the model once at the level of the TestClass object)
# ----------------------------------------------------

@pytest.fixture(scope = "function")
def setup_popart_ibm_methods(request):
    """
    Called before each method is run
    """
    pass
    def fin():
        """
        At the end of each method (test)
        """
        pass
    request.addfinalizer(fin)


# "Session" (all files in folder) setup/teardown
@pytest.fixture(scope = "class", autouse = True)
def compile_popart_ibm(request):
    """
    Set up function: call the IBM so as to generate test output of 
    Annual_outputs, Timestep_output, and Cost_effectiveness_output files
    """
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
    
    # Make a temporary copy of the code (remove this temporary directory if it already exists)
    shutil.rmtree(c.IBM_DIR_TEST, ignore_errors = True)
    shutil.copytree(c.IBM_DIR, c.IBM_DIR_TEST)
    
    # Turn on the macros for the output files above
    macros_one = [
        "WRITE_CALIBRATION", 
        "PRINT_EACH_RUN_OUTPUT", 
        "PRINT_ALL_RUNS", 
        "WRITE_EVERYTIMESTEP",
        "WRITE_ART_STATUS_BY_AGE_SEX"]
    
    macro_file = join(c.IBM_DIR_TEST, "constants.h")
    [adjust_macros(macro_file, macro, "1", macro_file) for macro in macros_one]
    

    # Construct the compilation command and compile
    compile_command = "make clean; make all"
    completed_compilation = subprocess.run([compile_command], 
        shell = True, cwd = c.IBM_DIR_TEST, capture_output = True)
    
    # Call the model
    completed_run = subprocess.run([c.command, c.DATA_DIR_TEST, str(c.NRUNS)], 
        capture_output = True)
    
    def fin():
        """
        At the end of each method (test)
        """
        shutil.rmtree(c.IBM_DIR_TEST, ignore_errors = True)
        shutil.rmtree(c.DATA_DIR_TEST, ignore_errors = True)
    request.addfinalizer(fin)


class TestClass(object):
    ############################################################################
    ################ Consistency of ART_status_age_sex #########################
    ############################################################################
    @pytest.mark.slow
    def test_art_status_number_men_overall(self, extended_model_output):
        # Check that the total population of men is consistent with timestep output
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "NM" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['N_m']

        np.testing.assert_array_equal(total.values, total_timestep.values)
    
    @pytest.mark.slow
    def test_art_status_number_women_overall(self, extended_model_output):
        # Check that the total population of women is consistent with timestep output
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "NF" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['N_f']

        np.testing.assert_array_equal(total.values, total_timestep.values)
    
    @pytest.mark.slow
    def test_art_status_number_positive_men(self, extended_model_output):
        # Check that the total number of HIV+ men is consistent with timestep outputs
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "HIV0" not in c if "NM" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['NPos_m']

        np.testing.assert_array_equal(total.values, total_timestep.values)
    
    @pytest.mark.slow
    def test_art_status_number_positive_women(self, extended_model_output):
        # Check that the total number of HIV+ women is consistent with timestep outputs
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "HIV0" not in c if "NF" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['NPos_f']

        np.testing.assert_array_equal(total.values, total_timestep.values)

    @pytest.mark.slow
    def test_vs_men(self, extended_model_output):
        # Check the total number of virally suppressed men is consistent with timestep outputs
        # Long-term viral suppression is coded as ART_status 2 (see constants.h)
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "ART2" in c if "NM" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['NVS_m']

        np.testing.assert_array_equal(total.values, total_timestep.values)

    @pytest.mark.slow
    def test_vs_women(self, extended_model_output):
        # Check the total number of virally suppressed women is consistent with timestep outputs
        # Long-term viral suppression is coded as ART_status 2 (see constants.h)
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "ART2" in c if "NF" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['NVS_f']

        np.testing.assert_array_equal(total.values, total_timestep.values)
    
    @pytest.mark.slow
    def test_art_men(self, extended_model_output):
        # Check the total number of men on ART is consistent with timestep outputs
        # "On ART" is coded as ART_status 1, 2, or 3 (see constants.h)
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "ART1" in c or "ART2" in c or "ART3" in c if "NM" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['NART_m']

        np.testing.assert_array_equal(total.values, total_timestep.values)

    @pytest.mark.slow
    def test_art_women(self, extended_model_output):
        # Check the total number of women on ART is consistent with timestep outputs
        # "On ART" is coded as ART_status 1, 2, or 3 (see constants.h)
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "ART1" in c or "ART2" in c or "ART3" in c if "NF" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['NART_f']

        np.testing.assert_array_equal(total.values, total_timestep.values)

    @pytest.mark.slow
    def test_knowledge_men(self, extended_model_output):
        # Check the total number of men who know their status is consistent with timestep outputs
        # "On ART" is coded as ART_status 1, 2, or 3 (see constants.h)
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "ART-1" not in c if "NM" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['N_knowpos_m']

        np.testing.assert_array_equal(total.values, total_timestep.values)
    
    @pytest.mark.slow
    def test_knowledge_women(self, extended_model_output):
        # Check the total number of women who know their status is consistent with timestep outputs
        # "On ART" is coded as ART_status 1, 2, or 3 (see constants.h)
        
        df_timestep, df_art_status = extended_model_output
        
        cols = [c for c in df_art_status.columns if "ART-1" not in c if "NF" in c]
        
        total = df_art_status[cols].sum(axis = 1)
        total_timestep = df_timestep['N_knowpos_f']

        np.testing.assert_array_equal(total.values, total_timestep.values)
        
    # if ((individual->ART_status==ARTNAIVE) || (individual->ART_status==EARLYART) || (individual->ART_status==LTART_VS) || (individual->ART_status==LTART_VU) || (individual->ART_status==ARTDROPOUT) || (individual->ART_status==CASCADEDROPOUT))

