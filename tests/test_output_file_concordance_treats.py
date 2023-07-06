#!/usr/bin/env python3
"""
Test internal validity of the annual outputs file using pytest. 

Several checks are performed to make sure incidence, prevalence, and other measures have internal
validity within the `Annual_outputs*.csv` files.  A single file is checked.  

W. Probert, 2018
"""

import subprocess, shutil, os, numpy as np, pandas as pd, pytest, sys
from os.path import join

sys.path.append("./src/popart_ibm")
from parameters import ParameterSet, SamplingSchedule

from utils import adjust_macros
import constants as c


@pytest.fixture
def treats_model_output():
    df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))
    df_timestep = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_timestep))
    df_cost_effectiveness = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_cost_effectiveness))
    df_treats = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_treats))
    
    return df_annual, df_timestep, df_cost_effectiveness, df_treats


######################################################
# Override compile_popart_ibm() in conftest.py 
# (since the setup/teardown protocol in setup_popart_ibm_methods() creates a new data directory for
# each method; here, we compile and run the model once at the level of the TestClass object)
# ----------------------------------------------------

@pytest.fixture(scope = "class")
def compile_popart_ibm(request):
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
        "WRITE_COST_EFFECTIVENESS_OUTPUT", 
        "PRINT_EACH_RUN_OUTPUT", 
        "PRINT_ALL_RUNS", 
        "WRITE_EVERYTIMESTEP",
        "WRITE_TREATS_OUTPUT"]
    
    macro_file = join(c.IBM_DIR_TEST, "constants.h")
    [adjust_macros(macro_file, macro, "1", macro_file) for macro in macros_one]
    

    # Construct the compilation command and compile
    compile_command = "make clean; make all"
    completed_compilation = subprocess.run([compile_command], 
        shell = True, cwd = c.IBM_DIR_TEST, capture_output = True)
    
    # Construct the executable command
    command = join(c.IBM_DIR_TEST, c.EXE)
    
    # Call the model
    completed_run = subprocess.run([command, c.DATA_DIR_TEST, str(c.NRUNS)], capture_output = True)
    
    def fin():
        """Remove the temporary code directory"""
        shutil.rmtree(c.IBM_DIR_TEST, ignore_errors = True)
        shutil.rmtree(c.DATA_DIR_TEST, ignore_errors = True)
    request.addfinalizer(fin)


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


class TestClass(object):
    
    ###################################################################################
    ############### Consistency of Annual_ and TREATS outputs #########################
    ###################################################################################
    
    def test_treats_number_overall(self, treats_model_output):
        # Check that the total population is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if "Ntot" in c]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual['TotalPopulation']

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_men_overall(self, treats_model_output):
        # Check that the total population of men is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if "Ntot_M" in c]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual['PopulationM']

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_women_overall(self, treats_model_output):
        # Check that the total population of women is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if "Ntot_F" in c]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual['PopulationF']

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_men_13_30(self, treats_model_output):
        # Check that the total population of men aged 13-30 is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if ("Ntot_M_13" in c) or
            ("Ntot_M_15" in c) or ("Ntot_M_20" in c) or ("Ntot_M_25" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual[['NMage13-18', 'NMage18-23', 'NMage23-30']].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_women_13_30(self, treats_model_output):
        # Check that the total population of women aged 13-30 is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if ("Ntot_F_13" in c) or
            ("Ntot_F_15" in c) or ("Ntot_F_20" in c) or ("Ntot_F_25" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual[['NFage13-18', 'NFage18-23', 'NFage23-30']].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)
    
    
    def test_treats_number_men_30_60(self, treats_model_output):
        # Check that the total population of men aged 30-60 is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if ("Ntot_M_30" in c) or
            ("Ntot_M_35" in c) or ("Ntot_M_40" in c) or ("Ntot_M_45" in c) or
            ("Ntot_M_50" in c) or ("Ntot_M_55" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual[['NMage30-40', 'NMage40-50', 'NMage50-60']].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_women_30_60(self, treats_model_output):
        # Check that the total population of women aged 30-60 is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if ("Ntot_F_30" in c) or
            ("Ntot_F_35" in c) or ("Ntot_F_40" in c) or ("Ntot_F_45" in c) or
            ("Ntot_F_50" in c) or ("Ntot_F_55" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual[['NFage30-40', 'NFage40-50', 'NFage50-60']].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_men_60_80(self, treats_model_output):
        # Check that the total population of men aged 60+ is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if ("Ntot_M_60" in c) or
            ("Ntot_M_65" in c) or ("Ntot_M_70" in c) or ("Ntot_M_75" in c) or ("Ntot_M_80" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual[['NMage60-80']].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_women_60_80(self, treats_model_output):
        # Check that the total population of women aged 60+ is consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if ("Ntot_F_60" in c) or
            ("Ntot_F_65" in c) or ("Ntot_F_70" in c) or ("Ntot_F_75" in c) or ("Ntot_F_80" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual[['NFage60-80']].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_positive_overall_cd4(self, treats_model_output):
        # Check that the number positive using CD4 category is correct
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if (("CD41" in c) or ("CD42" in c) or
            ("CD43" in c) or ("CD44" in c)) and ("Ntot" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual['NumberPositive']

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_positive_overall_spvl(self, treats_model_output):
        # Check that the number positive using SVPL category is correct
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if (("SPVL1" in c) or ("SPVL2" in c) or
            ("SPVL3" in c) or ("SPVL4" in c)) and ("Ntot" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual['NumberPositive']

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_positive_men_cd4(self, treats_model_output):
        # Check that the number positive men using CD4 category is correct
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if (("CD41" in c) or ("CD42" in c) or
            ("CD43" in c) or ("CD44" in c)) and ("Ntot_M" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual['NumberPositiveM']

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_number_positive_women_cd4(self, treats_model_output):
        # Check that the number positive women using CD4 category is correct
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if (("CD41" in c) or ("CD42" in c) or
            ("CD43" in c) or ("CD44" in c)) and ("Ntot_F" in c)]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual['NumberPositiveF']

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_incident_cases_against_NewCasesThisYear(self, treats_model_output):
        # Check that incident cases are consistent with NewCasesThisYear (which includes seed cases)
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if "Ninc" in c]
        treats = df_treats[cols].sum(axis = 1)
        annual = df_annual['NewCasesThisYear']

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_incident_cases_men(self, treats_model_output):
        # Check that incident cases in men are consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if "Ninc_M" in c]
        treats = df_treats[cols].sum(axis = 1)
        cols = [c for c in df_annual.columns if "IncMage" in c]
        annual = df_annual[cols].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_incident_cases_women(self, treats_model_output):
        # Check that incident cases in women are consistent
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if "Ninc_F" in c]
        treats = df_treats[cols].sum(axis = 1)
        cols = [c for c in df_annual.columns if "IncFage" in c]
        annual = df_annual[cols].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)

    
    def test_treats_died_from_HIV(self, treats_model_output):
        # Check that deaths from HIV-related causes, when stratified by gender, age, CD4 count, and 
        # SPVL, in TREATS file are the same as what is reported in Annual outputs file
        # NB: Annual outputs file is cumulative deaths
        
        df_annual, df_timestep, df_cost_effectiveness, df_treats = treats_model_output
        
        cols = [c for c in df_treats.columns if "NDied_from_HIV_" in c]
        treats = np.cumsum(df_treats[cols].sum(axis = 1))
        cols = [c for c in df_annual.columns if "NDied_from_HIV" == c]
        annual = df_annual[cols].sum(axis = 1)

        np.testing.assert_array_equal(treats.values, annual.values)

