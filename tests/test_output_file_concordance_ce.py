#!/usr/bin/env python3
"""
Test concordance between Annual_outputs* file and cost_effectiveness_* file

W. Probert, 2018
"""
import subprocess, shutil, os, pytest, sys, numpy as np, pandas as pd
from os.path import join

sys.path.append("./src/popart_ibm")
from parameters import ParameterSet, SamplingSchedule

from utils import adjust_macros
import constants as c


@pytest.fixture
def ce_model_output():
    df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))
    df_cost_effectiveness = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_cost_effectiveness))
    
    return df_annual, df_cost_effectiveness

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
        "WRITE_COST_EFFECTIVENESS_OUTPUT", 
        "PRINT_EACH_RUN_OUTPUT", 
        "PRINT_ALL_RUNS", 
        "WRITE_EVERYTIMESTEP"]
    
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
    def test_ce_Incidence(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['Incidence']
        ce = df_cost_effectiveness['Incidence']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_Prevalence(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['Prevalence']
        ce = df_cost_effectiveness['Prevalence']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_NumberPositive(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['NumberPositive']
        ce = df_cost_effectiveness['NumberPositive']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_NAnnual(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['NAnnual']
        ce = df_cost_effectiveness['NAnnual']
        np.testing.assert_array_equal(annual.values, ce.values)

    def test_ce_TotalPopulation(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['TotalPopulation']
        ce = df_cost_effectiveness['TotalPopulation']
        np.testing.assert_array_equal(annual.values, ce.values)
        
    def test_ce_NumberPositiveM(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['NumberPositiveM']
        ce = df_cost_effectiveness['NumberPositiveM']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_PopulationM(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['PopulationM']
        ce = df_cost_effectiveness['PopulationM']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_NumberPositiveF(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['NumberPositiveF']
        ce = df_cost_effectiveness['NumberPositiveF']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_PopulationF(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['PopulationF']
        ce = df_cost_effectiveness['PopulationF']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_AnnualNonPopartHIVtests(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = np.diff(df_annual['CumulativeNonPopartHIVtests'])
        ce = df_cost_effectiveness['AnnualNonPopartHIVtests']
        np.testing.assert_array_equal(annual, ce.values[1:])
    
    def test_ce_AnnualPopartHIVtests(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = np.diff(df_annual['CumulativePopartHIVtests'])
        ce = df_cost_effectiveness['AnnualPopartHIVtests']
        np.testing.assert_array_equal(annual, ce.values[1:])
    
    def test_ce_AnnualPopartHIVtests_split(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        ce_total = df_cost_effectiveness['AnnualPopartHIVtests']
        ce_positive = df_cost_effectiveness['AnnualPopartHIVtests_positive']
        ce_negative = df_cost_effectiveness['AnnualPopartHIVtests_negative']
        np.testing.assert_array_equal(ce_total.values, ce_positive.values + ce_negative.values)
    
    def test_ce_AnnualNonPopartCD4tests(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = np.diff(df_annual['CumulativeNonPopartCD4tests'])
        ce = df_cost_effectiveness['AnnualNonPopartCD4tests']
        np.testing.assert_array_equal(annual, ce.values[1:])
    
    def test_ce_AnnualPopartCD4tests(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = np.diff(df_annual['CumulativePopartCD4tests'])
        ce = df_cost_effectiveness['AnnualPopartCD4tests']
        np.testing.assert_array_equal(annual, ce.values[1:])
    
    def test_ce_NOnARTM(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['NOnARTM']
        ce = df_cost_effectiveness['NOnARTM']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_NOnARTF(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['NOnARTF']
        ce = df_cost_effectiveness['NOnARTF']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_NDied_from_HIV(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        annual = df_annual['NDied_from_HIV']
        ce = df_cost_effectiveness['NDied_from_HIV']
        np.testing.assert_array_equal(annual.values, ce.values)
    
    def test_ce_NDied_from_HIV_age_sex(self, ce_model_output):
        # Check that the number of individuals that died from HIV by age and sex equals the total
        df_annual, df_cost_effectiveness = ce_model_output
        
        colsf = [c for c in df_cost_effectiveness.columns if "NDied_from_HIVF" in c]
        colsm = [c for c in df_cost_effectiveness.columns if "NDied_from_HIVM" in c]
        
        ce_age_sex = df_cost_effectiveness[colsf].sum(axis = 1).cumsum() + \
            df_cost_effectiveness[colsm].sum(axis = 1).cumsum()
        ce = df_cost_effectiveness['NDied_from_HIV']
        
        np.testing.assert_array_equal(ce_age_sex.values, ce.values)

    def test_ce_NumberPositiveM_CD4(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        cols = ['NumberPositiveM_CD4_1', 'NumberPositiveM_CD4_2', 
            'NumberPositiveM_CD4_3', 'NumberPositiveM_CD4_4']
        
        ce_cd4 = df_cost_effectiveness[cols].sum(axis = 1)
        ce = df_cost_effectiveness['NumberPositiveM']
        np.testing.assert_array_equal(ce_cd4.values, ce.values)
    
    def test_ce_NumberPositiveF_CD4(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        cols = ['NumberPositiveF_CD4_1', 'NumberPositiveF_CD4_2', 
            'NumberPositiveF_CD4_3', 'NumberPositiveF_CD4_4']
        
        ce_cd4 = df_cost_effectiveness[cols].sum(axis = 1)
        ce = df_cost_effectiveness['NumberPositiveF']
        np.testing.assert_array_equal(ce_cd4.values, ce.values)
    
    def test_ce_NumberPositiveOnARTM_CD4(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        cols = ['NumberPositiveOnARTM_CD4_1', 'NumberPositiveOnARTM_CD4_2',
            'NumberPositiveOnARTM_CD4_3', 'NumberPositiveOnARTM_CD4_4']
        
        ce_cd4 = df_cost_effectiveness[cols].sum(axis = 1)
        ce = df_cost_effectiveness['NOnARTM']
        np.testing.assert_array_equal(ce_cd4.values, ce.values)
    
    def test_ce_NumberPositiveOnARTF_CD4(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        cols = ['NumberPositiveOnARTF_CD4_1', 'NumberPositiveOnARTF_CD4_2',
            'NumberPositiveOnARTF_CD4_3', 'NumberPositiveOnARTF_CD4_4']
        
        ce_cd4 = df_cost_effectiveness[cols].sum(axis = 1)
        ce = df_cost_effectiveness['NOnARTF']
        np.testing.assert_array_equal(ce_cd4.values, ce.values)
    
    def test_ce_NumberPositiveNotOnARTM_CD4(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        cols = ['NumberPositiveNotOnARTM_CD4_1', 'NumberPositiveNotOnARTM_CD4_2',
            'NumberPositiveNotOnARTM_CD4_3', 'NumberPositiveNotOnARTM_CD4_4']
        
        ce_cd4 = df_cost_effectiveness[cols].sum(axis = 1)
        ce = df_cost_effectiveness['NumberPositiveM'] - df_cost_effectiveness['NOnARTM']
        np.testing.assert_array_equal(ce_cd4.values[-20:], ce.values[-20:])
    
    def test_ce_NumberPositiveNotOnARTF_CD4(self, ce_model_output):
        df_annual, df_cost_effectiveness = ce_model_output
        cols = ['NumberPositiveNotOnARTF_CD4_1', 'NumberPositiveNotOnARTF_CD4_2',
            'NumberPositiveNotOnARTF_CD4_3', 'NumberPositiveNotOnARTF_CD4_4']
        
        ce_cd4 = df_cost_effectiveness[cols].sum(axis = 1)
        ce = df_cost_effectiveness['NumberPositiveF'] - df_cost_effectiveness['NOnARTF']
        np.testing.assert_array_equal(ce_cd4.values[-20:], ce.values[-20:])

