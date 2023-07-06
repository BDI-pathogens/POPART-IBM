#!/usr/bin/env python3
"""
Test internal validity of the annual outputs file, and concordance with timestep and transmission
tree files (made up of the individual and transmission file).  

Several checks are performed to make sure incidence, prevalence, and other measures have internal
validity within the `Annual_outputs*.csv` files.  A single file is checked.  

W. Probert, 2018
"""
import subprocess, shutil, os, pytest, sys, numpy as np, pandas as pd
from os.path import join

sys.path.append("./src/popart_ibm")
from parameters import ParameterSet, SamplingSchedule

from utils import adjust_macros
import constants as c

# Range of years over which we want to ensure consistency across files (arbitrary)
years = range(1976, 2020)

@pytest.fixture
def extended_model_output():
    df_annual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_annual))
    df_timestep = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_timestep))
    df_individual = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_individual))
    df_transmission = pd.read_csv(join(c.DATA_DIR_TEST, "Output", c.f_transmission))
    
    return df_annual, df_timestep, df_individual, df_transmission

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
        "WRITE_PHYLOGENETICS_OUTPUT"]
    
    macro_file = join(c.IBM_DIR_TEST, "constants.h")
    [adjust_macros(macro_file, macro, "1", macro_file) for macro in macros_one]
    
    # Generate output for phylo patch 0 (the "inside"/PopART patch)
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
    ####################################################################################
    ###################### Internal consistency of Annual_outputs ######################
    ####################################################################################
    
    def test_incidence(self, extended_model_output):
        # Check that incidence is calculated correctly
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        incidence = df_annual['NAnnual']/(df_annual['TotalPopulation'] - \
            df_annual['NumberPositive'])
        
        np.testing.assert_allclose(incidence.values, \
            df_annual['Incidence'].values, atol = 1e-5)
    
    def test_hiv_positive_population_by_sex(self, extended_model_output):
        # Check that total positive population size is sum of male and female parts
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        population = df_annual['NumberPositiveM'] + df_annual['NumberPositiveF']
        
        np.testing.assert_array_equal(population.values, df_annual['NumberPositive'].values)
    
    def test_number_positive_by_new_cases(self, extended_model_output):
        # Check that number positive is sum of incident cases
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        numberpositive = (df_annual['NewCasesThisYear'].cumsum() - \
            df_annual['NHIV_pos_dead'])
        
        np.testing.assert_array_equal(numberpositive.values, 
            df_annual['NumberPositive'].values)
    
    def test_total_population_by_sex(self, extended_model_output):
        # Check that total population size is sum of male and female populations
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        population = df_annual['PopulationM'] + df_annual['PopulationF']
        np.testing.assert_array_equal(population.values, df_annual['TotalPopulation'].values)
    
    def test_prevalence(self, extended_model_output):
        # Check that prevalance is calculated correctly
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        prevalence = df_annual['NumberPositive']/df_annual['TotalPopulation']
        np.testing.assert_allclose(prevalence.values, df_annual['Prevalence'].values, 
            atol = 1e-5)
    
    def test_new_cases_by_age(self, extended_model_output):
        # Check that new cases this year is sum across all age groups
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        inc_cols = [c for c in df_annual.columns if ("Inc" in c) & ("age" in c)]
        newcases = df_annual[inc_cols].sum(axis = 1)
        np.testing.assert_array_equal(newcases.values, df_annual['NewCasesThisYear'].values)
    
    def test_number_positive_by_age(self, extended_model_output):
        # Check that number positive is sum across all age groups
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        npos_cols = [c for c in df_annual.columns if ("NPos" in c) & ("age" in c)]
        numberpos = df_annual[npos_cols].sum(axis = 1)
        np.testing.assert_array_equal(numberpos.values, df_annual['NumberPositive'].values)
    
    def test_number_positive_men_by_age(self, extended_model_output):
        # Check number positive men is same across all age groups
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        ntotm_cols = [c for c in df_annual.columns if ("NMage" in c)]
        ntotm = df_annual[ntotm_cols].sum(axis = 1)
        np.testing.assert_array_equal(ntotm.values, df_annual['PopulationM'].values)
    
    def test_number_positive_women_by_age(self, extended_model_output):
        # Check number positive women is same across all age groups
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        ntotf_cols = [c for c in df_annual.columns if ("NFage" in c)]
        ntotf = df_annual[ntotf_cols].sum(axis = 1)
        np.testing.assert_array_equal(ntotf.values, df_annual['PopulationF'].values)
    
    def test_proportion_on_art(self, extended_model_output):
        # Check that proportion on ART is correct when looking across all genders
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        pronart = (df_annual['NOnARTM'] + \
            df_annual['NOnARTF'])/df_annual['NumberPositive']
        pronart.replace(np.nan, 0, inplace = True)
        
        np.testing.assert_allclose(pronart.values, df_annual['PropHIVPosONART'].values, 
            atol = 1e-5)
    
    def test_new_cases_this_year_by_risk(self, extended_model_output):
        # Check that the new cases this year, stratified by risk, adds to total
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        new_cases = df_annual['NewCasesThisYear_riskLow'] + \
                    df_annual['NewCasesThisYear_riskMed'] + \
                    df_annual['NewCasesThisYear_riskHigh']
        
        np.testing.assert_array_equal(new_cases.values, \
            df_annual['NewCasesThisYear'].values)
    
    def test_ndied_from_hiv_by_risk(self, extended_model_output):
        # Check that the number died from HIV, stratified by risk, adds to total
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        new_cases = df_annual['NDied_from_HIV_riskLow'].cumsum() + \
                    df_annual['NDied_from_HIV_riskMed'].cumsum() + \
                    df_annual['NDied_from_HIV_riskHigh'].cumsum()
        
        np.testing.assert_array_equal(new_cases.values, df_annual['NDied_from_HIV'].values)
    
    def test_prevalence_by_risk(self, extended_model_output):
        # Check that prevalence for each risk group (and prop in that risk group) checks out with 
        # overall prevalence (to about 4 d.p.)
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        prevalence = df_annual['Prop_riskLow'] * df_annual['Prevalence_riskLow'] + \
                     df_annual['Prop_riskMed'] * df_annual['Prevalence_riskMed'] + \
                     df_annual['Prop_riskHigh'] * df_annual['Prevalence_riskHigh']
        
        np.testing.assert_allclose(prevalence.values, \
            df_annual['Prevalence'].values, atol = 1e-4)
    
    
    def test_total_men_by_age(self, extended_model_output):
        # Check that the total number of men are consistent when summed across all ages
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        cols = [c for c in df_annual.columns if "NMage" in c]
        total_men = df_annual[cols].sum(axis = 1)

        np.testing.assert_array_equal(total_men.values, df_annual['PopulationM'].values)

    def test_total_women_by_age(self, extended_model_output):
        # Check that the total number of women are consistent when summed across all ages
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        cols = [c for c in df_annual.columns if "NFage" in c]
        total_women = df_annual[cols].sum(axis = 1)

        np.testing.assert_array_equal(total_women.values, df_annual['PopulationF'].values)

    def test_total_positive_men_by_age(self, extended_model_output):
        # Check that the total number of men are consistent when summed across all ages
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        cols = [c for c in df_annual.columns if "NPosMage" in c]
        total_positive_men = df_annual[cols].sum(axis = 1)

        np.testing.assert_array_equal(
            total_positive_men.values,
            df_annual['NumberPositiveM'].values)

    def test_total_positive_women_by_age(self, extended_model_output):
        # Check that the total number of women are consistent when summed across all ages
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        cols = [c for c in df_annual.columns if "NPosFage" in c]
        total_positive_women = df_annual[cols].sum(axis = 1)

        np.testing.assert_array_equal(
            total_positive_women.values, 
            df_annual['NumberPositiveF'].values)

    def test_total_cases_by_age_by_sex(self, extended_model_output):
        # Check that the total number incident cases stratified by age and sex is consistent with
        # total number of incident cases
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        cols = [c for c in df_annual.columns if "IncFage" in c]
        total_cases_women = df_annual[cols].sum(axis = 1)

        cols = [c for c in df_annual.columns if "IncMage" in c]
        total_cases_men = df_annual[cols].sum(axis = 1)

        total_cases = total_cases_men.values + total_cases_women.values

        np.testing.assert_array_equal(total_cases, df_annual['NewCasesThisYear'].values)
    
    ####################################################################################
    ###################### Internal consistency of Timestep_outputs ####################
    ####################################################################################
    
    def test_timestep_number_positive_women_by_status(self, extended_model_output):
        # Check that the number of women that know/do not know their status is equal to the 
        # total number of HIV positive women
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        total_positive_women = df_timestep['NNotKnowStatus_f'] + df_timestep['N_knowpos_f']
        
        np.testing.assert_array_equal(total_positive_women.values, df_timestep['NPos_f'].values)

    def test_timestep_number_positive_men_by_status(self, extended_model_output):
        # Check that the number of men that know/do not know their status is equal to the 
        # total number of HIV positive men
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        total_positive_men = df_timestep['NNotKnowStatus_m'] + df_timestep['N_knowpos_m']
        
        np.testing.assert_array_equal(total_positive_men.values, df_timestep['NPos_m'].values)

    ####################################################################################
    ############### Consistency of Timestep_ and complete transmission network #########
    ####################################################################################
    
    
    
    ####################################################################################
    ############### Consistency of Annual_ and complete transmission network ###########
    ####################################################################################
    
    ##########################
    # Cross-section measures #
    ##########################
    
    def test_total_number_annual_transmission(self, extended_model_output):
        # Check the total number of individuals (aged 14+) is the same in the complete 
        # transmission network and the annual outputs file on 31st December of each year
        # (or start of following year)
        
        # Calculate total population in complete network file for all years listed in
        # "Annual_output" file
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        total_population_annual = []
        for year in df_annual.Year.values:
            
            total_population = np.sum((df_individual.DoB + (c.AGE_ADULT + 1) <= year+1) & \
                ((df_individual.DoD >= year+1) | (df_individual.DoD == -1)))
            
            total_population_annual.append(total_population)
        
        total_population_annual = np.array(total_population_annual)
        
        np.testing.assert_array_equal(total_population_annual, df_annual['TotalPopulation'])

    def test_total_number_men_annual_transmission(self, extended_model_output):
        # Check the total number of men (aged 14+) is the same in the complete 
        # transmission network and the annual outputs file on 31st December of each year
        # (or start of following year)
        
        # Calculate total male population in complete network file 
        # for all years listed in "Annual_output" file
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        total_population_annual = []
        for year in df_annual.Year.values:
            
            total_population = np.sum((df_individual.DoB + (c.AGE_ADULT + 1) <= year+1) & \
                ((df_individual.DoD >= year+1) | (df_individual.DoD == -1)) & \
                (df_individual.Sex == "M"))
            
            total_population_annual.append(total_population)
        
        total_population_annual = np.array(total_population_annual)
        
        np.testing.assert_array_equal(total_population_annual, df_annual['PopulationM'])

    def test_total_number_women_annual_transmission(self, extended_model_output):
        # Check the total number of men (aged 14+) is the same in the complete 
        # transmission network and the annual outputs file on 31st December of each year
        
        # Calculate total male population in complete network file 
        # for all years listed in "Annual_output" file
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        total_population_annual = []
        for year in df_annual.Year.values:
            
            total_population = np.sum((df_individual.DoB + (c.AGE_ADULT + 1) <= year+1) & \
                ((df_individual.DoD >= year+1) | (df_individual.DoD == -1)) & \
                (df_individual.Sex == "F"))
            
            total_population_annual.append(total_population)
        
        total_population_annual = np.array(total_population_annual)
        
        np.testing.assert_array_equal(total_population_annual, df_annual['PopulationF'])

    def test_total_number_positive_men_annual_transmission(self, extended_model_output):
        # Check the total number of HIV-positive men (aged 14+) is the same in the complete
        # transmission network and the annual outputs file on 31st December of each year        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge "TimeOfInfection" into individual file
        df_individual_merged = pd.merge(df_individual, \
            df_transmission[["IdInfected", "TimeOfInfection"]], \
            left_on = "Id", right_on = "IdInfected", how = "left")
        
        # Calculate total HIV+ male population in complete network file
        # for all years listed in "Annual_output" file
        
        total_positive_men = []
        for year in df_annual.Year.values:
            
            # Check conditions that 1) was >14 yo in year in question, 2) not dead, 
            # 3) HIV+, 4) a man
            conditions = (df_individual_merged.DoB + (c.AGE_ADULT + 1) <= year+1) & \
            ((df_individual_merged.DoD >= year+1) | (df_individual_merged.DoD == -1)) & \
            (df_individual_merged.HIV_pos == 1) & \
            (df_individual_merged.TimeOfInfection < year+1) & \
            (df_individual_merged.Sex == "M")

            total_positive_men.append(np.sum(conditions))
        
        # Convert to numpy array
        total_positive_men = np.array(total_positive_men)
        
        np.testing.assert_array_equal(total_positive_men, df_annual['NumberPositiveM'])

    def test_total_number_positive_women_annual_transmission(self, extended_model_output):
        # Check the total number of HIV-positive women (aged 14+) is the same in the complete
        # transmission network and the annual outputs file on 31st December of each year
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge "TimeOfInfection" into individual file
        df_individual_merged = pd.merge(df_individual, \
            df_transmission[["IdInfected", "TimeOfInfection"]], \
            left_on = "Id", right_on = "IdInfected", how = "left")
        
        # Calculate total HIV+ female population in complete network file
        # for all years listed in "Annual_output" file
        
        total_positive_women = []
        for year in df_annual.Year.values:
            
            # Check conditions that 1) was >14 yo in year in question, 2) not dead, 
            # 3) HIV+, 4) a women
            conditions = (df_individual_merged.DoB + (c.AGE_ADULT + 1) <= year+1) & \
            ((df_individual_merged.DoD >= year+1) | (df_individual_merged.DoD == -1)) & \
            (df_individual_merged.HIV_pos == 1) & (df_individual_merged.TimeOfInfection < year+1) &\
            (df_individual_merged.Sex == "F")
            
            total_positive_women.append(np.sum(conditions))
        
        # Convert to numpy array
        total_positive_women = np.array(total_positive_women)
        
        np.testing.assert_array_equal(total_positive_women, df_annual['NumberPositiveF'])
    
    ##################################
    # Measures over an annual period #
    ##################################
    
    def test_ndead_annual_transmission(self, extended_model_output):
        # Calculate total number died in population in complete network file 
        # for all years listed in "Annual_output" file.  Compare with N_dead in Annual_output file.
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        ndied_complete = []
        for year in df_annual.Year.values:
            
            ndied = np.sum((df_individual.DoD >= year) & (df_individual.DoD < year+1))
            ndied_complete.append(ndied)
        
        ndied_complete = np.array(ndied_complete).cumsum()
        
        np.testing.assert_array_equal(ndied_complete, df_annual['N_dead'])

    def test_newcasesthisyear_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network
        # for all years listed in "Annual_output" file.  
        # Compare with NewCasesThisYear in Annual_output file.
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:
            
            newcasesthisyear = np.sum((df_transmission.TimeOfInfection >= year) & \
                (df_transmission.TimeOfInfection < year + 1))
            
            newcasesthisyear_complete.append(newcasesthisyear)
        
        newcasesthisyear_complete = np.array(newcasesthisyear_complete)
        
        np.testing.assert_array_equal(newcasesthisyear_complete, df_annual['NewCasesThisYear'])

    def test_newcasesthisyearfromoutside_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network caused by an 
        # individual in the outside patch (relative to the individual being infected)
        # for all years listed in "Annual_output" file.  
        # Compare with NewCasesThisYearFromOutside in Annual_output file.
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:
            
            newcasesthisyear = np.sum((df_transmission.TimeOfInfection >= year) & \
                (df_transmission.TimeOfInfection < year + 1) & \
                (df_transmission.IsInfectorOutsidePatch == 1))
            
            newcasesthisyear_complete.append(newcasesthisyear)
        
        newcasesthisyear_complete = np.array(newcasesthisyear_complete)
        
        np.testing.assert_array_equal(newcasesthisyear_complete,
            df_annual['NewCasesThisYearFromOutside'])

    def test_newcasesthisyearfromacute_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network caused by an 
        # individual in the acute phase of HIV infection
        # for all years listed in "Annual_output" file.  
        # Compare with NewCasesThisYearFromAcute in Annual_output file.
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:
            
            newcasesthisyear = np.sum((df_transmission.TimeOfInfection >= year) & \
                (df_transmission.TimeOfInfection < year + 1) & \
                (df_transmission.IsInfectorAcute == 1))
            
            newcasesthisyear_complete.append(newcasesthisyear)
        
        newcasesthisyear_complete = np.array(newcasesthisyear_complete)
        
        np.testing.assert_array_equal(newcasesthisyear_complete, 
            df_annual['NewCasesThisYearFromAcute'])

    def test_newcasesthisyear_13_18_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network for men ages
        # 13-18 (at start of the year) and for all years listed in "Annual_output" file.
        # Compare with IncMage13-18 in Annual_output file.
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge date of birth (DoB) and sex of infected individual into transmission file
        df_transmission_merged = pd.merge(df_transmission, \
            df_individual[["Id", "Sex", "DoB"]], \
            left_on = "IdInfected", right_on = "Id", how = "left")
        
        # Calculate age of the infected individual at start of the year of infection
        df_transmission_merged["AgeInfected"] = np.floor(df_transmission_merged.TimeOfInfection) -\
            df_transmission_merged.DoB
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:

            conditions = (df_transmission_merged.TimeOfInfection >= year) & \
                (df_transmission_merged.TimeOfInfection < year+1) & \
                (df_transmission_merged.Sex == "M") & \
                (df_transmission_merged.AgeInfected >= c.AGE_ADULT) & \
                (df_transmission_merged.AgeInfected < 18)

            newcasesthisyear_complete.append(np.sum(conditions))

        newcasesthisyear_complete = np.array(newcasesthisyear_complete)

        np.testing.assert_array_equal(newcasesthisyear_complete, df_annual['IncMage13-18'])

    def test_newcasesthisyear_18_23_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network for men ages
        # 18-23 (at start of the year) and for all years listed in "Annual_output" file.
        # Compare with IncMage18-23 in Annual_output file.
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge date of birth (DoB) and sex of infected individual into transmission file
        df_transmission_merged = pd.merge(df_transmission, \
            df_individual[["Id", "Sex", "DoB"]], \
            left_on = "IdInfected", right_on = "Id", how = "left")
        
        # Calculate age of the infected individual at start of the year of infection
        df_transmission_merged["AgeInfected"] = np.floor(df_transmission_merged.TimeOfInfection) -\
            df_transmission_merged.DoB
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:

            conditions = (df_transmission_merged.TimeOfInfection >= year) & \
                (df_transmission_merged.TimeOfInfection < year+1) & \
                (df_transmission_merged.Sex == "M") & \
                (df_transmission_merged.AgeInfected >= 18) & \
                (df_transmission_merged.AgeInfected < 23)

            newcasesthisyear_complete.append(np.sum(conditions))

        newcasesthisyear_complete = np.array(newcasesthisyear_complete)

        np.testing.assert_array_equal(newcasesthisyear_complete, df_annual['IncMage18-23'])

    def test_newcasesthisyear_23_30_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network for men ages
        # 23-30 (at start of the year) and for all years listed in "Annual_output" file.
        # Compare with IncMage23-30 in Annual_output file.
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge date of birth (DoB) and sex of infected individual into transmission file
        df_transmission_merged = pd.merge(df_transmission, \
            df_individual[["Id", "Sex", "DoB"]], \
            left_on = "IdInfected", right_on = "Id", how = "left")
        
        # Calculate age of the infected individual at start of the year of infection
        df_transmission_merged["AgeInfected"] = np.floor(df_transmission_merged.TimeOfInfection) -\
            df_transmission_merged.DoB
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:

            conditions = (df_transmission_merged.TimeOfInfection >= year) & \
                (df_transmission_merged.TimeOfInfection < year+1) & \
                (df_transmission_merged.Sex == "M") & \
                (df_transmission_merged.AgeInfected >= 23) & \
                (df_transmission_merged.AgeInfected < 30)

            newcasesthisyear_complete.append(np.sum(conditions))

        newcasesthisyear_complete = np.array(newcasesthisyear_complete)

        np.testing.assert_array_equal(newcasesthisyear_complete, df_annual['IncMage23-30'])

    def test_newcasesthisyear_30_40_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network for men ages
        # 30-40 (at start of the year) and for all years listed in "Annual_output" file.
        # Compare with IncMage30-40 in Annual_output file.
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge date of birth (DoB) and sex of infected individual into transmission file
        df_transmission_merged = pd.merge(df_transmission, \
            df_individual[["Id", "Sex", "DoB"]], \
            left_on = "IdInfected", right_on = "Id", how = "left")
        
        # Calculate age of the infected individual at start of the year of infection
        df_transmission_merged["AgeInfected"] = np.floor(df_transmission_merged.TimeOfInfection) -\
            df_transmission_merged.DoB
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:

            conditions = (df_transmission_merged.TimeOfInfection >= year) & \
                (df_transmission_merged.TimeOfInfection < year+1) & \
                (df_transmission_merged.Sex == "M") & \
                (df_transmission_merged.AgeInfected >= 30) & \
                (df_transmission_merged.AgeInfected < 40)

            newcasesthisyear_complete.append(np.sum(conditions))

        newcasesthisyear_complete = np.array(newcasesthisyear_complete)

        np.testing.assert_array_equal(newcasesthisyear_complete, df_annual['IncMage30-40'])

    def test_newcasesthisyear_40_50_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network for men ages
        # 40-50 (at start of the year) and for all years listed in "Annual_output" file.
        # Compare with IncMage40-50 in Annual_output file.
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge date of birth (DoB) and sex of infected individual into transmission file
        df_transmission_merged = pd.merge(df_transmission, \
            df_individual[["Id", "Sex", "DoB"]], \
            left_on = "IdInfected", right_on = "Id", how = "left")
        
        # Calculate age of the infected individual at start of the year of infection
        df_transmission_merged["AgeInfected"] = np.floor(df_transmission_merged.TimeOfInfection) -\
            df_transmission_merged.DoB
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:

            conditions = (df_transmission_merged.TimeOfInfection >= year) & \
                (df_transmission_merged.TimeOfInfection < year+1) & \
                (df_transmission_merged.Sex == "M") & \
                (df_transmission_merged.AgeInfected >= 40) & \
                (df_transmission_merged.AgeInfected < 50)

            newcasesthisyear_complete.append(np.sum(conditions))

        newcasesthisyear_complete = np.array(newcasesthisyear_complete)

        np.testing.assert_array_equal(newcasesthisyear_complete, df_annual['IncMage40-50'])

    def test_newcasesthisyear_50_60_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network for men ages
        # 50-60 (at start of the year) and for all years listed in "Annual_output" file.
        # Compare with IncMage50-60 in Annual_output file.
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge date of birth (DoB) and sex of infected individual into transmission file
        df_transmission_merged = pd.merge(df_transmission, \
            df_individual[["Id", "Sex", "DoB"]], \
            left_on = "IdInfected", right_on = "Id", how = "left")
        
        # Calculate age of the infected individual at start of the year of infection
        df_transmission_merged["AgeInfected"] = np.floor(df_transmission_merged.TimeOfInfection) -\
            df_transmission_merged.DoB
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:

            conditions = (df_transmission_merged.TimeOfInfection >= year) & \
                (df_transmission_merged.TimeOfInfection < year+1) & \
                (df_transmission_merged.Sex == "M") & \
                (df_transmission_merged.AgeInfected >= 50) & \
                (df_transmission_merged.AgeInfected < 60)
                
            newcasesthisyear_complete.append(np.sum(conditions))

        newcasesthisyear_complete = np.array(newcasesthisyear_complete)

        np.testing.assert_array_equal(newcasesthisyear_complete, df_annual['IncMage50-60'])

    def test_newcasesthisyear_60_annual_transmission(self, extended_model_output):
        # Calculate total number of cases in a year from the transmission network for men ages
        # 60+ (at start of the year) and for all years listed in "Annual_output" file.
        # Compare with IncMage60-80 in Annual_output file.
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Merge date of birth (DoB) and sex of infected individual into transmission file
        df_transmission_merged = pd.merge(df_transmission, \
            df_individual[["Id", "Sex", "DoB"]], \
            left_on = "IdInfected", right_on = "Id", how = "left")
        
        # Calculate age of the infected individual at start of the year of infection
        df_transmission_merged["AgeInfected"] = np.floor(df_transmission_merged.TimeOfInfection) -\
            df_transmission_merged.DoB
        
        newcasesthisyear_complete = []
        for year in df_annual.Year.values:

            conditions = (df_transmission_merged.TimeOfInfection >= year) & \
                (df_transmission_merged.TimeOfInfection < year+1) & \
                (df_transmission_merged.Sex == "M") & \
                (df_transmission_merged.AgeInfected >= 60)

            newcasesthisyear_complete.append(np.sum(conditions))

        newcasesthisyear_complete = np.array(newcasesthisyear_complete)

        np.testing.assert_array_equal(newcasesthisyear_complete, df_annual['IncMage60-80'])

    ####################################################################################
    ############### Consistency of Annual_ and Timestep_ outputs #######################
    ####################################################################################
    
    ##########################
    # Cross-section measures #
    ##########################
    
    def test_number_women(self, extended_model_output):
        # Check that the number of women is the same in the listed years in both files
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        timestep = df_timestep[df_timestep.Time.isin(years)]['N_f']
        annual = df_annual[df_annual.Year.isin(years)]['PopulationF']
    
        np.testing.assert_array_equal(timestep.values[1:], annual.values[:-1])
    
    def test_number_men(self, extended_model_output):
        # Check that the number of men is the same in the listed years in both files
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        timestep = df_timestep[df_timestep.Time.isin(years)]['N_m']
        annual = df_annual[df_annual.Year.isin(years)]['PopulationM']
        np.testing.assert_array_equal(timestep.values[1:], annual.values[:-1])
    
    def test_number_positive_women(self, extended_model_output):
        # Check that the number of HIV positive women is the same in the listed years in both files
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        timestep = df_timestep[df_timestep.Time.isin(years)]['NPos_f']
        annual = df_annual[df_annual.Year.isin(years)]['NumberPositiveF']
        np.testing.assert_array_equal(timestep.values[1:], annual.values[:-1])
    
    def test_number_positive_men(self, extended_model_output):
        # Check that the number of HIV positive men is the same in the listed years in both files
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        timestep = df_timestep[df_timestep.Time.isin(years)]['NPos_m']
        annual = df_annual[df_annual.Year.isin(years)]['NumberPositiveM']
        np.testing.assert_array_equal(timestep.values[1:], annual.values[:-1])
    
    def test_number_art_women(self, extended_model_output):
        # Check that the number of women on ART between both files
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        timestep = df_timestep[df_timestep.Time.isin(years)]['NART_f']
        annual = df_annual[df_annual.Year.isin(years)]['NOnARTF']
        np.testing.assert_array_equal(timestep.values[1:], annual.values[:-1])
    
    def test_number_art_men(self, extended_model_output):
        # Check that the number of men on ART between both files
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        timestep = df_timestep[df_timestep.Time.isin(years)]['NART_m']
        annual = df_annual[df_annual.Year.isin(years)]['NOnARTM']
        np.testing.assert_array_equal(timestep.values[1:], annual.values[:-1])
    
    def test_prop_men_circ(self, extended_model_output):
        # Check that the number of men on ART between both files
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        timestep = df_timestep[df_timestep.Time.isin(years)]['PropMenCirc']
        annual = df_annual[df_annual.Year.isin(years)]['PropMenCirc']
        np.testing.assert_array_equal(timestep.values[1:], annual.values[:-1])
    
    ##################################
    # Measures over an annual period #
    ##################################
    
    def test_incident_cases(self, extended_model_output):
        # Check annual incident cases are calculated correctly between the two files
        
        df_annual, df_timestep, df_individual, df_transmission = extended_model_output
        
        # Isolate cumulative cases in timestep outputs
        timestepf = df_timestep[df_timestep.Time.isin(years)]['Cumulative_Infected_f']
        timestepm = df_timestep[df_timestep.Time.isin(years)]['Cumulative_Infected_m']
        
        # Find total new cases from cumulative cases by sex
        timestep = timestepf.diff() + timestepm.diff()
        annual = df_annual[df_annual.Year.isin(years)]['NewCasesThisYear']
        
        # Remove the first value (the diff command means the first value is an NaN).
        np.testing.assert_array_equal(timestep[1:].values, annual[:-1].values)
