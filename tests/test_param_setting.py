#!/usr/bin/python3
"""
Tests of the setting/getting of parameter values
"""

import sys, pytest, numpy as np, pandas as pd, shutil, os
from os.path import join
from pathlib import Path

sys.path.append("./src/popart_ibm")
from parameters import ParameterSet, SamplingSchedule
import constants as c


# class-level setup/teardown
@pytest.fixture(scope = "class", autouse = True)
def compile_popart_ibm(request):
    pass

# class-level setup/teardown
@pytest.fixture(scope = "class", autouse = True)
def make_test_dir(request):
    os.mkdir(c.DATA_DIR_TEST)
    def fin():
        shutil.rmtree(c.DATA_DIR_TEST)

# Method ("function") setup/teardown
@pytest.fixture(scope = "function", autouse = True)
def setup_popart_ibm_methods(request):
    pass

@pytest.fixture(scope = "function")
def basic_model_output():
    pass

def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = sorted(funcarglist[0])
    metafunc.parametrize(
        argnames, [[funcargs[name] for name in argnames] for funcargs in funcarglist]
    )


class TestClass(object):
    params = {
        "test_get_param": [
            dict(
                file_type = "HIV",
                df = pd.DataFrame({"assortativity": [0.75]}),
                param_name = "assortativity"
            ),
            dict(
                file_type = "cascade",
                df = pd.DataFrame({"p_child_circ": [0.1532]}),
                param_name = "p_child_circ"
            ),
        ],
        "test_set_param": [
            dict(
                file_type = "HIV",
                df = pd.DataFrame({"assortativity": [0.95]}),
                param_name = "assortativity"
            ),
            dict(
                file_type = "cascade",
                df = pd.DataFrame({"p_child_circ": [0.1254]}),
                param_name = "p_child_circ"
            ),
        ],
        "test_read_patchinfo": [dict()],
        "test_sampling_schedule" : [
        dict(
            n_timestep_per_year = 1,
            starting_times = [[2013, 0]],
            ending_times = [[2013, 0]],
            final_coverage = 0.8,
            ages = [19],
            df_expected = pd.DataFrame({
                            "time": [2013.0], 
                            "year" : [2013], 
                            "t_step" : [0], 
                            "NM19": [0.8],
                            "NF19": [0.8]})
        ), 
        dict(
            n_timestep_per_year = 1,
            starting_times = [[2013, 0]],
            ending_times = [[2014, 0]],
            final_coverage = 0.8,
            ages = [19, 20],
            df_expected = pd.DataFrame({
                            "time": [2013.0, 2014.0], 
                            "year" : [2013, 2014], 
                            "t_step" : [0, 0], 
                            "NM19": [0.4, 0.4],
                            "NM20": [0.4, 0.4],
                            "NF19": [0.4, 0.4],
                            "NF20": [0.4, 0.4]})
        ),
        dict(
            n_timestep_per_year = 2,
            starting_times = [[2013, 0]],
            ending_times = [[2014, 0]],
            final_coverage = 0.6,
            ages = [19, 20],
            df_expected = pd.DataFrame({
                            "time": [2013.0, 2013.5, 2014.0], 
                            "year" : [2013, 2013, 2014], 
                            "t_step" : [0, 1, 0], 
                            "NM19": [0.2, 0.2, 0.2],
                            "NM20": [0.2, 0.2, 0.2],
                            "NF19": [0.2, 0.2, 0.2],
                            "NF20": [0.2, 0.2, 0.2]})
        ),
        dict(
            n_timestep_per_year = 4,
            starting_times = [[2013, 0]],
            ending_times = [[2013, 3]],
            final_coverage = 0.9,
            ages = [52, 41],
            df_expected = pd.DataFrame({
                            "time": [2013.0, 2013.25, 2013.5, 2013.75], 
                            "year" : [2013, 2013, 2013, 2013], 
                            "t_step" : [0, 1, 2, 3], 
                            "NM52": [0.225, 0.225, 0.225, 0.225],
                            "NM41": [0.225, 0.225, 0.225, 0.225],
                            "NF52": [0.225, 0.225, 0.225, 0.225],
                            "NF41": [0.225, 0.225, 0.225, 0.225]})
        )],
        "test_read_test_params_param_name": [
            dict(param_name = "p_vu_becomes_virally_suppressed"), # cascade
            dict(param_name = "sex_ratio"), # demographics
            dict(param_name = "p_child_circ"), # HIV
            dict(param_name = "initial_adult_population_size"), # init
            dict(param_name = "prop_compromise_from_males"), # partnerships
            dict(param_name = "PC_Retention_Round2"), # PC
            dict(param_name = "COUNTRY_HIV_TEST_START") # times
        ],
        "test_read_test_params": [dict()],
        "test_write_param_set": [dict()],
        "test_read_param_set": [dict()]
        
    }
    def test_get_param(self, file_type, df, param_name):
        """
        Test the get_param method of the ParameterSet class
        """
        
        test_filename = join(c.DATA_DIR_TEST, f"test_{file_type}.csv")
        df.to_csv(test_filename, index = False)
        
        params = ParameterSet()
        params.read_param_file(param_file = test_filename)
    
        np.testing.assert_equal(
            params.get_param(param_name)[0], 
            df[param_name].values[0])
        
        np.testing.assert_array_equal(params.params[file_type].to_numpy(), df.to_numpy())

    def test_set_param(self, file_type, df, param_name):
        """
        Test the set_param method of the ParameterSet class
        """
        
        test_filename = join(c.DATA_DIR_TEST, f"test_{file_type}.csv")
        df.to_csv(test_filename, index = False)
        
        params = ParameterSet()
        params.read_param_file(param_file = test_filename)
        
        new_val = np.random.uniform()
        params.set_param(param_name, new_val)
        
        # Check the value has changed
        np.testing.assert_equal(
            params.get_param(param_name)[0] == df[param_name].values[0],
            False)
        
        np.testing.assert_equal(params.get_param(param_name)[0], new_val)
    
    def test_read_patchinfo(self):
        """
        Test the ParameterSet class can read the patchinfo file correctly
        """
        
        p = ParameterSet()
        p.read_patchinfo(c.patchinfo_file)
        
        np.testing.assert_array_equal( 
            np.array(p.patchinfo), 
            np.array([[5, 4], [1, 0]])
        )

    def test_sampling_schedule(self, 
        n_timestep_per_year, 
        starting_times, 
        ending_times, 
        final_coverage, 
        ages, 
        df_expected):
        """
        Test the creation of sampling schedules as used for CHiPs and PC sampling schedules
        """
        S = SamplingSchedule(
            n_timestep_per_year = n_timestep_per_year, 
            starting_times = starting_times,
            ending_times = ending_times, 
            column_name_prefixes = ["N"], 
            ages = ages, 
            sexes = ['M', 'F'])
        
        df_schedule = S.create_uniform_schedule(
            S.starting_times[0], 
            S.ending_times[0], 
            final_coverage = final_coverage)
        
        pd.testing.assert_frame_equal(df_schedule, df_expected)
    
    def test_read_test_params_param_name(self, param_name):
        
        df_params_trans = pd.read_csv(c.parameter_transpose_file)
        val = df_params_trans.loc[df_params_trans.Name == param_name, "Value"].values[0]
        
        p = ParameterSet()
        p.read_test_parameters(
            parameter_transpose_file = c.parameter_transpose_file,
            fertility_file = c.fertility_file, 
            mortality_file = c.mortality_file, 
            patchinfo_file = c.patchinfo_file, 
            python_seed_file = c.python_seed_file)
        
        # Test the parameter is the same btw transpose file and that in the param object
        np.testing.assert_equal(p.get_param(param_name)[0], val)
    
    def test_read_test_params(self):
        
        p = ParameterSet()
        p.read_test_parameters(
            parameter_transpose_file = c.parameter_transpose_file,
            fertility_file = c.fertility_file, 
            mortality_file = c.mortality_file, 
            patchinfo_file = c.patchinfo_file, 
            python_seed_file = c.python_seed_file)
        
        # Test we can set/get params from these test params
        p.set_param("prop_compromise_from_males", 0.9)
        np.testing.assert_equal(p.get_param("prop_compromise_from_males")[0], 0.9)
        
        # Check the patchinfo file is read-in correctly
        np.testing.assert_array_equal( np.array(p.patchinfo), np.array([[5, 4], [1, 0]]) )
        
        # Check final coverage of the test CHiPs schedules
        np.testing.assert_equal( p.params["chipsuptake_round1"].sum(axis = 0)[10], 0.6)
        np.testing.assert_equal( p.params["chipsuptake_round2"].sum(axis = 0)[10], 0.7)
        np.testing.assert_equal( p.params["chipsuptake_round3"].sum(axis = 0)[10], 0.8)
        
        # Check dimensions of the post-trial CHiPs schedule
        np.testing.assert_equal( 
            p.params["chipsuptake_roundposttrial"].shape, 
            (c.NTIMESTEPS, len(range(18, 81))*2 + 3))
        
        # Check CHiPs schedules are created correctly
        np.testing.assert_equal( p.params["chipsuptake_round1"].sum(axis = 0)[10], 0.6)
        np.testing.assert_equal( p.params["chipsuptake_round2"].sum(axis = 0)[10], 0.7)
        np.testing.assert_equal( p.params["chipsuptake_round3"].sum(axis = 0)[10], 0.8)
        
        np.testing.assert_equal( 
            p.params["chipsuptake_roundposttrial"].shape, 
            (c.NTIMESTEPS, len(range(18, 81))*2 + 3))
        
        # Test PC schedules are created correctly
        np.testing.assert_equal( p.params["PC0"].shape[1], len(np.arange(18, 45))*2*3 + 3)
        np.testing.assert_equal( p.params["PC1"].shape[1], len(np.arange(18, 45))*2*3 + 3)
        np.testing.assert_equal( p.params["PC2"].shape[1], len(np.arange(18, 45))*2*3 + 3)
        np.testing.assert_equal( p.params["PC3"].shape[1], len(np.arange(18, 45))*2*3 + 3)
        np.testing.assert_equal( p.params["PC0"].sum(axis = 0)[10], 0)
        np.testing.assert_equal( p.params["PC1"].sum(axis = 0)[10], 0)
        np.testing.assert_equal( p.params["PC2"].sum(axis = 0)[10], 0)
        np.testing.assert_equal( p.params["PC3"].sum(axis = 0)[10], 0)
    
    def test_write_param_set(self):
        """
        Test we can write a whole parameter set to file
        """
        
        p = ParameterSet()
        p.read_test_parameters(
            parameter_transpose_file = c.parameter_transpose_file,
            fertility_file = c.fertility_file, 
            mortality_file = c.mortality_file, 
            patchinfo_file = c.patchinfo_file, 
            python_seed_file = c.python_seed_file)
        
        p.write_param_set(output_dir = join(c.DATA_DIR_TEST, "new_set"))
        
        np.testing.assert_equal(os.path.isdir(join(c.DATA_DIR_TEST, "new_set")), True)
        
        # Remove the folder
        shutil.rmtree(join(c.DATA_DIR_TEST, "new_set"))
    
    def test_read_param_set(self):
        
        q = ParameterSet()
        q.read_test_parameters(
            parameter_transpose_file = c.parameter_transpose_file,
            fertility_file = c.fertility_file, 
            mortality_file = c.mortality_file, 
            patchinfo_file = c.patchinfo_file, 
            python_seed_file = c.python_seed_file)
        
        q.write_param_set(output_dir = join(c.DATA_DIR_TEST, "new_set"))
        
        # Read the parameters from the written dir
        p = ParameterSet(param_input = join(c.DATA_DIR_TEST, "new_set"))
        
        # Test we can set/get params from these test params
        p.set_param("prop_compromise_from_males", 0.9)
        np.testing.assert_equal(p.get_param("prop_compromise_from_males")[0], 0.9)
        
        # Check the patchinfo file is read-in correctly
        np.testing.assert_array_equal( np.array(p.patchinfo), np.array([[5, 4], [1, 0]]) )
        
        # Check final coverage of the test CHiPs schedules
        np.testing.assert_equal( p.params["chipsuptake_round1"].sum(axis = 0)[10], 0.6)
        np.testing.assert_equal( p.params["chipsuptake_round2"].sum(axis = 0)[10], 0.7)
        np.testing.assert_equal( p.params["chipsuptake_round3"].sum(axis = 0)[10], 0.8)
        
        # Check dimensions of the post-trial CHiPs schedule
        np.testing.assert_equal( 
            p.params["chipsuptake_roundposttrial"].shape, 
            (c.NTIMESTEPS, len(range(18, 81))*2 + 3))
        
        # Check CHiPs schedules are created correctly
        np.testing.assert_equal( p.params["chipsuptake_round1"].sum(axis = 0)[10], 0.6)
        np.testing.assert_equal( p.params["chipsuptake_round2"].sum(axis = 0)[10], 0.7)
        np.testing.assert_equal( p.params["chipsuptake_round3"].sum(axis = 0)[10], 0.8)
    
        np.testing.assert_equal( 
            p.params["chipsuptake_roundposttrial"].shape, 
            (c.NTIMESTEPS, len(range(18, 81))*2 + 3))
        