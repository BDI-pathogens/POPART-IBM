#!/usr/bin/env python3
"""
Constants associated with tests of the PopART-IBM model

W. Probert, 2020
"""

from os.path import join

# Directories
#INPUTDIR = join("data", "SAMPLED_PARAMETERS", "PARAMS_COMMUNITY5")
NRUNS = 1
IBM_DIR = "src"
IBM_DIR_TEST = IBM_DIR + "_test"
DATA_DIR_TEST = "data_test"

TEST_PARAMS_DIR = "tests/data"

EXE = "popart-simul.exe"

# Construct the executable command
command = join(IBM_DIR_TEST, EXE)

# Filenames
f_timestep = "Timestep_outputs_CL05_Za_A_V1.2_patch0_Rand10_Run1_PCseed0_0.csv"
f_annual = "Annual_outputs_CL05_Za_A_V1.2_patch0_Rand10_Run1_PCseed0_0.csv"
f_cost_effectiveness = "Cost_effectiveness_CL05_Za_A_V1.2_patch0_Rand10_Run1_PCseed0_0.csv"
f_treats = "TREATS_CL05_Za_A_V1.2_patch0_Rand10_Run1_PCseed0_0.csv"
f_individual = "phylogenetic_individualdata_CL05_Za_A_V1.2_patch0_Rand10_Run1_PCseed0_0.csv"
f_transmission = "phylogenetic_transmission_CL05_Za_A_V1.2_patch0_Rand10_Run1_PCseed0_0.csv"
f_art_status = "ART_status_by_age_sex_CL05_Za_A_V1.2_patch0_Rand10_Run1_PCseed0_0.csv"

# Age at which individuals are created in the model
AGE_ADULT = 13

# Number of timesteps per year
NTIMESTEPS = 48



parameter_transpose_file = "tests/data/parameters_transpose.csv"
fertility_file = "tests/data/parameters_fertility.csv"
mortality_file = "tests/data/parameters_mortality.csv"
patchinfo_file = "tests/data/patchinfo.txt"
python_seed_file = "tests/data/python_seed.txt"
