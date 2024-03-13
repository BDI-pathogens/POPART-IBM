#!/usr/bin/python3
"""
ParameterSet class for organising parameter files as part of the individual-based model, 
PopART-IBM, as part of the HPTN 071 (PopART) trial.  

Usage
-----

from parameters import ParameterSet

##############################
# Reading a full parameter set
# ----------------------------

# Read a parameter set from a folder
params = ParameterSet(param_input = "PARAMS_COMMUNITY5")

# Get a parameter value (for both patches)
params.get_param("assortativity")
> 0.5

# Set a parameter value for both patches
params.set_param("assortativity", 0.921)

# Print the assortativity parameter within the partnerships file: 
params.partnerships.assortativity
> 0.921

# Print the complete "fertility" parameter file as a pandas dataframe
params.fertility

# Can also be accessed via the dictionary
params.params['fertility']

# Print the CHiPs uptake file for round 1 as a pandas dataframe
params.chipsuptake_round1

# Write a full parameter set to a new directory "PARAM_SET"
params.write_param_set(output_dir = "PARAM_SET")

#################################
# Reading a single parameter file
# -------------------------------

# Read a single parameter file
params = ParameterSet(param_input = "PARAMS_COMMUNITY5/param_processed_patchinfo.txt")

# Print the patchinfo file
params.params['param_processed_patchinfo.txt']

# Read the file with HIV-related parameters
params.read_param_file("PARAMS_COMMUNITY5/param_processed_patch0_HIV.csv")

# Print the file with HIV-related parameters (or simply params.HIV)
params.params['HIV']

# Change the efficacy of VMMC
params.set_param("eff_circ_vmmc", 0.5)

# Return the VMMC parameter
params.get_param("eff_circ_vmmc")
> 0.5

# Write the HIV file to current directory
params.write_param_file("HIV", output_dir = ".")

ls
> param_processed_patch0_HIV.csv param_processed_patch1_HIV.csv


Author: W. Probert
Created: Jan. 2019
"""

import copy, itertools, json, sys, os, re
from collections import OrderedDict
from os.path import join
from pathlib import Path
import pandas as pd, numpy as np

class AmbiguousFileName(Exception):
    pass

class AmbiguousParamName(Exception):
    pass


FILE_NAME_PREFIX = "param_processed"

# List of file names (or their suffixes)
FILE_NAMES_LINELIST = ["HIV", "cascade", "partnerships", "demographics", "times", "init", "PC"]
FILE_NAMES_DEMOGRAPHY = ["mortality", "fertility"]

FILE_NAMES_PATCH =["param_processed_patchinfo.txt"]
FILE_NAMES_MISC = ["python_seed.txt", "fitting_data_processed.txt"]

# Regular expressions for PC and CHiPs uptake files
regex_pc = re.compile(f".*{FILE_NAME_PREFIX}_patch._PC._community.*")
regex_chips = re.compile(f".*{FILE_NAME_PREFIX}_patch._chipsuptake.*")
regex_pc_schedule = re.compile("PC.")

class ParameterSet(object):
    """
    Class representing a parameter set(s) for the PopART-IBM model
    
    Arguments
    ---------
    param_input:
        path to a parameter file of the COVID19-IBM
    line_number :
        int of the line number of the parameter file to read (header line is 0)
    Methods
    -------
    
    read_param_dir
    
    
    get_param(param_name, file_name = None, line_number = None)
        List parameter value for parameter 'param_name'
        param_name: str
            Name of the parameter
        file_name: str
            Name of the file within which the parameter can be found (i.e. HIV, cascade, etc)
            Parameter files will be searched for "param_name" if file_name is not provided.  
        line_number: int
            Line number of the parameter file to present
    
    set_param(param_name, param_value)
        Set parameter value for 'param_name' to 'param_value'
    
    list_params()
        List ordered dictionary that stores parameters
    
    write_params(param_file)
        Write parameter set to csv file at location 'param_file'
    
    Optional keyword arguments
    --------------------------
    
    OUTPUT_SEP: str
        Separator for output files, currently a single whitespace (despite "CSV" file extensions).
    """
    def __init__(self, param_input = None, file_name = None, **kwargs):
        
        # Dictionaries for parameter files
        self.params = dict()
        self.misc = dict()
        
        self.FILE_NAMES_CHIPSUPTAKE = None
        self.FILE_NAMES_PC = None
        
        #############################
        # Optional keyword arguments
        # ---------------------------
        
        self.n_patches = kwargs.pop("n_patches", 2)
        self.n_chips_rounds = kwargs.pop("n_chips_rounds", 3)
        self.n_pc_rounds = kwargs.pop("n_pc_rounds", 4)
        
        self.OUTPUT_SEP = kwargs.pop("OUTPUT_SEP", " ")
        self.verbose = kwargs.pop("verbose", True)
        self.line_number = kwargs.pop("line_number", None)
        
        ########################
        # Read parameter file(s)
        # ----------------------
        
        if isinstance(param_input, str):
            # Read parameter directory
            if os.path.isdir(param_input):
                if self.verbose:
                    print(f"Reading parameter files from directory: {param_input}")
                self._param_dir = param_input
                self.read_param_dir(param_input)
        
            # Read parameter file
            if os.path.isfile(param_input):
                if self.verbose:
                    print(f"Reading parameter file: {param_input}")
                self.read_param_file(param_input)
                self._param_dir = os.path.dirname(param_input)
        
        # Read parameter file as a pandas dataframe or dict
        if isinstance(param_input, (pd.DataFrame, dict)):
            self.set_param_file(file_name, param_input)
        
    #####################
    # METHODS FOR READING
    # -------------------
    
    def read_param_dir(self, param_dir):
        files = [f for f in os.listdir(param_dir) if not f.startswith(".")]
        files = [f for f in files if not os.path.isdir(join(param_dir, f))]
        
        self.FILE_NAMES_CHIPSUPTAKE = [s for s in files if regex_chips.match(s)]
        self.FILE_NAMES_CHIPSUPTAKE = remove_popart_head_ext(self.FILE_NAMES_CHIPSUPTAKE)
        
        self.FILE_NAMES_PC = [s for s in files if regex_pc.match(s)]
        
        for f in files:
            self.read_param_file(join(param_dir, f))
        
        print("Read the following")
        self.summarise_params()
        print(self.misc.keys())
        print(self.params.keys())
    
    def read_param_file(self, param_file, file_name = None, patch = None):
        """Read a parameter file"""
        
        if not file_name:
            file_name = guess_file_name(param_file)
        
        if file_name in FILE_NAMES_LINELIST:
            self.read_file_whitespace(param_file, file_name, header = 0)
        
        if file_name in FILE_NAMES_DEMOGRAPHY:
            self.read_file_whitespace(param_file, file_name, header = None)
        
        if self.FILE_NAMES_CHIPSUPTAKE and (file_name in self.FILE_NAMES_CHIPSUPTAKE):
            self.read_file_whitespace(param_file, file_name, header = 0)
        
        if self.FILE_NAMES_PC and (file_name in self.FILE_NAMES_PC):
            self.read_file_raw(param_file, file_name)
        
        if file_name in FILE_NAMES_PATCH:
            self.read_patchinfo(param_file)
        
        if file_name in FILE_NAMES_MISC:
            self.read_file_raw(param_file, file_name)
    
    def read_test_parameters(self, parameter_transpose_file, mortality_file, 
        fertility_file, patchinfo_file, python_seed_file = None,
        chips_file = None, pc_file = None, **kwargs):
        """
        Read parameter test set as a consolidated set of files
        (as opposed to a processed dataset that the IBM can already use)
        """
        file_type_var = "File"
        param_name_var = "Name"
        param_value_var = "Value"

        df = pd.read_csv(parameter_transpose_file)
        files_grouped = df.groupby(file_type_var)
        
        for name, group in files_grouped:
            df_group = group.set_index(param_name_var)[[param_value_var]].T
            self.set_param_file(name, df_group)
        
        self.params["fertility"] = pd.read_csv(fertility_file, header = 0, index_col = 0, sep = ",")
        setattr(self, "fertility", self.params["fertility"])
        
        self.params["mortality"] = pd.read_csv(mortality_file, header = 0, index_col = 0, sep = ",")
        setattr(self, "mortality", self.params["mortality"])
        
        self.read_patchinfo(patchinfo_file)
        
        if python_seed_file:
            self.read_file_raw(python_seed_file, "python_seed.txt")
        
        # Create CHiPs schedules
        C = ChipsSchedule()
        
        self.FILE_NAMES_CHIPSUPTAKE = []
        for r in range(self.n_chips_rounds):
            rr = r + 1
            filename = f"chipsuptake_round{rr}"
            self.FILE_NAMES_CHIPSUPTAKE.append(filename)
            self.set_param_file(f"{filename}", C.schedule[f"chipsuptake_round{rr}"])
        
        # Add post-trial CHiPs round
        filename = "chipsuptake_roundposttrial"
        self.FILE_NAMES_CHIPSUPTAKE.append(filename)
        self.set_param_file(f"{filename}", C.schedule[f"chipsuptake_roundposttrial"])
        
        # Create PC schedule
        P = PCSchedule()
        
        self.FILE_NAMES_PC = []
        for r in range(self.n_pc_rounds):
            for p in range(self.n_patches):
                comm = self.patchinfo[0][p]
                filename = f"{FILE_NAME_PREFIX}_patch{p}_PC{r}_community{comm}.csv"
                self.FILE_NAMES_PC.append(filename)
                self.misc[f"{filename}"] = P.schedule[f"PC{r}"].to_string(index = False)
    
    def read_file_whitespace(self, param_file, file_name, header):
        self.params[file_name] = pd.read_csv(param_file, header = header, delim_whitespace = True)
        setattr(self, file_name, self.params[file_name])
    
    def read_file_raw(self, param_file, file_name):
        with open(param_file, 'r') as f:
            self.misc[file_name] = f.read()
            setattr(self, file_name, self.misc[file_name])
    
    def read_patchinfo(self, patchinfo_file):
        
        with open(patchinfo_file, 'r') as f:
            data = f.readlines()
        
        # Store list of ints in patch info file
        self.patchinfo = [[int(d) for d in l.split() ] for l in data]
    
    def find_pc_filename(self, p, filename):
        """
        Make PC filename based upon patch, and patch info
        """
        community_num = self.patchinfo[0][p]
        return f'{filename}_community{community_num}'
    
    #####################
    # METHODS FOR SETTING
    # -------------------
    def set_param(self, param_name, param_value, file_name = None, line_number = None):
        """
        Set parameter value or values within a file
        """
        if not file_name:
            file_names = self.find_param_name_file(param_name)
            if len(file_names) == 1:
                file_name = file_names[0]
            else:
                raise AmbiguousParamName(\
                    f"Pass file_name arg; Ambiguous or unknown parameter: {param_name}"\
                    )
        
        if line_number:
            self.params[file_name].at[line_number, param_name] = param_value
        else:
            self.params[file_name].loc[:, param_name] = param_value
    
    def set_param_file(self, file_name, df):
        """Set parameter file using input as dataframe or dict"""
        
        if isinstance(df, dict):
            df = pd.DataFrame(df)
        
        self.params[file_name] = df
        
        setattr(self, file_name, self.params[file_name])
    
    ###################################
    # METHODS FOR GETTING
    # ---------------------------------
    def get_param(self, param_name, file_name = None, line_number = None):
        """
        Returns a Series (of a particular type)
        """
        # Find which file this parameter name is from
        if not file_name:
            file_names = self.find_param_name_file(param_name)
            if len(file_names) == 1:
                file_name = file_names[0]
            else:
                raise AmbiguousParamName(\
                    f"Pass file_name arg; Ambiguous or unknown parameter: {param_name}"\
                    )
        
        if line_number:
            # FIXME: Perhaps check this is not beyond the end of the file
            return self.params[file_name][param_name].loc[[line_number]]
        else:
            return self.params[file_name][param_name]
    
    def get_param_file(self, file_name):
        return self.params[file_name]
    
    
    def summarise_params(self):
        """Print number of files and number of parameters per files"""
        files = sorted(self.list_files())
        n_files = len(files)
        
        total_cols = 0; total_rows = 0
        
        print(f"{n_files} files: {files}")
        
        print("\n")
        print('%(file_name)30s | %(rows)6s %(cols)6s' % \
            {'file_name': "File name", "rows": "# rows", "cols": "# cols"})
        for file_name in files:
            df = self.params[file_name]
            rows, cols = df.shape
            
            print('%(file_name)30s | %(rows)6d %(cols)6d' % \
                {'file_name': file_name, "rows": rows, "cols": cols})
            
            total_cols += cols; total_rows += rows
        
        print('%(total)30s | %(hline)12s' % {"total": "", "hline": "_"*13})
        print('%(total)30s | %(total_rows)6d %(total_cols)6d' % \
            {"total": "Total", "total_rows": total_rows, "total_cols": total_cols})
    
    
    def find_param_name_file(self, param_name):
        """Find the file name that a parameter name is associated with"""
        return [f for f in self.params.keys() if param_name in self.params[f].columns]
    
    def find_param_name(self, pattern):
        """Find pattern in parameter name"""
        return [p for f in self.params.keys() for p in self.params[f].columns if pattern in p]
    
    def dict_param_files(self):
        """Return a dict of parameter names"""
        return {k: list(v.columns.values) for (k, v) in self.params.items()}
    
    def list_files(self):
        """Return a list of parameter files names"""
        return list(self.params.keys())
    
    def list_params(self):
        """Return a flattened list of all parameters"""
        params_list = [self.params[file_name].columns.values for file_name in self.params.keys()]
        # Flatten list
        return [p for ps in params_list for p in ps]
    
    def list_params_linelist(self):
        """Return a flattened list of all parameters"""
        params_list = [self.params[file_name].columns.values for file_name in FILE_NAMES_LINELIST]
        # Flatten list
        return [p for ps in params_list for p in ps]
    
    #####################
    # METHODS FOR WRITING
    # -------------------
    def write_param_file(self, file_name, output_dir = ".", patches = None, **kwargs):
        
        if patches is None:
            patches = range(self.n_patches)
        
        if not hasattr(patches, "__iter__"):
            patches = list(patches)
        
        for patch in patches:
            
            if ( regex_pc_schedule.match(file_name) ):
                
                fname = self.find_pc_filename(patch, file_name)
                
                output_filename = join(output_dir,
                    f"{FILE_NAME_PREFIX}_patch{patch}_{fname}.csv")
            else:
                output_filename = join(output_dir,
                    f"{FILE_NAME_PREFIX}_patch{patch}_{file_name}.csv")
            
            self.params[file_name].to_csv(output_filename, \
                sep = self.OUTPUT_SEP, index = False, **kwargs)
    
    def write_param_file_raw(self, output, output_filename):
        """Write a file to disk with minimal processing"""
        
        with open(output_filename, 'w') as f:
            f.write(output)
    
    def write_patchinfo_file(self, output_filename):
        """Write patchinfo file to disk"""
        pd.DataFrame(self.patchinfo).to_csv(output_filename, 
            sep = " ", index = False, header = False)
    
    def write_param_set(self, output_dir = "."):
        """Write complete parameter set to a folder"""
        
        # Create output_dir if it doesn't already exist
        Path(output_dir).mkdir(parents = True, exist_ok = True)
        
        for f in FILE_NAMES_LINELIST:
            self.write_param_file(file_name = f, output_dir = output_dir)
        
        for f in FILE_NAMES_DEMOGRAPHY:
            self.write_param_file(file_name = f, output_dir = output_dir, header = False)
        
        for f in FILE_NAMES_PATCH:
            self.write_patchinfo_file(output_filename = join(output_dir, f))
        
        for f in self.FILE_NAMES_CHIPSUPTAKE:
            self.write_param_file(file_name = f, output_dir = output_dir)

        for f in self.FILE_NAMES_PC:
            self.write_param_file_raw(
                output = self.misc[f], 
                output_filename = join(output_dir, f))
        
        for f in FILE_NAMES_MISC:
            if f in self.misc:
                self.write_param_file_raw(
                    output = self.misc[f],
                    output_filename = join(output_dir, f))
    
    @property
    def param_dir(self):
        return self._param_dir


def remove_popart_head_ext(string):
    """
    Remove FILE_NAME_PREFIX associated with PopART-IBM input parameters and 
    file extension from string representation of a filename
    
    Arguments
    ---------
    string : str, list of str
        File name as a string (e.g. param_processed_patch1_demographics.csv)
    
    Returns
    -------
    
    
    """
    
    if not isinstance(string, list):
        string = [string]
    
    output_list = []
    for s in string:
        # Remove file extension
        sub = re.split(r"\.", s)[0]
        
        # Remove PopART-IBM related head
        sub = re.split(f'{FILE_NAME_PREFIX}_patch._', sub)[1]
        
        output_list.append(sub)
    
    if len(output_list) == 1:
        output_list = output_list[0]
    
    return output_list


def guess_file_name(param_file):
    """Attempt to guess file name if one is not provided"""
    basename = os.path.basename(param_file)
    
    file_name_matches_patch = [f for f in FILE_NAMES_PATCH if f in basename]
    if len(file_name_matches_patch) == 1:
        return file_name_matches_patch[0]
    
    file_name_matches_misc = [f for f in FILE_NAMES_MISC if f in basename]
    if len(file_name_matches_misc) == 1:
        return file_name_matches_misc[0]
    
    if regex_chips.match(basename):
        return remove_popart_head_ext(basename)
    
    if regex_pc.match(basename):
        return basename
    
    file_name_matches_demography = [f for f in FILE_NAMES_DEMOGRAPHY if f in basename]
    if len(file_name_matches_demography) == 1:
        return file_name_matches_demography[0]
    
    file_name_matches = [f for f in FILE_NAMES_LINELIST if f in basename]
    if len(file_name_matches) == 1:
        file_name = file_name_matches[0]
    else:
        raise AmbiguousFileName(\
            f"Pass file_name argument; Ambiguous or unknown file name within {param_file}"\
            )
    return file_name


####################################################
# Parameters related to creation of CHiPs schedules
# --------------------------------------------------

NTIMESTEP = 48
SEXES = ['M', 'F']

class SamplingSchedule(object):
    """
    Generic class representing an age/sex stratified sample scheme for the PopART-IBM
    """
    def __init__(self, **kwargs):
        self.n_rounds = kwargs.pop("n_rounds", None)
        self.n_timestep_per_year = kwargs.pop("n_timestep_per_year", NTIMESTEP)
        self.starting_times = kwargs.pop("starting_times", None)
        self.ending_times = kwargs.pop("ending_times", None)
        self.column_name_prefixes = kwargs.pop("column_name_prefixes", "")
        self.ages = kwargs.pop("ages", None)
        self.sexes = kwargs.pop("sexes", SEXES)
        
        self.timestep = 1./self.n_timestep_per_year
        
        # Dict to store the schedules
        self.schedule = dict()
    
    def create_uniform_schedule(self, 
        start_year_timestep = None, 
        end_year_timestep = None, 
        final_coverage = None, 
        starting_coverage = 0):
        """
        Create a uniform schedule for a single round of age/sex sampling for use with PopART-IBM
        
        Arguments
        ---------
        start_year_timestep : list of ints: [year, timestep]
            Ending year and ending timestep
        
        end_year_timestep : list of ints: [year, timestep]
            Ending year and ending timestep
        
        final_coverage: float
            final coverage across all age/sex
        
        starting_coverage
        """
        
        start_year, start_timestep = start_year_timestep
        end_year, end_timestep = end_year_timestep
        
        # Create an array of all time steps over the years in question
        years = np.arange(start_year, end_year + 1)
        arr = np.meshgrid(years, np.arange(self.n_timestep_per_year))
        
        # Rearrange into [year, timestep] form
        arr = np.vstack(np.array(arr).T)
        
        # Subset to just the timesteps/years of interest
        istart = start_timestep
        iend = ((len(years)-1)*self.n_timestep_per_year + end_timestep + 1)
        year_timestep = arr[istart:iend]
        
        # Create "decimal" year from year and timestep
        year = year_timestep[:, 0] + year_timestep[:, 1]/float(self.n_timestep_per_year)
        
        # Create headers
        cov_header = [pref + sex + str(age) \
            for pref in self.column_name_prefixes \
            for sex in self.sexes \
            for age in self.ages]
        
        ncol = len(cov_header)
        
        # Split coverage into a per-time step proportion
        coverage = np.diff(np.linspace(starting_coverage, final_coverage, len(year)+1))
        coverages = np.tile(coverage, (ncol, 1))
        
        # Stack output in required form
        output = np.hstack([year[:, None], year_timestep, coverages.T])
        
        time_header = ['time', 'year', 't_step']
        
        # Return dataframe
        df = pd.DataFrame(output)
        df.columns = time_header + cov_header
        
        # Convert data types
        df[["year", "t_step"]] = df[["year", "t_step"]].astype(int)
        
        return df


N_CHIPS_ROUNDS = 3
CHIPS_STARTING_TIMES = [[2013, 45], [2015, 24], [2016, 34]]
CHIPS_ENDING_TIMES = [[2015, 23], [2016, 31], [2017, 46]]
CHIPS_FINAL_COVERAGES = [0.6, 0.7, 0.8]
chips_ages = range(18, 81)

chips_column_name_prefixes = ["prop_hiv_status_known_3m"]
prefix = "chipsuptake_round"

class ChipsSchedule(SamplingSchedule):
    """
    Class representing a CHiPs sampling schedule
    CHiPs schedules are expressed as 'fraction of age/sex group visited in a time step'
    
    """
    def __init__(self, **kwargs):
        self.n_rounds = kwargs.pop("n_rounds", N_CHIPS_ROUNDS)
        self.n_timestep_per_year = kwargs.pop("n_timestep_per_year", NTIMESTEP)
        self.starting_times = kwargs.pop("starting_times", CHIPS_STARTING_TIMES)
        self.ending_times = kwargs.pop("ending_times", CHIPS_ENDING_TIMES)
        self.column_name_prefixes = kwargs.pop("column_name_prefixes", chips_column_name_prefixes)
        self.ages = kwargs.pop("ages", chips_ages)
        self.sexes = kwargs.pop("sexes", SEXES)
        
        self.timestep = 1./self.n_timestep_per_year
        
        # Dict to store the schedules
        self.schedule = dict()
        
        self.final_coverages = kwargs.pop("chips_final_coverages", CHIPS_FINAL_COVERAGES)
        self.create_chips_schedule()
    
    def create_chips_schedule(self):
        """
        Create a CHiPs schedule for all rounds + post trial rounds
        """
        
        # Loop through all rounds
        for r in range(self.n_rounds):
            start_year_timestep = self.starting_times[r]
            end_year_timestep = self.ending_times[r]
            
            df_chips_schedule = self.create_uniform_schedule(
                start_year_timestep, 
                end_year_timestep, 
                final_coverage = self.final_coverages[r])
            
            rr = r + 1
            self.schedule[f'{prefix}{rr}'] = df_chips_schedule
        
        # Find post-trial CHiPs coverage
        df_chips_schedule_posttrial = self.create_uniform_schedule(
            start_year_timestep = [2018, 0], 
            end_year_timestep = [2018, 47], 
            final_coverage = self.final_coverages[r])
        
        self.schedule[f'{prefix}posttrial'] = df_chips_schedule_posttrial


#################################################
# Parameters related to creation of PC schedules
# -----------------------------------------------

N_PC_ROUNDS = 4
PC_STARTING_TIMES = [[2014, 5], [2015, 27], [2016, 35], [2017, 24]]
PC_ENDING_TIMES = [[2015, 11], [2016, 22], [2017, 19], [2018, 23]]
PC_FINAL_COVERAGES = [0, 0, 0, 0]

pc_ages = np.arange(18, 45)

pc_column_name_prefixes = [
    "n_enrolled_in_PC0_with_HIV_test_HIVneg", 
    "n_enrolled_in_PC0_with_HIV_test_HIVposAware",
    "n_enrolled_in_PC0_with_HIV_test_HIVposUnaware"
]

class PCSchedule(SamplingSchedule):
    """
    Class representing a PC schedule for use within the PopART-IBM
    PC schedules are expressed as 'fraction of age/sex group visited in a time step'
    """
    def __init__(self, **kwargs):
        self.n_rounds = kwargs.pop("n_rounds", N_PC_ROUNDS)
        self.n_timestep_per_year = kwargs.pop("n_timestep_per_year", NTIMESTEP)
        self.starting_times = kwargs.pop("starting_times", PC_STARTING_TIMES)
        self.ending_times = kwargs.pop("ending_times", PC_ENDING_TIMES)
        self.column_name_prefixes = kwargs.pop("column_name_prefixes", pc_column_name_prefixes)
        self.ages = kwargs.pop("ages", pc_ages)
        self.sexes = kwargs.pop("sexes", SEXES)
        
        self.timestep = 1./self.n_timestep_per_year
        
        # Dict to store the schedules
        self.schedule = dict()
        
        self.final_coverages = kwargs.pop("pc_final_coverages", PC_FINAL_COVERAGES)
        self.create_pc_schedule()
    
    def create_pc_schedule(self):
        """
        Create a PC schedule for each round
        """
        
        # Loop through all rounds
        for r in range(self.n_rounds):
            
            start_year_timestep = self.starting_times[r]
            end_year_timestep = self.ending_times[r]
            
            df_pc_schedule = self.create_uniform_schedule(
                start_year_timestep, 
                end_year_timestep, 
                final_coverage = self.final_coverages[r])
            
            # Convert to integers
            cols = df_pc_schedule.columns
            cols = [c for c in cols if c != "time"]
            df_pc_schedule[cols] = df_pc_schedule[cols].astype(int)
            
            self.schedule[f'PC{r}'] = df_pc_schedule


