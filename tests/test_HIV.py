#!/usr/bin/env python3
"""
Tests of PopART-IBM associated with HIV using pytest. 

Several checks are performed by adjusting input parameters and comparing to expected outputs in
the `Annual_outputs*.csv` files.  A single simulation is checked (and only the intervention patch).

W. Probert, 2019
"""

import subprocess, shutil, os, pytest
from os.path import join
import numpy as np, pandas as pd

from utils import adjust_macros
import constants as c

class TestClass(object):
    def test_hiv_seeds(self, basic_model_output):
        """
        Test that HIV is seeded and a generalised (>10%) epidemic in 2010 is simulated
        """
        # Read model outputs
        df_annual, df_timestep = basic_model_output
        
        # Check if there were any HIV cases
        np.testing.assert_equal(df_annual.Prevalence[df_annual.Year == 2010].values > 0.1, 1)
