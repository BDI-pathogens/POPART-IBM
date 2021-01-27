#!/usr/bin/env python3
"""
Script to create markdown table from csv file
"""

import sys
from os.path import join
import pandas as pd, numpy as np
from utilities import create_markdown_from_df

if __name__ == "__main__":

    # Parse command line arguments
    output_file_csv = sys.argv[1]   # "output_file_dictionary_overview.csv"
    output_file = sys.argv[2]       # "output_file_dictionary_overview.md"

    df = pd.read_csv(output_file_csv, dtype = str)

    # Generate table for all parameters
    df_all = df.replace(np.nan, "-")
    markdown_table = create_markdown_from_df(df_all, 
        title = "Table: PopART-IBM output file dictionary")
    
    with open(output_file, 'w') as f:
        f.write(markdown_table)
