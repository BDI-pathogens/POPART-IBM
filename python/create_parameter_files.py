#!/usr/bin/env python3
"""
Script to create markdown tables of parameter files
from the test parameter set
"""

import sys, numpy as np, pandas as pd
from os.path import join
from utilities import create_markdown_from_df

input_file = "tests/data/parameters_transpose.csv"

file_type_var = "File"
param_name_var = "Name"
param_value_var = "Value"
param_file_stub = "param_processed"

n_patches = 2
SEP = " "

if __name__ == "__main__":
    
    ##########################
    # Make markdown files
    # for each parameter file
    # -----------------------

    df = pd.read_csv(input_file)

    files_grouped = df.groupby(file_type_var)

    # Make a csv file for each file type with column names as "Name"
    for name, group in files_grouped:
        group.Name = "`" + group.Name + "`"
        group = group.fillna("-")
        markdown_table = create_markdown_from_df(group, 
            title = f"Table: PopART-IBM parameter dictionary for {name} file ")
    
        output_file = f"doc/parameters/parameters_{name}.md"

        with open(output_file, 'w') as f:
            f.write(markdown_table)

    # Create markdown file of all parameter names
    df.Name = "`" + df.Name + "`"
    df = df.fillna("-")
    markdown_table = create_markdown_from_df(df, 
        title = "Table: PopART-IBM parameter dictionary")

    output_file = "doc/parameters/parameters.md"

    with open(output_file, 'w') as f:
        f.write(markdown_table)
