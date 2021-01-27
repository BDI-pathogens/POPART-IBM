#!/usr/bin/env python3
"""
Utilities associated with tests of the PopART-IBM model

W. Probert, 2020
"""

def adjust_macros(input_file, macro_name, macro_value, output_file):
    """
    Adjust macros within the constants.h file of the IBM to produce different output.  
    """
    f = open(input_file, 'r')
    data = f.readlines()
    f.close()
    
    # Split on white space
    d = [[i, line] for i, line in enumerate(data) if line.startswith('#define ' + macro_name)]
    print("Found", macro_name, "at line", d[0][0])
    line = d[0][1]
    line = line.split()
    new_line = " ".join([line[0], line[1], macro_value, "\n"])
    
    # Add the line into the output data ... 
    for i, line in enumerate(data):
        if i == d[0][0]:
            data[i] = new_line
            print(data[i])
    
    # Write to file ... 
    f = open(output_file, 'w')
    f.writelines(data)
    f.close()
