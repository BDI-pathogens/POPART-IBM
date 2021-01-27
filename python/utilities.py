#!/usr/bin/env python3
"""
Basic functions used in several Python files
"""

import os, re, sys

###########################################################################
# Basic debugging functions:
###########################################################################

def handle_error(err_messg):
    """
    Error-handling function
    
    If an error occurs, print error message and exit program.
    """
    print("---------")
    print("Called from: ", __file__)
    print(err_messg)
    print("---------")
    sys.exit(1)



def check_directory_exists(dir_to_check):
    """
    Checks if a directory exists or not. 
    
    If it does exist, do nothing. 
    If not, print an error message and exit.
    """
    if not(os.path.isdir(dir_to_check)):
        handle_error("Error: "+dir_to_check+" does not exist. Please check code. Exiting.")


def check_directory_exists_and_create(dir_to_check):
    """
    Checks if input directory exists. 
    If it does exist, do nothing, if not, create it.
    """
    if not(os.path.isdir(dir_to_check)):
        os.makedirs(dir_to_check)


########################################################
########### Generic functions used in code   ###########
########################################################


def merge_two_dicts(x, y):
    """
    Given two dicts, merge them into a new dict as a shallow copy.
    
    If y has any keys in common with x, then only the x values for those keys are kept.
    
    
    Parameters
    ----------
    x, y : dict
    
    
    Returns
    -------
    z : dict
    """
    z = x.copy()
    z.update(y)
    return z    


def add_comment(original_text, comment):
    """
    Append a C-style comment to a text string.
    
    
    Parameters
    ----------
    original_text : str
        Original text string
    
    comment : str
        Comment to be added after a double slash '//'
    
    
    Returns
    -------
    new_text : str
    """
    # Only need to add comment if comment is non-empty:
    if comment:
        new_text = original_text+ " // "+comment
    else:
        new_text = original_text
    return new_text


def run_sweave(sweave_filename,sweave_dir):
    """
    Run sweave to update data files, if needed.
    
    
    Parameters
    ----------
    sweave_filename : str
        Sweave filename (either within the current directory or sweave_dir)
    
    sweave_dir : str
        Directory in which the sweave file is located (if not in current directory)
    
    
    Returns
    -------
    Nothing returned.
    
    """
    currentdir = os.getcwd()

    # Move to the directory containing the sweave file (or print an error if that folder doesn't exist):
    try:
        os.chdir(sweave_dir)
    except:
        handle_error("Directory "+sweave_dir+" does not exist in function run_sweave(). Exiting\n")

    # Run Sweave (note this generates other files - e.g. the pdf - as well as the text files we need):
    os.system("R CMD Sweave --pdf "+sweave_filename)

    # Now go back to the original working directory:
    os.chdir(currentdir)


######### Functions which relate to reading in and processing files: ####################


def parse_file(f):
    """
    Read text file, split by line, return a list of lines.  
    
    
    Parameters
    ----------
    f : str
        Filename of the file to read
    
    
    Returns
    -------
    data : list
        Each element in the list is a line (as a str) of the input file
    """
    infile = open(f,"r")
    data = infile.read().rstrip().split("\n")
    infile.close()
    return data


def parse_line(l):
    """
    Split a string from a line of a parameter file into parameter, value, comment.
    
    Function takes a line from any parameter file of the form:
        "start_time_simul 1900   // Assumption - around 20 years for burn-in."
    and pulls out the parameter name (e.g. start_time_simul), the value (or range etc) - here 
    1900 - and any comments on that line (e.g. "// Assumption - around 20 years for burn-in.").
    Almost all param files have this format.  
    
    
    Parameters
    ----------
    l : str
        String of a line in a parameter file
    
    
    Returns
    -------
    A list of [paramname, paramvaluelist, comment]
    
    paramname : str
        The first thing in the input string is the parameter name
    paramvaluelist : str
        Second thing in the input string is the parameter value
    comment : str
        After the C-style comment (//) is the comment.  
    
    
    Example
    -------
    line = "eff_circ_vmmc 0.6  //  Cori 2013 Table S7"
    [paramname, paramvaluelist, comment] = parse_line(line)
    
    >> paramname
    >> 'eff_circ_vmmc'
    >> paramvaluelist
    >> ['0.6']
    >> comment
    >> '  Cori 2013 Table S7'
    
    line = "start_time_simul 1900   // Assumption - around 20 years for burn-in."
    [paramname, paramvaluelist, comment] = parse_line(line)
    
    >> paramname
    >> 'start_time_simul'
    >> paramvaluelist
    >> ['1900']
    >> comment
    >> ' Assumption - around 20 years for burn-in.'
    """
    # First remove any comments (ie anything starting with "//"
    lineparam = l.split("//")[0]
    
    # Now keep the comments
    try:
        # take everything after the first // and store as the string 'comment'. 
        # What I actually do is split the line every time there is a // into an array, 
        # then slice off the first element of the array and stick it back together with //s 
        # joining the elements of the array.
        comment = "//".join(l.split("//")[1:])
    except:
        comment = "" # If no "//" then there is no comment on this line.
        
    # First thing is the parameter name, followed by the value(s) of the parameter, with
    #  spaces/tabs (ie whitespace) separating items.
    # split() divides up a string by whitespace, making each non-whitespace thing a separate item
    # in a list  e.g. if s="a  bcd ef" then s.split()=["a","bcd","ef"]
    paramname = lineparam.split()[0] 
        
    # Divides up the line by whitespace, then make a list containing everything from the second
    # item onwards.  
    paramvaluelist = lineparam.split()[1:] 
    return [paramname,paramvaluelist,comment]


def remove_extra_whitespace(s):
    """
    Remove extra whitespace (spaces, tabs) from a string, return resultant string
    
    The re.sub() replaces the first argument with the second argument in s. ' +' means "one or more
    whitespace" so this replaces multiple whitespace with a single space.
    """
    s = s.rstrip() # Remove trailing whitespace.
    s = re.sub('\s\s+',' ',s)
    s = re.sub('\t+',' ',s)
    return s


def write_file(outfilename, outstring):
    """
    Write string `outstring` to file `outfilename` (removing trailing whitespace using rstrip()).
    """
    outfile = open(outfilename, "w")
    outfile.write(outstring.rstrip())
    outfile.close()


def create_markdown_from_df(df, title = None):
    """
    Create text string of markdown table from pandas dataframe.  
    Used in automated creation of markdown tables for model documentation.

    Arguments
    ---------
    df : pandas.DataFrame
        Dataframe for which one wants to generate a markdown table
    title : str
        Title string for header of the markdown file

    Returns
    -------
    str of markdown table of the form (columns of the table are column names of the pandas df):
    
    # <title>
    | Name | .... | Source |
    | ---- | ---- |  ----  |
    |  11  | .... |   1N   |
    |  21  | .... |   2N   |
    | .... | .... |  ....  |
    |  M1  | .... |   MN   |
    """
    
    NCOLS = df.shape[1]
    
    output = list()
    
    if title:
        output += ["# " + title]
    
    # Add header
    output += ["| " + " | ".join(df.columns) + " | "]
    output += ["| " + "".join([" ---- |" for i in range(NCOLS)])]
    
    for i, row in df.iterrows():
        output += ["| " + ' | '.join('{0}'.format(el) for el in row) + " |"]

    return("\n".join(output))
