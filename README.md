IBM for the PopART trial
========================

Description
-----------

PopART-IBM an individual-based model for simulating HIV epidemic in high-prevalence settings, as used in the HPTN 071 (PopART) trial.  A full description of the model is described in [Pickles et al., 2020](https://www.medrxiv.org/content/10.1101/2020.08.24.20181180v1).  A data dictionary describing all the output output files from the model is [here](doc/output_files/output_file_dictionary_overview.md) and a data dictionary describing the parameters in the model is described [here](doc/parameters/parameters.md).  


Compilation
----------------

For Mac and Unix-type systems, PopART-IBM requires a C compiler (such as gcc) and the [GSL](https://www.gnu.org/software/gsl/) libraries installed:

```bash
cd IBM_simul/src
make all
```

GSL can be downloaded from [here](ftp://ftp.gnu.org/gnu/gsl/).  


For Windows systems, please see [this](./doc/running_popartibm_on_windows.md) walkthrough.  

Usage
-----

```bash
cd IBM_simul/src
./popart-simul.exe <inputdir> <nruns>
```
 
 where
 
* `inputdir`: Directory where input parameter files ("param_processed*.csv") are located
* `nruns` : number of simulation runs in parameter files (num. of lines in parameter files to read in)

A basic [example](examples/example_101.py) illustrates how the parameter input files can be set up (steps 1 and 2) as is expected by the model.  

**Notes**

* Additional command-line arguments are described in [main.c](src/main.c).  
* The model will write all output files to the directory `inputdir/Output` (additional command-line arguments can adjust this).  
* The output files written will depending upon which macros are set to 1 within the file [constants.h](src/constants.h) (those beginning `WRITE_*`).  


Testing
-------

Tests are written [`pytest`](https://docs.pytest.org/en/stable/) using Python v3.6+, and run in the following manner: `python3 -m pytest`.  Some tests take a long time to run and so are not run by default, they can be invoked using the ` --runslow` option in `pytest` (`python3 -m pytest --runslow`).  

It is recommended that tests are run under a Python virtual environment.  The following will set up a Python virtual environment, install required modules, and run the tests: 

```
python3 -m venv venv
source venv/bin/activate
python3 -m pip install -r tests/requirements.txt
python3 -m pytest
deactivate
```

Contributing
------------

Contributions are most welcome.  Please see the documentation on [contributing](CONTRIBUTING.md) for further information.  Please contact the core team or raise an issue before making a pull request.  
