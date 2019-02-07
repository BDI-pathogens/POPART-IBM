IBM for the PopART trial
========================
Repository for the C code used to run the individual-based model for the PopART trial.  

Compilation
----------------
In order to run the model, the C code needs to first be compiled, linked, and an executable (`popart-simul.exe`) created.  

```
cd src
make clean; make all
cd ..
```


Running the model
-----------------

Once compiled, the model can be run in the following manner:

```bash
output_dir="./data/PARAMS"  # Folder where parameter files are kept
nruns=10                    # Number of lines in the parameter files
CF=0                        # Should a counterfactual be run (0) or not (1)

# Make a folder for output
mkdir -p $output_dir/Output

./src/popart-simul.exe $output_dir $nruns $CF
```


Output from the simulation
--------------------------

Output will be stored in `$output_dir/Output`.  Several files will be output depending upon which macros are switched 'on' (set to 1) within the file [constants.h](src/constants.h).  

