#!/bin/bash 

 
#SBATCH -A fraser.prj 
#SBATCH -J popart-ibm 

#SBATCH -o popart-ibm.out 

#SBATCH -e popart-ibm.err 
#SBATCH -p short
module unload GSL
module load GSL/2.6-GCC-8.3.0
module load intel
starting_run=1
number_of_runs=1000
number_of_sims=28
input_dir="/well/fraser/users/zkv138/PARAMETERS/POSTERIORS"

for community in {1..12}; do
        $(pwd)/popart-simul.exe ${input_dir}/PARAMS_COMMUNITY${community}_ACCEPTED $number_of_runs 0 $starting_run $number_of_sims $input_dir/Outputs
done

#get rid of anything non phylo_related

rm ${input_dir}/Outputs/CHIPS*
rm ${input_dir}/Outputs/Annual*
