/*  This file is part of the PopART IBM.

    The PopART IBM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The PopART IBM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the PopART IBM.  If not, see <http://www.gnu.org/licenses/>.
 */

/************************************************************************/
/******************************* Includes  ******************************/
/************************************************************************/

#include "constants.h"
#include "structures.h"
#include "input.h"
#include "init.h"
#include "utilities.h"
#include "partnership.h"
#include "checks.h"
#include "demographics.h"
#include "output.h"
#include "memory.h"
#include "fitting.h"
#include "debug.h"
#include "simul.h"

/************************************************************************/
/******************************* main function **************************/
/************************************************************************/

int main(int argc,char *argv[]){
    /*
    Main function
    
    Arguments
    ---------
    argc: int
        Number of arguments fed in the command line
    argv: char array
        Pointer to an array containing the arguments fed in the command line.  Currently the command
        line argument are the follwing:
        <inputdir> <nruns> <counterfactual> <i_startrun> <outputdir> <seed_offset> <PC_seed_offset>
        
        inputdir : str
            Directory where input parameter files ("param_processed*.csv") are located
        nruns : int
            Number of simulation runs in param files (num. of lines in parameter files to read in)
        counterfactual : int
            Should a counterfactual be run?  1 for yes, 0 for no
        i_startrun : int
            Simulation number from which to start from.  
        n_startrun : int
            Number of simulations to run (starting from i_startrun).  
        outputdir : str
            Directory where output files are stored (defaults to `inputdir`/Output/)
        seed_offset : int
            Offset to be used for the integer used to seed the GSL random number generator.  
        PC_seed_offset
            Offset to be used for the integer used to seed the PC sampling process.  
    */
    
    /* Check number of arguments */
    if(argc < 2){
        printf("Error. Execution syntax is the following:");
        printf("./popart-simul.exe <inputdir> <nruns> <counterfactual> <i_startrun> ");
        printf("<outputdir> <seed_offset> <PC_seed_offset>\n");
        printf("Exiting.\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if(VERBOSE_OUTPUT == 1){
        printf("-------------------------------------\n");
    }

    /*********************************************************/
    /*** Declaration of variables                          ***/
    /*********************************************************/

    int year, i, ir;
    
    int i_run; // Index for run number (ie we use parameter set 'i_run'), 
    /* i_run runs from 0..n_runs-1; so i_run is the index in allrunparameters[]
    referring to the current parameter set */
    
    int n_runs; // Total number of parameter sets used (read into memory)
    int i_startrun; // the first parameter set used (i_startrun); allows picking a single run
    int n_startrun; // the number of parameter sets used (i.e. run from i_startrun to n_startrun)
    int i_dhs_round; /* To let us know if/when we need to store data for DHS. */

    /* This stores the return value of carry_out_processes() - it is 1 if there was no fitting or
    if the run passed all the fits during the year, and 0 if there was a fit criterion which was
    not satisfied. This allows us to see if a run should be terminated early or not. */
    int fit_flag;
    int p; /* Index over patches. */
    int g, a, r;
    
    /* is_counterfactual = 1 if this is a counterfactual run (ie all patches set to be arm C) and
    is_counterfactual=0 if this is not a counterfactual run (ie it is a normal run). */
    int is_counterfactual;

    /* patch is an array of structs that each contain a single geographical unit (ie a "patch").
    
    A patch may be a trial community, district, region, or any other geographical unit.
    People in a single patch only interact with other patches in the following way:
     - external partnerships with people in other patches.
     - (eventually via migration between patches).
    */
    patch_struct *patch;
    all_partnerships *overall_partnerships;
    
    /* Contains pointers to all individuals in a cluster age AGE_ADULT or above, classified by
    gender, age and risk group */
    population *pop = NULL;

    /* Population: */
    /* This is the number of HIV+ individuals currently alive - use it to decide if HIV has died
    out in the population (for screening runs to stop early if have no HIV). */
    int n_infected_total[NPATCHES];
    
    /* At present we just set up this pointer. The memory for this array is allocated in
    load_fitting_data() separate to everything else. */
    fitting_data_struct *fitting_data[NPATCHES];

    /* Variables used in PANGEA comparison exercise: */
    long DEBUGNPARTNERSHIPS[NPATCHES];
    long DEBUGNSERODISCORDANTPARTNERSHIPS[NPATCHES];
    long DEBUGHIVNEGWITHHIVPOSPARTNERS[NPATCHES];
    long NACUTE[NPATCHES];
    long NCHRONIC[NPATCHES];
    
    /* Change the C random seed at the start of PopART if want to look at stochasticity with the
    same parameters. */
    int rng_seed_offset;
    
    /* Change the C random seed at the start of PC recruitment if want to look at stochasticity of
    PC recruitment with the same parameters. */
    int rng_seed_offset_PC;

    /* Read in command line arguments: */

    /* As argv[1] is a pointer, we set the pointer *input_file_directory to point at the same
    thing. This is a pre-existing thing that is being pointed at, so don't need to malloc (which
    would create a "blank" space for us to point at). */
    char *input_file_directory = argv[1];
    char *output_file_directory;
    
    /* Similarly to argv[1], it is easier to assign the pointer *output_file_directory here than in
    parse_command_line_arguments() below. */
    if(argc > 6){
        output_file_directory = argv[6];
    }else{
        /* Make empty if no argument (so current directory in linux/OS X). */
        output_file_directory = (char *)calloc(LONGSTRINGLENGTH, sizeof(char));
        /* Copy the output directory to output_filename. We take 1 from LONGSTRINGLENGTH as we add
        a slash after. */
        strcpy(output_file_directory, input_file_directory);
        add_slash(output_file_directory);
        join_strings_with_check(output_file_directory, "Output", LONGSTRINGLENGTH - 1,
            "'Output' and output_file_directory in main()");
    }
    
    /* Now parse any other command line arguments: */
    parse_command_line_arguments(argc, argv, &n_runs, &i_startrun, &n_startrun, &is_counterfactual,
        &rng_seed_offset, &rng_seed_offset_PC);
    
    /* Used in the output file name. */
    long python_rng_seed = get_python_seed(input_file_directory);

    /* Set up GSL variables for random numbers, and the values for n_runs and i_startrun */
    gsl_rng_env_setup();

    /* Set up gsl seed: */
    TYPE_RNG = gsl_rng_default;
    rng = gsl_rng_alloc (TYPE_RNG);

    /* Dummy variable used when storing random_seed in the output file string. */
    char tempstring[100];

    /* a pointer to a structure which will store items used for debuging */
    debug_struct *debug;

    /*********************************************************/
    /*** End of variable declarations                      ***/
    /*********************************************************/

    /*********************************************************/
    /* Parameters */
    /*********************************************************/
    /* We allow the input files to have multiple parameter sets. The array allrunparameters will
    contain all these parameter sets. Every time we want a new parameter set we take out the
    relevant array element from allrunparameters and copy it into param (actually we set the
    pointer for param to point at allrunparameters ie param = allrunparameters+i_run. */
    
    parameters *allrunparameters[NPATCHES];
    for(p = 0; p < NPATCHES; p++){
        
        allrunparameters[p] = malloc(n_runs * sizeof(parameters));
        
        if(allrunparameters[p] == NULL){
            printf("Unable to allocate allrunparameters in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        for(i_run = 0; i_run < n_runs; i_run++){
            
            allrunparameters[p][i_run].chips_params = malloc(sizeof(chips_param_struct));
            allrunparameters[p][i_run].PC_params = malloc(sizeof(PC_param_struct));
            allrunparameters[p][i_run].DHS_params = malloc(sizeof(DHS_param_struct));
            
            if(
            (allrunparameters[p][i_run].chips_params == NULL) ||
            (allrunparameters[p][i_run].PC_params == NULL) ||
            (allrunparameters[p][i_run].DHS_params == NULL)
            ){
                printf("Unable to allocate allrunparameters.chips_params/PC_params. ");
                printf("Execution aborted.");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }


    /*********************************************************/
    /*** Allocation of memory                              ***/
    /*********************************************************/

    int pc_enrolment_round, pc_size;
    patch = malloc(NPATCHES*sizeof(patch_struct));

    alloc_patch_memoryv2(patch);

    overall_partnerships = malloc(sizeof(all_partnerships));
    alloc_partnership_memoryv2(overall_partnerships);

    alloc_pop_memory(&pop, n_runs);

    /* Allocate PC memory. */
    for(pc_enrolment_round = 0; pc_enrolment_round < NPC_ENROLMENTS; pc_enrolment_round++){
        if(pc_enrolment_round == 0){
            
            /* Assume that PC0 is never more than 4000 in a community. */
            pc_size = 4000;
            
        }else if(pc_enrolment_round == 1){
            pc_size = 600;
        }else if(pc_enrolment_round == 2){
            pc_size = 600;
        }
        for(p = 0; p < NPATCHES; p++){
            alloc_pc_cohort_data(&patch[p].PC_cohort_data, pc_enrolment_round, pc_size);
        }
    }

    /* allocation of memory for debug variables */
    /* everything inside debug_struct is allocated statically so no need to allocate more memory */
    debug = malloc(sizeof(debug_struct));

    /*****************************************************/
    /*** READING PARAMETERS AND FITTING DATA           ***/
    /*****************************************************/

    /* Reads in all the parameters. Note that certain patch info (e.g. trial arm) is set here.
    In particular trial arm may be overwritten below if we are in a counterfactual scenario
    (is_counterfactual==1). */

    /* Read all the parameter sets - there should be n_runs of them. */
    read_param(input_file_directory, allrunparameters, n_runs, patch);
    
    for(i = 0; i < n_runs; i++){
        for(p = 0; p < NPATCHES; p++){
            for(ir = 1; ir < NCHIPSROUNDS; ir++){
                // Copy p_popart_to_cascade and p_circ_popart from first round to all rounds
                allrunparameters[p][i].p_popart_to_cascade[ir] =
                    allrunparameters[p][i].p_popart_to_cascade[0];
                
                allrunparameters[p][i].p_circ_popart[ir] = allrunparameters[p][i].p_circ_popart[0];
            }
        }
    }
    
    if(is_counterfactual == IS_COUNTERFACTUAL_RUN){
        for(p = 0; p < NPATCHES; p++){
            patch[p].trial_arm = ARM_C;
        }
        printf("Making this run a counterfactual\n");
    }

    /* Get country: country_setting values are 1 for Zambia and 2 for South Africa.
    There are constants defined in constants.h for this (ZAMBIA and SOUTH_AFRICA)
    This variable can be modified to be the cluster id (for example) if needed.
    Note we assume that the country setting does not vary across all the parameter sets.
    If this is not the case then we need to move this into the "for i_run" loop.
    country_setting determines demographics only (fertility+death rate). */
    get_setting(patch);

    file_label_struct *file_labels;
    file_labels = malloc(sizeof(file_label_struct));
    
    /* Mike's version of file_struct - same idea but with a file for each patch (rather than just
    patch0). */
    file_struct *file_data_store;
    file_data_store = malloc(sizeof(file_struct));

    /* This stores the prevalence, prevalence by gender and age group, incidence, number of HIV
    tests ever done, number of CD4 tests ever done, number on ART, number on ART and HIV+ etc. */
    output_struct *output;
    alloc_output_memory(&output);

    // Define the variable calibration_output_filename
    char *calibration_output_filename[NPATCHES];
    
    if(WRITE_CALIBRATION == 1){

        for(p = 0; p < NPATCHES; p++){
            calibration_output_filename[p] = (char *)calloc(LONGSTRINGLENGTH, sizeof(char));
        
            make_calibration_output_filename(calibration_output_filename[p], output_file_directory,
                python_rng_seed, patch, p, rng_seed_offset, rng_seed_offset_PC, is_counterfactual);
        
            /* Blank the calibration file. NOTE - this needs to be done outside of the i_run loop (so
            can't be done with  blank_debugging_files()) as only one of these files is generated for
            the whole set of runs. */
            blank_calibration_output_file(calibration_output_filename[p],
                allrunparameters[p][0].DHS_params->NDHSROUNDS);
        }
    }

    /* Set "partners outside community" file to just have the header (but be blank otherwise). */
    if(WRITE_PARTNERS_OUTSIDE_COMMUNITY == 1){
        initialise_partners_outside_community_file(output_file_directory, 0);
    }

    /****************************************************************************/
    /*   Loop over parameter sets from i_startrun...(i_startrun + n_startrun)   */
    /*   Remember C convention that start at 0 so index runs                    */
    /*   from (i_startrun-1)...(i_startrun-1 + n_startrun)                      */
    /****************************************************************************/

    /* SIMPLE_PARTNERSHIP_CHECK allows to either (=0) run the whole model normally or (=1) run a
    very simple partnership formation / dissolution which we used initially when designing
    partnership formation to check thigs work ok */
    if(SIMPLE_PARTNERSHIP_CHECK == 0){
        
        for(i_run = (i_startrun - 1); i_run < (i_startrun - 1 + n_startrun); i_run++){
            
            /* (re)initialise debug variables to be zero at the start of each run */
            initialise_debug_variables(debug);
            
            /* Reset PC cohort data to null. */
            for (pc_enrolment_round=0; pc_enrolment_round<NPC_ENROLMENTS;pc_enrolment_round++){
                set_to_null_pc_cohort_data(patch, 0, pc_enrolment_round, pc_size);
            }
            
            i_dhs_round = 0;  /* Reset the DHS round counter. */
            
            /* FIX - USE patch[0] FOR NOW AS THEY ARE ALL THE SAME. */
            gsl_rng_set(rng, (allrunparameters[0]+i_run)->rng_seed);
            
            /* Put the C random seed number into the saved data to make it easy to identify. */
            sprintf(tempstring, "%i,%ld,", i_run, 
                (allrunparameters[0]+i_run)->rng_seed + rng_seed_offset);
            
            if(WRITE_CALIBRATION == 1){
                for(p = 0; p < NPATCHES; p++){
                    /* This relates to the Calibration_output csv file. */
                    join_strings_with_check(output->calibration_outputs_combined_string[p], 
                        tempstring, SIZEOF_calibration_outputs, 
                        "tempstring and output->calibration_outputs_combined_string[p] in main()");
                }
            }

            /* This copies the address of the i_run th parameter set in allrunparameters into param.
            param now points to the same thing as allrunparameters[i_run] (rather than making param
            a copy of allrunparameters[i_run]).  Note that the way this is done, if we modify
            param->something it also modifies the element in allrunparameters[i_run] at the same
            time. */
            
            for(p = 0; p < NPATCHES; p++){
                patch[p].param = allrunparameters[p] + i_run; /* Use pointer arithmetic. */
                //print_param_struct(patch[p].param);      /* For debugging. */
                
                if(CHECKPARAMS == 1){
                    /* Check if parameters plausible (if not, then inputs and/or param ranges may
                    be wrong so check). */
                    check_if_parameters_plausible(patch[p].param);
                }else{
                    printf("The function check_if_parameters_plausible() is switched off for");
                    printf(" debugging. Please change CHECKPARAMS in constants.h to have value 1");
                    printf(" when finished. \n");
                }
                
                /* Set relevant counters to zero: */
                reinitialize_arrays_to_default(p, patch, overall_partnerships, output);

                /* Once first run is done, blank the individual_population key bits. */
                //// FOR DEBUGGING - MAYBE REMOVE FOR FINAL (FAST) CODE.
                if(i_run > (i_startrun - 1)){
                    blank_individual_array(patch[p].individual_population, patch[p].id_counter);
                }
                
                /* Set this to zero - determines where we are in the array fitting_data[]. */
                patch[p].i_fit = 0;
            } // END for(p = 0; p < NPATCHES; p++)
            
            if(VERBOSE_OUTPUT == 1){
                printf("Run number = %i, trial arm in patch 0 = %i\n", i_run, patch[0].trial_arm);
                fflush(stdout);
            }

            /* Make file labels of appropriate form. These are used in generating filenames. */
            make_output_label_struct(file_labels, python_rng_seed, i_run, rng_seed_offset,
                rng_seed_offset_PC, patch, is_counterfactual, 0);

            /* Generate the filenames based on these labels. */
            make_filenames_for_struct(file_labels, file_data_store, output_file_directory);

            /* make any files that use "a" (append) as the argument in fopen() blank at the start
            of the run.  We only need to do this if we are generating these files (ie if
            PRINT_EACH_RUN_OUTPUT==1). */
            if(PRINT_EACH_RUN_OUTPUT == 1){
                blank_debugging_files(file_data_store);
            }
            
            for(p = 0; p < NPATCHES; p++){
                /*********************************************************/
                /*** Initializing population characteristics ***/
                /*********************************************************/
                
                /* Initializing the population demographics (no partnerships). 
                Note that set_up_population sets id_counter to 0. */
                set_up_population(p, patch, pop);

                /* Initializing free partnerships based on how many free partnerships individuals
                have */
                init_available_partnerships(p, patch, overall_partnerships,pop);

                /* Set cumulative counters to zero: */
                init_cumulative_counters(patch[p].cumulative_outputs);
                
                
                
                /* Set calendar counters to zero: */
                init_calendar_counters(patch[p].calendar_outputs);
                
                // Reset counters of new infections and person-years each PC round.  
                for(g = 0; g < N_GENDER; g++){
                    for(a = 0; a < PC_AGE_RANGE_MAX; a++){
                        for(r = 0; r < NPC_ROUNDS; r++){
                            output->PC_ROUND_INFECTIONS[p][g][a][r] = 0;
                            output->PC_ROUND_PERSON_TIMESTEPS[p][g][a][r] = 0;
                        }
                    }
                }
            }
            
            /* Loop through multiple years of the simulation */
            for(year = patch[0].param->start_time_simul; 
            year < patch[0].param->end_time_simul; 
            year++){
                
                if(VERBOSE_OUTPUT == 1){
                    printf("Year %d\n",year);
                    fflush(stdout);
                }
                
                /* Save data on no. births/new adults/deaths if needed before it is overwritten. */
                if((WRITE_DEBUG_DEMOGRAPHICS_NBIRTHS_NEWADULTS_DEATHS == 1) &&
                    (PRINT_EACH_RUN_OUTPUT == 1)){
                    write_nbirths_nnewadults_ndeaths(file_data_store, patch, year);
                }

                /* Reset these annual counters to zero. */
                for(p = 0; p < NPATCHES; p++){
                    patch[p].DEBUG_NHIVPOSLASTYR = patch[p].DEBUG_NHIVPOS;
                    patch[p].PANGEA_N_ANNUALACUTEINFECTIONS = 0;
                    patch[p].PANGEA_N_ANNUALINFECTIONS = 0;
                    patch[p].DEBUG_NBIRTHS = 0;
                    patch[p].DEBUG_NNEWADULTS = 0;
                    patch[p].DEBUG_NDEATHS = 0;
                    
                    for(g = 0; g < N_GENDER; g++){
                        for(a = 0; a < (N_AGE_UNPD + 1); a++){
                            patch[p].py_died_from_HIV[g][a] = 0.0;
                            patch[p].n_died_from_HIV[g][a] = 0;
                        }
                    }
                }
                
                /**** FIX TO fit_flag = carry_out_processes(patch, year); ***/
                /* Carry out all weekly processes for a year. Fitting routines are also carried out
                in carry_out_processes, and if we don't fit some data then carry_out_processes()
                stops immediately and returns value 0 */
                fit_flag = -1; /* Default value - note that this should never be used except if we
                exit carry_out_processes() in some unspecified way. */
                
                /* Only write these outputs if we are not calibrating the model. */
                if(PRINT_EACH_RUN_OUTPUT == 1){
                    
                    /* Output number of people in every yearly age group in patch 0 to file: */
                    if(WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_ONEYEARINCKIDS == 1){
                        write_one_year_age_groups_including_kids(file_data_store, patch, 0, year);
                    }
                    if(WRITE_DEBUG_ART_STATE == 1){
                        write_art_states(patch, year, debug, file_data_store);
                    }
                    /* Only output this over a finite time (or every 10 years) as it is potentially
                    a large amount of data. */
                    if( (WRITE_DEBUG_HIV_STATES == 1) && 
                        (year <= DEBUG_MAX_HIV_STATE_OUTPUT_TIME || year%10 == 0)
                    ){
                        write_cd4_spvl_states(patch, year, file_data_store);
                    }
                }
                
                fit_flag = carry_out_processes(year, *fitting_data, patch, overall_partnerships,
                    output, rng_seed_offset, rng_seed_offset_PC, debug, file_data_store,
                    is_counterfactual);
                
                if(fit_flag == -1){
                    printf("Error.  Unexpected return from carry_out_processes(). Exiting.\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                
                if(fit_flag == 0 && PRINT_ALL_RUNS == 0){
                    break;
                }

                /* Age population by 1 year (update n_population, age_list, available partnerships)
                for each patch */
                for(p = 0; p < NPATCHES; p++){

                    /* Note: Update n_population first (as this uses age_list): */
                    update_n_population_ageing_by_one_year(patch, p);

                    /* Now update available partnerships (also uses age_list, so needs to be
                    updated before we update age_list): */
                    update_pop_available_partners_ageing_by_one_year(patch, p, 
                        overall_partnerships, (float) year + 1);

                    /* Now update age_list: */
                    age_population_by_one_year(patch[p].age_list);

                    /* Ageing by 1 year HIV positive people (update n_infected) */
                    age_population_size_one_year_age_by_one_year(patch[p].n_infected);

                    /* Age the 1-year age group population by 1 year: */
                        age_population_size_one_year_age_by_one_year(
                            patch[p].n_population_oneyearagegroups);

                    DEBUGNPARTNERSHIPS[p] = 0;
                    DEBUGNSERODISCORDANTPARTNERSHIPS[p] = 0;
                    DEBUGHIVNEGWITHHIVPOSPARTNERS[p] = 0;
                    NACUTE[p] = 0;
                    NCHRONIC[p] = 0;
                    
                    for(i = 0; i < patch[p].id_counter; i++){
                        if(patch[p].individual_population[i].cd4 != DEAD){
                            if(patch[p].individual_population[i].cd4 == DUMMYVALUE){
                                printf("Error. Trying to count partnerships of ");
                                printf("non-existent person. Exiting\n");
                                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                                fflush(stdout);
                                exit(1);
                            }
                            
                            DEBUGNPARTNERSHIPS[p] += patch[p].individual_population[i].n_partners;
                            if(patch[p].individual_population[i].HIV_status == UNINFECTED){
                                DEBUGNSERODISCORDANTPARTNERSHIPS[p] += 
                                    patch[p].individual_population[i].n_HIVpos_partners;
                                
                                if(patch[p].individual_population[i].n_HIVpos_partners > 0){
                                    DEBUGHIVNEGWITHHIVPOSPARTNERS[p]++;
                                }
                            }else if (patch[p].individual_population[i].HIV_status == ACUTE){
                                NACUTE[p]++;
                            }else if (patch[p].individual_population[i].HIV_status == CHRONIC){
                                NCHRONIC[p]++;
                            }
                        }
                    }
                    
                    /* store_annual_outputs is called at the end of the year, so the time at this
                    point is year+1. */
                    store_annual_outputs(patch, p, output, overall_partnerships,
                        n_infected_total + p, year, 0);
                    
                    store_annual_outputs(patch, p, output, overall_partnerships,
                        n_infected_total + p, year, 1);
                    
                    // Store the TREATS outputs
                    if(WRITE_TREATS_OUTPUT == 1){
                        store_treats_outputs(patch, p, output, overall_partnerships,
                            n_infected_total + p, year + 1, i_run + 1);
                    }
                    
                    /* Store the output associated with the cost-effectiveness analysis,
                    this is called every year.  */
                    if(WRITE_COST_EFFECTIVENESS_OUTPUT == 1){
                        store_cost_effectiveness_outputs(patch, p, output, overall_partnerships,
                            n_infected_total + p, year, i_run + 1);
                    }
                    
                    if(WRITE_ANNUAL_PARTNERSHIPS_OUTPUTS == 1){
                        store_annual_partnerships_outputs(patch, p, output, overall_partnerships,
                            n_infected_total + p, year + 1, 0);
                        store_annual_partnerships_outputs(patch, p, output, overall_partnerships,
                            n_infected_total + p, year+1, 1);
                        }

                    /* Print the sexual network every year for 5 years (used in DSMB 2016). */
                    if(
                        (year >= 1999 && year < 2005) && 
                        (WRITE_PARTNERSHIP_NETWORK_SNAPSHOT==1) && 
                        (PRINT_EACH_RUN_OUTPUT==1)
                    ){
                        print_partnership_network(file_data_store, output_file_directory,
                        file_labels, patch[p].individual_population, patch[p].id_counter, 
                        year + 1, p);
                    }
                    
                    if(
                        (p == 0) && 
                        (WRITE_PARTNERS_OUTSIDE_COMMUNITY == 1) && 
                        (PRINT_EACH_RUN_OUTPUT == 1)
                    ){
                        print_partners_outside_community(output_file_directory,
                            patch[p].individual_population, patch[p].id_counter, year + 1, p);
                    }
                    
                    /* Check if HIV has died out between the end of the introduction of HIV and the
                    start of POPART (this is CHIPS_START_YEAR[0]) - if it has, stop the simulation
                    and move on to the next parameter set. */
                    if(year>(patch[0].param->start_time_hiv+patch[0].param->n_years_HIV_seeding) && year<patch[0].param->CHIPS_START_YEAR[0]&& PRINT_ALL_RUNS==0){
                        
                        if(n_infected_total[p] == 0){
                            fit_flag = 0;
                            if(VERBOSE_OUTPUT == 1)
                                printf("No remaining HIV cases at year y=%i\n",year);
                        }
                    }

                    /* Store information needed for Calibration.csv file:
                    This part of the code writes a small set of information (prevalence, incidence)
                    into a file called Calibration...csv.  In that file is stored the prevalence
                    etc at certain timepoints for every run. We can then extract which runs meet
                    certain fitting charactersitcis and then rerun them to get the
                    Annual_output.csv and other more detailed files.
                    
                    This prevents us generating too much output (ie too many Gbs of data).
                    Note that store_annual_outputs is called at the end of the year, so the time at
                    this point is year+1. 
                    
                    We store DHS output at time point year-1 because we're taking a snapshot at the 
                    end of the year.  So if the time point for output of a DHS round is 2002 then we 
                    output the snapshot of the population at the END of 2001.  
                    */
                    if(WRITE_CALIBRATION == 1){
                        if(i_dhs_round < patch[p].param->DHS_params->NDHSROUNDS){
                            if(year == patch[p].param->DHS_params->DHS_YEAR[i_dhs_round] - 1){
                                if(VERBOSE_OUTPUT == 1){
                                    printf("Storing DHS data at time %i for DHS round %i\n", 
                                        year, i_dhs_round);
                                }
                                store_calibration_outputs_dhs(patch, p, output, year);
                            
                                /* Store everything after the final DHS: */
                                if(i_dhs_round == (patch[p].param->DHS_params->NDHSROUNDS - 1)){
                                    join_strings_with_check(
                                        output->calibration_outputs_combined_string[p],
                                        output->dhs_output_string[p], SIZEOF_calibration_outputs - 1,
                                        "output->dhs_output_string[p] and output->calibration_outputs_combined_string[p] in main()");
                                }
                                // Only increment i_dhs_round if we've looped through all patches    
                                if(p == NPATCHES - 1){
                                    i_dhs_round += 1;
                                }
                            }
                        }
                    }
                    /* Generates the Age_distribution_check files - these are used to validate the
                    model age distribution by gender against UNPD age distribution estimates. */
                    if(
                        WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER == 1 && 
                        PRINT_EACH_RUN_OUTPUT == 1
                    ){
                        write_demographics_byage_gender(patch, p, (float) year, file_data_store);
                    }
                }
                /* At least one patch has died out, and we are no longer seeding HIV. */
                if(fit_flag == 0 && PRINT_ALL_RUNS == 0){
                    printf("Stopping this run in year %i\n",year);
                    break;
                }
            } // for(year ... )
            
            /*********************************************************/
            /*** Printing stuff ***/
            /*********************************************************/
            /* Only output the full stuff if run successfully fitted data. */
            if(fit_flag == 1 || PRINT_ALL_RUNS == 1){
                /* We can switch off output if calibrating (ie set PRINT_EACH_RUN_OUTPUT to 0). */
                if(PRINT_EACH_RUN_OUTPUT == 1){
                    /* Only write this for patch p=0 as normally only have PopART intervention in
                    p=0. */
                    write_chips_data_visit(patch, 0, file_data_store, output);
                    
                    for(p = 0; p < NPATCHES; p++){
                        // Write the annual data stored in annual_outputs_string to a file
                        // called Annual_outputs.csv.
                        write_annual_outputs(file_data_store, output, p);
                        
                        if(WRITE_COST_EFFECTIVENESS_OUTPUT == 1){
                            write_cost_effectiveness_outputs(file_data_store, output, p);
                        }
                        
                        if(WRITE_ANNUAL_PARTNERSHIPS_OUTPUTS == 1){
                            write_annual_partnerships_outputs(file_data_store, output, p);
                        }
                        
                        if(WRITE_TREATS_OUTPUT == 1){
                            /* Write files for alignment with TREATS model. */
                            write_treats_outputs(file_data_store, output, p);
                        }
                        
                        if(WRITE_EVERYTIMESTEP == 1){
                            /* Write Timestep_outputs files (all ages). */
                            write_timestep_outputs(file_data_store, output, p, 0);
                            
                            /* Write Timestep_outputs PConly files (ages 18-44 only). */
                            write_timestep_outputs(file_data_store, output, p, 1);
                            
                            if(TIMESTEP_AGE == 1){
                                
                                /* Write Timestep_age_outputs files (all ages). */
                                write_timestep_age_outputs(file_data_store, output, p, 0);
                            
                                /* Write Timestep_age_outputs_PConly files (ages 18-44 only). */
                                write_timestep_age_outputs(file_data_store, output, p, 1);
                            }
                            
                            if(WRITE_ART_STATUS_BY_AGE_SEX == 1){
                                write_art_status_by_age_sex(file_data_store, output, p);
                            }
                        }
                    }
                }
            }
            
            if(fit_flag == 1 && WRITE_DEBUG_DEMOGRAPHICS_LIFE_EXPECTANCY ==1 && PRINT_EACH_RUN_OUTPUT==1){
                /* Output for patch 0 only. */
                output_life_expectancy(output_file_directory, patch, 0,  i_run+1);
            }

            /* Only output the full stuff if run successfully fitted data. */
            if(fit_flag == 1){
                if(WRITE_PHYLOGENETICS_OUTPUT >= 1 && PRINT_EACH_RUN_OUTPUT == 1){
                    /* call files to write individual data and transmission data to csv files: */
                    write_phylo_transmission_data(file_data_store, 
                        output->phylogenetics_output_string);
                    
                    write_phylo_individual_data(file_data_store, 
                        patch[PHYLO_PATCH].individual_population, patch[PHYLO_PATCH].id_counter);
                }
            }
            
            /* Write csv file of information to look at how long people live for when on ART. */
            if(WRITE_HIVSURVIVAL_OUTPUT == 1 && PRINT_EACH_RUN_OUTPUT == 1){
                write_hivpos_individual_data(file_data_store,
                    patch[PHYLO_PATCH].individual_population, patch[PHYLO_PATCH].id_counter);
            }

            if(WRITE_DEBUG_HIV_DURATION_KM == 1 && PRINT_EACH_RUN_OUTPUT == 1){
                write_hiv_duration_km_end_of_simulation(patch, 
                    patch[0].param->end_time_simul, file_data_store);
            }

            if(VERBOSE_OUTPUT == 1){
                printf("-------------------------------------\n");
                printf("Prevalent cases in patch 0:\n");
                print_population_from_one_year_data(patch[0].n_infected_wide_age_group,
                    patch[0].n_infected);
                
                printf("Incident cases in patch 0:\n");
                print_population_from_one_year_data(patch[0].n_newly_infected_wide_age_group,
                    patch[0].n_newly_infected);
                printf("-------------------------------------\n");
            }
            
            /* For prevalence it is 28 points at present (1990-2014, 2019, 2024, 2029) plus M/F
            data from PC0.
            * serostatus is 2000-2014, PC0 and 2019, 2024, 2029. */
            //int N_ANNUAL_PREVALENCE_PTS = 30;
            //int N_ANNUAL_SEROSTATUS_PTS = 20;
            
            if(WRITE_CALIBRATION == 1){
                
                // Write the CHIPS values for the calibration file to the character array
                // called `output->calibration_outputs_combined_string`.  
                for(p = 0; p < NPATCHES; p++){
                    store_calibration_outputs_chips(patch, p, output);
                    
                    // Combine PC "snapshot" outputs with PC "window" outputs
                    store_calibration_outputs_pc(patch, p, output);
                
                // Combine PC calibration values with calibration_outputs_combined_string
                join_strings_with_check(
                    output->calibration_outputs_combined_string[p],
                    output->pc_output_string[p], SIZEOF_calibration_outputs - 1,
                    "output->pc_output_string[p] and output->calibration_..._string[p] in main()");
                    
                    // Add a newline character to the output
                    strcat(output->calibration_outputs_combined_string[p], "\n");
                }
                
                /* Only want to write out to disk every NRUNSPERWRITETOFILE runs. So calculate i_run
                (mod NRUNSPERWRITETOFILE): */
                if((i_run % NRUNSPERWRITETOFILE == 0) || (i_run == n_runs - 1)){
                
                    for(p = 0; p < NPATCHES; p++){
                        write_calibration_outputs(calibration_output_filename[p],output, p);
                        
                        // Blank the string so we don't run out of memory
                        memset(output->calibration_outputs_combined_string[p], '\0',
                            SIZEOF_calibration_outputs*sizeof(char));
                        
                        memset(output->pc_output_string[p], '\0',
                            SIZEOF_calibration_outputs*sizeof(char));
                    }
                }
            }
            
            printf("-------------------------------------Simulated a total of: ");
            for(p = 0; p < NPATCHES; p++){
                printf("%ld ", patch[p].id_counter);
            }
            printf(" individuals per patch \n");
            
            /*****************************************************/
            /*** SOME MEMORY CHECKS                            ***/
            /*****************************************************/
            
            for(p = 0; p < NPATCHES; p++){
                if(patch[p].id_counter > MAX_POP_SIZE){
                    printf("Increase MAX_POP_SIZE (constants.h) to avoid memory issues\n");
                    fflush(stdout);
                }
            }

            if(VERBOSE_OUTPUT == 1){
                printf("-------------------------------------\n");
                printf("Simulated a total of: %ld partnerships\n",
                    overall_partnerships->n_partnerships[0]);
            }
            
            if(overall_partnerships->n_partnerships[0] > 
                MAX_POP_SIZE * MAX_PARTNERSHIPS_PER_INDIVIDUAL){
                printf("Allocate more memory to partnerships; no. partnerships currently greater");
                printf(" than MAX_POP_SIZE*MAX_PARTNERSHIPS_PER_INDIVIDUAL PARTNERHSIPS\n");
                fflush(stdout);
            }
            
            /* print matrix of age assortativity of partnerships at formation */
            if(CHECK_AGE_AND_RISK_ASSORTATIVITY == 1 && PRINT_EACH_RUN_OUTPUT == 1){
                print_assortativity(output_file_directory, debug, patch, i_run, file_data_store);
            }

        } // for(i_run = (i_startrun - 1); i_run < n_runs; i_run++)
    
    }else{ // else if SIMPLE_PARTNERSHIP_CHECK == 1

        /* SET ALL THE OTHER DEBUG SWITCH TO ZERO BECAUSE ALL THE OTHER DEBUG SWITCHES ARE RELATED
        TO THE OTHER MAIN CODE*/
        if(
        (SWEEP_THROUGH_TO_CHECK_LISTS != 0) ||
        (SWEEP_THROUGH_TO_CHECK_N_PARTNERS_OUTSIDE_AND_N_HIVPOS_PARTNERS_AND_N_HIVPOS_PARTNERS_OUTSIDE != 0) ||
        (CHECK_AGE_AND_RISK_ASSORTATIVITY != 0) ||
        (DEBUG_PARTNERSHIP_DURATION != 0 )||
        (WRITE_ANNUAL_PARTNERSHIPS_OUTPUTS != 0) ||
        (WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER!= 0) ||
        (WRITE_DEBUG_INITIAL_SPVL_DISTRIBUTION!= 0) ||
        (WRITE_DEBUG_CD4_AFTER_SEROCONVERSION!= 0) ||
        (WRITE_DEBUG_HIV_DURATION!= 0) ||
        (WRITE_DEBUG_HIV_DURATION_KM!= 0) ||
        (WRITE_DEBUG_HIV_STATES!= 0) ||
        (WRITE_DEBUG_ART_STATE!= 0) ||
        (WRITE_PARTNERSHIPS_AT_PC0 != 0)
        ){
            printf("Error. If SIMPLE_PARTNERSHIP_CHECK is 1 then all other swtches should be set");
            printf(" to zero, that is:\n");
            printf("- SWEAP_THROUGH_TO_CHECK_LISTS\n");
            printf("- SWEAP_THROUGH_TO_CHECK_N_PARTNERS_OUTSIDE_AND_N_HIVPOS_PARTNERS_AND_N_HIVPOS_PARTNERS_OUTSIDE\n");
            printf("- CHECK_AGE_AND_RISK_ASSORTATIVITY\n");
            printf("- DEBUG_PARTNERSHIP_DURATION\n");
            printf("- WRITE_ANNUAL_PARTNERSHIPS_OUTPUTS\n");
            printf("- WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER\n");
            printf("- WRITE_DEBUG_INITIAL_SPVL_DISTRIBUTION\n");
            printf("- WRITE_DEBUG_CD4_AFTER_SEROCONVERSION\n");
            printf("- WRITE_DEBUG_HIV_DURATION\n");
            printf("- WRITE_DEBUG_HIV_DURATION_KM\n");
            printf("- WRITE_DEBUG_HIV_STATES\n");
            printf("- WRITE_DEBUG_ART_STATE\n");
            printf("- WRITE_PARTNERSHIPS_AT_PC0\n");
            printf("Program terminating\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        /*********************************************************/
        /*** CHECKS AND DEBUGGING PARTNERSHIPS ***/
        /*********************************************************/

        patch[0].param = allrunparameters[0];
        
        /* The following functions or groups or functions are to be used ONE AT A TIME and without
        the main simulation working, otherwise individuals with same indexes will be created
        several times! */

        /***** Check 1 *****/
        /* This standalone function creates 2 women and 2 men, forms partnerships between them and
        prints output. Memory for these is allocated and freed inside the function. */
        
        reinitialize_arrays_to_default(0, patch, overall_partnerships, output);
        check_partnership_formation(overall_partnerships, allrunparameters[0], 
            debug, file_data_store);

        /***** Check 2 *****/
        /* Same as check_partnership_formation but with possible HIV transmission within
        partnerships. NOTE - only does patch p=0 at present. */
        reinitialize_arrays_to_default(0, patch, overall_partnerships, output);
        check_partnership_formation_and_HIV_acquisition(patch, 0, overall_partnerships, 
            output, debug, file_data_store);

        /* This standalone function creates 2 women and 2 men, forms partnerships, then dissolves
        some of them (at a time NOT given by the duration of the partnerships, so e.g. this is what
        would happen if one of the partners die). */
        reinitialize_arrays_to_default(0, patch, overall_partnerships, output);
        check_partnership_dissolution(overall_partnerships, allrunparameters[0], 
            debug, file_data_store);

        /* The functions below first creates an arbitrary population_size object with a certain
        distribution of the population and then calculates and prints the number of partnerships to
        be drawn between each gender/age/risk groups in one time step given this current population
        distribution (This allows checking that partnerships are drawn preferentially with similar
        age/risk groups) */
        reinitialize_arrays_to_default(0, patch, overall_partnerships, output);
        check_draw_number_partnership(patch, 0);

        /*********************************************************/
        /*** CHECKS AND DEBUGGING DEMOGRAPHICS ***/
        /*********************************************************/

        //print_demographics(individual_population, n_population, param->start_time_simul);
        //make_new_adults(cluster_hivneg_child_population,cluster_hivpos_child_population,
        //    individual_population, n_population, age_list, param->start_time_simul, param);
        //check_males_females(n_population,individual_population);
        //print_dob(n_population,individual_population);
        //print_population(n_population);
        //validate_ages_based_on_age_group(age_list, 14, param->start_time_simul+year+1);
        
        /* Print the details of this individual (12 is a random choice - can make anything from 0
        to 9999 with the current population size (10000) */
        //  print_individual(&individual_population[12]);
        
    } // else if SIMPLE_PARTNERSHIP_CHECK == 1
    
    
    /*****************************************************/
    /*** FREEING MEMORY                                ***/
    /*****************************************************/
    
    for(p = 0; p < NPATCHES; p++){
        for(i_run = 0; i_run < n_runs; i_run++){
            free(allrunparameters[p][i_run].chips_params);
            free(allrunparameters[p][i_run].PC_params);
            free(allrunparameters[p][i_run].DHS_params);
        }
    }
    
    free_pop_memory(pop,allrunparameters);
    free_output_memory(output);
    
    /* The fitting data array is allocated separately so free it separately: */
    for(p = 0; p < NPATCHES; p++){
	//free(fitting_data[p]); //this causes segmentation fault if uncommented (FDL)
        
        if(WRITE_CALIBRATION == 1){
            free(calibration_output_filename[p]);
        }
    }
    
    free(file_data_store);
    free(file_labels);
    
    /* Free patch memory last - need to free everything in patch first! .*/
    free_patch_memory(patch);
    
    free_partnership_memory(overall_partnerships);
    free(overall_partnerships);
    
    for(p = 0; p < NPATCHES; p++){
        free(allrunparameters[p]);
    }
    
    /* Only need to free memory if was allocated before - which only occurs if we don't explicitly
    specify the filename. */
    if(argc <= 4){
        free(output_file_directory);
    }
    
    free(debug);

    /***** Free GSL rng memory *****/
    gsl_rng_free(rng);

    if(VERBOSE_OUTPUT == 1){
        printf("-------------------------------------\n");
        printf("Simulation successfully finished running.\n");
        fflush(stdout);
    }
    return 0;
}
