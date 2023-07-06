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

/* Defines fixed constants and creates fixed-size global arrays for use throughout the code. */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "compat.h"

/* Standard libraries */
#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* calloc, exit, free */
#include <math.h>
#include <time.h>
#include <string.h>

/* GSL libraries */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

/************************************************************************/
/****** For switching on and off different parts of the code ************/
/************************************************************************/

#define VMMC_EFF_ZERO_AT_POPART_START 0 // Should the efficacy of VMMC be switched to zero after the start of PopART?  
#define SIMPLE_PARTNERSHIP_CHECK 0 /* if 0 then normal simulation is run, otherwise check functions are run (checks certain outputs - partnership formations etc). */
#define SWEEP_THROUGH_TO_CHECK_LISTS 0 /* if 1 then perform sweep through whole population once a year to check everyone is where they should be in list of susceptibles in serodiscordant partnerships and list of available partners
                                            NOTE THIS MAKES THE CODE VERY SLOW (because we sweep through every single every individual and their partners once a year for the whole simulation
                                            this doesn't produce any specific outputs, but stops with an error message if someone is not in a list it should belong to
                                            so if the code runs until the end with no error message it means the lists are updated as expected*/
#define SWEEP_THROUGH_TO_CHECK_N_PARTNERS_OUTSIDE_AND_N_HIVPOS_PARTNERS_AND_N_HIVPOS_PARTNERS_OUTSIDE 0 /* if 1 then perform sweep through whole population once a year to check the number of partners outside the community, the number of HIV+ partners and the number of HIV+ partners outside the community is recorded correctly for everyone
                                            this doesn't produce any specific outputs, but stops with an error message if someone is not in a list it should belong to
                                            so if the code runs until the end with no error message it means the numbers are updated as expected*/
#define CHECK_AGE_AND_RISK_ASSORTATIVITY 0 /* if 1 then performs checks, at partnership formation and once a year, of age and risk assortativity within partnerships
                                             outputs are generated in files
                                                 Age_assortativity_at_partnership_formation_runX.csv
                                                 Age_assortativity_cross_sectional_runX.csv
                                                 Risk_assortativity_at_partnership_formation_runX.csv
                                                 Risk_assortativity_cross_sectional_runX.csv
                                                 */
#define DEBUG_PARTNERSHIP_DURATION 0 /* if 1 then performs checks, at partnership formation, of partnership durations;
                                             outputs are generated in files
                                                 duration_partnership_between_high_high.csv
                                                 duration_partnership_between_low_low.csv
                                                 duration_partnership_between_med_med.csv
                                                 duration_partnership_within_high_high.csv
                                                 duration_partnership_within_low_low.csv
                                                 duration_partnership_within_med_med.csv
                                                 */

#define WRITE_ANNUAL_PARTNERSHIPS_OUTPUTS 0 // if 1 then writes files of the form Annual_partnerships_outputs_*.csv

#define WRITE_HIVSURVIVAL_OUTPUT 0 /*  Generates the files HIVsurvival_individualdata.csv - containing DoB, DoD, gender, date first on ART etc for all HIV+ in simulation. */

/************************************************************************/
/******************* Printing out more or less stuff - mostly for debugging - change eventually ******************/
/************************************************************************/

#define VERBOSE_OUTPUT 0  /* 1 if want to print normal things to screen, 0 for run on a cluster (where want to minimise output to screen). */
#define CHECKPARAMS 1 /* 1 means we run check_if_parameters_plausible() to check if parameters are in a suitable range. 0 if we are trying crazy params for debugging. */
#define PRINT_DEBUG_DEMOGRAPHICS 0      /* Prints extra demographic info to screen. */
//#define PRINT_DEBUG_DEMOGRAPHICS_NEWADULTS 0 /* Prints to screen the number of adults created in a timestep. Used to check e.g. OneYearAgeGp_CL01_Za_A_V1.2_patch0_Rand10_Run1_0.csv to check have the correct number of new adults (ie that the number of people aged 14 is equal to the number of new adults in the previous year - apart from deaths. */
#define WRITE_DEBUG_DEMOGRAPHICS_NBIRTHS_NEWADULTS_DEATHS 0           /* Generates the files BirthsNNewAdultsNdeaths.csv. For model validation. Outputs the number of births, new adults and deaths in a year. Compare with output of files produced by print_one_year_age_groups_including_kids() which allow us to see number of each of these in a year. */
#define WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER 0         /* ???? For model validation. Outputs the age distribution of the adult population using UNPD 5 year age groups. */
#define WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_ONEYEARINCKIDS 0    /* Generates the files OneYearAgeGp.csv.  For model validation. Outputs the age distribution of the population in one year age groups including kids. Note that this is not by gender. */
#define WRITE_DEBUG_DEMOGRAPHICS_LIFE_EXPECTANCY 0                    /* Generates the single file LifeExpectancy_Za.csv or LifeExpectancy_SA.csv ***ONLY IF*** the parameter end_time_simul is 2100 or later (so enough time for people to die). */
#define PRINT_DEBUG_INPUT 0             /* Prints extra stuff to screen from input.c */

#define WRITE_DEBUG_INITIAL_SPVL_DISTRIBUTION 0   /* Generates the files DEBUG_HIV_initialSPVLdistribution*.csv. */
#define WRITE_DEBUG_CD4_AFTER_SEROCONVERSION 0    /* Generates the files DEBUG_HIV_CD4_after_seroconversion*.csv. */
#define WRITE_DEBUG_HIV_DURATION 0                /* Generates the files DEBUG_HIVduration*.csv used to validate the time each person spends HIV+. */
#define WRITE_DEBUG_HIV_DURATION_KM 0             /* Generates the files DEBUG_HIVduration_KM*.csv used to do a Kaplan-Meier survival analysis of time spent HIV+. */
#define WRITE_DEBUG_HIV_STATES 0                  /* Generates the files  DEBUG_HIVstates_population*.csv used in cross-validation of model. */
#define DEBUG_MAX_HIV_STATE_OUTPUT_TIME 2016
#define WRITE_DEBUG_ART_STATE 0                   /* Generates the files DEBUG_ART_population_*.csv. These are then used by generate_art_distribution_files.py to make ART_distribution.csv and ART_transition_dist.csv files. */
#define WRITE_ART_STATUS_BY_AGE_SEX 0 // Write totals of individuals in each ART_status stratified by sex and year of age.  Write these for each time step.
#define WRITE_DEBUG_CHIPS_STATES 0                /* Generates the files CHIPS_outputs_annual*.csv containing the annual data on people when they are visited by CHiPs. */

#define WRITE_CALIBRATION 0 /* Write Calibration*.csv files to disk */
#define PRINT_ALL_RUNS 1 /* Use this if you want to print everything regardless of fitting. */
#define PRINT_EACH_RUN_OUTPUT 1 /* 0 if don't want to generate an output file for every run (when calibrating at present this is the case), or 1 if we do. */
#define WRITE_EVERYTIMESTEP 0 /* Generates the files Timestep_outputs*.csv */
#define TIMESTEP_AGE 0 /* Generates the files Timestep_age_outputs_*.csv */
#define WRITE_PHYLOGENETICS_OUTPUT 0 // 1    /* if 1 print phylo output to file, otherwise do not print */
#define WRITE_PARTNERSHIP_NETWORK_SNAPSHOT 0 /* if 1 then print out the sexual network at fixed times to allow network plots.  Writes the files Partnership_network_*.csv to disk.  The years at which partnerships are output are hard-coded in main.c */
#define WRITE_PARTNERS_OUTSIDE_COMMUNITY 0   /* if 1 then makes the file Partner_outside_inside_patch0.csv. */

#define WRITE_HAZARDS 0                      /* Generates the files Hazards_*.csv */


#define WRITE_PARTNERSHIPS_AT_PC0 0 /* Generates Distr_n_lifetime_partners and Distr_n_partners_lastyear csv files. NEEDED FOR ReadAnnualOutputs-knitr.Rnw.  */

#define FOLLOW_INDIVIDUAL -1 // 30295 // 28101 //  -1 // 1972 // 2727 // 267 // 4328  // if -1 then normal run, otherwise printing things and checking everything that happens to an individual with a certain ID

#define FOLLOW_PATCH 0 //1

#define WRITE_COST_EFFECTIVENESS_OUTPUT 0 /*  Generates a new file cost_effectiveness_$.csv */
#define WRITE_TREATS_OUTPUT 0 // Generates output for aligning the model with that used in the TREATS clinical trial.

/************************************************************************/
/************************** Random number variables *********************/
/************************************************************************/

const gsl_rng_type * TYPE_RNG;
gsl_rng * rng;

/************************************************************************/
/***************************** General variables ************************/
/************************************************************************/
#define FALSE 0
#define TRUE 1
#define DUMMYVALUE -99 /* This is a dummy value which will hopefully throw up exceptions if ever encountered in code. */

/************************************************************************/
/***************************** Time variables ***************************/
/************************************************************************/
/* time step is 1/4th of a month (approximately a week) */
#define N_TIME_STEP_PER_YEAR 48 /* This is the number of time steps in 1 year. */
#define TIME_STEP 1.0/N_TIME_STEP_PER_YEAR

//extern int MAX_N_YEARS;

#define MAX_N_YEARS 200 /* Maximum number of years the simulation will run for */

#define T_ROLLOUT_CHIPS_EVERYWHERE 2020 /* When we want post-popart CHiPs to roll out in contaminating patches. */
#define ROLL_OUT_CHIPS_INSIDE_PATCH 1
#define T_STOP_ROLLOUT_CHIPS_INSIDE_PATCH 2100 /* When to stop roll out of CHiPs to inside patch */

#define ALLOW_COUNTERFACTUAL_ROLLOUT 0 /* Should post-PopART rollout of CHiPs be allowed in counterfactual simulations?  Defaul is that it's switched off*/

/************************************************************************/
/***************************** Settings ***************************/
/************************************************************************/

/* Labels for countries (or clusters): */
#define ZAMBIA 1
#define SOUTH_AFRICA 2
/* First 12 clusters in M&E reports are always Zambia, so cluster numbers 1-12 for Zambia and >12 for South Africa: */
#define IS_ZAMBIA 12

#define NPATCHES 2

/* Trial arm: 0=ARM C, 1=ARM A, 2=ARM B. */
#define ARM_A 1
#define ARM_B 2
#define ARM_C 0

#define TIME_PC0 2014.5



/************************************************************************/
/***************************** Demographics ***************************/
/************************************************************************/

#define MAX_POP_SIZE 500000 /* the maximum population size for memory allocation */
#define MAX_N_PER_AGE_GROUP MAX_POP_SIZE/6   /* This is the maximum number of people in each adult age year group (ie 13, 14, 15...). The denominator is chosen to be really conservative (each age year will be <<5% of adult population)- as this is used for static memory allocation. It is also the maximum number of people who can die in each age group in a given timestep. */

//// BE CAREFUL, MAY NEED TO BE UPDATED IF VERY LONG PROJECTIONS ////
#define AGE_ADULT 13 /* age at which individuals enter the simulation */
#define N_AGE 7 /* number of age groups */
#define N_AGE_UNPD  14 /* Using UNPD 5 year age groups. */
#define AGE_CHIPS 18 /* Age at which CHiPs visit can occur. */
#define MAX_N_TIMESTEPS_PER_CHIPS_ROUND 96 /* Set this to be 2 years - no CHiPs round can therefore last >2 years (otherwise memory allocation issues). */
#define NCHIPSROUNDSFORFITTING 3 /* Number of CHiPs rounds we use to calibrate to. */

#define NDHSROUNDS_MAX 4 /* This is used for allocating memory to the array storing the times of the DHS rounds. In the code we use param->DHS_params->NDHSROUNDS once we know that, but that's only read in from param files at input time. */
#define AGE_DHS_MIN 15
#define AGE_DHS_MAX 59
#define DHS_AGE_RANGE_MAX 45  /* DHS runs from 15-59 so 45 age groups. */
/************************************************************************/
/***************************** PC constants *****************************/
/************************************************************************/
#define RUN_PC 0     /* If 1 then enroll and follow PC sample. */
#define NPC_ROUNDS 4 /* Number of rounds of PC visits (PC0, PC12, PC24, PC36). */
#define NPC_ENROLMENTS 1 /* Number of rounds of PC enrolment (PC0 at present only, will add PC12N+PC24N). */
#define AGE_PC_MIN 18 /* Min/Max ages at which PC enrollment can occur. */
#define AGE_PC_MAX 44
#define PC_AGE_RANGE_MAX 27  /* PC runs from 18-44 so 27 age groups. */

#define MAX_N_TIMESTEPS_PER_PC_ROUND 72 /* No PC round can last >1.5 years (in fact longest from data is 66 timesteps in community 9). */
#define N_PC_HIV_STRATA 3 /* This is the number of HIV-related categories we use for dividing up the PC sample - we want to include the right number of HIV-, HIV+ know status etc. */
#define MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP 200 /* Fixes size of list_ids_in_cohort[g][ap][i_pc_category][MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP]; */

extern const int AGE_GROUPS[N_AGE]; /* lower bounds of the age groups considered for partnership formation */
extern const int AGE_GROUPS_WITH_OLD[N_AGE+1]; /* lower bounds of the age groups considered for partnership formation */
extern const int AGE_GROUPS_UNPD[N_AGE_UNPD+1];

#define MAX_AGE 80 // Upper bound for age at start of simulation 

extern const int FIND_AGE_GROUPS_UNPD[MAX_AGE-AGE_ADULT+1];
extern const int FIND_AGE_GROUPS[MAX_AGE-AGE_ADULT+1]; /* Convert from (age-AGE_ADULT) to the AGE_GROUPS index. Each entry in the array is an AGE_GROUPS index. */

#define N_GENDER 2 /* Number of genders. */

#define MALE 0
#define FEMALE 1

#define DEAD -2 /* used as CD4 value to identify that people are dead */
#define DIEDBEFORECHIPSVISIT -2 /* Used to identify people who were scheduled to be visited by chips in current round but died beforehand. */


#define NSTEPS_GESTATION_TIME 36     /* This is timestep dependent (currently 36 timesteps = 9 months). */

#define N_UNPD_TIMEPOINTS 30 // the number of time periods for which fertility data is given by the UNPD
#define N_AGE_UNPD_FERTILITY 7 //no. of age groups in which fertility data is specified by the UNPD 
#define N_AGE_UNPD_MORTALITY 17 //no. of age groups in which mortality data is specified by the UNPD 
#define UNPD_START 1952.5 // UNPD estimates start in 1950-55. Note that because of the way they are calculated the actual point time is 1952.5.
#define UNPD_END 2097.5 // As above end period is 2095-2100 (note we are currently using medium fertility variant projections for future demographic parameters).
#define UNPD_FERTILITY_YOUNGEST_AGE 15 // UNPD fertility estimates only given for 15-49 year olds.
#define UNPD_FERTILITY_OLDEST_AGE 49

extern int POPART_SAMPLING_FRAME_ESTABLISHED;

/************************************************************************/
/***************************** Partnership ***************************/
/************************************************************************/

#define MAX_PARTNERSHIPS_PER_INDIVIDUAL 15 /* An individual can belong to up to MAX_PARTNERSHIPS_PER_INDIVIDUAL partnerships at any time point. */
#define MAX_N_PARTNERS_IN_OUTPUTS 100 /* the outputed distributions of partners will be written cumulatively from MAX_N_PARTNERS_IN_OUTPUTS onwards  */

#define MAX_BREAKUPS_PER_TIME_STEP MAX_PARTNERSHIPS_PER_INDIVIDUAL*MAX_POP_SIZE/2000 /* This is the maximum number of breakups that can happen in a given time step. Taken to be VERY conservative */

#define N_RISK 3

/* codes for sex riskiness of an individual. Note if we change these, then need to change RISK_GP_NAMES (defined in constants.c). */
#define LOW 0
#define MEDIUM 1
#define HIGH 2
extern const char RISK_GP_NAMES[N_RISK][5];

/************************************************************************/
/***************************** HIV related ***************************/
/************************************************************************/

/* Decide whether HIV test scheduling procedure happens for the whole population at fixed times,
 * or if each person has theirs scheduled sequentially. */
/* 1 = fixed times, 0 = individual. */
#define HIVTESTSCHEDULE 1

/* If 1, in hiv.c in hiv_acquisition() we boost by hand the high-high partnership transmission hazard by a factor of 2 and reduce the low-low partnership hazard by half.
 * That is based on the Cori 2013 PLOS One work. However it may lead to there appearing to be distinct HIV epidemics in low medium and high-risk.
 * We set the relative hazard by risk group to be 1 for now (5 OCtober 2016) as runs with no difference in hazard seem to have a good spread, and we therefore use Occam's razor - that if we don't need to make it more complex then we shouldn't. */
#define CHANGE_RR_BY_RISK_GROUP 0

#define DO_HIV_TESTING 1 /* Governs how background (non-CHiPs) HIV testing is carried out.
  0 = each person gets scheduled HIV tests sequentially.
  1 - we annually draw what percentage of pop gets tested and then draw a time to next test (which is uniform). */
#define RUN_POPART 1
#define YOUNGEST_AGE_SEED_WITH_HIV 18 // This gives the minimum age group index where we introduce HIV (1=18-22).
#define OLDEST_AGE_SEED_WITH_HIV 30 // This gives the maximum age group index where we introduce HIV (2=23-30).

#define ALLOW_EMERGENCY_ART 1 /* 1 - allow people to start emergency ART. 0 - stop this from happening. */

/* codes for HIV status */
#define UNINFECTED 0
#define ACUTE 1
#define CHRONIC 2


/* codes for indiv->ART_status - could merge with HIV status? Note that these are states and not processes. */
#define NARTEVENTS 8 // Currently runs from -1..6 so 8 events.
#define ARTNEG  -1 // If never tested HIV positive (note that this is tested, not serostatus).
#define ARTNAIVE 0 // Never been on ART.
#define EARLYART 1 // First few weeks/months before achieve viral suppression. Higher mortality and drop-out rate.
#define LTART_VS 2 // longer-term ART and Virally Suppressed (so low transmission, no CD4 progression or drug resistance).  
#define LTART_VU 3 // longer-term ART and Virally Unsuppressed (so higher transmission, could have (but not currently) CD4 progression and drug resistance).
#define ARTDROPOUT 4 // Has been on ART before but not currently.
#define CASCADEDROPOUT 5 // Dropped out of HIV care cascade prior to ever starting ART.
#define ARTDEATH 6 // Signals that the person needs to be removed as they die while on ART. 

/* Used as CD4 value to identify that people are not infected with HIV.
 * Note: CD4==-2 means the person is dead.
 * 0="CD4>500", 1="CD4 350-500", 2="CD4 200-350", 3="CD4 <200". */
#define CD4_UNINFECTED -1 

/* 4 set-point viral load categories (0="<4"; 1="4-4.5"; 2="4.5-5"; 3=">5") */
#define NSPVL 4
#define SPVL_INHERITANCE 0 /* 1 if inherit SPVL from infector, 0 if draw each one separately. */
/* 4 CD4 categories (0=">500", 1="350-500", 2= "200-350", 3="<200") */
#define NCD4 4

/* HIV progression events: */
#define NHIVEVENTS 4
#define HIVEVENT_ENDOFACUTE 0
#define HIVEVENT_CD4_PROGRESSION 1
#define HIVEVENT_STARTART_LOWCD4 2
#define HIVEVENT_AIDSDEATH 3
#define NOEVENT -1
#define EVENTAFTERENDSIMUL -2

/* HIV test/treatment cascade events indiv->next_cascade_event
 * Note that these are processes (in particular CD4 test), not states: */
/* Decides if a cascade event (e.g. an HIV test) is from PopART or not. 
 * If it is PopART then things happen faster (e.g. time to CD4 test is quicker), and
 * ART CD4 eligibility may be different. */
#define NOTPOPART 0
#define POPART 1

#define NCHIPSROUNDS 3 /* Number of rounds of CHiPS visits. */
#define CHIPSROUNDPOSTTRIAL -1 /* Indicates that we are post-trial. */
#define MAX_N_TIME_PERIODS_PER_ROUND 12 /* Numebr of time periods per round, used for ART initiation for which parameters are varied within a round over different time periods */
#define CHIPSNOTRUNNING -2 /* Output from get_chips_round() to show that we are before/between rounds. */

/* Note that these two numbers need to change if we add any states to the cascade: */
#define NCASCADEEVENTS_NONPOPART 7      /* We test if (EVENT<NCASCADEEVENTS_NONPOPART) to determine if the event is due to popart (=1) or not (=0). */
#define NCASCADEEVENTS 14

#define CASCADEEVENT_HIV_TEST_NONPOPART 0
#define CASCADEEVENT_CD4_TEST_NONPOPART 1
#define CASCADEEVENT_START_ART_NONPOPART 2
#define CASCADEEVENT_VS_NONPOPART 3
#define CASCADEEVENT_VU_NONPOPART 4
#define CASCADEEVENT_DROPOUT_NONPOPART 5
#define CASCADEEVENT_ARTDEATH_NONPOPART 6
#define CASCADEEVENT_HIV_TEST_POPART 7
#define CASCADEEVENT_CD4_TEST_POPART 8
#define CASCADEEVENT_START_ART_POPART 9
#define CASCADEEVENT_VS_POPART 10
#define CASCADEEVENT_VU_POPART 11
#define CASCADEEVENT_DROPOUT_POPART 12
#define CASCADEEVENT_ARTDEATH_POPART 13

#define NEVERHIVTESTED -1

/* Sensitivity and specificity of HIV tests. By default assume 1, but could potentially change these to be lower, or even time-varying. */
#define HIVTESTSENSITIVITY 1
#define HIVTESTSPECIFICITY 1


/* Determines how CD4 progression by SPVL is implemented - using 4 SPVL categories (CD4PROGRESSIONMODEL=0) or using the RR from the Cox PH model (CD4PROGRESSIONMODEL=1). */
#define CD4PROGRESSIONMODEL 1
#define USEFOURSPVLCATEGORIES 0
#define USECOXPH 1
#define COXMODELBASELINESPVL 4.0

/* Determines if this is a counterfactual run (i.e. set all arms to be arm C) or not. */
#define NOT_COUNTERFACTUAL_RUN 0
#define IS_COUNTERFACTUAL_RUN 1

/* Circumcision states for indiv->circ: */
#define UNCIRC 0
#define UNCIRC_WAITING_VMMC 1
#define VMMC 2
#define VMMC_HEALING 3
#define TRADITIONAL_MC 4


#define DEFAULT_N_HIV_PROGRESS_PER_TIME_STEP MAX_POP_SIZE/100   /* This is the default number of HIV+ people whose next progression event will happen in a given time step.
                                                               Taken to be moderately conservative - in a 50% prevalence population assume average of 1 HIV event per year. 
                                                                At some point in the future we will tune this to avoid too many reallocs(). */
#define DEFAULT_N_HIV_CASCADE_PER_TIME_STEP MAX_POP_SIZE/50   /* This is the default number of people who will have an HIV cascade event (HIV test, start ART etc) in a given time step.
                                                               Taken to be moderately conservative - in a 50% prevalence population assume average of 1 HIV event per year. 
                                                                At some point in the future we will tune this to avoid too many reallocs(). */

#define RESIZEMEM 100          /* If we run out of memory in DEFAULT_N_HIV_PROGRESS_PER_TIME_STEP, add this much extra. Note that the value is arbitrary so can change to optimise. */ 
#define RESIZEMEM_BREAKUP 100   /* If we run out of memory in planned_breakups[] array, add this much extra. Note that the value is arbitrary so can change to optimise. */

/* This is an array used by hiv_acquisition() to store the per-partnership hazard from all HIV+ partners of an individual. */
double PER_PARTNERSHIP_HAZARD_TEMPSTORE[MAX_PARTNERSHIPS_PER_INDIVIDUAL];

/* Defines amount that can be stored each year for a single patch. If we want to write very long lines we need to increase this. */
#define SIZEOF_annual_outputs_string 3000000
#define SIZEOF_annual_outputs_string_pconly 300000
/* This just is just the size for an array storing a single thing for a patch (prevalence, % know serostatus, % on ART). */
#define SIZEOF_annual_outputs_tempstore 100000
#define SIZEOF_calibration_outputs 3000000
#define SIZEOF_cost_effectiveness_outputs_string 3000000
/* To avoid excessive writing to disk, only write out these data every NRUNSPERWRITETOFILE runs. Don't set too big to avoid excessive memory usage. */
#define NRUNSPERWRITETOFILE 100
#define OUTPUTTIMESTEP 4 // Store data in timestep_outputs_string every OUTPUTTIMESTEP timestep to reduce size of file.
#define LONGSTRINGLENGTH 500 // This is the size of the arrays allocated for input/output filenames. Note that as we often specify full pathnames it should be quite big.

/* This stores the output we need for phylogenetics so that it is printed at the end of the simulation. 
 * It will contain the ID of infectee, infector, whether the infector was in acute infection etc. */
/* Estimate that each line in phylogenetics_output_string contains roughly 30 characters, and that we may have up to 
 * 50,000 transmission events per simulation. */

#define PHYLO_OUTPUT_STRING_LENGTH 400000 // This stores the output (all transmission info) for a run.
//#define PHYLO_OUTPUT_STRING_LENGTH 4000000 // This stores the output (all transmission info) for a run.
#define PHYLO_SCENARIO 0 /* used if we want to run several scenarios for PANGEA */

#define HAZARD_OUTPUT_STRING_LENGTH 20000 // This stores the hazard and associated factors.



/************************************************************************/
/* Define macros (bits of code that we always use).
 * Note - I found http://www.cprogramming.com/tutorial/cpreprocessor.html helpful to understand this! */
 /************************************************************************/

//#define AGE_INDEX(DoB,start_simul) ( (int) floor(DoB - start_simul) - (AGE_ADULT)) - youngest_age_group_index



#endif
