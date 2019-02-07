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

#ifndef STRUCTURE_H_
#define STRUCTURE_H_


/************************************************************************/
/******************************* Includes  ******************************/
/************************************************************************/

#include "constants.h"

/************************************************************************/
/*************************Define data structures*************************/
/************************************************************************/

/* Envisage desired number of partnerships for an individual to be a function of age, gender,
 * partnership risk (a variable describing some kind of preference for more or less partnerships
 * independent of age - including if in key population), possibly socio-economic, racial,  etc.
 * In addition, the partner picked will depend on geographical location.
 * We then use these to calculate the probability of each person i being picked */

/* This is a forward declaration of these structures -
 we need this because the two structures are mutually self-referencing. */
typedef struct individual individual;
typedef struct partnership partnership;


struct partnership{
    /* Structure partnership contains one pointer to a list of 2 persons, 
    and the time after which break-up occurs (in the absence of death): */
    individual* ptr[2];
    int duration_in_time_steps;
};


struct individual{
    long id; /* Unique identifier for each individual. */
    int gender; /* 0 = M, 1 = F. */
    double DoB;  /* Date of birth (stored as double). */
    double DoD;  /* Date of death (stored as double). */
    //////// think about comparative efficiency of DoB versus age ////////

    int patch_no;

    int time_to_delivery; /* -1 = not pregnant ; any positive int = number of time steps to delivery */
    //////// BE CAREFUL THIS IS TIME_STEP DEPENDENT, CHANGE IF TIME_STEP CHANGES ////////

    int HIV_status; /* 0 if uninfected, 1 if infected in acute infection, 2 if infected not in acute infection */
    int ART_status; /* -1 if never tested positive, 0 if positive but not yet on ART (or dropped out), 1 if on ART for <6 months, 2 if on ART for >=6 months and virally suppressed, 3 if on ART for >=6 months and not virally suppressed. */
    double t_sc;    /* time at which person seroconverts. Note this is used to check if the person is currently outside the window period of the given HIV test. */
    int cd4; /* Currently use -2: dead ; -1: uninfected ; 0: CD4>500, 1: CD4 350-500, 2: CD4 200-350, 3: CD4<200. */
    double SPVL_num_G; /* This is the genetic component of the log10(SPVL) - use for transmissibility and heritability. */
    double SPVL_num_E; /* This is the environmental component of the log10(SPVL) - use for transmissibility and heritability. SPVL = SPVL_G + SPVL_E. */
    double SPVL_infector; /* SPVL of infector (allows us to look at heritability) - the regression slope of viral loads V on viral loads of the infectors V’, for the last year. */
    int SPVL_cat; /* categorical variable. Derived from SPVL_num_G+SPVL_num_E in function get_spvl_cat(). 4 categories (0="<4"; 1="4-4.5"; 2="4.5-5"; 3=">5"). Use for CD4 progression.  */
    double DEBUGTOTALTIMEHIVPOS; /* For each person measure how long they are HIV+ for so can see population-level distribution of those who die by end of simulation. */
    double time_last_hiv_test;   /* Allows us to count proportion of population tested in last year, last 3 months, ever. */


    int next_HIV_event; /* -1 if not HIV+. Otherwise this stores the next HIV-biology related event to occur to this person (progression, AIDS death, starting ART because CD4<200). */
    long idx_hiv_pos_progression[2]; /* The indices which locate this individual in the hiv_pos_progression array. The first index is a function of the time to their next event (ie puts them in the group of people having an HIV event at some timestep dt) and the second is their location in this group. */
    /* for the above the origin of time is start_time_hiv */
    long debug_last_hiv_event_index; /* Stores the first index of idx_hiv_pos_progression for the last HIV progression event to happen - so can check that two HIV progression events do not happen to the same person in the same timestep. */

    int next_cascade_event; /* Stores the next testing/treatment cascade event to occur to this person (HIV test, CD4 test if in care but not on ART, start ART, stop ART, restart ART). First event is necessarily an HIV test. */ 
    long idx_cascade_event[2]; /* The indices which locate this individual in the cascade_event array. The first index is a function of the time to their next event (ie puts them in the group of people having a cascade event at some timestep dt) and the second is their location in this group. */    
    /* for the above the origin of time is COUNTRY_HIV_TEST_START */
    long debug_last_cascade_event_index; /* Stores the first index of idx_cascade_event for the last cascade event to happen - so can check that two cascade events do not happen to the same person in the same timestep. */

    int circ; /* 0 if uncircumcised, 1 if VMMC circumcised (and healed), 2 if VMMC during healing period, 3 if traditional circumcised (assumed at birth/youth). */

    long idx_vmmc_event[2];   /* The indices which locate this individual in the vmmc_event array. The first index is a function of the time to their next event (ie puts them in the group of people having a VMMC event at some timestep dt) and the second is their location in this group. */
    long debug_last_vmmc_event_index; /* Stores the first index of idx_vmmc_event for the last vmmc event to happen - so can check that two vmmc events do not happen to the same person in the same timestep. */

    // Pangea outputs for Olli                                                                                                     
    double PANGEA_t_prev_cd4stage; /* Time at which an individual last moved CD4 stage. Allows us to linearly estimate CD4 at ART initiation. */
    double PANGEA_t_next_cd4stage; /* Time at which individual will move to the next CD4 stage. Used with PANGEA_t_prev_cd4stage. */  
    double PANGEA_cd4atdiagnosis;  /* Based on t_next_cd4stage, estimate CD4 at diagnosis. */
    double PANGEA_cd4atfirstART;   /* Based on t_next_cd4stage, estimate CD4 at start of ART. */
    double PANGEA_t_diag;
    double PANGEA_date_firstARTstart;
    double PANGEA_date_startfirstVLsuppression;
    double PANGEA_date_endfirstVLsuppression;
    /* For debugging HIV cascade: */
    double DEBUG_cumulative_time_on_ART_VS;    /* Record the time someone spends on ART and VS. Zero if never VS. */
    double DEBUG_cumulative_time_on_ART_VU;    /* Record the time someone spends on ART and VU. Zero if never VU. */
    double DEBUG_cumulative_time_on_ART_early; /* Record the time someone spends in early ART. Zero if never starts ART. */
    double DEBUG_time_of_last_cascade_event;   /* Use to calculate time VS/VU/early. */

    int n_partners; /* current number of partners of this individual (who are in OR out of the community) */
    partnership* partner_pairs[MAX_PARTNERSHIPS_PER_INDIVIDUAL]; /* This is a list of the partnership pairs (with someone in the community) that this individual is in at a certain time. */

    /* The two following are only updated for HIV- individuals (who are in the community) */
    int n_HIVpos_partners; /* current number of HIV+ partners of this individual */
    int n_HIVpos_partners_outside; /* current number of HIV+ partners of this individual who are outside of the patch */
    partnership * partner_pairs_HIVpos[MAX_PARTNERSHIPS_PER_INDIVIDUAL]; /* This is a list of the partnership pairs (with someone in the community) with an HIV+ partner that this individual is in at a certain time. */

    int sex_risk; /* 0, 1, 2 for LOW, MEDIUM, HIGH risk */
    int max_n_partners; /* maximum number of partners this individual can have at a certain time point */
    int n_partners_outside; /* number of partners living outside the community at a certain time point */

    long idx_serodiscordant; // this is -1 if the individual is in no serodiscordant partnership, otherwise it is the index of this individual in the list susceptible_in_serodiscordant_partnership

    long idx_available_partner[MAX_PARTNERSHIPS_PER_INDIVIDUAL]; // this is filled in with -1 if the individual is not available for partnership (i.e. n_partners==max_n_partners), otherwise it is filled in from 0 to max_n+partners-n_partners with the index this individual is at in the list pop_available_partners (within its patch/gender/age/risk group)

    long n_lifetime_partners;             /* Counts the number of partners that someone has had in their lifetime. */
    long n_lifetime_partners_outside;             /* same but for partners outside the patch */
    long n_lifetimeminusoneyear_partners; /* Counts the number of partners that someone has had up to the start of the current year. Thus n_lifetime_partners-n_lifetimeminusoneyear_partners gives the number of partners in the most recent year. */ 
    long n_lifetimeminusoneyear_partners_outside; /* same but for partners outside the patch */

    long n_partnersminusoneyear; /* Counts the number of partners that someone has had at the start of the current year. Thus n_partnersminusoneyear + n_lifetime_partners-n_lifetimeminusoneyear_partners gives the number of partners in the most recent years (existing + new)*/

    int VISITED_BY_CHIPS_THISROUND;  /* Records whether a person has been visited by CHiPs (=1) or not (=0) this round. */
    int VISITEDBYCHIPS_TO_INIT_ART;  /* This is a flag used to generate survival curves for ART initiation. */
    int NCHIPSVISITS;                /* Records the number of visits a person has had by CHiPs throughout their life. */

    int PC_cohort_index;             /* -1 if not in a cohort. Otherwise if in PC0 this is a number in the range 0..(N_PC0-1).
                                        For PC12N we signify that it is PC12N by making it a number 100,000...(100,000+N_PC12N-1) as there are <100,000 in PC0 in each community.
                                        Similarly for PC24N it is a number 200,000...(200,000_NPC24N-1). */

    ////individual *age_list_ptr; /* This is a pointer to the place in the age list where the pointer to this person is (ie if you go to this address you will find a pointer to this structure). */
    ////////////////////////////////////////////////
    /* Do we include the following here:
       - member of household? (particularly pmtct)
       - socioeconomic_status, education/literacy (as a proxy of risk?);
       - clusternumber (includes which arm in, which country).
       - country; // 1 = Zambia, 2 = South Africa. (can ignore if have cluster number)
       - HSV;  1 = HSV2 +ve, 0 = HSV2 -ve.
       - healthcareutilization; // Includes VMMC, PMTCT, non-CHiPs HIV testing;
       - treatment
       - adherence; // Expected level of adherence to ART, loss-to-follow-up, etc.
       - double spatial_coords[2]; // Think of a simple measure (e.g. distance to travel to healthcare centre, time to travel to health-care centre) which we can get in practice.
     */
};


/* structure containing chips-related parameters: */
typedef struct{
    double prop_tested_by_chips_per_timestep[N_GENDER][MAX_AGE-AGE_CHIPS+1][MAX_N_TIMESTEPS_PER_CHIPS_ROUND][NCHIPSROUNDS];
    double prop_tested_by_chips_in_round[N_GENDER][MAX_AGE-AGE_CHIPS+1][NCHIPSROUNDS];
    int n_timesteps_per_round[NCHIPSROUNDS]; /* This is derived from the times e.g. (int) round((CHIPS_YEAR1_END-CHIPS_YEAR1_START)*N_TIME_STEP_PER_YEAR). */

    double prop_tested_by_chips_per_timestep_posttrial[N_GENDER][MAX_AGE-AGE_CHIPS+1][MAX_N_TIMESTEPS_PER_CHIPS_ROUND];
    double prop_tested_by_chips_in_round_posttrial[N_GENDER][MAX_AGE-AGE_CHIPS+1];
    int n_timesteps_per_round_posttrial; /* This is derived from the times e.g. (int) round((CHIPS_YEAR1_END-CHIPS_YEAR1_START)*N_TIME_STEP_PER_YEAR). */
} chips_param_struct;


/* structure containing PC-related parameters. Note that some of these parameters are read in by
 * read_pc_enrolment_params() and the rest by read_pc_future_params() */
typedef struct{
    int PC_START_TIMESTEP[NPC_ROUNDS]; /* This represents the timestep when the nth round of PC visits begins. */
    int PC_START_YEAR[NPC_ROUNDS];     /* This represents the  year when the nth round of PC visits begins.
                                  So the PC round starts at PC_START_TIMESTEP[n]+PC_START_TIMESTEP[n]*TIME_STEP. */
    int PC_END_TIMESTEP[NPC_ROUNDS]; /* This represents the timestep when the nth round of PC  ends. */
    int PC_END_YEAR[NPC_ROUNDS];     /* This represents the  year when the nth round of PC  ends. */

    int number_seen_by_PC_per_timestep[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA][MAX_N_TIMESTEPS_PER_PC_ROUND][NPC_ROUNDS];
    int number_enrolled_in_PC_round[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA][NPC_ENROLMENTS]; /* The number of people in each age group enrolled in PC in a single round. */
    int n_timesteps_per_round[NPC_ROUNDS]; /* This is derived from the times above. */

    double PC_retention[NPC_ROUNDS]; /* Retention in round. PC_retention[0] = 1, PC_retention[1] = (N in PC12)/(N in PC0) etc. */

    int cohort_size;      /* This is the total number of people in PC0. It equals the sum over age, gender and HIV strata of number_enrolled_in_PC_round[pc_round=0]. */
} PC_param_struct;


typedef struct{
    int NDHSROUNDS;
    int DHS_YEAR[NDHSROUNDS_MAX];
    /* Note: need to check that NDHSROUNDS_MAX is bigger than NDHSROUNDS when running code. */
} DHS_param_struct;


/* structure 'parameters': */
typedef struct {

    long rng_seed; /* This is the C random seed used for each run. */

    /********** demographics **********/
    double fertility_rate_by_age[N_AGE_UNPD_FERTILITY][N_UNPD_TIMEPOINTS]; /* UNPD outputs fertility for N_AGE_UNPD_FERTILITY=7 age groups at N_UNPD_TIMEPOINTS=30 time points (every 5 years from 1950-55 to 2095-2100). */
    /* These are based on regression of ln(UNPD estimates). UNPD gives estimates by 5 years age groups 0-4,...75-79 and then 80+. */
    double mortality_rate_by_gender_age_intercept[N_GENDER][N_AGE_UNPD_MORTALITY];
    double mortality_rate_by_gender_age_slope[N_GENDER][N_AGE_UNPD_MORTALITY];
    //double scale_fertility_param;
    double sex_ratio;     /* Proportion of new sexually active population who are male. */


    /********** times **********/
    double start_time_hiv;
    /* The two variables below are the year when HIV is seeded and the timestep in that year - so easier to check if HIV is seeded. */
    int start_time_hiv_discretised_year, start_time_hiv_discretised_timestep;
    int start_time_simul;
    int end_time_simul;
    double COUNTRY_HIV_TEST_START;
    double COUNTRY_ART_START;
    double COUNTRY_CD4_350_START;
    double COUNTRY_CD4_500_START;
    double COUNTRY_IMMEDIATE_ART_START;
    double COUNTRY_VMMC_START;

    int CHIPS_START_TIMESTEP[NCHIPSROUNDS]; /* This represents the timestep when the nth round of CHiPs visits begins. */
    int CHIPS_START_TIMESTEP_POSTTRIAL;     /* This is the timestep when each CHiPs round begins post-trial. */
    int CHIPS_START_YEAR[NCHIPSROUNDS];     /* This represents the  year when the nth round of CHiPs visits begins.
                                  So the round starts at CHIPS_START_YEAR[n]+CHIPS_START_TIMESTEP[n]*TIME_STEP. */
    int CHIPS_END_TIMESTEP[NCHIPSROUNDS]; /* This represents the timestep when the nth round of CHiPs visits ends. */
    int CHIPS_END_YEAR[NCHIPSROUNDS];     /* This represents the  year when the nth round of CHiPs visits ends.
                                  So the round finishes at CHIPS_END_YEAR[n]+CHIPS_END_TIMESTEP[n]*TIME_STEP.
                            (expect CHIPS_END_TIMESTEP[n] = CHIPS_END_TIMESTEP[n+1] -1 if the next CHIPS round starts immediately aftet the last one). */

    /********** hiv **********/
    double p_child_circ; /* probability of being traditionnally circumcised as a child */

    double eff_circ_vmmc;         /* Effectiveness of VMMC circumcision.  */
    double eff_circ_tmc;         /* Effectiveness of traditional male circumcision.  */
    double rr_circ_unhealed; /* Increased susceptibilty to HIV during VMMC healing period (ie just after operation).   */

    /* Time at which PMTCT programmes started */
    double t0_pmtct;
    double t50_pmtct;

    /* Average log viral load - this is used in hiv transmission and MTCT */
    double average_log_viral_load;

    /* Assuming average annual hazard of transmission with average viral load is 0.2. */
    double average_annual_hazard;

    /* Relative increase in transmission from acute infection, reference group is not acute untreated CD4>500 low SPVL */
    double RRacute_trans;
    double RRCD4[NCD4]; /* RRCD4[0] is 1 */
    //double RRSPVL;  /* RRSPVL compared to RRSPVL_baseline */
    //double RRSPVL_baseline;
    double SPVL_beta_k;
    double SPVL_beta_50;
    /* Relative reduction in infectivity on ART when virally suppressed and not virally suppressed (RR in infectivity =  1-effectiveness. Note in input file we use effectiveness cf HPTN052 96%) */
    double RR_ART_INITIAL;
    double RR_ART_VS;
    double RR_ART_VU;
    /* Relative increase for male to female transmission (relative to female-to-male). */
    double RRmale_to_female_trans;
    /* Range of duration (in years) of acute phase. Assume that we draw uniformly from this range for each person. */
    double min_dur_acute, max_dur_acute; /* measured in years */
    /* Currently these are simple (CD4 progression and AIDS death by VL). 
     * Will make more complex next... */
    double time_hiv_event[NCD4][NSPVL]; /* this contains the ATHENA derived CD4 duration times by SPVL category. */
    double CoxPH_SPVL_RR_CD4progression[NCD4]; /* This contains the ATHENA-derived Cox PH model estimates. */

    /* Note this is the increase in TIME (not rate) spent in each CD4 category when taking ART but virally unsuppressed. */
    double factor_for_slower_progression_ART_VU;
    /* Proportion of individuals having CD4 >500, 350-500, 200-350, <200 after end of acute phase. 
     * Note - this will depend on SPVL eventually.
     * We also calculate the cumulative totals (to make things computationally a bit quicker) ie p1, p1+p2, p1+p2+p3. Note that the sum of the 4 must be 1, so we don't store that.
     */
    double p_initial_cd4_gt500[NSPVL],p_initial_cd4_350_500[NSPVL],p_initial_cd4_200_350[NSPVL],p_initial_cd4_lt200[NSPVL];

    double cumulative_p_initial_cd4_gt500[NSPVL], cumulative_p_initial_cd4_350_500[NSPVL], cumulative_p_initial_cd4_200_350[NSPVL], cumulative_p_initial_cd4_lt200[NSPVL]; 

    //double p_initial_spvl_cat[NSPVL]; /* proportion of individuals in each SPVL category. Maybe rename p_initial to p? */
    double initial_SPVL_mu, initial_SPVL_sigma; /* proportion of individuals in each SPVL category. Maybe rename p_initial to p? */
    double SPVL_sigma_M, SPVL_sigma_E; /* Standard deviation for SPVL heritability. */

    double cumulative_p_initial_spvl_cat[NSPVL];

    double cumulative_p_misclassify_cd4[NCD4][NCD4]; /* the order ij is the same as in the ATHENA analysis */

    /********** partnerships **********/

    /* assortativity = proportion of contacts made within the same risk group ; 1-assortativity is the proportion of contacts made at random */
    double assortativity;

    /* parameter theta between 0 and 1 corresponding to the proportion of unbalanced partnership that will be balanced because of male compromise (theta=1 <-> males decide of everything) */
    double prop_compromise_from_males;

    /* Average number of new partners per year within the patch */
    double c_per_gender_within_patch[N_GENDER][N_AGE];
    double rel_rate_partnership_formation_between_patches;
    double c_per_gender_between_patches[N_GENDER][N_AGE]; // this is equal to c_per_gender_within_patch * rel_rate_partnership_formation_between_patches / (NPATCHES-1) (implicitly assuming that partnerships not within patch are split equally & randomly between other patches)

    /* relative riskiness of between patch partnership compared to within patch */
    double rr_hiv_between_vs_within_patch;

    double relative_number_partnerships_per_risk[N_RISK]; /* relative number of partnerships formed each year by individuals in a certain risk group, relative to a LOWER risk group person of the same sex and age group */
    /* relative_number_partnerships_per_risk[0] HAS TO BE 1 */

    /* Distribution of sexual partners over age  */
    double p_age_per_gender[N_GENDER][N_AGE][N_AGE]; /* constraint: sum_i p_age_per_gender[N_GENDER][N_AGE][i] = 1 */

    /* Maximum number of partners by risk group: */
    int max_n_part_noage[N_RISK];

    /* These are parameters used in time_to_partnership_dissolution() to determine when a partnership will break up.
     * At present I am assuming a Weibull distribution. */
    double breakup_scale_lambda_within_patch[N_RISK];
    double breakup_shape_k[N_RISK];

    double breakup_scale_lambda_between_patch[N_RISK];

    /********************************************************************/
    /* the following ones are not read in a file as their values change
     * as the population size and structure change ; they will be updated
     * at each time step */
    /********************************************************************/

    /* Distribution of sexual partners over risk  */
    /* They will be updated as a function of current population */
    double p_risk_per_gender[N_GENDER][N_RISK][N_RISK];

    /* Those will contain the relative number of new partners made by a female (resp male) in each risk group compared to an "average" female (resp male) */
    /* They will be updated as a function of current population */
    //double xi_per_gender[N_GENDER][N_RISK];

    double unbalanced_nb_f_to_m[N_AGE][N_RISK][N_AGE][N_RISK];
    double unbalanced_nb_m_to_f[N_AGE][N_RISK][N_AGE][N_RISK];

    long balanced_nb_f_to_m[N_AGE][N_RISK][N_AGE][N_RISK];

    /* Initial condition parameters: */
    double initial_population_size; /* Size of the 13+ population at the start of the simulation. */
    double initial_prop_age[N_AGE];    /* Proportion of sexually active population by age group at start of simulation. */
    double initial_prop_gender_risk[N_GENDER][N_RISK];  /* Proportion of individuals in each risk group (allow to be different for different genders?) */

    double initial_prop_infected_gender_risk[N_GENDER][N_RISK]; /* At time of HIV introduction, proportion of initially infected population in each gender/risk group */
    int n_years_HIV_seeding; /* Every year from start_time_hiv to start_time_hiv+n_years_HIV_seeding we seed the initial_prop_HIV_pos values above to be HIV positive (so annual introduction of HIV).
                                Set this to 0 if want a single introduction of HIV. */
    //// DO WE ALSO WANT A STRATIFICATION BY AGE?

    /* Cascade parameters: */

    /* Cascade probabilities: */

    double HIV_background_testing_rate_multiplier; /* Increase the background rate of HIV testing so we can fit the first 90 at time of CHiPs visit. */

    double HIV_rapid_test_sensitivity_CHIPS;       /* Represents the sensitivity of the rapid HIV test used by CHiPs. We might want this to be time-varying (to reflect better training/test kits used later in trial). */
    double p_collect_hiv_test_results_cd4_over200; /* given you've had an HIV test, probability that you get your results if you have CD4>200 */
    double p_collect_hiv_test_results_cd4_under200;

    //double p_collect_cd4_test_results_cd4_over200;
    //double p_collect_cd4_test_results_cd4_under200;

    /* This represents the probability that, given you've just collected the results of your HIV test and have found out you are HIV+, that you get your CD4 test results. Note that once you get your CD4 results you will (I think?) start ART with probability 1. */
    double p_collect_cd4_test_results_cd4_nonpopart;
    double p_collect_cd4_test_results_cd4_popartYEAR1;
    double p_collect_cd4_test_results_cd4_popartYEAR2onwards;


    /* Given you've just started ART several events can happen with the following probabilities: */
    double p_dies_earlyart_cd4[NCD4]; /* you die early */
    double p_leaves_earlyart_cd4_over200_if_not_die_early; /* drop out (high CD4) conditional on not dying */
    double p_leaves_earlyart_cd4_under200_if_not_die_early; /* drop out (low CD4) conditional on not dying*/
    double p_becomes_vs_after_earlyart_if_not_die_early_or_leave; /* become virally suppressed, conditional on other stuff above */
    /* or you remain virally unsuppressed */

    /* Given you've become virally suppressed, 3 possible events can happen with following probabilities */
    double p_stays_virally_suppressed; /* remain virally suppressed until you die */
    double p_stops_virally_suppressed; /* at some point you will become unsuppressed */
    /* or you will drop out */

    /* Given you've become virally unsuppressed, 2 possible events can happen with following probabilities */
    double p_vu_becomes_virally_suppressed; /* you become suppressed */
    /* or you drop out */

    /* Cascade times: */

    double time_to_background_HIVtestNOW;
    double time_to_background_HIVtest_maxval;
    double time_to_background_HIVtest_exponent;
    double time_to_background_HIVtest_midpoint;

    /* For each of these give the min and range, as we look to generate a RV of the form min+Uniform(0,1)*range. 
     * note that the max value is then t_dropout_min+t_dropout_range. 
     * The indices [2] reflect NOTPOPART (=0) and POPART(=1) */

    /* Time to drop out for someone who initiated ART but will drop out before making it to the VS/VU compartments. */ 
    double t_earlyart_dropout_min[2];
    double t_earlyart_dropout_range[2];  
    /* Time to drop out for someone who initiated ART but who dies before making it to the VS/VU compartments. */
    double t_dies_earlyart_min[2];
    double t_dies_earlyart_range[2];

    double t_end_early_art; /* Time to become virally suppressed when initiating ART. */

    double t_cd4_retest_min[2];
    double t_cd4_retest_range[2];

    double t_cd4_whenartfirstavail_min;
    double t_cd4_whenartfirstavail_range;

    /* Time between getting HIV test and getting CD4 test. */
    double t_delay_hivtest_to_cd4test_min[2];
    double t_delay_hivtest_to_cd4test_range[2];

    /* Delay in starting ART once eligible (ie once had CD4 test). */
    //double t_start_art_min[2]; // for the uniform version
    //double t_start_art_range[2]; // for the uniform version
    double t_start_art_mean_non_popart; // for the exponential version

    int n_time_periods_art_popart_per_round[NCHIPSROUNDS];

    // for the biexponential
    double t_start_art_mean_fast_popart[NCHIPSROUNDS][MAX_N_TIME_PERIODS_PER_ROUND];
    double t_start_art_mean_slow_popart[NCHIPSROUNDS][MAX_N_TIME_PERIODS_PER_ROUND];
    double p_start_art_mean_fast_popart[NCHIPSROUNDS][MAX_N_TIME_PERIODS_PER_ROUND];

    /* Time to go from viral suppression to non-suppression (if decided that will become VU). */
    double t_end_vs_becomevu_min[2];
    double t_end_vs_becomevu_range[2];

    /* Time to go from viral suppression to dropping out (if decided that will drop out). */
    double t_end_vs_dropout_min[2];
    double t_end_vs_dropout_range[2];

    /* Time to go from VU to VS (if decided that will become VS). */
    double t_end_vu_becomevs_min[2];
    double t_end_vu_becomevs_range[2];

    /* Time to go from VU to dropping out (if decided that will drop out). */
    double t_end_vu_dropout_min[2];
    double t_end_vu_dropout_range[2];

    /* Probability that someone who had dropped out of cascade is found by CHiPS team in a given year. */
    double p_popart_to_cascade[NCHIPSROUNDS];

    /* Put as cascade parameters at present - probability of getting VMMC if there is/is not popart 
     * (assuming VMMC is part of govt policy, ie after COUNTRY_VMMC_START). 
     * The [2] is again for non-popart/popart.*/
    double p_circ_nonpopart; /* probability of deciding to get circumcised (VMMC) following negative HIV test */
    double p_circ_popart[NCHIPSROUNDS]; /* Same as above but for popart. */
    double t_get_vmmc_min[2];
    double t_get_vmmc_range[2];
    double t_vmmc_healing; /* Time for VMMC wound to heal after op, during which time susceptibility may be higher. */

    //double prop_tested_by_chips[N_GENDER][NCHIPSROUNDS]; /* Proportion of population (by gender) visited by CHiPS each year. */
    chips_param_struct *chips_params;
    PC_param_struct *PC_params;
    DHS_param_struct *DHS_params;
} parameters;


// This contains the data for fitting a single timepoint. 
// We will make an array of n fitting_data_struct structures containing the n data points
// (ordered chronologically) that we fit to. */
typedef struct {
    int whatfitto; 
    /* At present these are the data I envisage we fit to:
    0 - overall prevalence, 1 - ANC prevalence (15-49), 2 = ANC prevalence (15-24),
    The following are based on PC data: 3-6 - male age gp 2-5 (18-22,23-30,31-40); 7-10 - female age gp 2-5. */

    double fit_time;
    int fit_year;
    int fit_timestep;
    double prevalence_point_est,prevalence_ll,prevalence_ul;
} fitting_data_struct;


/* Population per gender, age and risk group */
typedef struct {
    /* pop_per_gender_age_risk[g][a][r][k] points to the k^th individual 
    of gender g, age g and risk r */
    individual *pop_per_gender_age_risk[N_GENDER][N_AGE][N_RISK][MAX_N_PER_AGE_GROUP];
} population;


/* Population per gender, age and risk group */
typedef struct {
    /* pop_per_patch_gender_age_risk[p][g][a][r][k] points to the k^th individual of patch p, 
    gender g, age g and risk r */
    individual ******pop_per_patch_gender_age_risk;
} population_partners;


typedef struct {
    /* Population size per gender, age and risk group */
    long pop_size_per_gender_age_risk[N_GENDER][N_AGE][N_RISK];
} population_size;


typedef struct {
    /* Population size per gender, age and risk group */
    population_size pop_per_patch[NPATCHES];
} population_size_all_patches;


typedef struct {
    /* Population size per gender, age (one year) and risk group */

    long pop_size_per_gender_age1_risk[N_GENDER][MAX_AGE-AGE_ADULT][N_RISK];
    int youngest_age_group_index;

    /* We keep the oldest age group separate (this is anyone age MAX_AGE or over). */
    long pop_size_oldest_age_group_gender_risk[N_GENDER][N_RISK];
} population_size_one_year_age;


typedef struct {
    //These can be deduced from population_size
    long pop_size_per_gender_age[N_GENDER][N_AGE];
    long pop_size_per_gender_risk[N_GENDER][N_RISK];
    double prop_pop_per_gender_risk[N_GENDER][N_RISK];
    long total_pop_size_per_gender[N_GENDER];
    long total_pop_size;
} stratified_population_size;


/* Population size per gender, age and risk group */
typedef struct {
    double prop_size_per_gender_age_risk[N_GENDER][N_AGE][N_RISK];
    int youngest_age_group_index;
} proportion_population_size;


/* Population size per gender, age (one year) and risk group */
typedef struct {
    double prop_size_per_gender_age1_risk[N_GENDER][MAX_AGE-AGE_ADULT+1][N_RISK];
} proportion_population_size_one_year_age;


/* An array of pointers, where the ath pointer is an array of pointers to all adults in age group a (where a=0,1,2...MAX_AGE-AGE_ADULT) - note that these are not their actual ages.
   In order to make ageing easy we use a pointer which points to the youngest age group, and move this by 1 (mod n) each year (so that each group "moves" by 1). 
   One final complication is the oldest age group - we can kill them off at the end of each year or put them separate group which ages separately. */
typedef struct {
    /* This is a set of pointers to individuals stored in a matrix by age (for everybody aged AGE_ADULT to MAX_AGE), with pointers to all individuals in a given age group stored as a vector of pointers. */
    individual *age_group[MAX_AGE-AGE_ADULT][MAX_N_PER_AGE_GROUP];

    /* Index for the position of the vector of the youngest age group in age_group[][]. */ 
    ////1:individual **youngest_age_group; /* because of the way we coded ageing, the youngest age group is not always the 1st, it is the youngest_age_group^th */
    int youngest_age_group_index;

    /* We keep the oldest age group separate (this is anyone age MAX_AGE or over) and we kill them at a fast rate. */
    /* This is an array of pointers. */
    individual *oldest_age_group[MAX_N_PER_AGE_GROUP];
    /* This could equally well be a pointer to the last person in the oldest_age_group[] array. */
    long number_oldest_age_group;
    /////// DEBUGGING: Just reduced this from MAX_AGE-AGE_ADULT+1
    long number_per_age_group[MAX_AGE-AGE_ADULT]; //// SHOULD THIS BE A LONG? NOT SURE
    /* Points to the youngest element in number_per_age_group[]. */ 
    ////1:int *youngest_number_per_age_group_ptr;
    /* Index for the youngest element in number_per_age_group[]. This is used by update_age_list_death(). */
    ////1:int youngest_number_per_age_group_index;
}age_list_element;


typedef struct {
    age_list_element *age_list_by_gender[N_GENDER];
} age_list_struct;


/* This structure contains a list (n_child) of the number of children to be born at each timestep.
n_child is an array containing the number of children for each timestep, and
transition_to_adult_index_n_child is a pointer to the place in the list where the current adults 
are stored (rather than having to move everybody by one slot each timestep, it is easier to have a
pointer to the "end" of the list - ie kids about to turn AGE_ADULT.  At each time step, if n_t
people are born, then n_t is added to the "start" of the array and the pointer is incremented by 1
(so that it points to the "old" kids about to turn 13 at that timestep). */
typedef struct {
    /* Number of children (ie <=AGE_ADULT - as children turn to adults at
    AGE_ADULT+1-TIME_STEP=13.99 years) - this is an array where each element is the number of
    children who were born in a given timestep (so that in each timestep of the simulation we add
    more new people) */
    long n_child[(AGE_ADULT+1)*N_TIME_STEP_PER_YEAR];
    long *transition_to_adult_index_n_child;
    // A version of transition_to_adult_index_n_child, storing the index rather than a pointer
    int debug_tai;
} child_population_struct;


typedef struct{
    /* We divide the population into subgroups for the CHiPS teams to sample from each year. 
     * The same structure is then used to store the people who will get visited by CHiPs teams that
    year. */ 
    
    /* contains ID of people who may be visited by CHiPs in a round. The 100 is fairly arbitrary. */
    long list_ids_to_visit[N_GENDER][MAX_AGE-AGE_CHIPS+1][MAX_POP_SIZE/100];
    
    /* The number of people in each age group visited by CHiPs in a single round. */
    long number_to_visit[N_GENDER][MAX_AGE-AGE_CHIPS+1];
    
    /* This will be the number of people (m/f) we see each timestep in a given year
    * I use this form so it is easy to use the CHiPs reports to allow us to put in variation in
    * the number of people seen each week. 
    * Note that these arrays are only used in the variable which is the list of people to be seen.
    * It is not used in the sampling frame variable. */
    
    long number_to_see_per_timestep[N_GENDER][MAX_AGE-AGE_CHIPS+1][MAX_N_TIMESTEPS_PER_CHIPS_ROUND];
    long next_person_to_see[N_GENDER][MAX_AGE-AGE_CHIPS+1];
} chips_sample_struct;


/* This is the sampling frame we use for drawing the PC sample.
 * We divide the population into subgroups by gender, age and HIV status (HIV-, HIV+ unaware, HIV+ aware).
 * The same structure is then used to store the people who will get visited by PC teams that year.
 * Note that this structure gets recycled for each new enrolment (PC0, PC12N, PC24N).
 * Once an individual is *successfully* enrolled, they are transferred to PC_cohort_struct, which is the
 * structure that actually contains the PC data. */
typedef struct{
    /* contains the ID of people who *may* be enrolled by PC in a single round. In PC0 Makululu has 228 people aged 18 - this is the highest number from a single community in a single year age group. */
    long list_ids_potential_enrollees[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA][MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP];


    /* The number of people in each age group who could potentially be enrolled in PC. */
    long number_in_sample_including_reserves[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA];

    long next_person_to_see[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA];
    }PC_sample_struct;


    /* This stores people who have *successfully* been enrolled in a given PC round (PC0, PC12N, PC24N). */
    typedef struct{
    long list_ids_in_cohort[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA][MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP];
    /* This tells us how many in list_ids_in_cohort. */
    int number_in_cohort[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA];
    /* Index of where we are in list_ids_in_cohort at current timestep. This is analogous to next_person_to_see[] in PC_sample_struct, but in that one the next_person_to_visit included people who may not have been eligible in the cohort. */
    int next_person_to_visit[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA];

    /* This is people who leave the cohort after enrolment (migrated/not found) apart from death. We keep them so we can resample
     * if we want to allow people to come back into the cohort later on. */
    //long list_ids_in_dropped_out_from_cohort[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA][100];
    /* This tells us how many in list_ids_in_dropped_out_from_cohort. */
    //int number_in_dropped_out_from_cohort[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA];
} PC_cohort_struct;


/* Data that is saved for each person in PC cohort. */
typedef struct{
    long id; /* Unique identifier (ie indiv->id) allows access to other characteristics.  */
    double PC_visit_dates[NPC_ROUNDS]; /* Date of visits by PC. If the person was enrolled later (or dropped out) use -1 to signify that they were not visited in this round. */
    int HIV_status[NPC_ROUNDS];
    int ART_status[NPC_ROUNDS];
    int RETAINED_IN_COHORT[NPC_ROUNDS];            /* =1 if still in cohort, 0 if dropped out/died. */
    int gender;
    int ap;
    int serodiscordant_status[NPC_ROUNDS];                     /*Tracks whether person is HIV- and has an HIV+ partnership or not. */
} PC_cohort_individual_data_struct;


/* Structure will contain the PC data (actually pointers to the data).
Data is updated in PC_enroll_person() and PC_visit_person() in pc.c. */
typedef struct{
    PC_cohort_individual_data_struct *PC0_cohort_data;
    PC_cohort_individual_data_struct *PC12N_cohort_data;
    PC_cohort_individual_data_struct *PC24N_cohort_data;
    int PC_cohort_counter[NPC_ENROLMENTS];
} PC_cohort_data_struct;


typedef struct{
    long N_total_CD4_tests_nonpopart; /* Cum. num CD4 tests done (excluding those from PopART). */
    long N_total_HIV_tests_nonpopart; /* Cum. num HIV tests done (excluding those from PopART). */
    long N_total_CD4_tests_popart; /* Cum. num of CD4 tests done by PopART. */
    long N_total_HIV_tests_popart; /* Cum. num of HIV tests done by PopART. */
    
    /* Number of people who ever started ART (including those dead) for non-popart & popart: */
    long N_total_ever_started_ART_nonpopart;
    long N_total_ever_started_ART_popart;
} cumulative_outputs_struct;


/* 'patch' refers to a geographical unit - this may be trial community,
inside/outside popart, district, or any other geographical unit. */
typedef struct{
    age_list_struct *age_list;
    individual *individual_population;
    population_size *n_population;
    population_size_one_year_age *n_population_oneyearagegroups;
    stratified_population_size *n_population_stratified;
    child_population_struct *child_population;
    parameters *param;
    long *new_deaths;
    long *death_dummylist;

    population_size_one_year_age *n_infected;
    population_size_one_year_age *n_infected_cumulative;
    population_size_one_year_age *n_newly_infected;
    population_size *n_infected_wide_age_group;
    population_size *n_newly_infected_wide_age_group;

    long n_newly_infected_total;
    long n_newly_infected_total_from_outside;
    long n_newly_infected_total_from_acute;

    long n_newly_infected_total_pconly;
    long n_newly_infected_total_from_outside_pconly;
    long n_newly_infected_total_from_acute_pconly;

    int country_setting;
    individual ***hiv_pos_progression;
    long *n_hiv_pos_progression;
    long *size_hiv_pos_progression;
    individual ***cascade_events;
    long *n_cascade_events;
    long *size_cascade_events;
    individual ***vmmc_events;
    long *n_vmmc_events;
    long *size_vmmc_events;
    chips_sample_struct *chips_sample;
    PC_sample_struct *PC_sample; /* PC sample is the group of people (including reserves) before enrolment. Some people may not make it into the cohort (as they have moved out of the relevant subgroup). */
    PC_cohort_struct *PC_cohort; /* This is the actual PC cohort - everyone in this group was visited in PC0. */
    PC_cohort_data_struct *PC_cohort_data; /* This stores the PC data for each person (HIV status at each round, date of visit at each round etc). */
    cumulative_outputs_struct *cumulative_outputs;
    int n_fit;
    int i_fit;
    fitting_data_struct *fitting_data;
    /* Store number of births (of children who will survive to adulthood), deaths and new adults in a year for model validation. */
    long DEBUG_NBIRTHS;
    long DEBUG_NNEWADULTS;
    long DEBUG_NDEATHS;

    long DEBUG_NHIVPOS;
    long DEBUG_NHIVPOSLASTYR;
    long OUTPUT_NDIEDFROMHIV;
    long PANGEA_N_ANNUALACUTEINFECTIONS;
    long PANGEA_N_ANNUALINFECTIONS;
    /********* community/country/arm information *********/
    int community_id; /* Community number (1-21 in PopART). */

    /* Arm of the trial being modelled. Determines the HIV testing algorithm. */
    int trial_arm;

    long id_counter;  /* Unique ID of each person in the simulation.
     * It is also their position in the individual_population[] array. */

    long ****cross_sectional_distr_n_lifetime_partners;
    long ****cross_sectional_distr_n_partners_lastyear;
} patch_struct;


// Structure containing information on partnerships
typedef struct{
    partnership* partner_pairs;
    long *n_partnerships;
    population_partners* pop_available_partners;
    population_size_all_patches *n_pop_available_partners;
    long *new_partners_f_sorted;
    long *shuffled_idx;
    long *new_partners_f_non_matchable;
    long *new_partners_m;
    long *new_partners_m_sorted;
    long *partner_dummylist;

    partnership*** planned_breakups;
    long* n_planned_breakups;
    long* size_planned_breakups;

    individual** susceptible_in_serodiscordant_partnership;
    long *n_susceptible_in_serodiscordant_partnership;
} all_partnerships;


// Structure containing ART outputs used in debugging. 
// Housed in a single structure to make the debug_struct more readable. 
typedef struct{

    long n_start_emergency_art_fromuntested;
    long n_start_emergency_art_fromartnaive;
    long n_start_emergency_art_fromartdroupout;
    long n_start_emergency_art_fromcascadedropout;
    long n_start_emergency_art; /* Total to make sure we haven't missed any possible options. */

    long cascade_transitions[NARTEVENTS][NARTEVENTS];
} debug_art_vars;


// Structure containing miscellaneous stuff used for debugging.  
// These are not used unless debug mode is on 
typedef struct{

    long age_of_partners_at_partnership_formation[MAX_N_YEARS][N_AGE][N_AGE]; /* age_of_partners_at_partnership_formation[y][age_f][age_m] is the number of partnerships formed in year y (counted from start_time_simul) between females aged age_f and males aged age_m */
    long risk_of_partners_at_partnership_formation[MAX_N_YEARS][N_RISK][N_RISK]; /* risk_of_partners_at_partnership_formation[y][risk_f][risk_m] is the number of partnerships formed in year y (counted from start_time_simul) between females of risk group risk_f and males of risk group risk_m */
    double age_of_partners_cross_sectional[MAX_N_YEARS][N_AGE][N_AGE]; /* age_of_partners_cross_sectional[y][age_f][age_m] is the number of partnerships which exist at the end of year y (counted from start_time_simul) between females aged age_f and males aged age_m */
    double risk_of_partners_cross_sectional[MAX_N_YEARS][N_RISK][N_RISK]; /* risk_of_partners_cross_sectional[y][risk_f][risk_m] is the number of partnerships which exist at the end of year y  (counted from start_time_simul) between females of risk group risk_f and males of risk group risk_m */
    debug_art_vars art_vars[NPATCHES];

} debug_struct;


typedef struct{ /* structure which contains all the strings that are outputted */

    char *annual_outputs_string[NPATCHES];
    char *annual_outputs_string_pconly[NPATCHES];

    char *annual_partnerships_outputs_string[NPATCHES];
    char *annual_partnerships_outputs_string_pconly[NPATCHES];

    char *timestep_outputs_string[NPATCHES];
    char *timestep_outputs_string_PConly[NPATCHES];

    char *timestep_age_outputs_string[NPATCHES];
    char *timestep_age_outputs_string_PConly[NPATCHES];

    char *chips_output_string[NPATCHES];
    char *dhs_output_string[NPATCHES];
    char *calibration_outputs_combined_string[NPATCHES];

    char *phylogenetics_output_string;

    char *hazard_output_string; /* Stores hazard and other factors (e.g. age, risk gp, partner in/outside community) to allow us to examine whether the average hazard is credible. */

    long NCHIPS_VISITED[NPATCHES][N_GENDER][MAX_AGE-AGE_CHIPS+1][NCHIPSROUNDS];
    long NCHIPS_HIVPOS[NPATCHES][N_GENDER][MAX_AGE-AGE_CHIPS+1][NCHIPSROUNDS];
    long NCHIPS_HIVAWARE[NPATCHES][N_GENDER][MAX_AGE-AGE_CHIPS+1][NCHIPSROUNDS];
    long NCHIPS_ONART[NPATCHES][N_GENDER][MAX_AGE-AGE_CHIPS+1][NCHIPSROUNDS];
    long NCHIPS_VS[NPATCHES][N_GENDER][MAX_AGE-AGE_CHIPS+1][NCHIPSROUNDS];

} output_struct;


typedef struct{

    /* Label for files where data from each patch is stored in a separate file: */
    char filename_label_bypatch[NPATCHES][LONGSTRINGLENGTH];
    /* Label for files where  data from all patches is stored in a single file.
     * Used in demographics debugging NBirthsNNewAdultsNdeaths*.csv files */
    char filename_label_allpatches[LONGSTRINGLENGTH];

    /* Label for files where  data from all patches is stored in a single file.
     * Used in partnership debugging files */
    char filename_label_allpatches_witharm_communityno[LONGSTRINGLENGTH];

} file_label_struct;


typedef struct{

    /* Main outputs: */
    FILE *ANNUAL_OUTPUT_FILE[NPATCHES];
    char filename_annual_output[NPATCHES][LONGSTRINGLENGTH];

    FILE *TIMESTEP_OUTPUT_FILE[NPATCHES];
    /* These will use the same file pointer (as we close at the end of the relevant function we can recycle). */
    char filename_timestep_output[NPATCHES][LONGSTRINGLENGTH];
    char filename_timestep_output_PConly[NPATCHES][LONGSTRINGLENGTH];

    FILE *TIMESTEP_AGE_OUTPUT_FILE[NPATCHES];
    char filename_timestep_age_output[NPATCHES][LONGSTRINGLENGTH];
    char filename_timestep_age_output_PConly[NPATCHES][LONGSTRINGLENGTH];

    FILE *CALIBRATION_FILE[NPATCHES];
    char filename_calibration_output[NPATCHES][LONGSTRINGLENGTH];

    FILE *ANNUAL_PARTNERSHIP_OUTPUT_FILE[NPATCHES];
    char filename_annual_partnership_output[NPATCHES][LONGSTRINGLENGTH];

    /* Stores PC output: */
    FILE *PC_DATAFILE[NPATCHES];
    char filename_PC_data[NPATCHES][LONGSTRINGLENGTH];

    /* Stores ChiPs output at the end of each year: */
    FILE *ChipsAnnual_DATAFILE[NPATCHES];
    char filename_chipsannual_data[NPATCHES][LONGSTRINGLENGTH];

    /* Stores ChiPs output in the round when people are visited (ie as people are visited we first store their outputs then carry out visit (HIV test etc): */
    FILE *ChipsVisit_DATAFILE[NPATCHES];
    char filename_chipsvisit_data[NPATCHES][LONGSTRINGLENGTH];

    /* Partnership-related files: */
    FILE *DUR_BETWEEN_HIGHHIGH;
    FILE *DUR_BETWEEN_MEDMED;
    FILE *DUR_BETWEEN_LOWLOW;
    FILE *DUR_WITHIN_HIGHHIGH;
    FILE *DUR_WITHIN_MEDMED;
    FILE *DUR_WITHIN_LOWLOW;

    char filename_DUR_BETWEEN_HIGHHIGH[LONGSTRINGLENGTH];
    char filename_DUR_BETWEEN_MEDMED[LONGSTRINGLENGTH];
    char filename_DUR_BETWEEN_LOWLOW[LONGSTRINGLENGTH];
    char filename_DUR_WITHIN_HIGHHIGH[LONGSTRINGLENGTH];
    char filename_DUR_WITHIN_MEDMED[LONGSTRINGLENGTH];
    char filename_DUR_WITHIN_LOWLOW[LONGSTRINGLENGTH];

    FILE * distr_n_lifetime_partners[NPATCHES];
    FILE * distr_n_partners_lastyear[NPATCHES];

    char filename_distr_n_lifetime_partners[NPATCHES][LONGSTRINGLENGTH];
    char filename_distr_n_partners_lastyear[NPATCHES][LONGSTRINGLENGTH];

    FILE *age_assortativity;
    char filename_age_assortativity[LONGSTRINGLENGTH];

    FILE *age_assortativity_cross_sectional;
    char filename_age_assortativity_cross_sectional[LONGSTRINGLENGTH];

    FILE *risk_assortativity;
    char filename_risk_assortativity[LONGSTRINGLENGTH];

    FILE *risk_assortativity_cross_sectional;
    char filename_risk_assortativity_cross_sectional[LONGSTRINGLENGTH];

    /* Demographics-related files: */
    /* File outputs the age distribution of the adult population in 5 year age groups at fixed times (annually).
     * Generates files e.g. Age_distribution_check_CL01_Za_A_V1.2_patch0_Rand10_Run1_0.csv
     *   */
    FILE *AGEDISTRIBUTIONFILE[NPATCHES];
    char filename_debug_agedistribution[NPATCHES][LONGSTRINGLENGTH];

    FILE *NBIRTHS_NNEWADULTS_NDEATHS_FILE;
    char filename_debug_nnewadults_ndeaths_file[LONGSTRINGLENGTH];

    FILE *ONEYEARAGEDISTRIBUTIONFILE[NPATCHES];
    char filename_debug_one_yearage_dist_includingkids[NPATCHES][LONGSTRINGLENGTH];


    /* HIV-related files: */
    /* File outputs the total duration of any HIV-positive person who dies from AIDS (so not people who are HIV+ but die from
     * natural causes). Allows us to check that distribution of duration of HIV is realistic (nobody dying quickly/slowly).  */
    FILE *HIVDURATIONFILE[NPATCHES];
    char filename_debughivduration[NPATCHES][LONGSTRINGLENGTH];

    FILE *HIVDURATIONFILE_KM[NPATCHES];
    char filename_debughivduration_km[NPATCHES][LONGSTRINGLENGTH];

    FILE *HIVCD4_AFTER_SEROCONVERSION[NPATCHES];
    char filename_debughivcd4_after_seroconversion[NPATCHES][LONGSTRINGLENGTH];

    FILE *HIV_INITIAL_SPVL_DISTRIBUTION[NPATCHES];
    char filename_debuginitial_spvl_distribution[NPATCHES][LONGSTRINGLENGTH];

    /* File contains the gender, CD4 stage, SPVL and ART status at given times (yearly) of each individual HIV+ person. Allows us to check HIV biology. */
    FILE *HIVSTATEPOPULATIONFILE[NPATCHES];
    char filename_debug_hivpopulation[NPATCHES][LONGSTRINGLENGTH];

    /* File contains the number of HIV+ people in each ART stage at given times (yearly). Allows us to check we're not leaking people to e.g. drop-out. */
    FILE *ARTPOPULATIONFILE[NPATCHES];
    char filename_debug_artpopulation[NPATCHES][LONGSTRINGLENGTH];

    /* Hazard outputs - ie the hazard for every person in a s/d partnership over a given period: */
    FILE *HAZARD_FILE;
    char filename_hazard_output[LONGSTRINGLENGTH];


    /* Phylogenetic outputs: */
    FILE *PHYLOGENETIC_TRANSMISSION_FILE;
    char filename_phylogenetic_transmission[LONGSTRINGLENGTH];

    FILE *PHYLOGENETIC_INDIVIDUALDATA_FILE;
    char filename_phylogenetic_individualdata[LONGSTRINGLENGTH];

    FILE *HIVSURVIVAL_INDIVIDUALDATA_FILE;
    char filename_hivsurvival_individualdata[LONGSTRINGLENGTH];

} file_struct;

#endif
