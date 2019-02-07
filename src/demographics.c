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

/* Demographic processes for the PopART model:
per_woman_fertility_rate()
    Rate at which a given woman gets pregnant (or rate at which children born per capita).  
natural_death_rate()
    Natural death rate (as a function of age and calendar time at present). 
draw_sex_risk()
    Decides an individual's sexual risk group when they become adults (ie enter population at
    AGE_ADULT).
create_new_individual()
    Sets up a new individual structure when someone is born, initializing all parameters.
initialize_first_cascade_event_for_new_individual()
    Schedules first cascade event for a new adult if HIV testting has started (otherwise no need as
    everyone gets a cascase event assigned when HIV testing begins.) 
update_population_size_new_adult()
    Updates the population size variable when new adult enters adult population.
update_population_size_death()
    Updates population size variable when someone dies.
update_age_list_new_adult()
    Update age list structure when new adult enters adult population.
update_age_list_death()
    Update age list structure when someone dies.
get_age_index()
    Given DoB gives the index in the array age_list
get_age_group()
    Given time and DoB gives the age group index a person belongs to (age groups 13-17,18-22, etc).

 *** The following functions are related to the annual ageing of the whole population by 1 year 
(ie many lists divide the population into yearly age groups. At the end of each year these lists
are "moved" by 1 to represent ageing).

update_n_population_ageing_by_one_year()
    Update the n_population structure as people age by 1 year.
update_pop_available_partners_ageing_by_one_year()
    Updates list of available partners as people age by 1 year. 
age_population_by_one_year()
    Update the age_list structure as people age by 1 year.
age_population_size_one_year_age_by_one_year()
    Update the population

 *** The following functions remove people from the lists of scheduled events/available people.
remove_dead_person_from_susceptible_in_serodiscordant_partnership()
remove_dead_person_from_list_available_partners()
remove_dead_persons_partners()
remove_from_hiv_pos_progression()
    note this removes a person from HIV+ progression when on effective ART as well as death.
remove_from_cascade_events()
    Removes someone from cascade_events either due to death, or  we need to schedule an earlier
    event because of popart
remove_from_vmmc_events()
    removes man from scheduled VMMC event when dying.
deaths_natural_causes()
    For each age group calculate the probability of dying p_x per timestep, draws people to die 
    assuming a binomial process with prob p_x, then for each of them calls all the remove_()
    functions above to delete that individual from any lists of events/available people.
make_new_adults()
    Looks up how many kids reach adulthood at each timestep, and then makes them.
add_new_kids()
    This is a DUMMY function to add newly born babies to the child_population structures so that
    there are always new individuals to reach adulthood as time goes on
make_pop_from_age_list()
    Update the struct pop (counts of people by age/gender etc) based on the (always up-to-date)
    struct age_list.
individual_death_AIDS()
    For an indiv who dies of age this function calls all the processes to remove them from all
    appropriate lists.
 */

#include "structures.h"
#include "constants.h"
#include "utilities.h"
#include "demographics.h"
#include "init.h"
#include "hiv.h"
#include "debug.h"


double per_woman_fertility_rate(int age, parameters *param, int y0, double f){
    /* Calculate per-woman fertility rate based on age using UNPD rates
    
    Return the rate at which any one woman of age `age` gets pregnant and has an offspring that will
    survive until AGE_ADULT.  Note, in the simulation it is also checked that the women who get
    pregnant have at least one partner at that time.  
    
    Interpolate fertility rate over time and age.  UNPD fertility data is in 5-year age groups
    (e.g. 15-19, 20-24, 25-29) which are converted to an in index (e.g. 0, 1, 2).  
    
    Potential changes: fertility may depend on HIV status and time.  total_fertility_rate should be
    able to vary over time and country.
    
    
    Arguments
    ---------
    age : int
        Individual's age in years
    
    param : pointer to parameters strucuture
    
    y0 : int
        Time index for the array fertility_rate_by_age[][]
    
    f : double
        Interpolation coefficient over time
    
    Returns
    -------
    Per-year probability that a woman this age gets pregnant.
    
    */
    
    double result;
    if(age > UNPD_FERTILITY_OLDEST_AGE || age < UNPD_FERTILITY_YOUNGEST_AGE){
        
        printf("ERROR: in per_woman_fertility_rate() age %i lies outside fertile ages %i-%i\n",
            age, UNPD_FERTILITY_YOUNGEST_AGE, UNPD_FERTILITY_OLDEST_AGE);
        
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
        
    }
    
    // Calculate the index for the age in question (compared to the UNPD age categories)
    int a_unpd = (age - UNPD_FERTILITY_YOUNGEST_AGE) / 5;
    
    if(y0 == (N_UNPD_TIMEPOINTS - 1)){
        // No interpolation at the far end as 
        // param->fertility_rate_by_age[a_unpd][N_UNPD_TIMEPOINTS] not defined.
        result = param->fertility_rate_by_age[a_unpd][y0];
    }else{
        result = (1 - f) * param->fertility_rate_by_age[a_unpd][y0] +
            f * param->fertility_rate_by_age[a_unpd][y0 + 1];
    }
    return result;
}


void get_unpd_time_indices(double t, int *y0, double *f){
   /* Calculate index for arrays of UNPD fertility parameters and fraction of time through period
    
    This function is used when interpolating UNPD fertility rates over time.  Given the current
    time t, this function calculates the corresponding array index y0 (and y0+1) which are the
    array indices for fertility_rate_by_age that we interpolate over, and the fraction f so the
    interpolation over time is f*fertility_rate_by_age[][y0] + (1-f)*fertility_rate_by_age[][y0+1].
    
    Note that we do not interpolate over age here. 
    
    Arguments
    ---------
    
    t : double
        Current time
    y0 : pointer to an int
        
    f : pointer to a double
        
    
    Returns
    -------
    Nothing; the variables y0 and f are populated.  
    
    */
    if(t <= UNPD_START){
        *y0 = 0;
        *f = 0.0;
    }else if(t >= UNPD_END){
        *y0 = N_UNPD_TIMEPOINTS - 1;
        *f = 1.0;
    }else{
        *y0 = (int) floor((t - UNPD_START) / 5.0);
        *f = (t - UNPD_START) / 5.0 - *y0;
    }
    return;
}


/* Calculate the total probability of dying before reaching adulthood AGE_ADULT.
 * In this calculation we assume that the probability of dying aged 10-AGE_ADULT is negligible (which it is in Zambia/South Africa).
 * We therefore use P(dying in childhood) = 1 - P(don't die age 0-4)*P(don't die age 5-9).
 * Note that the mortality rates for e.g. dying age 0-4 are rates per year, so e.g. P(don't die age 0-4) = pow(1-mortality_rate_under5,5).*/
double childhood_mortality(parameters *param, double t){
    double mortality_rate_under5 = 0.0;
    double mortality_rate_5to10 = 0.0;
    double mortality_rate_childhood_total;
    int g;


    for (g=0;g<N_GENDER;g++){

        /* We by default assume that mortality in under 5 is mostly perinatal mortality (so occurs at time t) and that mortality in 5-10 year olds occurs uniformly over that age so on average at time t+7.5.
         * However we need to adjust these times for the fact that we only have data from 1950-2100. */
        if (t<1950){
            /* Average mortality over genders. The [0] and [1] indices refer to age groups 0-4 and 5-9: */
            mortality_rate_under5 += exp(param->mortality_rate_by_gender_age_intercept[g][0] + param->mortality_rate_by_gender_age_slope[g][0]*1950);
            mortality_rate_5to10 += exp(param->mortality_rate_by_gender_age_intercept[g][1] + param->mortality_rate_by_gender_age_slope[g][1]*1957.5);
        }
        else if (t>=2100){
            mortality_rate_under5 += exp(param->mortality_rate_by_gender_age_intercept[g][0] + param->mortality_rate_by_gender_age_slope[g][0]*2100);
            mortality_rate_5to10 += exp(param->mortality_rate_by_gender_age_intercept[g][1] + param->mortality_rate_by_gender_age_slope[g][1]*2100);
        }
        /* for times 2092.5-2100 */
        else if (t>=2092.5){
            mortality_rate_under5 += exp(param->mortality_rate_by_gender_age_intercept[g][0] + param->mortality_rate_by_gender_age_slope[g][0]*t);
            mortality_rate_5to10 += exp(param->mortality_rate_by_gender_age_intercept[g][1] + param->mortality_rate_by_gender_age_slope[g][1]*2100);
        }
        /* for times 1950-2092.5: */
        else{
            mortality_rate_under5 += exp(param->mortality_rate_by_gender_age_intercept[g][0] + param->mortality_rate_by_gender_age_slope[g][0]*t);
            mortality_rate_5to10 += exp(param->mortality_rate_by_gender_age_intercept[g][1] + param->mortality_rate_by_gender_age_slope[g][1]*(t+7.5));
        }
    }
    /* We take the average over genders for each rate: */
    mortality_rate_under5 = mortality_rate_under5/(N_GENDER*1.0);
    mortality_rate_5to10 = mortality_rate_5to10/(N_GENDER*1.0);

    /* Now combine to get an overall probability of mortality between birth and reaching adulthood. Each UNPD age group is 5 years, hence the power of 5. */
    mortality_rate_childhood_total = 1 - pow(1-mortality_rate_under5,5) * pow(1-mortality_rate_5to10,5);
    //printf("At time %6.4lf mortality_rate_childhood_total=%6.4lf\n",t,mortality_rate_childhood_total);
    return mortality_rate_childhood_total;
}


double natural_death_rate(int age, int g, parameters *param, double t){
   /* Return age- and gender- specific mortality rate for a given year
    
    Mortality depends upon age and gender of the group in question and also upon the year in 
    question.  Mortality rates are generated from UNPD mortality statistics.  These are adjusted
    to generate a natural (non-HIV related) mortality rate.  This function calculates the natural
    death rate (ie non-HIV related) using the  intercept and slope stored in the
    `mortality_rate_by_gender_age_intercept` array of the parameters structure.  
    
    UNPD mortality data is in 5 year age groups 0-4, 5-9, 10-14, 15-19, 20-24, 25-29, ... 
    which we index as 0, 1, 2.  
    
    
    Returns
    -------
    Probability of dying in the next calendar year (therefore need to multiply by timestep).
    
    Arguments
    ---------
    age : int
        Age of the individual
    g : int
        Gender of the individual
    param : pointer to param structure
        Parameter structure that includes mortality intercept and slope parameters
    t : double
        Time in decimal years
    */
    
    double mortality_rate;
    
    // The formula below gives the indexing for UNPD age-groups using integer division:
    int a_unpd = age / 5;
    
    
    // Check UNPD age-category is within a suitable range
    if((a_unpd > 16) || (a_unpd < 0)){
        printf("Error: Index for UNPD age-categories either too large or too small.\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    if(t < 1950){
        mortality_rate = exp(param->mortality_rate_by_gender_age_intercept[g][a_unpd] +
            param->mortality_rate_by_gender_age_slope[g][a_unpd] * 1950);
    }else if(t > 2100){
        mortality_rate = exp(param->mortality_rate_by_gender_age_intercept[g][a_unpd] +
            param->mortality_rate_by_gender_age_slope[g][a_unpd] * 2100);
    }else{
        mortality_rate = exp(param->mortality_rate_by_gender_age_intercept[g][a_unpd] +
            param->mortality_rate_by_gender_age_slope[g][a_unpd] * t);
    }
    return mortality_rate;
}


/* Function does: Assign new individuals to a specific risk group. As it stands the number of new individuals in a given risk class is fixed over time.
 * Function arguments: gender and pointer to param structure.
 * Function returns: the index of the risk group (currently 0,1,2). */
int draw_sex_risk(int gender, parameters *param){
    double x;
    x = gsl_rng_uniform (rng);
    if (x<=(param->initial_prop_gender_risk[gender][LOW]))
        return LOW;
    else if (x<=(param->initial_prop_gender_risk[gender][LOW]+param->initial_prop_gender_risk[gender][MEDIUM]))
        return MEDIUM;
    else 
        return HIGH;
}

/* Function does: creates entries for everything related to that new_person.
 * Function arguments: pointer to new person to be created, current time (for generating a DoB) and a pointer to the param structure (to get probabilities such as gender, MMC, etc).
 * Function returns: nothing. */
void create_new_individual(individual *new_adult, double t, parameters *param, int hivstatus, population_size_one_year_age *n_infected, patch_struct *patch, int p, all_partnerships *overall_partnerships){
    int i;

    new_adult->id = patch[p].id_counter;        /* Set the id to be the value of patch[p].id_counter. */

    new_adult->patch_no = p;

    // For debugging:
    if(new_adult->id==FOLLOW_INDIVIDUAL && new_adult->patch_no==FOLLOW_PATCH){
        printf("Creation of adult %ld in patch %d with hivstatus %i at t=%6.2f\n",new_adult->id,new_adult->patch_no,hivstatus,t);
        fflush(stdout);
    }

    /* Determine gender based on sex ratio parameter (which is kept fixed for all time): */
    if (gsl_ran_bernoulli(rng,(param->sex_ratio))==1)   /* Assume that the M/F sex ratio is unchanged over time. */
        new_adult->gender = MALE;
    else
        new_adult->gender = FEMALE;

    /* Assign a date of birth. Note: we currently assume that people do not enter the adult population aged 13.0, but instead 13.99, otherwise there are problems when ageing. 
     * As it is we ensure this way that someone who enters the population at the last timestep is aged 14.0 when they are aged to the next year-group one timestep later. */
    new_adult->DoB = t - AGE_ADULT - (N_TIME_STEP_PER_YEAR-1)/(1.0*N_TIME_STEP_PER_YEAR);      
    new_adult->DoD = -1;
    new_adult->DEBUGTOTALTIMEHIVPOS = 0;
    /* Assign a sex risk group: */
    new_adult->sex_risk = draw_sex_risk(new_adult->gender,param);  
    new_adult->n_lifetime_partners = 0;
    new_adult->n_lifetimeminusoneyear_partners = 0;
    new_adult->n_lifetime_partners_outside = 0;
    new_adult->n_lifetimeminusoneyear_partners_outside = 0;
    new_adult->n_partnersminusoneyear = 0;
    // For debugging:
    if(new_adult->id==FOLLOW_INDIVIDUAL && new_adult->patch_no==FOLLOW_PATCH)
        printf("New adult DoB = %f %li\n",new_adult->DoB,new_adult->id); 

    new_adult->time_to_delivery = -1;  /* Not pregnant when enters population. */

    new_adult->VISITEDBYCHIPS_TO_INIT_ART = FALSE;
    new_adult->VISITED_BY_CHIPS_THISROUND = FALSE;
    new_adult->NCHIPSVISITS = 0;

    new_adult->PC_cohort_index = -1; /* Not in PC cohort (for now). */

    /* Assign HIV status, allowing for the fact that some children may have had perinatal transmission (children are divided into HIV+/- at birth). 
     * Note that CHRONIC is 2 so need an if statement here. */
    if (hivstatus==0){
        new_adult->HIV_status = UNINFECTED;
        new_adult->ART_status = ARTNEG;
        new_adult->next_HIV_event = NOEVENT; /* Initialize at dummy value. */
        new_adult->next_cascade_event = NOEVENT; /* Initialize at dummy value. */
        new_adult->SPVL_num_G = 0;                  /* Initialize at dummy value. */
        new_adult->SPVL_num_E = 0;                  /* Initialize at dummy value. */
        new_adult->SPVL_infector = 0;                /* Initialize at dummy value. */
        new_adult->cd4 = CD4_UNINFECTED;                 /* Initialize at dummy value, here -1 */
        new_adult->SPVL_cat = -1;                            /* Initialize at dummy value. */
        new_adult->time_last_hiv_test = NEVERHIVTESTED;  /* Assume never previously tested. */
        new_adult->t_sc = -1;                            /* Initialize at dummy value. */
        new_adult->idx_hiv_pos_progression[0] = -1;     /* Initialize at dummy value. */
        new_adult->idx_hiv_pos_progression[1] = -1;     /* Initialize at dummy value. */
        new_adult->debug_last_hiv_event_index = -1;     /* Initialize at dummy value. */
        /* Note these two are also set/overwritten in initialize_first_cascade_event_for_new_individual().
         * However that function is only called if HIV testing has started, and we need to set them to dummy values (if doing several runs this is part of resetting the memory). */
        new_adult->idx_cascade_event[0] = -1;           /* Initialize at dummy value. */
        new_adult->idx_cascade_event[1] = -1;           /* Initialize at dummy value. */
        new_adult->debug_last_cascade_event_index = -1;     /* Initialize at dummy value. */
        new_adult->idx_vmmc_event[0] = -1;         /* Initialize at dummy value. */
        new_adult->idx_vmmc_event[1] = -1;      
        new_adult->debug_last_vmmc_event_index = -1;     /* Initialize at dummy value. */

        /* PANGEA stuff: */
        new_adult->PANGEA_t_prev_cd4stage = -1.0;
        new_adult->PANGEA_t_next_cd4stage = -1.0;
        new_adult->PANGEA_cd4atdiagnosis = -1.0;
        new_adult->PANGEA_cd4atfirstART = -1.0;
        new_adult->PANGEA_t_diag = -1.0;
        new_adult->PANGEA_date_firstARTstart = -1.0;
        new_adult->PANGEA_date_startfirstVLsuppression = -1.0;
        new_adult->PANGEA_date_endfirstVLsuppression = -1.0;

        /* Variables store cumulative amount of time a person spends on ART: */
        new_adult->DEBUG_cumulative_time_on_ART_VS = 0;
        new_adult->DEBUG_cumulative_time_on_ART_VU = 0;
        new_adult->DEBUG_cumulative_time_on_ART_early = 0;
        new_adult->DEBUG_time_of_last_cascade_event = -1; /* Dummy value. */
    }
    else{
        (n_infected->pop_size_per_gender_age1_risk[new_adult->gender][n_infected->youngest_age_group_index][new_adult->sex_risk]) += 1;
        printf("+++ One new HIV+ (new adult) \n");
        fflush(stdout);

        new_adult->HIV_status = CHRONIC;
        new_adult->ART_status = LTART_VS;                // Assume any new adult who has made it this far is successfully on ART 
        new_adult->SPVL_num_G = 0;                          // WRONG!!!
        new_adult->SPVL_num_E = 0;                          // WRONG!!!
        new_adult->SPVL_infector = 0;                       // WRONG!!!
        new_adult->cd4 = CD4_UNINFECTED;                 /// WRONG!!!
        new_adult->SPVL_cat = -1;                            /// WRONG!!!
        new_adult->time_last_hiv_test = t - AGE_ADULT;  ///  Assume that someone who survived this long was tested at birth. Not important anyway as we would be interested in ADULT testing in last year or less. 
        new_adult->t_sc = t - AGE_ADULT;                /* Seroconverted at birth. */
        new_adult->idx_vmmc_event[0] = -1;         /* Initialize at dummy value. */
        new_adult->idx_vmmc_event[1] = -1;

        new_adult->PANGEA_t_prev_cd4stage = -1.0;
        new_adult->PANGEA_t_next_cd4stage = -1.0;
        new_adult->PANGEA_cd4atdiagnosis = -1.0;
        new_adult->PANGEA_cd4atfirstART = -1.0;           /// WRONG?
        new_adult->PANGEA_t_diag = new_adult->time_last_hiv_test;
        new_adult->PANGEA_date_firstARTstart = -1.0;    /// WRONG!!!
        new_adult->PANGEA_date_startfirstVLsuppression = -1.0;
        new_adult->PANGEA_date_endfirstVLsuppression = -1.0;

        /* Assume that this person was on ART and virally suppressed their whole life minus an initial early ART period (otherwise the debug stats may look weird if someone is on ART but never had early ART): */
        new_adult->DEBUG_cumulative_time_on_ART_VS = AGE_ADULT - param[p].t_end_early_art;
        new_adult->DEBUG_cumulative_time_on_ART_VU = 0;
        new_adult->DEBUG_cumulative_time_on_ART_early = param[p].t_end_early_art;
        new_adult->DEBUG_time_of_last_cascade_event = t; /* Start counting additional time on ART from current time as we assume that childhood was spent VS (otherwise survival lower). */

        printf("Need to work out what is happening with new adults turning 13 who are HIV+\n"); 
        printf("Also need to schedule HIV and ART events for them. \n");
    }


    /* Set up partnerships later on, so initialize to zero here:. */
    new_adult->n_partners = 0;         
    new_adult->n_HIVpos_partners = 0;



    // At present set_max_n_partners() does not actually use age group or gender.
    /* set_max_n_partners() depends on gender, age group and risk group. This is a new adult so age group is 0. */
    new_adult->max_n_partners = set_max_n_partners(new_adult->gender, 0, new_adult->sex_risk, param);  

    /* Number of sexual partners outside cluster: */
    new_adult->n_partners_outside = 0;
    new_adult->n_HIVpos_partners_outside = 0;

    /* If male, decide if circumcised (as a child, not by trial) here: */
    if (new_adult->gender==MALE){
        if (gsl_ran_bernoulli(rng,(param->p_child_circ))==1)
            new_adult->circ = TRADITIONAL_MC;
        else
            new_adult->circ = UNCIRC;
    }
    else
        new_adult->circ = 0;   /* Women - set to zero. */   


    new_adult->idx_serodiscordant = -1;  /* Not in a serodiscordant partnership */


    /* Add all the available partnerships (max_n_partners as they do not have any yet) to the list of available partnerships
     * and create references to these in the new_adult individual structure. */
    for(i=new_adult->n_partners ; i<new_adult->max_n_partners ; i++){
        /* Add to end of pop_available_partners array element: */ 
        new_adult->idx_available_partner[i] = overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[new_adult->gender][0][new_adult->sex_risk]; /* Not yet in the list of available partners */
        /* Note that age group = 0 as new adults. */
        overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[p][new_adult->gender][0][new_adult->sex_risk][overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[new_adult->gender][0][new_adult->sex_risk]] = new_adult;
        overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[new_adult->gender][0][new_adult->sex_risk]++;
    }

    /* Above max_n_partners set the rest of the entries to -1 - this is a checking mechanism to ensure we never give them more than max_n_partners. */
    for(i=new_adult->max_n_partners ; i<MAX_PARTNERSHIPS_PER_INDIVIDUAL ; i++){
        new_adult->idx_available_partner[i] = -1; /* Not in the list of available partners */
    }

    /*if(new_adult->id==5486)
    {
            print_here_string("000000000000000000000000000000000000000000000000",0);
            print_here_string("Created individual ",new_adult->id);
            print_here_string("Patch ",new_adult->patch_no);
            print_here_string("Gender ",new_adult->gender);
            print_here_string("Circ ",new_adult->circ);
            print_here_string("Risk ",new_adult->sex_risk);
            print_here_string("Yearob ",(int) new_adult->DoB);
            print_here_string("Max_n_partners ",new_adult->max_n_partners);
            print_here_string("Current_n_partners ",new_adult->n_partners);
            print_here_string("111111111111111111111111111111111111111111111111",1);
    }*/

}


void initialize_first_cascade_event_for_new_individual(individual *new_adult, double t, 
    parameters *param, individual ***cascade_events, long *n_cascade_events, 
    long *size_cascade_events){
    /* Schedule HIV cascade events for individuals that transition from childhood to adulthood
    
    There are two types of HIV testing schedules in the model (determined by the macro called 
    `HIVTESTSCHEDULE`).  The first situation (where HIVTESTSCHEDULE == 0) means that individuals
    are scheduled HIV tests sequentially, the second situation (where HIVTESTSCHEDULE == 1) means
    that test scheduling procedure is performed for the whole population at fixed times.  
    
    Arguments
    ---------
    new_adult : pointer to an individual structure
        Individual for whom the test is to be scheduled
    t : double
        Current time in years
    param : pointer to a parameters structure
        Structure that stores all parameter values
    cascade_events : pointer to a multidimensional array of individuals
    n_cascade_events : pointer to a long
    size_cascade_events : pointer to a long
    */
    
    
    // Check the new adult is uninfected with HIV
    if(new_adult->HIV_status == UNINFECTED){
        // If each individual schedules their HIV tests sequentially draw a time for this person.
        if(HIVTESTSCHEDULE == 0){
            schedule_new_hiv_test(new_adult, param, t, cascade_events, n_cascade_events,
                size_cascade_events);
        
        // Otherwise, the test scheduling procedure happens for the whole population at fixed times.
        }else{
            // For a new adult (about to turn 14) assume that they won't get an HIV test until the
            // next scheduled time.  
            new_adult->next_cascade_event = NOEVENT;
            new_adult->idx_cascade_event[0] = NOEVENT;
            new_adult->idx_cascade_event[1] = -1;
        }
    }else{
        printf("What is the next cascade event if born HIV+?\n");
        // This is an assumption, change this code when we decide out what happens to people who
        // have been HIV+ since birth.  
        new_adult->next_cascade_event = NOEVENT;
        new_adult->idx_cascade_event[0] = -1;
        new_adult->idx_cascade_event[1] = -1;
    }
}


/* Function arguments: pointer to an individual new_adult, pointer to a population_size n_population
 * Function does: updates the population_size according to the new_adult after his/her birth
 * Function returns: nothing. */
void update_population_size_new_adult(individual *new_adult, population_size *n_population, population_size_one_year_age *n_population_oneyearagegroups,
        stratified_population_size *n_population_stratified){
    /* Add to first age group ag=0. */
    (n_population->pop_size_per_gender_age_risk[new_adult->gender][0][new_adult->sex_risk])++;
    (n_population_stratified->pop_size_per_gender_age[new_adult->gender][0])++;
    (n_population_stratified->pop_size_per_gender_risk[new_adult->gender][new_adult->sex_risk])++;
    (n_population_stratified->total_pop_size_per_gender[new_adult->gender])++;

    /* Now overall population: */
    //(n_population_stratified->pop_size_per_age_risk[0][new_adult->sex_risk])++;
    //(n_population_stratified->pop_size_per_age[0])++;
    //(n_population_stratified->pop_size_per_risk[new_adult->sex_risk])++;
    (n_population_stratified->total_pop_size)++;

    n_population_oneyearagegroups->pop_size_per_gender_age1_risk[new_adult->gender][n_population_oneyearagegroups->youngest_age_group_index][new_adult->sex_risk] += 1;

    int r, g; //// arguably we will do this many times so maybe write as a separate inline function?
    for (g=0; g<N_GENDER; g++){
        for (r=0; r<N_RISK; r++){
            n_population_stratified->prop_pop_per_gender_risk[g][r] = n_population_stratified->pop_size_per_gender_risk[g][r]/n_population_stratified->total_pop_size_per_gender[g];
        }
    }


}

/* Function does: Updates population size structure n_population when an individual dies.  
 * Function arguments: pointer to the specific individual who dies, pointer to population_size structure, age group index (age groups 13-18, 19-22, 23-30, etc). 
 * Function returns: Nothing. */
void update_population_size_death(individual *individual, population_size *n_population, population_size_one_year_age *n_population_oneyearagegroups,
        population_size_one_year_age *n_infected, stratified_population_size *n_population_stratified, int aa, age_list_struct *age_list){

    int ag = FIND_AGE_GROUPS[aa];
    int ai;
    if (PRINT_DEBUG_DEMOGRAPHICS)
        printf("Dead adult: ID = %li DoB = %f gender = %i risk = %i age gp =%i\n",individual->id,individual->DoB,individual->gender,individual->sex_risk,ag);

    (n_population->pop_size_per_gender_age_risk[individual->gender][ag][individual->sex_risk])--;

    if (aa<MAX_AGE-AGE_ADULT){
        ai = n_population_oneyearagegroups->youngest_age_group_index + aa;
        while (ai>(MAX_AGE-AGE_ADULT-1))
            ai = ai - (MAX_AGE-AGE_ADULT);
        n_population_oneyearagegroups->pop_size_per_gender_age1_risk[individual->gender][ai][individual->sex_risk] -=1;
    }
    else{
        n_population_oneyearagegroups->pop_size_oldest_age_group_gender_risk[individual->gender][individual->sex_risk] -=1;
    }

    /* Remove from prevalent cases if HIV+. */
    if (individual->HIV_status>UNINFECTED){
        if (aa<MAX_AGE-AGE_ADULT){
            ai = n_infected->youngest_age_group_index + aa; /* ai is the index of the two arrays age_list->number_per_age_group and age_list->age_group */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);
            // FOR DEBUGGING ONLY:
            //check_age_group_index(age_list, individual->gender, individual->id,ai); /// could be removed as only checking things ///
            (n_infected->pop_size_per_gender_age1_risk[individual->gender][ai][individual->sex_risk]) -= 1;
            //printf("--- One death of HIV+ (age group %d, gender %d risk group %d)\n",ai, individual->gender, individual->sex_risk);
            //fflush(stdout);
        }
        else{
            (n_infected->pop_size_oldest_age_group_gender_risk[individual->gender][individual->sex_risk]) -= 1;
            //printf("--- One death of HIV+ (old, gender %d risk group %d)\n",individual->gender, individual->sex_risk);
            //fflush(stdout);
        }
    }

    (n_population_stratified->pop_size_per_gender_age[individual->gender][ag])--;
    (n_population_stratified->pop_size_per_gender_risk[individual->gender][individual->sex_risk])--;
    (n_population_stratified->total_pop_size_per_gender[individual->gender])--;

    /* Now overall population: */
    //(n_population_stratified->pop_size_per_age_risk[ag][individual->sex_risk])--;
    //(n_population_stratified->pop_size_per_age[ag])--;
    //(n_population_stratified->pop_size_per_risk[individual->sex_risk])--;
    (n_population_stratified->total_pop_size)--;

    int r, g; //// arguably we will do this many times so maybe write as a separate inline function?
    for (g=0; g<N_GENDER; g++){
        for (r=0; r<N_RISK; r++){
            if (n_population_stratified->total_pop_size_per_gender[g]>0)
                n_population_stratified->prop_pop_per_gender_risk[g][r] = n_population_stratified->pop_size_per_gender_risk[g][r]/n_population_stratified->total_pop_size_per_gender[g];
            else
                n_population_stratified->prop_pop_per_gender_risk[g][r] = 0;
        }
    }

}

/* Function arguments: pointer to age_list, pointer to the new individual who has entered adult population.
 * Function does: Updates age_list when a new adult enters the adult population (from the child population). */
void update_age_list_new_adult(age_list_struct *age_list, individual *individual_ptr){
    int g = individual_ptr->gender;

    int yi = age_list->age_list_by_gender[g]->youngest_age_group_index; /* Temporary store so we don't have to keep referring to this index the long way. */

    /* Add the pointer to the new adult to the youngest age group: */
    age_list->age_list_by_gender[g]->age_group[yi][age_list->age_list_by_gender[g]->number_per_age_group[yi]] = individual_ptr;

    /* Adds one to the count of youngest age group. */
    (age_list->age_list_by_gender[g]->number_per_age_group[yi])++;
}


/* Function arguments: pointer to age_list, a=age of individual to die,
 *                      new_death = the index of the person to die within its age group, t=time
 * Function does: Updates age_list when an adult leaves the adult population (dying from natural or HIV-related causes).
 * Function returns: Nothing. */
void update_age_list_death(age_list_struct *age_list, int g, int aa, long new_death, double t, int p){
    if(g<0||g>1){
        printf("ERROR: UNKNOWN GENDER!!!!!!\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* We have divided up the population in age_list so that individuals aged 13-79 are kept in arrays (by year) within
     * age_list->age_group. People aged 80+ are kept separately and dealt with in the next part of the if statement. */
    if (aa<(MAX_AGE-AGE_ADULT)){
        //if (PRINT_DEBUG_DEMOGRAPHICS)
        //  printf("Getting rid of: %f %li aa=%i, n=%li. ",(age_list->age_list_by_gender[g]->age_group)[aa][new_death]->DoB,(age_list->age_list_by_gender[g]->age_group)[aa][new_death]->id,aa,new_death);

        /* We always want to keep the age_list->age_group[aa] arrays ordered so that the first age_list->number_per_age_group[aa]
         * individuals are still alive (anything after this point can be a dead person or uninitialized as we should
         * never access beyond that). To do this we just swap the last person in the list (who is still alive) with
         * the person who just died. If the last person in the list is the dead person it does not matter because we 
         * also decrease  number_per_age_group by 1 so they are moved outside the list of alive people in any case. */
        (age_list->age_list_by_gender[g]->age_group)[aa][new_death] = (age_list->age_list_by_gender[g]->age_group)[aa][age_list->age_list_by_gender[g]->number_per_age_group[aa]-1];

        (age_list->age_list_by_gender[g]->number_per_age_group[aa])--;

        //if (PRINT_DEBUG_DEMOGRAPHICS && (new_death<(age_list->age_list_by_gender[g]->number_per_age_group[aa])))
        //  printf("Swapped to: %f %li\n",(age_list->age_list_by_gender[g]->age_group)[aa][new_death]->DoB,(age_list->age_list_by_gender[g]->age_group)[aa][new_death]->id);
    }
    /* Now deal with individuals who are aged 80+ - these are kept in a separate array (oldest_age_group): */
    else{
        if (PRINT_DEBUG_DEMOGRAPHICS)
            printf("Getting rid of: %f %li gender=%i. ",(age_list->age_list_by_gender[g]->oldest_age_group)[new_death]->DoB,(age_list->age_list_by_gender[g]->oldest_age_group)[new_death]->id,g);

        /* This is the same swap as above. */
        (age_list->age_list_by_gender[g]->oldest_age_group)[new_death] = (age_list->age_list_by_gender[g]->oldest_age_group)[age_list->age_list_by_gender[g]->number_oldest_age_group-1];
        (age_list->age_list_by_gender[g]->number_oldest_age_group)--;

        if (PRINT_DEBUG_DEMOGRAPHICS && (new_death<(age_list->age_list_by_gender[g]->number_oldest_age_group)))
            printf("Swapped to: %f %li %i\n",(age_list->age_list_by_gender[g]->oldest_age_group)[new_death]->DoB,(age_list->age_list_by_gender[g]->oldest_age_group)[new_death]->id,(age_list->age_list_by_gender[g]->oldest_age_group)[new_death]->gender);
    }

}



int get_age_index(double DoB, double start_simul){ /// Do we ever use this?

    int ai = ( (int) floor(start_simul - DoB)) - AGE_ADULT;
    /* Here we MUST use a while loop instead of an if statement as if someone is born in 2100, then ai is still negative if we just do this once. */
    while (ai<0)
        ai += (MAX_AGE-AGE_ADULT);
    return ai;
}

//int get_age_indexv2(double DoB, double t){ /// Do we ever use this?
//
//  //int ai = ( (int) floor(start_simul - DoB)) - AGE_ADULT;
//  /* Here we MUST use a while loop instead of an if statement as if someone is born in 2100, then ai is still negative if we just do this once. */
//  while (aa<0)
//      aa += (MAX_AGE-AGE_ADULT);
//  return aa;
//}


int get_age_indexv2(double DoB, double t, int youngest_age_group_index){ /// Do we ever use this? Not currently used except as a check - may use as a 'standard method' later on. */
    int ai;
    int aa = (int) floor(floor(t) - DoB) - AGE_ADULT;

    if (aa<(MAX_AGE-AGE_ADULT)){
        ai = youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */

        while (ai>(MAX_AGE-AGE_ADULT-1))
            ai = ai - (MAX_AGE-AGE_ADULT);
    }
    else{
        //printf("OLD\n");
        ai=999;
    }
    return ai;
}


int get_age_group(double DoB, double t, const int age_groups[], int number_age_groups){
   /* Find index `ag` of AGE_GROUPS[] array to which someone with a DoB belongs at the beginning of
    the year of year t.  
    
    This function assumes nobody will be less than age_groups[0] since no individuals within the 
    IBM are generated with in that age group.  Note that age groups are generally 13-17, 18-22,
    etc, and defined within constants.c.  
    
    Arguments
    ---------
    DoB : double
        Date of birth of an individual (in decimal years)
    t : double
        Current time (in decimal years)
    age_groups : array of int
        Bins defining the age groups of interest
    number_age_groups : int
        The number of age groups (i.e. the length of the array age_groups)
    
    Returns
    -------
    ag : int
        Index of array age_groups[] in which the individual's (rounded down) age fits.  Any age
        greater than the largest age group will be counted in the final age group.  
    */
    double age = floor(t) - DoB;
    
    int ag = 0;
    
    if(age < age_groups[number_age_groups - 1]){
        while(age_groups[ag + 1] <= age){
            ag++;
        }
    }else{
        ag = number_age_groups - 1;
    }
    return ag;
}


int get_age_group_unpd(double DoB, double t){
    double age = floor(t) - DoB;
    int ag=0;
    //printf("age=%6.4lf %i\n",age,AGE_GROUPS_UNPD[N_AGE_UNPD]);
    if (age<AGE_GROUPS_UNPD[N_AGE_UNPD])
        while (AGE_GROUPS_UNPD[ag+1]<=age)
            ag++;
    else
        return  N_AGE_UNPD;
    return ag;
}


/* Function updates n_population as it ages by one year (ageing is by cohort). */
void update_n_population_ageing_by_one_year(patch_struct *patch, int p){
    /* aa+AGE_ADULT is the age of the person (so aa runs from 0..MAX_AGE-AGE_ADULT).
     * ai is the corresponding row index for them in age_list->age_group[ai][] and age_list->number_per_age_group[ai];
     * age_index is the index in the array AGE_GROUPS[] (which runs up to N_AGE). */
    int age_index, aa, ai, n;
    int n_age_ai;  /* Number of people with age index ai. */
    int g;

    for (g=0;g<N_GENDER;g++){
        /* Deliberately starting at age_index=1 - we are interested in ageing the population by one year, but those aged 12 turning 13 are dealt with separately as new adults. */
        for (age_index=1; age_index<N_AGE; age_index++){
            /* We are interested in those about to transition to the next age group - so take -1 as this is the age of those about to transition. */
            aa = AGE_GROUPS[age_index]-AGE_ADULT-1;
            ai = patch[p].age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);
            n_age_ai = patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai];
            /* Modify n_population for the given individuals of age aa+AGE_ADULT;
             * Note that these transitions are 17->18, 22->23, 30->31, 40->41, 50->51, 60->61. */
            for (n=0; n<n_age_ai; n++){
                patch[p].n_population->pop_size_per_gender_age_risk[g][age_index-1][patch[p].age_list->age_list_by_gender[g]->age_group[ai][n]->sex_risk]--;
                patch[p].n_population->pop_size_per_gender_age_risk[g][age_index][patch[p].age_list->age_list_by_gender[g]->age_group[ai][n]->sex_risk]++;
            }
        }
    }
}



/* Function updates pop_available_partners and n_pop_available_partners as population ages by one year (ageing is by cohort). 
 * It does this by going through each individual who is about to age (NOTE: it is important that this is called out before age_list is updated)
 * For each individual we go through each of their available partnerships and move each one to the new age group as needed. 
 * Code can probably be sped up - the issue with making pop_available_partners into 1 years age groups is that the number in each age group
 * will be small so more chance for having two identical partnerships formed (not clear how to prevent this). */
void update_pop_available_partners_ageing_by_one_year(patch_struct *patch, int p, all_partnerships *overall_partnerships, double t){
    /* aa+AGE_ADULT is the age of the person (so aa runs from 0..MAX_AGE-AGE_ADULT).
     * ai is the corresponding row index for them in age_list->age_group[ai][] and age_list->number_per_age_group[ai];
     * age_index is the index in the array AGE_GROUPS[] (which runs up to N_AGE). */
    int age_index, aa, ai, n, i, i2;
    int g, r;
    int n_age_ai;  /* Number of people with age index ai. */
    individual *this_person;
    individual *personB;

    //// Debug:
    //for 
    //printf("update_pop_available_partners_ageing_by_one_year: %li\n",n_pop_available_partners->pop_size_per_gender_age_risk[0][1][0]);
    //printf("Check this person is here: %li\n",pop_available_partners->pop_per_gender_age_risk[0][1][0][(n_pop_available_partners->pop_size_per_gender_age_risk[0][1][0])-1]->id);
    //fflush(stdout);
    for (g=0;g<N_GENDER;g++){
        /* Deliberately starting at age_index=1 - we are interested in ageing the population by one year, but those aged 12 turning 13 are dealt with separately as new adults. */
        for (age_index=1; age_index<N_AGE; age_index++){

            /* We are interested in those about to transition to the next age group - so take -1 as this is the age of those about to transition. */
            aa = AGE_GROUPS[age_index]-AGE_ADULT-1;

            ai = patch[p].age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);
            n_age_ai = patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai];
            if (PRINT_DEBUG_DEMOGRAPHICS){
                if (age_index==1){
                    printf("Check age is an 18 year old: %f\n",t-patch[p].age_list->age_list_by_gender[g]->age_group[ai][0]->DoB);
                    printf("Check age is an 18 year old: %f %f\n",t-patch[p].age_list->age_list_by_gender[g]->age_group[ai][patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai]-1]->DoB,patch[p].age_list->age_list_by_gender[g]->age_group[ai][patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai]-1]->DoB);
                }
            }
            /* Modify n_population for the given individuals of age aa+AGE_ADULT;
             * Note that these transitions are 17->18, 22->23, 30->31, 40->41, 50->51, 60->61. */

            for (n=0; n<n_age_ai; n++){
                /* Use the pointer as an alias for this person - makes code more readable + possibly a bit quicker? */
                this_person = patch[p].age_list->age_list_by_gender[g]->age_group[ai][n];
                /*if (this_person->id==10013)
                printf("10013 moving out of age group: %i with DoB %f\n",age_index-1,this_person->DoB);*/

                if(this_person->id==FOLLOW_INDIVIDUAL && this_person->patch_no==FOLLOW_PATCH)
                {
                    printf("Individual %ld gender %i from patch %d is aged: age group %d \n",this_person->id,g,this_person->patch_no,age_index-1);
                    fflush(stdout);
                }

                r = this_person->sex_risk;
                /*if(p!=this_person->patch_no)
                {
                    printf("Issue in update_pop_available_partners_ageing_by_one_year\n");
                    fflush(stdout);
                }*/

                //print_here_string("Check patch inside function",p);

                //printf("Current g a r = %i %i %i\n",g,age_index,r);
                //printf("XXupdate_pop_available_partners_ageing_by_one_year: %li\n",n_pop_available_partners->pop_size_per_gender_age_risk[0][1][0]);
                //printf("Check this person is here: ID=%li CD4=%i idx=%i\n",pop_available_partners->pop_per_gender_age_risk[0][1][0][(n_pop_available_partners->pop_size_per_gender_age_risk[0][1][0])-1]->id,pop_available_partners->pop_per_gender_age_risk[0][1][0][(n_pop_available_partners->pop_size_per_gender_age_risk[0][1][0])-1]->cd4,pop_available_partners->pop_per_gender_age_risk[0][1][0][(n_pop_available_partners->pop_size_per_gender_age_risk[0][1][0])-1]->idx_available_partner[0]);
                //fflush(stdout);
                /* We go over all available partnerships of this person. */
                i = 0;
                //printf("This person: %i %i %f\n",g,r,this_person->DoB);
                while ((this_person->idx_available_partner[i]>-1) && (i<(this_person->max_n_partners-this_person->n_partners)) && (overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][age_index-1][r]>0)){
                    /* Swap this person's available index with that of the last person in the array (we'll call them person B):
                     * Firstly swap the pointer with the pointer of person B: */
                    //printf("HEY2x %li %i %i %i\n",n_pop_available_partners->pop_size_per_gender_age_risk[g][age_index-1][r],g,age_index-1,r);
                    //printf("Last person: %i %i %f %li\n",pop_available_partners->pop_per_gender_age_risk[g][age_index-1][r][(n_pop_available_partners->pop_size_per_gender_age_risk[g][age_index-1][r])-1]->gender,pop_available_partners->pop_per_gender_age_risk[g][age_index-1][r][n_pop_available_partners->pop_size_per_gender_age_risk[g][age_index-1][r]-1]->sex_risk,pop_available_partners->pop_per_gender_age_risk[g][age_index-1][r][n_pop_available_partners->pop_size_per_gender_age_risk[g][age_index-1][r]-1]->DoB,pop_available_partners->pop_per_gender_age_risk[g][age_index-1][r][(n_pop_available_partners->pop_size_per_gender_age_risk[g][age_index-1][r])-1]->id);
                    personB = overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[p][g][age_index-1][r][overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][age_index-1][r]-1];

                    overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[p][g][age_index-1][r][this_person->idx_available_partner[i]] = personB;

                    /* Now adjust the idx_available_partner element of person B to reflect this change in pop_available_partners:
                     * Unfortunately we have to look through their partnerships to find the right one: */
                    i2 = personB->max_n_partners-personB->n_partners-1;
                    while ((personB->idx_available_partner[i2]!=overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][age_index-1][r]-1) && (personB->idx_available_partner[i2]>=0) && (i2>=0)){
                        //printf("HEY %i %i %li\n",i2,personB->idx_available_partner[i2],n_pop_available_partners->pop_size_per_gender_age_risk[g][age_index-1][r]-1);
                        if ((personB->idx_available_partner[i2]==-1)|| (i2<0)){
                            printf("Can't find person B's index in update_pop_available_partners_ageing_by_one_year\n");
                            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                            fflush(stdout);
                            exit(1);
                        }
                        i2--;
                    }
                    //  if (this_person->id==10013){
                    //  printf("Swapping 10013 idx = %i (%li) with %li's idx=%i (%li): \n",i, this_person->idx_available_partner[i],personB->id,i2,personB->idx_available_partner[i2]);
                    //  printf("Now 10013's indices are: %li %li %li %li %li %li %li %li %li \n",this_person->idx_available_partner[0],this_person->idx_available_partner[1],this_person->idx_available_partner[2],this_person->idx_available_partner[3],this_person->idx_available_partner[4],this_person->idx_available_partner[5],this_person->idx_available_partner[6],this_person->idx_available_partner[7],this_person->idx_available_partner[8]);
                    //}
                    ////* Note that in general we can swap indices from the same person. The only problem comes if we are trying to swap the same index (which should be because they are the last person).
                    /// * In that case we actually don't need to do anything.*/

                    //if (personB->idx_available_partner[i2]==this_person->idx_available_partner[i]){



                    /* Now we've got the correct index i2 set this to point to the new place in pop_available_partners (ie where this_person was): */
                    personB->idx_available_partner[i2] = this_person->idx_available_partner[i];

                    /* Next decrease the number of available partners in age_index-1: */
                    overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][age_index-1][r]--;

                    /* Add ageing person into the pop_available_partners for new age group (note that because we are adding to the end, no need for "-1" in the last index: */
                    overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[p][g][age_index][r][overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][age_index][r]] = this_person;
                    /* Modify idx_available_partner of this person: */
                    this_person->idx_available_partner[i] = overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][age_index][r];

                    /* Finally increase the number of available partners in this age group: */
                    overall_partnerships->n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][age_index][r]++;
                    /* Now go to next potential available partner: */
                    i++;
                }

            }
        }
    }
}


/* Function arguments: pointer to the age_list structure (which essentially contains lists of individuals in each group). 
 * Function does: moves the pointers for each age group by 1, moves MAX_AGE-1 age people into MAX_AGE group.
 * Function returns: nothing. */
void age_population_by_one_year(age_list_struct *age_list){
    int number_age_MAX_AGEminusone;
    int n;
    int g;


    for (g=0;g<N_GENDER;g++){
        /* If we have not reached the start of the array, move backwards to the previous element in the array. */
        if ((age_list->age_list_by_gender[g]->youngest_age_group_index) > 0){

            //printf("Oldest age group person:%f\n",age_list->age_list_by_gender[g]->oldest_age_group[(age_list->age_list_by_gender[g]->number_oldest_age_group)-1]->DoB);

            /* Move people aged MAX_AGE-1 into age_list->oldest_age_group[]. */
            ////1:/* Note: ((age_list->age_list_by_gender[g]->youngest_number_per_age_group_ptr)-1) should be a pointer to the previous group to the youngest age group (ie the group aged MAX_AGE-1). */
            number_age_MAX_AGEminusone = age_list->age_list_by_gender[g]->number_per_age_group[age_list->age_list_by_gender[g]->youngest_age_group_index-1];

            /* sending individuals aged 79 into 80+ */
            for (n=0; n<number_age_MAX_AGEminusone; n++){
                /* Copy pointers from one array of pointers to another. */
                /* What this does is it copies all the pointers of people turning AGE_MAX into the oldest_age_group array of pointers. */
                ////1:age_list->oldest_age_group[(age_list->number_oldest_age_group)+n] = (age_list->youngest_age_group-1)[n];
                age_list->age_list_by_gender[g]->oldest_age_group[(age_list->age_list_by_gender[g]->number_oldest_age_group)+n] = age_list->age_list_by_gender[g]->age_group[age_list->age_list_by_gender[g]->youngest_age_group_index-1][n];


                /* Make the pointer in the individual structure point to this place in the age_list. */
                //// Not necessary?:
                ////(age_list->oldest_age_group[(age_list->number_oldest_age_group)+n])->age_list_ptr = age_list->oldest_age_group[(age_list->number_oldest_age_group)+n];
                //if (PRINT_DEBUG_DEMOGRAPHICS)
                //  printf("Moving age group DoB = %f\n",(age_list->age_list_by_gender[g]->age_group[age_list->age_list_by_gender[g]->youngest_age_group_index-1][n])->DoB);
            }

            /* Update count in age_list->oldest_age_group[]. */
            age_list->age_list_by_gender[g]->number_oldest_age_group += number_age_MAX_AGEminusone;

            /* We have just removed everyone from the MAX_AGE-1 age group. This now becomes the counter for the youngest age group, so set to zero. */
            if (PRINT_DEBUG_DEMOGRAPHICS)
                printf("Number gender %i in youngest age group was %li, is %li\n",g, age_list->age_list_by_gender[g]->number_per_age_group[age_list->age_list_by_gender[g]->youngest_age_group_index],age_list->age_list_by_gender[g]->number_per_age_group[age_list->age_list_by_gender[g]->youngest_age_group_index-1]);

            ////1:*((age_list->youngest_number_per_age_group_ptr) -1) = 0;
            age_list->age_list_by_gender[g]->number_per_age_group[age_list->age_list_by_gender[g]->youngest_age_group_index-1] = 0;


            /* Move the pointer for the youngest age group to the start of the array. */
            (age_list->age_list_by_gender[g]->youngest_age_group_index)--;

            ////1:(age_list->youngest_age_group)--;
            ////1:(age_list->youngest_number_per_age_group_ptr)--;
            ////1:age_list->youngest_number_per_age_group_index--;
            /* Note: we probably ought to set the pointers in age_list for what used to be the MAX_AGE-1 group to NULL, but it's not essential as the count should tell us that there should be nobody there. */
            if (PRINT_DEBUG_DEMOGRAPHICS)
                if (age_list->age_list_by_gender[g]->number_per_age_group[age_list->age_list_by_gender[g]->youngest_age_group_index]>0)
                    printf("First entry in former young person list = %li %i %f\n",age_list->age_list_by_gender[g]->age_group[age_list->age_list_by_gender[g]->youngest_age_group_index][0]->id,g,age_list->age_list_by_gender[g]->age_group[age_list->age_list_by_gender[g]->youngest_age_group_index][0]->DoB);
            ////1:printf("First entry in former young person list = %i %i %f\n",        (age_list->youngest_age_group)[0]->id,(age_list->youngest_age_group)[0]->gender,(age_list->youngest_age_group)[0]->DoB);

        }
        else{
            /* Move people aged MAX_AGE-1 into age_list->oldest_age_group[]. */
            number_age_MAX_AGEminusone = age_list->age_list_by_gender[g]->number_per_age_group[MAX_AGE-AGE_ADULT-1];
            if (PRINT_DEBUG_DEMOGRAPHICS)
                printf("number to move to oldest age group = %i\n",number_age_MAX_AGEminusone);
            for (n=0; n<number_age_MAX_AGEminusone; n++){
                /* Copy pointers from one array of pointers to another. */
                if (PRINT_DEBUG_DEMOGRAPHICS)
                    printf("Moving age group DoB = %f\n",(age_list->age_list_by_gender[g]->age_group[MAX_AGE-AGE_ADULT-1][n]->DoB));
                age_list->age_list_by_gender[g]->oldest_age_group[(age_list->age_list_by_gender[g]->number_oldest_age_group)+n] = age_list->age_list_by_gender[g]->age_group[MAX_AGE-AGE_ADULT-1][n];
                /* Make the pointer in the individual structure point to this place in the age_list. */
                //// Not necessary?:
                ////(age_list->oldest_age_group[(age_list->number_oldest_age_group)+n])->age_list_ptr = age_list->oldest_age_group[(age_list->number_oldest_age_group)+n];
            }


            /* Update count in age_list->oldest_age_group[]. */
            age_list->age_list_by_gender[g]->number_oldest_age_group += age_list->age_list_by_gender[g]->number_per_age_group[MAX_AGE-AGE_ADULT-1];

            /* We have just removed everyone from the MAX_AGE-1 age group. This now becomes the counter for the youngest age group, so set to zero. */
            age_list->age_list_by_gender[g]->number_per_age_group[MAX_AGE-AGE_ADULT-1] = 0;

            /* Move the pointer for the youngest age group to the end of the array. */
            ////1:age_list->youngest_age_group = age_list->age_group[MAX_AGE-AGE_ADULT-1];
            age_list->age_list_by_gender[g]->youngest_age_group_index = MAX_AGE-AGE_ADULT-1;
            ////1:(age_list->youngest_number_per_age_group_ptr) = &(age_list->number_per_age_group[MAX_AGE-AGE_ADULT-1]);
            /* Note: we probably ought to set the pointers in age_list for what used to be the MAX_AGE-1 group to NULL, but it's not essential as the count should tell us that there should be nobody there. */
            ////1:age_list->youngest_number_per_age_group_index = MAX_AGE-AGE_ADULT-1;

        }
    }

}


//////////////////////
/* Function arguments: pointer to the population_size_one_year_age structure
 * Function does: moves the pointers for each age group by 1, moves MAX_AGE-1 age people into MAX_AGE group.
 * Function returns: nothing. */
void age_population_size_one_year_age_by_one_year(population_size_one_year_age *n_local_pop){
    int g,r;

    /* If we have not reached the start of the array, move backwards to the previous element in the array. */
    if ((n_local_pop->youngest_age_group_index) >= 1){
        for (g=0; g<N_GENDER; g++){
            for (r=0; r<N_RISK; r++){
                /* Merge people aged 79 (who are turning 80 now) into the 80+ year-age group: */
                n_local_pop->pop_size_oldest_age_group_gender_risk[g][r] += n_local_pop->pop_size_per_gender_age1_risk[g][n_local_pop->youngest_age_group_index-1][r]; 

                /* As everyone has aged by 1 year there are no people currently HIV+ aged 13: */
                n_local_pop->pop_size_per_gender_age1_risk[g][n_local_pop->youngest_age_group_index-1][r] = 0;
            }
        }
        /* Move the index for the youngest age group one back. */
        (n_local_pop->youngest_age_group_index)--;
    }
    else if ((n_local_pop->youngest_age_group_index)==0){

        for (g=0; g<N_GENDER; g++){
            for (r=0; r<N_RISK; r++){       
                /* Merge people aged 80+ into the 79 year-age group (so that the 79 group becomes the new 80+ group): */
                n_local_pop->pop_size_oldest_age_group_gender_risk[g][r] += n_local_pop->pop_size_per_gender_age1_risk[g][MAX_AGE-AGE_ADULT-1][r];

                /* As everyone has aged by 1 year there are no people currently HIV+ aged 13: */
                n_local_pop->pop_size_per_gender_age1_risk[g][MAX_AGE-AGE_ADULT-1][r] = 0;
            }
        }
        /* Move the index for the youngest age group to the right-hand end of the array: */
        (n_local_pop->youngest_age_group_index) = MAX_AGE-AGE_ADULT-1;
    }
    else{
        printf("Error: n_local_pop ageing process is not working!\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}



/* Function arguments: pointer to the person who died. 
 * Function does: removes a dead person from the list of susceptible_in_serodiscordant_partnership 
 * and n_susceptible_in_serodiscordant_partnership structs.
 * Only call this function if dead_person->idx_serodiscordant is >=0.
 * Function returns: nothing. */        
void remove_dead_person_from_susceptible_in_serodiscordant_partnership(individual *dead_person, 
        individual **susceptible_in_serodiscordant_partnership, long *n_susceptible_in_serodiscordant_partnership){
    int n,i;
    individual *a_partner;

    //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership line",731);

    /* We only need to update serodiscordant partnerships when the dead individual is seropositive. */ 
    if ((dead_person->HIV_status)>UNINFECTED){

        //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership in the if",0);

        if(dead_person->id==FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH){
            printf("Individual %ld from patch %d is dying - removing partners from susceptible in serodiscordant partnerships\n",dead_person->id,dead_person->patch_no);
            fflush(stdout);
        }

        //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership in the if",1);

        for (n=0; n<dead_person->n_partners; n++){

            //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership in the if in the for loop, n = ",n);

            //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership",0);

            a_partner = dead_person->partner_pairs[n]->ptr[1-dead_person->gender];

            //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership",1);

            /* Only need to deal with the partner if the partner is HIV-. */
            if (a_partner->HIV_status==UNINFECTED){

                //print_here_string("---",0);
                //print_here_string("remove from serodiscordant because HIV+ partner died",0);

                /* If they only have 1 seropositive partner (which should be the dead person), then remove them from the list susceptible_in_serodiscordant_partnership. */
                if (a_partner->n_HIVpos_partners==1){
                    /* Provided there is more than one HIV- in any serodiscordant partnership, do the following:
                     *  - swap out that person for the last person in the list.
                     *  - reduce the number of people in the list by 1.
                     *  - for the last person in the list, change their idx_serodiscordant index.
                     *  - for the partner, set their idx_serodiscordant index to -1 (as this was their only HIV+ partner). */
                    if ((*n_susceptible_in_serodiscordant_partnership)>1){
                        susceptible_in_serodiscordant_partnership[a_partner->idx_serodiscordant] = susceptible_in_serodiscordant_partnership[n_susceptible_in_serodiscordant_partnership[0] - 1];
                        (*n_susceptible_in_serodiscordant_partnership)--;
                        susceptible_in_serodiscordant_partnership[a_partner->idx_serodiscordant]->idx_serodiscordant = a_partner->idx_serodiscordant;
                        a_partner->idx_serodiscordant = -1;
                    }
                    /* If only one person, do the above apart from swapping (in this case set the pointer to NULL). */
                    else if ((*n_susceptible_in_serodiscordant_partnership)==1){
                        susceptible_in_serodiscordant_partnership[a_partner->idx_serodiscordant] = NULL;
                        (*n_susceptible_in_serodiscordant_partnership)--;           
                        a_partner->idx_serodiscordant = -1;
                    }
                    /* In either case their only seropositive partner has just died: */
                    a_partner->n_HIVpos_partners = 0;
                    if(a_partner->patch_no != dead_person->patch_no)
                        a_partner->n_HIVpos_partners_outside = 0;
                }
                /* Otherwise the partner is still in at least one serodiscordant partnership, so just need to reduce the number of seropositive partners by 1: */
                else{

                    //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership",20);
                    /* Firstly we need to find where in that list is the dead person. */

                    ////////////////////////////////////////////////////////////
                    ////// I AM COMMENTING THIS OUT AS I THINK THIS IS WRONG ////
                    //                  i=0;
                    //                  do{
                    //                      i++;
                    //                  } while ((a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender])!=dead_person);
                    ////////////////////////////////////////////////////////////////
                    ///// SHOULD BE INSTEAD (I THINK) ////
                    i=0;

                    if(dead_person->id == FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH)
                    {
                        printf("Person %li from patch %d has died with HIV status %d \n",dead_person->id,dead_person->patch_no, dead_person->HIV_status);
                        print_individual(dead_person);
                        if (VERBOSE_OUTPUT==1){
                            printf("Looking at HIV positive partner: %li in patch %d, who has %i HIV positive partners \n",a_partner->id,a_partner->patch_no, a_partner->n_HIVpos_partners);
                            print_individual(a_partner);
                        }
                        //printf("a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender] %p\n",a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender]);
                        //printf("dead_person %p\n",dead_person);

                        //printf("a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender]->id %li\n",a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender]->id);
                        //printf("a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender]->patch_no %d\n",a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender]->patch_no);
                        //printf("dead_person->id %li\n",dead_person->id);
                        //printf("dead_person->patch_no %d\n",dead_person->patch_no);

                    }

                    while ((a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender])!=dead_person){
                        //print_here_string("in the while before increment, i = ",i);
                        i++;
                        //print_here_string("in the while after increment, i = ",i);
                    }
                    ////////////////////////////////////////////////////////////////
                    //printf("CHECKME: %li %li\n", a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender]->id,dead_person->id);

                    //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership",21);

                    if(PRINT_DEBUG_DEMOGRAPHICS==1) printf("CHECKME: %li %li\n", a_partner->partner_pairs_HIVpos[i]->ptr[dead_person->gender]->id,dead_person->id);
                    /* Now swap out that partner (as they are dead) */
                    a_partner->partner_pairs_HIVpos[i] = a_partner->partner_pairs_HIVpos[a_partner->n_HIVpos_partners-1];

                    //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership",22);
                    /* Finally reduce number of seropositive partners by 1: */
                    (a_partner->n_HIVpos_partners)--;
                    if(a_partner->patch_no != dead_person->patch_no)
                        a_partner->n_HIVpos_partners_outside--;

                    //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership",23);
                }
            }
        }
    }
    /* If dead person is seronegative and in serodiscordant partnership and there is at least one other person in the same situation, swap the last person with the dead person in this list: */
    /* Provided there is more than one HIV- in any serodiscordant partnership, do the following:
     *  - swap out that person for the last person in the list.
     *  - reduce the number of people in the list by 1.
     *  - for the last person in the list, change their idx_serodiscordant index. */

    /* Otherwise the dead person is seronegative. In this case we only have to worry if the dead person 
     * has seropositive partners, in which case they are in susceptible_in_serodiscordant_partnership[] 
     * which needs updating. */ 
    else if (dead_person->idx_serodiscordant!=-1){

        //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership in the else",0);

        if(dead_person->id==FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH)
        {
            printf("Individual %ld from patch %d is dying- removing him/her from list of susceptibles in serodiscordant partnership\n",dead_person->id,dead_person->patch_no);
            fflush(stdout);
        }

        //print_here_string("remove_dead_person_from_susceptible_in_serodiscordant_partnership in the else",1);
        //print_here_string("---",0);
        //print_here_string("remove from serodiscordant because died whilst HIV-",0);

        if ((*n_susceptible_in_serodiscordant_partnership)>1){
            susceptible_in_serodiscordant_partnership[dead_person->idx_serodiscordant] = susceptible_in_serodiscordant_partnership[n_susceptible_in_serodiscordant_partnership[0] - 1];
            (*n_susceptible_in_serodiscordant_partnership)--;
            susceptible_in_serodiscordant_partnership[dead_person->idx_serodiscordant]->idx_serodiscordant = dead_person->idx_serodiscordant;
            dead_person->idx_serodiscordant = -1;
        }
        /* If only one person, do the above apart from swapping (in this case set the pointer to NULL). */
        else if ((*n_susceptible_in_serodiscordant_partnership)==1){
            susceptible_in_serodiscordant_partnership[dead_person->idx_serodiscordant] = NULL;
            (*n_susceptible_in_serodiscordant_partnership)--;
            dead_person->idx_serodiscordant = -1;
        }
    }
}

/* Function arguments: pointer to the person who died.
 * Function does: removes a dead person from the list of pop_available_partners
 * and n_pop_available_partners structs.
 * Function returns: nothing. */
void remove_dead_person_from_list_available_partners(double time_death, individual *dead_person, population_partners *pop_available_partners, population_size_all_patches *n_pop_available_partners){
    int n, g, ag, r, j, p;

    ag = get_age_group(dead_person->DoB,time_death, AGE_GROUPS, N_AGE);
    r = dead_person->sex_risk;
    g = dead_person->gender;
    p = dead_person->patch_no;

    if(dead_person->id==FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH)
    {
        printf("Individual %ld from patch %d is dying at time %lg - removing available partners. Npartners= %i MAX_Npartners=%i\n",dead_person->id,dead_person->patch_no,time_death,dead_person->n_partners,dead_person->max_n_partners);
        fflush(stdout);
    }


    ///// Debugging:
    //printf("Details to debug: %i %i %i, %li \n",g,ag,r,n_pop_available_partners->pop_size_per_gender_age_risk[g][ag][r] - 1);
    //printf("remove_dead_person_from_list_available_partners: %li\n",n_pop_available_partners->pop_size_per_gender_age_risk[0][6][2]);
    ///fflush(stdout);
    //printf("People: %li",n_pop_available_partners->pop_size_per_gender_age_risk[0][1][0]);
    //for (n=0; n<n_pop_available_partners->pop_size_per_gender_age_risk[0][6][2];n++){
    //  printf("%8li ",pop_available_partners->pop_per_gender_age_risk[0][6][2][n]->id);
    //}
    //printf("\n");
    //fflush(stdout);

    //printf(" %li\n",pop_available_partners->pop_per_gender_age_risk[g][ag][r][n_pop_available_partners->pop_size_per_gender_age_risk[g][ag][r] - 1]->id);
    ////    


    for (n=dead_person->n_partners; n<dead_person->max_n_partners; n++)
    {
        //if(dead_person->id==FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH)
        //printf("Hey %li %li\n",dead_person->idx_available_partner[n-dead_person->n_partners],n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r] - 1);
        if(dead_person->idx_available_partner[n-dead_person->n_partners]<n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r] - 1)
        {
            pop_available_partners->pop_per_patch_gender_age_risk[p][g][ag][r][dead_person->idx_available_partner[n-dead_person->n_partners]] = pop_available_partners->pop_per_patch_gender_age_risk[p][g][ag][r][n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r] - 1]; /* pointing to the last person instead of the current one */
            j = pop_available_partners->pop_per_patch_gender_age_risk[p][g][ag][r][dead_person->idx_available_partner[n-dead_person->n_partners]]->max_n_partners - pop_available_partners->pop_per_patch_gender_age_risk[p][g][ag][r][dead_person->idx_available_partner[n-dead_person->n_partners]]->n_partners - 1;
            while(pop_available_partners->pop_per_patch_gender_age_risk[p][g][ag][r][dead_person->idx_available_partner[n-dead_person->n_partners]]->idx_available_partner[j] != n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r] - 1)
            {
                j--;
            }
            pop_available_partners->pop_per_patch_gender_age_risk[p][g][ag][r][dead_person->idx_available_partner[n-dead_person->n_partners]]->idx_available_partner[j] = dead_person->idx_available_partner[n-dead_person->n_partners]; /* telling the person that has moved that they have. */
            /* switch idx_available partners of dead_person to -1 */
            //// This can probably be ignored but probably easier to identify issues etc. if it is done.
            dead_person->idx_available_partner[n-dead_person->n_partners] = -1;
        }
        n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r]--; /* decreasing the number of available females in that group by 1 */
    }
}

/* Function arguments: pointer to the person who died, the list of available partners (and the numbers of available partners) plus current time t.
 * Function does: removes the partnerships of people who have died (including serodiscordant partnerships), 
 * and add their partners back to "list of available partners".
 * NOTE: This function can be used for either death from natural causes or HIV-related death 
 * (or indeed anything that removes an individual from all partnerships such as permanent migration if this is ever implemented).
 * Function returns: nothing. */        
void remove_dead_persons_partners(individual *dead_person, population_partners *pop_available_partners, 
        population_size_all_patches *n_pop_available_partners, double t){

    int i,j,ag;
    /* All of these pointers will point to existing memory so no calls to malloc needed - they are there to make code readable. */
    long *n_ptr;
    individual *a_partner;
    partnership *a_partnership_ptr;

    if ( (PRINT_DEBUG_DEMOGRAPHICS==1) || (dead_person->id==FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH) ){
        if (dead_person->HIV_status>UNINFECTED)
            printf("Person %li from patch %d is HIV+ with %i partnerships\n",dead_person->id,dead_person->patch_no,dead_person->n_partners);
        else
            printf("Person %li from patch %d is HIV- with %i serodiscordant partnerships and %i partnerships\n",dead_person->id,dead_person->patch_no,dead_person->n_HIVpos_partners,dead_person->n_partners);
        if (dead_person->n_HIVpos_partners>dead_person->n_partners)
            printf("AAG\n");
        fflush(stdout);
    }

    /* For this we have to loop through all partnerships and serodiscordant partnerships. 
     * It is probably quicker to loop through them separately than to go through partnerships and then check if it is serodiscordant. */ 

    if (dead_person->id==FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH) {
        printf("\n------------ Dead individual characteristics:");
        print_individual(dead_person);
    }

    /**********************************/
    /* Loop through all partnerships: */
    /**********************************/
    for (i=0;i<dead_person->n_partners;i++){
        /* Pointer to the partnership: */
        a_partnership_ptr = dead_person->partner_pairs[i];
        /* This is a pointer to the partner: */
        a_partner = a_partnership_ptr->ptr[1-dead_person->gender];

        if (dead_person->id==FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH && VERBOSE_OUTPUT==1) {
            printf("\n------------ Partner %d of dead individual characteristics:",i);
            print_individual(a_partner);
        }

        if (a_partner->id==FOLLOW_INDIVIDUAL && a_partner->patch_no==FOLLOW_PATCH && VERBOSE_OUTPUT==1) {
            printf("\n------------ Individual %li from patch %d has a partner who just died (partner was %li from patch %d)\n",a_partner->id,a_partner->patch_no,dead_person->id,dead_person->patch_no);
            print_individual(a_partner);
        }

        //printf("Removing partner %li of dead person %li HIV %i %i\n",a_partner->id,dead_person->id,dead_person->HIV_status,a_partner->HIV_status);
        if (PRINT_DEBUG_DEMOGRAPHICS)
            printf("Removing partner %li of dead person %li\n",a_partner->id,dead_person->id);
        j=0;

        // For debugging:
        if (a_partner->n_partners <=0 ){
            printf("Error - partner has no partnerships!?\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        while ((j<a_partner->n_partners) && ((a_partner->partner_pairs[j])!=a_partnership_ptr))
            j++;
        // For debugging:
        if (j>=a_partner->n_partners){
            printf("Error - partnership not found between dead person id=%li from patch %d and apparent partner %li from patch %d\n",dead_person->id,dead_person->patch_no,a_partner->id,a_partner->patch_no);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        /* Move the last (ie n_partners-1) partnership to the jth partnership - note if j=n_partners-1 this does nothing but that's OK. */
        a_partner->partner_pairs[j] = a_partner->partner_pairs[a_partner->n_partners-1];
        /* Now reduce partnerships by 1. */
        a_partner->n_partners--;    
        if(a_partner->patch_no != dead_person->patch_no)
        {
            a_partner->n_partners_outside--;
        }
        /* Get the age group of this partner: */
        ag = get_age_group(a_partner->DoB,t, AGE_GROUPS, N_AGE);


        /* This is just a shorthand way to write. I create a pointer to the correct place (ie this is a reference so no need to malloc)
         *  - so that we can change the contents of the original place in the n_pop_available_partners struct. */ 
        n_ptr = &n_pop_available_partners->pop_per_patch[a_partner->patch_no].pop_size_per_gender_age_risk[a_partner->gender][ag][a_partner->sex_risk];
        /* Add this partner to the correct place in the pool of available partners. */
        pop_available_partners->pop_per_patch_gender_age_risk[a_partner->patch_no][a_partner->gender][ag][a_partner->sex_risk][*n_ptr] = a_partner;
        a_partner->idx_available_partner[a_partner->max_n_partners - a_partner->n_partners - 1] = *n_ptr;
        (*n_ptr)++;

    }

}

/* Removes someone from the hiv_pos_progression arrays. Note that normally this only happens for dead people.
 * However the same code is used when someone successfully starts ART, so have renamed to remove "dead_person".
 * The variable "reason" tells us whether this is due to non-AIDS death (reason=1), starting ART normally (reason=2), AIDS death (reason=3), or emergency ART (reason 4). */

void remove_from_hiv_pos_progression(individual *indiv, individual ***hiv_pos_progression, long *n_hiv_pos_progression, long *size_hiv_pos_progression, double t, parameters *param, int reason){
    if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        if (reason==1)
            printf("Removing individual %ld from HIV pos progression due to non-AIDS death at time %6.2f\n",indiv->id,t);
        else if (reason==2)
            printf("Removing individual %ld from HIV pos progression as starting ART/becoming VS (through normal cascade) at time %6.2f\n",indiv->id,t);
        else if (reason==3)
            printf("Removing individual %ld from HIV pos progression due to AIDS death at time %6.2f\n",indiv->id,t);
        else if (reason==4)
            printf("Removing individual %ld from HIV pos progression due to starting emergency ART at time %6.2f\n",indiv->id,t);
        else{
            printf("ERROR: Unknown reason for removing from HIV pos progression array Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fflush(stdout);
    }
    int array_index_for_hiv_event = (int) round((t-param->start_time_hiv)*N_TIME_STEP_PER_YEAR); /* index for current time in this array: hiv_pos_progression, only used for debugging */

    long i = indiv->idx_hiv_pos_progression[0]; /* index for hiv_pos_progression where the next hiv event for this individual is planned for */

    /* Person is removed from HIV positive progression - either due to death or early/VS ART - so no HIV event now scheduled. */
    indiv->next_HIV_event=NOEVENT;

    /* If no current event scheduled then stop: */
    if (i==NOEVENT||i==EVENTAFTERENDSIMUL){
        if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
            if (i==NOEVENT){
                printf("Nothing to remove for %li from hiv_pos_progression %li %li - no event scheduled\n",indiv->id,indiv->idx_hiv_pos_progression[0],indiv->idx_hiv_pos_progression[1]);
                fflush(stdout);
            }
            else if (i==EVENTAFTERENDSIMUL){
                printf("Nothing to remove for %li from hiv_pos_progression %li %li - event was after end of simulation\n",indiv->id,indiv->idx_hiv_pos_progression[0],indiv->idx_hiv_pos_progression[1]);
                fflush(stdout);
            }
        }

        //if (!((indiv->next_cascade_event==CASCADEEVENT_ARTDEATH_NONPOPART || indiv->next_cascade_event==CASCADEEVENT_ARTDEATH_POPART) && (indiv->ART_status==LTART_VS || indiv->ART_status==EARLYART))){
        if (!(indiv->ART_status==LTART_VS || indiv->ART_status==EARLYART || i==EVENTAFTERENDSIMUL)){
            printf("Error - no event to remove in remove_from_hiv_pos_progression() for %li. Exiting  ART status %i next_cascade_event %i %li %li,\n",indiv->id,indiv->ART_status,indiv->next_cascade_event,indiv->idx_cascade_event[0],indiv->idx_cascade_event[1]);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        /* Update the ART_status variable to indicate that this person died from AIDS. */
        if (reason==3)
            indiv->ART_status=ARTDEATH;
        return;
    }

    /* Within simul.c deaths_natural_causes() is called before we carry out HIV events at each timestep, so remove all people
     * who are in the current timestep or later (ie if reason==1). Also there is possibility that ART and HIV progression happen
     * in the same timestep (although as hiv event happens before cascade event this should NEVER happen).
     * For AIDS-death and emergency ART the person should not be removed from the list.
     * Note that we need to check if anyone has an HIV event scheduled which occurred
     * BEFORE the current time.
     * For AIDS death or */
    //if ((i>=array_index_for_hiv_event && reason<=2) || (i>array_index_for_hiv_event && (reason==3||reason==4))){
    if (i>=array_index_for_hiv_event && reason<=2){

        /* FOR DEBUGGING: */
        if (hiv_pos_progression[i][indiv->idx_hiv_pos_progression[1]]->id!=indiv->id){
            printf("ERROR: trying to swap out the wrong person in remove_from_hiv_pos_progression(). Trying to swap %li but in hiv_pos_progression[] the person is %li. Exiting\n",indiv->id,hiv_pos_progression[i][indiv->idx_hiv_pos_progression[1]]->id);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        /* We want to swap out the last person in the array hiv_pos_progression[i] for the indiv. */
        individual *person_to_move;
        person_to_move = hiv_pos_progression[i][n_hiv_pos_progression[i]-1];
        if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
            printf("Swapping %li out from HIV pos progression with %li n_hiv_pos_progression[%li]=%li\n",indiv->id,person_to_move->id,i,n_hiv_pos_progression[i]);
            fflush(stdout);
        }

        if(person_to_move->id==FOLLOW_INDIVIDUAL && person_to_move->patch_no==FOLLOW_PATCH){
            printf("Followed person being swapped %li out from HIV pos progression with %li n_hiv_pos_progression[%li]=%li\n",indiv->id,person_to_move->id,i,n_hiv_pos_progression[i]);
            fflush(stdout);
        }

        /* Now replace the indiv with the person_to_move in hiv_pos_progression: */
        hiv_pos_progression[i][indiv->idx_hiv_pos_progression[1]] = person_to_move;
        /* Update the details of person_to_move (note idx_hiv_pos_progression[0] remains the same): */
        person_to_move->idx_hiv_pos_progression[1] = indiv->idx_hiv_pos_progression[1];
        /* We have removed one person: */
        n_hiv_pos_progression[i]--;

    }

    /* for DEBUGGING: */
    //else if (!(i==array_index_for_hiv_event && (reason==3||reason==4))){
    else if (!(reason==3||reason==4)){
        printf("ERROR: ****Person %ld from patch %d in remove_from_hiv_pos_progression(), trying unsuccessfully to remove HIV event from past %ld %i\n",indiv->id,indiv->patch_no,i,array_index_for_hiv_event);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}



void remove_from_cascade_events(individual *indiv, individual ***cascade_events, long *n_cascade_events, long *size_cascade_events, double t, parameters *param){

    /* Don't need to do anything before start of HIV testing. Use this format so never have problems with index. */
    if (t<param->COUNTRY_HIV_TEST_START)
        return;
    if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        printf("Individual %ld at time %6.2f - removing from cascade events\n",indiv->id,t);
        fflush(stdout);
    }

    int array_index_for_cascade_event =  (int) (round((t - param->COUNTRY_HIV_TEST_START) * N_TIME_STEP_PER_YEAR));
    long i = indiv->idx_cascade_event[0];

    /* If no current event scheduled then stop: */
    if (i==NOEVENT){
        return;
    }


    /* Within simul.c deaths_natural_causes() is called before we carry out cascade events at each timestep, so 
     * remove all people who are in the current timestep or later. Note that we need to check if anyone has 
     * a cascade event scheduled which occurred BEFORE the current time. */
    if (i>=array_index_for_cascade_event){


        /* FOR DEBUGGING: */
        if (cascade_events[i][indiv->idx_cascade_event[1]]->id!=indiv->id){
            printf("ERROR: trying to swap out the wrong person in remove_from_cascade_events(). Trying to swap %li but in cascade_events[] the person is %li. Exiting\n",indiv->id,cascade_events[i][indiv->idx_cascade_event[1]]->id);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        /* We want to swap out the last person in the array cascade_events[i] for the indiv, if there is someone to swap. */
        if (n_cascade_events[i]>0){
            individual *person_to_move;

            person_to_move = cascade_events[i][n_cascade_events[i]-1];

            /* Now replace the indiv with the person_to_move in cascade_events: */

            cascade_events[i][indiv->idx_cascade_event[1]] = person_to_move;
            /* Update the details of person_to_move (note idx_cascade_event[0] remains the same): */
            person_to_move->idx_cascade_event[1] = indiv->idx_cascade_event[1];
            /* We have removed one person: */
            n_cascade_events[i]--;
        }
    }

    /* for DEBUGGING: */
    else{
        printf("ERROR: ****Person %ld from patch %d at t=%f cd4=%i, trying unsuccessfully to remove cascade event from past %ld %i\n",indiv->id,indiv->patch_no,t,indiv->cd4,i,array_index_for_cascade_event);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}


/* Function checks if we need to remove indiv from the list of scheduled VMMC events vmmc_events[] and 
 * removes them if necessary. */
void remove_from_vmmc_events(individual *indiv, individual ***vmmc_events, long *n_vmmc_events, long *size_vmmc_events, double t, parameters *param){

    /* Don't need to do anything if before start of VMMC.  */
    if (t<param->COUNTRY_VMMC_START)
        return;

    if (indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        printf("Individual %li is in remove_from_vmmc_events\n",indiv->id);
    }

    long i = indiv->idx_vmmc_event[0]; 

    /* If not currently scheduled for any VMMC events then return. */
    if (i==NOEVENT)
        return;


    /* FOR DEBUGGING: */
    if (vmmc_events[i][indiv->idx_vmmc_event[1]]->id!=indiv->id){
        printf("ERROR: trying to swap out the wrong person in remove_from_vmmc_events(). Trying to swap %li but in vmmc_events[] the person is %li. Exiting\n",indiv->id,vmmc_events[i][indiv->idx_vmmc_event[1]]->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }


    if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        printf("Individual %ld at time %6.2f - removing from VMMC events\n",indiv->id,t);
        fflush(stdout);
    }

    /* We want to swap out the last person in the array vmmc_events[i] for the indiv. */
    individual *person_to_move;
    person_to_move = vmmc_events[i][n_vmmc_events[i]-1];

    if (person_to_move->id==FOLLOW_INDIVIDUAL && person_to_move->patch_no==FOLLOW_PATCH){
        printf("Individual %li is person_to_move in remove_from_vmmc_events\n",person_to_move->id);
    }

    /* Now replace the indiv with the person_to_move in vmmc_events[]: */ 
    vmmc_events[i][indiv->idx_vmmc_event[1]] = person_to_move;

    /* Update the details of person_to_move (note idx_vmmc_event[0] remains the same): */
    person_to_move->idx_vmmc_event[1] = indiv->idx_vmmc_event[1];

    /* We have removed one person: */
    n_vmmc_events[i]--; 

}


void deaths_natural_causes(double t, patch_struct *patch, int p, 
    all_partnerships *overall_partnerships, file_struct *file_data_store){
   /* Age-specific deaths from natural causes; related calls to remove those from data structures
    
    This function determines the natural death rate (from natural_death_rate()) and then picks the 
    individuals in each age group who actually die.  Then other functions are called to delete
    these individuals from various lists within particular data structures (including sorting out
    partnerships and lists which the individual belonged to). 
    
    
    Arguments
    ---------
    t : double
        Current time in decimal years
    patch : Pointer to a patch_struct object
        Patch object
    p : int
        Patch identifier
    overall_partnerships : pointer to an all_partnerships structure
    file_data_store : pointer to a file_struct structure
    
    */

    int aa, ai, n_death_per_timestep;
    double p_death_per_timestep;
    int i,g, achecktemp;
    
    // Pointer to the person dying (so no need to malloc as pointing at pre-allocated memory) 
    individual *person_dying;

    // Note that we deal with deaths aged MAX_AGE separately - it is a separate array in age_list.
    
    // Loop over genders and age groups
    for(g = 0; g < N_GENDER; g++){
        for(aa = 0; aa < (MAX_AGE - AGE_ADULT); aa++){
            
            // ai is the index of the array age_list->number_per_age_group of the age group of
            // people you want to be dead 
            ai = patch[p].age_list->age_list_by_gender[g]->youngest_age_group_index + aa;
            
            while(ai > (MAX_AGE - AGE_ADULT - 1)){
                ai = ai - (MAX_AGE - AGE_ADULT);
            }
            
            if(PRINT_DEBUG_DEMOGRAPHICS == 1){
                if(patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai] > 0){
                    printf("Number of people[%i] age %i ", ai, aa + AGE_ADULT);
                    printf("gender %i = %li. DoB of first person is = %f\n",
                        g, patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai],
                        patch[p].age_list->age_list_by_gender[g]->age_group[ai][0]->DoB);
                }else{
                    printf("Number of people[%i] age %i gender %i = %li.\n", ai, aa + AGE_ADULT, g,
                        patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai]);
                }
            }
            
            p_death_per_timestep = 
                natural_death_rate(aa + AGE_ADULT, g, patch[p].param, t) * TIME_STEP;
            
            // This was checked against calculations from Excel file (and validated):
            //printf("For age %i at time %6.2f death rate per timestep = %f %f\n", 
            // age+AGE_ADULT,t,p_death_per_timestep,natural_death_rate(age+AGE_ADULT, t, param));
            
            // This is the number of people in age group `a` who die per timestep
            n_death_per_timestep =  gsl_ran_binomial(rng, p_death_per_timestep,
                patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai]);
            
            patch[p].DEBUG_NDEATHS = patch[p].DEBUG_NDEATHS + n_death_per_timestep;
            
            if(PRINT_DEBUG_DEMOGRAPHICS == 1){
                printf("1: Number of people age %i = %li. Number dying = %i\n", 
                    aa + AGE_ADULT, 
                    patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai],
                    n_death_per_timestep);
            }
            if(n_death_per_timestep > 0){
                // 
                gsl_ran_choose(rng, patch[p].new_deaths, n_death_per_timestep,
                    patch[p].death_dummylist, 
                    patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai], 
                    sizeof(long));
                
                // If ageing is switched off then we can accumulate too many individuals in
                // youngest age groups (as they don't die, they just accumulate). 
                
                // For DEBUGGING:
                if(n_death_per_timestep >
                    patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai]){
                    
                    printf("Too many people in age group aa = %i: Number = %li\n", 
                        ai, patch[p].age_list->age_list_by_gender[g]->number_per_age_group[ai]);
                    printf("Is ageing switched off?\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                
                for(i = n_death_per_timestep - 1; i >= 0; i--){
                    
                    person_dying = patch[p].age_list->age_list_by_gender[g]->age_group[ai][(int) patch[p].new_deaths[i]];
                    
                    // Now remove people who have died and to update their partnerships
                    achecktemp = get_age_index(person_dying->DoB, 
                        patch[p].param-> start_time_simul);
                    
                    if(person_dying->id == FOLLOW_INDIVIDUAL && p == FOLLOW_PATCH){
                        find_in_age_list(t, person_dying, patch[p].age_list, patch[p].param);
                    }
                    
                    // For debugging
                    if(person_dying->gender != g){
                        printf("Error - mismatch of gender in deaths_natural_causes(). Exiting\n");
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                    
                    if(ai != achecktemp){
                        printf("AAG %i %i t=%6.4f DoB=%6.4f", ai, achecktemp, t, person_dying->DoB);
                        printf(" age = %6.4f id=%li ", t-person_dying->DoB, person_dying->id);
                        printf("patch=%i param->start_time_simul=%i risk = %i gender = %i youngest_age_index=%i\n",p,patch[p].param-> start_time_simul,person_dying->sex_risk,g, patch[p].age_list->age_list_by_gender[g]->youngest_age_group_index);
                        
                        printf("Wrong id = %li\n",patch[p].age_list->age_list_by_gender[g]->age_group[achecktemp][(int) patch[p].new_deaths[i]]->id);
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                    
                    if(
                        (person_dying->id == FOLLOW_INDIVIDUAL) &&
                        (person_dying->patch_no == FOLLOW_PATCH)
                    ){
                        printf("Killing %li from patch %d ", 
                            person_dying->id, person_dying->patch_no);
                        printf("from natural causes at time %6.2f. Next HIV event was ", t);
                        printf("type= %i with indices %li %li\n\n", person_dying->next_HIV_event,
                            person_dying->idx_hiv_pos_progression[0],
                            person_dying->idx_hiv_pos_progression[1]);
                    }
                    remove_dead_person_from_susceptible_in_serodiscordant_partnership(person_dying,
                        overall_partnerships->susceptible_in_serodiscordant_partnership,
                        overall_partnerships->n_susceptible_in_serodiscordant_partnership);
                    
                    remove_dead_person_from_list_available_partners(t, person_dying,
 overall_partnerships->pop_available_partners,overall_partnerships->n_pop_available_partners);
                    
                    remove_dead_persons_partners(person_dying,
                        overall_partnerships->pop_available_partners,
                        overall_partnerships->n_pop_available_partners, t);
                        
                    if(person_dying->HIV_status > UNINFECTED){
                        // Note the final '1' argument means that the person is dying, 
                        // not starting ART.
                        remove_from_hiv_pos_progression(person_dying, patch[p].hiv_pos_progression,
                            patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression, t,
                            patch[p].param,1);
                    }
                    
                    if(PRINT_DEBUG_DEMOGRAPHICS == 1){
                        printf("Calling deaths_natural_causes() with %i partners\n",
                            person_dying->n_partners);
                    }
                    
                    remove_from_cascade_events(person_dying, patch[p].cascade_events,
                        patch[p].n_cascade_events, patch[p].size_cascade_events, t, 
                        patch[p].param);
                    
                    if(g == MALE){
                        remove_from_vmmc_events(person_dying, patch[p].vmmc_events,
                            patch[p].n_vmmc_events, patch[p].size_vmmc_events, t, patch[p].param);
                    }
                    
                    // Now update popn counts: n_population, n_infected, n_population_stratified
                    update_population_size_death(person_dying, patch[p].n_population,
                        patch[p].n_population_oneyearagegroups, patch[p].n_infected,
                        patch[p].n_population_stratified, aa, patch[p].age_list);

                    // Output time person was seropositive if HIV+ and not on ART.
                    // The final argument is reason for being removed from survival cohort. 1="AIDS
                    // death", 2="death from natural causes", 3="start ART". Note we don't bother
                    // with the end of the simulation for now.
                    if(WRITE_DEBUG_HIV_DURATION_KM == 1){
                        if (person_dying->HIV_status > UNINFECTED){
                            write_hiv_duration_km(person_dying, t, file_data_store, 2);
                        }
                    }
                    
                    // Assign dead person's CD4 count as being DEAD and assign date-of-death (DoD)
                    (person_dying)->cd4 = DEAD;
                    (person_dying)->DoD = t;
                    
                    update_age_list_death(patch[p].age_list, g, ai, 
                        (int) patch[p].new_deaths[i], t, p);
                }
            }
        }
        
        /******************** Now deal with oldest people: ********************/
        p_death_per_timestep = natural_death_rate(MAX_AGE, g, patch[p].param, t) * TIME_STEP;
        
        /* This is the number of people age >=MAX_AGE who die per timestep. */
        n_death_per_timestep =  gsl_ran_binomial(rng, p_death_per_timestep, 
            patch[p].age_list->age_list_by_gender[g]->number_oldest_age_group);

        /* Store number of deaths. */
        patch[p].DEBUG_NDEATHS = patch[p].DEBUG_NDEATHS + n_death_per_timestep;

        if(PRINT_DEBUG_DEMOGRAPHICS == 1){
            printf("2: Number of people age %i+ = %li. Number dying = %i\n", MAX_AGE,
                patch[p].age_list->age_list_by_gender[g]->number_oldest_age_group,
                n_death_per_timestep);
        }
        
        if (n_death_per_timestep>0){

            gsl_ran_choose(rng, patch[p].new_deaths, 
                n_death_per_timestep, patch[p].death_dummylist,
                patch[p].age_list->age_list_by_gender[g]->number_oldest_age_group, sizeof (long));
            
            // FOR DEBUGGING: To check we are pointing at the correct thing.
            if(PRINT_DEBUG_DEMOGRAPHICS == 1){
                printf("DoB of first adult age %i adults to be killed = %f\n",
                    MAX_AGE,
                    (patch[p].age_list->age_list_by_gender[g]->oldest_age_group[(int) patch[p].new_deaths[0]])->DoB);
            }
        }

        // Update their relationships
        for(i = n_death_per_timestep - 1; i >= 0; i--){
            person_dying = patch[p].age_list->age_list_by_gender[g]->oldest_age_group[(int) patch[p].new_deaths[i]];

            //printf("ID = %li Gender = %i %i %i %f\n",
            //(age_list->oldest_age_group[(int) new_deaths[i]])->id,
            // (age_list->oldest_age_group[(int) new_deaths[i]])->gender,
            // MAX_AGE,N_AGE-1,t-(age_list->oldest_age_group[(int) new_deaths[i]])->DoB);
            
            if(PRINT_DEBUG_DEMOGRAPHICS == 1){
                printf("ID = %li Gender = %i %i %i %f\n",
                    person_dying->id, g, MAX_AGE, N_AGE - 1, t - (person_dying->DoB));
            }
            
            /* Remove people who have died and to update their partnerships */
            remove_dead_person_from_susceptible_in_serodiscordant_partnership(person_dying, 
                overall_partnerships->susceptible_in_serodiscordant_partnership,
                overall_partnerships->n_susceptible_in_serodiscordant_partnership);
            
            remove_dead_person_from_list_available_partners(t, person_dying,
                overall_partnerships->pop_available_partners,
                overall_partnerships->n_pop_available_partners);
            
            remove_dead_persons_partners(person_dying,
                overall_partnerships->pop_available_partners,
                overall_partnerships->n_pop_available_partners, t);
            
            if(person_dying->HIV_status > UNINFECTED){
                /* Note the final '1' argument means that the person is dying, not starting ART. */
                remove_from_hiv_pos_progression(person_dying, patch[p].hiv_pos_progression,
                    patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression, t,
                    patch[p].param, 1);
                //patch[p].DEBUG_NHIVDEAD++;
            }
            
            remove_from_cascade_events(person_dying, patch[p].cascade_events,
                patch[p].n_cascade_events, patch[p].size_cascade_events,t, patch[p].param);
            
            
            if(g == MALE){
                remove_from_vmmc_events(person_dying, patch[p].vmmc_events, 
                    patch[p].n_vmmc_events, patch[p].size_vmmc_events, t, patch[p].param);
            }

            // Updates population counts
            update_population_size_death(person_dying, patch[p].n_population,
                patch[p].n_population_oneyearagegroups, patch[p].n_infected,
                patch[p].n_population_stratified, MAX_AGE-AGE_ADULT, patch[p].age_list);
            
            // Output time person was seropositive if HIV+
            // The final argument is reason for being removed from survival cohort. 1="AIDS death",
            // 2="death from natural causes", 3="start ART". Note we don't bother with the end of
            // the simulation for now. 
            if(WRITE_DEBUG_HIV_DURATION_KM == 1){
                if (person_dying->HIV_status > UNINFECTED)
                    write_hiv_duration_km(person_dying, t, file_data_store, 2);
            }
            
            // Assign CD4 count (of DEAD) to the dead person, and assign date-of-death
            (person_dying)->cd4 = DEAD;
            (person_dying)->DoD = t;
            // 
            update_age_list_death(patch[p].age_list, g, MAX_AGE-AGE_ADULT, 
                (int) patch[p].new_deaths[i], t, p);
        }
    }
}


/* Function does: Deals with transition to adulthood of children from child_population at each timestep. 
 * Children are assigned by hivstatus, but other characteristics (gender, risk, etc) assigned by create_new_individual() function.
 * Function arguments: pointers to child_population, the individual population, size of the population, age_list, params. Current time t.
 * Function returns: nothing. */
void make_new_adults(double t, patch_struct *patch, int p, all_partnerships *overall_partnerships){
    int hivstatus;
    if(PRINT_DEBUG_DEMOGRAPHICS == 1){
        //printf("Number of new HIV- (and HIV+) kids = %li %li\n",*(patch[p].child_population[0].transition_to_adult_index_n_child),*(patch[p].child_population[1].transition_to_adult_index_n_child));
        printf("Number of new HIV- (and HIV+) kids = %li %li\n",patch[p].child_population[0].n_child[patch[p].child_population[0].debug_tai],patch[p].child_population[1].n_child[patch[p].child_population[1].debug_tai]);
    }


    /* Store number of new adults for validation - compare with output of print_one_year_age_groups_including_kids() to see if the correct number of new adults are created. */
    //if (WRITE_DEBUG_DEMOGRAPHICS_NBIRTHS_NEWADULTS_DEATHS){
    //printf("patch[p].DEBUG_NNEWADULTS = %li",patch[p].DEBUG_NNEWADULTS);

    //if ((*(patch[p].child_population[0].transition_to_adult_index_n_child)!=patch[p].child_population[0].n_child[patch[p].child_population[0].debug_tai])||(*(patch[p].child_population[1].transition_to_adult_index_n_child)!=patch[p].child_population[1].n_child[patch[p].child_population[1].debug_tai]))
    //  printf("*(patch[p].child_pop) = %li %li  debug_tai = %li %li \n",*(patch[p].child_population[0].transition_to_adult_index_n_child),*(patch[p].child_population[1].transition_to_adult_index_n_child),patch[p].child_population[0].n_child[patch[p].child_population[0].debug_tai],patch[p].child_population[1].n_child[patch[p].child_population[1].debug_tai]);

    //patch[p].DEBUG_NNEWADULTS = patch[p].DEBUG_NNEWADULTS + *(patch[p].child_population[0].transition_to_adult_index_n_child)+*(patch[p].child_population[1].transition_to_adult_index_n_child);
    patch[p].DEBUG_NNEWADULTS = patch[p].DEBUG_NNEWADULTS + patch[p].child_population[0].n_child[patch[p].child_population[0].debug_tai]+patch[p].child_population[1].n_child[patch[p].child_population[1].debug_tai];
    //}

    //printf("Number of new HIV- (and HIV+) kids at time %6.2f in patch %i = %i %i\n",t,p,*(child_population[0].transition_to_adult_index_n_child),*(child_population[1].transition_to_adult_index_n_child));
    ////// Add HIV+ kids:
    /* Add all the kids for this  timestep: */
    for (hivstatus=0;hivstatus<=1;hivstatus++){
        //while (*(patch[p].child_population[hivstatus].transition_to_adult_index_n_child)>0){

        while (patch[p].child_population[hivstatus].n_child[patch[p].child_population[hivstatus].debug_tai]>0){
            /* This adds an individual (HIV-) to individual_population: */
            ////// DEBUG ERROR - this 
            // WRONG version: create_new_individual((individual_population+(n_population->total_pop_size)), t, param);
            /// Right version:

            create_new_individual((patch[p].individual_population+patch[p].id_counter), t, patch[p].param, hivstatus, patch[p].n_infected, patch, p, overall_partnerships);

            if (t>=patch[p].param->COUNTRY_HIV_TEST_START)
                initialize_first_cascade_event_for_new_individual((patch[p].individual_population+patch[p].id_counter), t, patch[p].param, patch[p].cascade_events, patch[p].n_cascade_events, patch[p].size_cascade_events);
            patch[p].id_counter++;

            if (patch[p].id_counter>MAX_POP_SIZE){
                printf("Too many adults in the simulation - exiting. Please increase MAX_POP_SIZE\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }



            /* This updates the n_population variable which counts number of people: */
            update_population_size_new_adult((patch[p].individual_population+patch[p].id_counter-1), patch[p].n_population, patch[p].n_population_oneyearagegroups, patch[p].n_population_stratified);

            update_age_list_new_adult(patch[p].age_list,(patch[p].individual_population+patch[p].id_counter-1));

            //printf("NEWID = %li Number of new kids left = %i, total pop = %li GENDER = %i\n",patch[p].id_counter,patch[p].child_population[hivstatus].n_child[patch[p].child_population[hivstatus].debug_tai],n_population_stratified->total_pop_size,(individual_population+patch[p].id_counter-1)->gender);
            if(PRINT_DEBUG_DEMOGRAPHICS == 1){
                //printf("NEWID = %li Number of new kids left = %li, total pop = %li GENDER = %i\n",patch[p].id_counter,*(patch[p].child_population[hivstatus].transition_to_adult_index_n_child),patch[p].n_population_stratified->total_pop_size,(patch[p].individual_population+patch[p].id_counter-1)->gender);
                printf("NEWID = %li Number of new kids left = %li, total pop = %li GENDER = %i\n",patch[p].id_counter,patch[p].child_population[hivstatus].n_child[patch[p].child_population[hivstatus].debug_tai],patch[p].n_population_stratified->total_pop_size,(patch[p].individual_population+patch[p].id_counter-1)->gender);
            }
            /* Have added a kid, so reduce the number we need to add by 1: */
            patch[p].child_population[hivstatus].n_child[patch[p].child_population[hivstatus].debug_tai]--;
            //(*(patch[p].child_population[hivstatus].transition_to_adult_index_n_child))--;


        }
        /* Note that the while loop set the number of kids in this slot to zero - this slot will now be used to store newborn kids. */
            /* Now we have added all the kids from this timestep, move the pointer to the place in the array for kids to add at the next time step. */
        //if ((patch[p].child_population[hivstatus].transition_to_adult_index_n_child)>&(patch[p].child_population[hivstatus].n_child[0]))
        //  (patch[p].child_population[hivstatus].transition_to_adult_index_n_child)--;
        //else
        //  patch[p].child_population[hivstatus].transition_to_adult_index_n_child = &(patch[p].child_population[hivstatus].n_child[(AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1]);

        if ((patch[p].child_population[hivstatus].debug_tai)>0){
            (patch[p].child_population[hivstatus].debug_tai)--;
        }else{
            patch[p].child_population[hivstatus].debug_tai = (AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1;
        }
    }
}


/* Function does: This is a DUMMY function to add newly born babies to the child_population structures so 
 * that there are always new individuals to reach adulthood as time goes on. The function has stochasticity 
 * so slight variation in number of births per timestep.  
 * NOTE: THIS FUNCTION CURRENTLY FORCES 10% OF CHILDREN TO BE HIV+.
 * Function arguments: pointers to child_population, the size of the population (to calculate the number 
 * of children to be born), params. 
 * Function returns: nothing. */
void add_new_kids(double t, patch_struct *patch, int p){
    int aa,ai;     /* index for age groups. */
    long n_births = 0;
    double age_group_fertility_rate_per_timestep = 0.0;
    /* Here we calculate the average per-woman fertility rate per timestep. Note that we ignore fertility in 65+ year olds! */

    /* Here we work out the interpolation index/fraction for this time (as this is the same for each age group). */
    int y0;
    double f;
    get_unpd_time_indices(t, &y0, &f);
    //printf("f=%lf y0=%i\n",f,y0);
    double childhood_mortality_rate = childhood_mortality(patch[p].param, t);
    //printf("t=%6.4lf\n",t);
    //print_here_string("HEY!!!",1);
    for (aa=(UNPD_FERTILITY_YOUNGEST_AGE-AGE_ADULT); aa<=(UNPD_FERTILITY_OLDEST_AGE-AGE_ADULT); aa++){
        //printf("aa=%i\n",aa);

        fflush(stdout);
        ai = aa + patch[p].age_list->age_list_by_gender[FEMALE]->youngest_age_group_index;
        while (ai>(MAX_AGE-AGE_ADULT-1))
            ai = ai - (MAX_AGE-AGE_ADULT);


        /* Get the fertility rate for this age group: */
        //printf("t=%lf age=%i y0=%i f=%lf peragefert = %lf n_age = %li\n",t,aa,y0,f,per_woman_fertility_rate(aa+AGE_ADULT, param, y0, f),age_list->age_list_by_gender[FEMALE]->number_per_age_group[ai]);

        /* We discount the fertility rate by the childhood mortality rate - so we only include children who will survive to adulthood. */
        age_group_fertility_rate_per_timestep = TIME_STEP*(1.0-childhood_mortality_rate)*per_woman_fertility_rate(aa+AGE_ADULT, patch[p].param, y0, f);      // * age_list->age_list_by_gender[FEMALE]->number_per_age_group[ai];

        //printf("t=%6.4lf age=%i afsr=%8.6lf childhood_mortality_rate=%6.4lf \n",t,aa,per_woman_fertility_rate(aa+AGE_ADULT, param, y0, f),childhood_mortality_rate);
        n_births += gsl_ran_binomial (rng, age_group_fertility_rate_per_timestep, patch[p].age_list->age_list_by_gender[FEMALE]->number_per_age_group[ai]);

    }
    /* Store number of new births for model validation/debugging: */
    patch[p].DEBUG_NBIRTHS = patch[p].DEBUG_NBIRTHS + n_births;

    //printf("t=%6.4lf births=%li\n",t,n_births);
    /* Now normalize to per-capita female population rate per timestep. */ 

//  long npop_check = 0;
//
//  for (aa=AGE_ADULT;aa<MAX_AGE;aa++){
//      ai = aa + age_list->age_list_by_gender[FEMALE]->youngest_age_group_index;
//      while (ai>(MAX_AGE-AGE_ADULT-1))
//          ai= ai - (MAX_AGE-AGE_ADULT);
//      npop_check += age_list->age_list_by_gender[FEMALE]->number_per_age_group[ai];
//  }
//  npop_check += age_list->age_list_by_gender[FEMALE]->number_oldest_age_group;

    //printf("npopf = %ld\n",npop_check);

    //total_population_fertility_rate *= (TIME_STEP /n_population_stratified->total_pop_size_per_gender[FEMALE]);


    //n_births = gsl_ran_binomial (rng, total_population_fertility_rate, n_population_stratified->total_pop_size_per_gender[FEMALE]);



    if (PRINT_DEBUG_DEMOGRAPHICS){
        //if ((patch[p].child_population[0].transition_to_adult_index_n_child)<&(patch[p].child_population[0].n_child[(AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1]))
        if ((patch[p].child_population[0].debug_tai)<((AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1))
            printf("BIRTHSx were: %li %li are: %i %i\n",patch[p].child_population[0].n_child[patch[p].child_population[0].debug_tai+1],patch[p].child_population[1].n_child[patch[p].child_population[1].debug_tai+1],(int) floor(n_births*1.0),(int) floor(n_births*0.0));
        else
            printf("BIRTHSy were: %li %li are: %i %i\n",patch[p].child_population[0].n_child[0],patch[p].child_population[1].n_child[0],(int) floor(n_births*1.0),(int) floor(n_births*0.0));
    }

    // This is a debugging routine for future use - assume that no children are HIV+ at this point.
    //if ((patch[p].child_population[0].transition_to_adult_index_n_child)<(&(patch[p].child_population[0].n_child[(AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1])))
    //  *(patch[p].child_population[0].transition_to_adult_index_n_child+1) = (int) floor(n_births*1.0);
    if ((patch[p].child_population[0].debug_tai) < ((AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1))
        patch[p].child_population[0].n_child[patch[p].child_population[0].debug_tai+1] = (int) floor(n_births*1.0);
    else
        (patch[p].child_population[0].n_child[0]) = (int) floor(n_births*1.0);


    //if ((patch[p].child_population[1].transition_to_adult_index_n_child)<(&(patch[p].child_population[1].n_child[(AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1])))
    //  *(patch[p].child_population[1].transition_to_adult_index_n_child+1) = (int) floor(n_births*0.0);
    if ((patch[p].child_population[1].debug_tai) < ((AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1))
        patch[p].child_population[1].n_child[patch[p].child_population[1].debug_tai+1] = (int) floor(n_births*0.0);
    else
        (patch[p].child_population[1].n_child[0]) = (int) floor(n_births*0.0);

}


/* Function arguments: pointers to pop, age_list, individual_population structure. 
 * Function does: goes through all the lists in age_list (by year) to get all currently alive individuals 
 * aged AGE_ADULT+ and sort them into the correct component of the pop structure.
 * This function can be used to get pop to seed HIV. 
 * Function returns: nothing. */
void make_pop_from_age_list(population *pop, age_list_struct *age_list, individual *individual_population){
    int aa, ai;
    int g, ag, r;
    long i;

    /* This is a (temporary) store for the index of the array pop for each gender x ag x risk group. 
     * It is only locally defined - this is fine as long as we only call this routine a few times. 
     * It is automatically defined - calloc() is not necessary as this is (always) a small array (independent of the 
     * size of the population) unless we stratify by hundreds of extra things. */
    long n_pop_temp[N_GENDER][N_AGE][N_RISK];
    /* Set array to zero: */
    for (g=0; g<N_GENDER; g++)
        for (ag=0; ag<N_AGE; ag++)
            for (r=0; r<N_RISK; r++)
                n_pop_temp[g][ag][r] = 0;

    /* First loop over all age groups by gender and year: */
    for (g=0; g<N_GENDER; g++){
        for (aa=0; aa<(MAX_AGE-AGE_ADULT); aa++){
            ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* a is the index of the two arrays age_list->number_per_age_group and age_list->age_group */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);

            /* Now loop over individuals in each year age group: */
            for (i=0; i<age_list->age_list_by_gender[g]->number_per_age_group[ai]; i++){
                ag = FIND_AGE_GROUPS[aa];
                r = age_list->age_list_by_gender[g]->age_group[ai][i]->sex_risk;
                pop->pop_per_gender_age_risk[g][ag][r][n_pop_temp[g][ag][r]] = age_list->age_list_by_gender[g]->age_group[ai][i];
                n_pop_temp[g][ag][r]++;
            }
        }

        /* Now for oldest individuals: */
        for (i=0; i<age_list->age_list_by_gender[g]->number_oldest_age_group; i++){
            ag = N_AGE-1;   /* Oldest age group. */
            r = age_list->age_list_by_gender[g]->oldest_age_group[i]->sex_risk;
            pop->pop_per_gender_age_risk[g][ag][r][n_pop_temp[g][ag][r]] = age_list->age_list_by_gender[g]->oldest_age_group[i];
            n_pop_temp[g][ag][r]++;
        }
    }
}




/* 
 * Routine which calls a subfunction to delete an individual dying of AIDS (including sorting out 
 * partnerships and lists which the individual belonged to). 
 * Note they are already dying of AIDS so we don't need to fix hiv_pos_progression.
 * Also note this function does not remove the individual from the cascade_events, this needs to be done separately by calling remove_from_cascade_events */
void individual_death_AIDS(age_list_struct *age_list, individual *dead_person, 
        population_size *n_population,  population_size_one_year_age *n_population_oneyearagegroups, population_size_one_year_age *n_infected,
        stratified_population_size *n_population_stratified, double t, parameters *param, 
        individual **susceptible_in_serodiscordant_partnership, 
        long *n_susceptible_in_serodiscordant_partnership, population_partners *pop_available_partners, 
        population_size_all_patches *n_pop_available_partners, individual ***cascade_events, long *n_cascade_events, long *size_cascade_events, patch_struct *patch, int p, file_struct *file_data_store){
    int aa, ai, age_list_index;
    int g = dead_person->gender;

    //print_here_string("individual_death_AIDS line",1466);
    //printf("Reached individual_death_AIDS() for individual %li in patch %i at t=%6.2f. Exiting\n",dead_person->id,dead_person->patch_no,t);
    //printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
    //fflush(stdout);
    //exit(1);

    if (dead_person->id==FOLLOW_INDIVIDUAL && dead_person->patch_no==FOLLOW_PATCH){
        printf("Killing %li from patch %d by HIV at time %6.2f\n",dead_person->id,dead_person->patch_no,t);
        fflush(stdout);
    }
    patch[p].OUTPUT_NDIEDFROMHIV++;

    if (WRITE_DEBUG_HIV_DURATION==1){
        /* Only consider people who are ART-naive for now: */
        //if (dead_person->ART_status==ARTNAIVE||dead_person->ART_status==ARTNEG)
        write_hiv_duration(dead_person, t, file_data_store);
        //      FILE *DEBUGHIVDURFILE;
        //      DEBUGHIVDURFILE = fopen("DEBUGHIVDURATION.csv","a");
        //      fprintf(DEBUGHIVDURFILE,"%6.4lf %d %8.6f %6.4lf %d %d %d\n",t,dead_person->patch_no,dead_person->DEBUGTOTALTIMEHIVPOS,dead_person->t_sc,dead_person->ART_status,dead_person->SPVL_cat,dead_person->cd4);
        //      fclose(DEBUGHIVDURFILE);
    }

    /* The final argument is reason for being removed from survival cohort. 1="AIDS death", 2="death from natural causes", 3="start ART". Note we don't bother with the end of the simulation for now. */
    if (WRITE_DEBUG_HIV_DURATION_KM==1){
    /* Only consider people who are ART-naive for now: */
        //if (dead_person->ART_status==ARTNAIVE||dead_person->ART_status==ARTNEG)
            write_hiv_duration_km(dead_person, t, file_data_store, 1);
    }

    //print_here_string("individual_death_AIDS line",1474);

    /* Note that we deal with deaths aged MAX_AGE separately - it is a separate array in age_list */
    aa = (int) floor(floor(t) - dead_person->DoB) - AGE_ADULT;
    if (aa<(MAX_AGE-AGE_ADULT)){

        //print_here_string("YOUNG",0);

        ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
        while (ai>(MAX_AGE-AGE_ADULT-1))
            ai = ai - (MAX_AGE-AGE_ADULT);

        //print_here_string("YOUNG",1);

        if (PRINT_DEBUG_DEMOGRAPHICS){
            if (age_list->age_list_by_gender[g]->number_per_age_group[ai]>0)
                printf("Number of people[%i] gender %i age %i = %li. DoB of dead person is = %f\n",ai,aa+AGE_ADULT,g,age_list->age_list_by_gender[g]->number_per_age_group[ai],dead_person->DoB);
            else
                printf("Number of people[%i] gender %i age %i = %li.\n",ai,aa+AGE_ADULT,g,age_list->age_list_by_gender[g]->number_per_age_group[ai]);
        }


        //// Remove this individual from individual_population?? - note this isn't easy.

        remove_dead_person_from_susceptible_in_serodiscordant_partnership(dead_person, susceptible_in_serodiscordant_partnership, n_susceptible_in_serodiscordant_partnership);

        remove_dead_person_from_list_available_partners(t, dead_person, pop_available_partners,n_pop_available_partners);

        remove_dead_persons_partners(dead_person, pop_available_partners, n_pop_available_partners, t);

        if (PRINT_DEBUG_DEMOGRAPHICS==1)
            printf("Calling deaths_natural_causes() with %i partners\n",dead_person->n_partners);

        /* Now update the n_population, n_infected and n_population_stratified counts. */
        update_population_size_death(dead_person, n_population, n_population_oneyearagegroups, n_infected, n_population_stratified, aa, age_list); /* Updates population counts. */


        //////// For DEBUGGING:
        dead_person->cd4 = DEAD;
        dead_person->DoD = t;
        //dead_person->DoB = -1;

        age_list_index = 0;
        while ((age_list_index<age_list->age_list_by_gender[g]->number_per_age_group[ai]) && (age_list->age_list_by_gender[g]->age_group[ai][age_list_index]->id!=dead_person->id)){
            age_list_index++;
        }

        if (age_list->age_list_by_gender[g]->age_group[ai][age_list_index]->id!=dead_person->id){
            printf("ERROR: Not sure why - didn't find dead person %li in patch %i. Exiting\n",dead_person->id,dead_person->patch_no);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        update_age_list_death(age_list, g, ai, age_list_index, t, p);


    }
    /******************** Now deal with oldest people: ********************/
    else{
        /* Do we bother to update their partnerships? YES!!! */

        if (PRINT_DEBUG_DEMOGRAPHICS)
            printf("ID = %li Gender = %i %i %i %f\n",dead_person->id,g,MAX_AGE,N_AGE-1,t-dead_person->DoB);

        /* Now call a function to remove people who have died and to update their partnerships. */
        remove_dead_person_from_susceptible_in_serodiscordant_partnership(dead_person, susceptible_in_serodiscordant_partnership, n_susceptible_in_serodiscordant_partnership);
        remove_dead_person_from_list_available_partners(t, dead_person, pop_available_partners,n_pop_available_partners);
        remove_dead_persons_partners(dead_person, pop_available_partners, n_pop_available_partners, t);

        // WRONG CODE: the following line is a call which shouldn't be here as we call this function outside individual_death_AIDS, so this is a repeat call.
        // Have had problems with trying to remove the same dead person twice.
        //remove_from_cascade_events(dead_person, cascade_events, n_cascade_events, size_cascade_events,t, param);

        update_population_size_death(dead_person, n_population, n_population_oneyearagegroups, n_infected, n_population_stratified, MAX_AGE-AGE_ADULT, age_list); /* Updates population counts. */

        dead_person->cd4 = DEAD;
        dead_person->DoD = t;


        age_list_index = 0;
        while ((age_list_index<age_list->age_list_by_gender[g]->number_oldest_age_group) && (age_list->age_list_by_gender[g]->oldest_age_group[age_list_index]->id!=dead_person->id)){
            age_list_index++;
        }
        update_age_list_death(age_list, g, MAX_AGE-AGE_ADULT, age_list_index, t, p);
    }


}


