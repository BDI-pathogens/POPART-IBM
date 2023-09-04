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

/* Contains all functions which initialize the population for the PopART model:
 * get_initial_population_distribution() - calculates the % of individuals in each age x risk group for each gender based 
 * initialize_child_population() - sets up the number of children (by timestep) at the start of the simulation - this will give the number of new adults at each timestep.
 * make_DoB() - gives an individual a DoB at the start of the simulation based on assuming that age is uniformly distributed in each age group and drawing a random number.
 * set_max_n_normal_partners() - gives an individual a maximum number of normal partners at any given time (exclude the number of clients when the individual is sexual worker).
 * set_max_n_clients() - gives a sexual worker a maximum number of clients at any given time.
 * set_population_count_zero() - initializes all the population counts in the "population" structure to be zero.
 * set_population_count_one_year_zero() - initializes all elements of a population_size_one_year_age struct to be zero. These are counts of prevalent and incident cases.
 * set_population_count_stratified_zero() - initializes n_population_stratified struct to be zero.
 * set_population_count_stratified() - based on the maximally stratified (ie gender/age/risk) population counts, generates other (less stratified e.g. by gender only) population counts needed by certain functions.
 * set_up_population() - sets up the "individual_population" - this is an array of param->initial_population_size "individual" structures - at the start of the simulation
 * init_available_partnerships() - sets up the available partnership lists.
 * init_cumulative_counters() - sets the cumulative counter variables in the struct cumulative_outputs to zero. This is called at the start of the simulation. 
 */


/************************************************************************/
/******************************* Includes  ******************************/
/************************************************************************/

#include "constants.h"
#include "init.h"
#include "demographics.h"
#include "utilities.h"
#include "debug.h"

/************************************************************************/
/******************************** functions *****************************/
/************************************************************************/

/* Function does: calculates the proportion of population in each risk group (gender x age x risk) at the start of 
 * the simulation from the initial values in the params_init.txt file and stores it in n_population[][][]. */
void get_initial_population_distribution(population_size *n_population, parameters *param){
    int r, ag;   /* Counters over risk group, age. */
    for (ag = 0; ag < N_AGE; ag++)
        for (r = 0; r < N_RISK; r++){
            n_population->pop_size_per_gender_age_risk[MALE][ag][r] = (long) floor(param->initial_population_size * param->sex_ratio * param->initial_prop_age[ag] * param->initial_prop_gender_risk[MALE][r]);
            n_population->pop_size_per_gender_age_risk[FEMALE][ag][r] = (long) floor(param->initial_population_size * (1.0 - param->sex_ratio) * param->initial_prop_age[ag] * param->initial_prop_gender_risk[FEMALE][r]);
            //printf("Pop size M/F =%li/%li age=%i risk=%i\n", n_population->pop_size_per_gender_age_risk[MALE][ag][r],n_population->pop_size_per_gender_age_risk[FEMALE][ag][r],ag,r);
        }
}

/* 
Function to generate initial population distribution that is limited to be the same size as the 
specified initial population size.  
*/
void get_initial_population_distribution_of_exact_init_size(population_size *n_population, parameters *param){
    int r, ag;   /* Counters over risk group, age. */
    
    int g;
    int genders[2] = {MALE, FEMALE};
    float sex_ratio, randn;
    
    int individual_counter = 0;
    
    generate_initial_individual: ;
    while(individual_counter < param->initial_population_size){
    
    float counter = 0.0;
    
    randn = gsl_rng_uniform_pos (rng);
    
    for (g = 0; g < 2; g++){
        for (ag = 0; ag < N_AGE; ag++){
            for (r = 0; r < N_RISK; r++){
                
                if (genders[g] == MALE){
                    sex_ratio = param->sex_ratio;
                }else{
                    sex_ratio = 1 - param->sex_ratio;
                }
                
                counter += sex_ratio * 
                    param->initial_prop_age[ag] * 
                    param->initial_prop_gender_risk[genders[g]][r];
                printf("counter: %f", counter);
                if (randn < counter){
                    printf("-%d   ", individual_counter);
                    // Add an individual to this strata
                    n_population->pop_size_per_gender_age_risk[genders[g]][ag][r] += 1;
                    individual_counter++;
                    goto generate_initial_individual;
                }
                }
            }
        }
        printf("counter: %f\n", counter);
    }
}


/* Function does: Initializes the child_population structure given fertility rate and number of women by age group.
 * Note: child_population is the number of children who are in each per-timestep age group. Thus at each timestep this
 * specifies the number of children reaching adulthood. */
void initialize_child_population(parameters *param, child_population_struct *child_population, 
        stratified_population_size *n_population_stratified, int country_setting, age_list_struct *age_list){

    int age_dt; /* This is the per-timestep age group of children, with range [0,AGE_ADULT) divided into units of 1 timestep. */
    int ai, aa;     /* index for age groups. */
    int n_births_per_timestep;
    double total_population_fertility_rate = 0.0;

    /* Here we calculate the the per-capita average population-level fertility rate per timestep. 
     * Note that we ignore fertility in 65+ year olds, so the upper bound is N_AGE-2. 
     * The first step is to calculate  the annual cumulative population-level fertility rate
     * (ie per-woman fertility rate per year summed over all women by age group): */
    //for (ag=0; ag< (N_AGE-2); ag++)


    /* Here we work out the interpolation index/fraction for this time
     * (as this is the same for each age group). Note that initialisation normally occurs
     * before the first available UNPD time period UNPD_START so this could in principle
     * be removed but safer to keep it in (and not much overhead). */
    int y0;
    double f;
    get_unpd_time_indices(param->start_time_simul, &y0, &f);

    for (aa=(UNPD_FERTILITY_YOUNGEST_AGE-AGE_ADULT); aa<=(UNPD_FERTILITY_OLDEST_AGE-AGE_ADULT); aa++){
        //printf("Initialising child pop: aa=%i\n",aa);
        //fflush(stdout);
        ai = aa + age_list->age_list_by_gender[FEMALE]->youngest_age_group_index;
        /* Get the fertility rate for this age group: */
        total_population_fertility_rate += per_woman_fertility_rate(aa+AGE_ADULT, param, y0, f) * age_list->age_list_by_gender[FEMALE]->number_per_age_group[ai];
    }
        //total_population_fertility_rate += per_woman_fertility_rate((AGE_GROUPS[ag]+AGE_GROUPS[ag+1])/2.0, param, param->start_time_simul, country_setting) * n_population_stratified->pop_size_per_gender_age[FEMALE][ag];

    /* Now normalize to per-capita female population rate per timestep - so multiply by timestep and divide by number of women: 
     * NB: add 0.0 to ensure always float (rather than integer) division. */ 
    total_population_fertility_rate *= (TIME_STEP/(n_population_stratified->total_pop_size_per_gender[FEMALE]+0.0));

    /* Draw number of children in each age_dt group using this fertility rate: */
    for (age_dt=0;age_dt<((AGE_ADULT+1)*N_TIME_STEP_PER_YEAR);age_dt++){ /* for each unit of age */

        n_births_per_timestep = gsl_ran_binomial (rng, total_population_fertility_rate, n_population_stratified->total_pop_size_per_gender[FEMALE]);
        // Here could draw by age group and then choose the mothers at random so that MTCT can be modelled

        /* All children are assumed HIV- at this point as HIV epidemic has not started. */
        child_population[0].n_child[age_dt] = n_births_per_timestep;
        child_population[1].n_child[age_dt] = 0;
    }
    /* These indices points to the oldest age group in each child_population[i] group: */
    child_population[0].transition_to_adult_index_n_child = ((child_population[0].n_child)+(AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1);
    child_population[1].transition_to_adult_index_n_child = ((child_population[1].n_child)+(AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1);
    child_population[0].debug_tai = (AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1;
    child_population[1].debug_tai = (AGE_ADULT+1)*N_TIME_STEP_PER_YEAR-1;
}

///////////// FIX THIS:
///////////// This is a placeholder - set max_n_partners correctly later
/* Function arguments: risk group of individual
 * Function does: This is a placeholder - should return a number of normal partners based on some data.
 * Function returns: The maximum number of normal partners that an individual can have in that risk group (exclude the number of clients of sexual worker)
 * (may also depend on age/gender?). */
int set_max_n_normal_partners(int g, int ag, int r, int sexual_worker_status, parameters *param){
    // assume the max normal partnership for a sexual worker is 1
    if (sexual_worker_status == SEXUAL_WORKER) {
        return 1;
    }
    else {
        return param->max_n_part_noage[r];
    }
}


/* Function arguments: sexual worker status of individual
 * Function does: return a number of clients based on power law distribution.
 * Function returns: The maximum number of clients that a sexual worker can have */
int set_max_n_clients(int sexual_worker_status) {
    double prob;
    int max_n_client = 101;

    // non-sexual worker has no client
    if (sexual_worker_status == NON_SEXUAL_WORKER) {
        return 0;
    }
    else {
        while (max_n_client > MAX_N_CLIENT) {
            prob = gsl_rng_uniform(rng);
            max_n_client = round(MIN_N_CLIENT * pow(1 - prob, -1 / (SCALING_POWER_LAW - 1)));
        }
        return max_n_client;
    }
}


/* Function arguments: age group of an individual, pointer to "parameter" structure, current time t, 
 * a pointer to age_year which is the age_group that the individual belongs to (we pass as a pointer so 
 * we can change it).
 * Function does: assumes that ages are uniformly distributed between the lower and upper limits of that 
 * age group, assigns an age randomly and then calculates the DoB. The variable age_group is also modified. 
 * Function returns: The DoB of that person. The variable age_group is also modified. */
double make_DoB(int ag, double t, int *aa_ptr){
    double age;
    /* AGE_GROUP[] gives the lower bound of each age group, so need to adapt if in last age group. 
     * Note that we initialise the population to be 14+ and then add in 13 year-olds over the course of the first year. 
     * Use function gsl_rng_uniform_pos as this samples on (0,1) and ensures we never have a DoB which is exactly
     * at the start of any year. */ 

    if (ag==0)
        age = (AGE_GROUPS[0]+1) + (AGE_GROUPS[1]-(AGE_GROUPS[0]+1))*gsl_rng_uniform_pos (rng);
    else if (ag<N_AGE-1)
        age = AGE_GROUPS[ag] + (AGE_GROUPS[ag+1]-AGE_GROUPS[ag])*gsl_rng_uniform_pos (rng);
    else
        age = AGE_GROUPS[ag] + (MAX_AGE-AGE_GROUPS[ag])*gsl_rng_uniform_pos (rng);  

    *aa_ptr = (int) (age-AGE_ADULT);

    /* This gives the DoB: */
    return t - age;
}


/* Function arguments: pointer to the "population" structure which contains information about the 
 * number of people in each risk, age and gender group and overall.
 * Function does: Initializes all of these as zero. 
 * Function is used to set n_pop_available_sexual_workers, n_pop_available_normal_partners to zero. */
void set_population_count_zero(population_size *n_pop){
    int g,ag,r;
    for (ag=0; ag<N_AGE; ag++){
        for (r=0; r<N_RISK; r++){
            for(g=0 ; g<N_GENDER ; g++){
                n_pop->pop_size_per_gender_age_risk[g][ag][r] = 0;
            }
        }
    }
}

/* Function does: Initializes n_population structure (which is number of people by gender x one-year age groups x risk group
 *  to zero and sets youngest_age_group_index as 0. */
void set_population_count_one_year_zero(population_size_one_year_age *n_population){
    int g,aa,r;
    n_population->youngest_age_group_index = 0;


    for (r=0; r<N_RISK; r++){
        for(g=0 ; g<N_GENDER ; g++){
            for (aa=0; aa<(MAX_AGE-AGE_ADULT); aa++){
                n_population->pop_size_per_gender_age1_risk[g][aa][r] = 0;
            }
            n_population->pop_size_oldest_age_group_gender_risk[g][r] = 0;
        }
    }
}


void set_population_count_stratified_zero(stratified_population_size *n_population_stratified){
    int g,ag,r;
    for (ag=0; ag<N_AGE; ag++){
        //n_population_stratified->pop_size_per_age[ag] = 0;
        for(g=0 ; g<N_GENDER ; g++){
            n_population_stratified->pop_size_per_gender_age[g][ag] = 0;
        }
        /*for (r=0; r<N_RISK; r++){
            n_population_stratified->pop_size_per_age_risk[ag][r] = 0;
        }*/
    }

    for (r=0; r<N_RISK; r++){
        //n_population_stratified->pop_size_per_risk[r] = 0;
        for(g=0 ; g<N_GENDER ; g++){
            n_population_stratified->pop_size_per_gender_risk[g][r] = 0;
            n_population_stratified->prop_pop_per_gender_risk[g][r] = 0.0;
        }
    }

    n_population_stratified->total_pop_size = 0;
    for(g=0 ; g<N_GENDER ; g++){
        n_population_stratified->total_pop_size_per_gender[g] = 0;
    }
}

void set_population_count_stratified(stratified_population_size *n_population_stratified, 
        population_size *pop){
    int g,ag,r;

    set_population_count_stratified_zero(n_population_stratified);           /* Sets to zero each element in this structure, which stores number of people in the population. */

    for (r=0; r<N_RISK; r++){
        for (ag=0; ag<N_AGE; ag++){
            for(g=0 ; g<N_GENDER ; g++)
            {
                n_population_stratified->pop_size_per_gender_age[g][ag] += pop->pop_size_per_gender_age_risk[g][ag][r];
                n_population_stratified->pop_size_per_gender_risk[g][r] += pop->pop_size_per_gender_age_risk[g][ag][r];
                n_population_stratified->total_pop_size_per_gender[g] += pop->pop_size_per_gender_age_risk[g][ag][r];
                //n_population_stratified->pop_size_per_age[ag] += pop->pop_size_per_gender_age_risk[g][ag][r];
                //n_population_stratified->pop_size_per_age_risk[ag][r] += pop->pop_size_per_gender_age_risk[g][ag][r];
                //n_population_stratified->pop_size_per_risk[r] += pop->pop_size_per_gender_age_risk[g][ag][r];
                n_population_stratified->total_pop_size += pop->pop_size_per_gender_age_risk[g][ag][r];
            }
        }
    }

    for (r=0; r<N_RISK; r++){
        for(g=0 ; g<N_GENDER ; g++)
        {
            /* Multiply denominator by 1.0 to make it a float (C will do integer division if both numerator and denominator are int). */
            n_population_stratified->prop_pop_per_gender_risk[g][r] = n_population_stratified->pop_size_per_gender_risk[g][r]/(1.0*n_population_stratified->total_pop_size_per_gender[g]);

        }
    }
}



/* Function arguments: pointer to the "population" structure which contains information about the number of 
 * people in each risk, age and gender group and overall; pointer to the "individual_population" array of 
 * individuals, pointer to the age-list structure; pointer to the "parameter" structure.
 * Function does: Assigns values to each "individual" structure within the "cluster population" - ie gives 
 * them gender, DoB, riskiness, etc.
 * NOTE: It assigns "-1" to variables such as SPVL and CD4 which will be initialized if the individual 
 * gets HIV, and it sets the numbers of current partners to zero as we initialize partnerships later on. */
void set_up_population(int p, patch_struct *patch, population *pop){

    int i;
    int g, r, ag;
    int aa, a_unpd;
    long i_x;

    /* Here we work out the interpolation index/fraction for this time
     * (as this is the same for each age group). Note that initialisation normally occurs
     * before the first available UNPD time period UNPD_START so this could in principle
     * be removed but safer to keep it in (and not much overhead). */
    int y0;
    double f;
    get_unpd_time_indices(patch[p].param->start_time_simul, &y0, &f);

    /* For debugging : */
    patch[p].DEBUG_NBIRTHS = 0;
    patch[p].DEBUG_NNEWADULTS = 0;
    patch[p].DEBUG_NDEATHS = 0;

    patch[p].DEBUG_NHIVPOS = 0;
    patch[p].DEBUG_NHIVPOSLASTYR = 0;
    /* Stores number of people who ever died from HIV-related illness. Can get e.g. yearly total by taking difference between two years*/
    patch[p].OUTPUT_NDIEDFROMHIV = 0;


    /* For each group we create a template person which has all characteristics in common with everyone in that group.
     * Specific individual characteristics are then assigned to each person after the template is used first to speed things up. */
    individual person_template;    

    int n_men_circumcised_as_children; /* For each group this is the number of men who are circumcised as children: */

    //// Not currently used as in core model children are created without specific mothers. However this may be needed for MTCT.
    double prop_pregnant_women;        /* For each group this is the proportion of women who are pregnant at the start of the simulation: */
    int n_pregnant_women;              /* For each group this is the number of women who are pregnant at the start of the simulation: */

    /* Work out how many people in each sub-population (age x risk x gender). */
    //get_initial_population_distribution_of_exact_init_size(patch[p].n_population,patch[p].param);
    get_initial_population_distribution(patch[p].n_population,patch[p].param);

    set_population_count_stratified(patch[p].n_population_stratified, patch[p].n_population);

    set_population_count_one_year_zero(patch[p].n_infected);
    set_population_count_one_year_zero(patch[p].n_art);
    set_population_count_one_year_zero(patch[p].n_virallysuppressed);
    set_population_count_one_year_zero(patch[p].n_cabo);
    set_population_count_one_year_zero(patch[p].n_newly_infected);
    set_population_count_one_year_zero(patch[p].n_infected_cumulative);


    /* age_list->age_list_by_gender[g]->youngest_age_group_index is the index of the youngest age group array in age_list->age_list_by_gender[g]->age_groups[].
     * As the population ages this index moves (so we don't have to move all the arrays). By default we start it at zero. */
    for (g=0;g<N_GENDER;g++)
        patch[p].age_list->age_list_by_gender[g]->youngest_age_group_index = 0;

    /* Similarly for n_population_oneyearagegroups. */
    patch[p].n_population_oneyearagegroups->youngest_age_group_index = 0;


    patch[p].id_counter = 0;     /* Initialize this variable in struct patch. */

    POPART_SAMPLING_FRAME_ESTABLISHED = 0; /* Initialize global variable to check that we have not yet done PopART sampling (to prevent errors where popart happens before the sampling is done). */

    /* This is a local counter over gender x age group x risk: */
    int id_counter_per_gender_age_risk[N_GENDER][N_AGE][N_RISK];
    for(g=0 ; g<N_GENDER ; g++){
        for(ag=0 ; ag<N_AGE ; ag++){
            for(r=0 ; r<N_RISK ; r++){
                id_counter_per_gender_age_risk[g][ag][r] = 0;
            }
        }
    }

    /* Create a template individual structure to copy. These are the parts which are always in common: */
    person_template.patch_no = p;
    person_template.HIV_status = UNINFECTED;
    person_template.ART_status = ARTNEG;
    person_template.t_sc = -1;                   /* Initialize at dummy value. */
    person_template.t_vmmc = -1;                   /* Initialize at dummy value. */
    person_template.SPVL_num_G = 0;              /* Initialize at dummy value. */
    person_template.SPVL_num_E = 0;              /* Initialize at dummy value. */
    person_template.SPVL_infector = 0;           /* Initialize at dummy value. */
    person_template.cd4 = CD4_UNINFECTED;        /* Initialize at dummy value. */
    person_template.SPVL_cat = -1;               /* Initialize at dummy value. */
    person_template.DEBUGTOTALTIMEHIVPOS = 0;
    person_template.next_HIV_event = NOEVENT;         /* Initialize at dummy value. */
    person_template.next_cascade_event = NOEVENT;        /* Initialize at dummy value. */
    person_template.time_last_hiv_test = NEVERHIVTESTED; /* No HIV testing at start of simulation as no HIV. */
    person_template.idx_hiv_pos_progression[0] = -1;     /* Initialize at dummy value. */
    person_template.idx_hiv_pos_progression[1] = -1;     /* Initialize at dummy value. */
    person_template.debug_last_hiv_event_index = -1;     /* Initialize at dummy value. */
    person_template.idx_cascade_event[0] = -1;           /* Initialize at dummy value. */
    person_template.idx_cascade_event[1] = -1;           /* Initialize at dummy value. */
    person_template.debug_last_cascade_event_index = -1;     /* Initialize at dummy value. */
    person_template.time_to_delivery = -1;                /* Not pregnant - allow a % of women to be pregnant later. */
    /* Set up partnerships later on. */
    person_template.n_clients = 0;
    person_template.n_normal_partners = 0;
    person_template.n_partners = 0;
    person_template.n_HIVpos_partners = 0;
    person_template.n_HIVpos_partners_outside = 0;
    person_template.n_partners_outside = 0;
    person_template.n_lifetime_partners = 0;
    person_template.n_lifetimeminusoneyear_partners = 0;
    person_template.n_lifetime_partners_outside = 0;
    person_template.n_lifetimeminusoneyear_partners_outside = 0;
    person_template.n_partnersminusoneyear = 0;
    person_template.idx_serodiscordant = -1;              /* Not in a serodiscordant partnership */

    person_template.circ = UNCIRC;                       /* Not circumcised - allow a % of men to be circumcised later. */
    person_template.idx_vmmc_event[0] = -1;         /* Initialize at dummy value. */
    person_template.idx_vmmc_event[1] = -1;
    person_template.debug_last_vmmc_event_index = -1;     /* Initialize at dummy value. */
    for (i=person_template.n_clients; i<MAX_N_CLIENT; i++){
        person_template.idx_available_sexual_worker[i] = -1; /* Not yet in the list of available sexual workers */
    }
    for (i=person_template.n_normal_partners; i<MAX_NORMAL_PARTNERSHIPS; i++){
        person_template.idx_available_normal_partner[i] = -1; /* Not yet in the list of available normal partners */
    }

    person_template.PANGEA_t_prev_cd4stage = -1.0;
    person_template.PANGEA_t_next_cd4stage = -1.0;
    person_template.PANGEA_cd4atdiagnosis = -1.0;
    person_template.PANGEA_cd4atfirstART = -1.0;
    person_template.PANGEA_t_diag = -1.0;
    person_template.PANGEA_date_firstARTstart = -1.0;
    person_template.PANGEA_date_startfirstVLsuppression = -1.0;
    person_template.PANGEA_date_endfirstVLsuppression = -1.0;

    person_template.VISITEDBYCHIPS_TO_INIT_ART = FALSE;
    person_template.VISITED_BY_CHIPS_THISROUND = FALSE;
    person_template.NCHIPSVISITS = 0;
    person_template.PC_cohort_index = -1; /* Not in PC cohort (for now). */

    /* Set this to zero here as we will sum over risk groups r later. */
    for (g = 0; g < N_GENDER; g++){
        for (aa = 0; aa < (MAX_AGE-AGE_ADULT); aa++)
            patch[p].age_list->age_list_by_gender[g]->number_per_age_group[aa] = 0;
        patch[p].age_list->age_list_by_gender[g]->number_oldest_age_group = 0;
    }

    for (r = 0; r < N_RISK; r++){
        for (g = 0; g < N_GENDER; g++){
            for (aa = 0; aa < (MAX_AGE-AGE_ADULT); aa++)
                patch[p].n_population_oneyearagegroups->pop_size_per_gender_age1_risk[g][aa][r] = 0;
            patch[p].n_population_oneyearagegroups->pop_size_oldest_age_group_gender_risk[g][r] = 0;
        }
    }

    double childhood_mortality_rate = childhood_mortality(patch[p].param, patch[p].param->start_time_simul);

    /* Now loop over gender x age x risk: */
    for (ag = 0; ag < N_AGE; ag++){
        /* For speed we draw the number of women in each age group who are pregnant. 
         * First term is to correct the annual rate of pregnancy for the fact that it only lasts 9 months 
         * - so only 3/4 of the annual fraction are pregnant at the start of the simulation. */
        prop_pregnant_women = 0;

        for (a_unpd=AGE_GROUPS[ag]; a_unpd<AGE_GROUPS[ag+1]; a_unpd++){
            if ((a_unpd>=UNPD_FERTILITY_YOUNGEST_AGE) && (a_unpd<=UNPD_FERTILITY_OLDEST_AGE)){
                //print_here_string("Fertility",a_unpd);
                /* We discount the fertility rate by the childhood mortality rate - so we only include children who will survive to adulthood. */
                prop_pregnant_women += (NSTEPS_GESTATION_TIME*TIME_STEP) * (1.0-childhood_mortality_rate)*per_woman_fertility_rate(a_unpd, patch[p].param, y0, f);
            }
        }

        /* We average over all age groups. Note that this isn't exactly the right distribution (generally more people in younger years than older) but it is close enough as we give the model some time to burn in demographics etc.
         * Note that we exclude the last age groups to avoid issues with indices - assumes the model ALWAYS models age groups beyond 49 years old. */
        if (ag<N_AGE-1)
            prop_pregnant_women = prop_pregnant_women/(a_unpd<(AGE_GROUPS[ag+1] - AGE_GROUPS[ag] - 1));
        else{
            if (prop_pregnant_women>0){ /* Just in case the model doesn't include enough age groups print an error message. */
                printf("Issue with normalising fertility rate in set_up_population(). Exiting\n");
                fflush(stdout);
                exit(1);
            }
        }
//      if (ag<N_AGE-1)
//          prop_pregnant_women = (NSTEPS_GESTATION_TIME*TIME_STEP) * per_woman_fertility_rate((AGE_GROUPS[ag]+AGE_GROUPS[ag+1])/2.0, param, param->start_time_simul, country_setting);
//      else
//          prop_pregnant_women = (NSTEPS_GESTATION_TIME*TIME_STEP) * per_woman_fertility_rate((AGE_GROUPS[ag]+MAX_AGE)/2.0, param, param->start_time_simul, country_setting);


        for (r = 0; r < N_RISK; r++){
            /* Calculate number of men circumcised: (we do it here so we don't need to calculate twice for M/F). */
            n_men_circumcised_as_children = (int) floor((patch[p].param->p_child_circ)*patch[p].n_population->pop_size_per_gender_age_risk[MALE][ag][r]);
            n_pregnant_women = (int) floor(prop_pregnant_women * patch[p].n_population->pop_size_per_gender_age_risk[FEMALE][ag][r]);
            for (g = 0; g < N_GENDER; g++){

                /* Set up parts of template which are specific to this age/sex/risk group: */
                person_template.gender = g;
                person_template.sex_risk = r;  /* 0, 1, 2 for LOW, MEDIUM, HIGH risk */
                person_template.sexual_worker_status = NON_SEXUAL_WORKER; 

                for (i_x = 0; i_x < patch[p].n_population->pop_size_per_gender_age_risk[g][ag][r]; i_x++){

                    /* Copy the template onto the person. Note that "id_counter" is a global variable which is 
                     * initialized at zero. Here we use it as both a counter and as the id number of the people 
                     * we are initializing - later on it is just used for the id number of new people. */
                    patch[p].individual_population[patch[p].id_counter] = person_template;

                    /* Now update the unique parts of this individual (ie the bits which aren't in common): */
                    patch[p].individual_population[patch[p].id_counter].id = patch[p].id_counter;

                    if (g==MALE){
                        /* Circumcise the first n_men_circumcised_as_children men in the group: */
                        if (n_men_circumcised_as_children>0){
                            patch[p].individual_population[patch[p].id_counter].circ = TRADITIONAL_MC;
                            n_men_circumcised_as_children--;
                        }
                    }
                    else{
                        /* First n_pregnant_women are pregnant in the group. */
                        if (n_pregnant_women>0){
                            // Again, in core model not using time_to_delivery although keep code to allow MTCT in extended versions of model.
                            patch[p].individual_population[patch[p].id_counter].time_to_delivery = gsl_rng_uniform_int (rng, NSTEPS_GESTATION_TIME);
                            n_pregnant_women--;
                        }
                    }

                    /* We can't assign ages increasing in the group because we already assign circumcision to the first k people in 
                     * the group. So use make_DoB() to assign uniform DoB consistent with age group ag: */                  

                    patch[p].individual_population[patch[p].id_counter].DoB = make_DoB(ag, patch[p].param->start_time_simul, &aa);
                    // Use code below to test different ways of getting age_group index:
                    //if (aa!=get_age_index(individual_population[patch[p].id_counter].DoB, param->start_time_simul)||aa!=get_age_indexv2(individual_population[patch[p].id_counter].DoB, patch[p].param->start_time_simul, patch[p].age_list->age_list_by_gender[g]->youngest_age_group_index))
                    //  printf("ERROR:::In make_DoB aa = %i, calc = %i v2 = %i\n",aa,get_age_index(individual_population[patch[p].id_counter].DoB, param->start_time_simul),get_age_indexv2(individual_population[patch[p].id_counter].DoB, param->start_time_simul, patch[p].age_list->age_list_by_gender[g]->youngest_age_group_index));
                    if (patch[p].individual_population[patch[p].id_counter].id==FOLLOW_INDIVIDUAL && p==FOLLOW_PATCH)
                        printf("In make_DoB DoB=%f start_time = %d aa = %i, calc = %i v2 = %i\n",patch[p].individual_population[patch[p].id_counter].DoB,patch[p].param->start_time_simul,aa,get_age_index(patch[p].individual_population[patch[p].id_counter].DoB, patch[p].param->start_time_simul),get_age_indexv2(patch[p].individual_population[patch[p].id_counter].DoB, patch[p].param->start_time_simul, patch[p].age_list->age_list_by_gender[g]->youngest_age_group_index));

                    patch[p].individual_population[patch[p].id_counter].DoD = -1;
                    if (patch[p].id_counter==FOLLOW_INDIVIDUAL && p==FOLLOW_PATCH)
                        printf("Making individual %li in patch %d at start of simulation with gender %d, age group ag=%i aa=%i DoB = %lf\n",patch[p].id_counter,p,g,ag,aa,patch[p].individual_population[patch[p].id_counter].DoB);

                    /* Now add this person to the correct age list, and then increment the number of people in that list by 1: */
                    if (aa<MAX_AGE-AGE_ADULT){
                        patch[p].age_list->age_list_by_gender[g]->age_group[aa][patch[p].age_list->age_list_by_gender[g]->number_per_age_group[aa]] = &patch[p].individual_population[patch[p].id_counter];
                        patch[p].age_list->age_list_by_gender[g]->number_per_age_group[aa] += 1;
                        if (patch[p].id_counter==FOLLOW_INDIVIDUAL && p==FOLLOW_PATCH){
                            printf("Individual pop id = %li %li %i\n",patch[p].individual_population[patch[p].id_counter].id,patch[p].age_list->age_list_by_gender[g]->number_per_age_group[aa],aa);
                            fflush(stdout);

                        }

                        patch[p].n_population_oneyearagegroups->pop_size_per_gender_age1_risk[g][aa][r] += 1;
                    }
                    else{
                        patch[p].age_list->age_list_by_gender[g]->oldest_age_group[patch[p].age_list->age_list_by_gender[g]->number_oldest_age_group] = &patch[p].individual_population[patch[p].id_counter];
                        patch[p].age_list->age_list_by_gender[g]->number_oldest_age_group += 1;

                        patch[p].n_population_oneyearagegroups->pop_size_oldest_age_group_gender_risk[g][r] += 1;
                    }
                
                    /* Draw max number of clients of a sexual worker */
                    patch[p].individual_population[patch[p].id_counter].max_n_clients = set_max_n_clients(person_template.sexual_worker_status);
                    /* Draw max number of normal partners of an individual */
                    // Note: set_max_n_normal_partners() is currently a placeholder function.
                    patch[p].individual_population[patch[p].id_counter].max_n_normal_partners = set_max_n_normal_partners(g,ag,r,person_template.sexual_worker_status,patch[p].param);
                    /* Draw max number of sexual partners based on gender, age, risk group: it includes the max number of clients a sexual worker can have */
                    patch[p].individual_population[patch[p].id_counter].max_n_partners = patch[p].individual_population[patch[p].id_counter].max_n_clients + patch[p].individual_population[patch[p].id_counter].max_n_normal_partners;

                    /* This is a list of pointers to all individuals of a given g, ag, r: */
                    (pop->pop_per_gender_age_risk[g][ag][r])[id_counter_per_gender_age_risk[g][ag][r]] = &patch[p].individual_population[patch[p].id_counter];


                    /*if(individual_population[patch[p].id_counter].id==5486)
                    {
                        print_here_string("000000000000000000000000000000000000000000000000",0);
                        print_here_string("Created individual in initial population with idx ",individual_population[patch[p].id_counter].id);
                        print_here_string("Patch ",individual_population[patch[p].id_counter].patch_no);
                        print_here_string("Gender ",individual_population[patch[p].id_counter].gender);
                        print_here_string("Circ ",individual_population[patch[p].id_counter].circ);
                        print_here_string("Risk ",individual_population[patch[p].id_counter].sex_risk);
                        print_here_string("Yearob ",(int) individual_population[patch[p].id_counter].DoB);
                        print_here_string("Max_n_partners ",individual_population[patch[p].id_counter].max_n_partners);
                        print_here_string("Current_n_partners ",individual_population[patch[p].id_counter].n_partners);
                        print_here_string("111111111111111111111111111111111111111111111111",1);
                    }*/

                    // For debugging:
                    if(patch[p].individual_population[patch[p].id_counter].id==FOLLOW_INDIVIDUAL && p==FOLLOW_PATCH){
                        printf("Creation of adult in initial population with id %ld in patch %d\n",patch[p].individual_population[patch[p].id_counter].id,patch[p].individual_population[patch[p].id_counter].patch_no);
                        fflush(stdout);
                    }
                    //printf("A:patch[%i].id_counter = %li\n",p,patch[p].id_counter);
                    /* In check_if_parameters_plausible() we check that initial_population_size is < MAX_POP_SIZE so don't need to check that here. */
                    patch[p].id_counter++;
                    //printf("B:patch[%i].id_counter = %li\n",p,patch[p].id_counter);
                    id_counter_per_gender_age_risk[g][ag][r]++;


                }
            }
        }
    }
    count_population_size_three_ways(patch, p, patch[p].param->start_time_simul);
    initialize_child_population(patch[p].param,patch[p].child_population,patch[p].n_population_stratified, patch[p].country_setting, patch[p].age_list);
}


/* This function fills in (*pop_available_sexual_workers) and (*n_pop_available_sexual_workers), (*pop_available_normal_partners) and (*n_pop_available_normal_partners) using (*pop) and (*n_population) 
 * it counts all the "free partnerships" and points to them in pop_available_sexual_workers, pop_available_normal_partners, with associated numbers of 
 * free partnerships stored in n_pop_available_sexual_workers, n_pop_available_normal_partners */
void init_available_partnerships(int p, patch_struct *patch, all_partnerships *overall_partnerships,population *pop){

    int g,ag,r;   /* Indices over gender, age group and risk group and patch. */
    long i;       /* Index to loop over number of people. */
    int j;        /* Index to loop over partnerships. */

    set_population_count_zero(&overall_partnerships->n_pop_available_sexual_workers->pop_per_patch[p]);
    set_population_count_zero(&overall_partnerships->n_pop_available_normal_partners->pop_per_patch[p]);
    for (g=0 ; g<N_GENDER ; g++){
        for (ag=0 ; ag<N_AGE ; ag++){
            for(r=0 ; r<N_RISK ; r++){
                /* Loop over each individual in this gender/age/risk category: */
                for(i=0 ; i<patch[p].n_population->pop_size_per_gender_age_risk[g][ag][r] ; i++){
                    /* for each free sexual worker partnership of this individual */
                    for(j=pop->pop_per_gender_age_risk[g][ag][r][i]->n_clients ; j<pop->pop_per_gender_age_risk[g][ag][r][i]->max_n_clients ; j++){
                        /* Adding this person as an available sexual worker in the list of available sexual workers */
                        overall_partnerships->pop_available_sexual_workers->pop_per_patch_gender_age_risk[p][g][ag][r][overall_partnerships->n_pop_available_sexual_workers->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r]] = pop->pop_per_gender_age_risk[g][ag][r][i];
                        /* Keeping track of where this person is in that list */
                        pop->pop_per_gender_age_risk[g][ag][r][i]->idx_available_sexual_worker[j] = overall_partnerships->n_pop_available_sexual_workers->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r];
                        overall_partnerships->n_pop_available_sexual_workers->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r]++;
                    }
                    /* for each free normal partnership of this individual */
                    for(j=pop->pop_per_gender_age_risk[g][ag][r][i]->n_normal_partners ; j<pop->pop_per_gender_age_risk[g][ag][r][i]->max_n_normal_partners ; j++){
                        /* Adding this person as an available normal partner in the list of available normal partners */
                        overall_partnerships->pop_available_normal_partners->pop_per_patch_gender_age_risk[p][g][ag][r][overall_partnerships->n_pop_available_normal_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r]] = pop->pop_per_gender_age_risk[g][ag][r][i];
                        /* Keeping track of where this person is in that list */
                        pop->pop_per_gender_age_risk[g][ag][r][i]->idx_available_normal_partner[j] = overall_partnerships->n_pop_available_normal_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r];
                        overall_partnerships->n_pop_available_normal_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r]++;
                    }
                }
            }
        }
    }
}

/* Set all cumulative counters to zero (number of HIV tests, CD4 tests etc. ) : */
void init_cumulative_counters(cumulative_outputs_struct *cumulative_outputs){
    cumulative_outputs->N_total_CD4_tests_nonpopart = 0;
    cumulative_outputs->N_total_HIV_tests_nonpopart = 0;
    cumulative_outputs->N_total_CD4_tests_popart = 0;
    cumulative_outputs->N_total_HIV_tests_popart = 0;
    cumulative_outputs->N_total_HIV_tests_popart_positive = 0;
    cumulative_outputs->N_total_HIV_tests_popart_negative = 0;
    cumulative_outputs->N_total_ever_started_ART_nonpopart = 0;
    cumulative_outputs->N_total_ever_started_ART_popart = 0;
}


/* init_calendar_counters()
Set all array elements of arrays in the calendar outputs to zero 
(number of HIV tests, CD4 tests etc. ).  Used in cost-effectiveness analysis.  

Arguments
---------
calendar_outputs : pointer to calendar_outputs_struct
    Patch-specific structure housing a number of arrays that count yearly events (index starts
    from first year of the simulation; i.e. a[0] is for year 1900 if the simulation starts 
    in the year 1900).  

Returns
-------
Nothing; sets the values of the calendar array elements to zero.  
*/
void init_calendar_counters(calendar_outputs_struct *calendar_outputs){
    
    int i, g, a, cd4, spvl;
    for (i = 0; i < MAX_N_YEARS; i++){
        calendar_outputs->N_calendar_CD4_tests_nonpopart[i] = 0;
        calendar_outputs->N_calendar_HIV_tests_nonpopart[i] = 0;
        calendar_outputs->N_calendar_CD4_tests_popart[i] = 0;
        calendar_outputs->N_calendar_HIV_tests_popart[i] = 0;
        calendar_outputs->N_calendar_HIV_tests_popart_positive[i] = 0;
        calendar_outputs->N_calendar_HIV_tests_popart_negative[i] = 0;
        calendar_outputs->N_calendar_started_ART_annual[i] = 0;
        calendar_outputs->N_calendar_dropout[i] = 0;
        calendar_outputs->N_calendar_CHIPS_visits[i] = 0;
        calendar_outputs->N_calendar_VMMC[i] = 0;
        
        // For use with TREATS output
        for(g = 0; g < N_GENDER; g++){
            for(a = 0; a < N_AGE_UNPD + 1; a++){
                calendar_outputs->N_calendar_infections[g][a][i] = 0;
                for(cd4 = 0; cd4 < NCD4 + 1; cd4++){
                    for(spvl = 0; spvl < NSPVL + 1; spvl++){
                        calendar_outputs->N_calendar_started_ART[g][a][cd4][spvl][i] = 0;
                        calendar_outputs->N_calendar_Died_from_HIV_OnART[g][a][cd4][spvl][i] = 0;
                        calendar_outputs->N_calendar_Died_from_HIV_ARTNaive[g][a][cd4][spvl][i] = 0;
                        calendar_outputs->N_calendar_AnnualDropoutOnART[g][a][cd4][spvl][i] = 0;
                    }
                }
            }
        }
    }
}


void initialise_debug_variables(debug_struct *debug){
    int year, age_f, age_m, risk_f, risk_m;
    int p, art_i, art_j;

    for(year=0 ; year<MAX_N_YEARS ; year++)
    {
        for(age_f=0 ; age_f<N_AGE ; age_f++)
        {
            for(age_m=0 ; age_m<N_AGE ; age_m++)
            {
                debug->age_of_partners_at_partnership_formation[year][age_f][age_m] = 0;
                debug->age_of_partners_cross_sectional[year][age_f][age_m] = 0;
            }
        }

        for(risk_f=0 ; risk_f<N_RISK ; risk_f++)
        {
            for(risk_m=0 ; risk_m<N_RISK ; risk_m++)
            {
                debug->risk_of_partners_at_partnership_formation[year][risk_f][risk_m] = 0;
                debug->risk_of_partners_cross_sectional[year][risk_f][risk_m] = 0;
            }
        }
    }

    for (p=0;p<NPATCHES;p++){
        debug->art_vars[p].n_start_emergency_art_fromuntested = 0;
        debug->art_vars[p].n_start_emergency_art_fromartnaive = 0;
        debug->art_vars[p].n_start_emergency_art_fromartdroupout = 0;
        debug->art_vars[p].n_start_emergency_art_fromcascadedropout = 0;
        debug->art_vars[p].n_start_emergency_art = 0;
        for (art_i=0;art_i<NARTEVENTS;art_i++){
            for (art_j=0;art_j<NARTEVENTS;art_j++){
                debug->art_vars[p].cascade_transitions[art_i][art_j] = 0;
            }
        }
    }
}


