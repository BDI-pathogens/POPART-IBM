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

/* 
Functions related to PopART/CHiPs:
create_popart_chips_samples()
    at the start of a PopART year, decide who will be visited by CHiPs and in what order.
schedule_chips_visits()
    decide how many of these people  will be visited each timestep.
carry_out_chips_visits_per_timestep()
    goes through the list of people being visited by CHiPs and sets up their next cascade event 
    accordingly (HIV test, CD4 test).
chips_visit_person()
    

Functions related to circumcision:
draw_if_VMMC()
    Decides if an individual will get VMMC, and if they do, calls schedule_vmmc() to schedule 
    circumcision for a future timestep in vmmc_events[].  
schedule_vmmc()
    Adds an event for circumcision to vmmc_events[] to schedule an individual to get circumcised in
    future.
schedule_vmmc_healing()
    Adds an event for someone to finish healing after a VMMC op to vmmc_events[].
finish_vmmc_healing()
    Updates a person's status when they have finished healing after a VMMC op. At this point they
    have no further VMMC events happening to them.
schedule_generic_vmmc_event()
    Does the actual adding of a VMMC event to vmmc_event[] - function is called by schedule_vmmc 
    and schedule_vmmc_healing.
carry_out_VMMC_events_per_timestep()
    Carry out any event associated with VMMC in the current time step.  
 */

/************************************************************************/
/******************************* Includes  ******************************/
/************************************************************************/

#include "interventions.h"
#include "structures.h"
#include "constants.h"
#include "hiv.h"

/************************************************************************/
/******************************** functions *****************************/
/************************************************************************/


void create_popart_chips_samples(age_list_struct *age_list, chips_sample_struct *chips_sample, 
    parameters *param, int chips_round, int p){
    
    /*
    This takes the population of currently alive people (using age_list) and firstly sub-divides
    them into a 'chips_sampling_frame', e.g. dividing up men and women, as CHiPs tends to visit
    more women than men.  We can also exclude certain people (e.g. <15 years old) as needed.  The
    sampling frame is then drawn from to pick chips_sample->size_n_m men and chips_sample->size_n_f
    women who will be the people visited in the year.  We then shuffle each of these lists, so the
    shuffled list will be the order in which people are visited.  Finally we call
    schedule_chips_visits() which sets the number of people to be visited in each timestep so that
    the total number of people visited in a year adds up to the correct total.
    
    NOTE: The way we divide up the population means that we may try to visit people who died during
    the year.  However, this should not be a big factor, and I think it may even mimic CHiPs in
    that people may move/die between enumeration/mapping and CHiPs visit. 
    
    Arguments
    ----------
    age_list : pointer to age_list_struct structure
        
    chips_sample : pointer to chips_sample_struct structure
        
    param : pointer to parameters structure
        Parameters structure
    chips_round : int
        CHiPs round of interest
    p : int
        Patch number
    
    Returns
    -------
    
    */
    int g;
    int aa, ai,i, ac;
    /* For use with FOLLOW_INDIVIDUAL - we store these the first time we find them so we can find 
    them easily next time: */
    int g_persontofollow = -1; /* Default value indicates that the FOLLOW_INDIVIDUAL did not turn 
        up when going through - this is because they are too young to be visited by CHiPs. */
    int ac_persontofollow = -1;
    
    if(chips_round < -1 || chips_round >= NCHIPSROUNDS){
        printf("ERROR: The calculated CHiPS round is %d, outside the range -1 to %d. Exiting\n",
            chips_round, NCHIPSROUNDS - 1);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    /* This is a temporary store of the sampling frame. Only called annually at present so should
    be OK as a local variable. */
    chips_sample_struct *csf; // csf = "chips sampling frame"
    csf = malloc(sizeof(chips_sample_struct)); 
    
    if(csf == NULL){ /* Check memory allocated successfully. */
        printf("Unable to allocate csf in create_popart_chips_samples().");
        printf(" Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    // Set counters to zero (we use these to index within each array).
    for(g = 0; g < N_GENDER; g++){
        /* ac is the var we use as the index for age group in the chips_sample_struct structure.*/
        for(ac = 0; ac < (MAX_AGE - AGE_CHIPS + 1); ac++){
            csf->number_to_visit[g][ac] = 0;
        }
    }
    
    for(g = 0; g < N_GENDER; g++){
        /* aa is age index for age_list by gender (not adjusting for youngest age group).
        Here it is chosen to correspond to ages from AGE_CHIPS to 79. */
        for(aa = (AGE_CHIPS - AGE_ADULT); aa < (MAX_AGE - AGE_ADULT); aa++){
            
            /* ai is the index of the array age_list->number_per_age_group of the age group of 
            people you want to be dead */
            ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa;
            while(ai > (MAX_AGE - AGE_ADULT - 1)){
                ai = ai - (MAX_AGE-AGE_ADULT);
            }
            
            for(i = 0; i < age_list->age_list_by_gender[g]->number_per_age_group[ai]; i++){
                /* This is the relationship between the index aa and ac (the index for the 
                chips_sample_struct structures). */
                ac = aa - (AGE_CHIPS - AGE_ADULT);
                
                /* For debugging; check if we're following an individual and a patch. */
                if(
                (age_list->age_list_by_gender[g]->age_group[ai][i]->id == FOLLOW_INDIVIDUAL) &&
                (age_list->age_list_by_gender[g]->age_group[ai][i]->patch_no == FOLLOW_PATCH)
                ){
                    printf("Possible CHiPs visit %ld %d %d in round %d \n",
                        age_list->age_list_by_gender[g]->age_group[ai][i]->id,ai,i,chips_round);
                    
                    fflush(stdout);
                    /* Now store their characteristics so it's easier to find them: */
                    g_persontofollow = g;
                    ac_persontofollow = ac;
                }
                /* Also for debugging. */
                if(
                (age_list->age_list_by_gender[g]->age_group[ai][i]->cd4 == DUMMYVALUE) ||
                (age_list->age_list_by_gender[g]->age_group[ai][i]->cd4 == DEAD)
                ){
                    printf("Error -trying to schedule CHiPs visit for dead/non-existent ");
                    printf("person %ld\n",
                        age_list->age_list_by_gender[g]->age_group[ai][i]->id);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                
                /* At the moment, just split by gender and age. 
                Can also add in some undersampling of high risk groups. */
                csf->list_ids_to_visit[g][ac][csf->number_to_visit[g][ac]] = age_list->age_list_by_gender[g]->age_group[ai][i]->id;
                
                csf->number_to_visit[g][ac]++;
            }
        }
        
        /* Now look at oldest age group: */
        for(i = 0; i < age_list->age_list_by_gender[g]->number_oldest_age_group; i++){
            
            if(
            (age_list->age_list_by_gender[g]->oldest_age_group[i]->id == FOLLOW_INDIVIDUAL) &&
            (age_list->age_list_by_gender[g]->oldest_age_group[i]->patch_no==FOLLOW_PATCH)
            ){
                
                printf("Possible CHiPs visit %ld oldest age gp i=%d\n",
                    age_list->age_list_by_gender[g]->oldest_age_group[i]->id,i);
                fflush(stdout);
                /* Now store their characteristics so it's easier to find them: */
                g_persontofollow = g;
                /* This is MAX_AGE-AGE_ADULT-(AGE_CHIPS-AGE_ADULT). */
                ac_persontofollow = MAX_AGE-AGE_CHIPS;
            }
            csf->list_ids_to_visit[g][MAX_AGE - AGE_CHIPS][csf->number_to_visit[g][MAX_AGE - AGE_CHIPS]] =
                    age_list->age_list_by_gender[g]->oldest_age_group[i]->id;
            csf->number_to_visit[g][MAX_AGE - AGE_CHIPS]++;
        }
    }
    
    /* These decide how many people are going to be visited by CHiPs. Can be implemented as a
    number rather than a percentage of population. */
    for(g = 0; g < N_GENDER; g++){
        /* Run from AGE_CHIPS to 80+ (note - this is different from previous loop). */
        for(ac = 0; ac < (MAX_AGE - AGE_CHIPS + 1); ac++){
            if(chips_round >= 0 && chips_round < NCHIPSROUNDS){
                chips_sample->number_to_visit[g][ac] = (int)  floor(csf->number_to_visit[g][ac]*param->chips_params->prop_tested_by_chips_in_round[g][ac][chips_round]);
            }else if (chips_round==CHIPSROUNDPOSTTRIAL){
                chips_sample->number_to_visit[g][ac] = (int)  floor(csf->number_to_visit[g][ac]*param->chips_params->prop_tested_by_chips_in_round_posttrial[g][ac]);
            }else{
                printf("ERROR: Unknown value of chips_round=%d. Exiting\n", chips_round);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            
            if (chips_sample->number_to_visit[g][ac]>0){
                
                /* Choose chips_sample->number_to_visit[g][ac] from the array csf->m. */
                gsl_ran_choose(rng, chips_sample->list_ids_to_visit[g][ac],
                    chips_sample->number_to_visit[g][ac],
                    csf->list_ids_to_visit[g][ac],
                    csf->number_to_visit[g][ac], sizeof (long));
                
                /* Randomise the order (as gsl_ran_choose maintains the order of 
                the original list). */
                gsl_ran_shuffle(rng, chips_sample->list_ids_to_visit[g][ac],
                    chips_sample->number_to_visit[g][ac], sizeof (long));
            }
        }
    }
    
    /* Check to see if this person was visited. Note that if g_persontofollow==-1 then they are too
    young to be in csf->list_ids_to_visit. */
    
    if(p == FOLLOW_PATCH && g_persontofollow > -1){
        for(i = 0; 
        i < csf->number_to_visit[g_persontofollow][ac_persontofollow]; 
        i++){
            
            if(csf->list_ids_to_visit[g_persontofollow][ac_persontofollow][i] == FOLLOW_INDIVIDUAL){
                
                printf("CHiPs visit for adult %ld, gender %i from patch %d now scheduled\n",
                    csf->list_ids_to_visit[g_persontofollow][ac_persontofollow][i],
                    g_persontofollow, p);
                fflush(stdout);
            }
        }
    }
    schedule_chips_visits(chips_sample, param, chips_round);
    free(csf);
}


void schedule_chips_visits(chips_sample_struct *chips_sample, parameters *param, int chips_round){
   /* Given a sample of people who are to be visited each year (currently chips_sample->m and
    chips_sample->f) schedule their visits in the arrays
    chips_sample->number_to_see_per_timestep_m/f.  These arrays contain the number of people to see
    at each timestep, and we run through e.g. chips_sample->m until we have seen that many people
    each timestep.  
    
    
    Arguments
    ---------
    chips_sample : pointer to a chips_sample_struct struct
    param : pointer to a parameters struct
        All the parameters of interest for the patch in question.  
    chips_round : int
        CHiPs round of interest
    
    Returns
    -------
    Nothing; 
    
    */
    
    
    
    int g,ac;
    double temp_chips_expected_cumulative_proportion_visited = 0;
    double temp_chips_expected_cumulative_number_visited = 0;
    long temp_chips_cumulative_number_scheduled;
    /* We use this to store the fraction of people who will be visited in a round who are visited
    in a given timestep. */
    double temp_fraction_visited_normalised; 

    /* We need to deal with rounding issues (described more below). To do this we ensure we are
    always within 1 of the expected cumulative number of visits.  However computer rounding means
    we need to adjust that very slightly, so subtract an arbitrary small amount. */
    double TOLERANCE = 1.0-1e-9;

    /* Initialise values in chips_sample->next_person_to_see[] so begin at the start of the list. */
    for(g = 0; g < N_GENDER; g++){
        
        /* Run from AGE_CHIPS to 80+. */
        for(ac = 0; ac <(MAX_AGE - AGE_CHIPS + 1); ac++){
            chips_sample->next_person_to_see[g][ac] = 0;
        }
    }

    // Determine the number of time steps for the chips round in question
    int n_steps_in_round;
    if(chips_round >= 0 && chips_round < NCHIPSROUNDS){
        n_steps_in_round = param->chips_params->n_timesteps_per_round[chips_round];
    }else if (chips_round == CHIPSROUNDPOSTTRIAL){
        n_steps_in_round = param->chips_params->n_timesteps_per_round_posttrial;
    }else{
        printf("ERROR: Unknown chips_round value =%d. Exiting\n", chips_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    int i;
    /* This is just a crude way to make sure we see everyone we're supposed to see that month: */ 
    for(g = 0; g < N_GENDER; g++){
        /* CHiPs visit everyone from from AGE_CHIPS to 80+. */
        for(ac = 0; ac < (MAX_AGE - AGE_CHIPS + 1); ac++){
            temp_chips_expected_cumulative_proportion_visited = 0;
            temp_chips_expected_cumulative_number_visited = 0;
            temp_chips_cumulative_number_scheduled = 0;
            
            for(i = 0; i < n_steps_in_round; i++){
                /* This is the fraction of people who will be visited in a round who are visited in
                a given timestep. */
                if(chips_round >= 0 && chips_round < NCHIPSROUNDS){
                    temp_fraction_visited_normalised = param->chips_params->prop_tested_by_chips_per_timestep[g][ac][i][chips_round]/param->chips_params->prop_tested_by_chips_in_round[g][ac][chips_round];
                }else if (chips_round == CHIPSROUNDPOSTTRIAL){
                    temp_fraction_visited_normalised = param->chips_params->prop_tested_by_chips_per_timestep_posttrial[g][ac][i]/param->chips_params->prop_tested_by_chips_in_round_posttrial[g][ac];
                }
                
                chips_sample->number_to_see_per_timestep[g][ac][i] = (int) round(
                    chips_sample->number_to_visit[g][ac] * temp_fraction_visited_normalised);
                
                /* Now deal with the fact that there are rounding issues - we may end up seeing
                more/less people than intended as a result. For example if we are to see 10 people
                over 9 timesteps equally at a rate of 1/9 then we would see 1 at each timestep and
                hence only 9 people at the end.  This code below will find a place to add the extra
                person, in such a way that the adjustments are distributed reasonably over the
                whole CHiPs round.     */
                
                /* This is the number of people we have scheduled so far: */
                temp_chips_cumulative_number_scheduled += chips_sample->number_to_see_per_timestep[g][ac][i];

                temp_chips_expected_cumulative_proportion_visited += temp_fraction_visited_normalised;
                /* This is the number of people we 'should' have visited so far. 
                We allow it to be non-integer. */
                temp_chips_expected_cumulative_number_visited = chips_sample->number_to_visit[g][ac] * temp_chips_expected_cumulative_proportion_visited;
                
                while(
                (temp_chips_cumulative_number_scheduled -
                temp_chips_expected_cumulative_number_visited) >= TOLERANCE &&
                (chips_sample->number_to_see_per_timestep[g][ac][i] >=1)
                ){
                    chips_sample->number_to_see_per_timestep[g][ac][i] -= 1;
                    temp_chips_cumulative_number_scheduled -= 1;
                }
                /* Print a warning message if the issue is not fixed: */
                if(
                (temp_chips_cumulative_number_scheduled - temp_chips_expected_cumulative_number_visited) >= TOLERANCE){
                    printf("Warning - may need to adjust schedule_chips_visits() at ");
                    printf("i=%i = %li %6.4lf - have scheduled incorrectly?\n", i,
                        temp_chips_cumulative_number_scheduled,
                        temp_chips_expected_cumulative_number_visited);
                }

                /* If we have scheduled too few people: */
                while((temp_chips_expected_cumulative_number_visited - temp_chips_cumulative_number_scheduled) >= TOLERANCE){
                    chips_sample->number_to_see_per_timestep[g][ac][i] += 1;
                    temp_chips_cumulative_number_scheduled += 1;
                }
                /* Print a warning message if the issue is not fixed: */
                if((temp_chips_expected_cumulative_number_visited - temp_chips_cumulative_number_scheduled) >= TOLERANCE){
                    printf("Warning - may need to adjust schedule_chips_visits() at ");
                    printf("i=%i = %li %6.4lf - have scheduled incorrectly?\n", i,
                        temp_chips_cumulative_number_scheduled,
                        temp_chips_expected_cumulative_number_visited);
                }
            }
        }
    }
}


/*************************** Functions which do the CHiPs visits: ************************/


void carry_out_chips_visits_per_timestep(int t0, int t_step, patch_struct *patch, int p, 
        int chips_round, debug_struct *debug, output_struct *output){
    /* Carry out the CHiPS visits for a given timestep at time t. 
    
    Arguments
    ---------
    t0 : int
    t_step : int
    patch : patch_struct struct
    p : int
        Patch identifier.  
    chips_round : int
        CHiPs round in question
    debug : pointer to a debug_struct struct
    output : pointer to an output_struct struct
    
    Returns
    -------
    Nothing; 
    */
    long i;  
    int g,ac;
    int n_steps_in_round;
    if(chips_round >= 0 && chips_round < NCHIPSROUNDS){
        n_steps_in_round = patch[p].param->chips_params->n_timesteps_per_round[chips_round];
    }else if(chips_round == CHIPSROUNDPOSTTRIAL){
        n_steps_in_round = patch[p].param->chips_params->n_timesteps_per_round_posttrial;
    }else{
        printf("ERROR: Unknown chips_round value =%d. Exiting\n", chips_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    /* t_i is the index of the timestep within chips_sample corresponding to the current time t.
    For example chips_sample->number_to_see_per_timestep_m[t_i] is no. men to visit at time t. */
    int t_i;
    
    /* If the trial has not ended yet: */
    if(chips_round > -1){
        t_i = (t0 - patch[p].param->CHIPS_START_YEAR[chips_round]) * N_TIME_STEP_PER_YEAR +
            (t_step-patch[p].param->CHIPS_START_TIMESTEP[chips_round]);
    }else{
        /* Assume that from end of trial onwards CHiPs rounds are ANNUAL, and that they start 1
        week after the last timestep of the previous round. */
        
        /* Copy the first timestep from the previous round. */
        t_i = t_step - patch[p].param->CHIPS_START_TIMESTEP_POSTTRIAL;
        if(t_i < 0){
            t_i = t_i + N_TIME_STEP_PER_YEAR;
        }
    }
    
    /* For debugging: */
    if(t_i >= MAX_N_TIMESTEPS_PER_CHIPS_ROUND){
        printf("Problem - MAX_N_TIMESTEPS_PER_CHIPS_ROUND is too small. ");
        printf("We are %i timesteps into CHiPs round %i. Exiting\n", t_i, chips_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    if(t_i < 0 || t_i >= n_steps_in_round){
        printf("Error: t_i=%i is outside range [0,%d]", t_i, n_steps_in_round);
        printf(" in CHiPs round %d at time %6.4lf. ", chips_round, t0 + t_step * TIME_STEP);
        printf("Exiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    /* Now CHiPs visits each sub-population (currently men/women by year age group ac) in turn: */
    /* Go through the list of id's of people to visit, and visit them: */
    for(g = 0; g < N_GENDER; g++){
        /* Run from AGE_CHIPS to 80+. */
        for (ac = 0; ac < (MAX_AGE - AGE_CHIPS + 1); ac++){
            
            for(i = patch[p].chips_sample->next_person_to_see[g][ac]; 
                (i < (patch[p].chips_sample->next_person_to_see[g][ac] +
                    patch[p].chips_sample->number_to_see_per_timestep[g][ac][t_i])); 
                 i++){
                    if((patch[p].chips_sample->next_person_to_see[g][ac] +
                        patch[p].chips_sample->number_to_see_per_timestep[g][ac][t_i]) > 
                    (patch[p].chips_sample->number_to_visit[g][ac])){
                    
                    /* Only print out error message if we are over 2 too many otherwise just skip
                    those people. */
                    if((patch[p].chips_sample->next_person_to_see[g][ac] + 
                        patch[p].chips_sample->number_to_see_per_timestep[g][ac][t_i]) >
                        (patch[p].chips_sample->number_to_visit[g][ac])){
                        
                        if(VERBOSE_OUTPUT == 1){
                        printf("Visited too many people in patch %i at ", p);
                        printf("time=%6.4f g=%i ac=%i Excess=%li Number = %li\n",
                            t0 + t_step * TIME_STEP, g, ac, 
                                                        patch[p].chips_sample->next_person_to_see[g][ac]+patch[p].chips_sample->number_to_see_per_timestep[g][ac][t_i] - patch[p].chips_sample->number_to_visit[g][ac],
                            patch[p].chips_sample->number_to_visit[g][ac]);
                        }
                    }
                    
                    }else{
                    
                    /* Send the address (ie pointer) to this person. */
                    chips_visit_person(&(patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]]), 
                        patch[p].cumulative_outputs,
                        patch[p].calendar_outputs,
                        t0 + t_step*TIME_STEP,
                        patch[p].cascade_events,
                        patch[p].n_cascade_events,
                        patch[p].size_cascade_events,
                        patch[p].hiv_pos_progression,
                        patch[p].n_hiv_pos_progression,
                        patch[p].size_hiv_pos_progression,
                        patch[p].param,
                        patch[p].vmmc_events,
                        patch[p].n_vmmc_events,
                        patch[p].size_vmmc_events, 
                        patch, p, chips_round, debug, output, g, ac);
                }
            }
            /* Update this index ready for the next timestep: */
            patch[p].chips_sample->next_person_to_see[g][ac] +=
                patch[p].chips_sample->number_to_see_per_timestep[g][ac][t_i];
        }
    }
}


void chips_visit_person(individual *indiv, cumulative_outputs_struct *cumulative_outputs,
    calendar_outputs_struct *calendar_outputs, double t, individual ***cascade_events, 
    long *n_cascade_events, long *size_cascade_events, individual ***hiv_pos_progression, 
    long *n_hiv_pos_progression, long *size_hiv_pos_progression, parameters *param, 
    individual ***vmmc_events, long *n_vmmc_events, long *size_vmmc_events, patch_struct *patch,
    int p, int chips_round, debug_struct *debug, output_struct *output, int g, int ac){
    /* 
    Because of the way we draw CHiPs visits at the beginning of the year, it is possible some people
    die before they are visited. If this is the case then do nothing more. 
    They are deleted from age_list so won't be in next year's sample. 
    
    
    Arguments
    ---------
    individual *indiv
    cumulative_outputs_struct *cumulative_outputs
    double tindividual ***cascade_events
    long *n_cascade_events
    long *size_cascade_events
    individual ***hiv_pos_progression
    long *n_hiv_pos_progression
    long *size_hiv_pos_progression
    parameters *param
    individual ***vmmc_events
    long *n_vmmc_events
    long *size_vmmc_events
    patch_struct *patch
    int p
    int chips_round
    debug_struct *debug
    output_struct *output
    int g
    int ac
    
    
    Returns
    -------
    Nothing;
    
    */
    
    /* We need a different variable as chips_round_including_end_trial is used as an array index. */
    int chips_round_including_end_trial;
    int year_idx = (int) floor(t) - param->start_time_simul;
    
    if(chips_round >= 0 && chips_round < NCHIPSROUNDS){
        chips_round_including_end_trial = chips_round;
    }else if(chips_round == CHIPSROUNDPOSTTRIAL){
        /* VMMC parameters should be like last round. */
        chips_round_including_end_trial = NCHIPSROUNDS - 1;
    }else{
        printf("ERROR: Unknown chips_round value =%d. Exiting\n", chips_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    if(indiv->cd4 == DUMMYVALUE){
        printf("Trying to CHiPS visit a non-existent person %d %ld !!! Exiting\n", p, indiv->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    if(indiv->cd4 == DEAD){
        indiv->VISITED_BY_CHIPS_THISROUND = DIEDBEFORECHIPSVISIT;
        return;
    }
    /* Save what was the old scheduled cascade event type in case needed: */
    int old_cascade_event = indiv->next_cascade_event;
    
    /* Update counters for whether visited this round and number of lifetime visits: */
    indiv->VISITED_BY_CHIPS_THISROUND = TRUE;
    
    /* We only record this during the trial. */
    if(chips_round > -1){
        output->NCHIPS_VISITED[p][g][ac][chips_round]++;
        
        if(indiv->HIV_status > UNINFECTED){
            output->NCHIPS_HIVPOS[p][g][ac][chips_round]++;
            if(
            (indiv->ART_status >= ARTNAIVE) && 
            (indiv->ART_status < ARTDEATH)
            ){
                output->NCHIPS_HIVAWARE[p][g][ac][chips_round]++;
                
                if(
                (indiv->ART_status == EARLYART) || 
                (indiv->ART_status == LTART_VS) || 
                (indiv->ART_status == LTART_VU)
                ){
                    output->NCHIPS_ONART[p][g][ac][chips_round]++;
                    
                    if(indiv->ART_status == LTART_VS){
                        output->NCHIPS_VS[p][g][ac][chips_round]++;
                    }
                }
            }
        }
    }
    indiv->NCHIPSVISITS++;
    
    // Record that there was a chips visit this year
    calendar_outputs->N_calendar_CHIPS_visits[year_idx]++;
    
    /* Are we following specific individuals or patches? */
    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("CHiPs visit for adult %ld from patch %d at time %lf with old_cascade_event %d\n",
            indiv->id, indiv->patch_no, t, old_cascade_event);
        fflush(stdout);
    }
    
    /* If their previous event was an HIV test, then at this step they are HIV tested by CHiPs. */
    if( (old_cascade_event == CASCADEEVENT_HIV_TEST_NONPOPART) 
            || (old_cascade_event == CASCADEEVENT_HIV_TEST_POPART)){
        
        /* Unschedule the current event from care cascade. */
        remove_from_cascade_events(indiv, cascade_events, n_cascade_events, 
                size_cascade_events, t, param);
        
        indiv->next_cascade_event = CASCADEEVENT_HIV_TEST_POPART;
        
        hiv_test_process(indiv, param, t, cascade_events, n_cascade_events, size_cascade_events, 
            hiv_pos_progression, n_hiv_pos_progression, size_hiv_pos_progression, 
            cumulative_outputs, calendar_outputs, vmmc_events, n_vmmc_events, size_vmmc_events,
            patch, p, debug);
        
        return;
    }else if(
        (old_cascade_event == CASCADEEVENT_CD4_TEST_NONPOPART) || 
        (old_cascade_event == CASCADEEVENT_CD4_TEST_POPART)
    ){
        /* Unschedule the current event from care cascade. */
        remove_from_cascade_events(indiv, cascade_events, n_cascade_events, 
            size_cascade_events, t, param);
        
        /* Eligibility for ART is determined by calendar time and trial arm.
        * If in arm A and has not dropped out then schedule start of ART. */
        if (patch[p].trial_arm == ARM_A){
            cumulative_outputs->N_total_CD4_tests_popart++;
            calendar_outputs->N_calendar_CD4_tests_popart[year_idx]++;
            indiv->ART_status = ARTNAIVE;
            indiv->next_cascade_event = CASCADEEVENT_START_ART_POPART;
            indiv->VISITEDBYCHIPS_TO_INIT_ART = 1;
            schedule_start_of_art(indiv, param, t, cascade_events, 
                n_cascade_events, size_cascade_events);
        
        /* In arm B but eligible for ART: */
        }else if (is_eligible_for_art(indiv, param, t, patch, p) > 0){
            
            if (patch[p].trial_arm != ARM_B){
                printf("ERROR: Not in arm B 1??? Exiting\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            
            indiv->ART_status = ARTNAIVE;
            cumulative_outputs->N_total_CD4_tests_popart++;
            calendar_outputs->N_calendar_CD4_tests_popart[year_idx]++;
            indiv->next_cascade_event = CASCADEEVENT_START_ART_POPART;
            indiv->VISITEDBYCHIPS_TO_INIT_ART = 1;
            schedule_start_of_art(indiv, param, t, cascade_events, 
                n_cascade_events, size_cascade_events);
        
        /* In arm B and not yet eligible for ART: */
        }else{
            if (patch[p].trial_arm != ARM_B){
                printf("ERROR: Not in arm B 2??? Exiting\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            
            /* Time to next CD4 test is again the sum of the time between 
            getting the HIV test and having the first CD4 test, and the time between consecutive
            CD4 tests. */
            cumulative_outputs->N_total_CD4_tests_popart++;
            calendar_outputs->N_calendar_CD4_tests_popart[year_idx]++;
            indiv->ART_status = ARTNAIVE;
            
            /* Note that we have two events happening here - time to get the CD4 test (to determine
            eligibility) and then time to retest CD4 again - we are scheduling the future CD4 test
             here. So I think this is correct. */
            
            double time_new_cd4 = t + param->t_delay_hivtest_to_cd4test_min[POPART] +
                param->t_delay_hivtest_to_cd4test_range[POPART] * gsl_rng_uniform (rng);
            
            indiv->next_cascade_event = CASCADEEVENT_CD4_TEST_POPART;
            schedule_generic_cascade_event(indiv, param, time_new_cd4, cascade_events,
                n_cascade_events, size_cascade_events,t);
        }
    }else if(old_cascade_event == NOEVENT){
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Adult %ld at time %lf is in old_cascade_event==NOEVENT if statement\n",
                indiv->id, t);
            fflush(stdout);
        }
        
        /* if already on ART (and has NOEVENT because of this - ie early ART or VS) then we don't
        need to do anything else. */
        if(
            (indiv->ART_status == EARLYART) || 
            (indiv->ART_status == LTART_VS)
        ){
            if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
                printf("Adult %ld at time %lf is already on ART so this is OK.\n", indiv->id, t);
                fflush(stdout);
            }
        }else{
            
            /* If individual has never tested positive this person needs to get an HIV test next.
             (note that this is tested, not serostatus). */
            if(indiv->ART_status == ARTNEG){
                
                // Assign next cascade event as PopART HIV test
                indiv->next_cascade_event = CASCADEEVENT_HIV_TEST_POPART;
                
                hiv_test_process(indiv, param, t, cascade_events, n_cascade_events,
                    size_cascade_events, hiv_pos_progression, n_hiv_pos_progression,
                    size_hiv_pos_progression, cumulative_outputs, calendar_outputs, 
                    vmmc_events, n_vmmc_events, size_vmmc_events, patch, p, debug);
            
            }else{
                
                /* Check that individual can only get here if previously dropped out. */
                if((indiv->ART_status != CASCADEDROPOUT) && (indiv->ART_status != ARTDROPOUT)){
                    
                    printf("ERROR: Incorrect cascade in chips_visit_person for ");
                    printf("ART status =%d for person %ld in patch %d. ",
                        indiv->ART_status, indiv->id, indiv->patch_no);
                    
                    printf("Exiting\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                /* Allow back into cascade if you previously dropped out at some given probability
                (note that this is currently 1): */
                if(
                gsl_ran_bernoulli(rng, 
                    param->p_popart_to_cascade[chips_round_including_end_trial]) == 1
                ){
                    /* If already know HIV+ then can go to clinic and get CD4 test etc. */
                    if(is_eligible_for_art(indiv, param, t, patch, p) > 0){
                        
                        cumulative_outputs->N_total_CD4_tests_popart++;
                        calendar_outputs->N_calendar_CD4_tests_popart[year_idx]++;
                        indiv->next_cascade_event = CASCADEEVENT_START_ART_POPART;
                        indiv->ART_status = ARTNAIVE;
                        indiv->VISITEDBYCHIPS_TO_INIT_ART = 1;
                        schedule_start_of_art(indiv, param, t, cascade_events, 
                            n_cascade_events, size_cascade_events);
                    }else{
                        /* 
                        Time to next CD4 test is again the sum of the time between
                        getting the HIV test and having the first CD4 test, and the time between
                        consecutive CD4 tests.
                        */
                        cumulative_outputs->N_total_CD4_tests_popart++;
                        calendar_outputs->N_calendar_CD4_tests_popart[year_idx]++;
                        indiv->ART_status = ARTNAIVE;
                        /*
                        Note that we have two events happening here - time to get the CD4 test 
                        (to determine eligibility) and then time to retest CD4 again - we are 
                        scheduling the future CD4 test here. So I think this is correct. 
                        */
                        double time_new_cd4 = t + param->t_delay_hivtest_to_cd4test_min[POPART] + 
                            param->t_delay_hivtest_to_cd4test_range[POPART] * gsl_rng_uniform (rng);
                        
                        indiv->next_cascade_event = CASCADEEVENT_CD4_TEST_POPART;
                        
                        schedule_generic_cascade_event(indiv, param, time_new_cd4, cascade_events, 
                            n_cascade_events, size_cascade_events, t);
                    }
                }
            }
        }
    }
}


/******************************************************************************************
 * These are VMMC intervention events - VMMC can either be through PopART or national 
 * policy/campaigns.
 ******************************************************************************************/
/* Determines if a man gets VMMC, and if so schedules the process: 
 * ASSUMPTION!!! - time is drawn with no data!!! */
void draw_if_VMMC(individual *indiv, parameters *param, individual ***vmmc_events, long *n_vmmc_events, long *size_vmmc_events, double t, int is_popart){
    double p_circ;
    int year, t_step, chips_round;

    /* For DEBUGGING: */
    if (indiv->gender==FEMALE||(indiv->circ!=UNCIRC)) 
    {
        printf("ERROR: This person %ld is gender=%d with circ=%d is in draw_if_VMMC. Exiting\n",
            indiv->id, indiv->gender, indiv->circ);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    /* Need to see if there is an HIV retest prior to VMMC. */
    if (t <= (param->COUNTRY_VMMC_START)){
        return;                       /* No VMMC so exit function. */
    }
    else{
        if (is_popart==POPART){
            year = (int) floor(t);
            t_step = (int) round((t - year)*N_TIME_STEP_PER_YEAR);
            chips_round = get_chips_round(param, year, t_step);

            /* Note that get_chips_round() returns -1 if after trial. We need to use chips_round as an array index so fix this: */
            if (chips_round==CHIPSROUNDPOSTTRIAL)
                chips_round = NCHIPSROUNDS-1; /* VMMC parameters should be like last round. */
            else if (chips_round<0 || chips_round>=NCHIPSROUNDS){
                printf("ERROR: Unknown chips_round value =%d in draw_if_VMMC() at t=%lf. Exiting\n",chips_round,t);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }

            //printf("Calling chips_sampling_frame() year=%d t_step=%d. chips_round=%d p=%lf\n",year,t_step,chips_round,param->p_circ_popart[chips_round]);
            p_circ= param->p_circ_popart[chips_round];

        }
        else
            p_circ= param->p_circ_nonpopart;
        /* Assume that probability of VMMC is higher with PopART (hence the index for p_circ): */

        if (gsl_ran_bernoulli(rng,p_circ)==1){
            //printf("Decided VMMC t=%lf is_popart=%d p=%lf\n",t,is_popart,p_circ);
            /* Schedule the person to get VMMC some time in the near future.
             * Note: we allow VMMC to happen quicker during PopART - hence pass is_popart to schedule_vmmc(): */
            schedule_vmmc(indiv, param, vmmc_events, n_vmmc_events, size_vmmc_events, t, is_popart);
        }
    }
}

/* Function is called when a person has just had a -ve HIV test, and decides to get VMMC at some time in the future.
 * The function determines when in the future the person will get VMMC and schedules this in vmmc_events[]. */
void schedule_vmmc(individual *indiv, parameters *param, individual ***vmmc_events, 
        long *n_vmmc_events, long *size_vmmc_events, double t, int is_popart){

    /* For DEBUGGING: */
    if (indiv->gender==FEMALE||(indiv->circ!=UNCIRC)) 
    {
        printf("ERROR: not sure why this person %ld gender=%d with circ=%d is in schedule_VMMC. Exiting\n",indiv->id,indiv->gender,indiv->circ);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    //printf("Individual %ld had been scheduled for VMMC at time %f, is_popart=%d\n",indiv->id,t,is_popart);
    double time_vmmc = t + param->t_get_vmmc_min[is_popart] + param->t_get_vmmc_range[is_popart]*gsl_rng_uniform (rng);
    //printf("Individual %ld had been scheduled for VMMC at t=%f for future time %f, is_popart=%d\n",indiv->id,t,time_vmmc,is_popart);
    /* Set status to waiting for VMMC: */
    indiv->circ = UNCIRC_WAITING_VMMC;
    schedule_generic_vmmc_event(indiv,param,vmmc_events,n_vmmc_events,size_vmmc_events,t,time_vmmc);
    /* Do we need to model people deciding more than once if they get VMMC? */
    return;
}


void schedule_vmmc_healing(individual *indiv, parameters *param, individual ***vmmc_events, 
    long *n_vmmc_events, long *size_vmmc_events, double t){
    /* Determine time between VMMC operation and healing, and schedules a healing event.
    
    This function is called when the VMMC operation is performed.   The healing event for the person
    in question is added to vmmc_events[], a time to healing is generated, and the individual's
    circ status is set to VMMC_HEALING.
    
    Arguments
    ---------
    indiv : pointer to an individual struct
    param : pointer to a parameter struct
    vmmc_events : multidimensional array of pointers to individual struct
    n_vmmc_events : array of long
        Array (of length N_TIME_STEP_PER_YEAR; see constants.h for def) of the number of VMMC events
        in a given time step for the current year.  
    size_vmmc_events : array of long
    t : double
        Current time in years.  
    
    Returns
    -------
    Nothing; a time to healing is drawn and schedule_generic_vmmc_event() is called.  
    */
    
    /* If the individual is female or they are not waiting for VMMC, throw and error */
    if(indiv->gender == FEMALE || (indiv->circ != UNCIRC_WAITING_VMMC)){
        printf("ERROR: not sure why this person %ld ", indiv->id);
        printf("gender=%d with circ=%d is in schedule_vmmc_healing.  ", indiv->gender, indiv->circ);
        printf("Exiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    double time_heal = t + param->t_vmmc_healing;
    
    // Change the individuals circumcision status to VMMC_HEALING
    indiv->circ = VMMC_HEALING;
    indiv->t_vmmc = t;
    
    //printf("Scheduling healing for %ld at time %lf for future time %lf\n",indiv->id,t,time_heal);
    schedule_generic_vmmc_event(indiv, param, vmmc_events, n_vmmc_events, size_vmmc_events, 
        t, time_heal);
    return;
}


/* Once someone has reached the end of the VMMC healing period, this function is called.
 * Function sets the individual so they no longer have any VMMC event, and their circ status is "VMMC". */
void finish_vmmc_healing(individual *indiv){
    
    /* For DEBUGGING: */
        if (indiv->gender==FEMALE||(indiv->circ!=VMMC_HEALING)) 
        {
            printf("ERROR: not sure why this person %ld gender=%d with circ=%d is in finish_vmmc_healing. Exiting\n",indiv->id,indiv->gender,indiv->circ);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
    indiv->circ = VMMC;
    /* Do not need to remove from vmmc_events[] array as this event is now in the past. */
    indiv->idx_vmmc_event[0] = -1;
    indiv->idx_vmmc_event[1] = -1;
    return;
}


void schedule_generic_vmmc_event(individual *indiv, parameters *param, individual ***vmmc_events,
    long *n_vmmc_events, long *size_vmmc_events, double t, double t_event){
    /*
    Add a VMMC event to the array vmmc_event[][] for an individual `indiv` at time `t_event`.  
    
    This function adds a VMMC event to the array vmmc_event for the time step at which it occurs in 
    the future, `t_step_vmmc_event`.  This function is called by schedule_vmmc() and
    schedule_vmmc_healing().
    
    Arguments
    ---------
    indiv : pointer to individual struct
        Individual for which the VMMC event is to be scheduled.  
    param : pointer to a parameters struct
    vmmc_events : 
    n_vmmc_events : array of long
    size_vmmc_events : array of long
    t : double
        Current time in years
    t_event : double
        Time in years at which the VMMC event is to take place.  
    
    Returns
    -------
    Nothing.  
    */
    
    /* Check the VMMC event doesn't happen before the end of the simulation */
    if (t_event <= param->end_time_simul){

        /* Note that if indiv->idx_vmmc_event[0] < 0 then the event is either "no event scheduled"
        (-1) or "an event was scheduled but only for after the end of the simulation" (-2).  
        Probably no way that -2 can ever happen but allow it as a possibility for now.  
         */
        if(
        (indiv->debug_last_vmmc_event_index == indiv->idx_vmmc_event[0]) &&
        (indiv->idx_vmmc_event[0] >= 0)
        ){
            printf("ERROR - trying to schedule a new vmmc event ");
            printf("(circumcision status=%d) in schedule_generic_vmmc_event() ", indiv->circ);
            printf("that occurs at the same time as the previous event ");
            printf("for person %ld in patch %d at time = %6.4f. ", indiv->id, indiv->patch_no, t);
            printf("Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        /* vmmc_events is an array of events occurring in the next 12 months only.  
        So the index at time t is calculated as follows:
        * Suppose VMMC starts in 2008.0. Then for the rest of 2008 we schedule as 
            (int) trunc((t-2008)*N_TIME_STEP_PER_YEAR).
        * However, after 2008 we want to reuse the same array. So in e.g. 2010.25 we take the
        decimal part (0.25) and then
        * (int) trunc(0.25*N_TIMESTEPS_PER_YEAR) gives the new array index. */  
        int t_step_vmmc_event = (int) (floor((t_event - floor(t_event)) * N_TIME_STEP_PER_YEAR));
        int t_step_current_time = (int) (floor((t - floor(t)) * N_TIME_STEP_PER_YEAR));
        
        /* Ensure that we never schedule a cascade event during the current timestep: */
        if (t_step_vmmc_event == t_step_current_time){
            if (t_step_vmmc_event < N_TIME_STEP_PER_YEAR){
                t_step_vmmc_event += 1;
            }else if (t_step_vmmc_event == N_TIME_STEP_PER_YEAR){
                t_step_vmmc_event = 0;
            }
        }
        
        if (t_step_vmmc_event < 0 || t_step_vmmc_event >= N_TIME_STEP_PER_YEAR){
            printf("Error - trying to schedule VMMC outside array indices %d\n",
                t_step_vmmc_event);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        if (indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Individual %ld is in schedule_generic_vmmc_event\n", indiv->id);
        }
        
        indiv->idx_vmmc_event[0] = t_step_vmmc_event;
        indiv->idx_vmmc_event[1] = n_vmmc_events[t_step_vmmc_event];
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("New generic VMMC event for adult %ld generated with array indices  %ld %ld\n",
                indiv->id, indiv->idx_vmmc_event[0], indiv->idx_vmmc_event[1]);
            fflush(stdout);
        }

        /* Check if we've run out of memory: */
        if (n_vmmc_events[t_step_vmmc_event] >= (size_vmmc_events[t_step_vmmc_event])){

            /* Note that realloc does not work (we need to pass a pointer to the pointer, which is
            really complicated as it propagates through several functions (so maybe make
            planned_breakups[time_breakup] **), so ditch this code for now and use the following
            lines: */
            printf("Unable to re-allocate vmmc_events[i]. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        /* Add a pointer to this individual to the vmmc_events array in the correct location */
        vmmc_events[t_step_vmmc_event][n_vmmc_events[t_step_vmmc_event]] = indiv;
        n_vmmc_events[t_step_vmmc_event]++;
    
    }else{ // If VMMC event lies after end of the simulation.  
        
        /* If next event scheduled for after the end of the simulation set to be dummy entries. */
        indiv->idx_vmmc_event[0] = -1;
        indiv->idx_vmmc_event[1] = -1;
        if (indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("No VMMC event scheduled for %ld as event lies after end of the simulation.\n",
                indiv->id);
        }
    }
}


void carry_out_VMMC_events_per_timestep(int t_step, double t, patch_struct *patch, int p){
    /*Carry out any event associated with VMMC in the current time step
    
    
    Arguments
    ---------
    t_step : int
        Current time step (used to index the patch[p].vmmc_events and patch[p].n_vmmc_events)
    t : double
        Current time in years.  
    patch : pointer to an array of patch_struct structures
        The array of patch_struct objects that house information on patches.  See structures.h for 
        a list of attributes that these objects have.  
    p : int
        Patch identifier (generally 0 or 1).  
    
    Returns
    -------
    Nothing; carries out VMMC events on individuals for which they are scheduled.  
    */
    
    /* For debugging: */
    if(t_step < 0 || t_step >= N_TIME_STEP_PER_YEAR){
        printf("ERROR: array index %d for vmmc event out of bounds", t_step);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    int n_events = patch[p].n_vmmc_events[t_step];
    individual *indiv;
    int n;
    //printf("Carrying out %d VMMC events at time t=%f\n",n_events,t);
    
    for(n = 0; n < n_events; n++){
        indiv = patch[p].vmmc_events[t_step][n];
        //printf("Person %ld with circ=%d is in vmmc_events.\n",indiv->id,indiv->circ);
        
        /* Throw an error if this individual is female */
        if (indiv->gender == FEMALE){
            printf("ERROR: There is a woman %ld in vmmc_events. Exiting\n",indiv->id);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        /* If this individual is dead, move on to the next person.  
        Note - we can set up a similar procedure to other lists to remove this person from this
        list but it is not necessary. As things stand, no VMMC event happens to the dead person and
        no new event is scheduled for them. */
        if(indiv->cd4 == DEAD){
            continue;
        }
        /* If uncircumcised but waiting for VMMC then at this timestep they get circumcised. */
        if (indiv->circ == UNCIRC_WAITING_VMMC){
            //printf("Person %ld with circ=%d is being scheduled for VMMC healing.\n",
            //      indiv->id,indiv->circ);
            schedule_vmmc_healing(indiv, patch[p].param, patch[p].vmmc_events,
                patch[p].n_vmmc_events, patch[p].size_vmmc_events, t);
            
            // Count the number of VMMC procedures in the current year by counting the 
            // time at which the VMMC procedure was performed.  
            int year_idx = (int) floor(t) - patch[p].param->start_time_simul;
            patch[p].calendar_outputs->N_calendar_VMMC[year_idx]++;
            
        }else if (indiv->circ == VMMC_HEALING){
            /* If current status is healing, then finish healing. Note that this is the last event
            in the VMMC process for this individual. */
            
            //printf("Person %ld with circ=%d is being scheduled to move to VMMC .\n",
            //    indiv->id,indiv->circ);
            finish_vmmc_healing(indiv);
            //printf("Person %ld with circ=%d now VMMC .\n",indiv->id,indiv->circ);
        }else{
            printf("ERROR: not sure why this person %ld with circ=%d is in vmmc_events. Exiting\n",
                indiv->id,indiv->circ);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    
    /* At this point we have carried out all the events stored in vmmc_events[t_step].
    
    We reuse the same array next year, so need to set n_vmmc_events[] to be zero.  Note that we do
    not need to set the elements in vmmc_events[t_step] to be blank
    as we overwrite any elements we use, and n_vmmc_events[] prevents us from accessing elements
    from previous years which have not been overwritten already. */
    patch[p].n_vmmc_events[t_step] = 0;
}
