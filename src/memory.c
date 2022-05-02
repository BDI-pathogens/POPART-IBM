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

/*Functions:
 * void reinitialize_arrays_to_default() - at the start of a new run set all necessary pointers to zero.
 *
 */


/************************************************************************/
/******************************* Includes  ******************************/
/************************************************************************/

#include "memory.h"
#include "utilities.h"
#include "constants.h"
#include "pc.h"
/************************************************************************/
/******************************** functions *****************************/
/************************************************************************/

/* Blanks all the used individual_population entries (up to id_counter). Note that this function MUST be called before set_up_population
 * as set_up_population resets id_counter to 0.
 */
void blank_individual_array(individual *individual_population, int id_counter){
    /* This is a blank template to make it easier to debug when we accidentally add people from previous runs who should not exist in the current run. */

    int i_id;
    individual blank_person_template;

    blank_person_template.HIV_status = DUMMYVALUE;
    blank_person_template.ART_status = DUMMYVALUE;
    blank_person_template.cd4 = DUMMYVALUE;
    blank_person_template.next_HIV_event = DUMMYVALUE;
    blank_person_template.next_cascade_event = DUMMYVALUE;
    blank_person_template.idx_hiv_pos_progression[0] = DUMMYVALUE;
    blank_person_template.idx_hiv_pos_progression[1] = DUMMYVALUE;
    blank_person_template.debug_last_hiv_event_index = DUMMYVALUE;
    blank_person_template.idx_cascade_event[0] = DUMMYVALUE;
    blank_person_template.idx_cascade_event[1] = DUMMYVALUE;
    blank_person_template.debug_last_cascade_event_index = DUMMYVALUE;
    blank_person_template.time_to_delivery = DUMMYVALUE;
    blank_person_template.n_partners = DUMMYVALUE;
    blank_person_template.circ = DUMMYVALUE;
    blank_person_template.idx_vmmc_event[0] = DUMMYVALUE;
    blank_person_template.idx_vmmc_event[1] = DUMMYVALUE;
    blank_person_template.debug_last_vmmc_event_index = DUMMYVALUE;

    blank_person_template.NCHIPSVISITS = 0;
    blank_person_template.VISITEDBYCHIPS_TO_INIT_ART = FALSE;
    blank_person_template.VISITED_BY_CHIPS_THISROUND = FALSE;

    /* Now blank out each individual: */
    for (i_id=0;i_id<id_counter;i_id++)
        individual_population[i_id] = blank_person_template;

}

/*reinitialize_arrays_to_default(patch[p].death_dummylist, overall_partnerships->partner_dummylist, overall_partnerships->n_partnerships,
                            overall_partnerships->n_susceptible_in_serodiscordant_partnership, patch[p].n_hiv_pos_progression,
                            patch[p].size_hiv_pos_progression, patch[p].n_cascade_events, patch[p].size_cascade_events,
                            patch[p].n_vmmc_events, patch[p].size_vmmc_events, overall_partnerships->n_planned_breakups,
                            overall_partnerships->size_planned_breakups, output->annual_outputs_string[p], output->annual_outputs_string_pconly[p], output->annual_partnerships_outputs_string[p], output->annual_partnerships_outputs_string_pconly[p], output->timestep_outputs_string[p], output->timestep_outputs_string_PConly[p],
                            output->annual_outputs_string_prevalence[p], output->annual_outputs_string_knowserostatus[p],
                            output->annual_outputs_string_knowserostatusandonart[p], output->phylogenetics_output_string);*/

void reinitialize_arrays_to_default(int p, patch_struct *patch, all_partnerships *overall_partnerships, output_struct *output)
{

    long i;
    int g, ac, a, chips_round, pc_round;
    for (i=0; i<MAX_N_PER_AGE_GROUP; i++)
        patch[p].death_dummylist[i] = i;         /* Initialize the dummy list. */

    for (i=0; i<MAX_N_PER_AGE_GROUP*MAX_PARTNERSHIPS_PER_INDIVIDUAL; i++)
        overall_partnerships->partner_dummylist[i] = i;         /* Initialize the dummy list. */
    *overall_partnerships->n_partnerships = 0;
    *overall_partnerships->n_susceptible_in_serodiscordant_partnership = 0;
    /* Initialise the number of people in each group to be zero (as no HIV at start of simulation): */
    for (i=0; i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR; i++)
        patch[p].n_hiv_pos_progression[i] = 0;

    /* Initialise the  size of the arrays to the default: */
    for (i=0; i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR; i++)
        patch[p].size_hiv_pos_progression[i] = DEFAULT_N_HIV_PROGRESS_PER_TIME_STEP;

    /* Initialise the number of people in each group to be zero (as no HIV at start of simulation): */
    for (i=0; i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR; i++)
        patch[p].n_cascade_events[i] = 0;

    /* Initialise the size of the arrays to the default: */
    for (i=0; i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR; i++)
        patch[p].size_cascade_events[i] = DEFAULT_N_HIV_PROGRESS_PER_TIME_STEP;


    /* Initialise the number of people in each group to be zero 
     * (as no VMMC at start of simulation - note traditional MC is dealt with separately): */
    for (i=0; i<N_TIME_STEP_PER_YEAR; i++)
        patch[p].n_vmmc_events[i] = 0;

    /* Initialise the size of the arrays to the default: */
    for (i=0; i<N_TIME_STEP_PER_YEAR; i++)
        patch[p].size_vmmc_events[i] = DEFAULT_N_HIV_PROGRESS_PER_TIME_STEP;

    for (i=0; i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR; i++)
        overall_partnerships->n_planned_breakups[i] = 0;

    for (i=0; i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR; i++)
        overall_partnerships->size_planned_breakups[i] = MAX_BREAKUPS_PER_TIME_STEP;

    for (i=0; i<NPC_ENROLMENTS; i++)
        patch[p].PC_cohort_data->PC_cohort_counter[i] = 0;

    for (g=0; g<N_GENDER; g++){
        for (ac=0; ac<(MAX_AGE-AGE_CHIPS+1); ac++){
            for (chips_round=0; chips_round<NCHIPSROUNDS; chips_round++){
                output->NCHIPS_VISITED[p][g][ac][chips_round] = 0;
                output->NCHIPS_HIVPOS[p][g][ac][chips_round] = 0;
                output->NCHIPS_HIVAWARE[p][g][ac][chips_round] = 0;
                output->NCHIPS_ONART[p][g][ac][chips_round] = 0;
                output->NCHIPS_VS[p][g][ac][chips_round] = 0;
            }
        }
    }
    
    // Reset counters of new infections and person-years each PC round.  
    for(g = 0; g < N_GENDER; g++){
        for(a = 0; a < PC_AGE_RANGE_MAX; a++){
            for(pc_round = 0; pc_round < NPC_ROUNDS; pc_round++){
                output->PC_ROUND_INFECTIONS[p][g][a][pc_round] = 0;
                output->PC_ROUND_PERSON_TIMESTEPS[p][g][a][pc_round] = 0;
                
                output->PC_NPOP[p][g][a][pc_round] = 0;
                output->PC_NPOSITIVE[p][g][a][pc_round] = 0;
                output->PC_NAWARE[p][g][a][pc_round] = 0;
                output->PC_NONART[p][g][a][pc_round] = 0;
                output->PC_NVS[p][g][a][pc_round] = 0;
            }
        }
    }

    /* Set annual outputs strings to be blank: */
    memset(output->annual_outputs_string[p], '\0', SIZEOF_annual_outputs_string*sizeof(char));
    memset(output->annual_outputs_string_pconly[p], '\0', SIZEOF_annual_outputs_string_pconly*sizeof(char));
    memset(output->annual_partnerships_outputs_string[p], '\0', SIZEOF_annual_outputs_string*sizeof(char));
    memset(output->annual_partnerships_outputs_string_pconly[p], '\0', SIZEOF_annual_outputs_string_pconly*sizeof(char));
    memset(output->timestep_outputs_string[p], '\0', (N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string*sizeof(char));
    memset(output->timestep_outputs_string_PConly[p], '\0', (N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string*sizeof(char));
    memset(output->timestep_age_outputs_string[p], '\0', (N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string*sizeof(char));
    memset(output->timestep_age_outputs_string_PConly[p], '\0', (N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string*sizeof(char));
    memset(output->dhs_output_string[p], '\0', SIZEOF_annual_outputs_tempstore*sizeof(char));
    memset(output->pc_output_string[p], '\0', SIZEOF_calibration_outputs*sizeof(char));
    //memset(output->annual_outputs_string_prevalence[p], '\0', SIZEOF_annual_outputs_tempstore*sizeof(char));
    //memset(output->annual_outputs_string_knowserostatus[p], '\0', SIZEOF_annual_outputs_tempstore*sizeof(char));
    //memset(output->annual_outputs_string_knowserostatusandonart[p], '\0', SIZEOF_annual_outputs_tempstore*sizeof(char));
    /* Note we only blank calibration_outputs_combined_string every NRUNSPERWRITETOFILE runs - this is done in main.c at present. */
    memset(output->phylogenetics_output_string[p], '\0', PHYLO_OUTPUT_STRING_LENGTH*sizeof(char));
    memset(output->hazard_output_string, '\0', HAZARD_OUTPUT_STRING_LENGTH*sizeof(char));
    memset(output->cost_effectiveness_outputs_string[p], '\0',
        SIZEOF_cost_effectiveness_outputs_string*sizeof(char));
    memset(output->treats_outputs_string[p], '\0',
        (N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string*sizeof(char));
    memset(output->art_status_by_age_sex_outputs_string[p], '\0',
        (N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_calibration_outputs*sizeof(char));
}

//void alloc_pop_memory(population **pop, parameters **allrunparameters, int n_runs){
void alloc_pop_memory(population **pop, int n_runs){
    *pop = malloc(sizeof(population));
    if(*pop==NULL)
    {
        printf("Unable to allocate pop in alloc_pop_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

}

void alloc_output_memory(output_struct **output)
{
    *output = malloc(sizeof(output_struct));
    int p;
    (*output)->hazard_output_string = (char *)calloc(HAZARD_OUTPUT_STRING_LENGTH,sizeof(char));
    for (p=0;p<NPATCHES;p++){
        (*output)->phylogenetics_output_string[p] = (char *)calloc(PHYLO_OUTPUT_STRING_LENGTH,sizeof(char));
        (*output)->annual_outputs_string[p] = (char *)calloc(SIZEOF_annual_outputs_string, sizeof(char));
        (*output)->annual_outputs_string_pconly[p] = (char *)calloc(SIZEOF_annual_outputs_string_pconly, sizeof(char));
        (*output)->annual_partnerships_outputs_string[p] = (char *)calloc(SIZEOF_annual_outputs_string, sizeof(char));
        (*output)->annual_partnerships_outputs_string_pconly[p] = (char *)calloc(SIZEOF_annual_outputs_string_pconly, sizeof(char));
        (*output)->timestep_outputs_string[p] = (char *)calloc((N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string, sizeof(char));
        (*output)->timestep_outputs_string_PConly[p] = (char *)calloc((N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string, sizeof(char));
        (*output)->timestep_age_outputs_string[p] = (char *)calloc((N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string, sizeof(char));
        (*output)->timestep_age_outputs_string_PConly[p] = (char *)calloc((N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string, sizeof(char));
        (*output)->chips_output_string[p] = (char *)calloc(SIZEOF_annual_outputs_string, sizeof(char));


        //(*output)->annual_outputs_string_prevalence[p] = (char *)calloc(SIZEOF_annual_outputs_tempstore, sizeof(char));
        //(*output)->annual_outputs_string_knowserostatus[p] = (char *)calloc(SIZEOF_annual_outputs_tempstore, sizeof(char));
        //(*output)->annual_outputs_string_knowserostatusandonart[p] = (char *)calloc(SIZEOF_annual_outputs_tempstore, sizeof(char));
        (*output)->dhs_output_string[p] = 
            (char *)calloc(SIZEOF_annual_outputs_tempstore, sizeof(char));
        (*output)->pc_output_string[p] = 
            (char *)calloc(SIZEOF_calibration_outputs, sizeof(char));
        (*output)->calibration_outputs_combined_string[p] = 
            (char *)calloc(SIZEOF_calibration_outputs, sizeof(char));
        (*output)->cost_effectiveness_outputs_string[p] = 
            (char *)calloc(SIZEOF_cost_effectiveness_outputs_string, sizeof(char));
        (*output)->treats_outputs_string[p] = (char *)calloc((N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_annual_outputs_string, sizeof(char));
        (*output)->art_status_by_age_sex_outputs_string[p] = (char *)calloc((N_TIME_STEP_PER_YEAR/OUTPUTTIMESTEP)*SIZEOF_calibration_outputs, sizeof(char));
    }
}

/* Allows us to blank cohort data for PC0, PC12N or PC24N at start of a run. */
void set_to_null_pc_cohort_data(patch_struct *patch, int p, int pc_enrolment_round, int PC_round_population_size){
    int i,pc_round;



    /* This is the counter for who we visit next during visits for each round. Note that we also need to reset between rounds. */
    reset_visit_counter(patch,  p);

    int g, ap, i_pc_category; /* Indices splitting up the cohort by gender, age, HIV status etc. */
    /* Reset the counter for visits: */
    for (g=0;g<N_GENDER;g++){
        /* Run from 18 to 44 (inclusive). */
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
//
//              patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category] = 0;
                /* This is the counter for how many people we need to visit. */
                patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category] = 0;
                for(i=0; i<MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP; i++)
                    patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i] = -1;
            }
        }
    }


    //printf("Setting PC cohort data counters to zero for round %i patch %i\n",pc_enrolment_round,p);
    if (pc_enrolment_round==0){
        //(*PC_cohort_data)->PC0_cohort_data = malloc(PC_round_population_size*sizeof(PC_cohort_individual_data_struct));
        for(i=0; i<PC_round_population_size; i++){
            patch[p].PC_cohort_data->PC0_cohort_data[i].gender = -1;
            patch[p].PC_cohort_data->PC0_cohort_data[i].ap = -1;
            for (pc_round=0; pc_round<NPC_ROUNDS; pc_round++){
                patch[p].PC_cohort_data->PC0_cohort_data[i].RETAINED_IN_COHORT[pc_round] = -1;
                patch[p].PC_cohort_data->PC0_cohort_data[i].HIV_status[pc_round] = -1;
                patch[p].PC_cohort_data->PC0_cohort_data[i].ART_status[pc_round] = -1;
                patch[p].PC_cohort_data->PC0_cohort_data[i].PC_visit_dates[pc_round] = -1;
            }
        }
    }
    else if(pc_enrolment_round==1){
        //PC_cohort_data->PC12N_cohort_data = malloc(PC_round_population_size*sizeof(PC_cohort_individual_data_struct));
        for(i=0; i<PC_round_population_size; i++){
            patch[p].PC_cohort_data->PC12N_cohort_data[i].gender = -1;
            patch[p].PC_cohort_data->PC12N_cohort_data[i].ap = -1;
            for (pc_round=0; pc_round<NPC_ROUNDS; pc_round++){
                patch[p].PC_cohort_data->PC12N_cohort_data[i].RETAINED_IN_COHORT[pc_round] = -1;
                patch[p].PC_cohort_data->PC12N_cohort_data[i].HIV_status[pc_round] = -1;
                patch[p].PC_cohort_data->PC12N_cohort_data[i].ART_status[pc_round] = -1;
                patch[p].PC_cohort_data->PC12N_cohort_data[i].PC_visit_dates[pc_round] = -1;
            }
        }
    }
    else if(pc_enrolment_round==2){
        //PC_cohort_data->PC24N_cohort_data = malloc(PC_round_population_size*sizeof(PC_cohort_individual_data_struct));
        for(i=0; i<PC_round_population_size; i++){
            patch[p].PC_cohort_data->PC24N_cohort_data[i].gender = -1;
            patch[p].PC_cohort_data->PC24N_cohort_data[i].ap = -1;
            for (pc_round=0; pc_round<NPC_ROUNDS; pc_round++){
                patch[p].PC_cohort_data->PC24N_cohort_data[i].RETAINED_IN_COHORT[pc_round] = -1;
                patch[p].PC_cohort_data->PC24N_cohort_data[i].HIV_status[pc_round] = -1;
                patch[p].PC_cohort_data->PC24N_cohort_data[i].ART_status[pc_round] = -1;
                patch[p].PC_cohort_data->PC24N_cohort_data[i].PC_visit_dates[pc_round] = -1;
            }
        }
    }
    else{
        printf("ERROR: Unknown PC round number. Exiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

}


void alloc_pc_cohort_data(PC_cohort_data_struct **PC_cohort_data, int pc_enrolment_round, int PC_round_population_size){
    int i,pc_round;


    //printf("Setting PC cohort data counters to zero for round %i \n",pc_enrolment_round);
    if (pc_enrolment_round==0){
        (*PC_cohort_data)->PC0_cohort_data = malloc(PC_round_population_size*sizeof(PC_cohort_individual_data_struct));
        for(i=0; i<PC_round_population_size; i++){
            (*PC_cohort_data)->PC0_cohort_data[i].gender = -1;
            (*PC_cohort_data)->PC0_cohort_data[i].ap = -1;
            for (pc_round=0; pc_round<NPC_ROUNDS; pc_round++){
                (*PC_cohort_data)->PC0_cohort_data[i].RETAINED_IN_COHORT[pc_round] = -1;
                (*PC_cohort_data)->PC0_cohort_data[i].HIV_status[pc_round] = -1;
                (*PC_cohort_data)->PC0_cohort_data[i].ART_status[pc_round] = -1;
                (*PC_cohort_data)->PC0_cohort_data[i].PC_visit_dates[pc_round] = -1;
            }
        }
    }
    else if(pc_enrolment_round==1){
        (*PC_cohort_data)->PC12N_cohort_data = malloc(PC_round_population_size*sizeof(PC_cohort_individual_data_struct));
        for(i=0; i<PC_round_population_size; i++){
            (*PC_cohort_data)->PC12N_cohort_data[i].gender = -1;
            (*PC_cohort_data)->PC12N_cohort_data[i].ap = -1;
            for (pc_round=0; pc_round<NPC_ROUNDS; pc_round++){
                (*PC_cohort_data)->PC12N_cohort_data[i].RETAINED_IN_COHORT[pc_round] = -1;
                (*PC_cohort_data)->PC12N_cohort_data[i].HIV_status[pc_round] = -1;
                (*PC_cohort_data)->PC12N_cohort_data[i].ART_status[pc_round] = -1;
                (*PC_cohort_data)->PC12N_cohort_data[i].PC_visit_dates[pc_round] = -1;
            }
        }
    }
    else if(pc_enrolment_round==2){
        (*PC_cohort_data)->PC24N_cohort_data = malloc(PC_round_population_size*sizeof(PC_cohort_individual_data_struct));
        for(i=0; i<PC_round_population_size; i++){
            (*PC_cohort_data)->PC24N_cohort_data[i].gender = -1;
            (*PC_cohort_data)->PC24N_cohort_data[i].ap = -1;
            for (pc_round=0; pc_round<NPC_ROUNDS; pc_round++){
                (*PC_cohort_data)->PC24N_cohort_data[i].RETAINED_IN_COHORT[pc_round] = -1;
                (*PC_cohort_data)->PC24N_cohort_data[i].HIV_status[pc_round] = -1;
                (*PC_cohort_data)->PC24N_cohort_data[i].ART_status[pc_round] = -1;
                (*PC_cohort_data)->PC24N_cohort_data[i].PC_visit_dates[pc_round] = -1;
            }
        }
    }
    else{
        printf("ERROR: Unknown PC round number. Exiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

}


void alloc_patch_memoryv2(patch_struct *patch){
    int p,g;
    /*
     *patch = malloc(NPATCHES*sizeof(patch_struct));
    for (p=0; p<NPATCHES; p++){
        patch[p] = malloc(sizeof(patch_struct));

        if(patch[p]==NULL){
            printf("Unable to allocate patch in alloc_patch_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
     */

    long i;
    int a, r;


    for (p=0; p<NPATCHES; p++){

        /* WE DO NOT NEED TO ALLOCATE MEMORY TO PARAM AS WE JUST USE *param AS A POINTER TO THE MEMORY ALLOCATED FOR allparameters. */
        //      patch[p].param = malloc(sizeof(parameters));
        //      if(patch[p].param==NULL)
        //      {
        //          printf("Unable to allocate param in alloc_all_memory. Execution aborted.");
        //          printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        //          fflush(stdout);
        //          exit(1);
        //      }

        patch[p].individual_population = malloc(MAX_POP_SIZE*sizeof(individual));  // at this stage individual_population[i].partner_pairs and partner_pairs_HIVpos point to random place where no space has been allocated for storing these partnership objects. This is done later on
        if(patch[p].individual_population==NULL)
        {
            printf("Unable to allocate individual_population in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }


        patch[p].n_population = malloc(sizeof(population_size));
        if(patch[p].n_population==NULL)
        {
            printf("Unable to allocate n_population in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].n_population_oneyearagegroups = malloc(sizeof(population_size_one_year_age));
        if(patch[p].n_population_oneyearagegroups==NULL)
        {
            printf("Unable to allocate n_population_oneyearagegroups in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].n_population_stratified = malloc(sizeof(stratified_population_size));
        if(patch[p].n_population_stratified==NULL)
        {
            printf("Unable to allocate n_population_stratified in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].n_infected = malloc(sizeof(population_size_one_year_age));
        if(patch[p].n_infected==NULL)
        {
            printf("Unable to allocate n_infected in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].n_newly_infected = malloc(sizeof(population_size_one_year_age));
        if(patch[p].n_newly_infected==NULL)
        {
            printf("Unable to allocate n_newly_infected in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].n_infected_cumulative = malloc(sizeof(population_size_one_year_age));
        if(patch[p].n_infected_cumulative==NULL)
        {
            printf("Unable to allocate n_infected_cumulative in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].n_infected_wide_age_group = malloc(sizeof(population_size));
        if(patch[p].n_infected_wide_age_group==NULL)
        {
            printf("Unable to allocate n_infected_wide_age_group in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].n_newly_infected_wide_age_group = malloc(sizeof(population_size));
        if(patch[p].n_newly_infected_wide_age_group==NULL)
        {
            printf("Unable to allocate n_newly_infected_wide_age_group in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].age_list = malloc(sizeof(age_list_struct));
        if(patch[p].age_list==NULL)
        {
            printf("Unable to allocate age_list in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        for (g=0;g<N_GENDER;g++){
            patch[p].age_list->age_list_by_gender[g] = malloc(sizeof(age_list_element));
            if(patch[p].age_list->age_list_by_gender[g]==NULL){
                printf("Unable to allocate patch[p].age_list->age_list_by_gender[%i] in alloc_all_memory. Execution aborted.",g);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }

        //
        //
        //
        //      patch[p].age_list->age_list_m = malloc(sizeof(age_list_element));
        //      if(patch[p].age_list->age_list_m==NULL)
        //      {
        //          printf("Unable to allocate patch[p].age_list->age_list_m in alloc_all_memory. Execution aborted.");
        //          printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        //          fflush(stdout);
        //          exit(1);
        //      }


        patch[p].child_population = malloc(2*sizeof(child_population_struct));  // The 2 is because we have HIV- and HIV+ lists.
        if(patch[p].child_population==NULL)
        {
            printf("Unable to allocate child_population in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].new_deaths = malloc(MAX_N_PER_AGE_GROUP*sizeof(long));  // There are <MAX_N_PER_AGE_GROUP people in an age group, so the number of people dying in that group has to be less than that number.
        if(patch[p].new_deaths==NULL)
        {
            printf("Unable to allocate new_deaths in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].death_dummylist = malloc(MAX_N_PER_AGE_GROUP*sizeof(long));  // This is just a dummy list of 0..MAX_N_PER_AGE_GROUP from which we pick the indices of people who will die in deaths_natural_causes() in demographics.c. This list of indices is stored in new_death[].
        if(patch[p].death_dummylist==NULL)
        {
            printf("Unable to allocate death_dummylist in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }


        patch[p].hiv_pos_progression = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(individual**));  // hiv_pos_progression[t] is a list of pointers to HIV+ individuals whose next progression event will happen at time step param.start_time_simul+t
        if(patch[p].hiv_pos_progression==NULL)
        {
            printf("Unable to allocate hiv_pos_progression in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        for(i=0 ; i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR ; i++){
            (patch[p].hiv_pos_progression)[i] = malloc(DEFAULT_N_HIV_PROGRESS_PER_TIME_STEP*sizeof(individual*));
            if((patch[p].hiv_pos_progression)[i]==NULL){
                printf("Unable to allocate hiv_pos_progression[i] in alloc_all_memory. Execution aborted.");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }

        patch[p].n_hiv_pos_progression = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(long));
        if(patch[p].n_hiv_pos_progression==NULL)
        {
            printf("Unable to allocate n_hiv_pos_progression in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].size_hiv_pos_progression = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(long));
        if(patch[p].size_hiv_pos_progression==NULL){
            printf("Unable to allocate size_hiv_pos_progression in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }


        /* These is the arrays that store the cascade events (HIV testing, starting/stopping ART): */
        patch[p].cascade_events = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(individual**));  // cascade_events[t] is a list of pointers to individuals whose next cascade event will happen at time step param.start_time_simul+t
        if(patch[p].cascade_events==NULL)
        {
            printf("Unable to allocate cascade_events in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        for(i=0 ; i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR ; i++){
            (patch[p].cascade_events)[i] = malloc(DEFAULT_N_HIV_CASCADE_PER_TIME_STEP*sizeof(individual*));
            if((patch[p].cascade_events)[i]==NULL){
                printf("Unable to allocate cascade_events[i] in alloc_all_memory. Execution aborted.");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }

        patch[p].n_cascade_events = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(long));
        if(patch[p].n_cascade_events==NULL)
        {
            printf("Unable to allocate n_cascade_events in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].size_cascade_events = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(long));
        if(patch[p].size_cascade_events==NULL){
            printf("Unable to allocate size_cascade_events in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        //////
        /* This schedules VMMC events. Note that unlike previous (HIV progression and cascade) we assume that VMMC events (scheduling VMMC, healing period) all take <1 year (very conservative - could make a few weeks if memory issues). */
        patch[p].vmmc_events = malloc(N_TIME_STEP_PER_YEAR*sizeof(individual**));  // vmmc_events[t] is a list of pointers to HIV+ individuals whose next progression event will happen at time step param->start_time_simul+t
        if(patch[p].vmmc_events==NULL)
        {
            printf("Unable to allocate vmmc_events in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        for(i=0 ; i<N_TIME_STEP_PER_YEAR ; i++){
            (patch[p].vmmc_events)[i] = malloc(DEFAULT_N_HIV_PROGRESS_PER_TIME_STEP*sizeof(individual*));
            if((patch[p].vmmc_events)[i]==NULL){
                printf("Unable to allocate vmmc_events[i] in alloc_all_memory. Execution aborted.");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }

        patch[p].n_vmmc_events = malloc(N_TIME_STEP_PER_YEAR*sizeof(long));
        if(patch[p].n_vmmc_events==NULL)
        {
            printf("Unable to allocate n_vmmc_events in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].size_vmmc_events = malloc(N_TIME_STEP_PER_YEAR*sizeof(long));
        if(patch[p].size_vmmc_events==NULL){
            printf("Unable to allocate size_vmmc_events in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }


        patch[p].chips_sample = malloc(sizeof(chips_sample_struct));
        if(patch[p].chips_sample==NULL)
        {
            printf("Unable to allocate chips_sample in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].PC_sample = malloc(sizeof(PC_sample_struct));
        if(patch[p].PC_sample==NULL)
        {
            printf("Unable to allocate PC_sample in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].PC_cohort = malloc(sizeof(PC_cohort_struct));
        if(patch[p].PC_cohort==NULL)
        {
            printf("Unable to allocate PC_cohort in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].PC_cohort_data = malloc(sizeof(PC_cohort_data_struct));
        if(patch[p].PC_cohort_data==NULL)
        {
            printf("Unable to allocate PC_cohort_data in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        patch[p].cumulative_outputs  = malloc(sizeof(chips_sample_struct));
        if(patch[p].cumulative_outputs==NULL)
        {
            printf("Unable to allocate cumulative_outputs in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        patch[p].calendar_outputs  = malloc(sizeof(calendar_outputs_struct));
        if(patch[p].calendar_outputs==NULL)
        {
            printf("Unable to allocate calendar_outputs in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }


        //      patch[p].i_fit = malloc(sizeof(int));
        //      if(patch[p].i_fit==NULL)
        //      {
        //          printf("Unable to allocate i_fit in alloc_all_memory. Execution aborted.");
        //          printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        //          fflush(stdout);
        //          exit(1);
        //      }


        patch[p].cross_sectional_distr_n_lifetime_partners = (long ****) malloc(N_GENDER*sizeof(long***));
        if(patch[p].cross_sectional_distr_n_lifetime_partners==NULL)
        {
            printf("Unable to allocate cross_sectional_distr_n_lifetime_partners in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        for(g=0 ; g<N_GENDER ; g++)
        {
            patch[p].cross_sectional_distr_n_lifetime_partners[g] = (long ***) malloc(N_AGE*sizeof(long**));
            if(patch[p].cross_sectional_distr_n_lifetime_partners[g]==NULL)
            {
                printf("Unable to allocate cross_sectional_distr_n_lifetime_partners[g] in alloc_all_memory. Execution aborted.");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            for(a=0 ; a<N_AGE ; a++)
            {
                patch[p].cross_sectional_distr_n_lifetime_partners[g][a] = (long **) malloc(N_RISK*sizeof(long*));
                if(patch[p].cross_sectional_distr_n_lifetime_partners[g][a]==NULL)
                {
                    printf("Unable to allocate cross_sectional_distr_n_lifetime_partners[g][a] in alloc_all_memory. Execution aborted.");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                for(r=0 ; r<N_RISK ; r++)
                {
                    patch[p].cross_sectional_distr_n_lifetime_partners[g][a][r] = (long *) malloc((MAX_N_PARTNERS_IN_OUTPUTS+1)*sizeof(long));
                    if(patch[p].cross_sectional_distr_n_lifetime_partners[g][a][r]==NULL)
                    {
                        printf("Unable to allocate cross_sectional_distr_n_lifetime_partners[g][a][r] in alloc_all_memory. Execution aborted.");
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                }
            }
        }

        patch[p].cross_sectional_distr_n_partners_lastyear = (long ****) malloc(N_GENDER*sizeof(long***));
        if(patch[p].cross_sectional_distr_n_partners_lastyear==NULL)
        {
            printf("Unable to allocate cross_sectional_distr_n_partners_lastyear in alloc_all_memory. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        for(g=0 ; g<N_GENDER ; g++)
        {
            patch[p].cross_sectional_distr_n_partners_lastyear[g] = (long ***) malloc(N_AGE*sizeof(long**));
            if(patch[p].cross_sectional_distr_n_partners_lastyear[g]==NULL)
            {
                printf("Unable to allocate cross_sectional_distr_n_partners_lastyear[g] in alloc_all_memory. Execution aborted.");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            for(a=0 ; a<N_AGE ; a++)
            {
                patch[p].cross_sectional_distr_n_partners_lastyear[g][a] = (long **) malloc(N_RISK*sizeof(long*));
                if(patch[p].cross_sectional_distr_n_partners_lastyear[g][a]==NULL)
                {
                    printf("Unable to allocate cross_sectional_distr_n_partners_lastyear[g][a] in alloc_all_memory. Execution aborted.");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                for(r=0 ; r<N_RISK ; r++)
                {
                    patch[p].cross_sectional_distr_n_partners_lastyear[g][a][r] = (long *) malloc((MAX_N_PARTNERS_IN_OUTPUTS+1)*sizeof(long));
                    if(patch[p].cross_sectional_distr_n_partners_lastyear[g][a][r]==NULL)
                    {
                        printf("Unable to allocate cross_sectional_distr_n_partners_lastyear[g][a][r] in alloc_all_memory. Execution aborted.");
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                }
            }
        }

    }
}

void alloc_partnership_memoryv2(all_partnerships *overall_partnerships){

    int i,j,k,l;
    long m;



    overall_partnerships->new_partners_f_sorted = malloc(MAX_N_PER_AGE_GROUP*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(long));
    if(overall_partnerships->new_partners_f_sorted==NULL)
    {
        printf("Unable to allocate new_partners_f_sorted in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    overall_partnerships->shuffled_idx = malloc(MAX_N_PER_AGE_GROUP*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(long));
    if(overall_partnerships->shuffled_idx==NULL)
    {
        printf("Unable to allocate shuffled_idx in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    overall_partnerships->new_partners_f_non_matchable = malloc(MAX_N_PER_AGE_GROUP*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(long));
    if(overall_partnerships->new_partners_f_non_matchable==NULL)
    {
        printf("Unable to allocate new_partners_f_non_matchable in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    overall_partnerships->new_partners_m = malloc(MAX_N_PER_AGE_GROUP*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(long));
    if(overall_partnerships->new_partners_m==NULL)
    {
        printf("Unable to allocate new_partners_m in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    overall_partnerships->new_partners_m_sorted = malloc(MAX_N_PER_AGE_GROUP*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(long));
    if(overall_partnerships->new_partners_m_sorted==NULL)
    {
        printf("Unable to allocate new_partners_m_sorted in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    overall_partnerships->partner_dummylist = malloc(MAX_N_PER_AGE_GROUP*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(long));
    if(overall_partnerships->partner_dummylist==NULL)
    {
        printf("Unable to allocate partner_dummylist in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }


    // THIS IS ***FAR*** TOO SMALL: should be something like (Duration of simulation) / ((Average age at death-AGE_ADULT)/(average # of partners)) * (MAX_POP_SIZE/2)
    overall_partnerships->partner_pairs = malloc(MAX_POP_SIZE*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(partnership)); // MAYBE NEEDS SOMETHING BIGGER BECAUSE OF DYNAMICS OF PARTNERSHIPS
    if(overall_partnerships->partner_pairs==NULL)
    {
        printf("Unable to allocate partner_pairs in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    overall_partnerships->n_partnerships = malloc(sizeof(long));  // only need space for one long
    if(overall_partnerships->n_partnerships==NULL)
    {
        printf("Unable to allocate n_partnerships in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    overall_partnerships->susceptible_in_serodiscordant_partnership = malloc(MAX_POP_SIZE*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(individual*));
    if(overall_partnerships->susceptible_in_serodiscordant_partnership==NULL)
    {
        printf("Unable to allocate susceptible_in_serodiscordant_partnership in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    overall_partnerships->n_susceptible_in_serodiscordant_partnership = malloc(sizeof(long));  // only need space for one long
    if(overall_partnerships->n_susceptible_in_serodiscordant_partnership==NULL)
    {
        printf("Unable to allocate n_susceptible_in_serodiscordant_partnership in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }


    overall_partnerships->pop_available_partners = malloc(sizeof(population_partners)); // Here only storing adresses to individuals who already exist and have been allocated memory for so should be ok
    if(overall_partnerships->pop_available_partners==NULL)
    {
        printf("Unable to allocate overall_partnerships->pop_available_partners in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /****************************************************************************************************/
    /**** Code to dynamically allocate memory for pop_available_partners->pop_per_patch_gender_age_risk. */
    /****************************************************************************************************/

    overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk = malloc(NPATCHES*sizeof(individual *****));
    if(overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk==NULL){
        printf("Unable to allocate memory for overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk.\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    for (i=0; i<NPATCHES; i++){
        overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i] = malloc(N_GENDER*sizeof(individual ****));
        if(overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i]==NULL){
            printf("Unable to allocate memory for overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i].\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        for (j=0; j<N_GENDER; j++){
            overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j] = malloc(N_AGE*sizeof(individual ***));
            if(overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j]==NULL){
                printf("Unable to allocate memory for overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j].\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            for (k=0; k<N_AGE; k++){
                overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j][k] = malloc(N_RISK*sizeof(individual **));
                if(overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j][k]==NULL){
                    printf("Unable to allocate memory for overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j][k].\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                for (l=0; l<N_RISK; l++){
                    overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j][k][l] = malloc(MAX_N_PER_AGE_GROUP*MAX_PARTNERSHIPS_PER_INDIVIDUAL*sizeof(individual *));
                    if(overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j][k]==NULL){
                        printf("Unable to allocate memory for overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[i][j][k][l].\n");
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                }
            }
        }
    }

    /****************************************************************************************************/
    /****************************************************************************************************/
    /****************************************************************************************************/






    overall_partnerships->n_pop_available_partners = malloc(sizeof(population_size_all_patches)); // This is ok as population_size contains only objects allocated statically
    if(overall_partnerships->n_pop_available_partners==NULL)
    {
        printf("Unable to allocate n_pop_available_partners in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }



    overall_partnerships->planned_breakups = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(partnership**)); // planned_breakups[t] is a list of pointers to partnerships planned to break up at time step param->start_time_simul+t
    if(overall_partnerships->planned_breakups==NULL)
    {
        printf("Unable to allocate planned_breakups in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    for(m=0 ; m<MAX_N_YEARS*N_TIME_STEP_PER_YEAR ; m++)
    {
        (overall_partnerships->planned_breakups)[m] = malloc(MAX_BREAKUPS_PER_TIME_STEP*sizeof(partnership*));
    }

    overall_partnerships->n_planned_breakups = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(long));
    if(overall_partnerships->n_planned_breakups==NULL)
    {
        printf("Unable to allocate n_planned_breakups in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }


    overall_partnerships->size_planned_breakups = malloc(MAX_N_YEARS*N_TIME_STEP_PER_YEAR*sizeof(long));
    if(overall_partnerships->size_planned_breakups==NULL)
    {
        printf("Unable to allocate size_planned_breakups in alloc_all_memory. Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

}


void free_all_patch_memory(parameters *param, individual *individual_population, population_size *n_population, population_size_one_year_age *n_population_oneyearagegroups, stratified_population_size *n_population_stratified, age_list_struct *age_list, child_population_struct *child_population,
        individual ***hiv_pos_progression, long *n_hiv_pos_progression, long *size_hiv_pos_progression, individual ***cascade_events, long *n_cascade_events, long *size_cascade_events, individual ***vmmc_events, long *n_vmmc_events, long *size_vmmc_events,
        long *new_deaths, long *death_dummylist,
        population_size_one_year_age *n_infected, population_size_one_year_age *n_newly_infected, population_size_one_year_age *n_infected_cumulative, population_size *n_infected_wide_age_group, population_size *n_newly_infected_wide_age_group,
        chips_sample_struct *chips_sample, cumulative_outputs_struct *cumulative_outputs, calendar_outputs_struct *calendar_outputs, long ****cross_sectional_distr_n_lifetime_partners, long ****cross_sectional_distr_n_partners_lastyear, PC_sample_struct *PC_sample, PC_cohort_struct *PC_cohort, PC_cohort_data_struct *PC_cohort_data)
{

    long i;
    int g, a, r;
    free(individual_population);
    free(n_population);
    free(n_population_oneyearagegroups);
    free(n_population_stratified);

    free(n_infected);
    free(n_newly_infected);
    free(n_infected_cumulative);
    free(n_infected_wide_age_group);
    free(n_newly_infected_wide_age_group);
    for (g=0;g<N_GENDER;g++)
        free(age_list->age_list_by_gender[g]);
    //  free(age_list->age_list_m);
    free(age_list);
    free(child_population);

    free(new_deaths);
    free(death_dummylist);

    /* Note we free each element within hiv_pos_progression in the function carry_out_HIV_events_per_timestep() in hiv.c. */
    for(i=0 ;i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR ;i++)
        free(hiv_pos_progression[i]);
    free(hiv_pos_progression);
    free(n_hiv_pos_progression);
    free(size_hiv_pos_progression);

    /* Note we free each element within cascade_events in the function carry_out_cascade_events_per_timestep() in hiv.c. */
    for(i=0 ;i<MAX_N_YEARS*N_TIME_STEP_PER_YEAR ;i++)
        free(cascade_events[i]);
    free(cascade_events);
    free(n_cascade_events);
    free(size_cascade_events);


    /* Didn't free vmmc events so free them here: */
    for(i=0; i<N_TIME_STEP_PER_YEAR; i++){
        free(vmmc_events[i]);
    }
    free(vmmc_events);
    free(n_vmmc_events);
    free(size_vmmc_events);

    free(chips_sample);
    free(PC_sample);
    free(PC_cohort);


    free_pc_cohort_data_memory(PC_cohort_data);  /* Firstly free memory inside: */
    free(PC_cohort_data);

    free(cumulative_outputs);
    free(calendar_outputs);
    
    for (g=0;g<N_GENDER;g++)
    {
        for (a=0;a<N_AGE;a++)
        {
            for (r=0;r<N_RISK;r++)
            {
                free(cross_sectional_distr_n_lifetime_partners[g][a][r]);
                free(cross_sectional_distr_n_partners_lastyear[g][a][r]);
            }
            free(cross_sectional_distr_n_lifetime_partners[g][a]);
            free(cross_sectional_distr_n_partners_lastyear[g][a]);

        }
        free(cross_sectional_distr_n_lifetime_partners[g]);
        free(cross_sectional_distr_n_partners_lastyear[g]);
    }
    free(cross_sectional_distr_n_lifetime_partners);
    free(cross_sectional_distr_n_partners_lastyear);


}

void free_all_partnership_memory(partnership *partner_pairs, long *n_partnerships, individual **susceptible_in_serodiscordant_partnership, long *n_susceptible_in_serodiscordant_partnership,
        population_partners* pop_available_partners, population_size_all_patches *n_pop_available_partners,
        partnership ***planned_breakups, long *n_planned_breakups, long *size_planned_breakups,
        long *new_partners_f_sorted, long *shuffled_idx, long *new_partners_f_non_matchable, long *new_partners_m, long *new_partners_m_sorted, long *partner_dummylist)
{
    int i,j,k,l;
    long m;

    /****************************************************************************************************/
    /**** Code to dynamically allocate memory for pop_available_partners->pop_per_patch_gender_age_risk. */
    /****************************************************************************************************/




    for (i=0; i<NPATCHES; i++){
        for (j=0; j<N_GENDER; j++){
            for (k=0; k<N_AGE; k++){
                for (l=0; l<N_RISK; l++)
                    free(pop_available_partners->pop_per_patch_gender_age_risk[i][j][k][l]);
                free(pop_available_partners->pop_per_patch_gender_age_risk[i][j][k]);
            }
            free(pop_available_partners->pop_per_patch_gender_age_risk[i][j]);
        }
        free(pop_available_partners->pop_per_patch_gender_age_risk[i]);
    }
    free(pop_available_partners->pop_per_patch_gender_age_risk);
    free(pop_available_partners);

    free(new_partners_f_sorted);
    free(shuffled_idx);
    free(new_partners_f_non_matchable);
    free(new_partners_m);
    free(new_partners_m_sorted);
    free(partner_dummylist);

    free(partner_pairs);
    free(n_partnerships);
    free(susceptible_in_serodiscordant_partnership);
    free(n_susceptible_in_serodiscordant_partnership);

    free(n_pop_available_partners);
    for(m=0 ; m<MAX_N_YEARS*N_TIME_STEP_PER_YEAR ; m++)
    {
        free(planned_breakups[m]);
    }
    free(planned_breakups);
    free(n_planned_breakups);
    free(size_planned_breakups);


}


void free_patch_memory(patch_struct *patch){
    int p;

    for (p=0; p<NPATCHES; p++){

        free_all_patch_memory(patch[p].param, patch[p].individual_population,
                patch[p].n_population, patch[p].n_population_oneyearagegroups, patch[p].n_population_stratified,
                patch[p].age_list, patch[p].child_population,
                patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression,
                patch[p].size_hiv_pos_progression, patch[p].cascade_events, patch[p].n_cascade_events,
                patch[p].size_cascade_events, patch[p].vmmc_events, patch[p].n_vmmc_events,
                patch[p].size_vmmc_events, patch[p].new_deaths, patch[p].death_dummylist,
                patch[p].n_infected,
                patch[p].n_newly_infected, patch[p].n_infected_cumulative, patch[p].n_infected_wide_age_group, patch[p].n_newly_infected_wide_age_group,
                patch[p].chips_sample, patch[p].cumulative_outputs, patch[p].calendar_outputs, patch[p].cross_sectional_distr_n_lifetime_partners, patch[p].cross_sectional_distr_n_partners_lastyear,
                patch[p].PC_sample, patch[p].PC_cohort, patch[p].PC_cohort_data);
    }
    free(patch);
}

void free_partnership_memory(all_partnerships *overall_partnerships){

    free_all_partnership_memory(overall_partnerships->partner_pairs,overall_partnerships->n_partnerships,
            overall_partnerships->susceptible_in_serodiscordant_partnership, overall_partnerships->n_susceptible_in_serodiscordant_partnership,
            overall_partnerships->pop_available_partners, overall_partnerships->n_pop_available_partners,
            overall_partnerships->planned_breakups, overall_partnerships->n_planned_breakups,
            overall_partnerships->size_planned_breakups,
            overall_partnerships->new_partners_f_sorted, overall_partnerships->shuffled_idx, overall_partnerships->new_partners_f_non_matchable,
            overall_partnerships->new_partners_m, overall_partnerships->new_partners_m_sorted, overall_partnerships->partner_dummylist);

}

void free_pop_memory(population *pop, parameters **allrunparameters){
    free(pop);
}

void allocate_fitting_data_memory(int n_fit, fitting_data_struct **fitting_data){
    *fitting_data = malloc((n_fit+1)*sizeof(fitting_data_struct));
    if(*fitting_data==NULL){
        printf("Unable to allocate fitting_data in allocate_fitting_data_memory(). Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}


void free_fitting_data_memory(fitting_data_struct *fitting_data){
    free(fitting_data);
}

void free_output_memory(output_struct *output){
    int p;
    for (p=0; p<NPATCHES; p++){

        free(output->annual_outputs_string[p]);
        free(output->annual_outputs_string_pconly[p]);
        free(output->annual_partnerships_outputs_string[p]);
        free(output->annual_partnerships_outputs_string_pconly[p]);
        free(output->timestep_outputs_string[p]);
        free(output->timestep_outputs_string_PConly[p]);
        free(output->timestep_age_outputs_string[p]);
        free(output->timestep_age_outputs_string_PConly[p]);
        free(output->chips_output_string[p]);
        free(output->dhs_output_string[p]);
        free(output->pc_output_string[p]);
        free(output->calibration_outputs_combined_string[p]);
        free(output->cost_effectiveness_outputs_string[p]);
        free(output->treats_outputs_string[p]);
        free(output->art_status_by_age_sex_outputs_string[p]);
        free(output->phylogenetics_output_string[p]);
    }
    free(output->hazard_output_string);

    free(output);
}

void free_pc_cohort_data_memory(PC_cohort_data_struct *PC_cohort_data){
    free(PC_cohort_data->PC0_cohort_data);
    if (NPC_ENROLMENTS>1)
        free(PC_cohort_data->PC12N_cohort_data);
    if (NPC_ENROLMENTS>2)
        free(PC_cohort_data->PC24N_cohort_data);
}
