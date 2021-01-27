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

/* Contains functions relevant to PC sampling and visits.
 * 1) Create PC list, and PC reserve list (of 5 people  of each age) for PC0 (note that PC12N and PC24N can draw from reserve list apart from new age 18).
 * 2) Add PC variables to indiv structure -  round of entry to PC (-1 not in PC, 0 = PC0, 1=PC12N, 2=PC24N), date of visits by PC (0, 12, 24, 36, set to -1 if not visited in that round), date of exit from PC.
 * 3) Function to update PC lists - count number of deaths and then draw n_dropout people to drop out (by age and gender) (n_dropout = number lost in data - number died in model). If one of those people died, then draw another person.
 * 3b) Update reserve list for PC12N and PC24N - maybe easier to create a new reserve list?
 * 4) write input functions for number of people visited at each timestep in PC0
 * 5) Pull out data from PC0 on number of people visited at each timestep.
 */
#include "pc.h"
#include "structures.h"
#include "constants.h"

/******************* Functions which schedule people for PC enrolment and visits: *******************/

/* This determines how we stratify PC participants' HIV status: */
int get_PC_HIV_stratum(individual *pc_participant){
    if (pc_participant->HIV_status==UNINFECTED)
        return 0;
    else if (pc_participant->HIV_status>UNINFECTED && (pc_participant->ART_status>ARTNEG) && (pc_participant->ART_status<ARTDEATH))
        return 1;  /* If HIV+ and aware: */
    else
        return 2;   /* If HIV+ but unaware: */
}

/* This takes the population of currently alive people (using age_list) and firstly sub-divides them into
 * a 'chips_sampling_frame', e.g. dividing up men and women, as CHiPs tends to visit more women than men.
 * We can also exclude certain people (e.g. <15 years old) as needed.
 * The sampling frame is then drawn from to pick chips_sample->size_n_m men and chips_sample->size_n_f women
 * who will be the people visited in the year.
 * We then shuffle each of these lists, so the shuffled list will be the order in which people are visited.
 * Finally we call schedule_chips_visits() which sets the number of people to be visited in each timestep
 * so that the total number of people visited in a year adds up to the correct total.
 * NOTE: The way we divide up the population means that we may try to visit people who died during the year.
 * However, this should not be a big factor, and I think it may even mimic CHiPs in that people may move/die
 * between enumeration/mapping and CHiPs visit.
 * The variable pc_round tells us what we are sampling: pc_round=0: PC0 entrant; pc_round=1: PC12N entrant; pc_round=2: PC24N entrant. */



void remove_extras_from_timestep_recruitment(patch_struct *patch, int p, int pc_enrolment_round, int g, int ap, int i_pc_category, int original_size, int number_to_remove){

    int *TEMP_SAMPLE_FRAME_TO_REMOVE;
    int *TEMP_LIST_TO_REMOVE;

    TEMP_SAMPLE_FRAME_TO_REMOVE = malloc(original_size*sizeof(int));
    TEMP_LIST_TO_REMOVE = malloc(number_to_remove*sizeof(int));
    if (TEMP_SAMPLE_FRAME_TO_REMOVE==NULL || TEMP_LIST_TO_REMOVE==NULL){
        printf("Unable to allocate TEMP_SAMPLE_FRAME_TO_REMOVE or TEMP_LIST_TO_REMOVE in remove_extras_from_timestep_recruitment(). Execution aborted.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    int n = 0;
    int i_dt,j;

    for (i_dt=0;i_dt<patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round]; i_dt++){
        for (j=0; j<patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]; j++){
            TEMP_SAMPLE_FRAME_TO_REMOVE[n] = i_dt;
            n++;
        }
    }
    if (n!=original_size){
        printf("Error - can't match sample size in remove_extras_from_timestep_recruitment(). Exiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }


    gsl_ran_choose(rng, TEMP_LIST_TO_REMOVE, number_to_remove, TEMP_SAMPLE_FRAME_TO_REMOVE, original_size, sizeof (int));



    for (n=0; n<number_to_remove; n++){
        /* We remove one person who was supposed to be seen from timestep TEMP_LIST_TO_REMOVE[n]. */
        i_dt = TEMP_LIST_TO_REMOVE[n];
        patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]--;
    }


    for (i_dt=0;i_dt<patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round]; i_dt++){
        if (patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]<0){
            printf("Error - have negative number of people to see per timestep in remove_extras_from_timestep_recruitment(). Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

    }



    free(TEMP_SAMPLE_FRAME_TO_REMOVE);
    free(TEMP_LIST_TO_REMOVE);

}

/* Draw up a list of people to be enrolled by PC (including reserves). */
void create_popart_pc_sample(patch_struct *patch, age_list_struct *age_list, PC_sample_struct *PC_sample, parameters *param, int pc_enrolment_round, int p){
    int g;
    int aa, ai,i, ap;
    int i_pc_category; /* Index splitting up the population by HIV status etc. */

    /* We are drawing our PC sample at the beginning of the round. Lots of things can happen in the mean time
     * (death, HIV infection, learning status) that make that person move outside the subpopulation - for example if someone was
     * uninfected when sampled but got HIV in the mean time they should no longer be in the 'uninfected' subpopulation.
     * So we need 'reserves' - people who are in the same population who can replace a person if needed. */
    int n_reserves;    /* Number of reserves to add for a given subpopulation. */
    double prop_reserves = 0.25; /* We want a minimum of 10% more people to account for deaths, movement to other subpopulations (e.g. HIV- gets infected). */
    int min_reserves = 15; /* Arbitrary minimum number of reserves. */


    /* For use with FOLLOW_INDIVIDUAL - we store these the first time we find them so we can find them easily next time: */
    int g_persontofollow = -1; /* Default value indicates that the FOLLOW_INDIVIDUAL did not turn up when going through - this is because they are too young to be visited by CHiPs. */
    int ap_persontofollow = -1;
    int i_pc_category_persontofollow = -1;

    individual *enrollee; /* Temporary pointer for the person we are enrolling at tthe time to make code more readable. As it points at an existing person, no need to malloc it. */
    int original_cohort_size = patch[p].param->PC_params->cohort_size;

    if (pc_enrolment_round<0 || pc_enrolment_round>2){
        printf("ERROR: pc_enrolment_round=%d takes range 0, 1 or 2 (PC0 entrant, PC12N entrant, PC24N entrant. Exiting\n",pc_enrolment_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    int max_number = 0;
    for (g=0;g<N_GENDER;g++){
        for (aa=(AGE_PC_MIN-AGE_ADULT); aa<=(AGE_PC_MAX-AGE_ADULT); aa++){
            ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);
            if (max_number<age_list->age_list_by_gender[g]->number_per_age_group[ai])
                max_number = age_list->age_list_by_gender[g]->number_per_age_group[ai];
        }
    }
    //printf("Maximum = %i\n",max_number);
    /* This acts as a temporary store for the PC sampling frame for a given subpopulation. We then draw the number of people
     * we want in the PC sample (+reserves) from this using gsl_ran_sample().
     * We assume that there are no more than 1/50th of the total population in a given age/gender/HIV status group.
     * Because of the way the sampling frame is constructed we need an array to store them. */
    //long TEMP_SAMPLE_FRAME[N_PC_HIV_STRATA][MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP];
    long *TEMP_SAMPLE_FRAME[N_PC_HIV_STRATA];
    for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
        TEMP_SAMPLE_FRAME[i_pc_category] = malloc(max_number*sizeof(long));
        if (TEMP_SAMPLE_FRAME[i_pc_category]==NULL){
            printf("Unable to allocate TEMP_SAMPLE_FRAME in create_popart_pc_sample(). Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    int n_sample[N_PC_HIV_STRATA];
    int number_not_recruited; /* Variable used if we can't make the full cohort size. */

    //  /* This is a temporary store of the sampling frame. Only called each round at present so should be OK as a local variable. */
    //  PC_sample_struct *PC_sampling_frame;
    //  PC_sampling_frame = malloc(sizeof(PC_sample_struct));
    //  if(PC_sampling_frame==NULL){ /* Check memory allocated successfully. */
    //      printf("Unable to allocate PC_sampling_frame in create_popart_pc(). Execution aborted.");
    //      printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
    //      fflush(stdout);
    //      exit(1);
    //  }

    /* Set counters to zero (we use these to index within each array). */
    //  for (g=0;g<N_GENDER;g++)
    //      /* ap is the variable we use as the index for age group in the PC_sample_struct structure.*/
    //      for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++)
    //          for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++)
    //              PC_sampling_frame->number_to_visit[g][ap][i_pc_category] = 0;


    for (g=0;g<N_GENDER;g++)
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++)
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++)
                PC_sample->next_person_to_see[g][ap][i_pc_category] = 0;


    for (g=0;g<N_GENDER;g++){
        /* aa is age index for age_list by gender (not adjusting for youngest age group).
         * Here it is chosen to correspond to ages from AGE_PC_MIN (=18) to AGE_PC_MAX (=44) inclusive. */
        for (aa=(AGE_PC_MIN-AGE_ADULT); aa<=(AGE_PC_MAX-AGE_ADULT); aa++){
            ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);

            ap = aa-(AGE_PC_MIN-AGE_ADULT); /* This is the relationship between the index aa and ap (the index for the PC_sample_struct structures). */
            //          printf("aa=%i ap=%i in create_popart_pc_sample()\n",aa,ap);
            //          fflush(stdout);

            /* TEMP_SAMPLE_FRAME[][] is a sampling frame for a given age group and gender.
             * So set counter to zero, and list to 'blank' (-1) here. */
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                n_sample[i_pc_category] = 0;
                for (i=0;i<MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP;i++)
                    TEMP_SAMPLE_FRAME[i_pc_category][i] = -1;
            }


            for (i=0;i<age_list->age_list_by_gender[g]->number_per_age_group[ai];i++){
                enrollee = age_list->age_list_by_gender[g]->age_group[ai][i];
                if (enrollee->PC_cohort_index==-1){  /* Only consider people who aren't already in PC. */
                    /* For this person work out what PC HIV stratum they belong to: */
                    i_pc_category = get_PC_HIV_stratum(enrollee);


                    //if (n_sample[i_pc_category]<MAX_NUMBER_PC_PARTICIPANTS_PER_GROUP){
                    /* Add the person to the sampling frame: */
                    TEMP_SAMPLE_FRAME[i_pc_category][n_sample[i_pc_category]] = enrollee->id;
                    n_sample[i_pc_category]++;

                    //                  if (g==1 && ap==15 && i_pc_category==2){
                    //                      printf("Adding %li to sample\n",enrollee->id);
                    //                  }

                    /*if(p==0 && enrollee.HIV_status>0)
                    {
                        printf("PCround %d HIVposIndivID %ld\n", 0, enrollee->id);
                        fflush(stdout);
                    }*/

                    /* For debugging. */
                    if(enrollee->id==FOLLOW_INDIVIDUAL && enrollee->patch_no==FOLLOW_PATCH){
                        printf("Possible PC enrolment %ld %d %d in round %d \n",enrollee->id,ai,i,pc_enrolment_round);
                        fflush(stdout);
                        /* Now store their characteristics so it's easier to find them: */
                        g_persontofollow = g;
                        ap_persontofollow = ap;
                        i_pc_category_persontofollow = i_pc_category;
                    }
                    /* Also for debugging. */
                    if(enrollee->cd4==DUMMYVALUE || enrollee->cd4==DEAD){
                        printf("Error -trying to schedule PC enrolment for dead/non-existent person %ld\n",enrollee->id);
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }


                    //}

                }
            }
            /* Now we have got the people of this age/gender to sample from, so we want to draw the people to visit. */
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                /* First calculate the number of reserves needed: */
                if(param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]>0){
                    n_reserves = (int) round(prop_reserves*param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]);
                    if (n_reserves<min_reserves)
                        n_reserves = min_reserves;

                }
                else
                    n_reserves = 0; /* No need for reserves if we're not seeing anyone. */
                //printf("param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] = %i n_sample=%i\n",param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round],n_sample[i_pc_category]);

                if (g==1 && ap==15 && i_pc_category==2){
                    printf("n_reserves = %i\n",n_reserves);
                    printf("param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] = %i n_sample=%i\n",param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round],n_sample[i_pc_category]);
                }

                /* Make sure no error in counting number of people: */
                if (n_sample[i_pc_category]<param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]){
                    /* Right now  param->PC_params->number_enrolled_in_PC_round[] is not actually the number of people who are recruited as we cannot exceed the sample size.
                     So we will reduce param->PC_params->number_enrolled_in_PC_round[]: */

                    number_not_recruited = param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] - n_sample[i_pc_category];

                    remove_extras_from_timestep_recruitment(patch, p, pc_enrolment_round, g, ap, i_pc_category, param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round], number_not_recruited);

                    param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round] = n_sample[i_pc_category];

                    /* Also update the cohort size to reflect non-recruitment: */
                    patch[p].param->PC_params->cohort_size = patch[p].param->PC_params->cohort_size - number_not_recruited;

                    //for (all_pc_rounds=pc_enrolment_round+1;all_pc_rounds<NPC_ROUNDS;all_pc_rounds++)
                    //  patch[p].param->PC_params->number_seen_in_PC_round[g][ap][i_pc_category][all_pc_rounds] = param->PC_params->number_seen_in_PC_round[g][ap][i_pc_category][pc_enrolment_round];

                    //printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    //fflush(stdout);
                    //exit(1);
                }

                /* Make sure that the total sample size including reserves doesn't exceed the number of people available. */
                PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category] = fmin(param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]+n_reserves,n_sample[i_pc_category]);
                if (g==1 && ap==15 && i_pc_category==2){
                    printf("PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category] = %li\n",PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category]);
                }

                //printf("number_in_sample_including_reserves = %ld\n",PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category]);
                //fflush(stdout);

                /* Choose the people who will be in PC in this subpopulation: */
                if (PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category]>0){
                    //                  if (g==1 && ap==15 && i_pc_category==2){
                    //                      printf("TEMP_SAMPLE_FRAME[i_pc_category] = %li\n",TEMP_SAMPLE_FRAME[i_pc_category][5]);
                    //                  }
                    gsl_ran_choose(rng, PC_sample->list_ids_potential_enrollees[g][ap][i_pc_category], PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category], TEMP_SAMPLE_FRAME[i_pc_category], n_sample[i_pc_category], sizeof (long));
                    /* Randomise the order (as gsl_ran_choose maintains the order of the original list). */
                    gsl_ran_shuffle(rng, PC_sample->list_ids_potential_enrollees[g][ap][i_pc_category], PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category], sizeof (long));
                }
            }
        }
    }

    if (original_cohort_size>patch[p].param->PC_params->cohort_size){
        printf("Warning: IBM population in certain strata was smaller than PC sample needed in pc_enrolment_round=%i. Cohort size reduced from %i to %i.\n",pc_enrolment_round,original_cohort_size,patch[p].param->PC_params->cohort_size);
        fflush(stdout);
    }


    /* Check to see if this person was visited. Note that if g_persontofollow==-1 then they are too young to be in chips_sampling_frame->list_ids_to_visit. */
    if (p==FOLLOW_PATCH && g_persontofollow>-1){
        for(i=0;i<param->PC_params->number_enrolled_in_PC_round[g_persontofollow][ap_persontofollow][i_pc_category_persontofollow][pc_enrolment_round];i++){
            if(PC_sample->list_ids_potential_enrollees[g_persontofollow][ap_persontofollow][i_pc_category_persontofollow][i]==FOLLOW_INDIVIDUAL){
                printf("PC enrolment for adult %ld, gender %d PC cat=%d from patch %d now scheduled\n",PC_sample->list_ids_potential_enrollees[g_persontofollow][ap_persontofollow][i_pc_category_persontofollow][i],g_persontofollow,i_pc_category_persontofollow,p);
                fflush(stdout);
            }
        }
    }


    for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++)
        free(TEMP_SAMPLE_FRAME[i_pc_category]);
}




/* Carry out the PC enrolment visits for a given timestep at time t. */
void carry_out_PC_enrolment_per_timestep(int t0, int t_step, patch_struct *patch, int p, int pc_enrolment_round){
    int g;
    int ap;
    int j;
    int i_pc_category; /* Index splitting up the population by HIV status etc. */
    int i_pc_thisperson;

    //int number_deleted = 0; /* Count how many people we didn't get to put in PC.  For debugging. */
    //printf("Carrying out enrolment for round %i at time %i + %i\n",pc_enrolment_round,t0,t_step);
    /* i_dt tells us the index of the timestep within PC_sample corresponding to the current time t so we can work out how many people to visit.*/
    int i_dt;

    int n_dropped_from_cohort = 0; /* Keep track of the number of people we couldn't enroll. */

    individual *enrollee; /* Temporary pointer for the person we are enrolling at tthe time to make code more readable. As it points at an existing person, no need to malloc it. */

    //printf("Number of timesteps = %i in round %i. Start = %i %i. End = %i %i\n",patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round],pc_enrolment_round, patch[p].param->PC_params->PC_START_YEAR[pc_enrolment_round], patch[p].param->PC_params->PC_START_TIMESTEP[pc_enrolment_round], patch[p].param->PC_params->PC_END_YEAR[pc_enrolment_round], patch[p].param->PC_params->PC_END_TIMESTEP[pc_enrolment_round]);

    i_dt = (t0 - patch[p].param->PC_params->PC_START_YEAR[pc_enrolment_round])*N_TIME_STEP_PER_YEAR + (t_step-patch[p].param->PC_params->PC_START_TIMESTEP[pc_enrolment_round]);

    /* For debugging: */
    if (i_dt>=MAX_N_TIMESTEPS_PER_PC_ROUND){
        printf("Problem - MAX_N_TIMESTEPS_PER_PC_ROUND is too small. We are %i timesteps into PC round %i. Exiting\n",i_dt,pc_enrolment_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if(i_dt<0||i_dt>=patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round])
    {
        printf("ERROR i_dt=%i is outside range [0,%i] in PC round %i at time %i+%i*TIME_STEP. Exiting\n",i_dt,patch[p].param->PC_params->n_timesteps_per_round[pc_enrolment_round],pc_enrolment_round,t0,t_step);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* Now PC visits each sub-population (currently men/women by year age group ap and HIV status stratum) in turn:
     * Go through the list of id's of people to visit, and visit them: */

    for (g=0;g<N_GENDER;g++){
        /* Run from 18 to 44. */
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                j = 0; /* This keeps track of how many people we have successfully enrolled of this subpopulation.
                We want to see patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round] people. */
                //              if (g==1 && ap==15 && i_pc_category==2)
                //                  printf("number_to_see=%i total_to_see = %i\n",patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round],patch[p].param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]);

                //              if (g==1 && ap==15 && i_pc_category==2){
                //                  printf("Want to enroll n= %i people in g=%i ac=%i i_pc=%i at dt=%i\n",patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round],g,ap,i_pc_category,i_dt);
                //              }


                while (j<(patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round])){

                    //                  printf("j= %i numberseen=%i\n",j,patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]);
                    //                  fflush(stdout);
                    /* Check if there are enough people in the PC sample (including reserves) - otherwise reduce the number of people to be seen and print a warning. */
                    if (patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category]<patch[p].PC_sample->number_in_sample_including_reserves[g][ap][i_pc_category]){
                        //                      printf("i= %i p=%i \n",i,p);
                        //                      fflush(stdout);
                        //                      printf("x= %li \n",patch[p].PC_sample->list_ids_potential_enrollees[g][ap][i_pc_category][i]);
                        //                      fflush(stdout);

                        enrollee = &(patch[p].individual_population[patch[p].PC_sample->list_ids_potential_enrollees[g][ap][i_pc_category][patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category]]]);
                        i_pc_thisperson = get_PC_HIV_stratum(enrollee);


                        //printf("Enrolling %li\n",enrollee.id);
                        //                      if (g==1 && ap==15 && i_pc_category==2){
                        //                          printf("Trying to enroll %li at dt=%i\n",enrollee->id,i_dt);
                        //                      }

                        /* We only keep this person if they are still in the right HIV stratum
                         * (so we get the correct number in the cohort who are HIV-, HIV+ unaware and HIV+ aware).
                         * Also make sure they are NOT already enrolled. */
                        //if (i_pc_thisperson==i_pc_category && enrollee->PC_cohort_index==-1){
                        if (i_pc_thisperson==i_pc_category){
                            //printf("Visiting %li\n",enrollee.id);
                            PC_enroll_person(enrollee, patch, t0+t_step*TIME_STEP, p, pc_enrolment_round, g, ap, i_pc_category);

                            //                          if (g==1 && ap==15 && i_pc_category==2){
                            //                              printf("Enrolling %li at dt=%i\n",enrollee->id,i_dt);
                            //                          }

                            patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category]++;  /* Move on to the next person on the list. */
                            j++;  /* We have successfully enrolled someone so increase j. */
                        }
                        else
                            patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category]++;  /* The current person was not eligible so move on to next person (will use reserves). Don't increment j as this was not successful. */
                    }
                    else{
                        patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round]--;
                        /* Also update cohort size as it will be smaller. */
                        patch[p].param->PC_params->cohort_size--;

                        /* Update how many people have been dropped out. */
                        n_dropped_from_cohort++;
                        //printf("Cohort_size reduced by one to = %i\n",patch[p].param->PC_params->cohort_size);

                        patch[p].param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]--;
                        //for (all_pc_rounds=0;all_pc_rounds<NPC_ROUNDS;all_pc_rounds++)
                        //patch[p].param->PC_params->number_seen_in_PC_round[g][ap][i_pc_category][all_pc_rounds]--;

                        //remove_a_visit(int t0, int t_step, patch_struct *patch, int p, int pc_enrolment_round, int g, int ap, int i_pc_category);



                        //number_deleted += 1;
                        //printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        //fflush(stdout);
                        //exit(1);
                        //                      if (g==1 && ap==15 && i_pc_category==2){
                        //                          printf("New target to enroll n= %i people at dt=%i\n",patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round],i_dt);
                        //                      }
                    }
                }

                ///* Update this index ready for the next timestep: */
                //patch[p].PC_sample->next_person_to_see[g][ap][i_pc_category] += patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_enrolment_round];
            }
        }
    }

    //printf("At enrolment timestep %i : Number in [1][15][2] = %i \n",i_dt,patch[p].param->PC_params->number_enrolled_in_PC_round[1][15][2][0]);
    if (n_dropped_from_cohort>0)
        printf("Warning: had to reduce PC sample by %i.\n",n_dropped_from_cohort);

    //printf("number not enrolled = %i\n",number_deleted);
}


//void remove_a_visit(int t0, int t_step, patch_struct *patch, int p, int pc_round, int g, int ap, int i_pc_category){
//  patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_round]
//
//}



void PC_enroll_person(individual *indiv, patch_struct *patch, double t, int p, int pc_enrolment_round, int g, int ap, int i_pc_category){
    /* Because of the way we draw CHiPs visits at the beginning of the year, it is possible some people
     * die before they are visited. If this is the case then do nothing more.
     * They are deleted from age_list so won't be in next year's sample. */


    if (indiv->cd4==DUMMYVALUE){
        printf("Trying to PC-enrol a non-existent person %d %ld !!! Exiting\n",p,indiv->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if(indiv->id==FOLLOW_INDIVIDUAL && p==FOLLOW_PATCH){
        printf("Enrolling %li at time %6.4lf g=%i ap=%i i_pc_category=%i\n",indiv->id,t,g,ap,i_pc_category);
        fflush(stdout);
    }

    //  if(indiv->id==25389 && p==0){
    //      printf("****Enrolling %li at time %6.4lf g=%i ap=%i i_pc_category=%i\n",indiv->id,t,g,ap,i_pc_category);
    //      fflush(stdout);
    //  }

    //  if (g==1 && ap==15 && i_pc_category==2){
    //      printf("Adding %li to sample. Number in that group = %i\n",indiv->id,patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category]);
    //  }

    /* Add this person to the cohort. */
    patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category]] = indiv->id;


    /* Increment the counter of number of people. */
    patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category]++;
    indiv->PC_cohort_index = patch[p].PC_cohort_data->PC_cohort_counter[pc_enrolment_round] + 100000*pc_enrolment_round; /* Person is in PC cohort. */
    //printf("XPC_cohort_index = %i\n",indiv->PC_cohort_index);
//  if(indiv->id==25389 && p==0){
//      printf("PC cohort index = %i. PC cohort size = %i\n",indiv->PC_cohort_index,patch[p].param->PC_params->cohort_size);
//  }

    /* PC_cohort_counter acts as the index for people in PCX_cohort_data.
     * Set up a temporary local variable (i.e. with shorter name) to make code more readable. */
    int pc_counter = patch[p].PC_cohort_data->PC_cohort_counter[pc_enrolment_round];

    /* Add this person's data to PC_cohort_data: */

    if (pc_enrolment_round==0){
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].gender = g;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].ap = ap;
        //printf("pc_counter = %i\n",pc_counter);
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].RETAINED_IN_COHORT[0] = 1;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].PC_visit_dates[0] = t;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].id =  indiv->id;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].HIV_status[0] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].ART_status[0] =  indiv->ART_status;
        if (indiv->HIV_status==UNINFECTED)
            patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].serodiscordant_status[0] = indiv->n_HIVpos_partners;
        else
            patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].serodiscordant_status[0] = 0;

        //printf("PC0 Enrolled  %li with pc_id=%i\n",patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].id,pc_counter);
    }
    else if (pc_enrolment_round==1){
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].gender = g;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].ap = ap;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].RETAINED_IN_COHORT[1] = 1;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].PC_visit_dates[1] = t;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].id =  indiv->id;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].HIV_status[1] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].ART_status[1] =  indiv->ART_status;
        //printf("*PC12 Enrolled  %li with pc_id=%i\n",patch[p].PC_cohort_data->PC12N_cohort_data[pc_counter].id,pc_counter);
    }
    else if (pc_enrolment_round==2){
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].gender = g;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].ap = ap;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].RETAINED_IN_COHORT[2] = 1;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].PC_visit_dates[2] = t;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].id =  indiv->id;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].HIV_status[2] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].ART_status[2] =  indiv->ART_status;
        //printf("*PC24 Enrolled  %li with pc_id=%i\n",patch[p].PC_cohort_data->PC24N_cohort_data[pc_counter].id,pc_counter);
    }
    //printf("id = %li\n",indiv->id);
    //if (pc_counter==53)
    //printf("***Enrolled  %li\n",patch[p].PC_cohort_data->PC0_cohort_data[pc_counter].id);

    /* Note this must be patch[p].PC_cohort_data->PC_cohort_counter, not the local variable pc_counter. */
    patch[p].PC_cohort_data->PC_cohort_counter[pc_enrolment_round]++;
    //printf("PC cohort counter = %i %i\n",patch[p].PC_cohort_data->PC_cohort_counter[pc_enrolment_round],indiv->PC_cohort_index);
}





// CURRENTLY UNIMPLEMENTED CODE.
///* Function takes the list of people who were in this round and removes a proportion who are not retained in the following round.
// * Function also resets visit counter for new round. */
//void PC_next_cohort_round(patch_struct *patch, int p, int pc_round){
//  int g, ap, i_pc_category; /* Indices splitting up the cohort by gender, age, HIV status etc. */
//  int n_remaining;
//
//  for (g=0;g<N_GENDER;g++){
//      /* Run from 18 to 44. */
//      for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
//          for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
//              /* How many people are left in this round? */
//              n_remaining = (int) round(patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category] * patch[p].param->PC_params->PC_retention[pc_round]);
//              if (n_remaining>0){
//                  gsl_ran_choose(rng, patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category], n_remaining, patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category], patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category], sizeof (long));
//              }
//              patch[p].PC_cohort->number_in_cohort[g][ap][i_pc_category] = n_remaining;
//          }
//      }
//  }
//}

/* Carry out the PC visits (AFTER enrolment round) for a given timestep at time t. */
void carry_out_PC_visits_per_timestep(int t0, int t_step, patch_struct *patch, int p, int pc_round, int pc_enrolment_round){
    int g;
    int ap;
    int i_pc_category; /* Index splitting up the population by HIV status etc. */

    int i_n; /* Temporary store to reduce verbosity of line. */

    /* i_dt tells us the index of the timestep within PC_sample corresponding to the current time t.*/
    int i_dt;

    i_dt = (t0 - patch[p].param->PC_params->PC_START_YEAR[pc_round])*N_TIME_STEP_PER_YEAR + (t_step-patch[p].param->PC_params->PC_START_TIMESTEP[pc_round]);

    /* For debugging: */
    if (i_dt>=MAX_N_TIMESTEPS_PER_PC_ROUND){
        printf("Problem - MAX_N_TIMESTEPS_PER_PC_ROUND is too small. We are %i timesteps into PC round %i. Exiting\n",i_dt,pc_round);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if(i_dt<0||i_dt>=patch[p].param->PC_params->n_timesteps_per_round[pc_round]){
        printf("ERROR i_dt=%i is outside range [0,%i] in PC round %i at time %6.4lf. Exiting\n",i_dt,patch[p].param->PC_params->n_timesteps_per_round[pc_round],pc_round,t0+t_step*TIME_STEP);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    //printf("In timestep dt=%i round %i we are looking to see %i people enrolled in round %i. Number in round = %i\n",i_dt,pc_round,patch[p].param->PC_params->number_seen_by_PC_per_timestep[1][15][2][i_dt][pc_round],pc_enrolment_round,patch[p].param->PC_params->number_enrolled_in_PC_round[1][15][2][pc_enrolment_round]);


    //**** HERE TRY THIS?
    //  printf("Hey2 %li\n",PC_sample->list_ids_potential_enrollees[0][14][2][0]);
    //printf("Hey2 id = %li\n",patch[0].PC_cohort->list_ids_in_cohort[0][14][2][0]);
    //  fflush(stdout);


    /* Now PC visits each sub-population (currently men/women by year age group ap and HIV status stratum) in turn:
     * Go through the list of id's of people to visit, and visit them: */
    for (g=0;g<N_GENDER;g++){
        /* Run from 18 to 44. */
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                //printf("Hey %i %i %i %i\n",patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category], patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category] + patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_round], patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category]<patch[p].param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round]);
                while ((patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category]<patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category] + patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_round])
                        && (patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category]<patch[p].param->PC_params->number_enrolled_in_PC_round[g][ap][i_pc_category][pc_enrolment_round])){

                    i_n = patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category];
                    //printf("In timestep dt=%i round %i we are looking to see %i people g=%i ap=%i i_pc=%i. Next person = %li\n",i_dt,pc_round,patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_round],g,ap,i_pc_category,patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]);
                    /* Visit this person. */
                    if (patch[p].individual_population[patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]].id<=0)
                        printf("Zero ids: id = %li, g=%i ap=%i i_pc = %i, i_dt=%i\n",patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n],g,ap,i_pc_category,i_dt);

                    //if (pc_round==3)
                    //  printf("Going to visit %li %li g=%i ap=%i i_pc_cat=%i n_to_visit=%i\n",patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n],patch[p].individual_population[patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]].id,g,ap,i_pc_category,patch[p].param->PC_params->number_seen_by_PC_per_timestep[g][ap][i_pc_category][i_dt][pc_round]);

                    //if (patch[p].individual_population[patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]].id==52067)
                    //printf("id = %li,g=%i ap=%i i_pc = %i, i_dt=%i\n",patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n],g,ap,i_pc_category,i_dt);
                    PC_visit_person(&(patch[p].individual_population[patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]]),
                            patch, t0+t_step*TIME_STEP, p, pc_round, pc_enrolment_round); /* Send the address (ie pointer) to this person. */

                    /*if(p==0 && patch[p].individual_population[patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]].HIV_status>0)
                    {
                        printf("PCround %d HIVposIndivID %ld\n", pc_round, (patch[p].individual_population[patch[p].PC_cohort->list_ids_in_cohort[g][ap][i_pc_category][i_n]]).id);
                        fflush(stdout);
                    }*/

                    patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category]++;
                }

            }
        }
    }
}


/* Store new data on person at PC round pc_round.
 * pc_enrolment_round reflects that people may have been recruited at different PC rounds (PC0, PC12N, PC24N). */
void PC_visit_person(individual *indiv, patch_struct *patch,  double t, int p, int pc_round, int pc_enrolment_round){
    /* Because of the way we draw CHiPs visits at the beginning of the year, it is possible some people
     * die before they are visited. If this is the case then do nothing more.
     * They are deleted from age_list so won't be in next year's sample. */
    if (indiv->cd4==DUMMYVALUE){
        printf("Trying to PC visit a non-existent person %d %ld !!! Exiting\n",p,indiv->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    int n = indiv->PC_cohort_index - 100000*pc_enrolment_round;
    //printf("HereRound%i %i\n",pc_round,n);
    if (n<0 || n>=patch[p].param->PC_params->cohort_size){
        printf("ERROR cohort index in PC_visit_person %i %i not defined when visiting %li - exiting\n",n,indiv->PC_cohort_index,indiv->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (p==0 && pc_enrolment_round==0 && pc_round==1){
        if (indiv->HIV_status>UNINFECTED){
            if (patch[p].PC_cohort_data->PC0_cohort_data[n].HIV_status[0]==UNINFECTED && VERBOSE_OUTPUT==1)
                printf("New PC12 infection id=%li\n",indiv->id);
        }
    }

    if (pc_enrolment_round==0){

        patch[p].PC_cohort_data->PC0_cohort_data[n].RETAINED_IN_COHORT[pc_round] = 1;
        patch[p].PC_cohort_data->PC0_cohort_data[n].PC_visit_dates[pc_round] = t;
        patch[p].PC_cohort_data->PC0_cohort_data[n].HIV_status[pc_round] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC0_cohort_data[n].ART_status[pc_round] =  indiv->ART_status;
        if (indiv->HIV_status==UNINFECTED)
            patch[p].PC_cohort_data->PC0_cohort_data[n].serodiscordant_status[pc_round] = indiv->n_HIVpos_partners;
        else
            patch[p].PC_cohort_data->PC0_cohort_data[n].serodiscordant_status[pc_round] = 0;
    }
    else if (pc_enrolment_round==1){
        patch[p].PC_cohort_data->PC12N_cohort_data[n].RETAINED_IN_COHORT[pc_round] = 1;
        patch[p].PC_cohort_data->PC12N_cohort_data[n].PC_visit_dates[pc_round] = t;
        patch[p].PC_cohort_data->PC12N_cohort_data[n].HIV_status[pc_round] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC12N_cohort_data[n].ART_status[pc_round] =  indiv->ART_status;
    }
    else if (pc_enrolment_round==2){
        patch[p].PC_cohort_data->PC24N_cohort_data[n].RETAINED_IN_COHORT[pc_round] = 1;
        patch[p].PC_cohort_data->PC24N_cohort_data[n].PC_visit_dates[pc_round] = t;
        patch[p].PC_cohort_data->PC24N_cohort_data[n].HIV_status[pc_round] =  indiv->HIV_status;
        patch[p].PC_cohort_data->PC24N_cohort_data[n].ART_status[pc_round] =  indiv->ART_status;
    }

    //  if (pc_enrolment_round==0)
    //      printf("***Visited  %li in round %i\n",patch[p].PC_cohort_data->PC0_cohort_data[n].id,pc_round);

}


/* get_pc_round

Return the index of the current PC round

Arguments
---------
t0 : int
    Current year
t_step : int
    Current time step
patch : pointer to a patch_struct
    Patch structure containing info on PC round times.  Assumed to have a  patch[0].param->PC_params->PC_END_TIMESTEP[NPC_ROUNDS]
    

Returns
-------
Index of the current PC round 0 = PC0, 1 = PC12 etc given the current (discrete) time.  A value of
-1 is returned if the input time is not a PC round.  
*/
int get_pc_round(int t0, int t_step, patch_struct *patch, int p){
    
    /* PC is only recruited in patch p=0. */
    if(p > 0){
        return -1;
    }

    /* If before/after PC years return -1. 
    Note that this is slightly more computationally efficient as we deal with most cases first. */
    if(t0 < patch[p].param->PC_params->PC_START_YEAR[0] ||
        t0 > patch[p].param->PC_params->PC_END_YEAR[NPC_ROUNDS - 1]){
        return -1;
    }

    int pc_round = 0;
    while(pc_round < NPC_ROUNDS){
        if((t0 == patch[p].param->PC_params->PC_START_YEAR[pc_round] &&
            t_step >= patch[p].param->PC_params->PC_START_TIMESTEP[pc_round]) ||
            (t0 > patch[p].param->PC_params->PC_START_YEAR[pc_round] &&
            t0 < patch[p].param->PC_params->PC_END_YEAR[pc_round]) ||
            (t0 == patch[p].param->PC_params->PC_END_YEAR[pc_round] &&
            t_step <= patch[p].param->PC_params->PC_END_TIMESTEP[pc_round])){
            return pc_round;
        }
        pc_round = pc_round +1;
    }
    return -1;
}


void reset_visit_counter(patch_struct *patch, int p){

    int g, ap, i_pc_category; /* Indices splitting up the cohort by gender, age, HIV status etc. */
    /* Reset the counter for visits: */
    for (g=0;g<N_GENDER;g++){
        /* Run from 18 to 44 (inclusive). */
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                /* This is the counter for who we visit next during visits for each round. */
                patch[p].PC_cohort->next_person_to_visit[g][ap][i_pc_category] = 0;

            }
        }
    }
}


/* Function used to debug PC sample - need to work out why we have enrolled more people than we are supposed to.
 * This function will check the number in each group and compare it against the input. */
void check_pc_after_enrollment(patch_struct *patch, int p){
    int Ncohort[N_GENDER][AGE_PC_MAX-AGE_PC_MIN+1][N_PC_HIV_STRATA];
    individual *person;
    int n_id;

    int g, ap, i_pc_category; /* Indices splitting up the cohort by gender, age, HIV status etc. */
    for (g=0; g<N_GENDER; g++){
        for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
            for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
                Ncohort[g][ap][i_pc_category] = 0;
            }
        }
    }


    for (n_id=0; n_id<patch[p].id_counter; n_id++){
        person = &(patch[p].individual_population[n_id]);
        if (person->PC_cohort_index>=0){
            g = patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].gender;
            ap = patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].ap;

            /* Check status is OK: */
            if (patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].HIV_status[0]<UNINFECTED|| patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].HIV_status[0]>CHRONIC){
                printf("PANIC! Exiting check_pc_after_enrollment()\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }

            /* Look up that person's status in PC0: */
            if (patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].HIV_status[0]==UNINFECTED)
                i_pc_category = 0;
            else if (patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].HIV_status[0]>UNINFECTED && patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].ART_status[0]>ARTNEG && patch[p].PC_cohort_data->PC0_cohort_data[person->PC_cohort_index].ART_status[0]<ARTDEATH)
                i_pc_category = 1;
            else
                i_pc_category = 2;

            Ncohort[g][ap][i_pc_category]++;

        }
    }

    //  for (g=0; g<N_GENDER; g++){
    //      for (ap=0; ap<(AGE_PC_MAX-AGE_PC_MIN+1); ap++){
    //          for (i_pc_category=0; i_pc_category<N_PC_HIV_STRATA; i_pc_category++){
    //              patch[0].param->PC_params->number_seen_in_PC_round[g][ap][i_pc_category][0]
    //          }
    //      }
    //  }


}
