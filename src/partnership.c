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

/* This file includes:
 * Partnership formation
 * Partnership dissolution
 */


/************************************************************************/
/******************************* Includes  ******************************/
/************************************************************************/

#include "partnership.h"
#include "utilities.h"
#include "output.h"

/************************************************************************/
/******************************** functions *****************************/
/************************************************************************/

/* Function does: fill in pair as a partnership between ind1 and ind2 formed at t_form_partnership
 *                  adds the pair to planned_breakups and n_planned_breakups
 *                  updates susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership accordingly
 *                  NOTE THAT UPDATE OF available_partners and n_available_partners is done outside of this function, afterwards
 *                  (because several partnerships are planned in advance based on the index of individuals in the list available_partnership,
 *                  it messes things up if we update the list after each partnership formation so we do it once they are all formed)
 * function arguments: pointer to the pair and the 2 individuals, flag to indicate if this partnership is related to sexual worker related activity, time of partnership formation, list of planned_breakups and n_planned_breakups, list of susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership, poiter to parameters.
 * Function returns: nothing */

void new_partnership(individual* ind1, individual* ind2, int sexual_worker_related_flag, double t_form_partnership,
        all_partnerships *overall_partnerships,
        parameters *param, debug_struct *debug,
        file_struct *file_data_store){
    /////// before calling this function one needs to check whether it is possible to form such a partnership
    /////// ie if ind1 and ind2 have some "free" partnerships,
    /////// are indeed of opposite sex,
    /////// and are not already in a partnership with each other!

    partnership *pair = &overall_partnerships->partner_pairs[*overall_partnerships->n_partnerships];
    int g;
    int age_f, age_m, risk_f, risk_m;

    /// FOR DEBUGGING, checking that these individuals are not dead!
    if(ind1->cd4==DEAD || ind2->cd4==DEAD)
    {
        printf("Problem: trying to make partnership with dead person...\n");
        if (ind1->cd4==DEAD)
            printf("Dead person ID = %li, other person = %li, dead person in patch %i\n",ind1->id,ind2->id,ind1->patch_no);
        else
            printf("Dead person ID = %li, other person = %li, dead person in patch %i\n",ind2->id,ind1->id,ind2->patch_no);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if(ind1->cd4==DUMMYVALUE || ind2->cd4==DUMMYVALUE)
    {
        printf("Problem: trying to make partnership with non-existent...\n");
        if (ind1->cd4==DUMMYVALUE)
            printf("Non-existent person ID = %li, other person = %li\n",ind1->id,ind2->id);
        else
            printf("Non-existent person ID = %li, other person = %li\n",ind2->id,ind1->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    //printf("Partnership formation involving individuals %ld and %ld\n",ind1->id,ind2->id);
    if((ind1->id==FOLLOW_INDIVIDUAL  && ind1->patch_no==FOLLOW_PATCH) || (ind2->id==FOLLOW_INDIVIDUAL  && ind2->patch_no==FOLLOW_PATCH))
    {
        //if (VERBOSE_OUTPUT==1){
        printf("-----------------------------------------------\n");
        printf("Starting partnership formation between individuals %ld from patch %d and %ld from patch %d\n",ind1->id,ind1->patch_no,ind2->id,ind2->patch_no);
        printf("Ind %ld from patch %d has %d partners and %d HIV+ partners (before this one)\n",ind1->id,ind1->patch_no,ind1->n_partners,ind1->n_HIVpos_partners);
        printf("Ind %ld from patch %d has %d partners and %d HIV+ partners (before this one)\n",ind2->id,ind2->patch_no,ind2->n_partners,ind2->n_HIVpos_partners);
        printf("----");
        print_individual(ind1);
        printf("----");
        print_individual(ind2);
        printf("----");
        fflush(stdout);
        //}
    }

    /* filling in a new partnership which points to the two individuals ind1 and ind2 and drawing a partnership duration */
    pair->ptr[ind1->gender] = ind1;
    pair->ptr[ind2->gender] = ind2;

    /* Filling in the age assortativity matrix (for debuging) */
    if(ind1->gender == FEMALE){
        age_f = get_age_group(ind1->DoB,t_form_partnership, AGE_GROUPS, N_AGE);
        age_m = get_age_group(ind2->DoB,t_form_partnership, AGE_GROUPS, N_AGE);

        risk_f = ind1->sex_risk;
        risk_m = ind2->sex_risk;

    }else
    {
        age_f = get_age_group(ind2->DoB,t_form_partnership, AGE_GROUPS, N_AGE);
        age_m = get_age_group(ind1->DoB,t_form_partnership, AGE_GROUPS, N_AGE);

        risk_f = ind2->sex_risk;
        risk_m = ind1->sex_risk;
    }
    pair->begin = t_form_partnership; /*current time in years*/

    if(CHECK_AGE_AND_RISK_ASSORTATIVITY == 1)
    {
        debug->age_of_partners_at_partnership_formation[(int)t_form_partnership - param->start_time_simul][age_f][age_m] ++;
        debug->risk_of_partners_at_partnership_formation[(int)t_form_partnership - param->start_time_simul][risk_f][risk_m] ++;
    }

    /* duration (in number of time steps) of the partnership */
    pair->duration_in_time_steps = time_to_partnership_dissolution(param,ind1->sex_risk,ind2->sex_risk, ind1->patch_no, ind2->patch_no); //// IF WE WANT SOMETHING ASYMETRICAL ACCORDING TO GENDER WILL NEED TO SPECIFY WHICH IS THE MALE AND WHICH IS THE FEMALE
    
    /* assign the partnership type */
    if (sexual_worker_related_flag == SEXUAL_WORKER_RELATED) {
        pair->sexual_worker_related = SEXUAL_WORKER_RELATED;
    }
    else {
        pair->sexual_worker_related = !SEXUAL_WORKER_RELATED;
    }

    if(DEBUG_PARTNERSHIP_DURATION ==1)
    {

        if (((ind1->patch_no+ind2->patch_no)==1) && (ind1->sex_risk==2) && (ind2->sex_risk==2)){
            file_data_store->DUR_BETWEEN_HIGHHIGH = fopen(file_data_store->filename_DUR_BETWEEN_HIGHHIGH,"a");
            fprintf(file_data_store->DUR_BETWEEN_HIGHHIGH,"%6.4f\n",(pair->duration_in_time_steps)*TIME_STEP);
            fclose(file_data_store->DUR_BETWEEN_HIGHHIGH);
        }
        else if(((ind1->patch_no+ind2->patch_no)==1) && (ind1->sex_risk==1) && (ind2->sex_risk==1)){
            file_data_store->DUR_BETWEEN_MEDMED = fopen(file_data_store->filename_DUR_BETWEEN_MEDMED,"a");
            fprintf(file_data_store->DUR_BETWEEN_MEDMED,"%6.4f\n",(pair->duration_in_time_steps)*TIME_STEP);
            fclose(file_data_store->DUR_BETWEEN_MEDMED);
        }
        else if(((ind1->patch_no+ind2->patch_no)==1) && (ind1->sex_risk==0) && (ind2->sex_risk==0)){
            file_data_store->DUR_BETWEEN_LOWLOW = fopen(file_data_store->filename_DUR_BETWEEN_LOWLOW,"a");
            fprintf(file_data_store->DUR_BETWEEN_LOWLOW,"%6.4f\n",(pair->duration_in_time_steps)*TIME_STEP);
            fclose(file_data_store->DUR_BETWEEN_LOWLOW);
        }
        else if (((ind1->patch_no+ind2->patch_no)==0) && (ind1->sex_risk==2) && (ind2->sex_risk==2)){
            file_data_store->DUR_WITHIN_HIGHHIGH = fopen(file_data_store->filename_DUR_WITHIN_HIGHHIGH,"a");
            fprintf(file_data_store->DUR_WITHIN_HIGHHIGH,"%6.4f\n",(pair->duration_in_time_steps)*TIME_STEP);
            fclose(file_data_store->DUR_WITHIN_HIGHHIGH);
        }
        else if(((ind1->patch_no+ind2->patch_no)==0) && (ind1->sex_risk==1) && (ind2->sex_risk==1)){
            file_data_store->DUR_WITHIN_MEDMED = fopen(file_data_store->filename_DUR_WITHIN_MEDMED,"a");
            fprintf(file_data_store->DUR_WITHIN_MEDMED,"%6.4f\n",(pair->duration_in_time_steps)*TIME_STEP);
            fclose(file_data_store->DUR_WITHIN_MEDMED);
        }
        else if(((ind1->patch_no+ind2->patch_no)==0) && (ind1->sex_risk==0) && (ind2->sex_risk==0)){
            file_data_store->DUR_WITHIN_LOWLOW = fopen(file_data_store->filename_DUR_WITHIN_LOWLOW,"a");
            fprintf(file_data_store->DUR_WITHIN_LOWLOW,"%6.4f\n",(pair->duration_in_time_steps)*TIME_STEP);
            fclose(file_data_store->DUR_WITHIN_LOWLOW);
        }
    }

    // This is for debug to check that all partnerships are broken up at some point
    /*if(t_form_partnership + pair->t_dissolve*TIME_STEP>param->end_time_simul)
    {
        printf("We will not see this partnership breakup (breakup planned for %g).\n",t_form_partnership + pair->t_dissolve*TIME_STEP);
        fflush(stdout);
        getchar();
    }*/

    /* updating the partners of ind1 and ind2 accordingly */
    ind1->n_partners++;
    if ((sexual_worker_related_flag == SEXUAL_WORKER_RELATED) && (ind1->sexual_worker_status == SEXUAL_WORKER)) {
        ind1->n_clients++;
    }
    ind1->partner_pairs[ind1->n_partners-1] = pair;
    if(ind2->HIV_status>0 && ind1->HIV_status==0) /* then we need to tell ind1 that he has a new HIV partner */
    {
        if(ind1->patch_no != ind2->patch_no)
        {
            ind1->n_HIVpos_partners_outside ++;
        }
        ind1->n_HIVpos_partners ++;
        ind1->partner_pairs_HIVpos[ind1->n_HIVpos_partners-1] = pair;
    }
    ind1->n_lifetime_partners++;

    ind2->n_partners++;
    if ((sexual_worker_related_flag == SEXUAL_WORKER_RELATED) && (ind2->sexual_worker_status == SEXUAL_WORKER)) {
        ind2->n_clients++;
    }
    ind2->partner_pairs[ind2->n_partners-1] = pair;
    if(ind1->HIV_status>0 && ind2->HIV_status==0)  /* then we need to tell ind2 that he has a new HIV partner */
    {
        if(ind1->patch_no != ind2->patch_no)
        {
            ind2->n_HIVpos_partners_outside ++;
        }
        ind2->n_HIVpos_partners++;
        ind2->partner_pairs_HIVpos[ind2->n_HIVpos_partners-1] = pair;
    }
    ind2->n_lifetime_partners++;

    if(ind1->patch_no != ind2->patch_no)
    {
        ind1->n_partners_outside++;
        ind2->n_partners_outside++;
        ind1->n_lifetime_partners_outside++;
        ind2->n_lifetime_partners_outside++;
    }

    /* adding this partnership to planned_breakups and n_planned_breakups */
    int time_breakup = N_TIME_STEP_PER_YEAR*(t_form_partnership- param->start_time_simul) + pair->duration_in_time_steps ;
    if(time_breakup<MAX_N_YEARS*N_TIME_STEP_PER_YEAR)
    {

        /* Check if we've run out of memory: */
        if (overall_partnerships->n_planned_breakups[time_breakup]>=(overall_partnerships->size_planned_breakups[time_breakup])){

            /* Note that realloc does not work (we need to pass a pointer to the pointer, which is really complicated as it propagates through several functions (so maybe make planned_breakups[time_breakup] **), so ditch this code for now and use the following lines: */
            printf("Unable to re-allocate planned_breakups[i]. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);


            /*
            size_planned_breakups[time_breakup] = size_planned_breakups[time_breakup] + RESIZEMEM_BREAKUP;
            printf("Reallocating stuff for planned_breakups[] now\n");
            fflush(stdout);

            planned_breakups[time_breakup] = realloc(planned_breakups[time_breakup],size_planned_breakups[time_breakup]*sizeof(partnership*));
            if (planned_breakups[time_breakup] == NULL){
                printf("Unable to re-allocate planned_breakups[i]. Execution aborted.");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
             */
        }

        overall_partnerships->planned_breakups[time_breakup][overall_partnerships->n_planned_breakups[time_breakup]] = pair;
        overall_partnerships->n_planned_breakups[time_breakup]++;
    }

    /* if the partnership is serodiscordant and the HIV- is not yet in the list of susceptible_in_serodiscordant_partnership,
     * adding the HIV- partner to susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership */
    if(is_serodiscordant(pair))
    {
        for(g=0 ; g<N_GENDER ; g++)
        {
            if(pair->ptr[g]->HIV_status==0) /* this is the susceptible individual in the serodiscordant partnership which is formed */
            {
                if(pair->ptr[g]->idx_serodiscordant==-1) /* this means this individual is not yet in the list susceptible_in_serodiscordant_partnership*/
                {
                    /* then add the susceptible from susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership */
                    //printf("n_susceptible_in_serodiscordant_partnership BEFORE - part formation %li\n",*n_susceptible_in_serodiscordant_partnership);
                    add_susceptible_to_list_serodiscordant_partnership(pair->ptr[g], overall_partnerships->susceptible_in_serodiscordant_partnership,  overall_partnerships->n_susceptible_in_serodiscordant_partnership);
                    //printf("n_susceptible_in_serodiscordant_partnership AFTER - part formation %li\n",*n_susceptible_in_serodiscordant_partnership);
                }
            }
        }
    }

    /*if(ind1->id==5486)
    {
        print_here_string("000000000000000000000000000000000000000000000000",0);
        print_here_string("Just formed a partnership involving individual ",ind1->id);
        print_here_string("From patch ",ind1->patch_no);
        print_here_string("Gender ",ind1->gender);
        print_here_string("Circ ",ind1->circ);
        print_here_string("Risk ",ind1->sex_risk);
        print_here_string("Yearob ",(int) ind1->DoB);
        print_here_string("Max_n_partners ",ind1->max_n_partners);
        print_here_string("Current_n_partners ",ind1->n_partners);
        print_here_string("111111111111111111111111111111111111111111111111",1);
    }

    if(ind2->id==5486)
    {
        print_here_string("000000000000000000000000000000000000000000000000",0);
        print_here_string("Just formed a partnership involving individual ",ind2->id);
        print_here_string("From patch ",ind2->patch_no);
        print_here_string("Gender ",ind2->gender);
        print_here_string("Circ ",ind2->circ);
        print_here_string("Risk ",ind2->sex_risk);
        print_here_string("Yearob ",(int) ind2->DoB);
        print_here_string("Max_n_partners ",ind2->max_n_partners);
        print_here_string("Current_n_partners ",ind2->n_partners);
        print_here_string("111111111111111111111111111111111111111111111111",1);
    }*/

    if((ind1->id==FOLLOW_INDIVIDUAL  && ind1->patch_no==FOLLOW_PATCH) || (ind2->id==FOLLOW_INDIVIDUAL  && ind2->patch_no==FOLLOW_PATCH))
    {
        //if (VERBOSE_OUTPUT==1){
        printf("Finishing partnership formation between individuals %ld from patch %d and %ld from patch %d\n",ind1->id,ind1->patch_no,ind2->id,ind2->patch_no);
        printf("Ind %ld from patch %d has %d partners and %d HIV+ partners (now that this one has formed)\n",ind1->id,ind1->patch_no,ind1->n_partners,ind1->n_HIVpos_partners);
        printf("Ind %ld from patch %d has %d partners and %d HIV+ partners (now that this one has formed)\n",ind2->id,ind2->patch_no,ind2->n_partners,ind2->n_HIVpos_partners);
        printf("----");
        print_individual(ind1);
        printf("----");
        print_individual(ind2);
        printf("----");
        printf("-----------------------------------------------\n");
        fflush(stdout);
        //}
    }

}

/* Function does: randomly draws a time to partnership dissolution as function of risk group of two partners
 * function arguments: two integers representing two risk groups of partners
 * Function returns: the number of time steps after which partnership will be dissolved
 * NOTE: currently no matter what the risk group of partners is, the duration of partnerships is set to a random number
 * drawn from a unique Exp distribution with mean 10 years, discretized so that it is exactly a multiple of time steps
 * the number returned is the number of time steps from formation to dissolution
 * this should be changed so that the function is realistic and different according to the risk group of the two partners */
int time_to_partnership_dissolution(parameters *param, int r1, int r2, int p, int q){ //// IF WE WANT SOMETHING ASYMETRICAL ACCORDING TO GENDER WILL NEED TO SPECIFY WHICH IS THE MALE AND WHICH IS THE FEMALE
    int overall_risk;

    /* This essentially does overall_risk=max(r1,r2) */
    if (r1==HIGH||r2==HIGH)
        overall_risk=HIGH;
    else if (r1==MEDIUM||r2==MEDIUM)
        overall_risk=MEDIUM;
    else
        overall_risk=LOW;

    //int t_partnership_dissolution = (int) (gsl_ran_exponential (rng, param->breakup_scale_lambda[overall_risk])* N_TIME_STEP_PER_YEAR);
    /* For Weibull(lambda,k) k=1 corresponds to the exponential dist. k<1 means prob of failure decreases over time, k>1 means increases over time. */
    //int t_partnership_dissolution = (int) (gsl_ran_weibull(rng, param->breakup_scale_lambda[overall_risk], param->breakup_shape_k[overall_risk])* N_TIME_STEP_PER_YEAR);
    //printf("Using Weibull\n");
    /* Using Gamma distribution: */
    /* gsl_ran_gamma(rng,a,b). a = shape k, b = scale theta. */
    /* corresponds to pdf p(x) dx = {1 \over \Gamma(a) b^a} x^{a-1} e^{-x/b} dx [https://www.gnu.org/software/gsl/manual/html_node/The-Gamma-Distribution.html] */
    /* so mean = ab and var = ab^2 */
    int t_partnership_dissolution;
    if(p==q) // within patch
    {
        t_partnership_dissolution = (int) (gsl_ran_gamma (rng, param->breakup_shape_k[overall_risk], param->breakup_scale_lambda_within_patch[overall_risk])* N_TIME_STEP_PER_YEAR);
        /*if (r1==HIGH && r2==HIGH)
        {
            printf("mean duration of partnership INSIDE should be (a=%lg*b=%lg) = %lg\n",param->breakup_shape_k[overall_risk], param->breakup_scale_lambda_within_patch[overall_risk], param->breakup_shape_k[overall_risk]* param->breakup_scale_lambda_within_patch[overall_risk]);
            fflush(stdout);
        }*/

    }else // between patch
    {
        t_partnership_dissolution = (int) (gsl_ran_gamma (rng, param->breakup_shape_k[overall_risk], param->breakup_scale_lambda_between_patch[overall_risk])* N_TIME_STEP_PER_YEAR);
        /*if (r1==HIGH && r2==HIGH)
        {
            printf("mean duration of partnership INSIDE should be (a=%lg*b=%lg) = %lg\n",param->breakup_shape_k[overall_risk], param->breakup_scale_lambda_between_patch[overall_risk], param->breakup_shape_k[overall_risk]* param->breakup_scale_lambda_between_patch[overall_risk]);
            fflush(stdout);
        }*/

    }
    //printf("pdf_gamma(%lg | a=%lg, b=%lg) = %lg\n",1.0, 1.0, 11.951, gsl_ran_gamma_pdf (1.0,1.0, 11.951));
    //printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
    //fflush(stdout);
    //exit(1);
    if(t_partnership_dissolution<1)
    {
        /* if the breakup time is less than a time step away then partnership will never be broken up;
         instead we make this duration 1 time step
         and these partnerships will be broken at next time step */
        t_partnership_dissolution = 1;
    }
    return t_partnership_dissolution;
}

/* Function does: breakup a given partnership and updates the lists of available partners and of serodiscordant couples accordingly
 * function arguments: the time at which partnership is broken up (double), a pointer to the partnership breaking up, pointers to the lists of pop_available_partners and n_pop_available_partners and susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership
 * Function returns: nothing */
void breakup(double time_breakup, partnership* breakup, all_partnerships *overall_partnerships)
{

    /* Only breakup partnership if the two individuals are alive - otherwise nothing happens. */
    if(breakup->ptr[0]->cd4 > DEAD && breakup->ptr[1]->cd4 > DEAD)
    {

        if((breakup->ptr[0]->id==FOLLOW_INDIVIDUAL  && breakup->ptr[0]->patch_no==FOLLOW_PATCH) || (breakup->ptr[1]->id==FOLLOW_INDIVIDUAL && breakup->ptr[1]->patch_no==FOLLOW_PATCH))
        {
            //if (VERBOSE_OUTPUT==1){
            printf("-----------------------------------------------\n");
            printf("Starting partnership breakup between individuals %ld from patch %d and %ld from patch %d\n",breakup->ptr[0]->id,breakup->ptr[0]->patch_no,breakup->ptr[1]->id,breakup->ptr[1]->patch_no);
            printf("Ind %ld from patch %d has %d partners and %d HIV+ partners (before this one breaks up) \n",breakup->ptr[0]->id,breakup->ptr[0]->patch_no,breakup->ptr[0]->n_partners,breakup->ptr[0]->n_HIVpos_partners);
            printf("Ind %ld from patch %d has %d partners and %d HIV+ partners (before this one breaks up) \n",breakup->ptr[1]->id,breakup->ptr[1]->patch_no,breakup->ptr[1]->n_partners,breakup->ptr[1]->n_HIVpos_partners);
            printf("----");
            print_individual(breakup->ptr[0]);
            printf("----");
            print_individual(breakup->ptr[1]);
            printf("----");
            fflush(stdout);
            //}
        }
        // static long int n_breakups = 0; // for debug to count the number of breakups and check that all partnerships are broken up at some point
        int g;
        long index_partner;
        long id_partners[N_GENDER]; /*  contains the id of the 2 partners */
        int patch_partners[N_GENDER]; /*  contains the patch of the 2 partners */
        for(g=0 ; g<N_GENDER ; g++)
        {
            id_partners[g] = breakup->ptr[g]->id;
            patch_partners[g] = breakup->ptr[g]->patch_no;
        }

        individual *male = breakup->ptr[0];
        individual *female = breakup->ptr[1];

        individual *who = male;
        individual *other = female;

        /* Printing id of partners */
        //printf("Partnership between male %ld and female %ld will be dissolved\n",id_partners[0],id_partners[1]);
        //fflush(stdout);

        /* update of lists needs to be done before partnership does not exist anymore */

        update_list_available_partners_breakup(time_breakup, breakup, overall_partnerships->pop_available_partners, overall_partnerships->n_pop_available_partners);

        update_list_susceptibles_in_serodiscordant_partnerships_breakup(breakup, overall_partnerships->susceptible_in_serodiscordant_partnership, overall_partnerships->n_susceptible_in_serodiscordant_partnership);

        /* removing the female from list of partners of the male and vice versa*/
        if(who->patch_no != other->patch_no)
        {
            who->n_partners_outside--;
            other->n_partners_outside--;
        }
        for(g=0 ; g<N_GENDER ; g++)
        {
            if(g==1)
            {
                who = female;
                other = male;
            }
            index_partner = 0;

            while(who->partner_pairs[index_partner]->ptr[1-g]->id!=id_partners[1-g] || who->partner_pairs[index_partner]->ptr[1-g]->patch_no!=patch_partners[1-g])
            {
                index_partner++;
            }

            /* decreasing the number of partners by 1 */
            who->n_partners--;
            /* If this partnership is sexual worker related and this individual is sexual worker, then decrease the number of clients by 1*/
            if ((breakup->sexual_worker_related == SEXUAL_WORKER_RELATED) && (who->sexual_worker_status == SEXUAL_WORKER)) {
                who->n_clients--;
            }
            /* copying the last partner in place of this one - NOTE FOR DEBUGGING: AC does this the opposite way to MP who firstly swaps with the n-1 entry and then decreases n. */
            who->partner_pairs[index_partner] = who->partner_pairs[who->n_partners];
            /* same within the list of HIV positive partners if the partner is HIV positive */

            if(who->HIV_status==0 && other->HIV_status>0) /* only bother if the couple is serodiscordant */
            {
                index_partner = 0;
                while(who->partner_pairs_HIVpos[index_partner]->ptr[1-g]->id!=id_partners[1-g] || who->partner_pairs_HIVpos[index_partner]->ptr[1-g]->patch_no!=patch_partners[1-g])
                {
                    index_partner++;
                }
                who->n_HIVpos_partners--;
                if(who->patch_no != other->patch_no)
                {
                    who->n_HIVpos_partners_outside--;
                }
                who->partner_pairs_HIVpos[index_partner] = who->partner_pairs_HIVpos[who->n_HIVpos_partners];
            }
        }

        // for debug to count the number of breakups and check that all partnerships are broken up at some point:
        /*n_breakups++;
    printf("Simulated a total of: %ld breakups\n",n_breakups);
    fflush(stdout);*/

        if((breakup->ptr[0]->id==FOLLOW_INDIVIDUAL  && breakup->ptr[0]->patch_no==FOLLOW_PATCH) || (breakup->ptr[1]->id==FOLLOW_INDIVIDUAL && breakup->ptr[1]->patch_no==FOLLOW_PATCH))
        {
            if (VERBOSE_OUTPUT==1){
                printf("Finishing partnership breakup between individuals %ld from patch %d and %ld from patch %d\n",breakup->ptr[0]->id,breakup->ptr[0]->patch_no,breakup->ptr[1]->id,breakup->ptr[1]->patch_no);
                printf("Ind %ld from patch %d has %d partners and %d HIV+ partners (after this one breaks up)\n",breakup->ptr[0]->id,breakup->ptr[0]->patch_no,breakup->ptr[0]->n_partners,breakup->ptr[0]->n_HIVpos_partners);
                printf("Ind %ld from patch %d has %d partners and %d HIV+ partners (after this one breaks up)\n",breakup->ptr[1]->id,breakup->ptr[1]->patch_no,breakup->ptr[1]->n_partners,breakup->ptr[1]->n_HIVpos_partners);
                printf("----");
                print_individual(breakup->ptr[0]);
                printf("----");
                print_individual(breakup->ptr[1]);
                printf("----");
                printf("-----------------------------------------------\n");
                fflush(stdout);
            }
        }

    }

}

/* Function does: updates the lists of available partners upon breakup of a partnership
 * function arguments: the time at which partnership is broken up (double), a pointer to the partnership breaking up, pointers to the lists of pop_available_partners and n_pop_available_partners
 * Function returns: nothing */
void update_list_available_partners_breakup(double time_breakup, partnership* breakup, population_partners* pop_available_partners, population_size_all_patches *n_pop_available_partners){
    int g, ag, r, p;

    /* when breaking up a partnership, it needs to be removed and the partners need to be put back in the list of available partners */
    for(g=0 ; g<N_GENDER ; g++)
    {
        ag = get_age_group(breakup->ptr[g]->DoB, time_breakup, AGE_GROUPS, N_AGE);
        r = breakup->ptr[g]->sex_risk;
        p = breakup->ptr[g]->patch_no;
        /* Adding this person as an available partner in the list of available partners */
        pop_available_partners->pop_per_patch_gender_age_risk[p][g][ag][r][n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r]] = breakup->ptr[g];
        /* Keeping track of where this person is in that list */
        breakup->ptr[g]->idx_available_partner[breakup->ptr[g]->max_n_partners - breakup->ptr[g]->n_partners] = n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r];
        n_pop_available_partners->pop_per_patch[p].pop_size_per_gender_age_risk[g][ag][r]++;
    }

}


void add_susceptible_to_list_serodiscordant_partnership(individual* indiv, 
    individual** susceptible_in_serodiscordant_partnership, 
    long *n_susceptible_in_serodiscordant_partnership){
    /* Add a susceptible to the list of serodiscordant partnerships
    
    This function is called following infection of an individual with HIV- partners or partnership
    formation with an HIV+ individual.  
    
    Arguments
    ---------
    indiv : pointer to an individual structure
        Pointer to the susceptible individual
    
    susceptible_in_serodiscordant_partnership : pointer to array of individual structures
    
    n_susceptible_in_serodiscordant_partnership : pointer to array of individual structures
    
    Returns
    -------
    Nothing; 
    */
    
    indiv->idx_serodiscordant = n_susceptible_in_serodiscordant_partnership[0];
    
    susceptible_in_serodiscordant_partnership[n_susceptible_in_serodiscordant_partnership[0]] =
        indiv;
    
    n_susceptible_in_serodiscordant_partnership[0]++;
}


void remove_susceptible_from_list_serodiscordant_partnership(individual* indiv, individual** susceptible_in_serodiscordant_partnership, long *n_susceptible_in_serodiscordant_partnership){
    
    /* Remove a susceptible to the list of serodiscordant partnership (can happen upon partnership breakup (including because death of partner) or seroconversion)
    * function arguments: a pointer to the susceptible individual, and pointers to the list of susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership
    
    Returns
    -------
    Nothing
    */
    
    /* removing this susceptible from susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership */
    //int k;
    if((n_susceptible_in_serodiscordant_partnership[0] > 0) && (indiv->idx_serodiscordant != -1)){
        //print_here_string("---",0);
        /*printf("Removing individual %li from patch %d from the list of susceptibles in serodiscordant partnerships\n",indiv->id,indiv->patch_no);

        printf("- Individual before removal is:\n");
        print_individual(indiv);

        printf("--- The list before removal is:\n");
        for(k=0 ; k<*n_susceptible_in_serodiscordant_partnership ; k++)
        {
            printf("        Indiv %li from patch %d,\n",susceptible_in_serodiscordant_partnership[k]->id,susceptible_in_serodiscordant_partnership[k]->patch_no);
        }
        fflush(stdout);*/

        susceptible_in_serodiscordant_partnership[indiv->idx_serodiscordant] = susceptible_in_serodiscordant_partnership[n_susceptible_in_serodiscordant_partnership[0] - 1]; /* pointing to the last susceptible instead of this one */
        /* changing the index of that partnership which has been "moved" from last */
        susceptible_in_serodiscordant_partnership[indiv->idx_serodiscordant]->idx_serodiscordant = indiv->idx_serodiscordant;
        n_susceptible_in_serodiscordant_partnership[0]--;

        indiv->idx_serodiscordant = -1; // putting back this index to -1

        /*printf("- Individual after removal is:\n");
        print_individual(indiv);

        printf("--- The list after removal is:\n");
        for(k=0 ; k<*n_susceptible_in_serodiscordant_partnership ; k++)
        {
            printf("        Indiv %li from patch %d,\n",susceptible_in_serodiscordant_partnership[k]->id,susceptible_in_serodiscordant_partnership[k]->patch_no);
        }
        fflush(stdout);*/
    }
}


/* Function does: updates the list of serodiscordant partnerships upon partnership breakup (the two individuals and their partners can enter or leave the list)
 * function arguments: a pointer to the partnership breaking up, and pointers to the list of susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership
 * Function returns: nothing */
void update_list_susceptibles_in_serodiscordant_partnerships_breakup(partnership* breakup, individual** susceptible_in_serodiscordant_partnership, long *n_susceptible_in_serodiscordant_partnership)
{
    /* if the partnership is serodiscordant and the susceptible has no other positive partners, removing the susceptible from susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership */
    int g;

    if(is_serodiscordant(breakup))
    {
        for(g=0 ; g<N_GENDER ; g++)
        {
            if(breakup->ptr[g]->HIV_status==0) /* this is the susceptible individual in the serodiscordant partnershipwhich is breaking up */
            {
                if(breakup->ptr[g]->n_HIVpos_partners==1) /* this means once this partnership has broken up the susceptible will have no more HIV+ partners */
                {
                    /* then remove the susceptible from susceptible_in_serodiscordant_partnership and n_susceptible_in_serodiscordant_partnership */
                    remove_susceptible_from_list_serodiscordant_partnership(breakup->ptr[g], susceptible_in_serodiscordant_partnership,  n_susceptible_in_serodiscordant_partnership);
                    //print_here_string("---",0);
                    //print_here_string("remove from serodiscordant because breakup partnership",0);
                }
            }
        }

    }
}


void draw_nb_new_partnerships(patch_struct *patch, parameters *param, int patch_f, int patch_m){
    
   /* Function does: fills in param->balanced_nb_f_to_m with the number of new partnerships between a female of age ag_f, risk r_f and a male of age ag_m, risk r_m over a time unit
    * function arguments: pointers to the population_size, the stratified_population_size, and the parameters
    * Function returns: nothing
    * NOTE: the functions used within this are defined in utilities.c */
    
    calcul_n_new_partners_f_to_m(patch, param, patch_f, patch_m);
    //printf("calcul_n_new_partners_f_to_m between patch_f %d and patch_m %d: %lg\n", patch_f, patch_m,param->unbalanced_nb_f_to_m[2][2][3][2]);
    //fflush(stdout);
    calcul_n_new_partners_m_to_f(patch, param, patch_f, patch_m);
    //printf("calcul_n_new_partners_m_to_f between patch_f %d and patch_m %d: %lg\n", patch_f, patch_m,param->unbalanced_nb_m_to_f[3][2][3][2]);
    //fflush(stdout);
    balance_contacts_arithmetic(param);
    //printf("balance_contacts_arithmetic between patch_f %d and patch_m %d: %ld\n", patch_f, patch_m,param->balanced_nb_f_to_m[2][2][3][2]);
    //fflush(stdout);
}


void draw_n_new_partnerships(double time, long n, long sexual_worker_related_flag, parameters *param, int ag_f, int r_f, int ag_m, 
    int r_m, int *n_non_matchable, all_partnerships *overall_partnerships, patch_struct *patch, 
    int patch_f, int patch_m, debug_struct *debug, file_struct *file_data_store){
    
    /*
    Form a certain number of partnerships (n) between given age/risk groups in men and women, 
    and updates all lists accordingly
    
    Notes
    -----
    Individuals selected for partnership formation are chosen from those who 
    have "available partnerships"
    
    Parameters
    ----------
    time : double
        the time of partnership formation
    
    n : long
        the number of partnerships to form

    sexual_worker_related_flag : long
        a flag to indicate whether this partnership is related to sexual worker activity (to replace the non-seuxal workers with sexual workers)
    
    ag_f, ag_m : int
        the age group of the female and male respectively
    
    r_f, r_m : int
        the risk group of the female and male respectively
    
    n_non_matchable : pointer to an int
        the number of partnership we failed to form
    
    overall_partnerships : pointer to an all_partnerships struct
    
    patch : pointer to a patch_struct struct
    
    patch_f, patch_m : int
    
    debug : a pointer to a debug_struct
    
    file_data_store : a pointer to a file_struct
    
    Returns
    -------
    Nothing
    */
    
    long k = 0;
    int i, j;
    long initial_selected_male;
    long initial_selected_female;
    int tmp;
    int idx_found;
    int flag;
    int current_n_normal_partnership;
    int maximum_n_normal_partnership;


    /* initialising n_non_matchable to zero: this will then be incremented to output the 
    number of partnerships (out of n) that failed to form */
    *n_non_matchable = 0;
    
    // Assign pointers to the population size structures for readability within the code
    population_size *popn_f;
    population_size *popn_m;
    popn_f = &overall_partnerships->n_pop_available_partners->pop_per_patch[patch_f];
    popn_m = &overall_partnerships->n_pop_available_partners->pop_per_patch[patch_m];
    
    individual ******popn_pgar;
    popn_pgar = overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk;

    //// ONCE WE HAVE DEBUGGED EVErYTHING SHOULD BE ABLE To GER RID OF idx_found and test on
    // positivity of j.  Leaving it for now as is useful for debugging
    
    /* Only run this if we draw 1 or more partnerships. */
    if(n > 0){
        
        /* draw the male partners */
        gsl_ran_choose(rng, 
            overall_partnerships->new_partners_m, 
            n, 
            overall_partnerships->partner_dummylist, 
            popn_m->pop_size_per_gender_age_risk[MALE][ag_m][r_m], 
            sizeof(long));
        
        // Loop through partnerships to be drawn
        for(k = 0; k < n; k++){

            // flag to indicate if this partnership successfully forms
            flag = 0;

            /* draw the female partner */
            gsl_ran_choose(rng, 
                overall_partnerships->new_partners_f_sorted, 
                1, 
                overall_partnerships->partner_dummylist, 
                popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f], 
                sizeof(long));

            
            // If it is not sexual worker related activity, then non-sexual worker involved partnership formation has priority
            if (sexual_worker_related_flag == !SEXUAL_WORKER_RELATED) {

                // Prefer to form normal female partnership first (replace selected female sexual workers with normal females if possible)
                if (popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]] -> sexual_worker_status == SEXUAL_WORKER) {
                    initial_selected_female = overall_partnerships->new_partners_f_sorted[0];
                    // Take the next one in the list that is not female sexual worker with free partnerships
                    if(overall_partnerships->new_partners_f_sorted[0] < popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1) {
                        do {
                            overall_partnerships->new_partners_f_sorted[0]++;
                        }while((overall_partnerships->new_partners_f_sorted[0] < popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1) && (popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]] -> sexual_worker_status == SEXUAL_WORKER));
                    }

                    // If not possible, starting from the first female until the initially selected one
                    if(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]] -> sexual_worker_status == SEXUAL_WORKER) {
                        overall_partnerships->new_partners_f_sorted[0] = 0;
                        while ((overall_partnerships->new_partners_f_sorted[0] < initial_selected_female) && (popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]] -> sexual_worker_status == SEXUAL_WORKER)) {
                            overall_partnerships->new_partners_f_sorted[0]++;
                        }
                    }

                    // If no more non-sexual workers are available, then consider to select sexual worker to form normal partnership
                    current_n_normal_partnership = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_clients;
                    maximum_n_normal_partnership = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_clients;
                    // If the selected sexual worker doesn't have enough normal partnership, loop through the next one
                    if ((overall_partnerships->new_partners_f_sorted[0] == initial_selected_female) && (current_n_normal_partnership == maximum_n_normal_partnership)) {
                        if(overall_partnerships->new_partners_f_sorted[0] < popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1) {
                            do {
                                overall_partnerships->new_partners_f_sorted[0]++;
                                current_n_normal_partnership = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_clients;
                                maximum_n_normal_partnership = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_clients;
                            }while((overall_partnerships->new_partners_f_sorted[0] < popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1) && (current_n_normal_partnership == maximum_n_normal_partnership));
                        }

                        // If not possible, starting from the first female until the initially selected one
                        if(current_n_normal_partnership == maximum_n_normal_partnership) {
                            overall_partnerships->new_partners_f_sorted[0] = 0;
                            current_n_normal_partnership = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_clients;
                            maximum_n_normal_partnership = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_clients;
                            while ((overall_partnerships->new_partners_f_sorted[0] < initial_selected_female) && (current_n_normal_partnership == maximum_n_normal_partnership)) {
                                overall_partnerships->new_partners_f_sorted[0]++;
                                current_n_normal_partnership = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_clients;
                                maximum_n_normal_partnership = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_clients;
                            }
                        }
                        // If can't find a available sexual worker with vacant normal partnership, skip this run
                        if (overall_partnerships->new_partners_f_sorted[0] == initial_selected_female) {
                            overall_partnerships->new_partners_f_non_matchable[*n_non_matchable] = k;
                            (*n_non_matchable)++;
                            flag = 0;
                            continue;
                        }
                    }
                }
            }

            // If it is sexual worker related activity, then selected female sexual worker whose n_cilents < max_cilents (also n_partners < max_partners)
            else {
                if (popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_clients == popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_clients) {
                    initial_selected_female = overall_partnerships->new_partners_f_sorted[0];
                    // Take the next one in the list that is female sexual worker with free sexual worker related partnerships
                    if(overall_partnerships->new_partners_f_sorted[0] < popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1) {
                        do {
                            overall_partnerships->new_partners_f_sorted[0]++;
                        }while((overall_partnerships->new_partners_f_sorted[0] < popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1) && (popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_clients == popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_clients));
                    }

                    // If not possible, starting from the first female until the initially selected one
                    if(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_clients == popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_clients) {
                        overall_partnerships->new_partners_f_sorted[0] = 0;
                        while ((overall_partnerships->new_partners_f_sorted[0] < initial_selected_female) && (popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_clients == popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_clients)) {
                            overall_partnerships->new_partners_f_sorted[0]++;
                        }
                    }

                    // If can't find a available sexual worker, skip this run
                    if (overall_partnerships->new_partners_f_sorted[0] == initial_selected_female) {
                        overall_partnerships->new_partners_f_non_matchable[*n_non_matchable] = k;
                        (*n_non_matchable)++;
                        flag = 0;
                        continue;
                    }
                }
            }

            /* Check that individuals are not already in a partnership with each other */
            if(
                are_in_partnership(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]], 
                    popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]) == 1){
                
                // If finding this male already as a partner, 
                // then choose a different one in the same age/risk group:
                initial_selected_male = overall_partnerships->new_partners_m[k];
                
                // Take the next one in the list that is not yet in a partnership with this female
                // and has free partnerships
                if(overall_partnerships->new_partners_m[k] < popn_m->pop_size_per_gender_age_risk[MALE][ag_m][r_m] - 1){
                    do
                    {
                        overall_partnerships->new_partners_m[k]++;
                    }while((overall_partnerships->new_partners_m[k] < popn_m->pop_size_per_gender_age_risk[MALE][ag_m][r_m] - 1) && ((are_in_partnership(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]], popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]) ==1) || (is_already_selected(overall_partnerships->new_partners_m, k, n)==1)));
                }

                /* If not possible, start back from male 0 in that group */
                if((are_in_partnership(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]], popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]) ==1) || (is_already_selected(overall_partnerships->new_partners_m, k, n)==1))
                {
                    overall_partnerships->new_partners_m[k] = 0;
                    /* if not possible, take the next one in the list that is not yet in a partnership with this female and has a free partnership, up to the initially selected one */
                    while( (overall_partnerships->new_partners_m[k] < initial_selected_male) && ((are_in_partnership(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]], popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]) ==1) || (is_already_selected(overall_partnerships->new_partners_m, k, n)==1)))
                    {
                        overall_partnerships->new_partners_m[k]++;
                    }
                }

                if(popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->id != popn_pgar[patch_m][MALE][ag_m][r_m][initial_selected_male]->id) /* if able to find an available male in this group which is not already in a partnership with this female */
                {
                    new_partnership(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]],
                            popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]], sexual_worker_related_flag,
                            time, overall_partnerships, param, debug, file_data_store);
                    
                    flag = 1;
                    (*overall_partnerships->n_partnerships) ++;

                    /* remove the corresponding idx_available_partnership right away, for females */
                    idx_found = 0;
                    for(i = 0; i <= popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners; i++){
                        if(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[i] == overall_partnerships->new_partners_f_sorted[0]){
                            idx_found = 1;
                            break;
                        }
                    }
                    if(idx_found == 0){
                        printf("problem a here, idx should have been found that has not (female indirectly chosen)\n");
                        fflush(stdout);
                    }
                    popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[i] = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners];
                    popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners] = -1;
                    /* and for males */
                    idx_found = 0;
                    for(i=0 ; i<=popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->max_n_partners - popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->n_partners ; i++)
                    {
                        if(popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->idx_available_partner[i] == overall_partnerships->new_partners_m[k])
                        {
                            idx_found = 1;
                            break;
                        }
                    }
                    if(idx_found==0)
                    {
                        printf("problem a here, idx should have been found that has not (male indirectly chosen)\n");
                        fflush(stdout);
                    }
                    popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->idx_available_partner[i] = popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->idx_available_partner[popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->max_n_partners - popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->n_partners];
                    popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->idx_available_partner[popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->max_n_partners - popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->n_partners] = -1;


                }
                else
                {
                    overall_partnerships->new_partners_f_non_matchable[*n_non_matchable] = k;
                    (*n_non_matchable)++;
                    flag = 0;
                }
            }else /* if the initially selected male was not already in a partnership with this female */
            {
                /* form a partnership */
                new_partnership(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]], popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]], sexual_worker_related_flag, time, overall_partnerships, param, debug, file_data_store);
                
                flag = 1;
                (*overall_partnerships->n_partnerships)++;
                /* remove the corresponding idx_available_partnership right away, for females */
                idx_found = 0;
                for(i=0 ; i<=popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners ; i++)
                {
                    if(popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[i] == overall_partnerships->new_partners_f_sorted[0])
                    {
                        idx_found = 1;
                        break;
                    }
                }
                if(idx_found==0)
                {
                    fflush(stdout);
                }
                popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[i] = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners];
                popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners] = -1;
                /* and for males */
                idx_found = 0;
                for(i=0 ; i<=popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->max_n_partners - popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->n_partners ; i++)
                {
                    if(popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->idx_available_partner[i] == overall_partnerships->new_partners_m[k])
                    {
                        idx_found = 1;
                        break;
                    }
                }
                if(idx_found==0)
                {
                    printf("problem here b, idx should have been found that has not (male directly chosen)\n");
                    fflush(stdout);
                }
                popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->idx_available_partner[i] = popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->idx_available_partner[popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->max_n_partners - popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->n_partners];
                popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->idx_available_partner[popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->max_n_partners - popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m[k]]->n_partners] = -1;

            }
            
            /* If we are able to match current selected female to a male, then remove this female */
            if (flag == 1) {

                /* removing this partnership with this female from the list of available partnerships */
                if(overall_partnerships->new_partners_f_sorted[0] < popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1)
                {
                    popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]] = popn_pgar[patch_f][FEMALE][ag_f][r_f][popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1]; /* pointing to the last female instead of the current one */
                    j = popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->max_n_partners - popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->n_partners - 1;
                    while(j>=0 && popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[j] != popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f] - 1)
                    {
                        j--;
                    }
                    if(j<0)
                    {
                        printf("Problem here j negative1\n");
                        fflush(stdout);
                    }
                    popn_pgar[patch_f][FEMALE][ag_f][r_f][overall_partnerships->new_partners_f_sorted[0]]->idx_available_partner[j] = overall_partnerships->new_partners_f_sorted[0]; /* telling the person that has moved that they have. */
                }
                popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f]--; /* decreasing the number of available females in that group by 1 */
            }
            // printf("after:%d\n", popn_f->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f]);
        }

        /* copies new_partners_m into new_partners_m_sorted */
        copy_array_long(overall_partnerships->new_partners_m_sorted, overall_partnerships->new_partners_m, n);

        /* sort new_partners_m_sorted (which can be not sorted if initial selected man was not selected in the end) */
        //qsort(new_partners_m_sorted, n, sizeof(long), (compfn) compare_longs);

        qsort(overall_partnerships->new_partners_m_sorted, n, sizeof(long), compare_longs);

        /* Remove individuals from the list of available partnerships, starting by the last ones, excluding those that we weren't able to match */
        for(k=n - 1 ; k>=0 ; k--)
        {
            tmp = 0;
            for(i=0 ; i<(*n_non_matchable) ; i++)
            {
                if(overall_partnerships->new_partners_m_sorted[k]==overall_partnerships->new_partners_m[overall_partnerships->new_partners_f_non_matchable[i]])
                {
                    tmp = 1;
                    break;
                }
            }
            if(tmp==0)
            {
                /* removing this male from the list of available partnerships */
                if(overall_partnerships->new_partners_m_sorted[k] < popn_m->pop_size_per_gender_age_risk[MALE][ag_m][r_m] - 1)
                {
                    popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m_sorted[k]] = popn_pgar[patch_m][MALE][ag_m][r_m][popn_m->pop_size_per_gender_age_risk[MALE][ag_m][r_m] - 1]; /* pointing to the last male instead of the current one */
                    j = popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m_sorted[k]]->max_n_partners - popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m_sorted[k]]->n_partners - 1;
                    while(j>=0 && popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m_sorted[k]]->idx_available_partner[j] != popn_m->pop_size_per_gender_age_risk[MALE][ag_m][r_m] - 1)
                    {
                        j--;
                    }
                    if(j<0)
                    {
                        printf("Problem here j negative2\n");
                        fflush(stdout);
                    }
                    popn_pgar[patch_m][MALE][ag_m][r_m][overall_partnerships->new_partners_m_sorted[k]]->idx_available_partner[j] = overall_partnerships->new_partners_m_sorted[k]; /* telling the person that has moved that they have. */
                }
                popn_m->pop_size_per_gender_age_risk[MALE][ag_m][r_m]--; /* decreasing the number of available males in that group by 1 */
            }
        }
    }
}


void draw_new_partnerships(double time, all_partnerships *overall_partnerships, patch_struct *patch,
    parameters *param, int patch_f, int patch_m, debug_struct *debug, 
    file_struct *file_data_store){
    // !!! here I kept param as an argument as we may want to generate a parameter set which is a
    // "mix" betwen the parameters of two patches if these are different

    /*
    Calculate the number of new partnerships to be drawn in a time step, form these
    partnerships and update all relevant lists

    Arguments
    ---------
    time : double
        the time of partnership formation in years
    overall_partnerships : pointer to all_partnerships struct
    
    patch : pointer to patch_struct struct
        the patch structure
    
    param : pointer to a parameters struct
        parameters of the model
    
    patch_f, patch_m : int
    
    debug : pointer to a debug_struct struct
    
    file_data_store: pointer to a file_struct struct
    
    
    Returns
    -------
    Nothing
    
    */

    
    int ag_f, ag_m, r_f, r_m; // indices
    int n_non_matchable;
    int tmp = 0;
    int index = 0;

    draw_nb_new_partnerships(patch, param, patch_f, patch_m);
    
    // Loop through all age/risk groups of both genders
    for(ag_f = 0; ag_f < N_AGE; ag_f++){
        for(r_f = 0; r_f < N_RISK; r_f++){
            for(ag_m = 0; ag_m < N_AGE; ag_m++){
                for(r_m = 0; r_m < N_RISK; r_m++){
                    index = 0;

                    //printf("ag_f = %d \t r_f = %d \t ag_m = %d \t r_m = %d \n", ag_f,r_f,ag_m,r_m);
                    //printf("Individual %ld's first index for available partner is: %ld\n",individual_population[10013].id,individual_population[10013].idx_available_partner[0]);
                    //printf("Number of partnerships to be drawn %ld\n", param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m]);
                    //printf("Number of available partnerships from females %ld\n", n_pop_available_partners->pop_size_per_age_risk_per_gender[FEMALE][ag_f][r_f]);
                    //printf("Number of available partnerships from males %ld\n", n_pop_available_partners->pop_size_per_age_risk_per_gender[MALE][ag_m][r_m]);
                    //fflush(stdout);

                    ////////////////////////////////////////
                    ///// IN CASE WE MAY DRAW A NUMBER OF PARTNERSHIPS LARGER THAN THE AVAILABILITY,
                    ///// FOR NOW ARTIFICIALLY REDUCING THE NUMBER OF PARTNERSHIPS TO BE DRAWN BUT THIS MAY NEED FURTHER THINKING
                    
                    if(param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] > overall_partnerships->n_pop_available_partners->pop_per_patch[patch_f].pop_size_per_gender_age_risk[FEMALE][ag_f][r_f]){
                        //getchar();
                        param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] = overall_partnerships->n_pop_available_partners->pop_per_patch[patch_f].pop_size_per_gender_age_risk[FEMALE][ag_f][r_f];
                    }
                    if(param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] > overall_partnerships->n_pop_available_partners->pop_per_patch[patch_m].pop_size_per_gender_age_risk[MALE][ag_m][r_m]){
                        //getchar();
                        param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] = overall_partnerships->n_pop_available_partners->pop_per_patch[patch_m].pop_size_per_gender_age_risk[MALE][ag_m][r_m];
                    }
                    ////////////////////////////////////////

                    draw_n_new_partnerships(time, param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m], !SEXUAL_WORKER_RELATED,
                        param, ag_f, r_f, ag_m, r_m, &n_non_matchable, overall_partnerships, patch, 
                        patch_f, patch_m, debug, file_data_store);
                    tmp = n_non_matchable;
                    /* if some of the pairs could not be formed because they were already in a partnership together and we couldn't find another suitable male for that female, we redraw the corresponding number of pairs */
                    while(tmp > 0 && index < 10){ /// HERE THE CONDITION index<10 is to avoid an infinite loop, think better about how to deal with the case where it's just impossible to form a new partnership (e.g. 1 available partner in each group but they are already in partnership)
                        index++;
                        draw_n_new_partnerships(time, n_non_matchable, !SEXUAL_WORKER_RELATED, param, ag_f, r_f, ag_m, r_m,
                            &tmp, overall_partnerships, patch, patch_f, patch_m, debug, 
                            file_data_store);
                        
                        n_non_matchable = tmp;
                    }
                    if(index == 10){
                        if(VERBOSE_OUTPUT == 1){
                            printf("Warning: some partnerships between males in ");
                            printf("age gp %d, risk gp %d and females in ", ag_m, r_m);
                            printf("age gp %d, risk gp %d couldn't be formed.\n", ag_f, r_f);
                            fflush(stdout);
                        }
                    }

                    // if allow female sexual workers in the simulation, the remaining partnerships are formed between female sexual workers
                    if (SEXUAL_WORKER_STRUCTURE == 1) {
                        int n_sexual_worker_non_matchable;
                        int n_relationship_sexual_worker = param->unbalanced_nb_m_to_f[ag_m][r_m][ag_f][r_f] - param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m];

                        index = 0;
                        // number of relationship between female sexual workers and males. Also have to ensure there are enough individuals to be drawn
                        if ((n_relationship_sexual_worker > 0) && (n_relationship_sexual_worker <= overall_partnerships->n_pop_available_partners->pop_per_patch[patch_m].pop_size_per_gender_age_risk[MALE][ag_m][r_m]) && (n_relationship_sexual_worker <= overall_partnerships->n_pop_available_partners->pop_per_patch[patch_f].pop_size_per_gender_age_risk[FEMALE][ag_f][r_f])) {
                            draw_n_new_partnerships(time, n_relationship_sexual_worker, SEXUAL_WORKER_RELATED, param, ag_f, r_f, ag_m, r_m,
                                &n_sexual_worker_non_matchable, overall_partnerships, patch, 
                                patch_f, patch_m, debug, file_data_store);
                            tmp = n_sexual_worker_non_matchable;
                            /* if some of the pairs could not be formed because they were already in a partnership together and we couldn't find another suitable male for that female sexual worker, we redraw the corresponding number of pairs */
                            while(tmp > 0 && index < 10  && (n_relationship_sexual_worker <= overall_partnerships->n_pop_available_partners->pop_per_patch[patch_m].pop_size_per_gender_age_risk[MALE][ag_m][r_m]) && (overall_partnerships->n_pop_available_partners->pop_per_patch[patch_f].pop_size_per_gender_age_risk[FEMALE][ag_f][r_f])){ /// HERE THE CONDITION index<10 is to avoid an infinite loop, think better about how to deal with the case where it's just impossible to form a new partnership (e.g. 1 available partner in each group but they are already in partnership)
                                index++;
                                draw_n_new_partnerships(time, n_sexual_worker_non_matchable, SEXUAL_WORKER_RELATED, param, ag_f, r_f, ag_m, r_m,
                                    &tmp, overall_partnerships, patch, patch_f, patch_m, debug, 
                                    file_data_store);
                                
                                n_sexual_worker_non_matchable = tmp;
                            }
                            if(index == 10){
                                if(VERBOSE_OUTPUT == 1){
                                    printf("Warning: some seuxal worker related partnerships between males in ");
                                    printf("age gp %d, risk gp %d and female sexual workers in ", ag_m, r_m);
                                    printf("age gp %d, risk gp %d couldn't be formed.\n", ag_f, r_f);
                                    fflush(stdout);
                                }
                            }
                        } 
                    }
                }
            }
        }
    }
}
