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

/* Useful functions for model
 *
 */


/************************************************************************/
/******************************* Includes  ******************************/
/************************************************************************/

#include "utilities.h"
#include "constants.h"
#include "init.h"


/* These are the functions in this file:
 * check_if_cannot_read_param() - prints a simple error message that we cannot find a given parameter.
 * print_here() - used as a simple debugging output, printing an integer.
 * print_here_string() - a slightly more sophisticated version of print_here() which can print a string too.
 * etc etc... to be filled in later on!
 * add_commas_to_calibration_output()
 *
 */

/************************************************************************/
/******************************** functions *****************************/
/************************************************************************/

/* Used in input.c and fitting.c to check the return value from fscanf. */
void check_if_cannot_read_param(int checkreadok, char *varname){
    if (checkreadok<1){
        printf("ERROR: Parameter %s not read. Exiting\n",varname);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}

void print_here(int n){
    printf("Here %i\n",n);
    fflush(stdout);
}

void print_here_string(char *s,int n){
    printf("Here %s %i\n",s,n);
    fflush(stdout);
}


void check_age_group_index(age_list_struct *age_list, int g, long person_id,int indextocheck){
    int i,j,aitrue;
    aitrue=-1; /* Assign dummy value. */
    for (i=0;i<(MAX_AGE-AGE_ADULT);i++){
        for (j=0;j<age_list->age_list_by_gender[g]->number_per_age_group[i];j++){
            if (age_list->age_list_by_gender[g]->age_group[i][j]->id==person_id)
                aitrue = i;
        }
    }
    if (aitrue==-1){
        printf("Person not found in check_age_group_index()\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (indextocheck!=aitrue){
        printf("Index for person %li does not match %i %i\n",person_id,indextocheck,aitrue);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

}


/* Function not currently used. */
void normalise_four_quantities(double *p1, double *p2, double *p3, double *p4){
    double normalization = *p1 + *p2 + *p3 + *p4;
    if (normalization>0){
        *p1 = *p1/normalization;
        *p2 = *p2/normalization;
        *p3 = *p3/normalization;
        *p4 = *p4/normalization;
    }
    // We don't use an else statement - hopefully
    else{
        printf("Not normalizing %f %f %f %f\n",*p1,*p2,*p3,*p4);
        *p1 = 0.25;
        *p2 = 0.25;
        *p3 = 0.25;
        *p4 = 0.25;
    }
}

/* Gives cumulative sums. Note we do not output a fourth var c4 as this is assumed ot be 1.
 * Not currently used. */
void cumulative_four_quantities(double p1, double p2, double p3, double p4, double *c1, double *c2, double *c3, double *c4){
    *c1 = p1;
    *c2 = *c1+p2;
    *c3 = *c2+p3;
    *c4 = *c3+p4;
    if (fabs(*c4-1.0)>1e-12){
        printf("ERROR: Cumulative total in cumulative_four_quantities() does not sum to 1. Exiting. \n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}


double hill_up(double x, double max_val, double exponent, double midpoint){
    /* Evaluate the upwards Hill function at point x.  
    
    This is only used in the MTCT transmission function so this is not currently used.  
    
    Arguments
    ---------
    x : double
        Point at which to evaluate the function
    max_val : double
        Max value parameters
    exponent : double
        Exponent parameter
    midpoint : double
        Midpoint parameter
    
    Returns
    -------
    Returns the Hill function (up) evaluated at point x (a double).  
    */
    double result = max_val * (pow(x, exponent))/((pow(x, exponent) + pow(midpoint, exponent)));
    return result;
}


double hill_down(double x, double max_val, double exponent, double midpoint){
    /* Evaluate the downwards Hill function at point x.  
    
    Used within schedule_new_hiv_test() to determine mean time until HIV test is scheduled.  
    
    
    Arguments
    ---------
    x : double
        Point at which to evaluate the function
    max_val : double
        Max value parameters
    exponent : double
        Exponent parameter
    midpoint : double
        Midpoint parameter
    
    Returns
    -------
    Returns the Hill function (down) evaluated at point x (a double).  
    
    */
    
    double result = max_val * (pow(midpoint, exponent))/((pow(x, exponent) + 
        pow(midpoint, exponent)));
    return result;
}


void calcul_population(population_size *pop, stratified_population_size *pop_strat){
    /* This function fills in pop_strat based on pop */

    int g, ag, r;

    /* summing populations */

    for(g=0 ; g<N_GENDER ; g++)
    {
        pop_strat->total_pop_size_per_gender[g] = 0;
        for(r=0 ; r<N_RISK ; r++)
        {
            pop_strat->pop_size_per_gender_risk[g][r] = 0;
        }
    }

    for(ag=0 ; ag<N_AGE ; ag++)
    {
        for(g=0 ; g<N_GENDER ; g++)
        {
            pop_strat->pop_size_per_gender_age[g][ag] = 0;
        }
        for(r=0 ; r<N_RISK ; r++)
        {
            //pop_strat->pop_size_per_age_risk[ag][r] = pop->pop_size_per_gender_age_risk[FEMALE][ag][r] + pop->pop_size_per_gender_age_risk[MALE][ag][r];
            for(g=0 ; g<N_GENDER ; g++)
            {
                pop_strat->pop_size_per_gender_age[g][ag] += pop->pop_size_per_gender_age_risk[g][ag][r];
                pop_strat->pop_size_per_gender_risk[g][r] += pop->pop_size_per_gender_age_risk[g][ag][r];
            }
        }
        //pop_strat->pop_size_per_age[ag] = pop_strat->pop_size_per_age_per_gender[FEMALE][ag] + pop_strat->pop_size_per_age_per_gender[MALE][ag];
        for(g=0 ; g<N_GENDER ; g++)
        {
            pop_strat->total_pop_size_per_gender[g] += pop_strat->pop_size_per_gender_age[g][ag];
        }
    }

    pop_strat->total_pop_size = pop_strat->total_pop_size_per_gender[FEMALE] + pop_strat->total_pop_size_per_gender[MALE];


    for(r=0 ; r<N_RISK ; r++)
    {
        //pop_strat->pop_size_per_risk[r] = pop_strat->pop_size_per_risk_per_gender[FEMALE][r] + pop_strat->pop_size_per_risk_per_gender[MALE][r];
        for(g=0 ; g<N_GENDER ; g++)
        {
            pop_strat->prop_pop_per_gender_risk[g][r] = (double)pop_strat->pop_size_per_gender_risk[g][r]/(double)pop_strat->total_pop_size_per_gender[g];
        }
    }

}

void calcul_pop_wider_age_groups(population_size *pop, population_size_one_year_age *pop_one_year){
    /* This function fills in pop based on pop_one_year */

    set_population_count_zero(pop);
    int total_pop = 0;
    int g, aa, ai, ag, r;
    for(ag=0 ; ag<N_AGE ; ag++)
    {
        for(aa=AGE_GROUPS_WITH_OLD[ag]- AGE_ADULT ; aa<AGE_GROUPS_WITH_OLD[ag+1]- AGE_ADULT ; aa++)
        {
            ai = pop_one_year->youngest_age_group_index + aa ; // ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);
            for(r=0 ; r<N_RISK ; r++)
            {
                for(g=0 ; g<N_GENDER ; g++)
                {
                    pop->pop_size_per_gender_age_risk[g][ag][r] += pop_one_year->pop_size_per_gender_age1_risk[g][ai][r];
                    total_pop+=pop_one_year->pop_size_per_gender_age1_risk[g][ai][r];
                }
            }

        }
    }

    for(r=0 ; r<N_RISK ; r++){
        for(g=0 ; g<N_GENDER ; g++){
            pop->pop_size_per_gender_age_risk[g][N_AGE-1][r] += pop_one_year->pop_size_oldest_age_group_gender_risk[g][r];
            total_pop+=pop_one_year->pop_size_oldest_age_group_gender_risk[g][r];
        }
    }
    //printf("total pop counted in calcul_pop_wider_age_groups() is %i\n",total_pop);
}

void calcul_prevalence(proportion_population_size *prevalence, population_size *pop, population_size *n_infected_wide_age_group, population_size_one_year_age *n_infected){
    /* This function fills in prevalence based on pop and n_infected (through calculation of n_infected_wide_age_group)*/

    calcul_pop_wider_age_groups(n_infected_wide_age_group, n_infected);

    int g, ag, r;
    for(g=0 ; g<N_GENDER ; g++)
    {
        for (ag=0; ag<N_AGE; ag++)
        {
            for(r=0 ; r<N_RISK ; r++)
            {
                prevalence->prop_size_per_gender_age_risk[g][ag][r] = n_infected_wide_age_group->pop_size_per_gender_age_risk[g][ag][r] / pop->pop_size_per_gender_age_risk[g][ag][r];
            }
        }
    }
}


void print_prevalence(population_size *pop, population_size *n_infected_wide_age_group, population_size_one_year_age *n_infected){
    /* This prints prevalence based on pop and n_infected*/

    calcul_pop_wider_age_groups(n_infected_wide_age_group, n_infected);

    int g, ag, r;
    for(g=0 ; g<N_GENDER ; g++)
    {
        if(g==0) printf("MALES\n"); else printf("FEMALES\n");
        for (ag=0; ag<N_AGE; ag++)
        {
            printf("age group %d\n",ag);
            for(r=0 ; r<N_RISK ; r++)
            {
                printf("risk group %d: prevalence = %lg\t", r, ((double) n_infected_wide_age_group->pop_size_per_gender_age_risk[g][ag][r]) / ( (double) pop->pop_size_per_gender_age_risk[g][ag][r]));
            }
            printf("\n");
        }
    }

    fflush(stdout);
}


void calcul_p_risk(int g, double p_risk[N_RISK][N_RISK], stratified_population_size *pop_strat, parameters *param){
    /* This fills in p_risk, so that p_risk[r,.] is the distribution over risk groups of sexual partners
     * of someone of gender g in risk group r.*/
    int i,j;
    //// I think it might be (slightly?) quicker to loop from j=0..N_RISK-1 and then add the extra param->assortativity outside the loop.
    for(i=0 ; i<N_RISK ; i++)
    {
        p_risk[i][i] = param->assortativity + (1-param->assortativity)*pop_strat->prop_pop_per_gender_risk[1-g][i];
        for(j=0 ; j<i ; j++)
        {
            p_risk[i][j] = (1-param->assortativity)*pop_strat->prop_pop_per_gender_risk[1-g][j];
        }
        for(j=i+1 ; j<N_RISK ; j++)
        {
            p_risk[i][j] = (1-param->assortativity)*pop_strat->prop_pop_per_gender_risk[1-g][j];
        }
    }
}


void calcul_n_new_partners_f_to_m(patch_struct *patch, parameters *param, int patch_f, int patch_m){
    int ag_f, ag_m ; /* indexes of ages for females and males */
    int r_f, r_m ; /* indexes of risk groups for females and males */

    /* Fills in p_risk_f */
    calcul_p_risk(FEMALE, param->p_risk_per_gender[FEMALE], patch[patch_f].n_population_stratified, param);


    /* Standardizes the relative number of new partners per risk group */
    //standardize_relative_number_partnerships_per_risk(FEMALE, param->xi_per_gender[FEMALE], patch[patch_f].n_population_stratified, param);

    /* Calculates the "desired" number of new partners for females per age and risk group */
    for(ag_f=0 ; ag_f<N_AGE ; ag_f++)
    {
        for(r_f=0 ; r_f<N_RISK ; r_f++)
        {
            for(ag_m=0 ; ag_m<N_AGE ; ag_m++)
            {
                for(r_m=0 ; r_m<N_RISK ; r_m++)
                {
                    //printf("ag_f = %d \t r_f = %d \t ag_m = %d \t r_m = %d \n", ag_f,r_f,ag_m,r_m);
                    //fflush(stdout);
                    if(patch_f==patch_m)
                    {
                        param->unbalanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] = TIME_STEP*patch[patch_f].n_population->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f]*param->c_per_gender_within_patch[FEMALE][ag_f]*param->relative_number_partnerships_per_risk[r_f]*param->p_age_per_gender[FEMALE][ag_f][ag_m]*param->p_risk_per_gender[FEMALE][r_f][r_m];
                    }else
                    {
                        param->unbalanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] = TIME_STEP*patch[patch_f].n_population->pop_size_per_gender_age_risk[FEMALE][ag_f][r_f]*param->c_per_gender_between_patches[FEMALE][ag_f]*param->relative_number_partnerships_per_risk[r_f]*param->p_age_per_gender[FEMALE][ag_f][ag_m]*param->p_risk_per_gender[FEMALE][r_f][r_m];
                    }
                    /* In English this means that the number of new partnerships within a time step between female ages ag_f risk r_f and males aged ag_m risk r_m,
                     * as desired by the females
                     * is the following product:
                     * TIME STEP x Number females in that group x rate of new partners for females in that age group x relative rate for females in that risk group
                     *           x proportion of new male partners within a certain age category ag_m
                     *           x proportion of new male partners within a certain risk category r_m */
                }
            }

        }
    }
}


void calcul_n_new_partners_m_to_f(patch_struct *patch, parameters *param, int patch_f, int patch_m){
    int ag_f, ag_m ; /* indexes of ages for females and males */
    int r_f, r_m ; /* indexes of risk groups for females and males */

    /* Fills in p_risk_m */
    calcul_p_risk(MALE, param->p_risk_per_gender[MALE], patch[patch_m].n_population_stratified, param);

    /* Standardizes the relative number of new partners per risk group */
    //standardize_relative_number_partnerships_per_risk(MALE, param->xi_per_gender[MALE], patch[patch_m].n_population_stratified, param);


    /* Calculates the "desired" number of new partners for males per age and risk group */
    for(ag_f=0 ; ag_f<N_AGE ; ag_f++)
    {
        for(r_f=0 ; r_f<N_RISK ; r_f++)
        {
            for(ag_m=0 ; ag_m<N_AGE ; ag_m++)
            {
                for(r_m=0 ; r_m<N_RISK ; r_m++)
                {
                    if(patch_f==patch_m)
                    {
                        param->unbalanced_nb_m_to_f[ag_m][r_m][ag_f][r_f] = TIME_STEP*patch[patch_m].n_population->pop_size_per_gender_age_risk[MALE][ag_m][r_m]*param->c_per_gender_within_patch[MALE][ag_m]*param->relative_number_partnerships_per_risk[r_m]*param->p_age_per_gender[MALE][ag_m][ag_f]*param->p_risk_per_gender[MALE][r_m][r_f];
                    }else
                    {
                        param->unbalanced_nb_m_to_f[ag_m][r_m][ag_f][r_f] = TIME_STEP*patch[patch_m].n_population->pop_size_per_gender_age_risk[MALE][ag_m][r_m]*param->c_per_gender_between_patches[MALE][ag_m]*param->relative_number_partnerships_per_risk[r_m]*param->p_age_per_gender[MALE][ag_m][ag_f]*param->p_risk_per_gender[MALE][r_m][r_f];
                    }
                }
            }

        }
    }
}

/* this fills in balanced_nb_f_to_m with the balanced number of
 * new partnerships between a female of age ag_f, risk r_f and a male of age ag_m, risk r_m over a time unit */
void balance_contacts_arithmetic(parameters *param)
{
    int ag_f, ag_m,r_f, r_m;
    double tmp;
    for(ag_f=0 ; ag_f<N_AGE ; ag_f++)
    {
        for(r_f=0 ; r_f<N_RISK ; r_f++)
        {
            for(ag_m=0 ; ag_m<N_AGE ; ag_m++)
            {
                for(r_m=0 ; r_m<N_RISK ; r_m++)
                {
                    tmp = (1-param->prop_compromise_from_males) * param->unbalanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] + param->prop_compromise_from_males * param->unbalanced_nb_m_to_f[ag_m][r_m][ag_f][r_f];
                    param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] = (long int) floor (tmp);
                    param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m] += gsl_ran_bernoulli (rng, tmp - floor(tmp)); /* Bernoulli trial with probability tmp - floor(tmp) */
                }
            }
        }
    }
}

int are_in_partnership(individual *indiv1, individual *indiv2)
{
    int res = 0;
    int i;

    for(i=0 ; i<indiv1->n_partners ; i++)
    {
        if(indiv1->partner_pairs[i]->ptr[1-indiv1->gender]->id==indiv2->id)
        {
            res = 1;
            break;
        }
    }

    return(res);
}

int has_free_partnership(individual *indiv)
{
    return(indiv->n_partners < indiv->max_n_partners);
}

/* is idx_new_partners_m[current_idx] equal to any of the other idx_new_partners_m[k]? */
int is_already_selected(long *idx_new_partners_m, long current_idx, long n_partners_m)
{
    long k;

    k = current_idx;
    if(k < n_partners_m - 1)
    {
        do
        {
            k++;
        }while((k < n_partners_m - 1) && (idx_new_partners_m[current_idx] != idx_new_partners_m[k]));
    }else
    {
        k = 0;
    }

    if(idx_new_partners_m[current_idx] != idx_new_partners_m[k])
    {
        if(current_idx >0)
        {
            k = current_idx - 1;
            while( (idx_new_partners_m[current_idx] != idx_new_partners_m[k]) && (k > 0))
            {
                k--;
            }
        }
    }

    if(idx_new_partners_m[current_idx] != idx_new_partners_m[k])
    {
        //printf("Answer = %d\n",0);
        //fflush(stdout);
        return(0);
    }else
    {
        //printf("Answer = %d\n",1);
        //fflush(stdout);
        return(1);
    }

}


void copy_array_long(long *dest, long *orig, long size){
    /* Copy an array of long integers
    
    Arguments
    ---------
    dest, orig : pointer to array of long integers
        orig is the array to be copied and dest is the location of the array in which to copy orig.
    
    size : long
        Length of the array orig to be copied
    */
    int k;
    
    for(k = 0; k < size ; k++){
        dest[k] = orig[k];
    }
}


int is_serodiscordant(partnership *pair){
    /* Determine if a partnership is serodiscordant or not
    
    This function returns 1 if partnership is serodiscordant (ie 1 HIV+ and 1 HIV- partner) and 0
    if seroconcordant (i.e. +/+ or -/-).  The check finds the sum and product of the two
    serostatuses.  If sum is >0 (so at least one is HIV+) and the product is zero (so at least one
    is HIV-) then the partnership is serodiscordant. 
    
    Arguments
    ---------
    pair : pointer to a partnership struct
        The partnership to test for serodiscordance.  
    
    Returns
    -------
    int
        1 if the partnership is serodiscordant; zero otherwise
    */
    return(
        ((pair->ptr[0]->HIV_status * pair->ptr[1]->HIV_status) == 0) && 
        ((pair->ptr[0]->HIV_status + pair->ptr[1]->HIV_status) > 0)
    );
}


int compare_longs (const void *a, const void *b){
    const long *da = (const long *) a;
    const long *db = (const long *) b;

    return (*da > *db) - (*da < *db);
}


void get_setting(patch_struct *patch){
   /* Convert community ID to country setting.  
    
    Query the community_id attribute of a patch structure and set the country_setting attribute
    according to where the community occurs (Zambia communities are <= 12; otherwise South Africa).
    
    The country settings (ZAMBIA and SOUTH_AFRICA) are integers defined within constants.h.  
    
    Arguments
    --------
    patch : pointer to array of patch_struct objects
    
    Returns
    -------
    Nothing; sets country_setting attribute of the patch object to either macro ZAMBIA or
    SOUTH_AFRICA (as defined within constants.h)
    */
    
    int p;
    for(p = 0; p < NPATCHES; p++){
        /* First 12 clusters are Zambia: */
        if(patch->community_id <= 12){
            
            if(VERBOSE_OUTPUT == 1){
                printf("Setting: Zambia\n");
            }
            patch[p].country_setting = ZAMBIA;
        }else{
            if(VERBOSE_OUTPUT == 1){
                printf("Setting: South Africa\n");
            }
            patch[p].country_setting = SOUTH_AFRICA;
        }
    }
}


/* Given the current time (discretised into year and timestep - ie t=year+t_step*TIME_STEP). */
int get_chips_round(parameters *param, int year, int t_step){
    //printf("Calling chips_sampling_frame() year=%i t_step=%i. Exiting\n",year,t_step);
    /* proposed modification:
    int round = 0;
    while (round < NCHIPSROUNDS){
        if ((year==param->CHIPS_START_YEAR[0] && t_step>=param->CHIPS_START_TIMESTEP[0]) || (year>param->CHIPS_START_YEAR[0] && year<param->CHIPS_START_YEAR[1])
                (year==param->CHIPS_START_YEAR[1] && t_step<param->CHIPS_START_TIMESTEP[1]))
            return r;
    }    */
    int chips_round;
    int saved_chips_round = CHIPSNOTRUNNING;   /* Default (-2) indicates we are outside of chips rounds - e.g. before trial, or during a break between rounds. */

    if ((year>param->CHIPS_END_YEAR[NCHIPSROUNDS-1])
            || (year==param->CHIPS_END_YEAR[NCHIPSROUNDS-1] && (t_step>param->CHIPS_END_TIMESTEP[NCHIPSROUNDS-1])))
        saved_chips_round = CHIPSROUNDPOSTTRIAL;  /* This indicates we are post-trial. */
    else{
        for (chips_round=0;chips_round<NCHIPSROUNDS;chips_round++){
            if ((year==param->CHIPS_START_YEAR[chips_round] && t_step>=param->CHIPS_START_TIMESTEP[chips_round] && (year<param->CHIPS_END_YEAR[chips_round] || (year==param->CHIPS_END_YEAR[chips_round] && t_step<=param->CHIPS_END_TIMESTEP[chips_round])))
                    || (year>param->CHIPS_START_YEAR[chips_round] && year<param->CHIPS_END_YEAR[chips_round])
                    || (year>param->CHIPS_START_YEAR[chips_round] && year==param->CHIPS_END_YEAR[chips_round] && t_step<=param->CHIPS_END_TIMESTEP[chips_round])){
                //printf("For round %i YEAR=%i %i T_STEP=%i %i\n",chips_round,param->CHIPS_START_YEAR[chips_round],param->CHIPS_END_YEAR[chips_round],param->CHIPS_START_TIMESTEP[chips_round],param->CHIPS_END_TIMESTEP[chips_round]);
                //printf("Chips round found to be %i at t=%i",chips_round,year);
                saved_chips_round = chips_round;
                break;
            }
        }
        if(saved_chips_round==CHIPSNOTRUNNING)
            if ((year>param->CHIPS_START_YEAR[0]) ||  ((year==param->CHIPS_START_YEAR[0]) && (t_step>=param->CHIPS_START_TIMESTEP[0])) ){
                printf("Warning - in get_chips_round() CHiPs does not seem to be running in year=%i, t_step=%i, CHIPS_START_YEAR=%i CHIPS_START_TIMESTEP=%i. Check that CHiPs were on between-rounds training.\n",year,t_step,param->CHIPS_START_YEAR[0],param->CHIPS_START_TIMESTEP[0]);
                fflush(stdout);
                //printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                //fflush(stdout);
                //exit(1);
            }

    }

    return saved_chips_round;
}

/* Function returns 1 if this is the start of a chips round, otherwise returns 0. 
Allows us to carry out chips sampling once per round. */
int is_start_of_chips_round(parameters *param, int year, int t_step, int trial_arm){
    
    int chips_round;
    int is_start_of_chips = 0;
    
    if(trial_arm == ARM_A || trial_arm == ARM_B){
        
        for(chips_round = 0; chips_round < NCHIPSROUNDS; chips_round++){
            if(
                (year == param->CHIPS_START_YEAR[chips_round] &&
                t_step == param->CHIPS_START_TIMESTEP[chips_round]) || 
                    
                (year > param->CHIPS_START_YEAR[NCHIPSROUNDS-1] &&
                t_step == param->CHIPS_START_TIMESTEP_POSTTRIAL) ||
                    
                 (year == param->CHIPS_END_YEAR[NCHIPSROUNDS-1] &&
                t_step == param->CHIPS_START_TIMESTEP_POSTTRIAL &&
                t_step > param->CHIPS_END_TIMESTEP[NCHIPSROUNDS-1])
            ){
                is_start_of_chips = 1;
                break;
            }
        }
    }else{
        /* The condition inside this loop is probably redundant (legacy code) but keep anyway*/
        if(year >= T_ROLLOUT_CHIPS_EVERYWHERE){
            // Allow rollout under CF scenarios
            if(ALLOW_COUNTERFACTUAL_ROLLOUT == 1){
                is_start_of_chips = 1;
            }else{
                if((year > param->CHIPS_END_YEAR[NCHIPSROUNDS-1] &&
                t_step == param->CHIPS_START_TIMESTEP_POSTTRIAL) ||
            
                (year == param->CHIPS_END_YEAR[NCHIPSROUNDS-1] &&
                t_step == param->CHIPS_START_TIMESTEP_POSTTRIAL &&
                t_step > param->CHIPS_END_TIMESTEP[NCHIPSROUNDS-1])
                ){
                    is_start_of_chips = 1;
                }
            }
        }
    }
    return is_start_of_chips;
}


/* Currently we have 7 options - directory for parameters,
 number of runs, start run number,
 * By setting rng_seed_offset to be an integer which is not 0, this changes the C seed from param->rng_seed to param->rng_seed+rng_seed_offset.
 * */
void parse_command_line_arguments(int argc, char **argv, int *n_runs, int *i_startrun, 
    int *n_startrun, int *is_counterfactual, int *rng_seed_offset, int *rng_seed_offset_PC){
    // Argc is the number of arguments. Note that the first argument is the filename so argc is always >=1.
    if (argc == 1 || argc > 9){
        //printf("ERROR: Arguments must include the directory of the parameter files. Optional second argument: number of runs (default=1000). Optional third argument: run numebr to start at (default is 1). Optional fourth argument: directory where results are written. \nExiting\n");
        printf("ERROR: Arguments must include the directory of the parameter files. \nOptional second argument: number of runs (default=1000). \nOptional third argument: is this run a counterfactual run (ie. all patches are assumed to have arm C testing/VMMC until 2020, after which all patches with patch number p>0 become popart-like while patch 0 remains without any UTT). Optional fourth argument: run numebr to start at (default is 1).\n Optional fifth argument: directory where results are written. Optional sixth argument: random seed offset (for stochasticity of impact). Optional seventh argument: random seed offset PC (for stochasticity of PC recruitment).\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (argc>2){
        /* The first argument of strtol() is the input string, the second is where it puts any characters in the first argument, the third is the base of the number to be returned. */
        *n_runs = strtol(argv[2],NULL,10);
    }
    else
        *n_runs = 1000;

    if (argc>3){
        *is_counterfactual = strtol(argv[3],NULL,10);
        /* Check that this only takes values 0 or 1: */
        if (!(*is_counterfactual==NOT_COUNTERFACTUAL_RUN || *is_counterfactual==IS_COUNTERFACTUAL_RUN)){
            printf("ERROR: 3rd argument (is_counterfactual) must be 0 (not counterfactual) or 1 (counterfactual) only.\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    else
        *is_counterfactual = NOT_COUNTERFACTUAL_RUN; /* Default value. */


    if (argc>4)
        *i_startrun = strtol(argv[4],NULL,10);
    else
        *i_startrun = 1;

    if (argc>5)
        *n_startrun = strtol(argv[5],NULL,10);
	else
        *n_startrun = (*n_runs) - (*i_startrun - 1);

    if (argc>7)
        *rng_seed_offset = strtol(argv[7],NULL,10);
    else
        *rng_seed_offset = 0;

    if (argc>8)
        *rng_seed_offset_PC = strtol(argv[8],NULL,10);
    else
        *rng_seed_offset_PC = 0;
    printf("Offset for PC is %i\n",*rng_seed_offset_PC);

    if (*i_startrun>*n_runs){
        printf("ERROR: 4th argument (i_startrun) must be <= 2nd argument (n_runs).\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}


void get_IBM_code_version(char *version, int stringlength){
   /* Return version number to the character string `version`.  
    
    Arguments
    ---------
    version : pointer to char array
        variable which will store the version number
    
    stringlength : int
        length of memory in `version` (ie max number of characters in the version number).
    
    Returns
    ------
    Nothing; appends version number to the input char array `version`
    */
    strncpy(version, "V1.2", stringlength);
    return;
}


void add_slash(char *filename){
    /* Add OS-depdendent slash ('/' or '\') to a character array.  
    
    Add a space that is OS-dependent, i.e. works for both Linux/OS X ("/") or Windows ("\\")
    respectively.
    
    If filename is empty then this corresponds to "current directory" so adding a slash would move
    to root directory instead - which would be very bad!!!  So only add "/" or "\" if the filename
    is not currently empty.
    
    Arguments
    ---------
    filename : pointer to a character array
        File name of interest that needs a slash added to it.  
    
    Returns
    -------
    Nothing; function adjusts the input character array
    */
    
    
    /* Check if filename is empty */
    if(filename[0] != '\0'){
        if(THISOS == 0){
            /* Add  \ when running on Windows OS. Need extra "\" to modify escape character. */
            strcat(filename,"\\"); 
        }else{
             /* Add / when running on Mac/Linux. */
            strcat(filename,"/");
        }
    }
}


void join_strings_with_check(char *dest, char *src, int dest_array_size, char *error_message){
    /* Join the string src to dest. ie dest = dest + source.
    
    Arguments
    ---------
    dest, src : pointer to char arrays
        Strings to be joined
    
    dest_array_size : int
        Size of the char array dest (so check memory).
    
    error_message : pointer to char
        Printed if the array dest is not big enough - should be of the form "SOURCE and DESTINATION
        in function f()."  Where SOURCE, DESTINATION and f() are changed as appropriate.  
    
    
    Note that the arguments dest and src are ordered in join_strings_with_check() to be as for
    strcat(dest, src).
    */
    
    /* -1 as the final character is '\0'. */
    //printf("Check %i, %lu\n",dest_array_size,sizeof(*dest));
    if((strlen(src) + strlen(dest)) < dest_array_size - 1){
        strcat(dest, src);
    }else{
        printf("Error joining arrays %s. Please make destination array larger.\n). Exiting\n",
            error_message);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
}


void make_output_label_struct(file_label_struct *file_labels, long python_rng_seed, int i_run, int rng_seed_offset, int rng_seed_offset_PC, patch_struct *patch, int is_counterfactual, int PCdata){
    int p;

    char version[21];  /* 20 is arbitrary but note we have to change the 20 in the call to get_IBM_code_version() below if we change this. */
    //char temp[LONGSTRINGLENGTH];

    char community_info[LONGSTRINGLENGTH];
    char clusternumber[20];    /* Again this is arbitrary - assume we never have stupidly long filenames! */
    char patchinfo[20];    /* Again this is arbitrary - assume we never have stupidly long filenames! */
    char run_info_ending[LONGSTRINGLENGTH];   /* Again this is arbitrary - assume we never have stupidly long filenames! */

    memset(clusternumber, '\0', sizeof(clusternumber));
    memset(community_info, '\0', sizeof(community_info));
    memset(patchinfo, '\0', sizeof(patchinfo));
    memset(run_info_ending, '\0', sizeof(run_info_ending));


    /* First sort out the parts of the filename which are independent of patch number: */
    get_IBM_code_version(version,20);

    if (is_counterfactual==NOT_COUNTERFACTUAL_RUN){
        if (PCdata==0)
            sprintf(run_info_ending,"_Rand%li_Run%i_PCseed%i_%i.csv",python_rng_seed,i_run+1,rng_seed_offset_PC,rng_seed_offset);
        else
            sprintf(run_info_ending,"_Rand%li_Run%i_PCseed%i_%i_PConly.csv",python_rng_seed,i_run+1,rng_seed_offset_PC,rng_seed_offset);
    }

    //Ann_out_CL01_SA_A_V1.2_patch0_RandX_RUn1_PCseed4_10.csv
    else if (is_counterfactual==IS_COUNTERFACTUAL_RUN){
        if (PCdata==0)
            sprintf(run_info_ending,"_Rand%li_Run%i_PCseed%i_%i_CF.csv",python_rng_seed,i_run+1,rng_seed_offset_PC,rng_seed_offset);
        else
            sprintf(run_info_ending,"_Rand%li_Run%i_PCseed%i_%i_PConly_CF.csv",python_rng_seed,i_run+1,rng_seed_offset_PC,rng_seed_offset);

        printf("Making this run a counterfactual with label=%s\n",run_info_ending);
    }
    else{
        printf("Unknown value for is_counterfactual=%i in make_output_label_struct(). Exiting\n",is_counterfactual);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    for (p=0;p<NPATCHES;p++){
        /* Add cluster number. Pad with a zero if needed: */
        if (patch[p].community_id<=9)
            sprintf(clusternumber,"_CL0%i",patch[p].community_id);
        else if ((patch[p].community_id<=1000000) && (patch[p].community_id>0))
            sprintf(clusternumber,"_CL%i",patch[p].community_id);
        else{
            printf("Error - community_id is too large %i. Trial communities go from 1-21. Exiting. \n",patch[p].community_id);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        /* community_info is the _COUNTRY_ARM_CODEVERSION_ data. */
        if (patch[p].country_setting==ZAMBIA)
            strcpy(community_info,"_Za_");
        else
            strcpy(community_info,"_SA_");
        if (patch[p].trial_arm==ARM_A)
            strcat(community_info,"A_");
        else if (patch[p].trial_arm==ARM_B)
            strcat(community_info,"B_");
        else
            strcat(community_info,"C_");
        strcat(community_info,version);

        sprintf(patchinfo,"_patch%i",p);


        if (strlen(clusternumber)>LONGSTRINGLENGTH-1){
            printf("Error - need to increase size of file_labels->filename_label_bypatch[p] as clusternumber[] in make_output_label_struct() is too long. Exiting\n");
        }
        strcpy(file_labels->filename_label_bypatch[p],clusternumber);

        join_strings_with_check(file_labels->filename_label_bypatch[p], community_info, LONGSTRINGLENGTH, "community_info and filename_label_bypatch[] in make_output_label_struct()");
        join_strings_with_check(file_labels->filename_label_bypatch[p], patchinfo, LONGSTRINGLENGTH, "patchinfo and filename_label_bypatch[] in make_output_label_struct()");
        join_strings_with_check(file_labels->filename_label_bypatch[p], run_info_ending, LONGSTRINGLENGTH, "run_info_ending and filename_label_bypatch[] in make_output_label_struct()");

        if (is_counterfactual==IS_COUNTERFACTUAL_RUN){
            printf("Making this run a counterfactual with label=%s\n",file_labels->filename_label_bypatch[p]);
        }
        /* These are just based on data from patch 0. */
        if (p==0){
            strcpy(file_labels->filename_label_allpatches,community_info);
            join_strings_with_check(file_labels->filename_label_allpatches, run_info_ending, LONGSTRINGLENGTH, "run_info_ending and filename_label_allpatches in make_output_label_struct()");

            /* Files e.g. duration_partnership_within_high_high_CL01_Za_A_V1.2_Rand10_Run2_0.csv. */
            strcpy(file_labels->filename_label_allpatches_witharm_communityno,clusternumber);
            join_strings_with_check(file_labels->filename_label_allpatches_witharm_communityno, community_info, LONGSTRINGLENGTH-1, "community_info and filename_label_allpatches_witharm_communityno in make_output_label_struct()");
            join_strings_with_check(file_labels->filename_label_allpatches_witharm_communityno, run_info_ending, LONGSTRINGLENGTH, "run_info_ending and filename_label_allpatches_witharm_communityno in make_output_label_struct()");

        }

    }

}


void make_filenames_for_struct(file_label_struct *file_labels, 
    file_struct *file_data_store, char *output_file_directory){
    /*
    Set up file names for output files.  
    
    This funciton has many calls to concatenate_filename() which concatenates the 2nd, 3rd, and 4th
    input arguments into a single string which is stored in the first argument (in the order
    2nd-4th-3rd args).  
    
    Arguments
    ---------
    file_labels : pointer to file_label_struct structure
        Structure that holds all the file names of files to be written.  There are three character 
        arrays within this struct:
        1) filename_label_bypatch, for patch-specific data (most commonly used file-type)
        2) filename_label_allpatches, for data that's summarised across both patches
        3) filename_label_allpatches_witharm_communityno, files used for debugging 
        (not entirely sure, WP)
    file_data_store : pointer to a file_struct structure
        Structure where the different file names will be stored.  Usually, when this function is 
        called these will all be empty and will be populated by calls to concatenate_filename().  
    output_file_directory : pointer to a char
        Output directory of interest, e.g. "./data/SAMPLED_PARAMETERS/PARAMS_COMMUNITY5/Output"
    
    Returns
    -------
    Nothing; file names are added to attributes of the `file_data_store` structure.  
    */
    
    int p;
    for(p = 0; p < NPATCHES; p++){
        concatenate_filename(file_data_store->filename_annual_output[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Annual_outputs");
        concatenate_filename(file_data_store->filename_timestep_output[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Timestep_outputs");
        concatenate_filename(file_data_store->filename_timestep_output_PConly[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Timestep_outputs_PConly");
        concatenate_filename(file_data_store->filename_timestep_age_output[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Timestep_age_outputs");
        concatenate_filename(file_data_store->filename_timestep_age_output_PConly[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Timestep_age_outputs_PConly");
        concatenate_filename(file_data_store->filename_annual_partnership_output[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Annual_partnerships_outputs");
        concatenate_filename(file_data_store->filename_PC_data[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "PC_outputs");
        concatenate_filename(file_data_store->filename_chipsannual_data[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "CHIPS_outputs_annual");
        concatenate_filename(file_data_store->filename_chipsvisit_data[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "CHIPS_outputs_visit");
        concatenate_filename(file_data_store->filename_debughivduration[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "DEBUG_HIVduration");
        concatenate_filename(file_data_store->filename_debughivduration_km[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "DEBUG_HIVduration_KM");
        concatenate_filename(file_data_store->filename_debughivcd4_after_seroconversion[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "DEBUG_HIV_CD4_after_seroconversion");
        concatenate_filename(file_data_store->filename_debuginitial_spvl_distribution[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "DEBUG_HIV_initialSPVLdistribution");
        concatenate_filename(file_data_store->filename_debug_artpopulation[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "DEBUG_ART_population");
        concatenate_filename(file_data_store->filename_debug_hivpopulation[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "DEBUG_HIVstates_population");
        concatenate_filename(file_data_store->filename_debug_agedistribution[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Age_distribution_check");
        concatenate_filename(file_data_store->filename_debug_one_yearage_dist_includingkids[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "OneYearAgeGp");

        concatenate_filename(file_data_store->filename_distr_n_lifetime_partners[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Distr_n_lifetime_partners");
        concatenate_filename(file_data_store->filename_distr_n_partners_lastyear[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Distr_n_partners_lastyear");
            
        if(WRITE_COST_EFFECTIVENESS_OUTPUT == 1){
        concatenate_filename(file_data_store->filename_cost_effectiveness_output[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "Cost_effectiveness");
        }
        
        if(WRITE_TREATS_OUTPUT == 1){
        concatenate_filename(file_data_store->filename_treats_output[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "TREATS");
        }
        
        if(WRITE_ART_STATUS_BY_AGE_SEX == 1){
        concatenate_filename(file_data_store->filename_art_status_by_age_sex_output[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            "ART_status_by_age_sex");
        }
    }
    /* Some filenames are special - only made for one run. */
    concatenate_filename(file_data_store->filename_debug_nnewadults_ndeaths_file,
        output_file_directory, file_labels->filename_label_bypatch[0],
        "NBirthsNNewAdultsNdeaths_Run");

    /* Set up partnership filenames: */
    concatenate_filename(file_data_store->filename_DUR_BETWEEN_HIGHHIGH,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "duration_partnership_between_high_high");
    concatenate_filename(file_data_store->filename_DUR_BETWEEN_MEDMED,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "duration_partnership_between_med_med");
    concatenate_filename(file_data_store->filename_DUR_BETWEEN_LOWLOW,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "duration_partnership_between_low_low");
    concatenate_filename(file_data_store->filename_DUR_WITHIN_HIGHHIGH,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "duration_partnership_within_high_high");
    concatenate_filename(file_data_store->filename_DUR_WITHIN_MEDMED,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "duration_partnership_within_med_med");
    concatenate_filename(file_data_store->filename_DUR_WITHIN_LOWLOW,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "duration_partnership_within_low_low");

    concatenate_filename(file_data_store->filename_age_assortativity,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "Age_assortativity_at_partnership_formation");
    concatenate_filename(file_data_store->filename_age_assortativity_cross_sectional,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "Age_assortativity_cross_sectional");
    concatenate_filename(file_data_store->filename_risk_assortativity,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "Risk_assortativity_at_partnership_formation");
    concatenate_filename(file_data_store->filename_risk_assortativity_cross_sectional,
        output_file_directory, file_labels->filename_label_allpatches_witharm_communityno,
        "Risk_assortativity_cross_sectional");


    /* Set up phylogenetic output filenames: */
    char phylo_trans_filename_temp[100];
    char phylo_indiv_filename_temp[100];
    sprintf(phylo_trans_filename_temp, "phylogenetic_transmission");
    sprintf(phylo_indiv_filename_temp, "phylogenetic_individualdata");
    char temp[20]; /* Size is arbitrary. */
    memset(temp, '\0', sizeof(temp));
    
    /* this adds the date of the simulation to the output filename - olli wanted this. */
    if (WRITE_PHYLOGENETICS_OUTPUT == 2){
        time_t tnow = time(NULL);
        struct tm tm = *localtime(&tnow);
        sprintf(temp,"_%02d%02d%d", tm.tm_mday, tm.tm_mon + 1, tm.tm_year - 100);
        join_strings_with_check(phylo_trans_filename_temp, temp, 100, 
            "phylo_trans_filename_temp and temp in make_filenames_for_struct()");
        join_strings_with_check(phylo_indiv_filename_temp, temp, 100, 
            "phylo_trans_filename_temp and temp in make_filenames_for_struct()");
    }
    for(p=0;p<NPATCHES;p++){
        concatenate_filename(file_data_store->filename_phylogenetic_transmission[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            phylo_trans_filename_temp);
        concatenate_filename(file_data_store->filename_phylogenetic_individualdata[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
            phylo_indiv_filename_temp);

        concatenate_filename(file_data_store->filename_hivsurvival_individualdata[p],
            output_file_directory, file_labels->filename_label_bypatch[p],
        "HIVsurvival_individualdata");
    }
    /* Output for hazards - note that we only want from patch 0. */
    concatenate_filename(file_data_store->filename_hazard_output,
        output_file_directory, file_labels->filename_label_bypatch[0],
        "Hazards");
}


void concatenate_filename(char *output_filename, char *output_file_directory, char *label_p, 
    char *file_tag){
   /* Join the file output directory, file_tag and label_p into a single string
    and store within `output_filename`.  This will eventually be used as the complete path 
    for writing files to disk.  
    
    The input argument `output_filename` is initially blanked out by this function.  
    
    "file_tag" is the start name of the file. Function *SHOULD* be designed to stop buffer overflow.
    
    Arguments
    ---------
    
    output_filename : pointer to a char
        Character array in which the complete file name will be stored (is initially blanked by this
        function).  This is usually an empty character string.  
    output_file_directory : pointer to a char
        Directory to where the file of interest should be saved.  
        e.g. "./data/SAMPLED_PARAMETERS/PARAMS_COMMUNITY5/Output"
    label_p : pointer to a char
        Final identifier for the file, e.g. "_CL05_Za_A_V1.2_patch0_Rand10_Run1_PCseed0_0.csv"
    file_tag : pointer to a char
        Tag for the file in question, e.g. "Annual_outputs"
    
    Returns
    -------
    Nothing; copies a string to `output_filename`.  
    */
    
    /* Ensure that output_filename is blank. */
    memset(output_filename, '\0', LONGSTRINGLENGTH*sizeof(char));
    
    /* If output_filename is not long enough then exit - so prevent buffer overflow.
     * We add 2 as we need a termination character '\0' and the slash from add_slash(). */
    if (LONGSTRINGLENGTH < strlen(output_file_directory) + strlen(label_p) +strlen(file_tag) + 2){
        
        printf("Inadequate memory allocated for file_tag=%s %lu %lu %lu %lu \nExiting\n",
            file_tag, sizeof(output_filename), strlen(output_file_directory), 
            strlen(label_p), strlen(file_tag));
        
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    /* Copy the output directory to output_filename. Note that while 'strcpy()' is regarded as 
    more prone to buffer overflow, we check this beforehand. */
    strcpy(output_filename, output_file_directory);
    
    /* Add a / or \ as needed if working in directory other than current local dir. */
    add_slash(output_filename);
    strcat(output_filename, file_tag);
    strcat(output_filename, label_p);
}


/* Note that this function can be called multiple times whenever we want to take snapshots at several times but want to
 * have flexibility as to the times when we store data - so we create several files each of which has the corresponding
 * year in the title. */
void make_filenames_for_snapshot(char *output_filename, char *output_file_directory, file_label_struct *file_labels, int year, int p, char *file_tag){
    char yearstring[6]; /* Each year has 4 digits + one char for "_" + one extra char for terminating '\0'. */
    sprintf(yearstring,"_%i",year);
    memset(output_filename, '\0', LONGSTRINGLENGTH*sizeof(char));            /* Ensure that output_filename is blank. */
    /* If output_filename is not long enough then exit - so prevent buffer overflow.
     * We add 2 as we need a termination character '\0' and the slash from add_slash(). */
    if (LONGSTRINGLENGTH < strlen(output_file_directory) + strlen(yearstring) + strlen(file_labels->filename_label_bypatch[p]) +strlen(file_tag) + 2){
        printf("Inadequate memory allocated for file_tag=%s %lu %lu %lu \nExiting\n",file_labels->filename_label_bypatch[p],sizeof(output_filename),strlen(output_file_directory),strlen(file_tag));
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* Copy the output directory to output_filename. Note that while 'strcpy()' is regarded as more prone to buffer overflow, we check this beforehand. */
    strcpy(output_filename,output_file_directory);
    add_slash(output_filename); /* Adds a / or \ as needed if working in directory other than current local dir. */
    strcat(output_filename,file_tag);
    strcat(output_filename,yearstring);
    strcat(output_filename,file_labels->filename_label_bypatch[p]);


}

/* Merge this with others above??? */
void make_calibration_output_filename(char *output_filename, char *output_file_directory, long python_rng_seed, patch_struct *patch, int p, int rng_seed_offset, int rng_seed_offset_PC, int is_counterfactual){

    char temp[LONGSTRINGLENGTH];     /*  This is arbitrary - assume we never have stupidly long filenames! */


    char version[21];     /* 20 is arbitrary but note we have to change the 20 in the call to get_IBM_code_version() below if we change this. */
    get_IBM_code_version(version,20);

    /* Ensure everything is blank: */
    memset(temp, '\0', sizeof(temp));

    strncpy(output_filename,output_file_directory,LONGSTRINGLENGTH);

    add_slash(output_filename); /* Adds a / or \ as needed if working in directory other than current local dir. */

    /* Add cluster number. Pad with a zero if needed: */
    if (patch[p].community_id<=9)
        sprintf(temp,"Calibration_output_CL0%i",patch[p].community_id);
    else
        sprintf(temp,"Calibration_output_CL%i",patch[p].community_id);

    strcat(output_filename,temp);

    if (patch[p].country_setting==ZAMBIA)
        strcat(output_filename,"_Za_");
    else
        strcat(output_filename,"_SA_");

    if (patch[p].trial_arm==ARM_A)
        strcat(output_filename,"A_");
    else if (patch[p].trial_arm==ARM_B)
        strcat(output_filename,"B_");
    else
        strcat(output_filename,"C_");

    strcat(output_filename,version);

    sprintf(temp,"_patch%i_",p);
    strcat(output_filename,temp);


    if (is_counterfactual==NOT_COUNTERFACTUAL_RUN)
        sprintf(temp,"Rand%li_PCseed%i_%i.csv",python_rng_seed,rng_seed_offset,rng_seed_offset_PC);
    else
        sprintf(temp,"Rand%li_PCseed%i_%i_CF.csv",python_rng_seed,rng_seed_offset,rng_seed_offset_PC);
    strcat(output_filename,temp);
    //printf("output filename = %s\n",output_filename);
}

/* Takes an existing string, counts the number of commas in it, and adds extra commas for padding
 * so that the total number of data points including blanks is NDATA (so NDATA-1 commas).
 * This is used for example for the calibration outputs
 * to ensure that the csv file is more readable. */
void add_commas_to_calibration_output(char *output_string,int NDATA){
    int i = 0;
    /* Find the first comma (if any) in the string: */
    char *pch=strchr(output_string,',');
    /* Count the rest of the commas, one by one. */
    while (pch!=NULL) {
        i++;
        pch=strchr(pch+1,',');
    }

    //int k;
    while(i<NDATA){
        //for(k=0;k<(NDATA-i-1);k++)
        strcat(output_string,",");
        i++;
    }
}


void print_param_struct(parameters *param){
    int g,ag,bg,icd4,spvl,r,i,y;
    int ac;

    printf("------------------------------------------------------------------------------------\n");
    printf("---------------------------- Printing param structure ------------------------------\n");
    printf("------------------------------------------------------------------------------------\n");

    for (y=0; y<N_UNPD_TIMEPOINTS; y++)
        for (ag = 0; ag<N_AGE_UNPD_FERTILITY; ag++)
            printf("param->fertility_rate_by_age[%i][%i]=%lg\n",ag,y,param->fertility_rate_by_age[ag][y]);
    for (g=0; g<N_GENDER;g++)
        for (ag = 0; ag<N_AGE_UNPD_MORTALITY; ag++)
            printf("param->mortality_rate_by_gender_age_intercept[%i][%i]=%lg\n",g,ag,param->mortality_rate_by_gender_age_intercept[g][ag]);
    for (g=0; g<N_GENDER;g++)
        for (ag = 0; ag<N_AGE_UNPD_MORTALITY; ag++)
            printf("param->mortality_rate_by_gender_age_slope[%i][%i]=%lg\n",g,ag,param->mortality_rate_by_gender_age_slope[g][ag]);
    printf("param->sex_ratio=%lg\n",param->sex_ratio);
    printf("param->p_child_circ=%lg\n",param->p_child_circ);
    printf("param->eff_circ_vmmc=%lg\n",param->eff_circ_vmmc);
    printf("param->eff_circ_tmc=%lg\n",param->eff_circ_tmc);
    printf("param->rr_circ_unhealed=%lg\n",param->rr_circ_unhealed);
    printf("param->t0_pmtct=%lg\n",param->t0_pmtct);
    printf("param->t50_pmtct=%lg\n",param->t50_pmtct);
    printf("param->average_log_viral_load=%lg\n",param->average_log_viral_load);
    printf("param->average_annual_hazard=%lg\n",param->average_annual_hazard);
    printf("param->RRacute_trans=%lg\n",param->RRacute_trans);
    printf("param->RRmale_to_female_trans=%lg\n",param->RRmale_to_female_trans);
    printf("param->RRCD4[0],param->RRCD4[1],param->RRCD4[2],param->RRCD4[3]=%lg %lg %lg %lg\n",param->RRCD4[0],param->RRCD4[1],param->RRCD4[2],param->RRCD4[3]);
    
    printf("param->SPVL_beta_k=%lg \n",param->SPVL_beta_k);
    printf("param->SPVL_beta_50=%lg \n",param->SPVL_beta_50);
    printf("param->RR_ART_INITIAL=%lg\n",param->RR_ART_INITIAL);
    printf("param->RR_ART_VS=%lg\n",param->RR_ART_VS);
    printf("param->RR_ART_VU=%lg\n",param->RR_ART_VU);
    printf("param->min_dur_acute,param->max_dur_acute=%lg %lg\n",param->min_dur_acute,param->max_dur_acute);
    for (spvl=0; spvl<NSPVL; spvl++)
        printf("param->p_initial_cd4_gt500[spvl]=%lg\n",param->p_initial_cd4_gt500[spvl]);
    for (spvl=0; spvl<NSPVL; spvl++)
        printf("param->p_initial_cd4_350_500[spvl]=%lg\n",param->p_initial_cd4_350_500[spvl]);
    for (spvl=0; spvl<NSPVL; spvl++)
        printf("param->p_initial_cd4_200_350[spvl]=%lg\n",param->p_initial_cd4_200_350[spvl]);
    for (spvl=0; spvl<NSPVL; spvl++)
        printf("param->p_initial_cd4_lt200[spvl]=%lg\n",param->p_initial_cd4_lt200[spvl]);

    printf("param->initial_SPVL_mu=%lg\n",param->initial_SPVL_mu);
    printf("param->initial_SPVL_sigma=%lg\n",param->initial_SPVL_sigma);
    printf("param->SPVL_sigma_M=%lg\n",param->SPVL_sigma_M);
    printf("param->SPVL_sigma_E=%lg\n",param->SPVL_sigma_E);

    for (icd4=0; icd4<NCD4; icd4++)
        for (spvl=0; spvl<NSPVL; spvl++)
            printf("param->time_hiv_event[icd4][spvl]= %lg\n",param->time_hiv_event[icd4][spvl]);
    printf("param->factor_for_slower_progression_ART_VU=%lg\n",param->factor_for_slower_progression_ART_VU);
    printf("param->assortativity=%lg\n",param->assortativity);
    printf("param->prop_compromise_from_males=%lg\n",param->prop_compromise_from_males);

    for (ag=0; ag<N_AGE; ag++)
        printf("param->c_per_gender_within_patch[FEMALE][ag]=%lg\n",param->c_per_gender_within_patch[FEMALE][ag]);
    for (ag=0; ag<N_AGE; ag++)
        printf("param->c_per_gender_within_patch[MALE][ag]=%lg\n",param->c_per_gender_within_patch[MALE][ag]);

    printf("param->rel_rate_partnership_formation_between_patches=%lg\n",param->rel_rate_partnership_formation_between_patches);

    for (r=0; r<N_RISK; r++)
        printf("param->relative_number_partnerships_per_risk[r]=%lg\n",param->relative_number_partnerships_per_risk[r]);
    for (g=0; g<N_GENDER; g++)
        for (ag=0; ag<N_AGE; ag++)
            for (bg=0; bg<N_AGE; bg++)
                printf("param->p_age_per_gender[g][ag][bg]=%lg\n",param->p_age_per_gender[g][ag][bg]);

    for(r=0 ; r<N_RISK ; r++)
        printf("param->max_n_part_noage[r]=%d\n",param->max_n_part_noage[r]);
    for(r=0 ; r<N_RISK ; r++)
        printf("param->breakup_scale_lambda_within_patch[r]=%lg\n",param->breakup_scale_lambda_within_patch[r]);
    for(r=0 ; r<N_RISK ; r++)
        printf("param->breakup_scale_lambda_between_patch[r]=%lg\n",param->breakup_scale_lambda_between_patch[r]);
    for(r=0 ; r<N_RISK ; r++)
        printf("param->breakup_shape_k[r]=%lg\n",param->breakup_shape_k[r]);

    printf("param->start_time_hiv=%lg\n",param->start_time_hiv);
    printf("param->start_time_simul=%i\n",param->start_time_simul);
    printf("param->end_time_simul=%i\n",param->end_time_simul);
    printf("param->COUNTRY_HIV_TEST_START=%lg\n",param->COUNTRY_HIV_TEST_START);
    printf("param->COUNTRY_ART_START=%lg\n",param->COUNTRY_ART_START);
    printf("param->COUNTRY_CD4_350_START=%lg\n",param->COUNTRY_CD4_350_START);
    printf("param->COUNTRY_CD4_500_START=%lg\n",param->COUNTRY_CD4_500_START);
    printf("param->COUNTRY_VMMC_START=%lg\n",param->COUNTRY_VMMC_START);
    for (i=0; i<NCHIPSROUNDS; i++){
        printf("param->CHIPS_START_YEAR[%i]=%i\n",i,param->CHIPS_START_YEAR[i]);
        printf("param->CHIPS_END_YEAR[%i]=%i\n",i,param->CHIPS_END_YEAR[i]);
        printf("param->CHIPS_START_TIMESTEP[%i]=%i\n",i,param->CHIPS_START_TIMESTEP[i]);
        printf("param->CHIPS_END_TIMESTEP[%i]=%i\n",i,param->CHIPS_END_TIMESTEP[i]);
    }

    printf("param->DHS_params->NDHSROUNDS=%i\n",param->DHS_params->NDHSROUNDS);
    for (i=0; i<param->DHS_params->NDHSROUNDS; i++)
        printf("param->DHS_params->DHS_YEAR[%i]=%i\n",i,param->DHS_params->DHS_YEAR[i]);


    printf("param->time_to_background_HIVtestNOW=%lg\n",param->time_to_background_HIVtestNOW);
    printf("param->time_to_background_HIVtest_maxval=%lg\n",param->time_to_background_HIVtest_maxval);
    printf("param->time_to_background_HIVtest_exponent=%lg\n",param->time_to_background_HIVtest_exponent);
    printf("param->time_to_background_HIVtest_midpoint=%lg\n",param->time_to_background_HIVtest_midpoint);

    printf("param->p_collect_hiv_test_results_cd4_over200=%lg\n",param->p_collect_hiv_test_results_cd4_over200);
    printf("param->p_collect_hiv_test_results_cd4_under200=%lg\n",param->p_collect_hiv_test_results_cd4_under200);

    printf("param->p_collect_cd4_test_results_cd4_nonpopart=%lg\n",param->p_collect_cd4_test_results_cd4_nonpopart);
    printf("param->p_collect_cd4_test_results_cd4_popartYEAR1=%lg\n",param->p_collect_cd4_test_results_cd4_popartYEAR1);
    printf("param->p_collect_cd4_test_results_cd4_popartYEAR2onwards=%lg\n",param->p_collect_cd4_test_results_cd4_popartYEAR2onwards);
    for (icd4=0; icd4<NCD4; icd4++)
        printf("param->p_dies_earlyart_cd4[icd4]=%lg\n",param->p_dies_earlyart_cd4[icd4]);
    printf("param->p_leaves_earlyart_cd4_over200_if_not_die_early=%lg\n",param->p_leaves_earlyart_cd4_over200_if_not_die_early);
    printf("param->p_leaves_earlyart_cd4_under200_if_not_die_early=%lg\n",param->p_leaves_earlyart_cd4_under200_if_not_die_early);
    printf("param->p_becomes_vs_after_earlyart_if_not_die_early_or_leave=%lg\n",param->p_becomes_vs_after_earlyart_if_not_die_early_or_leave);
    printf("param->p_stays_virally_suppressed=%lg\n",param->p_stays_virally_suppressed);
    printf("param->p_stops_virally_suppressed=%lg\n",param->p_stops_virally_suppressed);
    printf("param->p_vu_becomes_virally_suppressed=%lg\n",param->p_vu_becomes_virally_suppressed);
    printf("param->t_earlyart_dropout_min[NOTPOPART]=%lg\n",param->t_earlyart_dropout_min[NOTPOPART]);
    printf("param->t_earlyart_dropout_min[POPART]=%lg\n",param->t_earlyart_dropout_min[POPART]);
    printf("param->t_earlyart_dropout_range[NOTPOPART]=%lg\n",param->t_earlyart_dropout_range[NOTPOPART]);
    printf("param->t_earlyart_dropout_range[POPART]=%lg\n",param->t_earlyart_dropout_range[POPART]);
    printf("param->t_dies_earlyart_min[NOTPOPART]=%lg\n",param->t_dies_earlyart_min[NOTPOPART]);
    printf("param->t_dies_earlyart_min[POPART]=%lg\n",param->t_dies_earlyart_min[POPART]);
    printf("param->t_dies_earlyart_range[NOTPOPART]=%lg\n",param->t_dies_earlyart_range[NOTPOPART]);
    printf("param->t_dies_earlyart_range[POPART]=%lg\n",param->t_dies_earlyart_range[POPART]);
    printf("param->t_end_early_art=%lg\n",param->t_end_early_art);
    printf("param->t_cd4_retest_min[NOTPOPART]=%lg\n",param->t_cd4_retest_min[NOTPOPART]);
    printf("param->t_cd4_retest_min[POPART]=%lg\n",param->t_cd4_retest_min[POPART]);
    printf("param->t_cd4_retest_range[NOTPOPART]=%lg\n",param->t_cd4_retest_range[NOTPOPART]);
    printf("param->t_cd4_retest_range[POPART]=%lg\n",param->t_cd4_retest_range[POPART]);
    printf("param->t_cd4_whenartfirstavail_min=%lg\n",param->t_cd4_whenartfirstavail_min);
    printf("param->t_cd4_whenartfirstavail_range=%lg\n",param->t_cd4_whenartfirstavail_range);
    printf("param->t_delay_hivtest_to_cd4test_min[NOTPOPART]=%lg\n",param->t_delay_hivtest_to_cd4test_min[NOTPOPART]);
    printf("param->t_delay_hivtest_to_cd4test_min[POPART]=%lg\n",param->t_delay_hivtest_to_cd4test_min[POPART]);
    printf("param->t_delay_hivtest_to_cd4test_range[NOTPOPART]=%lg\n",param->t_delay_hivtest_to_cd4test_range[NOTPOPART]);
    printf("param->t_delay_hivtest_to_cd4test_range[POPART]=%lg\n",param->t_delay_hivtest_to_cd4test_range[POPART]);
    
    // Exponential version
    printf("param->t_start_art_mean_non_popart=%lg\n",param->t_start_art_mean_non_popart);
    for (i=0; i<NCHIPSROUNDS; i++)
        printf("param->n_time_periods_art_popart_per_round[%i]=%d\n",i,param->n_time_periods_art_popart_per_round[i]);
    
    int iround;
    for(iround = 0; iround < NCHIPSROUNDS; iround++){
        for(i = 0; i < param->n_time_periods_art_popart_per_round[iround]; i++){
            printf("param->t_start_art_mean_fast_popart[%d][%d]=%lg\n",
                iround, i, param->t_start_art_mean_fast_popart[iround][i]);
            printf("param->t_start_art_mean_slow_popart[%d][%d]=%lg\n",
                iround, i, param->t_start_art_mean_slow_popart[iround][i]);
            printf("param->p_start_art_mean_fast_popart[%d][%d]=%lg\n",
                iround, i, param->p_start_art_mean_fast_popart[iround][i]);
        }
    }
    
    printf("param->t_end_vs_becomevu_min[NOTPOPART]=%lg\n",param->t_end_vs_becomevu_min[NOTPOPART]);
    printf("param->t_end_vs_becomevu_min[POPART]=%lg\n",param->t_end_vs_becomevu_min[POPART]);
    printf("param->t_end_vs_becomevu_range[NOTPOPART]=%lg\n",param->t_end_vs_becomevu_range[NOTPOPART]);
    printf("param->t_end_vs_becomevu_range[POPART]=%lg\n",param->t_end_vs_becomevu_range[POPART]);
    printf("param->t_end_vs_dropout_min[NOTPOPART]=%lg\n",param->t_end_vs_dropout_min[NOTPOPART]);
    printf("param->t_end_vs_dropout_min[POPART]=%lg\n",param->t_end_vs_dropout_min[POPART]);
    printf("param->t_end_vs_dropout_range[NOTPOPART]=%lg\n",param->t_end_vs_dropout_range[NOTPOPART]);
    printf("param->t_end_vs_dropout_range[POPART]=%lg\n",param->t_end_vs_dropout_range[POPART]);
    printf("param->t_end_vu_becomevs_min[NOTPOPART]=%lg\n",param->t_end_vu_becomevs_min[NOTPOPART]);
    printf("param->t_end_vu_becomevs_min[POPART]=%lg\n",param->t_end_vu_becomevs_min[POPART]);
    printf("param->t_end_vu_becomevs_range[NOTPOPART]=%lg\n",param->t_end_vu_becomevs_range[NOTPOPART]);
    printf("param->t_end_vu_becomevs_range[POPART]=%lg\n",param->t_end_vu_becomevs_range[POPART]);
    printf("param->t_end_vu_dropout_min[NOTPOPART]=%lg\n",param->t_end_vu_dropout_min[NOTPOPART]);
    printf("param->t_end_vu_dropout_min[POPART]=%lg\n",param->t_end_vu_dropout_min[POPART]);
    printf("param->t_end_vu_dropout_range[NOTPOPART]=%lg\n",param->t_end_vu_dropout_range[NOTPOPART]);
    printf("param->t_end_vu_dropout_range[POPART]=%lg\n",param->t_end_vu_dropout_range[POPART]);
    for (i=0; i<NCHIPSROUNDS; i++)
        printf("param->p_popart_to_cascade[%i]=%lg\n",i,param->p_popart_to_cascade[i]);
    printf("param->p_circ_nonpopart=%lg\n",param->p_circ_nonpopart);
    for (i=0; i<NCHIPSROUNDS; i++)
        printf("param->p_circ_popart[%i]=%lg\n",i,param->p_circ_popart[i]);
    printf("param->t_get_vmmc_min[NOTPOPART]=%lg\n",param->t_get_vmmc_min[NOTPOPART]);
    printf("param->t_get_vmmc_min[POPART]=%lg\n",param->t_get_vmmc_min[POPART]);
    printf("param->t_get_vmmc_range[NOTPOPART]=%lg\n",param->t_get_vmmc_range[NOTPOPART]);
    printf("param->t_get_vmmc_range[POPART]=%lg\n",param->t_get_vmmc_range[POPART]);
    printf("param->t_vmmc_healing=%lg\n",param->t_vmmc_healing);

    for (g=0; g<N_GENDER;g++)
        for (ac = 0; ac<(MAX_AGE-AGE_CHIPS+1); ac++)
            for (i=0; i<NCHIPSROUNDS; i++)
                printf("param->prop_tested_by_chips_in_round[%i][%i][%i]=%lg\n",g,ac,i,param->chips_params->prop_tested_by_chips_in_round[g][ac][i]);

    printf("param->initial_population_size=%lg\n",param->initial_population_size);
    for (ag=0; ag<N_AGE; ag++)
        printf("param->initial_prop_age[ag]=%lg\n",param->initial_prop_age[ag]);
    for (g=0; g<N_GENDER; g++)
        for (r=0; r<N_RISK; r++)
            printf("param->initial_prop_gender_risk[g][r]=%lg\n",param->initial_prop_gender_risk[g][r]);
    for (g=0; g<N_GENDER; g++)
        for (r=0; r<N_RISK; r++)
            printf("param->initial_prop_infected_gender_risk[g][r]=%lg\n",param->initial_prop_infected_gender_risk[g][r]);

    printf("------------------------------------------------------------------------------------\n");
    printf("------------------------------------------------------------------------------------\n");
    printf("------------------------------------------------------------------------------------\n");

}


void check_if_parameters_plausible(parameters *param){
   /* Function goes through each input parameter and checks it against a pre-defined range. If any parameter doesn't
    * fit in the given range, then print an error message and exit.
    * Note that this function can prevent you from deliberately setting crazy values for debugging, so can disable by
    * setting CHECKPARAMS in constants.h to 0.
    */
    
    
    int g, ag, bg, icd4, jcd4, spvl, r, a_unpd, y, iquarter;
    int ac, chips_timestep, chips_round, dhs_round;
    double temp;
    
    for(y = 0; y < N_UNPD_TIMEPOINTS; y++){
        for (a_unpd = 0; a_unpd < N_AGE_UNPD_FERTILITY; a_unpd++){
            if (param->fertility_rate_by_age[a_unpd][y] < 1e-5 || param->fertility_rate_by_age[a_unpd][y] > 1.0){
                printf("Error: param->total_fertility_rate is outside expected range [0.1,8]\nExiting\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }

    if(param->sex_ratio < 0.4 || param->sex_ratio > 0.6){
        printf("Error: param->sex_ratio is outside expected range [0.4,0.6]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_child_circ<0 || param->p_child_circ>1){
        printf("Error: param->p_child_circ is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->eff_circ_vmmc<0.0 || param->eff_circ_vmmc>1.0){
        printf("Error: param->eff_circ_vmmc is outside expected range [0.0,1.0]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->eff_circ_tmc<0.0 || param->eff_circ_tmc>param->eff_circ_vmmc){
        printf("Error: param->eff_circ_tmc is outside expected range [0.0,param->eff_circ_vmmc]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->rr_circ_unhealed<0 || param->rr_circ_unhealed>3){
        printf("Error: param->rr_circ_unhealed is outside expected range [0,3]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->t0_pmtct<1990 || param->t0_pmtct>2015){
        printf("Error: param->t0_pmtct is outside expected range [1990,2015]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t50_pmtct<1990 || param->t50_pmtct>2015){
        printf("Error: param->t50_pmtct is outside expected range [1990,2015]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->t0_pmtct>param->t50_pmtct){
        printf("Error: param->t0_pmtct is bigger than param->t50_pmtct.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->average_log_viral_load<3 || param->average_log_viral_load>6){
        printf("Error: param->average_log_viral_load is outside expected range [3,6]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->average_annual_hazard<0.001 || param->average_annual_hazard>0.6){
        printf("Error: param->average_annual_hazard is outside expected range [0.01,0.6]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->RRacute_trans<1 || param->RRacute_trans>50){
        printf("Error: param->RRacute_trans is outside expected range [1,50]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->RRmale_to_female_trans<0.33 || param->RRmale_to_female_trans>5){
        printf("Error: param->RRmale_to_female_trans is outside expected range [0.33,5]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    for (icd4=0;icd4<NCD4;icd4++){
        if (param->RRCD4[icd4]<0 || param->RRCD4[icd4]>10){
            printf("Error: param->RRCD4[icd4] is outside expected range [0,10]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }

    if (param->SPVL_beta_k<0.9 || param->SPVL_beta_k>1.2){
        printf("Error: param->SPVL_beta_k =%lf is outside expected range [0.9,1.2]\nExiting\n",param->SPVL_beta_k);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->SPVL_beta_50<5000 || param->SPVL_beta_50>50000){
        printf("Error: param->SPVL_beta_50 is outside expected range [5000,50000]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->RR_ART_INITIAL<0 || param->RR_ART_INITIAL>1){
        printf("Error: param->RR_ART_INITIAL is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->RR_ART_VS<0 || param->RR_ART_VS>1){
        printf("Error: param->RR_ART_VS is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->RR_ART_VU<0 || param->RR_ART_VU>1){
        printf("Error: param->RR_ART_VU is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->RR_ART_INITIAL < param->RR_ART_VS){
        printf("Error: param->RR_ART_VS is bigger than param->RR_ART_INITIAL.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->RR_ART_VU < param->RR_ART_VS){
        printf("Error: param->RR_ART_VS is bigger than param->RR_ART_VU.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }


    if (param->min_dur_acute<0 || param->min_dur_acute>1){
        printf("Error: param->min_dur_acute is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->max_dur_acute<0 || param->max_dur_acute>2){
        printf("Error: param->max_dur_acute is outside expected range [0,2]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->min_dur_acute > param->max_dur_acute){
        printf("Error: param->min_dur_acute is bigger than param->max_dur_acute.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    for (spvl=0; spvl<NSPVL; spvl++){
        if (param->p_initial_cd4_gt500[spvl]<0 || param->p_initial_cd4_gt500[spvl]>1){
            printf("Error:param->p_initial_cd4_gt500[spvl] is outside expected range [0,1]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        if (param->p_initial_cd4_350_500[spvl]<0 || param->p_initial_cd4_350_500[spvl]>1){
            printf("Error:param->p_initial_cd4_350_500[spvl] is outside expected range [0,1]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        if (param->p_initial_cd4_200_350[spvl]<0 || param->p_initial_cd4_200_350[spvl]>1){
            printf("Error:param->p_initial_cd4_200_350[spvl] is outside expected range [0,1]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        if (param->p_initial_cd4_lt200[spvl]<0 || param->p_initial_cd4_lt200[spvl]>1){
            printf("Error:param->p_initial_cd4_lt200[spvl] is outside expected range [0,1]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }


    if (param->initial_SPVL_mu<0 || param->initial_SPVL_mu>5){
        printf("Error: param->initial_SPVL_mu is outside expected range [0,5]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->initial_SPVL_sigma<0 || param->initial_SPVL_sigma>2){
        printf("Error: param->initial_SPVL_sigma is outside expected range [0,2]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->SPVL_sigma_M<0.0 || param->SPVL_sigma_M>1.0){
        printf("Error: param->SPVL_sigma_M is outside expected range [0.01,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->SPVL_sigma_E<0.01 || param->SPVL_sigma_E>0.9){
        printf("Error: param->SPVL_sigma_E is outside expected range [0.01,0.9]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* Note that the first index icd4 is the true cd4, and the other index jcd4 is the measured cd4. */
    for (icd4=0; icd4<NCD4; icd4++){
        for (jcd4=0; jcd4<NCD4; jcd4++){
            /* temp converts back from cumulative probability to normal probability. */
            if (jcd4==0)
                temp = param->cumulative_p_misclassify_cd4[icd4][jcd4];
            else
                temp = param->cumulative_p_misclassify_cd4[icd4][jcd4] - param->cumulative_p_misclassify_cd4[icd4][jcd4-1];

            /* Probability of correct classification - should be high: */
            if (icd4==jcd4){
                if (temp<0.8||temp>1.0){
                    printf("Error:param->cumulative_p_misclassify_cd4[icd4][icd4] is outside expected range [0.8,1.0]\nExiting\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
            }
            /* Probability of incorrect classification - should be low: */
            else{
                if(temp < 0.0 || temp > 0.2){
                    printf("Error:param->cumulative_p_misclassify_cd4[icd4][jcd4] is outside expected range [0.0,0.2]\nExiting\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
            }
        }
    }


    for (icd4=0; icd4<NCD4; icd4++){
        for (spvl=0; spvl<NSPVL; spvl++){
            if (param->time_hiv_event[icd4][spvl]<0.0 || param->time_hiv_event[icd4][spvl]>18){
                printf("Error: param->time_hiv_event[icd4=%i][spvl=%i][0]=%f is outside expected range [0.25,18]\nExiting\n",icd4,spvl,param->time_hiv_event[icd4][spvl]);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }

    for (icd4=0; icd4<NCD4; icd4++){
        if (param->CoxPH_SPVL_RR_CD4progression[icd4]<1.0 || param->CoxPH_SPVL_RR_CD4progression[icd4]>3.5){
            printf("Error: param->CoxPH_SPVL_RR_CD4progression[icd4=%i]=%f is outside expected range [1.0,3.5]\nExiting\n",icd4,param->CoxPH_SPVL_RR_CD4progression[icd4]);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }


    if (param->factor_for_slower_progression_ART_VU<0 || param->factor_for_slower_progression_ART_VU>5){
        printf("Error: param->factor_for_slower_progression_ART_VU is outside expected range [0,5]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }



    if (param->prop_compromise_from_males<0 || param->prop_compromise_from_males>1){
        printf("Error:param->prop_compromise_from_males is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->assortativity<0 || param->assortativity>1){
        printf("Error:param->assortativity is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    for (g=0; g<N_GENDER; g++){
        for (ag=0; ag<N_AGE; ag++){
            for (bg=0; bg<N_AGE; bg++){
                if (param->p_age_per_gender[g][ag][bg]<0 || param->p_age_per_gender[g][ag][bg]>1){
                    printf("Error:param->p_age_per_gender[g][ag][bg]=%lf is outside expected range [0,1]\nExiting\n",param->p_age_per_gender[g][ag][bg]);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
            }
        }
    }


    if (param->time_to_background_HIVtestNOW<0.2 || param->time_to_background_HIVtestNOW>6.5){
        printf("Error:param->time_to_background_HIVtestNOW is outside expected range [0.2,6.5] years\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->time_to_background_HIVtest_maxval<1 || param->time_to_background_HIVtest_maxval>20){
        printf("Error:param->time_to_background_HIVtest_maxval is outside expected range [1,20] years\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->time_to_background_HIVtest_exponent<0 || param->time_to_background_HIVtest_exponent>10){
        printf("Error:param->time_to_background_HIVtest_exponent is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->time_to_background_HIVtest_midpoint<2 || param->time_to_background_HIVtest_midpoint>15){
        printf("Error:param->time_to_background_HIVtest_midpoint is outside expected range [2,15] years\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->HIV_rapid_test_sensitivity_CHIPS<0.5 || param->HIV_rapid_test_sensitivity_CHIPS>1){
        printf("Error:param->HIV_rapid_test_sensitivity_CHIPS is outside expected range [0.5,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->p_collect_hiv_test_results_cd4_over200<0 || param->p_collect_hiv_test_results_cd4_over200>1){
        printf("Error:param->p_collect_hiv_test_results_cd4_over200 is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_collect_hiv_test_results_cd4_under200<0 || param->p_collect_hiv_test_results_cd4_under200>1){
        printf("Error:param->p_collect_hiv_test_results_cd4_under200 is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_collect_hiv_test_results_cd4_over200 > param->p_collect_hiv_test_results_cd4_under200){
        printf("Error: param->p_collect_hiv_test_results_cd4_over200 is bigger than param->p_collect_hiv_test_results_cd4_under200.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->p_collect_cd4_test_results_cd4_nonpopart<0 || param->p_collect_cd4_test_results_cd4_nonpopart>1){
        printf("Error:param->p_collect_cd4_test_results_cd4_nonpopart is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_collect_cd4_test_results_cd4_popartYEAR1<0 || param->p_collect_cd4_test_results_cd4_popartYEAR1>1){
        printf("Error:param->p_collect_cd4_test_results_cd4_popartYEAR1 is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_collect_cd4_test_results_cd4_popartYEAR2onwards<0 || param->p_collect_cd4_test_results_cd4_popartYEAR2onwards>1){
        printf("Error:param->p_collect_cd4_test_results_cd4_popartYEAR2onwards is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->p_leaves_earlyart_cd4_over200_if_not_die_early<0 || param->p_leaves_earlyart_cd4_over200_if_not_die_early>1){
        printf("Error:param->p_leaves_earlyart_cd4_over200_if_not_die_early is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_leaves_earlyart_cd4_under200_if_not_die_early<0 || param->p_leaves_earlyart_cd4_under200_if_not_die_early>1){
        printf("Error:param->p_leaves_earlyart_cd4_under200_if_not_die_early is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->p_leaves_earlyart_cd4_under200_if_not_die_early > param->p_leaves_earlyart_cd4_over200_if_not_die_early){
        printf("Error: param->p_leaves_earlyart_cd4_under200_if_not_die_early is bigger than param->p_leaves_earlyart_cd4_over200_if_not_die_early.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->p_becomes_vs_after_earlyart_if_not_die_early_or_leave<0 || param->p_becomes_vs_after_earlyart_if_not_die_early_or_leave>1){
        printf("Error:param->p_becomes_vs_after_earlyart_if_not_die_early_or_leave is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_stays_virally_suppressed<0 || param->p_stays_virally_suppressed>1){
        printf("Error:param->p_stays_virally_suppressed is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_stops_virally_suppressed<0 || param->p_stops_virally_suppressed>1){
        printf("Error:param->p_stops_virally_suppressed is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->p_vu_becomes_virally_suppressed<0 || param->p_vu_becomes_virally_suppressed>1){
        printf("Error:param->p_vu_becomes_virally_suppressed is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    for (chips_round=0; chips_round<NCHIPSROUNDS; chips_round++){
        if (param->p_popart_to_cascade[chips_round]<0 || param->p_popart_to_cascade[chips_round]>1){
            printf("Error:param->p_popart_to_cascade[%i] is outside expected range [0,1]\nExiting\n",chips_round);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    if (param->p_circ_nonpopart<0 || param->p_circ_nonpopart>1){
        printf("Error:param->p_circ_nonpopart is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->chips_params->n_timesteps_per_round_posttrial<24 || param->chips_params->n_timesteps_per_round_posttrial>=96){
        printf("Error:param->chips_params->n_timesteps_per_round_posttrial =%d is outside expected range [24,95] weeks (0.5-2yrs)\n. Note that if you want a CHiPs round to last 2 years or more then you must increase the size of param->chips_params->prop_tested_by_chips_per_timestep. Exiting\n",param->chips_params->n_timesteps_per_round_posttrial);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    for (g=0; g<N_GENDER;g++){
        for (ac = 0; ac<(MAX_AGE-AGE_CHIPS+1); ac++){
            if (param->chips_params->prop_tested_by_chips_in_round_posttrial[g][ac]<0 || param->chips_params->prop_tested_by_chips_in_round_posttrial[g][ac]>1){
                printf("Error:param->chips_params->prop_tested_by_chips_in_round_posttrial[%i][%i] = %lf is outside expected range [0,1]\nExiting\n",g,ac,param->chips_params->prop_tested_by_chips_in_round_posttrial[g][ac]);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }

    for(chips_timestep=0; chips_timestep<param->chips_params->n_timesteps_per_round_posttrial;chips_timestep++){
        for (g=0; g<N_GENDER;g++){
            for (ac = 0; ac<(MAX_AGE-AGE_CHIPS+1); ac++){
                if (param->chips_params->prop_tested_by_chips_per_timestep_posttrial[g][ac][chips_timestep]<0 || param->chips_params->prop_tested_by_chips_per_timestep_posttrial[g][ac][chips_timestep]>1){
                    printf("Error:param->chips_params->prop_tested_by_chips_per_timestep_posttrial[%i][%i][%i] = %lf is outside expected range [0,1]\nExiting\n",g,ac,chips_timestep,param->chips_params->prop_tested_by_chips_per_timestep_posttrial[g][ac][chips_timestep]);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
            }
        }
    }


    for (chips_round=0; chips_round<NCHIPSROUNDS; chips_round++){
        if (param->chips_params->n_timesteps_per_round[chips_round]<1 || param->chips_params->n_timesteps_per_round[chips_round]>=96){
            printf("Error:param->chips_params->n_timesteps_per_round[%d] =%d is outside expected range [12,95] weeks (0.25-2yrs)\n. Note that if you want a CHiPs round to last 2 years or more then you must increase the size of param->chips_params->prop_tested_by_chips_per_timestep. Exiting\n",chips_round,param->chips_params->n_timesteps_per_round[chips_round]);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        if (param->p_circ_popart[chips_round]<0 || param->p_circ_popart[chips_round]>1){
            printf("Error:param->p_circ_popart[%i] is outside expected range [0,1]\nExiting\n",chips_round);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }


        for (g=0; g<N_GENDER;g++){
            for (ac = 0; ac<(MAX_AGE-AGE_CHIPS+1); ac++){
                //double prop_tested_by_chips_in_round[N_GENDER][MAX_AGE-AGE_CHIPS+1][NCHIPSROUNDS];
                if (param->chips_params->prop_tested_by_chips_in_round[g][ac][chips_round]<0 || param->chips_params->prop_tested_by_chips_in_round[g][ac][chips_round]>1){
                    printf("Error:param->chips_params->prop_tested_by_chips_in_round[%i][%i][%i] = %lf is outside expected range [0,1]\nExiting\n",g,ac,chips_round,param->chips_params->prop_tested_by_chips_in_round[g][ac][chips_round]);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
            }
        }

        for(chips_timestep=0; chips_timestep<param->chips_params->n_timesteps_per_round[chips_round];chips_timestep++){
            for (g=0; g<N_GENDER;g++){
                for (ac = 0; ac<(MAX_AGE-AGE_CHIPS+1); ac++){
                    if (param->chips_params->prop_tested_by_chips_per_timestep[g][ac][chips_timestep][chips_round]<0 || param->chips_params->prop_tested_by_chips_per_timestep[g][ac][chips_timestep][chips_round]>1){
                        printf("Error:param->chips_params->prop_tested_by_chips_per_timestep[%i][%i][%i][%i] = %lf is outside expected range [0,1]\nExiting\n",g,ac,chips_timestep,chips_round,param->chips_params->prop_tested_by_chips_per_timestep[g][ac][chips_timestep][chips_round]);
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                }
            }
        }
    }

    for (ag=0; ag<N_AGE; ag++){
        if (param->initial_prop_age[ag]<0 || param->initial_prop_age[ag]>1){
            printf("Error:param->initial_prop_age[ag] is outside expected range [0,1]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    for (g=0; g<N_GENDER; g++){
        for (r=0; r<N_RISK; r++){
            if (param->initial_prop_gender_risk[g][r]<0 || param->initial_prop_gender_risk[g][r]>1){
                printf("Error:param->initial_prop_gender_risk[g][r] is outside expected range [0,1]\nExiting\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }
    for (g=0; g<N_GENDER; g++){
        for (r=0; r<N_RISK; r++){
            if (param->initial_prop_infected_gender_risk[g][r]<0 || param->initial_prop_infected_gender_risk[g][r]>1){
                printf("Error:param->initial_prop_infected_gender_risk[g][r] is outside expected range [0,1]\nExiting\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }
    for (icd4=0; icd4<NCD4; icd4++){
        if (param->p_dies_earlyart_cd4[icd4]<0 || param->p_dies_earlyart_cd4[icd4]>1){
            printf("Error:param->p_dies_earlyart_cd4[icd4] is outside expected range [0,1]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }

    if (param->t_earlyart_dropout_min[NOTPOPART]<0 || param->t_earlyart_dropout_min[NOTPOPART]>1){
        printf("Error:param->t_earlyart_dropout_min[NOTPOPART] is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_earlyart_dropout_min[POPART]<0 || param->t_earlyart_dropout_min[POPART]>1){
        printf("Error:param->t_earlyart_dropout_min[POPART] is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_earlyart_dropout_range[NOTPOPART]<0 || param->t_earlyart_dropout_range[NOTPOPART]>1){
        printf("Error:param->t_earlyart_dropout_range[NOTPOPART] is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_earlyart_dropout_range[POPART]<0 || param->t_earlyart_dropout_range[POPART]>1){
        printf("Error:param->t_earlyart_dropout_range[POPART] is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_dies_earlyart_min[NOTPOPART]<0 || param->t_dies_earlyart_min[NOTPOPART]>1){
        printf("Error:param->t_dies_earlyart_min[NOTPOPART] is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_dies_earlyart_min[POPART]<0 || param->t_dies_earlyart_min[POPART]>1){
        printf("Error:param->t_dies_earlyart_min[POPART] is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_dies_earlyart_range[NOTPOPART]<0 || param->t_dies_earlyart_range[NOTPOPART]>1){
        printf("Error:param->t_dies_earlyart_range[NOTPOPART] is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_dies_earlyart_range[POPART]<0 || param->t_dies_earlyart_range[POPART]>1){
        printf("Error:param->t_dies_earlyart_range[POPART] is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_early_art<0 || param->t_end_early_art>1){
        printf("Error:param->t_end_early_art is outside expected range [0,1]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_cd4_retest_min[NOTPOPART]<0 || param->t_cd4_retest_min[NOTPOPART]>10){
        printf("Error:param->t_cd4_retest_min[NOTPOPART] is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_cd4_retest_min[POPART]<0 || param->t_cd4_retest_min[POPART]>10){
        printf("Error:param->t_cd4_retest_min[POPART] is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_cd4_retest_range[NOTPOPART]<0 || param->t_cd4_retest_range[NOTPOPART]>10){
        printf("Error:param->t_cd4_retest_range[NOTPOPART] is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_cd4_retest_range[POPART]<0 || param->t_cd4_retest_range[POPART]>10){
        printf("Error:param->t_cd4_retest_range[POPART] is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_cd4_whenartfirstavail_min<0 || param->t_cd4_whenartfirstavail_min>5){
        printf("Error:param->t_cd4_whenartfirstavail_min is outside expected range [0,5]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_cd4_whenartfirstavail_range<0 || param->t_cd4_whenartfirstavail_range>10){
        printf("Error:param->t_cd4_whenartfirstavail_range is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_delay_hivtest_to_cd4test_min[NOTPOPART]<0 || param->t_delay_hivtest_to_cd4test_min[NOTPOPART]>5){
        printf("Error:param->t_delay_hivtest_to_cd4test_min[NOTPOPART] is outside expected range [0,5]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_delay_hivtest_to_cd4test_min[POPART]<0 || param->t_delay_hivtest_to_cd4test_min[POPART]>5){
        printf("Error:param->t_delay_hivtest_to_cd4test_min[POPART] is outside expected range [0,5]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_delay_hivtest_to_cd4test_range[NOTPOPART]<0 || param->t_delay_hivtest_to_cd4test_range[NOTPOPART]>10){
        printf("Error:param->t_delay_hivtest_to_cd4test_range[NOTPOPART] is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_delay_hivtest_to_cd4test_range[POPART]<0 || param->t_delay_hivtest_to_cd4test_range[POPART]>10){
        printf("Error:param->t_delay_hivtest_to_cd4test_range[POPART] is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    // Exponential version

    if (param->t_start_art_mean_non_popart<0 || param->t_start_art_mean_non_popart>10){
        printf("Error:param->t_start_art_mean_non_popart is outside expected range [0,10]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    // Biexponential version

    for(chips_round=0 ; chips_round<NCHIPSROUNDS ; chips_round++){
        if (param->n_time_periods_art_popart_per_round[chips_round]<1 || param->n_time_periods_art_popart_per_round[chips_round]>MAX_N_TIME_PERIODS_PER_ROUND){
            printf("Error:param->p_start_art_mean_fast_popart is outside expected range [1,%d]\nExiting\n",MAX_N_TIME_PERIODS_PER_ROUND);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    
    int iround;
    for(iround = 0; iround < NCHIPSROUNDS; iround++){
        for(iquarter = 0; iquarter < param->n_time_periods_art_popart_per_round[iround]; iquarter++){

            if(param->t_start_art_mean_fast_popart[iround][iquarter]<0 || param->t_start_art_mean_fast_popart[iround][iquarter]>10){
                printf("Error:param->t_start_art_mean_fast_popart[%d][%d] is outside expected range [0,10]\nExiting\n", iround, iquarter);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            if(param->t_start_art_mean_slow_popart[iround][iquarter]<0 || param->t_start_art_mean_slow_popart[iround][iquarter]>16){
                printf("Error:param->t_start_art_mean_slow_popart[%d][%d] is outside expected range [0,10]\nExiting\n", iround, iquarter);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            if(param->t_start_art_mean_slow_popart[iround][iquarter] < param->t_start_art_mean_fast_popart[iround][iquarter]){
                printf("Error:param->t_start_art_mean_slow_popart[%d][%d] is shorter than param->t_start_art_mean_fast_popart[%d][%d] (%lg<%lg)\nExiting\n",iround, iquarter,iround, iquarter,param->t_start_art_mean_slow_popart[iround][iquarter],param->t_start_art_mean_fast_popart[iround][iquarter]);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            if (param->p_start_art_mean_fast_popart[iround][iquarter] < 0 || param->p_start_art_mean_fast_popart[iround][iquarter] > 1){
                printf("Error:param->p_start_art_mean_fast_popart[%d][%d] is outside expected range [0,1]\nExiting\n",iround, iquarter);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }
    
    if (param->t_end_vs_becomevu_min[NOTPOPART]<0 || param->t_end_vs_becomevu_min[NOTPOPART]>20){
        printf("Error:param->t_end_vs_becomevu_min[NOTPOPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vs_becomevu_min[POPART]<0 || param->t_end_vs_becomevu_min[POPART]>20){
        printf("Error:param->t_end_vs_becomevu_min[POPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vs_becomevu_range[NOTPOPART]<0 || param->t_end_vs_becomevu_range[NOTPOPART]>20){
        printf("Error:param->t_end_vs_becomevu_range[NOTPOPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vs_becomevu_range[POPART]<0 || param->t_end_vs_becomevu_range[POPART]>20){
        printf("Error:param->t_end_vs_becomevu_range[POPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vs_dropout_min[NOTPOPART]<0 || param->t_end_vs_dropout_min[NOTPOPART]>20){
        printf("Error:param->t_end_vs_dropout_min[NOTPOPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vs_dropout_min[POPART]<0 || param->t_end_vs_dropout_min[POPART]>20){
        printf("Error:param->t_end_vs_dropout_min[POPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vs_dropout_range[NOTPOPART]<0 || param->t_end_vs_dropout_range[NOTPOPART]>20){
        printf("Error:param->t_end_vs_dropout_range[NOTPOPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vs_dropout_range[POPART]<0 || param->t_end_vs_dropout_range[POPART]>100){
        printf("Error:param->t_end_vs_dropout_range[POPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vu_becomevs_min[NOTPOPART]<0 || param->t_end_vu_becomevs_min[NOTPOPART]>20){
        printf("Error:param->t_end_vu_becomevs_min[NOTPOPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vu_becomevs_min[POPART]<0 || param->t_end_vu_becomevs_min[POPART]>20){
        printf("Error:param->t_end_vu_becomevs_min[POPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vu_becomevs_range[NOTPOPART]<0 || param->t_end_vu_becomevs_range[NOTPOPART]>20){
        printf("Error:param->t_end_vu_becomevs_range[NOTPOPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vu_becomevs_range[POPART]<0 || param->t_end_vu_becomevs_range[POPART]>20){
        printf("Error:param->t_end_vu_becomevs_range[POPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vu_dropout_min[NOTPOPART]<0 || param->t_end_vu_dropout_min[NOTPOPART]>20){
        printf("Error:param->t_end_vu_dropout_min[NOTPOPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vu_dropout_min[POPART]<0 || param->t_end_vu_dropout_min[POPART]>20){
        printf("Error:param->t_end_vu_dropout_min[POPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vu_dropout_range[NOTPOPART]<0 || param->t_end_vu_dropout_range[NOTPOPART]>20){
        printf("Error:param->t_end_vu_dropout_range[NOTPOPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_end_vu_dropout_range[POPART]<0 || param->t_end_vu_dropout_range[POPART]>20){
        printf("Error:param->t_end_vu_dropout_range[POPART] is outside expected range [0,20]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_get_vmmc_min[NOTPOPART]<0 || param->t_get_vmmc_min[NOTPOPART]>100){
        printf("Error:param->t_get_vmmc_min[NOTPOPART]= %f is outside expected range [0,100]\nExiting\n",param->t_get_vmmc_min[NOTPOPART]);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_get_vmmc_min[POPART]<0 || param->t_get_vmmc_min[POPART]>100){
        printf("Error:param->t_get_vmmc_min[POPART]=%f is outside expected range [0,100] \nExiting\n",param->t_get_vmmc_min[POPART]);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_get_vmmc_range[NOTPOPART]<0 || param->t_get_vmmc_range[NOTPOPART]>5){
        printf("Error:param->t_get_vmmc_range[NOTPOPART]=%f is outside expected range [0,5]\nExiting\n",param->t_get_vmmc_range[NOTPOPART]);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_get_vmmc_range[POPART]<0 || param->t_get_vmmc_range[POPART]>5){
        printf("Error:param->t_get_vmmc_range[POPART] is outside expected range [0,5]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->t_vmmc_healing<0 || param->t_vmmc_healing>0.5){
        printf("Error:param->t_vmmc_healing is outside expected range [0,0.5]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->start_time_hiv<1950 || param->start_time_hiv>1990){
        printf("Error:param->start_time_hiv is outside expected range [1950,1990]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->start_time_simul<1900 || param->start_time_simul>1980){
        printf("Error:param->start_time_simul is outside expected range [1900,1980]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->start_time_simul>param->start_time_hiv){
        printf("Error: param->start_time_simul is bigger than param->start_time_hiv.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->end_time_simul<2000 || param->end_time_simul>2100){
        printf("Error:param->end_time_simul is outside expected range [2000,2100]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->COUNTRY_HIV_TEST_START<1970 || param->COUNTRY_HIV_TEST_START>2010){
        printf("Error:param->COUNTRY_HIV_TEST_START is outside expected range [1970,2010]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->COUNTRY_ART_START<1985 || param->COUNTRY_ART_START>2100){
        printf("Error:param->COUNTRY_ART_START is outside expected range [1985,2015]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->start_time_hiv>param->COUNTRY_HIV_TEST_START){
        printf("Error: param->start_time_hiv is bigger than param->COUNTRY_HIV_TEST_START.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->COUNTRY_HIV_TEST_START>param->COUNTRY_ART_START){
        printf("Error: param->COUNTRY_HIV_TEST_START is bigger than param->COUNTRY_ART_START.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->COUNTRY_CD4_350_START<param->COUNTRY_ART_START || param->COUNTRY_CD4_350_START>param->COUNTRY_CD4_500_START){
        printf("Error:param->COUNTRY_CD4_350_START is outside expected range [COUNTRY_ART_START,COUNTRY_CD4_500_START]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->COUNTRY_CD4_500_START>2100){
        printf("Error:param->COUNTRY_CD4_500_START is outside expected range between COUNTRY_CD4_350_START and 2020\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (param->COUNTRY_VMMC_START<2005 || param->COUNTRY_VMMC_START>2100){
        printf("Error:param->COUNTRY_VMMC_START is outside expected range [2005,2010]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->COUNTRY_HIV_TEST_START>param->COUNTRY_VMMC_START){
        printf("Error: param->COUNTRY_HIV_TEST_START is bigger than param->COUNTRY_VMMC_START.\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    for (chips_round=0; chips_round<NCHIPSROUNDS; chips_round++){
        if (param->CHIPS_START_YEAR[chips_round]<2013 || param->CHIPS_START_YEAR[chips_round]>2019){
            printf("Error:param->CHIPS_START_YEAR[chips_round]=%i is outside expected range [2013,2016]\nExiting\n",param->CHIPS_START_YEAR[chips_round]);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        //if (param->COUNTRY_ART_START>param->CHIPS_START_YEAR[chips_round]){
        //    printf("Error: param->COUNTRY_ART_START is bigger than param->CHIPS_START_YEAR[chips_round].\nExiting\n");
       //     printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
       //     fflush(stdout);
       //     exit(1);
       // }
        if (param->CHIPS_END_YEAR[chips_round]<2014 || param->CHIPS_END_YEAR[chips_round]>2020){
            printf("Error:param->CHIPS_END_YEAR[chips_round]=%i is outside expected range [2014,2020]\nExiting\n",param->CHIPS_END_YEAR[chips_round]);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        if (param->CHIPS_START_YEAR[chips_round]>param->CHIPS_END_YEAR[chips_round]){
            printf("Error: param->CHIPS_START_YEAR[chips_round] is bigger than param->CHIPS_END_YEAR[chips_round].\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        if (param->CHIPS_END_YEAR[chips_round]>param->end_time_simul){
            printf("Error: param->CHIPS_END_YEAR[chips_round] is bigger than param->end_time_simul.\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        if (param->CHIPS_START_TIMESTEP[chips_round]<0 || param->CHIPS_START_TIMESTEP[chips_round]>N_TIME_STEP_PER_YEAR){
            printf("Error:param->CHIPS_START_TIMESTEP[chips_round] is outside expected range [0,N_TIME_STEP_PER_YEAR]\nHas timestep changed?\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        if (param->CHIPS_END_TIMESTEP[chips_round]<0 || param->CHIPS_END_TIMESTEP[chips_round]>N_TIME_STEP_PER_YEAR){
            printf("Error:param->CHIPS_END_TIMESTEP[chips_round] is outside expected range [0,%i]\nHas timestep changed?\nExiting\n",N_TIME_STEP_PER_YEAR);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }


    if (param->chips_params->n_timesteps_per_round_posttrial<24 || param->chips_params->n_timesteps_per_round_posttrial>N_TIME_STEP_PER_YEAR){
        printf("Error:param->chips_params->n_timesteps_per_round_posttrial is outside expected range [24,%i]\n",N_TIME_STEP_PER_YEAR);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->DHS_params->NDHSROUNDS<2 || param->DHS_params->NDHSROUNDS>NDHSROUNDS_MAX){
        printf("Error:param->DHS_params->NDHSROUNDS is outside expected range [2,%i]\n",NDHSROUNDS_MAX);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    for (dhs_round=0; dhs_round<param->DHS_params->NDHSROUNDS; dhs_round++){
        if (param->DHS_params->DHS_YEAR[dhs_round]<2000 || param->DHS_params->DHS_YEAR[dhs_round]>2015){
            printf("Error:param->DHS_params->DHS_YEAR[dhs_round] is outside expected range [2000,2015]\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }



    if (param->initial_population_size<0 || param->initial_population_size>500000){
        printf("Error:param->initial_population_size is outside expected range [0,500000]\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if (param->initial_population_size>MAX_POP_SIZE){
        printf("Error:param->initial_population_size is bigger than MAX_POP_SIZE\nExiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    for(ag=0; ag<N_AGE; ag++){
        if (param->c_per_gender_within_patch[FEMALE][ag]<0 || param->c_per_gender_within_patch[FEMALE][ag]>20){
            printf("Error:param->c_per_gender_within_patch[FEMALE][ag] is outside expected range [0,20]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    for (ag=0; ag<N_AGE; ag++){
        if (param->c_per_gender_within_patch[MALE][ag]<0 || param->c_per_gender_within_patch[MALE][ag]>30){
            printf("Error:param->c_per_gender_within_patch[MALE][ag] is outside expected range [0,30]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    for (r=0; r<N_RISK; r++){
        if (param->relative_number_partnerships_per_risk[r]<0 || param->relative_number_partnerships_per_risk[r]>50){
            printf("Error:param->relative_number_partnerships_per_risk[r] is outside expected range [0,50]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    for(r=0 ; r<N_RISK ; r++){
        if (param->max_n_part_noage[r]<0 || param->max_n_part_noage[r]>20){
            printf("Error:param->max_n_part_noage[r] is outside expected range [0,20]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    for(r=0 ; r<N_RISK ; r++){
        if (param->breakup_scale_lambda_within_patch[r]<0 || param->breakup_scale_lambda_within_patch[r]>40){
            printf("Error:param->breakup_scale_lambda_within_patch[r=%i]=%6.4lf is outside expected range [0,40]\nExiting\n",r,param->breakup_scale_lambda_within_patch[r]);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    for(r = 0; r < N_RISK; r++){
        if(param->breakup_scale_lambda_between_patch[r] < 0 || param->breakup_scale_lambda_between_patch[r] > 30){
            printf("Error:param->breakup_scale_lambda_between_patch[r] is outside expected range [0,30]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
    for(r = 0; r < N_RISK; r++){
        if(param->breakup_shape_k[r] < 0 || param->breakup_shape_k[r] > 20){
            printf("Error:param->breakup_shape_k[r] is outside expected range [0,20]\nExiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
}
