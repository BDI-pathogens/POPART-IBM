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

/* Contains functions used only in debugging. The idea is it is easy to delete them afterwards.
 * */

#include "constants.h"
#include "output.h"
#include "debug.h"

/* Functions:
find_in_age_list()
print_age_list()
count_population_by_going_through_indiv()
count_population_by_going_through_age_list()
count_population_using_n_population()
count_population_size_three_ways()
is_in_patch_individual_population()
generate_demographics_byage_gender_file_header()
write_demographics_byage_gender()
blank_one_year_age_groups_including_kids()
write_one_year_age_groups_including_kids()
write_nbirths_nnewadults_ndeaths()
output_life_expectancy()
check_if_individual_should_be_in_list_susceptibles_in_serodiscordant_partnership()
check_if_individual_should_be_in_list_available_partners()
sweep_through_all_and_check_lists_serodiscordant_and_available_partners()
sweep_through_all_and_check_n_partners_outside_n_HIVpos_partners_and_n_HIVpos_partners_outside()
sweep_through_all_and_check_age_and_risk_of_partners()
blank_debugging_files()
write_hiv_duration()
 */


void find_in_age_list(double t,individual* person_to_find, age_list_struct *age_list, 
    parameters *param){

    int ai_age_calc; /* This is the calculated ai index. */
    /* Index for manual search through age_list looking for indvidiuals. */
    int ai_age_manual_search;
    int found_person = 0;
    int g = person_to_find->gender;
    int aa;
    aa = (int) floor(floor(t) - person_to_find->DoB) - AGE_ADULT;
    
    ai_age_calc = age_list->age_list_by_gender[g]->youngest_age_group_index + aa;
    
    while (ai_age_calc>(MAX_AGE-AGE_ADULT-1)){
        ai_age_calc = ai_age_calc - (MAX_AGE-AGE_ADULT);
    }
    
    int ai_age_calcv1 = get_age_index(person_to_find->DoB, param->start_time_simul);
    int ai_age_calcv2 = get_age_indexv2(person_to_find->DoB, t, 
        age_list->age_list_by_gender[g]->youngest_age_group_index);
    
    /* Looking for the person_to_find in the age_list --> 
    presumably only for debugging, could get rid of this in final code to speed up. 
    */
    
    int n;
    for(ai_age_manual_search = 0; 
        ai_age_manual_search < (MAX_AGE - AGE_ADULT); 
        ai_age_manual_search++
    ){
        n = 0;
        //printf("ai_age_manual_search=%i %li\n",ai_age_manual_search,age_list->number_per_age_group[ai_age_manual_search]);
        while(
            (n < age_list->age_list_by_gender[g]->number_per_age_group[ai_age_manual_search]) && 
            (age_list->age_list_by_gender[g]->age_group[ai_age_manual_search][n]->id != person_to_find->id)
        ){
            n++;
        }
        if((n < age_list->age_list_by_gender[g]->number_per_age_group[ai_age_manual_search]) && 
            (age_list->age_list_by_gender[g]->age_group[ai_age_manual_search][n]->id == person_to_find->id)
        ){
            printf("*********DEBUG DATA at time t=%f : ", t);
            printf("Individual %li Gender %i DoB = %f found with ai = %i. ", 
                    person_to_find->id, g, person_to_find->DoB, ai_age_manual_search);
            
            printf("Calc ai=%i, Calcv1 =%i Calcv2 =%i aa=%i, youngest_index=%i\n",
                ai_age_calc, ai_age_calcv1, ai_age_calcv2, aa, 
                age_list->age_list_by_gender[g]->youngest_age_group_index);
            
            found_person = 1;
            break;
        }
    }
    
    if(found_person == 0){
        n = 0;
        while((n < age_list->age_list_by_gender[g]->number_oldest_age_group) &&
            (age_list->age_list_by_gender[g]->oldest_age_group[n]->id != person_to_find->id)){
            n++;
        }
        if((n < age_list->age_list_by_gender[g]->number_oldest_age_group) && 
            (age_list->age_list_by_gender[g]->oldest_age_group[n]->id == person_to_find->id)){
            
            printf("*********DEBUG DATA: Individual %li DoB = %f gender %d at time t=%f found in oldest age go. Calc ai=%i, aa=%i, youngest_index=%i\n",person_to_find->id,person_to_find->DoB,g,t,ai_age_calc,aa,age_list->age_list_by_gender[g]->youngest_age_group_index);
            found_person = 1;
        }
    }
    if(found_person == 0){
        printf("**********PROBLEM2: person %li from patch %d gender %d not found when searching by hand in age_list %li\n",person_to_find->id,person_to_find->patch_no,g,age_list->age_list_by_gender[g]->number_per_age_group[2]);
        fflush(stdout);
        //printf("age_list->age_group[2][53]->id = %li %li\n",age_list->age_group[2][53]->id,person_to_find->id);

    }
    
    if(found_person == 1){
        printf("Found person %li from patch %d gender %d when searching by hand in age_list\n",
            person_to_find->id, g, person_to_find->patch_no);
    }
}


void print_age_list(age_list_struct *age_list){
    /* For each year, print the ID of each individual in that age.  
    
    This function is not currently used.  
    
    Parameters
    ----------
    age_list: pointer to an age_list_struct structure
    
    Returns
    -------
    Nothing; prints age lists
    
    */
    
    int ai;
    int n, g;
    
    // Loop through genders
    for(g = 0; g < N_GENDER; g++){
        printf("Gender = %i\n", g);
        
        // Loop through each year of age
        for(ai = 0; ai < (MAX_AGE - AGE_ADULT); ai++){
            n = 0;
            printf("ai=%i :", ai);
            while((n < age_list->age_list_by_gender[g]->number_per_age_group[ai])){
                printf("%li ", age_list->age_list_by_gender[g]->age_group[ai][n]->id);
                n++;
            }
        }
    }
}


void count_population_by_going_through_indiv(patch_struct *patch, long *n_m_indiv, long *n_f_indiv){
    /* Count the number of individuals in a patch, stratified by gender, and return.  
    
    The counts of number of males and females are returned in the objects n_m_indiv and n_f_indiv
    respectively.  This function loops through individuals in the `individual` array of the patch 
    object.  
    
    Parameters
    ----------
    patch : patch_struct structure
        Patch structure in the simulation
    n_m_indiv : pointer to a long
        Output long in which to store the number of males
    n_f_indiv: pointer to a long
        Output long in which to store the number of females
    
    Returns
    -------
    Nothing; but updates the input arguments n_m_indiv, n_f_indiv with counts of males, females.  
    */
    
    long i;
    *n_m_indiv = 0;
    *n_f_indiv = 0;
    
    // Loop through each individual in the patch in question
    for(i = 0; i < patch->id_counter; i++){
        // Check the individual is alive
        if(patch->individual_population[i].cd4 != DEAD){
            // Check the age of the individual and update the counters
            if(patch->individual_population[i].gender == MALE){
                *n_m_indiv = *n_m_indiv + 1;
            }else{
                *n_f_indiv = *n_f_indiv + 1;
            }
        }
    }
}


void count_population_by_going_through_age_list(patch_struct *patch, 
    long *n_m_agelist, long *n_f_agelist){
    /* Count the number of individuals in a patch, stratified by gender, and return.  
    
    In contrast to count_population_by_going_through_indiv(), this function counts the number of 
    individuals in each gender by traversing the arrays: 
        patch->age_list->age_list_by_gender[GENDER]->number_per_age_group[AGEYEAR]
        patch->age_list->age_list_by_gender[GENDER]->number_oldest_age_group.  
    
    Parameters
    ----------
    patch : patch_struct structure
        Patch structure in the simulation
    n_m_indiv : pointer to a long
        Output long in which to store the number of males
    n_f_indiv: pointer to a long
        Output long in which to store the number of females
    
    Returns
    -------
    Nothing; but updates the input arguments n_m_indiv, n_f_indiv with counts of males, females
    */
    
    int aa;
    *n_m_agelist = 0;
    *n_f_agelist = 0;
    for(aa = 0; aa < MAX_AGE - AGE_ADULT; aa++){
        *n_m_agelist = *n_m_agelist + 
            patch->age_list->age_list_by_gender[MALE]->number_per_age_group[aa];
        
        *n_f_agelist = *n_f_agelist + 
            patch->age_list->age_list_by_gender[FEMALE]->number_per_age_group[aa];
    }
    
    *n_m_agelist = *n_m_agelist + 
        patch->age_list->age_list_by_gender[MALE]->number_oldest_age_group;
    
    *n_f_agelist = *n_f_agelist + 
        patch->age_list->age_list_by_gender[FEMALE]->number_oldest_age_group;
}


void count_population_using_n_population(population_size *pop, long *n_m_pop, long *n_f_pop){
    int ag, r;

    *n_m_pop = 0;
    *n_f_pop = 0;


    for(ag=0 ; ag<N_AGE ; ag++){
        for(r=0 ; r<N_RISK ; r++){
            *n_m_pop = *n_m_pop + pop->pop_size_per_gender_age_risk[MALE][ag][r];
            *n_f_pop = *n_f_pop + pop->pop_size_per_gender_age_risk[FEMALE][ag][r];
        }
    }
}


void count_population_size_three_ways(patch_struct *patch, int p, double t){
    long n_m_indiv, n_f_indiv;
    count_population_by_going_through_indiv(&patch[p], &n_m_indiv, &n_f_indiv);

    long n_m_agelist, n_f_agelist;
    count_population_by_going_through_age_list(&patch[p], &n_m_agelist, &n_f_agelist);

    long n_m_pop, n_f_pop;
    count_population_using_n_population(patch[p].n_population, &n_m_pop, &n_f_pop);

    if ((n_m_indiv!=n_m_agelist) || (n_m_indiv!=n_m_pop) || (n_f_indiv!=n_f_agelist) || (n_f_indiv!=n_f_pop)){
        printf("Error - population sizes not matching in patch %i t=%6.4lf male: %ld %ld %ld female: %ld %ld %ld\n",p,t,n_m_indiv,n_m_agelist,n_m_pop,n_f_indiv,n_f_agelist,n_f_pop);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    //printf("Patch %i t=%6.4lf male: %ld %ld %ld female: %ld %ld %ld\n",p,t,n_m_indiv,n_m_agelist,n_m_pop,n_f_indiv,n_f_agelist,n_f_pop);
}


int is_in_patch_individual_population(long id, patch_struct *patch, int p)
{
    long k = 0;
    int res = 0;

    while((patch[p].individual_population[k].id!=id) && (k<(patch[p].id_counter-1)))
    {
        k++;
    }
    if(k<patch[p].id_counter-1)
    {
        res = 1;
        printf("Individual %ld was found in patch %d at position %ld in list\n",id,p,k);
        //print_individual(&patch[p].individual_population[k]);
        fflush(stdout);
    }else
    {
        printf("Individual %ld was not found in patch %d\n",id,p);
        fflush(stdout);
    }

    return(res);
}


/****************************************************************************************************************
 *****************************************************************************************************************
 * Generate output files for demographic validation:
 *****************************************************************************************************************
 ****************************************************************************************************************/

/* Function generates the text that goes in the header to the file Age_distribution_check_CL....csv. */
void generate_demographics_byage_gender_file_header(char *age_group_string, int size_age_group_string){
    int g,ag;
    char temp_string[100];

    sprintf(age_group_string,"Time,");
    for (g=0;g<N_GENDER;g++){
        for (ag=0; ag<N_AGE_UNPD; ag++){
            if(g==MALE)
                sprintf(temp_string,"M:%i-%i,",AGE_GROUPS_UNPD[ag],AGE_GROUPS_UNPD[ag+1]-1);
            else
                sprintf(temp_string,"F:%i-%i,",AGE_GROUPS_UNPD[ag],AGE_GROUPS_UNPD[ag+1]-1);
            join_strings_with_check(age_group_string, temp_string, size_age_group_string, "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
        }
        if(g==MALE)
            sprintf(temp_string,"M%i+,",AGE_GROUPS_UNPD[N_AGE_UNPD]);
        else
            sprintf(temp_string,"F:%i+,",AGE_GROUPS_UNPD[N_AGE_UNPD]);
        join_strings_with_check(age_group_string, temp_string, size_age_group_string, "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
    }

    for (g=0;g<N_GENDER;g++){
        for (ag=0; ag<N_AGE_UNPD; ag++){
            if(g==MALE)
                sprintf(temp_string,"DeadM:%i-%i,",AGE_GROUPS_UNPD[ag],AGE_GROUPS_UNPD[ag+1]-1);
            else
                sprintf(temp_string,"DeadF:%i-%i,",AGE_GROUPS_UNPD[ag],AGE_GROUPS_UNPD[ag+1]-1);
            join_strings_with_check(age_group_string, temp_string, size_age_group_string, "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
        }
        if(g==MALE)
            sprintf(temp_string,"DeadM%i+,",AGE_GROUPS_UNPD[N_AGE_UNPD]);
        else
            sprintf(temp_string,"DeadF:%i+,",AGE_GROUPS_UNPD[N_AGE_UNPD]);
        join_strings_with_check(age_group_string, temp_string, size_age_group_string, "temp_string and age_group_string in generate_demographics_byage_gender_file_header");
    }
    sprintf(temp_string,"CumulativeDead\n");
    join_strings_with_check(age_group_string, temp_string, size_age_group_string, "temp_string and age_group_string in generate_demographics_byage_gender_file_header");


}


/* Prints the number of people by age group and gender.
 * Function goes through individual_population in a given patch at a given time (currently called yearly in main.c).
 * Generates a csv file (for each run and patch) called  Age_distribution_check... .csv which gives the number of
 * adults in each 5 year UNPD age group by gender to compare against UNPD age distribution estimates for model validation.
 * Model also outputs number of dead people in each age group (their current age, not age ever).
 */

void write_demographics_byage_gender(patch_struct *patch, int p, double t, file_struct *file_data_store){
    long i;
    int ag;


    file_data_store->AGEDISTRIBUTIONFILE[p] = fopen(file_data_store->filename_debug_agedistribution[p],"a");
    if (file_data_store->AGEDISTRIBUTIONFILE[p]==NULL){
        printf("Cannot open output_file in write_demographics_byage_gender().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }


    /* Number of men/women in each age group (age groups correspond to UNPD 5 year age groups apart from 13-14 year olds. Note that we DO NOT want to count children as childhood mortality is taken to occur at birth so the number of children in the model is the number of children who will survive to age 13, not the number of children alive at any given time. */
    long  *num_age_people_m;
    long  *num_age_people_f;
    /* Number of deaths in the past year in men/women by age group : */
    long  *num_age_newdeaths_m;
    long  *num_age_newdeaths_f;
    long num_deaths_total = 0; /* This is a check that if we sum over num_age_newdeaths_m and num_age_newdeaths_f we get the new deaths per year. num_deaths_total is a cumulative measure so the difference between one year and the next is incident deaths. */

    /* Use calloc() so that these are initialized to zero. */
    num_age_people_m = (long*)calloc(N_AGE_UNPD+1, sizeof(long));
    num_age_people_f = (long*)calloc(N_AGE_UNPD+1, sizeof(long));
    num_age_newdeaths_m = (long*)calloc(N_AGE_UNPD+1, sizeof(long));
    num_age_newdeaths_f = (long*)calloc(N_AGE_UNPD+1, sizeof(long));

    for (i=0; i<patch[p].id_counter; i++){
        if (patch[p].individual_population[i].cd4>DEAD){
            //printf("i=%i dob=%6.4f t= %6.4f\n",get_age_group_unpd(individual_population[i].DoB,t),individual_population[i].DoB,t);
            if (patch[p].individual_population[i].gender==MALE)
                num_age_people_m[get_age_group_unpd(patch[p].individual_population[i].DoB,t)]++;
            else if (patch[p].individual_population[i].gender==FEMALE)
                num_age_people_f[get_age_group_unpd(patch[p].individual_population[i].DoB,t)]++;
            else{
                printf("ERROR: Unknown gender. Exiting\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
        else{
            num_deaths_total += 1;
            if (patch[p].individual_population[i].DoD<1800){
                printf("Error: someone has died before 1800. Exiting.\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            /* write_demographics_byage_gender() is called at the end of the year (so t is the start of the year), this counts anyone who died in the last year. */
            else if(patch[p].individual_population[i].DoD>=t){
                if (patch[p].individual_population[i].gender==MALE)
                    num_age_newdeaths_m[get_age_group_unpd(patch[p].individual_population[i].DoB,t)]++;
                else if (patch[p].individual_population[i].gender==FEMALE)
                    num_age_newdeaths_f[get_age_group_unpd(patch[p].individual_population[i].DoB,t)]++;
                else{
                    printf("ERROR: Unknown gender in dead person. Exiting\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
            }
        }
    }

    fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%6.4f,",t+1);
    for (ag=0; ag<=N_AGE_UNPD; ag++)
        fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%li,",num_age_people_m[ag]);
    for (ag=0; ag<=N_AGE_UNPD; ag++)
        fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%li,",num_age_people_f[ag]);

    for (ag=0; ag<=N_AGE_UNPD; ag++)
        fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%li,",num_age_newdeaths_m[ag]);
    for (ag=0; ag<=N_AGE_UNPD; ag++)
        fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%li,",num_age_newdeaths_f[ag]);


    fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%li\n",num_deaths_total);
    fclose(file_data_store->AGEDISTRIBUTIONFILE[p]);

    free(num_age_people_m);
    free(num_age_people_f);
    free(num_age_newdeaths_m);
    free(num_age_newdeaths_f);
}


void blank_one_year_age_groups_including_kids(file_struct *file_data_store){
    /* Write blank file and header for the file called OneYearAgeGp_*.csv
    
    For now we only use data from patch 0.  
    
    Arguments
    ---------
    file_data_store : pointer to file_struct structure
        This structure stores the file names of the files to write.  The file in question has its
        file name stored in the attribute called file_data_store->ONEYEARAGEDISTRIBUTIONFILE
    */
    
    // Open a connection to the file
    file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0] =
        fopen(file_data_store->filename_debug_one_yearage_dist_includingkids[0], "w");
    
    // Throw an error if we can't open the file.  
    if(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0] == NULL){
        printf("Cannot open file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0]\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    // Print header for year (time)
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0], "t,");
    
    int i;
    
    // Write headers for children (they don't have a gender in the model)
    for(i = 0; i < AGE_ADULT + 1; i++){
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0],"Nage%i-%i,", i, i + 1);
    }
    
    // Write headers for males
    for(i = AGE_ADULT + 1; i < MAX_AGE; i++){
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0],"NMage%i-%i,", i, i + 1);
    }
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0], "NMage80+,");
    
    // Write headers for females
    for(i = AGE_ADULT + 1; i < MAX_AGE; i++){
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0],"NFage%i-%i,", i, i + 1);
    }
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0], "NFage80+\n");
    
    fclose(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[0]);
    return;
}


void write_one_year_age_groups_including_kids(file_struct *file_data_store, patch_struct *patch,
    int p, double t){
    /* Write population size per yearly age groups to file for a single year (append to file)
    
    This function writes the totaol population size for children from ages 0-1 until 13-14 at
    yearly age gaps (children aren't split by gender in the simulation), population sizes are then 
    recorded from 14-15 until 80+ at yearly age groups split for males and females.  These
    values are taken from the structures `child_population` and `age_list->age_list_by_gender` 
    within the patch structure respectively.  The file in which these values are appended starts 
    with the prefix 'OneYearAgeGp_*.csv'.  This function is called from main.c and is only triggered
    if the macro (defined in constants.h) WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_ONEYEARINCKIDS is 1.  
    The function blank_one_year_age_groups_including_kids() sets up the file and header line that 
    this function writes to.  
    
    The `child_population` array that's part of the patch structure is indexed by HIV status of the
    children, 0 index for HIV- and 1 index for HIV+.  However, at this stage (Feb 2018) all children
    are HIV- so the second entry in the array is all zeros.  
    
    Arguments
    ---------
    file_data_store : pointer to a file_struct structure
    patch : pointer to a patch_struct structure
    p : int
    t : double
    
    Returns
    -------
    Nothing; writes a file to disk
    */
    
    
    int i_child, ai_m, ai_f, aa, hivstatus;
    
    // Open a connection to the file in question and append to it ("a")
    file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p] =
        fopen(file_data_store->filename_debug_one_yearage_dist_includingkids[p], "a");
    
    // Throw an error if the file can't be opened
    if(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p] == NULL){
        printf("Cannot open file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p]\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    // Append the year
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%6.4lf,", t);

    // The 2 is because we have HIV- and HIV+ children
    int child_start_index[2];
    
    // If the oldest child group is at the end, then the youngest age group is [0]
    for(hivstatus = 0; hivstatus < 2; hivstatus++){
        
        if(
        patch[p].child_population[hivstatus].debug_tai == (AGE_ADULT + 1) * N_TIME_STEP_PER_YEAR - 1
        ){
            child_start_index[hivstatus] = 0;
        }else{
            // Otherwise the youngest age group is the one after the oldest age group
            child_start_index[hivstatus] = patch[p].child_population[hivstatus].debug_tai + 1;
        }
    }

    // Write to disk the number of kids by 1 year age groups
    long nkids_by_year = 0;
    
    for(i_child = 0; i_child < ((AGE_ADULT + 1) * N_TIME_STEP_PER_YEAR); i_child++){
        
        // Sum up number of kids by yearly age group
        for(hivstatus = 0; hivstatus < 2; hivstatus++){
            nkids_by_year +=
                patch[p].child_population[hivstatus].n_child[child_start_index[hivstatus]];
        }
        
        /* Print out number of kids and then reset counter for next yearly age group: */
        if((i_child % N_TIME_STEP_PER_YEAR) == (N_TIME_STEP_PER_YEAR - 1)){
            fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li,", nkids_by_year);
            nkids_by_year = 0;
        }
        
        /* Now update pointer to loop over all kids: */
        for(hivstatus = 0; hivstatus < 2; hivstatus++){
            if(child_start_index[hivstatus] < (AGE_ADULT + 1) * N_TIME_STEP_PER_YEAR - 1){
                child_start_index[hivstatus]++;
            }else{
                child_start_index[hivstatus] = 0;
            }
        }
    }
    
    /* Write number of adults to disk.  
    Note we don't print the age 13 group here as we call this routine at the start of the year -
    so no adults aged 13 yet.  We count the age 13-14 group from the kids above.  Hence the range
    of aa is 0 to (MAX_AGE-AGE_ADULT-1) instead of (MAX_AGE-AGE_ADULT), and we add 1 to ai_m and
    ai_f
    */
    // Output count of individuals per age for males
    for(aa = 0; aa < (MAX_AGE - AGE_ADULT - 1); aa++){
        ai_m = aa + patch[p].age_list->age_list_by_gender[MALE]->youngest_age_group_index + 1;
        
        while(ai_m > (MAX_AGE - AGE_ADULT - 1)){
            ai_m = ai_m - (MAX_AGE - AGE_ADULT);
        }
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li,", 
            patch[p].age_list->age_list_by_gender[MALE]->number_per_age_group[ai_m]);
    }
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li,", 
        patch[p].age_list->age_list_by_gender[MALE]->number_oldest_age_group);
    
    // Output count of individuals per age for females
    for(aa = 0; aa < (MAX_AGE - AGE_ADULT - 1); aa++){
        ai_f = aa + patch[p].age_list->age_list_by_gender[FEMALE]->youngest_age_group_index + 1;
        
        while(ai_f > (MAX_AGE - AGE_ADULT - 1)){
            ai_f = ai_f - (MAX_AGE - AGE_ADULT);
        }
        fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li,", 
            patch[p].age_list->age_list_by_gender[FEMALE]->number_per_age_group[ai_f]);
    }
    fprintf(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p], "%li\n", 
        patch[p].age_list->age_list_by_gender[FEMALE]->number_oldest_age_group);
    
    // Close the connection to the file
    fclose(file_data_store->ONEYEARAGEDISTRIBUTIONFILE[p]);
}


void write_nbirths_nnewadults_ndeaths(file_struct *file_data_store, patch_struct *patch, int year){

    file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE = fopen(file_data_store->filename_debug_nnewadults_ndeaths_file ,"a");
    if (file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE==NULL){
        printf("Cannot open file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"%d,",year);
    int p;
    for(p=0;p<NPATCHES;p++)
        fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"%ld,%ld,%ld,",patch[p].DEBUG_NBIRTHS,patch[p].DEBUG_NNEWADULTS,patch[p].DEBUG_NDEATHS);
    fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"\n");
    fclose(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE);
}







/* Code not currently used. Would generate the duration of life of people in 5 year age cohorts. */
void output_life_expectancy(char *output_file_directory, patch_struct *patch, int p, int i_run){
    long i_id;
    int i_a;   /* Indexes which birth cohort you are in. For now make 5 year age cohorts (ie born 1900-1905, 1905-1910,...). */
    double age_at_death;
    /* If simulation run for too short a time this won't work: */
    if(patch[p].param->end_time_simul<2100){
        printf("Warning: need to make end_time_simul at least 2100 to run output_life_expectancy(). This function has not been run even though WRITE_DEBUG_DEMOGRAPHICS_LIFE_EXPECTANCY==1\n");
        return;
    }

    FILE *life_expectancy_file;
    char life_expectancy_filename[LONGSTRINGLENGTH];
    char templabel[LONGSTRINGLENGTH];     /* Length is arbitrary - but make sure we never have stupidly long filenames! */

    /* Make sure the arrays are blank: */
    memset(life_expectancy_filename, '\0', sizeof(life_expectancy_filename));
    memset(templabel, '\0', sizeof(templabel));

    /* Assume that country setting is same in all patches so use patch 0. */
    if (patch[0].country_setting==ZAMBIA)
        sprintf(templabel,"_Za.csv");
    else
        sprintf(templabel,"_SA.csv");



    strncpy(life_expectancy_filename,output_file_directory,LONGSTRINGLENGTH);
    add_slash(life_expectancy_filename); /* Adds a / or \ as needed if working in directory other than current local dir. */
    join_strings_with_check(life_expectancy_filename, "LifeExpectancy", LONGSTRINGLENGTH, "'LifeExpectancy' and life_expectancy_filename in output_life_expectancy()");
    join_strings_with_check(life_expectancy_filename, templabel, LONGSTRINGLENGTH, "templabel and life_expectancy_filename in output_life_expectancy()");
    printf("LE filename = %s\n",life_expectancy_filename);

    /* For simplicity we store all the data from every run in a single file. This means that we need to make the first run
     * (i_run=1) blank.
     */
    if (i_run==1){
        life_expectancy_file = fopen(life_expectancy_filename,"w");
        if (life_expectancy_file==NULL){
            printf("Cannot open life_expectancy_file\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        /* Print a header for the file: */
        fprintf(life_expectancy_file,"i_run,");
        for (i_a=0;i_a<20;i_a++){
            fprintf(life_expectancy_file,"%i-%i,",1900+i_a*5,1904+i_a*5);
        }
        fprintf(life_expectancy_file,"\n");
    }
    else{ /* if not the first run then just append to the existing file: */
        life_expectancy_file = fopen(life_expectancy_filename,"a");
        if (life_expectancy_file==NULL){
            printf("Cannot open life_expectancy_file\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }


    long n_by_age_cohort[20];            /* Store number of individual in each age cohort from 1900-2000 - so 20 age groups. */
    double cumulative_life_years[20];    /* Store number of life-years lived of people in cohort. */
    double adjusted_life_expectancy[20]; /* Life expectancy of each cohort adjusting for childhood mortality. */
    for (i_a=0;i_a<20;i_a++){
        n_by_age_cohort[i_a] = 0;
        cumulative_life_years[i_a] = 0;
    }

    /* Use this to estimate remaining life expectancy of someone still alive in 2000: */
    double annual_mortality_over80_in2000 = 0.5*(exp(patch[p].param->mortality_rate_by_gender_age_intercept[MALE][N_AGE_UNPD_MORTALITY-1] + patch[p].param->mortality_rate_by_gender_age_slope[MALE][N_AGE_UNPD_MORTALITY-1]*2000) + exp(patch[p].param->mortality_rate_by_gender_age_intercept[FEMALE][N_AGE_UNPD_MORTALITY-1] + patch[p].param->mortality_rate_by_gender_age_slope[FEMALE][N_AGE_UNPD_MORTALITY-1]*2000));
    //printf("annual_mortality_over80_in2000 = %lf\n",annual_mortality_over80_in2000); /* Checked in Zambia this was 0.16, so 6 years expected future life. */

    for (i_id=0;i_id<patch[p].id_counter;i_id++){
        /* Ignore anyone born before 1900 - ie the adults at initialisation - as we have incomplete mortality data for them.
         * Similarly assume that we may not have a good sample of those born after ~2000 as some may not die by the end of the simulation. */
        if ((patch[p].individual_population[i_id].DoB>1900) && (patch[p].individual_population[i_id].DoB<2000)){
            /* Is the person dead by the end of the simulation? */
            if (patch[p].individual_population[i_id].DoD>-1){
                age_at_death = patch[p].individual_population[i_id].DoD - patch[p].individual_population[i_id].DoB;
                i_a = (int) floor((patch[p].individual_population[i_id].DoB-1900)/5);
                n_by_age_cohort[i_a]++;                       /* Add one to counter in the age cohort i_a. */
                cumulative_life_years[i_a] += age_at_death;   /* Add the life years of that individual to age cohort i_a. */
            }
            else{
                i_a = (int) floor((patch[p].individual_population[i_id].DoB-1900)/5);
                //printf("Individual %li in age gp %i not dead yet DoB=%6.4f\n",patch[p].individual_population[i_id].id,i_a,patch[p].individual_population[i_id].DoB);
                n_by_age_cohort[i_a]++;                       /* Add one to counter in the age cohort i_a. */
                cumulative_life_years[i_a] += 100 + 1.0/annual_mortality_over80_in2000;   /* Add remaining expected number of life years of that individual to age cohort i_a. */
            }
        }
    }

    fprintf(life_expectancy_file,"Run%iUnadjusted,",i_run);
    for (i_a=0;i_a<20;i_a++){
        fprintf(life_expectancy_file,"%8.6lf,",cumulative_life_years[i_a]/(n_by_age_cohort[i_a]*1.0));
    }
    fprintf(life_expectancy_file,"\n");

    /* Now adjust estimate for childhood mortality: */
    fprintf(life_expectancy_file,"Run%iAdjusted,",i_run);
    double annual_mortality_rate_under5, mortality_rate_under5;
    double annual_mortality_rate_5to10,  mortality_rate_5to10;
    int g;
    double t;
    for (i_a=0;i_a<20;i_a++){
        t = 1902.5+i_a*5; /* Take mid-point of age cohort birth date and get corresponding mortality: */
        annual_mortality_rate_under5 = 0.0;
        annual_mortality_rate_5to10 = 0.0;
        for (g=0;g<N_GENDER;g++){

            /* We by default assume that mortality in under 5 is mostly perinatal mortality (so occurs at time t) and that mortality in 5-10 year olds occurs uniformly over that age so on average at time t+7.5.
             * However we need to adjust these times for the fact that we only have data from 1950-2100. */
            if (t<1950){
                /* Average mortality over genders. The [0] and [1] indices refer to age groups 0-4 and 5-9: */
                annual_mortality_rate_under5 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][0] + patch[p].param->mortality_rate_by_gender_age_slope[g][0]*1950);
                annual_mortality_rate_5to10 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][1] + patch[p].param->mortality_rate_by_gender_age_slope[g][1]*1957.5);
            }
            else if (t>=2100){ /* Note - should never need t>2000. This is copied from a function in demographics, and keep for completeness. */
                annual_mortality_rate_under5 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][0] + patch[p].param->mortality_rate_by_gender_age_slope[g][0]*2100);
                annual_mortality_rate_5to10 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][1] + patch[p].param->mortality_rate_by_gender_age_slope[g][1]*2100);
            }
            /* for times 2092.5-2100 */
            else if (t>=2092.5){
                annual_mortality_rate_under5 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][0] + patch[p].param->mortality_rate_by_gender_age_slope[g][0]*t);
                annual_mortality_rate_5to10 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][1] + patch[p].param->mortality_rate_by_gender_age_slope[g][1]*2100);
            }
            /* for times 1950-2092.5: */
            else{
                annual_mortality_rate_under5 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][0] + patch[p].param->mortality_rate_by_gender_age_slope[g][0]*t);
                annual_mortality_rate_5to10 += exp(patch[p].param->mortality_rate_by_gender_age_intercept[g][1] + patch[p].param->mortality_rate_by_gender_age_slope[g][1]*(t+7.5));
            }
        }
        /* We take the average of the annual_mortality_rates over gender: */
        mortality_rate_under5 = 1- pow(1-annual_mortality_rate_under5/(N_GENDER*1.0),5);
        mortality_rate_5to10 = 1- pow(1-annual_mortality_rate_5to10/(N_GENDER*1.0),5);
        adjusted_life_expectancy[i_a] = cumulative_life_years[i_a]/(n_by_age_cohort[i_a]*1.0)*(1-mortality_rate_under5)*(1-mortality_rate_5to10) + 7.5*(1-mortality_rate_under5)*mortality_rate_5to10 + 0*mortality_rate_under5;
        fprintf(life_expectancy_file,"%8.6lf,",adjusted_life_expectancy[i_a]);
    }
    fprintf(life_expectancy_file,"\n");

    fclose(life_expectancy_file);

}

/***************************************************************************************************/
/* End of demographic checks. */
/***************************************************************************************************/





/***************************************************************************************************/
/* Start of partnerships checks: */
/***************************************************************************************************/


// CHECK 1: check if individual is alive, HIV- and has HIV+ partners, and if so whether he is in the right place in the list of susceptibles in a serodiscordant partnership
void check_if_individual_should_be_in_list_susceptibles_in_serodiscordant_partnership(individual *temp_ind, all_partnerships * overall_partnerships){
    int i;
    int isInSerodiscordantCouple = 0;

    if(temp_ind->cd4 > DEAD && temp_ind->HIV_status==0)
    {
        if(temp_ind->n_partners>0)
        {
            for(i=0 ; i<temp_ind->n_partners ; i++) // loop over all partners in case there is something wrong with HIV+ partners tracking
            {
                // is partner HIV+?
                if(temp_ind->partner_pairs[i]->ptr[1-temp_ind->gender]->HIV_status>0)
                {
                    if(temp_ind->id==FOLLOW_INDIVIDUAL  && temp_ind->patch_no==FOLLOW_PATCH)
                    {
                        printf("seropositive partner: ");
                        print_individual(temp_ind->partner_pairs[i]->ptr[1-temp_ind->gender]);
                        fflush(stdout);
                    }
                    isInSerodiscordantCouple = 1;
                }
            }
        }
        if(isInSerodiscordantCouple>0)
        {
            // check this person is in the right place in the list of susceptibles in a serodiscordant partnership
            if(temp_ind->idx_serodiscordant<0 || temp_ind->idx_serodiscordant>=overall_partnerships->n_susceptible_in_serodiscordant_partnership[0])
            {
                printf("PROBLEM: individual %ld from patch %d is not at all in the list of susceptible individuals in a serodiscordant couple\n",temp_ind->id,temp_ind->patch_no);
                print_individual(temp_ind);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }else if(overall_partnerships->susceptible_in_serodiscordant_partnership[temp_ind->idx_serodiscordant]->id != temp_ind->id ||  overall_partnerships->susceptible_in_serodiscordant_partnership[temp_ind->idx_serodiscordant]->patch_no != temp_ind->patch_no)
            {
                printf("PROBLEM: individual %ld from patch %d is not found where should be in the list of susceptible individuals in a serodiscordant couple\n",temp_ind->id,temp_ind->patch_no);
                print_individual(temp_ind);
                printf("BUT the person pointed to in the list is at idx %li:\n",temp_ind->idx_serodiscordant);
                print_individual(overall_partnerships->susceptible_in_serodiscordant_partnership[temp_ind->idx_serodiscordant]);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }
}

// CHECK 2: check if individual is alive and has available partnerships, whether he is in the right place in the list of available partners
void check_if_individual_should_be_in_list_available_partners(individual *temp_ind, all_partnerships * overall_partnerships, int t0, int t_step){
    int i;
    long temp_id;
    int temp_patch;
    int ag;
    int n_times_found_in_list_available_partners;

    if(temp_ind->cd4 > DEAD && temp_ind->n_partners<temp_ind->max_n_partners)
    {
        n_times_found_in_list_available_partners = 0;
        // find ag the age group of temp_ind
        ag = get_age_group(temp_ind->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);

        for(i=0 ; i<overall_partnerships->n_pop_available_partners->pop_per_patch[temp_ind->patch_no].pop_size_per_gender_age_risk[temp_ind->gender][ag][temp_ind->sex_risk] ; i++)
        {
            temp_id = overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[temp_ind->patch_no][temp_ind->gender][ag][temp_ind->sex_risk][i]->id;
            temp_patch = overall_partnerships->pop_available_partners->pop_per_patch_gender_age_risk[temp_ind->patch_no][temp_ind->gender][ag][temp_ind->sex_risk][i]->patch_no;

            if(temp_id == temp_ind->id && temp_patch == temp_ind->patch_no)
            {
                n_times_found_in_list_available_partners++;
            }
        }
        if(n_times_found_in_list_available_partners != temp_ind->max_n_partners - temp_ind->n_partners)
        {
            printf("PROBLEM: individual %ld from patch %d is not found as many time as expected in the list of available partners\n",temp_ind->id,temp_ind->patch_no);
            print_individual(temp_ind);
            printf("Individual was found %d times in the list but has %d available partnerships\n",n_times_found_in_list_available_partners,temp_ind->max_n_partners - temp_ind->n_partners);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
}

void sweep_through_all_and_check_lists_serodiscordant_and_available_partners (patch_struct *patch, all_partnerships * overall_partnerships, int t0, int t_step)
{
    int p, k;
    individual temp_ind;

    for(p=0 ; p<NPATCHES; p++)
    {
        for(k=0 ; k<patch[p].id_counter; k++)
        {
            temp_ind = patch[p].individual_population[k];

            if(temp_ind.id==FOLLOW_INDIVIDUAL  && temp_ind.patch_no==FOLLOW_PATCH)
            {
                print_individual(&temp_ind);
            }

            // CHECK 1: check if individual is alive, HIV- and has HIV+ partners, and if so whether he is in the right place in the list of susceptibles in a serodiscordant partnership
            check_if_individual_should_be_in_list_susceptibles_in_serodiscordant_partnership(&temp_ind, overall_partnerships);

            // CHECK 2: check if individual is alive and has available partnerships, whether he is in the right place in the list of available partners
            check_if_individual_should_be_in_list_available_partners(&temp_ind, overall_partnerships, t0, t_step);
        }

    }
}


void sweep_through_all_and_check_n_partners_outside_n_HIVpos_partners_and_n_HIVpos_partners_outside (patch_struct *patch, all_partnerships * overall_partnerships, int t0, int t_step)
{

    int p, k, i;
    individual temp_ind;
    int temp_patch;
    int check_n_partners_outside, check_n_HIVpos_partners, check_n_HIVpos_partners_outside;

    for(p=0 ; p<NPATCHES; p++)
    {
        for(k=0 ; k<patch[p].id_counter; k++)
        {
            temp_ind = patch[p].individual_population[k];
            temp_patch = temp_ind.patch_no;

            if(temp_ind.cd4 > DEAD) // only do if person still alive
            {

                if(temp_ind.id==FOLLOW_INDIVIDUAL  && temp_patch==FOLLOW_PATCH)
                {
                    print_individual(&temp_ind);
                }

                // is n_partners_outside what it should be?
                check_n_HIVpos_partners = 0;
                check_n_partners_outside = 0;
                check_n_HIVpos_partners_outside = 0;
                if(temp_ind.n_partners>0)
                {
                    for(i=0 ; i<temp_ind.n_partners ; i++)
                    {
                        if(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]->HIV_status>0)
                        {
                            check_n_HIVpos_partners++;
                        }
                        if(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]->patch_no != temp_patch)
                        {
                            check_n_partners_outside++;
                            if(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]->HIV_status>0)
                            {
                                check_n_HIVpos_partners_outside++;
                            }
                        }
                    }
                }

                if(check_n_partners_outside != temp_ind.n_partners_outside)
                {
                    printf("PROBLEM: individual %ld from patch %d does not have the correct number of partners outside the patch\n",temp_ind.id,temp_patch);
                    print_individual(&temp_ind);
                    if(temp_ind.n_partners>0)
                    {
                        printf("--- PARTNERS OF THE PROBLEMATIC INDIVIDUAL ---\n");
                        fflush(stdout);
                        for(i=0 ; i<temp_ind.n_partners ; i++)
                        {
                            print_individual(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]);
                        }
                        printf("------\n");
                        fflush(stdout);
                    }
                    printf("Individual has %d partners outside his/her patch but n_partners_outside is %d \n",check_n_partners_outside,temp_ind.n_partners_outside);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                if(temp_ind.HIV_status==0) // we only keep track of HIV positive partners for the HIV negative individuals
                {
                    if(check_n_HIVpos_partners != temp_ind.n_HIVpos_partners)
                    {
                        printf("PROBLEM: individual %ld from patch %d does not have the correct number of HIV+ partners\n",temp_ind.id,temp_patch);
                        print_individual(&temp_ind);
                        if(temp_ind.n_partners>0)
                        {
                            printf("--- PARTNERS OF THE PROBLEMATIC INDIVIDUAL ---\n");
                            fflush(stdout);
                            for(i=0 ; i<temp_ind.n_partners ; i++)
                            {
                                print_individual(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]);
                            }
                            printf("------\n");
                            fflush(stdout);
                        }
                        printf("Individual has %d HIV+ partners but n_HIVpos_partners is %d \n",check_n_HIVpos_partners,temp_ind.n_HIVpos_partners);
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }

                    if(check_n_HIVpos_partners_outside != temp_ind.n_HIVpos_partners_outside)
                    {
                        printf("PROBLEM: individual %ld from patch %d does not have the correct number of HIV+ partners outside the patch\n",temp_ind.id,temp_patch);
                        print_individual(&temp_ind);
                        if(temp_ind.n_partners>0)
                        {
                            printf("--- PARTNERS OF THE PROBLEMATIC INDIVIDUAL ---\n");
                            fflush(stdout);
                            for(i=0 ; i<temp_ind.n_partners ; i++)
                            {
                                print_individual(temp_ind.partner_pairs[i]->ptr[1-temp_ind.gender]);
                            }
                            printf("------\n");
                            fflush(stdout);
                        }
                        printf("Individual has %d HIV+ partners outside his/her patch but n_HIVpos_partners_outside is %d \n",check_n_HIVpos_partners,temp_ind.n_HIVpos_partners_outside);
                        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                        fflush(stdout);
                        exit(1);
                    }
                }

            }

        }

    }
}

void sweep_through_all_and_check_age_and_risk_of_partners (patch_struct *patch, all_partnerships * overall_partnerships, int t0, int t_step, debug_struct *debug)
{

    int p, k, i;
    individual *ind1, *ind2;
    int age_f, age_m, risk_f, risk_m;

    for(age_f=0 ; age_f<N_AGE ; age_f++)
    {
        for(age_m=0 ; age_m<N_AGE ; age_m++)
        {
            debug->age_of_partners_cross_sectional[t0 - patch[0].param->start_time_simul][age_f][age_m] = 0.0;
        }
    }

    for(risk_f=0 ; risk_f<N_RISK ; risk_f++)
    {
        for(risk_m=0 ; risk_m<N_RISK ; risk_m++)
        {
            debug->risk_of_partners_cross_sectional[t0 - patch[0].param->start_time_simul][risk_f][risk_m] = 0.0;
        }
    }

    for(p=0 ; p<NPATCHES; p++)
    {

        for(k=0 ; k<patch[p].id_counter; k++)
        {
            ind1 = &patch[p].individual_population[k];

            if(ind1->id == FOLLOW_INDIVIDUAL && ind1->patch_no == FOLLOW_PATCH)
            {
                print_individual(ind1);
                fflush(stdout);
            }

            if(ind1->cd4 > DEAD && ind1->n_partners>0)
            {
                for(i=0 ; i<ind1->n_partners ; i++)
                {
                    ind2 = ind1->partner_pairs[i]->ptr[1-ind1->gender];

                    if(ind1->gender == FEMALE){
                        age_f = get_age_group(ind1->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);
                        age_m = get_age_group(ind2->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);

                        risk_f = ind1->sex_risk;
                        risk_m = ind2->sex_risk;
                    }else
                    {
                        age_f = get_age_group(ind2->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);
                        age_m = get_age_group(ind1->DoB,t0+t_step*TIME_STEP, AGE_GROUPS, N_AGE);

                        risk_f = ind2->sex_risk;
                        risk_m = ind1->sex_risk;
                    }

                    /*if ((ind1->id == FOLLOW_INDIVIDUAL && ind1->patch_no == FOLLOW_PATCH) || (ind2->id == FOLLOW_INDIVIDUAL && ind2->patch_no == FOLLOW_PATCH))
                    {
                        printf("Age groups of partners %d %d\n",age_f,age_m);
                        printf("Risk groups of partners %d %d\n",risk_f,risk_m);
                        fflush(stdout);
                    }*/

                    debug->age_of_partners_cross_sectional[t0 - patch[0].param->start_time_simul][age_f][age_m] += 0.5; // will count each partnership twice otherwise
                    debug->risk_of_partners_cross_sectional[t0 - patch[0].param->start_time_simul][risk_f][risk_m] += 0.5; // will count each partnership twice otherwise

                    /*printf("Partnerhip between %ld in patch %d and %ld in patch %d\n",ind1->id,ind1->patch_no,ind2->id,ind2->patch_no);
                    fflush(stdout);*/
                }
            }

        }

    }

}

/***************************************************************************************************/
/* End of partnerships checks: */
/***************************************************************************************************/


/***************************************************************************************************/
/* Start of HIV checks: */
/***************************************************************************************************/

/* From http://www.statssa.gov.za/?p=2973 (accessed 24 June 2016):
"....the number of AIDS related deaths is estimated to have decreased from 363 910 deaths in 2005 (51% of all deaths) to 171 733 deaths in 2014 (31% of all deaths). This can be associated with the increased rollout of antiretroviral therapy (ART)."
 */

/* Checks:
 *  - Number of HIV deaths
 *  - Duration of HIV without ART by SPVL
 *  - Duration of HIV with ART
 *  - Estimation of R0 and/or doubling time
 */




/* This blanks the debugging files at the beginning of the run. The reason we need to do this is because debugging files
 * are *inefficient* - we write to them continuously using the "a" (ie append) write option so that we don't have to
 * keep a lot of potentially large arrays that we do not need for non-debugging runs.
 * Essentially i/o (reading/writing from/to disk) is inefficient so we want to avoid it as much as possible, but for
 * debugging files we make an exception.
 * For files which save output and only write occasionally (e.g. Annual_output*.csv files) we don't need to blank the files
 * at the beginning of a run.

 */
void blank_debugging_files(file_struct *file_data_store){
    char age_group_string[1000];
    int p;

    for (p=0;p<NPATCHES;p++){

        if (WRITE_DEBUG_HIV_DURATION==1){
            file_data_store->HIVDURATIONFILE[p] = fopen(file_data_store->filename_debughivduration[p],"w");
            fprintf(file_data_store->HIVDURATIONFILE[p],"time,patch,id,Time_HIV+,Time_sc,ART_Status,SPVL_cat,CD4\n");
            fclose(file_data_store->HIVDURATIONFILE[p]);
        }
        if (WRITE_DEBUG_HIV_DURATION_KM==1){
            file_data_store->HIVDURATIONFILE_KM[p] = fopen(file_data_store->filename_debughivduration_km[p],"w");
            fprintf(file_data_store->HIVDURATIONFILE_KM[p],"time,time_pos,reason,gender,CD4,SPVLnum,SPVLcat\n");
            fclose(file_data_store->HIVDURATIONFILE_KM[p]);
        }

        if (WRITE_DEBUG_CD4_AFTER_SEROCONVERSION==1){
            file_data_store->HIVCD4_AFTER_SEROCONVERSION[p] = fopen(file_data_store->filename_debughivcd4_after_seroconversion[p],"w");
            fprintf(file_data_store->HIVCD4_AFTER_SEROCONVERSION[p],"id,CD4,SPVLcat\n");
            fclose(file_data_store->HIVCD4_AFTER_SEROCONVERSION[p]);
        }

        if (WRITE_DEBUG_INITIAL_SPVL_DISTRIBUTION==1){
            file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p] = fopen(file_data_store->filename_debuginitial_spvl_distribution[p],"w");
            fprintf(file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p],"SPVL_cat,SPVL_E,SPVL_G\n");
            fclose(file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p]);
        }

        if (WRITE_DEBUG_HIV_STATES==1){
            file_data_store->HIVSTATEPOPULATIONFILE[p] = fopen(file_data_store->filename_debug_hivpopulation[p],"w");
            fprintf(file_data_store->HIVSTATEPOPULATIONFILE[p],"time,cd4_status,spvl,art_status,cumulative_t_earlyART,cumulative_t_ARTVS,cumulative_t_ARTVU,npartners\n");
            fclose(file_data_store->HIVSTATEPOPULATIONFILE[p]);
        }
        /* If WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER==1 then we print age distribution at some specified times, so make sure that the file is initially blank. */
        if (WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_BY_GENDER==1){
            /* The "1000" refers to the size of the array age_group_string declared above. */
            generate_demographics_byage_gender_file_header(age_group_string, 1000);
            file_data_store->AGEDISTRIBUTIONFILE[p] = fopen(file_data_store->filename_debug_agedistribution[p],"w");
            fprintf(file_data_store->AGEDISTRIBUTIONFILE[p],"%s",age_group_string);
            fclose(file_data_store->AGEDISTRIBUTIONFILE[p]);
        }

        if(WRITE_DEBUG_ART_STATE==1){
            file_data_store->ARTPOPULATIONFILE[p] = fopen(file_data_store->filename_debug_artpopulation[p],"w");
            fprintf(file_data_store->ARTPOPULATIONFILE[p],"time,n_hivneg,n_hivpos_dontknowstatus,n_hivpos_knowposneverart,n_hivpos_earlyart,n_hivpos_artvs,n_hivpos_cabo,n_hivpos_artvu,n_hivpos_dropout,n_hivpos_cascadedropout,n_artdeath,");
            fprintf(file_data_store->ARTPOPULATIONFILE[p],"cumulative_n_start_emergency_art_fromuntested,cumulative_n_start_emergency_art_fromartnaive,cumulative_n_start_emergency_art_fromartdroupout,cumulative_n_start_emergency_art_fromcascadedropout,");
            fprintf(file_data_store->ARTPOPULATIONFILE[p],"cumulative_n_learnhivpos_fromuntested,cumulative_n_startART_fromuntested,cumulative_n_startART_fromartnaive,cumulative_n_startART_fromartdropout,cumulative_n_startART_fromcascadedropout,cumulative_n_becomeVS_fromearlyart,cumulative_n_becomeVS_fromartvu,cumulative_n_becomeCABO_fromearlyart,cumulative_n_becomeCABO_fromartvu,cumulative_n_becomeVU_fromearlyart,cumulative_n_becomeVU_fromartvs,cumulative_n_ARTdropout_fromearlyart,cumulative_n_ARTdropout_fromartvs,cumulative_n_ARTdropout_fromartvu,cumulative_n_cascadedropout_fromARTnaive,n_cascadedropout_fromARTneg,cumulative_n_aidsdeaths_fromuntested,cumulative_n_aidsdeaths_fromartnaive,cumulative_n_aidsdeaths_fromearlyart,cumulative_n_aidsdeaths_fromartvs,cumulative_n_aidsdeaths_fromcabo,cumulative_n_aidsdeaths_fromartvu,cumulative_n_aidsdeaths_fromartdropout,cumulative_n_aidsdeaths_fromcascadedropout\n");
            fclose(file_data_store->ARTPOPULATIONFILE[p]);
        }

    }


    if (WRITE_DEBUG_DEMOGRAPHICS_NBIRTHS_NEWADULTS_DEATHS==1){
        file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE = fopen(file_data_store->filename_debug_nnewadults_ndeaths_file,"w");
        fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"t,");
        for(p=0;p<NPATCHES;p++)
            fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"NBirthsPatch%i,NNewAdultsPatch%i,NDeaths%i,",p,p,p);
        fprintf(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE,"\n");
        fclose(file_data_store->NBIRTHS_NNEWADULTS_NDEATHS_FILE);
    }

    /* This header is quite complicated so put in a separate function. */
    if (WRITE_DEBUG_DEMOGRAPHICS_AGE_DISTRIBUTION_ONEYEARINCKIDS==1){
        blank_one_year_age_groups_including_kids(file_data_store);
    }

    if(WRITE_PHYLOGENETICS_OUTPUT==1)
        blank_phylo_transmission_data_file(file_data_store);

    if(DEBUG_PARTNERSHIP_DURATION ==1){
        file_data_store->DUR_BETWEEN_HIGHHIGH = fopen(file_data_store->filename_DUR_BETWEEN_HIGHHIGH,"w");
        if (file_data_store->DUR_BETWEEN_HIGHHIGH==NULL){
            printf("Cannot open DUR_BETWEEN_HIGHHIGH\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_BETWEEN_HIGHHIGH);

        file_data_store->DUR_BETWEEN_MEDMED = fopen(file_data_store->filename_DUR_BETWEEN_MEDMED,"w");
        if (file_data_store->DUR_BETWEEN_MEDMED==NULL){
            printf("Cannot open DUR_BETWEEN_MEDMED\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_BETWEEN_MEDMED);

        file_data_store->DUR_BETWEEN_LOWLOW = fopen(file_data_store->filename_DUR_BETWEEN_LOWLOW,"w");
        if (file_data_store->DUR_BETWEEN_LOWLOW==NULL){
            printf("Cannot open DUR_BETWEEN_LOWLOW\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_BETWEEN_LOWLOW);

        file_data_store->DUR_WITHIN_HIGHHIGH = fopen(file_data_store->filename_DUR_WITHIN_HIGHHIGH,"w");
        if (file_data_store->DUR_WITHIN_HIGHHIGH==NULL){
            printf("Cannot open DUR_WITHIN_HIGHHIGH\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_WITHIN_HIGHHIGH);

        file_data_store->DUR_WITHIN_MEDMED = fopen(file_data_store->filename_DUR_WITHIN_MEDMED,"w");
        if (file_data_store->DUR_WITHIN_MEDMED==NULL){
            printf("Cannot open DUR_WITHIN_MEDMED\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_WITHIN_MEDMED);

        file_data_store->DUR_WITHIN_LOWLOW = fopen(file_data_store->filename_DUR_WITHIN_LOWLOW,"w");
        if (file_data_store->DUR_WITHIN_LOWLOW==NULL){
            printf("Cannot open DUR_WITHIN_LOWLOW\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->DUR_WITHIN_LOWLOW);
    }

    if(CHECK_AGE_AND_RISK_ASSORTATIVITY ==1){
        file_data_store->age_assortativity = fopen(file_data_store->filename_age_assortativity,"w");
        if (file_data_store->age_assortativity==NULL){
            printf("Cannot open age_assortativity\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->age_assortativity);

        file_data_store->age_assortativity_cross_sectional = fopen(file_data_store->filename_age_assortativity_cross_sectional,"w");
        if (file_data_store->age_assortativity_cross_sectional==NULL){
            printf("Cannot open age_assortativity_cross_sectional\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->age_assortativity_cross_sectional);

        file_data_store->risk_assortativity = fopen(file_data_store->filename_risk_assortativity,"w");
        if (file_data_store->risk_assortativity==NULL){
            printf("Cannot open risk_assortativity\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->risk_assortativity);

        file_data_store->risk_assortativity_cross_sectional = fopen(file_data_store->filename_risk_assortativity_cross_sectional,"w");
        if (file_data_store->risk_assortativity_cross_sectional==NULL){
            printf("Cannot open risk_assortativity_cross_sectional\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        fclose(file_data_store->risk_assortativity_cross_sectional);
    }

    if (WRITE_HAZARDS==1)
        blank_hazard_file(file_data_store);
}

/* Prints how long each person who is HIV+ and dies of AIDS-related illness is alive for.
 * Allows us to check that HIV duration has the right distribution.
 * Can subset to ART-naive.
 */
void write_hiv_duration(individual *dead_person, double t, file_struct *file_data_store){
    int p = dead_person->patch_no;
    file_data_store->HIVDURATIONFILE[p] = fopen(file_data_store->filename_debughivduration[p],"a");
    if (file_data_store->HIVDURATIONFILE[p]==NULL){
        printf("Cannot open output file in write_hiv_duration().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->HIVDURATIONFILE[p],"%6.4lf,%d,%ld,%8.6f,%6.4lf,%d,%d,%d\n",t,dead_person->patch_no,dead_person->id,dead_person->DEBUGTOTALTIMEHIVPOS,dead_person->t_sc,dead_person->ART_status,dead_person->SPVL_cat,dead_person->cd4);
    fclose(file_data_store->HIVDURATIONFILE[p]);
}

/* The final argument is reason for being removed from survival cohort. 1="AIDS death", 2="death from natural causes", 3="start ART". Note we don't bother with the end of the simulation for now. */
void write_hiv_duration_km(individual *indiv, double t, file_struct *file_data_store, int reason){
    int p = indiv->patch_no;
    file_data_store->HIVDURATIONFILE_KM[p] = fopen(file_data_store->filename_debughivduration_km[p],"a");
    if (file_data_store->HIVDURATIONFILE_KM[p]==NULL){
        printf("Cannot open output file in write_hiv_duration_km().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->HIVDURATIONFILE_KM[p],"%6.4lf,%6.4lf,%d,%d,%d,%6.4lf,%d\n",t,t-indiv->t_sc,reason,indiv->gender,indiv->cd4,indiv->SPVL_num_E+indiv->SPVL_num_G,indiv->SPVL_cat);
    fclose(file_data_store->HIVDURATIONFILE_KM[p]);
}

/* The "4" is the "reason" (or censoring status) - 4 means reaching the end of the simulation before dying/starting ART. */
void write_hiv_duration_km_end_of_simulation(patch_struct *patch, double t, file_struct *file_data_store){
    int p;
    long n_id;
    for (p=0;p<NPATCHES;p++){
        printf("Opening file %s\n",file_data_store->filename_debughivduration_km[p]);
        file_data_store->HIVDURATIONFILE_KM[p] = fopen(file_data_store->filename_debughivduration_km[p],"a");
        if (file_data_store->HIVDURATIONFILE_KM[p]==NULL){
            printf("Cannot open output file in write_hiv_duration_km().\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        //fprintf(file_data_store->HIVDURATIONFILE_KM[p],"Hello world start\n");
        for (n_id=0;n_id<patch[p].id_counter;n_id++){
            if ((patch[p].individual_population[n_id].cd4>DEAD) && (patch[p].individual_population[n_id].HIV_status>UNINFECTED) && (patch[p].individual_population[n_id].ART_status==ARTNAIVE || patch[p].individual_population[n_id].ART_status==ARTNEG)){
                //fprintf(file_data_store->HIVDURATIONFILE_KM[p],"Hello world\n");
                fprintf(file_data_store->HIVDURATIONFILE_KM[p],"%6.4lf,%6.4lf,4,%d,%d,%6.4lf,%d\n",t,t-patch[p].individual_population[n_id].t_sc,patch[p].individual_population[n_id].gender,patch[p].individual_population[n_id].cd4,patch[p].individual_population[n_id].SPVL_num_E+patch[p].individual_population[n_id].SPVL_num_G,patch[p].individual_population[n_id].SPVL_cat);
                //printf("%6.4lf,%6.4lf,4,%d,%d,%6.4lf,%d\n",t,t-patch[p].individual_population[n_id].t_sc,patch[p].individual_population[n_id].gender,patch[p].individual_population[n_id].cd4,patch[p].individual_population[n_id].SPVL_num_E+patch[p].individual_population[n_id].SPVL_num_G,patch[p].individual_population[n_id].SPVL_cat);
            }
        }
        fclose(file_data_store->HIVDURATIONFILE_KM[p]);
    }
}

void write_cd4_at_seroconversion(individual *indiv, file_struct *file_data_store){
    int p = indiv->patch_no;
    file_data_store->HIVCD4_AFTER_SEROCONVERSION[p] = fopen(file_data_store->filename_debughivcd4_after_seroconversion[p],"a");
    if (file_data_store->HIVCD4_AFTER_SEROCONVERSION[p]==NULL){
        printf("Cannot open output file in write_hiv_duration_km().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->HIVCD4_AFTER_SEROCONVERSION[p],"%li,%i,%i\n",indiv->id,indiv->cd4,indiv->SPVL_cat);
    fclose(file_data_store->HIVCD4_AFTER_SEROCONVERSION[p]);
}


void write_initial_spvl_distribution(individual *seeded_infection, file_struct *file_data_store){
    int p = seeded_infection->patch_no;
    file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p] = fopen(file_data_store->filename_debuginitial_spvl_distribution[p],"a");
    if (file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p]==NULL){
        printf("Cannot open output file in write_hiv_duration_km().\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    fprintf(file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p],"%i,%6.4lf,%6.4lf\n",seeded_infection->SPVL_cat,seeded_infection->SPVL_num_E,seeded_infection->SPVL_num_G);
    fclose(file_data_store->HIV_INITIAL_SPVL_DISTRIBUTION[p]);


}

void write_cd4_spvl_states(patch_struct *patch, int year, file_struct *file_data_store){
    int p;
    long n_id;
    int art_status, cd4_status, npartners;
    double spvl, t_early_art, t_vs, t_vu;
    for (p=0;p<NPATCHES;p++){
        file_data_store->HIVSTATEPOPULATIONFILE[p] = fopen(file_data_store->filename_debug_hivpopulation[p],"a");
        for (n_id=0;n_id<patch[p].id_counter;n_id++){
            /* Only get data for HIV+ who are alive. */
            if ((patch[p].individual_population[n_id].cd4>DEAD) && (patch[p].individual_population[n_id].HIV_status>UNINFECTED)){
                art_status = patch[p].individual_population[n_id].ART_status;
                cd4_status = patch[p].individual_population[n_id].cd4;
                spvl = patch[p].individual_population[n_id].SPVL_num_E+patch[p].individual_population[n_id].SPVL_num_G;
                npartners = patch[p].individual_population[n_id].n_partners;
                t_early_art = patch[p].individual_population[n_id].DEBUG_cumulative_time_on_ART_early;
                t_vs = patch[p].individual_population[n_id].DEBUG_cumulative_time_on_ART_VS;
                t_vu = patch[p].individual_population[n_id].DEBUG_cumulative_time_on_ART_VU;
                fprintf(file_data_store->HIVSTATEPOPULATIONFILE[p],"%d,%d,%6.4lf,%d,%6.4lf,%6.4lf,%6.4lf,%i\n",year,cd4_status,spvl,art_status,t_early_art,t_vs,t_vu,npartners);
            }
        }
        fclose(file_data_store->HIVSTATEPOPULATIONFILE[p]);
    }
}


void print_debugerror_shouldbezero_exit(char *varname){
    printf("ERROR: Variable %s should be zero. Exiting\n",varname);
    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
    fflush(stdout);
    exit(1);
}

void write_art_states(patch_struct *patch, int year, debug_struct *debug, file_struct *file_data_store){
    int p;
    long n_id;
    int art_i;
    int art_status;

    long n_hivneg;
    long n_hivpos_dontknowstatus;
    long n_hivpos_knowposneverart;
    long n_hivpos_earlyart;
    long n_hivpos_artvs;
    long n_hivpos_cabo;
    long n_hivpos_artvu;
    long n_hivpos_dropout;
    long n_hivpos_cascadedropout;
    long n_artdeath;

    /****** Now look at transitions. ******/
    /* First AIDS deaths - should ONLY be possible if not ART VS. */
    long n_aidsdeaths_fromuntested,n_aidsdeaths_fromartnaive,n_aidsdeaths_fromearlyart,n_aidsdeaths_fromartvs,n_aidsdeaths_fromcabo,n_aidsdeaths_fromartvu,n_aidsdeaths_fromartdropout,n_aidsdeaths_fromcascadedropout;
    /* Next number of HIV+ who learn status (go from n_hivpos_dontknowstatus to n_hivpos_knowposneverart). */
    long n_learnhivpos_fromuntested;
    /* Number of people who get CD4 tested. */
    //long n_getcd4test_fromuntested,n_getcd4test_fromartnaive,n_getcd4test_fromartdropout,n_getcd4test_fromcascadedropout;
    /* Number of people who get ART tested: */
    long n_startART_fromuntested,n_startART_fromartnaive,n_startART_fromartdropout,n_startART_fromcascadedropout;
    /* Number of people becoming virally suppressed or cabo or virally unsuppressed: */
    long n_becomeVS_fromearlyart, n_becomeVS_fromartvu,n_becomeCABO_fromearlyart, n_becomeCABO_fromartvu, n_becomeVU_fromearlyart, n_becomeVU_fromartvs;
    /* Number of people stopping ART: */
    long n_ARTdropout_fromearlyart,n_ARTdropout_fromartvs,n_ARTdropout_fromartvu;
    /* Number of people dropping out of cascade before ART: */
    long n_cascadedropout_fromARTnaive, n_cascadedropout_fromARTneg;

    /* ART_status runs from -1 to 7 (ARTDEATH) so need array to go from 0 to ARTDEATH+2. */
    long temp_state_counter[ARTDEATH+2];
    long temp_hivnegstate_counter;




    for (p=0;p<NPATCHES;p++){
        for (art_i=0;art_i<(ARTDEATH+2);art_i++)
            temp_state_counter[art_i] = 0;
        temp_hivnegstate_counter = 0;
        for (n_id=0;n_id<patch[p].id_counter;n_id++){
            if (patch[p].individual_population[n_id].cd4>DEAD){
                art_status = patch[p].individual_population[n_id].ART_status;
                if (art_status<ARTNEG || art_status>=ARTDEATH){
                    printf("ERROR - undefined ART status for individual %ld in patch %d. Exiting\n",patch[p].individual_population[n_id].id,p);
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }

                temp_state_counter[art_status+1]++;
                if(patch[p].individual_population[n_id].HIV_status==UNINFECTED)
                    temp_hivnegstate_counter++;
            }
            else{ /* Check if dead person died while on ART: */
                if(patch[p].individual_population[n_id].ART_status==ARTDEATH)
                    temp_state_counter[ARTDEATH+1]++;
            }
        }
        /* Translate these into more readable variable names for clarity - this is not needed for coded. */
        n_hivneg = temp_hivnegstate_counter;
        n_hivpos_dontknowstatus = temp_state_counter[ARTNEG+1]-temp_hivnegstate_counter;
        n_hivpos_knowposneverart = temp_state_counter[ARTNAIVE+1];
        n_hivpos_earlyart = temp_state_counter[EARLYART+1];
        n_hivpos_artvs = temp_state_counter[LTART_VS+1];
        n_hivpos_cabo = temp_state_counter[CABO+1];
        n_hivpos_artvu = temp_state_counter[LTART_VU+1];
        n_hivpos_dropout = temp_state_counter[ARTDROPOUT+1];
        n_hivpos_cascadedropout = temp_state_counter[CASCADEDROPOUT+1];
        n_artdeath = temp_state_counter[ARTDEATH+1];

        /* Now look at transitions between states. Again, translate first for clarity: */
        if (debug->art_vars[p].n_start_emergency_art != (debug->art_vars[p].n_start_emergency_art_fromuntested + debug->art_vars[p].n_start_emergency_art_fromartnaive+ debug->art_vars[p].n_start_emergency_art_fromartdroupout+ debug->art_vars[p].n_start_emergency_art_fromcascadedropout)){
            printf("Error - number starting emergency ART in a given timestep doesn't add up. Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        n_learnhivpos_fromuntested = debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTNAIVE+1];
        if (debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTNEG+1]>0) /* Impossible */
            print_debugerror_shouldbezero_exit("debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTNEG+1]");
        for (art_i=EARLYART;art_i<=ARTDEATH;art_i++){
            if(debug->art_vars[p].cascade_transitions[art_i+1][ARTNEG+1]!=0){
                char errorstring[50];
                sprintf(errorstring,"art_vars[p].cascade_transitions[%i][ARTNEG+1] %li",art_i,debug->art_vars[p].cascade_transitions[art_i+1][ARTNEG+1]);
                print_debugerror_shouldbezero_exit(errorstring);
            }
        }

        //      n_getcd4test_fromuntested = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTNEG+1];
        //      n_getcd4test_fromartnaive = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTNAIVE+1];
        //      n_getcd4test_fromartdropout = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDROPOUT+1];
        //      n_getcd4test_fromcascadedropout = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][CASCADEDROPOUT+1];
        //      /* Check for impossible transitions: */
        //      for (art_i=EARLYART+1;art_i<ARTDROPOUT+1;art_i++){
        //          if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][art_i+1]!=0){
        //              char errorstring[50];
        //              sprintf(errorstring,"art_vars[p].cascade_transitions[ARTNAIVE+1][%i]",art_i);
        //              print_debugerror_shouldbezero_exit(errorstring);
        //          }
        //      }
        //      if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDEATH+1]!=0)
        //          print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDEATH+1]");



        n_startART_fromuntested = debug->art_vars[p].cascade_transitions[ARTNEG+1][EARLYART+1];
        n_startART_fromartnaive = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][EARLYART+1];
        n_startART_fromartdropout = debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][EARLYART+1];
        n_startART_fromcascadedropout = debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][EARLYART+1];
        /* Check for impossible transitions: */
        for (art_i=EARLYART;art_i<ARTDROPOUT;art_i++){
            if(debug->art_vars[p].cascade_transitions[art_i+1][EARLYART+1]!=0){
                char errorstring[50];
                sprintf(errorstring,"art_vars[p].cascade_transitions[%i][EARLYART+1]",art_i);
                print_debugerror_shouldbezero_exit(errorstring);
            }
        }
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][EARLYART+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][EARLYART+1]");

        // Now becoming VS - should only be possible from early ART or VU:
        n_becomeVS_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][LTART_VS+1];
        n_becomeVS_fromartvu = debug->art_vars[p].cascade_transitions[LTART_VU+1][LTART_VS+1];
        if(debug->art_vars[p].cascade_transitions[ARTNEG+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNAIVE+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VS+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VS+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[CABO+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CABO+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][LTART_VS+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][LTART_VS+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][LTART_VS+1]");

        // Now becoming CABO - should only be possible from early ART or VU:
        n_becomeCABO_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][CABO+1];
        n_becomeCABO_fromartvu = debug->art_vars[p].cascade_transitions[LTART_VU+1][CABO+1];
        if(debug->art_vars[p].cascade_transitions[ARTNEG+1][CABO+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][CABO+1]");
        if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][CABO+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNAIVE+1][CABO+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VS+1][CABO+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VS+1][CABO+1]");
        if(debug->art_vars[p].cascade_transitions[CABO+1][CABO+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CABO+1][CABO+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][CABO+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][CABO+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][CABO+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][CABO+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][CABO+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][CABO+1]");

        // Now becoming VU - should only be possible from early ART or VS:
        n_becomeVU_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][LTART_VU+1];
        n_becomeVU_fromartvs = debug->art_vars[p].cascade_transitions[LTART_VS+1][LTART_VU+1];
        if(debug->art_vars[p].cascade_transitions[ARTNEG+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNAIVE+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[CABO+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CABO+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VU+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VU+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][LTART_VU+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][LTART_VU+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][LTART_VU+1]");

        // Now dropping out from ART - should only be possible from early ART, VU or VS:
        n_ARTdropout_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][ARTDROPOUT+1];
        n_ARTdropout_fromartvs = debug->art_vars[p].cascade_transitions[LTART_VS+1][ARTDROPOUT+1];
        n_ARTdropout_fromartvu = debug->art_vars[p].cascade_transitions[LTART_VU+1][ARTDROPOUT+1];
        if(debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[CABO+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CABO+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][ARTDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][ARTDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][ARTDROPOUT+1]");

        // Now dropping out from cascade before ART - should only be possible from ARTNAIVE or ARTNEG:
        n_cascadedropout_fromARTnaive = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][CASCADEDROPOUT+1];
        n_cascadedropout_fromARTneg = debug->art_vars[p].cascade_transitions[ARTNEG+1][CASCADEDROPOUT+1];
        //if(debug->art_vars[p].cascade_transitions[ARTNEG+1][CASCADEDROPOUT+1]!=0)
        //  print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTNEG+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[EARLYART+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[EARLYART+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VS+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VS+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[CABO+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CABO+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[LTART_VU+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[LTART_VU+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDROPOUT+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[CASCADEDROPOUT+1][CASCADEDROPOUT+1]");
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][CASCADEDROPOUT+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][CASCADEDROPOUT+1]");

        n_aidsdeaths_fromuntested = debug->art_vars[p].cascade_transitions[ARTNAIVE+1][ARTDEATH+1];
        n_aidsdeaths_fromartnaive = debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTDEATH+1];
        n_aidsdeaths_fromearlyart = debug->art_vars[p].cascade_transitions[EARLYART+1][ARTDEATH+1];
        n_aidsdeaths_fromartvs = debug->art_vars[p].cascade_transitions[LTART_VS+1][ARTDEATH+1];
        n_aidsdeaths_fromcabo = debug->art_vars[p].cascade_transitions[CABO+1][ARTDEATH+1];
        n_aidsdeaths_fromartvu = debug->art_vars[p].cascade_transitions[LTART_VU+1][ARTDEATH+1];
        n_aidsdeaths_fromcascadedropout = debug->art_vars[p].cascade_transitions[CASCADEDROPOUT+1][ARTDEATH+1];
        n_aidsdeaths_fromartdropout = debug->art_vars[p].cascade_transitions[ARTDROPOUT+1][ARTDEATH+1];
        if(debug->art_vars[p].cascade_transitions[ARTDEATH+1][ARTDEATH+1]!=0)
            print_debugerror_shouldbezero_exit("art_vars[p].cascade_transitions[ARTDEATH+1][ARTDEATH+1]");

        file_data_store->ARTPOPULATIONFILE[p] = fopen(file_data_store->filename_debug_artpopulation[p],"a");
        fprintf(file_data_store->ARTPOPULATIONFILE[p],"%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,",year,n_hivneg,n_hivpos_dontknowstatus,n_hivpos_knowposneverart,n_hivpos_earlyart,n_hivpos_artvs,n_hivpos_cabo,n_hivpos_artvu,n_hivpos_dropout,n_hivpos_cascadedropout,n_artdeath);
        fprintf(file_data_store->ARTPOPULATIONFILE[p],"%ld,%ld,%ld,%ld,",debug->art_vars[p].n_start_emergency_art_fromuntested , debug->art_vars[p].n_start_emergency_art_fromartnaive, debug->art_vars[p].n_start_emergency_art_fromartdroupout, debug->art_vars[p].n_start_emergency_art_fromcascadedropout);
        fprintf(file_data_store->ARTPOPULATIONFILE[p],"%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",n_learnhivpos_fromuntested,n_startART_fromuntested,n_startART_fromartnaive,n_startART_fromartdropout,n_startART_fromcascadedropout,n_becomeVS_fromearlyart, n_becomeVS_fromartvu,n_becomeCABO_fromearlyart,n_becomeCABO_fromartvu,n_becomeVU_fromearlyart, n_becomeVU_fromartvs,n_ARTdropout_fromearlyart,n_ARTdropout_fromartvs,n_ARTdropout_fromartvu,n_cascadedropout_fromARTnaive,n_cascadedropout_fromARTneg, n_aidsdeaths_fromuntested,n_aidsdeaths_fromartnaive,n_aidsdeaths_fromearlyart,n_aidsdeaths_fromartvs,n_aidsdeaths_fromcabo,n_aidsdeaths_fromartvu,n_aidsdeaths_fromartdropout,n_aidsdeaths_fromcascadedropout);
        fclose(file_data_store->ARTPOPULATIONFILE[p]);

    }

}


/****************************************************************************************************************
 * Stuff related to CHiPs visits:
 ****************************************************************************************************************/

void reset_annual_chips_visit_counter(age_list_struct *age_list){
    int g,ai,n;
    for (g=0;g<N_GENDER;g++){
        for(ai=0;ai<(MAX_AGE-AGE_ADULT); ai++)
            for(n=0;n<age_list->age_list_by_gender[g]->number_per_age_group[ai];n++)
                age_list->age_list_by_gender[g]->age_group[ai][n]->VISITED_BY_CHIPS_THISROUND = FALSE;
        for(n=0;n<age_list->age_list_by_gender[g]->number_oldest_age_group;n++)
            age_list->age_list_by_gender[g]->oldest_age_group[n]->VISITED_BY_CHIPS_THISROUND = FALSE;
    }
}





/* This one only knows about people who are alive at the end of the round. */
void print_chips_statistics_using_age_list(age_list_struct *age_list, double t){
    int g,aa,ai,i, v;
    long n_visited_by_chips_this_round = 0;
    long total_chips_visits = 0;

    /* Fairly arbitrarily, we store the number of people who have been visited 0, 1, 2, 3, 4, 5+ times.The upper limit should
     * be bigger than the 4 visits scheduled in PopART, but in principle this could be more than 5 visits - you just need
     * to change the size of MAXVISITSRECORDED below. */
    int MAXVISITSRECORDED = 5;
    long *nvisits_distribution[N_GENDER];
    long n_chips_start_art = 0;
    long denom_chips_visits[N_GENDER]; /* Number of M/F eligible to be visited by ChiPs this round. */
    long denom_chips_visits_total = 0; /* Sum by gender of denom_chips_visits. */
    individual *person;            /* Temporary pointer to person currently being used, to make code more readable. As pointing to existing memory no malloc used. */
    for (g=0;g<N_GENDER;g++){
        nvisits_distribution[g] = malloc((MAXVISITSRECORDED+1)*sizeof(long)); /* Note that we go from 0..MAXVISITSRECORDED. */
        denom_chips_visits[g] = 0;
    }

    for (g=0;g<N_GENDER;g++)
        /* Again, we go from 0..MAXVISITSRECORDED - hence "<=" rather than "<". */
        for (v=0;v<=MAXVISITSRECORDED;v++)
            nvisits_distribution[g][v] = 0;

    int NDIED = 0;

    for (g=0;g<N_GENDER;g++){
        for (aa=(AGE_CHIPS-AGE_ADULT); aa<(MAX_AGE-AGE_ADULT); aa++){
            ai = age_list->age_list_by_gender[g]->youngest_age_group_index + aa; /* ai is the index of the array age_list->number_per_age_group of the age group of people you want to be dead */
            while (ai>(MAX_AGE-AGE_ADULT-1))
                ai = ai - (MAX_AGE-AGE_ADULT);
            for(i=0;i<age_list->age_list_by_gender[g]->number_per_age_group[ai];i++){
                person = age_list->age_list_by_gender[g]->age_group[ai][i];

                /* Note that there is a slight issue that people may be visited successfully by chips, then die (or be scheduled to receive a visit but die beforehand). To get around this use a special value for VISITED_BY_CHIPS_THISROUND. */
                if (person->VISITED_BY_CHIPS_THISROUND>DIEDBEFORECHIPSVISIT){
                    denom_chips_visits[g]++;
                    /* Since VISITED_BY_CHIPS_THISROUND=0 if not visited, and 1 if visited, can sum over this: */
                    if (person->VISITED_BY_CHIPS_THISROUND!=0 && person->VISITED_BY_CHIPS_THISROUND!=1)
                        printf("Error: at t=%6.4lf person.VISITED_BY_CHIPS_THISROUND=%i person.id=%li person.NCHIPSVISITS = %i\n",t,person->VISITED_BY_CHIPS_THISROUND,person->id,person->NCHIPSVISITS);
                    n_visited_by_chips_this_round += person->VISITED_BY_CHIPS_THISROUND;
                    /* v is the lifetime number of visits this individual has had by CHiPs. */
                    v = person->NCHIPSVISITS;
                    total_chips_visits += v; /* Use this to calculate mean number of visits per person in CHiPs denominator. */
                    if (v>=MAXVISITSRECORDED)
                        nvisits_distribution[g][MAXVISITSRECORDED]++;
                    else
                        nvisits_distribution[g][v]++;
                    if(person->VISITEDBYCHIPS_TO_INIT_ART)
                        n_chips_start_art++;
                }
                else{
                    NDIED++;
                }
            }
        }
    }

    printf("Number died = %i\n",NDIED);

    printf("%6.4lf n_thisround=%ld ",t,n_visited_by_chips_this_round);
    for (g=0;g<N_GENDER;g++){
        denom_chips_visits_total += denom_chips_visits[g];
    //  printf("%ld ",denom_chips_visits[g]);
    }
    for (v=0;v<=MAXVISITSRECORDED;v++)
        for (g=0;g<N_GENDER;g++)
            printf("%ld ",nvisits_distribution[g][v]);
    printf("\n");

    /* Free memory. Note - it's really important to NOT free person (as this memory is in use and was not locally allocated). */
    for (g=0;g<N_GENDER;g++)
        free(nvisits_distribution[g]);
}


void print_chips_statistics_using_chipsonly(patch_struct *patch, int p, double t){
    int g,ac,i, v;
    long n_visited_by_chips_this_round = 0;
    long total_chips_visits = 0;

    /* Fairly arbitrarily, we store the number of people who have been visited 0, 1, 2, 3, 4, 5+ times.The upper limit should
     * be bigger than the 4 visits scheduled in PopART, but in principle this could be more than 5 visits - you just need
     * to change the size of MAXVISITSRECORDED below. */
    int MAXVISITSRECORDED = 5;
    long *nvisits_distribution[N_GENDER];
    long n_chips_start_art = 0;
    long denom_chips_visits[N_GENDER]; /* Number of M/F eligible to be visited by ChiPs this round. */
    //long denom_chips_visits_total = 0; /* Sum by gender of denom_chips_visits. */
    for (g=0;g<N_GENDER;g++){
        nvisits_distribution[g] = malloc((MAXVISITSRECORDED+1)*sizeof(long)); /* Note that we go from 0..MAXVISITSRECORDED. */
        denom_chips_visits[g] = 0;
    }


    for (g=0;g<N_GENDER;g++)
        /* Again, we go from 0..MAXVISITSRECORDED - hence "<=" rather than "<". */
        for (v=0;v<=MAXVISITSRECORDED;v++)
            nvisits_distribution[g][v] = 0;

    for (g=0;g<N_GENDER;g++){
        /* Run from AGE_CHIPS to 80+. */
        for (ac=0; ac<(MAX_AGE-AGE_CHIPS+1); ac++){
            for (i=0;i<patch[p].chips_sample->number_to_visit[g][ac];i++){
                //
                denom_chips_visits[g]++;

                //              if (patch[p].chips_sample->list_ids_to_visit[g][ac][i]==9326)
                //                  printf("Checking %li. Status = %i\n",patch[p].chips_sample->list_ids_to_visit[g][ac][i],patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].VISITED_BY_CHIPS_THISROUND);

                //              if (patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].VISITED_BY_CHIPS_THISROUND!=0 && patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].VISITED_BY_CHIPS_THISROUND!=1)
                //                  printf("Error: at t=%6.4lf person.VISITED_BY_CHIPS_THISROUND=%i person.id=%li person.NCHIPSVISITS = %i\n",t,patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].VISITED_BY_CHIPS_THISROUND,patch[p].chips_sample->list_ids_to_visit[g][ac][i],patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].NCHIPSVISITS);

                /* Since VISITED_BY_CHIPS_THISROUND=0 if not visited, and 1 if visited, can sum over this: */
                n_visited_by_chips_this_round += patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].VISITED_BY_CHIPS_THISROUND;
                /* v is the lifetime number of visits this individual has had by CHiPs. */
                v = patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].NCHIPSVISITS;
                total_chips_visits += v; /* Use this to calculate mean number of visits per person in CHiPs denominator. */
                if (v>=MAXVISITSRECORDED)
                    nvisits_distribution[g][MAXVISITSRECORDED]++;
                else
                    nvisits_distribution[g][v]++;
                if(patch[p].individual_population[patch[p].chips_sample->list_ids_to_visit[g][ac][i]].VISITEDBYCHIPS_TO_INIT_ART)
                    n_chips_start_art++;
            }
        }
    }

    printf("%6.4lf n_thisround=%ld ",t,n_visited_by_chips_this_round);
    //for (g=0;g<N_GENDER;g++)
    //  printf("%ld ",denom_chips_visits[g]);
    for (v=0;v<=MAXVISITSRECORDED;v++)
        for (g=0;g<N_GENDER;g++)
            printf("%ld ",nvisits_distribution[g][v]);
    printf("\n");

    /* Free memory. Note - it's really important to NOT free person (as this memory is in use and was not locally allocated). */
    for (g=0;g<N_GENDER;g++)
        free(nvisits_distribution[g]);

}



void output_hazard_over_time_period(double t,double hazard, individual *susceptible, individual *pos_partner, file_struct *file_data_store, output_struct *output){
    int infector_hiv_cd4_acute;
    char temp_string_hazard[100];

    /* Merge acute and CD4 to reduce output: */
    if (pos_partner->HIV_status==ACUTE)
        infector_hiv_cd4_acute = -1;
    else
        infector_hiv_cd4_acute = pos_partner->cd4;

    sprintf(temp_string_hazard,"%lf,%8.6lf,%i,%i,%i,%i,%i,%i,%i,%8.6lf\n",hazard,t,susceptible->circ,pos_partner->gender,infector_hiv_cd4_acute,  pos_partner->ART_status, pos_partner->patch_no,  susceptible->sex_risk,  pos_partner->sex_risk,pos_partner->SPVL_num_E+pos_partner->SPVL_num_G);
    //printf("%s",temp_string_hazard);

    /* The -2 is because in C the last character in any string of length n is "\0" - so we only have n-1 characters in the array we can write to. Make -2 instead of -1 to be a bit more sure! */
    if((strlen(output->hazard_output_string)+strlen(temp_string_hazard))>(HAZARD_OUTPUT_STRING_LENGTH-2)){
        write_hazard_data(file_data_store, output->hazard_output_string);
        /* Empty output->phylogenetics_output_string. To do this we just need to set the first character to be '\0'. */
        (output->hazard_output_string)[0] = '\0';
        // DEBUG - checks that string length is reset to zero.
        //printf("New length = %lu\n",strlen(output->hazard_output_string));
        //printf("Error - need to increase size of HAZARD_OUTPUT_STRING_LENGTH. Exiting\n");
        //printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        //fflush(stdout);
        //exit(1);
    }

    /* Add to existing hazard output string. */
    strcat(output->hazard_output_string,temp_string_hazard);
}


void write_hazard_data(file_struct *file_data_store, char *hazard_output_string){
    file_data_store->HAZARD_FILE = fopen(file_data_store->filename_hazard_output,"a");
    fprintf(file_data_store->HAZARD_FILE,"%s",hazard_output_string);
    fclose(file_data_store->HAZARD_FILE);
}

void blank_hazard_file(file_struct *file_data_store){
    file_data_store->HAZARD_FILE = fopen(file_data_store->filename_hazard_output,"w");
    fprintf(file_data_store->HAZARD_FILE,"Hazard,t,SusceptibleCircStatus,PartnerGender,PartnerHIVstatus,PartnerARTstatus,PartnerPatchNumber,SusceptibleRiskGp,PartnerRiskGp,SPVL\n");
    fclose(file_data_store->HAZARD_FILE);
}


