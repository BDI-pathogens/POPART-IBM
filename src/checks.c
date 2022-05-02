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

#include "checks.h"
#include "utilities.h"
#include "partnership.h"
#include "demographics.h"
#include "hiv.h"

// Function: check_partnership_formation()
// 
// Creates 2 women and 2 men, forms partnerships between them and prints output.
// Memory for these is allocated and freed inside the function.
//
// Parameters
// ----------
//
// overall_partnerships: pointer to all_partnerships structure
// param: pointer to parameters structure
// debug : pointer to debug_struct structure
// file_data_store : pointer to a file_struct structure
//
// Returns
// -------
// Nothing; performs a test

void check_partnership_formation(all_partnerships *overall_partnerships, parameters *param, 
    debug_struct *debug, file_struct *file_data_store){
    
    printf("-------------------------------------\n");
    printf("Check partnership formation:\n");

    individual *indiv1, *indiv2, *indiv3, *indiv4;
    indiv1 = malloc(sizeof(individual));
    indiv2 = malloc(sizeof(individual));
    indiv3 = malloc(sizeof(individual));
    indiv4 = malloc(sizeof(individual));

    indiv1->id = 1;
    indiv2->id = 2;
    indiv3->id = 3;
    indiv4->id = 4;

    indiv1->idx_serodiscordant = -1;
    indiv2->idx_serodiscordant = -1;
    indiv3->idx_serodiscordant = -1;
    indiv4->idx_serodiscordant = -1;

    indiv1->gender = 1;
    indiv2->gender = 1;
    indiv3->gender = 0;
    indiv4->gender = 0;

    indiv1->HIV_status = 2;
    indiv2->HIV_status = 0;
    indiv3->HIV_status = 0;
    indiv4->HIV_status = 0;

    indiv1->n_partners = 0;
    indiv2->n_partners = 0;
    indiv3->n_partners = 0;
    indiv4->n_partners = 0;

    indiv1->n_HIVpos_partners = 0;
    indiv2->n_HIVpos_partners = 0;
    indiv3->n_HIVpos_partners = 0;
    indiv4->n_HIVpos_partners = 0;

    printf("- Before partnership formation:\n");

    print_partners(indiv1);
    print_partners(indiv2);
    print_partners(indiv3);
    print_partners(indiv4);

    /* forming partnership between them and checking things are OK */
    new_partnership( indiv1, indiv3,param->start_time_simul, 
        overall_partnerships, param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;
    new_partnership( indiv1, indiv4, param->start_time_simul, 
        overall_partnerships, param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;
    new_partnership( indiv2, indiv4, param->start_time_simul, 
        overall_partnerships, param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;

    printf("- After partnership formation:\n");

    print_partners(indiv1);
    print_partners(indiv2);
    print_partners(indiv3);
    print_partners(indiv4);

    free(indiv1);
    free(indiv2);
    free(indiv3);
    free(indiv4);
}


// Function: check_partnership_formation_and_HIV_acquisition()
// Same as check_partnership_formation() but with possible HIV transmission
//
// Creates 2 women (indiv1 and indiv2) and 2 men (indiv3 and indiv4)
//
// Parameters
// ----------
// patch_struct : pointer to a patch structure
// p : int
// all_partnerships : pointer to an overall_partnerships structure
// output_struct : pointer to an output structure
// debug_struct : pointer to a debug structure
// file_struct : pointer to a file_data_store structure
//
//
// Returns
// -------
// Nothing; performs a check

void check_partnership_formation_and_HIV_acquisition(patch_struct *patch, int p, 
    all_partnerships *overall_partnerships, output_struct *output, debug_struct *debug, 
    file_struct *file_data_store){

    int i;

    printf("-------------------------------------\n");
    printf("Check partnership formation and HIV acquisition:\n");

    individual *indiv1, *indiv2, *indiv3, *indiv4;
    indiv1 = malloc(sizeof(individual));
    indiv2 = malloc(sizeof(individual));
    indiv3 = malloc(sizeof(individual));
    indiv4 = malloc(sizeof(individual));

    indiv1->id = 1;
    indiv2->id = 2;
    indiv3->id = 3;
    indiv4->id = 4;

    indiv1->idx_serodiscordant = -1;
    indiv2->idx_serodiscordant = -1;
    indiv3->idx_serodiscordant = -1;
    indiv4->idx_serodiscordant = -1;

    indiv1->gender = 1;
    indiv2->gender = 1;
    indiv3->gender = 0;
    indiv4->gender = 0;

    indiv1->HIV_status = 1;
    indiv2->HIV_status = 0;
    indiv3->HIV_status = 0;
    indiv4->HIV_status = 0;

    indiv1->cd4 = 1;
    indiv2->cd4 = CD4_UNINFECTED;
    indiv3->cd4 = CD4_UNINFECTED;
    indiv4->cd4 = CD4_UNINFECTED;

    /* Genetic component of  log10(SPVL). Take a dummy value so that 
    total SPVL (=SPVL_num_G+SPVL_num_E) is 4.0 and 33% of SPVL is due to genetics. */
    indiv1->SPVL_num_G = 1.33;
    indiv1->SPVL_num_E = 2.67; /* Environmental component of  log10(SPVL). Take a dummy value */
    indiv1->SPVL_infector = -1; /* Dummy value to signify that seeded infection. */
    indiv1->SPVL_cat = get_spvl_cat(indiv1->SPVL_num_E + indiv1->SPVL_num_G);
    indiv1->cd4 = 3;

    indiv1->DEBUG_cumulative_time_on_ART_VS = 0;
    indiv1->DEBUG_cumulative_time_on_ART_VU = 0;
    indiv1->DEBUG_cumulative_time_on_ART_early = 0;

    indiv1->n_partners = 0;
    indiv2->n_partners = 0;
    indiv3->n_partners = 0;
    indiv4->n_partners = 0;

    indiv1->n_HIVpos_partners = 0;
    indiv2->n_HIVpos_partners = 0;
    indiv3->n_HIVpos_partners = 0;
    indiv4->n_HIVpos_partners = 0;

    indiv1->patch_no = 0;
    indiv2->patch_no = 0;
    indiv3->patch_no = 0;
    indiv4->patch_no = 0;

    printf("-------------------------------------\n");
    printf("Initial state\n");
    printf("-------------------------------------\n");

    print_HIV_status(indiv1);
    print_HIV_status(indiv2);
    print_HIV_status(indiv3);
    print_HIV_status(indiv4);

    /* forming partnership between them and checking things are OK */
    new_partnership(indiv1, indiv3, patch[p].param->start_time_simul, 
        overall_partnerships, patch[p].param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;
    new_partnership(indiv1, indiv4, patch[p].param->start_time_simul, 
        overall_partnerships, patch[p].param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;
    new_partnership(indiv2, indiv4, patch[p].param->start_time_simul, 
        overall_partnerships, patch[p].param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;

    print_partnership(&overall_partnerships->partner_pairs[0]);
    print_partnership(&overall_partnerships->partner_pairs[1]);
    print_partnership(&overall_partnerships->partner_pairs[2]);

    print_partners(indiv1);
    print_partners(indiv2);
    print_partners(indiv3);
    print_partners(indiv4);


    printf("List of serodiscordant partnerships:\n");
    for(i = 0; i < *overall_partnerships->n_partnerships; i++){
        if(is_serodiscordant(&overall_partnerships->partner_pairs[i])){
            print_partnership(&overall_partnerships->partner_pairs[i]);
        }
    }

    /* HIV acquisition processes and then checking whether some people have become infected */
    /* doing a long loop to check that by chance one of the susceptibles ends up being infected */
    for(i = 0; i < 100; i++){
        if(indiv2->HIV_status==0 && indiv2->n_HIVpos_partners>0){
            hiv_acquisition(indiv2, patch[p].param->start_time_hiv, 
                patch, p, overall_partnerships, output, debug, file_data_store, 0, 0);
        }
        if(indiv3->HIV_status==0 && indiv3->n_HIVpos_partners>0){
            hiv_acquisition(indiv3, patch[p].param->start_time_hiv, 
                patch, p, overall_partnerships, output, debug, file_data_store, 0, 0);
        }
        if(indiv4->HIV_status==0 && indiv4->n_HIVpos_partners>0){
            hiv_acquisition(indiv4, patch[p].param->start_time_hiv, 
                patch, p, overall_partnerships, output, debug, file_data_store, 0, 0);
        }
    }

    printf("-------------------------------------\n");
    printf("Final state\n");
    printf("-------------------------------------\n");

    print_partnership(&overall_partnerships->partner_pairs[0]);
    print_partnership(&overall_partnerships->partner_pairs[1]);
    print_partnership(&overall_partnerships->partner_pairs[2]);

    print_HIV_status(indiv1);
    print_HIV_status(indiv2);
    print_HIV_status(indiv3);
    print_HIV_status(indiv4);

    print_partners(indiv1);
    print_partners(indiv2);
    print_partners(indiv3);
    print_partners(indiv4);

    printf("List of serodiscordant partnerships:\n");
    for(i = 0; i < *overall_partnerships->n_partnerships; i++){
        if(is_serodiscordant(&overall_partnerships->partner_pairs[i])){
            print_partnership(&overall_partnerships->partner_pairs[i]);
        }
    }

    printf("-------------------------------------\n");
    printf("Phylogenetic output:\n");
    printf("IdInfector IdInfected IsInfectorAcute\n");
    for(int p=0; p<NPATCHES; p++){
        printf("%s\n", output->phylogenetics_output_string[p]);
    }
    printf("-------------------------------------\n");

    free(indiv1);
    free(indiv2);
    free(indiv3);
    free(indiv4);
}


// Function: check_partnership_dissolution()
// 
// Create 2 women and 2 men, forms partnerships, then dissolves some of them (at a time NOT
// given by the duration of the partnerships, so e.g. this is what would happen if one of the
// partners die).
// 
// 

void check_partnership_dissolution(all_partnerships *overall_partnerships, parameters *param, 
    debug_struct *debug, file_struct *file_data_store){
    
    long initial_n_partnerships = *overall_partnerships->n_partnerships;
    int i;
    /* creating 2 women (indiv1 and indiv2) and 2 men (indiv3 and indiv4) */

    printf("-------------------------------------\n");
    printf("Check partnership dissolution:\n");

    individual *indiv1, *indiv2, *indiv3, *indiv4;
    indiv1 = malloc(sizeof(individual));
    indiv2 = malloc(sizeof(individual));
    indiv3 = malloc(sizeof(individual));
    indiv4 = malloc(sizeof(individual));

    indiv1->id = 1;
    indiv2->id = 2;
    indiv3->id = 3;
    indiv4->id = 4;

    indiv1->idx_serodiscordant = -1;
    indiv2->idx_serodiscordant = -1;
    indiv3->idx_serodiscordant = -1;
    indiv4->idx_serodiscordant = -1;

    indiv1->gender = 1;
    indiv2->gender = 1;
    indiv3->gender = 0;
    indiv4->gender = 0;

    indiv1->HIV_status = 2;
    indiv2->HIV_status = 0;
    indiv3->HIV_status = 0;
    indiv4->HIV_status = 0;

    indiv1->n_partners = 0;
    indiv2->n_partners = 0;
    indiv3->n_partners = 0;
    indiv4->n_partners = 0;

    indiv1->n_HIVpos_partners = 0;
    indiv2->n_HIVpos_partners = 0;
    indiv3->n_HIVpos_partners = 0;
    indiv4->n_HIVpos_partners = 0;

    /* forming partnership between them and checking things are OK */
    new_partnership(indiv1, indiv3, param->start_time_simul, 
        overall_partnerships, param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;
    new_partnership(indiv1, indiv4, param->start_time_simul, 
        overall_partnerships, param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;
    new_partnership(indiv2, indiv4, param->start_time_simul, 
        overall_partnerships, param, debug, file_data_store);
    (*overall_partnerships->n_partnerships) ++;

    printf("- Before partnership dissolution:\n");

    print_partners(indiv1);
    print_partners(indiv2);
    print_partners(indiv3);
    print_partners(indiv4);

    printf("-----\n");
    print_individual(indiv1);
    printf("-----\n");
    print_individual(indiv2);
    printf("-----\n");
    print_individual(indiv3);
    printf("-----\n");
    print_individual(indiv4);
    printf("-----\n");

    printf("List of susceptibles in a serodiscordant partnerships:\n");
    for(i = 0; i < overall_partnerships->n_susceptible_in_serodiscordant_partnership[0] ; i++){
        print_individual(overall_partnerships->susceptible_in_serodiscordant_partnership[i]);
    }

    /* dissolving partnership between them and checking things are OK */
    breakup(param->start_time_simul, 
        &overall_partnerships->partner_pairs[initial_n_partnerships], overall_partnerships);

    breakup(param->start_time_simul, 
        &overall_partnerships->partner_pairs[*overall_partnerships->n_partnerships-1], 
        overall_partnerships);

    printf("- After partnership dissolution:\n");

    print_partners(indiv1);
    print_partners(indiv2);
    print_partners(indiv3);
    print_partners(indiv4);

    printf("-----\n");
    print_individual(indiv1);
    printf("-----\n");
    print_individual(indiv2);
    printf("-----\n");
    print_individual(indiv3);
    printf("-----\n");
    print_individual(indiv4);
    printf("-----\n");

    printf("List of susceptibles in a serodiscordant partnerships:\n");
    for(i = 0; i < overall_partnerships->n_susceptible_in_serodiscordant_partnership[0] ; i++){
        print_individual(overall_partnerships->susceptible_in_serodiscordant_partnership[i]);
    }

    free(indiv1);
    free(indiv2);
    free(indiv3);
    free(indiv4);
}


// Function: make_fake_population()
// The combination of the two functions below (to be executed together in this order)
// * first create an arbitrary population_size object with a certain distribution of the population
// * and then calculates and prints the number of partnerships to be drawn between each gender/age/risk groups in one time step given this current population distribution
// * (This allows checking that partnerships are drawn preferentially with similar age/risk groups)
// make_fake_population(n_population);
// check_draw_number_partnership(param);
//
// Parameters
// ----------
// n_pop: a pointer to a structure of type population_size, to be filled in by this function
//
// stratified_population_size *n_pop_strat
//

void make_fake_population(population_size *n_pop, stratified_population_size *n_pop_strat){

    int g, ag, r;

    for(g = 0; g < N_GENDER; g++){
        for(ag = 0; ag < N_AGE; ag++){
            for(r = 0; r < N_RISK; r++){
                n_pop->pop_size_per_gender_age_risk[g][ag][r] = gsl_ran_poisson (rng, 1200.0);
            }
        }
    }
    calcul_population(n_pop, n_pop_strat);
}


// Function : check_draw_number_partnership()
//
// Parameters
// ----------
//
// patch_struct : pointer to a patch structure
// p : int
//  Patch number
//
// Returns
// -------

void check_draw_number_partnership(patch_struct *patch, int p){

    /* filling in a population at random */
    make_fake_population(patch[p].n_population, patch[p].n_population_stratified);

    print_stratified_population(patch[p].n_population_stratified);
    print_population(patch[p].n_population);

    /* drawing number of partnerships between subgroups of this population */
    int ag_f, r_f, ag_m, r_m;

    draw_nb_new_partnerships(patch, patch[p].param,0,0);

    for(ag_f = 0; ag_f < N_AGE; ag_f++){
        for(r_f = 0; r_f < N_RISK; r_f++){
            for(ag_m = 0; ag_m < N_AGE; ag_m++){
                for(r_m = 0; r_m < N_RISK; r_m++){
                    printf("nb_partnerships_f_to_m[a_f=%d][r_f=%d][a_m=%d][r_m=%d]: %ld\n",
                        ag_f, r_f, ag_m, r_m, 
                        patch[p].param->balanced_nb_f_to_m[ag_f][r_f][ag_m][r_m]);
                    fflush(stdout);
                }
            }
        }
    }
}


// Function: check_available_partnerships()
// for debugging loops over all free partners but does not do super interesting stuff
//
//
// Parameters
// ----------
//
//
//
// Returns
// -------
// Nothing;

void check_available_partnerships(population_partners *pop_available_partners, 
    population_size_all_patches *n_pop_available_partners){
    
    int g, ag, r;
    long i, N;
    
    for(g = 0; g < N_GENDER; g++){
        for(ag = 0; ag < N_AGE; ag++){
            for(r = 0; r < N_RISK; r++){
                
                /* Loop through for each individual in this gender/age/risk category */
                N=n_pop_available_partners->pop_per_patch[0].pop_size_per_gender_age_risk[g][ag][r];
                
                for(i = 0; i < N; i++){
                    printf("g = %d \t a = %d \t r = %d \t i = %ld out of %ld\n", g, ag, r, i, N);
                    fflush(stdout);
                }
            }
        }
    }
}


// Function: check_males_females()
// Print the total population size, count and print number of males and females in the popn.
//
// Parameters
// ----------
// n_population_stratified: pointer to a stratified_population_size structure
// individual_population : pointer to an individual structure (or array)
//
// Returns
// -------
// Nothing; simply performs a check and prints output.

void check_males_females(stratified_population_size *n_population_stratified, 
    individual *individual_population){
    
    long n, n_m = 0, n_f = 0;
    
    printf("Total population size: %li\n", n_population_stratified->total_pop_size);
    
    for(n = 0; n < (n_population_stratified->total_pop_size); n++){
        if(individual_population[n].gender == MALE) n_m++; else n_f++;
    }
    printf("Total number of men: %li ; Total number of women: %li\n", n_m, n_f);
}


// Function: print_dob()
// Print the date of birth of every individual in the population
//
// Parameters
// ----------
// n_population_stratified : pointer to a stratified_population_size structure
// individual_population : pointer to an individual structure
//
// Returns
// -------
// Nothing; simply prints output of date-of-births
//
// Example
// -------
// Can be used within main()
// print_dob(n_population,individual_population);

void print_dob(stratified_population_size *n_population_stratified, 
    individual *individual_population){
    
    long i;
    
    for(i = 0; i < (n_population_stratified->total_pop_size); i++){
        printf("%f ", individual_population[i].DoB);
    }
}


// Function: validate_ages_based_on_age_group()
// Print out the age of everyone in age "age" at time t (calculated via DoB).
//
// This function is a validation test to ensure age groups haven't been miscalculated.
// Not currently used but could be used within the main() function.
//
//
// Parameters
// ----------
// age_list_struct : pointer to age_list structure
// age : int
// t : double
// 
// Returns
// -------
//
// Example
// -------

void validate_ages_based_on_age_group(age_list_struct *age_list, int age, double t){
    
    int aa, g;
    long i;
    double age_derived;

    for(g = 0; g < N_GENDER; g++){
        if(age < MAX_AGE){
            /* We need to adjust aa so that it is the correct index for age_list->age_group[]. */
            aa = (age - AGE_ADULT) + age_list->age_list_by_gender[g]->youngest_age_group_index;
            if(aa > (MAX_AGE - AGE_ADULT - 1)){
                aa = aa - (MAX_AGE - AGE_ADULT);
            }
            
            printf("Check: Number aged %i = %li\n", age, 
                age_list->age_list_by_gender[g]->number_per_age_group[aa]);
            
            for(i = 0; i < age_list->age_list_by_gender[g]->number_per_age_group[aa]; i++){
                age_derived = t - (age_list->age_list_by_gender[g]->age_group[aa][i])->DoB;
                printf("Derived age = %6.4f\n", age_derived);
            }
        }else{
            printf("Check: Number aged %i+ = %li\n", MAX_AGE, 
                age_list->age_list_by_gender[g]->number_oldest_age_group);
            
            for(i = 0; i < age_list->age_list_by_gender[g]->number_oldest_age_group; i++){
                age_derived = t - (age_list->age_list_by_gender[g]->oldest_age_group[i])->DoB;
                printf("Derived age = %6.4f\n", age_derived);
            }
        }
    }
}
