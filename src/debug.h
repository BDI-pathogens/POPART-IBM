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

#ifndef DEBUG_H_
#define DEBUG_H_

#include "constants.h"
#include "structures.h"
#include "debug.h"
#include "demographics.h"

void find_in_age_list(double ,individual* , age_list_struct *, parameters *);
void print_age_list(age_list_struct *);
void count_population_by_going_through_indiv(patch_struct *, long *, long *);
void count_population_by_going_through_age_list(patch_struct *, long *, long *);
void count_population_using_n_population(population_size *, long *, long *);
void count_population_size_three_ways(patch_struct *, int , double );
int is_in_patch_individual_population(long, patch_struct *, int);

void generate_demographics_byage_gender_file_header(char *, int);
void write_demographics_byage_gender(patch_struct *, int , double , file_struct *);
void blank_one_year_age_groups_including_kids(file_struct *);
void write_one_year_age_groups_including_kids(file_struct *, patch_struct *, int , double );
void write_nbirths_nnewadults_ndeaths(file_struct *, patch_struct *, int );
void output_life_expectancy(char *, patch_struct *, int , int );

void check_if_individual_should_be_in_list_susceptibles_in_serodiscordant_partnership(individual *,
    all_partnerships * );
void check_if_individual_should_be_in_list_available_partners(individual *, all_partnerships *,
    int , int );
void sweep_through_all_and_check_lists_serodiscordant_and_available_partners (patch_struct *,
    all_partnerships * , int , int );
void sweep_through_all_and_check_n_partners_outside_n_HIVpos_partners_and_n_HIVpos_partners_outside (patch_struct *, all_partnerships * , int , int );
void sweep_through_all_and_check_age_and_risk_of_partners (patch_struct *, all_partnerships *,
    int , int , debug_struct *);

void blank_debugging_files(file_struct *);
void write_hiv_duration(individual *, double , file_struct *);
void write_hiv_duration_km(individual *, double , file_struct *, int );
void write_hiv_duration_km_end_of_simulation(patch_struct *, double , file_struct *);
void write_cd4_at_seroconversion(individual *, file_struct *);
void write_initial_spvl_distribution(individual *, file_struct *);
void write_cd4_spvl_states(patch_struct *, int , file_struct *);
void print_debugerror_shouldbezero_exit(char *);
void write_art_states(patch_struct *, int , debug_struct *, file_struct *);
void reset_annual_chips_visit_counter(age_list_struct *);
void print_chips_statistics_using_age_list(age_list_struct *, double );
void print_chips_statistics_using_chipsonly(patch_struct *, int , double );

void output_hazard_over_time_period(double ,double , individual *, individual *, 
    file_struct *, output_struct *);
void write_hazard_data(file_struct *, char *);
void blank_hazard_file(file_struct *);

#endif /* DEBUG_H_ */
