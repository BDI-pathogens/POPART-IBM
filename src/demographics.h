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

#ifndef DEMOGRAPHICS_H_
#define DEMOGRAPHICS_H_

#include "structures.h"

double per_woman_fertility_rate(int , parameters *, int, double);
void get_unpd_time_indices(double , int *, double *);
double childhood_mortality(parameters *, double );
double natural_death_rate(int , int , parameters *, double );
int draw_sex_risk(int, int, parameters *);
void create_new_individual(individual *, double , parameters *, int, 
    population_size_one_year_age *, patch_struct *, int, all_partnerships *);
void update_population_size_new_adult(individual *, population_size *, 
    population_size_one_year_age *, stratified_population_size *);
void update_population_size_death(individual *, population_size *, 
    population_size_one_year_age *, population_size_one_year_age *,
    population_size_one_year_age *,population_size_one_year_age *,
    population_size_one_year_age *, stratified_population_size *,
    int, age_list_struct *);
void initialize_first_cascade_event_for_new_individual(individual *, double, 
    parameters *, individual ***, long *, long *);
void update_age_list_new_adult(age_list_struct *, individual *);
void update_age_list_death(age_list_struct *, int, int, long, double , int);
int get_age_index(double , double );
int get_age_indexv2(double , double , int);
int get_age_group(double , double , const int [], int);
int get_age_group_unpd(double , double );
void update_n_population_ageing_by_one_year(patch_struct *patch, int p);
void age_population_by_one_year(age_list_struct *);
void update_pop_available_partners_ageing_by_one_year(patch_struct *, int, 
    all_partnerships *, double );
void age_population_size_one_year_age_by_one_year(population_size_one_year_age *);
void remove_dead_person_from_susceptible_in_serodiscordant_partnership(individual *, 
    individual **, long *);
void remove_dead_person_from_list_available_partners(double, individual *,population_partners *,
    population_size_all_patches *);
void remove_dead_persons_partners(individual *, population_partners *, 
    population_size_all_patches *, double );
void remove_from_hiv_pos_progression(individual *, individual ***, long *, long *, 
    double, parameters *, int);
void remove_from_cascade_events(individual *, individual ***, long *, long *, double, parameters *);
void remove_from_vmmc_events(individual *, individual ***, long *, long *, double , parameters *);
void increase_population_count_art_or_virallysuppressed_or_cabo(population_size_one_year_age *, individual *, double);
void decrease_population_count_art_or_virallysuppressed_or_cabo(population_size_one_year_age *, individual *, double);
void deaths_natural_causes(double, patch_struct *, int , all_partnerships *,  file_struct *);
void make_new_adults(double, patch_struct *, int , all_partnerships *);
void add_new_kids(double , patch_struct *, int );
void make_pop_from_age_list(population *, age_list_struct *, individual *);
void individual_death_AIDS(age_list_struct *, individual *, population_size *, 
    population_size_one_year_age *, population_size_one_year_age *, 
    population_size_one_year_age *, population_size_one_year_age *, population_size_one_year_age *,
    stratified_population_size *, double , parameters *, individual **, 
    long *, population_partners *, population_size_all_patches *, individual ***, 
    long *, long *, patch_struct *, int , file_struct *);

#endif /* DEMOGRAPHICS_H_ */