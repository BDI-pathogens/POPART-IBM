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

#ifndef INIT_H_
#define INIT_H_



#include "demographics.h"
#include "constants.h"
#include "structures.h"
#include "demographics.h"

void get_initial_population_distribution(population_size *, parameters *);
int set_max_n_partners(int , int, int, parameters *);
double make_DoB(int , double, int *);
void set_population_count_zero(population_size*);
void set_population_count_one_year_zero(population_size_one_year_age *n_population);
void set_population_count_stratified(stratified_population_size*, population_size*);
void set_population_count_stratified_zero(stratified_population_size*);
void initialize_child_population(parameters *, child_population_struct *, stratified_population_size *, int, age_list_struct *);
void set_up_population(int, patch_struct *, population *);
void init_available_partnerships(int , patch_struct *, all_partnerships *,population *);
void init_cumulative_counters(cumulative_outputs_struct *);
void init_calendar_counters(calendar_outputs_struct *);
void initialise_debug_variables(debug_struct *);

#endif /* INIT_H_ */
