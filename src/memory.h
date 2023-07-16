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

#ifndef MEMORY_H_
#define MEMORY_H_

#include "structures.h"

void blank_individual_array(individual *, int);
void reinitialize_arrays_to_default(int ,patch_struct  *, all_partnerships *, output_struct *);
void alloc_pop_memory(population **, int);
void alloc_output_memory(output_struct **);
void set_to_null_pc_cohort_data(patch_struct *, int , int , int );
void alloc_pc_cohort_data(PC_cohort_data_struct **, int , int );
void alloc_patch_memoryv2(patch_struct *);
void alloc_partnership_memoryv2(all_partnerships *);
void free_all_partnership_memory(partnership *, long *, individual **, long *,
        population_partners* , population_size_all_patches *,
        partnership ***, long *, long *,
        long *, long *, long *, long *, long *, long *);
void free_all_patch_memory(parameters *, individual *, population_size *, 
        population_size_one_year_age *, stratified_population_size *, age_list_struct *, 
        child_population_struct *, individual ***, long *, long *, individual ***, 
        long *, long *, individual ***, long *, long *, long *, long *,
        population_size_one_year_age *, population_size_one_year_age *, 
        population_size_one_year_age *,population_size_one_year_age *,
        population_size_one_year_age *, population_size *, population_size *,
        chips_sample_struct *, cumulative_outputs_struct *, calendar_outputs_struct *, 
        long ****, long ****, PC_sample_struct *, PC_cohort_struct *, PC_cohort_data_struct *);
void free_partnership_memory(all_partnerships *);
void free_patch_memory(patch_struct *);
void free_pop_memory(population *, parameters **);
void allocate_fitting_data_memory(int , fitting_data_struct **);
void free_fitting_data_memory(fitting_data_struct *);
void free_output_memory(output_struct *);
void free_pc_cohort_data_memory(PC_cohort_data_struct *);

#endif /** MEMORY_H_ **/
