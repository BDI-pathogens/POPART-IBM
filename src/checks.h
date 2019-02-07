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

#ifndef CHECKS_H_
#define CHECKS_H_

#include "structures.h"

void check_partnership_formation(all_partnerships *, parameters *, debug_struct *, file_struct *);
void check_partnership_formation_and_HIV_acquisition(patch_struct *, int , all_partnerships *,
    output_struct *, debug_struct *, file_struct *);
void check_partnership_dissolution(all_partnerships *, parameters *, debug_struct *, file_struct *);
void make_fake_population(population_size *, stratified_population_size *);
void check_draw_number_partnership(patch_struct *, int);
void check_available_partnerships(population_partners *, population_size_all_patches *);
void check_males_females(stratified_population_size *,individual *);
void print_dob(stratified_population_size *,individual *);
void validate_ages_based_on_age_group(age_list_struct *, int , double );

#endif /* CHECKS_H_ */
