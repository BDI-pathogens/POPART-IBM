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


#ifndef PARTNERSHIP_H_
#define PARTNERSHIP_H_

#include "structures.h"

void new_partnership(individual* , individual* , int, double , all_partnerships *, parameters *, 
    debug_struct *, file_struct *);
int time_to_partnership_dissolution(parameters *, int r_m, int r_f, int p_m, int p_f);
void breakup(double, partnership*, all_partnerships *);
void update_list_available_partners_breakup(double , partnership* , population_partners*, 
    population_size_all_patches *);
void add_susceptible_to_list_serodiscordant_partnership(individual* , individual** , long *);
void remove_susceptible_from_list_serodiscordant_partnership(individual* , individual** , long *);
void update_list_susceptibles_in_serodiscordant_partnerships_breakup(partnership* , 
    individual** , long *);
void draw_nb_new_partnerships(patch_struct *, parameters *, int, int);
void draw_n_new_partnerships(double , long, long, parameters *, int , int , int , int , int *,
        all_partnerships *, patch_struct *, int , int , debug_struct *, file_struct *);
void draw_new_partnerships(double , all_partnerships *, patch_struct *, parameters *, int , int , 
    debug_struct *, file_struct *);

#endif /* PARTNERSHIP_H_ */
