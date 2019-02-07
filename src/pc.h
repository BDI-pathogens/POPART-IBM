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

#ifndef PC_H_
#define PC_H_


#include "constants.h"
#include "structures.h"
#include "utilities.h"
#include "output.h"
#include "debug.h"

int get_PC_HIV_stratum(individual *);

void remove_extras_from_timestep_recruitment(patch_struct *, int , int , int , int , 
    int , int , int );
void create_popart_pc_sample(patch_struct *, age_list_struct *, PC_sample_struct *, 
    parameters *, int , int );
void carry_out_PC_enrolment_per_timestep(int , int , patch_struct *, int , int );
void PC_enroll_person(individual *, patch_struct *, double , int , int , int , int , int );
void PC_next_cohort_round(patch_struct *, int , int );
void carry_out_PC_visits_per_timestep(int , int , patch_struct *, int , int , int );
void PC_visit_person(individual *, patch_struct *,  double , int , int , int );
void reset_visit_counter(patch_struct *, int );
int get_pc_round(int ,int , patch_struct *, int );
#endif /* PC_H_ */
