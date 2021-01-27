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

#ifndef INTERVENTIONS_H_
#define INTERVENTIONS_H_

#include "structures.h"
/* PopART CHiPs-specific processes */
void create_popart_chips_samples(age_list_struct *, chips_sample_struct *, parameters *, int, int);
void schedule_chips_visits(chips_sample_struct *, parameters *, int);
void carry_out_chips_visits_per_timestep(int, int , patch_struct *, int , int , debug_struct *, output_struct *);
void chips_visit_person(individual *, cumulative_outputs_struct *, calendar_outputs_struct *, 
    double ,individual ***, long *, long *, individual ***, long *, long *, parameters *, 
    individual ***, long *, long *, patch_struct *, int , int , debug_struct *, output_struct *,
    int, int);

/* VMMC intervention (both popart and non-popart). */
void draw_if_VMMC(individual *, parameters *, individual ***, long *, long *, double , int );
void schedule_vmmc(individual *,  parameters *, individual ***, long *, long *, double , int );
void schedule_vmmc_healing(individual *, parameters *, individual ***, long *, long *, double );
void schedule_generic_vmmc_event(individual *, parameters *, individual ***, long *, long *, double, double );
void carry_out_VMMC_events_per_timestep(int , double , patch_struct *, int );

#endif /* INTERVENTIONS_H_ */
