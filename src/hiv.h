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


#ifndef HIV_H_
#define HIV_H_


#include "constants.h"
#include "structures.h"
#include "utilities.h"
#include "output.h"
#include "debug.h"


int get_spvl_cat(double );
void draw_inital_SPVL(individual *, parameters *);
void inherit_spvl(individual *, individual *, parameters *);
double get_mean_time_hiv_progression(parameters *, individual *);
double get_RR_SPVL(double , parameters *);
double hiv_transmission_probability(individual *, parameters *);
void hiv_acquisition(individual*, double , patch_struct *, int, all_partnerships *, 
    output_struct *, debug_struct *, file_struct *, int, int);
int find_who_infected(int, double);
void inform_partners_of_seroconversion_and_update_list_serodiscordant_partnerships(individual *,
    individual **, long *);
void new_infection(double, int, individual *, individual *, population_size_one_year_age *,
    population_size_one_year_age *, age_list_struct *,  parameters *, individual ***, long *,
    long *, population_size_one_year_age *, file_struct *);
void draw_initial_infection(double, individual* , patch_struct *, int, all_partnerships *,
    output_struct *, file_struct *, int, int);
void next_hiv_event(individual *, individual ***, long *, long *, parameters *, double ,
    cumulative_outputs_struct *, calendar_outputs_struct *);
void carry_out_HIV_events_per_timestep(double , patch_struct *, int , all_partnerships *,
    debug_struct *, file_struct *);
int get_window_result(double ,double , patch_struct *);
int joins_preart_care(individual* , parameters *, double , cumulative_outputs_struct *,
    calendar_outputs_struct *);
int remains_in_cascade(individual* , parameters *, int);
int measured_cd4_cat(parameters *, int );
int art_cd4_eligibility_group(parameters *, double);
int is_eligible_for_art(individual* , parameters *, double , patch_struct *, int );
double get_time_emergency_start_ART(individual *, parameters *, double );
void start_ART_process(individual* , parameters *, double , individual ***, long *, long *,
    individual ***, long *, long *, int , file_struct *, calendar_outputs_struct *);
void draw_initial_hiv_tests(parameters *, age_list_struct *, double, individual ***, 
    long *, long *);
void draw_hiv_tests(parameters *, age_list_struct *, int , individual ***, long *, long *, int );
void schedule_generic_cascade_event(individual* , parameters *, double , individual ***, 
    long *, long *, double);
void schedule_new_hiv_test(individual *, parameters *, double, individual ***, long *, long *);
void probability_get_hiv_test_in_next_window(double *, double *, int , int , int , parameters *);
void schedule_hiv_test_fixedtime(individual* , parameters *, int , individual ***, long *, 
    long *, double , int , double *);
void hiv_test_process(individual* , parameters *, double, individual ***, long *, long *, 
    individual ***, long *, long *, cumulative_outputs_struct *, calendar_outputs_struct *, 
    individual ***, long *, long *, patch_struct *, int , debug_struct *);
void schedule_start_of_art(individual* , parameters *, double , individual ***, long *, long *);
void cd4_test_process(individual* , parameters *, double , individual ***, long *, long *,
    individual ***, long *, long *, cumulative_outputs_struct *, calendar_outputs_struct *,
    patch_struct *, int );
void virally_suppressed_process(individual* , parameters *, double , individual ***, long *, 
    long *, individual ***, long *, long *);
void virally_unsuppressed_process(individual* , parameters *, double , individual ***, long *,
    long *, individual ***, long *, long *, cumulative_outputs_struct *, calendar_outputs_struct *);
void dropout_process(individual* , parameters *, double , individual ***, long *, long *,
    individual ***, long *, long *, cumulative_outputs_struct *, calendar_outputs_struct *);
void carry_out_cascade_events_per_timestep(double , patch_struct *, int , all_partnerships *,
    debug_struct *, file_struct *);
int is_cabo_target_group(individual* , double t);
double PANGEA_get_cd4(individual* , double );
#endif /* HIV_H_ */
