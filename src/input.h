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

#ifndef INPUT_H_
#define INPUT_H_

#include "structures.h"


void read_param(char *, parameters **, int, patch_struct *);
void read_patch_info(char *, patch_struct *);
void read_demographic_params(char *, parameters *, int);
void read_hiv_params(char *, parameters *, int, int);
void read_partnership_params(char *, parameters *, int);
void read_time_params(char *, parameters *, int, int);
void read_cascade_params(char *, parameters *, int);
void read_chips_uptake_params(char *, parameters *);
void copy_chips_params(parameters **, int );
void read_pc0_enrolment_params(char *, int , parameters *, int , int );
void copy_pc_params( parameters **, int );
void read_pc_future_params(char *, parameters *, int );
void read_initial_params(char *, parameters *, int);
long get_python_seed(char *);
//void get_uptake_scenarios(char *, uptake_scenario_struct *, patch_struct *, int );
/* void change_popart_y3_params(parameters *, uptake_scenario_struct *); */
/* void force_y3_params(parameters *, uptake_scenario_struct *); */
#endif /* INPUT_H_ */
