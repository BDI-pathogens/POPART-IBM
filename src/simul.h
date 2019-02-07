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


#ifndef SIMUL_H_
#define SIMUL_H_

#include "structures.h"

void carry_out_partnership_processes_by_time_step(int , int , patch_struct *, all_partnerships *, 
    output_struct *, debug_struct *, file_struct *);
int carry_out_processes(int, fitting_data_struct *, patch_struct *, all_partnerships *, 
    output_struct *, int, int, debug_struct *, file_struct *, int);
int carry_out_processes_by_patch_by_time_step(int , int , fitting_data_struct *, patch_struct *, 
    int , all_partnerships *, output_struct *, int, int, debug_struct *, file_struct *, int);

#endif /* SIMUL_H_ */
