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

#ifndef FITTING_H_
#define FITTING_H_

int load_fitting_data_n(char *);
void load_fitting_data(char *, fitting_data_struct *, int);
void check_fitting_data(fitting_data_struct *, int );
int fit_data(int , int , fitting_data_struct *, patch_struct *, int );
double perform_target_fit(fitting_data_struct *, double );

#endif /* FITTING_H_ */
