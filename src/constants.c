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

#include "constants.h"

const int AGE_GROUPS[N_AGE] = {13,18,23,30,40,50,60};

// The 80 and over age group is only here for the ageing process. 
// We don't differentiate 61-79 and 80 and over for other processes.  
const int AGE_GROUPS_WITH_OLD[N_AGE+1] = {13,18,23,30,40,50,60,80};

const int AGE_GROUPS_UNPD[N_AGE_UNPD+1] = {13,15,20,25,30,35,40,45,50,55,60,65,70,75,80};

const int FIND_AGE_GROUPS_UNPD[MAX_AGE - AGE_ADULT + 1] = {
    0,0, // 13, 14
    1,1,1,1,1, // 15, 16, 17, 18, 19
    2,2,2,2,2, // 20, 21, 22, 23, 24
    3,3,3,3,3, // 25, 26, 27, 28, 29
    4,4,4,4,4, // 30, 31, 32, 33, 34
    5,5,5,5,5, // 35, 36, 37, 38, 39
    6,6,6,6,6, // 40, 41, 42, 43, 44
    7,7,7,7,7, // 45, 46, 47, 48, 49
    8,8,8,8,8, // 50, 51, 52, 53, 54
    9,9,9,9,9, // 55, 56, 57, 58, 59
    10,10,10,10,10, // 60, 61, 62, 63, 64
    11,11,11,11,11, // 65, 66, 67, 68, 69
    12,12,12,12,12, // 70, 71, 72, 73, 74
    13,13,13,13,13, // 75, 76, 77, 78, 79
    14 // 80 +
};


const int FIND_AGE_GROUPS[MAX_AGE-AGE_ADULT+1] = {
    0,0,0,0,0,
    1,1,1,1,1,
    2,2,2,2,2,2,2,
    3,3,3,3,3,3,3,3,3,3,
    4,4,4,4,4,4,4,4,4,4,
    5,5,5,5,5,5,5,5,5,5,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6
};

// Define names of the risk groups. 
// The 5 is the max length of each name +1 (ie length("High")+1)
const char RISK_GP_NAMES[N_RISK][5] = {"Low","Med","High"};

int POPART_SAMPLING_FRAME_ESTABLISHED;
