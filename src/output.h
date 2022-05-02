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

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "structures.h"
#include "demographics.h"
#include "hiv.h"

void print_individual(individual *);
void print_population(population_size *);
void print_population_from_one_year_data(population_size *, population_size_one_year_age *);
void print_population_one_year(population_size_one_year_age *);
void print_population_one_year_only_old(population_size_one_year_age *);
void print_stratified_population(stratified_population_size *);
void print_partners(individual *);
void print_partnership(partnership *);
void print_HIV_status(individual *);
void print_demographics(individual *, population_size *, stratified_population_size *, double );
void print_number_by_age(age_list_struct *);
void print_IDs_by_age(age_list_struct *);
void print_specific_IDs_by_age(long , age_list_struct *, int);
void print_number_by_age_grouped(age_list_struct *,population_size *, stratified_population_size*);
void update_outputs_gender_veryshort(individual *, long *, long *, long *, long *, long *, long *);
//void update_outputs_gender_veryshort(individual *, long *, long *, long *, long *, long *);
void update_annual_outputs_gender(individual *, long *, long *, long *, long *, int );
void update_annual_outputs_gender_cd4(individual *, long *, long *, long *, long *, 
    long [NCD4], long [NCD4], long [NCD4], int );
void update_calibration_outputs_gender(individual *, long *, long *, long *, long *);
void update_annual_outputs_riskgp(individual *, long *, long *, long *, long *, long *, long *, 
    long *, long *, long *, long *);
void update_annual_outputs_npartnersbygenderagerisk(individual *, long [N_GENDER][N_RISK][N_AGE+1],
    long [N_GENDER][N_RISK][N_AGE+1], long [N_GENDER][N_RISK][N_AGE+1], 
    long [N_GENDER][N_RISK][N_AGE+1], long [N_GENDER][N_RISK][N_AGE+1], 
    long [N_GENDER][N_RISK][N_AGE+1], long [N_GENDER][N_RISK][N_AGE+1], int , double );
void store_annual_outputs(patch_struct *, int , output_struct *, all_partnerships *, 
    int *, int , int );
void store_annual_partnerships_outputs(patch_struct *, int , output_struct *, 
    all_partnerships *, int *, int , int );
void store_timestep_outputs(patch_struct *, int , double , output_struct *, int, int, int);
void store_timestep_age_outputs(patch_struct *, int , double , output_struct *, int);
void store_calibration_outputs_dhs(patch_struct *, int , output_struct *, int );
void store_cost_effectiveness_outputs(patch_struct *, int , output_struct *, all_partnerships *,
    int *, int, int);
void store_treats_outputs(patch_struct *, int, output_struct *, all_partnerships *, 
        int *, int, int);
void store_art_status_by_age_sex(patch_struct *, int , double , output_struct *);
void save_calibration_outputs_pc(patch_struct *, int , output_struct *, int, int );
void save_person_timesteps_pc(patch_struct *, int, output_struct *, int , int );
void store_calibration_outputs_pc(patch_struct *, int , output_struct *);
void store_calibration_outputs_chips(patch_struct *, int ,  output_struct *);
void write_annual_outputs(file_struct *, output_struct *, int);
void write_annual_partnerships_outputs(file_struct *, output_struct *, int);
void write_timestep_outputs(file_struct *, output_struct *, int, int);
void write_timestep_age_outputs(file_struct *, output_struct *, int, int);
void blank_calibration_output_file(char *, int );
void write_calibration_outputs(char *, output_struct *, int);
void store_phylogenetic_transmission_output(output_struct *, double, individual*, 
    individual*, file_struct *, int, int);
void store_phylogenetic_transmission_initial_cases(output_struct *, parameters *, 
    individual*, file_struct *, int, int);
void blank_phylo_transmission_data_file(file_struct *);
void write_phylo_transmission_data(file_struct *, char *, int);
void write_phylo_individual_data(file_struct *, individual *,  long, int);
void write_hivpos_individual_data(file_struct *, individual *, long, int);
void print_partnership_network(file_struct *, char *, file_label_struct *, patch_struct *, int , int );
void print_partners_outside_community(char *, individual *, long , int , int );
void initialise_partners_outside_community_file(char *, int );
void print_assortativity(char *, debug_struct *, patch_struct *, int, file_struct *);
void sweep_through_all_and_compute_distribution_lifetime_and_lastyear_partners(patch_struct *, 
    all_partnerships * , int , int , output_struct *);
void write_distr_n_lifetime_partners_and_n_partners_lastyear(patch_struct *, file_struct *);
void write_pc_data(patch_struct *, int , file_struct *);
void write_chips_data_visit(patch_struct *, int , file_struct *, output_struct *);
void write_chips_data_annual(patch_struct *, int , int , int , int , file_struct *);
void write_cost_effectiveness_outputs(file_struct *, output_struct *, int);
void write_treats_outputs(file_struct *, output_struct *, int);
void write_art_status_by_age_sex(file_struct *, output_struct *, int);
#endif /* OUTPUT_H_ */
