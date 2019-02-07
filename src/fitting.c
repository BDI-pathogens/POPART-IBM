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

/* Contains potential fitting routines for PopART IBM.
 */

#include "structures.h"
#include "constants.h"
#include "fitting.h"
#include "memory.h"
#include "utilities.h"

/* Functions in this file:
    load_fitting_data()
        loads in fitting data from a file called "fitting_data_processed.txt" in the parameter 
        directory.
    check_fitting_data()
        checks that the fitting data lies within some fairly loose ranges.
 */

int load_fitting_data_n(char *file_directory){
   /* Load n time points of fitting data, allocate memory for fitting_data[n], store fitting data.
    
    Parameters
    ----------
    file_directory : pointer to a char array
        
    
    Returns
    -------
    n_fit : int
        Number of time points of fitting data.
    */
    
    /* Number of time points of data we fit to. */
    int n_fit; 
    FILE *fitting_file;
    char fitting_file_name[LONGSTRINGLENGTH];
    int checkreadok;

    // Add path before file name
    strncpy(fitting_file_name,file_directory,LONGSTRINGLENGTH);
    add_slash(fitting_file_name);
    join_strings_with_check(fitting_file_name, "fitting_data_processed.txt", LONGSTRINGLENGTH, 
        "'fitting_data_processed.txt' and fitting_file_name in load_fitting_data_n()");
    
    // Open parameter file
    if ((fitting_file = fopen(fitting_file_name, "r")) == NULL){
        printf("Cannot open fitting_data_processed.txt");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Fitting data read from: %s:\n",fitting_file_name);
        }
    }
    
    checkreadok = fscanf(fitting_file, "%i", &n_fit);
    check_if_cannot_read_param(checkreadok, "number of data points to fit to");
    fclose(fitting_file);
    return(n_fit);
}


void load_fitting_data(char *file_directory, fitting_data_struct *fitting_data, int n_fit){
    /*
    
    
    Parameters
    ----------
    file_directory : pointer to a char array
    fitting_data : pointer to a fitting_data_struct structure
    n_fit : int
        Number of fitted data points
    
    Returns
    -------
    Nothing
    
    */
    
    
    int i;
    FILE *fitting_file;
    char fitting_file_name[LONGSTRINGLENGTH];
     /* We reopen the file so need somewhere to stash the n_fit value for now. */
    int dummy_n_fit;
    int checkreadok;

    // Add path before file name
    strncpy(fitting_file_name, file_directory, LONGSTRINGLENGTH);
    add_slash(fitting_file_name); 
    join_strings_with_check(fitting_file_name, "fitting_data_processed.txt", LONGSTRINGLENGTH, 
        "'fitting_data_processed.txt' and fitting_file_name in load_fitting_data()");
    
    // Open parameter file
    if ((fitting_file = fopen(fitting_file_name, "r")) == NULL){
        printf("Cannot open fitting_data_processed.txt");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }else{
        if(VERBOSE_OUTPUT == 1){
            printf("Fitting data read from: %s:\n", fitting_file_name);
        }
    }
    
    checkreadok = fscanf(fitting_file, "%i", &dummy_n_fit);
    printf("Fitting to %i points\n", n_fit);

    for(i = 0; i < n_fit; i++){

        checkreadok = fscanf(fitting_file, "%i", &(fitting_data[i].whatfitto));
        check_if_cannot_read_param(checkreadok, "type of data fitted to");

        checkreadok = fscanf(fitting_file, "%lg", &(fitting_data[i].fit_time));
        check_if_cannot_read_param(checkreadok, "fit time");

        checkreadok = fscanf(fitting_file, "%lg", &(fitting_data[i].prevalence_point_est));
        check_if_cannot_read_param(checkreadok, "prevalence point est");

        checkreadok = fscanf(fitting_file, "%lg", &(fitting_data[i].prevalence_ll));
        check_if_cannot_read_param(checkreadok, "prevalence ll");

        checkreadok = fscanf(fitting_file, "%lg", &(fitting_data[i].prevalence_ul));
        check_if_cannot_read_param(checkreadok, "prevalence ul");
    }
    fclose(fitting_file);

    /* Run checks on prevalence to make sure plausible: */
    check_fitting_data(fitting_data, n_fit);

    /* Now get the (discrete) time for each data point: */
    for(i = 0; i < n_fit; i++){
        fitting_data[i].fit_year = (int) fitting_data[i].fit_time;
        fitting_data[i].fit_timestep = (int) ((fitting_data[i].fit_time - 
                fitting_data[i].fit_year)*N_TIME_STEP_PER_YEAR);
        
        if(fitting_data[i].fit_timestep >= N_TIME_STEP_PER_YEAR){
            printf("fitting_data[i].fit_timestep=%i is too big - exiting\n",
                fitting_data[i].fit_timestep);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
}


void check_fitting_data(fitting_data_struct *fitting_data, int n_fit){
    
    /* 
    Parameters
    ----------
    
    fitting_data : pointer to a fitting_data_struct stucture
    n_fit : int
    
    */
    
    int i;
    for(i = 0; i < n_fit; i++){
        
        /* Fitting data should be between 1960 and 2020 - note that this may change over time. */
        if((fitting_data[i].fit_time) < 1960 || (fitting_data[i].fit_time > 2020)){
            printf("ERROR: Fit time %lg (i=%i) is outside expected range 1960-2020\n. Exiting\n",
                fitting_data[i].fit_time,i);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        /* Prevalence has to lie in range [0,1]: */
        if (fitting_data[i].prevalence_point_est < 0 || fitting_data[i].prevalence_point_est > 1){
            printf("ERROR: Prevalence %lg has to lie in range [0,1]. Exiting\n",
                fitting_data[i].prevalence_point_est);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        // Check lower-limit < point estimate < upper limit 
        if((fitting_data[i].prevalence_ll > fitting_data[i].prevalence_point_est) || 
            (fitting_data[i].prevalence_ul < fitting_data[i].prevalence_point_est)){
            
            printf("ERROR: Point estimate must lie in range [LL,UL]: %lg [%lg,%lg]\nExiting\n",
                fitting_data[i].prevalence_point_est, fitting_data[i].prevalence_ll, 
                fitting_data[i].prevalence_ul);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }
}


int fit_data(int t0, int t_step, fitting_data_struct *fitting_data, patch_struct *patch, int p){
    
   /* Function does the following:
    * 1. Checks if there are any fitting criteria left
    * 2. If there are checks if there are any at the current timestep (there can be more than 1 condition per timestep - e.g. PC data may be fitted by age x sex etc.
    * NB AT PRESENT THIS IS SET UP JUST TO WORK WITH OVERALL PREVALENCE - CAN BE EASILY TWEAKED LATER FOR OTHER GROUPS.
    * 3. Calculate overall HIV prevalence
    * 4. Pass this to the function perform_target_fit() to actually do the fitting - at the moment it sits in a target-fitting like routine but can be modified.
    * 5. Return values: If all (target-fitting) conditions are passed (if there are any) then return 1 else exit when the first non-fitted condition is encountered, and return 0.
    */
    
    
    double this_fit;
    int g, r, aa, ai_npopulation, ai_positive;
    int npositive, npop, minage, maxage;
    double thisprev;

    // Check if there are any fitting criteria left.  
    if (patch[p].i_fit < patch[p].n_fit){
        
        while(
            (t0 == fitting_data[patch[p].i_fit].fit_year) && 
            (t_step == fitting_data[patch[p].i_fit].fit_timestep)
        ){

            /* Fit to prevalence in 15-49 year olds- so need to calculate prevalence first: */
            if (fitting_data[patch[p].i_fit].whatfitto == 0){
                
                npositive = 0;
                npop = 0;
                /* Fit to 15-49 year olds: */
                for(g = 0; g < N_GENDER; g++){
                    for(aa = (15 - AGE_ADULT); aa < (49 - AGE_ADULT); aa++){
                        npop += patch[p].age_list->age_list_by_gender[g]->number_per_age_group[aa];
                    }
                }
                
                for(g = 0; g < N_GENDER; g++){
                    for(r = 0; r < N_RISK; r++){
                        for(aa = (15 - AGE_ADULT); aa < (49 - AGE_ADULT); aa++){
                            npositive += 
                                patch[p].n_infected->pop_size_per_gender_age1_risk[g][aa][r];
                        }
                    }
                }
                if (npop>0){
                    thisprev = npositive/(npop + 0.0);
                }else{
                    thisprev = 0.0;
                }
            }else if(
                (fitting_data[patch[p].i_fit].whatfitto >= 1) && 
                (fitting_data[patch[p].i_fit].whatfitto <= 14)
            ){
                /* fitting_data[i_fit].whatfitto indices that we fit to:
                1 - male 15-49  - ie full DHS agegroup
                2 - male 15-19
                3 - male 20-24
                4 - male 25-29
                5 - male 30-34
                6 - male 35-39
                7 - male 40-44
                8 - female 15-49  - ie full DHS agegroup
                9 - female 15-19
                10 - female 20-24
                11 - female 25-29
                12 - female 30-34
                13 - female 35-39
                14 - female 40-44 */
                
                npositive = 0;
                npop = 0;
                /* MALE 0, FEMALE 1 */
                /* 1-7 is male, 8-14 is female. Use integer division e.g. 6/7 = 0, 7/7=1... */
                g = (fitting_data[patch[p].i_fit].whatfitto - 1)/7;
                
                if(g < 0 || g > 1){
                    printf("Error - gender undefined in fitting routine. Exiting\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }

                if(
                    fitting_data[patch[p].i_fit].whatfitto == 1 || 
                    fitting_data[patch[p].i_fit].whatfitto == 8){
                    minage = 15;
                    maxage = 49;
                }else if(
                    fitting_data[patch[p].i_fit].whatfitto == 2 || 
                    fitting_data[patch[p].i_fit].whatfitto == 9){
                    minage = 15;
                    maxage = 19;
                }else if(
                    fitting_data[patch[p].i_fit].whatfitto == 3 || 
                    fitting_data[patch[p].i_fit].whatfitto == 10){
                    minage = 20;
                    maxage = 24;
                }else if(
                    fitting_data[patch[p].i_fit].whatfitto == 4 || 
                    fitting_data[patch[p].i_fit].whatfitto == 11){
                    minage = 25;
                    maxage = 29;
                }else if(
                    fitting_data[patch[p].i_fit].whatfitto == 5 || 
                    fitting_data[patch[p].i_fit].whatfitto == 12){
                    minage = 30;
                    maxage = 30;
                }
                else if(
                    fitting_data[patch[p].i_fit].whatfitto == 6 || 
                    fitting_data[patch[p].i_fit].whatfitto == 13){
                    minage = 35;
                    maxage = 39;
                }else if(
                    fitting_data[patch[p].i_fit].whatfitto == 7 || 
                    fitting_data[patch[p].i_fit].whatfitto == 14){
                    minage = 40;
                    maxage = 44;
                }
                
                for(r = 0; r < N_RISK; r++){
                    for(aa = (minage - AGE_ADULT); aa < (maxage - AGE_ADULT); aa++){
                        ai_npopulation = aa + 
                            patch[p].n_population_oneyearagegroups->youngest_age_group_index;
                        
                        while(ai_npopulation > (MAX_AGE - AGE_ADULT - 1)){
                            ai_npopulation = ai_npopulation - (MAX_AGE-AGE_ADULT);
                        }
                        npop += patch[p].n_population_oneyearagegroups->pop_size_per_gender_age1_risk[g][ai_npopulation][r];
                        
                        ai_positive = aa + patch[p].n_infected->youngest_age_group_index;
                        while(ai_positive > (MAX_AGE - AGE_ADULT - 1)){
                            ai_positive = ai_positive - (MAX_AGE-AGE_ADULT);
                        }
                        npositive += patch[p].n_infected->pop_size_per_gender_age1_risk[g][ai_positive][r];
                    }
                }
                
                if(npop > 0){
                    thisprev = npositive/(npop+0.0);
                }else{
                    thisprev = 0.0;
                }
            }else{
                printf("Not currently set up to fit to other types of prevalence data. Exiting\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            
            this_fit = perform_target_fit(&(fitting_data[patch[p].i_fit]), thisprev);
            
            // This is a target fit for 1 CI width - change from 1 to something else to 
            // make target wider/narrower.
            
            if(this_fit > 2){
                //printf("HIV prev. = %f lies outside target prev. [%lf,%lf] at time %f\n",
                //thisprev, fitting_data[i].prevalence_ll,
                //fitting_data[i].prevalence_ul, t0 + t_step*TIME_STEP);
                return(0);
            }
            
            /* Move to the next fitting condition if there is one, else exit the while loop. */
            if(patch[p].i_fit < (patch[p].n_fit - 1)){
                patch[p].i_fit = patch[p].i_fit + 1;
            }else{
                break;
            }
        }
    }
    /* If we get here, then we fitted any conditions there are for this timestep, so return 1. */
    return(1);
}


double perform_target_fit(fitting_data_struct *thisfittingdata, double model_prev){
   /* Check how close prevalence is to a point estimate given some upper and lower limits.
    
    Function returns the number of multiples of the lower (if below the point est) or upper CI that
    you are away.
    
    Parameters
    ----------
    thisfittingdata : pointer to a fitting_data_struct structure
    model_prev : double
    
    Returns
    -------
    
    */
    
    if(model_prev < thisfittingdata->prevalence_point_est){
        double sigma_ll = thisfittingdata->prevalence_point_est - thisfittingdata->prevalence_ll;
        double epsilon = thisfittingdata->prevalence_point_est - model_prev;
        
        if(sigma_ll > 0){
            //printf("Prevalence fitted within %lf\n", epsilon/sigma_ll);
            return(epsilon/sigma_ll);
        }else{
            printf("ERROR: cannot define fitting range - please check. Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }else if(model_prev > thisfittingdata->prevalence_point_est){
        double sigma_ul = thisfittingdata->prevalence_ul - thisfittingdata->prevalence_point_est;
        double epsilon = model_prev - thisfittingdata->prevalence_point_est;
        if(sigma_ul > 0){
            //printf("Prevalence fitted within %lf\n", epsilon/sigma_ul);
            return(epsilon/sigma_ul);
        }else{
            printf("ERROR: cannot define fitting range - please check. Exiting\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
    }else{
        /* if not > or <, then model_prev==thisfittingdata->prevalence_point_est. */
        printf("Fitting with 0??? %f %f\n", model_prev, thisfittingdata->prevalence_point_est);
        return(0);
    }
}
