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


/* HIV infection processes for the PopART model:
 * hiv_transmission_probability(): For a HIV+ partner determines their per-timestep P(trans) given their CD4 stage, SPVL, ART status, gender.
 * hiv_acquisition():   For a given individual, determines if HIV acquisition occurs through sexual contact with any infected partner (ie in partner_pairs_HIVpos[]) at a given timestep - using hiv_transmission_probability() to get per-individual transmission probabilities and assuming they sum up to total transmission probability.
 * find_who_infected(): Decides which HIV+ partner was responsible for seroconversion (takes into account heterogeneities in infectivity from each partner).
 * inform_partners_of_seroconversion(): Updates the "individual" structure of each partner so that they know the seroconvertor is now HIV+ (so if the partner is HIV-, this partnership is added to the serodiscordant partnership list of that partner, etc).
 * new_infection(): Assign SPVL, CD4 values for new infection in "individual" structure.
 * draw_initial_infection(): Bernoulli trial to decide if a given individual is HIV+ at introduction of HIV into population.
 * next_hiv_event(): For HIV+ people, schedule the next HIV event to occur in hiv_pos_progression[] (function called when they seroconvert, at the time of each new event, and when effective ART is stopped).  
 * carry_out_HIV_events_per_timestep():At each timestep this function is called, and carries out the HIV-related events scheduled in hiv_pos_progression[i][] at that timestep i.
 * get_window_result(): Decides if someone who was recently infeted tests positive corresponding to the window period of the test used (function can include different tests used over time/different countries). Note there are some notes on this in HIVtestnotesNov2014.pptx 
 * joins_preart_care(): Decides if someone who just tested positive joins pre-ART care or not. At present people chose to leave care independent of CD4 (ie not caring if they are eligible or not) - ASSUMPTION!
 * remains_in_cascade():Analogue of joins_preart_care() for people who are already in pre-ART care - decides if they remain in care until next CD4 test. 
 * measured_cd4_cat(): Given the 'true' biological CD4 category, gives the measured CD4 based on ATHENA data. 
 * art_cd4_eligibility_group(): gives the current  CD4 eligibility group for starting ART given the time and country.
 * is_eligible_for_art(): Outputs whether someone is eligible for ART at time t given their CD4 and current guidelines in that country.
 * get_time_emergency_start_ART(): For someone who is symptomatic (currently CD4<200) draws a time at which they start ART. At present if this is less than the time to AIDS death then they start on emergency ART at that time.
 * start_ART_process(): This represents the first X months of ART. Function sets up an individual starting ART (ie sets indiv->ART_status, stops CD4 progression), and schedules next cascade event (HIV death, becoming virally suppressed/not suppressed, drop out of care). 
 * draw_initial_hiv_tests(): When HIV testing starts in the country (assumed it can start before ART) this function schedules an HIV test for everyone. Currently people may have a test far in the future (so far that they never have one).
 * schedule_generic_cascade_event(): given a time t_event and a person indiv, schedules an event in the array cascade_event[] for that person. Note that this function does not know what type of event it is, so that is set separately in the function calling schedule_generic_cascade_event().
 * schedule_new_hiv_test(): function that sets up a HIV test in the future for someone. Called by draw_initial_hiv_tests() at the start of HIV testing in that country and by hiv_test_process() for people who just tested negative.
 * hiv_test_process(): Carries out the processes of HIV testing for someone and sets up their next cascade event according to the test result. In principle this can have sensitivity, specificity, window period and other complications to the test (only window period currently included). 
 * schedule_start_of_art(): Function which calculates the time to start ART in someone who just tested their CD4 and is eligible for ART, and then schedules event in cascade_event[]. Function called by hiv_test_process() and cd4_test_process(). 
 * cd4_test_process(): Carries out the processes of repeat CD4 testing for someone and sets up their next cascade event according to the test result. This is for people who previously tested HIV and CD4 and were not eligible, but remained in care. In principle we could also include the mis-match between true CD4 and measured CD4 category using the table from the ATHENA paper.
 * virally_suppressed_process(): Carries out the processes when become virally suppressed after starting ART. 
 * virally_unsuppressed_process(): Carries out the processes when become virally unsuppressed after starting ART. This includes restarting HIV progression events (althoguh HIV progression can be slower than in untreated people).
 * dropout_process(): Carries out the processes when drop out of care (restarting HIV progression if needed).
 * carry_out_cascade_events_per_timestep(): Function is the equivalent to carry_out_HIV_events_per_timestep() but for cascade events. At each timestep this function is called, and carries out the cascade-related events scheduled in cascade_events[i][] at that timestep i.
 */


#include "hiv.h"
#include "constants.h"
#include "output.h"
#include "utilities.h"
#include "partnership.h"
#include "interventions.h"
#include "debug.h"
#include "pc.h"

/* Function: get_spvl_cat()

Return set-point viral load (SPVL) category given SPVL.  SPVL categories are as follows:
0: <4
1: [4, 4.5)
2:[4.5,5)
3: >=5
The units are log10(SPVL) per mL.

Arguments
---------
spvl_num : double
     An individual's set-point viral load

Returns
-------
An individual's set-point viral load category. 
*/

int get_spvl_cat(double spvl_num){
    if(spvl_num < 4.0)
        return 0;
    else if(spvl_num < 4.5)
        return 1;
    else if(spvl_num < 5.0)
        return 2;
    else
        return 3;
}


/* Function: draw_inital_SPVL()
Generate set-point viral load (SPVL) for seeded infections based on log normal distribution:
 
gsl_ran_lognormal (const gsl_rng * r, double zeta, double sigma) returns a random variate from
p(x) dx = {1 \over x \sqrt{2 \pi \sigma^2} } \exp(-(\ln(x) - \zeta)^2/2 \sigma^2) dx
 
Arguments
---------
seroconverter : pointer to an individual structure
param : pointer to a parameters structure
    The parameters of the simulation, needs to have the attributes initial_SPVL_mu, 
    initial_SPVL_sigma, SPVL_sigma_E.  

Returns
-------
Nothing; sets the SPVL_num_E and SPVL_num_G attributes of the seroconverter individual
*/

void draw_inital_SPVL(individual *seroconverter, parameters *param){
    
    double SPVL;
    if(SPVL_INHERITANCE == 1){
        SPVL = param->initial_SPVL_mu + gsl_ran_gaussian(rng, param->initial_SPVL_sigma);
        seroconverter->SPVL_num_E = gsl_ran_gaussian(rng, param->SPVL_sigma_E);
        seroconverter->SPVL_num_G = SPVL - seroconverter->SPVL_num_E;
    }else if(SPVL_INHERITANCE == 0){
        seroconverter->SPVL_num_G = param->initial_SPVL_mu;
        seroconverter->SPVL_num_E = gsl_ran_gaussian(rng, param->SPVL_sigma_E);
    }
}


/* Function: inherit_spvl()
Determine set-point viral load (SPVL) of seroconverter given SPVL of infector

Determines an individual's set-point viral load after infection.  There is an environmental 
component (E) to an individual's SPVL, drawn from a normal distribution with mean
infector->SPVL_num_G and std dev SPVL_sigma_M, and a genetic component (G) that is inherited from
the infector, drawn from a Normal distribution with mean 0 and standard deviation of SPVL_sigma_E.

Arguments
---------
seroconverter : pointer to an individual structure
    The individual that is seroconverting
infector : pointer to an individual structure
    The individual that is the infector
param : pointer to a parameters structure
    Structure housing all parameters of the simulation

Returns
-------
Nothing; assigns SPVL_num_G and SPVL_num_E attributes to the seroconverter individual
*/

void inherit_spvl(individual *seroconverter, individual *infector, parameters *param){
    
    // Genetic component of log10(SPVL)
    seroconverter->SPVL_num_G = infector->SPVL_num_G + gsl_ran_gaussian(rng, param->SPVL_sigma_M);
    
    //Environmental conponent to log10(SPVL)
    seroconverter->SPVL_num_E = gsl_ran_gaussian(rng, param->SPVL_sigma_E);
}


/* Function: get_mean_time_hiv_progression()
Calculate mean time to next CD4 progression event for an individual, adjusting for SPVL, 
and whether progression is by SPVL category or using a Cox PH model for SPVL. 

Arguments
---------
param : pointer to a parameters structure
    Params of the IBM; should have an attribute called param.time_hiv_event[ncd4_cat][nSPVL_cat]
    which has dimensions of number of CD4 categories (ncd4_cat) and SPVL categores (nSPVL_cat).     
indiv : pointer to an individual structure
    Individual for which to calculate the mean time to next CD4-progression event.  

Returns
-------
mean_time_to_next_event : double
    Mean time (in years) until the next CD4 progression event for the individual `individual`.  
*/

double get_mean_time_hiv_progression(parameters *param, individual *indiv){
    
    double mean_time_to_next_event;
    
    if(CD4PROGRESSIONMODEL == USEFOURSPVLCATEGORIES){
        mean_time_to_next_event = param->time_hiv_event[indiv->cd4][indiv->SPVL_cat];
    }else{
        
        /* Get the log10 SPVL of the individual: */
        double logSPVL = indiv->SPVL_num_E + indiv->SPVL_num_G;
        
        if(logSPVL < COXMODELBASELINESPVL){
            mean_time_to_next_event = param->time_hiv_event[indiv->cd4][0];
        }else{
            // Each increase of 1.0 in log10(SPVL) increases the hazard of progression by
            // CoxPH_SPVL_RR_CD4progression[indiv->cd4].
            double CoxPHmultiplier = pow(param->CoxPH_SPVL_RR_CD4progression[indiv->cd4],
                logSPVL-COXMODELBASELINESPVL);
            // FOR DEBUGGING
            if(CoxPHmultiplier < 0){
                printf("Error: CoxPHmultiplier=%lf. Exiting\n", CoxPHmultiplier);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            
            mean_time_to_next_event = param->time_hiv_event[indiv->cd4][0]/CoxPHmultiplier;
            
            if (mean_time_to_next_event < 0.0 || mean_time_to_next_event > 20){
                printf("Error: mean_time_to_next_event=%lf. %lf %lf Exiting\n",
                    mean_time_to_next_event, param->time_hiv_event[indiv->cd4][0], CoxPHmultiplier);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
        }
    }
    return mean_time_to_next_event;
}


/* Function: get_RR_SPVL()
Determine relative rate of increased infectiousness as a function of set-point viral load.

Increasing Hill function of infectiousness as function of viral load.  Returns the following:
V^beta_k / (V^beta_k + beta_50^beta_k).  See Fraser et al., (2007) PNAS.  

Arguments
---------
SPVL : double
    log10 set point viral load
param : pointer to a parameters structure
    parameters of the simulation (should have attributes SPVL_beta_k and SPVL_beta_50)
    SPVL_beta_50: viral load at which infectiousness is half its maximum
    SPVL_beta_k: steepness of the increase in infectiousness as a function of viral load

Returns
-------
double: V^beta_k / (V^beta_k + beta_50^beta_k)
*/

double get_RR_SPVL(double SPVL, parameters *param){
    // Get the SPVL from the log10(SPVL)
    double v = pow(10, SPVL);
    
    return pow(v, param->SPVL_beta_k) / 
        (pow(v, param->SPVL_beta_k) + pow(param->SPVL_beta_50, param->SPVL_beta_k));
}


/* Function: hiv_transmission_probability()
Determine per-timestep probability of transmission from an HIV+ partner 

Calculates a per-timestep probability of transmission to an HIV- partner from an HIV+ partner 
according to characteristics of the HIV+ partner such as CD4 stage, SPVL, ART status, gender, 
stage of HIV (i.e. acute).  


Arguments
---------
HIVpos_partner : pointer to an individual structure
    A pointer to the HIV+ partner
param : pointer to a parameters structure
    A pointer to a structure housing the parameters of the simulation.  

Returns
-------
hazard : double
    The per-timestep probabiliy of transmission to HIV- partner from HIV+ partner


Notes
-----
This function used to have code associated with the PANGEA project that used the function
intervention_effectiveness() to look at intervention scenarios that had different efficacies.  
The code was used in the following such manner:  
```c
if (HIVpos_partner->HIV_status == ACUTE){
    hazard = hazard*param->RRacute_trans;
}else{
...
    if (PHYLO_SCENARIO==2)
        hazard = hazard*param->RRCD4[HIVpos_partner->cd4] * 
            intervention_effectiveness(0.6, 2004, 2014, t);
```
Using `git show 649265ba` (a commit on 31 August 2018) will show the details of this code.  
*/

double hiv_transmission_probability(individual* HIVpos_partner, parameters *param){
    
    double hazard;
    
    // Assign hazard as baseline hazard
    hazard = param->average_annual_hazard;
    
    // Adjust hazard according to whether the partner is in the acute phase of HIV (else according 
    // to their CD4 category)
    if (HIVpos_partner->HIV_status == ACUTE){
        hazard *= param->RRacute_trans;
    }else{
        hazard *= param->RRCD4[HIVpos_partner->cd4];
    }
    // Adjust hazard according to SPVL
    hazard *= get_RR_SPVL(HIVpos_partner->SPVL_num_G + HIVpos_partner->SPVL_num_E, param);
    
    // Reduce infectivity if on ART
    if(HIVpos_partner->ART_status == EARLYART){
        hazard *= param->RR_ART_INITIAL;
    }else if(HIVpos_partner->ART_status == LTART_VS){
        hazard *= param->RR_ART_VS;
    }else if(HIVpos_partner->ART_status == LTART_VU){
        hazard *= param->RR_ART_VU;
    }
    
    // Adjust hazard according to gender of the partner
    if(HIVpos_partner->gender == MALE){
        hazard *= param->RRmale_to_female_trans;
    }
    return hazard;
}


/* Function: hiv_acquisition()
Determine whether a single individual in one (or more) serodiscordant partnerships gets
infected in a timestep.  If the individual gets infected then update all relevant lists (they 
are removed from list of HIV- people in serodiscordant partnerships).  (Not currently 
implemented - they will be added to list of people who are HIV+ with a time to next HIV 
progression event scheduled.  If infected addition output is generated (update incident and 
prevalent case counts, output phylogenetic data of interest). 

Arguments
---------
susceptible : pointer to an individual strucutre
    Pointer to the susceptible individual
time_infect : double
    Time (in decimal years) of infection
patch : pointer to a patch_struct structure
    Patch structure including information on all patches
p : int
    Patch number where the susceptible individual resides
overall_partnerships : pointer to an all_partnerships structure
output : pointer to an output_struct structure
debug : pointer to a debug_struct structure
file_data_store : pointer to a file_struct structure

Returns
-------
*/

void hiv_acquisition(individual* susceptible, double time_infect, patch_struct *patch, int p,
    all_partnerships *overall_partnerships, output_struct *output, debug_struct *debug, 
    file_struct *file_data_store, int t0, int t_step){

    double partner_SPVL_G; /* Environmental component of log-10 SPVL. */
    double partner_SPVL_E; /* Genetic component of log-10 SPVL. */
    individual *infector; /* Pointer to the person who infects. Makes code more readable. */
    
    double temp_hazard;   /* Store hazard including reduced susceptibility from VMMC if needed - used when we want to output - print to file - hazards from everyone over time. */
    
    if(susceptible->id == FOLLOW_INDIVIDUAL && susceptible->patch_no == FOLLOW_PATCH){
        printf("checking HIV acquisition for adult %ld from patch %d at time %6.4f\n",
            susceptible->id, susceptible->patch_no, time_infect);
        fflush(stdout);
    }
    if((susceptible->HIV_status > UNINFECTED) || (susceptible->cd4 != CD4_UNINFECTED)){
        printf("ERROR: Trying to infect an HIV+ person %li in patch %i\n",
            susceptible->id, susceptible->patch_no);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* This is the total (annual) hazard from all HIV+ partners, ignoring reduced susceptibility 
    if the individual is a circumcised man. */
    double total_hazard_ignore_circ = 0.0;

    /* This variable is total_hazard_ignore_circ converted to the hazard per timestep, 
     * and adjusted for circumcision if needed. */
    double total_hazard_per_timestep;

    /* Characteristics of the seropositive partner which may influence transmission.
     * Assumed to be heterosexual partners always: */
    int partner_gender = 1 - susceptible->gender;

    int i; /* Index for summing over partners. */

    int infector_index; /* Index (in partner_pairs_HIVpos[]) of partner who infected the seroconverter. */

    /* The number of HIV+ partners of this susceptible. */
    int npos = susceptible->n_HIVpos_partners;

    /* FOR DEBUGGING */
    if(npos < 1){
        printf("Problem: trying to infect someone who is in no serodiscordant partnership\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    if(susceptible->HIV_status!=UNINFECTED){
        printf("Problem: trying to infect someone who is already infected. HIV status = %i\n",
            susceptible->HIV_status);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* This will create an alias for each seropositive partner to save typing. */
    individual* temp_HIVpos_partner;

    /* Now sum hazard over each HIV+ partner: */
    for(i = 0; i < npos; i++){
        /* This is a pointer to the ith HIV+ partner of the susceptible: */
        temp_HIVpos_partner = susceptible->partner_pairs_HIVpos[i]->ptr[partner_gender];
        //partner_log_vl = log10( temp_HIVpos_partner -> viral_load);
        /* This is 0 if partner not in acute phase, 1 if they are. */
        //partner_acute = 2 - (temp_HIVpos_partner -> HIV_status);

        /* Store the hazard from each partnership so we can work out who infected whom: */
        //PER_PARTNERSHIP_HAZARD_TEMPSTORE[i] = hill_up(partner_log_vl,param->average_annual_hazard,1.0,param->average_log_viral_load)*(1+param->RRacute_trans*partner_acute);
        PER_PARTNERSHIP_HAZARD_TEMPSTORE[i] = hiv_transmission_probability(temp_HIVpos_partner,
            patch[p].param);

        // Add a multiplying factor for assortative risk mixing: high-high is twice as risky 
        // and low-low half as risky as other combinations 
        if(CHANGE_RR_BY_RISK_GROUP == 1){
            if(temp_HIVpos_partner->sex_risk == HIGH && susceptible->sex_risk == HIGH){
                PER_PARTNERSHIP_HAZARD_TEMPSTORE[i] *= 2.0;
            }
            if(temp_HIVpos_partner->sex_risk == LOW && susceptible->sex_risk == LOW){
                PER_PARTNERSHIP_HAZARD_TEMPSTORE[i] *= 0.5;
            }
        }

        // Add a multiplying factor for assortative community mizing: within community is more
        // risky than between community, reflecting higher frequency of sex act and lower condom
        // use 
        if(temp_HIVpos_partner->patch_no != susceptible->patch_no){
            PER_PARTNERSHIP_HAZARD_TEMPSTORE[i] *= patch[p].param->rr_hiv_between_vs_within_patch;
        }

        //printf("PER_PARTNERSHIP_HAZARD_TEMPSTORE[%i]=%lf r=%i %i\n",i,PER_PARTNERSHIP_HAZARD_TEMPSTORE[i],temp_HIVpos_partner->sex_risk,susceptible->sex_risk);
        total_hazard_ignore_circ += PER_PARTNERSHIP_HAZARD_TEMPSTORE[i];

        if(WRITE_HAZARDS == 1 && PRINT_EACH_RUN_OUTPUT == 1){
            /* Note that "time_infect" is actually time at which hiv_acquisition() is called - so represents time at which we are looking at hazard. */
            /* Only output over fixed time. */
            if(time_infect >= 2014.0 && time_infect < 2015 && p == 0){
                // Adjust hazard based upon circumcision status
                if((susceptible->circ) == UNCIRC || (susceptible->circ) == UNCIRC_WAITING_VMMC){
                    temp_hazard = PER_PARTNERSHIP_HAZARD_TEMPSTORE[i];
                }else if((susceptible->circ) == VMMC){
                    temp_hazard = PER_PARTNERSHIP_HAZARD_TEMPSTORE[i] * 
                            (1.0-patch[p].param->eff_circ_vmmc);
                }else if((susceptible->circ) == TRADITIONAL_MC){
                    temp_hazard = PER_PARTNERSHIP_HAZARD_TEMPSTORE[i] * 
                            (1.0-patch[p].param->eff_circ_tmc);
                }else if((susceptible->circ) == VMMC_HEALING){
                    /* Increased susceptibility if in healing period: */
                    temp_hazard = PER_PARTNERSHIP_HAZARD_TEMPSTORE[i] * 
                            patch[p].param->rr_circ_unhealed;
                }
                output_hazard_over_time_period(time_infect, temp_hazard, susceptible, 
                    temp_HIVpos_partner, file_data_store, output);
            }
        }
    }

    /* Adjust according to the circumcision status of the susceptible: */
    if((susceptible->circ) == UNCIRC || (susceptible->circ) == UNCIRC_WAITING_VMMC){
        total_hazard_per_timestep = total_hazard_ignore_circ * TIME_STEP;
    }else if((susceptible->circ) == VMMC){
        total_hazard_per_timestep = total_hazard_ignore_circ * 
                (1.0 - patch[p].param->eff_circ_vmmc) * TIME_STEP;
    }else if((susceptible->circ) == TRADITIONAL_MC){
        total_hazard_per_timestep = total_hazard_ignore_circ * 
                (1.0 - patch[p].param->eff_circ_tmc) * TIME_STEP;
    /* Increased susceptibility if in healing period: */
    }else if((susceptible->circ) == VMMC_HEALING){
        total_hazard_per_timestep = total_hazard_ignore_circ * 
            patch[p].param->rr_circ_unhealed * TIME_STEP;
    }else{
        fprintf(stderr,"ERROR: unknown circumcision status!!! Exiting\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    //printf("Individual %ld is subject to infection hazard %lg\n", susceptible->id, total_hazard_per_timestep);

    /* Now see if transmission occurs (Bernoulli trial - could replace by gsl_ran_bernoulli(rng,total_hazard_per_timestep): */ 
    double x = gsl_rng_uniform (rng);
    if(x <= total_hazard_per_timestep){

        /* Determine which of the seropositive partners was responsible for infecting.
         * so MMC status drops out when working out who the infector is (simply because VMMC reduces susceptibility,
         * To do this, pass total_hazard_ignore_circ - neither this nor PER_PARTNERSHIP_HAZARD_TEMPSTORE[] are scaled by eff_circ,
         * not infectivity). It is a bit faster to do it this way.
         * The output is the index in array partner_pairs_HIVpos of the individual person who infected. */
        /* Note if there was only one partner then that was the person responsible (so we don't draw as many random numbers). */
        if(npos == 1){
            infector_index = 0;
        }else{
            infector_index = find_who_infected(npos, total_hazard_ignore_circ);
        }
        
        /* Transmit HIV to this individual. */
        if(susceptible->id == FOLLOW_INDIVIDUAL && susceptible->patch_no == FOLLOW_PATCH){
            printf("Transmitting HIV infection to adult %ld in patch %d at time %6.4f\n",
                susceptible->id, susceptible->patch_no, time_infect);
            fflush(stdout);
        }
        infector = susceptible->partner_pairs_HIVpos[infector_index]->ptr[partner_gender];

        partner_SPVL_G = infector->SPVL_num_G;
        partner_SPVL_E = infector->SPVL_num_E;
        
        /* Store the partner's SPVL to look at heritability. */
        susceptible->SPVL_infector = partner_SPVL_G + partner_SPVL_E;

        new_infection(time_infect, FALSE, susceptible, infector , patch[p].n_infected, 
            patch[p].n_newly_infected, patch[p].age_list, patch[p].param, 
            patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression, 
            patch[p].size_hiv_pos_progression, patch[p].n_infected_cumulative, file_data_store);
        
        // Increment counters
        patch[p].n_newly_infected_total++;
        patch[p].n_newly_infected_total_by_risk[susceptible->sex_risk]++;
        
        if(susceptible->patch_no != infector->patch_no){
            patch[p].n_newly_infected_total_from_outside++;
        }
        if(infector->HIV_status == ACUTE){
            patch[p].n_newly_infected_total_from_acute++;
        }
        
        // Increment PConly outputs
        
        // Find the age of the newly infected individual
        int age = (int) floor(time_infect - susceptible->DoB);
        if((age >= AGE_PC_MIN) && (age <= AGE_PC_MAX)){
            patch[p].n_newly_infected_total_pconly++;
            patch[p].n_newly_infected_total_by_risk_pconly[susceptible->sex_risk]++;
            if(susceptible->patch_no != infector->patch_no){
                patch[p].n_newly_infected_total_from_outside_pconly++;
            }
            if(infector->HIV_status == ACUTE){
                patch[p].n_newly_infected_total_from_acute_pconly++;
            }
            
            // Does this infection take place during a PC round?  
            // If so, increment the incident infections in the PC age range.  
            int pc_round = get_pc_round(t0, t_step, patch, 0);
            
            if(pc_round != -1){
                output->PC_ROUND_INFECTIONS[p][susceptible->gender][age - AGE_PC_MIN][pc_round]++;
            }
        }

        patch[p].DEBUG_NHIVPOS++;
        patch[p].PANGEA_N_ANNUALINFECTIONS++;
        
        // Increment incident infections
        int g = susceptible->gender;
        if(age >= MAX_AGE){
            age = MAX_AGE;
        }
        int age_idx = FIND_AGE_GROUPS_UNPD[age - AGE_ADULT];
        int year_idx = (int) floor(time_infect - patch[p].param->start_time_simul);
        
        patch[p].calendar_outputs->N_calendar_infections[g][age_idx][year_idx]++;

        /* Update each of their HIV- partners so that they know that they have a serodiscordant partner and update list of serodiscrodant partnerships */
        inform_partners_of_seroconversion_and_update_list_serodiscordant_partnerships(susceptible, 
            overall_partnerships->susceptible_in_serodiscordant_partnership, 
            overall_partnerships->n_susceptible_in_serodiscordant_partnership);

        // Outputs of interest (note that PER_PARTNERSHIP_HAZARD_TEMPSTORE[] is a global array 
        //so we do not need to pass it): */

        if(susceptible->partner_pairs_HIVpos[infector_index]->ptr[partner_gender]->HIV_status == ACUTE){
            patch[p].PANGEA_N_ANNUALACUTEINFECTIONS++;
        }

        /*printf("-- time_infect %lg\n",time_infect);
        printf("-- susceptible->id %ld\n",susceptible->id);
        printf("-- infector->id %ld\n",susceptible -> partner_pairs_HIVpos[infector_index] -> ptr[partner_gender]->id);*/

        // Phylogenetic transmission outputs
        if(WRITE_PHYLOGENETICS_OUTPUT == 1){
            store_phylogenetic_transmission_output(output, time_infect, susceptible, 
                susceptible->partner_pairs_HIVpos[infector_index]->ptr[partner_gender], 
                file_data_store, t0, t_step);
        }
        //printf("Individual %d is infected by individual %d\n",
        //  susceptible->id,
        //    susceptible->partner_pairs_HIVpos[infector_index]->ptr[partner_gender]->id);
    }
}


/* Function: find_who_infected()
Determine which partner transmitted HIV infection if someone has more than one HIV+ partner.

Draw an index for the partner that is responsible weighted by the hazard of infection from 
that partner (stored in a global array called PER_PARTNERSHIP_HAZARD_TEMPSTORE[]).  


Arguments
---------
numberpos_partners : int
    Number of HIV+ partners
total_hazard : double
    Total probability of seroconversion in the timestep of interest.  Contributions to this 
    hazard are stored in the (global) array PER_PARTNERSHIP_HAZARD_TEMPSTORE[].  

Returns
-------
Index of the partner who transmitted the infection which can be used to reference the array 
susceptible->partner_pairs_HIVpos[] where `susceptible` is the individual structure of the 
susceptible partner.  
*/

int find_who_infected(int numberpos_partners, double total_hazard){
    
    double temp_hazard = 0.0;
    int partner_i = 0;

    // Pick a random number uniformly in [0, total_hazard)
    double x = gsl_rng_uniform(rng) * total_hazard;
    
    temp_hazard = PER_PARTNERSHIP_HAZARD_TEMPSTORE[partner_i];
    
    while((temp_hazard < x) && (partner_i < numberpos_partners - 1)){
        partner_i += 1;
        temp_hazard += PER_PARTNERSHIP_HAZARD_TEMPSTORE[partner_i];
    }
    return partner_i;
}


/* Function: inform_partners_of_seroconversion_and_update_list_serodiscordant_partnerships()
Update lists of serodiscorant partnerships following a seroconversion.  

Function firstly sets the n_HIVpos_partners attribute of the seroconverter to zero (note:
this variable is only >0 if the individual is HIV-, as well as has HIV+ partners, and it is -1
for someone who is HIV+).  Then the function goes through each partner of the seroconverter in
turn to see which are HIV-.
For each HIV- partner they are informed that they now have a seropositive partner (ie the
seroconverter), and if this HIV- partner did not have any seropositive partners before they
were not in the susceptible_in_serodiscordant_partnership list, and hence are now added to this
list. 

Note: n_partners can be zero for imported cases (i.e. those in which we seed the infection).  
This should be accounted for in this function.  

Arguments
---------
seroconverter : pointer to an individual
    The individual who has just seroconverted
susceptible_in_serodiscordant_partnership : pointer to an array of individuals
    
n_susceptible_in_serodiscordant_partnership : int

Returns
-------
Nothing; adjusts the partnerships attributes of the seroconverter and each of the 
seroconverter's partners, and 
*/

void inform_partners_of_seroconversion_and_update_list_serodiscordant_partnerships(
    individual* seroconverter, individual** susceptible_in_serodiscordant_partnership, 
    long *n_susceptible_in_serodiscordant_partnership){
    
    int i;
    
    // The gender of (all) partners of this seroconverter (since the model only has heterosexuals)
    int partner_gender = 1 - seroconverter->gender;
    
    // Set up a pointer to each partner in turn; an alias to save on code
    individual* temp_partner;
    
    // Update the list of serodiscordant and seropositive partners.
    // This person has just seroconverted, so set the count of HIV+ and serodiscordant partners to
    // zero as we no longer care about them. The list of indices does not need changing as this
    // counter tells us everything we need to know.  
    seroconverter->n_HIVpos_partners = 0;
    
    // Loop through the partners of the seroconverter
    for(i = 0; i < seroconverter->n_partners; i++){
        
        // Find pointer to the ith partner of the seroconverter
        // (it points to an already allocated person so no malloc needed)
        temp_partner = seroconverter->partner_pairs[i]->ptr[partner_gender]; 
        
        // Check if the partner is currently seronegative (do nothing if they're HIV+)
        if(temp_partner->HIV_status == UNINFECTED){
            // Note: temp_partner->n_HIVpos_partners is the number of HIV+ partners
            // of this partner (prior to the current seroconversion).

            // Add the current serconverter to the array of HIVpos partners of temp_partner
            temp_partner->partner_pairs_HIVpos[temp_partner->n_HIVpos_partners] =
                seroconverter->partner_pairs[i];
            
            // Increment this partner's number of HIV+ partners (n_HIVpos_partners)
            temp_partner->n_HIVpos_partners++;
            
            // Record if this is a between-patch partnership
            if(temp_partner->patch_no != seroconverter->patch_no){
                temp_partner->n_HIVpos_partners_outside++;
            }
            
            // This partnership has become serodiscordant so, unless the HIV- temp_partner was
            // already in the list of susceptible_in_serodiscordant_partnership, they have to be
            // added there
            
            // This means before seroconversion of seroconverter, 
            // this individual was in no serodiscordant partnerships
            if(temp_partner->idx_serodiscordant == -1 ){
                add_susceptible_to_list_serodiscordant_partnership(temp_partner,
                    susceptible_in_serodiscordant_partnership,
                    n_susceptible_in_serodiscordant_partnership);
            }
        } // We do nothing if the partner is currently seropositive
    }
    // Remove seroconverter from list of suscepts in serodiscordant partnerships, if appropriate
    remove_susceptible_from_list_serodiscordant_partnership(seroconverter,
        susceptible_in_serodiscordant_partnership, n_susceptible_in_serodiscordant_partnership);
}


/* Function: new_infection()
Set variables related to HIV in the individual structure for a newly infected individual. 

Allocates what the next hiv-related event to happen (after acute infection if relevant) is: CD4 progression or AIDS death/emergency start ART
Also updates the counts of incident and prevalent cases. 

*/
void new_infection(double time_infect, int SEEDEDINFECTION, individual* seroconverter, 
    individual *infector, population_size_one_year_age *n_infected, 
    population_size_one_year_age *n_newly_infected, age_list_struct *age_list, 
    parameters *param, individual ***hiv_pos_progression, long *n_hiv_pos_progression, 
    long *size_hiv_pos_progression, population_size_one_year_age *n_infected_cumulative,
    file_struct *file_data_store){

    /* SEEDEDINFECTION tells us whether this is one of the initial infections seeded at time (start_time_hiv...start_time_hiv+n_years_HIV_seeding), or an infection acquired since then.
     * For seeded infections we start them in chronic infection (to avoid a wave of initial acutes). */

    int g = seroconverter->gender;
    long ncheck; 
    int aa;
    double x;
    double time_to_next_event;
    int idx_next_event, idx_current_time;

    if(seroconverter->id==FOLLOW_INDIVIDUAL && seroconverter->patch_no==FOLLOW_PATCH){
        printf("New HIV infection of adult %ld from patch %d and gender %d at time %6.4f\n",seroconverter->id,seroconverter->patch_no,g,time_infect);
        fflush(stdout);
    }


    if(seroconverter->id==FOLLOW_INDIVIDUAL && seroconverter->patch_no==FOLLOW_PATCH)
    {
        //print_here_string("Start new_infection function",0);
        print_specific_IDs_by_age(FOLLOW_INDIVIDUAL,age_list,FOLLOW_PATCH);
    }

    // Store the time at which they got infected
    seroconverter->t_sc = time_infect;

    // Draw SPVL
    if(SEEDEDINFECTION == TRUE)
        // If a seeded infection (ie no modelled infector) draw from a set distribution.
        draw_inital_SPVL(seroconverter, param);
    else{
         // If infection came from another modelled individual draw the seroconverter's 
        // SPVL based on that of the infector, else assume as a seeded infection.  
        if(SPVL_INHERITANCE == 1){
            inherit_spvl(seroconverter, infector, param);
        }else{
            draw_inital_SPVL(seroconverter, param);
        }
    }
    double logSPVL = seroconverter->SPVL_num_G + seroconverter->SPVL_num_E;
    seroconverter->SPVL_cat = get_spvl_cat(logSPVL);

    // Output the SPVL of seeded individuals. 
    // Works best with quite a large seeded HIV+ population.
    if(SEEDEDINFECTION == TRUE && WRITE_DEBUG_INITIAL_SPVL_DISTRIBUTION == 1){
        write_initial_spvl_distribution(seroconverter, file_data_store);
    }

    // Draw CD4 category for end of acute phase (if not seeded infection), 
    // otherwise where a seeded infection starts:
    x = gsl_rng_uniform(rng);
    if(x < param->cumulative_p_initial_cd4_gt500[seroconverter->SPVL_cat]){
        seroconverter->cd4 = 0;
    }else if(x < param->cumulative_p_initial_cd4_350_500[seroconverter->SPVL_cat]){
        seroconverter->cd4 = 1;
    }else if(x < param->cumulative_p_initial_cd4_200_350[seroconverter->SPVL_cat]){
        seroconverter->cd4 = 2;
    }else{
        seroconverter->cd4 = 3;
    }
    
    if(WRITE_DEBUG_CD4_AFTER_SEROCONVERSION == 1){
        write_cd4_at_seroconversion(seroconverter, file_data_store);
    }

    if(seroconverter->id == FOLLOW_INDIVIDUAL && seroconverter->patch_no == FOLLOW_PATCH){
        printf("Initial CD4 for adult %ld from patch %d at time %6.4f is %i\n",
            seroconverter->id, seroconverter->patch_no, time_infect, seroconverter->cd4);
        fflush(stdout);
    }
    
    // For non-seed cases
    if(SEEDEDINFECTION == FALSE){
        // Person is now in acute stage.
        seroconverter->HIV_status = ACUTE;
        
        // Draw time to end of acute phase
        time_to_next_event = param->min_dur_acute + 
                gsl_rng_uniform(rng) * (param->max_dur_acute - param->min_dur_acute);
        
        if(seroconverter->id == FOLLOW_INDIVIDUAL && seroconverter->patch_no == FOLLOW_PATCH){
            printf("Individual %ld from patch %d", seroconverter->id, seroconverter->patch_no);
            printf(" is scheduled to progress from acute to chronic at ");
            printf("t=%6.4f. Currently t=%6.4f\n", time_to_next_event + time_infect, time_infect);
            fflush(stdout);
        }
        
        // In the absence of PopART assume that progression is just to next CD4 stage.
        seroconverter->next_HIV_event = HIVEVENT_ENDOFACUTE;
    }else{ 
        // For seeded infections
        
        // Assume person is in chronic stage
        seroconverter->HIV_status = CHRONIC;

        // Calculate the mean time to next event, adjusting for SPVL.
        double mean_time_to_next_event;
        mean_time_to_next_event = get_mean_time_hiv_progression(param, seroconverter);

        // Draw time to next CD4 stage
        time_to_next_event = gsl_ran_exponential(rng, mean_time_to_next_event);
        
        if(seroconverter->id == FOLLOW_INDIVIDUAL && seroconverter->patch_no == FOLLOW_PATCH){
            printf("Individual %ld from patch %d", seroconverter->id, seroconverter->patch_no);
            printf(" is seeded HIV+ at time %6.4f, ", time_infect);
            printf("CD4 stage = %i. Scheduled to progress to next CD4 stage at t = %6.4f\n",
                seroconverter->cd4, time_to_next_event + time_infect);
            fflush(stdout);
        }
        // These are seeded infections so assume there is NO ART at this point. 
        // So progression is just to next CD4 stage
        if(seroconverter->cd4 < 3){
            seroconverter->next_HIV_event = HIVEVENT_CD4_PROGRESSION;
        }else{
            seroconverter->next_HIV_event = HIVEVENT_AIDSDEATH;
        }
    }

    idx_next_event = (int) round(((time_infect-param->start_time_hiv) +
            time_to_next_event) * N_TIME_STEP_PER_YEAR);
    
    idx_current_time = (int) round(((time_infect-param->start_time_hiv)) *
        N_TIME_STEP_PER_YEAR);

    // Make sure that an event is not scheduled to occur in the current time-step
    if(idx_next_event == idx_current_time){
        idx_next_event++;
    }
    
    seroconverter->DEBUGTOTALTIMEHIVPOS += time_to_next_event;
    seroconverter->PANGEA_t_prev_cd4stage = time_infect;
    seroconverter->PANGEA_t_next_cd4stage = time_infect + time_to_next_event;
    
    if(time_infect+time_to_next_event < param->end_time_simul){
        seroconverter->idx_hiv_pos_progression[0] = idx_next_event;
        seroconverter->idx_hiv_pos_progression[1] =
            n_hiv_pos_progression[idx_next_event];

        // Check if we've run out of memory:
        if(n_hiv_pos_progression[idx_next_event] >=
            (size_hiv_pos_progression[idx_next_event])){
            
            // Note that realloc does not work (we need to pass a pointer to the pointer, 
            // which is really complicated as it propagates through several functions (so maybe
            // make planned_breakups[time_breakup] **), so ditch this code for now and use the
            // following lines: 
            printf("Unable to re-allocate hiv_pos_progression[i]. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        hiv_pos_progression[idx_next_event][n_hiv_pos_progression[idx_next_event]] = seroconverter;
        n_hiv_pos_progression[idx_next_event]++;
    }else{
        seroconverter->idx_hiv_pos_progression[0] = EVENTAFTERENDSIMUL;
        seroconverter->idx_hiv_pos_progression[1] = -1;

        if(seroconverter->id==FOLLOW_INDIVIDUAL && seroconverter->patch_no == FOLLOW_PATCH){
            printf("Not scheduling first HIV progression event after seroconversion ");
            printf("for individual %ld from patch %d as after end of simulation at t=%6.4f\n", 
                seroconverter->id, seroconverter->patch_no, time_infect + time_to_next_event);
            fflush(stdout);
        }
    }

    // Add seroconverter to lists of infected individuals, in the correct age category
    aa = (int) floor(floor(time_infect) - seroconverter->DoB) - AGE_ADULT;
    
    /* DEBUGGING - CAN REMOVE ALL THE DIFFERENT AI HERE */
    if(aa < (MAX_AGE - AGE_ADULT)){

        /* ai is the age index of the array n_infected->pop_size_per_gender_age1_risk[g][ai][r] for the person with DoB as above at t_infect. */

        // Indices for the prevalence, incidence and ???age???
        int ai_prev = n_infected->youngest_age_group_index + aa;
        while (ai_prev>(MAX_AGE-AGE_ADULT-1))
            ai_prev = ai_prev - (MAX_AGE-AGE_ADULT);

        int ai_inc = n_newly_infected->youngest_age_group_index + aa;
        while (ai_inc>(MAX_AGE-AGE_ADULT-1))
            ai_inc = ai_inc - (MAX_AGE-AGE_ADULT);

        int ai_inc_c = n_infected_cumulative->youngest_age_group_index + aa;
        while (ai_inc_c>(MAX_AGE-AGE_ADULT-1))
            ai_inc_c = ai_inc_c - (MAX_AGE-AGE_ADULT);

        int ai_age = age_list->age_list_by_gender[g]->youngest_age_group_index + aa;
        while (ai_age>(MAX_AGE-AGE_ADULT-1))
            ai_age = ai_age - (MAX_AGE-AGE_ADULT);
        
        /* looking for the seroconverter in the age_list --> presumably only for debugging, could get rid of this in final code to speed up. */
        ncheck = 0;
        while ((ncheck<age_list->age_list_by_gender[g]->number_per_age_group[ai_age]) && (age_list->age_list_by_gender[g]->age_group[ai_age][ncheck]->id!=seroconverter->id)){
            ncheck++;
        }

        if((ncheck>=age_list->age_list_by_gender[g]->number_per_age_group[ai_age])){
            printf("PROBLEM in new_infection() at time %f: Person not found %li from patch %d: should be in age group %i\n",time_infect,seroconverter->id,seroconverter->patch_no,ai_age);
            printf("time = %6.4f, DoB= %6.4f, aa = %i\n",time_infect,seroconverter->DoB,aa);

            /* For debugging: */
            find_in_age_list(time_infect,seroconverter,age_list,param);
            printf("Patch %d: ",seroconverter->patch_no);
            fflush(stdout);
            print_specific_IDs_by_age(FOLLOW_INDIVIDUAL,age_list,FOLLOW_PATCH);

            print_individual(seroconverter);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        // Add the seroconverter to n_infected in the correct age group
        (n_infected->pop_size_per_gender_age1_risk[seroconverter->gender][ai_prev][seroconverter->sex_risk])++;

        /* adding the seroconverter to n_newly_infected in the right age group */
        (n_newly_infected->pop_size_per_gender_age1_risk[seroconverter->gender][ai_inc][seroconverter->sex_risk])++;

        (n_infected_cumulative->pop_size_per_gender_age1_risk[seroconverter->gender][ai_inc_c][seroconverter->sex_risk])++;
        //printf("+++ One new HIV+ \n");
        //fflush(stdout);
    }else{

        /* looking for the seroconverter in the age_list --> presumably only for debugging, could get rid of this? */
        ncheck = 0;
        while ((ncheck<age_list->age_list_by_gender[g]->number_oldest_age_group) && (age_list->age_list_by_gender[g]->oldest_age_group[ncheck]->id!=seroconverter->id)){
            ncheck++;
        }
        if ((ncheck>=age_list->age_list_by_gender[g]->number_oldest_age_group))
            printf("PROBLEM: Person not found in oldest age group: %li in patch %d gender %d\n",seroconverter->id,seroconverter->patch_no,g);

        /* adding the seroconverter to n_infected in the right age group */
        (n_infected->pop_size_oldest_age_group_gender_risk[seroconverter->gender][seroconverter->sex_risk])++;

        /* adding the seroconverter to n_newly_infected in the right age group */
        (n_newly_infected->pop_size_oldest_age_group_gender_risk[seroconverter->gender][seroconverter->sex_risk])++;
        (n_infected_cumulative->pop_size_oldest_age_group_gender_risk[seroconverter->gender][seroconverter->sex_risk])++;
    }
}

/* Function: draw_initial_infection()
Seed initial HIV infection.

Determine if an individual from the population will be a seeding HIV infection. The proportion of 
individuals in the population that are converted to seeding cases (for each year or HIV seeding) is
stored in the array: patch[p].param->initial_prop_infected_gender_risk[gender][risk].  This function
is called within simul.c when HIV seeding is started.  For each individual in the population, the
test of whether or not they are a seed cases is a Bernoulli trial, infecting individual indiv with
probability param->initial_prop_infected_gender_risk[g][r] where g and r are the gender and risk of
indiv. 

Arguments
---------
t : double
    Year
indiv : pointer to an individual structure
    An individual in the population, a potential seed case
patch : pointer to a patch structure
p : int
    Patch number
overall_partnerships : pointer to an all_partnerships structure
output : pointer to a output_struct structure
    Used for the called functions that might write something to file
file_data_store : pointer to a file_struct structure
    Used for the called functions that might write something to file

Returns
-------
Nothing; sets an input individual (indiv) as a seeding case according to a Bernoulli trial.  If the
individual is a seeding case this function calls new_infection(), increments the number of infected
cases, sets the individual's SPVL to -1, and updates lists of serodiscorant partnerships following
a seroconversion.  
*/

void draw_initial_infection(double t, individual* indiv, patch_struct *patch, int p,
    all_partnerships *overall_partnerships, output_struct *output, file_struct *file_data_store, 
	int t0, int t_step){
    
    // Run a couple of checks to begin with
    if(indiv->cd4 == DUMMYVALUE){
        printf("Error. Using an uninitialised person.\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    if(indiv->cd4 == DEAD){
        printf("Trying to make a dead person acquire HIV at HIV introduction.\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    // No need to do anything if already HIV+.
    if(indiv->HIV_status > UNINFECTED){
        return;
    }
    
    // Draw a random uniform number and determine if this individual is a seed case or not
    double random = gsl_rng_uniform(rng);
    if(random < patch[p].param->initial_prop_infected_gender_risk[indiv->gender][indiv->sex_risk]){
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Seeding HIV infection of adult %ld in patch %d at time %6.4f\n", 
                indiv->id, indiv->patch_no,t);
            fflush(stdout);
        }
        
        // Generate a new infection
        // 'NULL' reflects that this is a seeded infection (normally the pointer to the infector).
        new_infection(t, TRUE, indiv, NULL, patch[p].n_infected, patch[p].n_newly_infected,
            patch[p].age_list, patch[p].param, patch[p].hiv_pos_progression,
            patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression,
            patch[p].n_infected_cumulative, file_data_store);
        
        // Increment the number of newly infected
        patch[p].n_newly_infected_total++;
        patch[p].n_newly_infected_total_by_risk[indiv->sex_risk]++;
        
        // PConly outputs
        if((t - indiv->DoB) > AGE_PC_MIN && (t - indiv->DoB) < (AGE_PC_MAX + 1)){
            patch[p].n_newly_infected_total_pconly++;
            patch[p].n_newly_infected_total_by_risk_pconly[indiv->sex_risk]++;
        }

        // SPVL_infector stores the SPVL of the person doing the infecting. 
        // As seeded infection set to -1 to indicate that this is a seeded infection. */
        indiv->SPVL_infector = -1;

        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Finished new_infection()\n");
            fflush(stdout);
        }
        // Update lists of serodiscorant partnerships following a seroconversion.  
        inform_partners_of_seroconversion_and_update_list_serodiscordant_partnerships(indiv,
            overall_partnerships->susceptible_in_serodiscordant_partnership,
            overall_partnerships->n_susceptible_in_serodiscordant_partnership);
        
        // Phylogenetics transmission outputs for seeded cases
        // (these are dummy outputs signifying that they are seed infections): */
        if(WRITE_PHYLOGENETICS_OUTPUT == 1){
        store_phylogenetic_transmission_initial_cases(output, patch[p].param, 
            indiv, file_data_store, t0, t_step);
        }
        patch[p].DEBUG_NHIVPOS++;
        
        // Increment incident infections
        int age = (int) floor(t - indiv->DoB);
        int g = indiv->gender;
        if(age >= MAX_AGE){
            age = MAX_AGE;
        }
        int age_idx = FIND_AGE_GROUPS_UNPD[age - AGE_ADULT];
        int year_idx = (int) floor(t - patch[p].param->start_time_simul);
        
        patch[p].calendar_outputs->N_calendar_infections[g][age_idx][year_idx]++;
    }
}


/* Function: next_hiv_event()
Determine the next HIV-related event for an individual given their current state (CD4, SPVL, etc). 

Events are CD4 progression, AIDS death, or starting ART because of CD4 low (bypassing testing -
these are people who have either never tested but turn up to clinic with AIDS symptoms, or who
know they are HIV+ but previously left the care cascade. It also potentially increases the rate
at which someone in the cascade can start ART when their CD4 gets low.


Arguments
---------
indiv : pointer to an individual structure
    The individual in question (for whom these events are to be scheduled).  
hiv_pos_progression : pointer to a multi-dimensional array of individuals
n_hiv_pos_progression : pointer to long
size_hiv_pos_progression : pointer to long
param : pointer to a parameters structure
    Parameters structure of the simulation
t : double
    Current time in years
cumulative_outputs : pointer to a cumulative_outputs_struct
    Structure that records the cumulative number of different events (such as CD4, HIV tests)

Returns
-------
Nothing;
*/

void next_hiv_event(individual* indiv, individual ***hiv_pos_progression, 
    long *n_hiv_pos_progression, long *size_hiv_pos_progression, parameters *param, double t,
    cumulative_outputs_struct *cumulative_outputs, calendar_outputs_struct *calendar_outputs){
    
    double time_to_next_event;
    int t_step_event;
    
    // Find the current year as an index
    int year_idx = (int) floor(t) - param->start_time_simul;
    
    if((indiv->cd4 < 0) || (indiv->cd4 > 3)){
        printf("ERROR: Unrecognised cd4  %i %li\n",indiv->cd4,indiv->id);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    /* Note that if indiv->idx_hiv_pos_progression[0]<0 then the event is either 
    "no event scheduled" (-1) or "an event was scheduled but only for after the end of the
    simulation" (-2).  What can happen is that an individual has an HIV progression event scheduled
    after the end of the simulation, but then they start ART and become virally unsuppressed. At
    that point they get a new HIV event but the old event is set to -1 (it was -2 but gets reset in
    the restart HIV progression procedure).
     */
    if(
        (indiv->debug_last_hiv_event_index == indiv->idx_hiv_pos_progression[0]) &&
        (indiv->idx_hiv_pos_progression[0] >= 0)
    ){
        printf("ERROR - trying to schedule a new hiv event (type=%d)", indiv->next_HIV_event);
        printf("in next_hiv_event() that occurs at the same time as the previous event for person");
        printf("%ld in patch %d at time = %6.4f. Exiting\n",indiv->id, indiv->patch_no, t);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Individual %ld in patch %d time = %6.4f in next_hiv_event()\n",
            indiv->id, indiv->patch_no, t);
        fflush(stdout);
    }
    
    /* We first calculate the mean time to next event, adjusting for SPVL, VS/VU, etc. */
    double mean_time_to_next_event;
    mean_time_to_next_event = get_mean_time_hiv_progression(param, indiv);

    /* If on ART but not virally suppressed, allow CD4 progression to be slower. */
    if (indiv->ART_status == LTART_VU){
        mean_time_to_next_event =
            mean_time_to_next_event * param->factor_for_slower_progression_ART_VU;
    }
    /* Now draw time to next CD4 stage: */
    time_to_next_event = gsl_ran_exponential(rng, mean_time_to_next_event);
    
    /* Ensure we never have an event in the current timestep. */
    if(time_to_next_event < TIME_STEP){
        time_to_next_event = TIME_STEP;
    }
    
    indiv->PANGEA_t_prev_cd4stage = indiv->PANGEA_t_next_cd4stage;
    indiv->PANGEA_t_next_cd4stage = t + time_to_next_event;
    
    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Individual %ld in patch %d ", indiv->id, indiv->patch_no);
        printf("time = %6.4f has cd4=%i", t, indiv->cd4);
        printf("Time to next event = %6.4f in function next_hiv_event() drawn from exp(%lf)\n",
            time_to_next_event,mean_time_to_next_event);
        fflush(stdout);
    }

    /* Note: hiv event progression is mostly linear. The only branching possibility is that 
    if currently CD4<200, allow the person to start ART as an emergency rather than dying before
    starting ART.  Note that we are assuming they are not already on ART. */
    if(indiv->cd4 == 3){
        if(
            (t >= param->COUNTRY_ART_START) && 
            (indiv->ART_status != LTART_VU) && 
            (DO_HIV_TESTING == 1) && 
            (ALLOW_EMERGENCY_ART == 1)
        ){
            double time_emergency_start_ART = get_time_emergency_start_ART(indiv,param,t);
            
            /* Decide if individual will start ART because of low CD4 (and symptoms) - at which
            point they may die quickly on ART - or die without starting ART: */
            
            if (time_emergency_start_ART < time_to_next_event){
                
                time_to_next_event = time_emergency_start_ART; 
                // Assume always a non-popart test.
                cumulative_outputs->N_total_CD4_tests_nonpopart++;
                calendar_outputs->N_calendar_CD4_tests_nonpopart[year_idx]++;
                indiv->next_HIV_event = HIVEVENT_STARTART_LOWCD4;
                
                if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
                    printf("Individual %ld in patch %d ",indiv->id, indiv->patch_no);
                    printf("time = %6.4f is scheduled to start emergency ART ", t);
                    printf("at time = %6.4f in function next_hiv_event()\n",t + time_to_next_event);
                    fflush(stdout);
                }
            }
            else{
                indiv->next_HIV_event = HIVEVENT_AIDSDEATH;
            }
        }else{
            /* If ART is not available yet, or if the person is on ART but virally unsuppressed, 
            and CD4<200, then their next HIV event is AIDS death. */
            indiv->next_HIV_event = HIVEVENT_AIDSDEATH;
        }
    }else{
        /* Next event is CD4 progression: */
        indiv->next_HIV_event = HIVEVENT_CD4_PROGRESSION;
    }

    indiv->DEBUGTOTALTIMEHIVPOS += time_to_next_event;  
    
    /* This gets the index for hiv_pos_progression
    (and n_hiv_pos_progression and size_hiv_pos_progression) arrays. */
    t_step_event = (int) round(((t - param->start_time_hiv) + 
        time_to_next_event) * N_TIME_STEP_PER_YEAR); 
    
    /* For debugging: */
    if(t_step_event == indiv->idx_hiv_pos_progression[0]){
        printf("ERROR - trying to schedule a new hiv event in next_hiv_event() that occurs ");
        printf("at the same time as the current event for person %ld ", indiv->id);
        printf("in patch %d at time = %6.4f. Exiting\n",indiv->patch_no, t);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    /* Only add an event if this happens before the end of the simulation: */
    if(t_step_event <= (param->end_time_simul - param->start_time_hiv) * N_TIME_STEP_PER_YEAR){
        
        indiv->debug_last_hiv_event_index = indiv->idx_hiv_pos_progression[0];
        indiv->idx_hiv_pos_progression[0] = t_step_event;
        indiv->idx_hiv_pos_progression[1] = n_hiv_pos_progression[t_step_event];
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Adding HIV progression event for %ld in patch %d ", indiv->id, indiv->patch_no);
            printf("at time = %6.4f i=%i. Indiv indices %ld %ld\n",
                t + time_to_next_event, t_step_event,
                indiv->idx_hiv_pos_progression[0], indiv->idx_hiv_pos_progression[1]);
            
            fflush(stdout);
        }

        /* Check if we've run out of memory: */
        if(n_hiv_pos_progression[t_step_event] >= (size_hiv_pos_progression[t_step_event])){
            
            /* Note that realloc does not work (we need to pass a pointer to the pointer, which is
            really complicated as it propagates through several functions (so maybe make
            planned_breakups[time_breakup] **), so ditch this code for now and use the following
            lines: */
            printf("Unable to re-allocate hiv_pos_progression[i]. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        hiv_pos_progression[t_step_event][n_hiv_pos_progression[t_step_event]] = indiv;
        n_hiv_pos_progression[t_step_event]++;

    }else{
        /* Next event happens after end of simulation so set to no event. */
        indiv->idx_hiv_pos_progression[0] = EVENTAFTERENDSIMUL;
        indiv->idx_hiv_pos_progression[1] = -1;

        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Not scheduling HIV progression event for %ld ", indiv->id);
            printf("in patch %d as after end of simulation at t=%6.4f\n",
                indiv->patch_no, t + time_to_next_event);
            fflush(stdout);
        }
    }
    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Next HIV progression event for %ld in patch %d ", indiv->id, indiv->patch_no);
        printf("is scheduled to be %i at t=%6.4f. Indiv indices are %li %li\n",
            indiv->next_HIV_event, t + time_to_next_event, 
            indiv->idx_hiv_pos_progression[0], indiv->idx_hiv_pos_progression[1]);
        
        fflush(stdout);
    }
}


/* Function: carry_out_HIV_events_per_timestep()
Go through the list of scheduled HIV events for HIV+ people for this timestep

For each person to whom some HIV event happens, draw their next HIV-based event (via 
indiv->next_HIV_event) unless they die from AIDS at this timestep.  The list of scheduled HIV
events are stored in the array patch[p].hiv_pos_progression and the list of possible events are
defined in constants.h.  This function checks the individual's attribute next_HIV_event and performs
the following activities based upon what this attribute is:

HIVEVENT_AIDSDEATH
    Remove the individual from list of care cascade events.  Call the fn individual_death_AIDS().  
HIVEVENT_ENDOFACUTE
    Change the individual's HIV_status attribute to CHRONIC.  
HIVEVENT_CD4_PROGRESSION
    Increment an individual's cd4 attribute by 1 (progress their CD4 category).  
HIVEVENT_STARTART_LOWCD4
    Start emergency ART.  Count the number of individuals transitioning to emergency ART from 
    other stages of the cascade.  Call start_ART_process().  
NOEVENT
    Nothing happens.  

Finally, next_HIV_event() is called to determine the next event for an individual.  

Arguments
---------

t : double
    Year in question
patch : pointer to a patch_struct structure
p : int
    Patch number
overall_partnerships : pointer to an all_partnerships structure
debug : pointer to a debug_struct structure
file_data_store : pointer to a file_struct structure

Returns
-------
Nothing; carries out HIV events for a particular timestep for individuals with a scheduled event and
schedules new events for those individuals.  Various lists and individual attributes are updated.  
*/

void carry_out_HIV_events_per_timestep(double t, patch_struct *patch, int p, 
    all_partnerships *overall_partnerships, debug_struct *debug, file_struct *file_data_store){

    int array_index_for_hiv_event = 
        (int) round((t - patch[p].param->start_time_hiv) * N_TIME_STEP_PER_YEAR);
    
    int n, n_events = patch[p].n_hiv_pos_progression[array_index_for_hiv_event];
    individual *indiv;
    
    // Loop through all HIV events this timestep
    for(n = 0; n < n_events; n++){
        // Assign pointer to the individual in question
        indiv = patch[p].hiv_pos_progression[array_index_for_hiv_event][n];
        
        if(indiv->cd4 == DUMMYVALUE){
            printf("Error. Using uninitialised person.\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Calling carry_out_HIV_events_per_timestep() for ");
            printf("%ld in patch %d at t=%6.4f, index = %i %i. Indiv indices are %li %li\n",
                indiv->id, indiv->patch_no, t, array_index_for_hiv_event, n,
                indiv->idx_hiv_pos_progression[0], indiv->idx_hiv_pos_progression[1]);
            fflush(stdout);
        }
        
        if(indiv->cd4 == DEAD){
            if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
                printf("Individual %ld from patch %d is dead, ", indiv->id, indiv->patch_no);
                printf("so not carrying out further HIV events for this person at t=%6.4f\n",t);
                fflush(stdout);
            }
            // Move to the next person. Note - we can set up a similar procedure to other lists to
            // remove this person from this list but it is not necessary. As things stand, no hiv
            // event happens to the dead person and no new event is scheduled for them.
            continue;
        }
        
        // If next event is AIDS death
        if(indiv->next_HIV_event == HIVEVENT_AIDSDEATH){
            
            if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
                printf("Individual %ld from patch %d dies from AIDS at t=%6.4f\n", 
                    indiv->id, indiv->patch_no, t);
                fflush(stdout);
            }
            
            // Person died so add 1 to death counter
            patch[p].DEBUG_NDEATHS = patch[p].DEBUG_NDEATHS + 1;
            
            // Remove individual from cascade events
            remove_from_cascade_events(indiv, patch[p].cascade_events, patch[p].n_cascade_events,
                patch[p].size_cascade_events,t, patch[p].param);
            
            // Function removes person from everything except the cascade event list:
            individual_death_AIDS(patch[p].age_list, indiv, patch[p].n_population,
                patch[p].n_population_oneyearagegroups, patch[p].n_infected,
                patch[p].n_population_stratified, t, patch[p].param,
                overall_partnerships->susceptible_in_serodiscordant_partnership,
                overall_partnerships->n_susceptible_in_serodiscordant_partnership,
                overall_partnerships->pop_available_partners,
                overall_partnerships->n_pop_available_partners, 
                patch[p].cascade_events, patch[p].n_cascade_events, 
                patch[p].size_cascade_events, patch, p, file_data_store);
        }else{
            // If the next event is not AIDS death then check other events that may occur
            // If next event is end of acute phase
            if(indiv->next_HIV_event == HIVEVENT_ENDOFACUTE){
                
                // Change individual's HIV_status to CHRONIC infection
                // Note that CD4 category has already been assigned upon seroconversion.
                indiv->HIV_status = CHRONIC;
                
                if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
                    printf("Individual %ld from patch %d is ", indiv->id, indiv->patch_no);
                    printf("progressing from acute to chronic at t=%6.4f\n",t);
                    fflush(stdout);
                }
            
            // If next event is CD4 progression
            }else if(indiv->next_HIV_event == HIVEVENT_CD4_PROGRESSION){
                
                if(indiv->cd4 < 0){
                    printf("Error. Trying to advance CD4 state=%i", indiv->cd4);
                    printf("for id=%ld in patch %d ", indiv->id, indiv->patch_no);
                    printf("when negative (-1 = uninfected, -2=dead). Exiting!\n");
                    printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                    fflush(stdout);
                    exit(1);
                }
                
                if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
                    printf("Individual %ld from patch %d ", indiv->id, indiv->patch_no);
                    printf("is changing CD4 category from %i to %i at t=%6.4f\n", 
                        indiv->cd4, indiv->cd4 + 1, t);
                    fflush(stdout);
                }
                
                // Increase CD4 category by 1
                indiv->cd4++;
                
            // If next event is emergency ART due to low CD4
            }else if(indiv->next_HIV_event == HIVEVENT_STARTART_LOWCD4){
                
                if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
                    printf("A: Starting emergency ART for adult %ld ", indiv->id);
                    printf("in patch %d at t=%6.4f. Index of HIV events = %li %li\n",
                        indiv->patch_no, t, indiv->idx_hiv_pos_progression[0],
                        indiv->idx_hiv_pos_progression[1]);
                    fflush(stdout);
                }
                
                // Assume that HIV testing still occurs at this point
                indiv->time_last_hiv_test = t;

                // PANGEA stuff - get CD4 at time of emergency ART:
                indiv->PANGEA_t_diag = t;
                indiv->PANGEA_cd4atdiagnosis = PANGEA_get_cd4(indiv, t);

                // Assume that HIV test happens at same time as start of ART:
                if(indiv->PANGEA_date_firstARTstart < 0){
                    indiv->PANGEA_date_firstARTstart = t;
                    indiv->PANGEA_cd4atfirstART = PANGEA_get_cd4(indiv, t);
                }
                
                // Update debug counters regarding transitions to emergency ART
                if(indiv->ART_status == ARTNEG){
                    debug->art_vars[p].n_start_emergency_art_fromartnaive++;
                }else if(indiv->ART_status == ARTNAIVE){
                    debug->art_vars[p].n_start_emergency_art_fromartnaive++;
                }else if(indiv->ART_status == ARTDROPOUT){
                    debug->art_vars[p].n_start_emergency_art_fromartdroupout++;
                }else if(indiv->ART_status == CASCADEDROPOUT){
                    debug->art_vars[p].n_start_emergency_art_fromcascadedropout++;
                }
                // Count totals to make sure have captured all possible ways to start emergency ART
                debug->art_vars[p].n_start_emergency_art++;

                // Start indiv on ART at low CD4.  
                // The second-to-last parameter in start_ART_process() states that this 
                // individual started ART because of AIDS symptoms (set to 1).  
                start_ART_process(indiv,patch[p].param, t, patch[p].cascade_events,
                    patch[p].n_cascade_events, patch[p].size_cascade_events,
                    patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression,
                    patch[p].size_hiv_pos_progression,1, file_data_store,
                    patch[p].calendar_outputs);
            }else{
                printf("ERROR: Unknown HIV event %i for ", indiv->next_HIV_event);
                printf("id=%li in patch %d with indices %i %i %li %li. Exiting.\n",
                    indiv->id, indiv->patch_no, array_index_for_hiv_event, n,
                    indiv->idx_hiv_pos_progression[0], indiv->idx_hiv_pos_progression[1]);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }
            
            // Schedule next HIV event. If they start ART this is dealt with in 
            // cascade_events (including AIDS death), so do not schedule a new HIV-related event
            // unless they stop ART. Note that start_ART_process sets indiv->next_HIV_event to
            // NOEVENT. 
            if(indiv->next_HIV_event != NOEVENT){
                if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
                    printf("Individual %ld in patch %d ", indiv->id, indiv->patch_no);
                    printf("time = %6.4f going from carry_out_HIV_events_per_timestep()", t);
                    printf(" to next_hiv_event()\n");
                    fflush(stdout);
                }
                
                next_hiv_event(indiv, patch[p].hiv_pos_progression, 
                    patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression,
                    patch[p].param, t, patch[p].cumulative_outputs, patch[p].calendar_outputs);
            }
        }
    }
}


/* Function: get_window_result()
Adjust HIV test to allow for the fact that antibody/antigen tests do not work for first few weeks.

Individuals who test during very early infection may not get a positive test result.  The window
period for the test is taken to be 60 days until 2006 (or 2010 in South Africa), roughly the window
period for first or second-generation HIV antibody tests, and 30 days from then onwards, with the
latter corresponding to the window period of more recent antigen/antibody tests.

This function returns 0 (meaning person tests negative) or 1 (meaning tests positive). 

See: M. S. Cohen, C. L. Gay, M. P. Busch, and F. M. Hecht. The detection of acute HIV infection. J. Infect. Dis., 202 Suppl 2:S270–277, Oct 2010.

Arguments
---------
time_since_sc : double
    Time since seroconversion
t : double
    Current year
patch : pointer to a patch_struct structure

Returns
-------
int of negative test result (0) or positive test result (1).  
*/

int get_window_result(double time_since_sc, double t, patch_struct *patch){
    
    if(patch->community_id <= IS_ZAMBIA){
        if(t <= 2006){
            if(time_since_sc >= 60.0/365)
                return 1;
            else
                return 0;
        }else{
            if(time_since_sc >= 30.0/365)
                return 1;
            else
                return 0;
        }
    }else{
        if(t <= 2010){
            if(time_since_sc >= 60.0/365)
                return 1;
            else
                return 0;
        }else{
            if(time_since_sc >= 30.0/365)
                return 1;
            else
                return 0;
        }
    }
}


/* Function: joins_preart_care()
Determine if a person will join pre-ART care or dropout as a function of time, their CD4, etc. 

Checks if an individual is in PopART or not (by checking the indiv->next_cascade_event attribute) 
and then determines if an individual collects their HIV test results.  For background care cascade
if the individual collects their HIV test results (the probability of which is stratifed by the 
individual's CD4 category) then it is assumed an individual has a CD4 test.  This function returns 
the result of whether that individual collects their CD4 test results (1) or drops out of care (0).
For individuals in PopART it is assumed that everyone collects their HIV test results.  This 
function then also returns whether individual's collect their CD4 test results (1) or drop out of 
care (0).  Note that the probability of collecting CD4 test results is really "1- Pr(do not stay in
care or start ART)".  

Arguments
---------
indiv : pointer to an individual structure
    Individual in question that might join pre-ART care or drop out
param : pointer to a parameters structure
    Structure of simulation parameters
t : double
    Current time (year)
cumulative_outputs : pointer to a cumulative_outputs_struct structure
    Structure containing cumulative outputs (for counting CD4 tests etc).  

Returns
-------
integer : 0 means the individual drops out, 1 means the individual stays in pre-ART care
*/

int joins_preart_care(individual* indiv, parameters *param, double t, 
    cumulative_outputs_struct *cumulative_outputs, calendar_outputs_struct *calendar_outputs){
        
    /* This tells us if the cd4 test is due to PopART (is_popart=1) or not (is_popart=0). */
    int is_popart = (indiv->next_cascade_event >= NCASCADEEVENTS_NONPOPART);
    int year_idx = (int) floor(t) - param->start_time_simul;
    double p_collects_hiv_test; 
    
    // Check if this person is part of the background care cascade or not
    if(is_popart == NOTPOPART){
        
        // Divide individuals into those with CD4>200 and those with CD4<=200
        if(indiv->cd4 < 3){
            p_collects_hiv_test = param->p_collect_hiv_test_results_cd4_over200;
        }else{
            p_collects_hiv_test = param->p_collect_hiv_test_results_cd4_under200;
        }
        
        // Check if the individual collects HIV test results
        if(gsl_ran_bernoulli(rng, p_collects_hiv_test) == 1){
            
            // If an individual collects their HIV test results then they also 
            // have their first CD4 test
            cumulative_outputs->N_total_CD4_tests_nonpopart++;
            calendar_outputs->N_calendar_CD4_tests_nonpopart[year_idx]++;
            
            // returns whether they drop out (0) or stays in cascade (1)
            return gsl_ran_bernoulli(rng, param->p_collect_cd4_test_results_cd4_nonpopart);
        }else{
            // If the individual does not collect HIV test results then assume they
            // drop out of the care cascade (return 0)
            return 0;
        }
    }else{
        // Individual is in PopART.  
        // Assume everyone in PopART gets their HIV test and has a CD4 test
        cumulative_outputs->N_total_CD4_tests_popart++;
        calendar_outputs->N_calendar_CD4_tests_popart[year_idx]++;
        
        // Divide probability of collecting CD4 results into CHiPs round 1 and later rounds. 
        // CHIPS_START_YEAR[1] is the start of round 2 (CHIPS_START_YEAR[0] is the start of R1):
        if(t < param->CHIPS_START_YEAR[1]){
            // returns whether they drop out (0) or stays in cascade (1)
            return gsl_ran_bernoulli(rng,param->p_collect_cd4_test_results_cd4_popartYEAR1);
        }else{
            // returns whether they drop out (0) or stays in cascade (1)
            return gsl_ran_bernoulli(rng,param->p_collect_cd4_test_results_cd4_popartYEAR2onwards);
        }
    }
}


/* Function: remains_in_cascade()
Test that an HIV+ person remains in the background care cascade until the next CD4 test when they've
just had a CD4 test and were not eligible. 

For PopART p_collect_cd4_test_results_cd4_popartYEAR1 and
p_collect_cd4_test_results_cd4_popartYEAR2onwards includes dropout at this step already so
don't double-count. 

Arguments
---------
indiv : pointer to an individual structure
    Pointer to the structure of the individual in question
param: pointer to parameters structure
    Parameters of the simulation
is_popart : int
    Indicator of whether this individual is in the PopART (1) community or not (0)
t : double
    Current time in decimal years

Returns
-------
integer: whether individual drops out (0) or remains (1) in cascade
*/

int remains_in_cascade(individual* indiv, parameters *param, int is_popart){
    
    // Assume individuals in PopART always remain in the cascade
    if(is_popart == POPART){
        return 1;
    }
    // Gets CD4 tested, determines of individual drops out (0) or stays in the cascade (1)
    return gsl_ran_bernoulli(rng, param->p_collect_cd4_test_results_cd4_nonpopart);
}


/* Function: measured_cd4_cat()
Given a 'real' CD4 category, draws the measured CD4

The parameter structure stores a multi-dimensional array of misclassification probabilities 
within the array `cumulative_p_misclassify_cd4`.  These are cumulative probabilities of 
misclassification of a 'real' CD4 category into any other CD4 category.  Such rates of 
misclassification could be based upon those in Cori et al (2015) AIDS (table 1) which is 
based up on data from the ATHENA cohort.  In PopART pre-unblinding projections these 
misclassification probabilities were all set to 0 (i.e. no chance of misclassification).  

Arguments
---------
param : pointer to a parameters structure
    Parameter structure should have a multi-dimensional array called 
    cumulative_p_misclassify_cd4 which gives the cumulative probabilities of misclassifying 
    a CD4 category into any other category.  
real_cd4_cat : int
    The real CD4 category of the individual in question.  

Returns
-------
The measured CD4 cat (int).  
*/ 

int measured_cd4_cat(parameters *param, int real_cd4_cat){
    
    double x = gsl_rng_uniform(rng);
    int i = 0;
    
    while(
        (i < NCD4) && 
        (param->cumulative_p_misclassify_cd4[real_cd4_cat][i] < x)
    ){
        // The final value of i will be the measured CD4 cat
        i++;
    }
    return i;
}


/* Function: art_cd4_eligibility_group()
Return the highest CD4 group that is eligible for ART at time t in a given setting.

Function returns the index corresponding to that group. E.g. if the function returns 2, then
CD4=2 and CD4=3 are eligible for ART - in other words CD4<350.

Arguments
---------
param : pointer to a parameters struct
    Parameter structure with information on changes in a country's ART guidelines.  
t : double
    Year in question (can be decimal)

Returns
-------
Integer giving the highest CD4 group that is eligible for ART at current time.  
*/

int art_cd4_eligibility_group(parameters *param, double t){

    if(t < param->COUNTRY_ART_START){
        return 4; /* Nobody eligible before start of ART. */
    }else if(t < param->COUNTRY_CD4_350_START){
        return 3;
    }else if(t < param->COUNTRY_CD4_500_START){
        return 2;
    }else if(t < param->COUNTRY_IMMEDIATE_ART_START){
        return 1;
    }else{
        return 0;
    }
}


/* Function: is_eligible_for_art()
Determin if individual is eligible for ART
 
Returned value depends on time at which test takes place (implicitly, the country in question),
may include misclassification of CD4 categories (see function measured_cd4_cat()), and depends
on whether this is a PopART arm A community (all HIV+ individuals are eligible if so).  

Arguments
---------
indiv : pointer to an individual structure
    Individual in question having their eligibility for ART tested.  
param : pointer to a parameters structure
    Structure holding all parameters
t : double
    Current time (in years)
patch : pointer to a patch_struct structure
    Patch structure with patch-related simulation parameters
p : int
    Patch number

Returns
-------
Integer denoting eligibility of an individual for ART; eligible (1) or non-eligible (0)  
*/

int is_eligible_for_art(individual* indiv, parameters *param, double t, patch_struct *patch, int p){
    
    /* If before ART in country then end. */
    if(t <= param->COUNTRY_ART_START){
        return 0;
    }

    // If Popart has started (ie CHiPs round 1 has started) and we are in arm A then any HIV+ 
    // individual is eligible, regardless of CD4 category.  
    if(
        (t >= (patch[p].param->CHIPS_START_YEAR[0] +
            patch[p].param->CHIPS_START_TIMESTEP[0]*TIME_STEP)) && 
        (patch[p].trial_arm == ARM_A)
    ){
        return 1;
    }
    
    // Draw CD4 category
    int measuredcd4 = measured_cd4_cat(param, indiv->cd4);
    
    // If the measured CD4 category is greater than or equal to the eligibility category then 
    // return 1, otherwise 0.  
    if(measuredcd4 >= art_cd4_eligibility_group(param, t)){
        return 1;
    }else{
        return 0;
    }
}


/* Function: get_time_emergency_start_ART()
Return length of time it takes to start ART once an individual's CD4 goes below 200.

At present the time is drawn from an identical distribution to that of the time it takes to go from
CD4 of 200 to death.  Thus 50% of people will die before starting emergency ART.  Note that an 
individual's CD4 category (i.e. indiv->cd4) should always be 3 within this function.  

Arguments
---------
indiv : pointer to an individual structure
    Individual in question (for whom a time to start emergency ART needs to be drawn).  
param : pointer to a parameters structure
    Parameters of the simulation
t : double
    Time in question (year in decimal)

Returns
-------
double : the duration of time until an individual starts emergency ART (or end time of the
simulation) if this time is beyond the end of the simulation.  
*/ 

double get_time_emergency_start_ART(individual* indiv, parameters *param, double t){
    
    double mean_time_to_next_event, time_to_start_emergency_art;
    
    // Calculate the mean time to next event, adjusting for SPVL, VS/VU, etc.
    mean_time_to_next_event = get_mean_time_hiv_progression(param, indiv);
    
    // If on ART but not virally suppressed, CD4 progression will be slower.
    if (indiv->ART_status == LTART_VU){
        mean_time_to_next_event =
            mean_time_to_next_event * param->factor_for_slower_progression_ART_VU;
    }

    // Draw time to start emergency ART from exponential distribution
    time_to_start_emergency_art = gsl_ran_exponential(rng, mean_time_to_next_event);
    
    if(time_to_start_emergency_art < TIME_STEP){
        time_to_start_emergency_art = TIME_STEP;
    }
    
    if(time_to_start_emergency_art < (param->end_time_simul - t)){
        return time_to_start_emergency_art;
    }else{
        return param->end_time_simul;
    }
}


/* Function: start_ART_process()
Set up state when beginning ART (including emergency ART).  

The crux of this function simply changes the `next_cascade_event` attribute of the 
individual in question and determines the time until this event (`time_to_next_event`).  These
are then passed to schedule_generic_cascade_event() which determines the next event.  

This state represents the first few months of ART when the individual is becoming virally
suppressed and may be more likely to drop out or die. 

This function also determines what happens in the next cascade step after someone starts ART 
(i.e. viral suppression, on ART but not virally suppressed, drops out, dies due to AIDS-related
illness).

Note: If is_emergency == 1, this means they started because of AIDS symptoms at low CD4 (so can
increase mortality rate)

No further CD4 progression (or other hiv_pos_progression[] events happens to them when they
start ART. Note that AIDS death is handled in cascade_events in this state (representing things
such as immune reconstitution inflammatory syndrome, where AIDS mortality is potentially
higher).  By comparison for somoneone not on ART they can die from AIDS through events
scheduled in hiv_pos_progression[]. 

Arguments
---------
indiv : pointer to an individual structure
    Individual in question that may be starting the ART process
param : pointer to a parameters structure
t : double
    Year in question (in decimal)
cascade_events : pointer to multidimenional array of individual structures
n_cascade_events : pointer to array of long
size_cascade_events : pointer to array of long
hiv_pos_progression: multidimensional array of pointers to individual structures
n_hiv_pos_progression : array of pointer to longs
size_hiv_pos_progression : array of pointer to long
is_emergency : int
    Is this emergency (1) ART or not (0).  
file_data_store : pointer to a file_struct structure.  

Returns
-------
Nothing; adjusts attributions of the individual in question and other lists of events if certain
processes are triggered.  
*/

void start_ART_process(individual* indiv, parameters *param, double t, 
    individual ***cascade_events, long *n_cascade_events, long *size_cascade_events, 
    individual ***hiv_pos_progression, long *n_hiv_pos_progression, long *size_hiv_pos_progression,
    int is_emergency, file_struct *file_data_store, calendar_outputs_struct *calendar_outputs){
    
    /* index for current time in this array: hiv_pos_progression, only used for debugging */
    int array_index_for_hiv_event = (int) round((t - param->start_time_hiv)*N_TIME_STEP_PER_YEAR);
    
    // Check if we are following a particular individual
    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Removing %li from HIV pos progression as starting ART. ", indiv->id);
        printf("Old (to be removed) array event = %li %li. Current index = %i \n", 
            indiv->idx_hiv_pos_progression[0],
            indiv->idx_hiv_pos_progression[1],
            array_index_for_hiv_event);
        fflush(stdout);
    }
    
    // Is this emergency ART?  
    if(is_emergency == 0){
        /* Note the final '2' argument means that the person is starting ART, not dying. */
        remove_from_hiv_pos_progression(indiv,  hiv_pos_progression, 
            n_hiv_pos_progression, size_hiv_pos_progression,t, param, 2);
    }else{
        /* The final '4' denotes that this is emergency ART. */
        remove_from_hiv_pos_progression(indiv,  hiv_pos_progression, 
            n_hiv_pos_progression, size_hiv_pos_progression,t, param, 4);
    }

    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Finished removing %li from HIV pos progression as starting ART. ", indiv->id);
        printf("Old (to be removed) array event = %li %li \n",
            indiv->idx_hiv_pos_progression[0],
            indiv->idx_hiv_pos_progression[1]);
        fflush(stdout);
    }
    
    indiv->next_HIV_event = NOEVENT;
    indiv->idx_hiv_pos_progression[0] = -1;
    indiv->idx_hiv_pos_progression[1] = -1;

    // Individuals who start emergency ART due to low CD4 need to have any cascade events removed
    if(is_emergency == 1){
        remove_from_cascade_events(indiv, cascade_events, n_cascade_events, 
            size_cascade_events, t, param);
        indiv->next_cascade_event = NOEVENT;
        indiv->idx_cascade_event[0] = -1;
        indiv->idx_cascade_event[1] = -1;
    }
    
    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Finished removing %li from HIV pos progression as starting ART. ", indiv->id);
        printf("Removed array event type = %i indices %li %li \n",
            indiv->next_HIV_event,
            indiv->idx_hiv_pos_progression[0],
            indiv->idx_hiv_pos_progression[1]);
        fflush(stdout);
    }
    
    /* Output time person was seropositive if HIV+ and not on ART.
    The final argument is reason for being removed from survival cohort:
    1="AIDS death", 2="death from natural causes", 3="start ART". 
    Note we don't bother with the end of the simulation for now. */
    if(WRITE_DEBUG_HIV_DURATION_KM == 1){
        /* Only consider people who are ART-naive for now: */
        if (indiv->ART_status == ARTNAIVE || indiv->ART_status == ARTNEG){
            write_hiv_duration_km(indiv, t, file_data_store, 3);
        }
    }
    
    indiv->ART_status = EARLYART;
    
    // Count the number of ART initiations
    int year_idx = (int) floor(t) - param->start_time_simul;
    calendar_outputs->N_calendar_started_ART_annual[year_idx]++;
    
    // Record the ART initiation for the TREATS output
    if(WRITE_TREATS_OUTPUT){
        int age = (int) floor(t - indiv->DoB);
        if(age >= MAX_AGE){
            age = MAX_AGE;
        }
    
        int age_idx = FIND_AGE_GROUPS_UNPD[age - AGE_ADULT];
        int cd4 = indiv->cd4;
        int spvl = indiv->SPVL_cat;
    
        // For uninfected individuals (SPVL_cat == -1), make sure indexing is correct
        spvl += 1;
        if(spvl < 0){
            printf("Error, SPVL_cat < 0\n");
            fflush(stdout);
            exit(1);
        }
    
        // HIV uninfected inidividuals (CD4 == -1) have index of 0.  
        cd4 += 1;
        if(cd4 < 0){
            printf("Error, cd4 category < 0\n");
            fflush(stdout);
            exit(1);
        }
        calendar_outputs->N_calendar_started_ART[indiv->gender][age_idx][cd4][spvl][year_idx]++;
    }
    
    /* This tells us if the cd4 test is due to PopART (is_popart=1) or not (is_popart=0). */
    int is_popart = (indiv->next_cascade_event >= NCASCADEEVENTS_NONPOPART);

    /* Now schedule what the next cascade event will be for them: */
    double x = gsl_rng_uniform (rng);
    double p_dropout, time_to_next_event;

    /* Determine probability individual drops out of ART care: */
    if(indiv->cd4 < 3){ /* check if CD4>200 */
        p_dropout = param->p_leaves_earlyart_cd4_over200_if_not_die_early*
            (1.0 - param->p_dies_earlyart_cd4[indiv->cd4]);
    }else{
        p_dropout = param->p_leaves_earlyart_cd4_under200_if_not_die_early*
            (1.0 - param->p_dies_earlyart_cd4[indiv->cd4]);
    }
    
    // Determine if individual drops out of the cascade
    // Note: Could rearrange if statement to be more efficient by putting highest prob event first
    // ( i.e. continuing on ART and becoming VS as first check in if statement)
    if(x < p_dropout){
        if(is_popart == NOTPOPART){
            indiv->next_cascade_event = CASCADEEVENT_DROPOUT_NONPOPART;         
            time_to_next_event = t + param->t_earlyart_dropout_min[NOTPOPART] + 
                gsl_rng_uniform (rng)*param->t_earlyart_dropout_range[NOTPOPART]+TIME_STEP;
        }else{
            indiv->next_cascade_event = CASCADEEVENT_DROPOUT_POPART;
            // assumed takes 0-6 months!!!
            time_to_next_event = t + param->t_earlyart_dropout_min[POPART] + 
                gsl_rng_uniform (rng)*param->t_earlyart_dropout_range[POPART]+TIME_STEP;
        }
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Scheduling next cascade event after initiating ART for adult %ld: dropout \n",
                indiv->id);
            fflush(stdout);
        }
        
        indiv->DEBUG_cumulative_time_on_ART_early += (time_to_next_event - t);
        
    /* Check if individual dies while on ART */
    }else if(x < (p_dropout + param->p_dies_earlyart_cd4[indiv->cd4])){
        if(is_popart == NOTPOPART){
            indiv->next_cascade_event = CASCADEEVENT_ARTDEATH_NONPOPART;
            time_to_next_event = t + param->t_dies_earlyart_min[NOTPOPART] + 
                gsl_rng_uniform(rng)*param->t_dies_earlyart_range[NOTPOPART]+TIME_STEP;
        }else{
            indiv->next_cascade_event = CASCADEEVENT_ARTDEATH_POPART;
            time_to_next_event = t + param->t_dies_earlyart_min[POPART] + 
                gsl_rng_uniform(rng)*param->t_dies_earlyart_range[POPART]+TIME_STEP;
        }
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Scheduling next cascade event after initiating ART ");
            printf("for adult %ld: dies on ART \n", indiv->id);
            fflush(stdout);
        }
        indiv->DEBUG_cumulative_time_on_ART_early += (time_to_next_event - t);
    
    /* Check if individual becomes virally suppressed on ART. */
    }else if(x < (p_dropout + param->p_dies_earlyart_cd4[indiv->cd4] +
                param->p_becomes_vs_after_earlyart_if_not_die_early_or_leave*
                (1.0 - p_dropout-param->p_dies_earlyart_cd4[indiv->cd4]) )
    ){
        if(is_popart == NOTPOPART){
            indiv->next_cascade_event = CASCADEEVENT_VS_NONPOPART;
            time_to_next_event = t + gsl_rng_uniform (rng)*param->t_end_early_art + TIME_STEP;
        }else{
            indiv->next_cascade_event = CASCADEEVENT_VS_POPART;
            time_to_next_event = t + gsl_rng_uniform (rng)*param->t_end_early_art + TIME_STEP;
        }
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Scheduling next cascade event after initiating ART for adult %ld: VS \n", 
                indiv->id);
            fflush(stdout);
        }
        
        indiv->DEBUG_cumulative_time_on_ART_early += (time_to_next_event - t);
    
    /* Otherwise, individual continues on ART but not fully virally suppressed. */
    }else{
        if(is_popart == NOTPOPART){
            indiv->next_cascade_event = CASCADEEVENT_VU_NONPOPART;
            time_to_next_event = t + gsl_rng_uniform (rng)*param->t_end_early_art + TIME_STEP;
        }else{
            indiv->next_cascade_event = CASCADEEVENT_VU_POPART;
            time_to_next_event = t + gsl_rng_uniform (rng)*param->t_end_early_art + TIME_STEP;
        }
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("Scheduling next cascade event after initiating ART");
            printf(" for adult %ld in patch %i: VU \n", indiv->id, indiv->patch_no);
            fflush(stdout);
        }
    }
    
    schedule_generic_cascade_event(indiv, param, time_to_next_event, 
        cascade_events, n_cascade_events, size_cascade_events,t);

    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Finished scheduling next cascade event after initiating ART ");
        printf("for adult %ld event = %i at array indices %ld %ld\n",
            indiv->id, indiv->next_cascade_event,
            indiv->idx_cascade_event[0], indiv->idx_cascade_event[1]);
        fflush(stdout);
    }
    return;
}


/* Function: draw_initial_hiv_tests()

Schedule HIV tests for individuals in the population.  When HIV testing first begins in the country
(or after PopART ends), if we switch back to clinic-based testing this function is called from
simul.c.  This function goes through every individual currently alive (using the array age_list) 
and schedules an HIV test for each person in the array cascade_events[].  

Arguments
---------
param : pointer to a parameters structure
    The parameters governing the simulation
age_list : pointer to an age_list_struct structure
t : double
    Current time (in decimal years)
cascade_events : multidimensional array of pointers to individual structures
n_cascade_events : pointer to array of long
    Saves the number of events within each slot of the cascade_events array
size_cascade_events : pointer to an array of long

Returns
-------
Nothing; HIV test are scheduled.  
*/

void draw_initial_hiv_tests(parameters *param, age_list_struct *age_list, double t, 
    individual ***cascade_events, long *n_cascade_events, long *size_cascade_events){
    
    int aa, k, g;
    
    for(g = 0; g < N_GENDER; g++){
        // for all but the age group 80+ (which is in a separate part of the age_list struct).
        for(aa = 0; aa < MAX_AGE - AGE_ADULT; aa++){
            // for each individual in that annual age group, schedule their first HIV test:
            for(k = 0; k < age_list->age_list_by_gender[g]->number_per_age_group[aa] ; k++){
                
                schedule_new_hiv_test(age_list->age_list_by_gender[g]->age_group[aa][k], 
                    param, t, cascade_events, n_cascade_events, size_cascade_events);
                
                if(age_list->age_list_by_gender[g]->age_group[aa][k]->id == FOLLOW_INDIVIDUAL &&
                    age_list->age_list_by_gender[g]->age_group[aa][k]->patch_no == FOLLOW_PATCH){
                    printf("Scheduling first HIV test for adult ");
                    printf("%ld ", age_list->age_list_by_gender[g]->age_group[aa][k]->id);
                    printf("(event type %i)at time index %li, array index = %li\n", 
                        age_list->age_list_by_gender[g]->age_group[aa][k]->next_cascade_event,
                        age_list->age_list_by_gender[g]->age_group[aa][k]->idx_cascade_event[0],
                        age_list->age_list_by_gender[g]->age_group[aa][k]->idx_cascade_event[1]);
                    fflush(stdout);
                }
            }
        }
        
        // Schedule HIV tests for each individual in the oldest age group
        for(k = 0; k < age_list->age_list_by_gender[g]->number_oldest_age_group; k++){
            
            schedule_new_hiv_test(age_list->age_list_by_gender[g]->oldest_age_group[k], 
                param, t, cascade_events, n_cascade_events, size_cascade_events);
            
            if(age_list->age_list_by_gender[g]->oldest_age_group[k]->id == FOLLOW_INDIVIDUAL &&
                 age_list->age_list_by_gender[g]->oldest_age_group[k]->patch_no == FOLLOW_PATCH){
                     printf("Scheduling first HIV test (at current time %6.4f) for OLDEST ", t);
                     printf("adult %ld (event type %i)at time index %li, array index = %li\n",
                         age_list->age_list_by_gender[g]->oldest_age_group[k]->id,
                         age_list->age_list_by_gender[g]->oldest_age_group[k]->next_cascade_event,
                         age_list->age_list_by_gender[g]->oldest_age_group[k]->idx_cascade_event[0],
                         age_list->age_list_by_gender[g]->oldest_age_group[k]->idx_cascade_event[1]
                    );
                fflush(stdout);
            }
        }
    }
}


/* Function: draw_hiv_tests()
Schedule HIV test for the whole population at fixed times.  

This is the function which does the HIV test scheduling.  Function based on
draw_initial_hiv_tests(). Allows us to draw HIV testing after if we switch back to clinic-based
testing this function is called from simul.c.  Function goes through every individual currently
alive (using the array age_list) and schedules an HIV test for each person in the array
cascade_events[].  Only schedule tests for people who are not "HIV+ and aware of serostatus".  

Arguments
---------
param : pointer to a parameters structure
    All parameters of the simulation
age_list : pointer to a age_list_struct structure
    Multidimensional array
year : int
    Current year
cascade_events : 3D array of pointers to individual structures
    
n_cascade_events : array of pointer to long
size_cascade_events : array of pointers to long
country_setting : int
    Either ZAMBIA (1) or SOUTH_AFRICA (2) (see constants.h for these definition)

Returns
-------
Nothing; schedules HIV tests for the population
*/

void draw_hiv_tests(parameters *param, age_list_struct *age_list, int year, 
    individual ***cascade_events, long *n_cascade_events, long *size_cascade_events, 
    int country_setting){
    
    int aa,k,g;
    
    // Array for storing probability of an indiv getting tested in the next window. 
    // Can be made more complex to allow differences by age, past history of testing. 
    double *p_test;
    
    p_test = malloc(N_GENDER*sizeof(double));
    double t_gap;

    // Draw probability that an indiv is tested in the next testing window,store it in p_test
    probability_get_hiv_test_in_next_window(p_test, &t_gap, country_setting, year, 
        param->COUNTRY_HIV_TEST_START, param);
    
    for(g = 0; g < N_GENDER; g++){
        // For all but the age group 80+ (which is in a separate part of the age_list struct)
        for(aa = 0; aa < MAX_AGE - AGE_ADULT; aa++){
            // For each individual in that annual age group, schedule their first HIV test
            for(k = 0; k < age_list->age_list_by_gender[g]->number_per_age_group[aa]; k++){
                
                // Only schedule tests for people who are not "HIV+ and aware of serostatus"
                if(age_list->age_list_by_gender[g]->age_group[aa][k]->ART_status == ARTNEG){

                    schedule_hiv_test_fixedtime(age_list->age_list_by_gender[g]->age_group[aa][k],
                        param, year, cascade_events, n_cascade_events, size_cascade_events, 
                        t_gap, country_setting, p_test);
                    
                    if(
                    (age_list->age_list_by_gender[g]->age_group[aa][k]->id == FOLLOW_INDIVIDUAL) &&
                    (age_list->age_list_by_gender[g]->age_group[aa][k]->patch_no == FOLLOW_PATCH)
                    ){
                        
                        printf("Scheduling new HIV test for adult %ld, ", 
                            age_list->age_list_by_gender[g]->age_group[aa][k]->id);
                        
                        printf("gender = %d (event type %i) at time index %li array index = %li\n", 
                        g, age_list->age_list_by_gender[g]->age_group[aa][k]->next_cascade_event,
                        age_list->age_list_by_gender[g]->age_group[aa][k]->idx_cascade_event[0],
                        age_list->age_list_by_gender[g]->age_group[aa][k]->idx_cascade_event[1]);
                        
                        fflush(stdout);
                    }
                }
            }
        }
        // For those in the last age group, loop through all individuals
        // Note - we may decide to switch this off
        for(k = 0; k < age_list->age_list_by_gender[g]->number_oldest_age_group; k++){ 
            
            // Only schedule tests for people who are not "HIV+ and aware of serostatus"
            if(age_list->age_list_by_gender[g]->oldest_age_group[k]->ART_status == ARTNEG){
                
                schedule_hiv_test_fixedtime(age_list->age_list_by_gender[g]->oldest_age_group[k],
                    param, year, cascade_events, n_cascade_events, size_cascade_events, t_gap,
                    country_setting, p_test);
                
                if(
                (age_list->age_list_by_gender[g]->oldest_age_group[k]->id == FOLLOW_INDIVIDUAL) &&
                (age_list->age_list_by_gender[g]->oldest_age_group[k]->patch_no == FOLLOW_PATCH)
                ){
                    printf("Scheduling first HIV test  (at current year %i) ", year);
                    printf("for OLDEST adult %ld, ", 
                        age_list->age_list_by_gender[g]->oldest_age_group[k]->id);
                    printf("gender = %d (event type %i)at time index %li, array index = %li\n",
                    g, age_list->age_list_by_gender[g]->oldest_age_group[k]->next_cascade_event, 
                    age_list->age_list_by_gender[g]->oldest_age_group[k]->idx_cascade_event[0], 
                    age_list->age_list_by_gender[g]->oldest_age_group[k]->idx_cascade_event[1]);
                    
                    fflush(stdout);
                }
            }
        }
    }
    free(p_test);
}


/* Function: schedule_generic_cascade_event()
Schedule the generic care cascade events in the array cascade_events[].

Each cascade event (HIV test, CD4 test, start/interrupt ART) has the same format for scheduling events, this single function is used for scheduling all the different cascade events.  


Arguments
---------
individual* indiv
parameters *param
double t_event
individual ***cascade_events
long *n_cascade_events
long *size_cascade_events
double t_now

Returns
-------
Nothing; 
*/

void schedule_generic_cascade_event(individual* indiv, parameters *param, double t_event,
    individual ***cascade_events, long *n_cascade_events, long *size_cascade_events, double t_now){
    
    // This is the index for cascade_events (and n_cascade_events and size_cascade_events) arrays
    int idx_event = (int) (round((t_event - param->COUNTRY_HIV_TEST_START) *
        N_TIME_STEP_PER_YEAR));
    
    // Ensure that we never schedule a cascade event during the current timestep:
    int idx_current_time = (int) (round((t_now - param->COUNTRY_HIV_TEST_START) * 
        N_TIME_STEP_PER_YEAR));
    
    // Make sure event is not scheduled for the current time
    if(idx_event == idx_current_time){
        idx_event += 1;
    }
    
    if(idx_event < idx_current_time){
        printf("Error. Scheduled event in the past.");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }
    
    // 
    if(idx_event <= (param->end_time_simul - param->COUNTRY_HIV_TEST_START)*N_TIME_STEP_PER_YEAR){
        
        indiv->idx_cascade_event[0] = idx_event;
        indiv->idx_cascade_event[1] = n_cascade_events[idx_event];
        
        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("New generic cascade event call for adult ");
            printf("%ld from patch %d generated at time t=%6.4f with array indices  %ld %ld\n",
                indiv->id, indiv->patch_no, t_event, 
                indiv->idx_cascade_event[0], indiv->idx_cascade_event[1]);
            fflush(stdout);
        }
        
        // Check if we've run out of memory
        if(n_cascade_events[idx_event] >= (size_cascade_events[idx_event])){

            // Note that realloc does not work (we need to pass a pointer to the pointer, which is
            // really complicated as it propagates through several functions (so maybe make
            // planned_breakups[time_breakup] **), so ditch this code for now and use the following
            // lines
            printf("Unable to re-allocate cascade_events[i]. Execution aborted.");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        cascade_events[idx_event][n_cascade_events[idx_event]] = indiv;
        n_cascade_events[idx_event]++;
        
    }else{
        // If next event scheduled for after the end of the simulation set all to be dummy entries
        indiv->idx_cascade_event[0] = NOEVENT;
        indiv->idx_cascade_event[1] = -1;
        
        if(VERBOSE_OUTPUT == 1){
            printf("No event scheduled for %li, ", indiv->id);
            printf(" as the next event lies after the end of the simulation.\n");
        }
    }
}


void schedule_new_hiv_test(individual* indiv, parameters *param, double t, 
    individual ***cascade_events, long *n_cascade_events, long *size_cascade_events){
   /* For a given individual, draw when they will next have an HIV test.  
        
    This is used for HIVTESTSCHEDULE=0 (each individual schedules their own tests sequentially).
    Only people who have never been tested or previously received a negative test are scheduled to
    be tested. 
    
    See file HIV testing function.xlsx to see what this looks like.
    note - format is hill_down(x,max_val,exponent,midpoint), where midpoint is time (since t=0)
    at which midpoint occurs. hill_down->0 as x->infty

    Previous version (May 2015):
    double mean_time_to_test = 1.0 + hill_down(t - param->COUNTRY_HIV_TEST_START, 12.0, 4, 8.0);
    */
    
    /* This tells us if the cd4 test is due to PopART (is_popart=1) or not (is_popart=0). */
    int is_popart = (indiv->next_cascade_event >= NCASCADEEVENTS_NONPOPART);
    
    //hill_down(double x, double max_val, double exponent, double midpoint);
    double mean_time_to_test = param->time_to_background_HIVtestNOW + 
        hill_down(t - param->COUNTRY_HIV_TEST_START, param->time_to_background_HIVtest_maxval,
            param->time_to_background_HIVtest_exponent, param->time_to_background_HIVtest_midpoint);
    
    double x = gsl_ran_exponential(rng, mean_time_to_test);
    
    if(x <= TIME_STEP){
        x = TIME_STEP;
    }
    
    double time_hiv_test = t + x;
    
    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("t=%6.4f, index=%i. ", t, 
            (int) (round((t - param->COUNTRY_HIV_TEST_START) * N_TIME_STEP_PER_YEAR)));
        
        printf("New HIV test for adult %ld scheduled at time %6.4f. Mean time to test = %f\n", 
            indiv->id,time_hiv_test,mean_time_to_test);
        fflush(stdout);
    }
    
    // Now call schedule_generic_cascade_event() to actually schedule the event
    if(time_hiv_test < (param->end_time_simul)){
        if(is_popart == NOTPOPART){
            indiv->next_cascade_event = CASCADEEVENT_HIV_TEST_NONPOPART;
        }else{
            indiv->next_cascade_event = CASCADEEVENT_HIV_TEST_POPART;
        }
        
        schedule_generic_cascade_event(indiv, param, time_hiv_test, cascade_events,
            n_cascade_events, size_cascade_events, t);
    }else{
        // Some people never test (as exponential distribution we truncate any tests occuring after
        // the end of the simulation).  
        // Note that we do not need to unschedule this indiv from the cascade_event list as the
        // event is this one and we are moving on in time
        indiv->next_cascade_event = NOEVENT;
        indiv->idx_cascade_event[0] = NOEVENT;
        indiv->idx_cascade_event[1] = -1;
    }
}


void probability_get_hiv_test_in_next_window(double *p_test, double *t_gap, int country_setting,
    int year, int COUNTRY_HIV_TEST_START, parameters *param){
    /*
    Calculate probability an individual has an HIV test in the next window (t_gap).  
    
    Arguments
    ---------
    p_test : pointer to an array of doubles
    t_gap : pointer to a double
        time period in which to schedule the HIV test
    country_setting : int
        An identifier for the country in question (see constants.h for the macro definitions for
        each country).  
    year : int
        Year (as an integer) in which the HIV test is taking place.  
    COUNTRY_HIV_TEST_START : int
        The year in which HIV testing started in the country in question.  
    param : pointer to a parameters struct
        The struct holding all the parameters of the model.  
    
    Returns
    -------
    Nothing; the variables which `p_test` and `t_gap` point to are populated.  
    */
    
    // Check that HIV testing has started
    if(year < COUNTRY_HIV_TEST_START){
        printf("probability_get_hiv_test_in_next_window() ");
        printf("called before start of HIV testing.\n");
    }
    
    if(year == COUNTRY_HIV_TEST_START){
        /* This refers to the period [COUNTRY_HIV_TEST_START, 2006]. */
        p_test[MALE] = param->p_HIV_background_testing_female_pre2006 * param->RR_HIV_background_testing_male;
        p_test[FEMALE] = param->p_HIV_background_testing_female_pre2006;
        *t_gap = 2006 - COUNTRY_HIV_TEST_START;
    }else{
        p_test[MALE] = param->p_HIV_background_testing_female_current*param->RR_HIV_background_testing_male;
        p_test[FEMALE] = param->p_HIV_background_testing_female_current;
        *t_gap = 1;
    }
}


void schedule_hiv_test_fixedtime(individual* indiv, parameters *param, int year, 
    individual ***cascade_events, long *n_cascade_events, long *size_cascade_events, double t_gap,
    int country_setting, double *p_test){
    /* For a given individual draw when they will next have an HIV test
    
    
    The probability of having a test is gender-specific and listed in the array `p_test`.  The 
    period of time in which the HIV test has to occur is between `year` and `year + t_gap`.  The 
    time until the test is drawn as a uniform random variable between these periods.  
    
    This func is used for HIVTESTSCHEDULE = 1 (HIV test scheduling procedure happens for the whole
    population at fixed times).  Only people who have never been tested or previously received a
    negative test are scheduled to be tested.  
    
    
    Arguments
    ---------
    indiv : pointer to an individual struct
    param : pointer to a parameters struct
    year : int
        Year in question
    individual ***cascade_events : pointer to a pointer to a pointer to an individual object
        Used because `cascade_events` is a multi-dimensional array
    n_cascade_events : pointer to a long
    size_cascade_events : pointer to a long
    t_gap : double
        Period of time (in years) in which the test should occur.  
    country_setting : int
    p_test : pointer to an array of doubles
        Array of the gender-specific probability of having an HIV test.  
    
    Returns
    -------
    Nothing; attributes of the individual are updated and schedule_generic_cascade_event is called. 
    */
    
    /* This tells us if the cd4 test is due to PopART (is_popart=1) or not (is_popart=0). */
    int is_popart = (indiv->next_cascade_event >= NCASCADEEVENTS_NONPOPART);

    // Draw to see if this person will test or not
    double x_test = gsl_rng_uniform(rng);
    
    if(x_test > p_test[indiv->gender]){
        /* No test this time, so make sure they are unscheduled: */
        indiv->next_cascade_event = NOEVENT;
        indiv->idx_cascade_event[0] = NOEVENT;
        indiv->idx_cascade_event[1] = -1;
    }else{ 
        // Schedule a test

        // Determine when the test will be held (time_hiv_test). 
        // Draw a uniform random variable from 0 to t_gap 
        // This has to happen before the end of the current testing window, hence the "-TIME_STEP"
        double x_time = gsl_rng_uniform(rng) * (t_gap - TIME_STEP);

        // Ensure that this always happens in the future (ie. at least TIME_STEP later).  
        if(x_time <= TIME_STEP){
            x_time = TIME_STEP;
        }
        
        double time_hiv_test = year + x_time;

        if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
            printf("At year=%i, new HIV test for adult %ld scheduled at future time %6.4f.\n",
                year,indiv->id, time_hiv_test);
            fflush(stdout);
        }

        // Schedule event if occurs before end of simulation. This should always be true, 
        // but keep for now just to be sure. */
        if(time_hiv_test < (param->end_time_simul)){
            // Now call schedule_generic_cascade_event() to actually schedule the event
            if(is_popart == NOTPOPART){
                indiv->next_cascade_event = CASCADEEVENT_HIV_TEST_NONPOPART;
            }else{
                indiv->next_cascade_event = CASCADEEVENT_HIV_TEST_POPART;
            }
            schedule_generic_cascade_event(indiv, param, time_hiv_test, cascade_events,
                n_cascade_events, size_cascade_events, (double) year);
        }else{
            // Note that we do not need to unschedule this indiv from the cascade_event list 
            // as the event is this one
            // and we are moving on in time. 
            indiv->next_cascade_event = NOEVENT;
            indiv->idx_cascade_event[0] = NOEVENT;
            indiv->idx_cascade_event[1] = -1;
        }
    }
}


/* This carries out the processes of someone testing, and assigns them a next event based on the test result
 * (ie start ART, wait in care, drop out of care). 
 * is_popart takes the values 0 (not popart) and 1 (popart) to indicate whether testing is PopART-run (CHiPs) or not.
 * PopART testing 
 * */
void hiv_test_process(individual* indiv, parameters *param, double t, individual ***cascade_events, 
        long *n_cascade_events, long *size_cascade_events, individual ***hiv_pos_progression, 
        long *n_hiv_pos_progression, long *size_hiv_pos_progression, 
        cumulative_outputs_struct *cumulative_outputs, calendar_outputs_struct *calendar_outputs, 
        individual ***vmmc_events, long *n_vmmc_events, long *size_vmmc_events, 
        patch_struct *patch, int p, debug_struct *debug){

    int year_idx = (int) floor(t) - param->start_time_simul;
    double time_new_cd4;
    double x_sensitivity;
    /* This tells us if the cd4 test is due to PopART (is_popart=1) or not (is_popart=0). */
    int is_popart = (indiv->next_cascade_event>=NCASCADEEVENTS_NONPOPART);

    //printf("hiv_test_process is_popart= %i\n",is_popart);

    int WINDOWNEGATIVE = 0; /* 0 if HIV+ tests negative (because recent infection), 1 if tests positive. 
                                Assign default value so always defined (otherwise may behave weirdly?) */

    int TRUE_POS = 1;    /* 1 if a HIV+ person tests +ve. */
    /* Record the time of their most recent test (so we can count % of people who hve tested in last X months): */
    indiv->time_last_hiv_test = t;

    /* Add one to the cumulative total of HIV tests carried out: */
    if(is_popart == NOTPOPART){
        cumulative_outputs->N_total_HIV_tests_nonpopart++; // Test is not from CHiPs (ie this is a clinic-based test).
        calendar_outputs->N_calendar_HIV_tests_nonpopart[year_idx]++;
    }else{
        cumulative_outputs->N_total_HIV_tests_popart++; // Home-based CHiPs test - note we distinguish as allows economists to distinguish home/facility-based testing costs
        calendar_outputs->N_calendar_HIV_tests_popart[year_idx]++;
    }
    /* HIV test window for someone who IS infected (NOTE - this is coded as >UNINFECTED):
     * NOTE: we call this first as what we want to know is their HIV test result (not HIV status). */
    if (indiv->HIV_status > UNINFECTED){
        WINDOWNEGATIVE = get_window_result(t-(indiv->t_sc), t, patch);
        /* Only draw rng if need to (ie <100% sensitivity). */
        if (is_popart == POPART){
            if (WINDOWNEGATIVE==0 && param->HIV_rapid_test_sensitivity_CHIPS<1.0){
                x_sensitivity = gsl_rng_uniform (rng);
                if (x_sensitivity>(param->HIV_rapid_test_sensitivity_CHIPS))
                    TRUE_POS = 0; /* Test false negative. */
            }
        }
    }



    if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        printf("HIV test result for adult %ld at time %6.4f is %d. True HIV status is %i. ART status is %i\n",indiv->id,t,(indiv->HIV_status!=UNINFECTED && WINDOWNEGATIVE!=0),indiv->HIV_status,indiv->ART_status);
        fflush(stdout);
    }

    /************************************************** 
     * Now look at whether person tests HIV+ or not: 
     ************************************************** */
    /* If test HIV-, then schedule a new HIV test (and VMMC if applicable): */
    if (indiv->HIV_status==UNINFECTED || WINDOWNEGATIVE==0 || TRUE_POS==0){
        /* If before PopART, schedule a new HIV test in the future: */
        if (is_popart==NOTPOPART){
            /* if each person schedules their own HIV tests sequentially. */
            if (HIVTESTSCHEDULE==0){
                indiv->next_cascade_event = CASCADEEVENT_HIV_TEST_NONPOPART;
                schedule_new_hiv_test(indiv, param, t, cascade_events, n_cascade_events, size_cascade_events);
            }
            /* if the whole population schedules their next HIV test at fixed times. */
            else{
                /* Do not schedule a new event at present - the next cascade event is assumed to be a new test at the next fixed time when HIV tests are scheduled . */
                indiv->next_cascade_event = NOEVENT;
                indiv->idx_cascade_event[0] = NOEVENT;
                indiv->idx_cascade_event[1] = -1;
            }

        }else{
            /* For POPART, do not schedule a new event at present - the next cascade event is assumed to be a new annual test at the next round, which is done through the function carry_out_chips_visits_per_timestep() rather than carry_out_cascade_events_per_timestep(). */ 
            indiv->next_cascade_event = NOEVENT;
            indiv->idx_cascade_event[0] = NOEVENT;
            indiv->idx_cascade_event[1] = -1;
            
            cumulative_outputs->N_total_HIV_tests_popart_negative++;
            calendar_outputs->N_calendar_HIV_tests_popart_negative[year_idx]++;
        }

        if (indiv->gender==MALE)
            if (indiv->circ==UNCIRC)   /* Only if not already circumcised (and not waiting for VMMC): */
                draw_if_VMMC(indiv,param,vmmc_events,n_vmmc_events,size_vmmc_events,t,is_popart);
        /* Rest of the code in this function is if test HIV+, so return. */
        return;
    }

    /* if get to here the individual tested HIV+ */

    if(is_popart == POPART){
        cumulative_outputs->N_total_HIV_tests_popart_positive++;
        calendar_outputs->N_calendar_HIV_tests_popart_positive[year_idx]++;
    }

    if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        printf("Adult %ld tested positive at time %6.4f\n",indiv->id,t);
        fflush(stdout);
    }

    /* PANGEA stuff: get the CD4 at diagnosis for HIV+ person: */
    indiv->PANGEA_cd4atdiagnosis = PANGEA_get_cd4(indiv, t); 
    indiv->PANGEA_t_diag = t;


    indiv->ART_status = ARTNAIVE;  /* Status changes as now known positive. */

    /* If ART has started then there are 3 possibilities for the next cascade event
     *  - drops out, waits until eligible, starts ART. */
    if (t>=param->COUNTRY_ART_START){       
        /* Only HIV+ diagnosed individuals remain in the function now. Three possible outcomes: */
        //printf("Need to schedule stuff for %ld\n",indiv->id);
        /* 1. Refuses to enter treatment/care. Next event will be determined by CD4 (e.g. wait until CD4<200), so add in code to update this as CD4 progression occurs. */
        if (joins_preart_care(indiv, param, t, cumulative_outputs, calendar_outputs) < 1){
            debug->art_vars[p].cascade_transitions[ARTNEG+1][CASCADEDROPOUT+1]++;
            /* Dropping out of the cascade is immediate once the HIV test is done. */
            dropout_process(indiv,param,t, cascade_events, n_cascade_events, size_cascade_events,
                hiv_pos_progression, n_hiv_pos_progression, size_hiv_pos_progression,
                cumulative_outputs, calendar_outputs);
        }
        /* 2. Eligible for ART and starts (after a delay) - next event is starting ART. 
         * Note that delays between HIV testing and starting ART (including getting CD4, picking up results, drugs)
         *  all built in to this function. 
         * Eligibility for ART is determined by calendar time and trial arm. */
        else if (is_eligible_for_art(indiv,param,t, patch, p)>0){
            debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTNAIVE+1]++;
            schedule_start_of_art(indiv,param,t, cascade_events, n_cascade_events, size_cascade_events);
        }
        /* 3. Not eligible for ART, so next event is a new CD4 test. */
        else{
            /* People test t_cd4_retest_min to t_cd4_retest_min+t_cd4_retest_range months after their last CD4 test (or first HIV test). 
             * Note we separate out the first and subsequent tests as the first CD4 test also includes 
             * the time for HIV testing etc so is longer.*/
            /* Takes an extra time due to the delay in between getting the HIV test and the CD4 test for the first CD4 test:
             * compared to between the nth and n+1th CD4 test. 
             * This is reflected by the parameters::
             * t_delay_hivtest_to_cd4test_min to t_delay_hivtest_to_cd4test_min+t_delay_hivtest_to_cd4test_range
             */
            if (is_popart==NOTPOPART){
                /* Time to next CD4 test is the sum of the time between getting the HIV test and having the first CD4 test, 
                 * and the time between consecutive CD4 tests. */
                time_new_cd4 = t + param->t_delay_hivtest_to_cd4test_min[NOTPOPART] + param->t_delay_hivtest_to_cd4test_range[NOTPOPART] * gsl_rng_uniform (rng);
                //+ param->t_cd4_retest_min[NOTPOPART] + param->t_cd4_retest_range[NOTPOPART] * gsl_rng_uniform (rng);

                indiv->next_cascade_event = CASCADEEVENT_CD4_TEST_NONPOPART;
            }
            /* For people in arm B: */
            else{

                /* Allow CD4 retesting to be more frequent. Time to next CD4 test is again the sum of the time between 
                 * getting the HIV test and having the first CD4 test, and the time between consecutive CD4 tests.*/
                /* DSMBBBBBBBBB - I think this is wrong!. */
                time_new_cd4 = t + param->t_delay_hivtest_to_cd4test_min[POPART] + param->t_delay_hivtest_to_cd4test_range[POPART] * gsl_rng_uniform (rng);
                //+ param->t_cd4_retest_min[POPART]    + param->t_cd4_retest_range[POPART] * gsl_rng_uniform (rng);

                indiv->next_cascade_event = CASCADEEVENT_CD4_TEST_POPART;
            }
            debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTNAIVE+1]++;
            schedule_generic_cascade_event(indiv, param, time_new_cd4, cascade_events, n_cascade_events, size_cascade_events,t);
        }
    }
    /* This is for what happens during the period when people were testing for HIV but ART was not widely available:
     * If test HIV+ before ART becomes available then the next cascade event is:
     *  either get a CD4 test in the future, or drop out: */
    else{
        /* 1. Drops out Next event will be determined by CD4 (e.g. wait until CD4<200), so add in code to update this as CD4 progression occurs. */
        if (!joins_preart_care(indiv,param,t,cumulative_outputs,calendar_outputs)){ /// SAME COMMENT AS BEFORE
            debug->art_vars[p].cascade_transitions[ARTNEG+1][CASCADEDROPOUT+1]++;
            dropout_process(indiv,param,t, cascade_events, n_cascade_events, size_cascade_events, hiv_pos_progression, n_hiv_pos_progression, size_hiv_pos_progression, cumulative_outputs,calendar_outputs);
        }else{
            /* 2. Has a CD4 test after ART becomes available. */
            /* Has CD4 test param->t_cd4_whenartfirstavail_min to param->t_cd4_whenartfirstavail_min+param->t_cd4_whenartfirstavail_range 
             * yrs after ART becomes available in the country. */
            double time_cd4;
            if (is_popart==NOTPOPART){
                time_cd4 = param->COUNTRY_ART_START + param->t_cd4_whenartfirstavail_min + param->t_cd4_whenartfirstavail_range*gsl_rng_uniform (rng);
                indiv->next_cascade_event = CASCADEEVENT_CD4_TEST_NONPOPART;
            }
            else{
                /* This should probably not be called as PopART should never happen before CD4 testing is available. */ 
                printf("ERROR: Not sure why here in hiv_test_process() when CD4 testing unavailable but PopART is!!!\n");
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
                //double time_cd4 = param->COUNTRY_ART_START + gsl_rng_uniform (rng);
                //indiv->next_cascade_event = CASCADEEVENT_CD4_TEST_POPART;
            }
            debug->art_vars[p].cascade_transitions[ARTNEG+1][ARTNAIVE+1]++;
            schedule_generic_cascade_event(indiv, param, time_cd4, cascade_events, n_cascade_events, size_cascade_events,t);
        }
    }
}


void schedule_start_of_art(individual* indiv, parameters *param, double t, 
    individual ***cascade_events, long *n_cascade_events, long *size_cascade_events){
   /* Give delay in starting ART after testing HIV+ (for those immediately eligible).
    
    This represents all the delays in collecting the HIV test, getting the CD4 test and results,
    and getting drugs.  
    
    Arguments
    ---------
    indiv : pointer to an individual structure
        The individual in question
    param : pointer to a parameters structure
        All parameters of the simulation
    t : double
        Current time (decimal year)
    cascade_events : pointer to a multidimensional array of individual structures
        
    n_cascade_events : pointer to an array of long
        
    size_cascade_events : pointer to an array of long
    
    Returns
    -------
    Nothing; the individual's next_cascade_event is updated, and the time until they start ART is
    drawn from a bi-exponential distribution, time_start_art, before calling the function
    schedule_generic_cascade_event() to actually schedule the event.  
    */
    
    double tmp, t_delay, time_start_art;
    double p_fast = 0, t_fast = 0, t_slow = 0;
    int idx_round = 0, period_in_round = 0;
    int year, t_step, timesteps_in_round_so_far; /* Discrete time if needed. */
    
    // Is CD4 test due to PopART (is_popart=1) or not (is_popart=0)?  
    int is_popart = (indiv->next_cascade_event >= NCASCADEEVENTS_NONPOPART);
    
    // Determine time until starting ART for those in the background process
    if(is_popart == NOTPOPART){
        
        // Uniform version:
        //time_start_art = t + param->t_start_art_min[NOTPOPART] + 
        //    gsl_rng_uniform (rng)*param->t_start_art_range[NOTPOPART];
        
        // Exponential version:
        t_delay = gsl_ran_exponential (rng, param->t_start_art_mean_non_popart);
        indiv->next_cascade_event = CASCADEEVENT_START_ART_NONPOPART;
        //printf("Add in counter for CD4 here!!!\n");
    }else{ 
        // Determine time until starting ART for those in PopART
        
        // Uniform version:
        //time_start_art = t + param->t_start_art_min[POPART] + 
        //    gsl_rng_uniform (rng)*param->t_start_art_range[POPART];
        
        // Exponential version:
        // time_start_art = t + gsl_ran_exponential (rng, param->t_start_art_mean[POPART]);
        
        // Biexponential version:
        tmp = gsl_rng_uniform(rng);

        while(t > param->CHIPS_END_YEAR[idx_round] + 
            param->CHIPS_END_TIMESTEP[idx_round]*TIME_STEP && idx_round < NCHIPSROUNDS - 2){
            idx_round++; // which round we're in
        }
        if(idx_round == 0){
            year = (int) t;
            t_step = (int) round((t - year)*N_TIME_STEP_PER_YEAR);
            
            timesteps_in_round_so_far =
                (year-param->CHIPS_START_YEAR[idx_round])*N_TIME_STEP_PER_YEAR +
                (t_step - param->CHIPS_START_TIMESTEP[idx_round]);
            
            period_in_round = (int) floor( (param->n_time_periods_art_popart_per_round[idx_round]*timesteps_in_round_so_far) / param->chips_params->n_timesteps_per_round[idx_round]); // which period in the round we're in
        }

        if(period_in_round >= 6 || period_in_round < 0){
            printf("Error - period in CHiPs round=%i outside range at time %6.4lf. Exiting\n",
                period_in_round,t);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }
        
        p_fast = param->p_start_art_mean_fast_popart[idx_round][period_in_round];
        t_fast = param->t_start_art_mean_fast_popart[idx_round][period_in_round];
        t_slow = param->t_start_art_mean_slow_popart[idx_round][period_in_round];
        
        if(tmp < p_fast){  // THIS IS NOT WORKING?
            t_delay = gsl_ran_exponential (rng, t_fast);
        }else{
            t_delay = gsl_ran_exponential (rng, t_slow);
        }
        indiv->next_cascade_event = CASCADEEVENT_START_ART_POPART;
        //printf("Add in counter for CD4 here!!!\n");
        //printf("time_start_art = %lf\n",time_start_art); // THIS ONE
        //printf("t=%lf\n",t);
    }

    /*Ensure that event happens at least 1 timestep in the future. */
    if(t_delay < TIME_STEP){
        t_delay = TIME_STEP;
    }
    time_start_art = t + t_delay;
    indiv->t_start_art = time_start_art;

    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Scheduling start of ART for adult %ld at time %6.4f\n", indiv->id, time_start_art);
        fflush(stdout);
    }
    schedule_generic_cascade_event(indiv, param, time_start_art, 
        cascade_events, n_cascade_events, size_cascade_events,t);
}


void cd4_test_process(individual* indiv, parameters *param, double t, individual ***cascade_events,
    long *n_cascade_events, long *size_cascade_events, individual ***hiv_pos_progression, 
    long *n_hiv_pos_progression, long *size_hiv_pos_progression, 
    cumulative_outputs_struct *cumulative_outputs, calendar_outputs_struct *calendar_outputs,
    patch_struct *patch, int p){
    /* 
    This is the process of repeated CD4 testing for someone who wasn't eligible for ART at
    their last CD4 test.  Note that the first CD4 test is carried out implicitly as part of
    hiv_test_process.  This function updates everything and schedules the next cascade event that
    will happen.  
    
    There are 3 possible events (as for when someone gets an HIV+ test): 
    
    1. Refuses to enter treatment/care. 
        In this case the next event will be determined by CD4 (e.g. wait until CD4<200), so add in
        code to update this as CD4 progression occurs.  Note that leaving is immediate, so we don't
        schedule an event.  
    
    2. Eligible for ART and starts (after a delay)
        In this case the next event is starting ART.  Note that delays between HIV testing and
        starting ART (including getting CD4, picking up results, drugs) all built in to this
        function.  
        
    3. Not eligible for ART.
        In this case the next event is a new CD4 test.  
    
    Assume people test 6-12 months after their last CD4 test (or first HIV test).  Note we separate
    out the first and subsequent tests as the first CD4 test also includes the time for HIV testing
    etc so is longer.  Make sure this is similar to the function for the time to the first repeated
    CD4 test.  Also note that testing is assumed to happen - people may then drop out after their
    test - so count all tests in cumulative outputs.  
    
    
    Arguments
    ---------
    individual* indiv
    parameters *param
    double t
    individual ***cascade_events
    long *n_cascade_events
    long *size_cascade_events
    individual ***hiv_pos_progression
    long *n_hiv_pos_progression
    long *size_hiv_pos_progression
    cumulative_outputs_struct *cumulative_outputs
    patch_struct *patch
    int p
    
    Returns
    -------
    Nothing; other events are scheduled and updated.  
    */
    
    // This tells us if the cd4 test is due to PopART (is_popart=1) or not (is_popart=0). 
    int is_popart = (indiv->next_cascade_event >= NCASCADEEVENTS_NONPOPART);
    int year_idx = (int) floor(t) - param->start_time_simul;
    
    // Count cumulative number of CD4 tests.  
    if(is_popart == NOTPOPART){
        cumulative_outputs->N_total_CD4_tests_nonpopart++;
        calendar_outputs->N_calendar_CD4_tests_nonpopart[year_idx]++;
    }else{
        cumulative_outputs->N_total_CD4_tests_popart++;
        calendar_outputs->N_calendar_CD4_tests_popart[year_idx]++;
    }
    
    // Print extra output if we're following a specific individual
    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("CD4 test process for adult %ld carried out at time %6.4f\n", indiv->id, t);
        fflush(stdout);
    }
    
    // Case 1: Refuses to enter treatment/care - next event will be determined by CD4 (see above)
    if(remains_in_cascade(indiv, param, is_popart) == 0){ 
        // If individual drops out, do nothing - wait until CD4 drops below 200.
        // Note that we do not need to unschedule this indiv from the cascade_event[] list as the
        // event is this one and we are moving on in time. 
        dropout_process(indiv, param, t, cascade_events, n_cascade_events, size_cascade_events,
            hiv_pos_progression, n_hiv_pos_progression, size_hiv_pos_progression,
            cumulative_outputs, calendar_outputs);
        return;
    
    // Case 2: Eligible for ART and starts (after a delay) - next event is starting ART. 
    }else if(is_eligible_for_art(indiv, param, t, patch, p) == 1){
        indiv->ART_status = ARTNAIVE;
        
        schedule_start_of_art(indiv, param,t, cascade_events, 
            n_cascade_events, size_cascade_events);
    
    // Case 3: Not eligible for ART - next event is a new CD4 test (see above).
    }else{
        indiv->ART_status = ARTNAIVE;
        
        double time_new_cd4;
        if(is_popart == NOTPOPART){
            time_new_cd4 = t + param->t_cd4_retest_min[NOTPOPART] +
                param->t_cd4_retest_range[NOTPOPART] * gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_CD4_TEST_NONPOPART;
        }else{
            time_new_cd4 = t + param->t_cd4_retest_min[POPART] + 
                param->t_cd4_retest_range[POPART] * gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_CD4_TEST_POPART;
        }
        schedule_generic_cascade_event(indiv, param, time_new_cd4, cascade_events,
            n_cascade_events, size_cascade_events, t);
    }
}


/* This sets up everything when someone becomes virally suppressed on ART and determines their next 
 * cascade event. */
void virally_suppressed_process(individual* indiv, parameters *param, double t, individual ***cascade_events, long *n_cascade_events, long *size_cascade_events,individual ***hiv_pos_progression, long *n_hiv_pos_progression, long *size_hiv_pos_progression){

    /* Only need to remove HIV progression event if becoming VS after being VU (those who go from early ART to VS are
     * assumed to already not be progressing).
     * Note that EVENTAFTERENDSIMUL also does not need to be removed. */
    if (indiv->next_HIV_event!=NOEVENT){
        /* Do not need to remove from hiv_pos_progression if there is no event scheduled as it would have occurred after the end of the simulation. */
        if (indiv->next_HIV_event>EVENTAFTERENDSIMUL)
            /* Note the final '2' argument means that the person is virally suppressed on ART, not dying. */
            remove_from_hiv_pos_progression(indiv,  hiv_pos_progression, n_hiv_pos_progression, size_hiv_pos_progression,t, param,2);
        indiv->idx_hiv_pos_progression[0] = -1;
        indiv->idx_hiv_pos_progression[1] = -1;
    }

    /* This tells us if the cd4 test is due to PopART (is_popart=1) or not (is_popart=0). */
    int is_popart = (indiv->next_cascade_event>=NCASCADEEVENTS_NONPOPART);

    //printf("virally_suppressed_process is_popart= %i\n",is_popart);

    indiv->ART_status = LTART_VS;


    if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        printf("Adult %ld has become VS at time %6.4f. Next HIV event is now %i\n",indiv->id,t,indiv->next_HIV_event);
        fflush(stdout);
    }
    /* Now decide what happens next. 3 possible states: continues being virally suppressed, becomes unsuppressed, or drops out of ART entirely. */
    double x = gsl_rng_uniform (rng);
    double p_stays_vs;
    if(indiv->gender == MALE){
        p_stays_vs = param->p_stays_virally_suppressed * param->p_stays_virally_suppressed_male;
    }else{
        p_stays_vs = param->p_stays_virally_suppressed;
    }
    
    if(x < p_stays_vs){
        /* Patient stays virally suppressed, no future cascade events for now: */
        indiv->next_cascade_event= NOEVENT;
        /* Set to dummy values */
        /* Note that we do not need to unschedule this indiv from the cascade_event[] list as the event is this one
         * and we are moving on in time. */
        indiv->idx_cascade_event[0] = NOEVENT; 
        indiv->idx_cascade_event[1] = -1;

        if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
            printf("Adult %ld will remain VS permanently with HIV status %i, ART status %i \n",indiv->id,indiv->HIV_status,indiv->ART_status);
            fflush(stdout);
        }
    }else if(x < ((p_stays_vs) + (param->p_stops_virally_suppressed))){
        /* Patient eventually ceases being virally suppressed. Schedule a transition to being virally unsuppressed. */

        /* Assume: randomly picked this to be 3-6 years later. */
        double time_end_vs;
        if (is_popart==NOTPOPART){
            time_end_vs = t + param->t_end_vs_becomevu_min[NOTPOPART] + param->t_end_vs_becomevu_range[NOTPOPART] * gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_VU_NONPOPART;
        }
        else{
            time_end_vs = t + param->t_end_vs_becomevu_min[POPART] + param->t_end_vs_becomevu_range[POPART] * gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_VU_POPART;
        }
        schedule_generic_cascade_event(indiv, param, time_end_vs, cascade_events, n_cascade_events, size_cascade_events,t);

        /* PANGEA stuff: record end of viral suppression. Note that this is in the future, so 
         * could be after someone dies, so need to fix in outputs to ensure this never happens. */
        indiv->PANGEA_date_endfirstVLsuppression = time_end_vs;
        if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
            printf("Adult %ld is scheduled to become VU at time %6.4f\n",indiv->id,time_end_vs);
            fflush(stdout);
        }
    }
    else{
        /* ASSUMPTION! Randomly picked this to be 3-6 years later. */
        double time_dropout;
        if (is_popart==NOTPOPART){
            time_dropout = t + param->t_end_vs_dropout_min[NOTPOPART] + param->t_end_vs_dropout_range[NOTPOPART] * gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_DROPOUT_NONPOPART;
        }
        else{
            time_dropout = t + param->t_end_vs_dropout_min[POPART] + param->t_end_vs_dropout_range[POPART] * gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_DROPOUT_POPART;
        }
        if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
            printf("Adult %ld is scheduled to drop out at time %6.4f\n",indiv->id,time_dropout);
            fflush(stdout);
        }
        schedule_generic_cascade_event(indiv, param, time_dropout, cascade_events, n_cascade_events, size_cascade_events,t);
    }
}

/* This sets up everything when someone becomes virally unsuppressed and determines their next cascade event. */
void virally_unsuppressed_process(individual* indiv, parameters *param, double t, individual ***cascade_events, long *n_cascade_events, long *size_cascade_events, individual ***hiv_pos_progression, long *n_hiv_pos_progression, long *size_hiv_pos_progression, cumulative_outputs_struct *cumulative_outputs, calendar_outputs_struct *calendar_outputs){

    /* This tells us if the cd4 test is due to PopART (is_popart=1) or not (is_popart=0). */
    int is_popart = (indiv->next_cascade_event>=NCASCADEEVENTS_NONPOPART);

    //printf("virally_unsuppressed_process is_popart= %i\n",is_popart);


    indiv->ART_status = LTART_VU;


    /* Need to allow CD4 progression again. */
    next_hiv_event(indiv, hiv_pos_progression, n_hiv_pos_progression, size_hiv_pos_progression, param, t, cumulative_outputs, calendar_outputs);

    if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        printf("Adult %ld has become VU at time %6.4f. next_hiv_event = %i indices=%li %li \n",indiv->id,t,indiv->next_HIV_event,indiv->idx_hiv_pos_progression[0],indiv->idx_hiv_pos_progression[1]);
        fflush(stdout);
    }

    /* Now decide what cascade event happens next. 3 possible states: continues being virally unsuppressed,
     * becomes suppressed, or drops out of ART entirely. */
    double x = gsl_rng_uniform (rng);
    if (x<(param->p_vu_becomes_virally_suppressed)){
        /* Next event is that patient becomes virally suppressed */
        /* Assume takes 2-3 years*/
        double time_become_vs;
        if (is_popart==NOTPOPART){
            time_become_vs = t + param->t_end_vu_becomevs_min[NOTPOPART]  + param->t_end_vu_becomevs_range[NOTPOPART]*gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_VS_NONPOPART;
        }
        else{
            time_become_vs = t + param->t_end_vu_becomevs_min[POPART]     + param->t_end_vu_becomevs_range[POPART]*gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_VS_POPART;
        }
        schedule_generic_cascade_event(indiv, param, time_become_vs, cascade_events, n_cascade_events, size_cascade_events,t);
        if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
            printf("Next cascade event for adult %ld is becoming VS\n",indiv->id);
            fflush(stdout);
        }
    }
    else{
        /* Assume: Randomly picked this to be 1-2 years later. */
        double time_dropout;
        if (is_popart==NOTPOPART){
            time_dropout = t + param->t_end_vu_dropout_min[NOTPOPART] + param->t_end_vu_dropout_range[NOTPOPART] * gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_DROPOUT_NONPOPART;
        }
        else{
            time_dropout = t + param->t_end_vu_dropout_min[POPART]    + param->t_end_vu_dropout_range[POPART] * gsl_rng_uniform (rng);
            indiv->next_cascade_event = CASCADEEVENT_DROPOUT_POPART;
        }
        schedule_generic_cascade_event(indiv, param, time_dropout, cascade_events, n_cascade_events, size_cascade_events,t);
        if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
            int array_index_for_now = (int) (round((t - param->COUNTRY_HIV_TEST_START) * N_TIME_STEP_PER_YEAR));
            int array_index_for_future_event = (int) (round((t - param->COUNTRY_HIV_TEST_START) * N_TIME_STEP_PER_YEAR));
            printf("Next cascade event for adult %ld is dropping out\n",indiv->id);

            printf("Current cascade indices at time %6.4f are %li %li. Dropout scheduled to occur at %f Indices are %i %i\n",t,indiv->idx_cascade_event[0],indiv->idx_cascade_event[1],time_dropout,array_index_for_now,array_index_for_future_event);
            fflush(stdout);
        }
    }
}


void dropout_process(individual* indiv, parameters *param, double t, individual ***cascade_events,
    long *n_cascade_events, long *size_cascade_events, individual ***hiv_pos_progression, 
    long *n_hiv_pos_progression, long *size_hiv_pos_progression, 
    cumulative_outputs_struct *cumulative_outputs, calendar_outputs_struct *calendar_outputs){
    /* Process events when someone drops out of care (including restarting their CD4
    progression if needed) and determining their next cascade event. 
    
    Arguments
    ---------
    indiv : pointer to an individual structure
    param : pointer to a parameter structure
    t : double
        Time out in question
    cascade_events : pointer (x3) to an individual (pointer to a multi-dimensional array of indivs)
    n_cascade_events : pointer to a long
    size_cascade_events : pointer to a long
    hiv_pos_progression : pointer (x3) to an individual structure
    cumulative_outputs : pointer to a cumulative_outputs_struct
        Structure that records cumultive numbers of events in the simulation. Dropout events need to
        be counted (used within next_hiv_event())
    
    Returns
    -------
    Nothing; updates different attributes of the individual structure.  
    
    */
    
    /* This tells us if the person dropped out during PopART (is_popart=1) or not (is_popart=0). */
    //int is_popart = (indiv->next_cascade_event>=NCASCADEEVENTS_NONPOPART);
    
    if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
        printf("Adult %ld has dropped out at time %6.4f. ", indiv->id, t);
        printf("Before dropping out next HIV event was %i ", indiv->next_HIV_event);
        printf("ART stat=%i next cascade event=%i\n",
            indiv->ART_status, indiv->next_cascade_event);
        fflush(stdout);
    }
    
    int g = indiv->gender, cd4 = indiv->cd4, spvl = indiv->SPVL_cat;
    // Find age of the dead person when they died
    int age = (int) floor(t - indiv->DoB);
    int year_idx = (int) floor(t) - param->start_time_simul;
    
    // Find the age index of this person (>=80 is its own category)
    // truncate >=80 to 80 so that indexing of FIND_AGE_GROUPS_UNPD works
    if(age >= MAX_AGE){
        age = MAX_AGE;
    }
    int age_idx = FIND_AGE_GROUPS_UNPD[age - AGE_ADULT];
    
    if( ( indiv->ART_status == EARLYART ) || 
        ( indiv->ART_status == LTART_VS ) || 
        ( indiv->ART_status == LTART_VU )
    ){
        calendar_outputs->N_calendar_AnnualDropoutOnART[g][age_idx][cd4][spvl][year_idx]++;
    }
    
    if(WRITE_COST_EFFECTIVENESS_OUTPUT == 1){
        calendar_outputs->N_calendar_dropout[year_idx]++;
    }
    
    indiv->ART_status = CASCADEDROPOUT;
    
    /* Need to allow CD4 progression again if they don't currently have CD4 progression event 
    (ie if they were on ART and virally suppressed). */
    /* CHANGE 8/1/15. */
    if(indiv->next_HIV_event == NOEVENT){
        next_hiv_event(indiv, hiv_pos_progression, n_hiv_pos_progression, 
            size_hiv_pos_progression, param, t, cumulative_outputs, calendar_outputs);
    }
    
    /* 
    No further cascade event happens. The two possibilities at present are that the person 
    dies from AIDS-related illness or starts emergency ART once their CD4 goes below 200.
    Both these possibilities are scheduled through the hiv_pos_progression array. 
    */ 

    /* Note that we do not need to unschedule this indiv from the cascade_event[] list as the event
    is this one and we are moving on in time. */
    indiv->next_cascade_event = NOEVENT;
    indiv->idx_cascade_event[0] = NOEVENT;
    indiv->idx_cascade_event[1] = -1;
    if(indiv->id == FOLLOW_INDIVIDUAL && indiv->patch_no == FOLLOW_PATCH){
        printf("Adult %ld has dropped out at time %6.4f. ", indiv->id, t);
        printf("After dropping out next HIV event is %i. ", indiv->next_HIV_event);
        printf("hiv_pos_progression indices %li %li\n", 
            indiv->idx_hiv_pos_progression[0], indiv->idx_hiv_pos_progression[1]);
        fflush(stdout);
    }
}


/* Go through the list of scheduled cascade events for this timestep (stored in hiv_pos_progression). 
 * For each person to whom some HIV event happens, draw their next HIV-based event (via next_hiv_event) unless they die from AIDS at this timestep.
 * */
void carry_out_cascade_events_per_timestep(double t, patch_struct *patch, int p, all_partnerships *overall_partnerships, debug_struct *debug, file_struct *file_data_store){

    int array_index_for_cascade_event = (int) (round((t - patch[p].param->COUNTRY_HIV_TEST_START) * N_TIME_STEP_PER_YEAR));
    int n_events = patch[p].n_cascade_events[array_index_for_cascade_event];
    individual *indiv;
    int n;
    //printf("Calling carry_out_cascade_events at time t=%f\n",t);
    for (n=0; n<n_events; n++){
        indiv = patch[p].cascade_events[array_index_for_cascade_event][n];

        if(indiv->id==FOLLOW_INDIVIDUAL && indiv->patch_no==FOLLOW_PATCH){
            printf("Adult %ld from patch %d is entering function carry_out_cascade_events_per_timestep at time %f, array_index_for_cascade_event is %i\n",indiv->id,indiv->patch_no, t, array_index_for_cascade_event);
            fflush(stdout);
        }

        /* Note that if indiv->idx_cascade_event[0]<0 then the event is either "no event scheduled" (-1) or "an event was scheduled but only for after the end of the simulation" (-2).
         * What can happen is that an individual dropped out or has a cascade event scheduled after the end of the simulation, but then the intervention acts to give them a new cascade event.
         */
        if (indiv->debug_last_cascade_event_index==indiv->idx_cascade_event[0] && indiv->idx_cascade_event[0]>=0){
            printf("ERROR - trying to schedule a new cascade event (type=%d) in carry_out_cascade_events_per_timestep() that occurs at the same time as the previous event for person %ld in patch %d at time = %6.4f. Exiting\n",indiv->next_cascade_event,indiv->id,indiv->patch_no,t);
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }

        if (indiv->cd4==DEAD){
            /* Move on to the next person. Note - we can set up a similar procedure to other lists to remove this person from this list but
               it is not necessary. As things stand, no hiv event happens to the dead person and no new event is scheduled for them. */
            continue;
        }
        if(indiv->cd4 == DUMMYVALUE){
            printf("Error. Using uninitialised person.\n");
            printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
            fflush(stdout);
            exit(1);
        }


        /* Decide if we kill this person: */

        if ((indiv->next_cascade_event==CASCADEEVENT_ARTDEATH_NONPOPART)||(indiv->next_cascade_event==CASCADEEVENT_ARTDEATH_POPART)){
            /* Update counters of cascade transitions first as these depend on the current state (they get modified by the subsequent function): */
            debug->art_vars[p].cascade_transitions[indiv->ART_status+1][ARTDEATH+1]++;
            if (indiv->ART_status==LTART_VU)
                indiv->DEBUG_cumulative_time_on_ART_VU += t - indiv->DEBUG_time_of_last_cascade_event;
            /* Note that the last argument '3' indicates that this is an AIDS-related death. */
            remove_from_hiv_pos_progression(indiv,  patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression,t, patch[p].param,3);
            /* Function removes person from everything except the cascade event list: */
            individual_death_AIDS(patch[p].age_list, indiv, patch[p].n_population, patch[p].n_population_oneyearagegroups, patch[p].n_infected, patch[p].n_population_stratified, t, patch[p].param,
                    overall_partnerships->susceptible_in_serodiscordant_partnership, overall_partnerships->n_susceptible_in_serodiscordant_partnership, overall_partnerships->pop_available_partners, overall_partnerships->n_pop_available_partners, patch[p].cascade_events, patch[p].n_cascade_events, patch[p].size_cascade_events, patch, p, file_data_store);
        }

        else{
            /* At this point "indiv->next_cascade_event" is what happens to them now. We then draw a new "next event" after. */
            if((indiv->next_cascade_event==CASCADEEVENT_HIV_TEST_NONPOPART)||(indiv->next_cascade_event==CASCADEEVENT_HIV_TEST_POPART)){
                /* Counters get modified in the function hiv_test_process. */
                //debug->art_vars[p].cascade_transitions[indiv->ART_status+1][ARTNEG+1]++;
                if (indiv->ART_status==ARTNAIVE)
                    printf("Individual %li transitioning illegally in patch %i\n.",indiv->id,p);
                /* Get their HIV test results and schedule the next cascade event accordingly. Note
                 * that the NOTPOPART indicates that this happens outside PopART. */
                hiv_test_process(indiv, patch[p].param, t, patch[p].cascade_events, patch[p].n_cascade_events, patch[p].size_cascade_events, patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression, patch[p].cumulative_outputs, patch[p].calendar_outputs, patch[p].vmmc_events, patch[p].n_vmmc_events, patch[p].size_vmmc_events, patch, p, debug);

            }
            else if ((indiv->next_cascade_event==CASCADEEVENT_CD4_TEST_NONPOPART)||(indiv->next_cascade_event==CASCADEEVENT_CD4_TEST_POPART)){
                //debug->art_vars[p].cascade_transitions[indiv->ART_status+1][ARTNAIVE+1]++;
                cd4_test_process(indiv, patch[p].param, t, patch[p].cascade_events, patch[p].n_cascade_events, patch[p].size_cascade_events, patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression, patch[p].cumulative_outputs, patch[p].calendar_outputs, patch, p);
            }
            else if ((indiv->next_cascade_event==CASCADEEVENT_START_ART_NONPOPART)||(indiv->next_cascade_event==CASCADEEVENT_START_ART_POPART)){
                if (indiv->PANGEA_date_firstARTstart<0){
                    indiv->PANGEA_date_firstARTstart = t;
                    indiv->PANGEA_cd4atfirstART = PANGEA_get_cd4(indiv, t);
                }
                debug->art_vars[p].cascade_transitions[indiv->ART_status+1][EARLYART+1]++;
                start_ART_process(indiv, patch[p].param, t, patch[p].cascade_events, patch[p].n_cascade_events, patch[p].size_cascade_events, patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression,0, file_data_store, patch[p].calendar_outputs);
            }
            else if ((indiv->next_cascade_event==CASCADEEVENT_VS_NONPOPART)||(indiv->next_cascade_event==CASCADEEVENT_VS_POPART)){
                /* Only record if this is the first date of viral suppression. */
                if (indiv->PANGEA_date_startfirstVLsuppression<0)
                    indiv->PANGEA_date_startfirstVLsuppression = t;
                debug->art_vars[p].cascade_transitions[indiv->ART_status+1][LTART_VS+1]++;
                if (indiv->ART_status==LTART_VU)
                    indiv->DEBUG_cumulative_time_on_ART_VU += (t-indiv->DEBUG_time_of_last_cascade_event);
                indiv->DEBUG_time_of_last_cascade_event = t;
                virally_suppressed_process(indiv,patch[p].param,t, patch[p].cascade_events, patch[p].n_cascade_events, patch[p].size_cascade_events, patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression);
            }
            else if ((indiv->next_cascade_event==CASCADEEVENT_VU_NONPOPART)||(indiv->next_cascade_event==CASCADEEVENT_VU_POPART)){
                debug->art_vars[p].cascade_transitions[indiv->ART_status+1][LTART_VU+1]++;
                if (indiv->ART_status==LTART_VS)
                    indiv->DEBUG_cumulative_time_on_ART_VS += (t-indiv->DEBUG_time_of_last_cascade_event);
                indiv->DEBUG_time_of_last_cascade_event = t;
                virally_unsuppressed_process(indiv,patch[p].param,t, patch[p].cascade_events, patch[p].n_cascade_events, patch[p].size_cascade_events, patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression, patch[p].cumulative_outputs, patch[p].calendar_outputs);
            }
            //get_virally_unsuppressed(indiv,param,t, cascade_events, n_cascade_events, size_cascade_events);
            else if ((indiv->next_cascade_event==CASCADEEVENT_DROPOUT_NONPOPART)||(indiv->next_cascade_event==CASCADEEVENT_DROPOUT_POPART)){
                debug->art_vars[p].cascade_transitions[indiv->ART_status+1][ARTDROPOUT+1]++;
                if (indiv->ART_status==LTART_VS)
                    indiv->DEBUG_cumulative_time_on_ART_VS += (t-indiv->DEBUG_time_of_last_cascade_event);
                else if (indiv->ART_status==LTART_VU)
                    indiv->DEBUG_cumulative_time_on_ART_VU += (t-indiv->DEBUG_time_of_last_cascade_event);
                dropout_process(indiv,patch[p].param,t, patch[p].cascade_events, patch[p].n_cascade_events, patch[p].size_cascade_events, patch[p].hiv_pos_progression, patch[p].n_hiv_pos_progression, patch[p].size_hiv_pos_progression, patch[p].cumulative_outputs, patch[p].calendar_outputs);
            }
            else{
                printf("ERROR: Unknown cascade event %i for id=%li in patch %i. Exiting.\n",indiv->next_cascade_event,indiv->id,p);
                printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
                fflush(stdout);
                exit(1);
            }

        }
    }
    /* We don't need this any more so free the memory: */
    //free(cascade_events[array_index_for_cascade_event]);

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Functions for PANGEA:
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PANGEA_get_cd4(individual* indiv, double t){
    if ((indiv->PANGEA_t_prev_cd4stage<0) || (indiv->PANGEA_t_next_cd4stage<0) || (indiv->HIV_status==UNINFECTED)){
        printf("ERROR: SHould not be in get_cd4 %li %f %f %i\n",indiv->id,indiv->PANGEA_t_prev_cd4stage,indiv->PANGEA_t_next_cd4stage,indiv->HIV_status==UNINFECTED);
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }


    double currentcd4,cd4start,cd4end;

    /* Get cd4 range for current stage of CD4: */
    if (indiv->cd4==0){
        cd4start = 800;
        cd4end = 500;
    }
    else if (indiv->cd4==1){
        cd4start = 500;
        cd4end = 350;
    }
    else if (indiv->cd4==2){
        cd4start = 350;
        cd4end = 200;
    }
    else if (indiv->cd4==3){
        cd4start = 200;
        cd4end = 0;
    }
    else{
        printf("ERROR: Unknown CD4!\n");
        printf("LINE %d; FILE %s\n", __LINE__, __FILE__);
        fflush(stdout);
        exit(1);
    }

    if (indiv->HIV_status==CHRONIC){

        double f = (t - indiv->PANGEA_t_prev_cd4stage) / (indiv->PANGEA_t_next_cd4stage - indiv->PANGEA_t_prev_cd4stage);

        currentcd4 = cd4start - (cd4start-cd4end) * f;
        //printf("cd4start=%f cd4end=%f f=%f currentcd4=%f t=%f %f %f\n",cd4start,cd4end,f,currentcd4,t,indiv->PANGEA_t_prev_cd4stage,indiv->PANGEA_t_next_cd4stage);
    }
    else{
        /* Assume if in acute that current CD4 is the upper end of the given CD4 compartment. */
        currentcd4 = cd4start;
    }

    return currentcd4;
}


