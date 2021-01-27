# Table: PopART-IBM parameter dictionary for HIV file 
| Name | Value | Description | File | Source | 
|  ---- | ---- | ---- | ---- | ---- |
| `p_child_circ` | 0.0618 | Proportion of male children who have undergone traditional circumcision | HIV | - |
| `eff_circ_vmmc` | 0.6 | Effectiveness of voluntary male medical circumcision in reducing susceptibility to HIV infection | HIV | - |
| `eff_circ_tmc` | 0.0 | Effectiveness of traditional male circumcision in reducing susceptibility to HIV infection | HIV | - |
| `rr_circ_unhealed` | 0.33 | Relative risk (compared to uncircumcised state) of acquiring HIV during VMMC healing phase | HIV | - |
| `t0_pmtct` | 1990.0 | Legacy parameter | HIV | - |
| `t50_pmtct` | 2000.0 | Legacy parameter | HIV | - |
| `average_log_viral_load` | 4.0 | Legacy parameter | HIV | - |
| `average_annual_hazard` | 0.12 | Annual HIV female-to-male transmission hazard in individuals with maximal SPVL and in CD4 stage >500 | HIV | - |
| `RRacute_trans` | 5.3 | Relative infectivity of acute and early HIV infection (AEHI) stage (compared to chronic HIV infection with CD4>500) | HIV | - |
| `RRmale_to_female_trans` | 2.0 | Relative infectivity of male-to-female HIV transmission (compared to female-to-male) | HIV | - |
| `RRCD4_1` | 1.0 | Relative infectivity in CD4 stage >500 (compared to chronic HIV infection with CD4>500) | HIV | - |
| `RRCD4_2` | 1.0 | Relative infectivity in CD4 stage 350-500 (compared to chronic HIV infection with CD4>500) | HIV | - |
| `RRCD4_3` | 1.0 | Relative infectivity in CD4 stage 200-350 (compared to chronic HIV infection with CD4>500) | HIV | - |
| `RRCD4_4` | 2.34 | Relative infectivity in CD4 stage <=200 (compared to chronic HIV infection with CD4>500) | HIV | - |
| `SPVL_beta_k` | 1.02 | Hill coefficient (the steepness of the increase in infectiousness as a function of set-point viral load) | HIV | - |
| `SPVL_beta_50` | 13938.0 | Set-point viral load at which infectiousness is half its maximum value | HIV | - |
| `ART_VS_INITIAL` | 0.5 | Relative infectivity of early ART stage (compared to no ART) | HIV | - |
| `ART_VS_EFFECT` | 1.0 | Relative infectivity of virally suppressed ART stage (compared to no ART) | HIV | - |
| `ART_VU_EFFECT` | 0.3 | Relative infectivity of virally unsuppressed ART stage (compared to no ART) | HIV | - |
| `Dur_acute_range_min` | 0.08 | Minimum duration of AEHI stage in years | HIV | - |
| `Dur_acute_range_max` | 0.25 | Maximum duration of AEHI stage in years | HIV | - |
| `p_initial_cd4_gt500_spvl0` | 0.8640000000000001 | Proportion of individuals with log10 SPVL <=4.0 who have initial CD4 >500 following HIV infection | HIV | - |
| `p_initial_cd4_350_500_spvl0` | 0.113 | Proportion of individuals with log10 SPVL <=4.0 who have initial CD4 350-500 following HIV infection | HIV | - |
| `p_initial_cd4_200_350_spvl0` | 0.023 | Proportion of individuals with log10 SPVL <=4.0 who have initial CD4 200-350 following HIV infection | HIV | - |
| `p_initial_cd4_lt200_spvl0` | 0.0 | Proportion of individuals with log10 SPVL <=4.0 who have initial CD4 <=200 following HIV infection | HIV | - |
| `p_initial_cd4_gt500_spvl1` | 0.78 | Proportion of individuals with log10 SPVL 4.0-4.5 who have initial CD4 >500 following HIV infection | HIV | - |
| `p_initial_cd4_350_500_spvl1` | 0.19 | Proportion of individuals with log10 SPVL 4.0-4.5 who have initial CD4 350-500 following HIV infection | HIV | - |
| `p_initial_cd4_200_350_spvl1` | 0.03 | Proportion of individuals with log10 SPVL 4.0-4.5 who have initial CD4 200-350 following HIV infection | HIV | - |
| `p_initial_cd4_lt200_spvl1` | 0.0 | Proportion of individuals with log10 SPVL 4.0-4.5 who have initial CD4 <=200 following HIV infection | HIV | - |
| `p_initial_cd4_gt500_spvl2` | 0.74 | Proportion of individuals with log10 SPVL 4.5-5.0 who have initial CD4 >500 following HIV infection | HIV | - |
| `p_initial_cd4_350_500_spvl2` | 0.21 | Proportion of individuals with log10 SPVL 4.5-5.0 who have initial CD4 350-500 following HIV infection | HIV | - |
| `p_initial_cd4_200_350_spvl2` | 0.05 | Proportion of individuals with log10 SPVL 4.5-5.0 who have initial CD4 200-350 following HIV infection | HIV | - |
| `p_initial_cd4_lt200_spvl2` | 0.0 | Proportion of individuals with log10 SPVL 4.5-5.0 who have initial CD4 <=200 following HIV infection | HIV | - |
| `p_initial_cd4_gt500_spvl3` | 0.71 | Proportion of individuals with log10 SPVL >5.0 who have initial CD4 >500 following HIV infection | HIV | - |
| `p_initial_cd4_350_500_spvl3` | 0.25 | Proportion of individuals with log10 SPVL >5.0 who have initial CD4 350-500 following HIV infection | HIV | - |
| `p_initial_cd4_200_350_spvl3` | 0.04 | Proportion of individuals with log10 SPVL >5.0 who have initial CD4 200-350 following HIV infection | HIV | - |
| `p_initial_cd4_lt200_spvl3` | 0.0 | Proportion of individuals with log10 SPVL >5.0 who have initial CD4 <=200 following HIV infection | HIV | - |
| `initial_SPVL_mu` | 4.74 | Mean SPVL in seeded HIV infections | HIV | - |
| `initial_SPVL_sigma` | 0.61 | Standard deviation of SPVL in seeded HIV infections | HIV | - |
| `SPVL_sigma_M` | 0.0 | Parameter used when SPVL inheritance is switched on (macro SPVL_INHERITANCE=1). By default SPVL inheritance is off | HIV | - |
| `SPVL_sigma_E` | 0.61 | Standard deviation of SPVL in non-seeded infections. | HIV | - |
| `p_misclassify_cd4_0_0` | 1.0 | Probability of classifying actual CD4 category >500 as >500 | HIV | - |
| `p_misclassify_cd4_0_1` | 0.0 | Probability of classifying actual CD4 category >500 as 350-500 | HIV | - |
| `p_misclassify_cd4_0_2` | 0.0 | Probability of classifying actual CD4 category >500 as 200-350 | HIV | - |
| `p_misclassify_cd4_0_3` | 0.0 | Probability of classifying actual CD4 category >500 as <=200 | HIV | - |
| `p_misclassify_cd4_1_0` | 0.0 | Probability of classifying actual CD4 category 350-500 as >500 | HIV | - |
| `p_misclassify_cd4_1_1` | 1.0 | Probability of classifying actual CD4 category 350-500 as 350-500 | HIV | - |
| `p_misclassify_cd4_1_2` | 0.0 | Probability of classifying actual CD4 category 350-500 as 200-350 | HIV | - |
| `p_misclassify_cd4_1_3` | 0.0 | Probability of classifying actual CD4 category 350-500 as <=200 | HIV | - |
| `p_misclassify_cd4_2_0` | 0.0 | Probability of classifying actual CD4 category 200-350 as >500 | HIV | - |
| `p_misclassify_cd4_2_1` | 0.0 | Probability of classifying actual CD4 category 200-350 as 350-500 | HIV | - |
| `p_misclassify_cd4_2_2` | 1.0 | Probability of classifying actual CD4 category 200-350 as 200-350 | HIV | - |
| `p_misclassify_cd4_2_3` | 0.0 | Probability of classifying actual CD4 category 200-350 as <=200 | HIV | - |
| `p_misclassify_cd4_3_0` | 0.0 | Probability of classifying actual CD4 category <=200 as >500 | HIV | - |
| `p_misclassify_cd4_3_1` | 0.0 | Probability of classifying actual CD4 category <=200 as 350-500 | HIV | - |
| `p_misclassify_cd4_3_2` | 0.0 | Probability of classifying actual CD4 category <=200 as 200-350 | HIV | - |
| `p_misclassify_cd4_3_3` | 1.0 | Probability of classifying actual CD4 category <=200 as <=200 | HIV | - |
| `time_gt500_to_500_spvl0` | 5.35 | Time (in years) spent in CD4 category >500 for individuals with log10 SPVL <=4.0 when not on ART | HIV | - |
| `time_350to500_to_350_spvl0` | 3.66 | Time (in years) spent in CD4 category 350-500 for individuals with log10 SPVL <=4.0 not on ART | HIV | - |
| `time_200to350_to_200_spvl0` | 7.62 | Time (in years) spent in CD4 category 200-350 for individuals with log10 SPVL <=4.0 not on ART | HIV | - |
| `time_lt200_to_death_spvl0` | 2.3 | Time (in years) spent in CD4 category <=200 until death for individuals with log10 SPVL <=4.0 not on ART | HIV | - |
| `time_gt500_to_500_spvl1` | 3.12 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category >500 for individuals with log10 SPVL 4-4.5 when not on ART | HIV | - |
| `time_350to500_to_350_spvl1` | 3.09 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category 350-500 for individuals with log10 SPVL 4-4.5 when not on ART | HIV | - |
| `time_200to350_to_200_spvl1` | 8.39 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category 200-350 for individuals with log10 SPVL 4-4.5 when not on ART | HIV | - |
| `time_lt200_to_death_spvl1` | 2.3 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category <=200 for individuals with log10 SPVL 4-4.5 when not on ART | HIV | - |
| `time_gt500_to_500_spvl2` | 2.35 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category >500 for individuals with log10 SPVL 4.5-5 when not on ART | HIV | - |
| `time_350to500_to_350_spvl2` | 2.32 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category 350-500 for individuals with log10 SPVL 4.5-5 when not on ART | HIV | - |
| `time_200to350_to_200_spvl2` | 6.57 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category 200-350 for individuals with log10 SPVL 4.5-5 when not on ART | HIV | - |
| `time_lt200_to_death_spvl2` | 2.3 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category <=200 for individuals with log10 SPVL 4.5-5 when not on ART | HIV | - |
| `time_gt500_to_500_spvl3` | 1.51 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category >500 for individuals with log10 SPVL >5 when not on ART | HIV | - |
| `time_350to500_to_350_spvl3` | 1.44 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category 350-500 for individuals with log10 SPVL >5 when not on ART | HIV | - |
| `time_200to350_to_200_spvl3` | 2.93 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category 200-350 for individuals with log10 SPVL >5 when not on ART | HIV | - |
| `time_lt200_to_death_spvl3` | 2.3 | Used when macro USEFOURSPVLCATEGORIES=1 for time spent (in years) in CD4 category <=200 for individuals with log10 SPVL >5 when not on ART | HIV | - |
| `cox_RRspvl_gt500_to_500` | 2.17 | Used when macro USEFOURSPVLCATEGORIES=0 (default setting). Factor by which time in CD4 category >500 is decreased, per 10-fold increase in SPVL | HIV | - |
| `cox_RRspvl_500_to_350` | 1.88 | Used when macro USEFOURSPVLCATEGORIES=0 (default setting). Factor by which time in CD4 category 350-500 is decreased, per 10-fold increase in SPVL | HIV | - |
| `cox_RRspvl_350_to_200` | 1.96 | Used when macro USEFOURSPVLCATEGORIES=0 (default setting). Factor by which time in CD4 category 200-350 is decreased, per 10-fold increase in SPVL | HIV | - |
| `cox_RRspvl_200_to_death` | 1.63 | Used when macro USEFOURSPVLCATEGORIES=0 (default setting). Factor by which time in CD4 category <=200 is decreased, per 10-fold increase in SPVL | HIV | - |
| `factor_for_slower_progression_ART_VU` | 1.0 | Multiplier for increased duration in each CD4 stage when on ART but virally unsuppressed (compared to being virally suppressed) | HIV | - |