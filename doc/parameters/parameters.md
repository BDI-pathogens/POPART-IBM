# Table: PopART-IBM parameter dictionary
| Name | Value | Description | File | Source | 
|  ---- | ---- | ---- | ---- | ---- |
| `time_to_background_HIVtestNOW` | 3.0 | Legacy parameter for different background testing function | cascade | - |
| `time_to_background_HIVtest_maxval` | 8.0 | Legacy parameter for different background testing function | cascade | - |
| `time_to_background_HIVtest_exponent` | 6.0 | Legacy parameter for different background testing function | cascade | - |
| `time_to_background_HIVtest_midpoint` | 6.5 | Legacy parameter for different background testing function | cascade | - |
| `p_HIV_background_testing_female_pre2006` | 0.5 | Baseline probability of a female having an HIV test in period from when HIV tests were first offered until 2006 in the background cascade | cascade | - |
| `p_HIV_background_testing_female_current` | 0.5 | Annual baseline probability of a female having an HIV test from 2006 onwards in the background cascade | cascade | - |
| `RR_HIV_background_testing_male` | 0.5 | Relative probability of a man having a HIV test in the background cascade compared to a woman | cascade | - |
| `HIV_rapid_test_sensitivity_CHIPS` | 1.0 | Sensitivity of the rapid HIV test used by CHiPs | cascade | - |
| `p_collect_hiv_test_results_cd4_over200` | 0.97 | Probability collect background HIV test results when CD4>200 or HIV-negative | cascade | - |
| `p_collect_hiv_test_results_cd4_under200` | 1.0 | Probability collect background HIV test results when HIV positive with CD4<=200 | cascade | - |
| `p_collect_cd4_test_results_cd4_nonpopart` | 0.85 | Probability collect background CD4 test results | cascade | - |
| `p_collect_cd4_test_results_cd4_popartYEAR1` | 1.0 | Probability collect CD4 test results following CHiPs HIV test for CHiPs round 1 | cascade | - |
| `p_collect_cd4_test_results_cd4_popartYEAR2onwards` | 1.0 | Probability collect CD4 test results following CHiPs HIV test for CHiPs round 2 onwards | cascade | - |
| `p_dies_earlyart_cd4_over500` | 0.0 | Probability die while on early ART and CD4>500 | cascade | - |
| `p_dies_earlyart_cd4_350_500` | 0.0 | Probability die while on early ART and CD4 350-500 | cascade | - |
| `p_dies_earlyart_cd4_200_350` | 0.0 | Probability die while on early ART and CD4 200-350 | cascade | - |
| `p_dies_earlyart_cd4_under200` | 0.08 | Probability die while on early ART and CD4<=200 | cascade | - |
| `p_leaves_earlyart_cd4_over200_if_not_die_early` | 0.0 | Probability drop out from early ART if CD4>200, given that do not die on early ART | cascade | - |
| `p_leaves_earlyart_cd4_under200_if_not_die_early` | 0.0 | Probability drop out from early ART if CD4<=200, given that do not die on early ART | cascade | - |
| `p_becomes_vs_after_earlyart_if_not_die_early_or_leave` | 0.9 | Probability become virally suppressed after early ART, given that do not die on early ART | cascade | - |
| `p_stays_virally_suppressed` | 0.6 | Probability that a woman who becomes virally suppressed remains virally suppressed for life | cascade | - |
| `p_stays_virally_suppressed_male` | 0.5 | Relative probability of a man who becomes virally suppressed remaining virally suppressed for life (compared to a woman) | cascade | - |
| `p_stops_virally_suppressed` | 0.1 | Probability that someone who becomes virally suppressed will eventually become virally unsuppressed (but remains on ART) | cascade | - |
| `p_stays_cabo` | 0.1 | Probability that targeted group individuals who becomes virally suppressed but wants to drop out will take cabotegravir medicine (no transmission, no drop-out) | cascade | - |
| `p_vu_becomes_virally_suppressed` | 0.0 | Probability that someone who becomes virally unsuppressed after early ART eventually becomes virally suppressed | cascade | - |
| `t_earlyart_dropout_nonpopart_min` | 0.0 | Minimum time (in years) from starting early ART (from background HIV test) to dropping out, for those who drop out | cascade | - |
| `t_earlyart_dropout_nonpopart_max` | 0.17 | Maximum time (in years) from starting early ART (from background HIV test) to dropping out, for those who drop out | cascade | - |
| `t_earlyart_dropout_popart_min` | 0.0 | Minimum time (in years) from starting early ART (from CHiPs HIV test) to dropping out, for those who drop out | cascade | - |
| `t_earlyart_dropout_popart_max` | 0.17 | Maximum time (in years) from starting early ART (from CHiPs HIV test) to dropping out, for those who drop out | cascade | - |
| `t_dies_earlyart_nonpopart_min` | 0.0 | Minimum time (in years) from starting early ART (from background HIV test) to dying, for those who die during early ART | cascade | - |
| `t_dies_earlyart_nonpopart_max` | 0.17 | Maximum time (in years) from starting early ART (from background HIV test) to dying, for those who die during early ART | cascade | - |
| `t_dies_earlyart_popart_min` | 0.0 | Minimum time (in years) from starting early ART (from CHiPs HIV test) to dying, for those who die during early ART | cascade | - |
| `t_dies_earlyart_popart_max` | 0.17 | Maximum time (in years) from starting early ART (from CHiPs HIV test) to dying, for those who die during early ART | cascade | - |
| `t_end_early_art` | 0.17 | Maximum duration (in years) of early ART stage (for those who will remain on ART (whether virally suppressed or unsuppressed) | cascade | - |
| `t_cd4_retest_nonpopart_min` | 0.9 | Minimum time between successive CD4 tests when not initially eligible for ART following positive background HIV test result | cascade | - |
| `t_cd4_retest_nonpopart_max` | 1.1 | Maximum time between successive CD4 tests when not initially eligible for ART following positive background HIV test result | cascade | - |
| `t_cd4_retest_popart_min` | 0.9 | Minimum time between successive CD4 tests when not initially eligible for ART following positive CHiPs HIV test result (arm B) | cascade | - |
| `t_cd4_retest_popart_max` | 1.1 | Maximum time between successive CD4 tests when not initially eligible for ART following positive CHiPs HIV test result (arm B) | cascade | - |
| `t_cd4_whenartfirstavail_min` | 0.0 | Minimum time (in years) until someone who is HIV-positive aware (and wants to start ART) has a CD4 test to determine eligibility for ART, when ART first becomes available | cascade | - |
| `t_cd4_whenartfirstavail_max` | 2.0 | Maximum time (in years) until someone who is HIV-positive aware (and wants to start ART) has a CD4 test to determine eligibility for ART, when ART first becomes available | cascade | - |
| `t_delay_hivtest_to_cd4test_nonpopart_min` | 0.02083333 | Minimum time in years from (positive) background HIV test to CD4 test | cascade | - |
| `t_delay_hivtest_to_cd4test_nonpopart_max` | 0.5 | Maximum time in years from (positive) background HIV test to CD4 test | cascade | - |
| `t_delay_hivtest_to_cd4test_popart_min` | 0.0208333 | Minimum time in years from (positive) CHiPs HIV test to CD4 test | cascade | - |
| `t_delay_hivtest_to_cd4test_popart_max` | 0.08 | Maximum time in years from (positive) CHiPs HIV test to CD4 test | cascade | - |
| `t_start_art_nonpopart_mean` | 0.55 | Mean time in years to start ART through background HIV testing (of those who decide to start ART) | cascade | - |
| `n_time_periods_art_popart_round_1` | 6.0 | Number of time periods in R1 with different distributions of time to start ART following a positive CHiPs test (see Supplement Table S4.5) | cascade | - |
| `n_time_periods_art_popart_round_2` | 1.0 | Number of time periods in R2 with different distributions of time to start ART following a positive CHiPs test | cascade | - |
| `n_time_periods_art_popart_round_3` | 1.0 | Number of time periods in R3 with different distributions of time to start ART following a positive CHiPs test | cascade | - |
| `t_start_art_popart_round1_mean_fast_1` | 0.198 | Mean time to start ART fast (years) in Round 1 Period 1 | cascade | - |
| `t_start_art_popart_round1_mean_fast_2` | 0.161 | Mean time to start ART fast (years) in Round 1 Period 2 | cascade | - |
| `t_start_art_popart_round1_mean_fast_3` | 0.11800000000000001 | Mean time to start ART fast (years) in Round 1 Period 3 | cascade | - |
| `t_start_art_popart_round1_mean_fast_4` | 0.135 | Mean time to start ART fast (years) in Round 1 Period 4 | cascade | - |
| `t_start_art_popart_round1_mean_fast_5` | 0.081 | Mean time to start ART fast (years) in Round 1 Period 5 | cascade | - |
| `t_start_art_popart_round1_mean_fast_6` | 0.027999999999999997 | Mean time to start ART fast (years) in Round 1 Period 6 | cascade | - |
| `t_start_art_popart_round1_mean_slow_1` | 2.7110000000000003 | Mean time to start ART slow (years) in Round 1 Period 1 | cascade | - |
| `t_start_art_popart_round1_mean_slow_2` | 2.378 | Mean time to start ART slow (years) in Round 1 Period 2 | cascade | - |
| `t_start_art_popart_round1_mean_slow_3` | 2.086 | Mean time to start ART slow (years) in Round 1 Period 3 | cascade | - |
| `t_start_art_popart_round1_mean_slow_4` | 1.62 | Mean time to start ART slow (years) in Round 1 Period 4 | cascade | - |
| `t_start_art_popart_round1_mean_slow_5` | 1.2819999999999998 | Mean time to start ART slow (years) in Round 1 Period 5 | cascade | - |
| `t_start_art_popart_round1_mean_slow_6` | 0.8540000000000001 | Mean time to start ART slow (years) in Round 1 Period 6 | cascade | - |
| `p_start_art_popart_round1_mean_fast_1` | 0.215 | Probability of being a fast starter in Round 1 Period 1 | cascade | - |
| `p_start_art_popart_round1_mean_fast_2` | 0.18600000000000003 | Probability of being a fast starter in Round 1 Period 2 | cascade | - |
| `p_start_art_popart_round1_mean_fast_3` | 0.18 | Probability of being a fast starter in Round 1 Period 3 | cascade | - |
| `p_start_art_popart_round1_mean_fast_4` | 0.2 | Probability of being a fast starter in Round 1 Period 4 | cascade | - |
| `p_start_art_popart_round1_mean_fast_5` | 0.313 | Probability of being a fast starter in Round 1 Period 5 | cascade | - |
| `p_start_art_popart_round1_mean_fast_6` | 0.175 | Probability of being a fast starter in Round 1 Period 6 | cascade | - |
| `t_start_art_popart_round2_mean_fast` | 0.075 | Mean time to start ART fast (years) in Round 2 | cascade | - |
| `t_start_art_popart_round2_mean_slow` | 1.262 | Mean time to start ART slow (years) in Round 2 | cascade | - |
| `p_start_art_popart_round2_mean_fast` | 0.226 | Probability of being a fast starter in Round 2 | cascade | - |
| `t_start_art_popart_round3_mean_fast` | 0.03 | Mean time to start ART fast (years) in Round 3 onwards | cascade | - |
| `t_start_art_popart_round3_mean_slow` | 0.795 | Mean time to start ART slow (years) in Round 3 onwards | cascade | - |
| `p_start_art_popart_round3_mean_fast` | 0.284 | Probability of being a fast starter in Round 3 onwards | cascade | - |
| `t_end_vs_becomevu_nonpopart_min` | 0.01 | Minimum time taken (in years) for someone VS and who started ART through background cascade (and who will eventually become VU) to become VU | cascade | - |
| `t_end_vs_becomevu_nonpopart_max` | 6.0 | Maximum time taken (in years) for someone VS and who started ART through background cascade (and who will eventually become VU) to become VU | cascade | - |
| `t_end_vs_becomevu_popart_min` | 0.01 | Minimum time taken (in years) for someone VS and who started ART through CHiPs testing (and who will eventually become VU) to become VU | cascade | - |
| `t_end_vs_becomevu_popart_max` | 6.0 | Maximum time taken (in years) for someone VS and who started ART through CHiPs testing (and who will eventually become VU) to become VU | cascade | - |
| `t_end_vs_dropout_nonpopart_min` | 0.01 | Minimum time taken (in years) for someone VS and who started ART through background cascade (and who will eventually drop out) to drop out | cascade | - |
| `t_end_vs_dropout_nonpopart_max` | 6.0 | Maximum time taken (in years) for someone VS and who started ART through background cascade (and who will eventually drop out) to drop out | cascade | - |
| `t_end_vs_dropout_popart_min` | 0.01 | Minimum time taken (in years) for someone VS and who started ART through CHiPs testing (and who will eventually drop out) to drop out | cascade | - |
| `t_end_vs_dropout_popart_max` | 6.0 | Maximum time taken (in years) for someone VS and who started ART through CHiPs testing (and who will eventually drop out) to drop out | cascade | - |
| `t_end_vu_becomevs_nonpopart_min` | 0.01 | Minimum time taken (in years) for someone VU and who started ART through background cascade (and who will eventually become VS) to become VS. Currently zero probability, as determined by p_vu_becomes_virally_suppressed | cascade | - |
| `t_end_vu_becomevs_nonpopart_max` | 6.0 | Maximum time taken (in years) for someone VU and who started ART through background cascade (and who will eventually become VS) to become VS. Currently zero probability, as determined by p_vu_becomes_virally_suppressed | cascade | - |
| `t_end_vu_becomevs_popart_min` | 0.01 | Minimum time taken (in years) for someone VU and who started ART through CHiPs testing (and who will eventually become VS) to become VS. Currently zero probability, as determined by p_vu_becomes_virally_suppressed | cascade | - |
| `t_end_vu_becomevs_popart_max` | 6.0 | Maximum time taken (in years) for someone VU and who started ART through CHiPs testing (and who will eventually become VS) to become VS. Currently zero probability, as determined by p_vu_becomes_virally_suppressed | cascade | - |
| `t_end_vu_dropout_nonpopart_min` | 0.01 | Minimum time taken (in years) for someone VU and who started ART through background cascade (and who does not eventually become VS) to drop out | cascade | - |
| `t_end_vu_dropout_nonpopart_max` | 6.0 | Maximum time taken (in years) for someone VU and who started ART through background cascade (and who does not eventually become VS) to drop out | cascade | - |
| `t_end_vu_dropout_popart_min` | 0.01 | Minimum time taken (in years) for someone VU and who started ART through CHiPs testing (and who does not eventually become VS) to drop out | cascade | - |
| `t_end_vu_dropout_popart_max` | 6.0 | Maximum time taken (in years) for someone VU and who started ART through CHiPs testing (and who does not eventually become VS) to drop out | cascade | - |
| `p_popart_to_cascade` | 1.0 | Probability that someone who had previously dropped out (either tested HIV+ but never linked to care, or was on ART and dropped out) will agree to PopART | cascade | - |
| `p_circ_nopopart` | 0.4 | Probability of man being circumcised following negative background HIV test result | cascade | - |
| `p_circ_popart` | 0.4 | Probability of man being circumcised following negative CHiPs HIV test result | cascade | - |
| `t_get_vmmc_nonpopart_min` | 0.25 | Minimum time (in years) between background HIV- test and getting VMMC | cascade | - |
| `t_get_vmmc_nonpopart_max` | 1.0 | Maximum time (in years) between background HIV- test and getting VMMC | cascade | - |
| `t_get_vmmc_popart_min` | 0.08333 | Minimum time (in years) between CHiPs HIV- test and getting VMMC | cascade | - |
| `t_get_vmmc_popart_max` | 1.0 | Maximum time (in years) between CHiPs HIV- test and getting VMMC | cascade | - |
| `t_vmmc_healing` | 0.038 | Time (in years) for VMMC wound to heal | cascade | - |
| `sex_ratio` | 0.5073892 | Fraction of new adults that are male | demographics | - |
| `sexual_worker_ratio` | 0.00663 | Fraction of new adults that are sexual worker | demographics | - |
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
| `rng_seed` | 2020.0 | Random starting seed | init | - |
| `initial_adult_population_size` | 2050.0 | Size of the population aged 14+ at the start of the simulation | init | - |
| `initial_prop_14_17` | 0.1583 | Proportion of the modelled population at the start of the simulation that are aged 14-17 | init | - |
| `initial_prop_18_22` | 0.1613 | Proportion of the modelled population at the start of the simulation that are aged 18-22 | init | - |
| `initial_prop_23_30` | 0.1885 | Proportion of the modelled population at the start of the simulation that are aged 23-30 | init | - |
| `initial_prop_31_40` | 0.1983 | Proportion of the modelled population at the start of the simulation that are aged 31-40 | init | - |
| `initial_prop_41_50` | 0.1353 | Proportion of the modelled population at the start of the simulation that are aged 41-50 | init | - |
| `initial_prop_51_60` | 0.0828 | Proportion of the modelled population at the start of the simulation that are aged 51-60 | init | - |
| `initial_prop_61_and_older` | 0.0755 | Proportion of the modelled population at the start of the simulation that are aged 61+ | init | - |
| `initial_low_risk_male` | 0.45 | Proportion of men who are in the low-risk activity group when entering adult population | init | - |
| `initial_low_risk_female` | 0.45 | Proportion of women who are in the low-risk activity group when entering adult population | init | - |
| `initial_med_not_low_risk_male` | 0.745 | Of those men who did not enter the low risk group, proportion of men who are in the medium-risk activity group when entering adult population | init | - |
| `initial_med_not_low_risk_female` | 0.745 | Of those women who did not enter the low risk group, proportion of women who are in the medium-risk activity group when entering adult population | init | - |
| `initial_prop_HIV_pos_low_risk_male` | 2e-05 | Unscaled % of low activity level male population seeded HIV+ each year | init | - |
| `initial_prop_HIV_pos_low_risk_female` | 2e-05 | Unscaled % of low activity level female population seeded HIV+ each year | init | - |
| `initial_prop_HIV_pos_med_risk_male` | 5e-05 | Unscaled % of medium activity level male population seeded HIV+ each year | init | - |
| `initial_prop_HIV_pos_med_risk_female` | 5e-05 | Unscaled % of medium activity level female population seeded HIV+ each year | init | - |
| `initial_prop_HIV_pos_high_risk_male` | 8e-05 | Unscaled % of high activity level male population seeded HIV+ each year | init | - |
| `initial_prop_HIV_pos_high_risk_female` | 8e-05 | Unscaled % of high activity level female population seeded HIV+ each year | init | - |
| `log_seed_multiplier` | 1.0 | log10 of factor multiplying seeded % at the start of the HIV epidemic | init | - |
| `n_years_HIV_seeding` | 5.0 | Number of years of HIV seeding after start of HIV epidemic | init | - |
| `assortativity` | 0.5 | Risk assortativity | partnerships | - |
| `prop_compromise_from_males` | 0.5 | theta, the proportion of compromise in number of sexual partners from males | partnerships | - |
| `c_f_1` | 0.011822581 | Within-community annual partnership formation rates for women in age group 14-17 | partnerships | - |
| `c_f_2` | 0.023645162 | Within-community annual partnership formation rates for women in age group 18-22 | partnerships | - |
| `c_f_3` | 0.018226193 | Within-community annual partnership formation rates for women in age group 23-30 | partnerships | - |
| `c_f_4` | 0.012009144 | Within-community annual partnership formation rates for women in age group 31-40 | partnerships | - |
| `c_f_5` | 0.009834833000000001 | Within-community annual partnership formation rates for women in age group 41-50 | partnerships | - |
| `c_f_6` | 0.004917416 | Within-community annual partnership formation rates for women in age group 51-60 | partnerships | - |
| `c_f_7` | 0.002458708 | Within-community annual partnership formation rates for women in age group 61+ | partnerships | - |
| `c_m_1` | 0.021872336000000003 | Within-community annual partnership formation rates for men in age group 14-17 | partnerships | - |
| `c_m_2` | 0.043744672000000005 | Within-community annual partnership formation rates for men in age group 18-22 | partnerships | - |
| `c_m_3` | 0.039466657 | Within-community annual partnership formation rates for men in age group 23-30 | partnerships | - |
| `c_m_4` | 0.023920296 | Within-community annual partnership formation rates for men in age group 31-40 | partnerships | - |
| `c_m_5` | 0.018824023 | Within-community annual partnership formation rates for men in age group 41-50 | partnerships | - |
| `c_m_6` | 0.006004572 | Within-community annual partnership formation rates for men in age group 51-60 | partnerships | - |
| `c_m_7` | 0.003002286 | Within-community annual partnership formation rates for men in age group 61+ | partnerships | - |
| `c_multiplier` | 2.25 | Multiplier to  for mis-reporting of number of sexual partners | partnerships | - |
| `rel_rate_partnership_formation_between_patches` | 0.45 | Relative rate of formation of partnerships between patches compared to within patches | partnerships | - |
| `rr_hiv_between_vs_within_patch` | 0.3 | Relative HIV transmission risk for partnerships between patches, compared to partnerships within same patch | partnerships | - |
| `relative_number_partnerships_per_risk_low` | 1.0 | Relative number of partnerships for low activity level group (compared to low activity level group) | partnerships | - |
| `relative_number_partnerships_per_risk_med` | 7.7 | Relative number of partnerships for medium activity level group (compared to low activity level group) | partnerships | - |
| `relative_number_partnerships_per_risk_high` | 18.6 | Relative number of partnerships for high activity level group (compared to low activity level group) | partnerships | - |
| `p_age_m_1_1` | 0.935028249 | Age mixing matrix for men | partnerships | - |
| `p_age_m_1_2` | 1.0 | - | partnerships | - |
| `p_age_m_1_3` | 0.0 | - | partnerships | - |
| `p_age_m_1_4` | 0.0 | - | partnerships | - |
| `p_age_m_1_5` | 0.0 | - | partnerships | - |
| `p_age_m_1_6` | 0.0 | - | partnerships | - |
| `p_age_m_1_7` | 0.0 | - | partnerships | - |
| `p_age_m_2_1` | 0.26859504100000003 | - | partnerships | - |
| `p_age_m_2_2` | 0.935028249 | - | partnerships | - |
| `p_age_m_2_3` | 1.0 | - | partnerships | - |
| `p_age_m_2_4` | 0.0 | - | partnerships | - |
| `p_age_m_2_5` | 0.0 | - | partnerships | - |
| `p_age_m_2_6` | 0.0 | - | partnerships | - |
| `p_age_m_2_7` | 0.0 | - | partnerships | - |
| `p_age_m_3_1` | 0.057416268 | - | partnerships | - |
| `p_age_m_3_2` | 0.640862944 | - | partnerships | - |
| `p_age_m_3_3` | 0.9434628979999999 | - | partnerships | - |
| `p_age_m_3_4` | 0.875 | - | partnerships | - |
| `p_age_m_3_5` | 1.0 | - | partnerships | - |
| `p_age_m_3_6` | 0.0 | - | partnerships | - |
| `p_age_m_3_7` | 0.0 | - | partnerships | - |
| `p_age_m_4_1` | 0.012030075 | - | partnerships | - |
| `p_age_m_4_2` | 0.112633181 | - | partnerships | - |
| `p_age_m_4_3` | 0.615780446 | - | partnerships | - |
| `p_age_m_4_4` | 0.9776785709999999 | - | partnerships | - |
| `p_age_m_4_5` | 1.0 | - | partnerships | - |
| `p_age_m_4_6` | 0.0 | - | partnerships | - |
| `p_age_m_4_7` | 0.0 | - | partnerships | - |
| `p_age_m_5_1` | 0.0 | - | partnerships | - |
| `p_age_m_5_2` | 0.023972603 | - | partnerships | - |
| `p_age_m_5_3` | 0.09473684199999999 | - | partnerships | - |
| `p_age_m_5_4` | 0.713178295 | - | partnerships | - |
| `p_age_m_5_5` | 0.932432432 | - | partnerships | - |
| `p_age_m_5_6` | 0.8 | - | partnerships | - |
| `p_age_m_5_7` | 1.0 | - | partnerships | - |
| `p_age_m_6_1` | 0.0 | - | partnerships | - |
| `p_age_m_6_2` | 0.0 | - | partnerships | - |
| `p_age_m_6_3` | 0.023972603 | - | partnerships | - |
| `p_age_m_6_4` | 0.09473684199999999 | - | partnerships | - |
| `p_age_m_6_5` | 0.713178295 | - | partnerships | - |
| `p_age_m_6_6` | 0.932432432 | - | partnerships | - |
| `p_age_m_6_7` | 1.0 | - | partnerships | - |
| `p_age_m_7_1` | 0.0 | - | partnerships | - |
| `p_age_m_7_2` | 0.0 | - | partnerships | - |
| `p_age_m_7_3` | 0.0 | - | partnerships | - |
| `p_age_m_7_4` | 0.023972603 | - | partnerships | - |
| `p_age_m_7_5` | 0.09473684199999999 | - | partnerships | - |
| `p_age_m_7_6` | 0.713178295 | - | partnerships | - |
| `p_age_m_7_7` | 1.0 | - | partnerships | - |
| `p_age_f_1_1` | 0.113636364 | Age mixing matrix for women | partnerships | - |
| `p_age_f_1_2` | 0.778329198 | - | partnerships | - |
| `p_age_f_1_3` | 0.891791045 | - | partnerships | - |
| `p_age_f_1_4` | 0.931034483 | - | partnerships | - |
| `p_age_f_1_5` | 1.0 | - | partnerships | - |
| `p_age_f_1_6` | 0.0 | - | partnerships | - |
| `p_age_f_1_7` | 0.0 | - | partnerships | - |
| `p_age_f_2_1` | 0.0036523009999999997 | - | partnerships | - |
| `p_age_f_2_2` | 0.113636364 | - | partnerships | - |
| `p_age_f_2_3` | 0.778329198 | - | partnerships | - |
| `p_age_f_2_4` | 0.891791045 | - | partnerships | - |
| `p_age_f_2_5` | 0.931034483 | - | partnerships | - |
| `p_age_f_2_6` | 1.0 | - | partnerships | - |
| `p_age_f_2_7` | 0.0 | - | partnerships | - |
| `p_age_f_3_1` | 0.001697073 | - | partnerships | - |
| `p_age_f_3_2` | 0.005949851 | - | partnerships | - |
| `p_age_f_3_3` | 0.27746900399999996 | - | partnerships | - |
| `p_age_f_3_4` | 0.889349112 | - | partnerships | - |
| `p_age_f_3_5` | 0.930481283 | - | partnerships | - |
| `p_age_f_3_6` | 0.7692307690000001 | - | partnerships | - |
| `p_age_f_3_7` | 1.0 | - | partnerships | - |
| `p_age_f_4_1` | 0.001708672 | - | partnerships | - |
| `p_age_f_4_2` | 0.002139495 | - | partnerships | - |
| `p_age_f_4_3` | 0.00728988 | - | partnerships | - |
| `p_age_f_4_4` | 0.432829374 | - | partnerships | - |
| `p_age_f_4_5` | 0.8979436409999999 | - | partnerships | - |
| `p_age_f_4_6` | 0.9104477609999999 | - | partnerships | - |
| `p_age_f_4_7` | 1.0 | - | partnerships | - |
| `p_age_f_5_1` | 0.0027662520000000003 | - | partnerships | - |
| `p_age_f_5_2` | 0.002773925 | - | partnerships | - |
| `p_age_f_5_3` | 0.0 | - | partnerships | - |
| `p_age_f_5_4` | 0.015299026 | - | partnerships | - |
| `p_age_f_5_5` | 0.535310734 | - | partnerships | - |
| `p_age_f_5_6` | 0.8632218840000001 | - | partnerships | - |
| `p_age_f_5_7` | 1.0 | - | partnerships | - |
| `p_age_f_6_1` | 0.0 | - | partnerships | - |
| `p_age_f_6_2` | 0.0027662520000000003 | - | partnerships | - |
| `p_age_f_6_3` | 0.002773925 | - | partnerships | - |
| `p_age_f_6_4` | 0.0 | - | partnerships | - |
| `p_age_f_6_5` | 0.015299026 | - | partnerships | - |
| `p_age_f_6_6` | 0.535310734 | - | partnerships | - |
| `p_age_f_6_7` | 1.0 | - | partnerships | - |
| `p_age_f_7_1` | 0.0 | - | partnerships | - |
| `p_age_f_7_2` | 0.0 | - | partnerships | - |
| `p_age_f_7_3` | 0.0027662520000000003 | - | partnerships | - |
| `p_age_f_7_4` | 0.002773925 | - | partnerships | - |
| `p_age_f_7_5` | 0.0 | - | partnerships | - |
| `p_age_f_7_6` | 0.015299026 | - | partnerships | - |
| `p_age_f_7_7` | 1.0 | - | partnerships | - |
| `max_n_part_noage_low` | 1.0 | Maximum number of concurrent partners for individuals in low activity level | partnerships | - |
| `max_n_part_noage_med` | 3.0 | Maximum number of concurrent partners for individuals in medium activity level | partnerships | - |
| `max_n_part_noage_high` | 10.0 | Maximum number of concurrent partners for individuals in high activity level | partnerships | - |
| `breakup_scale_lambda_low_within_patch` | 13.7 | Mean duration (in years) of low activity level partnerships | partnerships | - |
| `breakup_scale_lambda_med_within_patch` | 6.9 | Mean duration (in years) of medium activity level partnerships | partnerships | - |
| `breakup_scale_lambda_high_within_patch` | 4.4 | Mean duration (in years) of high activity level partnerships | partnerships | - |
| `breakup_shape_k_low` | 1.0 | Weibull scale factor chosen to be 1 so that duration of low activity level partnerships is exponential | partnerships | - |
| `breakup_shape_k_med` | 1.0 | Weibull scale factor chosen to be 1 so that duration of medium activity level partnerships is exponential | partnerships | - |
| `breakup_shape_k_high` | 1.0 | Weibull scale factor chosen to be 1 so that duration of high activity level partnerships is exponential | partnerships | - |
| `breakup_scale_multiplier_overall` | 1.5 | Multiplier scaling duration of all partnerships | partnerships | - |
| `breakup_scale_multiplier_between_vs_within_patch` | 0.4 | Multiplier scaling duration of partnerships between patches compared to within patches | partnerships | - |
| `PC_Retention_Round1` | 1.0 | Legacy parameter for population cohort retention at R1 (100%) | PC | - |
| `PC_Retention_Round2` | 0.9 | Legacy parameter for population cohort retention at R2 compared to R1 | PC | - |
| `PC_Retention_Round3` | 0.9 | Legacy parameter for population cohort retention at R3 compared to R2 | PC | - |
| `PC_Retention_Round4` | 0.9 | Legacy parameter for population cohort retention at R4 compared to R3 | PC | - |
| `start_time_hiv` | 1975.0 | Start of the HIV epidemic in the country | times | - |
| `start_time_simul` | 1900.0 | Start of simulation | times | - |
| `end_time_simul` | 2020.0 | End of simulation | times | - |
| `COUNTRY_HIV_TEST_START` | 2000.0 | Time when background HIV testing begins | times | - |
| `COUNTRY_ART_START` | 2004.0 | Time when ART first available | times | - |
| `COUNTRY_CD4_350_START` | 2011.0 | Time when ART guidelines changed to CD4<350 | times | - |
| `COUNTRY_CD4_500_START` | 2014.5 | Time when ART guidelines changed to CD4<500 (excluding PopART arm A communities) | times | - |
| `COUNTRY_IMMEDIATE_ART_START` | 2016.3333 | Time when ART guidelines changed to immediate treatment (excluding PopART arm A communities) | times | - |
| `COUNTRY_VMMC_START` | 2010.0 | Time when VMMC first available | times | - |
| `CHIPS_YEAR1_START` | 2013.9375 | Start of CHiPs round 1 | times | - |
| `CHIPS_YEAR1_END` | 2015.5 | End of CHiPs round 1 | times | - |
| `CHIPS_YEAR2_START` | 2015.5 | Start of CHiPs round 2 | times | - |
| `CHIPS_YEAR2_END` | 2016.6666670000002 | End of CHiPs round 2 | times | - |
| `CHIPS_YEAR3_START` | 2016.7083329999998 | Start of CHiPs round 3 | times | - |
| `CHIPS_YEAR3_END` | 2017.9791670000002 | End of CHiPs round 3 | times | - |
| `DURATION_POST_TRIAL_ROUND` | 48.0 | Number of timesteps in a post-trial CHiPs round (1 year) | times | - |
| `NDHSROUNDS` | 3.0 | Number of DHS rounds conducted in country | times | - |
| `DHSROUND1` | 2002.0 | Time of first DHS round | times | - |
| `DHSROUND2` | 2007.0 | Time of second DHS round | times | - |
| `DHSROUND3` | 2013.0 | Time of third DHS round | times | - |