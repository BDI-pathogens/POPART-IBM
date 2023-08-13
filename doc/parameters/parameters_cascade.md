# Table: PopART-IBM parameter dictionary for cascade file 
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