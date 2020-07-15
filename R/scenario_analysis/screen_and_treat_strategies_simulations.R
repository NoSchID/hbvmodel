# Strategies for screening and treatment simulation

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Load the calibrated parmsets
load(here("calibration", "input", "accepted_parmsets_123_180520.Rdata")) # params_mat_targets5

# ASSUMPTION A SIMULATIONS

# Prepare the parameter sets that don't change
scenario_a_parms <- parameter_list
scenario_a_parms$screening_coverage <- 0.9
scenario_a_parms$min_age_to_screen <- 30
scenario_a_parms$max_age_to_screen <- 70
scenario_a_parms$prop_to_vaccinate <- 0
scenario_a_parms$link_to_care_prob <- 0.8
scenario_a_parms$treatment_initiation_prob <- 1
scenario_a_parms$monitoring_prob <- 0.8
scenario_a_parms$apply_bdvacc <- 0
scenario_a_parms$apply_treat_it <- 0
#save(scenario_a_parms, file= here("analysis_input", "scenario_a_parms.Rdata"))

# ASSUMPTION B SIMULATIONS (Optimal with PMTCT)

scenario_b_parms <- parameter_list
scenario_b_parms$screening_coverage <- 0.9
scenario_b_parms$min_age_to_screen <- 30
scenario_b_parms$max_age_to_screen <- 70
scenario_b_parms$prop_to_vaccinate <- 0
scenario_b_parms$link_to_care_prob <- 0.8
scenario_b_parms$treatment_initiation_prob <- 1
scenario_b_parms$monitoring_prob <- 0.8
scenario_b_parms$apply_treat_it <- 0
scenario_b_parms$apply_bdvacc <- 1
scenario_b_parms$bdvacc_introtime <- 2020
scenario_b_parms$apply_bdvacc_linear_scale_up <- 1
#save(scenario_b_parms, file= here("analysis_input", "scenario_b_parms.Rdata"))

# ASSUMPTION D1 SIMULATIONS (Optimal with younger eligibility)
scenario_d1_parms <- parameter_list
scenario_d1_parms$screening_coverage <- 0.9
scenario_d1_parms$min_age_to_screen <- 15
scenario_d1_parms$max_age_to_screen <- 65
scenario_d1_parms$prop_to_vaccinate <- 0
scenario_d1_parms$link_to_care_prob <- 0.8
scenario_d1_parms$treatment_initiation_prob <- 1
scenario_d1_parms$monitoring_prob <- 0.8
scenario_d1_parms$apply_bdvacc <- 0
scenario_d1_parms$apply_treat_it <- 0
#save(scenario_d1_parms, file= here("analysis_input", "scenario_d1_parms.Rdata"))

# ASSUMPTION D2 SIMULATIONS (Optimal with different age groups)
scenario_d2_parms <- parameter_list
scenario_d2_parms$screening_coverage <- 0.9
scenario_d2_parms$min_age_to_screen <- 45
scenario_d2_parms$max_age_to_screen <- 70
scenario_d2_parms$prop_to_vaccinate <- 0
scenario_d2_parms$link_to_care_prob <- 0.8
scenario_d2_parms$treatment_initiation_prob <- 1
scenario_d2_parms$monitoring_prob <- 0.8
scenario_d2_parms$apply_bdvacc <- 0
scenario_d2_parms$apply_treat_it <- 0
#save(scenario_d2_parms, file= here("analysis_input", "scenario_d2_parms.Rdata"))

# ASSUMPTION D3 SIMULATIONS (Optimal with different age groups)
scenario_d3_parms <- parameter_list
scenario_d3_parms$screening_coverage <- 0.9
scenario_d3_parms$min_age_to_screen <- 15
scenario_d3_parms$max_age_to_screen <- 45
scenario_d3_parms$prop_to_vaccinate <- 0
scenario_d3_parms$link_to_care_prob <- 0.8
scenario_d3_parms$treatment_initiation_prob <- 1
scenario_d3_parms$monitoring_prob <- 0.8
scenario_d3_parms$apply_bdvacc <- 0
scenario_d3_parms$apply_treat_it <- 0
#save(scenario_d3_parms, file= here("analysis_input", "scenario_d3_parms.Rdata"))

# Simulate ----

# Run status quo scenario in parallel

# Sim0: counterfactual - status quo, but with cohort
#cohort_scenario_a_parms <- scenario_a_parms
#cohort_scenario_a_parms$apply_screen_not_treat <- 1
#sim0 <- run_one_screening_scenario(default_parameter_list = cohort_scenario_a_parms,
#                                   calibrated_parameter_sets = params_mat_targets5,
#                                   years_of_test = c(2020), monitoring_rate = 0,
#                                   drop_timesteps_before = 1960,
#                                   label = "status_quo_cohort")
#saveRDS(sim0, here("a_sim0_status_quo_040320.rds"))


# Sim1: counterfactual - status quo

#tic()
#library(parallel)
#cl <- makeCluster(2)
#clusterEvalQ(cl, {library(dplyr); library(tidyr); library(deSolve); library(binom)})
#clusterExport(cl, ls())
#sim1 <- run_one_scenario_parallel(default_parameter_list = scenario_a_parms,
#                                   calibrated_parameter_sets = params_mat_targets5,
#                                   drop_timesteps_before = 1960,
#                                   scenario = "vacc")
#stopCluster(cl)
#toc()

#saveRDS(sim1, here("a_sim1_status_quo_290220.rds"))

# sim2 to do

# Sim3: one-off screen, no monitoring
#sim3 <- run_one_screening_scenario(default_parameter_list = scenario_a_parms,
#                           calibrated_parameter_sets = params_mat,
#                           years_of_test = c(2020), monitoring_rate = 0,
#                            drop_timesteps_before = 1960,
#                            label = "screen_2020_monit_0")

# Sim4: one-off screen, monitor every 10 years
#sim4 <- run_one_screening_scenario(default_parameter_list = scenario_a_parms,
#                                   calibrated_parameter_sets = params_mat_targets5,
#                                   years_of_test = c(2020), monitoring_rate = 1/10,
#                                   drop_timesteps_before = 1960,
#                                   label = "screen_2020_monit_10")
#saveRDS(sim4, here("a_sim4_screen_2020_monit_10_030320.rds"))


# Sim5: one-off screen, monitor every 5 years
#sim5 <- run_one_screening_scenario(default_parameter_list = scenario_a_parms,
#                                   calibrated_parameter_sets = params_mat_targets5,
#                                   years_of_test = c(2020), monitoring_rate = 1/5,
#                                   drop_timesteps_before = 1960,
#                                   label = "screen_2020_monit_5")
#saveRDS(sim5, here("a_sim5_screen_2020_monit_5_020320.rds"))

# Sim6: one-off screen, monitor every year
#sim6 <- run_one_screening_scenario(default_parameter_list = scenario_a_parms,
#                                   calibrated_parameter_sets = params_mat_targets5,
#                                   years_of_test = c(2020), monitoring_rate = 1,
#                                   drop_timesteps_before = 1960,
#                                   label = "screen_2020_monit_1")

#saveRDS(sim6, here("a_sim6_screen_2020_monit_1_020320.rds"))


# Calculate some outcome examples
cum_hbv_deaths1 <- extract_cumulative_hbv_deaths(out1, scenario_label = "status_quo",
                                                 from_year = 2020, by_year = 2050)
cum_hbv_deaths3 <- extract_cumulative_hbv_deaths(out3, scenario_label = "screen1_monitor0",
                                                 from_year = 2020, by_year = 2050)
cum_hbv_deaths4 <- extract_cumulative_hbv_deaths(out4, scenario_label = "screen1_monitor10",
                                                 from_year = 2020, by_year = 2050)
cum_hbv_deaths5 <- extract_cumulative_hbv_deaths(out5, scenario_label = "screen1_monitor5",
                                                 from_year = 2020, by_year = 2050)
cum_hbv_deaths6 <- extract_cumulative_hbv_deaths(out6, scenario_label = "screen1_monitor1",
                                                 from_year = 2020, by_year = 2050)

ly3_c <- extract_cohort_life_years_lived(out3, "3")
ly4_c <- extract_cohort_life_years_lived(out4, "4")

calculate_cohort_number_averted(counterfactual_metric = ly3_c,
                         scenario_metric = ly4_c, summarise = FALSE)

lyl3 <- extract_life_years_lived(out3, scenario_label = "3",
                                 from_year = 2020, by_year = 2120)
lyl4 <- extract_life_years_lived(out4, scenario_label = "4",
                                 from_year = 2020, by_year = 2120)
calculate_number_averted(counterfactual_metric = lyl3,
                         scenario_metric = lyl4, summarise = FALSE)

cum_chronic_inf1 <- extract_cumulative_chronic_infections(out1, scenario_label = "status_quo",
                                                           from_year = 2020, by_year = 2050)
cum_chronic_inf3 <- extract_cumulative_chronic_infections(out3, scenario_label = "screen1_monitor0",
                                                          from_year = 2020, by_year = 2050)
calculate_number_averted(counterfactual_metric = cum_chronic_inf1,
                         scenario_metric = cum_chronic_inf3, summarise = FALSE)

lyl1 <- extract_life_years_lived(out1, scenario_label = "status_quo",
                                 from_year = 2020, by_year = 2050)
lyl3 <- extract_life_years_lived(out3, scenario_label = "screen1_monitor0",
                                 from_year = 2020, by_year = 2050)
calculate_number_averted(counterfactual_metric = lyl1,
                         scenario_metric = lyl3, summarise = FALSE)

summarise_average_age_at_death(out3,scenario_label = "screen1_monitor0")

# Extract outcomes (generic)----

# Load outcomes 0-6

#out <- readRDS(paste0(out_path, "a_extracted_output_0to6_050320.rds"))

#out0_cohort <- out$out0_cohort
#out1 <- out$out1
#out3 <- out$out3
#out4 <- out$out4
#out5 <- out$out5
#out6 <- out$out6

label <- "status_quo"

out_path <- "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/screen_and_treat_strategies/Scenario A with new parmsets/"

#out <- readRDS(paste0(out_path, "a_sim0_status_quo_040320.rds"))
#out <- readRDS(paste0(out_path, "a_sim1_status_quo_290220.rds"))
#out <- readRDS(paste0(out_path, "a_sim3_screen_2020_monit_0_280220.rds"))
#out <- readRDS(paste0(out_path, "a_sim4_screen_2020_monit_10_030320.rds"))
#out <- readRDS(paste0(out_path, "a_sim5_screen_2020_monit_5_020320.rds"))
#out <- readRDS(paste0(out_path, "a_sim6_screen_2020_monit_1_020320.rds"))
#out_lsoda <- readRDS(paste0(out_path, "sim6_lsoda.rds"))

out <- readRDS(paste0(out_path, "a_sim2_status_quo_220620.rds"))
out <- out[[1]]

# Cohort outcomes
cohort_age_at_death <- summarise_cohort_average_age_at_death(out,scenario_label = label)
cohort_cum_hbv_deaths <- extract_cohort_cumulative_hbv_deaths(out, label)
cohort_ly <- extract_cohort_life_years_lived(out,label)
cohort_size <- extract_cohort_size(out, label)

# Population outcomes
cum_hbv_deaths_2030 <- extract_cumulative_hbv_deaths(out, scenario_label = label,
                                                          from_year = 2020, by_year = 2030)
cum_hbv_deaths_2050 <- extract_cumulative_hbv_deaths(out, scenario_label = label,
                                                          from_year = 2020, by_year = 2050)
cum_hbv_deaths_2100 <- extract_cumulative_hbv_deaths(out, scenario_label = label,
                                                          from_year = 2020, by_year = 2100)
ly_2030 <- extract_life_years_lived(out, scenario_label = label,
                                         from_year = 2020, by_year = 2030)
ly_2050 <- extract_life_years_lived(out, scenario_label = label,
                                         from_year = 2020, by_year = 2050)
ly_2100 <- extract_life_years_lived(out, scenario_label = label,
                                         from_year = 2020, by_year = 2100)
interactions_2030 <- summarise_healthcare_interactions(out, from_year = 2020,
                                                            by_year = 2030, scenario_label = label)
interactions_2050 <- summarise_healthcare_interactions(out, from_year = 2020,
                                                            by_year = 2050, scenario_label = label)
interactions_2100 <- summarise_healthcare_interactions(out, from_year = 2020,
                                                            by_year = 2100, scenario_label = label)

# Change object names
#out4 <- list(cohort_age_at_death = cohort_age_at_death,
#             cohort_cum_hbv_deaths = cohort_cum_hbv_deaths,
#             cohort_ly = cohort_ly,
#             cohort_size = cohort_size,
#             cum_hbv_deaths_2030 = cum_hbv_deaths_2030,
#             cum_hbv_deaths_2050 = cum_hbv_deaths_2050,
#             cum_hbv_deaths_2100 = cum_hbv_deaths_2100,
#             ly_2030 = ly_2030,
#             ly_2050 = ly_2050,
#             ly_2100 = ly_2100,
#             interactions_2030 = interactions_2030,  # NA for no treatment
#             interactions_2050 = interactions_2050,  # NA for no treatment
#             interactions_2100 = interactions_2100)  # NA for no treatment

out2_ts <- summarise_time_series(out, scenario_label = label, summarise_percentiles = FALSE)

# Outcomes for the status quo scenario
out2 <- list(cum_hbv_deaths_2030 = cum_hbv_deaths_2030,
             cum_hbv_deaths_2050 = cum_hbv_deaths_2050,
             cum_hbv_deaths_2100 = cum_hbv_deaths_2100,
             ly_2030 = ly_2030,
             ly_2050 = ly_2050,
             ly_2100 = ly_2100,
             timeseries = out2_ts)
out2 <- list(out2)
names(out2) <- "status_quo"
saveRDS(out2, paste0(out_path, "a_out2_status_quo_250620.rds"))

rm(out)
gc()

out0to6 <- list(out0_cohort = out0_cohort,
                out1 = out1,
                out3 = out3,
                out4 = out4,
                out5 = out5,
                out6 = out6,
                out0_cohort_timeseries = out0_cohort_ts,
                out1_timeseries = out1_ts,
                out3_timeseries = out3_ts,
                out4_timeseries = out4_ts,
                out5_timeseries = out5_ts,
                out6_timeseries = out6_ts)

#saveRDS(out0to6, here("a_extracted_output_0to6_050320.rds"))

# Compare with lsoda
range(out6$cohort_age_at_death[-1]-out6_lsoda$cohort_age_at_death[-1])
range(out6$cohort_cum_hbv_deaths[-1]-out6_lsoda$cohort_cum_hbv_deaths[-1])
range(out6$cohort_ly[-1]-out6_lsoda$cohort_ly[-1])
range(out6$cohort_size[-1]/out6_lsoda$cohort_size[-1])
range(out6$cum_hbv_deaths_2030[-c(1:3)]-out6_lsoda$cum_hbv_deaths_2030[-c(1:3)])
range(out6$cum_hbv_deaths_2050[-c(1:3)]-out6_lsoda$cum_hbv_deaths_2050[-c(1:3)])
range(out6$cum_hbv_deaths_2100[-c(1:3)]-out6_lsoda$cum_hbv_deaths_2100[-c(1:3)])
range(out6$ly_2030[-c(1:3)]-out6_lsoda$ly_2030[-c(1:3)])
range(out6$ly_2050[-c(1:3)]-out6_lsoda$ly_2050[-c(1:3)])
range(out6$ly_2100[-c(1:3)]-out6_lsoda$ly_2100[-c(1:3)])
range(out6$interactions_2030$total_interactions[-c(1:3)]-out6_lsoda$interactions_2030$total_interactions[-c(1:3)])
range(out6$interactions_2050$total_interactions[-c(1:3)]-out6_lsoda$interactions_2050$total_interactions[-c(1:3)])
range(out6$interactions_2100$total_interactions[-c(1:3)]-out6_lsoda$interactions_2100$total_interactions[-c(1:3)])
