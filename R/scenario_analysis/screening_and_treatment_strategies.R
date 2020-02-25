# Strategies for screening and treatment simulation

require(here)  # for setting working directory
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Load the calibrated parmsets
load(here("calibration", "input", "mock_parmsets_210220.Rdata"))  # params_mat_mock
params_mat <- params_mat_mock  # change to what is read in

# ASSUMPTION A SIMULATIONS

# Prepare the parameter sets that don't change
scen_a <- parameter_list
scen_a$screening_coverage <- 0.9
scen_a$min_age_to_screen <- 30
scen_a$max_age_to_screen <- 70
scen_a$prop_to_vaccinate <- 0
scen_a$link_to_care_prob <- 0.8
scen_a$treatment_initiation_prob <- 1
scen_a$monitoring_prob <- 0.8
scen_a$apply_bdvacc <- 0
scen_a$apply_treat_it <- 0

# Simulate

# Sim1: counterfactual - status quo
sim1 <- apply(params_mat,1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = scen_a,
                          parms_to_change = list(b1 = as.list(x)$b1,b2 = as.list(x)$b2,b3 = as.list(x)$b3,
                                                 mtct_prob_s = as.list(x)$mtct_prob_s,mtct_prob_e = as.list(x)$mtct_prob_e,
                                                 alpha = as.list(x)$alpha, p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                                 pr_it_ir = as.list(x)$pr_it_ir,pr_ir_ic = as.list(x)$pr_ir_ic,
                                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,sag_loss_slope = as.list(x)$sag_loss_slope,
                                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                                 hccr_it = as.list(x)$hccr_it, hccr_ir = as.list(x)$hccr_ir,
                                                 hccr_enchb = as.list(x)$hccr_enchb, hccr_cc = as.list(x)$hccr_cc,
                                                 hccr_dcc = as.list(x)$hccr_dcc,mu_cc = as.list(x)$mu_cc,
                                                 mu_dcc = as.list(x)$mu_dcc,mu_hcc = as.list(x)$mu_hcc,
                                                 vacc_eff = as.list(x)$vacc_eff),
                          scenario = "vacc"))
out1 <- lapply(sim1, code_model_output)

# sim2 to do

# Sim3: one-off screen, no monitoring
sim3 <- apply(params_mat,1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = scen_a,
                          parms_to_change = list(b1 = as.list(x)$b1,b2 = as.list(x)$b2,b3 = as.list(x)$b3,
                                                 mtct_prob_s = as.list(x)$mtct_prob_s,mtct_prob_e = as.list(x)$mtct_prob_e,
                                                 alpha = as.list(x)$alpha, p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                                 pr_it_ir = as.list(x)$pr_it_ir,pr_ir_ic = as.list(x)$pr_ir_ic,
                                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,sag_loss_slope = as.list(x)$sag_loss_slope,
                                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                                 hccr_it = as.list(x)$hccr_it, hccr_ir = as.list(x)$hccr_ir,
                                                 hccr_enchb = as.list(x)$hccr_enchb, hccr_cc = as.list(x)$hccr_cc,
                                                 hccr_dcc = as.list(x)$hccr_dcc,mu_cc = as.list(x)$mu_cc,
                                                 mu_dcc = as.list(x)$mu_dcc,mu_hcc = as.list(x)$mu_hcc,
                                                 vacc_eff = as.list(x)$vacc_eff,
                                                 screening_years = c(2020),
                                                 monitoring_rate = 0),
                          scenario = "vacc_screen"))
out3 <- lapply(sim3, code_model_output)

# Sim4: one-off screen, monitor every 10 years
sim4 <- apply(params_mat,1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = scen_a,
                          parms_to_change = list(b1 = as.list(x)$b1,b2 = as.list(x)$b2,b3 = as.list(x)$b3,
                                                 mtct_prob_s = as.list(x)$mtct_prob_s,mtct_prob_e = as.list(x)$mtct_prob_e,
                                                 alpha = as.list(x)$alpha, p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                                 pr_it_ir = as.list(x)$pr_it_ir,pr_ir_ic = as.list(x)$pr_ir_ic,
                                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,sag_loss_slope = as.list(x)$sag_loss_slope,
                                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                                 hccr_it = as.list(x)$hccr_it, hccr_ir = as.list(x)$hccr_ir,
                                                 hccr_enchb = as.list(x)$hccr_enchb, hccr_cc = as.list(x)$hccr_cc,
                                                 hccr_dcc = as.list(x)$hccr_dcc,mu_cc = as.list(x)$mu_cc,
                                                 mu_dcc = as.list(x)$mu_dcc,mu_hcc = as.list(x)$mu_hcc,
                                                 vacc_eff = as.list(x)$vacc_eff,
                                                 screening_years = c(2020),
                                                 monitoring_rate = 1/10),
                          scenario = "vacc_screen"))
out4 <- lapply(sim4, code_model_output)

# Sim5: one-off screen, monitor every 5 years
sim5 <- apply(params_mat,1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = scen_a,
                          parms_to_change = list(b1 = as.list(x)$b1,b2 = as.list(x)$b2,b3 = as.list(x)$b3,
                                                 mtct_prob_s = as.list(x)$mtct_prob_s,mtct_prob_e = as.list(x)$mtct_prob_e,
                                                 alpha = as.list(x)$alpha, p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                                 pr_it_ir = as.list(x)$pr_it_ir,pr_ir_ic = as.list(x)$pr_ir_ic,
                                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,sag_loss_slope = as.list(x)$sag_loss_slope,
                                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                                 hccr_it = as.list(x)$hccr_it, hccr_ir = as.list(x)$hccr_ir,
                                                 hccr_enchb = as.list(x)$hccr_enchb, hccr_cc = as.list(x)$hccr_cc,
                                                 hccr_dcc = as.list(x)$hccr_dcc,mu_cc = as.list(x)$mu_cc,
                                                 mu_dcc = as.list(x)$mu_dcc,mu_hcc = as.list(x)$mu_hcc,
                                                 vacc_eff = as.list(x)$vacc_eff,
                                                 screening_years = c(2020),
                                                 monitoring_rate = 1/5),
                          scenario = "vacc_screen"))
out5 <- lapply(sim5, code_model_output)

# Sim6: one-off screen, monitor every year
sim6 <- apply(params_mat,1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = scen_a,
                          parms_to_change = list(b1 = as.list(x)$b1,b2 = as.list(x)$b2,b3 = as.list(x)$b3,
                                                 mtct_prob_s = as.list(x)$mtct_prob_s,mtct_prob_e = as.list(x)$mtct_prob_e,
                                                 alpha = as.list(x)$alpha, p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                                 pr_it_ir = as.list(x)$pr_it_ir,pr_ir_ic = as.list(x)$pr_ir_ic,
                                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,sag_loss_slope = as.list(x)$sag_loss_slope,
                                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                                 hccr_it = as.list(x)$hccr_it, hccr_ir = as.list(x)$hccr_ir,
                                                 hccr_enchb = as.list(x)$hccr_enchb, hccr_cc = as.list(x)$hccr_cc,
                                                 hccr_dcc = as.list(x)$hccr_dcc,mu_cc = as.list(x)$mu_cc,
                                                 mu_dcc = as.list(x)$mu_dcc,mu_hcc = as.list(x)$mu_hcc,
                                                 vacc_eff = as.list(x)$vacc_eff,
                                                 screening_years = c(2020),
                                                 monitoring_rate = 1),
                          scenario = "vacc_screen"))
out6 <- lapply(sim6, code_model_output)

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

hbv_deaths_averted <- rbind(calculate_number_averted(counterfactual_metric = cum_hbv_deaths1,
                                                     scenario_metric = cum_hbv_deaths3, summarise = FALSE),
                            calculate_number_averted(counterfactual_metric = cum_hbv_deaths1,
                                                     scenario_metric = cum_hbv_deaths4, summarise = FALSE),
                            calculate_number_averted(counterfactual_metric = cum_hbv_deaths1,
                                                     scenario_metric = cum_hbv_deaths5, summarise = FALSE),
                            calculate_number_averted(counterfactual_metric = cum_hbv_deaths1,
                                                     scenario_metric = cum_hbv_deaths6, summarise = FALSE),
                            calculate_number_averted(counterfactual_metric = cum_hbv_deaths3,
                                                     scenario_metric = cum_hbv_deaths4, summarise = FALSE),
                            calculate_number_averted(counterfactual_metric = cum_hbv_deaths3,
                                                     scenario_metric = cum_hbv_deaths5, summarise = FALSE),
                            calculate_number_averted(counterfactual_metric = cum_hbv_deaths3,
                                                     scenario_metric = cum_hbv_deaths6, summarise = FALSE))



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
