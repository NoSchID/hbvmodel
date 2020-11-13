# Explore the effect of treatment depending on which age group is targeted in screening in 2020

require(here)  # for setting working directory
source(here("R/imperial_model_interventions.R"))
load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans

# Simulations
# sim2 = no treatment, population
# simy = screen and treat 15-30 year olds in 2020
# simo = screen and treat 45-65 year olds in 2020
# sim1y = no treatment, isolate cohort of 15-30 year olds
# sim10 = no treatment, isolate cohort of 45-65 year olds
# all without monitoring

# On simulation 111
# sim2 = no treatment, population ----
sim2 <- apply(params_mat_accepted_kmeans[111,],1,
             function(x)
               run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                         parms_to_change =
                           list(b1 = as.list(x)$b1,
                                b2 = as.list(x)$b2,
                                b3 = as.list(x)$b3,
                                mtct_prob_s = as.list(x)$mtct_prob_s,
                                mtct_prob_e = as.list(x)$mtct_prob_e,
                                alpha = as.list(x)$alpha,
                                p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                pr_it_ir = as.list(x)$pr_it_ir,
                                pr_ir_ic = as.list(x)$pr_ir_ic,
                                eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                sag_loss_slope = as.list(x)$sag_loss_slope,
                                pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                hccr_it = as.list(x)$hccr_it,
                                hccr_ir = as.list(x)$hccr_ir,
                                hccr_enchb = as.list(x)$hccr_enchb,
                                hccr_cc = as.list(x)$hccr_cc,
                                hccr_dcc = as.list(x)$hccr_dcc,
                                mu_cc = as.list(x)$mu_cc,
                                mu_dcc = as.list(x)$mu_dcc,
                                mu_hcc = as.list(x)$mu_hcc,
                                vacc_eff = as.list(x)$vacc_eff,
                                screening_years = c(2020),
                                screening_coverage = 0.9,
                                apply_treat_it = 0,
                                prop_negative_to_remove_from_rescreening = 0,
                                apply_screen_not_treat = 0,
                                monitoring_rate = 0,
                                apply_repeat_screen = 0,
                                min_age_to_screen = 15,
                                max_age_to_screen = 30-da,
                                min_age_to_repeat_screen = 15,
                                max_age_to_repeat_screen = 60,
                                repeat_screening_years = seq(2030,2100, by = 10)),
                         drop_timesteps_before = 1960,
                         scenario = "vacc"))

out2 <- code_model_output(sim2[[1]])

# simy = screen and treat 15-30 year olds in 2020 ----
simy <- apply(params_mat_accepted_kmeans[111,],1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                          parms_to_change =
                            list(b1 = as.list(x)$b1,
                                 b2 = as.list(x)$b2,
                                 b3 = as.list(x)$b3,
                                 mtct_prob_s = as.list(x)$mtct_prob_s,
                                 mtct_prob_e = as.list(x)$mtct_prob_e,
                                 alpha = as.list(x)$alpha,
                                 p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                 pr_it_ir = as.list(x)$pr_it_ir,
                                 pr_ir_ic = as.list(x)$pr_ir_ic,
                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                 pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                 sag_loss_slope = as.list(x)$sag_loss_slope,
                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                 hccr_it = as.list(x)$hccr_it,
                                 hccr_ir = as.list(x)$hccr_ir,
                                 hccr_enchb = as.list(x)$hccr_enchb,
                                 hccr_cc = as.list(x)$hccr_cc,
                                 hccr_dcc = as.list(x)$hccr_dcc,
                                 mu_cc = as.list(x)$mu_cc,
                                 mu_dcc = as.list(x)$mu_dcc,
                                 mu_hcc = as.list(x)$mu_hcc,
                                 vacc_eff = as.list(x)$vacc_eff,
                                 screening_years = c(2020),
                                 screening_coverage = 0.9,
                                 apply_treat_it = 0,
                                 prop_negative_to_remove_from_rescreening = 0,
                                 apply_screen_not_treat = 0,
                                 monitoring_rate = 0,
                                 apply_repeat_screen = 0,
                                 min_age_to_screen = 15,
                                 max_age_to_screen = 30-da,
                                 min_age_to_repeat_screen = 15,
                                 max_age_to_repeat_screen = 60,
                                 repeat_screening_years = seq(2030,2100, by = 10)),
                          drop_timesteps_before = 1960,
                          scenario = "vacc_screen"))

outy <- code_model_output(simy[[1]])

# simo = screen and treat 45-65 year olds in 2020 ----
simo <- apply(params_mat_accepted_kmeans[111,],1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                          parms_to_change =
                            list(b1 = as.list(x)$b1,
                                 b2 = as.list(x)$b2,
                                 b3 = as.list(x)$b3,
                                 mtct_prob_s = as.list(x)$mtct_prob_s,
                                 mtct_prob_e = as.list(x)$mtct_prob_e,
                                 alpha = as.list(x)$alpha,
                                 p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                 pr_it_ir = as.list(x)$pr_it_ir,
                                 pr_ir_ic = as.list(x)$pr_ir_ic,
                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                 pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                 sag_loss_slope = as.list(x)$sag_loss_slope,
                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                 hccr_it = as.list(x)$hccr_it,
                                 hccr_ir = as.list(x)$hccr_ir,
                                 hccr_enchb = as.list(x)$hccr_enchb,
                                 hccr_cc = as.list(x)$hccr_cc,
                                 hccr_dcc = as.list(x)$hccr_dcc,
                                 mu_cc = as.list(x)$mu_cc,
                                 mu_dcc = as.list(x)$mu_dcc,
                                 mu_hcc = as.list(x)$mu_hcc,
                                 vacc_eff = as.list(x)$vacc_eff,
                                 screening_years = c(2020),
                                 screening_coverage = 0.9,
                                 apply_treat_it = 0,
                                 prop_negative_to_remove_from_rescreening = 0,
                                 apply_screen_not_treat = 0,
                                 monitoring_rate = 0,
                                 apply_repeat_screen = 0,
                                 min_age_to_screen = 45,
                                 max_age_to_screen = 65-da,
                                 min_age_to_repeat_screen = 15,
                                 max_age_to_repeat_screen = 60,
                                 repeat_screening_years = seq(2030,2100, by = 10)),
                          drop_timesteps_before = 1960,
                          scenario = "vacc_screen"))

outo <- code_model_output(simo[[1]])

# sim1y = no treatment, isolate cohort of 15-30 year olds ----
sim1y <- apply(params_mat_accepted_kmeans[111,],1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                          parms_to_change =
                            list(b1 = as.list(x)$b1,
                                 b2 = as.list(x)$b2,
                                 b3 = as.list(x)$b3,
                                 mtct_prob_s = as.list(x)$mtct_prob_s,
                                 mtct_prob_e = as.list(x)$mtct_prob_e,
                                 alpha = as.list(x)$alpha,
                                 p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                 pr_it_ir = as.list(x)$pr_it_ir,
                                 pr_ir_ic = as.list(x)$pr_ir_ic,
                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                 pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                 sag_loss_slope = as.list(x)$sag_loss_slope,
                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                 hccr_it = as.list(x)$hccr_it,
                                 hccr_ir = as.list(x)$hccr_ir,
                                 hccr_enchb = as.list(x)$hccr_enchb,
                                 hccr_cc = as.list(x)$hccr_cc,
                                 hccr_dcc = as.list(x)$hccr_dcc,
                                 mu_cc = as.list(x)$mu_cc,
                                 mu_dcc = as.list(x)$mu_dcc,
                                 mu_hcc = as.list(x)$mu_hcc,
                                 vacc_eff = as.list(x)$vacc_eff,
                                 screening_years = c(2020),
                                 screening_coverage = 0.9,
                                 apply_treat_it = 0,
                                 prop_negative_to_remove_from_rescreening = 0,
                                 apply_screen_not_treat = 1,
                                 monitoring_rate = 0,
                                 apply_repeat_screen = 0,
                                 min_age_to_screen = 15,
                                 max_age_to_screen = 30-da,
                                 min_age_to_repeat_screen = 15,
                                 max_age_to_repeat_screen = 60,
                                 repeat_screening_years = seq(2030,2100, by = 10)),
                          drop_timesteps_before = 1960,
                          scenario = "vacc_screen"))

out1y <- code_model_output(sim1y[[1]])

# sim1o = no treatment, isolate cohort of 45-65 year olds ----
sim1o <- apply(params_mat_accepted_kmeans[111,],1,
              function(x)
                run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                          parms_to_change =
                            list(b1 = as.list(x)$b1,
                                 b2 = as.list(x)$b2,
                                 b3 = as.list(x)$b3,
                                 mtct_prob_s = as.list(x)$mtct_prob_s,
                                 mtct_prob_e = as.list(x)$mtct_prob_e,
                                 alpha = as.list(x)$alpha,
                                 p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                 p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                 p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                 pr_it_ir = as.list(x)$pr_it_ir,
                                 pr_ir_ic = as.list(x)$pr_ir_ic,
                                 eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                 pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                 pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                 pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                 pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                 sag_loss_slope = as.list(x)$sag_loss_slope,
                                 pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                 cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                 pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                 cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                 cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                 cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                 hccr_it = as.list(x)$hccr_it,
                                 hccr_ir = as.list(x)$hccr_ir,
                                 hccr_enchb = as.list(x)$hccr_enchb,
                                 hccr_cc = as.list(x)$hccr_cc,
                                 hccr_dcc = as.list(x)$hccr_dcc,
                                 mu_cc = as.list(x)$mu_cc,
                                 mu_dcc = as.list(x)$mu_dcc,
                                 mu_hcc = as.list(x)$mu_hcc,
                                 vacc_eff = as.list(x)$vacc_eff,
                                 screening_years = c(2020),
                                 screening_coverage = 0.9,
                                 apply_treat_it = 0,
                                 prop_negative_to_remove_from_rescreening = 0,
                                 apply_screen_not_treat = 1,
                                 monitoring_rate = 0,
                                 apply_repeat_screen = 0,
                                 min_age_to_screen = 45,
                                 max_age_to_screen = 65-da,
                                 min_age_to_repeat_screen = 15,
                                 max_age_to_repeat_screen = 60,
                                 repeat_screening_years = seq(2030,2100, by = 10)),
                          drop_timesteps_before = 1960,
                          scenario = "vacc_screen"))

out1o <- code_model_output(sim1o[[1]])




# Extract HBV-related deaths by age
# Extract HCC deaths by age
# Extract life-years by age
# Calculate HBV deaths, HCC deaths averted and LY saved in the population
# Calculate HBV deaths, HCC deaths averted and LY saved in the screened+treated cohort


extract_outcomes_by_age <- function(sim, out) {

  # HBV deaths by age 2020-2100
  cum_hbv_deaths <- sim[,grepl("^cum_hbv_deathsf.",names(sim))]+
                sim[,grepl("^cum_hbv_deathsm.",names(sim))]+
                sim[,grepl("^cum_screened_hbv_deathsf.",names(sim))]+
                sim[,grepl("^cum_screened_hbv_deathsm.",names(sim))]+
                sim[,grepl("^cum_treated_hbv_deathsf.",names(sim))]+
                sim[,grepl("^cum_treated_hbv_deathsm.",names(sim))]

  # Incident deaths in each age group at each timestep
  # Incidence BY timestep shown in row
  inc_hbv_deaths <- cbind(time = sim$time, calculate_incident_numbers(cum_hbv_deaths))
  inc_hbv_deaths <- inc_hbv_deaths[which(sim$time > 2020 & sim$time <=2100),]

  # HCC deaths by age 2020-2100
  cum_hcc_deaths <- sim[,grepl("^cum_hcc_deathsf.",names(sim))]+
    sim[,grepl("^cum_hcc_deathsm.",names(sim))]+
    sim[,grepl("^cum_screened_hcc_deathsf.",names(sim))]+
    sim[,grepl("^cum_screened_hcc_deathsm.",names(sim))]+
    sim[,grepl("^cum_treated_hcc_deathsf.",names(sim))]+
    sim[,grepl("^cum_treated_hcc_deathsm.",names(sim))]

  # Incident deaths in each age group at each timestep
  # Incidence BY timestep shown in row
  inc_hcc_deaths <- cbind(time = sim$time, calculate_incident_numbers(cum_hcc_deaths))
  inc_hcc_deaths <- inc_hcc_deaths[which(sim$time > 2020 & sim$time <=2100),]

  # Cohort HBV deaths by age 2020-2100
  cum_hbv_deaths_cohort <- sim[,grepl("^cum_screened_hbv_deathsf.",names(sim))]+
    sim[,grepl("^cum_screened_hbv_deathsm.",names(sim))]+
    sim[,grepl("^cum_treated_hbv_deathsf.",names(sim))]+
    sim[,grepl("^cum_treated_hbv_deathsm.",names(sim))]

  # Incident deaths in each age group at each timestep
  # Incidence BY timestep shown in row
  inc_hbv_deaths_cohort <- cbind(time = sim$time, calculate_incident_numbers(cum_hbv_deaths_cohort))
  inc_hbv_deaths_cohort <- inc_hbv_deaths_cohort[which(sim$time > 2020 & sim$time <=2100),]

  # Cohort HCC deaths by age 2020-2100
  cum_hcc_deaths_cohort <- sim[,grepl("^cum_screened_hcc_deathsf.",names(sim))]+
    sim[,grepl("^cum_screened_hcc_deathsm.",names(sim))]+
    sim[,grepl("^cum_treated_hcc_deathsf.",names(sim))]+
    sim[,grepl("^cum_treated_hcc_deathsm.",names(sim))]

  # Incident deaths in each age group at each timestep
  # Incidence BY timestep shown in row
  inc_hcc_deaths_cohort <- cbind(time = sim$time, calculate_incident_numbers(cum_hcc_deaths_cohort))
  inc_hcc_deaths_cohort <- inc_hcc_deaths_cohort[which(sim$time > 2020 & sim$time <=2100),]

  # Life-years spent in each age group at each timestep 2020-2100
  ly <- cbind(time = out$time,
                      out$pop_female+out$pop_male)
  ly <- ly[ly$time >=2020 & ly$time <2100,]

  ly_cohort <- cbind(time = out$time,
                      out$screened_pop_female+out$screened_pop_male+
                        out$treated_pop_female+out$treated_pop_male)
  ly_cohort <- ly_cohort[ly_cohort$time >=2020 & ly_cohort$time <2100,]

  return(list(inc_hbv_deaths=inc_hbv_deaths,
              inc_hcc_deaths=inc_hcc_deaths,
              inc_hbv_deaths_cohort=inc_hbv_deaths_cohort,
              inc_hcc_deaths_cohort=inc_hcc_deaths_cohort,
              ly = ly,
              ly_cohort = ly_cohort,
              inc_hcc_deaths_treated = inc_hcc_deaths_treated))

}

sim2_out <- extract_outcomes_by_age(sim2[[1]]$out, out2)
simy_out <- extract_outcomes_by_age(simy[[1]]$out, outy)
simo_out <- extract_outcomes_by_age(simo[[1]]$out, outo)
sim1y_out <- extract_outcomes_by_age(sim1y[[1]]$out, out1y)
sim1o_out <- extract_outcomes_by_age(sim1o[[1]]$out, out1o)

# Averted outcomes by treatment by targeting young ages
hbv_deaths_avertedy <- sim2_out$inc_hbv_deaths[,-1]-simy_out$inc_hbv_deaths[,-1]
hcc_deaths_avertedy <- sim2_out$inc_hcc_deaths[,-1]-simy_out$inc_hcc_deaths[,-1]
ly_savedy <- simy_out$ly[,-1]-sim2_out$ly[,-1]
hbv_deaths_averted_cohorty <- sim1y_out$inc_hbv_deaths_cohort[,-1]-simy_out$inc_hbv_deaths_cohort[,-1]
hcc_deaths_averted_cohorty <- sim1y_out$inc_hcc_deaths_cohort[,-1]-simy_out$inc_hcc_deaths_cohort[,-1]
ly_saved_cohorty <- simy_out$ly_cohort[,-1]-sim1y_out$ly_cohort[,-1]
# Taking colSums calculates cumulative outcome in each age group between 2020-2100

# Averted outcomes by treatment by targeting old ages
hbv_deaths_avertedo <- sim2_out$inc_hbv_deaths[,-1]-simo_out$inc_hbv_deaths[,-1]
hcc_deaths_avertedo <- sim2_out$inc_hcc_deaths[,-1]-simo_out$inc_hcc_deaths[,-1]
ly_savedo <- simo_out$ly[,-1]-sim2_out$ly[,-1]
hbv_deaths_averted_cohorto <- sim1o_out$inc_hbv_deaths_cohort[,-1]-simo_out$inc_hbv_deaths_cohort[,-1]
hcc_deaths_averted_cohorto <- sim1o_out$inc_hcc_deaths_cohort[,-1]-simo_out$inc_hcc_deaths_cohort[,-1]
ly_saved_cohorto <- simo_out$ly_cohort[,-1]-sim1o_out$ly_cohort[,-1]

# Plots

# Only look at LY in the cohort!

# Black = averted by targeting young, blue = averted by targeting old

# Incident deaths averted over time:

# All deaths are averted in the cohort, so hbv_deaths_avertedy=hbv_deaths_averted_cohorty:
plot(x=seq(2020.5,2100,by=0.5), y = rowSums(hbv_deaths_avertedo), type = "l", ylim = c(-2.4,41), col = "blue")
abline(h=0, col = "grey")
lines(x=seq(2020.5,2100,by=0.5), y = rowSums(hbv_deaths_averted_cohorto), col = "blue")
lines(x=seq(2020.5,2100,by=0.5), y = rowSums(hcc_deaths_avertedo), lty = "dashed", col = "blue")
lines(x=seq(2020.5,2100,by=0.5), y = rowSums(hbv_deaths_avertedy))
lines(x=seq(2020.5,2100,by=0.5), y = rowSums(hcc_deaths_avertedy), lty = "dashed")
# HCC deaths are a subset of HBV deaths
# Screening and treating 15-30 year olds averts HCC deaths at first but experiences
# MORE deaths than without treatment later (postpones them)
plot(x=seq(2020,2100-da,by=0.5), y = rowSums(ly_saved_cohorty), type = "l")
lines(x=seq(2020,2100-da,by=0.5), y = rowSums(ly_saved_cohorto), col = "blue")

# Cumulative (over time) deaths averted by age:
# By targeting 15-30:
plot(x=ages, y = colSums(hbv_deaths_avertedy), type = "l", ylim = c(-2.4,25))
abline(h=0, col = "grey")
lines(x=ages, y = colSums(hcc_deaths_avertedy),  lty = "dashed")
abline(v=15, col = "grey")
abline(v=30, col = "grey")
# By targeting 45-60:
plot(x=ages, y = colSums(hbv_deaths_avertedo), type = "l", ylim = c(-2.4,25), col = "blue")
abline(h=0, col = "grey")
lines(x=ages, y = colSums(hcc_deaths_avertedo),  lty = "dashed", col = "blue")
abline(v=45, col = "grey")
abline(v=65, col = "grey")

# Is average age of HCC death changed?
# Plot cumulative HCC and HBV deaths in each age group - can't really see much
plot(x=ages, y = colSums(sim2_out$inc_hcc_deaths[,-1]), type = "l")
lines(x=ages, y = colSums(simy_out$inc_hcc_deaths[,-1]), col = "red")
lines(x=ages, y = colSums(simo_out$inc_hcc_deaths[,-1]), col = "blue")
plot(x=ages, y = colSums(sim2_out$inc_hbv_deaths[,-1]), type = "l")
lines(x=ages, y = colSums(simy_out$inc_hbv_deaths[,-1]), col = "red")
lines(x=ages, y = colSums(simo_out$inc_hbv_deaths[,-1]), col = "blue")

# Shift in pattern of HCC deaths by age in the screened and treated cohort:
plot(x=ages, y = colSums(sim1y_out$inc_hcc_deaths_cohort[,-1]), type = "l", ylim = c(0,20))
lines(x=ages, y = colSums(simy_out$inc_hcc_deaths_cohort[,-1]), col = "red")
plot(x=ages, y = colSums(sim1o_out$inc_hcc_deaths_cohort[,-1]), type = "l")
lines(x=ages, y = colSums(simo_out$inc_hcc_deaths_cohort[,-1]), col = "blue")

###
# Shift in pattern of HCC deaths by age in the treated cohort only:
# HCC deaths are postponed to older age!

# What would happen in those who need treatment but would not get it?
# Following the screened and treated cohort shows this, but gives the total in both
# Need to substract from this the deaths that occur in those screened, to get the deaths in
# treated untreated only.

x1 <- sim1y[[1]]$out[,grepl("^cum_screened_hcc_deathsf.",names(sim1y[[1]]$out))]+
  sim1y[[1]]$out[,grepl("^cum_screened_hcc_deathsm.",names(sim1y[[1]]$out))]+
  sim1y[[1]]$out[,grepl("^cum_treated_hcc_deathsf.",names(sim1y[[1]]$out))]+
  sim1y[[1]]$out[,grepl("^cum_treated_hcc_deathsm.",names(sim1y[[1]]$out))]
# x1 = deaths that would occur in the screened and treated people in total
# if those who need treatment don't get it
x2 <- simy[[1]]$out[,grepl("^cum_screened_hcc_deathsf.",names(simy[[1]]$out))]+
  simy[[1]]$out[,grepl("^cum_screened_hcc_deathsm.",names(simy[[1]]$out))]
# x2 = deaths in those screened but not treated if those who need treatment receive it
x3 <- x1-x2
# x3 = deaths that would occur in those who need treatment but don't receive it
x4 <- simy[[1]]$out[,grepl("^cum_treated_hcc_deathsf.",names(simy[[1]]$out))]+
simy[[1]]$out[,grepl("^cum_treated_hcc_deathsm.",names(simy[[1]]$out))]
# x4 = deaths that occur with treatment in those treated
# Get cumulative values in 2100

plot(x=ages, y = x3[which(simy[[1]]$out$time == 2100),], type = "l") #no treatment
lines(x=ages, y = x4[which(simy[[1]]$out$time == 2100),], col = "red")  #with treatment

# In the older age groups, this is not the case:
y1 <- sim1o[[1]]$out[,grepl("^cum_screened_hcc_deathsf.",names(sim1o[[1]]$out))]+
  sim1o[[1]]$out[,grepl("^cum_screened_hcc_deathsm.",names(sim1o[[1]]$out))]+
  sim1o[[1]]$out[,grepl("^cum_treated_hcc_deathsf.",names(sim1o[[1]]$out))]+
  sim1o[[1]]$out[,grepl("^cum_treated_hcc_deathsm.",names(sim1o[[1]]$out))]
# x1 = deaths that would occur in the screened and treated people in total
# if those who need treatment don't get it
y2 <- simo[[1]]$out[,grepl("^cum_screened_hcc_deathsf.",names(simo[[1]]$out))]+
  simo[[1]]$out[,grepl("^cum_screened_hcc_deathsm.",names(simo[[1]]$out))]
# x2 = deaths in those screened but not treated if those who need treatment receive it
y3 <- y1-y2
# x3 = deaths that would occur in those who need treatment but don't receive it
y4 <- simo[[1]]$out[,grepl("^cum_treated_hcc_deathsf.",names(simo[[1]]$out))]+
  simo[[1]]$out[,grepl("^cum_treated_hcc_deathsm.",names(simo[[1]]$out))]
# x4 = deaths that occur with treatment in those treated
# Get cumulative values in 2100

plot(x=ages, y = y3[which(simy[[1]]$out$time == 2100),], type = "l") #no treatment
lines(x=ages, y = y4[which(simy[[1]]$out$time == 2100),], col = "red")  #with treatment

# Could give the average age at HCC/HBV death for the treated cohort without monitoring
# (with monitoring it is not a single cohort anymore)
# or plot th cumulative deaths by age overall over time, showing a later decrease in
# deaths averted.

###


# Calculate average age at HCC death in the population:
sum((colSums(sim2_out$inc_hcc_deaths[,-1])*ages))/sum(colSums(sim2_out$inc_hcc_deaths[,-1]))
sum((colSums(simo_out$inc_hcc_deaths[,-1])*ages))/sum(colSums(simo_out$inc_hcc_deaths[,-1]))
sum((colSums(simy_out$inc_hcc_deaths[,-1])*ages))/sum(colSums(simy_out$inc_hcc_deaths[,-1]))

# Extension in average age at HCC death by treatment in the young cohort:
sum((colSums(simy_out$inc_hcc_deaths[,-1])*ages))/sum(colSums(simy_out$inc_hcc_deaths[,-1]))-
  sum((colSums(sim1y_out$inc_hcc_deaths[,-1])*ages))/sum(colSums(sim1y_out$inc_hcc_deaths[,-1]))
# 0.4 years
# In HBV deaths overall:
sum((colSums(simy_out$inc_hbv_deaths[,-1])*ages))/sum(colSums(simy_out$inc_hbv_deaths[,-1]))-
  sum((colSums(sim1y_out$inc_hbv_deaths[,-1])*ages))/sum(colSums(sim1y_out$inc_hbv_deaths[,-1]))
# 0.8 years

# Extension in average age at HCC death by treatment in the old cohort:
sum((colSums(simo_out$inc_hcc_deaths[,-1])*ages))/sum(colSums(simo_out$inc_hcc_deaths[,-1]))-
sum((colSums(sim1o_out$inc_hcc_deaths[,-1])*ages))/sum(colSums(sim1o_out$inc_hcc_deaths[,-1]))
# -0.2 years - reduction in average age at HCC deaths
# HBV deaths overall
sum((colSums(simo_out$inc_hbv_deaths[,-1])*ages))/sum(colSums(simo_out$inc_hbv_deaths[,-1]))-
  sum((colSums(sim1o_out$inc_hbv_deaths[,-1])*ages))/sum(colSums(sim1o_out$inc_hbv_deaths[,-1]))
# -0.5 years - reduction in average age at HBV death


# Exploration of LY saved difference between pop and cohort: ----
# Screening and treating 15-30 year olds averts HBV/HCC deaths at first but experiences
# MORE deaths than without treatment later (postpones them)
plot(x=seq(2020,2100-da,by=0.5), y = rowSums(ly_savedy), type = "l")
lines(x=seq(2020,2100-da,by=0.5), y = rowSums(ly_saved_cohorty), col = "red")
# LY saved in the cohort and in the population are NOT the same - why?? Pattern in cohort also makes
# more sense. Though LY saved in the population are higher, overall LY look almost identical:
plot(x=seq(2020,2100-da,by=0.5), y = rowSums(sim2_out$ly), type = "l")
lines(x=seq(2020,2100-da,by=0.5), y = rowSums(simy_out$ly), col = "red")

# In the older age groups, the LY saved is the same in the cohort as in the population
plot(x=seq(2020,2100-da,by=0.5), y = rowSums(ly_savedo), type = "l")
lines(x=seq(2020,2100-da,by=0.5), y = rowSums(ly_saved_cohorto), col = "red")

plot(x=ages, y = colSums(ly_savedy), type = "l")
lines(x=ages, y = colSums(ly_saved_cohorty), col = "red")
abline(v=15)
# This plot suggests LY are saved in age groups that are not affected by the treatment for the 15-30 years
plot(x=ages, y = colSums(ly_savedo), type = "l")
lines(x=ages, y = colSums(ly_saved_cohorto), col = "red")

# MORE infections with treatment:

out2$incident_chronic_infections$horizontal_chronic_infections[
  out2$incident_chronic_infections$time>=2020 & out2$incident_chronic_infections$time<=2100]-
  outy$incident_chronic_infections$horizontal_chronic_infections[
    outy$incident_chronic_infections$time>=2020 & outy$incident_chronic_infections$time<=2100]

out2$incident_chronic_infections$chronic_births[
  out2$incident_chronic_infections$time>=2020 & out2$incident_chronic_infections$time<=2100]-
  outy$incident_chronic_infections$chronic_births[
    outy$incident_chronic_infections$time>=2020 & outy$incident_chronic_infections$time<=2100]

plot(x=ages, y = ly_savedy[which(seq(2020,2100,by=0.5)==2021),], type = "l")
lines(x=ages, y = ly_saved_cohorty[which(seq(2020,2100,by=0.5)==2021),], col = "red")
plot(x=ages, y = ly_savedy[which(seq(2020,2100,by=0.5)==2025),], type = "l")
lines(x=ages, y = ly_saved_cohorty[which(seq(2020,2100,by=0.5)==2025),], col = "red")
plot(x=ages, y = ly_savedy[which(seq(2020,2100,by=0.5)==2050),], type = "l")
lines(x=ages, y = ly_saved_cohorty[which(seq(2020,2100,by=0.5)==2050),], col = "red")


View(cbind(sim2$`615035`$out$time,
           sim2$`615035`$out$cum_births,
           simy$`615035`$out$cum_births,
           sim2$`615035`$out$cum_births-simy$`615035`$out$cum_births))
# Screening and treating leads to a loss of total births
sim1y$`615035`$out$cum_births-simy$`615035`$out$cum_births
sim2$`615035`$out$cum_births-simo$`615035`$out$cum_births

## NO! Births are the same, just the cumulative number resets at 2020.5 once the treatment is in!!


diff(sim2$`615035`$out$cum_births[
  which(sim2$`615035`$out$time==2020.5):which(sim2$`615035`$out$time==2100)])-
  diff(simy$`615035`$out$cum_births[
    which(simy$`615035`$out$time==2020.5):which(simy$`615035`$out$time==2100)])
# On treatment, get more births overall, presumably due to improved survival. With screen but not treat
# these numbers are minimal so it is not due to rounding error.

diff(sim2$`615035`$out$cum_chronic_births[
  which(sim2$`615035`$out$time==2020.5):which(sim2$`615035`$out$time==2100)])-
  diff(simy$`615035`$out$cum_chronic_births[
    which(simy$`615035`$out$time==2020.5):which(simy$`615035`$out$time==2100)])
# Except for the first few years, also get more chronic births overall

# The proportion of chronic births out of total births should NOT be higher with treatment
(diff(sim2$`615035`$out$cum_chronic_births[
  which(sim2$`615035`$out$time==2020.5):which(sim2$`615035`$out$time==2100)])/
    diff(sim2$`615035`$out$cum_births[
      which(sim2$`615035`$out$time==2020.5):which(sim2$`615035`$out$time==2100)]))*100-
  (diff(simy$`615035`$out$cum_chronic_births[
    which(simy$`615035`$out$time==2020.5):which(simy$`615035`$out$time==2100)])/
     diff(simy$`615035`$out$cum_births[
       which(simy$`615035`$out$time==2020.5):which(simy$`615035`$out$time==2100)]))*100
# It seems to be approximately the same (tiny percentages)

sum(ly_savedy)/sum(ly_savedo)
sum(ly_saved_cohorty)/sum(ly_saved_cohorto)

