# Strategies for screening and treatment simulation

require(here)  # for setting working directory
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Load the calibrated parmsets
load(here("calibration", "input", "mock_parmsets_210220.Rdata"))  # params_mat_mock
load(here("calibration", "input", "accepted_parmsets_119_060120.Rdata")) # params_mat_targets5
#params_mat <- params_mat_mock  # change to what is read in

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

# Simulate ----

# Run status quo scenario in parallel
# Sim1: counterfactual - status quo

#tic()
#library(parallel)
#cl <- makeCluster(2)
#clusterEvalQ(cl, {library(dplyr); library(tidyr); library(deSolve); library(binom)})
#clusterExport(cl, ls())
sim1 <- run_one_scenario_parallel(default_parameter_list = scenario_a_parms,
                                   calibrated_parameter_sets = params_mat_targets5,
                                   drop_timesteps_before = 1960,
                                   scenario = "vacc")
#stopCluster(cl)
#toc()

#saveRDS(sim1, here("a_sim1_status_quo_290220.rds"))

# sim2 to do

# Sim3: one-off screen, no monitoring
sim3 <- run_one_screening_scenario(default_parameter_list = scenario_a_parms,
                           calibrated_parameter_sets = params_mat,
                           years_of_test = c(2020), monitoring_rate = 0,
                            drop_timesteps_before = 1960,
                            label = "screen_2020_monit_0")

# Sim4: one-off screen, monitor every 10 years
sim4 <- run_one_screening_scenario(default_parameter_list = scenario_a_parms,
                                   calibrated_parameter_sets = params_mat_targets5,
                                   years_of_test = c(2020), monitoring_rate = 1/10,
                                   drop_timesteps_before = 1960,
                                   label = "screen_2020_monit_10")
#saveRDS(sim4, here("a_sim4_screen_2020_monit_10_030320.rds"))


# Sim5: one-off screen, monitor every 5 years
sim5 <- run_one_screening_scenario(default_parameter_list = scenario_a_parms,
                                   calibrated_parameter_sets = params_mat_targets5,
                                   years_of_test = c(2020), monitoring_rate = 1/5,
                                   drop_timesteps_before = 1960,
                                   label = "screen_2020_monit_5")
#saveRDS(sim5, here("a_sim5_screen_2020_monit_5_020320.rds"))

# Sim6: one-off screen, monitor every year
sim6 <- run_one_screening_scenario(default_parameter_list = scenario_a_parms,
                                   calibrated_parameter_sets = params_mat_targets5,
                                   years_of_test = c(2020), monitoring_rate = 1,
                                   drop_timesteps_before = 1960,
                                   label = "screen_2020_monit_1")

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

### Extract outcomes (generic)----
out_path <- "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/screen_and_treat_strategies/Scenario A/"

label <- "screen_2020_monit_1"

#out <- readRDS(paste0(out_path, "a_sim1_status_quo_290220.rds"))
#out <- readRDS(paste0(out_path, "a_sim3_screen_2020_monit_0_280220.rds"))
#out <- readRDS(paste0(out_path, "a_sim4_screen_2020_monit_10_030320.rds"))
#out <- readRDS(paste0(out_path, "a_sim5_screen_2020_monit_5_020320.rds"))
#out <- readRDS(paste0(out_path, "a_sim6_screen_2020_monit_1_020320.rds"))

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
out6 <- list(cohort_age_at_death = cohort_age_at_death,
             cohort_cum_hbv_deaths = cohort_cum_hbv_deaths,
             cohort_ly = cohort_ly,
             cohort_size = cohort_size,
             cum_hbv_deaths_2030 = cum_hbv_deaths_2030,
             cum_hbv_deaths_2050 = cum_hbv_deaths_2050,
             cum_hbv_deaths_2100 = cum_hbv_deaths_2100,
             ly_2030 = ly_2030,
             ly_2050 = ly_2050,
             ly_2100 = ly_2100,
             interactions_2030 = interactions_2030,  # NA for no treatment
             interactions_2050 = interactions_2050,  # NA for no treatment
             interactions_2100 = interactions_2100)  # NA for no treatment

out6_ts <- summarise_time_series(out, scenario_label = label, summarise_percentiles = FALSE)

rm(out)
gc()

### Cohort comparison ----

# COUNTERFACTUAL = No monitoring (out3)

# Compare median age at death
age_at_death <- data.frame(rbind(out3$cohort_age_at_death, out4$cohort_age_at_death, out5$cohort_age_at_death,
                                 out6$cohort_age_at_death))
age_at_death <- t(age_at_death[,-1])
colnames(age_at_death) <- c("No monitoring", "Every 10 years", "Every 5 years", "Every year")
boxplot(age_at_death, ylim =c(70,73), ylab = "Mean age at death (years)")

# Extension in age at death compared to no monitoring
age_at_death_ext <- data.frame(cbind(age_at_death[,2]-age_at_death[,1],
                                     age_at_death[,3]-age_at_death[,1],
                                     age_at_death[,4]-age_at_death[,1]))
colnames(age_at_death_ext) <- c("Every 10 years", "Every 5 years", "Every year")
boxplot(age_at_death_ext, ylab = "Extension in mean age at death (years) compared to no monitoring",
        xlab = "Monitoring frequency", ylim = c(0,1.4))

# Compare cohort number of HBV deaths averted compared to no monitoring
cohort_deaths_averted <- rbind(
  calculate_cohort_number_averted(out3$cohort_cum_hbv_deaths, out4$cohort_cum_hbv_deaths, summarise = FALSE),
  calculate_cohort_number_averted(out3$cohort_cum_hbv_deaths, out5$cohort_cum_hbv_deaths, summarise = FALSE),
  calculate_cohort_number_averted(out3$cohort_cum_hbv_deaths, out6$cohort_cum_hbv_deaths, summarise = FALSE))
cohort_deaths_averted_long <- gather(cohort_deaths_averted, key = "sim", value = "value",
                                     -counterfactual, -scenario, -type)

ggplot(cohort_deaths_averted_long[cohort_deaths_averted_long$type == "number_averted",]) +
  geom_boxplot(aes(scenario, value)) +
  ylab("Cumulative number of HBV-related deaths averted by monitoring") +
  labs(title = "Cohort") +
  xlab("Monitoring frequency")

ggplot(cohort_deaths_averted_long[cohort_deaths_averted_long$type == "proportion_averted",]) +
  geom_boxplot(aes(scenario, value)) +
  ylab("Fraction of HBV-related deaths averted by monitoring") +
  xlab("Monitoring frequency") +
  labs(title = "Cohort") +
  ylim(0,1)

# Compare cohort number of life years gained compared to no monitoring
cohort_ly_gained <- rbind(
  calculate_cohort_number_averted(out4$cohort_ly, out3$cohort_ly, summarise = FALSE),
  calculate_cohort_number_averted(out5$cohort_ly, out3$cohort_ly, summarise = FALSE),
  calculate_cohort_number_averted(out6$cohort_ly, out3$cohort_ly, summarise = FALSE))
cohort_ly_gained_long <- gather(cohort_ly_gained, key = "sim", value = "value",
                                     -counterfactual, -scenario, -type)

ggplot(cohort_ly_gained_long[cohort_ly_gained_long$type == "number_averted",]) +
  geom_boxplot(aes(counterfactual, value)) +
  ylab("Life years saved by monitoring") +
  labs(title = "Cohort") +
  xlab("Monitoring frequency")

ggplot(cohort_ly_gained_long[cohort_ly_gained_long$type == "proportion_averted",]) +
  geom_boxplot(aes(counterfactual, value)) +
  ylab("Fraction of life years saved by monitoring") +
  xlab("Monitoring frequency") +
  labs(title = "Cohort")

### Compare population outcomes ----

# Compare population level number of HBV deaths averted
# Counterfactual = status quo
deaths_averted_sq <- rbind(calculate_number_averted(out1$cum_hbv_deaths_2030, out3$cum_hbv_deaths_2030),
                        calculate_number_averted(out1$cum_hbv_deaths_2050, out3$cum_hbv_deaths_2050),
                        calculate_number_averted(out1$cum_hbv_deaths_2100, out3$cum_hbv_deaths_2100),
                        calculate_number_averted(out1$cum_hbv_deaths_2030, out5$cum_hbv_deaths_2030),
                        calculate_number_averted(out1$cum_hbv_deaths_2050, out5$cum_hbv_deaths_2050),
                        calculate_number_averted(out1$cum_hbv_deaths_2100, out5$cum_hbv_deaths_2100))

ggplot(deaths_averted_sq[deaths_averted_sq$type == "number_averted",]) +
  geom_point(aes(x = as.factor(by_year), y = median, group = scenario, colour = scenario),
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(x = as.factor(by_year), ymin = lower, ymax = upper,
                    group = scenario, colour = scenario), position=position_dodge(width = 1)) +
  ylim(0,8000) +
  ylab("Cumulative number of HBV-related deaths averted compared to status quo")

ggplot(deaths_averted_sq[deaths_averted_sq$type == "proportion_averted",]) +
  geom_point(aes(x = as.factor(by_year), y = median, group = scenario, colour = scenario),
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(x = as.factor(by_year), ymin = lower, ymax = upper,
                    group = scenario, colour = scenario), position=position_dodge(width = 1)) +
  ylim(0,1) +
  ylab("Proportion of cumulative HBV-related deaths averted compared to status quo")
#  geom_point(data = cohort_deaths_averted[deaths_averted$type == "proportion_averted",],
#             aes(x = "cohort", y = median)) +
#  geom_errorbar(data = cohort_deaths_averted[deaths_averted$type == "proportion_averted",],
#                aes(x = "cohort",ymin = lower, ymax = upper))


# Compare population level number of HBV deaths averted
# Counterfactual = one-off screening but no monitoring
deaths_averted <- rbind(calculate_number_averted(out3$cum_hbv_deaths_2030, out5$cum_hbv_deaths_2030),
                           calculate_number_averted(out3$cum_hbv_deaths_2050, out5$cum_hbv_deaths_2050),
                           calculate_number_averted(out3$cum_hbv_deaths_2100, out5$cum_hbv_deaths_2100))

ggplot(deaths_averted[deaths_averted$type == "number_averted",]) +
  geom_point(aes(x = as.factor(by_year), y = median)) +
  geom_errorbar(aes(x = as.factor(by_year), ymin = lower, ymax = upper)) +
  ylab("Cumulative number of HBV-related deaths averted by monitoring") +
  geom_point(data = cohort_deaths_averted[deaths_averted$type == "number_averted",],
             aes(x = "cohort", y = median)) +
  geom_errorbar(data = cohort_deaths_averted[deaths_averted$type == "number_averted",],
               aes(x = "cohort",ymin = lower, ymax = upper))

ggplot(deaths_averted[deaths_averted$type == "proportion_averted",]) +
  geom_point(aes(x = as.factor(by_year), y = median)) +
  geom_errorbar(aes(x = as.factor(by_year), ymin = lower, ymax = upper)) +
  ylim(0,1) +
  ylab("Proportion of cumulative HBV-related deaths averted by monitoring") +
  geom_point(data = cohort_deaths_averted[deaths_averted$type == "proportion_averted",],
             aes(x = "cohort", y = median)) +
  geom_errorbar(data = cohort_deaths_averted[deaths_averted$type == "proportion_averted",],
                aes(x = "cohort",ymin = lower, ymax = upper))


# need to somehow record number of people in the cohort! E.g. deaths per population!

# Compare population level number of  life years gained
# Counterfactual = status_quo
ly_gained_sq <- rbind(calculate_number_averted(out3$ly_2030, out1$ly_2030),
                      calculate_number_averted(out3$ly_2050, out1$ly_2050),
                      calculate_number_averted(out3$ly_2100, out1$ly_2100),
                      calculate_number_averted(out5$ly_2030, out1$ly_2030),
                      calculate_number_averted(out5$ly_2050, out1$ly_2050),
                      calculate_number_averted(out5$ly_2100, out1$ly_2100))

ggplot(ly_gained_sq[ly_gained_sq$type == "number_averted",]) +
  geom_point(aes(x = as.factor(by_year), y = median, group = counterfactual, colour = counterfactual),
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(x = as.factor(by_year), ymin = lower, ymax = upper,
                    group = counterfactual, colour = counterfactual), position=position_dodge(width = 1)) +
  ylab("Life years gained compared to status quo")

ggplot(ly_gained_sq[ly_gained_sq$type == "proportion_averted",]) +
  geom_point(aes(x = as.factor(by_year), y = median, group =counterfactual, colour = counterfactual),
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(x = as.factor(by_year), ymin = lower, ymax = upper,
                    group = counterfactual, colour = counterfactual), position=position_dodge(width = 1)) +
  ylab("Fraction of life years gained compared to status quo")


# Compare population level number of  life years gained
# Counterfactual = one-off screening but no monitoring
ly_gained <- rbind(calculate_number_averted(out5$ly_2030, out3$ly_2030),
                   calculate_number_averted(out5$ly_2050, out3$ly_2050),
                   calculate_number_averted(out5$ly_2100, out3$ly_2100))

ggplot(ly_gained[ly_gained$type == "number_averted",]) +
  geom_point(aes(x = as.factor(by_year), y = median)) +
  geom_errorbar(aes(x = as.factor(by_year), ymin = lower, ymax = upper)) +
  geom_point(data = cohort_ly_gained[cohort_ly_gained$type == "number_averted",],
             aes(x = "cohort", y = median)) +
  geom_errorbar(data = cohort_ly_gained[cohort_ly_gained$type == "number_averted",],
                aes(x = "cohort",ymin = lower, ymax = upper)) +
  ylab("Number of life-years gained by monitoring")

ggplot(ly_gained[ly_gained$type == "proportion_averted",]) +
  geom_point(aes(x = as.factor(by_year), y = median)) +
  geom_errorbar(aes(x = as.factor(by_year), ymin = lower, ymax = upper)) +
  geom_point(data = cohort_ly_gained[cohort_ly_gained$type == "proportion_averted",],
             aes(x = "cohort", y = median)) +
  geom_errorbar(data = cohort_ly_gained[cohort_ly_gained$type == "proportion_averted",],
                aes(x = "cohort",ymin = lower, ymax = upper)) +
  ylab("Fraction of life-years gained by monitoring")

# Calculate outcomes per healthcare interaction

# HBV deaths averted per incremental treatment initiations - compared to status quo
# Number of HBV deaths averted compared to status quo
quantile(calculate_number_averted(out1$cum_hbv_deaths_2030, out3$cum_hbv_deaths_2030, summarise = FALSE)[1,-c(1:5)],
         prob = c(0.025,0.5,0.975))
quantile(calculate_number_averted(out1$cum_hbv_deaths_2030, out5$cum_hbv_deaths_2030, summarise = FALSE)[1,],
         prob = c(0.025,0.5,0.975))

# Number of HBV deaths averted per incremental treatment initiations compared to status quo
deaths_averted_per_int3 <- data.frame(rbind(
  c(by_year = 2030, unlist(calculate_number_averted(out1$cum_hbv_deaths_2030, out3$cum_hbv_deaths_2030, summarise = FALSE)[1,-c(1:5)]/
  out3$interactions_2030$total_treated[,-c(1:3)])),
  c(by_year = 2050, unlist(calculate_number_averted(out1$cum_hbv_deaths_2050, out3$cum_hbv_deaths_2050, summarise = FALSE)[1,-c(1:5)]/
                             out3$interactions_2050$total_treated[,-c(1:3)])),
  c(by_year = 2100, unlist(calculate_number_averted(out1$cum_hbv_deaths_2100, out3$cum_hbv_deaths_2100, summarise = FALSE)[1,-c(1:5)]/
                             out3$interactions_2100$total_treated[,-c(1:3)]))
  ))
deaths_averted_per_int3$scenario <- "screen1_monit0"

deaths_averted_per_int5 <- data.frame(rbind(
  c(by_year = 2030, unlist(calculate_number_averted(out1$cum_hbv_deaths_2030, out5$cum_hbv_deaths_2030, summarise = FALSE)[1,-c(1:5)]/
                             out5$interactions_2030$total_treated[,-c(1:3)])),
  c(by_year = 2050, unlist(calculate_number_averted(out1$cum_hbv_deaths_2050, out5$cum_hbv_deaths_2050, summarise = FALSE)[1,-c(1:5)]/
                             out5$interactions_2050$total_treated[,-c(1:3)])),
  c(by_year = 2100, unlist(calculate_number_averted(out1$cum_hbv_deaths_2100, out5$cum_hbv_deaths_2100, summarise = FALSE)[1,-c(1:5)]/
                             out5$interactions_2100$total_treated[,-c(1:3)]))
))
deaths_averted_per_int5$scenario <- "screen1_monit5"

deaths_averted_per_int <- rbind(
  gather(deaths_averted_per_int3, key = "sim", value = "value", -scenario, -by_year),
  gather(deaths_averted_per_int5, key = "sim", value = "value", -scenario, -by_year))

ggplot(data = deaths_averted_per_int, aes(as.factor(by_year), value)) +
  geom_boxplot(aes(fill = scenario)) +
  ylab("Cumulative HBV-related deaths averted per incremental treatment initiations")

# HBV deaths averted per incremental healthcare interactions - compared to status quo
# Number of HBV deaths averted compared to status quo
quantile(calculate_number_averted(out1$cum_hbv_deaths_2030, out3$cum_hbv_deaths_2030, summarise = FALSE)[1,-c(1:5)],
         prob = c(0.025,0.5,0.975))
quantile(calculate_number_averted(out1$cum_hbv_deaths_2030, out5$cum_hbv_deaths_2030, summarise = FALSE)[1,],
         prob = c(0.025,0.5,0.975))

# Number of HBV deaths averted per incremental assessment compared to status quo
deaths_averted_per_int3_2 <- data.frame(rbind(
  c(by_year = 2030, unlist(calculate_number_averted(out1$cum_hbv_deaths_2030, out3$cum_hbv_deaths_2030, summarise = FALSE)[1,-c(1:5)]/
                             out3$interactions_2030$total_interactions[,-c(1:3)])),
  c(by_year = 2050, unlist(calculate_number_averted(out1$cum_hbv_deaths_2050, out3$cum_hbv_deaths_2050, summarise = FALSE)[1,-c(1:5)]/
                             out3$interactions_2050$total_interactions[,-c(1:3)])),
  c(by_year = 2100, unlist(calculate_number_averted(out1$cum_hbv_deaths_2100, out3$cum_hbv_deaths_2100, summarise = FALSE)[1,-c(1:5)]/
                             out3$interactions_2100$total_interactions[,-c(1:3)]))
))
deaths_averted_per_int3_2$scenario <- "screen1_monit0"

deaths_averted_per_int5_2 <- data.frame(rbind(
  c(by_year = 2030, unlist(calculate_number_averted(out1$cum_hbv_deaths_2030, out5$cum_hbv_deaths_2030, summarise = FALSE)[1,-c(1:5)]/
                             out5$interactions_2030$total_interactions[,-c(1:3)])),
  c(by_year = 2050, unlist(calculate_number_averted(out1$cum_hbv_deaths_2050, out5$cum_hbv_deaths_2050, summarise = FALSE)[1,-c(1:5)]/
                             out5$interactions_2050$total_interactions[,-c(1:3)])),
  c(by_year = 2100, unlist(calculate_number_averted(out1$cum_hbv_deaths_2100, out5$cum_hbv_deaths_2100, summarise = FALSE)[1,-c(1:5)]/
                             out5$interactions_2100$total_interactions[,-c(1:3)]))
))
deaths_averted_per_int5_2$scenario <- "screen1_monit5"

deaths_averted_per_int_2 <- rbind(
  gather(deaths_averted_per_int3_2, key = "sim", value = "value", -scenario, -by_year),
  gather(deaths_averted_per_int5_2, key = "sim", value = "value", -scenario, -by_year))

ggplot(data = deaths_averted_per_int_2, aes(as.factor(by_year), value)) +
  geom_boxplot(aes(fill = scenario)) +
  ylab("Cumulative HBV-related deaths averted per incremental healthcare interactions")

# Need to think about what I would like to represent:
# compared to monitoring, not monitoring averts more deaths per incremental clinical assessments and treatment initiations,
# and averts a similar number of deaths per incremental healthcare interaction
# but Maud has said than screening people is relatively easy whereas the clinical assessment is difficult

# Life years gained per incremental healthcare interactions - compared to status quo
# Number of HBV deaths averted compared to status quo
quantile(calculate_number_averted(out1$cum_hbv_deaths_2030, out3$cum_hbv_deaths_2030, summarise = FALSE)[1,-c(1:5)],
         prob = c(0.025,0.5,0.975))
quantile(calculate_number_averted(out1$cum_hbv_deaths_2030, out5$cum_hbv_deaths_2030, summarise = FALSE)[1,],
         prob = c(0.025,0.5,0.975))

ly_gained_per_int3_2 <- data.frame(rbind(
  c(by_year = 2030, unlist(calculate_number_averted(out1$ly_2030, out3$ly_2030, summarise = FALSE)[1,-c(1:5)]*-1/
                             out3$interactions_2030$total_interactions[,-c(1:3)])),
  c(by_year = 2050, unlist(calculate_number_averted(out1$ly_2050, out3$ly_2050, summarise = FALSE)[1,-c(1:5)]*-1/
                             out3$interactions_2050$total_interactions[,-c(1:3)])),
  c(by_year = 2100, unlist(calculate_number_averted(out1$ly_2100, out3$ly_2100, summarise = FALSE)[1,-c(1:5)]*-1/
                             out3$interactions_2100$total_interactions[,-c(1:3)]))
))
ly_gained_per_int3_2$scenario <- "screen1_monit0"

ly_gained_per_int5_2 <- data.frame(rbind(
  c(by_year = 2030, unlist(calculate_number_averted(out1$ly_2030, out5$ly_2030, summarise = FALSE)[1,-c(1:5)]*-1/
                             out5$interactions_2030$total_interactions[,-c(1:3)])),
  c(by_year = 2050, unlist(calculate_number_averted(out1$ly_2050, out5$ly_2050, summarise = FALSE)[1,-c(1:5)]*-1/
                             out5$interactions_2050$total_interactions[,-c(1:3)])),
  c(by_year = 2100, unlist(calculate_number_averted(out1$ly_2100, out5$ly_2100, summarise = FALSE)[1,-c(1:5)]*-1/
                             out5$interactions_2100$total_interactions[,-c(1:3)]))
))
ly_gained_per_int5_2$scenario <- "screen1_monit5"

ly_gained_per_int_2 <- rbind(
  gather(ly_gained_per_int3_2, key = "sim", value = "value", -scenario, -by_year),
  gather(ly_gained_per_int5_2, key = "sim", value = "value", -scenario, -by_year))

ggplot(data = ly_gained_per_int_2, aes(as.factor(by_year), value)) +
  geom_boxplot(aes(fill = scenario)) +
  ylab("Life-years gained per incremental healthcare interactions")

# Another way of showing this is deaths averted on the y axis and median number of healthcare interactions

# Median values for no monitoring, monitoring every 5 years
test <- data.frame(deaths_averted_by_2050 = c(2597.0421888, 3277.3994132),
                   incremental_interactions_by_2050 = c(576298, 644428),
                   incremental_assessments_by_2050 = c(52498, 216173))

plot(x = test$incremental_interactions_by_2050, y = test$deaths_averted_by_2050, ylim = c(0,3300), xlim = c(50000,650000))
points(x = test$incremental_assessments_by_2050, y = test$deaths_averted_by_2050, ylim = c(0,3300), xlim = c(50000,650000), col = "red")

###

# Notes on cohort:
# The cohort we are following are those screened individuals who immediately engage with care
