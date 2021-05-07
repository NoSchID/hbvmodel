# Sensitivity analysis on cost-effectiveness of screening and monitoring
# Get code and load functions from monitoring_frequency_by_age_incremental_analysis

require(here)  # for setting working directory
require(ggplot2)
#require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
library(BCEA)
library(sensitivity)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Load data ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/sensitivity_analysis/"

out_path_monit <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"


# Vary progression to HCC on treatment (thccr_it, thccr_chb, thccr_cc, thccr_dcc)

# Status quo is the same
out2 <- readRDS(paste0(out_path_monit, "out2_status_quo_301120.rds"))
out2 <- out2[[1]]
out1_it <- readRDS(paste0(out_path_monit, "a1_it_out1_status_quo_cohort_200121.rds"))
out1_it <- out1_it[[1]]

# Default simulations for comparison:
out3_it <- readRDS(paste0(out_path_monit, "a1_it_out3_screen_2020_monit_0_180121.rds"))
out3_it <- out3_it[[1]]   # No monitoring
monit_out7 <- readRDS(paste0(out_path_monit, "a1_it_screen_2020_monit_out7_050321.rds"))
monit_out7 <- monit_out7[[1]]
out5_it <- readRDS(paste0(out_path_monit, "a1_it_out6_screen_2020_monit_1_240221.rds"))
out5_it <- out5_it[[1]]

# For PRCC:
out6_it <- readRDS(paste0(out_path_monit, "a1_it_out6_screen_2020_monit_1_130121.rds"))
out6_it <- out6_it[[1]]  # 1 year

# Option 1: increase all thccr parameters to 0.4
out3_s1 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_treatment_sensitivity_upper_030221.rds"))
out3_s1 <- out3_s1[[1]]
monit_out7_s1 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_treatment_sensitivity_upper_030221.rds"))
monit_out7_s1 <- monit_out7_s1[[1]]

# Option 2: increase thccr_cc and thccr_dcc to 0.8
out3_s2 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_treatment_sensitivity_vary_cirrhosis_upper_030221.rds"))
out3_s2 <- out3_s2[[1]]
monit_out7_s2 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_treatment_sensitivity_vary_cirrhosis_upper_030221.rds"))
monit_out7_s2 <- monit_out7_s2[[1]]

# Option 3: increase thccr_it and thccr_chb to 0.8
out3_s3 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_treatment_sensitivity_vary_chb_upper_030221.rds"))
out3_s3 <- out3_s3[[1]]
monit_out7_s3 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_treatment_sensitivity_vary_chb_upper_030221.rds"))
monit_out7_s3 <- monit_out7_s3[[1]]

# One-way (individual parameter variation)

# Vary tmu_dcc to 0.08 and 0.35
out3_tmu_dcc_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_tmu_dcc_lower_040221.rds"))
out3_tmu_dcc_lower <- out3_tmu_dcc_lower[[1]]
out3_tmu_dcc_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_tmu_dcc_upper_040221.rds"))
out3_tmu_dcc_upper <- out3_tmu_dcc_upper[[1]]
monit_out7_tmu_dcc_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_tmu_dcc_lower_040221.rds"))
monit_out7_tmu_dcc_lower <- monit_out7_tmu_dcc_lower[[1]]
monit_out7_tmu_dcc_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_tmu_dcc_upper_040221.rds"))
monit_out7_tmu_dcc_upper <- monit_out7_tmu_dcc_upper[[1]]
# Vary thccr_chb to 0.1 and 0.8
out3_thccr_chb_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_thccr_chb_lower_040221.rds"))
out3_thccr_chb_lower <- out3_thccr_chb_lower[[1]]
out3_thccr_chb_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_thccr_chb_upper_040221.rds"))
out3_thccr_chb_upper <- out3_thccr_chb_upper[[1]]
monit_out7_thccr_chb_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_thccr_chb_lower_040221.rds"))
monit_out7_thccr_chb_lower <- monit_out7_thccr_chb_lower[[1]]
monit_out7_thccr_chb_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_thccr_chb_upper_040221.rds"))
monit_out7_thccr_chb_upper <- monit_out7_thccr_chb_upper[[1]]
# Vary thccr_it to 0.1 and 0.8
out3_thccr_it_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_thccr_it_lower_080221.rds"))
out3_thccr_it_lower <- out3_thccr_it_lower[[1]]
out3_thccr_it_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_thccr_it_upper_080221.rds"))
out3_thccr_it_upper <- out3_thccr_it_upper[[1]]
monit_out7_thccr_it_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_thccr_it_lower_080221.rds"))
monit_out7_thccr_it_lower <- monit_out7_thccr_it_lower[[1]]
monit_out7_thccr_it_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_thccr_it_upper_090221.rds"))
monit_out7_thccr_it_upper <- monit_out7_thccr_it_upper[[1]]
# Vary thccr_cc to 0.1 and 0.8
out3_thccr_cc_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_thccr_cc_lower_080221.rds"))
out3_thccr_cc_lower <- out3_thccr_cc_lower[[1]]
out3_thccr_cc_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_thccr_cc_upper_080221.rds"))
out3_thccr_cc_upper <- out3_thccr_cc_upper[[1]]
monit_out7_thccr_cc_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_thccr_cc_lower_090221.rds"))
monit_out7_thccr_cc_lower <- monit_out7_thccr_cc_lower[[1]]
monit_out7_thccr_cc_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_thccr_cc_upper_090221.rds"))
monit_out7_thccr_cc_upper <- monit_out7_thccr_cc_upper[[1]]
# Vary thccr_dcc to 0.1 and 0.8
out3_thccr_dcc_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_thccr_dcc_lower_120221.rds"))
out3_thccr_dcc_lower <- out3_thccr_dcc_lower[[1]]
out3_thccr_dcc_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_thccr_dcc_upper_080221.rds"))
out3_thccr_dcc_upper <- out3_thccr_dcc_upper[[1]]
monit_out7_thccr_dcc_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_thccr_dcc_lower_080221.rds"))
monit_out7_thccr_dcc_lower <- monit_out7_thccr_dcc_lower[[1]]
monit_out7_thccr_dcc_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_thccr_dcc_upper_080221.rds"))
monit_out7_thccr_dcc_upper <- monit_out7_thccr_dcc_upper[[1]]

# Vary screening coverage to 0.5 and 1
out3_screening_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_screening_coverage_lower_050421.rds"))
out3_screening_lower <- out3_screening_lower[[1]]
out3_screening_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_screening_coverage_higher_010421.rds"))
out3_screening_upper <- out3_screening_upper[[1]]
monit_out7_screening_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_screening_coverage_lower_010421.rds"))
monit_out7_screening_lower <- monit_out7_screening_lower[[1]]
monit_out7_screening_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_screening_coverage_higher_310321.rds"))
monit_out7_screening_upper <- monit_out7_screening_upper[[1]]
# Vary linkage to care to 0.5 and 1
out3_assessment_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_link_to_care_prob_lower_310321.rds"))
out3_assessment_lower <- out3_assessment_lower[[1]]
out3_assessment_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_link_to_care_prob_higher_310321.rds"))
out3_assessment_upper <- out3_assessment_upper[[1]]
monit_out7_assessment_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_link_to_care_prob_lower_310321.rds"))
monit_out7_assessment_lower <- monit_out7_assessment_lower[[1]]
monit_out7_assessment_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_link_to_care_prob_higher_310321.rds"))
monit_out7_assessment_upper <- monit_out7_assessment_upper[[1]]
# Vary treatment initiation prob to 0.5
out3_treatment_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_treatment_prob_lower_310321.rds"))
out3_treatment_lower <- out3_treatment_lower[[1]]
monit_out7_treatment_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_treatment_prob_lower_310321.rds"))
monit_out7_treatment_lower <- monit_out7_treatment_lower[[1]]
# Vary monitoring prob to 0.5 and 1
monit_out7_monitoring_lower <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_monitoring_prob_lower_310321.rds"))
monit_out7_monitoring_lower <- monit_out7_monitoring_lower[[1]]
monit_out7_monitoring_upper <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_monitoring_prob_higher_310321.rds"))
monit_out7_monitoring_upper <- monit_out7_monitoring_upper[[1]]


# One-way sensitivity analysis of cost (median ICERs and uncertainty bounds)
cost_sensitivity_icer <- read.csv(file=paste0(out_path, "cost_sensitivity_results_170321.csv"))

# Technically out3 is dominated, but get ICER compared to SQ anyway for comparison purposes
# The ICER for monit_sim7 compared to monit_sim6 in the other file is: 338 (161-844)
# If monit_sim6 was excluded it would be: 271 (163-519)

# Sensitivity analysis for NO INFECTIVITY from treated carriers:
out3_no_inf <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_0_no_treated_infectivity_100321.rds"))
out3_no_inf <- out3_no_inf[[1]]

# Sensitivity analysis for impact of vaccination compared to no vacc:
out_no_vacc <- readRDS(paste0(out_path, "no_vacc_scenario_170321.rds"))
out_no_vacc <- out_no_vacc[[1]]

# PRCC on default simulation (need individual ICERS) ----
annual_discounting_rate <- 0.03

object_list <- list(out3_it, monit_out7,
                    out5_it)

# Extract interactions and person-years on treatment, and HBV-related deaths
# ad DALYs averted, in a for loop
age_interactions <- list()
age_interactions_py_on_treatment <- list()
age_hbv_deaths_averted <- list()
age_dalys_averted <- list()

for (i in 1:length(object_list)) {
  age_interactions[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          assemble_discounted_interactions_for_monitoring_frequencies(object_list[[i]],
                                                                      no_monitoring_object = out3_it))
  age_interactions_py_on_treatment[[i]] <-
    data.frame(scenario = object_list[[i]]$cohort_age_at_death$scenario,
               discount_outcome_2020_to_2100(scenario_object=object_list[[i]],
                                             object_to_subtract=NULL,
                                             outcome="py_on_treatment",
                                             yearly_discount_rate=annual_discounting_rate))
  age_hbv_deaths_averted[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          discount_outcome_2020_to_2100(scenario_object=out2,
                                        object_to_subtract=object_list[[i]],
                                        outcome="cum_hbv_deaths",
                                        yearly_discount_rate=annual_discounting_rate))
  age_dalys_averted[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          discount_outcome_2020_to_2100(scenario_object=out2,
                                        object_to_subtract=object_list[[i]],
                                        outcome="dalys",
                                        yearly_discount_rate=annual_discounting_rate))
}
age_interactions <- do.call("rbind", age_interactions)
age_interactions_py_on_treatment <- do.call("rbind", age_interactions_py_on_treatment)
age_hbv_deaths_averted <- do.call("rbind", age_hbv_deaths_averted)
age_dalys_averted <- do.call("rbind", age_dalys_averted)

age_interactions_py_on_treatment$sim <- gsub("[^0-9]", "", age_interactions_py_on_treatment$sim)
age_hbv_deaths_averted$sim <- gsub("[^0-9]", "", age_hbv_deaths_averted$sim)
colnames(age_hbv_deaths_averted)[colnames(age_hbv_deaths_averted) == "cum_hbv_deaths"] <-
  "value"

age_dalys_averted$sim <- gsub("[^0-9]", "", age_dalys_averted$sim)
colnames(age_dalys_averted)[colnames(age_dalys_averted) == "dalys"] <-
  "value"

age_df <- create_incremental_plot_df(interactions_df=age_interactions,
                                     py_on_treatment_df=age_interactions_py_on_treatment,
                                     deaths_averted_df=age_hbv_deaths_averted,
                                     ly_saved_df = age_dalys_averted, # replace LY by DALYs
                                     hbsag_test_cost = 8.3,
                                     clinical_assessment_cost = 33, #84.4,
                                     monitoring_assessment_cost = 25.5, #40.1,
                                     treatment_py_cost = 66.5,#60,
                                     #scenario_labels_obj = scenario_labels,
                                     ref_label = "No treatment")
colnames(age_df)[colnames(age_df)=="ly_saved"] <- "dalys_averted"

icer_list <- list()

for(i in 1:183) {
  print(i)
  icer_list[[i]] <- age_df[which(age_df$sim==
                                   unique(age_df$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result

# WORK ON icer_df
# Extract sim-specific ICERs for out3:
icer_monit_0 <- icer_df[icer_df$scenario=="screen_2020_monit_0",]
icer_monit_sim7 <- icer_df[icer_df$scenario=="screen_2020_monit_sim7",]
icer_monit_5 <- icer_df[icer_df$scenario=="screen_2020_monit_5",]



# Load parameter sets chosen based on kmeans clustering
load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans

# Check order of sims is the same:
all.equal(icer_monit_0$sim,rownames(params_mat_accepted_kmeans))
all.equal(icer_monit_sim7$sim,rownames(params_mat_accepted_kmeans))

all.equal(icer_monit_5$sim,rownames(params_mat_accepted_kmeans))

library(epiR)
prcc_icer_monit_0 <- epi.prcc(cbind(params_mat_accepted_kmeans, icer_monit_0$icer),
                              sided.test = 2, conf.level = 0.95)
# Checked that this gives the same est as pcc package
prcc_icer_monit_sim7 <- epi.prcc(cbind(params_mat_accepted_kmeans, icer_monit_sim7$icer),
                              sided.test = 2, conf.level = 0.95)

prcc_dalys_monit_0 <- epi.prcc(cbind(params_mat_accepted_kmeans, icer_monit_0$dalys_averted),
                              sided.test = 2, conf.level = 0.95)
prcc_dalys_monit_sim7 <- epi.prcc(cbind(params_mat_accepted_kmeans, icer_monit_sim7$dalys_averted-icer_monit_0$dalys_averted),
                                 sided.test = 2, conf.level = 0.95)

prcc_df <- data.frame(parameter= colnames(params_mat_accepted_kmeans),
                      out3_dalys_prcc = prcc_dalys_monit_0$est,
                      out3_dalys_p_value = prcc_dalys_monit_0$p.value,
          out3_icer_prcc = prcc_icer_monit_0$est,
           out3_icer_p_value = prcc_icer_monit_0$p.value,
          monit_sim7_dalys_prcc = prcc_dalys_monit_sim7$est,
          monit_sim7_dalys_p_value = prcc_dalys_monit_sim7$p.value,
          monit_sim7_icer_prcc = prcc_icer_monit_sim7$est,
          monit_sim7_icer_p_value = prcc_icer_monit_sim7$p.value)
prcc_df <- arrange(prcc_df, -abs(out3_icer_prcc))

parameter_names <- list(
  # Transmission parameters
  "beta1" = "b1",
  "beta2" = "b2",
  "beta3" = "b3",
  "Relative infectiousness in HBeAg+" = "alpha",
  "MTCT risk from HBeAg+ mother" = "mtct_prob_e",
  "MTCT risk from HBeAg- mother" = "mtct_prob_s",
  "Risk of chronic carriage at birth" = "p_chronic_in_mtct",
  "Coefficient for risk of chronic carriage (cr)" = "p_chronic_function_r",
  "Coefficient for risk of chronic carriage (cs)" = "p_chronic_function_s",
  "Vaccine efficacy" = "vacc_eff",
  # Natural history
  "Rate from HBeAg+ infection to CHB at age 0" = "pr_it_ir",
  "Rate from HBeAg+ CHB to HBeAg- infection at age 0" = "pr_ir_ic",
  "Coefficient for progression through HBeAg+ compartments" = "eag_prog_function_rate",
  "Rate from HBeAg+ to HBeAg- CHB" = "pr_ir_enchb",
  "Parameter for HBsAg loss" = "sag_loss_slope",
  "Rate from HBeAg- infection to CHB" = "pr_ic_enchb",
  # Liver disease
  "Rate from HBeAg+ CHB to CC in women" = "pr_ir_cc_female",
  "Minimum age for cirrhosis (HBeAg+)" = "pr_ir_cc_age_threshold",
  "Rate from HBeAg- CHB to CC in women" = "pr_enchb_cc_female",
  "Rate ratio for cirrhosis in men" = "cirrhosis_male_cofactor",
  "Rate of decompensation" = "pr_cc_dcc",
  "Coefficient for progression to HCC in women" = "cancer_prog_coefficient_female",
  "Minimum age for HCC" = "cancer_age_threshold",
  "Rate ratio for HCC in men" = "cancer_male_cofactor",
  "Rate ratio for HCC in HBeAg+ infection" = "hccr_it",
  "Rate ratio for HCC in HBeAg+ CHB" = "hccr_ir",
  "Rate ratio for HCC in HBeAg- CHB" = "hccr_enchb",
  "Rate ratio for HCC in CC" = "hccr_cc",
  "Rate from DCC to HCC" = "hccr_dcc",
  "Mortality rate from CC" = "mu_cc",
  "Mortality rate from DCC" = "mu_dcc",
  "Mortality rate from HCC" = "mu_hcc")

prcc_df$parameter <- factor(prcc_df$parameter)
levels(prcc_df$parameter) <- parameter_names

#write.csv(prcc_df, file = "prcc_estimates_030521.csv")

abs(prcc_df$monit_sim7_icer_prcc[rank(abs(prcc_df$monit_sim7_icer_prcc))==23])

# PREVIOUS ANALYSIS

# Run PRCC
prcc_icer_monit_0 <- pcc(params_mat_accepted_kmeans,
                         icer_monit_0$icer,
                         rank = TRUE, nboot = 100)   # checked and there was no difference between 100 and 1000 nboot
plot(prcc_icer_monit_0)
abline(h=0)

# Subset only significant parameters (95% CI not including 0):
icer_monit_0_significant_parms <- data.frame(parm = rownames(prcc_icer_monit_0$PRCC)[
  which(sign(prcc_icer_monit_0$PRCC$`max. c.i.`)==sign(prcc_icer_monit_0$PRCC$`min. c.i.`))],
  prcc_mean = prcc_icer_monit_0$PRCC$original[
    which(sign(prcc_icer_monit_0$PRCC$`max. c.i.`)==sign(prcc_icer_monit_0$PRCC$`min. c.i.`))],
  prcc_ci_lower = prcc_icer_monit_0$PRCC$`min. c.i.`[
    which(sign(prcc_icer_monit_0$PRCC$`max. c.i.`)==sign(prcc_icer_monit_0$PRCC$`min. c.i.`))],
  prcc_ci_upper = prcc_icer_monit_0$PRCC$`max. c.i.`[
    which(sign(prcc_icer_monit_0$PRCC$`max. c.i.`)==sign(prcc_icer_monit_0$PRCC$`min. c.i.`))])
icer_monit_0_significant_parms <- arrange(icer_monit_0_significant_parms, -abs(prcc_mean))

icer_monit_0_significant_parms$sign <- "Positive"
icer_monit_0_significant_parms$sign[which(icer_monit_0_significant_parms$prcc_mean<0)] <-
  "Negative"

# Run PRCC for monit_sim7 (this is ICER compared to no monitoring, not base case)
prcc_icer_monit_sim7 <- pcc(params_mat_accepted_kmeans,
                            icer_monit_sim7$icer,
                            rank = TRUE, nboot = 100)   # checked and there was no difference between 100 and 1000 nboot
plot(prcc_icer_monit_sim7)
abline(h=0)

# Subset only significant parameters (95% CI not including 0):
icer_monit_sim7_significant_parms <- data.frame(parm = rownames(prcc_icer_monit_sim7$PRCC)[
  which(sign(prcc_icer_monit_sim7$PRCC$`max. c.i.`)==sign(prcc_icer_monit_sim7$PRCC$`min. c.i.`))],
  prcc_mean = prcc_icer_monit_sim7$PRCC$original[
    which(sign(prcc_icer_monit_sim7$PRCC$`max. c.i.`)==sign(prcc_icer_monit_sim7$PRCC$`min. c.i.`))],
  prcc_ci_lower = prcc_icer_monit_sim7$PRCC$`min. c.i.`[
    which(sign(prcc_icer_monit_sim7$PRCC$`max. c.i.`)==sign(prcc_icer_monit_sim7$PRCC$`min. c.i.`))],
  prcc_ci_upper = prcc_icer_monit_sim7$PRCC$`max. c.i.`[
    which(sign(prcc_icer_monit_sim7$PRCC$`max. c.i.`)==sign(prcc_icer_monit_sim7$PRCC$`min. c.i.`))])
icer_monit_sim7_significant_parms <- arrange(icer_monit_sim7_significant_parms, -abs(prcc_mean))

icer_monit_sim7_significant_parms$sign <- "Positive"
icer_monit_sim7_significant_parms$sign[which(icer_monit_sim7_significant_parms$prcc_mean<0)] <-
  "Negative"

# IGNORE!
# Monitor all ages every 5 years compared to only youngest ages
prcc_icer_monit_5 <- pcc(params_mat_accepted_kmeans,
                         icer_monit_5$icer,
                         rank = TRUE, nboot = 100)   # checked and there was no difference between 100 and 1000 nboot
plot(prcc_icer_monit_5)
abline(h=0)

# Subset only significant parameters (95% CI not including 0):
icer_monit_5_significant_parms <- data.frame(parm = rownames(prcc_icer_monit_5$PRCC)[
  which(sign(prcc_icer_monit_5$PRCC$`max. c.i.`)==sign(prcc_icer_monit_5$PRCC$`min. c.i.`))],
  prcc_mean = prcc_icer_monit_5$PRCC$original[
    which(sign(prcc_icer_monit_5$PRCC$`max. c.i.`)==sign(prcc_icer_monit_5$PRCC$`min. c.i.`))],
  prcc_ci_lower = prcc_icer_monit_5$PRCC$`min. c.i.`[
    which(sign(prcc_icer_monit_5$PRCC$`max. c.i.`)==sign(prcc_icer_monit_5$PRCC$`min. c.i.`))],
  prcc_ci_upper = prcc_icer_monit_5$PRCC$`max. c.i.`[
    which(sign(prcc_icer_monit_5$PRCC$`max. c.i.`)==sign(prcc_icer_monit_5$PRCC$`min. c.i.`))])
icer_monit_5_significant_parms <- arrange(icer_monit_5_significant_parms, -abs(prcc_mean))

# Not many parameters are significantly associated with the all age vs age-specific monitoring
# ICER:
# pr_ic_enchb,  pr_enchb_cc_female, mu_cc, mu_hcc and cancer_prog_coefficient_female
# are over 25%. Afte that: cirrhosis_male_cofactor, hccr_enchb, hccr_cc, hccr_dcc

# I also checked for ICER of 5-yearly monitoring across all ages vs no monitoring directly,
# in which case the influential parameters are almost the same for the age-specific monitoring ICER.



## PLOTS

# Plot 10 most influential parameters for ICER of basic programme:
ggplot(subset(icer_monit_0_significant_parms, rank(abs(prcc_mean))>10)) +
  geom_col(aes(x=reorder(parm, abs(prcc_mean)), y = abs(prcc_mean), fill = sign)) +
  coord_flip() +
  ylab("Partial rank correlation coefficient for CER\nof 2020 screen & treat without monitoring") +
  xlab("Parameter or transition") +
  scale_x_discrete(labels = c("pr_ic_enchb" = "Progression from inactive carrier\nto HBeAg-negative CHB",
                              "pr_enchb_cc_female" = "Progression from\nHBeAg-negative CHB to CC",
                              "pr_ir_cc_female" = "Progression from\nHBeAg-positive CHB to CC",
                              "pr_ir_enchb" = "Progression from HBeAg-positive CHB\nto HBeAg-negative ENCHB",
                              "mtct_prob_e" = "MTCT risk from\nHBeAg-positive mother",
                              "pr_ir_ic" = "Progression from HBeAg-positive CHB\nto inactive carrier",
                              "p_chronic_function_r" = "Parameter of age-specific\nrisk of chronic carriage function",
                              "cirrhosis_male_cofactor" = "Rate ratio for progression\nto cirrhosis in men vs. women",
                              "b1" = "Horizontal transmission coefficient for\ntransmission among children aged <5 years",
                              "p_chronic_in_mtct" = "Risk of chronic carriage\nat birth")) +
  labs(fill="Sign", title = "10 parameters with highest PRCC") +
  theme_classic()

# Plot 10 most influential parameters for ICER of optimal monitoring programme:
ggplot(subset(icer_monit_sim7_significant_parms, rank(abs(prcc_mean))>7)) +
  geom_col(aes(x=reorder(parm, abs(prcc_mean)), y = abs(prcc_mean), fill = sign)) +
  coord_flip() +
  ylab("Partial rank correlation coefficient for ICER\nof 2020 screen & treat with monitoring every 5 years\nin <45 year olds vs. no monitoring") +
  xlab("Parameter or transition") +
  scale_x_discrete(labels = c("pr_ic_enchb" = "Progression from inactive carrier\nto HBeAg-negative CHB",
                              "pr_enchb_cc_female" = "Progression from\nHBeAg-negative CHB to CC",
                              "pr_ir_cc_female" = "Progression from\nHBeAg-positive CHB to CC",
                              "pr_ir_enchb" = "Progression from HBeAg-positive CHB\nto HBeAg-negative ENCHB",
                              "mtct_prob_e" = "MTCT risk from\nHBeAg-positive mother",
                              "mu_cc" = "Mortality rate from CC",
                              "pr_cc_dcc" = "Rate of decompensation",
                              "cirrhosis_male_cofactor" = "Rate ratio for progression\nto cirrhosis in men vs. women",
                              "pr_it_ir" = "Progression from HBeAg-positive infection\nto HBeAg-positive CHB",
                              "sag_loss_slope" = "Rate of HBsAg loss")) +
  labs(fill="Sign", title = "10 parameters with highest PRCC") +
  theme_classic()

# Not in monit ICER: pr_ir_ic, transmission-related parameters
# Additionally in monit ICER: mu_cc, pr_cc_dcc, pr_it_ir, sag_loss_slope

# Result:
# in both cases, parameters informing progression to treatment eligibility (pr_ic_enchb)
# is the most influential for ICER (negatively correlated)
# That being said, for the no monitoring programme almost all calibrated values of pr_ic_enchb
# lead to ICER < 518:
plot(x=params_mat_accepted_kmeans$pr_ic_enchb, y = icer_monit_0$icer,
     ylim = c(0,650))
abline(h=391)
abline(h=518)

plot(x=params_mat_accepted_kmeans$pr_ic_enchb, y = icer_monit_sim7$icer,
     ylim = c(0,1000))
abline(h=391)
abline(h=518)
# ICERs > WTP only occur if pr_ic_enchb < 0.004

# Other influential parameters:
# pr_enchb_cc_female, pr_ir_cc_female, pr_ir_enchb, mtct_prob_e, pr_ir_ic (no monitoring)
# pr_enchb_cc_female, cirrhosis_male_cofactor, pr_ir_cc_female, mu_cc, pr_ir_enchb

# For MTCT prob e the correlation is not obvious (probably depends on other parameters):
plot(x=params_mat_accepted_kmeans$mtct_prob_e, y = icer_monit_0$icer,
     ylim = c(0,650))
abline(h=391)
abline(h=518)

# It might be more useful to look qualitatively at ICER less or over the WTP!
# Also need to check monotonic relationship for all parameters

plot(x=params_mat_accepted_kmeans$pr_ic_enchb, y = icer_monit_sim7$icer,
     ylim = c(0,max(icer_monit_sim7$icer)))
abline(h=391)
abline(h=518)

# Associations with the monitoring ICER but not the no-monitoring CER:
plot(x=params_mat_accepted_kmeans$mu_cc, y = icer_monit_sim7$icer,
     ylim = c(0,max(icer_monit_sim7$icer)))
abline(h=391)
abline(h=518)

plot(x=params_mat_accepted_kmeans$pr_cc_dcc, y = icer_monit_sim7$icer,
     ylim = c(0,max(icer_monit_sim7$icer)))
abline(h=391)
abline(h=518)

plot(x=params_mat_accepted_kmeans$pr_it_ir, y = icer_monit_sim7$icer,
     ylim = c(0,max(icer_monit_sim7$icer)))
abline(h=391)
abline(h=518)
# monitoring programme is less likely to be cost-effective if
# rate of HBsAg loss is low

plot(x=params_mat_accepted_kmeans$pr_it_ir, y = icer_monit_sim7$icer,
     ylim = c(0,max(icer_monit_sim7$icer)))
abline(h=391)
abline(h=518)

# Looks like the monitoring programme is less likely to be cost-effective if:
# rate of decompensation is lower
# rate of HBsAg loss is low
# progression from IT to IR is lower

# It's notable that treatment ICER is mainly related to progressions to cirrhosis and not HCC,
# even though the HCC progression parameters are very uncertain and drive the whole age-specific
# result.

# 2-way sensitivity analysis with most influential parameters
two_way_df <- cbind(params_mat_accepted_kmeans,
                    data.frame(icer_monit_0=icer_monit_0$icer,
                               icer_monit_sim7=icer_monit_sim7$icer))

# Example with 518 as the switch where the colour diverges
ggplot(two_way_df,aes(x=pr_ic_enchb, y=pr_enchb_cc_female, z = icer_monit_sim7)) +
  stat_summary_2d() +
  # geom_point(shape = 1, col = 'white') +
  scale_fill_continuous_diverging(mid=518, palette = "Blue-Red 3",
                                  l1 = 30, l2 = 100, p1 = .9, p2 = 1.2) +
  theme_classic()
# ICER only >518 for low values of both parameters pr_ic_enchb and pr_enchb_cc_female
# Though there is not example of pr_enchb_cc_female being high with
# pr_ic_enchb being small
# Same pattern for 391 cut-off with monit_0

# Can't really do a proper 2-way sensitivity analysis with these sets because I don't
# have continuous vector for all combinations

ggplot(two_way_df,aes(x=pr_ic_enchb, y=pr_ir_cc_female, z = icer_monit_sim7)) +
  stat_summary_2d() +
  # geom_point(shape = 1, col = 'white') +
  scale_fill_continuous_diverging(mid=518, palette = "Blue-Red 3",
                                  l1 = 30, l2 = 100, p1 = .9, p2 = 1.2) +
  theme_classic()
# Less clear relationship for these 2 parameters or most others I've tried

ggplot(two_way_df,aes(x=pr_ic_enchb, y=cirrhosis_male_cofactor, z = icer_monit_sim7)) +
  stat_summary_2d() +
  # geom_point(shape = 1, col = 'white') +
  scale_fill_continuous_diverging(mid=518, palette = "Blue-Red 3",
                                  l1 = 30, l2 = 100, p1 = .9, p2 = 1.2) +
  theme_classic()
# Again ICER >518 for combination of low pr_ic_enchb and low cirrhosis_male_cofactor
# but again also no combination of high cofactor with low pr_ic_enchb

# So maybe the overall trend is actually just low pr_ic_enchb that drives the ICER,
# like for mtct_prob_e:
ggplot(two_way_df,aes(x=pr_ic_enchb, y=mtct_prob_e, z = icer_monit_sim7)) +
  stat_summary_2d() +
  # geom_point(shape = 1, col = 'white') +
  scale_fill_continuous_diverging(mid=518, palette = "Blue-Red 3",
                                  l1 = 30, l2 = 100, p1 = .9, p2 = 1.2) +
  theme_classic()


# PRCC for the incremental impact of annual vs 5-yearly monitoring ----
dalys_averted_cohort <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out5_it,
                                 scenario_objects = list(out6_it),
                                 outcome_to_avert = "cohort_dalys",
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "no monitoring")

outcome <- dalys_averted_cohort$value[dalys_averted_cohort$type=="proportion_averted"]

# Load parameter sets chosen based on kmeans clustering
load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans

# Run PRCC
prcc_diff_outcome <- pcc(params_mat_accepted_kmeans,
                         outcome,
                         rank = TRUE, nboot = 100)
plot(prcc_diff_outcome)
abline(h=0)

# Subset only significant parameters (95% CI not including 0):
diff_outcome_significant_parms <- data.frame(parm = rownames(prcc_diff_outcome$PRCC)[
  which(sign(prcc_diff_outcome$PRCC$`max. c.i.`)==sign(prcc_diff_outcome$PRCC$`min. c.i.`))],
  prcc_mean = prcc_diff_outcome$PRCC$original[
    which(sign(prcc_diff_outcome$PRCC$`max. c.i.`)==sign(prcc_diff_outcome$PRCC$`min. c.i.`))],
  prcc_ci_lower = prcc_diff_outcome$PRCC$`min. c.i.`[
    which(sign(prcc_diff_outcome$PRCC$`max. c.i.`)==sign(prcc_diff_outcome$PRCC$`min. c.i.`))],
  prcc_ci_upper = prcc_diff_outcome$PRCC$`max. c.i.`[
    which(sign(prcc_diff_outcome$PRCC$`max. c.i.`)==sign(prcc_diff_outcome$PRCC$`min. c.i.`))])
diff_outcome_significant_parms <- arrange(diff_outcome_significant_parms, -abs(prcc_mean))

# Want to know: why are returns with increasing monitoring frequency diminishing so quickly?

# This suggests that the impact of annual monitoring compared to 5-yearly would be
# higher if progression on the cirrhosis pathway was higher (pr_enchb_cc_female and pr_ic_enchb),
# which is stopped by treatment, and if progression on the HCC pathway was lower (cancer_male_cofactor,
# cancer_prog_coefficient_female)
# Good reminder that the impact of treatment does not just depend on how quickly
# people progress to disease, but also how treatment works to prevent that!
# Interestingly, the HCC parameters do not feature if looking at difference in proportion
# averted in the POPULATION instead of the cohort. But here we are really interested
# in the cohort effects.



# Create dataframes with cost (load functions from other script) ----
annual_discounting_rate <- 0.03

# CHANGE HERE FOR OBJECTS OF INTEREST
#object_list <- list(out3_it, monit_out7)
# Qualitative variations:
#object_list <- list(out3_s1, monit_out7_s1)
#object_list <- list(out3_s2, monit_out7_s2)
#object_list <- list(out3_s3, monit_out7_s3)
# One-way:
# Treatment parms:
#object_list <- list(out3_tmu_dcc_lower, monit_out7_tmu_dcc_lower)
#object_list <- list(out3_tmu_dcc_upper, monit_out7_tmu_dcc_upper)
#object_list <- list(out3_thccr_chb_lower, monit_out7_thccr_chb_lower)
#object_list <- list(out3_thccr_chb_upper, monit_out7_thccr_chb_upper)
#object_list <- list(out3_thccr_it_lower, monit_out7_thccr_it_lower)
#object_list <- list(out3_thccr_it_upper, monit_out7_thccr_it_upper)
#object_list <- list(out3_thccr_cc_lower, monit_out7_thccr_cc_lower)
#object_list <- list(out3_thccr_cc_upper, monit_out7_thccr_cc_upper)
#object_list <- list(out3_thccr_dcc_lower, monit_out7_thccr_dcc_lower)
#object_list <- list(out3_thccr_dcc_upper, monit_out7_thccr_dcc_upper)
# Coverage parms:
#object_list <- list(out3_assessment_lower, monit_out7_assessment_lower)
#object_list <- list(out3_assessment_upper, monit_out7_assessment_upper)
#object_list <- list(out3_treatment_lower, monit_out7_treatment_lower)
#object_list <- list(out3_it, monit_out7_monitoring_lower)
#object_list <- list(out3_it, monit_out7_monitoring_upper)
object_list <- list(out3_screening_lower, monit_out7_screening_lower)
#object_list <- list(out3_screening_upper, monit_out7_screening_upper)

# Extract interactions and person-years on treatment, and HBV-related deaths
# ad DALYs averted, in a for loop
age_interactions <- list()
age_interactions_py_on_treatment <- list()
age_hbv_deaths_averted <- list()
age_dalys_averted <- list()

for (i in 1:length(object_list)) {
  age_interactions[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          assemble_discounted_interactions_for_monitoring_frequencies(object_list[[i]],
                                                                      no_monitoring_object = out3_it))
  age_interactions_py_on_treatment[[i]] <-
    data.frame(scenario = object_list[[i]]$cohort_age_at_death$scenario,
               discount_outcome_2020_to_2100(scenario_object=object_list[[i]],
                                             object_to_subtract=NULL,
                                             outcome="py_on_treatment",
                                             yearly_discount_rate=annual_discounting_rate))
  age_hbv_deaths_averted[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          discount_outcome_2020_to_2100(scenario_object=out2,
                                        object_to_subtract=object_list[[i]],
                                        outcome="cum_hbv_deaths",
                                        yearly_discount_rate=annual_discounting_rate))
  age_dalys_averted[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          discount_outcome_2020_to_2100(scenario_object=out2,
                                        object_to_subtract=object_list[[i]],
                                        outcome="dalys",
                                        yearly_discount_rate=annual_discounting_rate))
}
age_interactions <- do.call("rbind", age_interactions)
age_interactions_py_on_treatment <- do.call("rbind", age_interactions_py_on_treatment)
age_hbv_deaths_averted <- do.call("rbind", age_hbv_deaths_averted)
age_dalys_averted <- do.call("rbind", age_dalys_averted)

age_interactions_py_on_treatment$sim <- gsub("[^0-9]", "", age_interactions_py_on_treatment$sim)
age_hbv_deaths_averted$sim <- gsub("[^0-9]", "", age_hbv_deaths_averted$sim)
colnames(age_hbv_deaths_averted)[colnames(age_hbv_deaths_averted) == "cum_hbv_deaths"] <-
  "value"

age_dalys_averted$sim <- gsub("[^0-9]", "", age_dalys_averted$sim)
colnames(age_dalys_averted)[colnames(age_dalys_averted) == "dalys"] <-
  "value"

age_df <- create_incremental_plot_df(interactions_df=age_interactions,
                                     py_on_treatment_df=age_interactions_py_on_treatment,
                                     deaths_averted_df=age_hbv_deaths_averted,
                                     ly_saved_df = age_dalys_averted, # replace LY by DALYs
                                     hbsag_test_cost = 8.3,
                                     clinical_assessment_cost = 33, #84.4,
                                     monitoring_assessment_cost = 25.5, #40.1,
                                     treatment_py_cost = 66.5,#60,
                                     #scenario_labels_obj = scenario_labels,
                                     ref_label = "No treatment")
colnames(age_df)[colnames(age_df)=="ly_saved"] <- "dalys_averted"

age_df %>% group_by(scenario) %>%
  summarise(cost=median(total_cost),
            dalys = median(dalys_averted))

icer_list <- list()

for(i in 1:183) {
  print(i)
  icer_list[[i]] <- age_df[which(age_df$sim==
                                    unique(age_df$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result



# Tornado plot for treatment effect ----
# ICERS are:
# Default simulations: out3 = 264 (163-493), monit_out7 = 321 (155-767)
# All thccr parms to 0.4: out3 = 286 (170-547), monit_out7 = 343 (163-894)
# Cirrhosis thccr parms to 0.8: out3 = 301 (178-587), monit_out7 = 345 (166-806)
# CHB thccr parms to 0.8: out3 = 290 (170-635), monit_out7 = 368 (164-1377)
# On average, these strategies all remain cost-effective despite increases in the ICER
# One-way sensitivity:
# tmu_dcc_lower:  out3 = 264 (162-491), monit_out7 = 321 (154-766)
# tmu_dcc_upper:  out3 = 266 (163-494), monit_out7 = 321 (155-768)
# thccr_chb_lower:  out3 = 260 (160-471), monit_out7 = 308 (152-705)
# thccr_chb_upper:  out3 = 286 (168-609), monit_out7 = 361 (163-1084)
# thccr_it_lower: out3 = 264 (163-490), monit_out7 = 320(154-761)
# thccr_it_upper:  out3 = 269 (164-509), monit_out7 = 329 (156-914)
# thccr_cc_lower: out3 = 256 (158-474), monit_out7 = 316(152-757)
# thccr_cc_upper:  out3 = 300 (178-585), monit_out7 = 345 (166-805)
# thccr_dcc_lower: out3 = 264 (163-492), monit_out7 = 321 (155-767)
# thccr_dcc_upper: out3 = 264 (163-493), monit_out7 = 321 (155-768)

# For tornado plot, take median ICER for default, and lower and upper value for each parameter.
# Have separate columns for lower and upper, one of which is negative from median.
# Add coord_flip and 2 separete layers for lower and upper

# Median ICERS:
tornado_df_out3 <- data.frame(outcome = "monit_0",
                              parm = c("thccr_it","thccr_chb","thccr_cc","thccr_dcc","tmu_dcc"),
                              lower_icer = c(264,260,256,264,264),
                              default_icer = c(264,264,264,264,264),
                              upper_icer = c(269,286,300,264,266))
tornado_df_out3$parm <- factor(tornado_df_out3$parm,
                               levels = c("thccr_it","thccr_chb","thccr_cc","thccr_dcc","tmu_dcc"))

tornado_df_monit_sim7 <- data.frame(outcome = "monit_sim7",
                                    parm = c("thccr_it","thccr_chb","thccr_cc","thccr_dcc","tmu_dcc"),
                                    lower_icer = c(320,308,316,321,321),
                                    default_icer = c(321,321,321,321,321),
                                    upper_icer = c(329,361,345,321,321))
tornado_df_monit_sim7$parm <- factor(tornado_df_monit_sim7$parm,
                               levels = c("thccr_it","thccr_chb","thccr_cc","thccr_dcc","tmu_dcc"))

lower_wtp <- 404

ggplot(tornado_df_out3) +
  geom_col(aes(x=parm, y = upper_icer-default_icer), fill = "blue") +
  geom_col(aes(x=parm, y = lower_icer-default_icer), fill = "red") +
  geom_text(aes(label = c(0.8,0.8,0.8,0.8,0.08), x = parm, y = upper_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = -0.5) +
  geom_text(aes(label = c(0.1,0.1,0.1,0.1,0.35), x = parm, y = lower_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = 1.5) +
  geom_hline(yintercept=lower_wtp-264, lty="dashed") +
  geom_hline(yintercept=0) +
  scale_x_discrete("Treatment effect parameter",
                   labels = c("thccr_it" = "RR for progression to HCC\nfrom HBeAg+ infection (0.19)",
                              "thccr_chb" = "RR for progression to HCC\nfrom CHB (0.27)",
                              "thccr_cc" = "RR for progression to HCC\nfrom CC (0.23)",
                              "thccr_dcc" = "RR for progression to HCC\nfrom DCC (0.13)",
                              "tmu_dcc" = "Mortality rate from\ntreated DCC (0.18)"),
                   limits = rev)+
  scale_y_continuous(breaks = c(200-264,0,lower_wtp-264),
                     labels = c(200,264,lower_wtp),
                     limits=c(200-264,(lower_wtp-264+10))) +
  ylab("Median incremental cost per averted DALY") +
  labs(title = "2020 screening and treatment without monitoring\nvs. status quo of no treatment") +
  theme_classic() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12)) +
  coord_flip()


ggplot(tornado_df_monit_sim7) +
  geom_col(aes(x=parm, y = upper_icer-default_icer), fill = "blue") +
  geom_col(aes(x=parm, y = lower_icer-default_icer), fill = "red") +
  geom_text(aes(label = c(0.8,0.8,0.8,0.8,0.08), x = parm, y = upper_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = -0.5) +
  geom_text(aes(label = c(0.1,0.1,0.1,0.1,0.35), x = parm, y = lower_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = 1.5) +
  geom_hline(yintercept=lower_wtp-321, lty="dashed") +
  geom_hline(yintercept=0) +
  scale_x_discrete("Treatment effect parameter",
                   labels = c("thccr_it" = "RR for progression to HCC\nfrom HBeAg+ infection (0.19)",
                              "thccr_chb" = "RR for progression to HCC\nfrom CHB (0.27)",
                              "thccr_cc" = "RR for progression to HCC\nfrom CC (0.23)",
                              "thccr_dcc" = "RR for progression to HCC\nfrom DCC (0.13)",
                              "tmu_dcc" = "Mortality rate from\ntreated DCC (0.18)"),
                   limits = rev)+
  scale_y_continuous(breaks = c(200-321,0,lower_wtp-321),
                     labels = c(200,321,lower_wtp),
                     limits=c(200-321,(lower_wtp-321+10))) +
  ylab("Median incremental cost per averted DALY") +
  labs(title = "2020 screening and treatment with optimal* monitoring\nscenario vs. treatment programme without monitoring") +
  theme_classic() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12)) +
  coord_flip()

# One-way sensitivity analysis of cost on ICERs (bar charts) ----

cost_sensitivity_icer$sensitivity_scenario <-
  factor(cost_sensitivity_icer$sensitivity_scenario,
         levels = c("screening_lower", "screening_higher",
                    "assessment_lower", "assessment_higher",
                    "treatment_lower", "treatment_higher",
                    "monitoring_lower", "monitoring_higher",
                    "assessment_monitoring_lower","assessment_monitoring_higher"))

cost_sensitivity_icer$scenario <-
  factor(cost_sensitivity_icer$scenario,
         levels = c("screen_2020_monit_0", "screen_2020_monit_sim6",  "screen_2020_monit_sim7",
                    "screen_2020_monit_sim2c", "screen_2020_monit_5","screen_2020_monit_4",
                    "screen_2020_monit_3","screen_2020_monit_2","screen_2020_monit_1"))

facet_labels_costs <- c("Screening = $4", "Screening = $17",
                        "Assessment = $16.5", "Assessment = $200",
                        "Treatment = $51", "Treatment = $155",
                        "Monitoring = $12.75", "Monitoring = $51",
                        "Assessment = $16.5\nMonitoring = $12.75",
                        "Assessment = $200\nMonitoring = $155")
names(facet_labels_costs) <- c("screening_lower", "screening_higher",
                               "assessment_lower", "assessment_higher",
                               "treatment_lower", "treatment_higher",
                               "monitoring_lower", "monitoring_higher",
                               "assessment_monitoring_lower",
                               "assessment_monitoring_higher")

# Every 2 and every 1 year frequencies are not shown to keep scales manageable
# Monitoring is in all ages unless otherwise indicated
# Changes in screening and assessment cost does not change ICER between incremental monitoring
# strategies.

# Only showing non-dominated strategies - Blank spaces are if dominated
cost_bars1 <- ggplot(subset(cost_sensitivity_icer, scenario != "screen_2020_monit_1" &
                scenario != "screen_2020_monit_2" &
                sensitivity_scenario != "assessment_monitoring_lower" &
                sensitivity_scenario != "assessment_monitoring_higher")) +
  geom_col(aes(x=scenario, y = icer_median, fill = scenario),
           width=0.95, col="black") +
  geom_errorbar(aes(x=scenario, ymin=icer_lower,
                    ymax=icer_upper), width=0.15) +
  geom_hline(yintercept=404, linetype="dashed", colour = "grey40") +
  geom_hline(yintercept=537, linetype="dashed", colour = "grey40") +
  geom_hline(yintercept=778, linetype="dashed", colour = "grey40") +
  scale_fill_manual(labels=c("screen_2020_monit_0" = "No monitoring",
                             "screen_2020_monit_sim6" = "5-yearly (15-30)",
                             "screen_2020_monit_sim7"= "5-yearly (15-45)",
                             "screen_2020_monit_sim2c"="4-yearly (15-45)",
                             "screen_2020_monit_5"="5-yearly",
                             "screen_2020_monit_4"="4-yearly",
                             "screen_2020_monit_3"="3-yearly"),
                    values=rev(brewer.pal(7, "Blues"))) +
  scale_x_discrete(labels=c("screen_2020_monit_0" = "No\nmonitoring",
                             "screen_2020_monit_sim6" = "5 (<30)",
                             "screen_2020_monit_sim7"= "5 (<45)",
                             "screen_2020_monit_sim2c"="4 (<45)",
                             "screen_2020_monit_5"="5",
                             "screen_2020_monit_4"="4",
                             "screen_2020_monit_3"="3")) +
  ylab("ICER") +
  xlab("Monitoring frequency (age group)") +
  facet_wrap(~sensitivity_scenario, scales="free_y", ncol=2,
             labeller = labeller(sensitivity_scenario = facet_labels_costs)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16))


# Show joint variation in assessment and monitoring separately (key message)
cost_bars2 <- ggplot(subset(cost_sensitivity_icer, (scenario != "screen_2020_monit_1" &
                scenario != "screen_2020_monit_2" &
                scenario != "screen_2020_monit_3") &
                (sensitivity_scenario == "assessment_monitoring_lower" |
                sensitivity_scenario == "assessment_monitoring_higher"))) +
  geom_col(aes(x=scenario, y = icer_median, fill = scenario),
           width=0.95, col="black") +
  geom_errorbar(aes(x=scenario, ymin=icer_lower,
                    ymax=icer_upper), width=0.15) +
  geom_hline(yintercept=404, linetype="dashed", colour = "grey40") +
  geom_hline(yintercept=537, linetype="dashed", colour = "grey40") +
  geom_hline(yintercept=778, linetype="dashed", colour = "grey40") +
  scale_fill_manual(labels=c("screen_2020_monit_0" = "No monitoring",
                             "screen_2020_monit_sim6" = "5-yearly (15-30)",
                             "screen_2020_monit_sim7"= "5-yearly (15-45)",
                             "screen_2020_monit_sim2c"="4-yearly (15-45)",
                             "screen_2020_monit_5"="5-yearly",
                             "screen_2020_monit_4"="4-yearly"),
                    values=rev(brewer.pal(6, "Blues"))) +
  scale_x_discrete(labels=c("screen_2020_monit_0" = "No\nmonitoring",
                            "screen_2020_monit_sim6" = "5 (<30)",
                            "screen_2020_monit_sim7"= "5 (<45)",
                            "screen_2020_monit_sim2c"="5 (<45)",
                            "screen_2020_monit_5"="5",
                            "screen_2020_monit_4"="4")) +
  ylab("ICER") +
  xlab("Monitoring frequency (age group)") +
  facet_wrap(~sensitivity_scenario, scales="free_y", ncol=2,
             labeller = labeller(sensitivity_scenario = facet_labels_costs)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16))

cost_bars1a <- arrangeGrob(cost_bars1, top = textGrob("A", x = unit(0.01, "npc"),
                                                      y   = unit(1, "npc"), just=c("left","top"),
                                                      gp=gpar(col="black", fontsize=22)))
cost_bars2b <- arrangeGrob(cost_bars2, top = textGrob("B", x = unit(0.01, "npc"),
                                                      y   = unit(1, "npc"), just=c("left","top"),
                                                      gp=gpar(col="black", fontsize=22)))

#png(file = "sensitivity_costs_barchart.png", width=350, height=300, units = "mm", res=300, pointsize = 0.99)
grid.arrange(cost_bars1a, cost_bars2b, nrow=2, heights=c(3,1))
#dev.off()

# Tornado plot for costs and coverage parameters ----

# ICERS for coverage:
# screening_lower: out3 = 271 (166-501), monit_out7 = 321 (155-767)
# screening_upper: out3 =  264 (162-492), monit_out7 =  321 (155-767)
# assessment_lower: out3 = 338 (203-642), monit_out7 =  321 (155-767)
# Monitoring probably gets more favourable when assessment is lower
# assessment_upper: out3 = 242 (149-449), monit_out7 =  321 (155-767)
# treatment_lower: out3 =  409 (243-784), monit_out7 =  409 (197-920)
# monitoring_lower: monit_out7 = 289 (143-679)
# monitoring_upper: monit_out7 = 342 (163-826)
# ICER gets more unfavourable with more monitoring because it just increases frequency

tornado_cov_out3 <- data.frame(outcome = "monit_0",
                                parm = c("screening","assessment","treatment","monitoring"),
                                lower_icer = c(271,338,409,264),
                                default_icer = c(264,264,264,264),
                                upper_icer = c(264,242,264,264))
tornado_cov_out3$parm <- factor(tornado_cov_out3$parm,
                                 levels = c("screening","assessment","treatment","monitoring"))

tornado_cov_monit_sim7 <- data.frame(outcome = "monit_sim7",
                                      parm = c("screening","assessment","treatment","monitoring"),
                                      lower_icer = c(321,321,409,289),
                                      default_icer = c(321,321,321,321),
                                      upper_icer = c(321,321,321,342))
tornado_cov_monit_sim7$parm <- factor(tornado_cov_monit_sim7$parm,
                                       levels = c("screening","assessment","treatment","monitoring"))



lower_wtp <- 404

# Note screening bar refers to reduction to 50% whereas going up to 100% did not affect median
ggplot(tornado_cov_out3) +
  geom_col(aes(x=parm, y = upper_icer-default_icer), fill = "blue") +
  geom_col(aes(x=parm, y = lower_icer-default_icer), fill = "red") +
  geom_text(aes(label = c(100,100,100,100), x = parm, y = upper_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = -0.5) +
  geom_text(aes(label = c(50,50,50,50), x = parm, y = lower_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = 1.5) +
  geom_hline(yintercept=lower_wtp-264, lty="dashed") +
  geom_hline(yintercept=0) +
  scale_x_discrete("Coverage (%)",
                   labels = c("screening" = "Screening (90)",
                              "assessment" = "Linkage to care (80)",
                              "treatment" = "Treatment initiation (100)",
                              "monitoring" = "Monitoring (80)"),
                   limits = rev)+
  scale_y_continuous(breaks = c(200-264,0,lower_wtp-264),
                     labels = c(200,264,lower_wtp),
                     limits=c(200-264,(420-264+10))) +
  ylab("Median incremental cost per averted DALY") +
  labs(title = "2020 screening and treatment without monitoring\nvs. status quo of no treatment") +
  theme_classic() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12)) +
  coord_flip()

# Note that assessment probability does not affect the CER of monit_sim7 compared to out3,
# but it makes monitoring more favourable compared to no monitoring.
ggplot(tornado_cov_monit_sim7) +
  geom_col(aes(x=parm, y = upper_icer-default_icer), fill = "blue") +
  geom_col(aes(x=parm, y = lower_icer-default_icer), fill = "red") +
  geom_text(aes(label = c(100,100,100,100), x = parm, y = upper_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = -0.5) +
  geom_text(aes(label = c(50,50,50,50), x = parm, y = lower_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = 1.5) +
  geom_hline(yintercept=lower_wtp-321, lty="dashed") +
  geom_hline(yintercept=0) +
  scale_x_discrete("Coverage (%)",
                   labels = c("screening" = "Screening (90)",
                              "assessment" = "Linkage to care (80)",
                              "treatment" = "Treatment initiation (100)",
                              "monitoring" = "Monitoring (80)"),
                   limits = rev)+
  scale_y_continuous(breaks = c(200-321,0,lower_wtp-321),
                     labels = c(200,321,lower_wtp),
                     limits=c(200-321,(420-321+10))) +
  ylab("Median incremental cost per averted DALY") +
  labs(title = "2020 screening and treatment with optimal* monitoring\nscenario vs. treatment programme without monitoring") +
  theme_classic() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12)) +
  coord_flip()


## COSTS

# Needed to run median ICERS separately here for:
# no monitoring vs no treatment and
# 5-yearly in <45 vs no monitoring

# Median ICERS:
subset(cost_sensitivity_icer, scenario == "screen_2020_monit_0")
# Treatment lower median ICER = 238
# Monitoring lower median ICER =264
# Screening higher median ICER = 379
# Assessment higher median ICER = 408

tornado_cost_out3 <- data.frame(outcome = "monit_0",
                              parm = c("screening","assessment","treatment","monitoring"),
                              lower_icer = c(209,252,238,264),
                              default_icer = c(264,264,264,264),
                              upper_icer = c(379,408,432,264))
tornado_cost_out3$parm <- factor(tornado_cost_out3$parm,
                               levels = c("screening","assessment","treatment","monitoring"))

# sim7 vs. no monitoring median ICERs:
# Screening lower = 321
# Assessment lower = 321
# Treatment lower =280
# Monitoring lower =243
# Screening upper =321
# Assessment upper =321
# Treatment upper =540
# Monitoring upper =476
tornado_cost_monit_sim7 <- data.frame(outcome = "monit_sim7",
                                parm = c("screening","assessment","treatment","monitoring"),
                                lower_icer = c(321,321,280,243),
                                default_icer = c(321,321,321,321),
                                upper_icer = c(321,321,540,476))
tornado_cost_monit_sim7$parm <- factor(tornado_cost_monit_sim7$parm,
                                 levels = c("screening","assessment","treatment","monitoring"))


lower_wtp <- 404

ggplot(tornado_cost_out3) +
  geom_col(aes(x=parm, y = upper_icer-default_icer), fill = "blue") +
  geom_col(aes(x=parm, y = lower_icer-default_icer), fill = "red") +
  geom_text(aes(label = c(17,200,155,51), x = parm, y = upper_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = -0.5) +
  geom_text(aes(label = c(4,16.5,51,12.75), x = parm, y = lower_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = 1.5 ) +
  geom_hline(yintercept=lower_wtp-264, lty="dashed") +
  geom_hline(yintercept=0) +
  scale_x_discrete("Cost per person (2020 US$)",
                   labels = c("screening" = "Screening (8.3)",
                              "assessment" = "Initial assessment (33)",
                              "treatment" = "Treatment\nper year (66.5)",
                              "monitoring" = "Monitoring assessment\nper visit (25.5)"),
                   limits = rev)+
  scale_y_continuous(breaks = c(200-264,0,lower_wtp-264,440-264),
                     labels = c(200,264,lower_wtp,440),
                     limits=c(200-264,(440-264+10))) +
  ylab("Median incremental cost per averted DALY") +
  labs(title = "2020 screening and treatment without monitoring\nvs. status quo of no treatment") +
  theme_classic() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12)) +
  coord_flip()


ggplot(tornado_cost_monit_sim7) +
  geom_col(aes(x=parm, y = upper_icer-default_icer), fill = "blue") +
  geom_col(aes(x=parm, y = lower_icer-default_icer), fill = "red") +
  geom_text(aes(label = c(17,200,155,51), x = parm, y = upper_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = -0.5) +
  geom_text(aes(label = c(4,16.5,51,12.75), x = parm, y = lower_icer-default_icer),
            position = position_dodge(width = 0.8), hjust = 1.5) +
  geom_hline(yintercept=lower_wtp-321, lty="dashed") +
  geom_hline(yintercept=0) +
  scale_x_discrete("Cost per person (2020 US$)",
                   labels = c("screening" = "Screening (8.3)",
                              "assessment" = "Initial assessment (33)",
                              "treatment" = "Treatment\nper year (66.5)",
                              "monitoring" = "Monitoring assessment\nper visit (25.5)"),
                   limits = rev)+
  scale_y_continuous(breaks = c(200-321,0,lower_wtp-321,580-321),
                     labels = c(200,321,lower_wtp,580),
                     limits=c(200-321,(580-321+10))) +
  ylab("Median incremental cost per averted DALY") +
  labs(title = "2020 screening and treatment with optimal* monitoring\nscenario vs. treatment programme without monitoring") +
  theme_classic() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12)) +
  coord_flip()



# Complete tornado plot (treatment effect, coverage, cost) ----

# No monitoring compared to no treatment
# Median ICERS

# Need to run previous sections!
levels(tornado_cost_out3$parm) <- list("screening_cost" = "screening",
                                       "assessment_cost" = "assessment",
                                       "treatment_cost" = "treatment",
                                       "monitoring_cost" = "monitoring")
levels(tornado_cost_monit_sim7$parm) <- list("screening_cost" = "screening",
                                       "assessment_cost" = "assessment",
                                       "treatment_cost" = "treatment",
                                       "monitoring_cost" = "monitoring")

tornado_out3 <- rbind(tornado_df_out3, tornado_cov_out3, tornado_cost_out3)
tornado_monit_sim7 <- rbind(tornado_df_monit_sim7, tornado_cov_monit_sim7, tornado_cost_monit_sim7)

lower_wtp <- 404
upper_wtp <- 537
# Ranges:
# thccr_it, thccr_chb, thccr_cc, thccr_dcc: 0.8, 0.1
# tmu_dcc: 0.08, 0.35
# screening, assessment, treatment, monitoring: 100, 50
# screening_cost: 17, 4
# assessment_cost: 200, 16.5
# treatment_cost 155, 51
# monitoring_cost: 51, 12.75

tornado_p1 <- ggplot(subset(tornado_out3, parm != "monitoring" & parm != "monitoring_cost")) +
  geom_col(aes(x=reorder(parm, -pmax(upper_icer, lower_icer)),
               y = upper_icer-default_icer), fill = "#377EB8", colour = "black", width = 0.6) +
  geom_col(aes(x=reorder(parm, -pmax(upper_icer, lower_icer)),
               y = lower_icer-default_icer), fill = "#E41A1C", colour = "black", width = 0.6) +
  #geom_text(aes(label = upper_parm, x = parm, y = upper_icer-default_icer),
 #           position = position_dodge(width = 0.8), hjust = -0.5) +
 # geom_text(aes(label = lower_parm, x = parm, y = lower_icer-default_icer),
  #          position = position_dodge(width = 0.8), hjust = 1.5) +
  geom_hline(yintercept=lower_wtp-264, lty="dashed") +
  geom_hline(yintercept=upper_wtp-264, lty="dashed") +
  geom_hline(yintercept=0) +
  scale_x_discrete("",
                   labels = c("screening_cost" = "Screening cost ($8.3, 4-17)",
                              "assessment_cost" = "Initial assessment cost ($33, 16.5-200)",
                              "treatment_cost" = "Treatment cost per year ($66.5, 51-155)",
                              "monitoring_cost" = "Monitoring cost per visit ($25.5, 12.75-51)",
                              "thccr_it" = "Rate ratio for progression to HCC from\nHBeAg+ infection (0.19, 0.1-0.8)",
                              "thccr_chb" = "Rate ratio for progression to HCC\nfrom CHB (0.27, 0.1-0.8)",
                              "thccr_cc" = "Rate ratio for progression to HCC\nfrom CC (0.23, 0.1-0.8)",
                              "thccr_dcc" = "Rate ratio for progression to HCC\nfrom DCC (0.13, 0.1-0.8)",
                              "tmu_dcc" = "Mortality rate from treated\nDCC (0.18, 0.08-0.35)",
                              "screening" = "Screening coverage (90%, 50-100)",
                              "assessment" = "Linkage to care (80%, 50-100)",
                              "treatment" = "Treatment initiation (100%, 50-100)",
                              "monitoring" = "Monitoring uptake (80%, 50-100)"),
                   limits = rev)+
  scale_y_continuous(breaks = c(200-264,0,lower_wtp-264,upper_wtp-264),
                     labels = c(200,264,lower_wtp,upper_wtp),
                     limits=c(200-264,(550-264))) +
  ylab("Median incremental cost (US$) per DALY averted") +
  theme_classic() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12)) +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
# Costs are per person, rate is in per person-year
# Monitoring parameters not shown because not relevant for scenario
# Blue represents an increase in the parameter value, red a decrease
# Mention that all median ICERs were well below higher estimates of the CE threshold (537)

# Monitoring compared to no monitoring
tornado_p2 <- ggplot(subset(tornado_monit_sim7, !(parm %in% c("screening_cost", "screening",
                                                    "assessment_cost", "assessment")))) +
  geom_col(aes(x=reorder(parm, -pmax(upper_icer, lower_icer)),
               y = upper_icer-default_icer), fill = "#377EB8", colour = "black", width = 0.6) +
  geom_col(aes(x=reorder(parm, -pmax(upper_icer, lower_icer)),
               y = lower_icer-default_icer), fill = "#E41A1C", colour = "black", width = 0.6) +
  #geom_text(aes(label = upper_parm, x = parm, y = upper_icer-default_icer),
  #           position = position_dodge(width = 0.8), hjust = -0.5) +
  # geom_text(aes(label = lower_parm, x = parm, y = lower_icer-default_icer),
  #          position = position_dodge(width = 0.8), hjust = 1.5) +
  geom_hline(yintercept=lower_wtp-321, lty="dashed") +
  geom_hline(yintercept=upper_wtp-321, lty="dashed") +
  geom_hline(yintercept=0) +
  scale_x_discrete("",
                   labels = c("screening_cost" = "Screening cost ($8.3, 4-17)",
                              "assessment_cost" = "Initial assessment cost ($33, 16.5-200)",
                              "treatment_cost" = "Treatment cost per year ($66.5, 51-155)",
                              "monitoring_cost" = "Monitoring cost per visit ($25.5, 12.75-51)",
                              "thccr_it" = "Rate ratio for progression to HCC from\nHBeAg+ infection (0.19, 0.1-0.8)",
                              "thccr_chb" = "Rate ratio for progression to HCC\nfrom CHB (0.27, 0.1-0.8)",
                              "thccr_cc" = "Rate ratio for progression to HCC\nfrom CC (0.23, 0.1-0.8)",
                              "thccr_dcc" = "Rate ratio for progression to HCC\nfrom DCC (0.13, 0.1-0.8)",
                              "tmu_dcc" = "Mortality rate from treated\nDCC (0.18, 0.08-0.35)",
                              "screening" = "Screening coverage (90%, 50-100)",
                              "assessment" = "Linkage to care (80%, 50-100)",
                              "treatment" = "Treatment initiation (100%, 50-100)",
                              "monitoring" = "Monitoring uptake (80%, 50-100)"),
                   limits = rev)+
  scale_y_continuous(breaks = c(200-321,0,lower_wtp-321,upper_wtp-321),
                     labels = c(200,321,lower_wtp,upper_wtp),
                     limits=c(200-321,(550-321))) +
  ylab("Median incremental cost (US$) per DALY averted") +
  theme_classic() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12)) +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
# Screening and initial assessment parameters not shown cause they don't influence going from
# no monitoring to monitoring

# Combine into 1 plot
library(grid)
tornado_p1a <- arrangeGrob(tornado_p1, top = textGrob("A", x = unit(0.01, "npc"),
                                                  y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=18)))
tornado_p2b <- arrangeGrob(tornado_p2, top = textGrob("B", x = unit(0.01, "npc"),
                                                       y   = unit(1, "npc"), just=c("left","top"),
                                                       gp=gpar(col="black", fontsize=18)))

#tiff(file = "sensitivity_tornado_plot.tiff", width=300, height=250, units = "mm", res=300, pointsize = 0.99)
grid.arrange(tornado_p1a, tornado_p2b, nrow =2)
#dev.off()

# Coverage variation and impact and CER plots ----
dalys_averted_cov <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(monit_out7, monit_out7_screening_lower,
                                                  monit_out7_assessment_lower,
                                                  monit_out7_treatment_lower,
                                                  monit_out7_monitoring_lower,
                                                  out3_it, out3_screening_lower,
                                                  out3_assessment_lower,
                                                  out3_treatment_lower),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
dalys_averted_cov <- dalys_averted_cov %>%
  filter(by_year==2100) %>%
  group_by(scenario, type) %>%
  summarise(median = median(value),
            cri_lower = quantile(value, 0.025),
            cri_upper = quantile(value, 0.975))
dalys_averted_cov$monitoring <- "Yes"
dalys_averted_cov$monitoring[dalys_averted_cov$scenario %in%
                               c("screen_2020_monit_0",
                                 "screen_2020_monit_0_screening_coverage_lower",
                                 "screen_2020_monit_0_link_to_care_prob_lower",
                                 "screen_2020_monit_0_treatment_prob_lower")] <- "No"

dalys_averted_cov$cov <- "Ambitious"
dalys_averted_cov$cov[dalys_averted_cov$scenario %in%
                        c("screen_2020_monit_sim7_screening_coverage_lower",
  "screen_2020_monit_0_screening_coverage_lower")] <- "Reduced screening"
dalys_averted_cov$cov[dalys_averted_cov$scenario %in% c(
  "screen_2020_monit_sim7_link_to_care_prob_lower",
  "screen_2020_monit_0_link_to_care_prob_lower")] <- "Reduced linkage\nto care"
dalys_averted_cov$cov[dalys_averted_cov$scenario %in% c(
  "screen_2020_monit_sim7_treatment_prob_lower",
  "screen_2020_monit_0_treatment_prob_lower")] <- "Reduced treatment initiation"
dalys_averted_cov$cov[dalys_averted_cov$scenario %in% c(
  "screen_2020_monit_sim7_monitoring_prob_lower")] <- "Reduced monitoring uptake"

ggplot(subset(dalys_averted_cov, type=="proportion_averted"))+
  geom_col(aes(x=cov, y = median, fill = monitoring)) +
  geom_errorbar(aes(x=cov, ymin=cri_lower, ymax=cri_upper,
                    group = monitoring), width = 0.15) +
  facet_wrap(~monitoring, scales="free_x") +
  theme_classic()

# Plot of median cost against median DALYs averted
# Costs here are discounted
coverage_cost_effect <- read.csv(file=paste0(out_path, "coverage_sensitivity_results_310321.csv"))

coverage_cost_effect$sensitivity_scenario <- factor(coverage_cost_effect$sensitivity_scenario)
levels(coverage_cost_effect$sensitivity_scenario) <- list("Ambitious uptake" = "default",
                                                       "Reduced screening" = "screening_lower",
                                                       "Reduced linkage\nto care" = "assessment_lower",
                                                       "Reduced treatment\ninitiation" = "treatment_lower",
                                                       "Reduced monitoring" = "monitoring_lower")

#tiff(file = "sensitivity_coverage_cer_plane.tiff", width=300, height=225, units = "mm", res=250, pointsize = 0.99)
ggplot(subset(coverage_cost_effect, !(scenario == "screen_2020_monit_0" &
                                        sensitivity_scenario == "Reduced monitoring"))) +
  geom_point(aes(x=dalys_averted/1000, y = total_cost/1000000,
                 shape=scenario, colour=sensitivity_scenario), size=8) +
  scale_colour_brewer("", palette="Accent") +
  scale_shape_discrete("Scenario",
                     labels = c("screen_2020_monit_0" = "No monitoring",
                                           "screen_2020_monit_sim7" =
                                             "Monitoring 5-yearly\nin <45 year olds")) +
  guides(colour = guide_legend(order = 2),
           shape = guide_legend(order = 1)) +
  theme_classic() +
  ylab("Cost (millions)") +
  xlab("DALYs averted (thousands)") +
  #facet_wrap(~sensitivity_scenario) +
  expand_limits(y = 0, x = 0) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
#dev.off()

# Qualitative sensitivity analysis of treatment impact on infections ----
# Compare out3_it (default) and out3_no_inf, where treated carriers
# do not transmit HBV at all

# How many more DALYs are averted by out3_no_inf?
quantile((out2$dalys[[16]][,-c(1:3)]-out3_no_inf$dalys[[16]][,-c(1:3)])-
(out2$dalys[[16]][,-c(1:3)]-out3_it$dalys[[16]][,-c(1:3)]), c(0.5,0.025,0.975))
# 3372 (487-13120) by 2100

quantile(((out2$dalys[[16]][,-c(1:3)]-out3_no_inf$dalys[[16]][,-c(1:3)])-
           (out2$dalys[[16]][,-c(1:3)]-out3_it$dalys[[16]][,-c(1:3)]))/
           (out2$dalys[[16]][,-c(1:3)]-out3_it$dalys[[16]][,-c(1:3)]), c(0.5,0.025,0.975))
# In % that is 3% (0.6-11%) more than what is achieved in the default

# What % of DALYS are averted in the treated cohort?
# With infectivity:
quantile((out1_it$cohort_dalys[,-1]-out3_it$cohort_dalys[,-1])/
  (out2$dalys[[16]][,-c(1:3)]-out3_it$dalys[[16]][,-c(1:3)]), c(0.5,0.025,0.975))
# 99% (95-104%)

# Without treated infectivity:
quantile((out1_it$cohort_dalys[,-1]-out3_no_inf$cohort_dalys[,-1])/
           (out2$dalys[[16]][,-c(1:3)]-out3_no_inf$dalys[[16]][,-c(1:3)]), c(0.5,0.025,0.975))
# 96% (88-99.9%)

# Infections averted with treated infectivity by 2100
quantile(apply(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time>=2020 &
                                                   out2$timeseries$total_chronic_infections$time<2100,-c(1,2)]-
  out3_it$timeseries$total_chronic_infections[out3_it$timeseries$total_chronic_infections$time>=2020 &
                                                   out3_it$timeseries$total_chronic_infections$time<2100,-c(1,2)],2,sum),
  c(0.5,0.025,0.975))
# 228 (-111-1303)


# Infections averted without treated infectivity by 2100
quantile(apply(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time>=2020 &
                                                          out2$timeseries$total_chronic_infections$time<2100,-c(1,2)]-
                 out3_no_inf$timeseries$total_chronic_infections[out3_no_inf$timeseries$total_chronic_infections$time>=2020 &
                                                               out3_no_inf$timeseries$total_chronic_infections$time<2100,-c(1,2)],2,sum),
         c(0.5,0.025,0.975))
# 973 (258-3217)

# What % of all new infections does this represent?
quantile(apply(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time>=2020 &
                                                          out2$timeseries$total_chronic_infections$time<2100,-c(1,2)]-
                 out3_no_inf$timeseries$total_chronic_infections[out3_no_inf$timeseries$total_chronic_infections$time>=2020 &
                                                                   out3_no_inf$timeseries$total_chronic_infections$time<2100,-c(1,2)],
      2,sum)/
  apply(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time>=2020 &
                                                   out2$timeseries$total_chronic_infections$time<2100,-c(1,2)],
        2,sum), c(0.5,0.025,0.975))
# 7% (5-13%)

# What about by 2030?
# With infectivity
quantile(apply(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time>=2020 &
                                                          out2$timeseries$total_chronic_infections$time<2030,-c(1,2)]-
                 out3_it$timeseries$total_chronic_infections[out3_it$timeseries$total_chronic_infections$time>=2020 &
                                                                   out3_it$timeseries$total_chronic_infections$time<2030,-c(1,2)],
               2,sum)/
           apply(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time>=2020 &
                                                            out2$timeseries$total_chronic_infections$time<2030,-c(1,2)],
                 2,sum), c(0.5,0.025,0.975))
# 4% (0.02-14%)

# No infectivity:
quantile(apply(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time>=2020 &
                                                          out2$timeseries$total_chronic_infections$time<2030,-c(1,2)]-
                 out3_no_inf$timeseries$total_chronic_infections[out3_no_inf$timeseries$total_chronic_infections$time>=2020 &
                                                                   out3_no_inf$timeseries$total_chronic_infections$time<2030,-c(1,2)],
               2,sum)/
           apply(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time>=2020 &
                                                            out2$timeseries$total_chronic_infections$time<2030,-c(1,2)],
                 2,sum), c(0.5,0.025,0.975))
# 10% (6-17%)

# Chronic infection incidence rate over time
# Shows that averting of infections happens in the short term (by 2030)
hbv_inc_over_time <-
  rbind(gather(out2$timeseries$total_chronic_infections_rate, key="sim", value = "value", -time,-scenario),
        gather(out3_it$timeseries$total_chronic_infections_rate, key="sim", value = "value", -time,-scenario),
        gather(out3_no_inf$timeseries$total_chronic_infections_rate, key="sim", value = "value", -time,-scenario))

hbv_inc_rate_reduction <- gather((out2$timeseries$total_chronic_infections_rate[,-c(1,2)]-
                                    out3_no_inf$timeseries$total_chronic_infections_rate[,-c(1,2)])/
    out2$timeseries$total_chronic_infections_rate[,-c(1,2)], key="sim", value = "value")
hbv_inc_rate_reduction$time <- rep(out2$timeseries$total_chronic_infections_rate$time, 183)

quantile(hbv_inc_rate_reduction[hbv_inc_rate_reduction$time==2030,]$value, c(0.5,0.025,0.975))
quantile(hbv_inc_rate_reduction[hbv_inc_rate_reduction$time==2050,]$value, c(0.5,0.025,0.975))

hbv_inc_over_time <- hbv_inc_over_time %>%
  group_by(time, scenario) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value,0.975))

# Incidence plot:
ggplot() +
  geom_ribbon(data=subset(hbv_inc_over_time, scenario == "status_quo"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5,fill = scenario, colour=scenario),
              linetype = "solid", alpha = 0) +
  geom_line(data=subset(hbv_inc_over_time, scenario == "status_quo"),
            aes(x=time, y = median/0.5), colour="#B2182B",
            linetype = "solid", size=0.75) +
  geom_ribbon(data=subset(hbv_inc_over_time, scenario == "screen_2020_monit_0_no_treated_infectivity"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, fill = scenario, colour=scenario),
              linetype = "solid",alpha = 0) +
  geom_line(data=subset(hbv_inc_over_time, scenario == "screen_2020_monit_0_no_treated_infectivity"),
            aes(x=time, y = median/0.5), colour="#2166AC",
            linetype = "solid", size=0.75) +
#  scale_fill_manual("Scenario",
#                    labels = c("status_quo" = "Infant vaccination only",
#                               "screen_2020_monit_sim7" = "Screening & treatment in 2020\n(monitor every 5 years\nin <45 year olds)"),
#                    values=c("status_quo" = "#D6604D",
#                             "screen_2020_monit_sim7" = "#92C5DE")) +
#  scale_colour_manual("Scenario",
#                      labels = c("status_quo" = "Infant vaccination only",
#                                 "screen_2020_monit_sim7" = "Screening & treatment in 2020\n(monitor every 5 years\nin <45 year olds)"),
#                      values=c("status_quo" = "#B2182B",
#                               "screen_2020_monit_sim7" = "#2166AC")) +
 # guides(fill=FALSE, colour = FALSE) +
  geom_vline(xintercept=2019.5, linetype="dashed") +
  ylab("New chronic HBV infections per person per year") +
  ylim(0,0.0015) +
  xlim(2015,2040) +
  theme_classic() +
  theme(legend.position=c(.78,.8),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13))

# Probably at some point it's also a question of built up immunity in the population.

# For appendix: comparison of status quo with no-vaccination scenario ----
hbv_deaths_over_time <-
  rbind(gather(out2$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario),
        gather(out_no_vacc$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario))

hbv_deaths_over_time <- hbv_deaths_over_time %>%
  group_by(time, scenario) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value,0.975))

hbv_inc_over_time <-
  rbind(gather(out2$timeseries$total_chronic_infections, key="sim", value = "value", -time,-scenario),
        gather(out_no_vacc$timeseries$total_chronic_infections, key="sim", value = "value", -time,-scenario))

hbv_inc_over_time <- hbv_inc_over_time %>%
  group_by(time, scenario) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value,0.975))


inc_plot <- ggplot() +
  geom_ribbon(data=subset(hbv_inc_over_time, scenario == "no_vacc"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, fill = scenario, colour=scenario),
              linetype = "solid",alpha = 0.8) +
  geom_line(data=subset(hbv_inc_over_time, scenario == "no_vacc"),
            aes(x=time, y = median/0.5), colour="#2166AC",
            linetype = "solid", size=0.75) +
  geom_ribbon(data=subset(hbv_inc_over_time, scenario == "status_quo"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5,fill = scenario, colour=scenario),
              linetype = "solid", alpha = 0.5) +
  geom_line(data=subset(hbv_inc_over_time, scenario == "status_quo"),
            aes(x=time, y = median/0.5), colour="#DC663A",
            linetype = "solid", size=0.75) +
  scale_fill_manual("Scenario",
                    labels = c("status_quo" = "Base case (infant vaccination)",
                               "no_vacc" = "No historical intervention"),
                    values=c("status_quo" = "#DC663A",
                             "no_vacc" = "#92C5DE")) +
  scale_colour_manual("Scenario",
                      labels = c("status_quo" = "Base case (infant vaccination)",
                                 "no_vacc" = "No historical intervention"),
                      values=c("status_quo" = "#DC663A",
                               "no_vacc" = "#2166AC")) +
  geom_vline(xintercept=1990, linetype="dashed", colour="grey30") +
  # guides(fill=FALSE, colour = FALSE) +
  ylab("Annual incident chronic HBV infections") +
  xlab("Year") +
  xlim(1985,2080) +
  theme_classic() +
  theme(legend.position=c(.2,.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_blank())


deaths_plot <- ggplot() +
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "no_vacc"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, fill = scenario, colour=scenario),
              linetype = "solid",alpha = 0.8) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "no_vacc"),
            aes(x=time, y = median/0.5), colour="#2166AC",
            linetype = "solid", size=0.75) +
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "status_quo"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5,fill = scenario, colour=scenario),
              linetype = "solid", alpha = 0.5) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "status_quo"),
            aes(x=time, y = median/0.5), colour="#DC663A",
            linetype = "solid", size=0.75) +
  scale_fill_manual("Scenario",
                    labels = c("status_quo" = "Continued infant vaccination\n(base case)",
                               "no_vacc" = "No historical intervention"),
                    values=c("status_quo" = "#DC663A",
                             "no_vacc" = "#92C5DE")) +
  scale_colour_manual("Scenario",
                      labels = c("status_quo" = "Continued infant vaccination\n(base case)",
                                 "no_vacc" = "No historical intervention"),
                      values=c("status_quo" = "#DC663A",
                               "no_vacc" = "#2166AC")) +
  geom_vline(xintercept=1990, linetype="dashed", colour="grey30") +
  # guides(fill=FALSE, colour = FALSE) +
  ylab("Annual HBV-related deaths") +
  xlab("Year") +
  xlim(1985,2080) +
  theme_classic() +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_blank())

library(grid)
inc_plot_a <- arrangeGrob(inc_plot, top = textGrob("A", x = unit(0.01, "npc"),
                                                   y   = unit(1, "npc"), just=c("left","top"),
                                                   gp=gpar(col="black", fontsize=20)))
deaths_plot_b <- arrangeGrob(deaths_plot, top = textGrob("B", x = unit(0.01, "npc"),
                                                         y   = unit(1, "npc"), just=c("left","top"),
                                                         gp=gpar(col="black", fontsize=20)))

#png(file = "no_historical_intervention_plot.png", width=300, height=260, units = "mm", res=300, pointsize = 0.99)
grid.arrange(inc_plot_a, deaths_plot_b, nrow=2)
#dev.off()
