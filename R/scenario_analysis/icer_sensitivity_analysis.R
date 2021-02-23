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

# Default simulations for comparison:
out3_it <- readRDS(paste0(out_path_monit, "a1_it_out3_screen_2020_monit_0_180121.rds"))
out3_it <- out3_it[[1]]   # No monitoring
monit_out7 <- readRDS(paste0(out_path_monit, "a1_it_monit_out7_161220.rds"))
monit_out7 <- monit_out7[[1]]
out5_it <- readRDS(paste0(out_path_monit, "a1_it_out5_screen_2020_monit_5_161220.rds"))
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

# Technically out3 is dominated, but get ICER compared to SQ anyway for comparison purposes
# The ICER for monit_sim7 compared to monit_sim6 in the other file is: 338 (161-844)
# If monit_sim6 was excluded it would be: 271 (163-519)

# Create dataframes with cost (load functions from other script) ----
annual_discounting_rate <- 0.03

# CHANGE HERE FOR OBJECTS OF INTEREST
#object_list <- list(out3_it, monit_out7)
# Qualitative variations:
#object_list <- list(out3_s1, monit_out7_s1)
#object_list <- list(out3_s2, monit_out7_s2)
#object_list <- list(out3_s3, monit_out7_s3)
# One-way:
#object_list <- list(out3_tmu_dcc_lower, monit_out7_tmu_dcc_lower)
#object_list <- list(out3_tmu_dcc_upper, monit_out7_tmu_dcc_upper)
#object_list <- list(out3_thccr_chb_lower, monit_out7_thccr_chb_lower)
#object_list <- list(out3_thccr_chb_upper, monit_out7_thccr_chb_upper)
#object_list <- list(out3_thccr_it_lower, monit_out7_thccr_it_lower)
#object_list <- list(out3_thccr_it_upper, monit_out7_thccr_it_upper)
#object_list <- list(out3_thccr_cc_lower, monit_out7_thccr_cc_lower)
#object_list <- list(out3_thccr_cc_upper, monit_out7_thccr_cc_upper)
#object_list <- list(out3_thccr_dcc_lower, monit_out7_thccr_dcc_lower)
object_list <- list(out3_thccr_dcc_upper, monit_out7_thccr_dcc_upper)

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

# Tornado plot ----
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


