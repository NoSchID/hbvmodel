# Results for IVHEM poster
require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
library(BCEA)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))


# Load data
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

# Status quo
out1 <- readRDS(paste0(out_path, "a1_out1_status_quo_cohort_240920.rds"))
out1 <- out1[[1]]
out2 <- readRDS(paste0(out_path, "out2_status_quo_180820.rds"))
out2 <- out2[[1]]

# No monitoring by age
# A4 (15-30) A5 (30-45) A2 (45-60)
a4_out3 <- readRDS(paste0(out_path, "a4_out3_screen_2020_monit_0_231020.rds"))
a4_out3 <- a4_out3[[1]]
a5_out3 <- readRDS(paste0(out_path, "a5_out3_screen_2020_monit_0_231020.rds"))
a5_out3 <- a5_out3[[1]]
a2_out3 <- readRDS(paste0(out_path, "a2_out3_screen_2020_monit_0_161020.rds"))
a2_out3 <- a2_out3[[1]]

# Age group analysis

# Number of deaths averted in each age group
deaths_averted_by_age1 <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                                 scenario_objects = list(a4_out3),
                                 counterfactual_label = "treatment programme without monitoring")
deaths_averted_by_age2 <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a5_out3),
                          counterfactual_label = "treatment programme without monitoring")
deaths_averted_by_age3 <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a2_out3),
                          counterfactual_label = "treatment programme without monitoring")


deaths_averted_by_age <- rbind(
  data.frame(age_group = "15-30",
                                  filter(deaths_averted_by_age1, by_year == 2100 &
                                    type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age2, by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-60",
             filter(deaths_averted_by_age3, by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value))
  )

ggplot(deaths_averted_by_age) +
  geom_boxplot(aes(x=age_group, y=value)) +
  ylim(0,4500)

# Need to get carriers by age in 2020

deaths_per_int1 <- plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(a4_out3),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "no treatment")
deaths_per_int2 <- plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                                      scenario_objects = list(a5_out3),
                                                                      interaction_type = "total_interactions",
                                                                      counterfactual_label = "no treatment")
deaths_per_int3 <- plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                                      scenario_objects = list(a2_out3),
                                                                      interaction_type = "total_interactions",
                                                                      counterfactual_label = "no treatment")

deaths_per_int <- rbind(
  data.frame(age_group = "15-30",
             filter(deaths_per_int1, by_year == 2100) %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_per_int2, by_year == 2100) %>% select(scenario, sim, value)),
  data.frame(age_group = "45-60",
             filter(deaths_per_int3, by_year == 2100) %>% select(scenario, sim, value))
)

ggplot(deaths_per_int) +
  geom_boxplot(aes(x=age_group, y=value*10000)) +
  ylab("HBV-related deaths averted per 10,000 interactions")
