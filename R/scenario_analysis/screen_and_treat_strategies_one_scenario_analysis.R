# Strategies for screening and treatment simulation - 1 set of assumptions (A, B, C etc)

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

## Load files (check manually) ----
# Assumption A set ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/screen_and_treat_strategies/Scenario A with new parmsets/"

# Status quo
out1 <- readRDS(paste0(out_path, "a_out1_status_quo_cohort_110620.rds"))   # was out0_cohort
out1 <- out1[[1]]
out2 <- readRDS(paste0(out_path, "a_out2_status_quo_080720.rds"))          # was out1
out2 <- out2[[1]]

# Monitoring
out3 <- readRDS(paste0(out_path, "a_out3_screen_2020_monit_0_080620.rds"))
out3 <- out3[[1]]
out4 <- readRDS(paste0(out_path, "a_out4_screen_2020_monit_10_080620.rds"))
out4 <- out4[[1]]
out5 <- readRDS(paste0(out_path, "a_out5_screen_2020_monit_5_080620.rds"))
out5 <- out5[[1]]
out6 <- readRDS(paste0(out_path, "a_out6_screen_2020_monit_1_080620.rds"))
out6 <- out6[[1]]

# Screening
out7 <- readRDS(paste0(out_path, "a_out7_monit_0_screen_20_090620.rds"))
out7 <- out7[[1]]
out8 <- readRDS(paste0(out_path, "a_out8_monit_0_screen_10_090620.rds"))
out8 <- out8[[1]]
out9 <- readRDS(paste0(out_path, "a_out9_monit_0_screen_5_090620.rds"))
out9 <- out9[[1]]
out10 <- readRDS(paste0(out_path, "a_out10_monit_0_screen_1_090620.rds"))
out10 <- out10[[1]]

# Combination of monitoring and screening
out11 <- readRDS(paste0(out_path, "a_out11_monit_5_screen_5_150720.rds"))
out11 <- out11[[1]]

# Assumption B set ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/screen_and_treat_strategies/Scenario B/"

# Status quo
out1 <- readRDS(paste0(out_path, "b_out1_bd_scale_up_cohort_030720.rds"))
out1 <- out1[[1]]
out2 <- readRDS(paste0(out_path, "b_out2_bd_scale_up_080720.rds"))
out2 <- out2[[1]]

# Monitoring
out3 <- readRDS(paste0(out_path, "b_out3_screen_2020_monit_0_030720.rds"))
out3 <- out3[[1]]
out4 <- readRDS(paste0(out_path, "b_out4_screen_2020_monit_10_030720.rds"))
out4 <- out4[[1]]
out5 <- readRDS(paste0(out_path, "b_out5_screen_2020_monit_5_030720.rds"))
out5 <- out5[[1]]
out6 <- readRDS(paste0(out_path, "b_out6_screen_2020_monit_1_030720.rds"))
out6 <- out6[[1]]

# Screening
out7 <- readRDS(paste0(out_path, "b_out7_monit_0_screen_20_050720.rds"))
out7 <- out7[[1]]
out8 <- readRDS(paste0(out_path, "b_out8_monit_0_screen_10_050720.rds"))
out8 <- out8[[1]]
out9 <- readRDS(paste0(out_path, "b_out9_monit_0_screen_5_050720.rds"))
out9 <- out9[[1]]
out10 <- readRDS(paste0(out_path, "b_out10_monit_0_screen_1_050720.rds"))
out10 <- out10[[1]]

# Combination of monitoring and screening
out11 <- readRDS(paste0(out_path, "b_out11_monit_5_screen_5_150720.rds"))
out11 <- out11[[1]]

# Assumption D1 set ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/screen_and_treat_strategies/Scenario D1/"

# Status quo
out1 <- readRDS(paste0(out_path, "out1_status_quo_cohort_110620.rds"))
out1 <- out1[[1]]
out2 <- readRDS(paste0(out_path, "out2_status_quo_080720.rds"))
out2 <- out2[[1]]

# Monitoring
out3 <- readRDS(paste0(out_path, "d1_out3_screen_2020_monit_0_090720.rds"))
out3 <- out3[[1]]
out4 <- readRDS(paste0(out_path, "d1_out4_screen_2020_monit_10_090720.rds"))
out4 <- out4[[1]]
out5 <- readRDS(paste0(out_path, "d1_out5_screen_2020_monit_5_090720.rds"))
out5 <- out5[[1]]
out6 <- readRDS(paste0(out_path, "d1_out6_screen_2020_monit_1_090720.rds"))
out6 <- out6[[1]]

# Screening
out7 <- readRDS(paste0(out_path, "d1_out7_monit_0_screen_20_090720.rds"))
out7 <- out7[[1]]
out8 <- readRDS(paste0(out_path, "d1_out8_monit_0_screen_10_090720.rds"))
out8 <- out8[[1]]
out9 <- readRDS(paste0(out_path, "d1_out9_monit_0_screen_5_090720.rds"))
out9 <- out9[[1]]
out10 <- readRDS(paste0(out_path, "d1_out10_monit_0_screen_1_090720.rds"))
out10 <- out10[[1]]

# Combination of monitoring and screening
out11 <- readRDS(paste0(out_path, "d1_out11_monit_5_screen_5_150720.rds"))
out11 <- out11[[1]]


## IMPACT OF MONITORING ----
### Cohort outcomes of monitoring (uses functions from calculate_outcomes.R) ----

# Note that, for out3, all cohorts have died by 2102 (cohort size < 1)
# The cohort we are following are those screened individuals who immediately engage with care

# Summary of average age at death for final table
cohort_age_at_death <- data.frame(rbind(out1$cohort_age_at_death,
                                    out3$cohort_age_at_death,
                                    out4$cohort_age_at_death,
                                    out5$cohort_age_at_death,
                                    out6$cohort_age_at_death))

cohort_age_at_death_long <- gather(cohort_age_at_death, key = "sim", value = "value", -scenario)

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

age_at_death_ext_long <- gather(age_at_death_ext, key = "scenario", value = "value")

age_at_death_ext_summary <- age_at_death_ext_long %>%
  group_by(scenario) %>%
  summarise(median = median(value),
            lower = quantile(value, prob = 0.025),
            upper = quantile(value, prob = 0.975))

# Boxplot
ggplot(age_at_death_ext_long) +
  geom_boxplot(aes(x = scenario, y = value), fill = "#F8766D", width = 0.25) +
  ylab("Extension in average age at death (years)") +
  xlab("Monitoring frequency") +
  labs(title= "Cohort impact compared to treatment programme without monitoring") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))

# Credible intervals
ggplot(age_at_death_ext_summary) +
  geom_point(aes(x = scenario, y = median*12), size = 5) +
  geom_errorbar(aes(x = scenario, ymin = lower*12, ymax = upper*12), width = 0.2) +
  ylab("Extension in average age at death (months)") +
  xlab("Monitoring frequency") +
  labs(title= "Cohort impact compared to treatment programme without monitoring") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  ylim(0,15)

# Compare cohort number of HBV deaths averted compared to no monitoring
cohort_deaths_averted_long <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
                                      scenario_objects = list(out4, out5, out6),
                                      counterfactual_label = "treatment programme without monitoring")

# Compare cohort number of life years gained compared to no monitoring
cohort_ly_gained_long <-
  plot_ly_gained_cohort(counterfactual_object = out3,
                        scenario_objects = list(out4, out5, out6),
                        counterfactual_label = "treatment programme without monitoring")

# COUNTERFACTUAL = No treatment - status quo (out0)

# Compare median age at death
age_at_death_sq <- data.frame(rbind(out1$cohort_age_at_death,
                                 out3$cohort_age_at_death,
                                 out4$cohort_age_at_death,
                                 out5$cohort_age_at_death,
                                 out6$cohort_age_at_death))
age_at_death_sq <- t(age_at_death_sq[,-1])
colnames(age_at_death_sq) <- c("No treatment", "No monitoring", "Every 10 years", "Every 5 years", "Every year")
boxplot(age_at_death_sq, ylim =c(68,73), ylab = "Mean age at death (years)", main = "Assessed cohort")

# Extension in age at death compared to no treatment (status quo)
age_at_death_ext_sq <- data.frame(cbind(age_at_death_sq[,2]-age_at_death_sq[,1],
                                     age_at_death_sq[,3]-age_at_death_sq[,1],
                                     age_at_death_sq[,4]-age_at_death_sq[,1],
                                     age_at_death_sq[,5]-age_at_death_sq[,1]))
colnames(age_at_death_ext_sq) <- c("No monitoring", "Every 10 years", "Every 5 years", "Every year")
boxplot(age_at_death_ext_sq, ylab = "Extension in mean age at death (years)",
        main = "Compared to no treatment (infant vaccine only)",
        xlab = "Monitoring frequency", ylim = c(0,5))

age_at_death_ext_sq_long <- gather(age_at_death_ext_sq, key = "scenario", value = "value")
age_at_death_ext_sq_long$scenario <- factor(age_at_death_ext_sq_long$scenario, levels =
                                             c("No monitoring", "Every 10 years", "Every 5 years", "Every year"))

# Boxplot
ggplot(age_at_death_ext_sq_long) +
  geom_boxplot(aes(x = scenario, y = value), fill = "#F8766D", width = 0.25) +
  ylab("Extension in average age at death (years)") +
  xlab("Monitoring frequency") +
  labs(title= "Cohort impact compared to no treatment (status quo scenario)") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15)) +
  ylim(0,5)

# Compare cohort number of HBV deaths averted compared to no treatment (status quo)
cohort_deaths_averted_sq_long <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out1,
                                      scenario_objects = list(out3,out4, out5, out6),
                                      counterfactual_label = "no treatment")

# Compare cohort number of life years gained compared to no monitoring
cohort_ly_gained_sq_long <-
  plot_ly_gained_cohort(counterfactual_object = out1,
                        scenario_objects = list(out3,out4, out5, out6),
                        counterfactual_label = "no treatment")

### Population outcomes (uses functions from calculate_outcomes.R) ----

# Labels for time periods on plot
period_labs <- c("2020-2030", "2020-2050", "2020-2100")
names(period_labs) <- c("2030", "2050", "2100")

# HBV DEATHS AVERTED

# COUNTERFACTUAL = STATUS QUO
deaths_averted_sq_long <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                scenario_objects = list(out3, out4, out5, out6),
                                                counterfactual_label = "no treatment programme")
# This plot suggests that maximum monitoring (yearly) has little additional effect in preventing deaths by 2030,
# but there is a bigger difference between yearly vs. no monitoring by 2050, and a medium difference by 2100
# In other words the immediate gains from treatment come from the initial treatment, in the medium term
# monitoring can prevent some deaths, but in the long term this additional gain diminishes because
# there is lack of further screening

# Plot proportion of deaths averted by year on x axis, for a no monitoring and a frequent monitoring scenario
prop_deaths_averted_by_year <- deaths_averted_sq_long[deaths_averted_sq_long$type == "proportion_averted" &
                                                        (deaths_averted_sq_long$scenario %in%
                                                           c("screen_2020_monit_0", "screen_2020_monit_1")),]
prop_deaths_averted_by_year$scenario_label[prop_deaths_averted_by_year$scenario == "screen_2020_monit_0"] <- "No monitoring"
prop_deaths_averted_by_year$scenario_label[prop_deaths_averted_by_year$scenario == "screen_2020_monit_1"] <-
  "Yearly monitoring"

ggplot(prop_deaths_averted_by_year) +
  geom_boxplot(aes(x = as.factor(by_year), y = value), fill = "#00BFC4") +
  facet_wrap(~ scenario_label) +
  scale_x_discrete(labels = c("2020-2030", "2020-2050", "2020-2100")) +
  xlab("Period") +
  ylab("Proportion of HBV-related deaths averted") +
  labs(title= "Population impact compared to no treatment (status quo scenario)") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15)) +
  ylim(0,0.55)
# Impact of treatment programme in preventing HBV-related deaths is immediate - largest proportion
# of HBV deaths averted by 2030. Monitoring reduces the gap between the immediate vs. long-term impact
# of the treatment programme a little


# COUNTERFACTUAL = ONE-OFF SCREENING BUT NO MONITORING (out3)
deaths_averted_long <- plot_hbv_deaths_averted(counterfactual_object = out3,
                                               scenario_objects = list(out4, out5, out6),
                                               counterfactual_label = "treatment programme without monitoring",
                                               outcome_to_plot = "number_averted")
# Gain from monitoring is most visible in the medium term (2050)

# LIFE YEARS GAINED

# COUNTERFACTUAL = STATUS QUO
# Population-level effect of screening/treatment/monitoring in the short and long term
ly_gained_sq_long <- plot_ly_gained(counterfactual_object = out2,
                                       scenario_objects = list(out3, out4, out5, out6),
                                       counterfactual_label = "treatment programme without monitoring")

# Plot proportion of life-years saved by year on x axis, for a no monitoring and a frequent monitoring scenario
ly_saved_by_year <- ly_gained_sq_long[ly_gained_sq_long$type == "proportion_averted" &
                                                        (ly_gained_sq_long$counterfactual %in%
                                                           c("screen_2020_monit_0", "screen_2020_monit_1")),]
ly_saved_by_year$scenario_label[ly_saved_by_year$counterfactual == "screen_2020_monit_0"] <- "No monitoring"
ly_saved_by_year$scenario_label[ly_saved_by_year$counterfactual == "screen_2020_monit_1"] <-
  "Yearly monitoring"

ggplot(ly_saved_by_year) +
  geom_boxplot(aes(x = as.factor(by_year), y = value), fill = "#00BFC4") +
  facet_wrap(~ scenario_label) +
  scale_x_discrete(labels = c("2020-2030", "2020-2050", "2020-2100")) +
  xlab("Period") +
  ylab("Proportion of life-years saved") +
  labs(title= "Population impact compared to no treatment (status quo scenario)") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))

# COUNTERFACTUAL = ONE-OFF SCREENING BUT NO MONITORING (out3)
ly_gained_long <- plot_ly_gained(counterfactual_object = out3,
                                 scenario_objects = list(out4, out5, out6),
                                 counterfactual_label = "treatment programme without monitoring")
# Gain from monitoring is most visible in the medium term (2050) (pattern (ratio) seems to be
# slightly stronger than for HBV deaths)

### Cohort-population level comparison (uses objects from previous 2 sections) ----

# 1) COUNTERFACTUAL = STATUS QUO
# This shows more in general the impact of this screening and treatment programme
# Rather than effect of monitoring specifically

# OUTCOME = HBV DEATHS AVERTED

# Combine cohort and population outcome by 2100 into the same dataframe to plot
deaths_averted_sq_pop <- deaths_averted_sq_long[deaths_averted_sq_long$by_year == 2100,-c(1,2)]
deaths_averted_sq_pop$level <- "Population-level"
deaths_averted_sq_cohort <- cohort_deaths_averted_sq_long
deaths_averted_sq_cohort$level <- "Cohort-level"
deaths_averted_sq_pop_and_cohort <- rbind(deaths_averted_sq_pop, deaths_averted_sq_cohort)

# Compare population vs cohort-level effect of treatment
ggplot(deaths_averted_sq_pop_and_cohort[deaths_averted_sq_pop_and_cohort$type == "proportion_averted",]) +
  geom_boxplot(aes(x = scenario, y = value, fill = level)) +
  scale_x_discrete(labels = c("No monitoring", "10 years", "5 years", "1 year")) +
  xlab("Monitoring frequency") +
  ylab("Proportion of HBV-related deaths averted") +
  labs(title= "Impact by 2100 compared to no treatment (status quo scenario)",
       fill = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15)) +
  ylim(0,1)

# Number of deaths averted in the cohort is the same as in the entire population by 2100:
# (when the whole cohort has died)
ggplot(deaths_averted_sq_pop_and_cohort[deaths_averted_sq_pop_and_cohort$type == "number_averted",]) +
  geom_boxplot(aes(x = scenario, y = value, fill = level)) +
  scale_x_discrete(labels = c("No monitoring", "10 years", "5 years", "1 year")) +
  xlab("Monitoring frequency") +
  ylab("Cumulative number of HBV-related deaths averted by treatment") +
  labs(title= "HBV deaths averted by treatment by 2100 (number)") +
  ylim(0,11000)

# OUTCOME = LIFE YEARS SAVED

# Combine cohort and population outcome by 2100 into the same dataframe to plot
ly_gained_sq_pop <- ly_gained_sq_long[ly_gained_sq_long$by_year == 2100,-c(1,2)]
ly_gained_sq_pop$level <- "pop"
ly_gained_sq_cohort <- cohort_ly_gained_sq_long
ly_gained_sq_cohort$level <- "cohort"
ly_gained_sq_pop_and_cohort <- rbind(ly_gained_sq_pop, ly_gained_sq_cohort)

ggplot(ly_gained_sq_pop_and_cohort[ly_gained_sq_pop_and_cohort$type == "proportion_averted",]) +
  geom_boxplot(aes(x = counterfactual, y = value, fill = level)) +
  scale_x_discrete(labels = c("No monitoring", "10 years", "5 years", "1 year")) +
  xlab("Monitoring frequency") +
  ylab("Fraction of life-years saved compared to status quo") +
  labs(title= "Life-years saved by treatment by 2100 (proportion)")

ggplot(ly_gained_sq_pop_and_cohort[ly_gained_sq_pop_and_cohort$type == "number_averted",]) +
  geom_boxplot(aes(x = counterfactual, y = value, fill = level)) +
  scale_x_discrete(labels = c("No monitoring", "10 years", "5 years", "1 year")) +
  xlab("Monitoring frequency") +
  ylab("Life-years saved compared to status quo") +
  labs(title= "Life-years saved by treatment by 2100 (number)") +
  ylim(0,300000)

# 2) COUNTERFACTUAL = ONE-OFF SCREENING BUT NO MONITORING (out3)
# Shows the effect of monitoring vs no monitoring of treatment-ineligibles

# OUTCOME = HBV DEATHS AVERTED

# Combine cohort and population outcome by 2100 into the same dataframe to plot
deaths_averted_pop <- deaths_averted_long[deaths_averted_long$by_year == 2100,-c(1,2)]
deaths_averted_pop$level <- "Population-level"
deaths_averted_cohort <- cohort_deaths_averted_long
deaths_averted_cohort$level <- "Cohort-level"
deaths_averted_pop_and_cohort <- rbind(deaths_averted_pop, deaths_averted_cohort)

# Compare population vs cohort-level effect of monitoring
ggplot(deaths_averted_pop_and_cohort[deaths_averted_pop_and_cohort$type == "proportion_averted",]) +
  geom_boxplot(aes(x = scenario, y = value, fill = level)) +
  scale_x_discrete(labels = c("10 years", "5 years", "1 year")) +
  xlab("Monitoring frequency") +
  ylab("Proportion of HBV-related deaths averted") +
  labs(title= "Impact by 2100 compared to treatment programme without monitoring",
       fill = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15)) +
  ylim(0,1)
# Monitoring is beneficial in the cohort but has less population-level benefit

# OUTCOME = LIFE YEARS SAVED

# Combine cohort and population outcome by 2100 into the same dataframe to plot
ly_gained_pop <- ly_gained_long[ly_gained_long$by_year == 2100,-c(1,2)]
ly_gained_pop$level <- "Population-level"
ly_gained_cohort <- cohort_ly_gained_long
ly_gained_cohort$level <- "Cohort-level"
ly_gained_pop_and_cohort <- rbind(ly_gained_pop, ly_gained_cohort)

# Compare population vs cohort-level effect of monitoring
ggplot(ly_gained_pop_and_cohort[ly_gained_pop_and_cohort$type == "proportion_averted",]) +
  geom_boxplot(aes(x = counterfactual, y = value, fill = level)) +
  scale_x_discrete(labels = c("10 years", "5 years", "1 year")) +
  xlab("Monitoring frequency") +
  ylab("Proportion of life-years saved") +
  labs(title= "Impact by 2100 compared to treatment programme without monitoring",
       fill = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))
# Monitoring is beneficial in the cohort but has less population-level benefit

### Population outcomes per healthcare interaction (uses functions from calculate_outcomes.R) ----

# COUNTERFACTUAL = STATUS QUO

# OUTCOME = HBV DEATHS AVERTED PER INCREMENTAL TREATMENT INITIATION
# compared to status quo, any healthcare interactions are incremental
deaths_averted_per_treated_sq_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                        scenario_objects = list(out3, out4,out5, out6),
                                                        interaction_type = "total_treated",
                                                        counterfactual_label = "no treatment programme")

# OUTCOME = HBV DEATHS AVERTED PER INCREMENTAL INTERACTION
deaths_averted_per_interaction_sq_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                        scenario_objects = list(out3, out4,out5, out6),
                                                        interaction_type = "total_interactions",
                                                        counterfactual_label = "no treatment programme")

# OUTCOME = LIFE YEARS GAINED PER INCREMENTAL HEALTHCARE INTERACTION

# One-off screen monitoring every 10 years
ly_gained_per_interaction_sq_long <-
  plot_ly_gained_per_healthcare_interaction(counterfactual_object = out2,
                                            scenario_objects = list(out3, out4,out5, out6),
                                            interaction_type = "total_interactions",
                                            counterfactual_label = "no treatment programme")

# COUNTERFACTUAL = TREATMENT BUT NO MONITORING (out3)

# OUTCOME = HBV DEATHS AVERTED PER INCREMENTAL HEALTHCARE INTERACTION
deaths_averted_per_interaction_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects = list(out4,out5, out6),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "treatment programme without monitoring")

# OUTCOME = LIFE YEARS GAINED PER INCREMENTAL HEALTHCARE INTERACTION
ly_gained_per_interaction_long <-
  plot_ly_gained_per_healthcare_interaction(counterfactual_object = out3,
                                            scenario_objects = list(out4,out5, out6),
                                            interaction_type = "total_interactions",
                                            counterfactual_label = "treatment programme without monitoring")


# Need to think about what I would like to represent:
# compared to monitoring, not monitoring averts more deaths per incremental clinical assessments and
# treatment initiations,
# and averts a similar number of deaths per incremental healthcare interaction
# but Maud has said than screening people is relatively easy whereas the clinical assessment is difficult

## IMPACT OF SCREENING ----
# Comparing out7-10 (screening frequencies) to out3 (one-off screening) and out2 (no treatment)

### Population outcomes ----

# HBV DEATHS AVERTED

# COUNTERFACTUAL = STATUS QUO
deaths_averted_sq_screen_long <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                     scenario_objects = list(out3, out7, out8, out9, out10),
                                                     counterfactual_label = "no treatment programme",
                                                     x_axis = "screening")
# Repeated screening, even if only every 20 years, affects proportion of HBV deaths averted
# particularly in the long-term. In the long term not much difference between 10, 5 or yearly screening,
# although yearly screening appears substantially better than every 5 years to avert deaths by 2030

# Plot proportion of deaths averted by year on x axis, for a no monitoring and a frequent monitoring scenario
prop_deaths_averted_by_year_screen <- deaths_averted_sq_screen_long[deaths_averted_sq_screen_long$type == "proportion_averted" &
                                                        (deaths_averted_sq_screen_long$scenario %in%
                                                           c("screen_2020_monit_0", "monit_0_screen_1")),]
prop_deaths_averted_by_year_screen$scenario_label[prop_deaths_averted_by_year_screen$scenario == "screen_2020_monit_0"] <- "One-off screening"
prop_deaths_averted_by_year_screen$scenario_label[prop_deaths_averted_by_year_screen$scenario == "monit_0_screen_1"] <-
  "Yearly screening"

ggplot(prop_deaths_averted_by_year_screen) +
  geom_boxplot(aes(x = as.factor(by_year), y = value), fill = "#00BFC4") +
  facet_wrap(~ scenario_label) +
  scale_x_discrete(labels = c("2020-2030", "2020-2050", "2020-2100")) +
  xlab("Period") +
  ylab("Proportion of HBV-related deaths averted") +
  labs(title= "Population impact compared to no treatment (status quo scenario)") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15)) +
  ylim(0,0.75)
# Less useful figure. Impact of yearly screening is already visible in the short term,
# but seems similar over time.


# COUNTERFACTUAL = ONE-OFF SCREENING BUT NO MONITORING (out3)
deaths_averted_screen_long <- plot_hbv_deaths_averted(counterfactual_object = out3,
                                               scenario_objects = list(out7, out8, out9, out10),
                                               counterfactual_label = "treatment programme with one-off screening",
                                               x_axis = "screening")

# Gain from yearly screening is similar at all timepoints. Gain from 10 and 5 yearly screening
# increases slightly over time and impact of different frequencies becomes more similar.
# Hypothesis: this could be because of increased impact of prevention over time.
# Order of magnitude of deaths averted by 2050 and 2100 is approx between 10 and 30%.

# LIFE YEARS GAINED

# COUNTERFACTUAL = STATUS QUO
ly_gained_sq_screen_long <- plot_ly_gained(counterfactual_object = out2,
                                        scenario_objects = list(out3, out7, out8, out9, out10),
                                        counterfactual_label = "no treatment programme",
                                        x_axis = "screening")

# Plot proportion of life-years saved by year on x axis, for a one-off and a frequent screening scenario
ly_saved_by_year_screen <- ly_gained_sq_screen_long[ly_gained_sq_screen_long$type == "proportion_averted" &
                                        (ly_gained_sq_screen_long$counterfactual %in%
                                           c("screen_2020_monit_0", "monit_0_screen_1")),]
ly_saved_by_year_screen$scenario_label[ly_saved_by_year_screen$counterfactual == "screen_2020_monit_0"] <- "One-off screening"
ly_saved_by_year_screen$scenario_label[ly_saved_by_year_screen$counterfactual == "monit_0_screen_1"] <-
  "Yearly screening"

ggplot(ly_saved_by_year_screen) +
  geom_boxplot(aes(x = as.factor(by_year), y = value), fill = "#00BFC4") +
  facet_wrap(~ scenario_label) +
  scale_x_discrete(labels = c("2020-2030", "2020-2050", "2020-2100")) +
  xlab("Period") +
  ylab("Proportion of life-years saved") +
  labs(title= "Population impact compared to no treatment (status quo scenario)") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))
# One of screen has very little effect on saving life years by 2100. Yearly screening improves that.

# COUNTERFACTUAL = ONE-OFF SCREENING AND NO MONITORING (out3)
ly_gained_screen_long <- plot_ly_gained(counterfactual_object = out3,
                                 scenario_objects = list(out7, out8, out9, out10),
                                 counterfactual_label = "treatment programme with one-off screening",
                                 x_axis = "screening")
# Gain from 10-yearly screening appears similar in the medium and long term, which are both larger
# than by 2030. This is a different pattern from HBV deaths, where the yearly screening had
# similar impact on all timescales.


# OUTCOME = HBV DEATHS AVERTED PER HEALTHCARE INTERACTION

# COUNTERFACTUAL = ONE-OFF SCREEN
deaths_averted_per_interaction_screen_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects = list(out7,out8, out9, out10),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "treatment programme with one-off screening",
                                                     x_axis = "screening")
# Note dropped variables (NA) are for screen every 20 and 10 years by 2030 (0 incremental interactions)
# For total interactions, 1 year requires by far the most incremental interactions per death averted,
# whereas the other screening frequencies are fairly similar. This is entirely due to HBsAg tests.
# For treatment initiations and assessments, this is actually higher for 5 years/20 years by 2030/2050
# than the other frequencies (all the same by 2100).
# By 2030, about 7,000 extra tests required to avert 1 death.

# COUNTERFACTUAL = NO TREATMENT PROGRAMME

# COUNTERFACTUAL = ONE-OFF SCREEN
deaths_averted_per_interaction_sq_screen_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(out3, out7,out8, out9, out10),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "no treatment programme",
                                                     x_axis = "screening")
# 20 and 10 years are fairly similar to one-off, compared to no treatment programme. 1 year frequency
# requires far more HBsAg tests per HBV death averted.

# OUTCOME = LIFE YEARS SAVED PER HEALTHCARE INTERACTION

# COUNTERFACTUAL = ONE-OFF SCREEN
ly_gained_per_interaction_screen_long <-
  plot_ly_gained_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects = list(out7,out8, out9, out10),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "treatment programme with one-off screening",
                                                     x_axis = "screening")
# Here in the medium term, 10 and 5 year screening frequency seems to require fewer extra tests
# per life-year saved than 20 year screening, but otherwise same pattern as for HBV deaths averted.


## IMPACT OF COMBINED MONITORING AND SCREENING  ----
### Population outcomes ----

# HBV DEATHS AVERTED

# COUNTERFACTUAL = STATUS QUO
deaths_averted_sq_combi_long <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                         scenario_objects = list(out3, out11),
                                                         counterfactual_label = "no treatment programme",
                                                         x_axis = "screening")

# COUNTERFACTUAL = ONE-OFF SCREENING BUT NO MONITORING (out3)
deaths_averted_combi_long <- plot_hbv_deaths_averted(counterfactual_object = out3,
                                                      scenario_objects = list(out11),
                                                      counterfactual_label = "treatment programme with one-off screening",
                                                      x_axis = "screening")

# LIFE YEARS GAINED

# COUNTERFACTUAL = STATUS QUO
ly_gained_sq_combi_long <- plot_ly_gained(counterfactual_object = out2,
                                           scenario_objects = list(out3, out11),
                                           counterfactual_label = "no treatment programme",
                                           x_axis = "screening")


# COUNTERFACTUAL = ONE-OFF SCREENING AND NO MONITORING (out3)
ly_gained_combi_long <- plot_ly_gained(counterfactual_object = out3,
                                        scenario_objects = list(out11),
                                        counterfactual_label = "treatment programme with one-off screening",
                                        x_axis = "screening")


# OUTCOME = HBV DEATHS AVERTED PER HEALTHCARE INTERACTION

# COUNTERFACTUAL = ONE-OFF SCREEN
deaths_averted_per_interaction_combi_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects = list(out11),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "treatment programme with one-off screening",
                                                     x_axis = "screening")

# COUNTERFACTUAL = NO TREATMENT PROGRAMME

# COUNTERFACTUAL = ONE-OFF SCREEN
deaths_averted_per_interaction_sq_combi_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(out3,out11),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "no treatment programme",
                                                     x_axis = "screening")

# OUTCOME = LIFE YEARS SAVED PER HEALTHCARE INTERACTION

# COUNTERFACTUAL = ONE-OFF SCREEN
ly_gained_per_interaction_combi_long <-
  plot_ly_gained_per_healthcare_interaction(counterfactual_object = out3,
                                            scenario_objects = list(out11),
                                            interaction_type = "total_interactions",
                                            counterfactual_label = "treatment programme with one-off screening",
                                            x_axis = "screening")

## TABLE OF ALL KEY OUTCOMES FOR 1 ASSUMPTION SET ----

# Missing age at death

scenario_d1_full_results <-
  list(
    # Monitoring analysis
    cohort_deaths_averted_long = cohort_deaths_averted_long,
    cohort_deaths_averted_sq_long = cohort_deaths_averted_sq_long,
    cohort_ly_gained_long = cohort_ly_gained_long,
    cohort_ly_gained_sq_long = cohort_ly_gained_sq_long,
    cohort_average_age_at_death_long = cohort_age_at_death_long,
    deaths_averted_long = deaths_averted_long,
    deaths_averted_sq_long = deaths_averted_sq_long,
    ly_gained_long = ly_gained_long,
    ly_gained_sq_long = ly_gained_sq_long,
    deaths_averted_per_interaction_long = deaths_averted_per_interaction_long,
    deaths_averted_per_interaction_sq_long = deaths_averted_per_interaction_sq_long,
    ly_gained_per_interaction_long = ly_gained_per_interaction_long,
    # Screening analysis
    deaths_averted_screen_long = deaths_averted_screen_long,
    deaths_averted_sq_screen_long = deaths_averted_sq_screen_long,
    ly_gained_screen_long = ly_gained_screen_long,
    ly_gained_sq_screen_long = ly_gained_sq_screen_long,
    deaths_averted_per_interaction_screen_long = deaths_averted_per_interaction_screen_long,
    deaths_averted_per_interaction_sq_screen_long = deaths_averted_per_interaction_sq_screen_long,
    ly_gained_per_interaction_screen_long = ly_gained_per_interaction_screen_long,
    # Combination of monitoring and screening
    deaths_averted_combi_long = deaths_averted_combi_long,
    deaths_averted_sq_combi_long = deaths_averted_sq_combi_long,
    ly_gained_combi_long = ly_gained_combi_long,
    ly_gained_sq_combi_long = ly_gained_sq_combi_long,
    deaths_averted_per_interaction_combi_long = deaths_averted_per_interaction_combi_long,
    deaths_averted_per_interaction_sq_combi_long = deaths_averted_per_interaction_sq_combi_long,
    ly_gained_per_interaction_combi_long = ly_gained_per_interaction_combi_long
  )
#saveRDS(scenario_d1_full_results, here("output", "screen_and_treat_results", "scenario_d1_full_results.rds"))

# Summary values (median, 2.5th and 97.5th percentile)
scenario_d1_summary_results <- list(
  # MONITORING ANALYSIS
  cohort_deaths_averted_long = (group_by(cohort_deaths_averted_long, counterfactual, scenario, type) %>%
                                  summarise(median = median(value),
                                            cri_lower = quantile(value, prob = 0.025),
                                            cri_upper = quantile(value, prob = 0.975))),
  cohort_deaths_averted_sq_long = (group_by(cohort_deaths_averted_sq_long, counterfactual, scenario, type) %>%
                                     summarise(median = median(value),
                                               cri_lower = quantile(value, prob = 0.025),
                                               cri_upper = quantile(value, prob = 0.975))),
  cohort_ly_gained_long = (group_by(cohort_ly_gained_long, counterfactual, scenario, type) %>%
                             summarise(median = median(value),
                                       cri_lower = quantile(value, prob = 0.025),
                                       cri_upper = quantile(value, prob = 0.975))),
  cohort_ly_gained_sq_long = (group_by(cohort_ly_gained_sq_long, counterfactual, scenario, type) %>%
                                summarise(median = median(value),
                                          cri_lower = quantile(value, prob = 0.025),
                                          cri_upper = quantile(value, prob = 0.975))),
  cohort_average_age_at_death_long = (group_by(cohort_age_at_death_long, scenario) %>%
                                        summarise(median = median(value),
                                                  cri_lower = quantile(value, prob = 0.025),
                                                  cri_upper = quantile(value, prob = 0.975))),
  deaths_averted_long = (group_by(deaths_averted_long, by_year,
                                  counterfactual, scenario, type) %>%
                           summarise(median = median(value),
                                     cri_lower = quantile(value, prob = 0.025),
                                     cri_upper = quantile(value, prob = 0.975))),
  deaths_averted_sq_long = (group_by(deaths_averted_sq_long, by_year,
                                     counterfactual, scenario, type) %>%
                              summarise(median = median(value),
                                        cri_lower = quantile(value, prob = 0.025),
                                        cri_upper = quantile(value, prob = 0.975))),
  ly_gained_long = (group_by(ly_gained_long, by_year,
                             counterfactual, scenario, type) %>%
                      summarise(median = median(value),
                                cri_lower = quantile(value, prob = 0.025),
                                cri_upper = quantile(value, prob = 0.975))),
  ly_gained_sq_long = (group_by(ly_gained_sq_long, by_year,
                                counterfactual, scenario, type) %>%
                         summarise(median = median(value),
                                   cri_lower = quantile(value, prob = 0.025),
                                   cri_upper = quantile(value, prob = 0.975))),
  interactions_per_death_averted_long =(group_by(deaths_averted_per_interaction_long,
                                                 by_year,scenario) %>%
                                          summarise(median = median(1/value),
                                                    cri_lower = quantile(1/value, prob = 0.025),
                                                    cri_upper = quantile(1/value, prob = 0.975))),
  interactions_per_death_averted_sq_long =(group_by(deaths_averted_per_interaction_sq_long,
                                                    by_year,scenario) %>%
                                             summarise(median = median(1/value),
                                                       cri_lower = quantile(1/value, prob = 0.025),
                                                       cri_upper = quantile(1/value, prob = 0.975))),
  interactions_per_ly_gained_long =(group_by(ly_gained_per_interaction_long,
                                            by_year,scenario) %>%
                                     summarise(median = median(1/value),
                                               cri_lower = quantile(1/value, prob = 0.025),
                                               cri_upper = quantile(1/value, prob = 0.975))),
  # SCREENING ANALYSIS
  deaths_averted_screen_long = (group_by(deaths_averted_screen_long, by_year,
                                         counterfactual, scenario, type) %>%
                                  summarise(median = median(value),
                                            cri_lower = quantile(value, prob = 0.025),
                                            cri_upper = quantile(value, prob = 0.975))),
  deaths_averted_sq_screen_long = (group_by(deaths_averted_sq_screen_long, by_year,
                                            counterfactual, scenario, type) %>%
                                     summarise(median = median(value),
                                               cri_lower = quantile(value, prob = 0.025),
                                               cri_upper = quantile(value, prob = 0.975))),
  ly_gained_screen_long = (group_by(ly_gained_screen_long, by_year,
                                    counterfactual, scenario, type) %>%
                             summarise(median = median(value),
                                       cri_lower = quantile(value, prob = 0.025),
                                       cri_upper = quantile(value, prob = 0.975))),
  ly_gained_sq_screen_long = (group_by(ly_gained_sq_screen_long, by_year,
                                       counterfactual, scenario, type) %>%
                                summarise(median = median(value),
                                          cri_lower = quantile(value, prob = 0.025),
                                          cri_upper = quantile(value, prob = 0.975))),
  interactions_per_death_averted_screen_long =(group_by(deaths_averted_per_interaction_long,
                                                        by_year,scenario) %>%
                                                 summarise(median = median(1/value),
                                                           cri_lower = quantile(1/value, prob = 0.025),
                                                           cri_upper = quantile(1/value, prob = 0.975))),
  interactions_per_death_averted_sq_screen_long =(group_by(deaths_averted_per_interaction_sq_long,
                                                           by_year,scenario) %>%
                                                    summarise(median = median(1/value),
                                                              cri_lower = quantile(1/value, prob = 0.025),
                                                              cri_upper = quantile(1/value, prob = 0.975))),
  interactions_per_ly_gained_screen_long =(group_by(ly_gained_per_interaction_long,
                                                   by_year,scenario) %>%
                                            summarise(median = median(1/value),
                                                      cri_lower = quantile(1/value, prob = 0.025),
                                                      cri_upper = quantile(1/value, prob = 0.975))),
  # MONITORING AND SCREENING COMBINATION
  deaths_averted_combi_long = (group_by(deaths_averted_combi_long, by_year,
                                         counterfactual, scenario, type) %>%
                                  summarise(median = median(value),
                                            cri_lower = quantile(value, prob = 0.025),
                                            cri_upper = quantile(value, prob = 0.975))),
  deaths_averted_sq_combi_long = (group_by(deaths_averted_sq_combi_long, by_year,
                                            counterfactual, scenario, type) %>%
                                     summarise(median = median(value),
                                               cri_lower = quantile(value, prob = 0.025),
                                               cri_upper = quantile(value, prob = 0.975))),
  ly_gained_combi_long = (group_by(ly_gained_combi_long, by_year,
                                    counterfactual, scenario, type) %>%
                             summarise(median = median(value),
                                       cri_lower = quantile(value, prob = 0.025),
                                       cri_upper = quantile(value, prob = 0.975))),
  ly_gained_sq_combi_long = (group_by(ly_gained_sq_combi_long, by_year,
                                       counterfactual, scenario, type) %>%
                                summarise(median = median(value),
                                          cri_lower = quantile(value, prob = 0.025),
                                          cri_upper = quantile(value, prob = 0.975))),
  interactions_per_death_averted_combi_long =(group_by(deaths_averted_per_interaction_long,
                                                        by_year,scenario) %>%
                                                 summarise(median = median(1/value),
                                                           cri_lower = quantile(1/value, prob = 0.025),
                                                           cri_upper = quantile(1/value, prob = 0.975))),
  interactions_per_death_averted_sq_combi_long =(group_by(deaths_averted_per_interaction_sq_long,
                                                           by_year,scenario) %>%
                                                    summarise(median = median(1/value),
                                                              cri_lower = quantile(1/value, prob = 0.025),
                                                              cri_upper = quantile(1/value, prob = 0.975))),
  interactions_per_ly_gained_combi_long =(group_by(ly_gained_per_interaction_long,
                                                    by_year,scenario) %>%
                                             summarise(median = median(1/value),
                                                       cri_lower = quantile(1/value, prob = 0.025),
                                                       cri_upper = quantile(1/value, prob = 0.975)))
)
#saveRDS(scenario_d1_summary_results, here("output", "screen_and_treat_results", "scenario_d1_summary_results.rds"))



## TIMESERIES PLOTS ----

hbv_deaths_rate <- rbind(cbind(out2$timeseries$total_hbv_deaths_rate[,c(1:2)],
                               median = apply(out2$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, median),
                               cri_lower = apply(out2$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, quantile, prob = 0.025),
                               cri_upper = apply(out2$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, quantile, prob = 0.975)),
                         cbind(out3$timeseries$total_hbv_deaths_rate[,c(1:2)],
                               median= apply(out3$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, median),
                               cri_lower = apply(out3$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, quantile, prob = 0.025),
                               cri_upper = apply(out3$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, quantile, prob = 0.975)),
                         cbind(out5$timeseries$total_hbv_deaths_rate[,c(1:2)],
                               median= apply(out5$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, median),
                               cri_lower = apply(out5$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, quantile, prob = 0.025),
                               cri_upper = apply(out5$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, quantile, prob = 0.975)),
                         cbind(out9$timeseries$total_hbv_deaths_rate[,c(1:2)],
                               median = apply(out9$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, median),
                               cri_lower = apply(out9$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, quantile, prob = 0.025),
                               cri_upper = apply(out9$timeseries$total_hbv_deaths_rate[,-c(1:2)], 1, quantile, prob = 0.975)))

hbv_deaths<- rbind(cbind(out2$timeseries$total_hbv_deaths[,c(1:2)],
                         median = apply(out2$timeseries$total_hbv_deaths[,-c(1:2)], 1, median),
                         cri_lower = apply(out2$timeseries$total_hbv_deaths[,-c(1:2)], 1, quantile, prob = 0.025),
                         cri_upper = apply(out2$timeseries$total_hbv_deaths[,-c(1:2)], 1, quantile, prob = 0.975)),
                   cbind(out3$timeseries$total_hbv_deaths_rate[,c(1:2)],
                         median= apply(out3$timeseries$total_hbv_deaths[,-c(1:2)], 1, median),
                         cri_lower = apply(out3$timeseries$total_hbv_deaths[,-c(1:2)], 1, quantile, prob = 0.025),
                         cri_upper = apply(out3$timeseries$total_hbv_deaths[,-c(1:2)], 1, quantile, prob = 0.975)),
                   cbind(out5$timeseries$total_hbv_deaths_rate[,c(1:2)],
                         median= apply(out5$timeseries$total_hbv_deaths[,-c(1:2)], 1, median),
                         cri_lower = apply(out5$timeseries$total_hbv_deaths[,-c(1:2)], 1, quantile, prob = 0.025),
                         cri_upper = apply(out5$timeseries$total_hbv_deaths[,-c(1:2)], 1, quantile, prob = 0.975)),
                   cbind(out9$timeseries$total_hbv_deaths_rate[,c(1:2)],
                         median = apply(out9$timeseries$total_hbv_deaths[,-c(1:2)], 1, median),
                         cri_lower = apply(out9$timeseries$total_hbv_deaths[,-c(1:2)], 1, quantile, prob = 0.025),
                         cri_upper = apply(out9$timeseries$total_hbv_deaths[,-c(1:2)], 1, quantile, prob = 0.975)))


# HBV-related deaths
ggplot(hbv_deaths_rate[hbv_deaths_rate$scenario %in% c("status_quo", "monit_0_screen_5"),]) +
  geom_line(aes(x=time, y = median*10000/0.5, group = scenario, colour = scenario), size =1)+
  geom_ribbon(aes(x=time, ymin=cri_lower*10000/0.5, ymax=cri_upper*10000/0.5, group = scenario, colour = scenario, fill = scenario),
              linetype = "dashed", alpha = 0)+
  labs(title = "HBV-related deaths",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "Status quo: historical infant vaccine coverage since 1990 and maintaining 93% coverage after 2018\n
       One-off screen+treat: infant vaccine and screening once in 2020\n
       Repeat screen+treat: Infant vaccine and screening every 5 years starting in 2020") +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(2015,2050)) +
  scale_color_manual(limits = c("status_quo", "screen_2020_monit_0", "screen_2020_monit_5", "monit_0_screen_5"),
                     labels = c("status_quo"="Status quo",
                                "screen_2020_monit_0"="One-off screen+treat",
                                "screen_2020_monit_5"="Monitor every 5 years",
                                "monit_0_screen_5" = "Screen every 5 years"),
                     values = c("screen_2020_monit_0"="steelblue", "status_quo"="orange", "monit_0_screen_5"="deeppink", "screen_2020_monit_5" = "black")) +
  #scale_fill_manual(limits = c("status_quo", "screen_once", "screen_5"),
  #                  labels = c("status_quo"="Status quo",
  #                             "screen_once"="One-off screen+treat",
  #                             "screen_5"="Repeat screen+treat"),
  #                  values = c("screen_once"="steelblue", "status_quo"="orange", "screen_5"="deeppink")) +
  ylab("HBV-related death rate  per 10,000 person-years")+
  xlab("Year")+
  ylim(0,5) +
  theme_classic()

ggplot(hbv_deaths) +
  geom_line(aes(x=time, y = median/0.5, group = scenario, colour = scenario), size =1)+
  geom_ribbon(aes(x=time, ymin=cri_lower/0.5, ymax=cri_upper/0.5, group = scenario, colour = scenario, fill = scenario),
              linetype = "dashed", alpha = 0.05)+
  labs(title = "HBV-related deaths",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "Status quo: historical infant vaccine coverage since 1990 and maintaining 93% coverage after 2018\n
       One-off screen+treat: infant vaccine and screening once in 2020\n
       Repeat screen+treat: Infant vaccine and screening every 5 years starting in 2020") +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(2015,2100)) +
  ylab("Annual number of HBV-related deaths")+
  xlab("Year")+
  theme_classic()

# Reduction in HBV related mortality rate TO DO
quantile((out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2019.5,-c(1,2)]-
  out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2030,-c(1,2)])/
  (out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2019.5,-c(1,2)]), prob =
    c(0.025,0.5,0.975))
quantile((out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2020,-c(1,2)]-
            out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2050,-c(1,2)])/
           (out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2020,-c(1,2)]), prob =
           c(0.025,0.5,0.975))

quantile((out3$timeseries$total_hbv_deaths_rate[out3$timeseries$total_hbv_deaths_rate$time == 2019.5,-c(1,2)]-
            out3$timeseries$total_hbv_deaths_rate[out3$timeseries$total_hbv_deaths_rate$time == 2030,-c(1,2)])/
           (out3$timeseries$total_hbv_deaths_rate[out3$timeseries$total_hbv_deaths_rate$time == 2019.5,-c(1,2)]), prob =
           c(0.025,0.5,0.975))
quantile((out3$timeseries$total_hbv_deaths_rate[out3$timeseries$total_hbv_deaths_rate$time == 2020,-c(1,2)]-
            out3$timeseries$total_hbv_deaths_rate[out3$timeseries$total_hbv_deaths_rate$time == 2050,-c(1,2)])/
           (out3$timeseries$total_hbv_deaths_rate[out3$timeseries$total_hbv_deaths_rate$time == 2020,-c(1,2)]), prob =
           c(0.025,0.5,0.975))

quantile((out5$timeseries$total_hbv_deaths_rate[out5$timeseries$total_hbv_deaths_rate$time == 2020,-c(1,2)]-
            out5$timeseries$total_hbv_deaths_rate[out5$timeseries$total_hbv_deaths_rate$time == 2030,-c(1,2)])/
           (out5$timeseries$total_hbv_deaths_rate[out5$timeseries$total_hbv_deaths_rate$time == 2020,-c(1,2)]), prob =
           c(0.025,0.5,0.975))
quantile((out5$timeseries$total_hbv_deaths_rate[out5$timeseries$total_hbv_deaths_rate$time == 2020,-c(1,2)]-
            out5$timeseries$total_hbv_deaths_rate[out5$timeseries$total_hbv_deaths_rate$time == 2050,-c(1,2)])/
           (out5$timeseries$total_hbv_deaths_rate[out5$timeseries$total_hbv_deaths_rate$time == 2020,-c(1,2)]), prob =
           c(0.025,0.5,0.975))


## REDUCTIONS FOR ELIMINATION ----

red_inf <-
  (out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time == 2015,-c(1:2)]-
     out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time == 2030,-c(1:2)])/
  out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time == 2015,-c(1:2)]

quantile(red_inf, prob = c(0.025,0.5,0.975))

red_inf_rate <-
  (out2$timeseries$total_chronic_infections_rate[out2$timeseries$total_chronic_infections_rate$time == 2015,-c(1:2)]-
     out2$timeseries$total_chronic_infections_rate[out2$timeseries$total_chronic_infections_rate$time == 2030,-c(1:2)])/
  out2$timeseries$total_chronic_infections_rate[out2$timeseries$total_chronic_infections_rate$time == 2015,-c(1:2)]

quantile(red_inf_rate, prob = c(0.025,0.5,0.975))

red_mort <-
  (out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time == 2015,-c(1:2)]-
     out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time == 2030,-c(1:2)])/
  out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time == 2015,-c(1:2)]

quantile(red_mort, prob = c(0.025,0.5,0.975))

red_mort_rate <-
  (out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2015,-c(1:2)]-
     out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2030,-c(1:2)])/
  out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2015,-c(1:2)]

quantile(red_mort_rate, prob = c(0.025,0.5,0.975))

## PREVIOUS ATTEMPTS ----
### Monitoring analysis: Plot population outcomes (y) by healthcare interactions (x) ----

# COUNTERFACTUAL = STATUS QUO

# OUTCOME = HBV DEATHS AVERTED
# EXPOSURE = MEDIAN INCREMENTAL INTERACTIONS

deaths_averted_sq_summary <- rbind(calculate_number_averted(out2$cum_hbv_deaths_2030,
                                                            out3$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2030)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2050,
                                                            out3$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2050)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2100,
                                                            out3$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2100)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2030,
                                                            out4$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2030)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2050,
                                                            out4$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2050)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2100,
                                                            out4$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2100)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2030,
                                                            out5$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2030)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2050,
                                                            out5$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2050)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2100,
                                                            out5$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2100)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2030,
                                                            out6$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2030)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2050,
                                                            out6$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2050)]]),
                                   calculate_number_averted(out2$cum_hbv_deaths_2100,
                                                            out6$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==2100)]]))

deaths_averted_sq_summary <- deaths_averted_sq_summary[deaths_averted_sq_summary$type == "number_averted",]

median_interactions_sq <- rbind(data.frame(by_year = 2030,
                                           scenario = c("screen_2020_monit_0", "screen_2020_monit_10", "screen_2020_monit_5", "screen_2020_monit_1"),
                                           treatments = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_treated[,-c(1:3)]))),
                                           assessments = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_assessed[,-c(1:3)]))),
                                           interactions = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2030)]]$total_interactions[,-c(1:3)])))),
                                data.frame(by_year = 2050,
                                           scenario = c("screen_2020_monit_0", "screen_2020_monit_10", "screen_2020_monit_5", "screen_2020_monit_1"),
                                           treatments = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_treated[,-c(1:3)]))),
                                           assessments = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_assessed[,-c(1:3)]))),
                                           interactions = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2050)]]$total_interactions[,-c(1:3)])))),
                                data.frame(by_year = 2100,
                                           scenario = c("screen_2020_monit_0", "screen_2020_monit_10", "screen_2020_monit_5", "screen_2020_monit_1"),
                                           treatments = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_treated[,-c(1:3)])),
                                                          median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_treated[,-c(1:3)]))),
                                           assessments = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_assessed[,-c(1:3)])),
                                                           median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_assessed[,-c(1:3)]))),
                                           interactions = c(median(unlist(out3$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out4$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out5$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_interactions[,-c(1:3)])),
                                                            median(unlist(out6$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_interactions[,-c(1:3)]))))
)

deaths_averted_by_interactions_sq <- left_join(deaths_averted_sq_summary,
                                               median_interactions_sq, by = c("by_year", "scenario"))


# Plots

# Note these don't show the uncertainty in the number of healthcare interactions,
# which is also substantial (overlapping between no monitoring, every 10 years and every 5 years)

# Healthcare interactions
# 2030
ggplot(data = deaths_averted_by_interactions_sq[deaths_averted_by_interactions_sq$by_year == 2030,],
       aes(x = interactions/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("HBV-related deaths averted/1000") +
  xlab("Median incremental healthcare interactions/1000") +
  labs(title = "By 2030, compared to status quo (no treatment)") +
  ylim(0,3.5) +
  xlim(0,1000)+
  coord_flip()

# 2050
ggplot(data = deaths_averted_by_interactions_sq[deaths_averted_by_interactions_sq$by_year == 2050,],
       aes(x = interactions/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("HBV-related deaths averted/1000") +
  xlab("Median incremental healthcare interactions/1000") +
  labs(title = "By 2050, compared to status quo (no treatment)") +
  ylim(0,7.3) +
  xlim(0,1400)+
  coord_flip()

#2100
ggplot(data = deaths_averted_by_interactions_sq[deaths_averted_by_interactions_sq$by_year == 2100,],
       aes(x = interactions/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("HBV-related deaths averted/1000") +
  xlab("Median incremental healthcare interactions/1000") +
  labs(title = "By 2100, compared to status quo (no treatment)") +
  ylim(0,8.6) +
  xlim(0,1600)+
  coord_flip()

# Treatment initiations
# 2030
ggplot(data = deaths_averted_by_interactions_sq[deaths_averted_by_interactions_sq$by_year == 2030,],
       aes(x = treatments/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("HBV-related deaths averted/1000") +
  xlab("Median incremental treatment initiations/1000") +
  labs(title = "By 2030, compared to status quo (no treatment)") +
  ylim(0,3.1) +
  xlim(0,10)+
  coord_flip()

# 2050
ggplot(data = deaths_averted_by_interactions_sq[deaths_averted_by_interactions_sq$by_year == 2050,],
       aes(x = treatments/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("HBV-related deaths averted/1000") +
  xlab("Median incremental treatment initiations/1000") +
  labs(title = "By 2050, compared to status quo (no treatment)") +
  ylim(0,7) +
  xlim(0,10)+
  coord_flip()

# Clinical assessments
# 2030
ggplot(data = deaths_averted_by_interactions_sq[deaths_averted_by_interactions_sq$by_year == 2030,],
       aes(x = assessments/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("HBV-related deaths averted/1000") +
  xlab("Median incremental clinical assessments/1000") +
  labs(title = "By 2030, compared to status quo (no treatment)") +
  ylim(0,3.1) +
  xlim(0,400)+
  coord_flip()

# 2050
ggplot(data = deaths_averted_by_interactions_sq[deaths_averted_by_interactions_sq$by_year == 2050,],
       aes(x = assessments/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("HBV-related deaths averted/1000") +
  xlab("Median incremental clinical assessments/1000") +
  labs(title = "By 2030, compared to status quo (no treatment)") +
  ylim(0,7) +
  xlim(0,900)+
  coord_flip()

# OUTCOME = LIFE YEARS SAVED

ly_gained_sq_summary <- rbind(calculate_number_averted(out3$ly[[which(seq(2025,2100, by = 5)==2030)]],
                                                       out2$ly_2030),
                              calculate_number_averted(out3$ly[[which(seq(2025,2100, by = 5)==2050)]],
                                                       out2$ly_2050),
                              calculate_number_averted(out3$ly[[which(seq(2025,2100, by = 5)==2100)]],
                                                       out2$ly_2100),
                              calculate_number_averted(out4$ly[[which(seq(2025,2100, by = 5)==2030)]],
                                                       out2$ly_2030),
                              calculate_number_averted(out4$ly[[which(seq(2025,2100, by = 5)==2050)]],
                                                       out2$ly_2050),
                              calculate_number_averted(out4$ly[[which(seq(2025,2100, by = 5)==2100)]],
                                                       out2$ly_2100),
                              calculate_number_averted(out5$ly[[which(seq(2025,2100, by = 5)==2030)]],
                                                       out2$ly_2030),
                              calculate_number_averted(out5$ly[[which(seq(2025,2100, by = 5)==2050)]],
                                                       out2$ly_2050),
                              calculate_number_averted(out5$ly[[which(seq(2025,2100, by = 5)==2100)]],
                                                       out2$ly_2100),
                              calculate_number_averted(out6$ly[[which(seq(2025,2100, by = 5)==2030)]],
                                                       out2$ly_2030),
                              calculate_number_averted(out6$ly[[which(seq(2025,2100, by = 5)==2050)]],
                                                       out2$ly_2050),
                              calculate_number_averted(out6$ly[[which(seq(2025,2100, by = 5)==2100)]],
                                                       out2$ly_2100))

ly_gained_sq_summary$scenario <- ly_gained_sq_summary$counterfactual

ly_gained_by_interactions_sq <- left_join(ly_gained_sq_summary,
                                          median_interactions_sq, by = c("by_year", "scenario"))
ly_gained_by_interactions_sq <- ly_gained_by_interactions_sq[ly_gained_by_interactions_sq$type == "number_averted",]


# Healthcare interactions
# 2030
ggplot(data = ly_gained_by_interactions_sq[ly_gained_by_interactions_sq$by_year == 2030,],
       aes(x = interactions/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("Life-years gained/1000") +
  xlab("Median incremental healthcare interactions/1000") +
  labs(title = "By 2030, compared to status quo (no treatment)") +
  ylim(0,17) +
  xlim(0,1000)+
  coord_flip()

# 2050
ggplot(data = ly_gained_by_interactions_sq[ly_gained_by_interactions_sq$by_year == 2050,],
       aes(x = interactions/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("Life-years gained/1000") +
  xlab("Median incremental healthcare interactions/1000") +
  labs(title = "By 2050, compared to status quo (no treatment)") +
  ylim(0,120) +
  xlim(0,1500)+
  coord_flip()

#2100
ggplot(data = ly_gained_by_interactions_sq[ly_gained_by_interactions_sq$by_year == 2100,],
       aes(x = interactions/1000)) +
  geom_point(aes(y = median/1000, colour = scenario), size = 5) +
  geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000, colour = scenario)) +
  ylab("Life-years gained/1000") +
  xlab("Median incremental healthcare interactions/1000") +
  labs(title = "By 2100, compared to status quo (no treatment)") +
  ylim(0,300) +
  xlim(0,1600)+
  coord_flip()

# COUNTERFACTUAL = ONE-OFF SCREENING BUT NO MONITORING (out3)

# OUTCOME = HBV DEATHS AVERTED



