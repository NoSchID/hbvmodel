# What is the optimal monitoring frequency for different age groups?
# Exploring this from a cohort perspective, assuming 15-60 year olds have been screened (A1)

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Function to plot boxplot whiskers as 95% percentile
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# Load files (A1) ----

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

# Status quo
out1 <- readRDS(paste0(out_path, "a1_out1_status_quo_cohort_240920.rds"))
out1 <- out1[[1]]
out2 <- readRDS(paste0(out_path, "out2_status_quo_180820.rds"))
out2 <- out2[[1]]

# No monitoring
out3 <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_201020.rds"))
out3 <- out3[[1]]

# Monitoring all age groups
out4 <- readRDS(paste0(out_path, "a1_out4_screen_2020_monit_10_290920.rds"))
out4 <- out4[[1]]
out5 <- readRDS(paste0(out_path, "a1_out5_screen_2020_monit_5_201020.rds"))
out5 <- out5[[1]]
out5a <- readRDS(paste0(out_path, "a1_out5a_screen_2020_monit_6_290920.rds"))
out5a <- out5a[[1]]
out5b <- readRDS(paste0(out_path, "a1_out5b_screen_2020_monit_7_290920.rds"))
out5b <- out5b[[1]]
out5c <- readRDS(paste0(out_path, "a1_out5c_screen_2020_monit_8_290920.rds"))
out5c <- out5c[[1]]
out5d <- readRDS(paste0(out_path, "a1_out5d_screen_2020_monit_9_290920.rds"))
out5d <- out5d[[1]]
out6 <- readRDS(paste0(out_path, "a1_out6_screen_2020_monit_1_201020.rds"))
out6 <- out6[[1]]
out6a <- readRDS(paste0(out_path, "a1_out6a_screen_2020_monit_2_280920.rds"))
out6a <- out6a[[1]]
out6b <- readRDS(paste0(out_path, "a1_out6b_screen_2020_monit_3_280920.rds"))
out6b <- out6b[[1]]
out6c <- readRDS(paste0(out_path, "a1_out6c_screen_2020_monit_4_280920.rds"))
out6c <- out6c[[1]]

# Monitoring different age groups until they age out (yearly)
monit_out1 <- readRDS(paste0(out_path, "a1_monit_out1_201020.rds"))
monit_out1 <- monit_out1[[1]]
monit_out2 <- readRDS(paste0(out_path, "a1_monit_out2_201020.rds"))
monit_out2 <- monit_out2[[1]]
monit_out3 <- readRDS(paste0(out_path, "a1_monit_out3_240920.rds"))
monit_out3 <- monit_out3[[1]]
monit_out4 <- readRDS(paste0(out_path, "a1_monit_out4_240920.rds"))
monit_out4 <- monit_out4[[1]]
monit_out5 <- readRDS(paste0(out_path, "a1_monit_out5_240920.rds"))
monit_out5 <- monit_out5[[1]]

# Monitoring different age groups until they age out (every 5 years)
monit_out6 <- readRDS(paste0(out_path, "a1_monit_out6_201020.rds"))
monit_out6 <- monit_out6[[1]]
monit_out7 <- readRDS(paste0(out_path, "a1_monit_out7_201020.rds"))
monit_out7 <- monit_out7[[1]]
monit_out8 <- readRDS(paste0(out_path, "a1_monit_out8_240920.rds"))
monit_out8 <- monit_out8[[1]]
monit_out9 <- readRDS(paste0(out_path, "a1_monit_out9_240920.rds"))
monit_out9 <- monit_out9[[1]]
monit_out10 <- readRDS(paste0(out_path, "a1_monit_out10_240920.rds"))
monit_out10 <- monit_out10[[1]]

# Labels of age groups being monitored
# sim1, sim6 = 15-30
# sim2, sim7 = 15-45
# sim3, sim8 = 30-45
# sim4, sim 9 = 30+
# sim5, sim10 = 45+

# Monitoring different age cohorts until they die (every year)
# A2 = screen and treat 45-60 year olds. A3 = screen and treat 30-60 year olds
a2_out3 <- readRDS(paste0(out_path, "a2_out3_screen_2020_monit_0_161020.rds"))
a2_out3 <- a2_out3[[1]]
a2_out5 <- readRDS(paste0(out_path, "a2_out5_screen_2020_monit_5_161020.rds"))
a2_out5 <- a2_out5[[1]]
a2_out6 <- readRDS(paste0(out_path, "a2_out6_screen_2020_monit_1_161020.rds"))
a2_out6 <- a2_out6[[1]]
a3_out3 <- readRDS(paste0(out_path, "a3_out3_screen_2020_monit_0_161020.rds"))
a3_out3 <- a3_out3[[1]]
a3_out5 <- readRDS(paste0(out_path, "a3_out5_screen_2020_monit_5_161020.rds"))
a3_out5 <- a3_out5[[1]]
a3_out6 <- readRDS(paste0(out_path, "a3_out6_screen_2020_monit_1_161020.rds"))
a3_out6 <- a3_out6[[1]]

scenario_labels <- list("No treatment" = "status_quo_cohort",
                        "No monitoring" = "screen_2020_monit_0",
                        "5-yearly all ages"="screen_2020_monit_5",
                        "Yearly all ages"="screen_2020_monit_1",
                        "Yearly 15-30"="screen_2020_monit_sim1",
                        "Yearly 15-45"="screen_2020_monit_sim2",
                        "Yearly 30-45"="screen_2020_monit_sim3",
                        "Yearly 30+"="screen_2020_monit_sim4",
                        "Yearly 45+"="screen_2020_monit_sim5",
                        "5-yearly 15-30"="screen_2020_monit_sim6",
                        "5-yearly 15-45"= "screen_2020_monit_sim7",
                        "5-yearly 30-45"="screen_2020_monit_sim8",
                        "5-yearly 30+"="screen_2020_monit_sim9",
                        "5-yearly 45+"="screen_2020_monit_sim10")

# Subset yearly frequencies
sub_yearly <- names(scenario_labels)[4:9]
# Subset 5-yearly frequencies
sub_5yearly <- names(scenario_labels)[c(3,10:14)]
# Selection for plot
sub_mixed <- names(scenario_labels)[c(3,4,5,8,10,13)]
# Since the best are 5-yearly 45+, 5-yearly 30+ and 5-yearly all ages, chose here <30, >30 and all

# Age group comparisons
sub_age_groups1 <- names(scenario_labels)[c(5,7,9,10,12,14)]
# 15-30, 30-45, 45+
sub_age_groups2 <- names(scenario_labels)[c(5,8,10,13)]
# 15-30, 30+
sub_age_groups3 <- names(scenario_labels)[c(6,9,11,14)]
# 15-45, 45+
sub_age_groups_all <- names(scenario_labels)[c(3,4)]

# 1) Comparing monitoring to no treatment ----
# Average age at death ----

# Summary of average age at death for final table
cohort_age_at_death <- data.frame(rbind(out1$cohort_age_at_death,
                                        out3$cohort_age_at_death,
                                        out5$cohort_age_at_death,
                                        out6$cohort_age_at_death,
                                        monit_out1$cohort_age_at_death,
                                        monit_out2$cohort_age_at_death,
                                        monit_out3$cohort_age_at_death,
                                        monit_out4$cohort_age_at_death,
                                        monit_out5$cohort_age_at_death,
                                        monit_out6$cohort_age_at_death,
                                        monit_out7$cohort_age_at_death,
                                        monit_out8$cohort_age_at_death,
                                        monit_out9$cohort_age_at_death,
                                        monit_out10$cohort_age_at_death))
levels(cohort_age_at_death$scenario) <- scenario_labels

cohort_age_at_death_long <- gather(cohort_age_at_death, key = "sim", value = "value", -scenario)

ggplot(cohort_age_at_death_long) +
  geom_boxplot(aes(x = reorder(scenario, value), y = value)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# monit_sim4 nearly as good as monitoring all yearly - this is 30+ year olds, yearly,
# followed by monit_sim5 (45+ year olds, yearly)
# worst are monit_sim6 and monit_sim1 (15-30, 5-yearly and yearly)
# But note differences are very small.

ggplot(subset(cohort_age_at_death_long, scenario %in%
                c("No treatment", "No monitoring", "Yearly all ages", "5-yearly all ages")),
       aes(x = reorder(scenario, value), y = value)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  scale_x_discrete("Scenario", labels = c("5-yearly all ages" = "Monitor\nevery 5 years",
                                          "Yearly all ages" = "Monitor\nevery year")) +
  ylab("Average age at death (years)\nin the screened HBV carrier cohort") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))

# Deaths averted ----
cohort_deaths_averted_sq_long <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out1,
                                 scenario_objects = list(out3,
                                                         out5,
                                                         out6),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "no treatment")

ggplot(cohort_deaths_averted_sq_long[cohort_deaths_averted_sq_long$type == "proportion_averted",],
       aes(scenario, value*100)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  ylab("Percentage of HBV-related deaths averted\nin the cohort") +
  scale_x_discrete("Monitoring frequency", labels =
                     c("screen_2020_monit_0" = "No monitoring",
                       "screen_2020_monit_5" = "Monitor\nevery 5 years",
                       "screen_2020_monit_1" = "Monitor\nevery year")) +
  theme_bw() +
  ylim(0,100) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))

# LY saved ----
cohort_ly_gained_sq_long <-
  plot_ly_gained_cohort(counterfactual_object = out1,
                        scenario_objects = list(out3,
                                                out5,
                                                out6),
                                 counterfactual_label = "no treatment")

ggplot(cohort_ly_gained_sq_long[cohort_ly_gained_sq_long$type == "proportion_averted",],
       aes(counterfactual, value*100)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  ylab("Percentage of life-years saved\nin the cohort") +
  scale_x_discrete("Monitoring frequency", labels =
                     c("screen_2020_monit_0" = "No monitoring",
                       "screen_2020_monit_5" = "Monitor\nevery 5 years",
                       "screen_2020_monit_1" = "Monitor\nevery year")) +
  theme_bw() +
  ylim(0,15) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))

# Deaths averted per clinical assessments and interactions (monitor all ages) ----
deaths_averted_per_assessment_sq_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(out3, out5, out6),
                                                     interaction_type = "total_assessed",
                                                     counterfactual_label = "no treatment")
#levels(deaths_averted_per_assessment_sq_long$scenario) <- scenario_labels
#deaths_averted_per_assessment_long$freq <- "Every year"
#deaths_averted_per_assessment_long$freq[deaths_averted_per_assessment_long$scenario %in% sub_5yearly] <-
#  "Every 5 years"

ggplot(data =deaths_averted_per_assessment_sq_long,
       aes(x=scenario, y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  ylab("HBV-related deaths averted\nper 10,000 treatment eligibility assessments") +
  xlab("Monitoring frequency") +
  scale_x_discrete(labels = c("screen_2020_monit_0" = "No monitoring",
                              "screen_2020_monit_5" = "Monitor\nevery 5 years",
                              "screen_2020_monit_1" = "Monitor\nevery year")) +
#  ylim(0,100) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))

deaths_averted_per_interactions_sq_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(out3, out5, out6),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "no treatment")
#levels(deaths_averted_per_assessment_sq_long$scenario) <- scenario_labels
#deaths_averted_per_assessment_long$freq <- "Every year"
#deaths_averted_per_assessment_long$freq[deaths_averted_per_assessment_long$scenario %in% sub_5yearly] <-
#  "Every 5 years"

ggplot(data =deaths_averted_per_interactions_sq_long,
       aes(x=scenario, y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  ylab("HBV-related deaths averted\nper 10,000 clinical interactions") +
  xlab("Monitoring frequency") +
  scale_x_discrete(labels = c("screen_2020_monit_0" = "No monitoring",
                              "screen_2020_monit_5" = "Monitor\nevery 5 years",
                              "screen_2020_monit_1" = "Monitor\nevery year")) +
  ylim(0,70) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))



# Barchart of total median interactions with/without monitoring ----
# Over the lifetime of the cohort = by 2100
interactions_bar_df <- data.frame(rbind(
  cbind(out3$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_screened,
        interaction_type = "total_screened"),
  cbind(out3$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_assessed,
        interaction_type = "total_assessed"),
  cbind(out3$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_treated,
        interaction_type = "total_treated"),
  cbind(out5$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_screened,
        interaction_type = "total_screened"),
  cbind(out5$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_assessed,
        interaction_type = "total_assessed"),
  cbind(out5$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_treated,
        interaction_type = "total_treated"),
  cbind(out6$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_screened,
        interaction_type = "total_screened"),
  cbind(out6$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_assessed,
        interaction_type = "total_assessed"),
  cbind(out6$interactions[[which(seq(2025,2100, by = 5)==2100)]]$total_treated,
        interaction_type = "total_treated")))

interactions_bar_df <- interactions_bar_df %>%
  gather(key = "sim", value = "value", -from_year, -by_year, -scenario, -interaction_type) %>%
  group_by(from_year, by_year, scenario, interaction_type) %>%
  summarise(median = median(value))

library(viridis)
ggplot(interactions_bar_df) +
  geom_bar(aes(x = scenario, y = median, fill = interaction_type), position="stack", width = 0.5, stat = "identity") +
  scale_fill_viridis(discrete = TRUE, direction =-1,
                     labels = c("total_screened" = "HBsAg tests",
                                "total_assessed" = "Treatment eligibility\nassessments",
                                "total_treated" = "Treatment initiations")) +
  xlab("Monitoring frequency") +
  scale_x_discrete(labels = c("screen_2020_monit_0" = "No monitoring",
                              "screen_2020_monit_5" = "Monitor\nevery 5 years",
                              "screen_2020_monit_1" = "Monitor\nevery year")) +
  ylab("Total clinical interactions 2020-2100") +
  labs(fill = "Resource type") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))

# 2) Yearly monitoring of which age group is most effective, compared to no monitoring? ----
# Average age at death ----
# Extension in average age at death compared to no monitoring
age_at_death_ext <- data.frame(rbind(cohort_age_at_death[3,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[4,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[5,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[6,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[7,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[8,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[9,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[10,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[11,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[12,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[13,]-cohort_age_at_death[2,],
                                     cohort_age_at_death[14,]-cohort_age_at_death[2,]))
age_at_death_ext$scenario <- cohort_age_at_death$scenario[-c(1:2)]

age_at_death_ext_long <- gather(age_at_death_ext, key = "sim", value = "value", -scenario)
age_at_death_ext_long$freq <- "Every year"
age_at_death_ext_long$freq[age_at_death_ext_long$scenario %in% sub_5yearly] <-
  "Every 5 years"

ggplot(age_at_death_ext_long, aes(x = reorder(scenario, value), y = value)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Extension in average age at death (years)") +
  xlab("Monitoring frequency and age group") +
  labs(title= "Cohort impact compared to treatment programme without monitoring") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))
# Extension in age is nearly the same for yearly monitoring in 30+ year olds
# as in all ages (best). This is followed by the same for 5-yearly monitoring in these age groups.
# This means 5-yearly monitoring in 30+ year olds is better than yearly monitoring in other age groups.
# Monitoring only 15-30 year olds is the worst strategy and barely worth it for this outcome.

# Subset plot for yearly frequencies *PLOT*
ggplot(subset(age_at_death_ext_long, scenario %in% sub_yearly),
       aes(x = reorder(gsub("Yearly ","", scenario), value), y = value, fill = scenario)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Extension in average age at death\nof cohort (years)") +
  xlab("Age group") +
  labs(title= "Yearly monitoring compared to no monitoring") +
  scale_x_discrete(labels = c("all ages" = "15+ (all)")) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Age group 1 plot
age_at_death_ext_long2 <- subset(age_at_death_ext_long,
                                                scenario %in% sub_age_groups1)
age_at_death_ext_long2$scenario <- gsub("5-yearly ", "",
                                                     gsub("Yearly ","",
                                                          age_at_death_ext_long2$scenario))

ggplot(age_at_death_ext_long2,
       aes(x = scenario, y = value, fill = freq)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge2", width = 0.5) +
  ylab("Extension in average age at death\nof cohort (years)") +
  xlab("Age group monitored (years)") +
  labs(fill = "Monitoring\nfrequency", title = paste0("Effect of monitoring compared to no monitoring")) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))

# Maybe check extension in age per 10,000 assessments

# Deaths averted ----
# Compare cohort number of HBV deaths averted compared to no monitoring
cohort_deaths_averted_long <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
                                 scenario_objects = list(out5,
                                                         out6,
                                                         monit_out1,
                                                         monit_out2,
                                                         monit_out3,
                                                         monit_out4,
                                                         monit_out5,
                                                         monit_out6,
                                                         monit_out7,
                                                         monit_out8,
                                                         monit_out9,
                                                         monit_out10),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")
levels(cohort_deaths_averted_long$scenario) <- scenario_labels
cohort_deaths_averted_long$freq <- "Every year"
cohort_deaths_averted_long$freq[cohort_deaths_averted_long$scenario %in% sub_5yearly] <-
  "Every 5 years"


p1 <- ggplot(cohort_deaths_averted_long[cohort_deaths_averted_long$type == "number_averted",],
       aes(reorder(scenario, value), value, fill = scenario)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Number of HBV-related deaths averted") +
  labs(title = paste0("Cohort impact compared to counterfactual:\n", "no monitoring")) +
  xlab("Monitoring frequency") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")
# Again highest for yearly all ages and nearly the same for yearly 30+. For HBV deaths,
# now followed by yearly 45+ before 5-yearly all ages and 5-yearly 30+.
# Again very small number for 15-30 year olds.

# Subset plot for yearly monitoring *PLOT*
ggplot(subset(cohort_deaths_averted_long, type == "number_averted" & scenario %in% sub_yearly),
       aes(reorder(gsub("Yearly ","", scenario), value), value, fill = scenario)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Number of HBV-related deaths averted") +
  labs(title = paste0("Yearly monitoring compared to no monitoring")) +
  scale_x_discrete(labels = c("all ages" = "15+ (all)")) +
  scale_fill_brewer(palette = "Set2") +
  xlab("Age group") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

ggplot(cohort_deaths_averted_long[cohort_deaths_averted_long$type == "proportion_averted",],
       aes(reorder(scenario, value), value)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Proportion of HBV-related deaths averted") +
  labs(title = paste0("Cohort impact compared to counterfactual:\n", "no monitoring")) +
  xlab("Monitoring frequency") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))
# Same as number but less uncertainty

# Age group 1 plot
cohort_deaths_averted_long2 <- subset(cohort_deaths_averted_long,
                                 scenario %in% sub_age_groups1)
cohort_deaths_averted_long2$scenario <- gsub("5-yearly ", "",
                                        gsub("Yearly ","",
                                             cohort_deaths_averted_long2$scenario))

ggplot(subset(cohort_deaths_averted_long2,  type == "number_averted"),
       aes(x = scenario, y = value, fill = freq)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge2", width = 0.5) +
  ylab("Number of HBV-related deaths averted") +
  xlab("Age group monitored (years)") +
  labs(fill = "Monitoring\nfrequency", title = paste0("Effect of monitoring compared to no monitoring")) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))

# LY saved ----
# Compare cohort number of life years gained compared to no monitoring
cohort_ly_gained_long <-
  plot_ly_gained_cohort(counterfactual_object = out3,
                        scenario_objects = list(out5,
                                                out6,
                                                monit_out1,
                                                monit_out2,
                                                monit_out3,
                                                monit_out4,
                                                monit_out5,
                                                monit_out6,
                                                monit_out7,
                                                monit_out8,
                                                monit_out9,
                                                monit_out10),
                        counterfactual_label = "treatment programme without monitoring")

levels(cohort_ly_gained_long$counterfactual) <- scenario_labels
cohort_ly_gained_long$freq <- "Every year"
cohort_ly_gained_long$freq[cohort_ly_gained_long$counterfactual %in% sub_5yearly] <-
  "Every 5 years"


ggplot(cohort_ly_gained_long[cohort_ly_gained_long$type == "number_averted",],
       aes(reorder(counterfactual, value), value)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Number of life-years saved") +
  labs(title = paste0("Cohort impact compared to counterfactual:\n", "no monitoring")) +
  xlab("Monitoring frequency") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        title = element_text(size = 15))
# Now the order is more similar to that shown with average age at death
# 5-yearly in 30+ year olds slightly higher than yearly in 45+ year olds
# (less extension of life possible among the older age groups only).
# Nevertheless, number of life-years saved again very low when monitoring only 15-30 year olds.

# Subset plot for yearly monitoring *PLOT*
ggplot(subset(cohort_ly_gained_long,type == "number_averted" & counterfactual %in% sub_yearly),
       aes(reorder(gsub("Yearly ","", counterfactual), value), value, fill = counterfactual)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Number of life-years saved") +
  labs(title = paste0("Yearly monitoring compared to no monitoring")) +
  scale_x_discrete(labels = c("all ages" = "15+ (all)")) +
  xlab("Age group") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Age group 1 plot
cohort_ly_gained_long2 <- subset(cohort_ly_gained_long,
                                      counterfactual %in% sub_age_groups1)
cohort_ly_gained_long2$counterfactual <- gsub("5-yearly ", "",
                                             gsub("Yearly ","",
                                                  cohort_ly_gained_long2$counterfactual))

ggplot(subset(cohort_ly_gained_long2,  type == "number_averted"),
       aes(x = counterfactual, y = value, fill = freq)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge2", width = 0.5) +
  ylab("Number of life-years saved") +
  xlab("Age group monitored (years)") +
  labs(fill = "Monitoring\nfrequency", title = paste0("Effect of monitoring compared to no monitoring")) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))

# Deaths averted per resource use ----
# HBV DEATHS AVERTED PER INCREMENTAL ASSESSMENT
deaths_averted_per_assessment_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects = list(out5,
                                                                             out6,
                                                                             monit_out1,
                                                                             monit_out2,
                                                                             monit_out3,
                                                                             monit_out4,
                                                                             monit_out5,
                                                                             monit_out6,
                                                                             monit_out7,
                                                                             monit_out8,
                                                                             monit_out9,
                                                                             monit_out10),
                                                     interaction_type = "total_assessed",
                                                     counterfactual_label = "treatment programme without monitoring")
levels(deaths_averted_per_assessment_long$scenario) <- scenario_labels
deaths_averted_per_assessment_long$freq <- "Every year"
deaths_averted_per_assessment_long$freq[deaths_averted_per_assessment_long$scenario %in% sub_5yearly] <-
  "Every 5 years"

timepoints <- c(2030,2050,2100)
period_labs <- c(paste0("2020-",timepoints[1]), paste0("2020-",timepoints[2]), paste0("2020-",timepoints[3]))
names(period_labs) <- c(as.character(timepoints[1]), as.character(timepoints[2]),
                        as.character(timepoints[3]))


p2 <- ggplot(data = subset(deaths_averted_per_assessment_long, by_year == 2100),
       aes(x=reorder(scenario,value), y=value*10000, fill = scenario)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(~by_year, ncol = 3, labeller=labeller(by_year = period_labs),scales = "free_y") +
#  xlab(x_axis_label) +
  ylab("HBV-related deaths averted per 10,000 assessments") +
#  labs(title = paste0("Population impact compared to counterfactual:\n", counterfactual_label)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle =90, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")
# The pattern for number of deaths averted PER 10,000 monitoring interactions is interestingly
# the opposite of the overall impact. Regardless of age group, 5-yearly monitoring here is
# always more resource-effective than yearly. Monitoring in <45 year olds is now more
# attractive than in 45+ year olds. Yearly monitoring in all ages and in 30+ year olds is the worst.

# Subset plot for yearly monitoring *PLOT*
ggplot(data = subset(deaths_averted_per_assessment_long, by_year == 2100 &  scenario %in% sub_yearly),
       aes(x=reorder(gsub("Yearly ","", scenario),value), y=value*10000, fill = scenario)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("HBV-related deaths averted\nper 10,000 monitoring assessments") +
  labs(title = paste0("Yearly monitoring compared to no monitoring")) +
  scale_x_discrete(labels = c("all ages" = "15+ (all)")) +
  scale_fill_brewer(palette = "Set2") +
  xlab("Age group") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle =90, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Subset plot by age group comparison
deaths_averted_per_assessment_long2 <- subset(deaths_averted_per_assessment_long, by_year == 2100 &
                                                scenario %in% sub_age_groups1)
deaths_averted_per_assessment_long2$scenario <- gsub("5-yearly ", "",
                                                     gsub("Yearly ","",
                                                          deaths_averted_per_assessment_long2$scenario))

ggplot(data =deaths_averted_per_assessment_long2,
       aes(x=scenario, y=value*10000, fill = freq)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("HBV-related deaths averted\nper 10,000 monitoring assessments") +
  labs(fill = "Monitoring\nfrequency", title = paste0("Effect of monitoring compared to no monitoring")) +
  xlab("Age group monitored (years)") +
  ylim(0,100) +
  theme_bw() +
#  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))

# LY saved per resource use ----
# HBV DEATHS AVERTED PER INCREMENTAL INTERACTION
# The difference here should be assessments + treatment initiations
# Basically identical to assessments only
deaths_averted_per_interaction_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects =list(out5,
                                                                            out6,
                                                                            monit_out1,
                                                                            monit_out2,
                                                                            monit_out3,
                                                                            monit_out4,
                                                                            monit_out5,
                                                                            monit_out6,
                                                                            monit_out7,
                                                                            monit_out8,
                                                                            monit_out9,
                                                                            monit_out10),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "no treatment programme")
levels(deaths_averted_per_interaction_long$scenario) <- scenario_labels

ggplot(data = deaths_averted_per_interaction_long,
       aes(x=reorder(scenario,value), y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(~by_year, ncol = 3, labeller=labeller(by_year = period_labs),scales = "free_y") +
  #  xlab(x_axis_label) +
  ylab("HBV-related deaths averted per 10,000 interactions") +
  #  labs(title = paste0("Population impact compared to counterfactual:\n", counterfactual_label)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle =90, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))

# LIFE YEARS GAINED PER INCREMENTAL HEALTHCARE INTERACTION
ly_gained_per_interaction_long <-
  plot_ly_gained_per_healthcare_interaction(counterfactual_object = out3,
                                            scenario_objects = list(out5,
                                                                    out6,
                                                                    monit_out1,
                                                                    monit_out2,
                                                                    monit_out3,
                                                                    monit_out4,
                                                                    monit_out5,
                                                                    monit_out6,
                                                                    monit_out7,
                                                                    monit_out8,
                                                                    monit_out9,
                                                                    monit_out10),
                                            interaction_type = "total_interactions",
                                            counterfactual_label = "treatment programme without monitoring")

ly_gained_per_assessment_long <-
  plot_ly_gained_per_healthcare_interaction(counterfactual_object = out3,
                                            scenario_objects = list(out5,
                                                                    out6,
                                                                    monit_out1,
                                                                    monit_out2,
                                                                    monit_out3,
                                                                    monit_out4,
                                                                    monit_out5,
                                                                    monit_out6,
                                                                    monit_out7,
                                                                    monit_out8,
                                                                    monit_out9,
                                                                    monit_out10),
                                            interaction_type = "total_assessed",
                                            counterfactual_label = "treatment programme without monitoring")
levels(ly_gained_per_assessment_long$scenario) <- scenario_labels
ly_gained_per_assessment_long$freq <- "Every year"
ly_gained_per_assessment_long$freq[ly_gained_per_assessment_long$scenario %in% sub_5yearly] <-
  "Every 5 years"

timepoints <- c(2030,2050,2100)
period_labs <- c(paste0("2020-",timepoints[1]), paste0("2020-",timepoints[2]), paste0("2020-",timepoints[3]))
names(period_labs) <- c(as.character(timepoints[1]), as.character(timepoints[2]),
                        as.character(timepoints[3]))

ggplot(data = subset(ly_gained_per_assessment_long, by_year == 2100),
       aes(x=reorder(scenario,value), y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(~by_year, ncol = 3, labeller=labeller(by_year = period_labs),scales = "free_y") +
  #  xlab(x_axis_label) +
  ylab("Life-years saved per 10,000 assessments") +
  #  labs(title = paste0("Population impact compared to counterfactual:\n", counterfactual_label)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle =90, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))
# 5-yearly monitoring in 15-30 year olds is best (and yearly in 15-30 year olds is also better
# than 5-yearly in 30+ and 45+ year olds). When looking at life-years (but also at deaths averted),
# in relation to resource use, focus on younger people is more favourable.

# Subset plot for yearly monitoring *PLOT*
ggplot(data = subset(ly_gained_per_assessment_long, by_year == 2100 &  scenario %in% sub_yearly),
       aes(x=reorder(gsub("Yearly ","", scenario),value), y=value*10000, fill = scenario)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Life-years saved\nper 10,000 monitoring assessments") +
  labs(title = paste0("Yearly monitoring compared to no monitoring")) +
  scale_x_discrete(labels = c("all ages" = "15+ (all)")) +
  scale_fill_brewer(palette = "Set2") +
  xlab("Age group") +
  theme_bw() +
#  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle =90, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Subset plot by age group comparison
ly_gained_per_assessment_long2 <- subset(ly_gained_per_assessment_long, by_year == 2100 &
                                                scenario %in% sub_age_groups1)
ly_gained_per_assessment_long2$scenario <- gsub("5-yearly ", "",
                                                     gsub("Yearly ","",
                                                          ly_gained_per_assessment_long2$scenario))

ggplot(data =ly_gained_per_assessment_long2,
       aes(x=scenario, y=value*10000, fill = freq)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  ylab("Life-years saved\nper 10,000 monitoring assessments") +
  labs(fill = "Monitoring\nfrequency", title = paste0("Effect of monitoring compared to no monitoring")) +
  xlab("Age group monitored (years)") +
#  ylim(0,100) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))

# Both for deaths averted and life-years saved per 10,000 assessments,
# focusing on specific age groups is better than yearly or 5-yearly monitoring in all ages.
# 5-yearly monitoring in given age groups is also better than yearly monitoring in the same age group.
# The best strategies are all in <45 year olds: 15-30, 30-45 or 15-45 depending on the outcome of interest.
# These 3 all have similar numbers of deaths averted per assessment, but the life-years saved is
# higher for the 15-30 year olds than other <45 year groups.
# Yet the actual number of deaths averted and life-years saved with this strategy will be minimal!

# What if we compare:
# yearly 30+ vs yearly all ages? => pretty much same per assessment
# 5-yearly all ages vs yearly 45+, yearly 30+ and yearly all ages and 5-yearly 30+
# => 5-yearly 30+ and all ages are pretty much the same per assessment, but
# are far better than the yearly strategies.
# These all have a good amount of deaths averted.

# This result comes from the fact that very few people in the screened cohort are under 30 years old,
# as a result of the lower prevalence in this age group.

grid.arrange(p1,p2, ncol = 1)
# 5-yearly monitoring in all ages, 30+ and 45+ year olds have the same relative positions for
# impact and resource use.

# 3) What is the optimal monitoring frequency if monitoring all age groups? ----
monitoring_freq_allages_label <- list("10"="screen_2020_monit_10",
                                      "9"="screen_2020_monit_9",
                                      "8"="screen_2020_monit_8",
                                      "7"="screen_2020_monit_7",
                                      "6"="screen_2020_monit_6",
                                      "5"="screen_2020_monit_5",
                                        "4"="screen_2020_monit_4",
                                        "3"="screen_2020_monit_3",
                                        "2"="screen_2020_monit_2",
                                        "1"="screen_2020_monit_1")

# Impact per resource use ----

deaths_averted_per_interaction_allages_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects = list(out4,
                                                                             out5,
                                                                             out5a,
                                                                             out5b,
                                                                             out5c,
                                                                             out5d,
                                                                             out6,
                                                                             out6a,
                                                                             out6b,
                                                                             out6c),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "treatment programme without monitoring")
levels(deaths_averted_per_interaction_allages_long$scenario) <- monitoring_freq_allages_label


ggplot(data = subset(deaths_averted_per_interaction_allages_long, by_year == 2100),
       aes(x=scenario, y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "#00BFC4") +
  ylab("HBV-related deaths averted\nper 10,000 incremental clinical interactions") +
  xlab("Monitoring frequency (years)") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

ly_gained_per_interaction_allages_long <-
  plot_ly_gained_per_healthcare_interaction(counterfactual_object = out3,
                                            scenario_objects = list(out4,
                                                                    out5,
                                                                    out5a,
                                                                    out5b,
                                                                    out5c,
                                                                    out5d,
                                                                    out6,
                                                                    out6a,
                                                                    out6b,
                                                                    out6c),
                                            interaction_type = "total_interactions",
                                            counterfactual_label = "treatment programme without monitoring")


levels(ly_gained_per_interaction_allages_long$scenario) <- monitoring_freq_allages_label

ggplot(data = subset(ly_gained_per_interaction_allages_long, by_year == 2100),
       aes(x=scenario, y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "#00BFC4") +
  ylab("Life-years saved\nper 10,000 incremental clinical interactions") +
  xlab("Monitoring frequency (years)") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Absolute deaths averted ----
cohort_deaths_averted_allages_long <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
                                 scenario_objects = list(out4,
                                                         out5,
                                                         out5a,
                                                         out5b,
                                                         out5c,
                                                         out5d,
                                                         out6,
                                                         out6a,
                                                         out6b,
                                                         out6c),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")

levels(cohort_deaths_averted_allages_long$scenario) <- monitoring_freq_allages_label

ggplot(subset(cohort_deaths_averted_allages_long, type == "proportion_averted"),
       aes(x =scenario, y= value*100)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "#00BFC4") +
  ylab("Percentage of HBV-related deaths averted") +
  labs(title = paste0("Effect of monitoring in the cohort compared to no monitoring")) +
  xlab("Monitoring frequency (years)") +
  ylim(0,100) +
  theme_bw() +
#  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Absolute LY saved ----
cohort_ly_gained_allages_long <-
  plot_ly_gained_cohort(counterfactual_object = out3,
                        scenario_objects = list(out4,
                                                out5,
                                                out5a,
                                                out5b,
                                                out5c,
                                                out5d,
                                                out6,
                                                out6a,
                                                out6b,
                                                out6c),
                        counterfactual_label = "treatment programme without monitoring")

levels(cohort_ly_gained_allages_long$counterfactual) <- monitoring_freq_allages_label

ggplot(subset(cohort_ly_gained_allages_long, type == "proportion_averted"),
       aes(x=counterfactual, y=value*100)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "#00BFC4") +
  ylab("Percentage of life-years saved") +
  labs(title = paste0("Effect of monitoring in the cohort compared to no monitoring")) +
  xlab("Monitoring frequency (years)") +
  ylim(0,5) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Maybe compare different frequencies to yearly monitoring here
#deaths_averted_per_assessment_long$frequency <- "Every 5 years"
#deaths_averted_per_assessment_long$frequency[deaths_averted_per_assessment_long$scenario %in% sub_yearly] <- "Yearly"
#deaths_averted_per_assessment_long$age_group <- gsub("Yearly ","", deaths_averted_per_assessment_long$scenario)
#deaths_averted_per_assessment_long$age_group <- gsub("5-yearly ","", deaths_averted_per_assessment_long$age_group)

# Incremental deaths averted vs incremental resource use ("cost-effectiveness" plane)----
# Including status quo
#out3 = 0, out4 = 10, out5d = 9, out5c = 8, out5b = 7, out5a = 6, out5 = 5,
#out6c = 4, out6b = 3, out6a = 2, out6 =1
freq_interactions <- rbind(
  cbind(scenario = "No monitoring",
        out3$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 10 years",
        out4$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 9 years",
        out5d$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 8 years",
        out5c$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 7 years",
        out5b$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 6 years",
        out5a$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 5 years",
        out5$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 4 years",
        out6c$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 3 years",
        out6b$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 2 years",
        out6a$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 1 year",
        out6$interactions[[16]]$total_interactions[-c(1:3)]))

freq_interactions_long <- gather(freq_interactions, key = "sim",
                                 value = "interactions", - scenario)

freq_cohort_deaths_averted <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out1,
                                 scenario_objects = list(out3,
                                                         out4,
                                                         out5,
                                                         out5a,
                                                         out5b,
                                                         out5c,
                                                         out5d,
                                                         out6,
                                                         out6a,
                                                         out6b,
                                                         out6c),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")
monitoring_freq_label <- list("No monitoring"="screen_2020_monit_0",
                                      "Every 10 years"="screen_2020_monit_10",
                                      "Every 9 years"="screen_2020_monit_9",
                                      "Every 8 years"="screen_2020_monit_8",
                                      "Every 7 years"="screen_2020_monit_7",
                                      "Every 6 years"="screen_2020_monit_6",
                                      "Every 5 years"="screen_2020_monit_5",
                                      "Every 4 years"="screen_2020_monit_4",
                                      "Every 3 years"="screen_2020_monit_3",
                                      "Every 2 years"="screen_2020_monit_2",
                                      "Every 1 year"="screen_2020_monit_1")
levels(freq_cohort_deaths_averted$scenario) <- monitoring_freq_label

freq_cohort_deaths_averted <- filter(freq_cohort_deaths_averted, type == "number_averted") %>%
  select(scenario, sim, value)
freq_cohort_deaths_averted$sim <- paste0("X", as.numeric(gsub("[^0-9]", "", freq_cohort_deaths_averted$sim)))
# Check simulations match:
unique(freq_cohort_deaths_averted$sim) == colnames(freq_interactions[-1])

incremental_df_freq <- left_join(freq_interactions_long,
                                 freq_cohort_deaths_averted, by = c("scenario", "sim"))
colnames(incremental_df_freq)[4] <- "deaths_averted"

# Add status quo point at 0 (reference)
incremental_df_freq <- rbind(data.frame(scenario = rep("Status quo", 183),
                                   sim = unique(incremental_df_freq$sim),
                                   interactions = rep(0,183),
                                   deaths_averted = rep(0,183)),
                             incremental_df_freq)

incremental_df_freq_summary <- group_by(incremental_df_freq, scenario) %>%
  summarise(median_interactions = median(interactions),
            median_deaths_averted = median(deaths_averted))

# CORRECT PLOT showing all cumulative deaths averted and incremental interactions
# compared to status quo of no treatment
# The x and y value of a point shows the cumulative number of deaths averted and interactions
# invested compared to no treatment. The slope of the lines show the increment in this
# between consecutive strategies.
# However in interpreting this need to keep in mind the type of interactions varies between
# different strategies
ggplot(data = incremental_df_freq) +
  geom_line(aes(y = deaths_averted, x = interactions, group = sim),
            colour = "grey", alpha = 0.3) +
  geom_point(aes(y = deaths_averted, x= interactions,
                 group = scenario, colour = scenario), alpha = 0.15) +
  stat_ellipse(data = incremental_df_freq[incremental_df_freq$scenario != "Status quo",],
               aes(y=deaths_averted,x=interactions,
                   group = scenario, fill= scenario),
               geom = "polygon",
               alpha = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_line(data = incremental_df_freq_summary, aes(x = median_interactions,
                                               y = median_deaths_averted), size = 1) +
  geom_point(data = incremental_df_freq_summary, aes(x = median_interactions,
                                                y = median_deaths_averted,
                                                group = scenario, colour = scenario),
             size = 5) +
  geom_point(data = incremental_df_freq_summary, aes(x = median_interactions,
                                                     y = median_deaths_averted,
                                                     group = scenario),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual(values = brewer.pal(11,"RdYlBu")) +
  scale_colour_manual("Monitoring frequency",
                      values = c("black", brewer.pal(11,"RdYlBu")),
                      labels = c("Status quo" = "No treatment")) +
  guides(fill=FALSE) +
  xlab("Incremental total clinical interactions") +
  ylab("Incremental HBV-related deaths averted") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  coord_flip()
# Need to use coord_flip to get scenario lines connected in the right order because of some negative values

# Same for life-years
freq_cohort_ly_saved <-
  plot_ly_gained_cohort(counterfactual_object = out1,
                                 scenario_objects = list(out3,
                                                         out4,
                                                         out5,
                                                         out5a,
                                                         out5b,
                                                         out5c,
                                                         out5d,
                                                         out6,
                                                         out6a,
                                                         out6b,
                                                         out6c),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")
colnames(freq_cohort_ly_saved)[1:2] <- c("scenario", "counterfactual")
levels(freq_cohort_ly_saved$scenario) <- monitoring_freq_label

freq_cohort_ly_saved <- filter(freq_cohort_ly_saved, type == "number_averted") %>%
  select(scenario, sim, value)
freq_cohort_ly_saved$sim <- paste0("X", as.numeric(gsub("[^0-9]", "", freq_cohort_ly_saved$sim)))
# Check simulations match:
unique(freq_cohort_ly_saved$sim) == colnames(freq_interactions[-1])

incremental_df_ly_freq <- left_join(freq_interactions_long,
                                    freq_cohort_ly_saved, by = c("scenario", "sim"))
colnames(incremental_df_ly_freq)[4] <- "ly_saved"

# Add status quo point at 0 (reference)
incremental_df_ly_freq <- rbind(data.frame(scenario = rep("Status quo", 183),
                                        sim = unique(incremental_df_ly_freq$sim),
                                        interactions = rep(0,183),
                                        ly_saved = rep(0,183)),
                                incremental_df_ly_freq)

incremental_df_ly_freq_summary <- group_by(incremental_df_ly_freq, scenario) %>%
  summarise(median_interactions = median(interactions),
            median_ly_saved = median(ly_saved))

# CORRECT PLOT showing all cumulative LY gained and incremental interactions
# compared to status quo of no treatment
# However in interpreting this need to keep in mind the type of interactions varies between
# different strategies
ggplot(data = incremental_df_ly_freq) +
  geom_line(aes(y = ly_saved, x = interactions, group = sim),
            colour = "grey", alpha = 0.3) +
  geom_point(aes(y = ly_saved, x= interactions,
                 group = scenario, colour = scenario), alpha = 0.15) +
  stat_ellipse(data = incremental_df_ly_freq[incremental_df_ly_freq$scenario != "Status quo",],
               aes(y=ly_saved,x=interactions,
                   group = scenario, fill= scenario),
               geom = "polygon",
               alpha = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_line(data = incremental_df_ly_freq_summary, aes(x = median_interactions,
                                                    y = median_ly_saved), size = 1) +
  geom_point(data = incremental_df_ly_freq_summary, aes(x = median_interactions,
                                                     y = median_ly_saved,
                                                     group = scenario, colour = scenario),
             size = 5) +
  geom_point(data = incremental_df_ly_freq_summary, aes(x = median_interactions,
                                                     y = median_ly_saved,
                                                     group = scenario),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual(values = brewer.pal(11,"RdYlBu")) +
  scale_colour_manual("Monitoring frequency",
                      values = c("black", brewer.pal(11,"RdYlBu")),
                      labels = c("Status quo" = "No treatment")) +
  guides(fill=FALSE) +
  xlab("Incremental total clinical interactions") +
  ylab("Incremental life-years saved") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  coord_flip()

# 4) Incremental benefits and interactions between a set of strategies (monitoring by age) ----

## Approach 1: monitoring given age group from entry until aging out ----
# Incrementally intensive strategies are: 15-30 years, 15-45 years, 15+ (all ages)

# INCREMENTAL INTERACTIONS IN EACH SCENARIO BY 2100 COMPARED TO NO MONITORING
# out3, monit_out6, monit_out7, out5 (5 yearly)

# Order could be:
# 5-yearly in 15-30 => yearly in 15-30 (better or worse?) => if better:
# yearly in 15-30 + 5-yearly in 30-45 => yearly in 15-30 + yearly in 30-45 etc
# Right now, going from Yearly 15-30 to 5-yearly 15-45, which makes no sense
# Can I derive the increment for this from my current simulations?
# However, interesting observation that going from yearly monitoring in 15-45
# to monitoring 5-yearly across all ages leads to more deaths averted at a lower cost!

# Create data frame with all interactions and outcomes of interest
interactions1 <- rbind(
  cbind(scenario = "screen_2020_monit_5",
        left_join(gather(out5$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(out5$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_1",
        left_join(gather(out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(out6$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim1",
        left_join(gather(monit_out1$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out1$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim2",
        left_join(gather(monit_out2$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out2$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim6",
        left_join(gather(monit_out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out6$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim7",
        left_join(gather(monit_out7$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out7$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim"))
)
interactions1$total_interactions <- interactions1$monitoring_assessments+
  interactions1$treatment_initiations
levels(interactions1$scenario) <- scenario_labels
interactions1$sim <- gsub("[^0-9]", "", interactions1$sim)

# Extract person-years on treatment
interactions1_py_on_treatment <- rbind(
  data.frame(scenario = "screen_2020_monit_5",
             sim = names(out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_1",
             sim = names(out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment= out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim1",
             sim = names(monit_out1$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out1$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim2",
             sim = names(monit_out2$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out2$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim6",
             sim = names(monit_out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim7",
             sim = names(monit_out7$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out7$py_on_treatment[[16]]-out3$py_on_treatment[[16]]))
levels(interactions1_py_on_treatment$scenario) <- scenario_labels

# Combine interactions and assign costs based on Shevanthi's Gambia paper:
# $15.77 per monitoring assessment and $84.88 per person-year of treatment
interactions1 <- left_join(interactions1, interactions1_py_on_treatment,
                           by=c("scenario", "sim")) %>%
  mutate(monitoring_cost = monitoring_assessments*15.77,
         treatment_cost = py_on_treatment*84.88,
         total_cost = monitoring_cost+treatment_cost)
# Check person-years on treatment decreases if older age groups are included:
group_by(interactions1, scenario) %>%
  summarise(median(py_on_treatment/treatment_initiations))

# Outcome 1: HBV related deaths averted
cohort_deaths_averted_long <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
                                 scenario_objects = list(out5,
                                                         out6,
                                                         monit_out1,
                                                         monit_out2,
                                                         monit_out6,
                                                         monit_out7),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")
levels(cohort_deaths_averted_long$scenario) <- scenario_labels

df1 <- subset(cohort_deaths_averted_long, type == "number_averted") %>%
  select(scenario, sim, value)
df1$sim <- gsub("[^0-9]", "", df1$sim)
df1 <- df1 %>%
  left_join(interactions1, by = c("scenario", "sim")) %>%
  rename(deaths_averted = value) %>%
  drop_na()

# Outcome 2: Life-years saved
cohort_ly_gained_long <-
  plot_ly_gained_cohort(counterfactual_object = out3,
                        scenario_objects = list(out5,
                                                out6,
                                                monit_out1,
                                                monit_out2,
                                                monit_out6,
                                                monit_out7),
                        outcome_to_plot = "number_averted",
                        counterfactual_label = "treatment programme without monitoring")
colnames(cohort_ly_gained_long)[1:2] <- c("scenario", "counterfactual")
levels(cohort_ly_gained_long$scenario) <- scenario_labels

df1_ly <- subset(cohort_ly_gained_long, type == "number_averted") %>%
  select(scenario, sim, value)
df1_ly$sim <- gsub("[^0-9]", "", df1_ly$sim)
df1 <- df1_ly %>%
  left_join(df1, by = c("scenario", "sim")) %>%
  rename(ly_saved = value) %>%
  drop_na()

# Split by monitoring frequency
df1$frequency <- "Every 5 years"
df1$frequency[df1$scenario %in% c("Yearly all ages", "Yearly 15-45", "Yearly 15-30")] <- "Every 1 year"
# Split by age group
df1$age_group <- "15+ (all ages)"
df1$age_group[df1$scenario %in% c("Yearly 15-45", "5-yearly 15-45")] <- "15-45"
df1$age_group[df1$scenario %in% c("Yearly 15-30", "5-yearly 15-30")] <- "15-30"

# Add reference no monitoring point at 0
df1$scenario <- as.character(df1$scenario)
df1 <- rbind(c("No monitoring", "X", rep(0,(ncol(df1)-4)), "Every 5 years", "15-30"),
             c("No monitoring", "X", rep(0,(ncol(df1)-4)), "Every 5 years", "15-45"),
             c("No monitoring", "X", rep(0,(ncol(df1)-4)), "Every 5 years", "15+ (all ages)"),
             c("No monitoring", "X", rep(0,(ncol(df1)-4)), "Every 1 year", "15-30"),
             c("No monitoring", "X", rep(0,(ncol(df1)-4)), "Every 1 year", "15-45"),
             c("No monitoring", "X", rep(0,(ncol(df1)-4)), "Every 1 year", "15+ (all ages)"),
             df1)
df1[,-c(1:2,ncol(df1)-1,ncol(df1))] <- apply(df1[,-c(1:2,ncol(df1)-1,ncol(df1))], 2, as.numeric)
df1$scenario <- factor(df1$scenario)

# Use BCEA package to find dominated and extended dominated strategies
library(BCEA)
find_dominated_strategies <- function(df, exposure, outcome) {
plot <- ceef.plot(bcea(e=cbind(matrix(rep(0, 183)),
                       as.matrix(select(df, scenario, sim, outcome) %>%
                                   filter(sim != "X") %>%
                                   spread(key="scenario", value = outcome)%>%
                                   select(-sim))),
               c=cbind(matrix(rep(0, 183)),
                       as.matrix(select(df, scenario, sim, exposure) %>%
                                   filter(sim != "X") %>%
                                   spread(key="scenario", value = exposure)%>%
                                   select(-sim))),
               ref=1,
               interventions=c("No monitoring", colnames(select(df, scenario, sim, outcome) %>%
                 filter(sim != "X") %>%
                 spread(key="scenario", value = outcome)%>%
                 select(-sim))),
               Kmax=50000000,
               plot=FALSE),
          graph="base", relative = FALSE)
return(plot)
}
find_dominated_strategies(df1, "total_interactions", "deaths_averted")
find_dominated_strategies(df1,"total_cost", "deaths_averted")
find_dominated_strategies(df1,"total_interactions", "ly_saved")
find_dominated_strategies(df1, "total_cost", "ly_saved")
# Deaths averted and interactions: Yearly 15-45 (A), Yearly 15-30 (E)
# Deaths averted and cost: Yearly 15-45 (A), Yearly 15-30 (E), 5-yearly 15-45 (E), 5-yearly 15-30 (E),
# LY saved and interactions: Yearly 15-45 (A), Yearly 15-30 (E)
# LY saved and cost: Yearly 15-45 (A), Yearly 15-30 (E), 5-yearly 15-30 (E)

# Add labels for this to dataframe
df1$frontier_interactions <- "Include"
df1$frontier_interactions[df1$scenario %in% c("Yearly 15-45", "Yearly 15-30")] <- "Dominated"
df1$frontier_deaths_averted_cost <- "Include"
df1$frontier_deaths_averted_cost[df1$scenario %in% c("Yearly 15-45", "Yearly 15-30",
                                                     "5-yearly 15-45", "5-yearly 15-30")] <- "Dominated"
df1$frontier_ly_saved_cost <- "Include"
df1$frontier_ly_saved_cost[df1$scenario %in% c("Yearly 15-45", "Yearly 15-30", "5-yearly 15-30")] <-
  "Dominated"

df1_summary <- df1 %>%
  group_by(scenario, frequency, age_group, frontier_interactions, frontier_deaths_averted_cost, frontier_ly_saved_cost) %>%
  summarise(median_deaths_averted = median(deaths_averted),
            median_ly_saved = median(ly_saved),
            median_interactions = median(total_interactions),
            median_cost = median(total_cost))

# Plots ----
# Plot: interactions vs deaths averted
ggplot(df1) +
  geom_line(data= subset(df1, frontier_interactions== "Include"),
            aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  stat_ellipse(data=subset(df1, scenario != "No monitoring"),
                aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted),
                   fill= reorder(scenario, deaths_averted)),
               geom = "polygon", na.rm = FALSE, alpha = 0.3) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  # Overlay median
  geom_line(data = subset(df1_summary, frontier_interactions == "Include"),
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual("Monitoring strategies\nstarting at entry\ninto cohort",
                    values=brewer.pal(6,"RdYlBu")) +
  scale_colour_manual("Monitoring strategies\nstarting at entry\ninto cohort",
                      values=c("black", brewer.pal(6,"RdYlBu"))) +
  guides(fill=FALSE) +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs deaths averted (separated by frequency)
ggplot(df1) +
  geom_line(aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = df1_summary,
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~frequency) +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs deaths averted
ggplot(df1) +
  geom_line(data= subset(df1, frontier_deaths_averted_cost == "Include"),
            aes(x = deaths_averted, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_cost,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_cost,
                   group = reorder(scenario, deaths_averted),
                   fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df1_summary, frontier_deaths_averted_cost == "Include"),
            aes(y = median_cost,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_cost,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_cost,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental cost (USD)") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs LY saved
ggplot(df1) +
  geom_line(data= subset(df1, frontier_interactions == "Include"),
            aes(x = ly_saved, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_interactions,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_interactions,
                   group = reorder(scenario, ly_saved), fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df1_summary, frontier_interactions == "Include"),
            aes(y = median_interactions,
                x = median_ly_saved), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental number of clinical interactions") +
#  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs LY saved
ggplot(df1) +
  geom_line(data= subset(df1, frontier_ly_saved_cost == "Include"),
            aes(x = ly_saved, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_cost,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_cost,
                   group = reorder(scenario, ly_saved),
                   fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df1_summary, frontier_ly_saved_cost == "Include"),
            aes(y = median_cost,
                x = median_ly_saved), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
#  geom_abline(intercept =0, slope = 240, linetype = "dashed") +   # WTP: need to compare to status quo?
#  geom_abline(intercept =0, slope = 1460, linetype = "dashed") +  # WTP: need to compare to status quo?
#  geom_abline(intercept =0, slope = 487, linetype = "dashed") +   # WTP: need to compare to status quo?
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental cost (USD)") +
 # xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

## Approach 2: monitoring given age cohort from entry until death ----
# Incrementally intensive strategies are: 45+ years, 30+ years, 15+ (all ages)
# Start with 5-yearly frequency
# Increment between no monitoring and status quo: out3-out1
# Increment between monitoring 45+ year olds and no monitoring in any age group: a2_out5-a2_out3
# Increment between adding 30-45 year old to monitoring vs only 45+ year olds:
# (a3_out5-a3_out3)-(a2_out5-a2_out3)
# Increment between adding 30-45 year old to monitoring vs only 30+ year olds:
#(a1_out5-a1_out3)-(a3_out5-a3_out3)

# Create data frame with all interactions and outcomes of interest
interactions2 <- rbind(
  cbind(scenario = "5-yearly 15+ (all ages)",
        left_join(gather(out5$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(out5$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Yearly 15+ (all ages)",
        left_join(gather(out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(out6$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "5-yearly 45+",
        left_join(gather(a2_out5$interactions[[16]]$total_assessed[-c(1:3)]-a2_out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(a2_out5$interactions[[16]]$total_treated[-c(1:3)]-a2_out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "5-yearly 30+",
        left_join(gather(a3_out5$interactions[[16]]$total_assessed[-c(1:3)]-a3_out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(a3_out5$interactions[[16]]$total_treated[-c(1:3)]-a3_out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Yearly 45+",
        left_join(gather(a2_out6$interactions[[16]]$total_assessed[-c(1:3)]-a2_out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(a2_out6$interactions[[16]]$total_treated[-c(1:3)]-a2_out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Yearly 30+",
        left_join(gather(a3_out6$interactions[[16]]$total_assessed[-c(1:3)]-a3_out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(a3_out6$interactions[[16]]$total_treated[-c(1:3)]-a3_out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim"))
)
interactions2$total_interactions <- interactions2$monitoring_assessments+
  interactions2$treatment_initiations
#levels(interactions2$scenario) <- scenario_labels
interactions2$sim <- gsub("[^0-9]", "", interactions2$sim)

# Extract person-years on treatment
interactions2_py_on_treatment <- rbind(
  data.frame(scenario = "5-yearly 15+ (all ages)",
             sim = names(out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "Yearly 15+ (all ages)",
             sim = names(out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment= out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "5-yearly 45+",
             sim = names(a2_out5$py_on_treatment[[16]]-a2_out3$py_on_treatment[[16]]),
             py_on_treatment = a2_out5$py_on_treatment[[16]]-a2_out3$py_on_treatment[[16]]),
  data.frame(scenario = "5-yearly 30+",
             sim = names(a3_out5$py_on_treatment[[16]]-a3_out3$py_on_treatment[[16]]),
             py_on_treatment = a3_out5$py_on_treatment[[16]]-a3_out3$py_on_treatment[[16]]),
  data.frame(scenario = "Yearly 45+",
             sim = names(a2_out6$py_on_treatment[[16]]-a2_out3$py_on_treatment[[16]]),
             py_on_treatment = a2_out6$py_on_treatment[[16]]-a2_out3$py_on_treatment[[16]]),
  data.frame(scenario = "Yearly 30+",
             sim = names(a3_out6$py_on_treatment[[16]]-a3_out3$py_on_treatment[[16]]),
             py_on_treatment = a3_out6$py_on_treatment[[16]]-a3_out3$py_on_treatment[[16]]))

# Combine interactions and assign costs based on Shevanthi's Gambia paper:
# $15.77 per monitoring assessment and $84.88 per person-year of treatment
interactions2 <- left_join(interactions2, interactions2_py_on_treatment,
                           by=c("scenario", "sim")) %>%
  mutate(monitoring_cost = monitoring_assessments*15.77,
         treatment_cost = py_on_treatment*84.88,
         total_cost = monitoring_cost+treatment_cost)
# Check person-years on treatment increases if younger age groups are included:
group_by(interactions2, scenario) %>%
  summarise(median(py_on_treatment/treatment_initiations))

# Outcome 1: HBV related deaths averted
cohort_deaths_averted2 <- rbind(
  cbind(scenario = "5-yearly 45+",
        a2_out3$cum_hbv_deaths[[16]][-c(1:3)]-a2_out5$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "5-yearly 30+",
        a3_out3$cum_hbv_deaths[[16]][-c(1:3)]-a3_out5$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "5-yearly 15+ (all ages)",
        out3$cum_hbv_deaths[[16]][-c(1:3)]-out5$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 45+",
        a2_out3$cum_hbv_deaths[[16]][-c(1:3)]-a2_out6$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 30+",
        a3_out3$cum_hbv_deaths[[16]][-c(1:3)]-a3_out6$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 15+ (all ages)",
        out3$cum_hbv_deaths[[16]][-c(1:3)]-out6$cum_hbv_deaths[[16]][-c(1:3)]))
colnames(cohort_deaths_averted2)[-1] <- as.numeric(gsub("[^0-9]", "",
                                                           colnames(cohort_deaths_averted2)[-1]))
cohort_deaths_averted2_long <- gather(cohort_deaths_averted2, key = "sim",
                                         value = "deaths_averted", - scenario)

df2 <- cohort_deaths_averted2_long %>%
  left_join(interactions2, by = c("scenario", "sim")) %>%
  drop_na()

# Outcome 2: Life-years saved

# Life-years saved
cohort_ly_gained2 <- rbind(
  cbind(scenario = "5-yearly 45+",
        a2_out5$ly[[16]][-c(1:3)]-a2_out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "5-yearly 30+",
        a3_out5$ly[[16]][-c(1:3)]-a3_out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "5-yearly 15+ (all ages)",
        out5$ly[[16]][-c(1:3)]-out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 45+",
        a2_out6$ly[[16]][-c(1:3)]-a2_out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 30+",
        a3_out6$ly[[16]][-c(1:3)]-a3_out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 15+ (all ages)",
        out6$ly[[16]][-c(1:3)]-out3$ly[[16]][-c(1:3)]))
colnames(cohort_ly_gained2)[-1] <- as.numeric(gsub("[^0-9]", "",
                                                      colnames(cohort_ly_gained2)[-1]))
cohort_ly_gained2_long <- gather(cohort_ly_gained2, key = "sim",
                                    value = "ly_saved", - scenario)

df2 <- cohort_ly_gained2_long %>%
  left_join(df2, by = c("scenario", "sim")) %>%
  drop_na()

# Split by monitoring frequency
df2$frequency <- "Every 5 years"
df2$frequency[df2$scenario %in% c("Yearly 15+ (all ages)", "Yearly 45+", "Yearly 30+")] <- "Every 1 year"
# Split by birth cohort age group
df2$age_group <- "15+ (all ages)"
df2$age_group[df2$scenario %in% c("Yearly 45+", "5-yearly 45+")] <- "45+"
df2$age_group[df2$scenario %in% c("Yearly 30+", "5-yearly 30+")] <- "30+"

# Add reference no monitoring point at 0
df2$scenario <- as.character(df2$scenario)
df2 <- rbind(c("No monitoring", "X", rep(0,(ncol(df2)-4)), "Every 5 years", "45+"),
             c("No monitoring", "X", rep(0,(ncol(df2)-4)), "Every 5 years", "30+"),
             c("No monitoring", "X", rep(0,(ncol(df2)-4)), "Every 5 years", "15+ (all ages)"),
             c("No monitoring", "X", rep(0,(ncol(df2)-4)), "Every 1 year", "45+"),
             c("No monitoring", "X", rep(0,(ncol(df2)-4)), "Every 1 year", "30+"),
             c("No monitoring", "X", rep(0,(ncol(df2)-4)), "Every 1 year", "15+ (all ages)"),
             df2)
df2[,-c(1:2,ncol(df2)-1,ncol(df2))] <- apply(df2[,-c(1:2,ncol(df2)-1,ncol(df2))], 2, as.numeric)
df2$scenario <- factor(df2$scenario)

# Find dominated strategies for each combination of exposure and outcome using BCEA package
find_dominated_strategies(df2, "total_interactions", "deaths_averted")
find_dominated_strategies(df2,"total_cost", "deaths_averted")
find_dominated_strategies(df2,"total_interactions", "ly_saved")
find_dominated_strategies(df2, "total_cost", "ly_saved")
# Deaths averted & interactions, LY saved & any exp.:
# Yearly 30+, Yearly 45+, 5-yearly 30+, 5-yearly 45+
# Deaths averted & cost: Yearly 30+, 5-yearly 45+, Yearly 45+

# Add labels to this dataframe
df2$frontier_deaths_averted_interactions <- "Include"
df2$frontier_deaths_averted_interactions[df2$scenario %in% c("Yearly 30+", "5-yearly 30+",
                                                "5-yearly 45+", "Yearly 45+")] <- "Dominated"
df2$frontier_deaths_averted_cost <- "Include"
df2$frontier_deaths_averted_cost[df2$scenario %in% c("Yearly 30+", "5-yearly 45+", "Yearly 45+")] <-
  "Dominated"
df2$frontier_ly_saved <- "Include"
df2$frontier_ly_saved[df2$scenario %in% c("Yearly 30+", "5-yearly 30+",
                                          "5-yearly 45+", "Yearly 45+")] <- "Dominated"

df2_summary <- df2 %>%
  group_by(scenario, frequency, age_group, frontier_deaths_averted_interactions,
           frontier_deaths_averted_cost, frontier_ly_saved) %>%
  summarise(median_deaths_averted = median(deaths_averted),
            median_ly_saved = median(ly_saved),
            median_interactions = median(total_interactions),
            median_cost = median(total_cost))

# Plots ----
# Plot: interactions vs deaths averted (combined)
ggplot(df2) +
  geom_line(data= subset(df2, frontier_deaths_averted_interactions == "Include"),
            aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)),
               alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df2_summary, frontier_deaths_averted_interactions == "Include"),
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  #xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs deaths averted (separated by birth cohort)
ggplot(df2) +
  geom_line(data= df2,
            aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = df2_summary,
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~age_group) +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  #xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs deaths averted (separate by frequency)
ggplot(df2) +
  geom_line(data= df2,
            aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = df2_summary,
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~frequency) +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs deaths averted
ggplot(df2) +
  geom_line(data= subset(df2, frontier_cost == "Include"),
            aes(x = deaths_averted, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_cost,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_cost,
                   group = reorder(scenario, deaths_averted),
                   fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df2_summary, frontier_cost == "Include"),
            aes(y = median_cost,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental cost (USD)") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs LY saved
ggplot(df2) +
  geom_line(data= subset(df2, frontier_interactions == "Include"),
            aes(x = ly_saved, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_interactions,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_interactions,
                   group = reorder(scenario, ly_saved), fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df2_summary, frontier_interactions == "Include"),
            aes(y = median_interactions,
                x = median_ly_saved), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental number of clinical interactions") +
  #  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs LY saved
ggplot(df2) +
  geom_line(data= subset(df2, frontier_cost == "Include"),
            aes(x = ly_saved, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_cost,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_cost,
                   group = reorder(scenario, ly_saved),
                   fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df2_summary, frontier_cost == "Include"),
            aes(y = median_cost,
                x = median_ly_saved), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  geom_abline(intercept =0, slope = 240, linetype = "dashed") +
  geom_abline(intercept =0, slope = 1460, linetype = "dashed") +
  geom_abline(intercept =0, slope = 487, linetype = "dashed") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental cost (USD)") +
  # xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs LY saved (separated by frequency - no dominated strategies)
ggplot(df2) +
  geom_line(data= df2,
            aes(x = ly_saved, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_cost,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_cost,
                   group = reorder(scenario, ly_saved),
                   fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = df2_summary,
            aes(y = median_cost,
                x = median_ly_saved), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~frequency) +
  geom_abline(intercept =0, slope = 240, linetype = "dashed") +
  geom_abline(intercept =0, slope = 1460, linetype = "dashed") +
  geom_abline(intercept =0, slope = 487, linetype = "dashed") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental cost (USD)") +
  # xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# This plot shows the increase in deaths averted and # total interactions between each 2 strategies,
# with "treatment without monitoring" being at 0 (status quo of no treatment is shown as negative to that).
# The assumption here is that everyone is monitored from the start, but it stops at some age (e.g. 30, 45).
# Would ideally exclude the lines/points outside the 95% interval

## Combine approaches 1 and 2 ----
df1 <- select(df, scenario, sim, interactions, value) %>%
  filter(scenario %in% c("No treatment", "No monitoring",
                         "5-yearly 15-30", "Yearly 15-30",
                         "5-yearly 15-45", "5-yearly all ages",
                         "Yearly all ages",  "Yearly 15-45"))
colnames(df1)[4] <- "deaths_averted"
df1$scenario <- factor(paste("By age:", df1$scenario))
# Remove duplicate scenarios (all ages) from approach 2:
df2 <- filter(age_cohort_df, scenario != "5-yearly 15+ (all ages)" & scenario != "Yearly 15+ (all ages)")
df2$scenario <- factor(paste("By birth cohort:", df2$scenario))
combined_df <- rbind(df1, df2)
# Label dominated strategies (based on median):
combined_df$scenario <- as.character(combined_df$scenario)
combined_df$scenario <- factor(combined_df$scenario,
                               levels = c("By birth cohort: No monitoring",
                                          "By age: 5-yearly 15-30",
                                          "By birth cohort: 5-yearly 45+",
                                          "By age: Yearly 15-30",
                                          "By age: 5-yearly 15-45",
                                          "By birth cohort: 5-yearly 30+",
                                          "By age: 5-yearly all ages",
                                          "By age: Yearly all ages",
                                          "By birth cohort: Yearly 45+",
                                          "By age: Yearly 15-45",
                                          "By birth cohort: Yearly 30+"))
levels(combined_df$scenario) <- list("No monitoring" = "By birth cohort: No monitoring",
                                     "By age: 15-30,\nevery 5 years" = "By age: 5-yearly 15-30",
                                     "By birth cohort: 45+,\nevery 5 years" = "By birth cohort: 5-yearly 45+",
                                     "By age: 15-30,\nevery 1 year" = "By age: Yearly 15-30",
                                     "By age: 15-45,\nevery 5 years" = "By age: 5-yearly 15-45",
                                     "By birth cohort:\n30+, every 5 years" = "By birth cohort: 5-yearly 30+",
                                     "All ages,\nevery 5 years" = "By age: 5-yearly all ages",
                                     "All ages,\nevery 1 year" = "By age: Yearly all ages",
                                     "By birth cohort: 45+,\nevery 1 year" = "By birth cohort: Yearly 45+",
                                     "By age: 15-45,\nevery 1 year" = "By age: Yearly 15-45",
                                     "By birth cohort: 30+,\nevery 1 year" = "By birth cohort: Yearly 30+")

combined_df$frontier <- "Include"
combined_df$frontier[combined_df$scenario %in% c("By birth cohort: 45+,\nevery 1 year",
                                                 "By age: 15-45,\nevery 1 year",
                                                 "By birth cohort: 30+,\nevery 1 year")] <- "Dominated"


# Summary:
combined_df_summary <- group_by(combined_df, scenario, frontier) %>%
  summarise(median_interactions = median(interactions),
            median_deaths_averted = median(deaths_averted))

# Plot:
ggplot(combined_df) +
  geom_point(aes(x = deaths_averted,interactions,
                 group = scenario, colour = scenario), alpha = 0.4) +
  geom_line(data=subset(combined_df, frontier == "Include"),
            aes(x = deaths_averted,interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  stat_ellipse(geom = "polygon",
               aes(deaths_averted,interactions,
                   group = scenario, fill= scenario), alpha = 0.3) +
  geom_line(data = subset(combined_df_summary, frontier == "Include"),
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = combined_df_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = scenario, colour = scenario), size = 5) +
  geom_point(data = combined_df_summary,
             aes(y = median_interactions,
             x = median_deaths_averted,
              group = scenario),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
#  xlim(-100,7100) +
#  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# This one might look quite different in life-years or if using costs (treatment time)

# Life-years saved
df_ly1 <- select(df_ly, scenario, sim, interactions, value) %>%
  filter(scenario %in% c("No treatment", "No monitoring",
                         "5-yearly 15-30", "Yearly 15-30",
                         "5-yearly 15-45", "5-yearly all ages",
                         "Yearly all ages",  "Yearly 15-45"))
colnames(df_ly1)[4] <- "ly_saved"
df_ly1$scenario <- factor(paste("By age:", df_ly1$scenario))
# Remove duplicate scenarios (all ages) from approach 2:
df_ly2 <- filter(age_cohort_ly_df, scenario != "5-yearly 15+ (all ages)" &
                   scenario != "Yearly 15+ (all ages)")
df_ly2$scenario <- factor(paste("By birth cohort:", df_ly2$scenario))
combined_df_ly <- rbind(df_ly1, df_ly2)
# Label dominated strategies (based on median):
combined_df_ly$scenario <- as.character(combined_df_ly$scenario)
combined_df_ly$scenario <- factor(combined_df_ly$scenario,
                               levels = c("By birth cohort: No monitoring",
                                          "By age: 5-yearly 15-30",
                                          "By birth cohort: 5-yearly 45+",
                                          "By age: Yearly 15-30",
                                          "By age: 5-yearly 15-45",
                                          "By birth cohort: 5-yearly 30+",
                                          "By age: 5-yearly all ages",
                                          "By age: Yearly all ages",
                                          "By birth cohort: Yearly 45+",
                                          "By age: Yearly 15-45",
                                          "By birth cohort: Yearly 30+"))
levels(combined_df_ly$scenario) <- list("No monitoring" = "By birth cohort: No monitoring",
                                     "By age: 15-30,\nevery 5 years" = "By age: 5-yearly 15-30",
                                     "By birth cohort: 45+,\nevery 5 years" = "By birth cohort: 5-yearly 45+",
                                     "By age: 15-30,\nevery 1 year" = "By age: Yearly 15-30",
                                     "By age: 15-45,\nevery 5 years" = "By age: 5-yearly 15-45",
                                     "By birth cohort:\n30+, every 5 years" = "By birth cohort: 5-yearly 30+",
                                     "All ages,\nevery 5 years" = "By age: 5-yearly all ages",
                                     "All ages,\nevery 1 year" = "By age: Yearly all ages",
                                     "By birth cohort: 45+,\nevery 1 year" = "By birth cohort: Yearly 45+",
                                     "By age: 15-45,\nevery 1 year" = "By age: Yearly 15-45",
                                     "By birth cohort: 30+,\nevery 1 year" = "By birth cohort: Yearly 30+")

combined_df_ly$frontier <- "Include"
combined_df_ly$frontier[combined_df_ly$scenario %in% c("By birth cohort: 45+,\nevery 1 year",
                                                 "By age: 15-45,\nevery 1 year",
                                                 "By birth cohort: 30+,\nevery 1 year",
                                                 "By birth cohort: 45+,\nevery 5 years",
                                                 "By birth cohort:\n30+, every 5 years")] <- "Dominated"


# Summary:
combined_df_ly_summary <- group_by(combined_df_ly, scenario, frontier) %>%
  summarise(median_interactions = median(interactions),
            median_ly_saved = median(ly_saved))

# Plot:
ggplot(combined_df_ly) +
  geom_point(aes(x = ly_saved,interactions,
                 group = scenario, colour = scenario), alpha = 0.4) +
  geom_line(data=subset(combined_df_ly, frontier == "Include"),
            aes(x = ly_saved,interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  stat_ellipse(geom = "polygon",
               aes(ly_saved,interactions,
                   group = scenario, fill= scenario), alpha = 0.3) +
  geom_line(data = subset(combined_df_ly_summary, frontier == "Include"),
            aes(y = median_interactions,
                x = median_ly_saved), size = 1) +
  geom_point(data = combined_df_ly_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = scenario, colour = scenario), size = 5) +
  geom_point(data = combined_df_ly_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = scenario),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental number of clinical interactions") +
  #  xlim(-100,7100) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))


## Add cost to approach 1 ----
# Compared to no monitoring, the cost comes from monitoring assessments (total_assessed)
# and treatment (person-years on treatment due to initiations after the initial screen)
interactions1_monitoring <- rbind(
  cbind(scenario = "screen_2020_monit_5",
        gather(out5$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "screen_2020_monit_1",
        gather(out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "screen_2020_monit_sim1",
        gather(monit_out1$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "screen_2020_monit_sim2",
        gather(monit_out2$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "screen_2020_monit_sim6",
        gather(monit_out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "screen_2020_monit_sim7",
        gather(monit_out7$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
               key = "sim", value = "value"))
)
levels(interactions1_monitoring$scenario) <- scenario_labels
colnames(interactions1_monitoring)[3] <- "monitoring_assessments"
interactions1_monitoring$sim <- gsub("[^0-9]", "", interactions1_monitoring$sim)

interactions1_py_on_treatment <- rbind(
  data.frame(scenario = "screen_2020_monit_5",
             sim = names(out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_1",
             sim = names(out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment= out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim1",
             sim = names(monit_out1$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out1$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim2",
             sim = names(monit_out2$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out2$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim6",
             sim = names(monit_out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim7",
             sim = names(monit_out7$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out7$py_on_treatment[[16]]-out3$py_on_treatment[[16]]))
levels(interactions1_py_on_treatment$scenario) <- scenario_labels

# Assign costs based on Shevanthi's Gambia paper:
# $15.77 per monitoring assessment and $84.88 per person-year of treatment
interactions1_cost <- left_join(interactions1_monitoring, interactions1_py_on_treatment,
                                by=c("scenario", "sim")) %>%
  mutate(monitoring_cost = monitoring_assessments*15.77,
         treatment_cost = py_on_treatment*84.88,
         total_cost = monitoring_cost+treatment_cost)

df1_cost <- subset(cohort_deaths_averted_long, type == "number_averted")
df1_cost$sim <- gsub("[^0-9]", "", df1_cost$sim)
df1_cost <- df1_cost %>%
  left_join(interactions1_cost, by = c("scenario", "sim")) %>%
  drop_na()

df1_cost_summary <- df1_cost %>%
  group_by(scenario) %>%
  summarise(median_deaths = median(value),
            median_cost = median(total_cost))

ggplot(df1_cost) +
  geom_line(aes(y = total_cost, x= value,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(y=total_cost, x= value,
                 group = scenario, colour = scenario), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(y=total_cost, x=value, group = scenario, fill= scenario), alpha = 0.3) +
  geom_point(data = df1_cost_summary,
             aes(y = median_cost, x = median_deaths, colour = scenario), size = 5) +
  geom_line(data = df1_cost_summary,
            aes(y = median_cost, x = median_deaths), colour = "black", size = 1) +
  scale_colour_brewer(palette = "Set1", direction = -1) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(colour = "Monitoring strategy", fill = "Monitoring strategy") +
  #  scale_fill_viridis_d() +
  #  scale_colour_viridis_d() +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental cost") +
  xlim(-100,7100) +
#  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))

# 5) Test: monitoring frequencies by separate age groups

# Yearly: sim1,sim3, sim5, 5-yearly: sim6,sim8,sim10

# Create data frame with all interactions and outcomes of interest
freq_by_age_interactions <- rbind(
  cbind(scenario = "screen_2020_monit_sim1",
        left_join(gather(monit_out1$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out1$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim3",
        left_join(gather(monit_out3$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out3$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim5",
        left_join(gather(monit_out5$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out5$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim6",
        left_join(gather(monit_out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out6$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim8",
        left_join(gather(monit_out8$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out8$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim10",
        left_join(gather(monit_out10$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out10$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim"))
)
freq_by_age_interactions$total_interactions <-freq_by_age_interactions$monitoring_assessments+
  freq_by_age_interactions$treatment_initiations
levels(freq_by_age_interactions$scenario) <- scenario_labels
freq_by_age_interactions$sim <- gsub("[^0-9]", "", freq_by_age_interactions$sim)

# Outcome 1: HBV related deaths averted
freq_by_age_deaths_averted_long <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
                                 scenario_objects = list(monit_out1,
                                                         monit_out3,
                                                         monit_out5,
                                                         monit_out6,
                                                         monit_out8,
                                                         monit_out10),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")
levels(freq_by_age_deaths_averted_long$scenario) <- scenario_labels

freq_by_age_df <- subset(freq_by_age_deaths_averted_long, type == "number_averted") %>%
  select(scenario, sim, value)
freq_by_age_df$sim <- gsub("[^0-9]", "", freq_by_age_df$sim)
freq_by_age_df <- freq_by_age_df %>%
  left_join(freq_by_age_interactions, by = c("scenario", "sim")) %>%
  rename(deaths_averted = value) %>%
  drop_na()

freq_by_age_df$age_group <- "15-30"
freq_by_age_df$age_group[freq_by_age_df$scenario %in% c("Yearly 30-45", "5-yearly 30-45")] <- "30-45"
freq_by_age_df$age_group[freq_by_age_df$scenario %in% c("Yearly 45+", "5-yearly 45+")] <- "45+"

# Add reference no monitoring point at 0
freq_by_age_df <- rbind(c("No monitoring", "X", rep(0,(ncol(freq_by_age_df)-3)), "15-30"),
                        c("No monitoring", "X", rep(0,(ncol(freq_by_age_df)-3)), "30-45"),
                        c("No monitoring", "X", rep(0,(ncol(freq_by_age_df)-3)), "45+"),
                        freq_by_age_df)
freq_by_age_df[,-c(1:2,ncol(freq_by_age_df))] <-
  apply(freq_by_age_df[,-c(1:2, ncol(freq_by_age_df))], 2, as.numeric)

# Manually find dominated strategies for each combination of exposure and outcome,
# based on median
# For deaths averted outcomes, it is Yearly 15-45
# For LY saved outcomes, it is Yearly 15-45 and Yearly 15-30
#df1$frontier_deaths_averted <- "Include"
#df1$frontier_deaths_averted[df1$scenario == "Yearly 15-45"] <- "Dominated"
#df1$frontier_ly_saved <- "Include"
#df1$frontier_ly_saved[df1$scenario %in% c("Yearly 15-45", "Yearly 15-30")] <- "Dominated"

freq_by_age_df_summary <- freq_by_age_df %>%
  group_by(scenario, age_group) %>%
  summarise(median_deaths_averted = median(deaths_averted),
            median_interactions = median(total_interactions))

# Plot: interactions vs deaths averted
ggplot(freq_by_age_df) +
  geom_line(aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = freq_by_age_df_summary,
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = freq_by_age_df_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = freq_by_age_df_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~age_group, scales = "free") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
#  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15),
        aspect.ratio = 3.5)


# BCEA package test ----
library(BCEA)
?bcea

interactions_test <- t(interactions_by_2100[,-1])
interactions_test <- cbind(interactions_test[,2]-interactions_test[,1],
                           interactions_test[,3]-interactions_test[,1],
                           interactions_test[,4]-interactions_test[,1])
interactions_test <- cbind(interactions_test, rep(0,183))
outcome_test <- t(deaths_averted_by_2100[,-1])
outcome_test <- cbind(outcome_test[,2]-outcome_test[,1],
                      outcome_test[,3]-outcome_test[,1],
                      outcome_test[,4]-outcome_test[,1])
outcome_test <- cbind(outcome_test, rep(0,183))


test <- bcea(e=outcome_test,c=interactions_test,
             ref=4,
             interventions=c(as.character(interactions_by_2100$scenario[-1]), "No monitoring"),
             Kmax=50000,
             plot=TRUE)

ceef.plot(test, graph="base", relative = FALSE)

#
ceef_interactions_by_2100 <- rbind(
cbind(scenario = "No monitoring",
      out3$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
cbind(scenario = "5-Yearly 15-30",
      monit_out6$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
cbind(scenario = "5-Yearly 15-45",
      monit_out7$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
cbind(scenario = "5-Yearly 15+ (all ages)",
      out5$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
cbind(scenario = "Yearly 15-30",
      monit_out1$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
cbind(scenario = "Yearly 15-45",
      monit_out2$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
cbind(scenario = "Yearly 15+ (all ages)",
      out6$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]))
ceef_interactions <- t(ceef_interactions_by_2100[,-1])

ceef_deaths_averted_by_2100 <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
                                 scenario_objects = list(out3, monit_out6, monit_out7, out5,
                                                         monit_out1, monit_out2, out6),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "no monitoring")

ceef_deaths_averted_by_2100<- filter(ceef_deaths_averted_by_2100, type == "number_averted") %>%
  select(scenario, sim, value) %>%
  spread(key = "scenario", value = "value") %>%
  select(-sim)

test2 <- bcea(e=as.matrix(ceef_deaths_averted_by_2100),c=ceef_interactions,
              ref=1,
              interventions=c(as.character(ceef_interactions_by_2100$scenario)),
              Kmax=50000,
              plot=TRUE)

ceef.plot(test2, graph="base", relative = FALSE)



