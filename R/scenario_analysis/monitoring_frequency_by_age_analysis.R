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

# Load files ----

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
