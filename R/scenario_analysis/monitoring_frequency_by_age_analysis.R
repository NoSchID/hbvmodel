# What is the optimal monitoring frequency for different age groups?
# Exploring this from a cohort perspective, assuming 15-60 year olds have been screened (A1)

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
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
out3 <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_240920.rds"))
out3 <- out3[[1]]

# Monitoring all age groups
out5 <- readRDS(paste0(out_path, "a1_out5_screen_2020_monit_5_240920.rds"))
out5 <- out5[[1]]
out6 <- readRDS(paste0(out_path, "a1_out6_screen_2020_monit_1_240920.rds"))
out6 <- out6[[1]]

# Monitoring different age groups (yearly)
monit_out1 <- readRDS(paste0(out_path, "a1_monit_out1_240920.rds"))
monit_out1 <- monit_out1[[1]]
monit_out2 <- readRDS(paste0(out_path, "a1_monit_out2_240920.rds"))
monit_out2 <- monit_out2[[1]]
monit_out3 <- readRDS(paste0(out_path, "a1_monit_out3_240920.rds"))
monit_out3 <- monit_out3[[1]]
monit_out4 <- readRDS(paste0(out_path, "a1_monit_out4_240920.rds"))
monit_out4 <- monit_out4[[1]]
monit_out5 <- readRDS(paste0(out_path, "a1_monit_out5_240920.rds"))
monit_out5 <- monit_out5[[1]]

# Monitoring different age groups (every 5 years)
monit_out6 <- readRDS(paste0(out_path, "a1_monit_out6_240920.rds"))
monit_out6 <- monit_out6[[1]]
monit_out7 <- readRDS(paste0(out_path, "a1_monit_out7_240920.rds"))
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

# 1) Yearly monitoring of which age group is most effective, compared to no monitoring? ----

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

# Maybe check extension in age per 10,000 assessments

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


# INCREMENTAL INTERACTIONS IN EACH SCENARIO BY 2100 COMPARED TO NO MONITORING
interactions <- rbind(
cbind(scenario = "screen_2020_monit_5",
      gather(out5$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
       key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_1",
      gather(out6$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim1",
      gather(monit_out1$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim2",
      gather(monit_out2$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim3",
      gather(monit_out3$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim4",
      gather(monit_out4$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim5",
      gather(monit_out5$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim6",
      gather(monit_out6$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim7",
      gather(monit_out7$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim8",
      gather(monit_out8$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim9",
      gather(monit_out9$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value")),
cbind(scenario = "screen_2020_monit_sim10",
      gather(monit_out10$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)],
             key = "sim", value = "value"))
)
levels(interactions$scenario) <- scenario_labels
interactions <- arrange(interactions, scenario)
colnames(interactions)[3] <- "interactions"
interactions$sim <- gsub("[^0-9]", "", interactions$sim)

df <- subset(cohort_deaths_averted_long, type == "number_averted")
df$sim <- gsub("[^0-9]", "", df$sim)
df <- df %>%
  left_join(interactions, by = c("scenario", "sim"))

df_summary <- df %>%
  group_by(scenario) %>%
  summarise(median_deaths = median(value),
            median_int = median(interactions))

# Subset yearly
sub_yearly <- names(scenario_labels)[4:9]
# Subset 5-yearly
sub_5yearly <- names(scenario_labels)[c(3,10:14)]
# Selection
sub_mixed <- names(scenario_labels)[c(3,4,5,8,10,13)]
# Since the best are 5-yearly 45+, 5-yearly 30+ and 5-yearly all ages, chose here <30, >30 and all

ggplot(subset(df, scenario %in% sub_mixed)) +
#  geom_point(aes(x = interactions, y = value, colour = scenario)) +
  stat_ellipse(geom = "polygon",
               mapping=aes(interactions, value, group = scenario,
                           fill= scenario), alpha = 0.3) +
  geom_point(data = subset(df_summary, scenario %in% sub_mixed),
             aes(x = median_int, y = median_deaths, colour = scenario), size = 5) +
#  scale_colour_brewer(palette = "Set3") +
#  scale_fill_brewer(palette = "Set3") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  ylab("Cumulative number of HBV-related deaths averted") +
  xlab("Incremental number of healthcare interactions") +
  theme_bw()
# 95% prediction ellipse (A prediction ellipse is a region for predicting the location of a new observation
# under the assumption that the population is bivariate normal)
# Need to change ellipse not to go below 0!



