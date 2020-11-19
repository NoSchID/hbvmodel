# Population-level impact of repeat screening

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

# Load files (A1/E1) ----

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/repeat_screening_analysis/"

# No monitoring:

# Status quo
out2 <- readRDS(paste0(out_path, "out2_status_quo_180820.rds"))
out2 <- out2[[1]]

# No repeat screening (2020 only)
out3 <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_240920.rds"))
out3 <- out3[[1]]

# With repeat screening in the same age group (15-60) and rescreening of HBsAg-negatives
out8 <- readRDS(paste0(out_path, "a1_out8_monit_0_screen_10_290920.rds"))
out8 <- out8[[1]]
out9 <- readRDS(paste0(out_path, "a1_out9_monit_0_screen_5_290920.rds"))
out9 <- out9[[1]]
out10 <- readRDS(paste0(out_path, "a1_out10_monit_0_screen_1_290920.rds"))
out10 <- out10[[1]]

# With repeat screening in the unscreened age group (depends on repeat screening frequency)
out8a <- readRDS(paste0(out_path, "a1_out8a_monit_0_screen_10a_021020.rds"))
out8a <- out8a[[1]]
out9a <- readRDS(paste0(out_path, "a1_out9a_monit_0_screen_5a_021020.rds"))
out9a <- out9a[[1]]
out10a <- readRDS(paste0(out_path, "a1_out10a_monit_0_screen_1a_021020.rds"))
out10a <- out10a[[1]]

# With repeat screening in the same age group (15-60) but without rescreening of HBsAg-negatives
out8b <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_061020.rds"))
out8b <- out8b[[1]]
out9b <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_061020.rds"))
out9b <- out9b[[1]]
out10b <- readRDS(paste0(out_path, "a1_out10b_monit_0_screen_1b_061020.rds"))
out10b <- out10b[[1]]

# Vary number of repeat screening events for a 10-year frequency and without rescreening
out8b_2030 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2030_071020.rds"))
out8b_2030 <- out8b_2030[[1]]
out8b_2040 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2040_071020.rds"))
out8b_2040 <- out8b_2040[[1]]
out8b_2050 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2050_071020.rds"))
out8b_2050 <- out8b_2050[[1]]
out8b_2060 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2060_071020.rds"))
out8b_2060 <- out8b_2060[[1]]

## Vary number of repeat screening events for a 10-year frequency and without rescreening
# 10% screening coverage!
out8b_2020_cov10 <- readRDS(paste0(out_path, "e1_out3_screen_2020_monit_0_191120.rds"))
out8b_2020_cov10 <- out8b_2020_cov10[[1]]
out8b_2030_cov10 <- readRDS(paste0(out_path, "e1_out8b_monit_0_screen_10b_2030_191120.rds"))
out8b_2030_cov10 <- out8b_2030_cov10[[1]]
out8b_2040_cov10 <- readRDS(paste0(out_path, "e1_out8b_monit_0_screen_10b_2040_191120.rds"))
out8b_2040_cov10 <- out8b_2040_cov10[[1]]
out8b_2050_cov10 <- readRDS(paste0(out_path, "e1_out8b_monit_0_screen_10b_2050_191120.rds"))
out8b_2050_cov10 <- out8b_2050_cov10[[1]]
out8b_2060_cov10 <- readRDS(paste0(out_path, "e1_out8b_monit_0_screen_10b_2060_191120.rds"))
out8b_2060_cov10 <- out8b_2060_cov10[[1]]


# Vary number of repeat screening events for a 5-year frequency and without rescreening
out9b_2030 <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_2030_121020.rds"))
out9b_2030 <- out9b_2030[[1]]
out9b_2040 <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_2040_121020.rds"))
out9b_2040 <- out9b_2040[[1]]
out9b_2050 <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_2050_121020.rds"))
out9b_2050 <- out9b_2050[[1]]
out9b_2060 <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_2060_121020.rds"))
out9b_2060 <- out9b_2060[[1]]

# Access channels analysis
# Workplace screening
out3_wpl <-  readRDS(paste0(out_path, "wpl1_out3_screen_2020_monit_0_301020.rds"))
out3_wpl <- out3_wpl[[1]]
# ANC screening
out3_anc <-  readRDS(paste0(out_path, "anc1_out3_screen_2020_monit_0_021120.rds"))
out3_anc <- out3_anc[[1]]


# Labels
timepoints <- c(2050,2100)
period_labs <- c(paste0("2020-",timepoints[1]), paste0("2020-",timepoints[2]))
names(period_labs) <- c(as.character(timepoints[1]), as.character(timepoints[2]))

## Repeat screening compared to status quo of no treatment ----
# HBV deaths averted
deaths_averted_sq_long <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                  scenario_objects = list(out3,
                                                                          out8,
                                                                          out9,
                                                                          out10,
                                                                          out8a,
                                                                          out9a,
                                                                          out10a,
                                                                          out8b,
                                                                          out9b,
                                                                          out10b),
                                                  counterfactual_label = "no treatment programme",
                                                  x_axis = "screening")
# Repeated screening in the full age group, even if only every 10 years, affects proportion of HBV deaths averted
# particularly in the long-term. In the long term not much difference between 10, 5 or yearly screening,
# although yearly screening appears substantially better than every 5 years to avert deaths by 2030.
# For repeat screening only in the previously untargeted age group, this does not appear to make
# a big difference compared to no repeat screening (much lower than the full age group rescreening,
# and same for all frequencies as you would expect).
# This is likely because of the lack of monitoring (since at this age not many people will be
# identified for treatment).
# No rescreening of HBsAg-negatives leads to the same effect as with rescreening.

# Proportion
ggplot(subset(deaths_averted_sq_long,
              scenario %in% c("screen_2020_monit_0", "monit_0_screen_10b",
                              "monit_0_screen_5b", "monit_0_screen_1b")
              & type=="proportion_averted" & by_year %in% c(2050,2100)),
       aes(scenario, value*100)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  facet_wrap(~by_year, labeller=labeller(by_year = period_labs)) +
  ylab("Percentage of HBV-related deaths averted") +
  scale_x_discrete("Screening frequency", labels =
                     c("screen_2020_monit_0" = "No repeat\nscreen",
                       "monit_0_screen_10b" = "10 years",
                       "monit_0_screen_5b" = "5 years",
                       "monit_0_screen_1b" = "Yearly")) +
  theme_bw() +
  ylim(0,100) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))

# Number
ggplot(subset(deaths_averted_sq_long,
              scenario %in% c("screen_2020_monit_0", "monit_0_screen_10b",
                              "monit_0_screen_5b", "monit_0_screen_1b")
              & type=="number_averted" & by_year %in% c(2050,2100)),
       aes(scenario, value)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  facet_wrap(~by_year, labeller=labeller(by_year = period_labs)) +
  ylab("Number of HBV-related deaths averted") +
  scale_x_discrete("Screening frequency", labels =
                     c("screen_2020_monit_0" = "No repeat\nscreen",
                       "monit_0_screen_10b" = "10 years",
                       "monit_0_screen_5b" = "5 years",
                       "monit_0_screen_1b" = "Yearly")) +
  theme_bw() +
  ylim(0,10000) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))


# LY saved
ly_gained_sq_screen_long <- plot_ly_gained(counterfactual_object = out2,
                                           scenario_objects = list(out3,
                                                                   out8,
                                                                   out9,
                                                                   out10,
                                                                   out8a,
                                                                   out9a,
                                                                   out10a,
                                                                   out8b,
                                                                   out9b,
                                                                   out10b),
                                           counterfactual_label = "no treatment programme",
                                           x_axis = "screening")

# Number
ggplot(subset(ly_gained_sq_screen_long,
              counterfactual %in% c("screen_2020_monit_0", "monit_0_screen_10b",
                              "monit_0_screen_5b", "monit_0_screen_1b")
              & type=="number_averted" & by_year %in% c(2050,2100)),
       aes(counterfactual, value)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  facet_wrap(~by_year, labeller=labeller(by_year = period_labs), scales = "free_y") +
  ylab("Number of life-years saved") +
  scale_x_discrete("Screening frequency", labels =
                     c("screen_2020_monit_0" = "No repeat\nscreen",
                       "monit_0_screen_10b" = "10 years",
                       "monit_0_screen_5b" = "5 years",
                       "monit_0_screen_1b" = "Yearly")) +
  theme_bw() +
 # ylim(0,10000) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))


# Deaths averted per resource use
deaths_averted_per_interaction_sq_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(out3,
                                                                             out8,
                                                                             out9,
                                                                             out10,
                                                                             out8a,
                                                                             out9a,
                                                                             out10a,
                                                                             out8b,
                                                                             out9b,
                                                                             out10b),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "no treatment programme",
                                                     x_axis = "screening")


ggplot(data = deaths_averted_per_interaction_sq_long,
       aes(x=scenario, y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(~by_year, scales = "free_y") +
  ylab("HBV-related deaths averted\nper 10,000 incremental interactions") +
  labs(fill = "Repeat screening \nfrequency", title = paste0("Effect of repeat screening compared to one-off screening")) +
  xlab("Screening frequency (years)") +
  #  ylim(0,100) +
  theme_bw() +
  #  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))
# Repeat screening in only the previously unscreened age groups averts more deaths per interaction
# than repeat screening in all the age groups, but again way fewer than no repeat screening.
# Here the all-age repeat screening with no re-screening of HBsAg-negatives is as good as the new age
# group screening.

ggplot(subset(deaths_averted_per_interaction_sq_long,
              scenario %in% c("screen_2020_monit_0", "monit_0_screen_10b",
                              "monit_0_screen_5b", "monit_0_screen_1b")
              & by_year %in% c(2050,2100)),
       aes(x=scenario, y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5, fill = "grey") +
  facet_wrap(~by_year, labeller=labeller(by_year = period_labs)) +
  ylab("HBV-related deaths averted\nper 10,000 interactions") +
  scale_x_discrete("Screening frequency", labels =
                     c("screen_2020_monit_0" = "No repeat\nscreen",
                       "monit_0_screen_10b" = "10 years",
                       "monit_0_screen_5b" = "5 years",
                       "monit_0_screen_1b" = "Yearly")) +
  theme_bw() +
#  ylim(0,10000) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15))


## Repeat screening compared to a one-off screen in 2020 ----
deaths_averted_screen_long <- plot_hbv_deaths_averted(counterfactual_object = out3,
                                                      scenario_objects = list(out8, out9, out10,
                                                                              out8a, out9a, out10a,
                                                                              out8b, out9b, out10b),
                                                      counterfactual_label = "treatment programme with one-off screening",
                                                      x_axis = "screening")
# By 2100, any repeat screening freq has similar impact. However at earlier timepoints,
# more frequent rescreening is more beneficial. It is mainly the impact of yearly
# screening that declines over time, whereas the 10 and 5 year frequencies maintain about
# the same % of deaths averted in 2050 and 2100.
# Hypothesis: Leveling of the different frequencies over time could be because of increased impact
# of prevention over time.
# Order of magnitude of deaths averted by 2050 and 2100 is approx between 10 and 30%.
deaths_averted_screen_long$scenario <- factor(deaths_averted_screen_long$scenario,
                                              levels =
                                                c("monit_0_screen_10",
                                                 "monit_0_screen_10b",
                                                 "monit_0_screen_10a",
                                                 "monit_0_screen_5",
                                                 "monit_0_screen_5b",
                                                 "monit_0_screen_5a",
                                                 "monit_0_screen_1",
                                                 "monit_0_screen_1b",
                                                 "monit_0_screen_1a"))

deaths_averted_screen_long$freq <- "10 years"
deaths_averted_screen_long$freq[deaths_averted_screen_long$scenario %in%
                                  c("monit_0_screen_5", "monit_0_screen_5a", "monit_0_screen_5b")] <- "5 years"
deaths_averted_screen_long$freq[deaths_averted_screen_long$scenario %in%
                                  c("monit_0_screen_1", "monit_0_screen_1a", "monit_0_screen_1b")] <- "1 year"
deaths_averted_screen_long$freq <- factor(deaths_averted_screen_long$freq, levels = list("10 years" = "10 years",
                                                                                         "5 years" = "5 years",
                                                                                         "1 year" = "1 year"))

# Proportion
ggplot(subset(deaths_averted_screen_long,
              type=="proportion_averted" & by_year %in% c(2050,2100)),
       aes(scenario, value*100, fill = scenario)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(by_year~freq, ncol=3, labeller=labeller(by_year = period_labs), scales = "free_x") +
  ylab("Percentage of HBV-related deaths averted") +
  scale_x_discrete("Repeat screening strategy", labels =
                     c("monit_0_screen_10b" = "B",
                       "monit_0_screen_5b" = "B",
                       "monit_0_screen_1b" = "B",
                       "monit_0_screen_10" = "A",
                       "monit_0_screen_5" = "A",
                       "monit_0_screen_1" = "A",
                       "monit_0_screen_10a" = "C",
                       "monit_0_screen_5a" = "C",
                       "monit_0_screen_1a" = "C")) +
  theme_bw() +
#  ylim(0,40) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Life-years
ly_gained_screen_long <- plot_ly_gained(counterfactual_object = out3,
                                        scenario_objects = list(out8, out9, out10,
                                                                out8a, out9a, out10a,
                                                                out8b, out9b, out10b),
                                        counterfactual_label = "treatment programme with one-off screening",
                                        x_axis = "screening")
ly_gained_screen_long$counterfactual <- factor(ly_gained_screen_long$counterfactual,
                                              levels =
                                                c("monit_0_screen_10",
                                                  "monit_0_screen_10b",
                                                  "monit_0_screen_10a",
                                                  "monit_0_screen_5",
                                                  "monit_0_screen_5b",
                                                  "monit_0_screen_5a",
                                                  "monit_0_screen_1",
                                                  "monit_0_screen_1b",
                                                  "monit_0_screen_1a"))

ly_gained_screen_long$freq <- "10 years"
ly_gained_screen_long$freq[ly_gained_screen_long$counterfactual %in%
                                  c("monit_0_screen_5", "monit_0_screen_5a", "monit_0_screen_5b")] <- "5 years"
ly_gained_screen_long$freq[ly_gained_screen_long$counterfactual %in%
                                  c("monit_0_screen_1", "monit_0_screen_1a", "monit_0_screen_1b")] <- "1 year"
ly_gained_screen_long$freq <- factor(ly_gained_screen_long$freq, levels = list("10 years" = "10 years",
                                                                                         "5 years" = "5 years",
                                                                                         "1 year" = "1 year"))


ggplot(subset(ly_gained_screen_long,
              type=="number_averted" & by_year %in% c(2050,2100)),
       aes(counterfactual, value, fill = counterfactual)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(by_year~freq, ncol=3, labeller=labeller(by_year = period_labs), scales = "free") +
  ylab("Number of life-years saved") +
  scale_x_discrete("Repeat screening strategy", labels =
                     c("monit_0_screen_10b" = "B",
                       "monit_0_screen_5b" = "B",
                       "monit_0_screen_1b" = "B",
                       "monit_0_screen_10" = "A",
                       "monit_0_screen_5" = "A",
                       "monit_0_screen_1" = "A",
                       "monit_0_screen_10a" = "C",
                       "monit_0_screen_5a" = "C",
                       "monit_0_screen_1a" = "C")) +
  theme_bw() +
  #  ylim(0,40) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")


# Deaths averted per resource use
deaths_averted_per_interaction_screen_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects = list(out8, out9, out10,
                                                                             out8a, out9a, out10a,
                                                                             out8b, out9b, out10b),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "treatment programme with one-off screening",
                                                     x_axis = "screening")
# Note dropped variables (NA) are for screen every 20 and 10 years by 2030 (0 incremental interactions)
# For total interactions, 1 year requires by far the most incremental interactions per death averted,
# whereas the other screening frequencies are fairly similar. This is entirely due to HBsAg tests.
deaths_averted_per_interaction_screen_long$scenario <- factor(deaths_averted_per_interaction_screen_long$scenario,
                                              levels =
                                                c("monit_0_screen_10",
                                                  "monit_0_screen_10b",
                                                  "monit_0_screen_10a",
                                                  "monit_0_screen_5",
                                                  "monit_0_screen_5b",
                                                  "monit_0_screen_5a",
                                                  "monit_0_screen_1",
                                                  "monit_0_screen_1b",
                                                  "monit_0_screen_1a"))

deaths_averted_per_interaction_screen_long$freq <- "10 years"
deaths_averted_per_interaction_screen_long$freq[deaths_averted_per_interaction_screen_long$scenario %in%
                                  c("monit_0_screen_5", "monit_0_screen_5a", "monit_0_screen_5b")] <- "5 years"
deaths_averted_per_interaction_screen_long$freq[deaths_averted_per_interaction_screen_long$scenario %in%
                                  c("monit_0_screen_1", "monit_0_screen_1a", "monit_0_screen_1b")] <- "1 year"
deaths_averted_per_interaction_screen_long$freq <- factor(deaths_averted_per_interaction_screen_long$freq, levels = list("10 years" = "10 years",
                                                                                         "5 years" = "5 years",
                                                                                         "1 year" = "1 year"))


ggplot(data = deaths_averted_per_interaction_screen_long,
       aes(x=scenario, y=value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(~by_year, scales = "free_y") +
  ylab("HBV-related deaths averted\nper 10,000 incremental interactions") +
  labs(fill = "Repeat screening \nfrequency", title = paste0("Effect of repeat screening compared to one-off screening")) +
  xlab("Screening frequency (years)") +
#  ylim(0,100) +
  theme_bw() +
  #  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text =element_text(size = 14),
        title = element_text(size = 15))
# Again, no rescreening of HBsAg-negatives has the same impact as with rescreening, but at far fewer
# HBsAg tests.

ggplot(subset(deaths_averted_per_interaction_screen_long,
              by_year %in% c(2050,2100)),
       aes(scenario, value*10000, fill = scenario)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(by_year~freq, ncol=3, labeller=labeller(by_year = period_labs), scales = "free_x") +
  ylab("HBV-related deaths averted\nper 10,000 incremental interactions") +
  scale_x_discrete("Repeat screening strategy", labels =
                     c("monit_0_screen_10b" = "B",
                       "monit_0_screen_5b" = "B",
                       "monit_0_screen_1b" = "B",
                       "monit_0_screen_10" = "A",
                       "monit_0_screen_5" = "A",
                       "monit_0_screen_1" = "A",
                       "monit_0_screen_10a" = "C",
                       "monit_0_screen_5a" = "C",
                       "monit_0_screen_1a" = "C")) +
  theme_bw() +
  #  ylim(0,40) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")


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

## HBV deaths averted and LY saved over time compared to one-off screen----
# Plot of cumulative number of deaths averted over time
deaths_averted_sq_long1 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                   scenario_objects = list(out3, out8b, out9b, out10b),
                                                   counterfactual_label = "no treatment programme",
                                                   timepoints = c(2025,2030,2035))
deaths_averted_sq_long2 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                   scenario_objects = list(out3, out8b, out9b, out10b),
                                                   counterfactual_label = "no treatment programme",
                                                   timepoints = c(2040,2045,2050))
deaths_averted_sq_long3 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                   scenario_objects = list(out3, out8b, out9b, out10b),
                                                   counterfactual_label = "no treatment programme",
                                                   timepoints = c(2055,2060,2065))
deaths_averted_sq_long4 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                   scenario_objects = list(out3, out8b, out9b, out10b),
                                                   counterfactual_label = "no treatment programme",
                                                   timepoints = c(2070,2075,2080))

deaths_averted_sq_long_time <- rbind(deaths_averted_sq_long1, deaths_averted_sq_long2,
                                     deaths_averted_sq_long3, deaths_averted_sq_long4)

deaths_averted_sq_long_time <- deaths_averted_sq_long_time %>%
  filter(type == "number_averted") %>%
  group_by(from_year, by_year, counterfactual, scenario) %>%
  summarise(median = median(value),
            cri_lower = quantile(value, 0.025),
            cri_upper = quantile(value, 0.975))

ly_gained_sq_long1 <- plot_ly_gained(counterfactual_object = out2,
                                     scenario_objects = list(out3, out8b, out9b, out10b),
                                     counterfactual_label = "no treatment programme",
                                     timepoints = c(2025,2030,2035))
ly_gained_sq_long2 <- plot_ly_gained(counterfactual_object = out2,
                                     scenario_objects = list(out3, out8b, out9b, out10b),
                                     counterfactual_label = "no treatment programme",
                                     timepoints = c(2040,2045,2050))
ly_gained_sq_long3 <- plot_ly_gained(counterfactual_object = out2,
                                     scenario_objects = list(out3, out8b, out9b, out10b),
                                     counterfactual_label = "no treatment programme",
                                     timepoints = c(2055,2060,2065))
ly_gained_sq_long4 <- plot_ly_gained(counterfactual_object = out2,
                                     scenario_objects = list(out3, out8b, out9b, out10b),
                                     counterfactual_label = "no treatment programme",
                                     timepoints = c(2070,2075,2080))
ly_gained_sq_long_time <- rbind(ly_gained_sq_long1, ly_gained_sq_long2,
                                ly_gained_sq_long3, ly_gained_sq_long4)

ly_gained_sq_long_time <- ly_gained_sq_long_time %>%
  filter(type == "number_averted") %>%
  group_by(from_year, by_year, counterfactual, scenario) %>%
  summarise(median = median(value),
            cri_lower = quantile(value, 0.025),
            cri_upper = quantile(value, 0.975))

# Plots
p1 <- ggplot(deaths_averted_sq_long_time) +
  geom_line(aes(x = by_year, y = median, group = scenario, colour = scenario)) +
  geom_ribbon(aes(x=by_year, ymin=cri_lower, ymax=cri_upper, group = scenario,
                  fill = scenario), alpha = 0.1)+
  theme_classic() +
  xlab("Time") + ylab("Cumulative number of\nHBV-related deaths averted") +
  theme(legend.position= "none")

p2 <- ggplot(ly_gained_sq_long_time) +
  geom_line(aes(x = by_year, y = median, group = counterfactual, colour = counterfactual)) +
  geom_ribbon(aes(x=by_year, ymin=cri_lower, ymax=cri_upper, group = counterfactual,
                  fill = counterfactual), alpha = 0.1)+
  theme_classic()+
  xlab("Time") + ylab("Cumulative number of\nlife-years saved") +
  theme(legend.position= "bottom")

grid.arrange(p1,p2,ncol =1)
# Would need to add earlier timesteps and indicate timing of screening in plot


## Compare number/duration of repeat screening events for 5/10 year frequency ----
deaths_averted_by_duration_sq_long <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                      scenario_objects = list(out3, out8b_2030, out8b_2040,
                                                                              out8b_2050, out8b_2060, out8b),
                                                      counterfactual_label = "no treatment programme",
                                                      x_axis = "screening")

deaths_averted_by_duration_long <- plot_hbv_deaths_averted(counterfactual_object = out3,
                                                              scenario_objects = list(out8b_2030, out8b_2040,
                                                                                      out8b_2050, out8b_2060, out8b),
                                                              counterfactual_label = "treatment programme with one-off screening",
                                                              x_axis = "screening")

ly_gained_by_duration_sq_long <- plot_ly_gained(counterfactual_object = out2,
                                           scenario_objects = list(out3, out8b_2030, out8b_2040,
                                                                   out8b_2050, out8b_2060, out8b),
                                           counterfactual_label = "no treatment programme",
                                           x_axis = "screening")

deaths_averted_per_interaction_by_duration_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(out3,out8b_2030, out8b_2040,
                                                                             out8b_2050, out8b_2060, out8b),
                                                       interaction_type = "total_interactions",
                                                     counterfactual_label = "treatment programme with one-off screening",
                                                     x_axis = "screening")


ggplot(subset(deaths_averted_per_interaction_by_duration_long,
              by_year %in% 2100),
       aes(scenario, value*10000)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(~by_year, ncol=3, labeller=labeller(by_year = period_labs), scales = "free_x") +
  ylab("HBV-related deaths averted\nper 10,000 incremental interactions") +
#  scale_x_discrete("Repeat screening strategy", labels =
#                     c("monit_0_screen_10b" = "B",
#                       "monit_0_screen_5b" = "B",
#                       "monit_0_screen_1b" = "B",
#                       "monit_0_screen_10" = "A",
#                       "monit_0_screen_5" = "A",
##                       "monit_0_screen_1" = "A",
#                       "monit_0_screen_10a" = "C",
#                       "monit_0_screen_5a" = "C",
#                       "monit_0_screen_1a" = "C")) +
  theme_bw() +
  #  ylim(0,40) +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle=90,hjust=1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

ggplot(subset(deaths_averted_by_duration_long,
              by_year %in% 2100 & type == "proportion_averted"),
       aes(scenario, value)) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge", width = 0.5) +
  facet_wrap(~by_year, ncol=3, labeller=labeller(by_year = period_labs), scales = "free_x") +
  ylab("HBV-related deaths averted\nper 10,000 incremental interactions") +
  #  scale_x_discrete("Repeat screening strategy", labels =
  #                     c("monit_0_screen_10b" = "B",
  #                       "monit_0_screen_5b" = "B",
  #                       "monit_0_screen_1b" = "B",
  #                       "monit_0_screen_10" = "A",
  #                       "monit_0_screen_5" = "A",
  ##                       "monit_0_screen_1" = "A",
  #                       "monit_0_screen_10a" = "C",
  #                       "monit_0_screen_5a" = "C",
  #                       "monit_0_screen_1a" = "C")) +
  theme_bw() +
  ylim(0,0.25) +
  theme(axis.text = element_text(size = 15),
        #axis.text.x = element_text(angle=90,hjust=1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 15),
        legend.position = "none")

# Plot deaths averted against interactions directly
interactions_sq <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        gather(out3$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2030",
        gather(out8b_2030$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2040",
        gather(out8b_2040$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2050",
        gather(out8b_2050$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2060",
        gather(out8b_2060$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b",
        gather(out8b$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")))
#levels(interactions$scenario) <- scenario_labels
#interactions <- arrange(interactions, scenario)
colnames(interactions_sq)[3] <- "interactions"
interactions_sq$sim <- gsub("[^0-9]", "", interactions_sq$sim)

df_deaths_by_interactions <- subset(deaths_averted_by_duration_sq_long, by_year == 2100 & type == "number_averted")
df_deaths_by_interactions$sim <- gsub("[^0-9]", "", df_deaths_by_interactions$sim)
df_deaths_by_interactions <- df_deaths_by_interactions %>%
  left_join(interactions_sq, by = c("scenario", "sim"))

#df_to_plot <- subset(df, scenario %in% sub_mixed)
#levels(df_to_plot$scenario) <-
#  list("Yearly All ages" ="Yearly all ages",
#       "Yearly >30 years"="Yearly 30+",
#       "Yearly <30 years"="Yearly 15-30",
#       "Every 5 years All ages"="5-yearly all ages",
#       "Every 5 years >30 years"="5-yearly 30+",
#       "Every 5 years <30 years"="5-yearly 15-30")

df_deaths_by_interactions_summary <- df_deaths_by_interactions %>%
  group_by(scenario) %>%
  summarise(median_deaths = median(value),
            median_int = median(interactions))

ggplot(df_deaths_by_interactions) +
  stat_ellipse(geom = "polygon",
               aes(interactions, value, group = scenario, fill= scenario), alpha = 0.2) +
  #geom_point(data = df,
  #           aes(x=interactions, y = value, group = scenario, colour= scenario), alpha = 0.5) +
  geom_point(data = df_deaths_by_interactions_summary,
             aes(x = median_int, y = median_deaths, colour = scenario), size = 6) +
  labs(colour = "Number of\nscreening\ncampaigns", fill = "Number of\nscreening\ncampaigns") +
  scale_fill_viridis_d(labels = c("screen_2020_monit_0" = "1 (2020)",
                                  "monit_0_screen_10b_2030" = "2 (2020-2030)",
                                  "monit_0_screen_10b_2040" = "3 (2020-2040)",
                                  "monit_0_screen_10b_2050" = "4 (2020-2050)",
                                  "monit_0_screen_10b_2060" = "5 (2020-2060)",
                                  "monit_0_screen_10b" = "9 (2020-2100)")) +
  scale_colour_viridis_d(labels = c("screen_2020_monit_0" = "1 (2020)",
                                    "monit_0_screen_10b_2030" = "2 (2020-2030)",
                                    "monit_0_screen_10b_2040" = "3 (2020-2040)",
                                    "monit_0_screen_10b_2050" = "4 (2020-2050)",
                                    "monit_0_screen_10b_2060" = "5 (2020-2060)",
                                    "monit_0_screen_10b" = "9 (2020-2100)")) +
  ylab("Cumulative number of\nHBV-related deaths averted") +
  xlab("Incremental number of clinical interactions") +
  labs(title = paste0("Repeat screening strategies (every 10 years)\ncompared to no treatment programme")) +
  theme_bw() +
  xlim(0,8000000) +
  ylim(0,11000) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        title = element_text(size = 15))
# 95% prediction ellipse (A prediction ellipse is a region for predicting the location of a new observation
# under the assumption that the population is bivariate normal)
# Need to change ellipse not to go below 0!

df_ly_by_interactions <- subset(ly_gained_by_duration_sq_long, by_year == 2100 & type == "number_averted")
df_ly_by_interactions$sim <- gsub("[^0-9]", "", df_ly_by_interactions$sim)
colnames(df_ly_by_interactions)[3:4] <- c("scenario", "counterfactual")
df_ly_by_interactions <- df_ly_by_interactions %>%
  left_join(interactions_sq, by = c("scenario", "sim"))

df_ly_by_interactions_summary <- df_ly_by_interactions %>%
  group_by(scenario) %>%
  summarise(median_deaths = median(value),
            median_int = median(interactions))

ggplot(df_ly_by_interactions) +
  stat_ellipse(geom = "polygon",
               aes(interactions, value, group =scenario, fill= scenario), alpha = 0.2) +
  #geom_point(data = df,
  #           aes(x=interactions, y = value, group = scenario, colour= scenario), alpha = 0.5) +
  geom_point(data = df_ly_by_interactions_summary,
             aes(x = median_int, y = median_deaths, colour = scenario), size = 6) +
  labs(colour = "Number of\nscreening\ncampaigns", fill = "Number of\nscreening\ncampaigns") +
  scale_fill_viridis_d(labels = c("screen_2020_monit_0" = "1 (2020)",
                                  "monit_0_screen_10b_2030" = "2 (2020-2030)",
                                  "monit_0_screen_10b_2040" = "3 (2020-2040)",
                                  "monit_0_screen_10b_2050" = "4 (2020-2050)",
                                  "monit_0_screen_10b_2060" = "5 (2020-2060)",
                                  "monit_0_screen_10b" = "9 (2020-2100)")) +
  scale_colour_viridis_d(labels = c("screen_2020_monit_0" = "1 (2020)",
                                    "monit_0_screen_10b_2030" = "2 (2020-2030)",
                                    "monit_0_screen_10b_2040" = "3 (2020-2040)",
                                    "monit_0_screen_10b_2050" = "4 (2020-2050)",
                                    "monit_0_screen_10b_2060" = "5 (2020-2060)",
                                    "monit_0_screen_10b" = "9 (2020-2100)")) +
  ylab("Cumulative number of\nlife-years saved") +
  xlab("Incremental number of clinical interactions") +
  labs(title = paste0("Repeat screening strategies (every 10 years)\ncompared to no treatment programme")) +
  theme_bw() +
  xlim(0,8000000) +
  ylim(0,500000) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        title = element_text(size = 15))

# Add 5 yearly screening frequency
deaths_averted_by_duration_sq_long2 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                              scenario_objects = list(out3, out8b_2030, out8b_2040,
                                                                                      out8b_2050, out8b_2060, out8b,
                                                                                      out9b_2030, out9b_2040,
                                                                                      out9b_2050, out9b_2060, out9b),
                                                              counterfactual_label = "no treatment programme",
                                                              x_axis = "screening")


# Plot deaths averted against interactions directly
interactions_sq2 <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        gather(out3$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2030",
        gather(out8b_2030$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2040",
        gather(out8b_2040$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2050",
        gather(out8b_2050$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2060",
        gather(out8b_2060$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b",
        gather(out8b$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_5b_2030",
        gather(out9b_2030$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_5b_2040",
        gather(out9b_2040$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_5b_2050",
        gather(out9b_2050$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_5b_2060",
        gather(out9b_2060$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_5b",
        gather(out9b$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")))

#levels(interactions$scenario) <- scenario_labels
#interactions <- arrange(interactions, scenario)
colnames(interactions_sq2)[3] <- "interactions"
interactions_sq2$sim <- gsub("[^0-9]", "", interactions_sq2$sim)

df_deaths_by_interactions2 <- subset(deaths_averted_by_duration_sq_long2, by_year == 2100 & type == "number_averted")
df_deaths_by_interactions2$sim <- gsub("[^0-9]", "", df_deaths_by_interactions2$sim)
df_deaths_by_interactions2 <- df_deaths_by_interactions2 %>%
  left_join(interactions_sq2, by = c("scenario", "sim"))

df_deaths_by_interactions_summary2 <- df_deaths_by_interactions2 %>%
  group_by(scenario) %>%
  summarise(median_deaths = median(value),
            median_int = median(interactions))

ggplot(df_deaths_by_interactions2) +
  stat_ellipse(geom = "polygon",
               aes(interactions, value, group = scenario, fill= scenario), alpha = 0.2) +
  #geom_point(data = df,
  #           aes(x=interactions, y = value, group = scenario, colour= scenario), alpha = 0.5) +
  geom_point(data = df_deaths_by_interactions_summary2,
             aes(x = median_int, y = median_deaths, colour = scenario), size = 6) +
  labs(colour = "Number of\nscreening\ncampaigns", fill = "Number of\nscreening\ncampaigns") +
  scale_fill_viridis_d(labels = c("screen_2020_monit_0" = "1 (2020)",
                                  "monit_0_screen_10b_2030" = "2 (2020-2030)",
                                  "monit_0_screen_10b_2040" = "3 (2020-2040)",
                                  "monit_0_screen_10b_2050" = "4 (2020-2050)",
                                  "monit_0_screen_10b_2060" = "5 (2020-2060)",
                                  "monit_0_screen_10b" = "9 (2020-2100)")) +
  scale_colour_viridis_d(labels = c("screen_2020_monit_0" = "1 (2020)",
                                    "monit_0_screen_10b_2030" = "2 (2020-2030)",
                                    "monit_0_screen_10b_2040" = "3 (2020-2040)",
                                    "monit_0_screen_10b_2050" = "4 (2020-2050)",
                                    "monit_0_screen_10b_2060" = "5 (2020-2060)",
                                    "monit_0_screen_10b" = "9 (2020-2100)")) +
  ylab("Cumulative number of\nHBV-related deaths averted") +
  xlab("Incremental number of clinical interactions") +
  labs(title = paste0("Repeat screening strategies (every 10 years)\ncompared to no treatment programme")) +
  theme_bw() +
  xlim(0,8000000) +
  ylim(0,11000) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        title = element_text(size = 15))
# 95% prediction ellipse (A prediction ellipse is a region for predicting the location of a new observation
# under the assumption that the population is bivariate normal)
# Need to change ellipse not to go below 0!


# Add 10% screening coverage (for 10 yearly frequency)
deaths_averted_by_duration_sq_long_cov <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                               scenario_objects = list(out3, out8b_2030, out8b_2040,
                                                                                       out8b_2050, out8b_2060,
                                                                                       out8b_2020_cov10, out8b_2030_cov10,
                                                                                       out8b_2040_cov10,
                                                                                       out8b_2050_cov10, out8b_2060_cov10),
                                                               counterfactual_label = "no treatment programme",
                                                               x_axis = "screening")


# Plot deaths averted against interactions directly
interactions_sq_cov <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        gather(out3$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2030",
        gather(out8b_2030$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2040",
        gather(out8b_2040$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2050",
        gather(out8b_2050$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2060",
        gather(out8b_2060$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "screen_2020_monit_0_cov10",
        gather(out8b_2020_cov10$interactions[[16]]$total_interactions[-c(1:3)],
               key = "sim", value = "value")),
  cbind(scenario = "monit_0_screen_10b_2030_cov10",
        data.frame(sim=rownames(out8b_2030_cov10$interactions[[16]]$total_interactions),
        value=out8b_2030_cov10$interactions[[16]]$total_interactions$`total_screened + total_assessed + total_treated`)),
  cbind(scenario = "monit_0_screen_10b_2040_cov10",
        data.frame(sim=rownames(out8b_2040_cov10$interactions[[16]]$total_interactions),
                   value=out8b_2040_cov10$interactions[[16]]$total_interactions$`total_screened + total_assessed + total_treated`)),
  cbind(scenario = "monit_0_screen_10b_2050_cov10",
        data.frame(sim=rownames(out8b_2050_cov10$interactions[[16]]$total_interactions),
                   value=out8b_2050_cov10$interactions[[16]]$total_interactions$`total_screened + total_assessed + total_treated`)),
  cbind(scenario = "monit_0_screen_10b_2060_cov10",
        data.frame(sim=rownames(out8b_2060_cov10$interactions[[16]]$total_interactions),
                   value=out8b_2060_cov10$interactions[[16]]$total_interactions$`total_screened + total_assessed + total_treated`))
  )

#levels(interactions$scenario) <- scenario_labels
#interactions <- arrange(interactions, scenario)
colnames(interactions_sq_cov)[3] <- "interactions"
interactions_sq_cov$sim <- gsub("[^0-9]", "", interactions_sq_cov$sim)

df_deaths_by_interactions_cov <- subset(deaths_averted_by_duration_sq_long_cov,
                                        by_year == 2100 & type == "number_averted")
df_deaths_by_interactions_cov$sim <- gsub("[^0-9]", "", df_deaths_by_interactions_cov$sim)
df_deaths_by_interactions_cov$scenario <- as.character(df_deaths_by_interactions_cov$scenario)
df_deaths_by_interactions_cov <- df_deaths_by_interactions_cov %>%
  left_join(interactions_sq_cov, by = c("scenario", "sim"))

df_deaths_by_interactions_summary_cov <- df_deaths_by_interactions_cov %>%
  group_by(scenario) %>%
  summarise(median_deaths = median(value),
            median_int = median(interactions))

# Separate by time and coverage
df_deaths_by_interactions_cov$year <- factor(df_deaths_by_interactions_cov$scenario)
levels(df_deaths_by_interactions_cov$year) <- list("2020 90%" = "screen_2020_monit_0",
                                                   "2020 10%"="screen_2020_monit_0_cov10",
                                                   "2030 90%"="monit_0_screen_10b_2030",
                                                   "2030 10%"="monit_0_screen_10b_2030_cov10",
                                                   "2040 90%"="monit_0_screen_10b_2040",
                                                   "2040 10%"= "monit_0_screen_10b_2040_cov10",
                                                   "2050 90%" = "monit_0_screen_10b_2050",
                                                   "2050 10%" = "monit_0_screen_10b_2050_cov10",
                                                   "2060 90%" = "monit_0_screen_10b_2060",
                                                   "2060 10%" = "monit_0_screen_10b_2060_cov10")
df_deaths_by_interactions_cov$year <- as.character(df_deaths_by_interactions_cov$year)
df_deaths_by_interactions_cov <- separate(df_deaths_by_interactions_cov,
                                          col = year, into = c("year", "cov"), sep = " ")
df_deaths_by_interactions_cov$deaths_averted_per_interaction <- df_deaths_by_interactions_cov$value/
  df_deaths_by_interactions_cov$interactions


ggplot(df_deaths_by_interactions_cov) +
  stat_ellipse(geom = "polygon",
               aes(interactions, value, group = scenario, fill= scenario), alpha = 0.2) +
  #geom_point(data = df,
  #           aes(x=interactions, y = value, group = scenario, colour= scenario), alpha = 0.5) +
  geom_point(data = df_deaths_by_interactions_summary_cov,
             aes(x = median_int, y = median_deaths, colour = scenario), size = 6) +
  labs(colour = "Number of\nscreening\ncampaigns", fill = "Number of\nscreening\ncampaigns") +
  scale_fill_viridis_d(labels = c("screen_2020_monit_0" = "1 (2020)",
                                  "monit_0_screen_10b_2030" = "2 (2020-2030)",
                                  "monit_0_screen_10b_2040" = "3 (2020-2040)",
                                  "monit_0_screen_10b_2050" = "4 (2020-2050)",
                                  "monit_0_screen_10b_2060" = "5 (2020-2060)",
                                  "monit_0_screen_10b" = "9 (2020-2100)")) +
  scale_colour_viridis_d(labels = c("screen_2020_monit_0" = "1 (2020)",
                                    "monit_0_screen_10b_2030" = "2 (2020-2030)",
                                    "monit_0_screen_10b_2040" = "3 (2020-2040)",
                                    "monit_0_screen_10b_2050" = "4 (2020-2050)",
                                    "monit_0_screen_10b_2060" = "5 (2020-2060)",
                                    "monit_0_screen_10b" = "9 (2020-2100)")) +
  ylab("Cumulative number of\nHBV-related deaths averted") +
  xlab("Incremental number of clinical interactions") +
  labs(title = paste0("Repeat screening strategies (every 10 years)\ncompared to no treatment programme")) +
  theme_bw() +
  xlim(0,8000000) +
  ylim(0,11000) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        title = element_text(size = 15))

# Plot against time
ggplot(df_deaths_by_interactions_cov) +
  stat_summary(aes(x = year, y = value, group =scenario, colour = cov),
               fun = "median", geom = "point", size = 6) +
  stat_summary(aes(x = year, y = value, group =scenario, colour = cov),
               fun.min = function(x) quantile(x, 0.025),
               fun.max = function(x) quantile(x, 0.975),
               geom = "errorbar", width = 0.5) +
  labs(colour = "Screening\ncoverage") +
  scale_colour_viridis_d(end = 0.9) +
  facet_wrap(~cov, scales = "free") +
  ylab("Cumulative number of\nHBV-related deaths averted") +
  xlab("End year of screening") +
  labs(title = paste0("Repeat screening strategies (every 10 years from 2020)\ncompared to no treatment programme")) +
  theme_bw() +
  expand_limits(y = 0) +
#  xlim(0,8000000) +
#  ylim(0,11000) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 14),
        strip.text = element_blank(),
        legend.text = element_text(size = 13),
        title = element_text(size = 15))

ggplot(df_deaths_by_interactions_cov) +
  stat_summary(aes(x = year, y = deaths_averted_per_interaction*10000, group =scenario, colour = cov),
               fun = "median", geom = "point", size = 6) +
  stat_summary(aes(x = year, y = deaths_averted_per_interaction*10000, group =scenario, colour = cov),
               fun.min = function(x) quantile(x, 0.025),
               fun.max = function(x) quantile(x, 0.975),
               geom = "errorbar", width = 0.5) +
  labs(colour = "Screening\ncoverage") +
  scale_colour_viridis_d(end = 0.9) +
  facet_wrap(~cov, scales = "free") +
  ylab("HBV-related deaths averted\nper 10,000 interactions") +
  xlab("End year of screening") +
  labs(title = paste0("Repeat screening strategies (every 10 years from 2020)\ncompared to no treatment programme")) +
  theme_bw() +
  expand_limits(y = 0) +
  #  xlim(0,8000000) +
  #  ylim(0,11000) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 14),
        strip.text = element_blank(),
        legend.text = element_text(size = 13),
        title = element_text(size = 15))


## Cumulative number of HBsAg tests over time ----
tests <- as.data.frame(sapply(lapply(out10b$interactions, "[[", "total_screened"), "[", 4:186))
tests <- data.frame(t(tests))
tests$year <- seq(2025,2100,by=5)
tests <- gather(tests, key= "sim", value = "value", -year)
tests$value <- as.numeric(unlist(tests$value))

tests_10 <- as.data.frame(sapply(lapply(out8b$interactions, "[[", "total_screened"), "[", 4:186))
tests_10 <- data.frame(t(tests_10))
tests_10$year <- seq(2025,2100,by=5)
tests_10 <- gather(tests_10, key= "sim", value = "value", -year)
tests_10$value <- as.numeric(unlist(tests_10$value))

ggplot() +
  geom_line(data = tests, aes(x= year, y = value, group = sim), col = "grey") +
  geom_line(data = tests_10, aes(x= year, y = value, group = sim), col = "red") +
  theme_bw() +
  ylim(0,8100000)

## Plot deaths averted against prevalence at time of test (testplot) ----
testdf <- data.frame(scenario = df_deaths_by_interactions_summary$scenario[1:5])
testdf$new_deaths_averted_median <- c(df_deaths_by_interactions_summary$median_deaths[1],
                                      diff(df_deaths_by_interactions_summary$median_deaths[1:5]))
testdf$x <- c(2020,2030,2040,2050,2060)
testdf$prev <- apply(out2$timeseries$number_infected[which(out2$timeseries$number_infected$time %in% c(2020,2030,2040,2050,2060)),-c(1:2)],
                     1,median)

plot(x=testdf$prev, y = testdf$new_deaths_averted_median, xlim = c(0,120000))


## Test: One-off population based vs workplace screening ----
plot_hbv_deaths_averted(counterfactual_object = out2,
                        scenario_objects = list(out3,
                                                out3_wpl),
                        counterfactual_label = "no treatment programme",
                        x_axis = "screening")

plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                   scenario_objects = list(out3,
                                                                           out3_wpl),
                                                   interaction_type = "total_interactions",
                                                   counterfactual_label = "no treatment programme",
                                                   x_axis = "screening")

plot_ly_gained(counterfactual_object = out2,
                        scenario_objects = list(out3,
                                                out3_wpl),
                        counterfactual_label = "no treatment programme",
                        x_axis = "screening")

plot_ly_gained_per_healthcare_interaction(counterfactual_object = out2,
                                                   scenario_objects = list(out3,
                                                                           out3_wpl),
                                                   interaction_type = "total_interactions",
                                                   counterfactual_label = "no treatment programme",
                                                   x_axis = "screening")



# Compare this to out5
d1 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                 scenario_objects = list(out3),
                                                 outcome_to_plot = "number_averted",
                                                 counterfactual_label = "no treatment")
d2 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                     scenario_objects = list(out3_wpl),
                                                     outcome_to_plot = "number_averted",
                                                     counterfactual_label = "no treatment")

i1 <- out3$interactions[[16]]$total_interactions[,-c(1:3)]
i2<- out3_wpl$interactions[[16]]$total_interactions[,-c(1:3)]

deaths_df <- cbind(d1$value[d1$by_year==2100 &
                              d1$type == "number_averted"],
                   d2$value[d2$by_year==2100 &
                              d2$type == "number_averted"],
                   rep(0,183))
int_df <- cbind(unlist(i1), unlist(i2),rep(0,183))


ceef.plot(bcea(e=deaths_df,
               c=int_df,ref=3,interventions=c("General population", "Workplace", "No treatment")),graph="base")


# TEST: 1 REPEAT SCREEN VS ALL AGE MONITORING AND ACCESS CHANNELS ----
out_path2 <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

# No monitoring and 5-yearly monitoring
out3 <- readRDS(paste0(out_path2, "a1_out3_screen_2020_monit_0_201020.rds"))
out3 <- out3[[1]]
out5 <- readRDS(paste0(out_path2, "a1_out5_screen_2020_monit_5_201020.rds"))
out5 <- out5[[1]]

# Compare this to out5
monit_scenario_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                 scenario_objects = list(out5),
                                                 outcome_to_plot = "number_averted",
                                                 counterfactual_label = "no treatment")
screen_scenario_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                  scenario_objects = list(out3),
                                                  outcome_to_plot = "number_averted",
                                                  counterfactual_label = "no treatment")
repscreen_scenario_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                     scenario_objects = list(out8b_2030),
                                                     outcome_to_plot = "number_averted",
                                                     counterfactual_label = "no treatment")

wpl_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                              scenario_objects = list(out3_wpl),
                              outcome_to_plot = "number_averted",
                              counterfactual_label = "no treatment")

anc_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                      scenario_objects = list(out3_anc),
                                      outcome_to_plot = "number_averted",
                                      counterfactual_label = "no treatment")


monit_scenario_int <- out5$interactions[[16]]$total_interactions[,-c(1:3)]
screen_scenario_int <- out3$interactions[[16]]$total_interactions[,-c(1:3)]
repscreen_scenario_int <- out8b_2030$interactions[[16]]$total_interactions[,-c(1:3)]
wpl_int <- out3_wpl$interactions[[16]]$total_interactions[,-c(1:3)]
anc_int <- out3_anc$interactions[[16]]$total_interactions[,-c(1:3)]

# NOTE monit_scenario is only 1 with monitoring

deaths_df <- cbind(monit_scenario_deaths$value[monit_scenario_deaths$by_year==2100 &
                                                 monit_scenario_deaths$type == "number_averted"],
                   screen_scenario_deaths$value[screen_scenario_deaths$by_year==2100 &
                                                  screen_scenario_deaths$type == "number_averted"],
                   repscreen_scenario_deaths$value[repscreen_scenario_deaths$by_year==2100 &
                                                     repscreen_scenario_deaths$type == "number_averted"],
                   wpl_deaths$value[wpl_deaths$by_year==2100 &
                                      wpl_deaths$type == "number_averted"],
                   anc_deaths$value[anc_deaths$by_year==2100 &
                                                     anc_deaths$type == "number_averted"],
                   rep(0,183))
int_df <- cbind(unlist(monit_scenario_int), unlist(screen_scenario_int), unlist(repscreen_scenario_int),
                unlist(wpl_int), unlist(anc_int),
                rep(0,183))

scenario_labels2 <- c("2020 population screen,\n5-yearly monitoring",
                     "2020 population screen,\nno monitoring",
                     "2020+2030 population screen,\nno monitoring",
                     "2020 workplace screen,\nno monitoring",
                     "2020 antenatal screen,\nno monitoring",
                     "No treatment")

ceef.plot.median(bcea(e=deaths_df,c=int_df,ref=ncol(deaths_df),
               interventions=scenario_labels2),
          graph="base")

deaths_df <- data.frame(deaths_df)
colnames(deaths_df) <- scenario_labels2
deaths_df$sim <- rownames(deaths_df)
deaths_df <- gather(deaths_df, key = "scenario", value = "deaths_averted", -sim)

int_df <- data.frame(int_df)
colnames(int_df) <- scenario_labels2
int_df$sim <- as.character(seq(1:183))
int_df <- gather(int_df, key = "scenario", value = "interactions", -sim)

combi <- left_join(deaths_df, int_df, by = c("scenario", "sim"))
# Checked that simulation order is the same

combi$frontier <- "Include"
combi$frontier[combi$scenario %in% c("2020 population screen,\nno monitoring",
                                     "2020+2030 population screen,\nno monitoring",
                                     "2020 antenatal screen,\nno monitoring")] <- "Dominated"

combi_summary <- combi %>%
  group_by(scenario, frontier) %>%
  summarise(median_deaths_averted = median(deaths_averted),
            median_interactions = median(interactions))

# Plot
ggplot(combi) +
  geom_line(data= subset(combi, frontier== "Include"),
            aes(x = deaths_averted, y= interactions,
                group = sim), colour = "grey", alpha = 0.3) +
  stat_ellipse(data=subset(combi, scenario != "No treatment"),
               aes(x=deaths_averted,y=interactions,
                   group = reorder(scenario, deaths_averted),
                   fill= reorder(scenario, deaths_averted)),
               geom = "polygon", na.rm = FALSE, alpha = 0.2) +
  geom_point(aes(x = deaths_averted, y = interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)),
             alpha = 0.15) +
  # Overlay median
  geom_line(data = subset(combi_summary, frontier == "Include"),
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = combi_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = combi_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual("Screening scenarios",
                    values=rev(brewer.pal(5,"RdYlBu"))) +
  scale_colour_manual("Screening scenarios",
                      values=c("black", rev(brewer.pal(5,"RdYlBu")))) +
  guides(fill=FALSE,
         colour=guide_legend(keywidth=0.1,keyheight=0.5,default.unit="inch")) +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  #  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

###

# In out5, average time on treatment is 22.3 years - assume the same for repeat screen

monit_scenario_cost <- out5$interactions[[16]]$total_screened[,-c(1:3)]*10.38+
  out3$interactions[[16]]$total_assessed[,-c(1:3)]*120+
  (out5$interactions[[16]]$total_assessed[,-c(1:3)]-out3$interactions[[16]]$total_assessed[,-c(1:3)])*15.77+
  out3$interactions[[16]]$total_treated[,-c(1:3)]*22.3*84.88

repscreen_scenario_cost <- out8b_2030$interactions[[16]]$total_screened[,-c(1:3)]*10.38+
  out8b_2030$interactions[[16]]$total_assessed[,-c(1:3)]*120+
  out8b_2030$interactions[[16]]$total_treated[,-c(1:3)]*22.3*84.88

cost_df <- cbind(unlist(monit_scenario_cost), unlist(repscreen_scenario_cost),rep(0,183))

ceef.plot(bcea(e=deaths_df,c=cost_df,ref=3,interventions=c("Monitoring", "Repeat screen", "No treatment")),graph="base")


