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

# Load files (A1) ----

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/repeat_screening_analysis/"

# No monitoring:

# Status quo
out2 <- readRDS(paste0(out_path, "out2_status_quo_180820.rds"))
out2 <- out2[[1]]

# No repeat screening (2020 only)
out3 <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_240920.rds"))
out3 <- out3[[1]]

# With repeat screening in the same age group (15-60)
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

# Re-screening of the same age group compared to status quo of no treatment ----
# HBV deaths averted
deaths_averted_sq_long <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                  scenario_objects = list(out3,
                                                                          out8,
                                                                          out9,
                                                                          out10,
                                                                          out8a,
                                                                          out9a,
                                                                          out10a),
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

# LY saved
ly_gained_sq_screen_long <- plot_ly_gained(counterfactual_object = out2,
                                           scenario_objects = list(out3,
                                                                   out8,
                                                                   out9,
                                                                   out10,
                                                                   out8a,
                                                                   out9a,
                                                                   out10a),
                                           counterfactual_label = "no treatment programme",
                                           x_axis = "screening")

# Deaths averted per resource use
deaths_averted_per_interaction_sq_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(out3,
                                                                             out8,
                                                                             out9,
                                                                             out10,
                                                                             out8a,
                                                                             out9a,
                                                                             out10a),
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

# Re-screening of the same age group compared to a one-off screen in 2020 ----
deaths_averted_screen_long <- plot_hbv_deaths_averted(counterfactual_object = out3,
                                                      scenario_objects = list(out8, out9, out10,
                                                                              out8a, out9a, out10a),
                                                      counterfactual_label = "treatment programme with one-off screening",
                                                      x_axis = "screening")
# By 2100, any repeat screening freq has similar impact. However at earlier timepoints,
# more frequent rescreening is more beneficial. It is mainly the impact of yearly
# screening that declines over time, whereas the 10 and 5 year frequencies maintain about
# the same % of deaths averted in 2050 and 2100.
# Hypothesis: Leveling of the different frequencies over time could be because of increased impact
# of prevention over time.
# Order of magnitude of deaths averted by 2050 and 2100 is approx between 10 and 30%.
ly_gained_screen_long <- plot_ly_gained(counterfactual_object = out3,
                                        scenario_objects = list(out8, out9, out10,
                                                                out8a, out9a, out10a),
                                        counterfactual_label = "treatment programme with one-off screening",
                                        x_axis = "screening")

# Deaths averted per resource use
deaths_averted_per_interaction_screen_long <-
  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
                                                     scenario_objects = list(out8, out9, out10,
                                                                             out8a, out9a, out10a),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "treatment programme with one-off screening",
                                                     x_axis = "screening")
# Note dropped variables (NA) are for screen every 20 and 10 years by 2030 (0 incremental interactions)
# For total interactions, 1 year requires by far the most incremental interactions per death averted,
# whereas the other screening frequencies are fairly similar. This is entirely due to HBsAg tests.

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
