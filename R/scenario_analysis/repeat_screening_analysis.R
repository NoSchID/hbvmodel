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

# Load files (A1/E1 (10% coverage)/EA1 (50% screening coverage)) ----

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/repeat_screening_analysis/"
out_path_monit <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

out_path_full_output <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/"

out2_it <- readRDS(paste0(out_path_full_output, "out2_status_quo_repeat_screens_comparison_280121.rds"))
out2_it <- out2_it[[1]]
out3_it <- readRDS(paste0(out_path_full_output, "a1_it_out3_screen_2020_monit_0_repeat_screens_comparison_280121.rds"))
out3_it <- out3_it[[1]]
# Varying number of repeat screening events for a 10-year frequency and without rescreening
# No monitoring:
out8b_it_2030 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_0_screen_10b_2030_repeat_screens_comparison_280121.rds"))
out8b_it_2030 <- out8b_it_2030[[1]]
out8b_it_2040 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_0_screen_10b_2040_repeat_screens_comparison_280121.rds"))
out8b_it_2040 <- out8b_it_2040[[1]]
out8b_it_2050 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_0_screen_10b_2050_repeat_screens_comparison_280121.rds"))
out8b_it_2050 <- out8b_it_2050[[1]]
out8b_it_2060 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_0_screen_10b_2060_repeat_screens_comparison_280121.rds"))
out8b_it_2060 <- out8b_it_2060[[1]]
# With monitoring every 5 years in <45 year olds
monit_out7_it <- readRDS(paste0(out_path_full_output, "a1_it_out3_screen_2020_monit_sim7_repeat_screens_comparison_290121.rds"))
monit_out7_it <- monit_out7_it[[1]]
out8b_it_2030_monit_sim7 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2030_repeat_screens_comparison_290121.rds"))
out8b_it_2030_monit_sim7 <- out8b_it_2030_monit_sim7[[1]]
out8b_it_2040_monit_sim7 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2040_repeat_screens_comparison_290121.rds"))
out8b_it_2040_monit_sim7 <- out8b_it_2040_monit_sim7[[1]]
out8b_it_2050_monit_sim7 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2050_repeat_screens_comparison_290121.rds"))
out8b_it_2050_monit_sim7 <- out8b_it_2050_monit_sim7[[1]]
out8b_it_2060_monit_sim7 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2060_repeat_screens_comparison_290121.rds"))
out8b_it_2060_monit_sim7 <- out8b_it_2060_monit_sim7[[1]]
# With 50% screening coverage and no monitoring
out3_it_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out3_screen_2020_monit_0_cov50_repeat_screens_comparison_290121.rds"))
out3_it_cov50 <- out3_it_cov50[[1]]
out8b_it_2030_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_0_screen_10b_2030_cov50_repeat_screens_comparison_290121.rds"))
out8b_it_2030_cov50 <- out8b_it_2030_cov50[[1]]
out8b_it_2040_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_0_screen_10b_2040_cov50_repeat_screens_comparison_290121.rds"))
out8b_it_2040_cov50 <- out8b_it_2040_cov50[[1]]
out8b_it_2050_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_0_screen_10b_2050_cov50_repeat_screens_comparison_290121.rds"))
out8b_it_2050_cov50 <- out8b_it_2050_cov50[[1]]
out8b_it_2060_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_0_screen_10b_2060_cov50_repeat_screens_comparison_290121.rds"))
out8b_it_2060_cov50 <- out8b_it_2060_cov50[[1]]
# With 50% screening coverage and monitoring every 5 years in <45 year olds
monit_out7_it_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out3_screen_2020_monit_sim7_cov50_repeat_screens_comparison_010221.rds"))
monit_out7_it_cov50 <- monit_out7_it_cov50[[1]]
out8b_it_2030_monit_sim7_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2030_cov50_repeat_screens_comparison_010221.rds"))
out8b_it_2030_monit_sim7_cov50 <- out8b_it_2030_monit_sim7_cov50[[1]]
out8b_it_2040_monit_sim7_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2040_cov50_repeat_screens_comparison_010221.rds"))
out8b_it_2040_monit_sim7_cov50 <- out8b_it_2040_monit_sim7_cov50[[1]]
out8b_it_2050_monit_sim7_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2050_cov50_repeat_screens_comparison_010221.rds"))
out8b_it_2050_monit_sim7_cov50 <- out8b_it_2050_monit_sim7_cov50[[1]]
out8b_it_2060_monit_sim7_cov50 <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2060_cov50_repeat_screens_comparison_010221.rds"))
out8b_it_2060_monit_sim7_cov50 <- out8b_it_2060_monit_sim7_cov50[[1]]

# OLD SIMULATIONS WITHOUT IT:
# No monitoring:

# Status quo
out2 <- readRDS(paste0(out_path, "out2_status_quo_180820.rds"))
out2 <- out2[[1]]

# Screened untreated different cohorts
# No repeat screening
# 90% cov
out1 <- readRDS(paste0(out_path_monit, "a1_out1_status_quo_cohort_240920.rds"))
out1 <- out1[[1]]
# 10% cov
out8b_2020_cov10_sq <- readRDS(paste0(out_path, "e1_out1_screen_2020_status_quo_cohort_cov10_231120.rds"))
out8b_2020_cov10_sq <- out8b_2020_cov10_sq[[1]]
# With repeat screening
# 90% cov
out8b_2030_sq <- readRDS(paste0(out_path, "a1_out1_screen_2030_status_quo_cohort_231120.rds"))
out8b_2030_sq <- out8b_2030_sq[[1]]
out8b_2040_sq <- readRDS(paste0(out_path, "a1_out1_screen_2040_status_quo_cohort_231120.rds"))
out8b_2040_sq <- out8b_2040_sq[[1]]
# 10% cov
out8b_2030_cov10_sq <- readRDS(paste0(out_path, "e1_out1_screen_2030_status_quo_cohort_cov10_231120.rds"))
out8b_2030_cov10_sq <- out8b_2030_cov10_sq[[1]]
out8b_2040_cov10_sq <- readRDS(paste0(out_path, "e1_out1_screen_2040_status_quo_cohort_cov10_231120.rds"))
out8b_2040_cov10_sq <- out8b_2040_cov10_sq[[1]]

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

### VARY NUMBER OF REPEAT SCREENING EVENTS:

# Vary number of repeat screening events for a 10-year frequency and without rescreening
out8b_2030 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2030_071020.rds"))
out8b_2030 <- out8b_2030[[1]]
out8b_2040 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2040_071020.rds"))
out8b_2040 <- out8b_2040[[1]]
out8b_2050 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2050_071020.rds"))
out8b_2050 <- out8b_2050[[1]]
out8b_2060 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2060_071020.rds"))
out8b_2060 <- out8b_2060[[1]]
# 10-yearly monitoring
out4 <- readRDS(paste0(out_path_monit, "a1_out4_screen_2020_monit_10_290920.rds"))  # One-off screen
out4 <- out4[[1]]
out8b_2030_monit10 <- readRDS(paste0(out_path, "a1_out8b_monit_10_screen_10b_2030_201120.rds"))
out8b_2030_monit10 <- out8b_2030_monit10[[1]]
out8b_2040_monit10 <- readRDS(paste0(out_path, "a1_out8b_monit_10_screen_10b_2040_201120.rds"))
out8b_2040_monit10 <- out8b_2040_monit10[[1]]
# 5 yearly monitoring
out5 <- readRDS(paste0(out_path_monit, "a1_out5_screen_2020_monit_5_201020.rds"))  # One-off screen
out5 <- out5[[1]]
out8b_2030_monit5 <- readRDS(paste0(out_path, "a1_out8b_monit_5_screen_10b_2030_201120.rds"))
out8b_2030_monit5 <- out8b_2030_monit5[[1]]
out8b_2040_monit5 <- readRDS(paste0(out_path, "a1_out8b_monit_5_screen_10b_2040_231120.rds"))
out8b_2040_monit5 <- out8b_2040_monit5[[1]]

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
# 10-yearly monitoring
out8b_2020_cov10_monit10 <- readRDS(paste0(out_path, "e1_out4_screen_2020_monit_10_201120.rds"))
out8b_2020_cov10_monit10 <- out8b_2020_cov10_monit10[[1]]
out8b_2030_cov10_monit10 <- readRDS(paste0(out_path, "e1_out8b_monit_10_screen_10b_2030_201120.rds"))
out8b_2030_cov10_monit10 <- out8b_2030_cov10_monit10[[1]]
out8b_2040_cov10_monit10 <- readRDS(paste0(out_path, "e1_out8b_monit_10_screen_10b_2040_201120.rds"))
out8b_2040_cov10_monit10 <- out8b_2040_cov10_monit10[[1]]
# 5-yearly monitoring
out8b_2020_cov10_monit5 <- readRDS(paste0(out_path, "e1_out5_screen_2020_monit_5_201120.rds"))
out8b_2020_cov10_monit5 <- out8b_2020_cov10_monit5[[1]]
out8b_2030_cov10_monit5 <- readRDS(paste0(out_path, "e1_out8b_monit_5_screen_10b_2030_201120.rds"))
out8b_2030_cov10_monit5 <- out8b_2030_cov10_monit5[[1]]
out8b_2040_cov10_monit5 <- readRDS(paste0(out_path, "e1_out8b_monit_5_screen_10b_2040_201120.rds"))
out8b_2040_cov10_monit5 <- out8b_2040_cov10_monit5[[1]]

# 50% screening coverage!
out8b_2020_cov50 <- readRDS(paste0(out_path, "ea1_out3_screen_2020_monit_0_251120.rds"))
out8b_2020_cov50 <- out8b_2020_cov50[[1]]
out8b_2030_cov50 <- readRDS(paste0(out_path, "ea1_out8b_monit_0_screen_10b_2030_251120.rds"))
out8b_2030_cov50 <- out8b_2030_cov50[[1]]
out8b_2040_cov50 <- readRDS(paste0(out_path, "ea1_out8b_monit_0_screen_10b_2040_251120.rds"))
out8b_2040_cov50 <- out8b_2040_cov50[[1]]
# 10-yearly monitoring
out8b_2020_cov50_monit10 <- readRDS(paste0(out_path, "ea1_out4_screen_2020_monit_10_251120.rds"))
out8b_2020_cov50_monit10 <- out8b_2020_cov50_monit10[[1]]
out8b_2030_cov50_monit10 <- readRDS(paste0(out_path, "ea1_out8b_monit_10_screen_10b_2030_251120.rds"))
out8b_2030_cov50_monit10 <- out8b_2030_cov50_monit10[[1]]
out8b_2040_cov50_monit10 <- readRDS(paste0(out_path, "ea1_out8b_monit_10_screen_10b_2040_251120.rds"))
out8b_2040_cov50_monit10 <- out8b_2040_cov50_monit10[[1]]
# 5-yearly monitoring
out8b_2020_cov50_monit5 <- readRDS(paste0(out_path, "ea1_out5_screen_2020_monit_5_251120.rds"))
out8b_2020_cov50_monit5 <- out8b_2020_cov50_monit5[[1]]
out8b_2030_cov50_monit5 <- readRDS(paste0(out_path, "ea1_out8b_monit_5_screen_10b_2030_251120.rds"))
out8b_2030_cov50_monit5 <- out8b_2030_cov50_monit5[[1]]
out8b_2040_cov50_monit5 <- readRDS(paste0(out_path, "ea1_out8b_monit_5_screen_10b_2040_251120.rds"))
out8b_2040_cov50_monit5 <- out8b_2040_cov50_monit5[[1]]

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



## Plot prevalence of untreated carriers and HBV deaths by age after treatment ----

# Full output for 2020 screen without monitoring
out3f <- readRDS(paste0("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/",
                        "a1_out3_screen_2020_monit_0_full_output_131220.rds"))
out3f <- out3f[[1]]

# Number of carriers in 2020
carriers_2020 <- do.call("rbind", lapply(out3f$carriers_female, function(x) x[(which(out3f$time==2020)),]))+
  do.call("rbind", lapply(out3f$carriers_male, function(x) x[(which(out3f$time==2020)),]))


carriers_2020_grouped <- data.frame(age_group = rep(c("Age 15-30", "Age 30-45", "Age 45-65"), each = 183),
                                  value = c(rowSums(carriers_2020[,which(ages==15):which(ages==30-da)]),
                                            rowSums(carriers_2020[,which(ages==30):which(ages==45-da)]),
                                            rowSums(carriers_2020[,which(ages==45):which(ages==65-da)])))

carriers_2020 <- gather(carriers_2020, key = "age", value = "value")
carriers_2020$age <- rep(seq(0,99.5,by=0.5), each = 183)

ggplot(carriers_2020) +
  stat_summary(aes(x=age, y = value), fun=median, geom = "bar")

ggplot(carriers_2020_grouped) +
  stat_summary(aes(x=age_group, y = value), fun=median, geom = "bar")

# Number of carriers in 2030
carriers_2030 <- do.call("rbind", lapply(out3f$carriers_female, function(x) x[(which(out3f$time==2030)),]))+
  do.call("rbind", lapply(out3f$carriers_male, function(x) x[(which(out3f$time==2030)),]))

treated_carriers_2030 <- do.call("rbind", lapply(out3f$treated_carriers_female, function(x) x[(which(out3f$time==2030)),]))+
  do.call("rbind", lapply(out3f$treated_carriers_male, function(x) x[(which(out3f$time==2030)),]))

carriers_2030_grouped <- data.frame(age_group = rep(c("Age 15-30", "Age 30-45", "Age 45-65"), each = 183),
                                    value = c(rowSums(carriers_2030[,which(ages==15):which(ages==30-da)]),
                                              rowSums(carriers_2030[,which(ages==30):which(ages==45-da)]),
                                              rowSums(carriers_2030[,which(ages==45):which(ages==65-da)])))
treated_carriers_2030_grouped <- data.frame(age_group = rep(c("Age 15-30", "Age 30-45", "Age 45-65"), each = 183),
                                    value = c(rowSums(treated_carriers_2030[,which(ages==15):which(ages==30-da)]),
                                              rowSums(treated_carriers_2030[,which(ages==30):which(ages==45-da)]),
                                              rowSums(treated_carriers_2030[,which(ages==45):which(ages==65-da)])))


carriers_2030 <- gather(carriers_2030, key = "age", value = "value")
carriers_2030$age <- rep(seq(0,99.5,by=0.5), each = 183)

treated_carriers_2030 <- gather(treated_carriers_2030, key = "age", value = "value")
treated_carriers_2030$age <- rep(seq(0,99.5,by=0.5), each = 183)

ggplot(carriers_2030) +
  stat_summary(aes(x=age, y = value), fun=median, geom = "bar") +
  stat_summary(data=treated_carriers_2030, aes(x=age, y = value), fun=median,
               geom = "bar", fill="red")

ggplot(carriers_2030_grouped) +
  stat_summary(aes(x=age_group, y = value), fun=median, geom = "bar") +
  stat_summary(data=treated_carriers_2030_grouped,
               aes(x=age_group, y = value), fun=median, geom = "bar", fill="red")

# Combined plot to see decline in numbers
comb_carriers <- rbind(cbind(carriers_2020_grouped, time = "2020"),
                       cbind(carriers_2030_grouped, time = "2030"))

ggplot(comb_carriers) +
  stat_summary(aes(x=age_group, y = value, fill = time), fun=median, geom = "bar",
               position = "dodge")
# Need to look at carriers who can still be targeted in screening (those not in cohort)

## THESIS PLOT ----
## NEW WITH IT TREATED: Plot of repeat screening impact compared to prevalence and treatment need ----

obj_list <- list(out3_it, out8b_it_2030, out8b_it_2040, out8b_it_2050, out8b_it_2060,
                 monit_out7_it, out8b_it_2030_monit_sim7, out8b_it_2040_monit_sim7,
                 out8b_it_2050_monit_sim7, out8b_it_2060_monit_sim7,
                 out3_it_cov50, out8b_it_2030_cov50, out8b_it_2040_cov50,
                 out8b_it_2050_cov50, out8b_it_2060_cov50,
                 monit_out7_it_cov50, out8b_it_2030_monit_sim7_cov50, out8b_it_2040_monit_sim7_cov50,
                 out8b_it_2050_monit_sim7_cov50, out8b_it_2060_monit_sim7_cov50)
hbv_deaths_averted_sq_n <- list()
hbv_deaths_averted_sq_p <- list()

for(i in 1:length(obj_list)) {
  hbv_deaths_averted_sq_n[[i]] <- data.frame(scenario = obj_list[[i]]$cum_hbv_deaths_by_2100$scenario,
                                             type = "number",
                                             out2_it$cum_hbv_deaths_by_2100[,-c(1:3)]-
                                               obj_list[[i]]$cum_hbv_deaths_by_2100[,-c(1:3)])
  hbv_deaths_averted_sq_p[[i]] <- data.frame(scenario = obj_list[[i]]$cum_hbv_deaths_by_2100$scenario,
                                             type = "proportion",
                                             (out2_it$cum_hbv_deaths_by_2100[,-c(1:3)]-
                                               obj_list[[i]]$cum_hbv_deaths_by_2100[,-c(1:3)])/
                                               out2_it$cum_hbv_deaths_by_2100[,-c(1:3)])

}

hbv_deaths_averted_sq <- rbind(do.call("rbind", hbv_deaths_averted_sq_n),
                               do.call("rbind", hbv_deaths_averted_sq_p))
# List of labels to match scenario to combine with df
categories_df <- data.frame(scenario = unique(hbv_deaths_averted_sq$scenario),
                            screening_end_year = rep(c("2020","2030","2040","2050","2060"), times= 4),
                            monitoring = rep(c("No", "5<45", "No", "5<45"), each = 5),
                            screening_coverage = c(rep("90%", 2*5),rep("50%",2*5)))
categories_df$screening_coverage_end_year <- paste(categories_df$screening_coverage,
                                                   categories_df$screening_end_year)
categories_df$screening_end_year_monitoring <- paste(categories_df$screening_end_year,
                                                     categories_df$monitoring)
categories_df$screening_coverage_monitoring <- paste(categories_df$screening_coverage,
                                                   categories_df$monitoring)

hbv_deaths_averted_sq <- left_join(hbv_deaths_averted_sq, categories_df,
                                   by = "scenario")

hbv_deaths_averted_sq <- gather(hbv_deaths_averted_sq, key = "sim", value = "value",
                                -scenario, -type, -screening_end_year,
                                -monitoring, -screening_coverage, -screening_coverage_end_year,
                                -screening_end_year_monitoring,-screening_coverage_monitoring)
hbv_deaths_averted_sq$sim <- gsub("[^0-9]", "", hbv_deaths_averted_sq$sim)

subset(hbv_deaths_averted_sq, type=="proportion") %>%
  group_by(scenario) %>%
  summarise(median=median(value),
            lower = quantile(value, 0.025),
            upper= quantile(value, 0.975))

ggplot(subset(hbv_deaths_averted_sq, type == "proportion")) +
  #geom_boxplot(aes(x=screening_coverage_end_year, y = value, fill = monitoring))
  stat_summary(aes(x=screening_coverage_end_year, y = value, fill = monitoring), fun = "median",
               geom= "col", position = "dodge")

# Extract remaining number of carriers needing treatment AT time of screen (not after)
# Remaining treatment need with treatment provided
timepoints <- rep(seq(2020,2060,10),4)
remaining_treatment_need_n <- list()
remaining_treatment_need_p <- list()

for(i in 1:length(obj_list)) {
    remaining_treatment_need_n[[i]] <- data.frame(scenario = obj_list[[i]]$cum_hbv_deaths_by_2100$scenario,
               at_time = timepoints[i], type ="number",
               remaining_treatment_need=
                 obj_list[[i]]$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[obj_list[[i]]$time==timepoints[i]]+
                 obj_list[[i]]$treatment_eligible_carriers_screened_over_time_15_to_65[obj_list[[i]]$time==timepoints[i]])
    # As a proportion of the total Gambian population
    remaining_treatment_need_p[[i]] <- data.frame(scenario = obj_list[[i]]$cum_hbv_deaths_by_2100$scenario,
                                                  at_time = timepoints[i], type ="proportion",
                                                  remaining_treatment_need=
                                                    (obj_list[[i]]$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[obj_list[[i]]$time==timepoints[i]]+
                                                    obj_list[[i]]$treatment_eligible_carriers_screened_over_time_15_to_65[obj_list[[i]]$time==timepoints[i]])/
                                                    obj_list[[i]]$total_pop_15_to_65_over_time[obj_list[[i]]$time==timepoints[i]])
}

# Remaining/total treatment need if there was no treatment programme (to see effect of vaccination alone)
sq_treatment_need <- rbind(
  data.frame(scenario = "sq", at_time = "2020",type ="number",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2020]),
  data.frame(scenario = "sq", at_time = "2030",type ="number",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2030]),
  data.frame(scenario = "sq", at_time = "2040",type ="number",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2040]),
  data.frame(scenario = "sq", at_time = "2050",type ="number",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2050]),
  data.frame(scenario = "sq", at_time = "2060",type ="number",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2060]),
  data.frame(scenario = "sq", at_time = "2020",type ="proportion",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2020]/
               out2_it$total_pop_15_to_65_over_time[out2_it$time==2020]),
  data.frame(scenario = "sq", at_time = "2030",type ="proportion",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2030]/
               out2_it$total_pop_15_to_65_over_time[out2_it$time==2030]),
  data.frame(scenario = "sq", at_time = "2040",type ="proportion",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2040]/
               out2_it$total_pop_15_to_65_over_time[out2_it$time==2040]),
  data.frame(scenario = "sq", at_time = "2050",type ="proportion",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2050]/
               out2_it$total_pop_15_to_65_over_time[out2_it$time==2050]),
  data.frame(scenario = "sq", at_time = "2060",type ="proportion",
             remaining_treatment_need=out2_it$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it$time==2060]/
               out2_it$total_pop_15_to_65_over_time[out2_it$time==2060])
)


remaining_treatment_need <- rbind(do.call("rbind", remaining_treatment_need_n),
                                  do.call("rbind", remaining_treatment_need_p))
remaining_treatment_need$sim <- as.character(rep(1:183))
sq_treatment_need$sim <- as.character(rep(1:183))

remaining_treatment_need <- left_join(remaining_treatment_need, categories_df,
                                   by = "scenario")

sq_treatment_need$screening_end_year <- sq_treatment_need$at_time
sq_treatment_need$monitoring <- "No monitoring"
sq_treatment_need$screening_coverage <- NA
sq_treatment_need$screening_coverage_end_year <- NA
sq_treatment_need$screening_end_year_monitoring <- sq_treatment_need$at_time
sq_treatment_need$screening_coverage_monitoring <- NA

remaining_treatment_need <- rbind(remaining_treatment_need, sq_treatment_need)

total_df <- full_join(hbv_deaths_averted_sq, remaining_treatment_need,
                      by = c("scenario", "sim", "type",
                             "screening_end_year", "monitoring","screening_coverage",
                             "screening_coverage_end_year", "screening_end_year_monitoring",
                             "screening_coverage_monitoring"))
colnames(total_df)[colnames(total_df)=="value"] <- "deaths_averted"
total_df <- gather(total_df, key = "outcome", value = "value", -scenario,
                   -type, -sim, -screening_end_year, -monitoring, -screening_coverage,
                   -screening_coverage_end_year,-at_time, -screening_end_year_monitoring,
                   -screening_coverage_monitoring)

ggplot(subset(total_df, type == "proportion"& outcome %in% c("deaths_averted",
                                                             "remaining_treatment_need") &
                scenario != "sq"),
       aes(x=screening_coverage_end_year, y = value*100, fill = monitoring)) +
  stat_summary(fun="median", geom="col", colour="black", position = "dodge")+
  #  stat_summary(fun.min= function(x) quantile(x,0.025),
  #               fun.max= function(x) quantile(x,0.975),
  #               geom="errorbar", width = 0.2)+
  facet_wrap(~outcome, scales="free") +
  scale_fill_viridis_d() +
  theme_classic()
# Maybe show total number/proportion needing treatment AND those (like here) who haven't received it yet
# Note how the total proportion of HBV deaths averted does not go over 35%
# This is to a large degree because of no monitoring, but also:
# imperfect treatment effect and imperfect uptake
# Plot shows: cumulative deaths averted by 2100 against duration of repeat screening,
# and remaining treatment need at time of each repeat screen.
# See Fig 1 or 7 for plots: https://www.nature.com/articles/s41467-020-17528-3.pdf
# Fill colours will be split by monitoring status & coverage
# Maybe add remaining treatment need without treatment
# See my other plot for inspiration

# Have lower coverage as a subset of bars
# Need 2 separate plots to overlay these because order varies by outcome
ggplot() +
  stat_summary(data=subset(total_df, type == "number"& outcome %in% c("deaths_averted",
                                                                 "remaining_treatment_need") &
                        screening_coverage=="90%" ),
               aes(x=screening_end_year, y = value, fill = monitoring),
               fun="median", geom="col",colour="black", position = "dodge")+
  stat_summary(data=subset(total_df, type == "number"& outcome %in% c("deaths_averted",
                                                                      "remaining_treatment_need") &
                             screening_coverage=="50%"),
               aes(x=screening_end_year, y = value, fill = monitoring),
               fun="median", geom="col", colour="black", position = "dodge")+
  facet_wrap(~outcome, scales="free") +
  scale_fill_viridis_d() +
  theme_classic()

# Combine with showing vaccination effect alone on treatment need (cross):
ggplot() +
  stat_summary(data=subset(total_df, type == "proportion" & outcome == "remaining_treatment_need" &
                             scenario == "sq"),
               aes(x=screening_end_year, y = value*100), shape = 4, size=3,
               fun="median", geom="point",colour="black")+
  stat_summary(data=subset(total_df, type == "proportion" & outcome == "remaining_treatment_need" &
                             scenario != "sq" & screening_coverage == "90%"),
               aes(x=screening_end_year, y = value*100, fill = monitoring),
               fun="median", geom="col",colour="black", position = "dodge")+
  scale_fill_viridis_d() +
  theme_classic()
# WHO are the unidentified eligible carriers at the end? Those who did not
# take up the assessment/treatment and also new MTCT cases.

# Note that I looked at this with carriers on treatment, but that is not the ideal plot to show that overall treatment need is declining due to
# vaccination, cause by showing carriers on treatment it already includes the effect
# of prolonged life as a result of treatment

# To add: proportion of all treatment eligible carriers on treatment? Need to calculate this
# with reference to SQ scenario.

# MAIN PLOT
# LEFT: Deaths averted - overlay 50% coverage (which should be lighter)

# Remove 2060 timepoint:
total_df <- subset(total_df, screening_end_year != 2060)

p1 <- ggplot() +
  stat_summary(data=subset(total_df, type == "proportion"& outcome =="deaths_averted" & scenario != "sq" &
                             screening_coverage=="90%"),
               aes(x=monitoring, y = value*100,
                   fill = screening_coverage_monitoring,colour="l1"),
               fun="median", geom="bar")+
  stat_summary(data=subset(total_df, type == "proportion"& outcome =="deaths_averted" & scenario != "sq" &
                             screening_coverage=="50%"),
               aes(x=monitoring, y = value*100,
                   fill = screening_coverage_monitoring,colour="l2"),
               fun="median", geom="bar")+
  facet_wrap(~screening_end_year, ncol = 5, strip.position="bottom") +
  scale_fill_manual("Monitoring strategy",
                    values = c("50% No" = "#A180A9",
                               "50% 5<45" = "#90C7C5",
                               "90% No" = "#440154",
                               "90% 5<45" = "#21908C"),
                    breaks = c("90% 5<45", "90% No"),
                    labels =  c("90% 5<45" = "Every 5 years\nin <45 year olds",
                                "90% No" = "No monitoring")) +
  scale_colour_manual("Screening coverage",
                      values = c("l1" = "black",
                                 "l2" = "black"),
                      labels = c("l1" = "50%",
                                 "l2" = "90%")) +
  guides(color = guide_legend(override.aes = list(fill = c("grey80", "black"),
                                                  colour = c("white", "white")),
                                                  order =2),
         fill = guide_legend(order=1)) +
  ylab("Cumulative HBV-related deaths averted (%)") + # averted by 2100
  xlab("End year of screening") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# RIGHT: Remaining treatment need
p2 <- ggplot() +
  stat_summary(data=subset(total_df, type == "proportion"& outcome =="remaining_treatment_need" &
                             scenario != "sq" &
                             screening_coverage=="50%"),
               aes(x=monitoring, y = value*100,
                   fill = screening_coverage_monitoring,colour="l2"),
               fun="median", geom="bar")+
  stat_summary(data=subset(total_df, type == "proportion"& outcome =="remaining_treatment_need" &
                             scenario != "sq" &
                             screening_coverage=="90%"),
               aes(x=monitoring, y = value*100,
                   fill = screening_coverage_monitoring,colour="l1"),
               fun="median", geom="bar") +
  stat_summary(data=subset(total_df, type == "proportion" &
                             outcome == "remaining_treatment_need" & scenario == "sq"),
            aes(x=1.5, y = value*100, shape = "Base case"), size=5,
               fun="median", geom="point", colour="black")+
  facet_wrap(~screening_end_year, ncol = 5, strip.position="bottom") +
  scale_fill_manual(values = c("50% No" = "#A180A9",
                               "50% 5<45" = "#90C7C5",
                               "90% No" = "#440154",
                               "90% 5<45" = "#21908C"),
                    breaks = c("90% 5<45", "90% No"),
                    labels =  c("90% 5<45" = "Screening and treatment,\nmonitor 5-yearly in <45 year olds",
                                "90% No" = "Screening and treatment,\nno monitoring")) +
  scale_colour_manual("Screening coverage",
                      values = c("l1" = "black",
                                 "l2" = "black"),
                      labels = c("l1" = "50%",
                                 "l2" = "90%")) +
  scale_shape_manual("", values=c("Base case" = 18)) +
  guides(color = guide_legend(override.aes = list(fill = c("grey80", "black"),
                                                  colour = c("white", "white")),
                              order =3),
         fill = guide_legend(title=NULL, order=2, reverse = T),
         shape =guide_legend(title=NULL, order=1)) +
  ylab("Unmet treatment need in total population (%)") + # Unmet need in total targeted population at time of screening
  xlab("Year of screening") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,1.15)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = c(0.73, 0.85),
        legend.spacing.y = unit(0, "pt"),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
# It is surprising how few of the eligible carriers are identified by the programme.
# Is this because of losses in the cascade? Or uncertainty?
# => to quite a large degree because of uncertainty in NUMBERS (demography has an
# effect here too). Though related to vaccination showing the PROPORTION has similar
# effect.

library(grid)
p1a <- arrangeGrob(p1, top = textGrob("A", x = unit(0.01, "npc"),
                                                  y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=18)))
p2b <- arrangeGrob(p2, top = textGrob("B", x = unit(0.01, "npc"),
                                      y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))

# THESIS PLOT
#png(file = "repeat_screen_effect.png", width=315, height=135, units = "mm", res=300, pointsize = 0.99)
grid.arrange(p1a,p2b,ncol=2)
#dev.off()
# Cross indicates treatment need decline with vaccination alone.

# With only 90% coverage
p1b <- ggplot() +
  stat_summary(data=subset(total_df, type == "proportion"& outcome =="deaths_averted" & scenario != "sq" &
                             screening_coverage=="90%"),
               aes(x=monitoring, y = value*100,
                   fill = screening_coverage_monitoring,colour="l1"),
               fun="median", geom="bar")+
  facet_wrap(~screening_end_year, ncol = 5, strip.position="bottom") +
  scale_fill_manual("Monitoring strategy",
                    values = c("50% 5<45" = "#A180A9",
                               "50% No" = "#90C7C5",
                               "90% 5<45" = "#440154",
                               "90% No" = "#21908C"),
                    breaks = c("90% 5<45", "90% No"),
                    labels =  c("90% 5<45" = "Every 5 years\nin <45 year olds",
                                "90% No" = "No monitoring")) +
  scale_colour_manual("Screening coverage",
                      values = c("l1" = "black",
                                 "l2" = "black"),
                      labels = c("l1" = "50%",
                                 "l2" = "90%")) +
  guides(color = guide_legend(override.aes = list(fill = c("grey80", "black"),
                                                  colour = c("white", "white")),
                              order =2),
         fill = guide_legend(order=1)) +
  ylab("Cumulative HBV-related deaths\naverted by 2100 (%)") +
  xlab("End year of screening") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

# RIGHT: Remaining treatment need
p2b <- ggplot() +
  stat_summary(data=subset(total_df, type == "proportion"& outcome =="remaining_treatment_need" &
                             scenario != "sq" &
                             screening_coverage=="90%"),
               aes(x=monitoring, y = value*100,
                   fill = screening_coverage_monitoring),
               fun="median", geom="bar", col="black") +
  stat_summary(data=subset(total_df, type == "proportion" & outcome == "remaining_treatment_need" &
                             scenario == "sq"),
               aes(x=1.5, y = value*100), shape = 4, size=3,
               fun="median", geom="point",colour="black")+
  facet_wrap(~screening_end_year, ncol = 5, strip.position="bottom") +
  scale_fill_manual("Monitoring strategy",
                    values = c("50% 5<45" = "#A180A9",
                               "50% No" = "#90C7C5",
                               "90% 5<45" = "#440154",
                               "90% No" = "#21908C"),
                    breaks = c("90% 5<45", "90% No"),
                    labels =  c("90% 5<45" = "Every 5 years\nin <45 year olds",
                                "90% No" = "No monitoring")) +
#  guides(color = guide_legend(override.aes = list(fill = c("grey80", "black"),
#                                                  colour = c("white", "white")),
#                              order =2),
#         fill = guide_legend(order=1)) +
  ylab("Unmet treatment need in total\ntargeted population at time of screening (%)") +
  xlab("Year of screening") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = c(0.75, 0.75),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

grid.arrange(p1b,p2b)

###

# What proportion of treatment eligible carriers are on treatment over time??
total_eligible <- out3_it$treatment_eligible_carriers_undiagnosed_over_time+
  out3_it$treatment_eligible_carriers_screened_over_time+
  out3_it$treated_carriers_over_time
prop_treated <- out3_it$treated_carriers_over_time/total_eligible

prop_treated <- gather(data.frame(prop_treated), key = "sim", value = "value")
prop_treated$time <- rep(out3_it$time, 183)

ggplot(prop_treated) +
  geom_line(aes(x=time, y = value, group = sim), col = "grey50") +
  stat_summary(aes(x=time, y = value), fun = "median", geom = "line", col = "red")

# What proportion of ALL carriers are on treatment over time?
prop_treated_all_carriers <- out3_it$treated_carriers_over_time/out3_it$total_carriers_over_time
# Of all carriers in the population only 8% would be treated in 2020, going down from there
prop_treated_all_carriers <- gather(data.frame(prop_treated_all_carriers), key = "sim", value = "value")
prop_treated_all_carriers$time <- rep(out3_it$time, 183)

quantile(prop_treated_all_carriers[prop_treated_all_carriers$time==2020.5,]$value,
         c(0.5,0.025,0.975))

ggplot(prop_treated_all_carriers) +
  geom_line(aes(x=time, y = value, group = sim), col = "grey50") +
  stat_summary(aes(x=time, y = value), fun = "median", geom = "line", col = "red")
# Actually the proportion of all carriers on treatment remains relatively similar over time!

