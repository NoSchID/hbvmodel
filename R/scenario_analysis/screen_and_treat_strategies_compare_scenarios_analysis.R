# Strategies for screening and treatment simulation - compare sets of assumptions (A, B, C etc)
require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

## Load files ----
scenario_a_full_results <- readRDS(here("output", "screen_and_treat_results",
                                        "scenario_a_full_results.rds"))

scenario_b_full_results <- readRDS(here("output", "screen_and_treat_results",
                                        "scenario_b_full_results.rds"))

scenario_d1_full_results <- readRDS(here("output", "screen_and_treat_results",
                                        "scenario_d1_full_results.rds"))


### Impact of different treatment programmes on the population level ----

# Full dataframe
hbv_deaths_averted_sq <- rbind(
  cbind(scenario_a_full_results$deaths_averted_sq_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_sq_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_sq_long, assumption = "d1"),
  cbind(scenario_a_full_results$deaths_averted_sq_screen_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_sq_screen_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_sq_screen_long, assumption = "d1"))
hbv_deaths_averted_sq$scenario <- gsub("b_", "", hbv_deaths_averted_sq$scenario)
hbv_deaths_averted_sq$scenario <- factor(hbv_deaths_averted_sq$scenario, levels =
                                           c("screen_2020_monit_0", "screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1",
                                             "monit_0_screen_20", "monit_0_screen_10",
                                             "monit_0_screen_5", "monit_0_screen_1"))

# Impact of the basic treatment programme in the different scenarios
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap( ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Shows that:
# - most of the basic treatment impact happens in the short-term.
# - not much difference between the infant vaccine only and PMTCT scale-up scenarios.
# - Including younger age groups (D1) is always more beneficial, especially in the short-term.

# Impact of the different test/treat strategies by scenario assumption
# Plots comparing the impact of: No monitoring or repeat screening, 5-yearly monitoring, 5-yearly screening
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Including younger age groups always slightly better except with 5-yearly screening by 2100, when all strategies avert the same proportion of HBV deaths
# (if you screen often there is no real need to include younger age groups)
# Plot looks similar if looking at yearly intervals.
# Check: Is the benefit of repeated screening reduced if younger age groups are included?

# Comparing repeat screening to monitoring
ggplot(data = subset(hbv_deaths_averted_sq,
                       type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# For A and B, the repeat screen generally looks better than basic programme or with monitoring.
# Once younger ages are included (D1), the impact of repeat screening is about the same as
# that of monitoring especially by 2100. This makes sense as there are fewer people.
# left to identify as carriers. Plot looks similar if looking at yearly intervals.

# Is the impact in terms of LY saved the same?

### Impact of monitoring and screening on the population level ----

# Full dataframe
hbv_deaths_averted <- rbind(
  cbind(scenario_a_full_results$deaths_averted_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_long, assumption = "d1"),
  cbind(scenario_a_full_results$deaths_averted_screen_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_screen_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_screen_long, assumption = "d1"))
hbv_deaths_averted$scenario <- gsub("b_", "", hbv_deaths_averted$scenario)
hbv_deaths_averted$scenario <- factor(hbv_deaths_averted$scenario, levels =
                                        c("screen_2020_monit_10",
                                          "screen_2020_monit_5", "screen_2020_monit_1",
                                          "monit_0_screen_20", "monit_0_screen_10",
                                          "monit_0_screen_5", "monit_0_screen_1"))

### MONITORING

# Impact of the different monitoring strategies by scenario assumption
ggplot(data = subset(hbv_deaths_averted,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Population-level monitoring is more beneficial if younger age groups are included.
# This is true for all intervals and becomes clearer in the medium and long-term.
# Monitoring is a medium to long-term strategy: it has little benefit by 2030 but more thereafter.

# Comparing monitoring frequencies
ggplot(data = subset(hbv_deaths_averted,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# This is not so easy to tell using this plot, as would have to compare e.g yearly
# to 5 yearly monitoring directly.
# For all scenarios + all timepoints, more frequent monitoring is always slightly better.
# It appears that this effect could be stronger for D1 than for A and B, reinforcing
# the notion that the population-level benefit of monitoring increases if younger age groups are included.

### SCREENING

# Impact of the different screening strategies by scenario assumption
ggplot(data = subset(hbv_deaths_averted,
                     type == "proportion_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Repeated screening every 5 or 1 year has similar impact in the short as in the long term.
# Scenarios are similar overall but very frequent screening is slightly less beneficial in the long term
# when younger ages are included (D1)

# Comparing screening frequencies
ggplot(data = subset(hbv_deaths_averted,
                     type == "proportion_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# In the short term, more frequent screening is always more beneficial, but this
# pattern levels off in the medium and long term.
# By 2030, yearly screening looks much better in effect than 5-yearly screening,
# By 2100 screening every 10, 5 and 1 year look very similar, which is clearest for D1.

### IMPACT PER RESOURCE USE

# Full dataframe
hbv_deaths_averted_per_interaction <- rbind(
  cbind(scenario_a_full_results$deaths_averted_per_interaction_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_per_interaction_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_per_interaction_long, assumption = "d1"),
  cbind(scenario_a_full_results$deaths_averted_per_interaction_screen_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_per_interaction_screen_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_per_interaction_screen_long, assumption = "d1"))
hbv_deaths_averted_per_interaction$scenario <- gsub("b_", "", hbv_deaths_averted_per_interaction$scenario)
hbv_deaths_averted_per_interaction$scenario <- factor(hbv_deaths_averted_per_interaction$scenario, levels =
                                        c("screen_2020_monit_10",
                                          "screen_2020_monit_5", "screen_2020_monit_1",
                                          "monit_0_screen_20", "monit_0_screen_10",
                                          "monit_0_screen_5", "monit_0_screen_1"))
# Note that it does not make sense to compare the interactions between monitoring and screening
# as they require very different levels of resources

# MONITORING

# Interactions required per death averted for different monitoring strategies
# Different plots for different timescales due to varying y axis
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                       by_year == 2030 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~scenario) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                     by_year == 2050 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~scenario) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                     by_year == 2100 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~scenario) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# At all timepoints the interactions per death averted are very similar across A, B and D1
# although it tends to be slightly lower for D1.
# Deaths averted per interaction always declines over time.

# Comparing monitoring frequencies
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                     by_year == 2030 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                     by_year == 2050 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                     by_year == 2100 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# In all scenarios and at all timepoints, yearly monitoring always requires far more
# interactions per death averted than 10 or 5 yearly monitoring, which are similar.
# Since this picture is so consistent might be worth looking at something intermediate
# betweem 5 and 1 year.

# SCREENING

# Interactions required per death averted for different screening strategies
# Note varying y scales
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(by_year~scenario, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# For all timeframes and screening frequencies, scenario D1 requires a higher number of interactions
# per death averted (but this is mostly screenings). Makes sense because screening many young people
# that will not die from HBV anytime soon.
# For screening, the number of interactions per death averted increases over time as opposed to
# for monitoring.

# Comparing screening frequencies
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                     by_year == 2030 &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                     by_year == 2050 &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(data = subset(hbv_deaths_averted_per_interaction,
                     by_year == 2100 &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# For all timescales and scenarios, screening every 20, 10 and 5 years is relatively
# similar whereas the number of interactions per death averted becomes much higher
# for yearly screening. The jump between 5 years and 1 year appears biggest for D1.


### Compare scenario results ----

hbv_deaths_averted_sq <- rbind(
  cbind(scenario_a_full_results$deaths_averted_sq_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_sq_long, assumption = "b"))
hbv_deaths_averted_sq$scenario <- gsub("b_", "", hbv_deaths_averted_sq$scenario)
hbv_deaths_averted_sq$scenario <- factor(hbv_deaths_averted_sq$scenario, levels =
                                           c("screen_2020_monit_0", "screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1"))

ggplot(hbv_deaths_averted_sq[hbv_deaths_averted_sq$type == "proportion_averted",]) +
  geom_boxplot(aes(x=scenario, y = value, colour = assumption)) +
  facet_wrap(~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# No difference by 2030 or 2050, but slightly higher proportion of deaths averted for all scenarios by 2100
# with additional BD

hbv_deaths_averted <- rbind(
  cbind(scenario_a_full_results$deaths_averted_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_long, assumption = "b"))
hbv_deaths_averted$scenario <- gsub("b_", "", hbv_deaths_averted$scenario)
hbv_deaths_averted$scenario <- factor(hbv_deaths_averted$scenario, levels =
                                           c("screen_2020_monit_0", "screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1"))

ggplot(hbv_deaths_averted[hbv_deaths_averted$type == "proportion_averted",]) +
  geom_boxplot(aes(x=scenario, y = value, colour = assumption)) +
  facet_wrap(~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# No difference by 2030 or 2050, but slightly higher proportion of deaths averted for all
# monitoring scenarios by 2100 with additional BD

# For HBV deaths averted, monitoring of any frequency seems to get slightly more beneficial in the
# long term if PMTCT is scaled-up

ly_gained_sq <- rbind(
  cbind(scenario_a_full_results$ly_gained_sq_long, assumption = "a"),
  cbind(scenario_b_full_results$ly_gained_sq_long, assumption = "b"))
ly_gained_sq$counterfactual <- gsub("b_", "", ly_gained_sq$counterfactual)
ly_gained_sq$counterfactual <- factor(ly_gained_sq$counterfactual, levels =
                                           c("screen_2020_monit_0", "screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1"))

ggplot(ly_gained_sq[ly_gained_sq$type == "proportion_averted",]) +
  geom_boxplot(aes(x=counterfactual, y = value, colour = assumption)) +
  facet_wrap(~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Impact of treatment with or without monitoring on saving life years appears identical between
# scenarios A and B.

# Screening
hbv_deaths_averted_screen_sq <- rbind(
  cbind(scenario_a_full_results$deaths_averted_sq_screen_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_sq_screen_long, assumption = "b"))
hbv_deaths_averted_screen_sq$scenario <- gsub("b_", "", hbv_deaths_averted_screen_sq$scenario)
hbv_deaths_averted_screen_sq$scenario <- factor(hbv_deaths_averted_screen_sq$scenario, levels =
                                           c("screen_2020_monit_0", "monit_0_screen_20",
                                             "monit_0_screen_10", "monit_0_screen_5",
                                             "monit_0_screen_1"))

ggplot(hbv_deaths_averted_screen_sq[hbv_deaths_averted_screen_sq$type == "proportion_averted",]) +
  geom_boxplot(aes(x=scenario, y = value, colour = assumption)) +
  facet_wrap(~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# No difference by 2030 or 2050, but slightly higher proportion of deaths averted for all scenarios by 2100
# with additional BD


### TO DO: LY, COHORT, ABSOLUTE NUMBERS ----
