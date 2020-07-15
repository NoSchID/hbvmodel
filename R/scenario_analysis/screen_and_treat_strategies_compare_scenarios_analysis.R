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


### Impact of different treatment programmes on the population level (basic, with monitoring, with screening) ----

# Full dataframe
hbv_deaths_averted_sq <- rbind(
  cbind(scenario_a_full_results$deaths_averted_sq_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_sq_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_sq_long, assumption = "d1"),
  cbind(scenario_a_full_results$deaths_averted_sq_screen_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_sq_screen_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_sq_screen_long, assumption = "d1"),
  cbind(scenario_a_full_results$deaths_averted_sq_combi_long, assumption = "a"),
  cbind(scenario_b_full_results$deaths_averted_sq_combi_long, assumption = "b"),
  cbind(scenario_d1_full_results$deaths_averted_sq_combi_long, assumption = "d1"))
hbv_deaths_averted_sq$scenario <- gsub("b_", "", hbv_deaths_averted_sq$scenario)
hbv_deaths_averted_sq$scenario <- factor(hbv_deaths_averted_sq$scenario, levels =
                                           c("screen_2020_monit_0", "screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1",
                                             "monit_0_screen_20", "monit_0_screen_10",
                                             "monit_0_screen_5", "monit_0_screen_1", "monit_5_screen_5"))

ly_gained_sq <- rbind(
  cbind(scenario_a_full_results$ly_gained_sq_long, assumption = "a"),
  cbind(scenario_b_full_results$ly_gained_sq_long, assumption = "b"),
  cbind(scenario_d1_full_results$ly_gained_sq_long, assumption = "d1"),
  cbind(scenario_a_full_results$ly_gained_sq_screen_long, assumption = "a"),
  cbind(scenario_b_full_results$ly_gained_sq_screen_long, assumption = "b"),
  cbind(scenario_d1_full_results$ly_gained_sq_screen_long, assumption = "d1"),
  cbind(scenario_a_full_results$ly_gained_sq_combi_long, assumption = "a"),
  cbind(scenario_b_full_results$ly_gained_sq_combi_long, assumption = "b"),
  cbind(scenario_d1_full_results$ly_gained_sq_combi_long, assumption = "d1"))
colnames(ly_gained_sq)[colnames(ly_gained_sq) %in% c("counterfactual", "scenario")] <- c("scenario", "counterfactual")
ly_gained_sq$scenario <- gsub("b_", "", ly_gained_sq$scenario)
ly_gained_sq$scenario <- factor(ly_gained_sq$scenario, levels =
                                           c("screen_2020_monit_0", "screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1",
                                             "monit_0_screen_20", "monit_0_screen_10",
                                             "monit_0_screen_5", "monit_0_screen_1", "monit_5_screen_5"))

# 1) Impact of the basic treatment programme in the different scenarios
# Deaths proportion
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap( ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Shows that:
# - most of the basic treatment impact in terms of % of deaths averted happens in the short-term.
# - not much difference between the infant vaccine only and PMTCT scale-up scenarios.
# - Including younger age groups (D1) is always more beneficial, especially in the short-term.
# Without monitoring or repeat screening, preferable to include younger age groups.

# Deaths number
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap( ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# In contrast to proportion averted, the number of deaths averted increases  until 2050 (then stays the same)
# D1 still higher but more overlap across the scenarios.

# LY proportion
ggplot(data = subset(ly_gained_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap( ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# - most of the basic treatment impact in terms of % ly saved happens in the medium term (2050)
# and is lowest in the long-term.
# - no difference between the infant vaccine only and PMTCT scale-up scenarios.
# - Including younger age groups (D1) is always more beneficial, especially in the medium-long-term
# (compared to short-term for HBV deaths)
# Without monitoring or repeat screening, preferable to include younger age groups.

# LY number
ggplot(data = subset(ly_gained_sq,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap( ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Number of life years saved keeps increasing over time (even though proportionally this becomes less
# after 2050). Always higher for D1, especially in the long-term.

# 2) Impact of the different test/treat strategies by scenario assumption
# Plots comparing the impact of: No monitoring or repeat screening, 5-yearly monitoring, 5-yearly screening
# Deaths proportion
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5", "monit_5_screen_5"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Including younger age groups always slightly better except with 5-yearly screening by 2100, when all strategies avert the same proportion of HBV deaths
# (if you screen often there is no real need to include younger age groups)
# Plot looks similar if looking at yearly intervals.

# Deaths number
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5", "monit_5_screen_5"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Same pattern as for proportion but with less distinction between the scenarios.

# LY proportion
ggplot(data = subset(ly_gained_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5", "monit_5_screen_5"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# D1 slightly better at all timepoints and for all strategies (slightly different to looking at
# deaths averted). Highest proportion of LY saved always by 2050.

# LY number
ggplot(data = subset(ly_gained_sq,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5", "monit_5_screen_5"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Number of life years saved always increases over time and is always slightly higher for D1.

# 3) Comparing repeat screening vs monitoring vs both
# Deaths proportion
ggplot(data = subset(hbv_deaths_averted_sq,
                       type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5", "monit_5_screen_5"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# For A and B, the repeat screen generally looks better than basic programme or with monitoring.
# Once younger ages are included (D1), the impact of repeat screening is about the same as
# that of monitoring especially by 2100. This makes sense as there are fewer people.
# left to identify as carriers. Plot looks similar if looking at yearly intervals.
# The combination of monitoring and screening looks substantially better by 2050 and especially by 2100
# in all scenarios, but this is most obvious for D1.

# Deaths number
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Broadly similar pattern

# LY proportion
ggplot(data = subset(ly_gained_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# All nearly identical by 2030. In 2050 and 2100, repeat screening is always slightly higher
# than monitoring for all scenarios including for D1. 5 year monitoring appears to have
# very little impact on LY saved for A and B (unless younger ages are included).

# LY number
ggplot(data = subset(ly_gained_sq,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_0", "screen_2020_monit_5",
                                       "monit_0_screen_5"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Same as LY proportion

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

ly_gained <- rbind(
  cbind(scenario_a_full_results$ly_gained_long, assumption = "a"),
  cbind(scenario_b_full_results$ly_gained_long, assumption = "b"),
  cbind(scenario_d1_full_results$ly_gained_long, assumption = "d1"),
  cbind(scenario_a_full_results$ly_gained_screen_long, assumption = "a"),
  cbind(scenario_b_full_results$ly_gained_screen_long, assumption = "b"),
  cbind(scenario_d1_full_results$ly_gained_screen_long, assumption = "d1"))
colnames(ly_gained)[colnames(ly_gained) %in% c("counterfactual", "scenario")] <- c("scenario", "counterfactual")
ly_gained$scenario <- gsub("b_", "", ly_gained$scenario)
ly_gained$scenario <- factor(ly_gained$scenario, levels =
                                        c("screen_2020_monit_10",
                                          "screen_2020_monit_5", "screen_2020_monit_1",
                                          "monit_0_screen_20", "monit_0_screen_10",
                                          "monit_0_screen_5", "monit_0_screen_1"))


### MONITORING

# Impact of the different monitoring strategies by scenario assumption
# Deaths proportion
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

# Deaths number
ggplot(data = subset(hbv_deaths_averted,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Higher for D1 for all intervals and becomes clearer in the medium and long-term.
# Very few extra deaths averted by monitoring by 2030 (around 200 for yearly monitoring).

# LY proportion
ggplot(data = subset(ly_gained,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Same conclusion as deaths proportion: always higher for D1, especially long-term.

# LY number
ggplot(data = subset(ly_gained,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Same conclusion: always higher for D1, especially long-term.

# Comparing monitoring frequencies
# Deaths proportion
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

# Deaths number
ggplot(data = subset(hbv_deaths_averted,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Same as proportion

# LY proportion
ggplot(data = subset(ly_gained,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Same as deaths: for all scenarios + all timepoints, more frequent monitoring is always slightly better.
# But difficult to tell how much that added benefit is in this comparison.

### SCREENING

# Impact of the different screening strategies by scenario assumption
# Deaths proportion
ggplot(data = subset(hbv_deaths_averted,
                     type == "proportion_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Repeated screening every 5 or 1 year has similar prop averted in the short as in the long term.
# Scenarios are similar overall but very frequent screening is slightly less beneficial in the long term
# when younger ages are included (D1)

# Deaths number
ggplot(data = subset(hbv_deaths_averted,
                     type == "number_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Number averted increases over time but otherwise same as prop.

# LY proportion
ggplot(data = subset(ly_gained,
                     type == "proportion_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Proportion of LY saved is largest in 2050 and 2100 and seems the same across scenario assumptions.
# This is different from HBV deaths which suggested that frequent screening is slightly less beneficial
# where younger age groups were included.

# LY number
ggplot(data = subset(ly_gained,
                     type == "number_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Increases over time but seems the same across scenarios (like proportion LY).

# Comparing screening frequencies
# Deaths proportion
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

# Deaths number
ggplot(data = subset(hbv_deaths_averted,
                     type == "number_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Broadly the same conclusion as for prop.

# LY proportion
ggplot(data = subset(ly_gained,
                     type == "proportion_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Like for deaths averted, in the short term more frequent screening is always more beneficial but this
# pattern levels off in the medium and long term. This is much more pronounced for D1.

# LY number
ggplot(data = subset(ly_gained,
                     type == "number_averted" &
                       scenario %in% c("monit_0_screen_20", "monit_0_screen_10",
                                       "monit_0_screen_5", "monit_0_screen_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(by_year ~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Broadly the same

### IMPACT PER RESOURCE USE ----

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

ly_gained_per_interaction <- rbind(
  cbind(scenario_a_full_results$ly_gained_per_interaction_long, assumption = "a"),
  cbind(scenario_b_full_results$ly_gained_per_interaction_long, assumption = "b"),
  cbind(scenario_d1_full_results$ly_gained_per_interaction_long, assumption = "d1"),
  cbind(scenario_a_full_results$ly_gained_per_interaction_screen_long, assumption = "a"),
  cbind(scenario_b_full_results$ly_gained_per_interaction_screen_long, assumption = "b"),
  cbind(scenario_d1_full_results$ly_gained_per_interaction_screen_long, assumption = "d1"))
ly_gained_per_interaction$scenario <- gsub("b_", "", ly_gained_per_interaction$scenario)
ly_gained_per_interaction$scenario <- factor(ly_gained_per_interaction$scenario, levels =
                                                        c("screen_2020_monit_10",
                                                          "screen_2020_monit_5", "screen_2020_monit_1",
                                                          "monit_0_screen_20", "monit_0_screen_10",
                                                          "monit_0_screen_5", "monit_0_screen_1"))

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

# LY GAINED PER INTERACTION
ggplot(data = subset(ly_gained_per_interaction,
                     by_year == 2030 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(data = subset(ly_gained_per_interaction,
                     by_year == 2050 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
ggplot(data = subset(ly_gained_per_interaction,
                     by_year == 2100 &
                       scenario %in% c("screen_2020_monit_10", "screen_2020_monit_5",
                                       "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~ assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))


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


### TO DO: COHORT ----
