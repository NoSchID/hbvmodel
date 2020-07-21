# Strategies for screening and treatment simulation - compare sets of assumptions with different
# age groups being screened and potential modifiers

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

## Load files ----
# A = optimal coverage, ages 30-70
scenario_a_full_results <- readRDS(here("output", "screen_and_treat_results",
                                        "scenario_a_full_results.rds"))
# B = A with PMTCT scale up
#scenario_b_full_results <- readRDS(here("output", "screen_and_treat_results",
#                                        "scenario_b_full_results.rds"))
# D1 = optimal coverage, ages 15-65
scenario_d1_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_d1_full_results.rds"))
# D2 = optimal coverage, ages 45-70
scenario_d2_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_d2_full_results.rds"))
# D3 = optimal coverage, ages 15-45
scenario_d3_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_d3_full_results.rds"))

# E = Low screening coverage, ages 30-70
scenario_e_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_e_full_results.rds"))

# E1 = Low screening coverage, ages 15-65
scenario_e1_full_results <- readRDS(here("output", "screen_and_treat_results",
                                        "scenario_e1_full_results.rds"))

# E2 = Low screening coverage, ages 45-70
scenario_e2_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_e2_full_results.rds"))

# E3 = Low screening coverage, ages 15-45
scenario_e3_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_e3_full_results.rds"))

## Population-level outcomes of the treatment programme (without/with monitoring) ----

# Full dataframes of HBV deaths averted and LY saved compared to infant vaccine only
hbv_deaths_averted_sq <- rbind(
  cbind(scenario_a_full_results$deaths_averted_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_sq_long, assumption = "d3"),
  cbind(scenario_e_full_results$deaths_averted_sq_long, assumption = "e"),
  cbind(scenario_e1_full_results$deaths_averted_sq_long, assumption = "e1"),
  cbind(scenario_e2_full_results$deaths_averted_sq_long, assumption = "e2"),
  cbind(scenario_e3_full_results$deaths_averted_sq_long, assumption = "e3"))
#hbv_deaths_averted_sq$scenario <- gsub("b_", "", hbv_deaths_averted_sq$scenario)
hbv_deaths_averted_sq$scenario <- factor(hbv_deaths_averted_sq$scenario, levels =
                                           c("screen_2020_monit_0", "screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1"))
hbv_deaths_averted_sq$screening_coverage <- "Optimal"
hbv_deaths_averted_sq$screening_coverage[hbv_deaths_averted_sq$assumption %in% c("e", "e1", "e2", "e3")] <-
  "Low"

ly_gained_sq <- rbind(
  cbind(scenario_a_full_results$ly_gained_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$ly_gained_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$ly_gained_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$ly_gained_sq_long, assumption = "d3"),
  cbind(scenario_e_full_results$ly_gained_sq_long, assumption = "e"),
  cbind(scenario_e1_full_results$ly_gained_sq_long, assumption = "e1"),
  cbind(scenario_e2_full_results$ly_gained_sq_long, assumption = "e2"),
  cbind(scenario_e3_full_results$ly_gained_sq_long, assumption = "e3"))
colnames(ly_gained_sq)[colnames(ly_gained_sq) %in% c("counterfactual", "scenario")] <- c("scenario", "counterfactual")
#ly_gained_sq$scenario <- gsub("b_", "", ly_gained_sq$scenario)
ly_gained_sq$scenario <- factor(ly_gained_sq$scenario, levels =
                                  c("screen_2020_monit_0", "screen_2020_monit_10",
                                    "screen_2020_monit_5", "screen_2020_monit_1"))
ly_gained_sq$screening_coverage <- "Optimal"
ly_gained_sq$screening_coverage[ly_gained_sq$assumption %in% c("e", "e1", "e2", "e3")] <-
  "Low"


# BASIC PROGRAMME (NO MONITORING)
# Deaths proportion
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(screening_coverage ~ by_year, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Most of the basic treatment impact in terms of % of deaths averted happens in the short-term.
# 15-65 years (D1) averts the largest proportion of deaths at every timepoint. 30-70 years (A) and
# 15-45 years (D3) are similar. Old age groups only - 45-70 (D2) - is far worse than all the others.

# In contrast to proportion averted, the number of deaths averted increases  until 2050 (then stays the same)
# Pattern between scenarios still the same but more overlap between bars.

# Looks exactly the same if the screening coverage is lowered to 10%, and the % or number averted appears to be
# about 1/10th of those with optimal coverage. Suggests that the impact of the basic programme
# on averting HBV deaths and saving LY behaves linearly in relation to screening coverage, but would
# need a few more datapoints.

# LY proportion
ggplot(data = subset(ly_gained_sq,
                     type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(screening_coverage ~ by_year, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Most of the basic treatment impact in terms of % ly saved happens in the medium term (2050)
# Pattern broadly the same except that by 2100, D3 is closer to D1 (young age groups only
# looking better than the 30-70 years).

# Number of life years saved keeps increasing over time (even though proportionally this becomes less
# after 2050). Pattern the same as for LY proportion (D3 and D1 very similar in 2100).

# Effect of the lower screening coverage on LY saved appears the same as for deaths averted.

# WITH MONITORING
# Deaths proportion
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "proportion_averted" &
                       assumption %in% c("a", "d1", "d2", "d3"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(assumption ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# In all scenarios by 2030, nearly all the benefit in averting HBV deaths appears to come
# from the initial programme. However in the long-term the monitoring
# appears to become more beneficial.

# LY proportion
ggplot(data = subset(ly_gained_sq,
                     type == "proportion_averted")) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(assumption ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# For proportion of LY saved, early all the benefit seems to come from the initial
# programme at all timepoints.

# TRADE OFF BETWEEN MONITORING AND AGE GROUPS SCREENED(e.g. screen all without monitoring,
# or screen only young with monitoring)
# Look at where monitoring is most useful to chose comparison
# Can do the same with coverage

# Prop deaths averted with yearly monitoring:
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "proportion_averted" &
                       ((scenario == "screen_2020_monit_0" & assumption == "d1") |
                       (scenario == "screen_2020_monit_1" & assumption == "a") |
                       (scenario == "screen_2020_monit_1" & assumption == "d3")))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(~by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Suggests that, by 2030, including a broader age range in the initial screen (15-65 years) and never monitoring
# is better than including a more limited age range (30-70 years or 15-45 years) and monitoring yearly.
# However, by 2050, this changes: yearly monitoring for the restricted age groups looks slightly
# better than the one-off broad screen (although fairly similar). Need to think which one here
# requires more resources - those involved in monitoring are more difficult than one-off screen,
# so the small extra benefit here may not warrant the difficult follow-up.
# In 2100, this benefit seems to increase and is now also higher for D3 than for A.
# D2 already excluded from plot due to extremely low benefit.

# Prop deaths averted with 5-yearly monitoring
ggplot(data = subset(hbv_deaths_averted_sq,
                     type == "proportion_averted" &
                       ((scenario == "screen_2020_monit_0" & assumption == "d1") |
                          (scenario == "screen_2020_monit_5" & assumption == "a") |
                          (scenario == "screen_2020_monit_5" & assumption == "d3")))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(~by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# If monitoring is only every 5 years rather than yearly, the 2050 prop of deaths averted
# is basically identical for all scenarios and a small advantage to monitoring only appears by 2100.

# Prop LY saved with yearly monitoring:
ggplot(data = subset(ly_gained_sq,
                     type == "proportion_averted" &
                       ((scenario == "screen_2020_monit_0" & assumption == "d1") |
                          (scenario == "screen_2020_monit_1" & assumption == "a") |
                          (scenario == "screen_2020_monit_1" & assumption == "d3")))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(~by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# For LY saved, 30-70 year age group with yearly monitoring always stays below 15-65 without monitoring.
# Only the 15-45 with yearly monitoring is slighlty better than D1 by 2100.

# Prop LY saved with 5-yearly monitoring:
ggplot(data = subset(ly_gained_sq,
                     type == "proportion_averted" &
                       ((scenario == "screen_2020_monit_0" & assumption == "d1") |
                          (scenario == "screen_2020_monit_5" & assumption == "a") |
                          (scenario == "screen_2020_monit_5" & assumption == "d3")))) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(~by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Effect for 5-yearly monitoring appears to be the same as yearly. Also very similar for 10-yearly.

# Results overall may suggest that focusing on the initial included age range (and uptake?)
# may be far more important than monitoring. But still want to show an ideal scenario for
# countries who can implement it.

# TRADE OFF BETWEEN SCREENING COVERAGE AND AGE GROUPS SCREENED?
# Probably to some degree but it would still have to be much higher than 10%.
# But could possibly use this to show whether screening uptake is more important than retention.

## Outcomes of basic treatment programme in relation to resource use ----

# Full dataframe of HBV deaths averted per healthcare interaction
hbv_deaths_averted_per_interaction_sq <- rbind(
  cbind(scenario_a_full_results$deaths_averted_per_interaction_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_per_interaction_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_per_interaction_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_per_interaction_sq_long, assumption = "d3"),
  cbind(scenario_e_full_results$deaths_averted_per_interaction_sq_long, assumption = "e"),
  cbind(scenario_e1_full_results$deaths_averted_per_interaction_sq_long, assumption = "e1"),
  cbind(scenario_e2_full_results$deaths_averted_per_interaction_sq_long, assumption = "e2"),
  cbind(scenario_e3_full_results$deaths_averted_per_interaction_sq_long, assumption = "e3"),
  cbind(scenario_a_full_results$deaths_averted_per_assessment_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_per_assessment_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_per_assessment_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_per_assessment_sq_long, assumption = "d3"),
  cbind(scenario_e_full_results$deaths_averted_per_assessment_sq_long, assumption = "e"),
  cbind(scenario_e1_full_results$deaths_averted_per_assessment_sq_long, assumption = "e1"),
  cbind(scenario_e2_full_results$deaths_averted_per_assessment_sq_long, assumption = "e2"),
  cbind(scenario_e3_full_results$deaths_averted_per_assessment_sq_long, assumption = "e3"),
  cbind(scenario_a_full_results$deaths_averted_per_treatment_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_per_treatment_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_per_treatment_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_per_treatment_sq_long, assumption = "d3"),
  cbind(scenario_e_full_results$deaths_averted_per_treatment_sq_long, assumption = "e"),
  cbind(scenario_e1_full_results$deaths_averted_per_treatment_sq_long, assumption = "e1"),
  cbind(scenario_e2_full_results$deaths_averted_per_treatment_sq_long, assumption = "e2"),
  cbind(scenario_e3_full_results$deaths_averted_per_treatment_sq_long, assumption = "e3"))
#hbv_deaths_averted_per_interaction$scenario <- gsub("b_", "", hbv_deaths_averted_per_interaction$scenario)
hbv_deaths_averted_per_interaction_sq$scenario <- factor(hbv_deaths_averted_per_interaction_sq$scenario, levels =
                                                        c("screen_2020_monit_0", "screen_2020_monit_10",
                                                          "screen_2020_monit_5", "screen_2020_monit_1"))
hbv_deaths_averted_per_interaction_sq$screening_coverage <- "Optimal"
hbv_deaths_averted_per_interaction_sq$screening_coverage[
  hbv_deaths_averted_per_interaction_sq$assumption %in% c("e", "e1", "e2", "e3")] <-
  "Low"

# Total interactions
ggplot(data = subset(hbv_deaths_averted_per_interaction_sq,
                     interaction_type == "total_interactions" &
                      scenario == c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# At all timepoints, interactions per death averted are similar for A and D2, and for D1 and D3.
# It's highest for those scenarios involving younger age groups, highest for D3 by 2030.
# Despite this, by 2050 they are fairly similar:
# Assumption A: 217 (110-409)
# Assumption D1: 338 (174-648)
# Judging this would depend on what type these interactions are.
# Effect of 10% screening coverage: interactions per death averted seem to be exactly the same
# irrespective of the screening coverage.

# Comparing this to overall impact:
# D1 best, D2 worst, but A and D3 were not far off from D1 (especially D3 for life years)
# For resource use, A and D2 are best, and D3 worst in the short-term only.
# Another question is, how many resources are used in total?

# Total assessments
ggplot(data = subset(hbv_deaths_averted_per_interaction_sq,
                     interaction_type == "total_assessed" &
                       scenario == c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Total assessments per HBV death averted a tiny bit higher for D3 and D1 (by 2030 at least).
# By 2050 looks nearly identical for all. Unaffected by lower screening coverage.

# Total treatment initiations
ggplot(data = subset(hbv_deaths_averted_per_interaction_sq,
                     interaction_type == "total_treated" &
                       scenario == c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# By 2050 and 2100, treatment initiations per death averted are a tiny bit higher for D2.
# Again unaffected by low screening coverage.

# NEED TO ADD LIFE YEARS SAVED!!!

# TRADE OFF BETWEEN MONITORING AND AGE GROUPS SCREENED
# Here also need to look at the interaction type
ggplot(data = subset(hbv_deaths_averted_per_interaction_sq,
                       ((scenario == "screen_2020_monit_0" & assumption == "d1") |
                          (scenario == "screen_2020_monit_1" & assumption == "a") |
                          (scenario == "screen_2020_monit_1" & assumption == "d3")))) +
  geom_boxplot(aes(x=scenario, y = 1/value, fill = assumption)) +
  facet_wrap(~by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Interactions per death averted is similar for D1 without monitoring as for A with yearly monitoring
# Higher for D3 with yearly monitoring. But also note interactions are not equally costly
# so D1 likely best both in terms of overall impact and resource use.
# If monitoring is only 5 yearly, D3 becomes more similar to D1 and A has slightly lower
# interactions required per life-year saved.


## To add: Total resource use! This will vary by coverage etc. ----
## Cohort-level outcomes of the basic treatment programme (no monitoring) ----

# Full dataframes of HBV deaths averted and LY saved compared to infant vaccine only
cohort_hbv_deaths_averted_sq <- rbind(
 cbind(scenario_a_full_results$cohort_deaths_averted_sq_long, assumption = "a"),
 cbind(scenario_d1_full_results$cohort_deaths_averted_sq_long, assumption = "d1"),
 cbind(scenario_d2_full_results$cohort_deaths_averted_sq_long, assumption = "d2"),
 cbind(scenario_d3_full_results$cohort_deaths_averted_sq_long, assumption = "d3"))
cohort_hbv_deaths_averted_sq$scenario <- factor(cohort_hbv_deaths_averted_sq$scenario, levels =
                                           c("screen_2020_monit_0", "screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1"))

cohort_ly_gained_sq <- rbind(
  cbind(scenario_a_full_results$cohort_ly_gained_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$cohort_ly_gained_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$cohort_ly_gained_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$cohort_ly_gained_sq_long, assumption = "d3"))
cohort_ly_gained_sq$scenario <- factor(cohort_ly_gained_sq$scenario, levels =
                                                  c("screen_2020_monit_0", "screen_2020_monit_10",
                                                    "screen_2020_monit_5", "screen_2020_monit_1"))

# Proportions on the cohort level should be identical for lower screening coverage (since it
# represents the same cohort, just smaller).

# Deaths proportion
ggplot(data = subset(cohort_hbv_deaths_averted_sq,
                    type == "proportion_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# % deaths averted on the cohort level by the basic programme is relatively similar
# no matter the age groups, but looks to be highest in D2 (old age groups only)
# and lowest in D3 (young age groups only). Makes sense that the individual-level benefit
# would be higher for older people than for younger people.

# Deaths number
ggplot(data = subset(cohort_hbv_deaths_averted_sq,
                     type == "number_averted" &
                       scenario %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Looking at the absolute number of deaths averted (which is not specifically a cohort outcome
# as this would be the same in the population), the pattern is reversed: most deaths averted in D2
# followed by D3, presumably cause there are far more young people.
# Number of deaths averted in the cohort should be same as number of deaths averted in the
# total population as the treatment has no significant benefit to those not treated (e.g. reduced transmission).

# LY proportion
ggplot(data = subset(cohort_ly_gained_sq,
                     type == "proportion_averted" &
                       counterfactual %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Proportion of LY saved in the cohort pretty much identical across age groups

# LY number
ggplot(data = subset(cohort_ly_gained_sq,
                     type == "number_averted" &
                       counterfactual %in% c("screen_2020_monit_0"))) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Absolute number of LY saved (again not a cohort-specific outcome) highest in D1 followed by D3
# (those including younger age groups) and far lower in D2.

## Population-level outcomes of monitoring (compared to no monitoring) ----

# Full dataframes of HBV deaths averted and LY saved compared to no monitoring programme
hbv_deaths_averted <- rbind(
  cbind(scenario_a_full_results$deaths_averted_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_long, assumption = "d3"))
hbv_deaths_averted$scenario <- factor(hbv_deaths_averted$scenario, levels =
                                           c("screen_2020_monit_10",
                                             "screen_2020_monit_5", "screen_2020_monit_1"))

ly_gained <- rbind(
  cbind(scenario_a_full_results$ly_gained_long, assumption = "a"),
  cbind(scenario_d1_full_results$ly_gained_long, assumption = "d1"),
  cbind(scenario_d2_full_results$ly_gained_long, assumption = "d2"),
  cbind(scenario_d3_full_results$ly_gained_long, assumption = "d3"))
colnames(ly_gained)[colnames(ly_gained) %in% c("counterfactual", "scenario")] <- c("scenario", "counterfactual")
ly_gained$scenario <- factor(ly_gained$scenario, levels =
                                  c("screen_2020_monit_10",
                                    "screen_2020_monit_5", "screen_2020_monit_1"))

# Is monitoring more beneficial depending on which age group was included?

# Deaths proportion
ggplot(data = subset(hbv_deaths_averted,
                     type == "proportion_averted")) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(scenario ~ by_year) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# For all frequencies and at all timepoints, monitoring appears most beneficial for D1 and D3 (involving
# younger age groups). It is much less beneficial if only older age groups are included (D2).
# In the long-term, the effect of monitoring is fairly similar for A and D3.

# Comparing monitoring frequencies
# Deaths proportion
ggplot(data = subset(hbv_deaths_averted,
                     type == "proportion_averted")) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
  facet_wrap(assumption ~ by_year, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# More monitoring always more beneficial to some degree, but this slightly levels off over time
# (e.g. 1 yearly looks far better than 5-yearly by 2030, but by 2050 this doesn't seem so different).
# Again trade-off between screening frequency and included age groups, e.g.:
# 5-yearly monitoring in D1 appears to have similar effect to yearly monitoring in D3.

## Cohort-level outcomes of monitoring (compared to no monitoring) ----

# Full dataframes of HBV deaths averted and LY saved compared to programme without monitoring
cohort_hbv_deaths_averted <- rbind(
  cbind(scenario_a_full_results$cohort_deaths_averted_long, assumption = "a"),
  cbind(scenario_d1_full_results$cohort_deaths_averted_long, assumption = "d1"),
  cbind(scenario_d2_full_results$cohort_deaths_averted_long, assumption = "d2"),
  cbind(scenario_d3_full_results$cohort_deaths_averted_long, assumption = "d3"))
cohort_hbv_deaths_averted$scenario <- factor(cohort_hbv_deaths_averted$scenario, levels =
                                               c("screen_2020_monit_10",
                                                 "screen_2020_monit_5", "screen_2020_monit_1"))

# Deaths proportion
ggplot(data = subset(cohort_hbv_deaths_averted,
                    type == "proportion_averted")) +
  geom_boxplot(aes(x=scenario, y = value, fill = assumption)) +
#  facet_wrap(~assumption) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# Impact of monitoring in D2 less than the others, but otherwise the difference doesn't
# seem to be huge between the scenarios. D3 highest followed by D1, then A.
# Monitoring gets slightly more beneficial for the cohort if screening is started at a younger age.
# Monitoring is fairly beneficial on the cohort level overall though (even 10 yearly can avert another 25% of deaths)

# ISSUE WITH A NEGATIVE NUMBER IN D3  - likely due to runtime not being long enough?
