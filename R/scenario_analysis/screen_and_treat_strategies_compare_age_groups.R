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

# F1 = Low assessment uptake, ages 15-65
scenario_f1_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_f1_basic_results.rds"))

# G1 = Low treatment uptake/retention, ages 15-65
scenario_g1_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_g1_basic_results.rds"))

# BX = vaccine introduction in 2004, ages 30-70
scenario_bx_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_bx_basic_results.rds"))

# BX1 = vaccine introduction in 2004, ages 15-65
scenario_bx1_full_results <- readRDS(here("output", "screen_and_treat_results",
                                         "scenario_bx1_basic_results.rds"))

# BX2 = vaccine introduction in 2004, ages 45-70
scenario_bx2_full_results <- readRDS(here("output", "screen_and_treat_results",
                                          "scenario_bx2_basic_results.rds"))

# BX3 = vaccine introduction in 2004, ages 15-45
scenario_bx3_full_results <- readRDS(here("output", "screen_and_treat_results",
                                          "scenario_bx3_basic_results.rds"))

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

## Plot of population- vs cohort effect of monitoring in the 15-65 year group ----

# HBV deaths averted (by 2100 in the population)
pop_cohort_comp_deaths <- rbind(
      scenario_d1_full_results$cohort_deaths_averted_sq_long,
      scenario_d1_full_results$deaths_averted_sq_long[
        scenario_d1_full_results$deaths_averted_sq_long$by_year == 2100,-c(1,2)])

# LY saved
pop_cohort_comp_ly <- rbind(
  scenario_d1_full_results$cohort_ly_gained_sq_long,
  scenario_d1_full_results$ly_gained_sq_long[
    scenario_d1_full_results$ly_gained_sq_long$by_year == 2100,-c(1,2)])
# Why is number of LY saved not the same in the population as in the cohort?

# HBV deaths averted plot

# Facet label names
cohort_pop_label <- c("Cohort-level", "Population-level")
names(cohort_pop_label) <- c("status_quo_cohort", "status_quo")

ggplot(subset(pop_cohort_comp_deaths, type == "proportion_averted" &
                scenario %in% c("screen_2020_monit_0", "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = counterfactual)) +
  facet_wrap(~counterfactual, labeller = labeller(counterfactual = cohort_pop_label)) +
  xlab("Treatment programme strategy") + ylab("Proportion of HBV-related deaths averted\ncompared to no treatment") +
  scale_x_discrete(labels=c("screen_2020_monit_0" = "No monitoring",
                            "screen_2020_monit_1" = "Yearly monitoring")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(0,1)

# What proportion of deaths averted in the yearly monitoring programme are from the
# initial screen+treat?
# Total number averted by basic programme compared to no treatment divided by
# total number averted monitoring programme compared to no treatment
quantile(subset(pop_cohort_comp_deaths, type == "number_averted" &
         scenario == "screen_2020_monit_0" & counterfactual == "status_quo_cohort")$value/
  subset(pop_cohort_comp_deaths, type == "number_averted" &
           scenario == "screen_2020_monit_1" & counterfactual == "status_quo_cohort")$value,
  prob = c(0.5,0.025,0.975))
# In the cohort, 61% (49-76%) of averted HBV deaths are averted through the initial screening and
# treatment.
quantile(subset(pop_cohort_comp_deaths, type == "number_averted" &
                  scenario == "screen_2020_monit_0" & counterfactual == "status_quo")$value/
           subset(pop_cohort_comp_deaths, type == "number_averted" &
                    scenario == "screen_2020_monit_1" & counterfactual == "status_quo")$value,
         prob = c(0.5,0.025,0.975))
# Pretty much same values in the population.
# What about LY?
quantile(subset(pop_cohort_comp_ly, type == "number_averted" &
                  counterfactual == "screen_2020_monit_0" & scenario == "status_quo_cohort")$value/
           subset(pop_cohort_comp_ly, type == "number_averted" &
                    counterfactual == "screen_2020_monit_1" & scenario == "status_quo_cohort")$value,
         prob = c(0.5,0.025,0.975))
# In cohort, 67% (56-81%) but slightly higher on average in population.

# How many extra deaths are averted by monitoring? Same in cohort and population.
quantile(subset(pop_cohort_comp_deaths, type == "number_averted" &
         scenario == "screen_2020_monit_1" & counterfactual == "status_quo_cohort")$value-
  subset(pop_cohort_comp_deaths, type == "number_averted" &
           scenario == "screen_2020_monit_0" & counterfactual == "status_quo_cohort")$value,
  prob = c(0.5,0.025,0.975))
# This number is the same as previously calculated when comparing monitoring to no monitoring:
median(subset(pop_cohort_comp_deaths_monit, type == "number_averted" &
         scenario == "screen_2020_monit_1" & x == "cohort")$value)
# What proportion of deaths does this represent compared to the deaths that would occur
# in the programme without monitoring? Need to look at original number of deaths, NOT
# deaths averted by the basic treatment programme.

# Yearly monitoring would avert an additional 2348 (734-5298) HBV-related deaths. This represents
# 50% (23-72%) of HBV-related deaths that would occur in the screened and treated cohort in the absence
# of monitoring, and 19% (8-28%) of HBV-related deaths occurring in the Gambian population without
# monitoring.
quantile(subset(pop_cohort_comp_deaths_monit, type == "proportion_averted" &
         scenario == "screen_2020_monit_1" & x == "pop")$value, prob = c(0.5,0.025,0.975))

# Effect of monitoring on HBV deaths averted (by 2100 in the population)
pop_cohort_comp_deaths_monit <- rbind(
  cbind(scenario_d1_full_results$cohort_deaths_averted_long, x="cohort"),
  cbind(scenario_d1_full_results$deaths_averted_long[
    scenario_d1_full_results$deaths_averted_long$by_year == 2100,-c(1,2)], x= "pop"))

# HBV deaths averted
ggplot(subset(pop_cohort_comp_deaths_monit, type == "proportion_averted" &
                scenario %in% c("screen_2020_monit_5", "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = x)) +
  facet_wrap(~x) +
  theme_bw() +
  ylim(0,1)
# Shows that monitoring averts more larger proportion of deaths compared to no monitoring in the
# cohort than in the population

ggplot(subset(scenario_d1_full_results$deaths_averted_long, type == "proportion_averted" &
                scenario %in% c("screen_2020_monit_5", "screen_2020_monit_1"))) +
  geom_boxplot(aes(x=scenario, y = value, fill = by_year)) +
  facet_wrap(~by_year) +
  theme_bw()

## One-way sensitivity analysis of the cascade of care (example) ----

# Currently for 1 age group (15-65) and basic programme only:
# compare D1 (optimal - 90% screening, 80% assessment, 100% treatment),
# E1 (10% screening coverage), F1 (40% assessment) and G1 (40% treatment)

# HBV deaths averted compared to infant vaccine only
cascade_deaths_averted_sq <- rbind(
  cbind(scenario_d1_full_results$deaths_averted_sq_long, assumption = "d1"),
  cbind(scenario_e1_full_results$deaths_averted_sq_long, assumption = "e1"),
  cbind(scenario_f1_full_results$deaths_averted_sq_long, assumption = "f1"),
  cbind(scenario_g1_full_results$deaths_averted_sq_long, assumption = "g1")
)
cascade_deaths_averted_sq$screening_coverage <- 0.9
cascade_deaths_averted_sq$screening_coverage[cascade_deaths_averted_sq$assumption == "e1"] <- 0.1
cascade_deaths_averted_sq$assessment_uptake <- 0.8
cascade_deaths_averted_sq$assessment_uptake[cascade_deaths_averted_sq$assumption == "f1"] <- 0.4
cascade_deaths_averted_sq$treatment_uptake <- 1
cascade_deaths_averted_sq$treatment_uptake[cascade_deaths_averted_sq$assumption == "g1"] <- 0.4

# Visualise what result could look like (but would plot median and 95% CrI with ribbon)

p1 <- ggplot(subset(cascade_deaths_averted_sq,
              assumption %in% c("d1", "e1") & type == "proportion_averted" & by_year == 2050),
             aes(x = screening_coverage, y = value, group = assumption)) +
  geom_boxplot(width = 0.1) +
  xlim(0,1.1) + ylim(0,1) +
  ylab("Proportion of HBV deaths averted") +
  stat_summary(fun=median, geom="smooth", aes(group=1), lwd=1) +
  theme_classic()

p2 <- ggplot(subset(cascade_deaths_averted_sq,
              assumption %in% c("d1", "f1") & type == "proportion_averted" & by_year == 2050),
             aes(x = assessment_uptake, y = value, group = assumption)) +
  geom_boxplot(width = 0.1) +
  xlim(0,1.1) + ylim(0,1) +
  ylab("") +
  stat_summary(fun=median, geom="smooth", aes(group=1), lwd=1) +
  theme_classic()

p3 <- ggplot(subset(cascade_deaths_averted_sq,
                    assumption %in% c("d1", "g1") & type == "proportion_averted" & by_year == 2050),
             aes(x = treatment_uptake, y = value, group = assumption)) +
  geom_boxplot(width = 0.1) +
  xlim(0,1.1) + ylim(0,1) +
  ylab("") +
  stat_summary(fun=median, geom="smooth", aes(group=1), lwd=1) +
  theme_classic()

library(gridExtra)
grid.arrange(p1,p2,p3, ncol = 3)
# Slope for assessment uptake appears steeper than for screening coverage and treatment uptake
# (this becomes more visible over time)
# Linear trend may not be true.

# HBV deaths averted per interaction
cascade_deaths_averted_per_interaction_sq <- rbind(
  cbind(scenario_d1_full_results$deaths_averted_per_interaction_sq_long, assumption = "d1"),
  cbind(scenario_e1_full_results$deaths_averted_per_interaction_sq_long, assumption = "e1"),
  cbind(scenario_f1_full_results$deaths_averted_per_interaction_sq_long, assumption = "f1"),
  cbind(scenario_g1_full_results$deaths_averted_per_interaction_sq_long, assumption = "g1")
)
cascade_deaths_averted_per_interaction_sq$screening_coverage <- 0.9
cascade_deaths_averted_per_interaction_sq$screening_coverage[cascade_deaths_averted_per_interaction_sq$assumption == "e1"] <- 0.1
cascade_deaths_averted_per_interaction_sq$assessment_uptake <- 0.8
cascade_deaths_averted_per_interaction_sq$assessment_uptake[cascade_deaths_averted_per_interaction_sq$assumption == "f1"] <- 0.4
cascade_deaths_averted_per_interaction_sq$treatment_uptake <- 1
cascade_deaths_averted_per_interaction_sq$treatment_uptake[cascade_deaths_averted_per_interaction_sq$assumption == "g1"] <- 0.4


p1a <- ggplot(subset(cascade_deaths_averted_sq,
                    assumption %in% c("d1", "e1") & type == "proportion_averted" & by_year == 2030),
             aes(x = screening_coverage, y = 1/value, group = assumption)) +
  geom_boxplot(width = 0.1) +
  xlim(0,1.1) + ylim(0,35) +
  ylab("Interactions per death averted") +
  stat_summary(fun=median, geom="smooth", aes(group=1), lwd=1) +
  theme_classic()

p2a <- ggplot(subset(cascade_deaths_averted_sq,
                    assumption %in% c("d1", "f1") & type == "proportion_averted" & by_year == 2030),
             aes(x = assessment_uptake, y = 1/value, group = assumption)) +
  geom_boxplot(width = 0.1) +
  xlim(0,1.1) + ylim(0,35) +
  ylab("") +
  stat_summary(fun=median, geom="smooth", aes(group=1), lwd=1) +
  theme_classic()

p3a <- ggplot(subset(cascade_deaths_averted_sq,
                    assumption %in% c("d1", "g1") & type == "proportion_averted" & by_year == 2030),
             aes(x = treatment_uptake, y = 1/value, group = assumption)) +
  geom_boxplot(width = 0.1) +
  xlim(0,1.1) + ylim(0,35) +
  ylab("") +
  stat_summary(fun=median, geom="smooth", aes(group=1), lwd=1) +
  theme_classic()

grid.arrange(p1a,p2a,p3a, ncol = 3)
# Higher uptake along the casecade lead to reductions in the number of healthcare
# interactions required to avert 1 HBV death. Reduction appears largest (steepest slope)
# for screening, but HBsAg tests are the easiest and cheapest intervention.
# Pattern is similar although less pronounced by 2050.
# Linear trend may not be true.

## Compare population-level impact and age groups if infant vaccine was introduced in 2004----
# Issue with cohort projections: not everyone has died in BX1 and BX3
# Overall conclusion from the population-level analysis:
# Later vaccine introduction favours inclusion/focus on younger age groups even more!

# Full dataframes of HBV deaths averted and LY saved compared to infant vaccine only
hbv_deaths_averted_sq_vaccintro <- rbind(
  cbind(scenario_a_full_results$deaths_averted_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$deaths_averted_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$deaths_averted_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$deaths_averted_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$deaths_averted_sq_long, assumption = "bx3"))
hbv_deaths_averted_sq_vaccintro$vaccintro <- "1990"
hbv_deaths_averted_sq_vaccintro$vaccintro[hbv_deaths_averted_sq_vaccintro$assumption
                                          %in% c("bx", "bx1", "bx2", "bx3")] <- "2004"

ly_gained_sq_vaccintro <- rbind(
  cbind(scenario_a_full_results$ly_gained_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$ly_gained_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$ly_gained_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$ly_gained_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$ly_gained_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$ly_gained_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$ly_gained_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$ly_gained_sq_long, assumption = "bx3"))
colnames(ly_gained_sq_vaccintro)[colnames(ly_gained_sq_vaccintro) %in% c("counterfactual", "scenario")] <- c("scenario", "counterfactual")
ly_gained_sq_vaccintro$vaccintro <- "1990"
ly_gained_sq_vaccintro$vaccintro[ly_gained_sq_vaccintro$assumption
                                          %in% c("bx", "bx1", "bx2", "bx3")] <- "2004"

# Deaths
ggplot(data = subset(hbv_deaths_averted_sq_vaccintro,
                     type == "proportion_averted" & by_year == 2030)) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(~ vaccintro, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# By 2030: Effect for BX1/D1 (15-65) is very similar, and also BX3/D3
# BX and BX2 are worse than A and D2. For vaccintro in 1990, the order
# from best to worst was: D1, A, D3, D2. For 2004 vaccintro, it is BX1, BX3, A, BX2.
# => Later vaccine introduction favours screening of younger age groups even more.
# Number averted is higher for 2004 than for 1990 (higher absolute impact of treatment
# if prevention is lagging behind).

ggplot(data = subset(hbv_deaths_averted_sq_vaccintro,
                     type == "proportion_averted" & by_year == 2050)) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(~ vaccintro, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# By 2050: Same except that the overall treatment impact % for the best age group (15-65)
# is lower (median 30%) than if vaccine was introduced in 1990 (median 40%).
# However, the absolute impact (number of deaths averted) is still higher, there
# would just be more HBV deaths overall.

ggplot(data = subset(hbv_deaths_averted_sq_vaccintro,
                     type == "proportion_averted" & by_year == 2100)) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(~ vaccintro, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# 2100 same pattern as 2050. BX3 nearly as good as BX1 (BX and BX2 much lower). But here for default vaccintro
# D3 also equivalent to A.

# LY saved
ggplot(data = subset(ly_gained_sq_vaccintro,
                     type == "proportion_averted" & by_year == 2030)) +
  geom_boxplot(aes(x=assumption, y = value, fill = assumption)) +
  facet_wrap(~ vaccintro, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# For all timepoints:
# Both number and proportion LY saved higher for 2004 than for 1990. Same age group order as for deaths.

# Overall conclusion of impact: If infant vaccine had been introduced later, treatment would
# avert a smaller proportion of HBV deaths but would save a larger proportion of life years.
# 15-65 years is still the best age group to screen for overall impact, but later vaccine introduction
# favours focus on the younger age groups. 30-70 year olds looking far worse and screening
# 45-70 year olds only has negligible impact on the population level.
# The 30-70 age group might be better for the early vaccine intro cause here it corresponds to
# those that were not vaccinated. For the 2004 intro, 16+ year olds would not be vaccinated.

# What about compared to resource use?

# Full dataframe of HBV deaths averted per healthcare interaction
hbv_deaths_averted_per_interaction_sq_vaccintro <- rbind(
  cbind(scenario_a_full_results$deaths_averted_per_interaction_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_per_interaction_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_per_interaction_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_per_interaction_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$deaths_averted_per_interaction_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$deaths_averted_per_interaction_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$deaths_averted_per_interaction_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$deaths_averted_per_interaction_sq_long, assumption = "bx3"),
  cbind(scenario_a_full_results$deaths_averted_per_test_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_per_test_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_per_test_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_per_test_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$deaths_averted_per_test_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$deaths_averted_per_test_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$deaths_averted_per_test_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$deaths_averted_per_test_sq_long, assumption = "bx3"),
  cbind(scenario_a_full_results$deaths_averted_per_assessment_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_per_assessment_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_per_assessment_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_per_assessment_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$deaths_averted_per_assessment_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$deaths_averted_per_assessment_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$deaths_averted_per_assessment_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$deaths_averted_per_assessment_sq_long, assumption = "bx3"),
  cbind(scenario_a_full_results$deaths_averted_per_treatment_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$deaths_averted_per_treatment_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$deaths_averted_per_treatment_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$deaths_averted_per_treatment_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$deaths_averted_per_treatment_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$deaths_averted_per_treatment_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$deaths_averted_per_treatment_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$deaths_averted_per_treatment_sq_long, assumption = "bx3"))
hbv_deaths_averted_per_interaction_sq_vaccintro$vaccintro <- "1990"
hbv_deaths_averted_per_interaction_sq_vaccintro$vaccintro[hbv_deaths_averted_per_interaction_sq_vaccintro$assumption
                                          %in% c("bx", "bx1", "bx2", "bx3")] <- "2004"

# LY gained per healthcare interaction
ly_gained_per_interaction_sq_vaccintro <- rbind(
  cbind(scenario_a_full_results$ly_gained_per_interaction_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$ly_gained_per_interaction_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$ly_gained_per_interaction_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$ly_gained_per_interaction_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$ly_gained_per_interaction_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$ly_gained_per_interaction_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$ly_gained_per_interaction_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$ly_gained_per_interaction_sq_long, assumption = "bx3"),
  cbind(scenario_a_full_results$ly_gained_per_assessment_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$ly_gained_per_assessment_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$ly_gained_per_assessment_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$ly_gained_per_assessment_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$ly_gained_per_assessment_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$ly_gained_per_assessment_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$ly_gained_per_assessment_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$ly_gained_per_assessment_sq_long, assumption = "bx3"),
  cbind(scenario_a_full_results$ly_gained_per_treatment_sq_long, assumption = "a"),
  cbind(scenario_d1_full_results$ly_gained_per_treatment_sq_long, assumption = "d1"),
  cbind(scenario_d2_full_results$ly_gained_per_treatment_sq_long, assumption = "d2"),
  cbind(scenario_d3_full_results$ly_gained_per_treatment_sq_long, assumption = "d3"),
  cbind(scenario_bx_full_results$ly_gained_per_treatment_sq_long, assumption = "bx"),
  cbind(scenario_bx1_full_results$ly_gained_per_treatment_sq_long, assumption = "bx1"),
  cbind(scenario_bx2_full_results$ly_gained_per_treatment_sq_long, assumption = "bx2"),
  cbind(scenario_bx3_full_results$ly_gained_per_treatment_sq_long, assumption = "bx3"))
ly_gained_per_interaction_sq_vaccintro$vaccintro <- "1990"
ly_gained_per_interaction_sq_vaccintro$vaccintro[ly_gained_per_interaction_sq_vaccintro$assumption
                                                          %in% c("bx", "bx1", "bx2", "bx3")] <- "2004"

# Number needed to interact, assess and treat
nnt_df <- group_by(hbv_deaths_averted_per_interaction_sq_vaccintro,
                   vaccintro, assumption, scenario, by_year, interaction_type) %>%
  summarise(median = median(1/value),
            cri_lower = quantile(1/value, prob = 0.025),
            cri_upper = quantile(1/value, prob = 0.975))

ggplot(data = subset(nnt_df,
                     assumption %in% c("a", "d1", "d2", "d3") & by_year == 2050 & vaccintro == "1990" &
                       scenario == "screen_2020_monit_0")) +
  geom_col(aes(x= reorder(interaction_type, -median), y = median, group = assumption, fill = assumption),
           position = "dodge") +
  geom_errorbar(aes(x= reorder(interaction_type, -median), ymin = cri_lower, ymax = cri_upper,
                    group = assumption), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE) +
  facet_wrap( ~ interaction_type, scales= "free") +
  theme_classic()

ggplot(data = subset(nnt_df,
                     assumption %in% c("a", "d1", "d2", "d3") & by_year == 2050 & vaccintro == "1990" &
                       scenario == "screen_2020_monit_0")) +
  geom_col(aes(x= reorder(interaction_type, -median), y = median, group = assumption, fill = assumption),
           position = "dodge") +
  geom_errorbar(aes(x= reorder(interaction_type, -median), ymin = cri_lower, ymax = cri_upper,
                    group = assumption), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE) +
  facet_wrap( ~ assumption) +
  theme_classic()

ggplot(data = subset(nnt_df,
                     assumption =="d1" & by_year == 2050 & vaccintro == "1990" &
                       scenario == "screen_2020_monit_0")) +
  geom_col(aes(x= reorder(interaction_type, -median), y = median, group = interaction_type, fill = interaction_type),
           position = "stack") +
  geom_errorbar(aes(x= reorder(interaction_type, -median), ymin = cri_lower, ymax = cri_upper,
                    group = assumption), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE) +
  theme_classic()

# Total interactions
ggplot(data = subset(hbv_deaths_averted_per_interaction_sq_vaccintro,
                     interaction_type == "total_interactions" & by_year == 2030)) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~vaccintro, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# By 2030: For the 1990 vaccintro, there is a bit of variation with higher interactions
# required per death averted in D1 and D3 than A and D2
# However in the 2004 vaccintro, they all look extremely similar on average with only BX
# slightly lower than the others.
# The number of interactions required to avert 1 death is substantially lower for the 2004 than the
# 1990 scenario.

ggplot(data = subset(hbv_deaths_averted_per_interaction_sq_vaccintro,
                     interaction_type == "total_interactions" & by_year == 2050)) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~vaccintro, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# By 2050 and 2100, the pattern for 1990 has not changed but for 2004 vaccintro only BX2
# now stands out as requiring more interactions per death averted.

# Assessments
ggplot(data = subset(hbv_deaths_averted_per_interaction_sq_vaccintro,
                     interaction_type == "total_assessed" & by_year == 2030)) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~vaccintro, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# At all timepoints:
# Looking at assessments for 1990, they are very similar but there is now some variation
# for 2004, which is higher for BX1 and BX3 (If vaccine is introduced later,
# screening younger people requires more liver disease assessments to avert 1 death than screening older age groups.

# Treatments
ggplot(data = subset(hbv_deaths_averted_per_interaction_sq_vaccintro,
                     interaction_type == "total_treated" & by_year == 2030)) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~vaccintro, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# By 2030, Treatment pattern by age is the same for either vaccine intro: slightly higher if including
# younger age groups.
# By 2100, relatively similar overall across age groups except for D2 and BX2 being higher.

# LY gained

# Total interactions
ggplot(data = subset(ly_gained_per_interaction_sq_vaccintro,
                     interaction_type == "total_interactions" & by_year == 2030)) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~vaccintro, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# By 2030, very similar to deaths, except that BX2 slightly higher than the other BX.
# By 2050, BX2 still the outlier for the 2004 vaccintro and BX1 and BX3 are the lowest.
# For 1990 vaccintro, D1-D3 are now very similar and only A is slightly lower.
# By 2100, D2 also requires far more interactions than the others.
# Probably shows that focus on older age groups has fewer long-term benefits and gets even more
# "expensive" compared to other screening strategies if vaccine was introduced later.
# Nevertheless the later vaccine intro means that fewer interactions are required overall
# to save 1 LY.

# Assessments
ggplot(data = subset(ly_gained_per_interaction_sq_vaccintro,
                     interaction_type == "total_assessed" & by_year == 2030)) +
  geom_boxplot(aes(x=assumption, y = 1/value, fill = assumption)) +
  facet_wrap(~vaccintro, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1))
# By 2030, all the same in 1990, minimally higher for BX1 and BX3 in 2004.
# By 2050 and 2100, similar trend to total interactions with BX2 (and D2) standing out as
# particularly high.




