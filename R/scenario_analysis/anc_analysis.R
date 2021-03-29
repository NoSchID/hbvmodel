# ANC analysis

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
library(BCEA)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Load functions in monitoring_frequency_by_age_incremental_analysis

# Load data ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/anc_analysis/"

# Status quo
out2 <- readRDS(paste0(out_path, "out2_status_quo_301120.rds"))
out2 <- out2[[1]]

# Population-based strategies with monitoring 5 years in <45 year olds
monit_out7 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_050321.rds"))
monit_out7 <- monit_out7[[1]]

out8b_2030_sim7 <- readRDS(paste0(out_path, "a1_it_out8b_monit_sim7_screen_10b_2030_290321.rds"))
out8b_2030_sim7 <- out8b_2030_sim7[[1]]

# Combination of population-based and ANC strategies with monitoring 5 years in <45 year olds
pop_2020_anc_2030_sim7 <- readRDS(paste0(out_path, "pop_2020_anc_2030_no_rescreen_monit_sim7_210121.rds"))
pop_2020_anc_2030_sim7 <- pop_2020_anc_2030_sim7[[1]]
pop_2020_anc_2040_sim7 <- readRDS(paste0(out_path, "pop_2020_anc_2040_no_rescreen_monit_sim7_220121.rds"))
pop_2020_anc_2040_sim7 <- pop_2020_anc_2040_sim7[[1]]

# For strategies involving population-based screening, need no monitoring
# simulations to extract monitoring assessments:
out3_it <- readRDS(paste0(out_path, "a1_it_out3_screen_2020_monit_0_180121.rds"))
out3_it <- out3_it[[1]]
out8b_2030_monit_0 <- readRDS(paste0(out_path, "a1_it_out8b_monit_0_screen_10b_2030_080121.rds"))
out8b_2030_monit_0 <- out8b_2030_monit_0[[1]]
pop_2020_anc_2030_monit_0 <- readRDS(paste0(out_path, "pop_2020_anc_2030_no_rescreen_monit_0_210121.rds"))
pop_2020_anc_2030_monit_0 <- pop_2020_anc_2030_monit_0[[1]]
pop_2020_anc_2040_monit_0 <- readRDS(paste0(out_path, "pop_2020_anc_2040_no_rescreen_monit_0_210121.rds"))
pop_2020_anc_2040_monit_0 <- pop_2020_anc_2040_monit_0[[1]]


# ANC strategies
# No monitoring
anc_2020_monit_0 <- readRDS(paste0(out_path, "anc1_out3_screen_2020_monit_0_150121.rds"))
anc_2020_monit_0 <- anc_2020_monit_0[[1]]
anc_2030_monit_0 <- readRDS(paste0(out_path, "anc1_2030_no_rescreen_monit_0_150121.rds"))
anc_2030_monit_0 <- anc_2030_monit_0[[1]]
anc_2040_monit_0 <- readRDS(paste0(out_path, "anc1_2040_no_rescreen_monit_0_150121.rds"))
anc_2040_monit_0 <- anc_2040_monit_0[[1]]

# Monitor 5-yearly in <45 year olds
anc_2020_monit_sim7 <- readRDS(paste0(out_path, "anc1_out3_screen_2020_monit_sim7_290321.rds"))
anc_2020_monit_sim7 <- anc_2020_monit_sim7[[1]]
anc_2030_monit_sim7 <- readRDS(paste0(out_path, "anc1_2030_no_rescreen_monit_sim7_290321.rds"))
anc_2030_monit_sim7 <- anc_2030_monit_sim7[[1]]
anc_2040_monit_sim7 <- readRDS(paste0(out_path, "anc1_2040_no_rescreen_monit_sim7_290321.rds"))
anc_2040_monit_sim7 <- anc_2040_monit_sim7[[1]]
anc_2100_monit_sim7 <- readRDS(paste0(out_path, "anc1_2100_no_rescreen_monit_sim7_290321.rds"))
anc_2100_monit_sim7 <- anc_2100_monit_sim7[[1]]



# 1) Cost-effectiveness of ANC compared to population-based screening (no PMTCT) ----
# Create data frame with discounted costs and outcomes ----
annual_discounting_rate <- 0.03
time_horizon_year <- 2100

assemble_discounted_interactions_for_screening_strategies <- function(scenario_object,
                                                                      discount_rate = annual_discounting_rate,
                                                                      assessment_object,
                                                                      time_horizon=time_horizon_year) {
  # Compares monitoring to out3 and everything else to status quo
  # at discount rate of 3%
  out <- left_join(
    left_join(
      left_join(discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_screened",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "hbsag_tests",
                                              end_year = time_horizon),
                discount_outcome_2020_to_2100(scenario_object=assessment_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_assessed",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "clinical_assessments",
                                              end_year = time_horizon)),
      discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                    object_to_subtract=assessment_object,
                                    outcome="interactions",
                                    interaction_outcome="total_assessed",
                                    yearly_discount_rate=discount_rate,
                                    interaction_colname = "monitoring_assessments",
                                    end_year = time_horizon)),
    discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                  object_to_subtract=NULL,
                                  outcome="interactions",
                                  interaction_outcome="total_treated",
                                  yearly_discount_rate=discount_rate,
                                  interaction_colname = "treatment_initiations",
                                  end_year = time_horizon))

  return(out)
}

# NEED TO CHANGE ASSESSMENT OBJECT FOR ALL STRATEGIES INVOLVING MONITORING
# Comparator objects are: out3, out8b_2030_monit_0, pop_2020_anc_2030_monit_0, pop_2020_anc_2040_monit_0

# Analysis objects are:
# monit_out7, out8b_2030_sim7, pop_2020_anc_2030_sim7, pop_2020_anc_2040_sim7,
# anc_2020_monit_0, anc_2030_monit_0, anc_2040_monit_0,
# anc_2020_monit_sim7, anc_2030_monit_sim7, anc_2040_monit_sim7
# For anc_2100_monit_sim7 need to run without monitoring first

anc_interactions <- rbind(
  # Population-based screening:
  cbind(scenario = "screen_2020_monit_0",   # delete this one later
        assemble_discounted_interactions_for_screening_strategies(out3_it,
                                                                  assessment_object = out3_it)),
  cbind(scenario = "screen_2020_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(monit_out7,
                                                                  assessment_object = out3_it)),
  cbind(scenario = "monit_sim7_screen_10b_2030",
        assemble_discounted_interactions_for_screening_strategies(out8b_2030_sim7,
                                                                  assessment_object = out8b_2030_monit_0)),
  # Combination of pop and ANC testing:
  cbind(scenario = "pop_2020_anc_2030_no_rescreen_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2030_sim7,
                                                                  assessment_object = pop_2020_anc_2030_monit_0)),
  cbind(scenario = "pop_2020_anc_2040_no_rescreen_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2040_sim7,
                                                                  assessment_object = pop_2020_anc_2040_monit_0)),  # ANC strategies
# ANC strategies:
  cbind(scenario = "screen_2020_anc_monit_0",
        assemble_discounted_interactions_for_screening_strategies(anc_2020_monit_0,
                                                                  assessment_object = anc_2020_monit_0)),
  cbind(scenario = "screen_2020_anc_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(anc_2020_monit_sim7,
                                                                  assessment_object = anc_2020_monit_0)),
cbind(scenario = "anc_2030_no_rescreen_monit_0",
      assemble_discounted_interactions_for_screening_strategies(anc_2030_monit_0,
                                                                assessment_object = anc_2030_monit_0)),
cbind(scenario = "anc_2030_no_rescreen_monit_sim7",
      assemble_discounted_interactions_for_screening_strategies(anc_2030_monit_sim7,
                                                                assessment_object = anc_2030_monit_0)),
cbind(scenario = "anc_2040_no_rescreen_monit_0",
      assemble_discounted_interactions_for_screening_strategies(anc_2040_monit_0,
                                                                assessment_object = anc_2040_monit_0)),
cbind(scenario = "anc_2040_no_rescreen_monit_sim7",
      assemble_discounted_interactions_for_screening_strategies(anc_2040_monit_sim7,
                                                                assessment_object = anc_2040_monit_0))

)


# PY on treatment, DALYS and deaths averted:
object_list <- list(out3_it, monit_out7, out8b_2030_sim7,
                    pop_2020_anc_2030_sim7, pop_2020_anc_2040_sim7,
                    anc_2020_monit_0, anc_2020_monit_sim7, anc_2030_monit_0, anc_2030_monit_sim7,
                    anc_2040_monit_0, anc_2040_monit_sim7)

# Extract interactions and person-years on treatment, and HBV-related deaths
# ad DALYs averted, in a for loop
anc_py_on_treatment <- list()
anc_hbv_deaths_averted <- list()
anc_dalys_averted <- list()

for (i in 1:length(object_list)) {
  anc_py_on_treatment[[i]] <-
    data.frame(scenario = object_list[[i]]$cohort_age_at_death$scenario,
               discount_outcome_2020_to_2100(scenario_object=object_list[[i]],
                                             object_to_subtract=NULL,
                                             outcome="py_on_treatment",
                                             yearly_discount_rate=annual_discounting_rate,
                                             end_year = time_horizon_year))
  anc_hbv_deaths_averted[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          discount_outcome_2020_to_2100(scenario_object=out2,
                                        object_to_subtract=object_list[[i]],
                                        outcome="cum_hbv_deaths",
                                        yearly_discount_rate=annual_discounting_rate,
                                        end_year = time_horizon_year))
  anc_dalys_averted[[i]] <-
    cbind(scenario = object_list[[i]]$cohort_age_at_death$scenario,
          discount_outcome_2020_to_2100(scenario_object=out2,
                                        object_to_subtract=object_list[[i]],
                                        outcome="dalys",
                                        yearly_discount_rate=annual_discounting_rate,
                                        end_year = time_horizon_year))
}
anc_py_on_treatment <- do.call("rbind", anc_py_on_treatment)
anc_hbv_deaths_averted <- do.call("rbind", anc_hbv_deaths_averted)
anc_dalys_averted <- do.call("rbind", anc_dalys_averted)

anc_py_on_treatment$sim <- gsub("[^0-9]", "", anc_py_on_treatment$sim)
anc_hbv_deaths_averted$sim <- gsub("[^0-9]", "", anc_hbv_deaths_averted$sim)
colnames(anc_hbv_deaths_averted)[colnames(anc_hbv_deaths_averted) == "cum_hbv_deaths"] <-
  "value"

anc_dalys_averted$sim <- gsub("[^0-9]", "", anc_dalys_averted$sim)
colnames(anc_dalys_averted)[colnames(anc_dalys_averted) == "dalys"] <-
  "value"

# Combine into full dataframe with costs ----
hbsag_test_cost <- 8.3  ## ensure hbsag_test_cost is the same as below
anc_df <- create_incremental_plot_df(interactions_df=anc_interactions,
                                     py_on_treatment_df=anc_py_on_treatment,
                                     deaths_averted_df=anc_hbv_deaths_averted,
                                     ly_saved_df = anc_dalys_averted, # replace LY by DALYs
                                     hbsag_test_cost = hbsag_test_cost,
                                     clinical_assessment_cost = 33, #84.4,
                                     monitoring_assessment_cost = 25.5, #40.1,
                                     treatment_py_cost = 66.5, #60,
                                     #scenario_labels_obj = scenario_labels,
                                     ref_label = "No treatment")
colnames(anc_df)[colnames(anc_df)=="ly_saved"] <- "dalys_averted"

# Correct the cost of HBsAg test in ANC to 1.7 per HBsAg test:

# Strategies with only ANC testing:
anc_only_strategies <- c("screen_2020_anc_monit_0", "screen_2020_anc_monit_sim7",
  "anc_2030_no_rescreen_monit_0", "anc_2030_no_rescreen_monit_sim7",
  "anc_2040_no_rescreen_monit_0", "anc_2030_no_rescreen_monit_sim7")

anc_df[
  anc_df$scenario %in% anc_only_strategies,]$screening_cost <-
  anc_df[
    anc_df$scenario %in% anc_only_strategies,]$hbsag_tests*1.7
anc_df[
  anc_df$scenario %in% anc_only_strategies,]$total_cost <-
  anc_df[
    anc_df$scenario %in% anc_only_strategies,]$screening_cost +
  anc_df[
    anc_df$scenario %in% anc_only_strategies,]$assessment_cost +
  anc_df[
    anc_df$scenario %in% anc_only_strategies,]$monitoring_cost+
  anc_df[
    anc_df$scenario %in% anc_only_strategies,]$treatment_cost

# Strategies combining population and ANC testing

# Need to distinguish between the screens in the population test and those conducted in ANC
# (compare to screening cost in one-time pop test)
# Normal cost from the pop-based test + reduced cost for subsequent ANC screens
pop_anc_combi_strategies <- c("pop_2020_anc_2030_no_rescreen_monit_sim7",
                              "pop_2020_anc_2040_no_rescreen_monit_sim7")

for (i in 1:length(pop_anc_combi_strategies)) {
  anc_df[
    anc_df$scenario == pop_anc_combi_strategies[i],]$screening_cost <-
    anc_df[anc_df$scenario ==
                              "screen_2020_monit_0",]$hbsag_tests*hbsag_test_cost +
    (anc_df[
      anc_df$scenario == pop_anc_combi_strategies[i],]$hbsag_tests-
       anc_df[
         anc_df$scenario == "screen_2020_monit_0",]$hbsag_tests)*1.7
  anc_df[
    anc_df$scenario == pop_anc_combi_strategies[i],]$total_cost <-
    anc_df[
      anc_df$scenario == pop_anc_combi_strategies[i],]$screening_cost +
    anc_df[
      anc_df$scenario == pop_anc_combi_strategies[i],]$assessment_cost +
    anc_df[
      anc_df$scenario == pop_anc_combi_strategies[i],]$monitoring_cost+
    anc_df[
      anc_df$scenario == pop_anc_combi_strategies[i],]$treatment_cost
}

# Delete pop 2020 no monitoring strategy for analysis:
anc_df <- subset(anc_df, scenario != "screen_2020_monit_0")


# Analysis: calculate dominance and ICERs ----
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- anc_df[which(anc_df$sim==
                                             unique(anc_df$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust =1))

# Calculate ICER by simulation on non-dominated strategies
anc_df2 <- subset(anc_df, scenario %in% c("monit_sim7_screen_10b_2030",
                                          "screen_2020_anc_monit_0",
                                          "pop_2020_anc_2040_no_rescreen_monit_sim7",
                                          "anc_2030_no_rescreen_monit_0",
                                          "anc_2040_no_rescreen_monit_0"))

icer_list <- list()

for(i in 1:183) {
  print(i)
  icer_list[[i]] <- anc_df2[which(anc_df2$sim==
                                    unique(anc_df2$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result

# Strategies up until combination of 2020 population-based screen + ANC until 2040
# are cost-effective. Repeat population-based screen in 2030 is not dominated
# but under current CE threshold not cost-effective compared to ANC strategies
