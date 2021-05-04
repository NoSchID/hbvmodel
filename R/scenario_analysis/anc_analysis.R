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
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/repeat_screening_anc_analysis/"

# Status quo
out2 <- readRDS(paste0(out_path, "out2_status_quo_301120.rds"))
out2 <- out2[[1]]

# Population-based strategies with monitoring 5 years in <45 year olds
monit_out7 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_050321.rds"))
monit_out7 <- monit_out7[[1]]

# Repeat screening strategies
out8b_2030_sim7 <- readRDS(paste0(out_path, "a1_it_out8b_monit_sim7_screen_10b_2030_290321.rds"))
out8b_2030_sim7 <- out8b_2030_sim7[[1]]
out8a_2030_sim7 <- readRDS(paste0(out_path, "a1_it_out8a_monit_sim7_screen_10a_2030_130421.rds"))
out8a_2030_sim7 <- out8a_2030_sim7[[1]]
out8b_2040_sim7 <- readRDS(paste0(out_path, "a1_it_out8b_monit_sim7_screen_10b_2040_130421.rds"))
out8b_2040_sim7 <- out8b_2040_sim7[[1]]
out8a_2040_sim7 <- readRDS(paste0(out_path, "a1_it_out8a_monit_sim7_screen_10a_2040_130421.rds"))
out8a_2040_sim7 <- out8a_2040_sim7[[1]]
out8b_2050_sim7 <- readRDS(paste0(out_path, "a1_it_out8b_monit_sim7_screen_10b_2050_130421.rds"))
out8b_2050_sim7 <- out8b_2050_sim7[[1]]
out8a_2050_sim7 <- readRDS(paste0(out_path, "a1_it_out8a_monit_sim7_screen_10a_2050_130421.rds"))
out8a_2050_sim7 <- out8a_2050_sim7[[1]]

# Combination of population-based and ANC strategies with monitoring 5 years in <45 year olds
pop_2020_anc_2030_sim7 <- readRDS(paste0(out_path, "pop_2020_anc_2030_no_rescreen_monit_sim7_210121.rds"))
pop_2020_anc_2030_sim7 <- pop_2020_anc_2030_sim7[[1]]
pop_2020_anc_2040_sim7 <- readRDS(paste0(out_path, "pop_2020_anc_2040_no_rescreen_monit_sim7_220121.rds"))
pop_2020_anc_2040_sim7 <- pop_2020_anc_2040_sim7[[1]]
pop_2020_anc_2050_sim7 <- readRDS(paste0(out_path, "pop_2020_anc_2050_no_rescreen_monit_sim7_300421.rds"))
pop_2020_anc_2050_sim7 <- pop_2020_anc_2050_sim7[[1]]
# With rescreens:
pop_2020_anc_2030_sim7a <- readRDS(paste0(out_path, "pop_2020_anc_2030_with_rescreen_monit_sim7_300421.rds"))
pop_2020_anc_2030_sim7a <- pop_2020_anc_2030_sim7a[[1]]
pop_2020_anc_2040_sim7a <- readRDS(paste0(out_path, "pop_2020_anc_2040_with_rescreen_monit_sim7_300421.rds"))
pop_2020_anc_2040_sim7a <- pop_2020_anc_2040_sim7a[[1]]
pop_2020_anc_2050_sim7a <- readRDS(paste0(out_path, "pop_2020_anc_2050_with_rescreen_monit_sim7_300421.rds"))
pop_2020_anc_2050_sim7a <- pop_2020_anc_2050_sim7a[[1]]


# THESE NEED NO MONITORING COMPARATOR!!

# For strategies involving population-based screening, need no monitoring
# simulations to extract monitoring assessments:
out3_it <- readRDS(paste0(out_path, "a1_it_out3_screen_2020_monit_0_180121.rds"))
out3_it <- out3_it[[1]]
out8b_2030_monit_0 <- readRDS(paste0(out_path, "a1_it_out8b_monit_0_screen_10b_2030_220421.rds"))
out8b_2030_monit_0 <- out8b_2030_monit_0[[1]]
out8a_2030_monit_0 <- readRDS(paste0(out_path, "a1_it_out8a_monit_0_screen_10a_2030_220421.rds"))
out8a_2030_monit_0 <- out8a_2030_monit_0[[1]]
out8b_2040_monit_0 <- readRDS(paste0(out_path, "a1_it_out8b_monit_0_screen_10b_2040_220421.rds"))
out8b_2040_monit_0 <- out8b_2040_monit_0[[1]]
out8a_2040_monit_0 <- readRDS(paste0(out_path, "a1_it_out8a_monit_0_screen_10a_2040_220421.rds"))
out8a_2040_monit_0 <- out8a_2040_monit_0[[1]]
out8b_2050_monit_0 <- readRDS(paste0(out_path, "a1_it_out8b_monit_0_screen_10b_2050_220421.rds"))
out8b_2050_monit_0 <- out8b_2050_monit_0[[1]]
out8a_2050_monit_0 <- readRDS(paste0(out_path, "a1_it_out8a_monit_0_screen_10a_2050_220421.rds"))
out8a_2050_monit_0 <- out8a_2050_monit_0[[1]]
pop_2020_anc_2030_monit_0 <- readRDS(paste0(out_path, "pop_2020_anc_2030_no_rescreen_monit_0_210121.rds"))
pop_2020_anc_2030_monit_0 <- pop_2020_anc_2030_monit_0[[1]]
pop_2020_anc_2040_monit_0 <- readRDS(paste0(out_path, "pop_2020_anc_2040_no_rescreen_monit_0_210121.rds"))
pop_2020_anc_2040_monit_0 <- pop_2020_anc_2040_monit_0[[1]]
pop_2020_anc_2050_monit_0 <- readRDS(paste0(out_path, "pop_2020_anc_2050_no_rescreen_monit_0_040521.rds"))
pop_2020_anc_2050_monit_0 <- pop_2020_anc_2050_monit_0[[1]]
# With rescreens:
#pop_2020_anc_2030_monit_0a <- readRDS(paste0(out_path, "pop_2020_anc_2030_with_rescreen_monit_0_050521.rds"))
#pop_2020_anc_2030_monit_0a <- pop_2020_anc_2030_monit_0a[[1]]
pop_2020_anc_2040_monit_0a <- readRDS(paste0(out_path, "pop_2020_anc_2040_with_rescreen_monit_0_040521.rds"))
pop_2020_anc_2040_monit_0a <- pop_2020_anc_2040_monit_0a[[1]]
pop_2020_anc_2050_monit_0a <- readRDS(paste0(out_path, "pop_2020_anc_2050_with_rescreen_monit_0_040521.rds"))
pop_2020_anc_2050_monit_0a <- pop_2020_anc_2050_monit_0a[[1]]

# ANC strategies
# No monitoring
anc_2020_monit_0 <- readRDS(paste0(out_path, "anc1_out3_screen_2020_monit_0_150121.rds"))
anc_2020_monit_0 <- anc_2020_monit_0[[1]]
anc_2030_monit_0 <- readRDS(paste0(out_path, "anc1_2030_no_rescreen_monit_0_150121.rds"))
anc_2030_monit_0 <- anc_2030_monit_0[[1]]
anc_2040_monit_0 <- readRDS(paste0(out_path, "anc1_2040_no_rescreen_monit_0_150121.rds"))
anc_2040_monit_0 <- anc_2040_monit_0[[1]]
# With rescreen:
anc_2030_monit_0a <- readRDS(paste0(out_path, "anc1_2030_with_rescreen_monit_0_140121.rds"))
anc_2030_monit_0a <- anc_2030_monit_0a[[1]]
anc_2040_monit_0a <- readRDS(paste0(out_path, "anc1_2040_with_rescreen_monit_0_140121.rds"))
anc_2040_monit_0a <- anc_2040_monit_0a[[1]]

# Monitor 5-yearly in <45 year olds
anc_2020_monit_sim7 <- readRDS(paste0(out_path, "anc1_out3_screen_2020_monit_sim7_290321.rds"))
anc_2020_monit_sim7 <- anc_2020_monit_sim7[[1]]
anc_2030_monit_sim7 <- readRDS(paste0(out_path, "anc1_2030_no_rescreen_monit_sim7_290321.rds"))
anc_2030_monit_sim7 <- anc_2030_monit_sim7[[1]]
anc_2040_monit_sim7 <- readRDS(paste0(out_path, "anc1_2040_no_rescreen_monit_sim7_290321.rds"))
anc_2040_monit_sim7 <- anc_2040_monit_sim7[[1]]

# PREVENTION SIMULATIONS
anc_2100_monit_sim7 <- readRDS(paste0(out_path, "anc1_2100_no_rescreen_monit_sim7_090421.rds"))
anc_2100_monit_sim7 <- anc_2100_monit_sim7[[1]]
anc_2100_monit_sim7_bd <- readRDS(paste0(out_path, "anc1_2100_no_rescreen_monit_sim7_bd_090421.rds"))
anc_2100_monit_sim7_bd <- anc_2100_monit_sim7_bd[[1]]
anc_2100_monit_sim7_ppt <- readRDS(paste0(out_path, "anc1_2100_no_rescreen_monit_sim7_bd_ppt_090421.rds"))
anc_2100_monit_sim7_ppt <- anc_2100_monit_sim7_ppt[[1]]

## 1) Cost-effectiveness of ANC compared to repeat population-based screening (no PMTCT) ----
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
  # Repeat population screening:
  cbind(scenario = "monit_sim7_screen_10b_2030",
        assemble_discounted_interactions_for_screening_strategies(out8b_2030_sim7,
                                                                  assessment_object = out8b_2030_monit_0)),
  cbind(scenario = "monit_sim7_screen_10a_2030",
        assemble_discounted_interactions_for_screening_strategies(out8a_2030_sim7,
                                                                  assessment_object = out8a_2030_monit_0)),
  cbind(scenario = "monit_sim7_screen_10b_2040",
        assemble_discounted_interactions_for_screening_strategies(out8b_2040_sim7,
                                                                  assessment_object = out8b_2040_monit_0)),
  cbind(scenario = "monit_sim7_screen_10a_2040",
        assemble_discounted_interactions_for_screening_strategies(out8a_2040_sim7,
                                                                  assessment_object = out8a_2040_monit_0)),
  cbind(scenario = "monit_sim7_screen_10b_2050",
        assemble_discounted_interactions_for_screening_strategies(out8b_2050_sim7,
                                                                  assessment_object = out8b_2050_monit_0)),
  cbind(scenario = "monit_sim7_screen_10a_2050",
        assemble_discounted_interactions_for_screening_strategies(out8a_2050_sim7,
                                                                  assessment_object = out8a_2050_monit_0)),
    # Combination of pop and ANC testing:
  cbind(scenario = "pop_2020_anc_2030_no_rescreen_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2030_sim7,
                                                                  assessment_object = pop_2020_anc_2030_monit_0)),
  cbind(scenario = "pop_2020_anc_2040_no_rescreen_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2040_sim7,
                                                                  assessment_object = pop_2020_anc_2040_monit_0)),
  cbind(scenario = "pop_2020_anc_2050_no_rescreen_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2050_sim7,
                                                                  assessment_object = pop_2020_anc_2050_monit_0)),
 #  cbind(scenario = "pop_2020_anc_2030_with_rescreen_monit_sim7",
 #       assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2030_sim7a,
 #                                                                 assessment_object = pop_2020_anc_2030_monit_0a)),
  cbind(scenario = "pop_2020_anc_2040_with_rescreen_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2040_sim7a,
                                                                  assessment_object = pop_2020_anc_2040_monit_0a)),
  cbind(scenario = "pop_2020_anc_2050_with_rescreen_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2050_sim7a,
                                                                  assessment_object = pop_2020_anc_2050_monit_0a)),
 # ANC strategies
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
cbind(scenario = "anc_2030_with_rescreen_monit_0",
      assemble_discounted_interactions_for_screening_strategies(anc_2030_monit_0a,
                                                                assessment_object = anc_2030_monit_0a)),
cbind(scenario = "anc_2040_no_rescreen_monit_0",
      assemble_discounted_interactions_for_screening_strategies(anc_2040_monit_0,
                                                                assessment_object = anc_2040_monit_0)),
cbind(scenario = "anc_2040_no_rescreen_monit_sim7",
      assemble_discounted_interactions_for_screening_strategies(anc_2040_monit_sim7,
                                                                assessment_object = anc_2040_monit_0)),
cbind(scenario = "anc_2040_with_rescreen_monit_0",
      assemble_discounted_interactions_for_screening_strategies(anc_2040_monit_0a,
                                                                assessment_object = anc_2040_monit_0a))
)


# PY on treatment, DALYS and deaths averted:
object_list <- list(out3_it, monit_out7, out8b_2030_sim7, out8a_2030_sim7,
                    out8b_2040_sim7,out8a_2040_sim7,out8b_2050_sim7,out8a_2050_sim7,
                    pop_2020_anc_2030_sim7, pop_2020_anc_2040_sim7, pop_2020_anc_2050_sim7,
                    #pop_2020_anc_2030_sim7a,
                    pop_2020_anc_2040_sim7a, pop_2020_anc_2050_sim7a,
                    anc_2020_monit_0, anc_2020_monit_sim7,
                    anc_2030_monit_0, anc_2030_monit_sim7, anc_2030_monit_0a,
                    anc_2040_monit_0, anc_2040_monit_sim7, anc_2040_monit_0a)

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
  "anc_2040_no_rescreen_monit_0", "anc_2030_no_rescreen_monit_sim7",
  "anc_2030_with_rescreen_monit_0", "anc_2040_with_rescreen_monit_0")

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
                              "pop_2020_anc_2040_no_rescreen_monit_sim7",
                              "pop_2020_anc_2050_no_rescreen_monit_sim7",
                              #"pop_2020_anc_2030_with_rescreen_monit_sim7",
                              "pop_2020_anc_2040_with_rescreen_monit_sim7",
                              "pop_2020_anc_2050_with_rescreen_monit_sim7")

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
anc_df2 <- subset(anc_df, scenario %in% c("monit_sim7_screen_10a_2050",
                                          "monit_sim7_screen_10b_2050",
                                          "monit_sim7_screen_10b_2040",
                                          "screen_2020_anc_monit_0",
                                          "anc_2030_no_rescreen_monit_0",
                                          "pop_2020_anc_2050_no_rescreen_monit_sim7",
                                          "anc_2040_no_rescreen_monit_0",
                                          "pop_2020_anc_2040_no_rescreen_monit_sim7"))


#anc_df2 <- subset(anc_df, scenario %in% c("monit_sim7_screen_10b_2030",
#                                          "screen_2020_anc_monit_0",
#                                          "pop_2020_anc_2040_no_rescreen_monit_sim7",
#                                          "anc_2030_no_rescreen_monit_0",
#                                          "anc_2040_no_rescreen_monit_0"))

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

# CE plane plot ----

# Group strategies by target population and rescreen (for shape)

pop_with_rescreen <- c("monit_sim7_screen_10a_2030",
                       "monit_sim7_screen_10a_2040",
                       "monit_sim7_screen_10a_2050")
anc_with_rescreen <- c("anc_2030_with_rescreen_monit_0", "anc_2040_with_rescreen_monit_0")
pop_anc_combi_with_rescreen <- c(
                              #"pop_2020_anc_2030_with_rescreen_monit_sim7",
                              "pop_2020_anc_2040_with_rescreen_monit_sim7",
                              "pop_2020_anc_2050_with_rescreen_monit_sim7")

anc_df$rescreen_group <- "no_rescreen"
anc_df$rescreen_group[anc_df$scenario %in% c(pop_with_rescreen,
                                             anc_with_rescreen,
                                             pop_anc_combi_with_rescreen )] <- "with_rescreen"

anc_df$scenario_without_rescreen_cat <- anc_df$scenario

levels(anc_df$scenario_without_rescreen_cat) <- list(
  "ANC 2020-2030, no monitoring" = "anc_2030_no_rescreen_monit_0",
  "ANC 2020-2030, with monitoring"="anc_2030_no_rescreen_monit_sim7",
  "ANC 2020-2030, no monitoring"="anc_2030_with_rescreen_monit_0",
  "ANC 2020-2040, no monitoring"= "anc_2040_no_rescreen_monit_0",
  "ANC 2020-2040, with monitoring"="anc_2040_no_rescreen_monit_sim7",
  "ANC 2020-2040, no monitoring"="anc_2040_with_rescreen_monit_0",
  "All adults 2020+2030, with monitoring" ="monit_sim7_screen_10a_2030" ,
  "All adults 2020+2030+2040, with monitoring" ="monit_sim7_screen_10a_2040" ,
  "All adults 2020+2030+2040+2050, with monitoring"="monit_sim7_screen_10a_2050" ,
  "All adults 2020+2030, with monitoring"="monit_sim7_screen_10b_2030" ,
  "All adults 2020+2030+2040, with monitoring"="monit_sim7_screen_10b_2040" ,
  "All adults 2020+2030+2040+2050, with monitoring" ="monit_sim7_screen_10b_2050" ,
  "No treatment"="No treatment" ,
  "All adults 2020 + ANC 2020-2030, with monitoring" ="pop_2020_anc_2030_no_rescreen_monit_sim7" ,
  "All adults 2020 + ANC 2020-2040, with monitoring"="pop_2020_anc_2040_no_rescreen_monit_sim7",
  "All adults 2020 + ANC 2020-2040, with monitoring"="pop_2020_anc_2040_with_rescreen_monit_sim7",
  "All adults 2020 + ANC 2020-2050, with monitoring"="pop_2020_anc_2050_no_rescreen_monit_sim7",
  "All adults 2020 + ANC 2020-2050, with monitoring"="pop_2020_anc_2050_with_rescreen_monit_sim7",
  "ANC 2020, no monitoring"="screen_2020_anc_monit_0",
  "ANC 2020, with monitoring"= "screen_2020_anc_monit_sim7"  ,
  "All adults 2020, no monitoring" ="screen_2020_monit_0" ,
  "All adults 2020, with monitoring"="screen_2020_monit_sim7")

anc_df_medians <- group_by(anc_df, rescreen_group,
                           scenario_without_rescreen_cat, scenario) %>%
            summarise(dalys_averted=median(dalys_averted),
            total_cost = median(total_cost)) %>%
  filter(scenario != "No treatment")



# THESIS PLOT
# Plot medians only with different shapes for ANC/repeat/combi strategies
ggplot() +
  geom_point(data=anc_df_medians,
             aes(x=dalys_averted, y = total_cost,
                 colour = reorder(scenario_without_rescreen_cat, total_cost)), size = 5) +
    geom_point(data=subset(anc_df_medians, rescreen_group == "with_rescreen"),
               aes(x=dalys_averted, y = total_cost, shape = rescreen_group), size = 5,
               #shape = 21,
               colour="black") +
  scale_colour_manual(values = c("black", brewer.pal(12,"Paired")))+
  scale_shape_manual(labels = c("with_rescreen" = "With individual re-testing"),
                     values=c("with_rescreen" = 21))+
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(shape = guide_legend(order = 0),
         colour = guide_legend(order = 1)) +
  theme_classic() +
  theme(legend.position = c(0.2, 0.75),
        #legend.box = "horizontal",
        legend.title = element_blank())

# Plot with full uncertainty
ggplot(subset(anc_df, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_point(data= group_by(anc_df, scenario) %>%
               summarise(dalys_averted=median(dalys_averted),
                         total_cost = median(total_cost)),
             aes(x=dalys_averted, y = total_cost, group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), size = 5) +
  scale_fill_manual(values = c(rev(brewer.pal(11,"RdYlBu")), brewer.pal(4, "Paired"))) +
  scale_colour_manual("Scenarios",
                      #                      labels = c("screen_2020_monit_0" = "2020 screen, no monitoring",
                      #                                 "screen_2020_monit_10" = "2020 screen, monitor every 10 years",
                      #                                 "monit_0_screen_10b_2030" = "2020+2030 screen, no monitoring",
                      #                                 "monit_0_screen_10b_2040" = "2020+2030+2040 screen, no monitoring",
                      #                                 "monit_10_screen_10b_2030" = "2020+2030 screen, monitor every 10 years",
                      #                                 "monit_10_screen_10b_2040" = "2020+2030+2040 screen, monitor every 10 years"),
                      values = c("black", rev(brewer.pal(11,"RdYlBu")), brewer.pal(4, "Paired"))) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
#  geom_abline(slope=391, linetype = "dashed") +
#  geom_abline(slope=518, linetype = "dashed") +
  theme_classic()

## 2) Elimination analysis with BD and PPT ----
# Reduce chronic infection incidence by 90% from 2015 ----

# New chronic infections in 2015
out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,
                                         -c(1:2)]
# Reduction by given year (status quo):
quantile((out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,
                                                   -c(1:2)]-
            out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2030,
                                                     -c(1:2)])/
           out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,
                                                    -c(1:2)],
         c(0.5,0.025,0.975))
# In 2030: 60% (41-74%)
# Median >=90% in 2054
# Lower percentile >=90% in year 2077

# Treatment in ANC without added prevention:
quantile((anc_2100_monit_sim7$timeseries$total_chronic_infections[
  anc_2100_monit_sim7$timeseries$total_chronic_infections$time==2015,
                                                   -c(1:2)]-
    anc_2100_monit_sim7$timeseries$total_chronic_infections[
      anc_2100_monit_sim7$timeseries$total_chronic_infections$time==2030,
                                                     -c(1:2)])/
    anc_2100_monit_sim7$timeseries$total_chronic_infections[
      anc_2100_monit_sim7$timeseries$total_chronic_infections$time==2015,
                                                    -c(1:2)],
         c(0.5,0.025,0.975))
# In 2030: 62% (42-77%)
# Median >=90% in 2053

# Treatment in ANC with added BD:
quantile((anc_2100_monit_sim7_bd$timeseries$total_chronic_infections[
  anc_2100_monit_sim7_bd$timeseries$total_chronic_infections$time==2015,
  -c(1:2)]-
    anc_2100_monit_sim7_bd$timeseries$total_chronic_infections[
      anc_2100_monit_sim7_bd$timeseries$total_chronic_infections$time==2039,
      -c(1:2)])/
    anc_2100_monit_sim7_bd$timeseries$total_chronic_infections[
      anc_2100_monit_sim7_bd$timeseries$total_chronic_infections$time==2015,
      -c(1:2)],
  c(0.5,0.025,0.975))
# In 2030: 83% (68-89%)
# Median >=90% in 2039

# Treatment in ANC with added BD and PPT:
quantile((anc_2100_monit_sim7_ppt$timeseries$total_chronic_infections[
  anc_2100_monit_sim7_ppt$timeseries$total_chronic_infections$time==2015,
  -c(1:2)]-
    anc_2100_monit_sim7_ppt$timeseries$total_chronic_infections[
      anc_2100_monit_sim7_ppt$timeseries$total_chronic_infections$time==2037,
      -c(1:2)])/
    anc_2100_monit_sim7_ppt$timeseries$total_chronic_infections[
      anc_2100_monit_sim7_ppt$timeseries$total_chronic_infections$time==2015,
      -c(1:2)],
  c(0.5,0.025,0.975))
# In 2030: 85% (69-90%)
# Median >=90% in 2037

# Achieve <=0.1% HBsAg prevalence among <5 year olds ----

# Add out2_carriers

# Treatment only
quantile(anc_2100_monit_sim7$timeseries$prev_under5[
  anc_2100_monit_sim7$timeseries$prev_under5$time==2049,-c(1,2)]*100,
  c(0.5,0.025,0.975))
# In 2030: 0.41% (0.08-1.55)
# Median <= 0.1% in 2049

# Treatment + BD
quantile(anc_2100_monit_sim7_bd$timeseries$prev_under5[
  anc_2100_monit_sim7_bd$timeseries$prev_under5$time==2032,-c(1,2)]*100,
  c(0.5,0.025,0.975))
# In 2030: 0.13% (0.03-0.49)
# Median <= 0.1% in 2032

# Treatment + BD + PPT
quantile(anc_2100_monit_sim7_ppt$timeseries$prev_under5[
  anc_2100_monit_sim7_ppt$timeseries$prev_under5$time==2031,-c(1,2)]*100,
  c(0.5,0.025,0.975))
# In 2030: 0.11% (0.02-0.40)
# Median <= 0.1% in 2031



# Reduce HBV deaths by 65% from 2015 ----

# Status quo:
quantile((out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,
                                           -c(1:2)]-
            out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2030,
                                             -c(1:2)])/
           out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,
                                            -c(1:2)],
         c(0.5,0.025,0.975))
# In 2030: 10% (-7-34%)
# Median >=65% in 2062

# One-off population treatment:
quantile((monit_out7$timeseries$total_hbv_deaths[
  monit_out7$timeseries$total_hbv_deaths$time==2015,-c(1:2)]-
    monit_out7$timeseries$total_hbv_deaths[
      monit_out7$timeseries$total_hbv_deaths$time==2030,-c(1:2)])/
    monit_out7$timeseries$total_hbv_deaths[
      monit_out7$timeseries$total_hbv_deaths$time==2015,-c(1:2)],
  c(0.5,0.025,0.975))
# In 2030: 49% (6-44%) - effect declines again between 2030 and 2035
# Median >=65% in 2060 => actually slightly later than continuous ANC and
# almost same as SQ

# 2 population treatments
quantile((out8b_2030_sim7$timeseries$total_hbv_deaths[
  out8b_2030_sim7$timeseries$total_hbv_deaths$time==2015,-c(1:2)]-
    out8b_2030_sim7$timeseries$total_hbv_deaths[
      out8b_2030_sim7$timeseries$total_hbv_deaths$time==2055,-c(1:2)])/
    out8b_2030_sim7$timeseries$total_hbv_deaths[
      out8b_2030_sim7$timeseries$total_hbv_deaths$time==2015,-c(1:2)],
  c(0.5,0.025,0.975))
# In 2030: 63% (49-74%)
# Median 64% in 2032 but increasing again after
# Median >=65% in 2056

# ANC Treatment only:
quantile((anc_2100_monit_sim7$timeseries$total_hbv_deaths[
  anc_2100_monit_sim7$timeseries$total_hbv_deaths$time==2015,-c(1:2)]-
    anc_2100_monit_sim7$timeseries$total_hbv_deaths[
      anc_2100_monit_sim7$timeseries$total_hbv_deaths$time==2056,-c(1:2)])/
    anc_2100_monit_sim7$timeseries$total_hbv_deaths[
      anc_2100_monit_sim7$timeseries$total_hbv_deaths$time==2015,-c(1:2)],
         c(0.5,0.025,0.975))
# In 2030: 26% (6-44%)
# Median >=65% in 2056

# With BD:
quantile((anc_2100_monit_sim7_bd$timeseries$total_hbv_deaths[
  anc_2100_monit_sim7_bd$timeseries$total_hbv_deaths$time==2015,-c(1:2)]-
    anc_2100_monit_sim7_bd$timeseries$total_hbv_deaths[
      anc_2100_monit_sim7_bd$timeseries$total_hbv_deaths$time==2054,-c(1:2)])/
    anc_2100_monit_sim7_bd$timeseries$total_hbv_deaths[
      anc_2100_monit_sim7_bd$timeseries$total_hbv_deaths$time==2015,-c(1:2)],
  c(0.5,0.025,0.975))
# In 2030: 26% (6-45%)
# Median >=65% in 2054

# With BD+PPT:
quantile((anc_2100_monit_sim7_ppt$timeseries$total_hbv_deaths[
  anc_2100_monit_sim7_ppt$timeseries$total_hbv_deaths$time==2015,-c(1:2)]-
    anc_2100_monit_sim7_ppt$timeseries$total_hbv_deaths[
      anc_2100_monit_sim7_ppt$timeseries$total_hbv_deaths$time==2054,-c(1:2)])/
    anc_2100_monit_sim7_ppt$timeseries$total_hbv_deaths[
      anc_2100_monit_sim7_ppt$timeseries$total_hbv_deaths$time==2015,-c(1:2)],
  c(0.5,0.025,0.975))
# In 2030: 26% (6-45%)
# Median >=65% in 2054




plot(x=anc_2100_monit_sim7$timeseries$total_hbv_deaths$time,
     y = anc_2100_monit_sim7$timeseries$total_hbv_deaths[,3])
plot(x=monit_out7$timeseries$total_hbv_deaths$time,
     y = monit_out7$timeseries$total_hbv_deaths[,3])
plot(x=out8b_2030_sim7$timeseries$total_hbv_deaths$time,
     y = out8b_2030_sim7$timeseries$total_hbv_deaths[,3])
plot(x=anc_2020_monit_0$timeseries$total_hbv_deaths$time,
     y = anc_2020_monit_0$timeseries$total_hbv_deaths[,3])
# Something seems off about the ANC - how can mortality be reduced as much as for the pop one?
# Deaths are still miscalculated here:
plot(x=anc_2040_monit_sim7$timeseries$total_hbv_deaths$time,
     y = anc_2040_monit_sim7$timeseries$total_hbv_deaths[,5])
plot(x=pop_2020_anc_2040_sim7$timeseries$total_hbv_deaths$time,
     y = pop_2020_anc_2040_sim7$timeseries$total_hbv_deaths[,3])

# Reduce HBV mortality to ≤5 per 100 000 (CDA absolute target) ----
quantile(out2$timeseries$total_hbv_deaths_rate[
  out2$timeseries$total_hbv_deaths_rate$time==2040,-c(1:2)]*100000,
         c(0.5,0.025,0.975))
# In 2030: 7.4 (3.7-13.0)
# Median <5 in 2040

# ANC Treatment only:
quantile(anc_2100_monit_sim7$timeseries$total_hbv_deaths_rate[
      anc_2100_monit_sim7$timeseries$total_hbv_deaths_rate$time==2035,-c(1:2)]*100000,
  c(0.5,0.025,0.975))
# In 2030: 6.3 (3.3-10.4)
# Median <5 in 2035

# One-off population treatment:
quantile(monit_out7$timeseries$total_hbv_deaths_rate[
  monit_out7$timeseries$total_hbv_deaths_rate$time==2022,-c(1:2)]*100000,
  c(0.5,0.025,0.975))
# Median <5 in 2022
