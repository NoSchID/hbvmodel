######################################################
### Source code for fitting the Imperial HBV model ###
######################################################

## @knitr part1
### Load packages ----
require(here)
require(truncnorm)
require(binom)  # for calculation of confidence intervals
### Load main script with model, data and other functions ----
source(here("R/imperial_model_main.R"))

### Load calibration data and assign weights ----
# Define my datapoints to fit to
calibration_data_path <- "data/calibration_data"
input_hbsag_dataset <- read.csv(here(calibration_data_path,
                                     "hbsag_prevalence.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)
input_antihbc_dataset <- read.csv(here(calibration_data_path,
                                       "antihbc_prevalence.csv"),
                                  header = TRUE, check.names = FALSE,
                                  stringsAsFactors = FALSE)
input_hbeag_dataset <- read.csv(here(calibration_data_path,
                                     "hbeag_prevalence.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)
input_natural_history_prev_dataset <- read.csv(here(calibration_data_path,
                                                    "natural_history_prevalence.csv"),
                                               header = TRUE, check.names = FALSE,
                                               stringsAsFactors = FALSE)

input_mtct_risk_dataset <- read.csv(here(calibration_data_path,
                                         "mtct_risk.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

input_progression_rates <- read.csv(here(calibration_data_path,
                                         "progression_rates.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

input_mortality_curves <- read.csv(here(calibration_data_path,
                                        "mortality_curves.csv"),
                                   header = TRUE, check.names = FALSE,
                                   stringsAsFactors = FALSE)

input_risk_of_chronic_carriage <- read.csv(here(calibration_data_path,
                                                "risk_of_chronic_carriage.csv"),
                                           header = TRUE, check.names = FALSE,
                                           stringsAsFactors = FALSE)

input_globocan_incidence_data <- read.csv(here(calibration_data_path,
                                               "globocan_incidence_data.csv"),
                                          header = TRUE, check.names = FALSE,
                                          stringsAsFactors = FALSE)

# Input GLOBOCAN incidence data is population-wide
# Need to multiply by PAF from GLCS to get HBV-related HCC only
paf_average <- 0.56
input_globocan_incidence_data$data_value <- input_globocan_incidence_data$data_value*paf_average

input_globocan_mortality_curve <- read.csv(here(calibration_data_path,
                                                "globocan_mortality_curve.csv"),
                                           header = TRUE, check.names = FALSE,
                                           stringsAsFactors = FALSE,
                                           fileEncoding="UTF-8-BOM")

input_odds_ratios <- read.csv(here(calibration_data_path,
                                   "odds_ratios.csv"),
                              header = TRUE, check.names = FALSE,
                              stringsAsFactors = FALSE)

input_gbd_cirrhosis_mortality <- read.csv(here(calibration_data_path,
                                               "gbd_cirrhosis_mortality.csv"),
                                          header = TRUE, check.names = FALSE,
                                          stringsAsFactors = FALSE)


input_liver_disease_demography <- read.csv(here(calibration_data_path,
                                                "liver_disease_demography.csv"),
                                           header = TRUE, check.names = FALSE,
                                           stringsAsFactors = FALSE)


# Calculate 95% CIs for proportions using Wilson method (better for small samples)
# For data points that don't have 95% CIs
# For rate datasets, exclude datasets where rate = 0
# Function:
calculate_95_ci <- function(input_dataset, data_type) {

  if(data_type == "proportion") {

    input_dataset[is.na(input_dataset$ci_lower),] <-
      filter(input_dataset,
             is.na(ci_lower)) %>%
      mutate(ci_lower = replace(ci_lower,
                                values = as.numeric(unlist(binom.confint(data_value*sample_size, sample_size,
                                                           methods = "wilson")["lower"])))) %>%
      mutate(ci_upper = replace(ci_upper,
                                values = as.numeric(unlist(binom.confint(data_value*sample_size, sample_size,
                                                       methods = "wilson")["upper"]))))

    return(input_dataset)

  } else if(data_type == "rate") {

    if(is.null(input_dataset$events_number)) {

      input_dataset[is.na(input_dataset$ci_lower) & input_dataset$data_value != 0,] <-
        filter(input_dataset, is.na(ci_lower) & data_value != 0) %>%
        mutate(ci_lower = replace(ci_lower,
                                  values = data_value/(exp(1.96/sqrt(data_value*py_at_risk))))) %>%
        mutate(ci_upper = replace(ci_upper,
                                  values = data_value*(exp(1.96/sqrt(data_value*py_at_risk)))))
      return(input_dataset)

    } else if(is.null(input_dataset$events_number) == FALSE) {

      input_dataset[is.na(input_dataset$ci_lower) & input_dataset$data_value != 0
                    & is.na(input_dataset$data_value) == FALSE,] <-
        filter(input_dataset, is.na(ci_lower) & data_value != 0 & is.na(data_value) == FALSE) %>%
        mutate(ci_lower = replace(ci_lower,
                                  values = data_value/(exp(1.96/sqrt(events_number))))) %>%
        mutate(ci_upper = replace(ci_upper,
                                  values = data_value*(exp(1.96/sqrt(events_number)))))
      return(input_dataset)

    }

  } else {
    print("data_type can be rate or proportion")
  }

}

# Apply to all proportion datasets
input_hbsag_dataset <- calculate_95_ci(input_hbsag_dataset, "proportion")
input_antihbc_dataset <- calculate_95_ci(input_antihbc_dataset, "proportion")
input_hbeag_dataset <- calculate_95_ci(input_hbeag_dataset, "proportion")
input_natural_history_prev_dataset <- calculate_95_ci(input_natural_history_prev_dataset,
                                                      "proportion")
input_mtct_risk_dataset <- calculate_95_ci(input_mtct_risk_dataset, "proportion")
input_risk_of_chronic_carriage <- calculate_95_ci(input_risk_of_chronic_carriage, "proportion")
input_liver_disease_demography[input_liver_disease_demography$outcome == "hcc_prop_male" |
                                 input_liver_disease_demography$outcome == "cirrhosis_prop_male",] <-
  calculate_95_ci(input_liver_disease_demography[input_liver_disease_demography$outcome == "hcc_prop_male" |
                                                   input_liver_disease_demography$outcome == "cirrhosis_prop_male",],
                  "proportion")

# Apply to rates datasets
input_progression_rates <- calculate_95_ci(input_progression_rates, "rate")  # this doesn't add any CIs
input_globocan_incidence_data <- calculate_95_ci(input_globocan_incidence_data, "rate")

# From HBeAg dataset, need to remove points in liver disease patients! These are already in nat hist prev
input_hbeag_dataset <- input_hbeag_dataset[input_hbeag_dataset$pop_group_clinical != "HBsAg_positive_HCC_patients_cirrhosis_unknown" &
                      input_hbeag_dataset$pop_group_clinical != "HBsAg_positive_cirrhosis_patients_no_HCC",]

# OPTION 1: WEIGHTING SCHEME BASED ON QUALITY

# Add quality weights to datapoints
# First add standard weight of 1 to all:
input_hbsag_dataset$quality_weight <- 1
input_antihbc_dataset$quality_weight <- 1
input_hbeag_dataset$quality_weight <- 1
input_natural_history_prev_dataset$quality_weight <- 1
input_mtct_risk_dataset$quality_weight <- 1
input_progression_rates$quality_weight <- 1
input_mortality_curves$quality_weight <- 1
input_globocan_incidence_data$quality_weight <- 1
input_globocan_mortality_curve$quality_weight <- 1
input_gbd_cirrhosis_mortality$quality_weight <- 1
input_risk_of_chronic_carriage$quality_weight <- 1
input_odds_ratios$quality_weight <- 1
input_liver_disease_demography$quality_weight <- 1

# Downweight specific datapoints according to expert opinion:

# All GBD cirrhosis mortality datapoints as these are modelled estimates
# and there is no data on cirrhosis cases/mortality
input_gbd_cirrhosis_mortality$quality_weight <- 0.5
# Survival curve data: Diarra (A3) and Shimakawa (A4) (feedback from Mark)
input_mortality_curves$quality_weight[input_mortality_curves$id_paper == "A3" |
                                        input_mortality_curves$id_paper == "A4"] <- 0.5
# Cirrhosis prevalence in HCC patients from GLCS (feedback from Maud, difficult to measure)
input_natural_history_prev_dataset$quality_weight[
  input_natural_history_prev_dataset$id_unique == "id_gmb2_1_1999_incident_hcc_cases_from_cc" |
    input_natural_history_prev_dataset$id_unique == "id_gmb2_1_1999_incident_hcc_cases_from_dcc"] <- 0.5

# Sample size/quality weighting strategy: upweight the best in each dataset
# try to represent all ages and time periods
input_hbsag_dataset$quality_weight[input_hbsag_dataset$sample_size >= 250 |
                                     input_hbsag_dataset$geographic_scope == "National"] <- 1.5
input_antihbc_dataset$quality_weight[input_antihbc_dataset$sample_size >= 250 |
                                       input_antihbc_dataset$geographic_scope == "National"] <- 1.5
input_hbeag_dataset$quality_weight[input_hbeag_dataset$sample_size >= 50] <- 1.5
input_natural_history_prev_dataset$quality_weight[input_natural_history_prev_dataset$sample_size >= 250] <- 1.5
# Risk of chronic carriage: upweight those where the age at infection was known and sample size > 25
input_risk_of_chronic_carriage$quality_weight[(input_risk_of_chronic_carriage$id_paper == "6" |
                                 input_risk_of_chronic_carriage$id_paper == "19" |
                                 (input_risk_of_chronic_carriage$id_paper == "GMB6" &
                                    input_risk_of_chronic_carriage$id_proc == "c")) &
                                 input_risk_of_chronic_carriage$sample_size > 25] <- 1.5
# No upweights on odds ratios
# MTCT risk: upweight Kodjoh cause they confirmed the maternal HBsAg persistence post-natally and had highest sample
input_mtct_risk_dataset$quality_weight[input_mtct_risk_dataset$id_paper == "PLUS3"] <- 1.5
# No upweights on liver disease demography or mortality curves
# Progression rates: upweight those with person-years information (only Shimakawa)
# FOI and infection incidence is setting-specific, and Olubuyide mortality rate is neither HBV only nor Gambia
input_progression_rates$quality_weight[is.na(input_progression_rates$py_at_risk) == FALSE] <- 1.5

# Test: give weights so that TRANSMISSION datapoints = NAT HIST datapoints
#input_hbsag_dataset$quality_weight <- 1.37
#input_antihbc_dataset$quality_weight <- 1.37
#input_progression_rates$quality_weight[12:14] <- 1.37
#input_natural_history_prev_dataset$quality_weight[input_natural_history_prev_dataset$id_unique ==
#                                     "id_1_1_1986_incident_chronic_births"] <- 1.37
#input_mtct_risk_dataset$quality_weight <- 1.37

# OPTION 2: WEIGHTING SCHEME BASED ON UPWEIGHTING TRANSMISSION

# Add quality weights to datapoints
# First add standard weight of 1 to all:
input_hbsag_dataset$quality_weight <- 1
input_antihbc_dataset$quality_weight <- 1
input_hbeag_dataset$quality_weight <- 1
input_natural_history_prev_dataset$quality_weight <- 1
input_mtct_risk_dataset$quality_weight <- 1
input_progression_rates$quality_weight <- 1
input_mortality_curves$quality_weight <- 1
input_globocan_incidence_data$quality_weight <- 1
input_globocan_mortality_curve$quality_weight <- 1
input_gbd_cirrhosis_mortality$quality_weight <- 1
input_risk_of_chronic_carriage$quality_weight <- 1
input_odds_ratios$quality_weight <- 1
input_liver_disease_demography$quality_weight <- 1

# Downweight specific datapoints according to expert opinion:

# All GBD cirrhosis mortality datapoints as these are modelled estimates
# and there is no data on cirrhosis cases/mortality
input_gbd_cirrhosis_mortality$quality_weight <- 0.5
# Survival curve data: Diarra (A3) and Shimakawa (A4) (feedback from Mark)
input_mortality_curves$quality_weight[input_mortality_curves$id_paper == "A3" |
                                        input_mortality_curves$id_paper == "A4"] <- 0.5
# Cirrhosis prevalence in HCC patients from GLCS (feedback from Maud, difficult to measure)
input_natural_history_prev_dataset$quality_weight[
  input_natural_history_prev_dataset$id_unique == "id_gmb2_1_1999_incident_hcc_cases_from_cc" |
    input_natural_history_prev_dataset$id_unique == "id_gmb2_1_1999_incident_hcc_cases_from_dcc"] <- 0.5

# Upweight HBsAg, anti-HBc and HBeAg prevalence to 10 and upweight Shimakawa rates to 5
input_hbsag_dataset$quality_weight <- 10
input_antihbc_dataset$quality_weight <- 10
input_hbeag_dataset$quality_weight <- 10
input_progression_rates$quality_weight[input_progression_rates$id_paper == "1"] <- 5
input_natural_history_prev_dataset$quality_weight[input_natural_history_prev_dataset$id_unique ==
                                     "id_1_1_1986_incident_chronic_births"] <- 10

# OPTION 3: BASIC WEIGHTING SCHEME FOR MEDIAN REL DIFF

# Add quality weights to datapoints
# First add standard weight of 1 to all:
input_hbsag_dataset$quality_weight <- 1
input_antihbc_dataset$quality_weight <- 1
input_hbeag_dataset$quality_weight <- 1
input_natural_history_prev_dataset$quality_weight <- 1
input_mtct_risk_dataset$quality_weight <- 1
input_progression_rates$quality_weight <- 1
input_mortality_curves$quality_weight <- 1
input_globocan_incidence_data$quality_weight <- 1
input_globocan_mortality_curve$quality_weight <- 1
input_gbd_cirrhosis_mortality$quality_weight <- 1
input_risk_of_chronic_carriage$quality_weight <- 1
input_odds_ratios$quality_weight <- 1
input_liver_disease_demography$quality_weight <- 1

# Downweight specific datapoints according to expert opinion:

# All GBD cirrhosis mortality datapoints as these are modelled estimates
# and there is no data on cirrhosis cases/mortality
input_gbd_cirrhosis_mortality$quality_weight <- 0.5
# Survival curve data: Diarra (A3) and Shimakawa (A4) (feedback from Mark)
input_mortality_curves$quality_weight[input_mortality_curves$id_paper == "A3" |
                                        input_mortality_curves$id_paper == "A4"] <- 0.5
# Cirrhosis prevalence in HCC patients from GLCS (feedback from Maud, difficult to measure)
input_natural_history_prev_dataset$quality_weight[
  input_natural_history_prev_dataset$id_unique == "id_gmb2_1_1999_incident_hcc_cases_from_cc" |
    input_natural_history_prev_dataset$id_unique == "id_gmb2_1_1999_incident_hcc_cases_from_dcc"] <- 0.5

# Need to change name of this list
calibration_datasets_list <- list(hbsag_prevalence = input_hbsag_dataset,
                                  antihbc_prevalence = input_antihbc_dataset,
                                  hbeag_prevalence = input_hbeag_dataset,
                                  natural_history_prevalence = input_natural_history_prev_dataset,
                                  mtct_risk = input_mtct_risk_dataset,
                                  progression_rates = input_progression_rates,
                                  mortality_curves = input_mortality_curves,
                                  globocan_incidence_data = input_globocan_incidence_data,
                                  globocan_mortality_curve = input_globocan_mortality_curve,
                                  gbd_cirrhosis_mortality = input_gbd_cirrhosis_mortality,
                                  p_chronic = input_risk_of_chronic_carriage,
                                  odds_ratios = input_odds_ratios,
                                  liver_disease_demography = input_liver_disease_demography)


### Load input datasets: to work on ----
#prior_vaccine_efficacy <- read.csv(here("data-raw", "input_infant_vaccine_efficacy.csv"),
#                                   header = TRUE, check.names = FALSE,
#                                   stringsAsFactors = FALSE)
#prior_mtct_risk <- read.csv(here("data-raw", "input_mtct_risk.csv"),
#                            header = TRUE, check.names = FALSE,
#                            stringsAsFactors = FALSE)
#prior_paf_liver_disease <- read.csv(here("data-raw", "input_paf_liver_disease.csv"),
#                                    header = TRUE, check.names = FALSE,
#                                    stringsAsFactors = FALSE)
#prior_progression_rates <- read.csv(here("data-raw", "input_progression_rates.csv"),
#                                    header = TRUE, check.names = FALSE,
#                                    stringsAsFactors = FALSE)

### Define functions for calibration ----

## Main calibration function:
# Runs the model (and shadow models), maps calibration datapoints to matching model output
# (summary statistics) and calculates distance metric between datapoints and model prediction.
# Distance function sums the absolute weighted normalised difference for each data-model pair.
fit_model <- function(..., default_parameter_list, parms_to_change = list(...),
                          scenario = "vacc", data_to_fit) {

  # Simulation with given parameter set ----

  # Simulation parameters for fitting procedure:
  # Simulation starts in 1880, runs for 140 years
  # This is because data starts in 1980, so the population can reach equilibrium with
  # given parameter set
  parameters_for_fit <- generate_parameters(default_parameter_list = default_parameter_list,
                                            parms_to_change = c(parms_to_change,
                                                                sim_starttime = 1880))

  sim <- run_model(sim_duration = 140, init_pop_vector = init_pop_sim,
                   default_parameter_list = parameters_for_fit,
                   parms_to_change = NULL,
                   scenario = "vacc")
  out <- code_model_output_summary(sim)

  # Save population distributions for shadow models
  model_pop1978 <- sim[which(out$time==1978),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1983 <- sim[which(out$time==1983),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1985 <- sim[which(out$time==1985),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1995 <- sim[which(out$time==1995),1:(2*n_infectioncat*n_agecat)+1]
  model_pop2005 <- sim[which(out$time==2005),1:(2*n_infectioncat*n_agecat)+1]
  model_pop2012 <- sim[which(out$time==2012),1:(2*n_infectioncat*n_agecat)+1]

  # Check: Is HBsAg prevalence in 1980 within reasonable limits? ----
  #prev_check <- sum(out$carriers[which(out$time == 1980),])/sum(out$pop[which(out$time == 1980),])

  #if (prev_check < 0.01 | prev_check > 0.5) {
  #  return(list(parameter_set = parameters_for_fit,
  #              error_term = NA))
  #}

  # Matching datasets to model ouput ----

  # Define my datapoints to fit to:
  data <- data_to_fit

  ## 1) GLOBOCAN age- and sex specific HBV-related HCC incidence in The Gambia  ----
  # In 2018
  globocan_hcc_incidence_2018 <- map_incidence_rates(rate_outcome = "hcc_incidence",
                                                     rate_num = "cum_incident_hcc",
                                                     rate_denom = "pop",
                                                     rate_timepoint = 2018,
                                                     rate_dataset = data_to_fit$globocan_incidence_data,
                                                     model_sim = sim, model_out = out)
  # In 1998
  globocan_hcc_incidence_1998 <- map_incidence_rates(rate_outcome = "hcc_incidence",
                                                     rate_num = "cum_incident_hcc",
                                                     rate_denom = "pop",
                                                     rate_timepoint = 1998,
                                                     rate_dataset = data_to_fit$globocan_incidence_data,
                                                     model_sim = sim, model_out = out)

  # In 1988
  globocan_hcc_incidence_1988 <- map_incidence_rates(rate_outcome = "hcc_incidence",
                                                     rate_num = "cum_incident_hcc",
                                                     rate_denom = "pop",
                                                     rate_timepoint = 1988,
                                                     rate_dataset = data_to_fit$globocan_incidence_data,
                                                     model_sim = sim, model_out = out)

  # GLOBOCAN age- and sex-specific HCC mortality rate in The Gambia in 2018:
  globocan_hcc_mortality_2018 <- map_incidence_rates(rate_outcome = "hcc_mortality",
                                                     rate_num = "cum_hcc_deaths",
                                                     rate_denom = "pop",
                                                     rate_timepoint = 2018,
                                                     rate_dataset = data_to_fit$globocan_incidence_data,
                                                     model_sim = sim, model_out = out)

  # Combine HCC incidence and mortality sets and map to GLOBOCAN input data
  globocan_incidence_output <- rbind(globocan_hcc_incidence_2018,
                                     globocan_hcc_incidence_1998,
                                     globocan_hcc_incidence_1988,
                                     globocan_hcc_mortality_2018)

  # For merging, transform factors to character objects
  globocan_incidence_output$outcome <- as.character(globocan_incidence_output$outcome)
  globocan_incidence_output$sex <- as.character(globocan_incidence_output$sex)

  mapped_globocan_incidence <- left_join(data_to_fit$globocan_incidence_data,
                                         globocan_incidence_output,
                                         by = c("outcome", "time", "sex", "age_min", "age_max"))

  ## 2) GBD age- and sex-specific HBV-related cirrhosis death rates in The Gambia ----
  # In 2017:
  gbd_cirrhosis_mortality_2017 <- map_incidence_rates(rate_outcome = "cirrhosis_mortality",
                                                      rate_num = "incident_cirrhosis_deaths",
                                                      rate_denom = "pop",
                                                      rate_timepoint = 2017,
                                                      rate_dataset = data_to_fit$gbd_cirrhosis_mortality,
                                                      model_sim = sim, model_out = out)

  gbd_cirrhosis_mortality_1990 <- map_incidence_rates(rate_outcome = "cirrhosis_mortality",
                                                      rate_num = "incident_cirrhosis_deaths",
                                                      rate_denom = "pop",
                                                      rate_timepoint = 1990,
                                                      rate_dataset = data_to_fit$gbd_cirrhosis_mortality,
                                                      model_sim = sim, model_out = out)

  # Combine cirrhosises mortality datasets and map to GBD input data
  gbd_cirrhosis_mortality_output <- rbind(gbd_cirrhosis_mortality_2017,
                                          gbd_cirrhosis_mortality_1990)

  # For merging, transform factors to character objects
  gbd_cirrhosis_mortality_output$outcome <- as.character(gbd_cirrhosis_mortality_output$outcome)
  gbd_cirrhosis_mortality_output$sex <- as.character(gbd_cirrhosis_mortality_output$sex)

  # Merge with calibration dataset
  mapped_gbd_cirrhosis_mortality <- left_join(data_to_fit$gbd_cirrhosis_mortality,
                                              gbd_cirrhosis_mortality_output,
                                              by = c("outcome", "time", "sex", "age_min", "age_max"))

  ## 3) RISK OF CHRONIC CARRIAGE BY AGE AT INFECTION ----
  model_p_chronic <- data.frame(outcome = "p_chronic",
                                age = ages,
                                model_value = unlist(head(sim[,grepl("^p_chronic_function.",names(sim))],1)))

  model_p_chronic$outcome <- as.character(model_p_chronic$outcome)

  mapped_p_chronic <- full_join(data_to_fit$p_chronic,
                                model_p_chronic, by = c("outcome", "age"))

  ## 4) AGE- AND SEX- SPECIFIC SEROMARKER PREVALENCE ----
  # Mapping matching model output to specific year and age

  # HBsAg prevalence in the population:
  mapped_output_hbsag <- map_seromarker_prev(seromarker_num = "carriers",
                                             seromarker_denom = "pop",
                                             prev_dataset = data_to_fit$hbsag_prevalence,
                                             model_output = out)
  mapped_output_hbsag$outcome <- "HBsAg_prevalence"

  # HBeAg prevalence in chronic carriers:
  mapped_output_hbeag <- map_seromarker_prev(seromarker_num = "eag_positive",
                                             seromarker_denom = "carriers",
                                             prev_dataset = data_to_fit$hbeag_prevalence,
                                             model_output = out)
  mapped_output_hbeag$outcome <- "HBeAg_prevalence"

  # anti-HBc prevalence in the population:
  mapped_output_antihbc <- map_seromarker_prev(seromarker_num = "ever_infected",
                                               seromarker_denom = "pop",
                                               prev_dataset = data_to_fit$antihbc_prevalence,
                                               model_output = out)
  mapped_output_antihbc$outcome <- "Anti_HBc_prevalence"

  # Combine into a single dataframe
  mapped_seromarker_prevalence <- rbind(mapped_output_hbsag, mapped_output_hbeag, mapped_output_antihbc)

  ## 5) NATURAL HISTORY-RELATED PREVALENCE ----
  # Various estimates in chronic carriers/liver disease subgroups
  # Outputs calculated manually and matched to specific agegroup and year

  # Prepare output dataframe:
  model_output_nat_hist <- data.frame(id_unique = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                      age_min = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                      age_max = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                      model_value = rep(0, nrow(data_to_fit$natural_history_prevalence)))

  # Prepare denominators/numerators in common to several studies:

  # Male carriers, aged 27-35.5 years, in 2013
  denom_gmb1_2_2013 <- (sum(out$carriers_male[which(out$time == 2013),(which(ages ==27):which(ages ==35.5))]))

  # Carriers without liver disease, aged 4.5 to 21.5 years, in 1986
  denom_1_1_1986 <- sum(sim[,grepl("^ITf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^IRf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^ICf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
                          sim[,grepl("^ITm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^IRm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^ICm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])

  # Carriers, aged 33 to 47 years, in 2012
  denom_gmb1_1_2012 <- (sum(out$carriers[which(out$time == 2012),(which(ages ==33):which(ages ==47))]))

  # Carriers without liver disease, aged 8 to 95.5 years, in 2013
  denom_1_1_2013 <- sum(sim[,grepl("^ITf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^IRf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^ICf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
                          sim[,grepl("^ITm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^IRm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^ICm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])

  # Carriers with significant liver fibrosis or cirrhosis in 2013
  num_1_1_2013 <- sim[,grepl("^IRf.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^CCf.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^DCCf.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^IRm.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^CCm.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2013),]

  # Incident CC cases in 1999
  denom_gmb15_2_1999 <-
    sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),]+
    sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),]+
    sim[,grepl("^cum_enchb_to_ccf.",names(sim))][which(sim$time == 1999),]+
    sim[,grepl("^cum_enchb_to_ccm.",names(sim))][which(sim$time == 1999),]-
    sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),]-
    sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),]-
    sim[,grepl("^cum_enchb_to_ccf.",names(sim))][which(sim$time == 1998),]-
    sim[,grepl("^cum_enchb_to_ccm.",names(sim))][which(sim$time == 1998),]

  # Incident HBeAg-positive non-cirrhotic HCC cases in 1990
  # (same intermediate timepoint used for 2 studies)
  num_gmb12_gmb15_1990 <-
    sim[,grepl("^cum_it_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ir_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1991),]-
    sim[,grepl("^cum_it_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ir_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1990),]

  # Incident non-cirrhotic HCC cases in 1990
  # (same intermediate timepoint used for 2 studies)
  denom_gmb12_gmb15_1990 <-
    sim[,grepl("^cum_it_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ir_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ic_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ic_to_hccm.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_enchb_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_enchb_to_hccm.",names(sim))][which(sim$time == 1991),] -
    sim[,grepl("^cum_it_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ir_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ic_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ic_to_hccm.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_enchb_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_enchb_to_hccm.",names(sim))][which(sim$time == 1990),]

  # Calculate each datapoint manually and distinguish by minimum and maximum age and
  # a unique ID as follows:
  # id_*study ID*_*group ID*_*datapoint time*_*numerator*
  model_output_nat_hist[1,] <-
    c("id_gmb1_2_2013_ic",
      age_min = 27,
      age_max = 35,
      (sum(sim[,grepl("^ICm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[2,] <-
    c("id_gmb1_2_2013_ir_enchb",
      age_min = 27,
      age_max = 35,
      (sum(sim[,grepl("^IRm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
         sum(sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[3,] <-
    c("id_gmb1_2_2013_cc_dcc",
      age_min = 27,
      age_max = 35,
      (sum(sim[,grepl("^CCm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
         sum(sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[4,] <-
    c("id_gmb1_2_2013_hcc",
      age_min = 27,
      age_max = 35,
      (sum(sim[,grepl("^HCCm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[5,] <-
    c("id_1_1_1986_it",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^ITf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^ITm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[6,] <-
    c("id_1_1_1986_ir",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^IRf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^IRm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[7,] <-
    c("id_1_1_1986_enchb",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[8,] <-
    c("id_1_1_1986_ic",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^ICf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^ICm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[9,] <-
    c("id_1_1_1986_hcc",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^HCCf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^HCCm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        (sum(out$carriers[which(out$time == 1986),(which(ages ==4.5):which(ages ==21.5))])))

  model_output_nat_hist[10,] <-
    c("id_gmb1_1_2012_ic",
      age_min = 33,
      age_max = 47,
      (sum(sim[,grepl("^ICf.",names(sim))][which(sim$time == 2012),which(ages ==33):which(ages ==47)]+
             sim[,grepl("^ICm.",names(sim))][which(sim$time == 2012),which(ages ==33):which(ages ==47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[11,] <-
    c("id_gmb1_1_2012_ir_enchb",
      age_min = 33,
      age_max = 47,
      (sum(sim[,grepl("^IRf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^IRm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[12,] <-
    c("id_gmb1_1_2012_cc_dcc",
      age_min = 33,
      age_max = 47,
      (sum(sim[,grepl("^CCf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^DCCf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^CCm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[13,] <-
    c("id_gmb1_1_2012_hcc",
      age_min = 33,
      age_max = 47,
      (sum(sim[,grepl("^HCCf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^HCCm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[14,] <-
    c("id_1_1_2013_it",
      age_min = 8,
      age_max = 95,
      sum(sim[,grepl("^ITf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            sim[,grepl("^ITm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[15,] <-
    c("id_1_1_2013_ir",
      age_min = 8,
      age_max = 95,
      sum(sim[,grepl("^IRf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            sim[,grepl("^IRm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[16,] <-
    c("id_1_1_2013_enchb",
      age_min = 8,
      age_max = 95,
      sum(sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[17,] <-
    c("id_1_1_2013_ic",
      age_min = 8,
      age_max = 95,
      sum(sim[,grepl("^ICf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            sim[,grepl("^ICm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[18,] <-
    c("id_1_1_2013_cc_dcc",
      age_min = 8,
      age_max = 95,
      (sum(sim[,grepl("^CCf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
             sim[,grepl("^DCCf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
             sim[,grepl("^CCm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
             sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==8):which(ages ==95.5))])))

  model_output_nat_hist[19,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 8,
      age_max = 29,
      (sum(num_1_1_2013[,which(ages ==8):which(ages == 29.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==8):which(ages ==29.5))])))

  model_output_nat_hist[20,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 30,
      age_max = 39,
      (sum(num_1_1_2013[,which(ages ==30):which(ages == 39.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==30):which(ages ==39.5))])))

  model_output_nat_hist[21,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 40,
      age_max = 49,
      (sum(num_1_1_2013[,which(ages ==40):which(ages == 49.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==40):which(ages ==49.5))])))

  model_output_nat_hist[22,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 50,
      age_max = 95,
      (sum(num_1_1_2013[,which(ages ==50):which(ages == 95.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==50):which(ages ==95.5))])))

  # CC prevalence in HCC
  # Approximate as proportion of incident HCC cases in 1999 originating from CC
  model_output_nat_hist[23,] <-
    c("id_gmb2_1_1999_incident_hcc_cases_from_cc",
      age_min = 15,
      age_max = 67,
      (sum(sim[,grepl("^cum_cc_to_hcc.",names(sim))][which(sim$time == 1999),])-
         sum(sim[,grepl("^cum_cc_to_hcc.",names(sim))][which(sim$time == 1998),]))/
        (sum(sim[,grepl("^cum_incident_hcc.",names(sim))][which(sim$time == 1999),])-
           sum(sim[,grepl("^cum_incident_hcc.",names(sim))][which(sim$time == 1998),])))

  # DCC prevalence in HCC, 1999
  # Approximate as proportion of incident HCC cases in 1999 originating from DCC
  model_output_nat_hist[24,] <-
    c("id_gmb2_1_1999_incident_hcc_cases_from_dcc",
      age_min = 15,
      age_max = 67,
      (sum(sim[,grepl("^cum_dcc_to_hcc.",names(sim))][which(sim$time == 1999),])-
         sum(sim[,grepl("^cum_dcc_to_hcc.",names(sim))][which(sim$time == 1998),]))/
        (sum(sim[,grepl("^cum_incident_hcc.",names(sim))][which(sim$time == 1999),])-
           sum(sim[,grepl("^cum_incident_hcc.",names(sim))][which(sim$time == 1998),])))

  # HBeAg prevalence in cirrhosis patients, 1999
  # Approximate as proportion of incident CC cases in 1999 originating from IR
  model_output_nat_hist[25,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 17,
      age_max = 34,
      (sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),which(ages == 17):which(ages == 34.5)])+
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 17):which(ages == 34.5)])-
         sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),which(ages == 17):which(ages == 34.5)])-
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),which(ages == 17):which(ages == 34.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 17):which(ages == 34.5)]))

  model_output_nat_hist[26,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 35,
      age_max = 44,
      (sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),which(ages == 35):which(ages == 44.5)])+
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 35):which(ages == 44.5)])-
         sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),which(ages == 35):which(ages == 44.5)])-
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),which(ages == 35):which(ages == 44.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 35):which(ages == 44.5)]))

  model_output_nat_hist[27,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 45,
      age_max = 54,
      (sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),which(ages == 45):which(ages == 54.5)])+
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 45):which(ages == 54.5)])-
         sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),which(ages == 45):which(ages == 54.5)])-
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),which(ages == 45):which(ages == 54.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 45):which(ages == 54.5)]))

  model_output_nat_hist[28,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 55,
      age_max = 64,
      (sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),which(ages == 55):which(ages == 64.5)])+
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 55):which(ages == 64.5)])-
         sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),which(ages == 55):which(ages == 64.5)])-
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),which(ages == 55):which(ages == 64.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 55):which(ages == 64.5)]))

  # HBeAg prevalence in HCC patients, 1982
  # Approximate as proportion of non-cirrhotic incident HCC cases in 1990
  # originating from IT and IR
  model_output_nat_hist[29,] <-
    c("id_gmb12_1_1982_hbeag_hcc",
      age_min = 15,
      age_max = 49,
      sum(num_gmb12_gmb15_1990[which(ages == 15):which(ages==49.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 15):which(ages == 49.5)]))

  model_output_nat_hist[30,] <-
    c("id_gmb12_1_1982_hbeag_hcc",
      age_min = 50,
      age_max = 72,
      sum(num_gmb12_gmb15_1990[which(ages == 50):which(ages==72.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 50):which(ages == 72.5)]))

  model_output_nat_hist[31,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 17,
      age_max = 34,
      sum(num_gmb12_gmb15_1990[which(ages == 17):which(ages == 34.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 17):which(ages == 34.5)]))

  model_output_nat_hist[32,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 35,
      age_max = 44,
      sum(num_gmb12_gmb15_1990[which(ages == 35):which(ages == 44.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 35):which(ages == 44.5)]))

  model_output_nat_hist[33,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 45,
      age_max = 54,
      sum(num_gmb12_gmb15_1990[which(ages == 45):which(ages == 54.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 45):which(ages == 54.5)]))

  model_output_nat_hist[34,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 55,
      age_max = 64,
      sum(num_gmb12_gmb15_1990[which(ages == 55):which(ages == 64.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 55):which(ages == 64.5)]))

  model_output_nat_hist[35,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 65,
      age_max = 87,
      sum(num_gmb12_gmb15_1990[which(ages == 65):which(ages == 87.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 65):which(ages == 87.5)]))

  # Proportion of chronic carriers attributable to mother-to-child transmission
  model_output_nat_hist[36,] <-
    c("id_1_1_1986_incident_chronic_births",
      age_min = 0,
      age_max = 99.5,
      sum(sim[which(sim$time == 1985),grepl("^cum_chronic_births",names(sim))])/
        (sum(sim[which(sim$time == 1985),grepl("^cum_chronic_births",names(sim))])+
           sum(sim[,grepl("^cum_chronic_infections.",names(sim))][which(sim$time == 1985),])))

  # Turn age into numeric format:
  model_output_nat_hist$age_min <- as.numeric(model_output_nat_hist$age_min)
  model_output_nat_hist$age_max <- as.numeric(model_output_nat_hist$age_max)
  model_output_nat_hist$model_value <- as.numeric(model_output_nat_hist$model_value)

  # Merge with the dataset to fit to:
  mapped_nat_hist_prevalence <- left_join(data_to_fit$natural_history_prevalence,
                                          model_output_nat_hist,
                                          by = c("id_unique", "age_min", "age_max"))

  ## 6) MOTHER-TO-CHILD TRANSMISSION RISK ----
  # Overall mother-to-child transmission risk in given year (irrespective of maternal HBeAg)
  mtct_risk <- data.frame(time = out$time[out$time %in% data_to_fit$mtct_risk$time],
                          model_value = rowSums(out$eag_positive_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba] * parameters_for_fit$mtct_prob_e +
                                                  (out$carriers_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]-
                                                     out$eag_positive_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]) * parameters_for_fit$mtct_prob_s)/
                            rowSums(out$carriers_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]))

  mapped_mtct_risk <- left_join(data_to_fit$mtct_risk, mtct_risk, by = "time")

  ## 7) NATURAL HISTORY PROGRESSION RATES AND MORTALITY CURVES ----
  # Some are calculated from full model output and others require a shadow model

  # Prepare output storage for progression rates and mortality curves:
  progression_rates <- data.frame(outcome = c("shadow1_eag_loss_m",
                                              "shadow1_eag_loss_f",
                                              "shadow1a_hcc_incidence_m",
                                              "shadow1b_hcc_incidence_m",
                                              "shadow1a_hcc_incidence_f",
                                              "shadow1b_hcc_incidence_f",
                                              "shadow1_dcc_incidence",
                                              "shadow1_mortality_m",
                                              "shadow1_mortality_f",
                                              "shadow3_mortality",
                                              "shadow2_sag_loss",
                                              "gmb6_1_a_foi",
                                              "gmb6_1_b_foi",
                                              "gmb7_1_chronic_infection_incidence"),
                                  model_value = 0,
                                  stringsAsFactors = FALSE)

  model_mort_curves <- data.frame(outcome = data_to_fit$mortality_curves$outcome,
                                  model_value = 0,
                                  stringsAsFactors = FALSE)

  ## 6a) TRANSMISSION-REALTED RATES

  # Force of infection (any infection) in 0.5-8.5 year olds (GMB6)
  # Numerator = cumulative incidence of chronic infections and transitions to immune compartment (= any horizontal infection) over follow-up
  # Denominator = person-time in susceptible compartment
  progression_rates[progression_rates$outcome=="gmb6_1_a_foi","model_value"] <-
    (sum(sim[,grepl("^cum_infectionsf.",names(sim))][which(sim$time == 1984),(which(ages == 0.5):which(ages == 8.5))]+
           sim[,grepl("^cum_infectionsm.",names(sim))][which(sim$time == 1984),(which(ages == 0.5):which(ages == 8.5))]) -
       sum(sim[,grepl("^cum_infectionsf.",names(sim))][which(sim$time == 1980),(which(ages == 0.5):which(ages == 8.5))]+
             sim[,grepl("^cum_infectionsm.",names(sim))][which(sim$time == 1980),(which(ages == 0.5):which(ages == 8.5))]))/
    ((sum(sim[,grepl("^Sf.",names(sim))][(which(sim$time == 1980):which(sim$time == 1983.5)),(which(ages == 0.5):which(ages == 8.5))]+
            sim[,grepl("^Sm.",names(sim))][(which(sim$time == 1980):which(sim$time == 1983.5)),(which(ages == 0.5):which(ages == 8.5))]))*dt)

  progression_rates[progression_rates$outcome=="gmb6_1_b_foi","model_value"] <-
    progression_rates[progression_rates$outcome=="gmb6_1_a_foi","model_value"]

  # Incidence rate of chronic infections in 0.5-7.5 year olds (GMB7)
  # Numerator = cumulative incidence of chronic infections over follow-up
  # Denominator = person-time in susceptible compartment
  progression_rates[progression_rates$outcome=="gmb7_1_chronic_infection_incidence","model_value"] <-
    (sum(sim[,grepl("^cum_chronic_infectionsf.",names(sim))][which(sim$time == 1982),(which(ages == 0.5):which(ages == 7.5))]+
           sim[,grepl("^cum_chronic_infectionsm.",names(sim))][which(sim$time == 1982),(which(ages == 0.5):which(ages == 7.5))]) -
       sum(sim[,grepl("^cum_chronic_infectionsf.",names(sim))][which(sim$time == 1981),(which(ages == 0.5):which(ages == 7.5))]+
             sim[,grepl("^cum_chronic_infectionsm.",names(sim))][which(sim$time == 1981),(which(ages == 0.5):which(ages == 7.5))]))/
    ((sum(sim[,grepl("^Sf.",names(sim))][(which(sim$time == 1981):which(sim$time == 1981.5)),(which(ages == 0.5):which(ages == 7.5))]+
            sim[,grepl("^Sm.",names(sim))][(which(sim$time == 1981):which(sim$time == 1981.5)),(which(ages == 0.5):which(ages == 7.5))]))*dt)

  ## 6b) SHADOW MODELS 1 a and b: SHIMAKAWA NATURAL HISTORY COHORT

  # Follow 2 cohorts of chronic carriers - 1 of 0-19 year olds (1a) and
  # one of 20-29 year olds (1b) for 28 years, starting in 1985
  # no one had HCC at baseline and we assume no one had DCC at baseline
  # Compartments: IT, IR, IC, ENCHB, CC
  # Switch off births, migation, betas and vaccination to prevent influx of new carriers
  shadow1a_sim <- run_shadow_model(init_age_from = 0, init_age_to = 19.5, init_sex = "both",
                                   init_compartment_from = 2, init_compartment_to = 6,
                                   shadow_default_parameter_list = parameters_for_fit,
                                   shadow_init = model_pop1985, shadowsim_duration = 29,
                                   shadow_parms_to_change = list(sim_starttime = 1985,
                                                                 births_on = 0,
                                                                 migration_on = 0,
                                                                 b1 = 0,
                                                                 b2 = 0,
                                                                 b3 = 0))
  # Total initial population size in shadow model 1a
  shadow1a_init_pop <- sum(shadow1a_sim[1,1:(2*n_infectioncat*n_agecat)])

  shadow1b_sim <- run_shadow_model(init_age_from = 20, init_age_to = 29.5, init_sex = "both",
                                   init_compartment_from = 2, init_compartment_to = 6,
                                   shadow_default_parameter_list = parameters_for_fit,
                                   shadow_init = model_pop1985, shadowsim_duration = 29,
                                   shadow_parms_to_change = list(sim_starttime = 1985,
                                                                 births_on = 0,
                                                                 migration_on = 0,
                                                                 b1 = 0,
                                                                 b2 = 0,
                                                                 b3 = 0))
  # Total initial population size in shadow model 1b
  shadow1b_init_pop <- sum(shadow1b_sim[1,1:(2*n_infectioncat*n_agecat)])

  # Calculate outputs to fit to:
  shadow1a_out <- code_model_output_summary(shadow1a_sim)
  shadow1b_out <- code_model_output_summary(shadow1b_sim)

  # Total HCC incidence rate per person-year over follow-up
  # Numerator = cumulative number of incident HCC cases over follow-up
  # Denominator = person-timestep at risk * dt (carrier compartments other than HCC)
  # MODEL 1a
  # In women:
  progression_rates[progression_rates$outcome=="shadow1a_hcc_incidence_f","model_value"] <-
    sum(tail(shadow1a_sim[,grepl("^cum_incident_hccf.",names(shadow1a_sim))],1))/
    ((sum(head(shadow1a_out$carriers_female,-1)) -
        sum(head(shadow1a_sim[,grepl("^HCCf.",names(shadow1a_sim))],-1)))*dt)

  # In men:
  progression_rates[progression_rates$outcome=="shadow1a_hcc_incidence_m","model_value"] <-
    sum(tail(shadow1a_sim[,grepl("^cum_incident_hccm.",names(shadow1a_sim))],1))/
    ((sum(head(shadow1a_out$carriers_male,-1)) -
        sum(head(shadow1a_sim[,grepl("^HCCm.",names(shadow1a_sim))],-1)))*dt)

  # MODEL 1b
  # In women:
  progression_rates[progression_rates$outcome=="shadow1b_hcc_incidence_f","model_value"] <-
    sum(tail(shadow1b_sim[,grepl("^cum_incident_hccf.",names(shadow1b_sim))],1))/
    ((sum(head(shadow1b_out$carriers_female,-1)) -
        sum(head(shadow1b_sim[,grepl("^HCCf.",names(shadow1b_sim))],-1)))*dt)

  # In men:
  progression_rates[progression_rates$outcome=="shadow1b_hcc_incidence_m","model_value"] <-
    sum(tail(shadow1b_sim[,grepl("^cum_incident_hccm.",names(shadow1b_sim))],1))/
    ((sum(head(shadow1b_out$carriers_male,-1)) -
        sum(head(shadow1b_sim[,grepl("^HCCm.",names(shadow1b_sim))],-1)))*dt)


  # Total incidence of non-malignant ESLD (DCC) per person-year (both sexes)
  # Numerator = cumulative number of incident DCC cases over follow-up (at last timestep)
  # - cumulative number of transitions from DCC to HCC
  # Denominator = person-timestep at risk * dt (carrier compartments other than DCC and HCC)
  # MODEL 1a:
  shadow1a_dcc_rate <- (sum(tail(shadow1a_sim[,grepl("^cum_incident_dcc.",names(shadow1a_sim))],1))-
                          sum(tail(shadow1a_sim[,grepl("^cum_dcc_to_hcc.",names(shadow1a_sim))],1)))/
    (sum(head(shadow1a_out$carriers,-1)) -
       (sum(head(shadow1a_sim[,grepl("^HCC.",names(shadow1a_sim))],-1))) -
       (sum(head(shadow1a_sim[,grepl("^DCC.",names(shadow1a_sim))],-1)))*dt)

  # MODEL 1b:
  shadow1b_dcc_rate <- (sum(tail(shadow1b_sim[,grepl("^cum_incident_dcc.",names(shadow1b_sim))],1))-
                          sum(tail(shadow1b_sim[,grepl("^cum_dcc_to_hcc.",names(shadow1b_sim))],1)))/
    (sum(head(shadow1b_out$carriers,-1)) -
       (sum(head(shadow1b_sim[,grepl("^HCC.",names(shadow1b_sim))],-1))) -
       (sum(head(shadow1b_sim[,grepl("^DCC.",names(shadow1b_sim))],-1)))*dt)

  # USE AVERAGE ACROSS AGE GROUPS:
  progression_rates[progression_rates$outcome=="shadow1_dcc_incidence","model_value"] <-
    weighted.mean(x = c(shadow1a_dcc_rate, shadow1b_dcc_rate),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  # Incidence rate of eAg loss per person-year
  # Numerator = cumulative number of cases of eAg loss at last timestep
  # Denominator = person-time spent in IT and IR compartments
  # MODEL 1a
  # In women:
  shadow1a_eag_loss_f <-
    (sum(tail(shadow1a_sim[,grepl("^cum_eag_lossf.",names(shadow1a_sim))],1)))/
    ((sum(head(shadow1a_sim[,grepl("^ITf.",names(shadow1a_sim))],-1))+
        sum(head(shadow1a_sim[,grepl("^IRf.",names(shadow1a_sim))],-1)))*dt)

  # In men:
  shadow1a_eag_loss_m <-
    (sum(tail(shadow1a_sim[,grepl("^cum_eag_lossm.",names(shadow1a_sim))],1)))/
    ((sum(head(shadow1a_sim[,grepl("^ITm.",names(shadow1a_sim))],-1))+
        sum(head(shadow1a_sim[,grepl("^IRm.",names(shadow1a_sim))],-1)))*dt)

  # MODEL 1b
  # In women:
  shadow1b_eag_loss_f <-
    (sum(tail(shadow1b_sim[,grepl("^cum_eag_lossf.",names(shadow1b_sim))],1)))/
    ((sum(head(shadow1b_sim[,grepl("^ITf.",names(shadow1b_sim))],-1))+
        sum(head(shadow1b_sim[,grepl("^IRf.",names(shadow1b_sim))],-1)))*dt)

  # In men:
  shadow1b_eag_loss_m <-
    (sum(tail(shadow1b_sim[,grepl("^cum_eag_lossm.",names(shadow1b_sim))],1)))/
    ((sum(head(shadow1b_sim[,grepl("^ITm.",names(shadow1b_sim))],-1))+
        sum(head(shadow1b_sim[,grepl("^IRm.",names(shadow1b_sim))],-1)))*dt)

  # USE AVERAGE ACROSS GROUPS
  progression_rates[progression_rates$outcome=="shadow1_eag_loss_f","model_value"] <-
    weighted.mean(x = c(shadow1a_eag_loss_f, shadow1b_eag_loss_f),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  progression_rates[progression_rates$outcome=="shadow1_eag_loss_m","model_value"] <-
    weighted.mean(x = c(shadow1a_eag_loss_m, shadow1b_eag_loss_m),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  # Mortality rate from any cause (HBV-related deaths + background mortality)
  # This was measured in the whole cohort (no matter where they progressed to)
  # MODEL 1a
  # In women:
  shadow1a_mortality_ratef <- (sum(tail(shadow1a_sim[,grepl("^cum_hbv_deathsf.",names(shadow1a_sim))],1)) +
                                 sum(tail(shadow1a_sim[,grepl("^cum_deathsf.",names(shadow1a_sim))],1)))/
    (sum(head(shadow1a_out$pop_female,-1))*dt)

  # In men:
  shadow1a_mortality_ratem <- (sum(tail(shadow1a_sim[,grepl("^cum_hbv_deathsm.",names(shadow1a_sim))],1)) +
                                 sum(tail(shadow1a_sim[,grepl("^cum_deathsm.",names(shadow1a_sim))],1)))/
    (sum(head(shadow1a_out$pop_male,-1))*dt)

  # MODEL 1b
  # In women:
  shadow1b_mortality_ratef <- (sum(tail(shadow1b_sim[,grepl("^cum_hbv_deathsf.",names(shadow1b_sim))],1)) +
                                 sum(tail(shadow1b_sim[,grepl("^cum_deathsf.",names(shadow1b_sim))],1)))/
    (sum(head(shadow1b_out$pop_female,-1))*dt)

  # In men:
  shadow1b_mortality_ratem <- (sum(tail(shadow1b_sim[,grepl("^cum_hbv_deathsm.",names(shadow1b_sim))],1)) +
                                 sum(tail(shadow1b_sim[,grepl("^cum_deathsm.",names(shadow1b_sim))],1)))/
    (sum(head(shadow1b_out$pop_male,-1))*dt)

  # USE AVERAGE ACROSS AGE GROUPS
  progression_rates[progression_rates$outcome=="shadow1_mortality_f","model_value"] <-
    weighted.mean(x = c(shadow1a_mortality_ratef, shadow1b_mortality_ratef),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))
  progression_rates[progression_rates$outcome=="shadow1_mortality_m","model_value"] <-
    weighted.mean(x = c(shadow1a_mortality_ratem, shadow1b_mortality_ratem),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  ## 6c) SHADOW MODEL 2: COURSAGET CHRONIC CARRIER COHORT

  # Follow a cohort of chronic carriers aged 0-2 years
  # starting in 1978, for 7 years
  # Assume no one had HCC or DCC at baseline based on their age
  # compartments: IT, IR, IC, ENCHB, CC
  # Only represent progression within chronic carriers, so switch off births,
  # migation, betas and vaccination
  shadow2_sim <- run_shadow_model(init_age_from = 0, init_age_to = 2.5, init_sex = "both",
                                  init_compartment_from = 2, init_compartment_to = 6,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1978, shadowsim_duration = 8,
                                  shadow_parms_to_change = list(sim_starttime = 1978,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Overall rate of sAg loss
  # Numerator = cumulative number of incident transitions from IC to immune
  # Denominator = person-time in IC compartment at risk
  progression_rates[progression_rates$outcome=="shadow2_sag_loss","model_value"] <-
    (sum(tail(shadow2_sim[,grepl("^cum_sag_loss.",names(shadow2_sim))],1)))/
    (sum(head(shadow2_sim[,grepl("^IC.",names(shadow2_sim))],-1))*dt)

  ## 6d) SHADOW MODEL 3: OLUBUYIDE

  # Olubuyide (A6) cohort of CC, DCC and HCC patients, followed from 1983 for 6 years
  # Scaling initial population to achieve distribution: 60% HCC, 20% CC, 20% DCC
  shadow3_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1983, shadowsim_duration = 7,
                                  shadow_parms_to_change = list(sim_starttime = 1983,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Output: Mortality rate in CC, DCC and HCC from any cause (HBV-related or background mortality)
  # Numerator = sum of cumulative number of incident deaths between 1989 and 1983
  # Denominator = person-time in compartments at risk (compensated cirrhosis, decompensated cirrhosis, HCC) * dt
  progression_rates[progression_rates$outcome=="shadow3_mortality","model_value"] <-
    (sum(tail(shadow3_sim[,grepl("^cum_hbv_deaths.",names(shadow3_sim))],1)) +
       sum(tail(shadow3_sim[,grepl("^cum_background_deaths_ld.",names(shadow3_sim))],1)))/
    ((sum(head(shadow3_sim[,grepl("^CC.",names(shadow3_sim))],-1)) +
        sum(head(shadow3_sim[,grepl("^DCC.",names(shadow3_sim))],-1)) +
        sum(head(shadow3_sim[,grepl("^HCC.",names(shadow3_sim))],-1)))*dt)

  ## 6e) SHADOW MODEL 4: SHIMAKAWA COMPENSATED CIRRHOSIS COHORT

  # To fit survival curve at time interval 0.5 years
  # Follow a cohort of compensated cirrhosis patients
  # starting in 2012, for 0.5 years
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow4_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 6,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop2012, shadowsim_duration = 1,
                                  shadow_parms_to_change = list(sim_starttime = 2012,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from HBV-related cause or background) after 0.5 years
  # Numerator = sum of incident deaths at next timestep
  # Can use background deaths from all LD patients because DCC and HCC compartments are empty
  # Denominator = number in compensated cirrhosis compartment at first timestep
  model_mort_curves[model_mort_curves$outcome=="shadow4_cum_mortality","model_value"] <-
    (sum(shadow4_sim[,grepl("^cum_hbv_deaths.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]) +
       sum(shadow4_sim[,grepl("^cum_background_deaths_ld.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]))/
    (sum(shadow4_sim[,grepl("^CC.",names(shadow4_sim))][which(shadow4_sim$time == 2012),]))

  # Proportion of deaths due to DCC and HCC
  mapped_nat_hist_prevalence$model_value[
    mapped_nat_hist_prevalence$id_unique == "id_a4_1_2014_shadow_incident_deaths"] <-
    (sum(shadow4_sim[,grepl("^cum_hcc_deaths.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]) +
       sum(shadow4_sim[,grepl("^cum_dcc_deaths.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]))/
    (sum(shadow4_sim[,grepl("^cum_hbv_deaths.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]) +
       sum(shadow4_sim[,grepl("^cum_background_deaths_ld.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]))

  ## 6f) SHADOW MODEL 5: YANG HCC COHORT

  # To fit survival curve at time intervals 0.5, 1 and 1.5 years
  # Follow a cohort of HCC patients, starting in 2012, for 1.5 years
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow5_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 8, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop2012, shadowsim_duration = 2,
                                  shadow_parms_to_change = list(sim_starttime = 2012,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from HCC or background) after 0.5, 1 and 1.5 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because DCC and HCC compartments are empty
  # Denominator = number in compensated cirrhosis compartment at t0

  for (i in 1:3) {  # i = timestep
    model_mort_curves[model_mort_curves$outcome=="shadow5_cum_mortality","model_value"][i] <-
      ((sum(shadow5_sim[,grepl("^cum_hcc_deaths.",names(shadow5_sim))][(i+1),])) +
         (sum(shadow5_sim[,grepl("^cum_background_deaths_ld.",names(shadow5_sim))][(i+1),])))/
      sum(shadow5_sim[,grepl("^HCC.",names(shadow5_sim))][which(shadow5_sim$time == 2012),])
  }

  ## 6g) SHADOW MODEL 6: DIARRA CIRRHOSIS COHORT

  # To fit survival curve at time intervals 0.5 and 1
  # Follow a cohort of CC and DCC patients, starting in 2005, for 1 year
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow6_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 7,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop2005, shadowsim_duration = 1.5,
                                  shadow_parms_to_change = list(sim_starttime = 2005,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from cirrhosis, HCC or background) after 0.5 and 1 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because HCC compartment is empty at baseline
  # Denominator = number in CC and DCC compartments at t0

  for (i in 1:2) {  # i = timestep
    model_mort_curves[model_mort_curves$outcome=="shadow6_cum_mortality","model_value"][i] <-
      ((sum(shadow6_sim[,grepl("^cum_hbv_deaths.",names(shadow6_sim))][(i+1),])) +
         (sum(shadow6_sim[,grepl("^cum_background_deaths_ld.",names(shadow6_sim))][(i+1),])))/
      (sum(shadow6_sim[,grepl("^CC.",names(shadow6_sim))][which(shadow6_sim$time == 2005),]) +
         sum(shadow6_sim[,grepl("^DCC.",names(shadow6_sim))][which(shadow6_sim$time == 2005),]))
  }


  # Cumulative HCC incidence after 0.5 and 1 years
  # Numerator = sum of icnident HCC cases at given timestep
  # Can use HCC cases coming from any compartment because there are no people in the other
  # compartments to transition
  # Denominator = number in CC and DCC compartments at t0

  for (i in 1:2) {  # i = timestep index, so 1 = t0
    model_mort_curves[model_mort_curves$outcome=="shadow6_cum_hcc_incidence","model_value"][i] <-
      (sum(shadow6_sim[,grepl("^cum_incident_hcc.",names(shadow6_sim))][(i+1),]))/
      (sum(shadow6_sim[,grepl("^CC.",names(shadow6_sim))][which(shadow6_sim$time == 2005),]) +
         sum(shadow6_sim[,grepl("^DCC.",names(shadow6_sim))][which(shadow6_sim$time == 2005),]))
  }


  ## 6h) SHADOW MODEL 7: GLOBOCAN SURVIVAL CURVE

  # To fit survival curve at time intervals 1, 3 and 5
  # Follow a cohort of HCC patients, starting in 1993, until 1997.5
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow7_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 8, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1995, shadowsim_duration = 5.5,
                                  shadow_parms_to_change = list(sim_starttime = 1995,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))


  # Calculate output for fitting:
  # Cumulative mortality probability (from HCC or background) after 1, 3 and 5 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because CC and DCC compartments are empty at baseline
  # Denominator = number in HCC compartment at t01
  mapped_globocan_mortality_curve <- data_to_fit$globocan_mortality_curve
  mapped_globocan_mortality_curve$model_value[1] <-
    ((sum(shadow7_sim[,grepl("^cum_hcc_deaths.",names(shadow7_sim))][which(shadow7_sim$time == 1996),])) +
       (sum(shadow7_sim[,grepl("^cum_background_deaths_ld.",names(shadow7_sim))][which(shadow7_sim$time == 1996),])))/
    (sum(shadow7_sim[,grepl("^HCC.",names(shadow7_sim))][which(shadow7_sim$time == 1995),]))

  mapped_globocan_mortality_curve$model_value[2] <-
    ((sum(shadow7_sim[,grepl("^cum_hcc_deaths.",names(shadow7_sim))][which(shadow7_sim$time == 1998),])) +
       (sum(shadow7_sim[,grepl("^cum_background_deaths_ld.",names(shadow7_sim))][which(shadow7_sim$time == 1998),])))/
    (sum(shadow7_sim[,grepl("^HCC.",names(shadow7_sim))][which(shadow7_sim$time == 1995),]))

  mapped_globocan_mortality_curve$model_value[3] <-
    ((sum(shadow7_sim[,grepl("^cum_hcc_deaths.",names(shadow7_sim))][which(shadow7_sim$time == 2000),])) +
       (sum(shadow7_sim[,grepl("^cum_background_deaths_ld.",names(shadow7_sim))][which(shadow7_sim$time == 2000),])))/
    (sum(shadow7_sim[,grepl("^HCC.",names(shadow7_sim))][which(shadow7_sim$time == 1995),]))

  ## Combine model predictions with input datasets
  mapped_progression_rates <- left_join(data_to_fit$progression_rates, progression_rates, by = "outcome")
  mapped_mortality_curves <- cbind(data_to_fit$mortality_curves, model_value = model_mort_curves$model_value)
  mapped_mortality_curves <- bind_rows(mapped_mortality_curves,
                                       mapped_globocan_mortality_curve) # add GLOBOCAN shadow model


  ## 8) ASSOCIATION DATA ----

  # Prepare dataframe
  association_output <- data.frame(outcome = c("odds_ratio_current_hbeag_positivity_and_hcc",
                                               "odds_ratio_current_hbeag_positivity_and_cirrhosis",
                                               "odds_ratio_male_sex_and_significant_liver_fibrosis_or_cirrhosis"))

  # GMB15 GLCS: Association of current HBeAg status and HCC
  # Outcome = HCC, exposure = concurrent HBeAg status
  # Cases = HCC compartment, controls = all other carrier compartments minus CC and DCC
  # because their HBeAg status is not known
  # Cases and controls were recruited in 1999 (midpoint) and were aged 15-83.5 years
  # 83% of participants were male so assume this represents the OR in males
  # Age range was derived from mean age +/- 2*SD
  # Use the proportion of incident non-cirrhotic HCC cases originating from IT and IR
  # as a proxy for the proportion HBeAg-positive among cases and multiply with the total cases

  assoc_hbeag_hcc <- data.frame(
    total_cases = sum(sim[,grepl("^HCCm.",names(sim))][which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]),
    # Exposed controls = IT and IR compartment, which are HCC-free by definition
    exposed_controls = sum(sim[,grepl("^ITm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                             sim[,grepl("^IRm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    # Unexposed controls = IC and ENCHB compartment
    unexposed_controls = sum(sim[,grepl("^ICm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                               sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    prop_exposed_cases = (sum(sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
                            sum(sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]))/
      (sum(sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(sim[,grepl("^cum_ic_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(sim[,grepl("^cum_enchb_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]))
  )

  # Exposed cases = HBeAg-positive HCC patients (approximated)
  assoc_hbeag_hcc$exposed_cases <- assoc_hbeag_hcc$total_cases *
    assoc_hbeag_hcc$prop_exposed_cases
  # Unexposed cases = HBeAg-negative HCC patients
  assoc_hbeag_hcc$unexposed_cases <- assoc_hbeag_hcc$total_cases -
    assoc_hbeag_hcc$exposed_cases

  hbeag_hcc_odds_ratio <- (
    assoc_hbeag_hcc$exposed_cases * assoc_hbeag_hcc$unexposed_controls)/
    (assoc_hbeag_hcc$unexposed_cases * assoc_hbeag_hcc$exposed_controls)

  # GMB15 GLCS: Association of current HBeAg status and cirrhosis
  # Outcome = CC or DCC, exposure = concurrent HBeAg status
  # Cases = CC+DCC compartment, controls = all other carrier compartments minus HCC
  # because their HBeAg status is not known
  # Cases and controls were recruited in 1999 (midpoint) and were aged 15-83.5 years
  # 80% of participants were male so assume this represents the OR in males
  # Age range was derived from mean age +/- 2*SD
  # Use the proportion of incident CC cases originating from IR as a proxy for
  # the proportion HBeAg-positive among cases and multiply this with the total cases

  assoc_hbeag_cirrhosis <- data.frame(
    total_cases = sum(sim[,grepl("^CCm.",names(sim))][which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]+
                        sim[,grepl("^DCCm.",names(sim))][which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]),
    total_pop = sum(out$carriers_male[which(out$time == 1999),which(ages == 15):which(ages == 83.5)])-
      sum(sim[,grepl("^HCCm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    # Exposed controls = IT and IR compartment, which are cirrhosis-free by definition
    exposed_controls = sum(sim[,grepl("^ITm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                             sim[,grepl("^IRm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    prop_exposed_cases = sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])/
      sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
            sim[,grepl("^cum_enchb_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])
  )

  # Exposed cases = HBeAg-positive CC or DCC patients (approximated)
  assoc_hbeag_cirrhosis$exposed_cases <- assoc_hbeag_cirrhosis$total_cases *
    assoc_hbeag_cirrhosis$prop_exposed_cases
  # Unexposed cases = HBeAg-negative CC or DCC patients
  assoc_hbeag_cirrhosis$unexposed_cases <- assoc_hbeag_cirrhosis$total_cases -
    assoc_hbeag_cirrhosis$exposed_cases
  assoc_hbeag_cirrhosis$total_controls <- assoc_hbeag_cirrhosis$total_pop -
    assoc_hbeag_cirrhosis$total_cases
  # Unexposed controls = HBeAg-negative non-cirrhotic carriers
  assoc_hbeag_cirrhosis$unexposed_controls <- assoc_hbeag_cirrhosis$total_controls -
    assoc_hbeag_cirrhosis$exposed_controls

  hbeag_cirrhosis_odds_ratio <- (
    assoc_hbeag_cirrhosis$exposed_cases * assoc_hbeag_cirrhosis$unexposed_controls)/
    (assoc_hbeag_cirrhosis$unexposed_cases * assoc_hbeag_cirrhosis$exposed_controls)

  # 1-1 Shimakawa: Association of sex and significant liver fibrosis or cirrhosis
  # Use outcome as a proxy for CC/DCC, exposure = male sex
  # Cases = CC+DCC compartment, controls = all other carrier compartments minus HCC
  # Cases and controls were recruited in 2013 and were aged 8-95.5 years
  assoc_sex_cirrhosis <- data.frame(
    # Exposed cases = males with CC or DCC
    exposed_cases = sum(sim[,grepl("^CCm.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]+
                          sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    # Unexposed cases = females with CC or DCC
    unexposed_cases = sum(sim[,grepl("^CCf.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]+
                            sim[,grepl("^DCCf.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    total_exposed = sum(out$carriers_male[which(out$time == 2013),which(ages == 8):which(ages == 95.5)])-
      sum(sim[,grepl("^HCCm.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    total_unexposed = sum(out$carriers_female[which(out$time == 2013),which(ages == 8):which(ages == 95.5)])-
      sum(sim[,grepl("^HCCf.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)])
  )

  # Exposed controls = males in IT, IR, IC or ENCHB compartment
  assoc_sex_cirrhosis$exposed_controls <- assoc_sex_cirrhosis$total_exposed -
    assoc_sex_cirrhosis$exposed_cases
  # Unexposed controls = females in IT, IR, IC or ENCHB compartment
  assoc_sex_cirrhosis$unexposed_controls <- assoc_sex_cirrhosis$total_unexposed -
    assoc_sex_cirrhosis$unexposed_cases

  sex_cirrhosis_odds_ratio <-
    (assoc_sex_cirrhosis$exposed_cases * assoc_sex_cirrhosis$unexposed_controls)/
    (assoc_sex_cirrhosis$unexposed_cases * assoc_sex_cirrhosis$exposed_controls)


  # Combine odds ratios in dataframe
  association_output$model_value <- c(hbeag_hcc_odds_ratio,
                                      hbeag_cirrhosis_odds_ratio,
                                      sex_cirrhosis_odds_ratio)

  # Merge with data to fit to (first transform factor to character vector)
  association_output$outcome <- as.character(association_output$outcome)
  mapped_odds_ratios <- left_join(data_to_fit$odds_ratios,
                                  association_output,
                                  by = "outcome")

  ## 9) LIVER DISEASE MEAN AGE AND SEX ----
  liver_disease_demography_output <-
    data.frame(outcome = c("hcc_prop_male",
                           "cirrhosis_prop_male",
                           "hcc_mean_age",
                           "cirrhosis_mean_age"),
               model_value = 0)


  # HCC proportion male
  liver_disease_demography_output$model_value[1] <-
    sum(sim[,grepl("^HCCm.",names(sim))][sim$time == 1999,])/
    sum(sim[,grepl("^HCC.",names(sim))][sim$time == 1999,])

  # Cirrhosis proportion male
  liver_disease_demography_output$model_value[2] <-
    (sum(sim[,grepl("^CCm.",names(sim))][sim$time == 1999,])+
       sum(sim[,grepl("^DCCm.",names(sim))][sim$time == 1999,]))/
    (sum(sim[,grepl("^CC.",names(sim))][sim$time == 1999,])+
       sum(sim[,grepl("^DCC.",names(sim))][sim$time == 1999,]))

  # HCC mean age
  liver_disease_demography_output$model_value[3] <-
    (sum(sim[,grepl("^HCCm.",names(sim))][sim$time == 1999,]*ages)+
       sum(sim[,grepl("^HCCf.",names(sim))][sim$time == 1999,]*ages))/
    sum(sim[,grepl("^HCC.",names(sim))][sim$time == 1999,])

  # Cirrhosis mean age
  liver_disease_demography_output$model_value[4] <-
    (sum(sim[,grepl("^CCm.",names(sim))][sim$time == 1999,]*ages)+
       sum(sim[,grepl("^CCf.",names(sim))][sim$time == 1999,]*ages)+
       sum(sim[,grepl("^DCCm.",names(sim))][sim$time == 1999,]*ages)+
       sum(sim[,grepl("^DCCf.",names(sim))][sim$time == 1999,]*ages))/
    (sum(sim[,grepl("^CC.",names(sim))][sim$time == 1999,])+
       sum(sim[,grepl("^DCC.",names(sim))][sim$time == 1999,]))

  # Map output to calibration dataset
  liver_disease_demography_output$outcome <- as.character(liver_disease_demography_output$outcome)
  mapped_liver_disease_demography <- left_join(data_to_fit$liver_disease_demography,
                                               liver_disease_demography_output,
                                               by = "outcome")


  # Combine all mapped outputs and calculate distance metric ----
  mapped_output_complete <- list(globocan_hcc_incidence = mapped_globocan_incidence,
                                 gbd_cirrhosis_mortality = mapped_gbd_cirrhosis_mortality,
                                 risk_of_chronic_carriage = mapped_p_chronic,
                                 seromarker_prevalence = mapped_seromarker_prevalence,
                                 nat_hist_prevalence = mapped_nat_hist_prevalence,
                                 mtct_risk = mapped_mtct_risk,
                                 progression_rates = mapped_progression_rates,
                                 mortality_curves = mapped_mortality_curves,
                                 odds_ratios = mapped_odds_ratios,
                                 mapped_liver_disease_demography = mapped_liver_disease_demography)

  # Calculate the distance metric/error term:
  #error_term <- calculate_distance_metrics(mapped_output_complete, metric = "mean_rel_diff")

  # Return relevant info (given parameter set, error term and the matched datapoints and outputs)
  #res <- list(parameter_set = parameters_for_fit,
  #            error_term = error_term,
  #            mapped_output = mapped_output_complete,
  #            full_output = out)

  # Remove mapped output where no data point is available (missing data_value)
  mapped_output_for_error <- lapply(mapped_output_complete, function(x) x[!is.na(x$data_value),])
  # Extract model prediction matching datapoints
  model_prediction <- unlist(sapply(mapped_output_for_error, function(x) x$model_value))
  params_matrix <- t(as.matrix(unlist(parameters_for_fit)))

  res <- list(parameter_set = params_matrix,
              model_value = model_prediction)

  return(res)

}

fit_model_full_output <- function(..., default_parameter_list, parms_to_change = list(...),
                      scenario = "vacc", data_to_fit) {

  # Simulation with given parameter set ----

  # Simulation parameters for fitting procedure:
  # Simulation starts in 1880, runs for 140 years
  # This is because data starts in 1980, so the population can reach equilibrium with
  # given parameter set
  parameters_for_fit <- generate_parameters(default_parameter_list = default_parameter_list,
                                            parms_to_change = c(parms_to_change,
                                                                sim_starttime = 1880))

  sim <- run_model(sim_duration = 140, init_pop_vector = init_pop_sim,
                   default_parameter_list = parameters_for_fit,
                   parms_to_change = NULL,
                   scenario = "vacc")
  out <- code_model_output_summary(sim)

  # Save population distributions for shadow models
  model_pop1978 <- sim[which(out$time==1978),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1983 <- sim[which(out$time==1983),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1985 <- sim[which(out$time==1985),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1995 <- sim[which(out$time==1995),1:(2*n_infectioncat*n_agecat)+1]
  model_pop2005 <- sim[which(out$time==2005),1:(2*n_infectioncat*n_agecat)+1]
  model_pop2012 <- sim[which(out$time==2012),1:(2*n_infectioncat*n_agecat)+1]

  # Check: Is HBsAg prevalence in 1980 within reasonable limits? ----
  #prev_check <- sum(out$carriers[which(out$time == 1980),])/sum(out$pop[which(out$time == 1980),])

  #if (prev_check < 0.01 | prev_check > 0.5) {
  #  return(list(parameter_set = parameters_for_fit,
  #              error_term = NA))
  #}

  # Matching datasets to model ouput ----

  # Define my datapoints to fit to:
  data <- data_to_fit

  ## 1) GLOBOCAN age- and sex specific HBV-related HCC incidence in The Gambia  ----
  # In 2018
  globocan_hcc_incidence_2018 <- map_incidence_rates(rate_outcome = "hcc_incidence",
                                                     rate_num = "cum_incident_hcc",
                                                     rate_denom = "pop",
                                                     rate_timepoint = 2018,
                                                     rate_dataset = data_to_fit$globocan_incidence_data,
                                                     model_sim = sim, model_out = out)
  # In 1998
  globocan_hcc_incidence_1998 <- map_incidence_rates(rate_outcome = "hcc_incidence",
                                                     rate_num = "cum_incident_hcc",
                                                     rate_denom = "pop",
                                                     rate_timepoint = 1998,
                                                     rate_dataset = data_to_fit$globocan_incidence_data,
                                                     model_sim = sim, model_out = out)

  # In 1988
  globocan_hcc_incidence_1988 <- map_incidence_rates(rate_outcome = "hcc_incidence",
                                                     rate_num = "cum_incident_hcc",
                                                     rate_denom = "pop",
                                                     rate_timepoint = 1988,
                                                     rate_dataset = data_to_fit$globocan_incidence_data,
                                                     model_sim = sim, model_out = out)

  # GLOBOCAN age- and sex-specific HCC mortality rate in The Gambia in 2018:
  globocan_hcc_mortality_2018 <- map_incidence_rates(rate_outcome = "hcc_mortality",
                                                     rate_num = "cum_hcc_deaths",
                                                     rate_denom = "pop",
                                                     rate_timepoint = 2018,
                                                     rate_dataset = data_to_fit$globocan_incidence_data,
                                                     model_sim = sim, model_out = out)

  # Combine HCC incidence and mortality sets and map to GLOBOCAN input data
  globocan_incidence_output <- rbind(globocan_hcc_incidence_2018,
                                     globocan_hcc_incidence_1998,
                                     globocan_hcc_incidence_1988,
                                     globocan_hcc_mortality_2018)

  # For merging, transform factors to character objects
  globocan_incidence_output$outcome <- as.character(globocan_incidence_output$outcome)
  globocan_incidence_output$sex <- as.character(globocan_incidence_output$sex)

  mapped_globocan_incidence <- left_join(data_to_fit$globocan_incidence_data,
                                         globocan_incidence_output,
                                         by = c("outcome", "time", "sex", "age_min", "age_max"))

  ## 2) GBD age- and sex-specific HBV-related cirrhosis death rates in The Gambia ----
  # In 2017:
  gbd_cirrhosis_mortality_2017 <- map_incidence_rates(rate_outcome = "cirrhosis_mortality",
                                                      rate_num = "incident_cirrhosis_deaths",
                                                      rate_denom = "pop",
                                                      rate_timepoint = 2017,
                                                      rate_dataset = data_to_fit$gbd_cirrhosis_mortality,
                                                      model_sim = sim, model_out = out)

  gbd_cirrhosis_mortality_1990 <- map_incidence_rates(rate_outcome = "cirrhosis_mortality",
                                                      rate_num = "incident_cirrhosis_deaths",
                                                      rate_denom = "pop",
                                                      rate_timepoint = 1990,
                                                      rate_dataset = data_to_fit$gbd_cirrhosis_mortality,
                                                      model_sim = sim, model_out = out)

  # Combine cirrhosises mortality datasets and map to GBD input data
  gbd_cirrhosis_mortality_output <- rbind(gbd_cirrhosis_mortality_2017,
                                          gbd_cirrhosis_mortality_1990)

  # For merging, transform factors to character objects
  gbd_cirrhosis_mortality_output$outcome <- as.character(gbd_cirrhosis_mortality_output$outcome)
  gbd_cirrhosis_mortality_output$sex <- as.character(gbd_cirrhosis_mortality_output$sex)

  # Merge with calibration dataset
  mapped_gbd_cirrhosis_mortality <- left_join(data_to_fit$gbd_cirrhosis_mortality,
                                              gbd_cirrhosis_mortality_output,
                                              by = c("outcome", "time", "sex", "age_min", "age_max"))

  ## 3) RISK OF CHRONIC CARRIAGE BY AGE AT INFECTION ----
  model_p_chronic <- data.frame(outcome = "p_chronic",
                                age = ages,
                                model_value = unlist(head(sim[,grepl("^p_chronic_function.",names(sim))],1)))

  model_p_chronic$outcome <- as.character(model_p_chronic$outcome)

  mapped_p_chronic <- full_join(data_to_fit$p_chronic,
                                model_p_chronic, by = c("outcome", "age"))

  ## 4) AGE- AND SEX- SPECIFIC SEROMARKER PREVALENCE ----
  # Mapping matching model output to specific year and age

  # HBsAg prevalence in the population:
  mapped_output_hbsag <- map_seromarker_prev(seromarker_num = "carriers",
                                             seromarker_denom = "pop",
                                             prev_dataset = data_to_fit$hbsag_prevalence,
                                             model_output = out)
  mapped_output_hbsag$outcome <- "HBsAg_prevalence"

  # HBeAg prevalence in chronic carriers:
  mapped_output_hbeag <- map_seromarker_prev(seromarker_num = "eag_positive",
                                             seromarker_denom = "carriers",
                                             prev_dataset = data_to_fit$hbeag_prevalence,
                                             model_output = out)
  mapped_output_hbeag$outcome <- "HBeAg_prevalence"

  # anti-HBc prevalence in the population:
  mapped_output_antihbc <- map_seromarker_prev(seromarker_num = "ever_infected",
                                               seromarker_denom = "pop",
                                               prev_dataset = data_to_fit$antihbc_prevalence,
                                               model_output = out)
  mapped_output_antihbc$outcome <- "Anti_HBc_prevalence"

  # Combine into a single dataframe
  mapped_seromarker_prevalence <- rbind(mapped_output_hbsag, mapped_output_hbeag, mapped_output_antihbc)

  ## 5) NATURAL HISTORY-RELATED PREVALENCE ----
  # Various estimates in chronic carriers/liver disease subgroups
  # Outputs calculated manually and matched to specific agegroup and year

  # Prepare output dataframe:
  model_output_nat_hist <- data.frame(id_unique = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                      age_min = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                      age_max = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                      model_value = rep(0, nrow(data_to_fit$natural_history_prevalence)))

  # Prepare denominators/numerators in common to several studies:

  # Male carriers, aged 27-35.5 years, in 2013
  denom_gmb1_2_2013 <- (sum(out$carriers_male[which(out$time == 2013),(which(ages ==27):which(ages ==35.5))]))

  # Carriers without liver disease, aged 4.5 to 21.5 years, in 1986
  denom_1_1_1986 <- sum(sim[,grepl("^ITf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^IRf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^ICf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
                          sim[,grepl("^ITm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^IRm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^ICm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])

  # Carriers, aged 33 to 47 years, in 2012
  denom_gmb1_1_2012 <- (sum(out$carriers[which(out$time == 2012),(which(ages ==33):which(ages ==47))]))

  # Carriers without liver disease, aged 8 to 95.5 years, in 2013
  denom_1_1_2013 <- sum(sim[,grepl("^ITf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^IRf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^ICf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
                          sim[,grepl("^ITm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^IRm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^ICm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])

  # Carriers with significant liver fibrosis or cirrhosis in 2013
  num_1_1_2013 <- sim[,grepl("^IRf.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^CCf.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^DCCf.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^IRm.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^CCm.",names(sim))][which(sim$time == 2013),]+
    sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2013),]

  # Incident CC cases in 1999
  denom_gmb15_2_1999 <-
    sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),]+
    sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),]+
    sim[,grepl("^cum_enchb_to_ccf.",names(sim))][which(sim$time == 1999),]+
    sim[,grepl("^cum_enchb_to_ccm.",names(sim))][which(sim$time == 1999),]-
    sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),]-
    sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),]-
    sim[,grepl("^cum_enchb_to_ccf.",names(sim))][which(sim$time == 1998),]-
    sim[,grepl("^cum_enchb_to_ccm.",names(sim))][which(sim$time == 1998),]

  # Incident HBeAg-positive non-cirrhotic HCC cases in 1990
  # (same intermediate timepoint used for 2 studies)
  num_gmb12_gmb15_1990 <-
    sim[,grepl("^cum_it_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ir_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1991),]-
    sim[,grepl("^cum_it_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ir_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1990),]

  # Incident non-cirrhotic HCC cases in 1990
  # (same intermediate timepoint used for 2 studies)
  denom_gmb12_gmb15_1990 <-
    sim[,grepl("^cum_it_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ir_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ic_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_ic_to_hccm.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_enchb_to_hccf.",names(sim))][which(sim$time == 1991),]+
    sim[,grepl("^cum_enchb_to_hccm.",names(sim))][which(sim$time == 1991),] -
    sim[,grepl("^cum_it_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ir_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ic_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_ic_to_hccm.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_enchb_to_hccf.",names(sim))][which(sim$time == 1990),]-
    sim[,grepl("^cum_enchb_to_hccm.",names(sim))][which(sim$time == 1990),]

  # Calculate each datapoint manually and distinguish by minimum and maximum age and
  # a unique ID as follows:
  # id_*study ID*_*group ID*_*datapoint time*_*numerator*
  model_output_nat_hist[1,] <-
    c("id_gmb1_2_2013_ic",
      age_min = 27,
      age_max = 35,
      (sum(sim[,grepl("^ICm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[2,] <-
    c("id_gmb1_2_2013_ir_enchb",
      age_min = 27,
      age_max = 35,
      (sum(sim[,grepl("^IRm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
         sum(sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[3,] <-
    c("id_gmb1_2_2013_cc_dcc",
      age_min = 27,
      age_max = 35,
      (sum(sim[,grepl("^CCm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
         sum(sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[4,] <-
    c("id_gmb1_2_2013_hcc",
      age_min = 27,
      age_max = 35,
      (sum(sim[,grepl("^HCCm.",names(sim))][which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[5,] <-
    c("id_1_1_1986_it",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^ITf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^ITm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[6,] <-
    c("id_1_1_1986_ir",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^IRf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^IRm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[7,] <-
    c("id_1_1_1986_enchb",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[8,] <-
    c("id_1_1_1986_ic",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^ICf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^ICm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[9,] <-
    c("id_1_1_1986_hcc",
      age_min = 4.5,
      age_max = 21.5,
      sum(sim[,grepl("^HCCf.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            sim[,grepl("^HCCm.",names(sim))][which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        (sum(out$carriers[which(out$time == 1986),(which(ages ==4.5):which(ages ==21.5))])))

  model_output_nat_hist[10,] <-
    c("id_gmb1_1_2012_ic",
      age_min = 33,
      age_max = 47,
      (sum(sim[,grepl("^ICf.",names(sim))][which(sim$time == 2012),which(ages ==33):which(ages ==47)]+
             sim[,grepl("^ICm.",names(sim))][which(sim$time == 2012),which(ages ==33):which(ages ==47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[11,] <-
    c("id_gmb1_1_2012_ir_enchb",
      age_min = 33,
      age_max = 47,
      (sum(sim[,grepl("^IRf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^IRm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[12,] <-
    c("id_gmb1_1_2012_cc_dcc",
      age_min = 33,
      age_max = 47,
      (sum(sim[,grepl("^CCf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^DCCf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^CCm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[13,] <-
    c("id_gmb1_1_2012_hcc",
      age_min = 33,
      age_max = 47,
      (sum(sim[,grepl("^HCCf.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             sim[,grepl("^HCCm.",names(sim))][which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[14,] <-
    c("id_1_1_2013_it",
      age_min = 8,
      age_max = 95,
      sum(sim[,grepl("^ITf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            sim[,grepl("^ITm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[15,] <-
    c("id_1_1_2013_ir",
      age_min = 8,
      age_max = 95,
      sum(sim[,grepl("^IRf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            sim[,grepl("^IRm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[16,] <-
    c("id_1_1_2013_enchb",
      age_min = 8,
      age_max = 95,
      sum(sim[,grepl("^ENCHBf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[17,] <-
    c("id_1_1_2013_ic",
      age_min = 8,
      age_max = 95,
      sum(sim[,grepl("^ICf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            sim[,grepl("^ICm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[18,] <-
    c("id_1_1_2013_cc_dcc",
      age_min = 8,
      age_max = 95,
      (sum(sim[,grepl("^CCf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
             sim[,grepl("^DCCf.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
             sim[,grepl("^CCm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
             sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==8):which(ages ==95.5))])))

  model_output_nat_hist[19,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 8,
      age_max = 29,
      (sum(num_1_1_2013[,which(ages ==8):which(ages == 29.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==8):which(ages ==29.5))])))

  model_output_nat_hist[20,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 30,
      age_max = 39,
      (sum(num_1_1_2013[,which(ages ==30):which(ages == 39.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==30):which(ages ==39.5))])))

  model_output_nat_hist[21,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 40,
      age_max = 49,
      (sum(num_1_1_2013[,which(ages ==40):which(ages == 49.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==40):which(ages ==49.5))])))

  model_output_nat_hist[22,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 50,
      age_max = 95,
      (sum(num_1_1_2013[,which(ages ==50):which(ages == 95.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==50):which(ages ==95.5))])))

  # CC prevalence in HCC
  # Approximate as proportion of incident HCC cases in 1999 originating from CC
  model_output_nat_hist[23,] <-
    c("id_gmb2_1_1999_incident_hcc_cases_from_cc",
      age_min = 15,
      age_max = 67,
      (sum(sim[,grepl("^cum_cc_to_hcc.",names(sim))][which(sim$time == 1999),])-
         sum(sim[,grepl("^cum_cc_to_hcc.",names(sim))][which(sim$time == 1998),]))/
        (sum(sim[,grepl("^cum_incident_hcc.",names(sim))][which(sim$time == 1999),])-
           sum(sim[,grepl("^cum_incident_hcc.",names(sim))][which(sim$time == 1998),])))

  # DCC prevalence in HCC, 1999
  # Approximate as proportion of incident HCC cases in 1999 originating from DCC
  model_output_nat_hist[24,] <-
    c("id_gmb2_1_1999_incident_hcc_cases_from_dcc",
      age_min = 15,
      age_max = 67,
      (sum(sim[,grepl("^cum_dcc_to_hcc.",names(sim))][which(sim$time == 1999),])-
         sum(sim[,grepl("^cum_dcc_to_hcc.",names(sim))][which(sim$time == 1998),]))/
        (sum(sim[,grepl("^cum_incident_hcc.",names(sim))][which(sim$time == 1999),])-
           sum(sim[,grepl("^cum_incident_hcc.",names(sim))][which(sim$time == 1998),])))

  # HBeAg prevalence in cirrhosis patients, 1999
  # Approximate as proportion of incident CC cases in 1999 originating from IR
  model_output_nat_hist[25,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 17,
      age_max = 34,
      (sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),which(ages == 17):which(ages == 34.5)])+
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 17):which(ages == 34.5)])-
         sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),which(ages == 17):which(ages == 34.5)])-
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),which(ages == 17):which(ages == 34.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 17):which(ages == 34.5)]))

  model_output_nat_hist[26,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 35,
      age_max = 44,
      (sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),which(ages == 35):which(ages == 44.5)])+
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 35):which(ages == 44.5)])-
         sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),which(ages == 35):which(ages == 44.5)])-
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),which(ages == 35):which(ages == 44.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 35):which(ages == 44.5)]))

  model_output_nat_hist[27,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 45,
      age_max = 54,
      (sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),which(ages == 45):which(ages == 54.5)])+
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 45):which(ages == 54.5)])-
         sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),which(ages == 45):which(ages == 54.5)])-
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),which(ages == 45):which(ages == 54.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 45):which(ages == 54.5)]))

  model_output_nat_hist[28,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 55,
      age_max = 64,
      (sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1999),which(ages == 55):which(ages == 64.5)])+
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 55):which(ages == 64.5)])-
         sum(sim[,grepl("^cum_ir_to_ccf.",names(sim))][which(sim$time == 1998),which(ages == 55):which(ages == 64.5)])-
         sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1998),which(ages == 55):which(ages == 64.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 55):which(ages == 64.5)]))

  # HBeAg prevalence in HCC patients, 1982
  # Approximate as proportion of non-cirrhotic incident HCC cases in 1990
  # originating from IT and IR
  model_output_nat_hist[29,] <-
    c("id_gmb12_1_1982_hbeag_hcc",
      age_min = 15,
      age_max = 49,
      sum(num_gmb12_gmb15_1990[which(ages == 15):which(ages==49.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 15):which(ages == 49.5)]))

  model_output_nat_hist[30,] <-
    c("id_gmb12_1_1982_hbeag_hcc",
      age_min = 50,
      age_max = 72,
      sum(num_gmb12_gmb15_1990[which(ages == 50):which(ages==72.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 50):which(ages == 72.5)]))

  model_output_nat_hist[31,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 17,
      age_max = 34,
      sum(num_gmb12_gmb15_1990[which(ages == 17):which(ages == 34.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 17):which(ages == 34.5)]))

  model_output_nat_hist[32,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 35,
      age_max = 44,
      sum(num_gmb12_gmb15_1990[which(ages == 35):which(ages == 44.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 35):which(ages == 44.5)]))

  model_output_nat_hist[33,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 45,
      age_max = 54,
      sum(num_gmb12_gmb15_1990[which(ages == 45):which(ages == 54.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 45):which(ages == 54.5)]))

  model_output_nat_hist[34,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 55,
      age_max = 64,
      sum(num_gmb12_gmb15_1990[which(ages == 55):which(ages == 64.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 55):which(ages == 64.5)]))

  model_output_nat_hist[35,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 65,
      age_max = 87,
      sum(num_gmb12_gmb15_1990[which(ages == 65):which(ages == 87.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 65):which(ages == 87.5)]))

  # Proportion of chronic carriers attributable to mother-to-child transmission
  model_output_nat_hist[36,] <-
    c("id_1_1_1986_incident_chronic_births",
      age_min = 0,
      age_max = 99.5,
      sum(sim[which(sim$time == 1985),grepl("^cum_chronic_births",names(sim))])/
        (sum(sim[which(sim$time == 1985),grepl("^cum_chronic_births",names(sim))])+
           sum(sim[,grepl("^cum_chronic_infections.",names(sim))][which(sim$time == 1985),])))

  # Turn age into numeric format:
  model_output_nat_hist$age_min <- as.numeric(model_output_nat_hist$age_min)
  model_output_nat_hist$age_max <- as.numeric(model_output_nat_hist$age_max)
  model_output_nat_hist$model_value <- as.numeric(model_output_nat_hist$model_value)

  # Merge with the dataset to fit to:
  mapped_nat_hist_prevalence <- left_join(data_to_fit$natural_history_prevalence,
                                          model_output_nat_hist,
                                          by = c("id_unique", "age_min", "age_max"))

  ## 6) MOTHER-TO-CHILD TRANSMISSION RISK ----
  # Overall mother-to-child transmission risk in given year (irrespective of maternal HBeAg)
  mtct_risk <- data.frame(time = out$time[out$time %in% data_to_fit$mtct_risk$time],
                          model_value = rowSums(out$eag_positive_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba] * parameters_for_fit$mtct_prob_e +
                                                  (out$carriers_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]-
                                                     out$eag_positive_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]) * parameters_for_fit$mtct_prob_s)/
                            rowSums(out$carriers_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]))

  mapped_mtct_risk <- left_join(data_to_fit$mtct_risk, mtct_risk, by = "time")

  ## 7) NATURAL HISTORY PROGRESSION RATES AND MORTALITY CURVES ----
  # Some are calculated from full model output and others require a shadow model

  # Prepare output storage for progression rates and mortality curves:
  progression_rates <- data.frame(outcome = c("shadow1_eag_loss_m",
                                              "shadow1_eag_loss_f",
                                              "shadow1a_hcc_incidence_m",
                                              "shadow1b_hcc_incidence_m",
                                              "shadow1a_hcc_incidence_f",
                                              "shadow1b_hcc_incidence_f",
                                              "shadow1_dcc_incidence",
                                              "shadow1_mortality_m",
                                              "shadow1_mortality_f",
                                              "shadow3_mortality",
                                              "shadow2_sag_loss",
                                              "gmb6_1_a_foi",
                                              "gmb6_1_b_foi",
                                              "gmb7_1_chronic_infection_incidence"),
                                  model_value = 0,
                                  stringsAsFactors = FALSE)

  model_mort_curves <- data.frame(outcome = data_to_fit$mortality_curves$outcome,
                                  model_value = 0,
                                  stringsAsFactors = FALSE)

  ## 6a) TRANSMISSION-REALTED RATES

  # Force of infection (any infection) in 0.5-8.5 year olds (GMB6)
  # Numerator = cumulative incidence of chronic infections and transitions to immune compartment (= any horizontal infection) over follow-up
  # Denominator = person-time in susceptible compartment
  progression_rates[progression_rates$outcome=="gmb6_1_a_foi","model_value"] <-
    (sum(sim[,grepl("^cum_infectionsf.",names(sim))][which(sim$time == 1984),(which(ages == 0.5):which(ages == 8.5))]+
           sim[,grepl("^cum_infectionsm.",names(sim))][which(sim$time == 1984),(which(ages == 0.5):which(ages == 8.5))]) -
       sum(sim[,grepl("^cum_infectionsf.",names(sim))][which(sim$time == 1980),(which(ages == 0.5):which(ages == 8.5))]+
             sim[,grepl("^cum_infectionsm.",names(sim))][which(sim$time == 1980),(which(ages == 0.5):which(ages == 8.5))]))/
    ((sum(sim[,grepl("^Sf.",names(sim))][(which(sim$time == 1980):which(sim$time == 1983.5)),(which(ages == 0.5):which(ages == 8.5))]+
            sim[,grepl("^Sm.",names(sim))][(which(sim$time == 1980):which(sim$time == 1983.5)),(which(ages == 0.5):which(ages == 8.5))]))*dt)

  progression_rates[progression_rates$outcome=="gmb6_1_b_foi","model_value"] <-
    progression_rates[progression_rates$outcome=="gmb6_1_a_foi","model_value"]

  # Incidence rate of chronic infections in 0.5-7.5 year olds (GMB7)
  # Numerator = cumulative incidence of chronic infections over follow-up
  # Denominator = person-time in susceptible compartment
  progression_rates[progression_rates$outcome=="gmb7_1_chronic_infection_incidence","model_value"] <-
    (sum(sim[,grepl("^cum_chronic_infectionsf.",names(sim))][which(sim$time == 1982),(which(ages == 0.5):which(ages == 7.5))]+
           sim[,grepl("^cum_chronic_infectionsm.",names(sim))][which(sim$time == 1982),(which(ages == 0.5):which(ages == 7.5))]) -
       sum(sim[,grepl("^cum_chronic_infectionsf.",names(sim))][which(sim$time == 1981),(which(ages == 0.5):which(ages == 7.5))]+
             sim[,grepl("^cum_chronic_infectionsm.",names(sim))][which(sim$time == 1981),(which(ages == 0.5):which(ages == 7.5))]))/
    ((sum(sim[,grepl("^Sf.",names(sim))][(which(sim$time == 1981):which(sim$time == 1981.5)),(which(ages == 0.5):which(ages == 7.5))]+
            sim[,grepl("^Sm.",names(sim))][(which(sim$time == 1981):which(sim$time == 1981.5)),(which(ages == 0.5):which(ages == 7.5))]))*dt)

  ## 6b) SHADOW MODELS 1 a and b: SHIMAKAWA NATURAL HISTORY COHORT

  # Follow 2 cohorts of chronic carriers - 1 of 0-19 year olds (1a) and
  # one of 20-29 year olds (1b) for 28 years, starting in 1985
  # no one had HCC at baseline and we assume no one had DCC at baseline
  # Compartments: IT, IR, IC, ENCHB, CC
  # Switch off births, migation, betas and vaccination to prevent influx of new carriers
  shadow1a_sim <- run_shadow_model(init_age_from = 0, init_age_to = 19.5, init_sex = "both",
                                   init_compartment_from = 2, init_compartment_to = 6,
                                   shadow_default_parameter_list = parameters_for_fit,
                                   shadow_init = model_pop1985, shadowsim_duration = 29,
                                   shadow_parms_to_change = list(sim_starttime = 1985,
                                                                 births_on = 0,
                                                                 migration_on = 0,
                                                                 b1 = 0,
                                                                 b2 = 0,
                                                                 b3 = 0))
  # Total initial population size in shadow model 1a
  shadow1a_init_pop <- sum(shadow1a_sim[1,1:(2*n_infectioncat*n_agecat)])

  shadow1b_sim <- run_shadow_model(init_age_from = 20, init_age_to = 29.5, init_sex = "both",
                                   init_compartment_from = 2, init_compartment_to = 6,
                                   shadow_default_parameter_list = parameters_for_fit,
                                   shadow_init = model_pop1985, shadowsim_duration = 29,
                                   shadow_parms_to_change = list(sim_starttime = 1985,
                                                                 births_on = 0,
                                                                 migration_on = 0,
                                                                 b1 = 0,
                                                                 b2 = 0,
                                                                 b3 = 0))
  # Total initial population size in shadow model 1b
  shadow1b_init_pop <- sum(shadow1b_sim[1,1:(2*n_infectioncat*n_agecat)])

  # Calculate outputs to fit to:
  shadow1a_out <- code_model_output_summary(shadow1a_sim)
  shadow1b_out <- code_model_output_summary(shadow1b_sim)

  # Total HCC incidence rate per person-year over follow-up
  # Numerator = cumulative number of incident HCC cases over follow-up
  # Denominator = person-timestep at risk * dt (carrier compartments other than HCC)
  # MODEL 1a
  # In women:
  progression_rates[progression_rates$outcome=="shadow1a_hcc_incidence_f","model_value"] <-
    sum(tail(shadow1a_sim[,grepl("^cum_incident_hccf.",names(shadow1a_sim))],1))/
    ((sum(head(shadow1a_out$carriers_female,-1)) -
        sum(head(shadow1a_sim[,grepl("^HCCf.",names(shadow1a_sim))],-1)))*dt)

  # In men:
  progression_rates[progression_rates$outcome=="shadow1a_hcc_incidence_m","model_value"] <-
    sum(tail(shadow1a_sim[,grepl("^cum_incident_hccm.",names(shadow1a_sim))],1))/
    ((sum(head(shadow1a_out$carriers_male,-1)) -
        sum(head(shadow1a_sim[,grepl("^HCCm.",names(shadow1a_sim))],-1)))*dt)

  # MODEL 1b
  # In women:
  progression_rates[progression_rates$outcome=="shadow1b_hcc_incidence_f","model_value"] <-
    sum(tail(shadow1b_sim[,grepl("^cum_incident_hccf.",names(shadow1b_sim))],1))/
    ((sum(head(shadow1b_out$carriers_female,-1)) -
        sum(head(shadow1b_sim[,grepl("^HCCf.",names(shadow1b_sim))],-1)))*dt)

  # In men:
  progression_rates[progression_rates$outcome=="shadow1b_hcc_incidence_m","model_value"] <-
    sum(tail(shadow1b_sim[,grepl("^cum_incident_hccm.",names(shadow1b_sim))],1))/
    ((sum(head(shadow1b_out$carriers_male,-1)) -
        sum(head(shadow1b_sim[,grepl("^HCCm.",names(shadow1b_sim))],-1)))*dt)


  # Total incidence of non-malignant ESLD (DCC) per person-year (both sexes)
  # Numerator = cumulative number of incident DCC cases over follow-up (at last timestep)
  # - cumulative number of transitions from DCC to HCC
  # Denominator = person-timestep at risk * dt (carrier compartments other than DCC and HCC)
  # MODEL 1a:
  shadow1a_dcc_rate <- (sum(tail(shadow1a_sim[,grepl("^cum_incident_dcc.",names(shadow1a_sim))],1))-
                          sum(tail(shadow1a_sim[,grepl("^cum_dcc_to_hcc.",names(shadow1a_sim))],1)))/
    (sum(head(shadow1a_out$carriers,-1)) -
       (sum(head(shadow1a_sim[,grepl("^HCC.",names(shadow1a_sim))],-1))) -
       (sum(head(shadow1a_sim[,grepl("^DCC.",names(shadow1a_sim))],-1)))*dt)

  # MODEL 1b:
  shadow1b_dcc_rate <- (sum(tail(shadow1b_sim[,grepl("^cum_incident_dcc.",names(shadow1b_sim))],1))-
                          sum(tail(shadow1b_sim[,grepl("^cum_dcc_to_hcc.",names(shadow1b_sim))],1)))/
    (sum(head(shadow1b_out$carriers,-1)) -
       (sum(head(shadow1b_sim[,grepl("^HCC.",names(shadow1b_sim))],-1))) -
       (sum(head(shadow1b_sim[,grepl("^DCC.",names(shadow1b_sim))],-1)))*dt)

  # USE AVERAGE ACROSS AGE GROUPS:
  progression_rates[progression_rates$outcome=="shadow1_dcc_incidence","model_value"] <-
    weighted.mean(x = c(shadow1a_dcc_rate, shadow1b_dcc_rate),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  # Incidence rate of eAg loss per person-year
  # Numerator = cumulative number of cases of eAg loss at last timestep
  # Denominator = person-time spent in IT and IR compartments
  # MODEL 1a
  # In women:
  shadow1a_eag_loss_f <-
    (sum(tail(shadow1a_sim[,grepl("^cum_eag_lossf.",names(shadow1a_sim))],1)))/
    ((sum(head(shadow1a_sim[,grepl("^ITf.",names(shadow1a_sim))],-1))+
        sum(head(shadow1a_sim[,grepl("^IRf.",names(shadow1a_sim))],-1)))*dt)

  # In men:
  shadow1a_eag_loss_m <-
    (sum(tail(shadow1a_sim[,grepl("^cum_eag_lossm.",names(shadow1a_sim))],1)))/
    ((sum(head(shadow1a_sim[,grepl("^ITm.",names(shadow1a_sim))],-1))+
        sum(head(shadow1a_sim[,grepl("^IRm.",names(shadow1a_sim))],-1)))*dt)

  # MODEL 1b
  # In women:
  shadow1b_eag_loss_f <-
    (sum(tail(shadow1b_sim[,grepl("^cum_eag_lossf.",names(shadow1b_sim))],1)))/
    ((sum(head(shadow1b_sim[,grepl("^ITf.",names(shadow1b_sim))],-1))+
        sum(head(shadow1b_sim[,grepl("^IRf.",names(shadow1b_sim))],-1)))*dt)

  # In men:
  shadow1b_eag_loss_m <-
    (sum(tail(shadow1b_sim[,grepl("^cum_eag_lossm.",names(shadow1b_sim))],1)))/
    ((sum(head(shadow1b_sim[,grepl("^ITm.",names(shadow1b_sim))],-1))+
        sum(head(shadow1b_sim[,grepl("^IRm.",names(shadow1b_sim))],-1)))*dt)

  # USE AVERAGE ACROSS GROUPS
  progression_rates[progression_rates$outcome=="shadow1_eag_loss_f","model_value"] <-
    weighted.mean(x = c(shadow1a_eag_loss_f, shadow1b_eag_loss_f),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  progression_rates[progression_rates$outcome=="shadow1_eag_loss_m","model_value"] <-
    weighted.mean(x = c(shadow1a_eag_loss_m, shadow1b_eag_loss_m),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  # Mortality rate from any cause (HBV-related deaths + background mortality)
  # This was measured in the whole cohort (no matter where they progressed to)
  # MODEL 1a
  # In women:
  shadow1a_mortality_ratef <- (sum(tail(shadow1a_sim[,grepl("^cum_hbv_deathsf.",names(shadow1a_sim))],1)) +
                                 sum(tail(shadow1a_sim[,grepl("^cum_deathsf.",names(shadow1a_sim))],1)))/
    (sum(head(shadow1a_out$pop_female,-1))*dt)

  # In men:
  shadow1a_mortality_ratem <- (sum(tail(shadow1a_sim[,grepl("^cum_hbv_deathsm.",names(shadow1a_sim))],1)) +
                                 sum(tail(shadow1a_sim[,grepl("^cum_deathsm.",names(shadow1a_sim))],1)))/
    (sum(head(shadow1a_out$pop_male,-1))*dt)

  # MODEL 1b
  # In women:
  shadow1b_mortality_ratef <- (sum(tail(shadow1b_sim[,grepl("^cum_hbv_deathsf.",names(shadow1b_sim))],1)) +
                                 sum(tail(shadow1b_sim[,grepl("^cum_deathsf.",names(shadow1b_sim))],1)))/
    (sum(head(shadow1b_out$pop_female,-1))*dt)

  # In men:
  shadow1b_mortality_ratem <- (sum(tail(shadow1b_sim[,grepl("^cum_hbv_deathsm.",names(shadow1b_sim))],1)) +
                                 sum(tail(shadow1b_sim[,grepl("^cum_deathsm.",names(shadow1b_sim))],1)))/
    (sum(head(shadow1b_out$pop_male,-1))*dt)

  # USE AVERAGE ACROSS AGE GROUPS
  progression_rates[progression_rates$outcome=="shadow1_mortality_f","model_value"] <-
    weighted.mean(x = c(shadow1a_mortality_ratef, shadow1b_mortality_ratef),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))
  progression_rates[progression_rates$outcome=="shadow1_mortality_m","model_value"] <-
    weighted.mean(x = c(shadow1a_mortality_ratem, shadow1b_mortality_ratem),
                  w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  ## 6c) SHADOW MODEL 2: COURSAGET CHRONIC CARRIER COHORT

  # Follow a cohort of chronic carriers aged 0-2 years
  # starting in 1978, for 7 years
  # Assume no one had HCC or DCC at baseline based on their age
  # compartments: IT, IR, IC, ENCHB, CC
  # Only represent progression within chronic carriers, so switch off births,
  # migation, betas and vaccination
  shadow2_sim <- run_shadow_model(init_age_from = 0, init_age_to = 2.5, init_sex = "both",
                                  init_compartment_from = 2, init_compartment_to = 6,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1978, shadowsim_duration = 8,
                                  shadow_parms_to_change = list(sim_starttime = 1978,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Overall rate of sAg loss
  # Numerator = cumulative number of incident transitions from IC to immune
  # Denominator = person-time in IC compartment at risk
  progression_rates[progression_rates$outcome=="shadow2_sag_loss","model_value"] <-
    (sum(tail(shadow2_sim[,grepl("^cum_sag_loss.",names(shadow2_sim))],1)))/
    (sum(head(shadow2_sim[,grepl("^IC.",names(shadow2_sim))],-1))*dt)

  ## 6d) SHADOW MODEL 3: OLUBUYIDE

  # Olubuyide (A6) cohort of CC, DCC and HCC patients, followed from 1983 for 6 years
  # Scaling initial population to achieve distribution: 60% HCC, 20% CC, 20% DCC
  shadow3_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1983, shadowsim_duration = 7,
                                  shadow_parms_to_change = list(sim_starttime = 1983,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Output: Mortality rate in CC, DCC and HCC from any cause (HBV-related or background mortality)
  # Numerator = sum of cumulative number of incident deaths between 1989 and 1983
  # Denominator = person-time in compartments at risk (compensated cirrhosis, decompensated cirrhosis, HCC) * dt
  progression_rates[progression_rates$outcome=="shadow3_mortality","model_value"] <-
    (sum(tail(shadow3_sim[,grepl("^cum_hbv_deaths.",names(shadow3_sim))],1)) +
       sum(tail(shadow3_sim[,grepl("^cum_background_deaths_ld.",names(shadow3_sim))],1)))/
    ((sum(head(shadow3_sim[,grepl("^CC.",names(shadow3_sim))],-1)) +
        sum(head(shadow3_sim[,grepl("^DCC.",names(shadow3_sim))],-1)) +
        sum(head(shadow3_sim[,grepl("^HCC.",names(shadow3_sim))],-1)))*dt)

  ## 6e) SHADOW MODEL 4: SHIMAKAWA COMPENSATED CIRRHOSIS COHORT

  # To fit survival curve at time interval 0.5 years
  # Follow a cohort of compensated cirrhosis patients
  # starting in 2012, for 0.5 years
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow4_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 6,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop2012, shadowsim_duration = 1,
                                  shadow_parms_to_change = list(sim_starttime = 2012,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from HBV-related cause or background) after 0.5 years
  # Numerator = sum of incident deaths at next timestep
  # Can use background deaths from all LD patients because DCC and HCC compartments are empty
  # Denominator = number in compensated cirrhosis compartment at first timestep
  model_mort_curves[model_mort_curves$outcome=="shadow4_cum_mortality","model_value"] <-
    (sum(shadow4_sim[,grepl("^cum_hbv_deaths.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]) +
       sum(shadow4_sim[,grepl("^cum_background_deaths_ld.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]))/
    (sum(shadow4_sim[,grepl("^CC.",names(shadow4_sim))][which(shadow4_sim$time == 2012),]))

  # Proportion of deaths due to DCC and HCC
  mapped_nat_hist_prevalence$model_value[
    mapped_nat_hist_prevalence$id_unique == "id_a4_1_2014_shadow_incident_deaths"] <-
    (sum(shadow4_sim[,grepl("^cum_hcc_deaths.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]) +
       sum(shadow4_sim[,grepl("^cum_dcc_deaths.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]))/
    (sum(shadow4_sim[,grepl("^cum_hbv_deaths.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]) +
       sum(shadow4_sim[,grepl("^cum_background_deaths_ld.",names(shadow4_sim))][which(shadow4_sim$time == 2012.5),]))

  ## 6f) SHADOW MODEL 5: YANG HCC COHORT

  # To fit survival curve at time intervals 0.5, 1 and 1.5 years
  # Follow a cohort of HCC patients, starting in 2012, for 1.5 years
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow5_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 8, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop2012, shadowsim_duration = 2,
                                  shadow_parms_to_change = list(sim_starttime = 2012,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from HCC or background) after 0.5, 1 and 1.5 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because DCC and HCC compartments are empty
  # Denominator = number in compensated cirrhosis compartment at t0

  for (i in 1:3) {  # i = timestep
    model_mort_curves[model_mort_curves$outcome=="shadow5_cum_mortality","model_value"][i] <-
      ((sum(shadow5_sim[,grepl("^cum_hcc_deaths.",names(shadow5_sim))][(i+1),])) +
         (sum(shadow5_sim[,grepl("^cum_background_deaths_ld.",names(shadow5_sim))][(i+1),])))/
      sum(shadow5_sim[,grepl("^HCC.",names(shadow5_sim))][which(shadow5_sim$time == 2012),])
  }

  ## 6g) SHADOW MODEL 6: DIARRA CIRRHOSIS COHORT

  # To fit survival curve at time intervals 0.5 and 1
  # Follow a cohort of CC and DCC patients, starting in 2005, for 1 year
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow6_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 7,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop2005, shadowsim_duration = 1.5,
                                  shadow_parms_to_change = list(sim_starttime = 2005,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from cirrhosis, HCC or background) after 0.5 and 1 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because HCC compartment is empty at baseline
  # Denominator = number in CC and DCC compartments at t0

  for (i in 1:2) {  # i = timestep
    model_mort_curves[model_mort_curves$outcome=="shadow6_cum_mortality","model_value"][i] <-
      ((sum(shadow6_sim[,grepl("^cum_hbv_deaths.",names(shadow6_sim))][(i+1),])) +
         (sum(shadow6_sim[,grepl("^cum_background_deaths_ld.",names(shadow6_sim))][(i+1),])))/
      (sum(shadow6_sim[,grepl("^CC.",names(shadow6_sim))][which(shadow6_sim$time == 2005),]) +
         sum(shadow6_sim[,grepl("^DCC.",names(shadow6_sim))][which(shadow6_sim$time == 2005),]))
  }


  # Cumulative HCC incidence after 0.5 and 1 years
  # Numerator = sum of icnident HCC cases at given timestep
  # Can use HCC cases coming from any compartment because there are no people in the other
  # compartments to transition
  # Denominator = number in CC and DCC compartments at t0

  for (i in 1:2) {  # i = timestep index, so 1 = t0
    model_mort_curves[model_mort_curves$outcome=="shadow6_cum_hcc_incidence","model_value"][i] <-
      (sum(shadow6_sim[,grepl("^cum_incident_hcc.",names(shadow6_sim))][(i+1),]))/
      (sum(shadow6_sim[,grepl("^CC.",names(shadow6_sim))][which(shadow6_sim$time == 2005),]) +
         sum(shadow6_sim[,grepl("^DCC.",names(shadow6_sim))][which(shadow6_sim$time == 2005),]))
  }


  ## 6h) SHADOW MODEL 7: GLOBOCAN SURVIVAL CURVE

  # To fit survival curve at time intervals 1, 3 and 5
  # Follow a cohort of HCC patients, starting in 1993, until 1997.5
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow7_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 8, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1995, shadowsim_duration = 5.5,
                                  shadow_parms_to_change = list(sim_starttime = 1995,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))


  # Calculate output for fitting:
  # Cumulative mortality probability (from HCC or background) after 1, 3 and 5 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because CC and DCC compartments are empty at baseline
  # Denominator = number in HCC compartment at t01
  mapped_globocan_mortality_curve <- data_to_fit$globocan_mortality_curve
  mapped_globocan_mortality_curve$model_value[1] <-
    ((sum(shadow7_sim[,grepl("^cum_hcc_deaths.",names(shadow7_sim))][which(shadow7_sim$time == 1996),])) +
       (sum(shadow7_sim[,grepl("^cum_background_deaths_ld.",names(shadow7_sim))][which(shadow7_sim$time == 1996),])))/
    (sum(shadow7_sim[,grepl("^HCC.",names(shadow7_sim))][which(shadow7_sim$time == 1995),]))

  mapped_globocan_mortality_curve$model_value[2] <-
    ((sum(shadow7_sim[,grepl("^cum_hcc_deaths.",names(shadow7_sim))][which(shadow7_sim$time == 1998),])) +
       (sum(shadow7_sim[,grepl("^cum_background_deaths_ld.",names(shadow7_sim))][which(shadow7_sim$time == 1998),])))/
    (sum(shadow7_sim[,grepl("^HCC.",names(shadow7_sim))][which(shadow7_sim$time == 1995),]))

  mapped_globocan_mortality_curve$model_value[3] <-
    ((sum(shadow7_sim[,grepl("^cum_hcc_deaths.",names(shadow7_sim))][which(shadow7_sim$time == 2000),])) +
       (sum(shadow7_sim[,grepl("^cum_background_deaths_ld.",names(shadow7_sim))][which(shadow7_sim$time == 2000),])))/
    (sum(shadow7_sim[,grepl("^HCC.",names(shadow7_sim))][which(shadow7_sim$time == 1995),]))

  ## Combine model predictions with input datasets
  mapped_progression_rates <- left_join(data_to_fit$progression_rates, progression_rates, by = "outcome")
  mapped_mortality_curves <- cbind(data_to_fit$mortality_curves, model_value = model_mort_curves$model_value)
  mapped_mortality_curves <- bind_rows(mapped_mortality_curves,
                                       mapped_globocan_mortality_curve) # add GLOBOCAN shadow model


  ## 8) ASSOCIATION DATA ----

  # Prepare dataframe
  association_output <- data.frame(outcome = c("odds_ratio_current_hbeag_positivity_and_hcc",
                                               "odds_ratio_current_hbeag_positivity_and_cirrhosis",
                                               "odds_ratio_male_sex_and_significant_liver_fibrosis_or_cirrhosis"))

  # GMB15 GLCS: Association of current HBeAg status and HCC
  # Outcome = HCC, exposure = concurrent HBeAg status
  # Cases = HCC compartment, controls = all other carrier compartments minus CC and DCC
  # because their HBeAg status is not known
  # Cases and controls were recruited in 1999 (midpoint) and were aged 15-83.5 years
  # 83% of participants were male so assume this represents the OR in males
  # Age range was derived from mean age +/- 2*SD
  # Use the proportion of incident non-cirrhotic HCC cases originating from IT and IR
  # as a proxy for the proportion HBeAg-positive among cases and multiply with the total cases

  assoc_hbeag_hcc <- data.frame(
    total_cases = sum(sim[,grepl("^HCCm.",names(sim))][which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]),
    # Exposed controls = IT and IR compartment, which are HCC-free by definition
    exposed_controls = sum(sim[,grepl("^ITm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                             sim[,grepl("^IRm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    # Unexposed controls = IC and ENCHB compartment
    unexposed_controls = sum(sim[,grepl("^ICm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                               sim[,grepl("^ENCHBm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    prop_exposed_cases = (sum(sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
                            sum(sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]))/
      (sum(sim[,grepl("^cum_it_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(sim[,grepl("^cum_ir_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(sim[,grepl("^cum_ic_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(sim[,grepl("^cum_enchb_to_hccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]))
  )

  # Exposed cases = HBeAg-positive HCC patients (approximated)
  assoc_hbeag_hcc$exposed_cases <- assoc_hbeag_hcc$total_cases *
    assoc_hbeag_hcc$prop_exposed_cases
  # Unexposed cases = HBeAg-negative HCC patients
  assoc_hbeag_hcc$unexposed_cases <- assoc_hbeag_hcc$total_cases -
    assoc_hbeag_hcc$exposed_cases

  hbeag_hcc_odds_ratio <- (
    assoc_hbeag_hcc$exposed_cases * assoc_hbeag_hcc$unexposed_controls)/
    (assoc_hbeag_hcc$unexposed_cases * assoc_hbeag_hcc$exposed_controls)

  # GMB15 GLCS: Association of current HBeAg status and cirrhosis
  # Outcome = CC or DCC, exposure = concurrent HBeAg status
  # Cases = CC+DCC compartment, controls = all other carrier compartments minus HCC
  # because their HBeAg status is not known
  # Cases and controls were recruited in 1999 (midpoint) and were aged 15-83.5 years
  # 80% of participants were male so assume this represents the OR in males
  # Age range was derived from mean age +/- 2*SD
  # Use the proportion of incident CC cases originating from IR as a proxy for
  # the proportion HBeAg-positive among cases and multiply this with the total cases

  assoc_hbeag_cirrhosis <- data.frame(
    total_cases = sum(sim[,grepl("^CCm.",names(sim))][which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]+
                        sim[,grepl("^DCCm.",names(sim))][which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]),
    total_pop = sum(out$carriers_male[which(out$time == 1999),which(ages == 15):which(ages == 83.5)])-
      sum(sim[,grepl("^HCCm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    # Exposed controls = IT and IR compartment, which are cirrhosis-free by definition
    exposed_controls = sum(sim[,grepl("^ITm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                             sim[,grepl("^IRm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    prop_exposed_cases = sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])/
      sum(sim[,grepl("^cum_ir_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
            sim[,grepl("^cum_enchb_to_ccm.",names(sim))][which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])
  )

  # Exposed cases = HBeAg-positive CC or DCC patients (approximated)
  assoc_hbeag_cirrhosis$exposed_cases <- assoc_hbeag_cirrhosis$total_cases *
    assoc_hbeag_cirrhosis$prop_exposed_cases
  # Unexposed cases = HBeAg-negative CC or DCC patients
  assoc_hbeag_cirrhosis$unexposed_cases <- assoc_hbeag_cirrhosis$total_cases -
    assoc_hbeag_cirrhosis$exposed_cases
  assoc_hbeag_cirrhosis$total_controls <- assoc_hbeag_cirrhosis$total_pop -
    assoc_hbeag_cirrhosis$total_cases
  # Unexposed controls = HBeAg-negative non-cirrhotic carriers
  assoc_hbeag_cirrhosis$unexposed_controls <- assoc_hbeag_cirrhosis$total_controls -
    assoc_hbeag_cirrhosis$exposed_controls

  hbeag_cirrhosis_odds_ratio <- (
    assoc_hbeag_cirrhosis$exposed_cases * assoc_hbeag_cirrhosis$unexposed_controls)/
    (assoc_hbeag_cirrhosis$unexposed_cases * assoc_hbeag_cirrhosis$exposed_controls)

  # 1-1 Shimakawa: Association of sex and significant liver fibrosis or cirrhosis
  # Use outcome as a proxy for CC/DCC, exposure = male sex
  # Cases = CC+DCC compartment, controls = all other carrier compartments minus HCC
  # Cases and controls were recruited in 2013 and were aged 8-95.5 years
  assoc_sex_cirrhosis <- data.frame(
    # Exposed cases = males with CC or DCC
    exposed_cases = sum(sim[,grepl("^CCm.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]+
                          sim[,grepl("^DCCm.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    # Unexposed cases = females with CC or DCC
    unexposed_cases = sum(sim[,grepl("^CCf.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]+
                            sim[,grepl("^DCCf.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    total_exposed = sum(out$carriers_male[which(out$time == 2013),which(ages == 8):which(ages == 95.5)])-
      sum(sim[,grepl("^HCCm.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    total_unexposed = sum(out$carriers_female[which(out$time == 2013),which(ages == 8):which(ages == 95.5)])-
      sum(sim[,grepl("^HCCf.",names(sim))][which(sim$time == 2013),which(ages == 8):which(ages == 95.5)])
  )

  # Exposed controls = males in IT, IR, IC or ENCHB compartment
  assoc_sex_cirrhosis$exposed_controls <- assoc_sex_cirrhosis$total_exposed -
    assoc_sex_cirrhosis$exposed_cases
  # Unexposed controls = females in IT, IR, IC or ENCHB compartment
  assoc_sex_cirrhosis$unexposed_controls <- assoc_sex_cirrhosis$total_unexposed -
    assoc_sex_cirrhosis$unexposed_cases

  sex_cirrhosis_odds_ratio <-
    (assoc_sex_cirrhosis$exposed_cases * assoc_sex_cirrhosis$unexposed_controls)/
    (assoc_sex_cirrhosis$unexposed_cases * assoc_sex_cirrhosis$exposed_controls)


  # Combine odds ratios in dataframe
  association_output$model_value <- c(hbeag_hcc_odds_ratio,
                                      hbeag_cirrhosis_odds_ratio,
                                      sex_cirrhosis_odds_ratio)

  # Merge with data to fit to (first transform factor to character vector)
  association_output$outcome <- as.character(association_output$outcome)
  mapped_odds_ratios <- left_join(data_to_fit$odds_ratios,
                                  association_output,
                                  by = "outcome")

  ## 9) LIVER DISEASE MEAN AGE AND SEX ----
  liver_disease_demography_output <-
    data.frame(outcome = c("hcc_prop_male",
                           "cirrhosis_prop_male",
                           "hcc_mean_age",
                           "cirrhosis_mean_age"),
               model_value = 0)


  # HCC proportion male
  liver_disease_demography_output$model_value[1] <-
    sum(sim[,grepl("^HCCm.",names(sim))][sim$time == 1999,])/
    sum(sim[,grepl("^HCC.",names(sim))][sim$time == 1999,])

  # Cirrhosis proportion male
  liver_disease_demography_output$model_value[2] <-
    (sum(sim[,grepl("^CCm.",names(sim))][sim$time == 1999,])+
       sum(sim[,grepl("^DCCm.",names(sim))][sim$time == 1999,]))/
    (sum(sim[,grepl("^CC.",names(sim))][sim$time == 1999,])+
       sum(sim[,grepl("^DCC.",names(sim))][sim$time == 1999,]))

  # HCC mean age
  liver_disease_demography_output$model_value[3] <-
    (sum(sim[,grepl("^HCCm.",names(sim))][sim$time == 1999,]*ages)+
       sum(sim[,grepl("^HCCf.",names(sim))][sim$time == 1999,]*ages))/
    sum(sim[,grepl("^HCC.",names(sim))][sim$time == 1999,])

  # Cirrhosis mean age
  liver_disease_demography_output$model_value[4] <-
    (sum(sim[,grepl("^CCm.",names(sim))][sim$time == 1999,]*ages)+
       sum(sim[,grepl("^CCf.",names(sim))][sim$time == 1999,]*ages)+
       sum(sim[,grepl("^DCCm.",names(sim))][sim$time == 1999,]*ages)+
       sum(sim[,grepl("^DCCf.",names(sim))][sim$time == 1999,]*ages))/
    (sum(sim[,grepl("^CC.",names(sim))][sim$time == 1999,])+
       sum(sim[,grepl("^DCC.",names(sim))][sim$time == 1999,]))

  # Map output to calibration dataset
  liver_disease_demography_output$outcome <- as.character(liver_disease_demography_output$outcome)
  mapped_liver_disease_demography <- left_join(data_to_fit$liver_disease_demography,
                                               liver_disease_demography_output,
                                               by = "outcome")


  # Combine all mapped outputs and calculate distance metric ----
  mapped_output_complete <- list(globocan_hcc_incidence = mapped_globocan_incidence,
                                 gbd_cirrhosis_mortality = mapped_gbd_cirrhosis_mortality,
                                 risk_of_chronic_carriage = mapped_p_chronic,
                                 seromarker_prevalence = mapped_seromarker_prevalence,
                                 nat_hist_prevalence = mapped_nat_hist_prevalence,
                                 mtct_risk = mapped_mtct_risk,
                                 progression_rates = mapped_progression_rates,
                                 mortality_curves = mapped_mortality_curves,
                                 odds_ratios = mapped_odds_ratios,
                                 mapped_liver_disease_demography = mapped_liver_disease_demography)

  # Calculate the distance metric/error term:
  error_term <- calculate_distance_metrics(mapped_output_complete, metric = "mean_rel_diff")

  # Return relevant info (given parameter set, error term and the matched datapoints and outputs)
  res <- list(parameter_set = parameters_for_fit,
              error_term = error_term,
              mapped_output = mapped_output_complete,
              full_output = out)

  # Remove mapped output where no data point is available (missing data_value)
  #mapped_output_for_error <- lapply(mapped_output_complete, function(x) x[!is.na(x$data_value),])
  # Extract model prediction matching datapoints
  #model_prediction <- unlist(sapply(mapped_output_for_error, function(x) x$model_value))
  #params_matrix <- t(as.matrix(unlist(parameters_for_fit)))

  #res <- list(parameter_set = params_matrix,
  #            model_value = model_prediction)

  return(res)

}

## Sub-functions called within main calibration function:

# Function to simulate shadow models within main calibration function
run_shadow_model <- function(init_age_from, init_age_to, init_sex,
                             init_compartment_from, init_compartment_to,
                             shadow_init, shadowsim_duration,
                             shadow_default_parameter_list,
                             shadow_parms_to_change) {

  # Input arguments explanation
  # The shadow model allows to follow a cohort of individuals in:
  # given compartments: init_compartment_from to init_compartment_to
  # > Compartment numbers are: 1 = Susceptible, 2 = IT, 3 = IR, 4 = IC, 5 = ENCHB,
  # > 6 = CC, 7 = DCC, 8 = HCC, 9 = Recovered
  # given age groups: init_age_from to init_age_to
  # of a given sex: init_sex takes "both", "female" or "male" as input
  # To mimic an empirical study, the cohort has a given starting year (defined in shadow_parms_to_change)
  # Define shadow_init is the population distribution simulated by the overall fitting algorithm in the starting year
  # shadowsim_duration is the length of follow-up defined in the study + 1
  # shadow_default_parameter_list should be the parameters used by overall fitting algorithm
  # shadow_parms_to_change needs to be adapted to turn off any influx into the cohort to observe
  # That involves turning off at least births and migration
  # Vaccination is automatically switched off

  # Set up cohort to follow:
  # Define index for groups of interest (infection compartments and age groups)
  if (init_sex == "both") {
    shadow_index <- c(t(mapply(seq, from = (which(ages %in% seq(init_age_from,init_age_to,da))+(init_compartment_from-1)*n_agecat),
                               to = n_agecat*(init_compartment_to), by =n_agecat)),  # for women
                      t(mapply(seq, from = (which(ages %in% seq(init_age_from,init_age_to,da))+(9+init_compartment_from-1)*n_agecat),
                               to= 2*n_agecat*(n_infectioncat-(n_infectioncat-init_compartment_to)/2), by =n_agecat)))  # for men
  } else if (init_sex == "female") {
    shadow_index <- c(t(mapply(seq, from = (which(ages %in% seq(init_age_from,init_age_to,da))+(init_compartment_from-1)*n_agecat),
                               to = n_agecat*(init_compartment_to), by =n_agecat)))
  } else if (init_sex == "male") {
    shadow_index <- c(t(mapply(seq, from = (which(ages %in% seq(init_age_from,init_age_to,da))+(9+init_compartment_from-1)*n_agecat),
                               to= 2*n_agecat*(n_infectioncat-(n_infectioncat-init_compartment_to)/2), by =n_agecat)))
  } else  {
    return(print("init_sex can be both, female or male"))
  }

  # Get init pop from full model output in a given year (defined in function call)
  shadow_init_pop <- c(shadow_init, output_storage)
  shadow_init_pop <- unlist(shadow_init_pop)

  # Set all age groups and compartments other than those to follow to nearly 0
  # (model does not seem to run if it is exactly 0)
  shadow_init_pop[-c(shadow_index,((2*n_agecat*n_infectioncat+1):length(shadow_init_pop)))] <-
    0.000000001

  # Special cases: representing a cohort of liver disease patients
  # In studies of cirrhosis, decompensated patients are overrepresented
  # To simulate these cohorts, increase DCC patients so the cohort is
  # 50% CC and 50% DCC patients
  # If the cohort includes CC, DCC and HCC, distribute as 20, 20 and 60% (Olubuyide study A6)
  if (init_compartment_from == 6 & init_compartment_to == 7) {

    cc_index <- grep("^CC", names(shadow_init_pop))  # index for CC compartments
    dcc_index <- grep("^DCC", names(shadow_init_pop))  # index for DCC compartments
    cc_to_dcc_ratio <- sum(shadow_init_pop[cc_index])/sum(shadow_init_pop[dcc_index])  # current ratio of CC to DCC patients
    shadow_init_pop[dcc_index] <- shadow_init_pop[dcc_index] * cc_to_dcc_ratio  # Increase DCC to obtain 50-50 ratio

  } else if (init_compartment_from == 6 & init_compartment_to == 8) {

    cc_index <- grep("^CC", names(shadow_init_pop))  # index for CC compartments
    dcc_index <- grep("^DCC", names(shadow_init_pop))  # index for DCC compartments
    hcc_index <- grep("^HCC", names(shadow_init_pop))  # index for HCC compartments
    cc_to_dcc_ratio <- sum(shadow_init_pop[cc_index])/sum(shadow_init_pop[dcc_index])  # current ratio of CC to DCC patients
    shadow_init_pop[dcc_index] <- shadow_init_pop[dcc_index] * cc_to_dcc_ratio  # Increase DCC to obtain 50-50 ratio
    dcc_to_hcc_ratio <- sum(shadow_init_pop[dcc_index])/sum(shadow_init_pop[hcc_index])  # current ratio of HCC to DCC patients
    shadow_init_pop[hcc_index] <- shadow_init_pop[hcc_index] * dcc_to_hcc_ratio * 3  # Increase HCC to obtain 60-40 ratio HCC to cirrhosis
  }

  # Set possible transitions within cohort to follow (no new additions) in function call

  # Simulate the cohort
  shadow_sim <- run_model(sim_duration = shadowsim_duration, init_pop_vector = shadow_init_pop,
                          default_parameter_list = shadow_default_parameter_list,
                          parms_to_change = shadow_parms_to_change,
                          scenario = "no_vacc")

  return(shadow_sim)

}

# Function to map age- and sex-specific seroprevalence data to
# matching model output by year, sex and age
# Used to fit to HBsAg, HBeAg and anti-HBc seroprevalence data
map_seromarker_prev <- function(seromarker_num, seromarker_denom, prev_dataset, model_output) {

  seromarker_num_female <- paste0(seromarker_num, "_female")
  seromarker_num_male <- paste0(seromarker_num, "_male")
  seromarker_denom_female <- paste0(seromarker_denom, "_female")
  seromarker_denom_male <- paste0(seromarker_denom, "_male")

  # Filter the output dataset by the year of interest and calculate prevalence for all ages

  # For data from both sexes:
  model_prev_subset_both <- data.frame(time = model_output$time[model_output$time %in%
                                                                  prev_dataset$time[prev_dataset$sex == "Mixed"]],
                                       sex = "Mixed",
                                       prev = model_output[[seromarker_num]][model_output$time %in%
                                                                               prev_dataset$time[prev_dataset$sex == "Mixed"],]/
                                         model_output[[seromarker_denom]][model_output$time %in% prev_dataset$time[prev_dataset$sex == "Mixed"],])


  # For women:
  model_prev_subset_female <- data.frame(time = model_output$time[model_output$time %in%
                                                                    prev_dataset$time[prev_dataset$sex == "Female"]],
                                         sex = "Female",
                                         prev = model_output[[seromarker_num_female]][model_output$time %in%
                                                                                        prev_dataset$time[prev_dataset$sex == "Female"],]/
                                           model_output[[seromarker_denom_female]][model_output$time %in% prev_dataset$time[prev_dataset$sex == "Female"],])


  # For men:
  model_prev_subset_male <- data.frame(time = model_output$time[model_output$time %in%
                                                                  prev_dataset$time[prev_dataset$sex == "Male"]],
                                       sex = "Male",
                                       prev = model_output[[seromarker_num_male]][model_output$time %in%
                                                                                    prev_dataset$time[prev_dataset$sex == "Male"],]/
                                         model_output[[seromarker_denom_male]][model_output$time %in% prev_dataset$time[prev_dataset$sex == "Male"],])


  # Assign all columns the same names to combine into 1 dataframe
  names(model_prev_subset_both) <- c("time", "sex", paste0("prev", index$ages_all))
  names(model_prev_subset_female) <- c("time", "sex", paste0("prev", index$ages_all))
  names(model_prev_subset_male) <- c("time", "sex", paste0("prev", index$ages_all))

  # Combine sex-specific dataframes and turn into long format
  model_prev_subset <- rbind(model_prev_subset_female,
                             model_prev_subset_male,
                             model_prev_subset_both)

  model_prev_subset <- gather(model_prev_subset,
                              key = "age", value = "model_value", -time, -sex)
  model_prev_subset$age <- ages[as.numeric(gsub("\\D", "", model_prev_subset$age))]  # Assign ages as column

  # Merge with the dataset to fit to (transform factor to character vector)
  model_prev_subset$sex <- as.character(model_prev_subset$sex)

  mapped_output_seromarker <- full_join(prev_dataset,
                                        model_prev_subset,
                                        by = c("sex", "time", "age"))

  return(mapped_output_seromarker)

}

# Function to fit to age-, sex- and time-specific incidence rates from GLOBOCAN/GBD
# by mapping model output to year, age and sex
map_incidence_rates <- function(rate_outcome, rate_num, rate_denom, rate_timepoint, rate_dataset,
                                model_sim, model_out) {

  # Extract numerator dataset

  if (rate_outcome == "cirrhosis_mortality") {

    # Cirrhosis deaths are all HBV-related deaths minus HCC-related deaths
    # For women:
    num_f <- model_sim[,grepl("^cum_hbv_deathsf.",names(model_sim))] -
      model_sim[,grepl("^cum_hcc_deathsf.",names(model_sim))]
    # For men:
    num_m <- model_sim[,grepl("^cum_hbv_deathsm.",names(model_sim))] -
      model_sim[,grepl("^cum_hcc_deathsm.",names(model_sim))]

  } else {

    # For women:
    num_f <- select(model_sim, starts_with(paste0(rate_num,"f")))
    # For men:
    num_m <- select(model_sim, starts_with(paste0(rate_num,"m")))
  }

  # For women:
  rate_f <- subset(rate_dataset, outcome == rate_outcome & sex == "Female"
                   & time == rate_timepoint)[,c("outcome", "time", "sex", "age_min", "age_max")]
  rate_f$model_value <- 0
  rate_f$model_events <- 0

  # Extract denominator  dataset
  denom_f_label <- paste0(rate_denom, "_female")
  denom_f <- model_out[[denom_f_label]]

  for (i in 1:nrow(rate_f)){
    rate_f$model_value[i] <-
      (sum(num_f[which(model_sim$time == (rate_f$time[i]+1)), which(ages == rate_f$age_min[i]):which(ages == rate_f$age_max[i])])-
         sum(num_f[which(model_sim$time == rate_f$time[i]), which(ages == rate_f$age_min[i]):which(ages == rate_f$age_max[i])]))/
      sum(denom_f[which(model_sim$time == (rate_f$time[i]+0.5)),which(ages == rate_f$age_min[i]):which(ages == rate_f$age_max[i])])
  }

  for (i in 1:nrow(rate_f)){
    rate_f$model_events[i] <-
      (sum(num_f[which(model_sim$time == (rate_f$time[i]+1)), which(ages == rate_f$age_min[i]):which(ages == rate_f$age_max[i])])-
         sum(num_f[which(model_sim$time == rate_f$time[i]), which(ages == rate_f$age_min[i]):which(ages == rate_f$age_max[i])]))
    }

  # For men:
  denom_m_label <- paste0(rate_denom, "_male")
  rate_m <- subset(rate_dataset, outcome == rate_outcome & sex == "Male"
                   & time == rate_timepoint)[,c("outcome", "time", "sex", "age_min", "age_max")]
  rate_m$model_value <- 0
  rate_m$model_events <- 0

  # Extract denominator  dataset
  denom_m_label <- paste0(rate_denom, "_male")
  denom_m <- model_out[[denom_m_label]]

  for (i in 1:nrow(rate_m)){
    rate_m$model_value[i] <-
      (sum(num_m[which(model_sim$time == (rate_m$time[i]+1)),
                 which(ages == rate_m$age_min[i]):which(ages == rate_m$age_max[i])])-
         sum(num_m[which(model_sim$time == rate_m$time[i]),
                   which(ages == rate_m$age_min[i]):which(ages == rate_m$age_max[i])]))/
      sum(denom_m[which(model_sim$time == (rate_m$time[i]+0.5)),which(ages == rate_m$age_min[i]):which(ages == rate_f$age_max[i])])
  }

  for (i in 1:nrow(rate_m)){
    rate_m$model_events[i] <-
      (sum(num_m[which(model_sim$time == (rate_m$time[i]+1)),
                 which(ages == rate_m$age_min[i]):which(ages == rate_m$age_max[i])])-
         sum(num_m[which(model_sim$time == rate_m$time[i]),
                   which(ages == rate_m$age_min[i]):which(ages == rate_m$age_max[i])]))
  }


  # Combine HCC incidence and mortality sets and map to GLOBOCAN input data
  incidence_output <- rbind(rate_f, rate_m)

  return(incidence_output)

}

# Function to add noise to simulated output
# This is required for fitting deterministic models using ABC

# Distance function: metric which calculates the distance between the model output
# and the empirical data
# Takes as input a dataframe with the datapoints and the matching model output
calculate_distance <- function(mapped_output) {

  # The distance metric is the sum of weighted normalised absolute differences
  # between each datapoint and matching model output

  # Remove mapped output where no data point is available (missing data_value)
  mapped_output_for_error <- lapply(mapped_output, function(x) x[!is.na(x$data_value),])

  # Extract datapoints, their assigned quality weights and the matching model prediction
  # into separate vectors
  datapoints <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$data_value)))
  quality_weights <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$quality_weight)))
  model_prediction <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$model_value)))

  # Calculate vector if differences between each datapoint and simulated output
  data_model_diff <- datapoints-model_prediction  # observation - prediction

  # Take absolute value of difference, normalise it by dividing by the datapoint
  # multiply by the respective quality weight and sum across the vector of differences
  # For datapoints that are 0, replace with a value of 1 upon dividing.
  # This automatically gives a lower weight to these points though.
  error_term <- sum(quality_weights * (abs(data_model_diff)/
                                         replace(datapoints, datapoints==0, 1)))
  # if datapoint = 0, divide by 1. This automatically gives a lower weight to these points though.

  return(error_term)

}
# Calculate_distance_metrics function allows switching between metrics
calculate_distance_metrics <- function(mapped_output, metric) {

  # Remove mapped output where no data point is available (missing data_value)
  mapped_output_for_error <- lapply(mapped_output, function(x) x[!is.na(x$data_value),])

  # Extract datapoints, their assigned quality weights and the matching model prediction
  # into separate vectors
  datapoints <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$data_value)))
  quality_weights <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$quality_weight)))
  model_prediction <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$model_value)))

  # Calculate vector of differences between each datapoint and simulated output
  data_model_diff <- datapoints-model_prediction  # observation - prediction

  # Metric 1: Mean weighted relative difference
  # Mean of the weighted normalised absolute differences
  # between each datapoint and matching model output

  if(metric == "sum_rel_diff") {

    # Take absolute value of difference, normalise it by dividing by the datapoint
    # multiply by the respective quality weight and sum across the vector of differences
    # For datapoints that are 0, replace with a value of 1 upon dividing.
    # This automatically gives a lower weight to these points though.
    error_term <- sum(quality_weights *
                        (abs(data_model_diff)/replace(datapoints, datapoints==0, 1)))
    # if datapoint = 0, divide by 1. This automatically gives a lower weight to these points though.

  }

  else if(metric == "mean_rel_diff") {
    error_term <- sum(quality_weights * (abs(data_model_diff)/replace(datapoints, datapoints==0, 1)))/
                   sum(quality_weights)
  }

  else if(metric == "median_rel_diff") {
    error_term <- median(quality_weights * (abs(data_model_diff)/replace(datapoints, datapoints==0, 1)))
  }

  else if(metric == "euclidian") {
    error_term <- sqrt(sum((quality_weights * data_model_diff)^2))
  }

  else if(metric == "chisq") {
    error_term <- sum(quality_weights *
                        (data_model_diff^2/replace(datapoints, datapoints==0, 1)))
  }

  else if(metric == "rel_diff_vect") {
    error_term <- quality_weights *(abs(data_model_diff)/replace(datapoints, datapoints==0, 1))
  }

  else if(metric == "median_rel_squares") {
    error_term <- median(quality_weights *
                           (data_model_diff^2/replace(datapoints, datapoints==0, 1)))
  }

  else if(metric == "rel_euclidian") {
    error_term <- sqrt(sum((quality_weights*data_model_diff/replace(datapoints, datapoints==0, 1))^2))
  }

  return(error_term)

}

# Algorithm for random sampling and simulation with those parameter sets
run_on_cluster <- function(n_sims, data) {

  # 1) Draw parameter sets randomly from prior distribution
  n_sims <- n_sims  # number of simulations/parameter sets
  # First sample parameters where prior distributions depend on each other
  b1 <- runif(n_sims, 0.03, 0.7)
  b2 <- runif(n_sims, 0, b1)
  b3 <- runif(n_sims, 0, b1)
  mtct_prob_s <- rbeta(n_sims, 1.5,13.5)
  mtct_prob_e <- runif(n_sims, mtct_prob_s, 0.9)
  hccr_it <- rtruncnorm(n_sims, a=1, b=Inf, mean=6, sd=3) # normal truncated at 1
  hccr_cc <- runif(n_sims, hccr_it, 100)
  hccr_enchb <- runif(n_sims, hccr_it, hccr_cc)
  hccr_ir <- runif(n_sims, hccr_enchb, hccr_cc)
  # Combine in dataframe and sample remaining parameters
  params_mat <- data.frame(b1 = b1,
                           b2 = b2,
                           b3 = b3,
                           mtct_prob_s = mtct_prob_s,
                           mtct_prob_e = mtct_prob_e,
                           alpha = runif(n_sims, 1.5,10),
                           p_chronic_in_mtct = rbeta(n_sims, 10.49,1.3),
                           p_chronic_function_r = rnorm(n_sims,0.65,0.1),
                           p_chronic_function_s = rnorm(n_sims,0.46,0.1),
                           pr_it_ir = rgamma(n_sims,3.63,26.27),
                           pr_ir_ic = runif(n_sims, 0,1),
                           eag_prog_function_rate = runif(n_sims,0,0.01),
                           pr_ir_enchb = rgamma(n_sims, 1.49, 97.58),
                           pr_ir_cc_female = runif(n_sims, 0.005, 0.05),
                           pr_ir_cc_age_threshold = sample(0:15, n_sims,replace = TRUE),
                           pr_ic_enchb = rgamma(n_sims, 3.12, 141.30),
                           sag_loss_slope = rnorm(n_sims, 0.0004106, 0.00005),
                           pr_enchb_cc_female = rgamma(n_sims, 2.3, 123.8),
                           cirrhosis_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                           pr_cc_dcc = rgamma(n_sims,17.94,423.61),
                           cancer_prog_coefficient_female = runif(n_sims, 0.0001, 0.0003),
                           cancer_age_threshold = sample(0:15, n_sims,replace = TRUE),
                           cancer_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                           hccr_it = hccr_it,
                           hccr_ir = hccr_ir,
                           hccr_enchb = hccr_enchb,
                           hccr_cc = hccr_cc,
                           hccr_dcc = rgamma(n_sims, 8.09, 101.33),
                           mu_cc = rgamma(n_sims, 4.25, 124.91),
                           mu_dcc = rgamma(n_sims, 1.49, 0.98),
                           mu_hcc = rgamma(n_sims, 1.49, 0.98),
                           vacc_eff = rbeta(n_sims, 7.07, 0.37))

  # 2) Run n_sims simulations with these parameter sets and
  # calculate the distance function to compare to data
  out_mat <- apply(params_mat,1,
                   function(x)
                     fit_model(default_parameter_list = parameter_list,
                               data_to_fit = data,
                               parms_to_change =
                                 list(b1 = as.list(x)$b1,
                                      b2 = as.list(x)$b2,
                                      b3 = as.list(x)$b3,
                                      mtct_prob_s = as.list(x)$mtct_prob_s,
                                      mtct_prob_e = as.list(x)$mtct_prob_e,
                                      alpha = as.list(x)$alpha,
                                      p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                      p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                      p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                      pr_it_ir = as.list(x)$pr_it_ir,
                                      pr_ir_ic = as.list(x)$pr_ir_ic,
                                      eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                      pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                      pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                      pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                      pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                      sag_loss_slope = as.list(x)$sag_loss_slope,
                                      pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                      cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                      pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                      cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                      cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                      cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                      hccr_it = as.list(x)$hccr_it,
                                      hccr_ir = as.list(x)$hccr_ir,
                                      hccr_enchb = as.list(x)$hccr_enchb,
                                      hccr_cc = as.list(x)$hccr_cc,
                                      hccr_dcc = as.list(x)$hccr_dcc,
                                      mu_cc = as.list(x)$mu_cc,
                                      mu_dcc = as.list(x)$mu_dcc,
                                      mu_hcc = as.list(x)$mu_hcc,
                                      vacc_eff = as.list(x)$vacc_eff
                                 )))

  # 3) Return matrix of parameter sets and matching error term
  #out_mat_subset <- sapply(out_mat, "[[", "error_term")
  #res_mat <- cbind(params_mat, error_term = out_mat_subset)

  return(out_mat)

}

# This function uses my first choice of priors
run_on_cluster_parallel <- function(n_sims, data) {

  # 1) Draw parameter sets randomly from prior distribution
  n_sims <- n_sims  # number of simulations/parameter sets
  # First sample parameters where prior distributions depend on each other
  b1 <- runif(n_sims, 0.03, 0.7)
  b2 <- runif(n_sims, 0, b1)
  b3 <- runif(n_sims, 0, b1)
  mtct_prob_s <- rbeta(n_sims, 1.5,13.5)
  mtct_prob_e <- runif(n_sims, mtct_prob_s, 0.9)
  hccr_it <- rtruncnorm(n_sims, a=1, b=Inf, mean=6, sd=3) # normal truncated at 1
  hccr_cc <- runif(n_sims, hccr_it, 100)
  hccr_enchb <- runif(n_sims, hccr_it, hccr_cc)
  hccr_ir <- runif(n_sims, hccr_enchb, hccr_cc)
  # Combine in dataframe and sample remaining parameters
  params_mat <- data.frame(b1 = b1,
                           b2 = b2,
                           b3 = b3,
                           mtct_prob_s = mtct_prob_s,
                           mtct_prob_e = mtct_prob_e,
                           alpha = runif(n_sims, 1.5,10),
                           p_chronic_in_mtct = rbeta(n_sims, 10.49,1.3),
                           p_chronic_function_r = rnorm(n_sims,0.65,0.1),
                           p_chronic_function_s = rnorm(n_sims,0.46,0.1),
                           pr_it_ir = rgamma(n_sims,3.63,26.27),
                           pr_ir_ic = runif(n_sims, 0,1),
                           eag_prog_function_rate = runif(n_sims,0,0.01),
                           pr_ir_enchb = rgamma(n_sims, 1.49, 97.58),
                           pr_ir_cc_female = runif(n_sims, 0.005, 0.05),
                           pr_ir_cc_age_threshold = sample(0:15,n_sims,replace=TRUE),
                           pr_ic_enchb = rgamma(n_sims, 3.12, 141.30),
                           sag_loss_slope = rnorm(n_sims, 0.0004106, 0.00005),
                           pr_enchb_cc_female = rgamma(n_sims, 2.3, 123.8),
                           cirrhosis_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                           pr_cc_dcc = rgamma(n_sims,17.94,423.61),
                           cancer_prog_coefficient_female = runif(n_sims, 0.0001, 0.0003),
                           cancer_age_threshold = sample(0:15,n_sims,replace=TRUE),
                           cancer_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                           hccr_it = hccr_it,
                           hccr_ir = hccr_ir,
                           hccr_enchb = hccr_enchb,
                           hccr_cc = hccr_cc,
                           hccr_dcc = rgamma(n_sims, 8.09, 101.33),
                           mu_cc = rgamma(n_sims, 4.25, 124.91),
                           mu_dcc = rgamma(n_sims, 1.49, 0.98),
                           mu_hcc = rgamma(n_sims, 1.49, 0.98),
                           vacc_eff = rbeta(n_sims, 7.07, 0.37))

  # 2) Run n_sims simulations with these parameter sets and
  # calculate the distance function to compare to data
  out_mat <- parApply(cl = NULL, params_mat, 1,
                      function(x) fit_model(default_parameter_list = parameter_list,
                                            data_to_fit = data,
                                            parms_to_change =
                                              list(b1 = as.list(x)$b1,
                                                   b2 = as.list(x)$b2,
                                                   b3 = as.list(x)$b3,
                                                   mtct_prob_s = as.list(x)$mtct_prob_s,
                                                   mtct_prob_e = as.list(x)$mtct_prob_e,
                                                   alpha = as.list(x)$alpha,
                                                   p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                                   p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                                   p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                                   pr_it_ir = as.list(x)$pr_it_ir,
                                                   pr_ir_ic = as.list(x)$pr_ir_ic,
                                                   eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                                   pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                                   pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                                   pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                                   pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                                   sag_loss_slope = as.list(x)$sag_loss_slope,
                                                   pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                                   cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                                   pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                                   cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                                   cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                                   cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                                   hccr_it = as.list(x)$hccr_it,
                                                   hccr_ir = as.list(x)$hccr_ir,
                                                   hccr_enchb = as.list(x)$hccr_enchb,
                                                   hccr_cc = as.list(x)$hccr_cc,
                                                   hccr_dcc = as.list(x)$hccr_dcc,
                                                   mu_cc = as.list(x)$mu_cc,
                                                   mu_dcc = as.list(x)$mu_dcc,
                                                   mu_hcc = as.list(x)$mu_hcc,
                                                   vacc_eff = as.list(x)$vacc_eff)))

  # 3) Return matrix of parameter sets and matching error term
  #out_mat_subset <- sapply(out_mat, "[[", "error_term")
  #res_mat <- cbind(params_mat, error_term = out_mat_subset)

  # TEST: only return full output if median rel diff is less than 0.5
  #best_fit_ids <- which(sapply(out_mat, "[[", "error_term") < 0.5)
  #res <- list(res_mat = res_mat,
  #            out_mat = out_mat[best_fit_ids])

  #return(res_mat)

  #return(res)

  return(out_mat)

}

# This function uses my updated choice of priors - less informative and more in line with African data
run_on_cluster_parallel_new_priors <- function(n_sims, data) {

  # 1) Draw parameter sets randomly from prior distribution
  n_sims <- n_sims  # number of simulations/parameter sets
  # First sample parameters where prior distributions depend on each other
  b1 <- runif(n_sims, 0.03, 0.7)
  b2 <- runif(n_sims, 0, b1)
  b3 <- runif(n_sims, 0, b1)
  mtct_prob_s <- rbeta(n_sims, 1.5,13.5)
  mtct_prob_e <- runif(n_sims, mtct_prob_s, 0.9)
  hccr_it <- rtruncnorm(n_sims, a=1, b=Inf, mean=6, sd=3) # normal truncated at 1
  hccr_cc <- runif(n_sims, hccr_it, 100)
  hccr_enchb <- runif(n_sims, hccr_it, hccr_cc)
  hccr_ir <- runif(n_sims, hccr_enchb, hccr_cc)
  # Combine in dataframe and sample remaining parameters
  params_mat <- data.frame(b1 = b1,
                           b2 = b2,
                           b3 = b3,
                           mtct_prob_s = mtct_prob_s,
                           mtct_prob_e = mtct_prob_e,
                           alpha = runif(n_sims, 1.5,10),
                           p_chronic_in_mtct = rbeta(n_sims, 10.49,1.3),
                           p_chronic_function_r = rnorm(n_sims,0.65,0.1),
                           p_chronic_function_s = rnorm(n_sims,0.46,0.1),
                           pr_it_ir = rgamma(n_sims,3.63,26.27),
                           pr_ir_ic = runif(n_sims, 0,1),
                           eag_prog_function_rate = runif(n_sims,0,0.01),
                           pr_ir_enchb = rgamma(n_sims, 1.22, 44.20),
                           pr_ir_cc_female = runif(n_sims, 0.005, 0.05),
                           pr_ir_cc_age_threshold = sample(0:15,n_sims,replace=TRUE),
                           pr_ic_enchb = rgamma(n_sims, 3.49, 83.05),
                           sag_loss_slope = rnorm(n_sims, 0.0004106, 0.00005),
                           pr_enchb_cc_female = rgamma(n_sims, 1.23, 22.33),
                           cirrhosis_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                           pr_cc_dcc = rgamma(n_sims,17.94,423.61),
                           cancer_prog_coefficient_female = runif(n_sims, 0.0001, 0.0003),
                           cancer_age_threshold = sample(0:15,n_sims,replace=TRUE),
                           cancer_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                           hccr_it = hccr_it,
                           hccr_ir = hccr_ir,
                           hccr_enchb = hccr_enchb,
                           hccr_cc = hccr_cc,
                           hccr_dcc = rgamma(n_sims, 3.08, 29.76),
                           mu_cc = rgamma(n_sims, 4.25, 124.91),
                           mu_dcc = rgamma(n_sims, 2.18, 1.18),
                           mu_hcc = rgamma(n_sims, 2.18, 1.18),
                           vacc_eff = rbeta(n_sims, 7.07, 0.37))

  # 2) Run n_sims simulations with these parameter sets and
  # calculate the distance function to compare to data
  out_mat <- parApply(cl = NULL, params_mat, 1,
                      function(x) fit_model(default_parameter_list = parameter_list,
                                            data_to_fit = data,
                                            parms_to_change =
                                              list(b1 = as.list(x)$b1,
                                                   b2 = as.list(x)$b2,
                                                   b3 = as.list(x)$b3,
                                                   mtct_prob_s = as.list(x)$mtct_prob_s,
                                                   mtct_prob_e = as.list(x)$mtct_prob_e,
                                                   alpha = as.list(x)$alpha,
                                                   p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                                   p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                                   p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                                   pr_it_ir = as.list(x)$pr_it_ir,
                                                   pr_ir_ic = as.list(x)$pr_ir_ic,
                                                   eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                                   pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                                   pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                                   pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                                   pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                                   sag_loss_slope = as.list(x)$sag_loss_slope,
                                                   pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                                   cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                                   pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                                   cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                                   cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                                   cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                                   hccr_it = as.list(x)$hccr_it,
                                                   hccr_ir = as.list(x)$hccr_ir,
                                                   hccr_enchb = as.list(x)$hccr_enchb,
                                                   hccr_cc = as.list(x)$hccr_cc,
                                                   hccr_dcc = as.list(x)$hccr_dcc,
                                                   mu_cc = as.list(x)$mu_cc,
                                                   mu_dcc = as.list(x)$mu_dcc,
                                                   mu_hcc = as.list(x)$mu_hcc,
                                                   vacc_eff = as.list(x)$vacc_eff)))

  # 3) Return matrix of parameter sets and matching error term
  #out_mat_subset <- sapply(out_mat, "[[", "error_term")
  #res_mat <- cbind(params_mat, error_term = out_mat_subset)

  # TEST: only return full output if median rel diff is less than 0.5
  #best_fit_ids <- which(sapply(out_mat, "[[", "error_term") < 0.5)
  #res <- list(res_mat = res_mat,
  #            out_mat = out_mat[best_fit_ids])

  #return(res_mat)

  #return(res)

  return(out_mat)

}


### Set up calibration ----

# For fitting, the model is run from 1880 to allow 100 years (1 generation)
# before the first datapoints

# Load initial population saved from previous model run
load(here("data/simulated_inits_1880.RData"))  # this is saved from previous model run
init_pop_sim <- c("Sf" = select(model_pop1880, starts_with("Sf")),
                  "ITf" = select(model_pop1880, starts_with("ITf")),
                  "IRf" = select(model_pop1880, starts_with("IRf")),
                  "ICf" = select(model_pop1880, starts_with("ICf")),
                  "ENCHBf" = select(model_pop1880, starts_with("ENCHBf")),
                  "CCf" = select(model_pop1880, starts_with("CCf")),
                  "DCCf" = select(model_pop1880, starts_with("DCCf")),
                  "HCCf" = select(model_pop1880, starts_with("HCCf")),
                  "Rf" = select(model_pop1880, starts_with("Rf")),
                  "Sm" = select(model_pop1880, starts_with("Sm")),
                  "ITm" = select(model_pop1880, starts_with("ITm")),
                  "IRm" = select(model_pop1880, starts_with("IRm")),
                  "ICm" = select(model_pop1880, starts_with("ICm")),
                  "ENCHBm" = select(model_pop1880, starts_with("ENCHBm")),
                  "CCm" = select(model_pop1880, starts_with("CCm")),
                  "DCCm" = select(model_pop1880, starts_with("DCCm")),
                  "HCCm" = select(model_pop1880, starts_with("HCCm")),
                  "Rm" = select(model_pop1880, starts_with("Rm")),
                  output_storage)
init_pop_sim <- unlist(init_pop_sim)

