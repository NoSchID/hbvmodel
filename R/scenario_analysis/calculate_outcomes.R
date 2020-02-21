# Functions to calculate analysis outcomes

# Calculate healthcare interactions (note IT interactions are recorded separately)
# Function can be applied to run_model_output (single simulation)
# Need to add record of vaccination of susceptibles!
calculate_healthcare_interactions <- function(sim, from_year_vector, by_year_vector, scenario_label) {

  # Timing of screening programme
  time_of_screening <- sim$input_parameters$screening_years

  # Number of HBsAg tests
  total_sag_tests <- unique(sim$out$total_screened_susceptible +
                              sim$out$total_screened_immune +
                              sim$out$total_screened_it +
                              sim$out$total_screened_chb +
                              sim$out$total_screened_cirrhosis +
                              sim$out$total_screened_ineligible)

  # Number of liver disease assessments if IT is NOT treated
  if(sim$input_parameters$apply_treat_it == 0) {

    total_identified_as_ineligible <- unique(sim$out$total_screened_ineligible + sim$out$total_screened_it) *
      sim$input_parameters$link_to_care_prob

    total_identified_as_eligible <-
      unique(sim$out$total_screened_chb +
               sim$out$total_screened_cirrhosis) * sim$input_parameters$link_to_care_prob

    # Number of liver disease assessments if IT is treated

  } else if(sim$input_parameters$apply_treat_it == 1) {
    total_identified_as_ineligible <- unique(sim$out$total_screened_ineligible) *
      sim$input_parameters$link_to_care_prob

    total_identified_as_eligible <-
      unique(sim$out$total_screened_it + sim$out$total_screened_chb +
               sim$out$total_screened_cirrhosis) * sim$input_parameters$link_to_care_prob

  }


  # Number of people immediately starting treatment
  total_immediate_treatment_initiations <- total_identified_as_eligible *
    sim$input_parameters$treatment_initiation_prob

  # Combine into dataframe
  res1 <- data.frame(screening_years = time_of_screening,
                    total_screened = total_sag_tests[total_sag_tests != 0],
                    total_identified_as_ineligible = total_identified_as_ineligible[total_identified_as_ineligible != 0],
                    total_identified_as_eligible = total_identified_as_eligible[total_identified_as_eligible != 0],
                    total_immediate_treatment_initiations = total_immediate_treatment_initiations[total_immediate_treatment_initiations != 0],
                    scenario = scenario_label)

  # Monitoring and treatment

  # Number of monitoring events after screening by compartment
  # Need to change ineligible depending on whether IT is included
  monitored_ineligible <- sim$out[,grepl("^cum_monitored_icf", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_icm", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_hccf", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_hccm", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_rf", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_rm", colnames(sim$out))]

  monitored_it <- sim$out[,grepl("^cum_monitored_itf", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_itm", colnames(sim$out))]
  monitored_ir <- sim$out[,grepl("^cum_monitored_irf", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_irm", colnames(sim$out))]
  monitored_enchb <- sim$out[,grepl("^cum_monitored_enchbf", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_enchbm", colnames(sim$out))]
  monitored_cc <- sim$out[,grepl("^cum_monitored_ccf", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_ccm", colnames(sim$out))]
  monitored_dcc <- sim$out[,grepl("^cum_monitored_dccf", colnames(sim$out))] +
    sim$out[,grepl("^cum_monitored_dccm", colnames(sim$out))]

  # Number of treatment initiations as a result of monitoring
  # Sum those monitoring events in treatment eligible compartments and multiply by treatment initiation probability
  cum_monitoring_treatment_initiations <-
    ((apply(monitored_ir[which(sim$out$time %in% by_year_vector),],1,sum)-
    apply(monitored_ir[which(sim$out$time %in% from_year_vector),],1,sum)) +
    (apply(monitored_enchb[which(sim$out$time %in% by_year_vector),],1,sum)-
       apply(monitored_enchb[which(sim$out$time %in% from_year_vector),],1,sum)) +
    (apply(monitored_cc[which(sim$out$time %in% by_year_vector),],1,sum)-
       apply(monitored_cc[which(sim$out$time %in% from_year_vector),],1,sum)) +
    (apply(monitored_dcc[which(sim$out$time %in% by_year_vector),],1,sum)-
       apply(monitored_dcc[which(sim$out$time %in% from_year_vector),],1,sum))) *
    sim$input_parameters$treatment_initiation_prob

  if (sim$input_parameters$apply_treat_it == 1) {
    cum_monitoring_treatment_initiations_it <-
      (apply(monitored_it[which(sim$out$time %in% by_year_vector),],1,sum)-
         apply(monitored_it[which(sim$out$time %in% from_year_vector),],1,sum)) *
      sim$input_parameters$treatment_initiation_prob
  } else if (sim$input_parameters$apply_treat_it == 0) {
    cum_monitoring_treatment_initiations_it <- 0
  }


  # Combine into dataframe
  res2 <- data.frame(from_year = from_year_vector,
                     by_year = by_year_vector,
                     cum_monitoring_events_it =
                       apply(monitored_it[which(sim$out$time %in% by_year_vector),],1,sum)-
                       apply(monitored_it[which(sim$out$time %in% from_year_vector),],1,sum),
                     cum_monitoring_events_ir =
                       apply(monitored_ir[which(sim$out$time %in% by_year_vector),],1,sum)-
                       apply(monitored_ir[which(sim$out$time %in% from_year_vector),],1,sum),
                     cum_monitoring_events_enchb =
                       apply(monitored_enchb[which(sim$out$time %in% by_year_vector),],1,sum)-
                       apply(monitored_enchb[which(sim$out$time %in% from_year_vector),],1,sum),
                     cum_monitoring_events_cc =
                       apply(monitored_cc[which(sim$out$time %in% by_year_vector),],1,sum)-
                       apply(monitored_cc[which(sim$out$time %in% from_year_vector),],1,sum),
                     cum_monitoring_events_dcc =
                       apply(monitored_dcc[which(sim$out$time %in% by_year_vector),],1,sum)-
                       apply(monitored_dcc[which(sim$out$time %in% from_year_vector),],1,sum),
                     cum_monitoring_events_ineligible =
                       apply(monitored_ineligible[which(sim$out$time %in% by_year_vector),],1,sum)-
                       apply(monitored_ineligible[which(sim$out$time %in% from_year_vector),],1,sum),
                     cum_monitoring_treatment_initiations,
                     cum_monitoring_treatment_initiations_it,
                     scenario = scenario_label)

  res3 <- ifelse(sim$input_parameters$apply_treat_it == 1, "IT >30 treated", "IT not treated")

  return(list(res1,res2, res3))

}

calculate_healthcare_interactions(sim[[1]], from_year_vector = 2020, by_year_vector = c(2030,2050),
                                  scenario_label = "test")

# Extract time series
# Function can be applied to run_model_output (single simulation)
extract_time_series <- function(output_file, scenario_label) {

  # Note: Calculating incidence and rates with the correct denominator:
  # annual incidence at time t / population at risk at time t-1
  # On a vector, achieve this by:
  # - removing the first element (tail(x,-1)) to get incidence IN (not by) the given timestep
  # - removing the last element (head(x,-1)) of the denominator and/or time

  # HBsAg prevalence
  total_prev <- data.frame(number_infected = head(output_file$infectioncat_total$carriers,-1),
                           prev = head(output_file$infectioncat_total$carriers,-1)/
                             head(output_file$pop_total$pop_total,-1))

  # New cases of chronic HBV carriage per timestep and per population
  total_chronic_incidence <- data.frame(total_chronic_infections = tail(output_file$incident_chronic_infections$horizontal_chronic_infections+
                                          output_file$incident_chronic_infections$chronic_births,-1),
                                        chronic_births = tail(output_file$incident_chronic_infections$chronic_births,-1))
  total_chronic_incidence$total_chronic_infections_rate <-
    total_chronic_incidence$total_chronic_infections/head(output_file$pop_total$pop_total,-1)

  # HBV-related deaths per timestep and per population
  total_hbv_deaths <- data.frame(total_hbv_deaths = tail(output_file$hbv_deaths$incident_number_total+
                                   output_file$screened_hbv_deaths$incident_number_total+
                                   output_file$treated_hbv_deaths$incident_number_total, -1),
                                 hbv_deaths_male = tail(output_file$hbv_deaths$incident_number_male+
                                   output_file$screened_hbv_deaths$incident_number_male+
                                   output_file$treated_hbv_deaths$incident_number_male, -1))
  total_hbv_deaths$total_hbv_deaths_rate <- total_hbv_deaths$total_hbv_deaths/
                                                head(output_file$pop_total$pop_total,-1)
  total_hbv_deaths$hbv_deaths_rate_male <-  total_hbv_deaths$hbv_deaths_male/
    head(output_file$pop_total$pop_male,-1)

  outcome_df <- cbind(time = head(output_file$time, -1),
                      scenario = scenario_label,
                      total_prev, total_chronic_incidence, total_hbv_deaths)


  #return(list(total_prev = total_prev,
  #            total_chronic_incidence = total_chronic_incidence,
  #            total_hbv_deaths = total_hbv_deaths))

  return(outcome_df)

}
# To add: age-standardised HBV related deaths, HCC incidence maybe

# Generic function to extract any 1 outcome from model output
extract_time_series_generic <- function(output_file, type, numerator1, numerator2, denominator1 = NULL, denominator2 = NULL,
                                        scenario_label) {

  # Note: Calculating incidence and rates with the correct denominator:
  # annual incidence at time t / population at risk at time t-1
  # On a vector, achieve this by:
  # - removing the first element (tail(x,-1)) to get incidence IN (not by) the given timestep
  # - removing the last element (head(x,-1)) of the denominator and/or time

  if (type == "prevalence" & is.null(denominator1)) {
    outcome <- data.frame(time = output_file$time,
                          outcome = output_file[[numerator1]][numerator2])
  } else  if (type == "prevalence" & !(is.null(denominator1))) {
    outcome <- data.frame(time = output_file$time,
                          outcome = output_file[[numerator1]][numerator2]/
                            output_file[[denominator1]][denominator2])
  } else if (type == "incidence" & is.null(denominator1)) {
    outcome <- data.frame(time = head(output_file$time,-1),
                          outcome = tail(output_file[[numerator1]][numerator2],-1))
  } else if (type == "incidence" & !(is.null(denominator1))) {
    outcome <- data.frame(time = head(output_file$time,-1),
                          outcome = tail(output_file[[numerator1]][numerator2],-1)/
                            head(output_file[[denominator1]][denominator2],-1))
  }

  return(outcome)

}

# This function applies extract_time_series to all sets and returns the median, 2.5th and 97.5th percentile
summarise_time_series <- function(output_file, scenario_label) {



}


test <- lapply(out,extract_time_series, "test")

x <- list()

for (i in 3:11) {
  x[[i-2]] <- as.data.frame(sapply(test, "[[", i))
}

# Need to name each list element (outcome), then summarise

# Add function for primary outcomes

