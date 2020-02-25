# Functions to calculate analysis outcomes

# These are all applied to the output from code_model_output
# If the argument is output_file, it is applied to 1 parmset,
# if if the argument is output_files, it is applied to a list of parmsets

# Calculate healthcare interactions (note IT interactions are recorded separately):
# 2 separate functions for the interactions that occur on screening vs. monitoring
calculate_screening_interactions <- function(output_file, scenario_label) {

  # Timing of screening programme
  time_of_screening <- output_file$input_parameters$screening_years

  # Number of HBsAg tests
  total_sag_tests <- unique(output_file$full_output$total_screened_susceptible +
                              output_file$full_output$total_screened_immune +
                              output_file$full_output$total_screened_it +
                              output_file$full_output$total_screened_chb +
                              output_file$full_output$total_screened_cirrhosis +
                              output_file$full_output$total_screened_ineligible)

  # Number of liver disease assessments if IT is NOT treated
  if(output_file$input_parameters$apply_treat_it == 0) {

    total_identified_as_ineligible <- unique(output_file$full_output$total_screened_ineligible +
                                               output_file$full_output$total_screened_it) *
      output_file$input_parameters$link_to_care_prob

    total_identified_as_eligible <-
      unique(output_file$full_output$total_screened_chb +
               output_file$full_output$total_screened_cirrhosis) * output_file$input_parameters$link_to_care_prob

    # Number of liver disease assessments if IT is treated

  } else if(output_file$input_parameters$apply_treat_it == 1) {
    total_identified_as_ineligible <- unique(output_file$full_output$total_screened_ineligible) *
      output_file$input_parameters$link_to_care_prob

    total_identified_as_eligible <-
      unique(output_file$full_output$total_screened_it + output_file$full_output$total_screened_chb +
               output_file$full_output$total_screened_cirrhosis) * output_file$input_parameters$link_to_care_prob

  }


  # Number of people immediately starting treatment
  total_immediate_treatment_initiations <- total_identified_as_eligible *
    output_file$input_parameters$treatment_initiation_prob

  # Combine into dataframe
  res1 <- data.frame(scenario = scenario_label,
                     screening_years = time_of_screening,
                     total_screened = total_sag_tests[total_sag_tests != 0],
                     total_identified_as_ineligible = total_identified_as_ineligible[total_identified_as_ineligible != 0],
                     total_identified_as_eligible = total_identified_as_eligible[total_identified_as_eligible != 0],
                     total_immediate_treatment_initiations = total_immediate_treatment_initiations[total_immediate_treatment_initiations != 0])


  return(res1)

}

# For monitoring, need to add record of vaccination of susceptibles!
calculate_monitoring_interactions <- function(output_file, from_year, by_year, scenario_label) {

  # Monitoring and treatment

  # Number of monitoring events after screening by compartment
  # Need to change ineligible depending on whether IT is included
  monitored_ineligible <- output_file$full_output[,grepl("^cum_monitored_icf", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_icm", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_hccf", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_hccm", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_rf", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_rm", colnames(output_file$full_output))]

  monitored_it <- output_file$full_output[,grepl("^cum_monitored_itf", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_itm", colnames(output_file$full_output))]
  monitored_ir <- output_file$full_output[,grepl("^cum_monitored_irf", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_irm", colnames(output_file$full_output))]
  monitored_enchb <- output_file$full_output[,grepl("^cum_monitored_enchbf", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_enchbm", colnames(output_file$full_output))]
  monitored_cc <- output_file$full_output[,grepl("^cum_monitored_ccf", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_ccm", colnames(output_file$full_output))]
  monitored_dcc <- output_file$full_output[,grepl("^cum_monitored_dccf", colnames(output_file$full_output))] +
    output_file$full_output[,grepl("^cum_monitored_dccm", colnames(output_file$full_output))]

  # Number of treatment initiations as a result of monitoring
  # Sum those monitoring events in treatment eligible compartments and multiply by treatment initiation probability
  cum_monitoring_treatment_initiations <-
    ((apply(monitored_ir[which(output_file$time == by_year),],1,sum)-
        apply(monitored_ir[which(output_file$time == from_year),],1,sum)) +
       (apply(monitored_enchb[which(output_file$time == by_year),],1,sum)-
          apply(monitored_enchb[which(output_file$time == from_year),],1,sum)) +
       (apply(monitored_cc[which(output_file$time == by_year),],1,sum)-
          apply(monitored_cc[which(output_file$time == from_year),],1,sum)) +
       (apply(monitored_dcc[which(output_file$time == by_year),],1,sum)-
          apply(monitored_dcc[which(output_file$time == from_year),],1,sum))) *
    output_file$input_parameters$treatment_initiation_prob

  if (output_file$input_parameters$apply_treat_it == 1) {
    cum_monitoring_treatment_initiations_it <-
      (apply(monitored_it[which(output_file$time == by_year),],1,sum)-
         apply(monitored_it[which(output_file$time == from_year),],1,sum)) *
      output_file$input_parameters$treatment_initiation_prob
  } else if (output_file$input_parameters$apply_treat_it == 0) {
    cum_monitoring_treatment_initiations_it <- 0
  }


  # Combine into dataframe
  res <- data.frame(scenario = scenario_label,
                    from_year = from_year,
                     by_year = by_year,
                     cum_monitoring_events_it =
                       apply(monitored_it[which(output_file$time == by_year),],1,sum)-
                       apply(monitored_it[which(output_file$time == from_year),],1,sum),
                     cum_monitoring_events_ir =
                       apply(monitored_ir[which(output_file$time == by_year),],1,sum)-
                       apply(monitored_ir[which(output_file$time == from_year),],1,sum),
                     cum_monitoring_events_enchb =
                       apply(monitored_enchb[which(output_file$time == by_year),],1,sum)-
                       apply(monitored_enchb[which(output_file$time == from_year),],1,sum),
                     cum_monitoring_events_cc =
                       apply(monitored_cc[which(output_file$time == by_year),],1,sum)-
                       apply(monitored_cc[which(output_file$time == from_year),],1,sum),
                     cum_monitoring_events_dcc =
                       apply(monitored_dcc[which(output_file$time == by_year),],1,sum)-
                       apply(monitored_dcc[which(output_file$time == from_year),],1,sum),
                     cum_monitoring_events_ineligible =
                       apply(monitored_ineligible[which(output_file$time == by_year),],1,sum)-
                       apply(monitored_ineligible[which(output_file$time == from_year),],1,sum),
                     cum_monitoring_treatment_initiations,
                     cum_monitoring_treatment_initiations_it)

  return(res)

}

# Function to summarise healthcare interactions
# Total number of healthcare interactions = HBsAg tests+liver assessments+treatment initiations,
# initially + through monitoring
summarise_healthcare_interactions <- function(output_files, from_year, by_year, scenario_label) {
  # Immediate interactions upon screening
  screening_interactions <- lapply(output_files, calculate_screening_interactions, scenario_label)

  total_screened <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                               sapply(screening_interactions, "[[", "total_screened"))

  total_assessed_immediately <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                                           sapply(screening_interactions, "[[", "total_identified_as_ineligible")+
                                             sapply(screening_interactions, "[[", "total_identified_as_eligible"))

  total_treated_immediately <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                                          sapply(screening_interactions, "[[", "total_immediate_treatment_initiations"))

  # Add interactions together for all years >= from_year and < by_year
  total_screened <- apply(total_screened[total_screened$screening_years>=from_year &
                                           total_screened$screening_years<by_year,-1],2,sum)

  total_assessed_immediately <- apply(total_assessed_immediately[
    total_assessed_immediately$screening_years>=from_year &
      total_assessed_immediately$screening_years<by_year,-1],2,sum)

  total_treated_immediately <- apply(total_treated_immediately[
    total_treated_immediately$screening_years >= from_year &
      total_treated_immediately$screening_years<by_year,-1],2,sum)

  # Interactions during and after monitoring
  monitoring_interactions <- data.frame(sapply(output_files, calculate_monitoring_interactions, from_year, by_year, scenario_label))
  monitoring_interactions <- as.data.frame(apply(monitoring_interactions,2,unlist))

  total_monitored <- monitoring_interactions["cum_monitoring_events_it",]+
    monitoring_interactions["cum_monitoring_events_ir",]+
    monitoring_interactions["cum_monitoring_events_enchb",]+
    monitoring_interactions["cum_monitoring_events_cc",]+
    monitoring_interactions["cum_monitoring_events_dcc",]+
    monitoring_interactions["cum_monitoring_events_ineligible",]   # ignore rowname

  total_treated_after_monitoring <- monitoring_interactions["cum_monitoring_treatment_initiations",]+
    monitoring_interactions["cum_monitoring_treatment_initiations_it",]

  total_assessed <- total_assessed_immediately+total_monitored
  rownames(total_assessed) <- NULL
  total_treated <- total_treated_immediately+total_treated_after_monitoring
  rownames(total_treated) <- NULL


  total_screened_res <- cbind(data.frame(from_year = from_year,
                                     by_year = by_year,
                                     scenario = scenario_label),
                          t(as.data.frame(total_screened)))
  total_assessed_res <- cbind(data.frame(from_year = from_year,
                                         by_year = by_year,
                                         scenario = scenario_label),
                          as.data.frame(total_assessed))
  total_treated_res <- cbind(data.frame(from_year = from_year,
                                        by_year = by_year,
                                        scenario = scenario_label),
                         as.data.frame(total_treated))


  total_interactions_res <- cbind(data.frame(from_year = from_year,
                                             by_year = by_year,
                                             scenario = scenario_label),
                                  t(as.data.frame(total_screened))+total_assessed+total_treated)


  res <- list(total_interactions = total_interactions_res,
              total_screened = total_screened_res,
              total_assessed = total_assessed_res,
              total_treated = total_treated_res)

  return(res)


}

# Outcomes: total screened (HBsAg)
# Total liver disease assessments
# Total treatment initiations

# Function to calculate age-standardised rate of HBV-related deaths
calculate_age_standardised_hbv_deaths_rate <- function(output_file) {
  # Age-standardised incidence of HBV-related deaths per 100000 per timestep

  # a) Calculate crude age-specific rates per person-year: need age-specific number of deaths and age-specific population size
 deaths_by_age <- output_file$full_output[,grepl("^cum_hbv_deathsf.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_hbv_deathsm.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_screened_hbv_deathsf.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_screened_hbv_deathsm.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_treated_hbv_deathsf.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_treated_hbv_deathsm.",names(output_file$full_output))]

  deaths_by_age <- calculate_incident_numbers(deaths_by_age)

  pop_by_age <- output_file$pop

  # Crude rate
  deaths_rate <- deaths_by_age[-1,]/pop_by_age[-nrow(pop_by_age),]
  deaths_rate[is.na(deaths_rate)|deaths_rate == Inf] <- 0   # where pop was 0

  # b) Multiply by reference population (Gambian pop in 2020)
  ref_pop <- pop_by_age[output_file$time==2020,]
  deaths_rate_ref_pop <- sweep(as.matrix(deaths_rate), MARGIN=2, as.matrix(ref_pop), "*")

  # c) Calculate total expected deaths (sum of age-specific numbers) and divide by total Gambian popsize in 2020
  age_standardised_rate <- rowSums(deaths_rate_ref_pop)/rowSums(ref_pop)

  return(age_standardised_rate)
}

# Extract time series
# Function can be applied to run_model_output+code_model_output (single simulation)
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
  total_hbv_deaths$total_hbv_deaths_age_standardised_rate <-
    calculate_age_standardised_hbv_deaths_rate(output_file)

  outcome_df <- cbind(time = head(output_file$time, -1),
                      scenario = scenario_label,
                      total_prev, total_chronic_incidence, total_hbv_deaths)


  #return(list(total_prev = total_prev,
  #            total_chronic_incidence = total_chronic_incidence,
  #            total_hbv_deaths = total_hbv_deaths))

  return(outcome_df)

}
# To add: age-standardised HBV related deaths, HCC incidence maybe

# This function applies extract_time_series to all sets and returns the median, 2.5th and 97.5th percentile
# Calls extract_time_series function
summarise_time_series <- function(output_files, scenario_label, summarise_percentiles = TRUE) {

  # Extract timeseries outcomes for each simulation individually
  timeseries <- lapply(output_files, extract_time_series, scenario_label)

  time <- lapply(timeseries, "[[", "time")
  scenario <- lapply(timeseries, "[[", "scenario")

  outcome_list <- list()

  # Except for first 2 columns (time and scenario), extract each outcome into a separate list element
  for (i in 3:ncol(timeseries[[1]])) {
    outcome_list[[i-2]] <- as.data.frame(sapply(timeseries, "[[", i))
  }

  outcome_list <- Map(cbind,time, scenario, outcome_list)
  outcome_list <- lapply(outcome_list, setNames, c("time", "scenario", colnames(outcome_list[[1]])[-c(1,2)]))

  # Rename list to outcome
  names(outcome_list) <- colnames(timeseries[[1]])[3:ncol(timeseries[[1]])]
  # outcome_list contains each outcome as a list element, with individual simulations being the columns
  # in each dataframe

  # Take median, 2.5th and 97.th percentile of all simulations
  if (summarise_percentiles == TRUE) {

    time <- lapply(outcome_list, "[[", "time")
    scenario <- lapply(outcome_list, "[[", "scenario")
    median <- lapply(outcome_list, function (x) apply(x[,-c(1,2)], 1, median))
    lower <- lapply(outcome_list, function (x) apply(x[,-c(1,2)], 1, quantile, prob = 0.025))
    upper <- lapply(outcome_list, function (x) apply(x[,-c(1,2)], 1, quantile, prob = 0.975))

    summary <- Map(cbind, time, scenario, median, lower, upper)
    summary <- lapply(summary, as.data.frame)
    summary <- lapply(summary, setNames, c("time", "scenario", "median", "lower", "upper"))
    summary <- lapply(summary, function(x) {x[,"scenario"] <- scenario_label ; x})  # reassign scenario as a character

    return(summary)

  } else if (summarise_percentiles == FALSE) {
    return(outcome_list)
  }
}

# Extract cumulative HBV-related deaths (for all simulations)
# Function automatically uses correct time step to match period starting at from up until (excluding) by.
# Only works for 1 time period, but can use rbind to combine results for several periods
extract_cumulative_hbv_deaths <- function(output_files, scenario_label, from_year, by_year) {

  # Extract absolute incident HBV-related deaths per timestep
  incident_hbv_deaths <- data.frame(time = head(output_files[[1]]$time,-1),
                                    deaths = tail(sapply(lapply(output_files,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                                                sapply(lapply(output_files,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                                                sapply(lapply(output_files,"[[", "treated_hbv_deaths"), "[[", "incident_number_total"),-1))
  colnames(incident_hbv_deaths)[1] <- "time"

  cum_hbv_deaths <- as.data.frame(t(apply(incident_hbv_deaths[which(incident_hbv_deaths$time == from_year):
                                                which(incident_hbv_deaths$time ==(by_year-da)),-1],2,sum)))

  res <- data.frame(from_year = from_year,
                    by_year = by_year,
                    scenario = scenario_label,
                    cum_hbv_deaths = cum_hbv_deaths)

  return(res)
}

# Extract chronic infection (for all simulations)
# Function automatically uses correct time step to match period starting at from up until (excluding) by.
# Only works for 1 time period, but can use rbind to combine results for several periods
extract_cumulative_chronic_infections <- function(output_files, scenario_label, from_year, by_year) {

  # Extract new chronic infections per timestep
  incident_chronic_infections <- data.frame(time = head(output_files[[1]]$time,-1),
                                            deaths = tail(sapply(lapply(output_files,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                                                            sapply(lapply(output_files,"[[", "incident_chronic_infections"), "[[", "chronic_births"),-1))
  colnames(incident_chronic_infections)[1] <- "time"

  cum_chronic_infections <- as.data.frame(t(apply(incident_chronic_infections[which(incident_chronic_infections$time == from_year):
                                                                which(incident_chronic_infections$time ==(by_year-da)),-1],2,sum)))

  res <- data.frame(from_year = from_year,
                    by_year = by_year,
                    scenario = scenario_label,
                    cum_chronic_infections = cum_chronic_infections)

  return(res)
}

# Extract life years lived (for all simulations)
extract_life_years_lived <- function(output_files, scenario_label, from_year, by_year, sex_to_return = "both") {

  # Extract male and female population size per timestep
  pop_male <- data.frame(time = head(output_files[[1]]$time,-1),
                         pop = head(sapply(lapply(output_files,"[[", "pop_total"), "[[", "pop_male"),-1))
  colnames(pop_male)[1] <- "time"

  pop_female <- data.frame(time = head(output_files[[1]]$time,-1),
                         pop = head(sapply(lapply(output_files,"[[", "pop_total"), "[[", "pop_female"),-1))
  colnames(pop_female)[1] <- "time"

  # Life years lived in given period = sum of population size at each timestep * da
  # Calculating life years lived up until (excluding) the defined by_year
  life_years_male <- as.data.frame(t(apply(pop_male[which(pop_male$time == from_year):
                                                           which(pop_male$time ==(by_year-da)),-1],2,sum)))*da
  life_years_female <- as.data.frame(t(apply(pop_female[which(pop_female$time == from_year):
                                                      which(pop_female$time ==(by_year-da)),-1],2,sum)))*da

  life_years_total <- life_years_male+life_years_female

  if (sex_to_return == "both") {
    res_total <- data.frame(from_year = from_year,
                            by_year = by_year,
                            scenario = scenario_label,
                            life_years_total = life_years_total)

    return(res_total)
  } else if (sex_to_return == "male") {
    res_male <- data.frame(from_year = from_year,
                           by_year = by_year,
                           scenario = scenario_label,
                           life_years_male = life_years_male)

    return(res_male)
  } else if (sex_to_return == "female") {
    res_female <- data.frame(from_year = from_year,
                           by_year = by_year,
                           scenario = scenario_label,
                           life_years_female = life_years_female)
    return(res_female)
  } else {
    print("sex_to_return has to be male, female or both")
  }


}

# Function to compare scenarios
# Can take as input the output from the following functions:
# extract_cumulative_hbv_deaths, extract_cumulative_chronic_infections, extract_life_years_lived
# Calculation is metric1-metric2, percentage is with metric1 as reference
# For life years, take absolute value to get the life years GAINED in the scenario
calculate_number_averted <- function(counterfactual_metric, scenario_metric, summarise = TRUE) {

  # Check that time period matches between the 2 scenarios
  if (counterfactual_metric$from_year == scenario_metric$from_year &
      counterfactual_metric$by_year == scenario_metric$by_year) {

    # Calculate number averted during time period by scenario compared to counterfactual
    n_averted <- counterfactual_metric[,-c(1:3)]-scenario_metric[,-c(1:3)]
    # Calculate number averted during time period by scenario compared to counterfactual proportional to counterfactual
    prop_averted <- (counterfactual_metric[,-c(1:3)]-scenario_metric[,-c(1:3)])/counterfactual_metric[,-c(1:3)]

    if (summarise == TRUE) {

      n_averted_res <- data.frame(from_year = counterfactual_metric$from_year,
                        by_year = counterfactual_metric$by_year,
                        counterfactual = counterfactual_metric$scenario,
                        scenario = scenario_metric$scenario,
                        type = "number_averted",
                        median = apply(n_averted,1,median),
                        lower = apply(n_averted,1,quantile, prob = 0.025),
                        upper = apply(n_averted,1,quantile, prob = 0.975))

      prop_averted_res <- data.frame(from_year = counterfactual_metric$from_year,
                                  by_year = counterfactual_metric$by_year,
                                  counterfactual = counterfactual_metric$scenario,
                                  scenario = scenario_metric$scenario,
                                  type = "proportion_averted",
                                  median = apply(prop_averted,1,median),
                                  lower = apply(prop_averted,1,quantile, prob = 0.025),
                                  upper = apply(prop_averted,1,quantile, prob = 0.975))

      return(rbind(n_averted_res, prop_averted_res))

    } else if (summarise == FALSE) {

      n_averted_res <- cbind(data.frame(from_year = counterfactual_metric$from_year,
                                        by_year = counterfactual_metric$by_year,
                                        counterfactual = counterfactual_metric$scenario,
                                        scenario = scenario_metric$scenario,
                                        type = "number_averted"),
                             n_averted)


      prop_averted_res <- cbind(data.frame(from_year = counterfactual_metric$from_year,
                                           by_year = counterfactual_metric$by_year,
                                           counterfactual = counterfactual_metric$scenario,
                                           scenario = scenario_metric$scenario,
                                           type = "proportion_averted"),
                                prop_averted)

      return(rbind(n_averted_res, prop_averted_res))

    }

  } else {
    print("Time period not matching.")
  }

}

# Function to calculate average age at death for one simulation
# This function is applied to sim
calculate_average_age_at_death <- function(output_file) {

  # This function calculates the median age at death of a screened+treated cohort
  # only valid if there is only one screening event

  # Extract total cumulative deaths at the last 2 timesteps
  last_timesteps <- c(length(output_file$time)-1, length(output_file$time))

  # Add all background deaths and HBV-related deaths by age occurring after screening and treatment
  # Note this also includes deaths of those who were susceptible at test and received catch-up vaccine
  total_deaths_by_age <- output_file$full_output[last_timesteps,grepl("^cum_screened_deathsf.",names(output_file$full_output))]+
    output_file$full_output[last_timesteps,grepl("^cum_screened_deathsm.",names(output_file$full_output))]+
    output_file$full_output[last_timesteps,grepl("^cum_treated_deathsf.",names(output_file$full_output))]+
    output_file$full_output[last_timesteps,grepl("^cum_treated_deathsm.",names(output_file$full_output))]+
    output_file$full_output[last_timesteps,grepl("^cum_screened_hbv_deathsf.",names(output_file$full_output))]+
    output_file$full_output[last_timesteps,grepl("^cum_screened_hbv_deathsm.",names(output_file$full_output))]+
    output_file$full_output[last_timesteps,grepl("^cum_treated_hbv_deathsf.",names(output_file$full_output))]+
    output_file$full_output[last_timesteps,grepl("^cum_treated_hbv_deathsm.",names(output_file$full_output))]

  # Check that by the end of the simulation everyone in the cohort has died
  # (new deaths at last timestep <0.5)
  total_deaths_sum <- apply(total_deaths_by_age,1,sum)

  if (total_deaths_sum[2]-total_deaths_sum[1] >= 0.5) {

    print("Not everyone has died.")
    print(paste(total_deaths_sum[2]-total_deaths_sum[1], "new deaths"))

  # Confirm there is only 1 screening event
  } else if (length(output_file$input_parameters$screening_years)>1L) {

    print("More than one screened cohort.")

  } else {

    # Cumulative number of age-specific deaths at last timestep
    total_cum_deaths <- total_deaths_by_age[nrow(total_deaths_by_age),]
    # Multiply cumulative number of deaths at each age by age at deaths
    # Sum and divide by total number of cumulative deaths
    median_age_at_death <- sum(total_cum_deaths*ages)/sum(total_cum_deaths)

    return(median_age_at_death)

  }

}
# Function to apply to multiple sims
summarise_average_age_at_death <- function(output_files, scenario_label) {

  median_age_at_death <- as.data.frame(t(sapply(output_files, calculate_average_age_at_death)))

  res <- cbind(data.frame(scenario = scenario_label),
               median_age_at_death)

  return(res)

}

# Examples of running the function

# Using my preset median, lower and upper mock set:
#t <- summarise_time_series(out, "test", summarise_percentiles = FALSE)
#t <- lapply(t, setNames, c("time", "scenario", "median", "lower", "upper"))

#monit <- extract_life_years_lived(out, scenario_label="with_monitoring", from_year=2020, by_year=2050)
#no_monit <- extract_life_years_lived(out2, scenario_label="no_monitoring", from_year=2020, by_year=2050)
# Using my preset median, lower and upper mock set:
#avert <- calculate_number_averted(no_monit, monit, summarise = FALSE)
#colnames(avert) <- c("from_year", "by_year", "counterfactual", "scenario", "type", "median", "lower", "upper")
