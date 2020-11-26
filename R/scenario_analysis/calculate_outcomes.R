### FUNCTIONS FOR ANALYSIS OF MODEL RESULTS ###

# Read in data for calculation of YLL/DALYS ----
life_expectancy_female <- read.csv(file = here("data", "remaining_life_expectancy_female.csv"))
life_expectancy_male <- read.csv(file = here("data", "remaining_life_expectancy_male.csv"))

# Functions to calculate analysis outcomes ----

# These are all applied to the output from code_model_output
# If the argument is output_file, it is applied to 1 parmset,
# if if the argument is output_files, it is applied to a list of parmsets
# They are already used within the simulation function that runs on the cluster

# Calculate healthcare interactions (note IT interactions are recorded separately):
# 2 separate functions for the interactions that occur on screening vs. monitoring
calculate_screening_interactions <- function(output_file, scenario_label) {

  # Timing of screening programme
  if (output_file$input_parameters$apply_repeat_screen == 0) {
    time_of_screening <- output_file$input_parameters$screening_years
  } else if (output_file$input_parameters$apply_repeat_screen == 1) {
    time_of_screening <- c(output_file$input_parameters$screening_years,
                           output_file$input_parameters$repeat_screening_years)
  }

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

  # Add interactions together for all years >= from_year and < by_year

  if (length(unique(screening_interactions[[1]]$screening_years))>1) {

  total_screened <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                                 sapply(screening_interactions, "[[", "total_screened"))

  total_assessed_immediately <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                                             sapply(screening_interactions, "[[", "total_identified_as_ineligible")+
                                                 sapply(screening_interactions, "[[", "total_identified_as_eligible"))

  total_treated_immediately <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                                            sapply(screening_interactions, "[[", "total_immediate_treatment_initiations"))

  total_screened <- apply(total_screened[total_screened$screening_years>=from_year &
                                           total_screened$screening_years<by_year,-1],2,sum)

  total_assessed_immediately <- apply(total_assessed_immediately[
    total_assessed_immediately$screening_years>=from_year &
      total_assessed_immediately$screening_years<by_year,-1],2,sum)

  total_treated_immediately <- apply(total_treated_immediately[
    total_treated_immediately$screening_years >= from_year &
      total_treated_immediately$screening_years<by_year,-1],2,sum)

  total_screened_res <- cbind(data.frame(from_year = from_year,
                                         by_year = by_year,
                                         scenario = scenario_label),
                              t(total_screened))

  # Interactions during and after monitoring
  monitoring_interactions <- data.frame(sapply(output_files, calculate_monitoring_interactions, from_year, by_year, scenario_label))
  monitoring_interactions <- as.data.frame(apply(monitoring_interactions,2,unlist))

  total_monitored <- as.numeric(monitoring_interactions["cum_monitoring_events_it",])+
    as.numeric(monitoring_interactions["cum_monitoring_events_ir",])+
    as.numeric(monitoring_interactions["cum_monitoring_events_enchb",])+
    as.numeric(monitoring_interactions["cum_monitoring_events_cc",])+
    as.numeric(monitoring_interactions["cum_monitoring_events_dcc",])+
    as.numeric(monitoring_interactions["cum_monitoring_events_ineligible",])   # ignore rowname

  total_treated_after_monitoring <- as.numeric(monitoring_interactions["cum_monitoring_treatment_initiations",])+
    as.numeric(monitoring_interactions["cum_monitoring_treatment_initiations_it",])

  total_assessed <- total_assessed_immediately+total_monitored
  rownames(total_assessed) <- NULL
  total_treated <- total_treated_immediately+total_treated_after_monitoring
  rownames(total_treated) <- NULL

  total_assessed_res <- cbind(data.frame(from_year = from_year,
                                         by_year = by_year,
                                         scenario = scenario_label),
                              t(total_assessed))
  total_treated_res <- cbind(data.frame(from_year = from_year,
                                        by_year = by_year,
                                        scenario = scenario_label),
                             t(total_treated))


  total_interactions_res <- cbind(data.frame(from_year = from_year,
                                             by_year = by_year,
                                             scenario = scenario_label),
                                  t(total_screened)+t(total_assessed)+t(total_treated))


  res <- list(total_interactions = total_interactions_res,
              total_screened = total_screened_res,
              total_assessed = total_assessed_res,
              total_treated = total_treated_res)

  } else {

    total_screened <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                                 t(sapply(screening_interactions, "[[", "total_screened")))

    total_assessed_immediately <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                                             t(sapply(screening_interactions, "[[", "total_identified_as_ineligible")+
                                                 sapply(screening_interactions, "[[", "total_identified_as_eligible")))

    total_treated_immediately <- data.frame(screening_years = screening_interactions[[1]]$screening_years,
                                            t(sapply(screening_interactions, "[[", "total_immediate_treatment_initiations")))


    total_screened <- total_screened[,-1]
    total_assessed_immediately <- total_assessed_immediately[,-1]
    total_treated_immediately <- total_treated_immediately [,-1]

    total_screened_res <- cbind(data.frame(from_year = from_year,
                                           by_year = by_year,
                                           scenario = scenario_label),
                                total_screened)

    # Interactions during and after monitoring
    monitoring_interactions <- data.frame(sapply(output_files, calculate_monitoring_interactions, from_year, by_year, scenario_label))
    monitoring_interactions <- as.data.frame(apply(monitoring_interactions,2,unlist))

    total_monitored <- as.numeric(monitoring_interactions["cum_monitoring_events_it",])+
      as.numeric(monitoring_interactions["cum_monitoring_events_ir",])+
      as.numeric(monitoring_interactions["cum_monitoring_events_enchb",])+
      as.numeric(monitoring_interactions["cum_monitoring_events_cc",])+
      as.numeric(monitoring_interactions["cum_monitoring_events_dcc",])+
      as.numeric(monitoring_interactions["cum_monitoring_events_ineligible",])   # ignore rowname

    total_treated_after_monitoring <- as.numeric(monitoring_interactions["cum_monitoring_treatment_initiations",])+
      as.numeric(monitoring_interactions["cum_monitoring_treatment_initiations_it",])

    total_assessed <- total_assessed_immediately+total_monitored
    rownames(total_assessed) <- NULL
    total_treated <- total_treated_immediately+total_treated_after_monitoring
    rownames(total_treated) <- NULL

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
                                    total_screened+total_assessed+total_treated)


    res <- list(total_interactions = total_interactions_res,
                total_screened = total_screened_res,
                total_assessed = total_assessed_res,
                total_treated = total_treated_res)

  }

  return(res)


}

# Outcomes: total screened (HBsAg)
# Total liver disease assessments
# Total treatment initiations

# Function to calculate age-standardised rate of HBV-related deaths
calculate_age_standardised_hbv_deaths_rate <- function(output_file) {
  # Age-standardised incidence of HBV-related deaths per 100000 per timestep

  # a) Calculate crude age-specific rates per person-year: need age-specific number of deaths and age-specific population size
 deaths_by_age <-
   output_file$full_output[,grepl("^cum_hbv_deathsf.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_hbv_deathsm.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_screened_hbv_deathsf.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_screened_hbv_deathsm.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_treated_hbv_deathsf.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_treated_hbv_deathsm.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_negative_hbv_deathsf.",names(output_file$full_output))]+
   output_file$full_output[,grepl("^cum_negative_hbv_deathsm.",names(output_file$full_output))]

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
  total_chronic_incidence <-
    data.frame(total_chronic_infections =
                 tail(output_file$incident_chronic_infections$horizontal_chronic_infections+
                        output_file$incident_chronic_infections$horizontal_negative_chronic_infections+
                        output_file$incident_chronic_infections$chronic_births,-1),
               chronic_births = tail(output_file$incident_chronic_infections$chronic_births,-1))
  total_chronic_incidence$total_chronic_infections_rate <-
    total_chronic_incidence$total_chronic_infections/head(output_file$pop_total$pop_total,-1)

  # HBV-related deaths per timestep and per population
  total_hbv_deaths <- data.frame(total_hbv_deaths = tail(output_file$hbv_deaths$incident_number_total+
                                   output_file$screened_hbv_deaths$incident_number_total+
                                   output_file$treated_hbv_deaths$incident_number_total+
                                     output_file$negative_hbv_deaths$incident_number_total, -1),
                                 hbv_deaths_male = tail(output_file$hbv_deaths$incident_number_male+
                                   output_file$screened_hbv_deaths$incident_number_male+
                                   output_file$treated_hbv_deaths$incident_number_male+
                                     output_file$negative_hbv_deaths$incident_number_male, -1))
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

  outcome_list <- list()

  # Except for first 2 columns (time and scenario), extract each outcome into a separate list element
  for (i in 3:ncol(timeseries[[1]])) {
    outcome_list[[i-2]] <- as.data.frame(sapply(timeseries, "[[", i))
  }

  time <- lapply(timeseries, "[[", "time")[1:length(outcome_list)]
  scenario <- lapply(timeseries, "[[", "scenario")[1:length(outcome_list)]

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
                                                sapply(lapply(output_files,"[[", "treated_hbv_deaths"), "[[", "incident_number_total")+
                                                sapply(lapply(output_files,"[[", "negative_hbv_deaths"), "[[", "incident_number_total"),-1))
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
  incident_chronic_infections <-
    data.frame(time = head(output_files[[1]]$time,-1),
               deaths = tail(sapply(lapply(output_files,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                               sapply(lapply(output_files,"[[", "incident_chronic_infections"), "[[", "horizontal_negative_chronic_infections")+
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

# Extract years of life lost YLL (for all simulations)
extract_yll_and_dalys <- function(output_files, scenario_label, from_year, by_year, sex_to_return = "both",
                          disability_weight_dcc=0.178, disability_weight_hcc=0.54) {

  # Extract full output from simulations
  out <- lapply(output_files, "[[", "full_output")
  timevec <- out[[1]]$time
  sim_names <- names(output_files)

  # Extract cumulative number of HBV-related deaths (from cirrhosis and HCC),
  # person-time in DCC and HCC compartments
  cum_hbv_deathsf <- list()
  cum_hbv_deathsm <- list()
  dcc_prevf <- list()
  dcc_prevm <- list()
  hcc_prevf <- list()
  hcc_prevm <- list()
  for (i in 1:length(out)) {
    cum_hbv_deathsf[[i]] <- out[[i]][,grepl("^cum_hbv_deathsf.",names(out[[i]]))] +
      out[[i]][,grepl("^cum_screened_hbv_deathsf.",names(out[[i]]))] +
      out[[i]][,grepl("^cum_treated_hbv_deathsf.",names(out[[i]]))] +
      out[[i]][,grepl("^cum_negative_hbv_deathsf.",names(out[[i]]))]
    cum_hbv_deathsm[[i]] <- out[[i]][,grepl("^cum_hbv_deathsm.",names(out[[i]]))] +
      out[[i]][,grepl("^cum_screened_hbv_deathsm.",names(out[[i]]))] +
      out[[i]][,grepl("^cum_treated_hbv_deathsm.",names(out[[i]]))] +
      out[[i]][,grepl("^cum_negative_hbv_deathsm.",names(out[[i]]))]
    dcc_prevf[[i]] <- out[[i]][,grepl("^DCCf.",names(out[[i]]))] +
      out[[i]][,grepl("^S_DCCf.",names(out[[i]]))] +
      out[[i]][,grepl("^T_DCCf.",names(out[[i]]))] +
      out[[i]][,grepl("^N_DCCf.",names(out[[i]]))]
    dcc_prevm[[i]] <- out[[i]][,grepl("^DCCm.",names(out[[i]]))] +
      out[[i]][,grepl("^S_DCCm.",names(out[[i]]))] +
      out[[i]][,grepl("^T_DCCm.",names(out[[i]]))] +
      out[[i]][,grepl("^N_DCCm.",names(out[[i]]))]
    hcc_prevf[[i]] <- out[[i]][,grepl("^HCCf.",names(out[[i]]))] +
      out[[i]][,grepl("^S_HCCf.",names(out[[i]]))] +
      out[[i]][,grepl("^T_HCCf.",names(out[[i]]))] +
      out[[i]][,grepl("^N_HCCf.",names(out[[i]]))]
    hcc_prevm[[i]] <- out[[i]][,grepl("^DCCm.",names(out[[i]]))] +
      out[[i]][,grepl("^S_HCCm.",names(out[[i]]))] +
      out[[i]][,grepl("^T_HCCm.",names(out[[i]]))] +
      out[[i]][,grepl("^N_HCCm.",names(out[[i]]))]
  }

  # Remove out
  rm(out)
  gc()

  # Turn into incidence: new cases since the last timestep (i.e. incidence BY timestep)
  inc_hbv_deathsf <- lapply(cum_hbv_deathsf, calculate_incident_numbers)
  inc_hbv_deathsm <- lapply(cum_hbv_deathsm, calculate_incident_numbers)

  # Calculate total years of life lost between from_year and by_year
  # Life expectancy already in units of years
  # Steps: extract HBV deaths in each year of interest, multiply by respective life expectency
  # and sum across ages and timesteps for each simulation
  yll_female <- sapply(lapply(lapply(inc_hbv_deathsf, function(x) x[which(timevec == from_year+da):which(timevec == by_year),]),
         `*`, life_expectancy_female[which(life_expectancy_female$time==from_year+da):
                                       which(life_expectancy_female$time==by_year),-1]), sum)
  yll_male <- sapply(lapply(lapply(inc_hbv_deathsm, function(x) x[which(timevec == from_year+da):which(timevec == by_year),]),
                  `*`, life_expectancy_male[which(life_expectancy_male$time==from_year+da):
                                                which(life_expectancy_male$time==by_year),-1]), sum)
  yll <- yll_female+yll_male

  yll_output_female <- as.data.frame(t(yll_female))
  colnames(yll_output_female) <- sim_names
  yll_output_male <- as.data.frame(t(yll_male))
  colnames(yll_output_male) <- sim_names
  yll_output <- as.data.frame(t(yll))
  colnames(yll_output) <- sim_names

  # Calculate disability-adjusted years of life (years of life lost to disability)
  # between from_year and by_year
  # Steps: extract number in each compartment at each timestep of interest,
  # sum across ages and timesteps, multiply by da to convert to years, multiply by respective
  # disability weight

  yld_female <- sapply(lapply(dcc_prevf, function(x) x[which(timevec == from_year):which(timevec == by_year-da),]),sum)*
    dt*disability_weight_dcc +
    sapply(lapply(hcc_prevf, function(x) x[which(timevec == from_year):which(timevec == by_year-da),]),sum)*
    dt*disability_weight_hcc

  yld_male <- sapply(lapply(dcc_prevm, function(x) x[which(timevec == from_year):which(timevec == by_year-da),]),sum)*
    dt*disability_weight_dcc +
    sapply(lapply(hcc_prevm, function(x) x[which(timevec == from_year):which(timevec == by_year-da),]),sum)*
    dt*disability_weight_hcc

  daly_female <- as.data.frame(t(yll_female+yld_female))
  colnames(daly_female) <- sim_names
  daly_male <- as.data.frame(t(yll_male+yld_male))
  colnames(daly_male) <- sim_names

  daly <- as.data.frame(t(yll_female+yld_female+yll_male+yld_male))
  colnames(daly) <- sim_names

  if (sex_to_return == "both") {
    res_total <- list(yll=data.frame(from_year = from_year,
                            by_year = by_year,
                            scenario = scenario_label,
                            yll = yll_output),
                      daly=data.frame(from_year = from_year,
                                 by_year = by_year,
                                 scenario = scenario_label,
                                 daly = daly))

    return(res_total)
  } else if (sex_to_return == "male") {
    res_male <- list(yll=data.frame(from_year = from_year,
                           by_year = by_year,
                           scenario = scenario_label,
                           yll_male = yll_output_male),
                     daly=data.frame(from_year = from_year,
                                by_year = by_year,
                                scenario = scenario_label,
                                daly_male = daly_male)
    )

    return(res_male)
  } else if (sex_to_return == "female") {
    res_female <- list(yll=data.frame(from_year = from_year,
                             by_year = by_year,
                             scenario = scenario_label,
                             yll_female = yll_output_female),
    daly=data.frame(from_year = from_year,
               by_year = by_year,
               scenario = scenario_label,
               daly_female =daly_female)
    )
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

  # Check that time period matches between the 2 scenarios (if it is a population-level metric)
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

# COHORT OUTCOMES (all evaluated in 2100)

# Function to calculate average age at death of cohort for one simulation
# This function is applied to sim
calculate_cohort_average_age_at_death <- function(output_file) {

  # This function calculates the median age at death of a screened+treated cohort
  # only valid if there is only one screening event

  # UPDATE 14/08/20: Calculate outcomes in 2100 even if not entire cohort has died at this point
  # First tried to do this before negative numbers in cohort occurs or at last timestep
  # (for each simulation individually) but actually it makes less sense to calculate this outcome
  # at different timesteps for each simulation.

  # Extract total cumulative deaths in 2100
    last_timestep <- which(output_file$time==2100)

  # Add all background deaths and HBV-related deaths by age occurring after screening and treatment
  # Note this also includes deaths of those who were susceptible at test and received catch-up vaccine
  total_deaths_by_age <-
    output_file$full_output[last_timestep,grepl("^cum_screened_deathsf.",names(output_file$full_output))]+
    output_file$full_output[last_timestep,grepl("^cum_screened_deathsm.",names(output_file$full_output))]+
    output_file$full_output[last_timestep,grepl("^cum_treated_deathsf.",names(output_file$full_output))]+
    output_file$full_output[last_timestep,grepl("^cum_treated_deathsm.",names(output_file$full_output))]+
    output_file$full_output[last_timestep,grepl("^cum_screened_hbv_deathsf.",names(output_file$full_output))]+
    output_file$full_output[last_timestep,grepl("^cum_screened_hbv_deathsm.",names(output_file$full_output))]+
    output_file$full_output[last_timestep,grepl("^cum_treated_hbv_deathsf.",names(output_file$full_output))]+
    output_file$full_output[last_timestep,grepl("^cum_treated_hbv_deathsm.",names(output_file$full_output))]

  ##########################################################################
  # No longer checking this within this function:
  # Check that by the end of the simulation everyone in the cohort has died
  # (new deaths at last timestep <0.5)
  #total_deaths_sum <- apply(total_deaths_by_age,1,sum)

  #if (total_deaths_sum[2]-total_deaths_sum[1] >= 0.5) {

  #  print("Not everyone has died.")
  #  print(paste(total_deaths_sum[2]-total_deaths_sum[1], "new deaths"))

  # Confirm there is only 1 screening event
  #} else
  ###########################################################################

  # Confirm there is only 1 screening event
  if (length(output_file$input_parameters$screening_years)>1L) {

    print("More than one screened cohort.")

  } else {

    # Cumulative number of age-specific deaths at previously defined last timestep
    total_cum_deaths <- total_deaths_by_age
    # Multiply cumulative number of deaths at each age by age at deaths
    # Sum and divide by total number of cumulative deaths
    median_age_at_death <- sum(total_cum_deaths*ages)/sum(total_cum_deaths)

    return(median_age_at_death)

  }

}
# Function to apply to multiple sims
summarise_cohort_average_age_at_death <- function(output_files, scenario_label) {

  median_age_at_death <- as.data.frame(t(sapply(output_files, calculate_cohort_average_age_at_death)))

  res <- cbind(data.frame(scenario = scenario_label),
               median_age_at_death)

  return(res)

}

# Function to calculate cumulative HBV deaths averted in a screened+treated cohort
extract_cohort_cumulative_hbv_deaths <- function(output_files, scenario_label) {

  # Confirm there is only 1 screening event
  if (length(output_files[[1]]$input_parameters$screening_years)>1L) {

    print("More than one screened cohort.")

  } else {

    # UPDATE 14/08/20: Calculate outcomes by 2100 even if not entire cohort has died at this point
    # First tried to do this before negative numbers in cohort occurs or at last timestep
    # (for each simulation individually) but actually it makes less sense to calculate this outcome
    # at different timesteps for each simulation.

    # Calculate sum of incident deaths until 2100
    last_timestep <- which(output_files[[1]]$time==2100)

    # Extract absolute incident HBV-related deaths per timestep
    incident_hbv_deaths <- data.frame(time = head(output_files[[1]]$time,-1),
                                      deaths = tail(sapply(lapply(output_files,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                                                      sapply(lapply(output_files,"[[", "treated_hbv_deaths"), "[[", "incident_number_total"),-1))
    colnames(incident_hbv_deaths)[1] <- "time"

    cum_hbv_deaths <- as.data.frame(t(apply(incident_hbv_deaths[1:last_timestep,-1],2,sum)))

    res <- data.frame(scenario = scenario_label,
                      cum_hbv_deaths = cum_hbv_deaths)

    return(res)
  }

}

# Function to calculate life-years lived in a screened+treated cohort
extract_cohort_life_years_lived <- function(output_files, scenario_label, sex_to_return = "both") {

  # Confirm there is only 1 screening event
  if (length(output_files[[1]]$input_parameters$screening_years)>1L) {

    print("More than one screened cohort.")

  } else {

    # UPDATE 14/08/20: Calculate outcomes by 2100 even if not entire cohort has died at this point
    # First tried to do this before negative numbers in cohort occurs or at last timestep
    # (for each simulation individually) but actually it makes less sense to calculate this outcome
    # at different timesteps for each simulation.

    # Calculate sum of life-years until 2100
    last_timestep <- which(output_files[[1]]$time==2100)

    cohort_pop_male <- sapply(lapply(output_files, "[[", "treated_pop_male"),rowSums)+
    sapply(lapply(output_files, "[[", "screened_pop_male"),rowSums)

    cohort_pop_female <- sapply(lapply(output_files, "[[", "treated_pop_female"),rowSums)+
    sapply(lapply(output_files, "[[", "screened_pop_female"),rowSums)

    cohort_pop <- cohort_pop_male+cohort_pop_female

    # Life years lived in given period = sum of population size at each timestep * da
    # Calculating life years lived up until (excluding) the defined by_year
    life_years_male <- as.data.frame(t(apply(cohort_pop_male[1:last_timestep,],2,sum)))*da
    life_years_female <- as.data.frame(t(apply(cohort_pop_female[1:last_timestep,],2,sum)))*da

    life_years_total <- life_years_male+life_years_female

    if (sex_to_return == "both") {
      res_total <- data.frame(scenario = scenario_label,
                              life_years_total = life_years_total)

      return(res_total)
    } else if (sex_to_return == "male") {
      res_male <- data.frame(scenario = scenario_label,
                             life_years_male = life_years_male)

      return(res_male)
    } else if (sex_to_return == "female") {
      res_female <- data.frame(scenario = scenario_label,
                               life_years_female = life_years_female)
      return(res_female)
    } else {
      print("sex_to_return has to be male, female or both")
    }

  }


}

# Function to extract initial cohort population size (upon screen)
# This cohort represents those who were screened+assessed+treated/not eligible, so these are carriers who
# have completed the liver disease assessment and received a treatment decision (either taken up treatment
# or not treatment eligible). It excludes carriers who refused liver disease assessment or treatment.
extract_cohort_size <- function(output_files, scenario_label, sex_to_return = "both") {

  screened_pop_female <- sapply(lapply(output_files, "[[", "screened_pop_female"), rowSums)
  screened_pop_male <- sapply(lapply(output_files, "[[", "screened_pop_male"), rowSums)
  treated_pop_female <- sapply(lapply(output_files, "[[", "treated_pop_female"), rowSums)
  treated_pop_male <- sapply(lapply(output_files, "[[", "treated_pop_male"), rowSums)

  # Confirm there is only 1 screening event
  if (length(output_files[[1]]$input_parameters$screening_years)>1L) {

    print("More than one screened cohort.")

  } else {

    first_timestep <- output_files[[1]]$input_parameters$screening_years+da

    cohort_size_total <- screened_pop_female[output_files[[1]]$time == first_timestep]+
      screened_pop_male[output_files[[1]]$time == first_timestep]+
      treated_pop_female[output_files[[1]]$time == first_timestep]+
      treated_pop_male[output_files[[1]]$time == first_timestep]

    cohort_size_male <- screened_pop_male[output_files[[1]]$time == first_timestep]+
      treated_pop_male[output_files[[1]]$time == first_timestep]

    cohort_size_female <- screened_pop_female[output_files[[1]]$time == first_timestep]+
      treated_pop_female[output_files[[1]]$time == first_timestep]

    if (sex_to_return == "both") {
      res_total <- data.frame(scenario = scenario_label,
                              cohort_size_total = data.frame(t(cohort_size_total)))

      return(res_total)
    } else if (sex_to_return == "male") {
      res_male <- data.frame(scenario = scenario_label,
                             cohort_size_male = data.frame(t(cohort_size_male)))

      return(res_male)
    } else if (sex_to_return == "female") {
      res_female <- data.frame(scenario = scenario_label,
                               cohort_size_female = data.frame(t(cohort_size_female)))
      return(res_female)
    } else {
      print("sex_to_return has to be male, female or both")
    }

  }


}

# Function to check the remaining cohort population size in 2100, at which cohort outcomes are evaluated
# (NOT when the entire cohort has died)
extract_cohort_size_at_outcome <- function(output_files, scenario_label) {

  last_timestep <- which(output_files[[1]]$time==2100)
  cohort_size <- 0

 for (i in 1:length(output_files)) {
   cohort_size[i] <- output_files[[i]]$infectioncat_total$screened_pop[last_timestep]+
     output_files[[i]]$infectioncat_total$treated_pop[last_timestep]
 }

  res <- cbind(data.frame(scenario = scenario_label),
               cohort_size)

  return(res)
}

# Compare scenarios for cohort
calculate_cohort_number_averted <- function(counterfactual_metric, scenario_metric, summarise = TRUE) {

  # Calculate number averted during time period by scenario compared to counterfactual
  n_averted <- counterfactual_metric[,which(colnames(counterfactual_metric)!="scenario")]-
    scenario_metric[,which(colnames(scenario_metric)!="scenario")]
  # Calculate number averted during time period by scenario compared to counterfactual proportional to counterfactual
  prop_averted <- (counterfactual_metric[,which(colnames(counterfactual_metric)!="scenario")]-
                     scenario_metric[,which(colnames(scenario_metric)!="scenario")])/
    counterfactual_metric[,which(colnames(counterfactual_metric)!="scenario")]

  if (summarise == TRUE) {

    n_averted_res <- data.frame(counterfactual = counterfactual_metric$scenario,
                                scenario = scenario_metric$scenario,
                                type = "number_averted",
                                median = apply(n_averted,1,median),
                                lower = apply(n_averted,1,quantile, prob = 0.025),
                                upper = apply(n_averted,1,quantile, prob = 0.975))

    prop_averted_res <- data.frame(counterfactual = counterfactual_metric$scenario,
                                   scenario = scenario_metric$scenario,
                                   type = "proportion_averted",
                                   median = apply(prop_averted,1,median),
                                   lower = apply(prop_averted,1,quantile, prob = 0.025),
                                   upper = apply(prop_averted,1,quantile, prob = 0.975))

    return(rbind(n_averted_res, prop_averted_res))

  } else if (summarise == FALSE) {

    n_averted_res <- cbind(data.frame(counterfactual = counterfactual_metric$scenario,
                                      scenario = scenario_metric$scenario,
                                      type = "number_averted"),
                           n_averted)


    prop_averted_res <- cbind(data.frame(counterfactual = counterfactual_metric$scenario,
                                         scenario = scenario_metric$scenario,
                                         type = "proportion_averted"),
                              prop_averted)

    return(rbind(n_averted_res, prop_averted_res))

  }

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


# Functions to conduct analysis/compare scenarios ----
## Cohort functions (for monitoring analysis only) ----

# Calculate HBV deaths averted in the cohort
# Arguments are 1 counterfactual, list of scenarios to compare to it
# Optional: specify counterfactual label for plot title
# Plotted outcome is "proportion_averted" by default, can be switched to "number_averted"
plot_hbv_deaths_averted_cohort <- function(counterfactual_object, scenario_objects,
                                           counterfactual_label = "",
                                           outcome_to_plot = "proportion_averted") {

  cohort_deaths_averted <- list()

  for (i in 1:length(scenario_objects)) {
    cohort_deaths_averted[[i]] <- calculate_cohort_number_averted(counterfactual_object$cohort_cum_hbv_deaths,
                                                                  scenario_objects[[i]]$cohort_cum_hbv_deaths, summarise = FALSE)
  }

  cohort_deaths_averted <- do.call("rbind", cohort_deaths_averted)

  cohort_deaths_averted_long_original <- gather(cohort_deaths_averted, key = "sim", value = "value",
                                                -counterfactual, -scenario, -type)

  cohort_deaths_averted_long <- cohort_deaths_averted_long_original

  # Relabel scenarios for plot
  levels(cohort_deaths_averted_long$scenario)[
    grepl("screen_2020_monit_0$", levels(cohort_deaths_averted_long$scenario))] <- "Never"
  levels(cohort_deaths_averted_long$scenario)[
    grepl("screen_2020_monit_10$", levels(cohort_deaths_averted_long$scenario))] <- "10 years"
  levels(cohort_deaths_averted_long$scenario)[
    grepl("screen_2020_monit_5$", levels(cohort_deaths_averted_long$scenario))] <- "5 years"
  levels(cohort_deaths_averted_long$scenario)[
    grepl("screen_2020_monit_1$", levels(cohort_deaths_averted_long$scenario))] <- "1 year"
  levels(cohort_deaths_averted_long$scenario)[
    grepl("monit_0_screen_20$", levels(cohort_deaths_averted_long$scenario))] <- "20 years"
  levels(cohort_deaths_averted_long$scenario)[
    grepl("monit_0_screen_10$", levels(cohort_deaths_averted_long$scenario))] <- "10 years"
  levels(cohort_deaths_averted_long$scenario)[
    grepl("monit_0_screen_5$", levels(cohort_deaths_averted_long$scenario))] <- "5 years"
  levels(cohort_deaths_averted_long$scenario)[
    grepl("monit_0_screen_1$", levels(cohort_deaths_averted_long$scenario))] <- "1 year"


  # Choose y axis label based on outcome to plot (proportion or number)
  if (outcome_to_plot == "proportion_averted") {
    y_axis_label <- "Proportion of HBV-related deaths averted"
  } else if (outcome_to_plot == "number_averted") {
    y_axis_label <- "Number of HBV-related deaths averted"
  }

  print(ggplot(cohort_deaths_averted_long[cohort_deaths_averted_long$type == outcome_to_plot,]) +
          geom_boxplot(aes(scenario, value), fill = "#F8766D", width = 0.25) +
          ylab(y_axis_label) +
          labs(title = paste0("Cohort impact compared to counterfactual:\n", counterfactual_label)) +
          xlab("Monitoring frequency") +
          #    scale_x_discrete(labels = c("10 years", "5 years", "1 year")) +
          theme_classic() +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
          theme(axis.text = element_text(size = 15),
                axis.title = element_text(size = 15),
                title = element_text(size = 15)))

  return(cohort_deaths_averted_long_original)

}

# Example applications
#cohort_deaths_averted_long <-
#  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
#                                 scenario_objects = list(out4, out5, out6),
#                                 counterfactual_label = "treatment programme without monitoring")
#cohort_deaths_averted_sq_long <-
#  plot_hbv_deaths_averted_cohort(counterfactual_object = out1,
#                                 scenario_objects = list(out3,out4, out5, out6),
#                                 counterfactual_label = "no treatment")

# Calculate life years saved in the cohort
plot_ly_gained_cohort <- function(counterfactual_object, scenario_objects,
                                  counterfactual_label = "",
                                  outcome_to_plot = "proportion_averted") {

  cohort_ly_gained <- list()

  for (i in 1:length(scenario_objects)) {
    cohort_ly_gained[[i]] <- calculate_cohort_number_averted(scenario_objects[[i]]$cohort_ly,
                                                             counterfactual_object$cohort_ly,
                                                             summarise = FALSE)
  }


  cohort_ly_gained <- do.call("rbind", cohort_ly_gained)

  cohort_ly_gained_long_original <- gather(cohort_ly_gained, key = "sim", value = "value",
                                           -counterfactual, -scenario, -type)

  cohort_ly_gained_long <- cohort_ly_gained_long_original

  # Relabel scenarios for plot
  levels(cohort_ly_gained_long$counterfactual)[
    grepl("screen_2020_monit_0$", levels(cohort_ly_gained_long$counterfactual))] <- "Never"
  levels(cohort_ly_gained_long$counterfactual)[
    grepl("screen_2020_monit_10$", levels(cohort_ly_gained_long$counterfactual))] <- "10 years"
  levels(cohort_ly_gained_long$counterfactual)[
    grepl("screen_2020_monit_5$", levels(cohort_ly_gained_long$counterfactual))] <- "5 years"
  levels(cohort_ly_gained_long$counterfactual)[
    grepl("screen_2020_monit_1$", levels(cohort_ly_gained_long$counterfactual))] <- "1 year"
  levels(cohort_ly_gained_long$counterfactual)[
    grepl("monit_0_screen_20$", levels(cohort_ly_gained_long$counterfactual))] <- "20 years"
  levels(cohort_ly_gained_long$counterfactual)[
    grepl("monit_0_screen_10$", levels(cohort_ly_gained_long$counterfactual))] <- "10 years"
  levels(cohort_ly_gained_long$counterfactual)[
    grepl("monit_0_screen_5$", levels(cohort_ly_gained_long$counterfactual))] <- "5 years"
  levels(cohort_ly_gained_long$counterfactual)[
    grepl("monit_0_screen_1$", levels(cohort_ly_gained_long$counterfactual))] <- "1 year"

  # Choose y axis label based on outcome to plot (proportion or number)
  if (outcome_to_plot == "proportion_averted") {
    y_axis_label <- "Proportion of life years saved"
  } else if (outcome_to_plot == "number_averted") {
    y_axis_label <- "Number of life years saved"
  }

  print(ggplot(cohort_ly_gained_long[cohort_ly_gained_long$type == outcome_to_plot,]) +
          geom_boxplot(aes(counterfactual, value), fill = "#F8766D", width = 0.25) +
          ylab(y_axis_label) +
          labs(title = paste0("Cohort impact compared to counterfactual:\n", counterfactual_label)) +
          xlab("Monitoring frequency") +
          theme_classic() +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
          theme(axis.text = element_text(size = 15),
                axis.title = element_text(size = 15),
                title = element_text(size = 15)))

  return(cohort_ly_gained_long_original)

}

# Example applications
#cohort_ly_gained_long <-
#  plot_ly_gained_cohort(counterfactual_object = out3,
#                        scenario_objects = list(out4, out5, out6),
#                        counterfactual_label = "treatment programme without monitoring")
#cohort_ly_gained_sq_long <-
#  plot_ly_gained_cohort(counterfactual_object = out1,
#                        scenario_objects = list(out3,out4, out5, out6),
#                        counterfactual_label = "no treatment")

## HBV deaths averted and LY saved ----
# Calculate HBV deaths averted on the population level
# Currently for 2030, 2050 and 2100 fixed

plot_hbv_deaths_averted <- function(counterfactual_object, scenario_objects,
                                    counterfactual_label = "",
                                    outcome_to_plot = "proportion_averted",
                                    x_axis = "monitoring", timepoints = c(2030,2050,2100)) {

  period_labs <- c(paste0("2020-",timepoints[1]), paste0("2020-",timepoints[2]), paste0("2020-",timepoints[3]))
  names(period_labs) <- c(as.character(timepoints[1]), as.character(timepoints[2]),
                          as.character(timepoints[3]))


  deaths_averted <- list()

  for (i in 1:length(scenario_objects)) {
    deaths_averted[[i]] <- rbind(calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                                          summarise = FALSE),
                                 calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                                          summarise = FALSE),
                                 calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                                          summarise = FALSE))

  }

  deaths_averted <- do.call("rbind", deaths_averted)

  deaths_averted_long_original <- gather(deaths_averted, key = "sim", value = "value", -from_year,
                                         -by_year, -counterfactual, -scenario, - type)

  deaths_averted_long_original$by_year <- factor(deaths_averted_long_original$by_year)

  deaths_averted_long <- deaths_averted_long_original

  # Relabel scenarios for plot
  levels(deaths_averted_long$scenario)[
    grepl("screen_2020_monit_0$", levels(deaths_averted_long$scenario))] <- "Never"
  levels(deaths_averted_long$scenario)[
    grepl("screen_2020_monit_10$", levels(deaths_averted_long$scenario))] <- "10 years"
  levels(deaths_averted_long$scenario)[
    grepl("screen_2020_monit_5$", levels(deaths_averted_long$scenario))] <- "5 years"
  levels(deaths_averted_long$scenario)[
    grepl("screen_2020_monit_1$", levels(deaths_averted_long$scenario))] <- "1 year"
  levels(deaths_averted_long$scenario)[
    grepl("monit_0_screen_20$", levels(deaths_averted_long$scenario))] <- "20 years"
  levels(deaths_averted_long$scenario)[
    grepl("monit_0_screen_10$", levels(deaths_averted_long$scenario))] <- "10 years"
  levels(deaths_averted_long$scenario)[
    grepl("monit_0_screen_5$", levels(deaths_averted_long$scenario))] <- "5 years"
  levels(deaths_averted_long$scenario)[
    grepl("monit_0_screen_1$", levels(deaths_averted_long$scenario))] <- "1 year"

  # Choose y axis label based on outcome to plot (proportion or number)
  if (outcome_to_plot == "proportion_averted") {
    y_axis_label <- "Proportion of HBV-related deaths averted"
  } else if (outcome_to_plot == "number_averted") {
    y_axis_label <- "Number of HBV-related deaths averted"
  }

  # Chose x axis label based on monitoring or screening impact analysis
  if (x_axis == "monitoring") {
    x_axis_label <- "Monitoring frequency"
  } else if (x_axis == "screening") {
    x_axis_label <- "Repeat screening frequency"
  }

  print(ggplot(deaths_averted_long[deaths_averted_long$type == outcome_to_plot,]) +
          geom_boxplot(aes(x = scenario, y = value), fill = "#00BFC4") +
          facet_wrap(~ by_year, ncol = 3, labeller=labeller(by_year = period_labs)) +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          labs(title = paste0("Population impact compared to counterfactual:\n", counterfactual_label)) +
          theme_classic() +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
          theme(axis.text = element_text(size = 15),
                axis.text.x = element_text(angle =45, hjust = 1),
                axis.title = element_text(size = 15),
                strip.text = element_text(size = 15),
                title = element_text(size = 15)))

  return(deaths_averted_long_original)

}
# Example application
#deaths_averted_long <- plot_hbv_deaths_averted(counterfactual_object = out3,
#                                               scenario_objects = list(out4, out5, out6),
#                                               counterfactual_label = "treatment programme without monitoring")

# Calculate life years saved on the population level
# Currently for 2030, 2050 and 2100 fixed
plot_ly_gained <- function(counterfactual_object, scenario_objects,
                           counterfactual_label = "",
                           outcome_to_plot = "proportion_averted",
                           x_axis = "monitoring", timepoints = c(2030,2050,2100)) {

  period_labs <- c(paste0("2020-",timepoints[1]), paste0("2020-",timepoints[2]), paste0("2020-",timepoints[3]))
  names(period_labs) <- c(as.character(timepoints[1]), as.character(timepoints[2]),
                          as.character(timepoints[3]))

  ly_gained <- list()

  for (i in 1:length(scenario_objects)) {
    ly_gained[[i]] <- rbind(calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                                     counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                                     summarise = FALSE),
                            calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                                     counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                                     summarise = FALSE),
                            calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                                     counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                                     summarise = FALSE))

  }

  ly_gained <- do.call("rbind", ly_gained)

  ly_gained_long_original <- gather(ly_gained, key = "sim", value = "value", -from_year,
                                    -by_year, -counterfactual, -scenario, - type)

  ly_gained_long_original$by_year <- factor(ly_gained_long_original$by_year)

  ly_gained_long <- ly_gained_long_original

  # Relabel coutnerfactuals for plot
  levels(ly_gained_long$counterfactual)[
    grepl("screen_2020_monit_0$", levels(ly_gained_long$counterfactual))] <- "Never"
  levels(ly_gained_long$counterfactual)[
    grepl("screen_2020_monit_10$", levels(ly_gained_long$counterfactual))] <- "10 years"
  levels(ly_gained_long$counterfactual)[
    grepl("screen_2020_monit_5$", levels(ly_gained_long$counterfactual))] <- "5 years"
  levels(ly_gained_long$counterfactual)[
    grepl("screen_2020_monit_1$", levels(ly_gained_long$counterfactual))] <- "1 year"
  levels(ly_gained_long$counterfactual)[
    grepl("monit_0_screen_20$", levels(ly_gained_long$counterfactual))] <- "20 years"
  levels(ly_gained_long$counterfactual)[
    grepl("monit_0_screen_10$", levels(ly_gained_long$counterfactual))] <- "10 years"
  levels(ly_gained_long$counterfactual)[
    grepl("monit_0_screen_5$", levels(ly_gained_long$counterfactual))] <- "5 years"
  levels(ly_gained_long$counterfactual)[
    grepl("monit_0_screen_1$", levels(ly_gained_long$counterfactual))] <- "1 year"



  # Choose y axis label based on outcome to plot (proportion or number)
  if (outcome_to_plot == "proportion_averted") {
    y_axis_label <- "Proportion of life years saved"
  } else if (outcome_to_plot == "number_averted") {
    y_axis_label <- "Number of life years saved"
  }

  # Chose x axis label based on monitoring or screening impact analysis
  if (x_axis == "monitoring") {
    x_axis_label <- "Monitoring frequency"
  } else if (x_axis == "screening") {
    x_axis_label <- "Repeat screening frequency"
  }

  print(ggplot(ly_gained_long[ly_gained_long$type == outcome_to_plot,]) +
          geom_boxplot(aes(x = counterfactual, y = value), fill = "#00BFC4") +
          facet_wrap(~ by_year, ncol = 3, labeller=labeller(by_year = period_labs)) +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          labs(title = paste0("Population impact compared to counterfactual:\n", counterfactual_label)) +
          theme_classic() +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
          theme(axis.text = element_text(size = 15),
                axis.text.x = element_text(angle =45, hjust = 1),
                axis.title = element_text(size = 15),
                strip.text = element_text(size = 15),
                title = element_text(size = 15)))

  return(ly_gained_long_original)

}

# Example application
#ly_gained_long <- plot_ly_gained(counterfactual_object = out3,
#                                 scenario_objects = list(out4, out5, out6),
#                                 counterfactual_label = "treatment programme without monitoring")

## HBV deaths averted and LY saved per healthcare interactions ----

# Function to calculate incremental healthcare interactions per HBV deaths averted (or opposite)
# for counterfactual not being current status quo
# Currently for 2030, 2050 and 2100 fixed
# Interactions seelcted with interaction_type: "total_interactions" (default),
# "total_screened" (incremental HBsAg tests only)
# "total_assessed" (incremental liver disease assessments only) and
# "total_treated" (incremental treatment initiations only)

plot_hbv_deaths_averted_per_healthcare_interaction <- function(counterfactual_object, scenario_objects,
                                                               interaction_type = "total_interactions",
                                                               counterfactual_label = "",
                                                               x_axis = "monitoring",
                                                               timepoints = c(2030,2050,2100)) {

  period_labs <- c(paste0("2020-",timepoints[1]), paste0("2020-",timepoints[2]), paste0("2020-",timepoints[3]))
  names(period_labs) <- c(as.character(timepoints[1]), as.character(timepoints[2]),
                          as.character(timepoints[3]))

  # Calculating HBV deaths averted per interaction but plotting the opposite
  deaths_averted_per_interaction <- list()

  # Check whether interactions are recorded in counterfactual_object (this would not be the case if no screening)

  if (is.na(counterfactual_object$interactions)==TRUE) {

    for (i in 1:length(scenario_objects)) {
      deaths_averted_per_interaction[[i]] <- data.frame(rbind(
        c(by_year = timepoints[1], scenario = as.character(scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[1])]]$scenario),
          unlist(calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[1])]][[interaction_type]][,-c(1:3)]))),
        c(by_year = timepoints[2], scenario = as.character(scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[2])]]$scenario),
          unlist(calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[2])]][[interaction_type]][,-c(1:3)]))),
        c(by_year = timepoints[3], scenario = as.character(scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[3])]]$scenario),
          unlist(calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[3])]][[interaction_type]][,-c(1:3)])))
      ))

    }

  } else {

    for (i in 1:length(scenario_objects)) {
      deaths_averted_per_interaction[[i]] <- data.frame(rbind(
        c(by_year = timepoints[1], scenario = as.character(scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[1])]]$scenario),
          unlist(calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[1])]][[interaction_type]][,-c(1:3)]-
                      counterfactual_object$interactions[[which(seq(2025,2100, by = 5)==timepoints[1])]][[interaction_type]][,-c(1:3)]))),
        c(by_year = timepoints[2], scenario = as.character(scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[2])]]$scenario),
          unlist(calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[2])]][[interaction_type]][,-c(1:3)]-
                      counterfactual_object$interactions[[which(seq(2025,2100, by = 5)==timepoints[2])]][[interaction_type]][,-c(1:3)]))),
        c(by_year = timepoints[3], scenario = as.character(scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[3])]]$scenario),
          unlist(calculate_number_averted(counterfactual_object$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                          scenario_objects[[i]]$cum_hbv_deaths[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[3])]][[interaction_type]][,-c(1:3)]-
                      counterfactual_object$interactions[[which(seq(2025,2100, by = 5)==timepoints[3])]][[interaction_type]][,-c(1:3)])))
      ))

    }

  }


  deaths_averted_per_interaction <- do.call("rbind", deaths_averted_per_interaction)

  deaths_averted_per_interaction_long_original <- gather(deaths_averted_per_interaction, key = "sim",
                                                         value = "value", -scenario, -by_year)
  deaths_averted_per_interaction_long_original$value <-
    as.numeric(deaths_averted_per_interaction_long_original$value)

  deaths_averted_per_interaction_long_original$by_year <- as.factor(deaths_averted_per_interaction_long_original$by_year)

  # Add column to indicate type of healthcare interaction
  deaths_averted_per_interaction_long_original$interaction_type <- interaction_type

  deaths_averted_per_interaction_long <- deaths_averted_per_interaction_long_original

  # Relabel scenarios for plot
  levels(deaths_averted_per_interaction_long$scenario)[
    grepl("screen_2020_monit_0$", levels(deaths_averted_per_interaction_long$scenario))] <- "Never"
  levels(deaths_averted_per_interaction_long$scenario)[
    grepl("screen_2020_monit_10$", levels(deaths_averted_per_interaction_long$scenario))] <- "10 years"
  levels(deaths_averted_per_interaction_long$scenario)[
    grepl("screen_2020_monit_5$", levels(deaths_averted_per_interaction_long$scenario))] <- "5 years"
  levels(deaths_averted_per_interaction_long$scenario)[
    grepl("screen_2020_monit_1$", levels(deaths_averted_per_interaction_long$scenario))] <- "1 year"
  levels(deaths_averted_per_interaction_long$scenario)[
    grepl("monit_0_screen_20$", levels(deaths_averted_per_interaction_long$scenario))] <- "20 years"
  levels(deaths_averted_per_interaction_long$scenario)[
    grepl("monit_0_screen_10$", levels(deaths_averted_per_interaction_long$scenario))] <- "10 years"
  levels(deaths_averted_per_interaction_long$scenario)[
    grepl("monit_0_screen_5$", levels(deaths_averted_per_interaction_long$scenario))] <- "5 years"
  levels(deaths_averted_per_interaction_long$scenario)[
    grepl("monit_0_screen_1$", levels(deaths_averted_per_interaction_long$scenario))] <- "1 year"


  # Choose y axis label based on type of interaction
  if (interaction_type == "total_interactions") {
    y_axis_label <- "Incremental healthcare interactions\nper HBV death averted"
  } else if (interaction_type == "total_screened") {
    y_axis_label <- "Incremental HBsAg tests (screening)\nper HBV death averted"
  } else if (interaction_type == "total_assessed") {
    y_axis_label <- "Incremental treatment eligibility assessments\nper HBV death averted"
  } else if (interaction_type == "total_treated") {
    y_axis_label <- "Incremental treatment initiations\nper HBV death averted"
  }

  # Chose x axis label based on monitoring or screening impact analysis
  if (x_axis == "monitoring") {
    x_axis_label <- "Monitoring frequency"
  } else if (x_axis == "screening") {
    x_axis_label <- "Repeat screening frequency"
  }

  print(ggplot(data = deaths_averted_per_interaction_long) +
          geom_boxplot(aes(x=scenario, y=1/value), fill = "#00BFC4") +
          facet_wrap(~by_year, ncol = 3, labeller=labeller(by_year = period_labs),scales = "free_y") +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          labs(title = paste0("Population impact compared to counterfactual:\n", counterfactual_label)) +
          theme_classic() +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
          theme(axis.text = element_text(size = 15),
                axis.text.x = element_text(angle =45, hjust = 1),
                axis.title = element_text(size = 15),
                strip.text = element_text(size = 15),
                title = element_text(size = 15)))

  return(deaths_averted_per_interaction_long_original)

}

# Example application
#deaths_averted_per_interaction_long <-
#  plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out3,
#                                                    scenario_objects = list(out4,out5, out6),
#                                                     interaction_type = "total_interactions",
#                                                     counterfactual_label = "treatment programme without monitoring")


# Function to calculate incremental healthcare interactions per life year saved (or opposite)
# for counterfactual not being current status quo
# Currently for 2030, 2050 and 2100 fixed
# Interactions seelcted with interaction_type: "total_interactions" (default),
# "total_screened" (incremental HBsAg tests only)
# "total_assessed" (incremental liver disease assessments only) and
# "total_treated" (incremental treatment initiations only)
plot_ly_gained_per_healthcare_interaction <- function(counterfactual_object, scenario_objects,
                                                      interaction_type = "total_interactions",
                                                      counterfactual_label = "",
                                                      x_axis = "monitoring",
                                                      timepoints = c(2030,2050,2100)) {

  period_labs <- c(paste0("2020-",timepoints[1]), paste0("2020-",timepoints[2]), paste0("2020-",timepoints[3]))
  names(period_labs) <- c(as.character(timepoints[1]), as.character(timepoints[2]),
                          as.character(timepoints[3]))

  ly_gained_per_interaction <- list()

  if (is.na(counterfactual_object$interactions)==TRUE) {

    for (i in 1:length(scenario_objects)) {
      ly_gained_per_interaction[[i]] <- data.frame(rbind(
        c(by_year = timepoints[1], scenario = as.character(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[1])]]$scenario),
          unlist(calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                          counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[1])]][[interaction_type]][,-c(1:3)]))),
        c(by_year = timepoints[2], scenario = as.character(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[2])]]$scenario),
          unlist(calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                          counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[2])]][[interaction_type]][,-c(1:3)]))),
        c(by_year = timepoints[3], scenario = as.character(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[3])]]$scenario),
          unlist(calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                          counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[3])]][[interaction_type]][,-c(1:3)])))
      ))

    }

  } else {

    for (i in 1:length(scenario_objects)) {
      ly_gained_per_interaction[[i]] <- data.frame(rbind(
        c(by_year = timepoints[1], scenario = as.character(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[1])]]$scenario),
          unlist(calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                          counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[1])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[1])]][[interaction_type]][,-c(1:3)]-
                      counterfactual_object$interactions[[which(seq(2025,2100, by = 5)==timepoints[1])]][[interaction_type]][,-c(1:3)]))),
        c(by_year = timepoints[2], scenario = as.character(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[2])]]$scenario),
          unlist(calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                          counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[2])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[2])]][[interaction_type]][,-c(1:3)]-
                      counterfactual_object$interactions[[which(seq(2025,2100, by = 5)==timepoints[2])]][[interaction_type]][,-c(1:3)]))),
        c(by_year = timepoints[3], scenario = as.character(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[3])]]$scenario),
          unlist(calculate_number_averted(scenario_objects[[i]]$ly[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                          counterfactual_object$ly[[which(seq(2025,2100, by = 5)==timepoints[3])]],
                                          summarise = FALSE)[1,-c(1:5)]/
                   (scenario_objects[[i]]$interactions[[which(seq(2025,2100, by = 5)==timepoints[3])]][[interaction_type]][,-c(1:3)]-
                      counterfactual_object$interactions[[which(seq(2025,2100, by = 5)==timepoints[3])]][[interaction_type]][,-c(1:3)])))
      ))

    }

  }

  ly_gained_per_interaction <- do.call("rbind", ly_gained_per_interaction)

  ly_gained_per_interaction_long_original <- gather(ly_gained_per_interaction, key = "sim",
                                                    value = "value", -scenario, -by_year)
  ly_gained_per_interaction_long_original$value <-
    as.numeric(ly_gained_per_interaction_long_original$value)

  ly_gained_per_interaction_long_original$by_year <-
    as.factor(ly_gained_per_interaction_long_original$by_year)

  # Add column to indicate type of healthcare interaction
  ly_gained_per_interaction_long_original$interaction_type <- interaction_type

  ly_gained_per_interaction_long <- ly_gained_per_interaction_long_original

  # Relabel scenarios for plot
  levels(ly_gained_per_interaction_long$scenario)[
    grepl("screen_2020_monit_0$", levels(ly_gained_per_interaction_long$scenario))] <- "Never"
  levels(ly_gained_per_interaction_long$scenario)[
    grepl("screen_2020_monit_10$", levels(ly_gained_per_interaction_long$scenario))] <- "10 years"
  levels(ly_gained_per_interaction_long$scenario)[
    grepl("screen_2020_monit_5$", levels(ly_gained_per_interaction_long$scenario))] <- "5 years"
  levels(ly_gained_per_interaction_long$scenario)[
    grepl("screen_2020_monit_1$", levels(ly_gained_per_interaction_long$scenario))] <- "1 year"
  levels(ly_gained_per_interaction_long$scenario)[
    grepl("monit_0_screen_20$", levels(ly_gained_per_interaction_long$scenario))] <- "20 years"
  levels(ly_gained_per_interaction_long$scenario)[
    grepl("monit_0_screen_10$", levels(ly_gained_per_interaction_long$scenario))] <- "10 years"
  levels(ly_gained_per_interaction_long$scenario)[
    grepl("monit_0_screen_5$", levels(ly_gained_per_interaction_long$scenario))] <- "5 years"
  levels(ly_gained_per_interaction_long$scenario)[
    grepl("monit_0_screen_1$", levels(ly_gained_per_interaction_long$scenario))] <- "1 year"


  # Choose y axis label based on type of interaction
  if (interaction_type == "total_interactions") {
    y_axis_label <- "Incremental healthcare interactions\nper life year saved"
  } else if (interaction_type == "total_screened") {
    y_axis_label <- "Incremental HBsAg tests (screening)\nper life year saved"
  } else if (interaction_type == "total_assessed") {
    y_axis_label <- "Incremental treatment eligibility assessments\nper life year saved"
  } else if (interaction_type == "total_treated") {
    y_axis_label <- "Incremental treatment initiations\nper life year saved"
  }

  # Chose x axis label based on monitoring or screening impact analysis
  if (x_axis == "monitoring") {
    x_axis_label <- "Monitoring frequency"
  } else if (x_axis == "screening") {
    x_axis_label <- "Repeat screening frequency"
  }


  print(ggplot(data = ly_gained_per_interaction_long) +
          geom_boxplot(aes(x=scenario, y=1/value), fill = "#00BFC4") +
          facet_wrap(~by_year, ncol = 3, labeller=labeller(by_year = period_labs), scales="free_y") +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          labs(title = paste0("Population impact compared to counterfactual:\n", counterfactual_label)) +
          theme_classic() +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
          theme(axis.text = element_text(size = 15),
                axis.text.x = element_text(angle =45, hjust = 1),
                axis.title = element_text(size = 15),
                strip.text = element_text(size = 15),
                title = element_text(size = 15)))

  return(ly_gained_per_interaction_long_original)

}

# Example application
#ly_gained_per_interaction_long <-
#  plot_ly_gained_per_healthcare_interaction(counterfactual_object = out3,
#                                            scenario_objects = list(out4,out5, out6),
#                                            interaction_type = "total_interactions",
#                                            counterfactual_label = "treatment programme without monitoring")


