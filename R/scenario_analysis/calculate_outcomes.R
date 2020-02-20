# Functions to calculate analysis outcomes

# Calculate healthcare interactions
# Function can be applied to run_model output
calculate_healthcare_interactions <- function(sim, from_year_vector, by_year_vector) {

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
                    total_immediate_treatment_initiations = total_immediate_treatment_initiations[total_immediate_treatment_initiations != 0])

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
                     cum_monitoring_treatment_initiations_it)

  return(list(res1,res2))

}

calculate_healthcare_interactions(sim[[1]], from_year_vector = 2020, by_year_vector = c(2030,2050))
