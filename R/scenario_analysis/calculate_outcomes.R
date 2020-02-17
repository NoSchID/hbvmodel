# Functions to calculate analysis outcomes

# Calculate healthcare interactions
# Function can be applied to run_model output
calculate_healthcare_interactions <- function(sim) {

  # Timing of screening programme
  time_of_screening <- sim$input_parameters$screening_years

  # Number of HBsAg tests
  total_sag_tests <- unique(sim$out$total_screened_uninfected + sim$out$total_screened_chb +
                              sim$out$total_screened_cirrhosis +
                              sim$out$total_screened_ineligible)
  # Number of liver disease assessments
  total_identified_as_ineligible <- unique(sim$out$total_screened_ineligible) *
    sim$input_parameters$link_to_care_prob

  total_identified_as_eligible <-
    unique(sim$out$total_screened_chb +
             sim$out$total_screened_cirrhosis) * sim$input_parameters$link_to_care_prob

  # Number of people immediately starting treatment
  total_immediate_treatment_initiations <- total_identified_as_eligible *
    sim$input_parameters$treatment_initiation_prob

  # Combine into dataframe
  res <- data.frame(screening_years = time_of_screening,
                    total_screened = total_sag_tests[total_sag_tests != 0],
                    total_identified_as_ineligible = total_identified_as_ineligible[total_identified_as_ineligible != 0],
                    total_identified_as_eligible = total_identified_as_eligible[total_identified_as_eligible != 0],
                    total_immediate_treatment_initiations = total_immediate_treatment_initiations[total_immediate_treatment_initiations != 0])

  return(res)

}

