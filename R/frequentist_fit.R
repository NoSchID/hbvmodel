# Fit the model in a frequentist framework using least squares

library(here)
source(here("R", "imperial_model_calibration.R"))

# Function to calculate relative sum of squares
calculate_sos <- function(mapped_output) {
  mapped_output_for_error <- lapply(mapped_output, function(x) x[!is.na(x$data_value),])
  datapoints <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$data_value)))
  #quality_weights <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$quality_weight)))
  model_prediction <- as.numeric(unlist(lapply(mapped_output_for_error, function(x) x$model_value)))

  # Calculate vector of differences between each datapoint and simulated output
  data_model_diff <- datapoints-model_prediction  # observation - prediction
  data_model_mean <- (datapoints+model_prediction)/2

  sos <- sum((data_model_diff/data_model_mean)^2)

  return(sos)
}


# Test fit for 1 parameter:
#least_squares_function <- function(pr_ic_enchb) {
#  temp <- fit_model_full_output(default_parameter_list = parameter_list,
#                                data_to_fit = calibration_datasets_list,
#                                parms_to_change = list(pr_ic_enchb = pr_ic_enchb))
#  sse <- calculate_sos(temp$mapped_output)
#  return(sse)
#}
#inits_test <- c(pr_ic_enchb = 0.05)
#optim(fn=least_squares_function, par = inits_test, method = "Brent", lower = 0, upper = 1)

## Vary all parameters

# Define prior ranges
range <- data.frame(b1 = c(0.03,0.7),
                    b2 = c(0,0.7),
                    b3 = c(0,0.7),
                    mtct_prob_s = c(0,0.3),
                    mtct_prob_e = c(0.3,0.9),
                    alpha = c(1,10),
                    p_chronic_in_mtct = c(0.5,1),
                    p_chronic_function_r = c(0.25,1),
                    p_chronic_function_s = c(0.05,0.9),
                    pr_it_ir = c(0.01,0.6),
                    pr_ir_ic = c(0,2),
                    eag_prog_function_rate = c(0,0.01),
                    pr_ir_enchb = c(0,0.2),
                    pr_ir_cc_female = c(0.005,0.05),
                    pr_ir_cc_age_threshold = c(0,15),
                    pr_ic_enchb = c(0.002,0.2),
                    sag_loss_slope = c(0.0002,0.0006),
                    pr_enchb_cc_female = c(0,0.5),
                    cirrhosis_male_cofactor = c(1,20),
                    pr_cc_dcc = c(0.01,0.09),
                    cancer_prog_coefficient_female = c(0,0.0003),
                    cancer_age_threshold = c(0,15),
                    cancer_male_cofactor = c(1,20),
                    hccr_it = c(1,19),
                    hccr_ir = c(1,100),
                    hccr_enchb = c(1,100),
                    hccr_cc = c(1,100),
                    hccr_dcc = c(0.003,0.6),
                    mu_cc = c(0,0.1),
                    mu_dcc = c(0,10),
                    mu_hcc = c(0,10),
                    vacc_eff = c(0.5,1))
range <- as.data.frame(t(range))
names(range) <- c("min", "max")
range$parm <- rownames(range)

least_squares_function_all <- function(all_parms) {
  b1 = all_parms["b1"]
  b2 = all_parms["b2"]
  b3 = all_parms["b3"]
  mtct_prob_s = all_parms["mtct_prob_s"]
  mtct_prob_e = all_parms["mtct_prob_e"]
  alpha = all_parms["alpha"]
  p_chronic_in_mtct = all_parms["p_chronic_in_mtct"]
  p_chronic_function_r = all_parms["p_chronic_function_r"]
  p_chronic_function_s = all_parms["p_chronic_function_s"]
  pr_it_ir = all_parms["pr_it_ir"]
  pr_ir_ic = all_parms["pr_ir_ic"]
  eag_prog_function_rate = all_parms["eag_prog_function_rate"]
  pr_ir_enchb = all_parms["pr_ir_enchb"]
  pr_ir_cc_female = all_parms["pr_ir_cc_female"]
  pr_ir_cc_age_threshold = all_parms["pr_ir_cc_age_threshold"]
  pr_ic_enchb = all_parms["pr_ic_enchb"]
  sag_loss_slope = all_parms["sag_loss_slope"]
  pr_enchb_cc_female = all_parms["pr_enchb_cc_female"]
  cirrhosis_male_cofactor = all_parms["cirrhosis_male_cofactor"]
  pr_cc_dcc = all_parms["pr_cc_dcc"]
  cancer_prog_coefficient_female = all_parms["cancer_prog_coefficient_female"]
  cancer_age_threshold = all_parms["cancer_age_threshold"]
  cancer_male_cofactor = all_parms["cancer_male_cofactor"]
  hccr_it = all_parms["hccr_it"]
  hccr_ir = all_parms["hccr_ir"]
  hccr_enchb = all_parms["hccr_enchb"]
  hccr_cc = all_parms["hccr_cc"]
  hccr_dcc = all_parms["hccr_dcc"]
  mu_cc = all_parms["mu_cc"]
  mu_dcc = all_parms["mu_dcc"]
  mu_hcc = all_parms["mu_hcc"]
  vacc_eff = all_parms["vacc_eff"]

  parms <- data.frame(all_parms)
  parms$parm <- rownames(parms)
  merge <- left_join(range, parms)

  if (any(merge$all_parms < merge$min) |
      any(merge$all_parms > merge$max) |
      b1 < b2 |
      b1 < b3 |
      mtct_prob_e < mtct_prob_s |
      hccr_ir > hccr_cc |
      hccr_enchb > hccr_ir |
      hccr_it > hccr_enchb) {
    sse <- 1000000000
    return(sse)
  } else {
  temp <- fit_model_full_output(default_parameter_list = parameter_list,
                                data_to_fit = calibration_datasets_list,
                                parms_to_change = list(b1 = b1,
                                                       b2 = b2,
                                                       b3 = b3,
                                                       mtct_prob_s = mtct_prob_s,
                                                       mtct_prob_e = mtct_prob_e,
                                                       alpha = alpha,
                                                       p_chronic_in_mtct = p_chronic_in_mtct,
                                                       p_chronic_function_r = p_chronic_function_r,
                                                       p_chronic_function_s = p_chronic_function_s,
                                                       pr_it_ir = pr_it_ir,
                                                       pr_ir_ic = pr_ir_ic,
                                                       eag_prog_function_rate = eag_prog_function_rate,
                                                       pr_ir_enchb = pr_ir_enchb,
                                                       pr_ir_cc_female = pr_ir_cc_female,
                                                       pr_ir_cc_age_threshold = pr_ir_cc_age_threshold,
                                                       pr_ic_enchb = pr_ic_enchb,
                                                       sag_loss_slope = sag_loss_slope,
                                                       pr_enchb_cc_female = pr_enchb_cc_female,
                                                       cirrhosis_male_cofactor = cirrhosis_male_cofactor,
                                                       pr_cc_dcc = pr_cc_dcc,
                                                       cancer_prog_coefficient_female = cancer_prog_coefficient_female,
                                                       cancer_age_threshold = cancer_age_threshold,
                                                       cancer_male_cofactor = cancer_male_cofactor,
                                                       hccr_it = hccr_it,
                                                       hccr_ir = hccr_ir,
                                                       hccr_enchb = hccr_enchb,
                                                       hccr_cc = hccr_cc,
                                                       hccr_dcc = hccr_dcc,
                                                       mu_cc = mu_cc,
                                                       mu_dcc = mu_dcc,
                                                       mu_hcc = mu_hcc,
                                                       vacc_eff = vacc_eff))
  sse <- calculate_sos(temp$mapped_output)
  return(sse)
  }
}

# Starting values
starting_values1 <- unlist(parameter_list[1:32])  # default parameter list: done
starting_values1["alpha"] <- 9
starting_values1["pr_ir_cc_age_threshold"] <- 10
starting_values2 <- c(b1 = 0.13, b2 = 0.04, b3 = 0.01, mtct_prob_s = 0.05, mtct_prob_e = 0.6,
                      alpha = 7, p_chronic_in_mtct = 0.89, p_chronic_function_r = 0.65,
                      p_chronic_function_s = 0.46, pr_it_ir = 0.1, pr_ir_ic = 0.8, eag_prog_function_rate = 0,
                      pr_ir_enchb = 0.005, pr_ir_cc_female = 0.01, pr_ir_cc_age_threshold = 0, pr_ic_enchb = 0.01,
                      sag_loss_slope = 0.000451, pr_enchb_cc_female = 0.008, cirrhosis_male_cofactor = 5,
                      pr_cc_dcc = 0.04, cancer_prog_coefficient_female = 0.00022, cancer_age_threshold = 0,
                      cancer_male_cofactor = 3, hccr_it = 5, hccr_ir = 15, hccr_enchb = 10, hccr_cc = 25,
                      hccr_dcc = 0.07, mu_cc = 0.005, mu_dcc = 0.8, mu_hcc = 1.5, vacc_eff = 0.95)  # manual fit
starting_values3 <- c(b1 = 0.6, b2 = 0.01, b3 = 0.1, mtct_prob_s = 0.2, mtct_prob_e = 0.8,
                      alpha = 2, p_chronic_in_mtct = 0.7, p_chronic_function_r = 0.3,
                      p_chronic_function_s = 0.7, pr_it_ir = 0.5, pr_ir_ic = 0.1, eag_prog_function_rate = 0.01,
                      pr_ir_enchb = 0.05, pr_ir_cc_female = 0.04, pr_ir_cc_age_threshold = 10, pr_ic_enchb = 0.005,
                      sag_loss_slope = 0.00025, pr_enchb_cc_female = 0.001, cirrhosis_male_cofactor = 2,
                      pr_cc_dcc = 0.09, cancer_prog_coefficient_female = 0.0003, cancer_age_threshold = 10,
                      cancer_male_cofactor = 12, hccr_it = 2, hccr_ir = 8, hccr_enchb = 4, hccr_cc = 15,
                      hccr_dcc = 0.5, mu_cc = 0, mu_dcc = 0.05, mu_hcc = 0.05, vacc_eff = 0.6) # manual variation
starting_values4 <- c(b1 = 0.6935, b2 = 0.4969, b3 = 0.0287, mtct_prob_s = 0.1361, mtct_prob_e =0.4,
                      alpha = 2.8248, p_chronic_in_mtct = 0.9320, p_chronic_function_r =  0.7110,
                      p_chronic_function_s = 0.5471, pr_it_ir = 0.4279, pr_ir_ic = 0.0160, eag_prog_function_rate = 0.0048,
                      pr_ir_enchb = 0.1216, pr_ir_cc_female = 0.0059, pr_ir_cc_age_threshold = 11, pr_ic_enchb = 0.0251,
                      sag_loss_slope = 0.0006, pr_enchb_cc_female = 0.0301, cirrhosis_male_cofactor = 10,
                      pr_cc_dcc = 0.0574, cancer_prog_coefficient_female = 0.0002, cancer_age_threshold = 1,
                      cancer_male_cofactor = 13.7530, hccr_it = 1.9665, hccr_ir = 5.0727, hccr_enchb = 3.0000, hccr_cc = 53.6816,
                      hccr_dcc = 0.5379, mu_cc = 0.01, mu_dcc = 2.9266, mu_hcc = 2.0062, vacc_eff =  0.6943)
# picked random number between min and max here

# Also need to try different distance functions!

# Optimise with different starting values: using least squares
#res <- optim(fn=least_squares_function_all, par = starting_values1)
#res2 <- optim(fn=least_squares_function_all, par = starting_values2)

run_frequentist_fit_cluster <- function(fit_function, starting_parms) {
  res <- optim(fn=fit_function, par = starting_parms)
  return(res)
}
