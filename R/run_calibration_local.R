#################################
### Run calibration on laptop ###
#################################

### Load packages and source file ----
require(here)
require(truncnorm)

# Load functions and data required for calibration
source(here("R/imperial_model_calibration.R"))

### Run the model with each parameter set: vary parameters manually ----
# Draw parameter sets from prior distribution using Latin Hypercube Sampling
n_sims <- 1  # number of simulations
n_parms_to_vary <- 3  # number of parameters to infer - this requires manual adaptations below
lhs_samples <- randomLHS(n_sims, n_parms_to_vary) # draw 100 samples from uniform distribution U(0,1) using a Latin Hypercube design
params_mat <- data.frame(b1 = lhs_samples[,1],
                         b2 = lhs_samples[,2],
                         mtct_prob_s = lhs_samples[,3])
params_mat$b1 <- 0 + (0.2-0) * params_mat$b1 # rescale U(0,1) to be U(0,0.2)
params_mat$b2 <- 0 + (0.01-0) * params_mat$b2 # rescale U(0,1) to be U(0,0.01)
params_mat$mtct_prob_s <- 0 + (0.5-0) * params_mat$mtct_prob_s # rescale U(0,1) to be U(0,0.5)

params_mat <- data.frame(b1 = 0.13, b2 = 0.04, mtct_prob_s = 0.05)

## @knitr part2

time1 <- proc.time()
out_mat <- apply(params_mat,1,
                 function(x)
                   fit_model(default_parameter_list = parameter_list,
                             data_to_fit = calibration_datasets_list,
                             parms_to_change =
                               list(b1 = as.list(x)$b1,
                                    b2 = as.list(x)$b2,
                                    mtct_prob_s = as.list(x)$mtct_prob_s,
                                    mtct_prob_e = 0.6,  # decrease
                                    alpha = 7,
                                    b3 = 0.01,
                                    eag_prog_function_rate = 0,
                                    pr_it_ir = 0.1,
                                    pr_ir_ic = 0.8,
                                    pr_ir_cc_female = 0.028,
                                    pr_ir_cc_age_threshold = 15,
                                    pr_ir_enchb = 0.005,
                                    pr_ic_enchb = 0.01,
                                    pr_enchb_cc_female = 0.005, # 0.005, 0.016
                                    hccr_dcc = 0.07,  # 5 times increase
                                    hccr_it = 5,
                                    hccr_ir = 15,
                                    hccr_enchb = 10,
                                    hccr_cc = 25,
                                    cirrhosis_male_cofactor = 5,  # increase, 20
                                    cancer_prog_coefficient_female = 0.00022,  # doubled 0.0002
                                    cancer_age_threshold = 10,
                                    cancer_male_cofactor = 3,
                                    mu_cc = 0.005, # decrease
                                    mu_hcc = 1.5,  # increase
                                    mu_dcc = 0.8  # 1
                               )))

sim_duration = proc.time() - time1
paste(sim_duration["elapsed"], "seconds")

# Profiling
#profvis(fit_model(default_parameter_list = parameter_list,
#                      data_to_fit = calibration_datasets_list))

# Matrix of parameter values and error term
out_mat_subset <- sapply(out_mat, "[[", "error_term")
res_mat <- cbind(params_mat, error_term = out_mat_subset)
res_mat[res_mat$error_term == min(res_mat$error_term),]

## @knitr part3
### Run the model with each parameter set: vary all parameters ----
# Option 2: Draw parameter sets randomly from prior distribution
n_sims <- 100  # number of simulations/parameter sets
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
                         pr_ir_cc_age_threshold = round(runif(n_sims, 0, 15),0),
                         pr_ic_enchb = rgamma(n_sims, 3.12, 141.30),
                         sag_loss_slope = rnorm(n_sims, 0.0004106, 0.00005),
                         pr_enchb_cc_female = rgamma(n_sims, 2.3, 123.8),
                         cirrhosis_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                         pr_cc_dcc = rgamma(n_sims,17.94,423.61),
                         cancer_prog_coefficient_female = runif(n_sims, 0.0001, 0.0003),
                         cancer_age_threshold = round(runif(n_sims, 0, 15),0),
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
# Parameters to fix: hccr_ic = 1 and cancer_prog_constant_female = 0

time1 <- proc.time()
out_mat <- apply(params_mat,1,
                 function(x)
                   fit_model(default_parameter_list = parameter_list,
                             data_to_fit = calibration_datasets_list,
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
                               )))  # increase

sim_duration = proc.time() - time1
paste(sim_duration["elapsed"], "seconds")

# Profiling
#profvis(fit_model(default_parameter_list = parameter_list,
#                      data_to_fit = calibration_datasets_list))

# Matrix of parameter values and error term
out_mat_subset <- sapply(out_mat, "[[", "error_term")
res_mat <- cbind(params_mat, error_term = out_mat_subset)
res_mat[res_mat$error_term == min(res_mat$error_term),]

### Output calibration plots ----
library(grid)
library(ggplot2)
library(gridExtra)
# Loop to create plot set for every parameter combination
pdf(file = here("output/random_fit_plots", "newdata_plots.pdf"), paper="a4r")
plot_list = list()
for (i in 1:length(out_mat)) {
#for (i in c(21,30,42,44,80,99)) {
  # Parameter set table and error
  p_parms <- grid.arrange(tableGrob(lapply(out_mat[[i]]$parameter_set[1:17], function(x) round(x,6)),
                                    rows = names(out_mat[[i]]$parameter_set[1:17]),
                                    cols = "Parameters",
                                    theme = ttheme_minimal(base_size = 8)),
                          tableGrob(lapply(out_mat[[i]]$parameter_set[18:34], function(x) round(x,6)),
                                    rows = names(out_mat[[i]]$parameter_set[18:34]),
                                    cols = "Parameters (cont.)", theme=ttheme_minimal(base_size = 8)),
                          tableGrob(out_mat[[i]]$error_term, cols = "Error term"),
                          nrow = 1)

  # OUTPUTS

  ## HBsAg prevalence by time and age

  # Define study labels
  hbsag_studies <- unique(data.frame(time = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                   outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$time,
                                     paper_first_author = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                                 outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                     paper_year = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                         outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                     study_link = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                         outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$study_link))
  years_with_several_studies <- hbsag_studies[duplicated(hbsag_studies$time),1]
  hbsag_studies_double <- data.frame(time = years_with_several_studies,
                                     paper_first_author = "Several studies",
                                     paper_year = "Several studies",
                                     study_link = "Several studies")
  hbsag_studies_double$label <- c("Thursz 1995, Bellamy 1998", "Whittle 1991, Whittle 1995", "Whittle 1995, Kirk 2004")
  hbsag_studies_unique <- hbsag_studies[!(hbsag_studies$time %in% years_with_several_studies),]
  hbsag_studies_unique$label <- paste(hbsag_studies_unique$paper_first_author, hbsag_studies_unique$paper_year)
  hbsag_study_labels <- rbind(hbsag_studies_unique, hbsag_studies_double)

  # Make plot
  p_hbsag1 <- print(ggplot(data = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                         outcome == "HBsAg_prevalence")) +
                      geom_line(aes(x = age, y = model_value, linetype = "Model",
                                    colour = sex)) +
                      geom_point(aes(x = age, y = data_value,
                                     fill = "Data", colour = sex),
                                 shape = 4, stroke = 1.5) +
                      scale_linetype_manual(name = NULL, values = c("Model" = "solid")) +
                      scale_fill_manual(name = NULL, values = c("Data" = "black")) +
                      geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
                      facet_wrap(~ time, ncol = 3) +
                      geom_text(size = 3, data = hbsag_study_labels,
                                mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
                      labs(title = "HBsAg prevalence over time and by age",
                           y = "Prevalence (proportion)", x = "Age (years)",
                           colour = "Sex",
                           caption = "Keneba Manduar cohort: Whittle studies, Van der Sande 2005 | GHIS: Chotard 1992, Fortuin 1993, Wild 1993 | GLCS: Kirk 2004 | PROLIFICA: Lemoine 2016") +
                      theme_bw() +
                      theme(plot.title = element_text(hjust = 0.5),
                            plot.caption = element_text(hjust = 0, size = 6),
                            legend.margin=margin(t = 0, unit="cm")) +
                      ylim(0,0.6))

  # Carrier prevalence over time
  p_hbsag2 <- print(ggplot() +
                      geom_line(aes(x = out_mat[[i]]$full_output$time,
                                    y = apply(out_mat[[i]]$full_output$carriers,1,sum)/
                                      apply(out_mat[[i]]$full_output$pop,1,sum))) +
                      labs(title = "Modelled HBsAg prevalence over time", y = "Prevalence (proportion)", x = "Time") +
                      theme_bw() +
                      theme(plot.title = element_text(hjust = 0.5)) +
                      ylim(0,0.6))

  ## Anti-HBc prevalence by time and age

  # Define study labels
  anti_hbc_studies <- unique(data.frame(time = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                      outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$time,
                                        paper_first_author = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                                    outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                        paper_year = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                            outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                        study_link = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                            outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$study_link))
  years_with_several_studies_antihbc <- anti_hbc_studies[duplicated(anti_hbc_studies$time),1]
  anti_hbc_studies_double <- data.frame(time = years_with_several_studies_antihbc,
                                        paper_first_author = "Several studies",
                                        paper_year = "Several studies",
                                        study_link = "Several studies")
  anti_hbc_studies_double$label <- "Thursz 1995, Bellamy 1998"
  anti_hbc_studies_unique <- anti_hbc_studies[!(anti_hbc_studies$time %in% years_with_several_studies_antihbc),]
  anti_hbc_studies_unique$label <- paste(anti_hbc_studies_unique$paper_first_author, anti_hbc_studies_unique$paper_year)
  antihbc_study_labels <- rbind(anti_hbc_studies_unique, anti_hbc_studies_double)

  # Make plot
  p_antihbc <- print(ggplot(data = out_mat[[i]]$mapped_output$seromarker_prevalence[
    out_mat[[i]]$mapped_output$seromarker_prevalence$outcome == "Anti_HBc_prevalence",]) +
      geom_line(aes(x = age, y = model_value, linetype = "Model", colour = sex)) +
      geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
                 shape = 4, stroke = 1.5) +
      geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
      scale_linetype_manual(name = NULL, values = c("Model" = "solid")) +
      scale_fill_manual(name = NULL, values = c("Data" = "black")) +
      facet_wrap(~ time, ncol = 3) +
      geom_text(size = 3, data = antihbc_study_labels,
                mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
      labs(title = "Anti-HBc prevalence over time and by age",
           y = "Prevalence (proportion)", x = "Age (years)",
           colour = "Sex",
           caption = "Keneba Manduar cohort: Whittle studies") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.caption = element_text(hjust = 0, size = 6),
            legend.margin=margin(t = 0, unit="cm")) +
      ylim(0,1))

  ## HBeAg prevalence by time and age

  # Define study labels
  hbeag_studies <- unique(data.frame(time = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                   outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$time,
                                     paper_first_author = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                                 outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                     paper_year = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                         outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                     study_link = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                         outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$study_link))
  years_with_several_studies_hbeag <- unique(hbeag_studies[duplicated(hbeag_studies$time),1])
  hbeag_studies_double <- data.frame(time = years_with_several_studies_hbeag,
                                     paper_first_author = "Several studies",
                                     paper_year = "Several studies",
                                     study_link = "Several studies")
  hbeag_studies_double$label <- c("Whittle 1995, Mendy 2008", "Van der Sande 2006, Mendy 2008")
  hbeag_studies_unique <- hbeag_studies[!(hbeag_studies$time %in% years_with_several_studies_hbeag),]
  hbeag_studies_unique$label <- paste(hbeag_studies_unique$paper_first_author, hbeag_studies_unique$paper_year)
  hbeag_study_labels <- rbind(hbeag_studies_unique, hbeag_studies_double)

  # Make plot
  p_hbeag <- print(ggplot(data = out_mat[[i]]$mapped_output$seromarker_prevalence[
    out_mat[[i]]$mapped_output$seromarker_prevalence$outcome == "HBeAg_prevalence",]) +
      geom_line(aes(x = age, y = model_value, linetype = "Model", colour = sex)) +
      geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
                 shape = 4, stroke = 1.5) +
      geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
      scale_linetype_manual(name = NULL, values = c("Model" = "solid")) +
      scale_fill_manual(name = NULL, values = c("Data" = "black")) +
      facet_wrap(~ time, ncol = 3) +
      geom_text(size = 3, data = hbeag_study_labels,
                mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
      labs(title = "HBeAg prevalence over time and by age",
           y = "Prevalence (proportion)", x = "Age (years)",
           colour = "Sex",
           caption = "Keneba Manduar cohort: Whittle studies, Mendy 1999 & 2008, Van der Sande 2006, Shimakawa 2016 |\nGHIS: Chotard 1992, Fortuin 1993, Whittle 1995, Mendy 1999, Peto 2014 | GLCS: Mendy 2010 | PROLIFICA: Lemoine 2016") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.caption = element_text(hjust = 0, size = 6),
            legend.margin=margin(t = 0, unit="cm")) +
      ylim(0,1))

  ## GLOBOCAN PAF-adjusted cancer incidence and mortality in 2018
  globocan_outcome_facet_labels <- c("HCC case incidence", "HCC mortality")
  names(globocan_outcome_facet_labels) <- c("hcc_incidence", "hcc_mortality")

  p_globocan1 <- print(ggplot(data = subset(out_mat[[i]]$mapped_output$globocan_hcc_incidence,
                                            time == 2018)) +
                         geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                         geom_point(aes(x = paste(age_min,"-",age_max), y = data_value*100000, colour = "Data"),
                                    shape = 4, stroke = 1.5) +
                         geom_errorbar(aes(x = paste(age_min,"-",age_max), ymax = ci_upper*100000, ymin = ci_lower*100000),
                                       col = "red", width = 0.2) +
                         scale_fill_manual("", values = c("Model" = "gray35")) +
                         scale_colour_manual("", values = c("Data" = "red")) +
                         facet_grid(outcome ~ sex,
                                    labeller = labeller(outcome =globocan_outcome_facet_labels)) +
                         theme_bw() +
                         labs(title = "GLOBOCAN HBV-related HCC incidence and mortality rates in 2018",
                              y = "Cases/deaths per 100000 PY", x = "Age (years)",
                              subtitle = "GLOBOCAN rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
                         theme(plot.title = element_text(hjust = 0.5),
                               plot.subtitle = element_text(hjust = 0.5, size = 10),
                               legend.margin=margin(t = 0, unit="cm")) +
                         ylim(0,100))

  ## Modelled age pattern in corresponding number of HCC cases
  p_hcc_pattern1 <- print(ggplot(data = subset(out_mat[[1]]$mapped_output$globocan_hcc_incidence,
                                               outcome == "hcc_incidence" & time == 2018)) +
                            geom_col(aes(x = paste(age_min,"-",age_max), y = model_events)) +
                            facet_grid(~ sex, scales = "free") +
                            theme_bw() +
                            labs(title = "Modelled age pattern in number of incident HBV-attributable HCC cases in 2018",
                                 y = "Number of cases", x = "Age (years)") +
                            theme_bw() +
                            theme(plot.title = element_text(hjust = 0.5),
                                  plot.subtitle = element_text(hjust = 0.5, size = 10),
                                  legend.margin=margin(t = 0, unit="cm")))

  ## GLOBOCAN PAF-adjusted cancer incidence in 1988 and 1998
  p_globocan2 <- print(ggplot(data = subset(out_mat[[i]]$mapped_output$globocan_hcc_incidence,
                                            time != 2018)) +
                         geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                         #  geom_errorbar(aes(x = paste(age_min,"-",age_max), ymax = ci_upper, ymin = ci_lower)) +
                         geom_point(aes(x = paste(age_min,"-",age_max), y = data_value*100000, colour = "Data"),
                                    shape = 4, stroke = 1.5) +
                         scale_fill_manual("", values = c("Model" = "gray35")) +
                         scale_colour_manual("", values = c("Data" = "red")) +
                         facet_grid(time ~ sex) +
                         labs(title = "GLOBOCAN HBV-related HCC incidence rates in 1988 and 1998",
                              y = "Cases per 100000 PY", x = "Age (years)",
                              subtitle = "GLOBOCAN rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
                         theme_bw() +
                         theme(plot.title = element_text(hjust = 0.5),
                               axis.text.x = element_text(angle = 90),
                               plot.subtitle = element_text(hjust = 0.5, size = 10),
                               legend.margin = margin(t = 0, unit="cm")) +
                         ylim(0,100))

  ## Modelled age pattern in corresponding number of HCC cases
  p_hcc_pattern2 <- print(ggplot(data = subset(out_mat[[1]]$mapped_output$globocan_hcc_incidence,
                                               outcome == "hcc_incidence" & time != 2018)) +
                            geom_col(aes(x = paste(age_min,"-",age_max), y = model_events)) +
                            facet_grid(time ~ sex, scales = "free") +
                            theme_bw() +
                            labs(title = "Modelled age pattern in number of incident HBV-attributable HCC cases\nin 1988 and 1998",
                                 y = "Number of cases", x = "Age (years)") +
                            theme_bw() +
                            theme(plot.title = element_text(hjust = 0.5),
                                  plot.subtitle = element_text(hjust = 0.5, size = 10),
                                  axis.text.x = element_text(angle = 90),
                                  legend.margin=margin(t = 0, unit="cm")))

  ## GBD HBV-related cirrhosis mortality rate
  p_gbd <- print(ggplot(data = out_mat[[i]]$mapped_output$gbd_cirrhosis_mortality) +
                   geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                   geom_point(aes(x = paste(age_min,"-",age_max), y = data_value*100000, colour = "Data"),
                              shape = 4, stroke = 1.5) +
                   geom_errorbar(aes(x = paste(age_min,"-",age_max), ymax = ci_upper*100000, ymin = ci_lower*100000),
                                 col = "red", width = 0.5) +
                   scale_fill_manual("", values = c("Model" = "gray35")) +
                   scale_colour_manual("", values = c("Data" = "red")) +
                   facet_grid(time ~ sex) +
                   labs(title = "Global Burden of Disease Study HBV-related cirrhosis mortality rates",
                        y = "Deaths per 100000 PY", x = "Age (years)") +
                   theme_bw() +
                   theme(plot.title = element_text(hjust = 0.5),
                         axis.text.x = element_text(angle = 90),
                         legend.margin = margin(t = 0, unit="cm")) +
                   ylim(0,300))

  # Demographic characteristics of HBV-related liver disease patients
  # Proportion male
  plot_ld_prop_male <-
    ggplot(data = subset(out_mat[[i]]$mapped_output$mapped_liver_disease_demography,
                         grepl("prop_male$",outcome))) +
    geom_col(aes(x = gsub("_.*$","",outcome), y = model_value)) +
    geom_point(aes(x = gsub("_.*$","",outcome), y = data_value, colour = "Data"),
               size = 5, shape = 4, stroke = 1.5) +
    geom_errorbar(aes(x = gsub("_.*$","",outcome), ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.2) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(y = "Proportion male", x = "") +
    theme_bw() +
    theme(legend.margin = margin(t = 0, unit="cm"),
          legend.position = "bottom",
          legend.justification = "right") +
    ylim(0,1)

  # Mean age
  plot_ld_mean_age <- ggplot(data = subset(out_mat[[i]]$mapped_output$mapped_liver_disease_demography,
                                           grepl("mean_age$",outcome))) +
    geom_col(aes(x = gsub("_.*$","",outcome), y = model_value, fill = "Model")) +
    geom_point(aes(x = gsub("_.*$","",outcome), y = data_value),
               col = "red", size = 5, shape = 4, stroke = 1.5) +
    geom_errorbar(aes(x = gsub("_.*$","",outcome), ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.2) +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    labs(y = "Mean age (years)", x = "") +
    theme_bw() +
    theme(legend.margin = margin(t = 0, unit="cm"),
          legend.position = "bottom",
          legend.justification = "left") +
    ylim(0,100)

  ## Combined liver disease demography
  p_ld_demog <- grid.arrange(plot_ld_prop_male, plot_ld_mean_age, nrow = 1,
                             top = "HBV-related liver disease patients: demographic characteristics\nGambia Liver Cancer Study (Mendy 2010)")

  ## Risk of chronic carriage: change to author and year and add caption which ones are from Edmunds
  p_p_chronic <- print(ggplot(data = out_mat[[i]]$mapped_output$risk_of_chronic_carriage) +
                         geom_line(aes(x = age, y = model_value, group = "Model", linetype = "Model")) +
                         geom_point(data = subset(out_mat[[i]]$mapped_output$risk_of_chronic_carriage,
                                                  is.na(data_value) == FALSE),
                                    aes(x = age, y = data_value, colour = paste(paper_first_author, paper_year)),
                                    shape = 4, stroke = 1.5) +
                         geom_errorbar(data = subset(out_mat[[i]]$mapped_output$risk_of_chronic_carriage,
                                                     is.na(data_value) == FALSE),
                                       aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = paste(paper_first_author, paper_year))) +
                         scale_linetype_manual(name = "", values = c("Model" = "solid")) +
                         labs(title = "Risk of chronic carriage by age at infection",
                              y = "Risk (proportion)", x = "Age (years)",
                              colour = "Data",
                              caption = "Gambian studies - Keneba Manduar cohort: Whittle 1990* | GHIS: Wild 1993, Fortuin 1993\nWest African studies - Senegal: Barin 1981, Marinier 1985*, Coursaget 1987* |  Liberia: Prince 1981\n* In Edmunds 1993 review") +
                         theme_bw() +
                         theme(plot.title = element_text(hjust = 0.5),
                               plot.caption = element_text(hjust = 0, size = 6),
                               legend.title = element_text(size = 9)) +
                         ylim(0,1) +
                         xlim(0,30))

  ## Mortality curves
  # Add articifial zeros at first timestep to allow plotting of step curves
  mortality_curves_zeros <- out_mat[[i]]$mapped_output$mortality_curves
  mortality_curves_zeros$time_interval_years <- 0
  mortality_curves_zeros$data_value <- 0
  mortality_curves_zeros$model_value <- 0
  mortality_curves_zeros$number_at_risk <- mortality_curves_zeros$sample_size
  mortality_curves_zeros <- unique(mortality_curves_zeros)

  # Add labels for panels with reference and study population
  mort_curves_labels <- c("Shimakawa 2016:\ncumulative mortality in\ncomp. cirrhosis patients\n(Gambia)",
                          "Yang 2017:\ncumulative mortality in\n HCC patients\n(sS Africa, not Gambia)",
                          "Diarra 2010:\ncumulative HCC incid. in\n mixed cirrhosis patients\n(Mali)",
                          "Diarra 2010:\ncumulative mortality in\nmixed cirrhosis patients\n(Mali)",
                          "Bah 2011 (IARC):\ncumulative mortality in\nHCC patients\n(Gambia)")
  names(mort_curves_labels) <- c("shadow4_cum_mortality", "shadow5_cum_mortality",
                                 "shadow6_cum_hcc_incidence", "shadow6_cum_mortality",
                                 "shadow7_cum_mortality")

  p_mort_curves <- print(ggplot(data = rbind(out_mat[[i]]$mapped_output$mortality_curves,
                                             mortality_curves_zeros)) +
                           geom_step(aes(x = time_interval_years, y = model_value, linetype = "Model")) +
                           geom_point(aes(x = time_interval_years, y = data_value, colour = "Data"),
                                      shape = 4, stroke = 1.5) +
                           facet_grid(~ outcome, scales = "free",
                                      labeller = labeller(outcome = mort_curves_labels)) +
                           scale_colour_manual(name = "", values = c("Data" = "red")) +
                           scale_linetype_manual(name = "", values = c("Model" = "solid")) +
                           labs(title = "Cumulative probability of death/HCC over time",
                                y = "Cumulative probability", x = "Follow-up time (years)") +
                           theme_bw() +
                           theme(plot.title = element_text(hjust = 0.5),
                                 legend.position = "bottom",
                                 strip.text.x = element_text(size = 8)))

  ## ORs
  p_or <- print(ggplot(data = out_mat[[i]]$mapped_output$odds_ratios) +
                  geom_col(aes(x = gsub("odds_ratio_", "", outcome), y = log(model_value),
                               fill = "Model")) +
                  geom_point(aes(x = gsub("odds_ratio_", "", outcome), y = log(data_value),
                                 colour = "Data"),
                             shape = 4, size = 3, stroke = 2) +
                  geom_errorbar(aes(x = gsub("odds_ratio_", "", outcome),
                                    ymax = log(ci_upper), ymin = log(ci_lower)), col= "red", width = 0.2) +
                  geom_hline(aes(yintercept=0), colour = "blue") +
                  geom_text(aes(3.5, 0.25, label = ">0 positive association", vjust = 0.25, hjust = 0),
                            size = 3, colour = "blue") +
                  geom_text(aes(3.5, -0.25, label = "<0 negative association", vjust = 0.25, hjust = 0),
                            size = 3, colour = "blue") +
                  coord_cartesian(xlim = c(0.75,3.25), clip = "off") +
                  labs(title = "Log odds ratios for liver disease outcomes",
                       subtitle = "Gambia Liver Cancer Study (Mendy 2010)\nKeneba Manduar chronic carrier cohort (Shimakawa 2016)",
                       y = "ln(OR)", x = "") +
                  scale_x_discrete(breaks=c("current_hbeag_positivity_and_cirrhosis",
                                            "current_hbeag_positivity_and_hcc",
                                            "male_sex_and_significant_liver_fibrosis_or_cirrhosis"),
                                   labels=c("Odds of cirrhosis in\ncurrent HBeAg-positives\nvs.\ncurrent HBeAg-negatives",
                                            "Odds of HCC in\ncurrent HBeAg-positives vs.\ncurrent HBeAg-negatives",
                                            "Odds of\nsignificant liver fibrosis\nor cirrhosis\nin males vs. females")) +
                  scale_fill_manual("", values = c("Model" = "gray35")) +
                  scale_colour_manual("", values = c("Data" = "red")) +
                  theme_bw() +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5, size = 8),
                        axis.text.x = element_text(size = 8),
                        plot.margin = unit(c(1,3,1,1), "lines")))

  # Natural history prevalence plots
  out_mat[[i]]$mapped_output$nat_hist_prevalence$model_num <-
    gsub(".*[[:digit:]]{4}_", "",out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique)

  # GMB1 PROLIFICA plots: infection phase in chronic carriers
  gmb1_facet_labels <- c("Male blood donors", "Community screening pop.")
  names(gmb1_facet_labels) <- c("Male", "Mixed")

  plot1_gmb1 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                     id_paper == "GMB1" &
                                       model_num != "cc_dcc" & model_num != "hcc"),
                       aes(x = model_num)) +
    geom_col(aes(y = model_value))+
    geom_point(aes(y = data_value), shape = 4, size = 1.5, stroke = 2, col = "red") +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    facet_grid(~sex, scales = "free", labeller = labeller(sex = gmb1_facet_labels)) +
    scale_x_discrete(breaks=c("ic", "ir_enchb", "it_ic"),
                     labels=c("IC", "IR or\nENCHB", "IT or IC")) +
    labs(title = "Infection phase in chronic carriers",
         subtitle = "PROLIFICA (Lemoine 2016)",
         y = "Prevalence (proportion)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    ylim(0,1)

  # GMB1 plots: liver disease in chronic carriers
  plot2_gmb1 <-ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "GMB1" &
                                      (model_num == "cc_dcc" | model_num == "hcc")),
                      aes(x = gsub("_", " or ", toupper(model_num)))) +
    geom_col(aes(y = model_value))+
    geom_point(aes(y = data_value), shape = 4, size = 1.5, stroke = 2, col = "red") +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    facet_grid(~sex, scales = "free", labeller = labeller(sex = gmb1_facet_labels)) +
    labs(title = "Liver disease in chronic carriers",
         subtitle = "PROLIFICA (Lemoine 2016)",
         y = "Prevalence (proportion)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    ylim(0,0.1)

  # 1-1 plots: infection phase in chronic carriers without liver disease
  study_1_facet_labels <- c("1986: median age 11 years", "2013: median age 38 years")
  names(study_1_facet_labels) <- c(1986, 2013)

  plot1_1 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                  grepl(".*it,_ir,_ic_and_enchb$", out_mat[[i]]$mapped_output$nat_hist_prevalence$outcome)),
                    aes(x = toupper(model_num))) +
    geom_col(aes(y = model_value))+
    geom_point(aes(y = data_value), shape = 4, size = 1.5, stroke = 2, col = "red") +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    facet_grid(~time, scales = "free", labeller = labeller(time = study_1_facet_labels)) +
    labs(title = "Infection phase in chronic carriers\nwithout liver disease",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Prevalence (proportion)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    ylim(0,1)

  # 1-1 plots: liver disease in chronic carriers
  plot2_1 <-ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                 id_paper == "1" &
                                   (outcome == "hcc_prevalence_in_chronic_carriers" |
                                      outcome == "cc_and_dcc_prevalence_in_chronic_carriers")),
                   aes(x = gsub("_", " or ", toupper(model_num)))) +
    geom_col(aes(y = model_value))+
    geom_point(aes(y = data_value), shape = 4, size = 1.5, stroke = 2, col = "red") +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    facet_grid(~time, scales = "free_x", labeller = labeller(time = study_1_facet_labels)) +
    labs(title = "Liver disease in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Prevalence (proportion)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    ylim(0,0.1)

  # 1-1 plots: chronic carriers by age
  plot3_1 <- ggplot(data = out_mat[[i]]$mapped_output$nat_hist_prevalence[
    out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_2013_ir_enchb_cc_dcc",]) +
    geom_col(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = model_value, fill = "Model"))+
    geom_point(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = data_value, colour = "Data"),
               shape = 4, size = 1.5, stroke = 2) +
    geom_errorbar(aes(x = reorder(paste(age_min,"-",age_max), age_min),
                      ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Significant liver fibrosis or cirrhosis in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Prevalence (proportion)", x = "Age (years)",
         caption = "\nIT = Immune tolerant, IR = Immune reactive, IC = Inactive carrier, ENCHB = HBeAg-negative chronic hepatitis B,\nCC = Compensated cirrhosis, DCC = Decompensated cirrhosis, HCC = Hepatocellular carcinoma") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8),
          legend.margin = margin(t = 0, unit="cm"),
          legend.position = "left",
          plot.caption = element_text(hjust = 0, size = 8)) +
    ylim(0,0.75)

  ## Natural history prevalence PLOT 1
  p_nat_hist_prev1 <- grid.arrange(plot1_gmb1, plot1_1, plot2_gmb1, plot2_1,
                                   plot3_1, nrow = 3,
                                   layout_matrix = rbind(c(1,2), c(3,4), c(5,5)),
                                   top = "Prevalence measures in chronic carriers")

  # A4 Proportion of deaths from DCC and HCC in comp. cirrhosis cohort
  plot_nat_hist_a4 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                           id_paper == "A4")) +
    geom_col(aes(x = model_num, y = model_value, fill = "Model"))+
    geom_point(aes(x = model_num, y = data_value, colour = "Data"),
               shape = 4, size = 3, stroke = 2) +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Proportion of deaths from\nDCC and HCC in\ncompensated cirrhosis cohort",
         subtitle = "(Shimakawa 2016)",
         y = "Proportion", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 9),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank(),
          legend.position = "left") +
    ylim(0,1)

  # GMB2 Cirrhosis prevalence in HBsAg-positive HCC patients (GLCS)
  plot_nat_hist_gmb2 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                             id_paper == "GMB2")) +
    geom_col(aes(x = model_num, y = model_value))+
    geom_point(aes(x = model_num, y = data_value),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = model_num,
                      ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.2) +
    scale_x_discrete(breaks=c("incident_hcc_cases_from_cc",
                              "incident_hcc_cases_from_dcc"),
                     labels=c("Compensated\ncirrhosis",
                              "Decompensated\ncirrhosis")) +
    labs(title = "Cirrhosis prevalence in\nHBsAg-positive HCC patients",
         subtitle = "Gambia Liver Cancer Study (Umoh 2011)",
         y = "Proportion", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(size = 7),
          plot.subtitle = element_text(hjust = 0.5, size = 8)) +
    ylim(0,1)

  # HBeAg prevalence in liver disease patients: GMB12 and GMB15
  nat_hist_glcs_facet_labels <- c("Ryder 1992: HCC patients",
                                  "Mendy 2010 (GLCS): HCC patients",
                                  "Mendy 2010 (GLCS): cirrhosis patients")
  names(nat_hist_glcs_facet_labels) <- c("id_gmb12_1_1982_hbeag_hcc",
                                         "id_gmb15_1_1999_hbeag_hcc",
                                         "id_gmb15_2_1999_hbeag_cirrhosis")

  plot_nat_hist_glcs <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                             outcome == "hbeag_prevalence_in_hcc" |
                                               outcome == "hbeag_prevalence_in_cirrhosis")) +
    geom_col(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = model_value))+
    geom_point(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = data_value),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = reorder(paste(age_min,"-",age_max),age_min),
                      ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.2) +
    facet_grid(~id_unique, scales = "free_x", labeller = labeller(id_unique = nat_hist_glcs_facet_labels)) +
    labs(title = "HBeAg prevalence in HBsAg-positive HCC/cirrhosis patients",
         y = "Proportion", x = "Age group (years)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(size = 8),
          strip.text.x = element_text(size = 8)) +
    ylim(0,1)

  ## Natural history prevalence PLOT 2
  p_nat_hist_prev2 <- grid.arrange(plot_nat_hist_a4,  plot_nat_hist_gmb2, plot_nat_hist_glcs, nrow = 2,
                                   layout_matrix = rbind(c(1,2),
                                                         c(3,3)),
                                   top = "Prevalence measures in liver disease patients")

  # Vertical transmission plots
  # 1-1 plots: chronic infections due to vertical transmission
  plot_nat_hist_1 <- ggplot(data = out_mat[[i]]$mapped_output$nat_hist_prevalence[
    out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_1986_incident_chronic_births",]) +
    geom_col(aes(x = model_num, y = model_value))+
    geom_point(aes(x = model_num, y = data_value),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = model_num, ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.1) +
    labs(title = "Proportion of chronic infection\ndue to vertical transmission\nin unvaccinated pop.",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Proportion", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank()) +
    ylim(0,1)

  # MTCT risk
  plot_mtct <- ggplot(data = out_mat[[i]]$mapped_output$mtct_risk,
                      aes(x = paste(paper_first_author, paper_year))) +
    geom_col(aes(y = model_value, fill = "Model"))+
    geom_point(aes(y = data_value, colour = "Data"),
               shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.1) +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Mother-to-child transmission risk",
         y = "Proportion", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,1)

  ## Progression rates
  out_mat[[i]]$mapped_output$progression_rates$type <-
    gsub("shadow[[:digit:]].{0,1}_", "", out_mat[[i]]$mapped_output$progression_rates$outcome)
  out_mat[[i]]$mapped_output$progression_rates$type <- gsub("_.$", "", out_mat[[i]]$mapped_output$progression_rates$type)

  # Study 1: HCC incidence
  plot_1_hcc_incidence <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                               id_paper == "1" & type == "hcc_incidence")) +
    geom_col(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = model_value*100000))+
    geom_point(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = data_value*100000),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), ymax = ci_upper*100000, ymin = ci_lower*100000),
                  col = "red", width = 0.1)  +
    facet_grid(~sex, scales = "free") +
    labs(title = "HCC incidence rate in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Cases per 100000 PY", x = "Baseline age group (years)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8))

  # Study 1: DCC incidence
  plot_1_dcc_incidence <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                               id_paper == "1" & type == "dcc_incidence")) +
    geom_col(aes(x = outcome, y = model_value*100000))+
    geom_point(aes(x = outcome, y = data_value*100000),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = outcome, ymax = ci_upper*100000, ymin = ci_lower*100000),
                  col = "red", width = 0.1)  +
    labs(title = "DCC incidence rate in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort\n(Shimakawa 2016)",
         y = "Cases per 100000 PY", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank()) +
    ylim(0,100)


  # Study 1: Mortality
  plot_1_mortality <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                           id_paper == "1" & type == "mortality")) +
    geom_col(aes(x = sex, y = model_value*100000))+
    geom_point(aes(x = sex, y = data_value*100000),
               col = "red", shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(x = sex, ymax = ci_upper*100000, ymin = ci_lower*100000),
                  col = "red", width = 0.1)  +
    labs(title = "All-cause mortality rate in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Deaths per 100000 PY", x = "Sex") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8)) +
    ylim(0,1000)

  # Study A6: Mortality in cirrhosis cohort
  plot_a6_mortality <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                            id_paper == "A6")) +
    geom_col(aes(x = outcome, y = model_value*100, fill = "Model"))+
    geom_point(aes(x = outcome, y = data_value*100, colour = "Data"),
               shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(x = outcome, ymax = ci_upper*100, ymin = ci_lower*100),
                  col = "red", width = 0.1)  +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Mortality rate in liver disease patients",
         subtitle = "Mixed cohort of Nigerian CC, DCC and HCC patients\n(Olubuyide 1996)",
         y = "Deaths per 100 PY", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank(),
          legend.position = "bottom") +
    ylim(0,60)

  ## Combined liver disease incidence and mortality rates
  p_prog_rates1 <- grid.arrange(plot_1_hcc_incidence, plot_1_dcc_incidence,
                                plot_1_mortality, plot_a6_mortality, nrow = 2, widths = 4:3,
                                top = "Liver disease-related rates")

  # Study 1: HBeAg loss
  plot_1_eag_loss <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                          id_paper == "1" & type == "eag_loss")) +
    geom_col(aes(x = sex, y = model_value*100))+
    geom_point(aes(x = sex, y = data_value*100),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = sex, ymax = ci_upper*100, ymin = ci_lower*100),
                  col = "red", width = 0.1)  +
    labs(title = "Rate of HBeAg loss in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Cases per 100 PY", x = "Sex") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8)) +
    ylim(0,12)

  # Study 6: HBsAg loss
  plot_6_sag_loss <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                          id_paper == "6")) +
    geom_col(aes(x = outcome, y = model_value*100, fill = "Model"))+
    geom_point(aes(x = outcome, y = data_value*100, colour = "Data"),
               shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(x = outcome, ymax = ci_upper*100, ymin = ci_lower*100),
                  col = "red", width = 0.1)  +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Rate of HBsAg loss in\nchronic carrier children",
         subtitle = "Coursaget 1987 (Senegal)",
         y = "Cases per 100 PY", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank()) +
    ylim(0,5)

  ## Combined seromarker loss rates
  p_prog_rates2 <- grid.arrange(plot_1_eag_loss, plot_6_sag_loss, nrow = 1, widths = 2:1,
                                top = "Seromarker clearance rates")

  # Transmission-related data from GMB6 and GMB7
  plot_horizontal_transmission <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                                       id_paper == "GMB6" | id_paper == "GMB7")) +
    geom_col(aes(x = gsub(".*_","",outcome), y = model_value))+
    geom_point(aes(x = gsub(".*_","",outcome), y = data_value),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = gsub(".*_","",outcome), ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.1)  +
    scale_x_discrete(breaks=c("foi",
                              "incidence"),
                     labels=c("Force of infection\nin children in\nKeneba and Manduar\n(Whittle 1990)",
                              "Chronic infection incidence\nin susceptible children\n(Ryder 1984)")) +
    labs(title = "Horizontal transmission-related rates",
         y = "Rate (per PY)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,2)


  ## Combined transmission-related plot
  p_transmission_rates <- grid.arrange(plot_mtct, plot_nat_hist_1, plot_horizontal_transmission,
                                       layout_matrix = rbind(c(1,1),
                                                             c(2,3)),
                                       top = "Transmission-related measures")

  # List of all plots
  plot_list[[i]] <- list(p_parms, p_hbsag1, p_hbsag2, p_antihbc, p_hbeag,
                         p_globocan1, p_globocan2, p_hcc_pattern1, p_hcc_pattern2,
                         p_gbd, p_ld_demog,
                         p_p_chronic, p_mort_curves, p_or,
                         p_nat_hist_prev1, p_nat_hist_prev2,
                         p_prog_rates1, p_prog_rates2, p_transmission_rates)

}
dev.off()

## @knitr part5
### ABC rejection algorithm ----
# Steps/principles of ABC rejection for parameter estimation (Toni):
# 1) Define prior distribution
# 2) Sample a parameter - a particle - from the prior distribition
# 3) Use the deterministic model to simulate an output - the trajectory
# 4) Add noise to the simulated outputs/trajectory to be mapped to the data
# 5) Define a distance function which compares the simulated output to the
#    experimental data
# 6) Define a tolerance level epsilon, representing the desired level of agreement
#    between the data and the model output
# 7) If the distance function is < epsilon, accept the particle sampled in 2).
#    If the distance function is > epsilon, reject the particle sampled in 2)
# 8) Define the desired number of particles/parameter sets N
# 9) Repeat steps 1-7) until N particles have been accepted.
#    The accepted particles represent an approximation of the posterior distribution.

# ABC rejection alogorithm from LSHTM fitting course
run_rejection_algorithm2 <- function(N, epsilon) {


  # Set up empty matrix to store results
  results <- data.frame(b1 = 0, b2 = 0, mtct_prob_s = 0, error = 0)

  # Initialise the loop with i=0
  i <- 0

  # while the length of the accepted values (result) is less than the desired length (N)
  while(i<N) {

    # A) DRAW PARAMETERS TO BE VARIED FROM THE PRIOR DISTRIBUTION

    # Draw a new theta from prior distributions
    b1_sample <- runif(1, min = 0.01, max = 0.4)
    b2_sample <- runif(1, min = 0, max = b1)
    mtct_prob_s_sample <- runif(1, min = 0, max = 0.3)

    # Combine parameters to be varied
    params_mat <- data.frame(b1 = b1_sample, b2 = b2_sample, mtct_prob_s = mtct_prob_s_sample)

    # B) RUN THE MODEL WITH THE GIVEN PARAMETER SET,
    # CALCUTE AND PERTURB OUTPUT, ANC CALCULATE DISTANCE TO DATA
    # Noise still needs to be added
    out_mat <- apply(params_mat,1,
                     function(x)
                       fit_model(default_parameter_list = parameter_list,
                                 data_to_fit = calibration_datasets_list,
                                 parms_to_change =
                                   list(b1 = as.list(x)$b1,
                                        b2 = as.list(x)$b2,
                                        mtct_prob_s = as.list(x)$mtct_prob_s,
                                        mtct_prob_e = 0.6,  # decrease
                                        alpha = 7,
                                        b3 = 0.01,
                                        eag_prog_function_intercept = 0.1,
                                        eag_prog_function_rate = 0,
                                        pr_it_ir = 1,  # fix
                                        pr_ir_ic = 8,
                                        pr_ir_cc_female = 0.1,
                                        pr_ir_cc_age_threshold = 30,  # increase
                                        pr_ir_enchb = 0.005,
                                        pr_enchb_cc_female = 0.005, # 0.005, 0.016
                                        hccr_dcc = 0.2,  # 5 times increase
                                        hccr_ir = 16,  # doubled
                                        hccr_enchb = 6,
                                        hccr_cc = 25,
                                        cirrhosis_male_cofactor = 5,  # increase, 20
                                        cancer_prog_coefficient_female = 0.00017,  # doubled 0.0002
                                        cancer_age_threshold = 15,
                                        cancer_male_cofactor = 5,
                                        mu_cc = 0.005, # decrease
                                        mu_hcc = 1.5,  # increase
                                        mu_dcc = 0.8  # 1
                                   )))  # increase

    dist <- out_mat[[1]]$error_term

    params_mat <- data.frame(params_mat, error = dist)

    # C) COMPARE THE DISTANCE TO THE EPSILON TOLERANCE
    # If the distance is within the epsilon window,
    # accept and store the parameter values.

    if(dist<=epsilon){

      results <- rbind(results,params_mat)

    }

    # Update i (dimension of results store)
    i <- dim(results)[1]-1

    # D) REPEAT PROCEDURE UNTIL N PARAMETER SETS HAVE BEEN ACCEPTED

  }
  # return the accepted values
  return(results[-1,])
}
# Run the algorithm
res <- run_rejection_algorithm2(N=10, epsilon = 300)

# Issue with this is that it uses a while loop and defines the required number of accepted particles
# This cannot be parallelised
# However it is possible to calculate with all the parameter sets before accepting/rejecting
# Only issue is need to figure out how many runs to make to get enough accepted parameter sets
run_rejection_algorithm <- function(n) {

  # A) RANDOMLY DRAW n PARAMETERS TO BE VARIED FROM THE PRIOR DISTRIBUTION

  # Draw a new theta from prior distributions
  b1_sample <- runif(n, min = 0.01, max = 0.4)
  b2_sample <- runif(n, min = 0, max = b1_sample)
  mtct_prob_s_sample <- runif(n, min = 0, max = 0.3)

  # Combine parameters to be varied
  params_mat <- data.frame(b1 = b1_sample, b2 = b2_sample, mtct_prob_s = mtct_prob_s_sample)

  # B) RUN THE MODEL FOR ALL PARAMETER SETS,
  # CALCUTE AND PERTURB OUTPUT, ANC CALCULATE DISTANCE TO DATA
  # Noise still needs to be added
  out_mat <- apply(params_mat,1,
                   function(x)
                     fit_model(default_parameter_list = parameter_list,
                               data_to_fit = calibration_datasets_list,
                               parms_to_change =
                                 list(b1 = as.list(x)$b1,
                                      b2 = as.list(x)$b2,
                                      mtct_prob_s = as.list(x)$mtct_prob_s,
                                      mtct_prob_e = 0.6,  # decrease
                                      alpha = 7,
                                      b3 = 0.01,
                                      eag_prog_function_intercept = 0.1,
                                      eag_prog_function_rate = 0,
                                      pr_it_ir = 1,  # fix
                                      pr_ir_ic = 8,
                                      pr_ir_cc_female = 0.1,
                                      pr_ir_cc_age_threshold = 30,  # increase
                                      pr_ir_enchb = 0.005,
                                      pr_enchb_cc_female = 0.005, # 0.005, 0.016
                                      hccr_dcc = 0.2,  # 5 times increase
                                      hccr_ir = 16,  # doubled
                                      hccr_enchb = 6,
                                      hccr_cc = 25,
                                      cirrhosis_male_cofactor = 5,  # increase, 20
                                      cancer_prog_coefficient_female = 0.00017,  # doubled 0.0002
                                      cancer_age_threshold = 15,
                                      cancer_male_cofactor = 5,
                                      mu_cc = 0.005, # decrease
                                      mu_hcc = 1.5,  # increase
                                      mu_dcc = 0.8  # 1
                                 )))  # increase

  # C) COMBINE INTO OUTPUT TABLE
  out_mat_subset <- sapply(out_mat, "[[", "error_term")
  res_mat <- cbind(params_mat, error_term = out_mat_subset)

  return(res_mat)

}
res <- run_rejection_algorithm(n=10)


### Parallelised code: vary all parameters ----
# Set up cluster
cl <- makeCluster(4)
clusterEvalQ(cl, {library(dplyr); library(tidyr); library(deSolve)})
clusterExport(cl, ls())

time1 <- proc.time()
out_mat <- parApply(cl = cl, params_mat,1,
                    function(x) fit_model(default_parameter_list = parameter_list,
                                          data_to_fit = calibration_datasets_list,
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
sim_duration = proc.time() - time1
sim_duration["elapsed"]/60

# Important: stop cluster!!
stopCluster(cl)
# 10 sims take 3.33 min, 100 sims take 35 min

# Matrix of parameter values and error term
out_mat_subset <- sapply(out_mat, "[[", "error_term")
res_mat <- cbind(params_mat, error_term = out_mat_subset)
res_mat[res_mat$error_term == min(res_mat$error_term),]

### Target fitting approach ----

# Fit to overall prevalence for minimum SSE
plot(x = seq(1960,2019.5, by = 0.5),
     y = out_mat[[which(res_mat$sse == min(res_mat$sse))]]$carrier_prev_total,
     type = "l", ylim = c(0,0.5))
points(x = hbsag_dataset$time,
       y = hbsag_dataset$prev,
       col = "red")

# Fit to age-specific prevalence for minimum SSE
#plot(x = ages,
#     y = out_mat[[which(res_mat$sse == min(res_mat$sse))]]$prev_by_age_1980,
#     type = "l", ylim = c(0,0.5))
#points(x = hbsag_by_age_dataset_1980$age,
#       y = hbsag_by_age_dataset_1980$prev,
#       col = "red")


# Target fitting approach / filtration
res_mat$fit <- 0

res_mat[(res_mat$prev_est_1980 >= hbsag_dataset$ci_lower[1]) &
          (res_mat$prev_est_1980 <= hbsag_dataset$ci_upper[1]) &
          (res_mat$prev_est_2015 >= hbsag_dataset$ci_lower[2]) &
          (res_mat$prev_est_2015 <= hbsag_dataset$ci_upper[2]),]$fit <- 1
table(res_mat$fit)

# Box plot of parameter estimates
boxplot(subset(res_mat,fit==1)$b1,subset(res_mat,fit==1)$b2,subset(res_mat,fit==1)$mtct_prob_s,
        names = c("b1", "b2", "mtct_prob_s"),ylim=c(0,0.5))

# Box plot of priors and posteriors
par(mfrow=c(1,3))
boxplot(res_mat$b1, subset(res_mat,fit==1)$b1, col=c(grey(0.6),2),ylim=c(0,0.6),
        names=c("Prior","Posterior"), main="LHS for 'b1':\nprior/posteriors distributions")

boxplot(res_mat$b2, subset(res_mat,fit==1)$b2, col=c(grey(0.6),2),ylim=c(0,0.6),
        names=c("Prior","Posterior"), main="LHS for 'b2':\nprior/posteriors distributions")

boxplot(res_mat$mtct_prob_s, subset(res_mat,fit==1)$mtct_prob_s, col=c(grey(0.6),2),
        ylim=c(0,0.6), names=c("Prior","Posterior"),
        main="LHS for 'mtct_prob_s':\nprior/posteriors distributions")
par(mfrow=c(1,1))








### Save a training dataset with known parameter values ----

# Save the training dataset
mapped_output <- unlist(lapply(out_mat, "[[", "mapped_output"), recursive = FALSE)
mapped_output2 <- lapply(mapped_output, function(x) {names(x)[names(x) == "model_value"] <- "training_data_value"; x })

training_hbsag_prevalence <- left_join(calibration_datasets_list$hbsag_prevalence,
                                       mapped_output2$seromarker_prevalence[mapped_output2$seromarker_prevalence$outcome == "HBsAg_prevalence",])

training_antihbc_prevalence <- left_join(calibration_datasets_list$antihbc_prevalence,
                                         mapped_output2$seromarker_prevalence[mapped_output2$seromarker_prevalence$outcome == "Anti_HBc_prevalence",])

training_hbeag_prevalence <- left_join(calibration_datasets_list$hbeag_prevalence,
                                       mapped_output2$seromarker_prevalence[mapped_output2$seromarker_prevalence$outcome == "HBeAg_prevalence",])

training_natural_history_prevalence <- left_join(calibration_datasets_list$natural_history_prevalence,
                                                 mapped_output2$nat_hist_prevalence)

training_mtct_risk <- left_join(calibration_datasets_list$mtct_risk,
                                mapped_output2$mtct_risk)

training_progression_rates <-  left_join(calibration_datasets_list$progression_rates,
                                         mapped_output2$progression_rates)

training_mortality_curves <- left_join(calibration_datasets_list$mortality_curves,
                                       mapped_output2$mortality_curves)

training_globocan_incidence_data <- left_join(calibration_datasets_list$globocan_incidence_data,
                                              mapped_output2$globocan_hcc_incidence)
training_globocan_incidence_data$model_events <- NULL

training_globocan_mortality_curve <- left_join(calibration_datasets_list$globocan_mortality_curve,
                                               mapped_output2$mortality_curves)
training_globocan_mortality_curve <- select(training_globocan_mortality_curve, c(names(calibration_datasets_list$globocan_mortality_curve), "training_data_value"))

training_gbd_cirrhosis_mortality <- left_join(calibration_datasets_list$gbd_cirrhosis_mortality,
                                              mapped_output2$gbd_cirrhosis_mortality)
training_gbd_cirrhosis_mortality$model_events <- NULL

training_p_chronic <- left_join(calibration_datasets_list$p_chronic,
                                mapped_output2$risk_of_chronic_carriage)

training_odds_ratios <- left_join(calibration_datasets_list$odds_ratios,
                                  mapped_output2$odds_ratios)

training_liver_disease_demography <- left_join(calibration_datasets_list$liver_disease_demography,
                                               mapped_output2$mapped_liver_disease_demography)

# Checks
training_hbsag_prevalence[,-ncol(training_hbsag_prevalence)] == calibration_datasets_list$hbsag_prevalence
dim(training_hbsag_prevalence)
dim(calibration_datasets_list$hbsag_prevalence)

training_antihbc_prevalence[,-ncol(training_antihbc_prevalence)] == calibration_datasets_list$antihbc_prevalence
dim(training_antihbc_prevalence)
dim(calibration_datasets_list$antihbc_prevalence)

training_hbeag_prevalence[,-ncol(training_hbeag_prevalence)] == calibration_datasets_list$hbeag_prevalence
dim(training_hbeag_prevalence)
dim(calibration_datasets_list$hbeag_prevalence)

training_natural_history_prevalence[,-ncol(training_natural_history_prevalence)] == calibration_datasets_list$natural_history_prevalence
dim(training_natural_history_prevalence)
dim(calibration_datasets_list$natural_history_prevalence)

training_mtct_risk[,-ncol(training_mtct_risk)] == calibration_datasets_list$mtct_risk

training_progression_rates[,-ncol(training_progression_rates)] == calibration_datasets_list$progression_rates
dim(training_progression_rates)
dim(calibration_datasets_list$progression_rates)

training_mortality_curves[,-ncol(training_mortality_curves)] == calibration_datasets_list$mortality_curves
dim(training_mortality_curves)
dim(calibration_datasets_list$mortality_curves)

training_globocan_incidence_data[,-ncol(training_globocan_incidence_data)] == calibration_datasets_list$globocan_incidence_data
dim(training_globocan_incidence_data)
dim(calibration_datasets_list$globocan_incidence_data)

training_globocan_mortality_curve[,-ncol(training_globocan_mortality_curve)] ==
  calibration_datasets_list$globocan_mortality_curve
dim(training_globocan_mortality_curve)
dim(calibration_datasets_list$globocan_mortality_curve)

training_gbd_cirrhosis_mortality[,-ncol(training_gbd_cirrhosis_mortality)] ==
  calibration_datasets_list$gbd_cirrhosis_mortality
dim(training_gbd_cirrhosis_mortality)
dim(calibration_datasets_list$gbd_cirrhosis_mortality)

training_p_chronic[,-ncol(training_p_chronic)] ==
  calibration_datasets_list$p_chronic

training_odds_ratios[,-ncol(training_odds_ratios)] ==
  calibration_datasets_list$odds_ratios

training_liver_disease_demography[,-ncol(training_liver_disease_demography)] ==
  calibration_datasets_list$liver_disease_demography

# Combine into list

training_datasets_list <- list(hbsag_prevalence = training_hbsag_prevalence,
                               antihbc_prevalence = training_antihbc_prevalence,
                               hbeag_prevalence = training_hbeag_prevalence,
                               natural_history_prevalence = training_natural_history_prevalence,
                               mtct_risk = training_mtct_risk,
                               progression_rates = training_progression_rates,
                               mortality_curves = training_mortality_curves,
                               globocan_incidence_data = training_globocan_incidence_data,
                               globocan_mortality_curve = training_globocan_mortality_curve,
                               gbd_cirrhosis_mortality = training_gbd_cirrhosis_mortality,
                               p_chronic = training_p_chronic,
                               odds_ratios = training_odds_ratios,
                               liver_disease_demography = training_liver_disease_demography)

#save(training_datasets_list, file = here("output", "distance_metric_trials", "training_datasets_list.Rdata"))

### Analyse cluster fits ----
#load(here("output", "fits", "cluster_fit_150719.Rdata"))
load(here("output", "fits", "cluster_fit_with_prev_check_170719.Rdata")) # loads out

out_mat <- out

# Dataframe of parameter values and error term
out_mat_subset <- sapply(out_mat, "[[", "error_term")
res_mat <- data.frame(t(sapply(out_mat, "[[", "parameter_set"))[,1:32], error_term = out_mat_subset)

# Distribution of error
quantile(res_mat$error_term, na.rm = TRUE)

best_fit_ids <- which(res_mat$error_term < 1000)
medium_fit_ids <- which(res_mat$error_term >= 1000 & res_mat$error_term < 2000)
worst_fit_ids <- which(res_mat$error_term > 6000)

# Packages for making plots
require(ggplot2)  # for calibration plots
require(gridExtra)  # for calibration plots
require(grid)  # for calibration plots

pdf(file = here("output/fits", "add_prev_check_170719_mediumfits_transmission_weights5.pdf"), paper="a4r")
plot_list = list()
#for (i in 1:length(out_mat)) {
for (i in medium_fit_ids) {
  # Parameter set table and error
  p_parms <- grid.arrange(tableGrob(lapply(out_mat[[i]]$parameter_set[1:17], function(x) round(x,6)),
                                    rows = names(out_mat[[i]]$parameter_set[1:17]),
                                    cols = "Parameters",
                                    theme = ttheme_minimal(base_size = 8)),
                          tableGrob(lapply(out_mat[[i]]$parameter_set[18:34], function(x) round(x,6)),
                                    rows = names(out_mat[[i]]$parameter_set[18:34]),
                                    cols = "Parameters (cont.)", theme=ttheme_minimal(base_size = 8)),
                          tableGrob(out_mat[[i]]$error_term, cols = "Error term"),
                          nrow = 1)

  # OUTPUTS

  ## HBsAg prevalence by time and age

  # Define study labels
  hbsag_studies <- unique(data.frame(time = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                   outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$time,
                                     paper_first_author = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                                 outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                     paper_year = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                         outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                     study_link = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                         outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$study_link))
  years_with_several_studies <- hbsag_studies[duplicated(hbsag_studies$time),1]
  hbsag_studies_double <- data.frame(time = years_with_several_studies,
                                     paper_first_author = "Several studies",
                                     paper_year = "Several studies",
                                     study_link = "Several studies")
  hbsag_studies_double$label <- c("Thursz 1995, Bellamy 1998", "Whittle 1991, Whittle 1995", "Whittle 1995, Kirk 2004")
  hbsag_studies_unique <- hbsag_studies[!(hbsag_studies$time %in% years_with_several_studies),]
  hbsag_studies_unique$label <- paste(hbsag_studies_unique$paper_first_author, hbsag_studies_unique$paper_year)
  hbsag_study_labels <- rbind(hbsag_studies_unique, hbsag_studies_double)

  # Make plot
  p_hbsag1 <- print(ggplot(data = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                         outcome == "HBsAg_prevalence")) +
                      geom_line(aes(x = age, y = model_value, linetype = "Model",
                                    colour = sex)) +
                      geom_point(aes(x = age, y = data_value,
                                     fill = "Data", colour = sex),
                                 shape = 4, stroke = 1.5) +
                      scale_linetype_manual(name = NULL, values = c("Model" = "solid")) +
                      scale_fill_manual(name = NULL, values = c("Data" = "black")) +
                      geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
                      facet_wrap(~ time, ncol = 3) +
                      geom_text(size = 3, data = hbsag_study_labels,
                                mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
                      labs(title = "HBsAg prevalence over time and by age",
                           y = "Prevalence (proportion)", x = "Age (years)",
                           colour = "Sex",
                           caption = "Keneba Manduar cohort: Whittle studies, Van der Sande 2005 | GHIS: Chotard 1992, Fortuin 1993, Wild 1993 | GLCS: Kirk 2004 | PROLIFICA: Lemoine 2016") +
                      theme_bw() +
                      theme(plot.title = element_text(hjust = 0.5),
                            plot.caption = element_text(hjust = 0, size = 6),
                            legend.margin=margin(t = 0, unit="cm")) +
                      ylim(0,0.6))

  # Carrier prevalence over time
  p_hbsag2 <- print(ggplot() +
                      geom_line(aes(x = out_mat[[i]]$full_output$time,
                                    y = apply(out_mat[[i]]$full_output$carriers,1,sum)/
                                      apply(out_mat[[i]]$full_output$pop,1,sum))) +
                      labs(title = "Modelled HBsAg prevalence over time", y = "Prevalence (proportion)", x = "Time") +
                      theme_bw() +
                      theme(plot.title = element_text(hjust = 0.5)) +
                      ylim(0,0.6))

  ## Anti-HBc prevalence by time and age

  # Define study labels
  anti_hbc_studies <- unique(data.frame(time = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                      outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$time,
                                        paper_first_author = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                                    outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                        paper_year = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                            outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                        study_link = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                            outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$study_link))
  years_with_several_studies_antihbc <- anti_hbc_studies[duplicated(anti_hbc_studies$time),1]
  anti_hbc_studies_double <- data.frame(time = years_with_several_studies_antihbc,
                                        paper_first_author = "Several studies",
                                        paper_year = "Several studies",
                                        study_link = "Several studies")
  anti_hbc_studies_double$label <- "Thursz 1995, Bellamy 1998"
  anti_hbc_studies_unique <- anti_hbc_studies[!(anti_hbc_studies$time %in% years_with_several_studies_antihbc),]
  anti_hbc_studies_unique$label <- paste(anti_hbc_studies_unique$paper_first_author, anti_hbc_studies_unique$paper_year)
  antihbc_study_labels <- rbind(anti_hbc_studies_unique, anti_hbc_studies_double)

  # Make plot
  p_antihbc <- print(ggplot(data = out_mat[[i]]$mapped_output$seromarker_prevalence[
    out_mat[[i]]$mapped_output$seromarker_prevalence$outcome == "Anti_HBc_prevalence",]) +
      geom_line(aes(x = age, y = model_value, linetype = "Model", colour = sex)) +
      geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
                 shape = 4, stroke = 1.5) +
      geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
      scale_linetype_manual(name = NULL, values = c("Model" = "solid")) +
      scale_fill_manual(name = NULL, values = c("Data" = "black")) +
      facet_wrap(~ time, ncol = 3) +
      geom_text(size = 3, data = antihbc_study_labels,
                mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
      labs(title = "Anti-HBc prevalence over time and by age",
           y = "Prevalence (proportion)", x = "Age (years)",
           colour = "Sex",
           caption = "Keneba Manduar cohort: Whittle studies") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.caption = element_text(hjust = 0, size = 6),
            legend.margin=margin(t = 0, unit="cm")) +
      ylim(0,1))

  ## HBeAg prevalence by time and age

  # Define study labels
  hbeag_studies <- unique(data.frame(time = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                   outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$time,
                                     paper_first_author = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                                 outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                     paper_year = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                         outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                     study_link = subset(out_mat[[i]]$mapped_output$seromarker_prevalence,
                                                         outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$study_link))
  years_with_several_studies_hbeag <- unique(hbeag_studies[duplicated(hbeag_studies$time),1])
  hbeag_studies_double <- data.frame(time = years_with_several_studies_hbeag,
                                     paper_first_author = "Several studies",
                                     paper_year = "Several studies",
                                     study_link = "Several studies")
  hbeag_studies_double$label <- c("Whittle 1995, Mendy 2008", "Van der Sande 2006, Mendy 2008")
  hbeag_studies_unique <- hbeag_studies[!(hbeag_studies$time %in% years_with_several_studies_hbeag),]
  hbeag_studies_unique$label <- paste(hbeag_studies_unique$paper_first_author, hbeag_studies_unique$paper_year)
  hbeag_study_labels <- rbind(hbeag_studies_unique, hbeag_studies_double)

  # Make plot
  p_hbeag <- print(ggplot(data = out_mat[[i]]$mapped_output$seromarker_prevalence[
    out_mat[[i]]$mapped_output$seromarker_prevalence$outcome == "HBeAg_prevalence",]) +
      geom_line(aes(x = age, y = model_value, linetype = "Model", colour = sex)) +
      geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
                 shape = 4, stroke = 1.5) +
      geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
      scale_linetype_manual(name = NULL, values = c("Model" = "solid")) +
      scale_fill_manual(name = NULL, values = c("Data" = "black")) +
      facet_wrap(~ time, ncol = 3) +
      geom_text(size = 3, data = hbeag_study_labels,
                mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
      labs(title = "HBeAg prevalence over time and by age",
           y = "Prevalence (proportion)", x = "Age (years)",
           colour = "Sex",
           caption = "Keneba Manduar cohort: Whittle studies, Mendy 1999 & 2008, Van der Sande 2006, Shimakawa 2016 |\nGHIS: Chotard 1992, Fortuin 1993, Whittle 1995, Mendy 1999, Peto 2014 | GLCS: Mendy 2010 | PROLIFICA: Lemoine 2016") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.caption = element_text(hjust = 0, size = 6),
            legend.margin=margin(t = 0, unit="cm")) +
      ylim(0,1))

  ## GLOBOCAN PAF-adjusted cancer incidence and mortality in 2018
  globocan_outcome_facet_labels <- c("HCC case incidence", "HCC mortality")
  names(globocan_outcome_facet_labels) <- c("hcc_incidence", "hcc_mortality")

  p_globocan1 <- print(ggplot(data = subset(out_mat[[i]]$mapped_output$globocan_hcc_incidence,
                                            time == 2018)) +
                         geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                         geom_point(aes(x = paste(age_min,"-",age_max), y = data_value*100000, colour = "Data"),
                                    shape = 4, stroke = 1.5) +
                         geom_errorbar(aes(x = paste(age_min,"-",age_max), ymax = ci_upper*100000, ymin = ci_lower*100000),
                                       col = "red", width = 0.2) +
                         scale_fill_manual("", values = c("Model" = "gray35")) +
                         scale_colour_manual("", values = c("Data" = "red")) +
                         facet_grid(outcome ~ sex,
                                    labeller = labeller(outcome =globocan_outcome_facet_labels)) +
                         theme_bw() +
                         labs(title = "GLOBOCAN HBV-related HCC incidence and mortality rates in 2018",
                              y = "Cases/deaths per 100000 PY", x = "Age (years)",
                              subtitle = "GLOBOCAN rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
                         theme(plot.title = element_text(hjust = 0.5),
                               plot.subtitle = element_text(hjust = 0.5, size = 10),
                               legend.margin=margin(t = 0, unit="cm")) +
                         ylim(0,100))

  ## Modelled age pattern in corresponding number of HCC cases
  p_hcc_pattern1 <- print(ggplot(data = subset(out_mat[[1]]$mapped_output$globocan_hcc_incidence,
                                               outcome == "hcc_incidence" & time == 2018)) +
                            geom_col(aes(x = paste(age_min,"-",age_max), y = model_events)) +
                            facet_grid(~ sex, scales = "free") +
                            theme_bw() +
                            labs(title = "Modelled age pattern in number of incident HBV-attributable HCC cases in 2018",
                                 y = "Number of cases", x = "Age (years)") +
                            theme_bw() +
                            theme(plot.title = element_text(hjust = 0.5),
                                  plot.subtitle = element_text(hjust = 0.5, size = 10),
                                  legend.margin=margin(t = 0, unit="cm")))

  ## GLOBOCAN PAF-adjusted cancer incidence in 1988 and 1998
  p_globocan2 <- print(ggplot(data = subset(out_mat[[i]]$mapped_output$globocan_hcc_incidence,
                                            time != 2018)) +
                         geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                         #  geom_errorbar(aes(x = paste(age_min,"-",age_max), ymax = ci_upper, ymin = ci_lower)) +
                         geom_point(aes(x = paste(age_min,"-",age_max), y = data_value*100000, colour = "Data"),
                                    shape = 4, stroke = 1.5) +
                         scale_fill_manual("", values = c("Model" = "gray35")) +
                         scale_colour_manual("", values = c("Data" = "red")) +
                         facet_grid(time ~ sex) +
                         labs(title = "GLOBOCAN HBV-related HCC incidence rates in 1988 and 1998",
                              y = "Cases per 100000 PY", x = "Age (years)",
                              subtitle = "GLOBOCAN rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
                         theme_bw() +
                         theme(plot.title = element_text(hjust = 0.5),
                               axis.text.x = element_text(angle = 90),
                               plot.subtitle = element_text(hjust = 0.5, size = 10),
                               legend.margin = margin(t = 0, unit="cm")) +
                         ylim(0,100))

  ## Modelled age pattern in corresponding number of HCC cases
  p_hcc_pattern2 <- print(ggplot(data = subset(out_mat[[1]]$mapped_output$globocan_hcc_incidence,
                                               outcome == "hcc_incidence" & time != 2018)) +
                            geom_col(aes(x = paste(age_min,"-",age_max), y = model_events)) +
                            facet_grid(time ~ sex, scales = "free") +
                            theme_bw() +
                            labs(title = "Modelled age pattern in number of incident HBV-attributable HCC cases\nin 1988 and 1998",
                                 y = "Number of cases", x = "Age (years)") +
                            theme_bw() +
                            theme(plot.title = element_text(hjust = 0.5),
                                  plot.subtitle = element_text(hjust = 0.5, size = 10),
                                  axis.text.x = element_text(angle = 90),
                                  legend.margin=margin(t = 0, unit="cm")))

  ## GBD HBV-related cirrhosis mortality rate
  p_gbd <- print(ggplot(data = out_mat[[i]]$mapped_output$gbd_cirrhosis_mortality) +
                   geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                   geom_point(aes(x = paste(age_min,"-",age_max), y = data_value*100000, colour = "Data"),
                              shape = 4, stroke = 1.5) +
                   geom_errorbar(aes(x = paste(age_min,"-",age_max), ymax = ci_upper*100000, ymin = ci_lower*100000),
                                 col = "red", width = 0.5) +
                   scale_fill_manual("", values = c("Model" = "gray35")) +
                   scale_colour_manual("", values = c("Data" = "red")) +
                   facet_grid(time ~ sex) +
                   labs(title = "Global Burden of Disease Study HBV-related cirrhosis mortality rates",
                        y = "Deaths per 100000 PY", x = "Age (years)") +
                   theme_bw() +
                   theme(plot.title = element_text(hjust = 0.5),
                         axis.text.x = element_text(angle = 90),
                         legend.margin = margin(t = 0, unit="cm")) +
                   ylim(0,300))

  # Demographic characteristics of HBV-related liver disease patients
  # Proportion male
  plot_ld_prop_male <-
    ggplot(data = subset(out_mat[[i]]$mapped_output$mapped_liver_disease_demography,
                         grepl("prop_male$",outcome))) +
    geom_col(aes(x = gsub("_.*$","",outcome), y = model_value)) +
    geom_point(aes(x = gsub("_.*$","",outcome), y = data_value, colour = "Data"),
               size = 5, shape = 4, stroke = 1.5) +
    geom_errorbar(aes(x = gsub("_.*$","",outcome), ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.2) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(y = "Proportion male", x = "") +
    theme_bw() +
    theme(legend.margin = margin(t = 0, unit="cm"),
          legend.position = "bottom",
          legend.justification = "right") +
    ylim(0,1)

  # Mean age
  plot_ld_mean_age <- ggplot(data = subset(out_mat[[i]]$mapped_output$mapped_liver_disease_demography,
                                           grepl("mean_age$",outcome))) +
    geom_col(aes(x = gsub("_.*$","",outcome), y = model_value, fill = "Model")) +
    geom_point(aes(x = gsub("_.*$","",outcome), y = data_value),
               col = "red", size = 5, shape = 4, stroke = 1.5) +
    geom_errorbar(aes(x = gsub("_.*$","",outcome), ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.2) +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    labs(y = "Mean age (years)", x = "") +
    theme_bw() +
    theme(legend.margin = margin(t = 0, unit="cm"),
          legend.position = "bottom",
          legend.justification = "left") +
    ylim(0,100)

  ## Combined liver disease demography
  p_ld_demog <- grid.arrange(plot_ld_prop_male, plot_ld_mean_age, nrow = 1,
                             top = "HBV-related liver disease patients: demographic characteristics\nGambia Liver Cancer Study (Mendy 2010)")

  ## Risk of chronic carriage: change to author and year and add caption which ones are from Edmunds
  p_p_chronic <- print(ggplot(data = out_mat[[i]]$mapped_output$risk_of_chronic_carriage) +
                         geom_line(aes(x = age, y = model_value, group = "Model", linetype = "Model")) +
                         geom_point(data = subset(out_mat[[i]]$mapped_output$risk_of_chronic_carriage,
                                                  is.na(data_value) == FALSE),
                                    aes(x = age, y = data_value, colour = paste(paper_first_author, paper_year)),
                                    shape = 4, stroke = 1.5) +
                         geom_errorbar(data = subset(out_mat[[i]]$mapped_output$risk_of_chronic_carriage,
                                                     is.na(data_value) == FALSE),
                                       aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = paste(paper_first_author, paper_year))) +
                         scale_linetype_manual(name = "", values = c("Model" = "solid")) +
                         labs(title = "Risk of chronic carriage by age at infection",
                              y = "Risk (proportion)", x = "Age (years)",
                              colour = "Data",
                              caption = "Gambian studies - Keneba Manduar cohort: Whittle 1990* | GHIS: Wild 1993, Fortuin 1993\nWest African studies - Senegal: Barin 1981, Marinier 1985*, Coursaget 1987* |  Liberia: Prince 1981\n* In Edmunds 1993 review") +
                         theme_bw() +
                         theme(plot.title = element_text(hjust = 0.5),
                               plot.caption = element_text(hjust = 0, size = 6),
                               legend.title = element_text(size = 9)) +
                         ylim(0,1) +
                         xlim(0,30))

  ## Mortality curves
  # Add articifial zeros at first timestep to allow plotting of step curves
  mortality_curves_zeros <- out_mat[[i]]$mapped_output$mortality_curves
  mortality_curves_zeros$time_interval_years <- 0
  mortality_curves_zeros$data_value <- 0
  mortality_curves_zeros$model_value <- 0
  mortality_curves_zeros$number_at_risk <- mortality_curves_zeros$sample_size
  mortality_curves_zeros <- unique(mortality_curves_zeros)

  # Add labels for panels with reference and study population
  mort_curves_labels <- c("Shimakawa 2016:\ncumulative mortality in\ncomp. cirrhosis patients\n(Gambia)",
                          "Yang 2017:\ncumulative mortality in\n HCC patients\n(sS Africa, not Gambia)",
                          "Diarra 2010:\ncumulative HCC incid. in\n mixed cirrhosis patients\n(Mali)",
                          "Diarra 2010:\ncumulative mortality in\nmixed cirrhosis patients\n(Mali)",
                          "Bah 2011 (IARC):\ncumulative mortality in\nHCC patients\n(Gambia)")
  names(mort_curves_labels) <- c("shadow4_cum_mortality", "shadow5_cum_mortality",
                                 "shadow6_cum_hcc_incidence", "shadow6_cum_mortality",
                                 "shadow7_cum_mortality")

  p_mort_curves <- print(ggplot(data = rbind(out_mat[[i]]$mapped_output$mortality_curves,
                                             mortality_curves_zeros)) +
                           geom_step(aes(x = time_interval_years, y = model_value, linetype = "Model")) +
                           geom_point(aes(x = time_interval_years, y = data_value, colour = "Data"),
                                      shape = 4, stroke = 1.5) +
                           facet_grid(~ outcome, scales = "free",
                                      labeller = labeller(outcome = mort_curves_labels)) +
                           scale_colour_manual(name = "", values = c("Data" = "red")) +
                           scale_linetype_manual(name = "", values = c("Model" = "solid")) +
                           labs(title = "Cumulative probability of death/HCC over time",
                                y = "Cumulative probability", x = "Follow-up time (years)") +
                           theme_bw() +
                           theme(plot.title = element_text(hjust = 0.5),
                                 legend.position = "bottom",
                                 strip.text.x = element_text(size = 8)))

  ## ORs
  p_or <- print(ggplot(data = out_mat[[i]]$mapped_output$odds_ratios) +
                  geom_col(aes(x = gsub("odds_ratio_", "", outcome), y = log(model_value),
                               fill = "Model")) +
                  geom_point(aes(x = gsub("odds_ratio_", "", outcome), y = log(data_value),
                                 colour = "Data"),
                             shape = 4, size = 3, stroke = 2) +
                  geom_errorbar(aes(x = gsub("odds_ratio_", "", outcome),
                                    ymax = log(ci_upper), ymin = log(ci_lower)), col= "red", width = 0.2) +
                  geom_hline(aes(yintercept=0), colour = "blue") +
                  geom_text(aes(3.5, 0.25, label = ">0 positive association", vjust = 0.25, hjust = 0),
                            size = 3, colour = "blue") +
                  geom_text(aes(3.5, -0.25, label = "<0 negative association", vjust = 0.25, hjust = 0),
                            size = 3, colour = "blue") +
                  coord_cartesian(xlim = c(0.75,3.25), clip = "off") +
                  labs(title = "Log odds ratios for liver disease outcomes",
                       subtitle = "Gambia Liver Cancer Study (Mendy 2010)\nKeneba Manduar chronic carrier cohort (Shimakawa 2016)",
                       y = "ln(OR)", x = "") +
                  scale_x_discrete(breaks=c("current_hbeag_positivity_and_cirrhosis",
                                            "current_hbeag_positivity_and_hcc",
                                            "male_sex_and_significant_liver_fibrosis_or_cirrhosis"),
                                   labels=c("Odds of cirrhosis in\ncurrent HBeAg-positives\nvs.\ncurrent HBeAg-negatives",
                                            "Odds of HCC in\ncurrent HBeAg-positives vs.\ncurrent HBeAg-negatives",
                                            "Odds of\nsignificant liver fibrosis\nor cirrhosis\nin males vs. females")) +
                  scale_fill_manual("", values = c("Model" = "gray35")) +
                  scale_colour_manual("", values = c("Data" = "red")) +
                  theme_bw() +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5, size = 8),
                        axis.text.x = element_text(size = 8),
                        plot.margin = unit(c(1,3,1,1), "lines")))

  # Natural history prevalence plots
  out_mat[[i]]$mapped_output$nat_hist_prevalence$model_num <-
    gsub(".*[[:digit:]]{4}_", "",out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique)

  # GMB1 PROLIFICA plots: infection phase in chronic carriers
  gmb1_facet_labels <- c("Male blood donors", "Community screening pop.")
  names(gmb1_facet_labels) <- c("Male", "Mixed")

  plot1_gmb1 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                     id_paper == "GMB1" &
                                       model_num != "cc_dcc" & model_num != "hcc"),
                       aes(x = model_num)) +
    geom_col(aes(y = model_value))+
    geom_point(aes(y = data_value), shape = 4, size = 1.5, stroke = 2, col = "red") +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    facet_grid(~sex, scales = "free", labeller = labeller(sex = gmb1_facet_labels)) +
    scale_x_discrete(breaks=c("ic", "ir_enchb", "it_ic"),
                     labels=c("IC", "IR or\nENCHB", "IT or IC")) +
    labs(title = "Infection phase in chronic carriers",
         subtitle = "PROLIFICA (Lemoine 2016)",
         y = "Prevalence (proportion)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    ylim(0,1)

  # GMB1 plots: liver disease in chronic carriers
  plot2_gmb1 <-ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "GMB1" &
                                      (model_num == "cc_dcc" | model_num == "hcc")),
                      aes(x = gsub("_", " or ", toupper(model_num)))) +
    geom_col(aes(y = model_value))+
    geom_point(aes(y = data_value), shape = 4, size = 1.5, stroke = 2, col = "red") +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    facet_grid(~sex, scales = "free", labeller = labeller(sex = gmb1_facet_labels)) +
    labs(title = "Liver disease in chronic carriers",
         subtitle = "PROLIFICA (Lemoine 2016)",
         y = "Prevalence (proportion)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    ylim(0,0.1)

  # 1-1 plots: infection phase in chronic carriers without liver disease
  study_1_facet_labels <- c("1986: median age 11 years", "2013: median age 38 years")
  names(study_1_facet_labels) <- c(1986, 2013)

  plot1_1 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                  grepl(".*it,_ir,_ic_and_enchb$", out_mat[[i]]$mapped_output$nat_hist_prevalence$outcome)),
                    aes(x = toupper(model_num))) +
    geom_col(aes(y = model_value))+
    geom_point(aes(y = data_value), shape = 4, size = 1.5, stroke = 2, col = "red") +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    facet_grid(~time, scales = "free", labeller = labeller(time = study_1_facet_labels)) +
    labs(title = "Infection phase in chronic carriers\nwithout liver disease",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Prevalence (proportion)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    ylim(0,1)

  # 1-1 plots: liver disease in chronic carriers
  plot2_1 <-ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                 id_paper == "1" &
                                   (outcome == "hcc_prevalence_in_chronic_carriers" |
                                      outcome == "cc_and_dcc_prevalence_in_chronic_carriers")),
                   aes(x = gsub("_", " or ", toupper(model_num)))) +
    geom_col(aes(y = model_value))+
    geom_point(aes(y = data_value), shape = 4, size = 1.5, stroke = 2, col = "red") +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    facet_grid(~time, scales = "free_x", labeller = labeller(time = study_1_facet_labels)) +
    labs(title = "Liver disease in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Prevalence (proportion)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    ylim(0,0.1)

  # 1-1 plots: chronic carriers by age
  plot3_1 <- ggplot(data = out_mat[[i]]$mapped_output$nat_hist_prevalence[
    out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_2013_ir_enchb_cc_dcc",]) +
    geom_col(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = model_value, fill = "Model"))+
    geom_point(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = data_value, colour = "Data"),
               shape = 4, size = 1.5, stroke = 2) +
    geom_errorbar(aes(x = reorder(paste(age_min,"-",age_max), age_min),
                      ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.4) +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Significant liver fibrosis or cirrhosis in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Prevalence (proportion)", x = "Age (years)",
         caption = "\nIT = Immune tolerant, IR = Immune reactive, IC = Inactive carrier, ENCHB = HBeAg-negative chronic hepatitis B,\nCC = Compensated cirrhosis, DCC = Decompensated cirrhosis, HCC = Hepatocellular carcinoma") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8),
          legend.margin = margin(t = 0, unit="cm"),
          legend.position = "left",
          plot.caption = element_text(hjust = 0, size = 8)) +
    ylim(0,0.75)

  ## Natural history prevalence PLOT 1
  p_nat_hist_prev1 <- grid.arrange(plot1_gmb1, plot1_1, plot2_gmb1, plot2_1,
                                   plot3_1, nrow = 3,
                                   layout_matrix = rbind(c(1,2), c(3,4), c(5,5)),
                                   top = "Prevalence measures in chronic carriers")

  # A4 Proportion of deaths from DCC and HCC in comp. cirrhosis cohort
  plot_nat_hist_a4 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                           id_paper == "A4")) +
    geom_col(aes(x = model_num, y = model_value, fill = "Model"))+
    geom_point(aes(x = model_num, y = data_value, colour = "Data"),
               shape = 4, size = 3, stroke = 2) +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Proportion of deaths from\nDCC and HCC in\ncompensated cirrhosis cohort",
         subtitle = "(Shimakawa 2016)",
         y = "Proportion", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 9),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank(),
          legend.position = "left") +
    ylim(0,1)

  # GMB2 Cirrhosis prevalence in HBsAg-positive HCC patients (GLCS)
  plot_nat_hist_gmb2 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                             id_paper == "GMB2")) +
    geom_col(aes(x = model_num, y = model_value))+
    geom_point(aes(x = model_num, y = data_value),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = model_num,
                      ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.2) +
    scale_x_discrete(breaks=c("incident_hcc_cases_from_cc",
                              "incident_hcc_cases_from_dcc"),
                     labels=c("Compensated\ncirrhosis",
                              "Decompensated\ncirrhosis")) +
    labs(title = "Cirrhosis prevalence in\nHBsAg-positive HCC patients",
         subtitle = "Gambia Liver Cancer Study (Umoh 2011)",
         y = "Proportion", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(size = 7),
          plot.subtitle = element_text(hjust = 0.5, size = 8)) +
    ylim(0,1)

  # HBeAg prevalence in liver disease patients: GMB12 and GMB15
  nat_hist_glcs_facet_labels <- c("Ryder 1992: HCC patients",
                                  "Mendy 2010 (GLCS): HCC patients",
                                  "Mendy 2010 (GLCS): cirrhosis patients")
  names(nat_hist_glcs_facet_labels) <- c("id_gmb12_1_1982_hbeag_hcc",
                                         "id_gmb15_1_1999_hbeag_hcc",
                                         "id_gmb15_2_1999_hbeag_cirrhosis")

  plot_nat_hist_glcs <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                             outcome == "hbeag_prevalence_in_hcc" |
                                               outcome == "hbeag_prevalence_in_cirrhosis")) +
    geom_col(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = model_value))+
    geom_point(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = data_value),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = reorder(paste(age_min,"-",age_max),age_min),
                      ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.2) +
    facet_grid(~id_unique, scales = "free_x", labeller = labeller(id_unique = nat_hist_glcs_facet_labels)) +
    labs(title = "HBeAg prevalence in HBsAg-positive HCC/cirrhosis patients",
         y = "Proportion", x = "Age group (years)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(size = 8),
          strip.text.x = element_text(size = 8)) +
    ylim(0,1)

  ## Natural history prevalence PLOT 2
  p_nat_hist_prev2 <- grid.arrange(plot_nat_hist_a4,  plot_nat_hist_gmb2, plot_nat_hist_glcs, nrow = 2,
                                   layout_matrix = rbind(c(1,2),
                                                         c(3,3)),
                                   top = "Prevalence measures in liver disease patients")

  # Vertical transmission plots
  # 1-1 plots: chronic infections due to vertical transmission
  plot_nat_hist_1 <- ggplot(data = out_mat[[i]]$mapped_output$nat_hist_prevalence[
    out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_1986_incident_chronic_births",]) +
    geom_col(aes(x = model_num, y = model_value))+
    geom_point(aes(x = model_num, y = data_value),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = model_num, ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.1) +
    labs(title = "Proportion of chronic infection\ndue to vertical transmission\nin unvaccinated pop.",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Proportion", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank()) +
    ylim(0,1)

  # MTCT risk
  plot_mtct <- ggplot(data = out_mat[[i]]$mapped_output$mtct_risk,
                      aes(x = paste(paper_first_author, paper_year))) +
    geom_col(aes(y = model_value, fill = "Model"))+
    geom_point(aes(y = data_value, colour = "Data"),
               shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.1) +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Mother-to-child transmission risk",
         y = "Proportion", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,1)

  ## Progression rates
  out_mat[[i]]$mapped_output$progression_rates$type <-
    gsub("shadow[[:digit:]].{0,1}_", "", out_mat[[i]]$mapped_output$progression_rates$outcome)
  out_mat[[i]]$mapped_output$progression_rates$type <- gsub("_.$", "", out_mat[[i]]$mapped_output$progression_rates$type)

  # Study 1: HCC incidence
  plot_1_hcc_incidence <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                               id_paper == "1" & type == "hcc_incidence")) +
    geom_col(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = model_value*100000))+
    geom_point(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = data_value*100000),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), ymax = ci_upper*100000, ymin = ci_lower*100000),
                  col = "red", width = 0.1)  +
    facet_grid(~sex, scales = "free") +
    labs(title = "HCC incidence rate in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Cases per 100000 PY", x = "Baseline age group (years)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8))

  # Study 1: DCC incidence
  plot_1_dcc_incidence <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                               id_paper == "1" & type == "dcc_incidence")) +
    geom_col(aes(x = outcome, y = model_value*100000))+
    geom_point(aes(x = outcome, y = data_value*100000),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = outcome, ymax = ci_upper*100000, ymin = ci_lower*100000),
                  col = "red", width = 0.1)  +
    labs(title = "DCC incidence rate in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort\n(Shimakawa 2016)",
         y = "Cases per 100000 PY", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank()) +
    ylim(0,100)


  # Study 1: Mortality
  plot_1_mortality <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                           id_paper == "1" & type == "mortality")) +
    geom_col(aes(x = sex, y = model_value*100000))+
    geom_point(aes(x = sex, y = data_value*100000),
               col = "red", shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(x = sex, ymax = ci_upper*100000, ymin = ci_lower*100000),
                  col = "red", width = 0.1)  +
    labs(title = "All-cause mortality rate in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Deaths per 100000 PY", x = "Sex") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8)) +
    ylim(0,1000)

  # Study A6: Mortality in cirrhosis cohort
  plot_a6_mortality <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                            id_paper == "A6")) +
    geom_col(aes(x = outcome, y = model_value*100, fill = "Model"))+
    geom_point(aes(x = outcome, y = data_value*100, colour = "Data"),
               shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(x = outcome, ymax = ci_upper*100, ymin = ci_lower*100),
                  col = "red", width = 0.1)  +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Mortality rate in liver disease patients",
         subtitle = "Mixed cohort of Nigerian CC, DCC and HCC patients\n(Olubuyide 1996)",
         y = "Deaths per 100 PY", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank(),
          legend.position = "bottom") +
    ylim(0,60)

  ## Combined liver disease incidence and mortality rates
  p_prog_rates1 <- grid.arrange(plot_1_hcc_incidence, plot_1_dcc_incidence,
                                plot_1_mortality, plot_a6_mortality, nrow = 2, widths = 4:3,
                                top = "Liver disease-related rates")

  # Study 1: HBeAg loss
  plot_1_eag_loss <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                          id_paper == "1" & type == "eag_loss")) +
    geom_col(aes(x = sex, y = model_value*100))+
    geom_point(aes(x = sex, y = data_value*100),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = sex, ymax = ci_upper*100, ymin = ci_lower*100),
                  col = "red", width = 0.1)  +
    labs(title = "Rate of HBeAg loss in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Cases per 100 PY", x = "Sex") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8)) +
    ylim(0,12)

  # Study 6: HBsAg loss
  plot_6_sag_loss <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                          id_paper == "6")) +
    geom_col(aes(x = outcome, y = model_value*100, fill = "Model"))+
    geom_point(aes(x = outcome, y = data_value*100, colour = "Data"),
               shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(x = outcome, ymax = ci_upper*100, ymin = ci_lower*100),
                  col = "red", width = 0.1)  +
    scale_fill_manual("", values = c("Model" = "gray35")) +
    scale_colour_manual("", values = c("Data" = "red")) +
    labs(title = "Rate of HBsAg loss in\nchronic carrier children",
         subtitle = "Coursaget 1987 (Senegal)",
         y = "Cases per 100 PY", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_blank()) +
    ylim(0,5)

  ## Combined seromarker loss rates
  p_prog_rates2 <- grid.arrange(plot_1_eag_loss, plot_6_sag_loss, nrow = 1, widths = 2:1,
                                top = "Seromarker clearance rates")

  # Transmission-related data from GMB6 and GMB7
  plot_horizontal_transmission <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                                       id_paper == "GMB6" | id_paper == "GMB7")) +
    geom_col(aes(x = gsub(".*_","",outcome), y = model_value))+
    geom_point(aes(x = gsub(".*_","",outcome), y = data_value),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = gsub(".*_","",outcome), ymax = ci_upper, ymin = ci_lower),
                  col = "red", width = 0.1)  +
    scale_x_discrete(breaks=c("foi",
                              "incidence"),
                     labels=c("Force of infection\nin children in\nKeneba and Manduar\n(Whittle 1990)",
                              "Chronic infection incidence\nin susceptible children\n(Ryder 1984)")) +
    labs(title = "Horizontal transmission-related rates",
         y = "Rate (per PY)", x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,2)


  ## Combined transmission-related plot
  p_transmission_rates <- grid.arrange(plot_mtct, plot_nat_hist_1, plot_horizontal_transmission,
                                       layout_matrix = rbind(c(1,1),
                                                             c(2,3)),
                                       top = "Transmission-related measures")

  # List of all plots
  plot_list[[i]] <- list(p_parms, p_hbsag1, p_hbsag2, p_antihbc, p_hbeag,
                         p_globocan1, p_globocan2, p_hcc_pattern1, p_hcc_pattern2,
                         p_gbd, p_ld_demog,
                         p_p_chronic, p_mort_curves, p_or,
                         p_nat_hist_prev1, p_nat_hist_prev2,
                         p_prog_rates1, p_prog_rates2, p_transmission_rates)

}
dev.off()

prev <- 0
for (i in 1:length(out_mat)) {
  prev[i] <- sum(out_mat[[i]]$full_output$carriers[201,])/sum(out_mat[[i]]$full_output$pop[201,])
}
prev[best_fit_ids]

median(res_mat[which(prev > 0.1 & prev < 0.25),"error_term"])

plot(out_mat[[1]]$full_output$time,
     (apply(out_mat[[90]]$full_output$carriers,1,sum)/apply(out_mat[[90]]$full_output$pop,1,sum)),
     ylim = c(0,0.5))
abline(h = 0.11)
