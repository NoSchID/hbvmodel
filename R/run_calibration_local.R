#################################
### Run calibration on laptop ###
#################################

### Load packages and source file ----
require(here)
require(truncnorm)
library(grid)
library(ggplot2)
library(gridExtra)

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
                   fit_model_full_output(default_parameter_list = parameter_list,
                             data_to_fit = calibration_datasets_list,
                             parms_to_change =
                               list(b1 = as.list(x)$b1,
                                    b2 = as.list(x)$b2,
                                    mtct_prob_s = as.list(x)$mtct_prob_s,
                                    mtct_prob_e = 0.6,
                                    alpha = 7.07812500,
                                    b3 = 0.01976563,
                                    p_chronic_function_r = 0.65595093,
                                    p_chronic_function_s = 0.46,
                                    eag_prog_function_rate = 0,
                                    pr_it_ir = 0.1,
                                    pr_ir_ic = 0.83906250,
                                    pr_ir_cc_female = 0.01,
                                    pr_ir_cc_age_threshold =0,
                                    pr_ir_enchb = 0.005,
                                    pr_ic_enchb = 0.01,
                                    sag_loss_slope = 0.000451,
                                    pr_enchb_cc_female =  0.01181470,
                                    pr_cc_dcc = 0.04,
                                    hccr_dcc = 0.08953125,
                                    hccr_it = 6.25,
                                    hccr_ir = 15,
                                    hccr_enchb = 11.25000000,
                                    hccr_cc = 25,
                                    cirrhosis_male_cofactor = 5.15625000,
                                    cancer_prog_coefficient_female = 0.00022,
                                    cancer_age_threshold = 0,
                                    cancer_male_cofactor = 3.625,
                                    mu_cc = 0.005,
                                    mu_hcc = 1.81250000,
                                    mu_dcc = 0.8)))

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
                         pr_ir_enchb = rgamma(n_sims, 1.22, 44.20),
                         pr_ir_cc_female = runif(n_sims, 0.005, 0.05),
                         pr_ir_cc_age_threshold = sample(0:15,n_sims,replace=TRUE),
                         pr_ic_enchb = rgamma(n_sims, 2.18, 118.16),
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

### Output calibration plots code (PDF) ----
library(grid)
library(ggplot2)
library(gridExtra)
# Loop to create plot set for every parameter combination
pdf(file = here("output/fits/Frequentist fit", "sse_scale_max_start_at_manual_fit_test_horizontal.pdf"), paper="a4r")
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


### Create 1,000,000 LHS samples ----
library(lhs)
library(truncnorm)

lhs_mat <- randomLHS(1000000, 32)

b1 <- qunif(lhs_mat[,1], 0.03, 0.7)
b2 <- qunif(lhs_mat[,2], 0, b1)
b3 <- qunif(lhs_mat[,3], 0, b1)
mtct_prob_s <- qbeta(lhs_mat[,4], 1.5,13.5)
mtct_prob_e <- qunif(lhs_mat[,5], mtct_prob_s, 0.9)
hccr_it <- qtruncnorm(lhs_mat[,6], a=1, b=Inf, mean=6, sd=3) # normal truncated at 1
hccr_cc <- qunif(lhs_mat[,7], hccr_it, 100)
hccr_enchb <- qunif(lhs_mat[,8], hccr_it, hccr_cc)
hccr_ir <- runif(lhs_mat[,9], hccr_enchb, hccr_cc)

params_mat <- data.frame(b1 = b1,
                         b2 = b2,
                         b3 = b3,
                         mtct_prob_s = mtct_prob_s,
                         mtct_prob_e = mtct_prob_e,
                         alpha = qunif(lhs_mat[,10], 1.5,10),
                         p_chronic_in_mtct = qbeta(lhs_mat[,11], 10.49,1.3),
                         p_chronic_function_r = qnorm(lhs_mat[,12],0.65,0.1),
                         p_chronic_function_s = qnorm(lhs_mat[,13],0.46,0.1),
                         pr_it_ir = qgamma(lhs_mat[,14],3.63,26.27),
                         pr_ir_ic = qunif(lhs_mat[,15], 0,1),
                         eag_prog_function_rate = qunif(lhs_mat[,16],0,0.01),
                         pr_ir_enchb = qgamma(lhs_mat[,17], 1.22, 44.20),
                         pr_ir_cc_female = qunif(lhs_mat[,18], 0.005, 0.05),
                         pr_ir_cc_age_threshold = floor(qunif(lhs_mat[,19],0,15)),
                         pr_ic_enchb = qgamma(lhs_mat[,20], 2.18, 118.16),
                         sag_loss_slope = qnorm(lhs_mat[,21], 0.0004106, 0.00005),
                         pr_enchb_cc_female = qgamma(lhs_mat[,22], 1.23, 22.33),
                         cirrhosis_male_cofactor = qtruncnorm(lhs_mat[,23], a = 1, mean = 3.5, sd = 4),
                         pr_cc_dcc = qgamma(lhs_mat[,24],17.94,423.61),
                         cancer_prog_coefficient_female = qunif(lhs_mat[,25], 0.0001, 0.0003),
                         cancer_age_threshold = floor(qunif(lhs_mat[,26],0,15)),
                         cancer_male_cofactor = qtruncnorm(lhs_mat[,27], a = 1, mean = 3.5, sd = 4),
                         hccr_it = hccr_it,
                         hccr_ir = hccr_ir,
                         hccr_enchb = hccr_enchb,
                         hccr_cc = hccr_cc,
                         hccr_dcc = qgamma(lhs_mat[,28], 3.08, 29.76),
                         mu_cc = qgamma(lhs_mat[,29], 4.25, 124.91),
                         mu_dcc = qgamma(lhs_mat[,30], 2.18, 1.18),
                         mu_hcc = qgamma(lhs_mat[,31], 2.18, 1.18),
                         vacc_eff = qbeta(lhs_mat[,32], 7.07, 0.37))

library(here)
#save(params_mat, file = here("calibration_input", "lhs_samples_1000000.Rdata"))

### Run with LHS calibration best parameters----
# Run simulations with the best-fit parameter sets

# From first calibration:
#load(here("calibration", "input", "accepted_parmsets_119_060120.Rdata")) # params_mat_targets5
#load(here("calibration", "input", "target_threshold_parms_328_060120.Rdata")) # params_mat_targets5_2
# params_mat_targets5: chosen so that 100% of accepted sets fall within targets
# params_mat_targets5_2: chosen so that 99% of accepted sets fall within targets

# From recalibration
load(here("calibration", "input", "accepted_parmsets_123_180520.Rdata")) # params_mat_accepted
params_mat_targets5 <- params_mat_accepted  # rename so I don't have to change code

### Visualise posteriors for different cutoffs
load(here("calibration", "input", "lhs_samples_1000000.Rdata"))

quantile(params_mat_targets5$pr_ic_enchb, prob = c(0.025,0.5,0.975))
quantile(params_mat$pr_ic_enchb, prob = c(0.025,0.5,0.975))

# Progression to treatment eligibility test
range((exp(params_mat_targets5$eag_prog_function_rate*ages[which(ages==33)])*params_mat_targets5$pr_it_ir)+
  params_mat_targets5$pr_ic_enchb)
mean((exp(params_mat_targets5$eag_prog_function_rate*ages[which(ages==33)])*params_mat_targets5$pr_it_ir)+
        params_mat_targets5$pr_ic_enchb)

# Generate table of all prior and posterior summaries
posterior_summary <- gather(params_mat_accepted, key = "parameter", value = "sim") %>%
  group_by(parameter) %>%
  summarise(post_mean = round(mean(sim),4),
            post_median = round(median(sim),4),
            post_cri_lower = round(quantile(sim, prob = 0.025),4),
            post_cri_upper = round(quantile(sim, prob = 0.975),4))

prior_summary <- gather(params_mat, key = "parameter", value = "sim") %>%
  group_by(parameter) %>%
  summarise(prior_mean = round(mean(sim),4),
            prior_median = round(median(sim),4),
            prior_cri_lower = round(quantile(sim, prob = 0.025),4),
            prior_cri_upper = round(quantile(sim, prob = 0.975),4))

prior_posterior_summary <- left_join(prior_summary, posterior_summary, by ="parameter")
#write.csv(prior_posterior_summary, file=here("calibration", "output", "prior_posterior_summary_table_150720.csv"), row.names = FALSE)

plot_prior_posterior <- function(parm) {
  plot(density(posterior[,parm]), xlim = c(min(min(prior[,parm]),min((posterior[,parm]))), max(max(prior[,parm]),max((posterior[,parm])))),
       ylim = c(0, max(max(density(prior[,parm])$y),max((density(posterior[,parm])$y)))), main= parm,
       lwd=3, col="red")
  lines(density(prior[,parm]), lwd=3, lty=2, col="blue")
  legend("bottomleft", legend=c("prior density","posterior density"),
         col=c("blue","red"), lty=c(3,1), lwd=c(3,3), cex = 1)
}

prior <- as.data.frame(params_mat)
posterior <- params_mat_targets5

par(mfrow=c(1,3))
plot_prior_posterior("b1")
plot_prior_posterior("b2")
plot_prior_posterior("b3")
plot_prior_posterior("alpha")
plot_prior_posterior("mtct_prob_s")
plot_prior_posterior("mtct_prob_e")
plot_prior_posterior("p_chronic_in_mtct")
plot_prior_posterior("p_chronic_function_r")
plot_prior_posterior("p_chronic_function_s")
plot_prior_posterior("pr_it_ir")
plot_prior_posterior("pr_ir_ic")
plot_prior_posterior("eag_prog_function_rate")
plot_prior_posterior("pr_ir_enchb")
plot_prior_posterior("pr_ir_cc_female")
plot_prior_posterior("pr_ir_cc_age_threshold")
plot_prior_posterior("pr_ic_enchb")
plot_prior_posterior("sag_loss_slope")
plot_prior_posterior("pr_enchb_cc_female")
plot_prior_posterior("cirrhosis_male_cofactor")
plot_prior_posterior("pr_cc_dcc")
plot_prior_posterior("cancer_prog_coefficient_female")
plot_prior_posterior("cancer_age_threshold")
plot_prior_posterior("cancer_male_cofactor")
plot_prior_posterior("hccr_it")
plot_prior_posterior("hccr_ir")
plot_prior_posterior("hccr_enchb")
plot_prior_posterior("hccr_cc")
plot_prior_posterior("hccr_dcc")
plot_prior_posterior("mu_cc")
plot_prior_posterior("mu_dcc")
plot_prior_posterior("mu_hcc")
plot_prior_posterior("vacc_eff")

library(HDInterval)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

comp_prior_post2 <- rbind(cbind(apply(prior, 2, Mode),
                                t(apply(prior, 2, hdi))),
                          cbind(apply(posterior, 2, Mode),
                                t(apply(posterior, 2, hdi))))

colnames(comp_prior_post2) <- c("mode", "lower", "upper")
comp_prior_post2 <- as.data.frame(comp_prior_post2)
comp_prior_post2$type <- c(rep("prior", 32), rep("post", 32))
comp_prior_post2$parm <- c(colnames(prior), colnames(posterior))

comp_prior_post <- rbind(t(apply(prior, 2, quantile, prob = c(0.025, 0.5, 0.975))),
                         t(apply(posterior, 2, quantile, prob = c(0.025, 0.5, 0.975))))
colnames(comp_prior_post) <- c("lower", "median", "upper")
comp_prior_post <- as.data.frame(comp_prior_post)
comp_prior_post$type <- c(rep("prior", 32), rep("post", 32))
comp_prior_post$parm <- c(colnames(prior), colnames(posterior))

ggplot(comp_prior_post) +
  geom_point(aes(x=type, y = median)) +
  geom_errorbar(aes(x=type, ymin = lower, ymax = upper)) +
  facet_wrap(~parm, scales = "free")

ggplot(comp_prior_post2) +
  geom_point(aes(x=type, y = mode)) +
  geom_errorbar(aes(x=type, ymin = lower, ymax = upper)) +
  facet_wrap(~parm, scales = "free")

ggplot(comp_prior_post[comp_prior_post$parm == "alpha",]) +
  geom_point(aes(x=type, y = median)) +
  geom_errorbar(aes(x=type, ymin = lower, ymax = upper))

ggplot(comp_prior_post2[comp_prior_post2$parm == "alpha",]) +
  geom_point(aes(x=type, y = mode)) +
  geom_errorbar(aes(x=type, ymin = lower, ymax = upper))


(((0.000238 * (ages - 9))^2)  * c(rep(0, times = 9/da),rep(1, times = n_agecat - 9/da)))[41:121]
# 7.341062e-05 female 45 year old IC, 0.003517926 IR, 0.002378231 ENCHB
# 0.0013 male 45 year old IC, 0.06229759 IR, 0.04211516 ENCHB
# In mixed female cohort assuming 90% IC, 5% IR and 5% ENCHB: 0.00036
# In mixed male cohort: 0.0064

# Not in parallel
out_mat <- apply(params_mat_targets5, 1,
                    function(x)
                      fit_model_full_output(default_parameter_list = parameter_list,
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
                                                   vacc_eff = as.list(x)$vacc_eff)))
save(out_mat, file = here("calibration", "output", "model_fit_output_123_180520.Rdata"))
out_mat_copy <- out_mat

#new_vacc_meth <- out_mat
#output1 <- new_vacc_meth2
#output2 <- new_vacc_meth
#new_vacc_meth <- output2

# CHECK VACCINE COVERAGE IN MODEL
# For vaccine coverage, have to substract earlier recovery % (ca 4 %)
vacc_cov_2020 <- 119L
vacc_cov_1997 <- 119L
vacc_cov_2010 <- 119L
rec_1990 <- 119L

for (i in 1:119) {
rec_1990[i] <- (new_vacc_meth[[i]]$full_output$ever_infected[221,which(ages==1)]-
                  new_vacc_meth[[i]]$full_output$carriers[221,which(ages==1)])/
  new_vacc_meth[[i]]$full_output$pop[221,which(ages==1)]

vacc_cov_1997[i] <- ((new_vacc_meth[[i]]$full_output$ever_infected[235,which(ages==1)]-
                        new_vacc_meth[[i]]$full_output$carriers[235,which(ages==1)])/
                       new_vacc_meth[[i]]$full_output$pop[235,which(ages==1)])-rec_1990[i]

vacc_cov_2010[i] <- ((new_vacc_meth[[i]]$full_output$ever_infected[261,which(ages==1)]-
                  new_vacc_meth[[i]]$full_output$carriers[261,which(ages==1)])/
  new_vacc_meth[[i]]$full_output$pop[261,which(ages==1)])-rec_1990[i]

vacc_cov_2020[i] <- ((new_vacc_meth[[i]]$full_output$ever_infected[280,which(ages==1)]-
                        new_vacc_meth[[i]]$full_output$carriers[280,which(ages==1)])/
                       new_vacc_meth[[i]]$full_output$pop[280,which(ages==1)])-rec_1990[i]
}

test <- as.data.frame(cbind(vacc_cov_2020,unlist(params_mat_targets5$vacc_eff)))
test$high_eff <- 1
test$high_eff[test$V2<0.97] <- 0

plot(x=1:119, y = vacc_cov_1997, ylim = c(0,1))
points(x=1:119, y = vacc_cov_2010, col = "blue")
points(x=1:119, y = vacc_cov_2020, col = "red")
points(x=1:119, y = params_mat_targets5$vacc_eff, col = "green")

boxplot(vacc_cov_2020~high_eff, test, ylim = c(0.7,1))

# In parallel
library(parallel)
cl <- makeCluster(4)
clusterEvalQ(cl, {library(dplyr); library(tidyr); library(deSolve); library(binom)})
clusterExport(cl, ls())
out_mat <- parApply(cl = cl, params_mat_targets5_2, 1,
                    function(x)
                      fit_model_full_output(default_parameter_list = parameter_list,
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
# Important: stop cluster!!
stopCluster(cl)
#save(out_mat, file = here("calibration", "output", "model_fit_output_328_070120.Rdata"))

my_output <- out_mat

out_mat <- output2

#Extract mapped output
best_fits_mapped_output <- lapply(out_mat, "[[", "mapped_output")
# Extract model values from mapped output
best_fits_model_values <- lapply(best_fits_mapped_output, function(x) sapply(x, '[', "model_value"))
# Turn model values into vectors within list
best_fits_model_values2 <- lapply(best_fits_model_values, function(x) unlist(x))
# Turn into matrix with columns = runs, rows = summary stats
best_fits_model_values_mat <- do.call("cbind", best_fits_model_values2)

# Calculate the median and 95% percentile for each summary statistic
median_model_values <- apply(best_fits_model_values_mat,1,median)
ci_lower_model_values <- apply(best_fits_model_values_mat,1,quantile, probs = 0.025)
ci_upper_model_values <- apply(best_fits_model_values_mat,1,quantile, probs = 0.975)


# Need to replace the model values in out_mat_median[[1]]$mapped_output etc with these average values!
#best_fits_500_sse[1,-ncol(best_fits_500_sse)]
out_mat_median <- apply(params_mat[1,],1,
                       function(x)
                       fit_model_full_output(default_parameter_list = parameter_list,
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
out_mat_ci_lower <- out_mat_median
out_mat_ci_upper <- out_mat_median

out_mat_median[[1]]$mapped_output$globocan_hcc_incidence$model_value <-
  median_model_values[grepl("^globocan_hcc_incidence.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$gbd_cirrhosis_mortality$model_value <-
  median_model_values[grepl("^gbd_cirrhosis_mortality.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$risk_of_chronic_carriage$model_value <-
  median_model_values[grepl("^risk_of_chronic_carriage.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$seromarker_prevalence$model_value <-
  median_model_values[grepl("^seromarker_prevalence.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$nat_hist_prevalence$model_value <-
  median_model_values[grepl("^nat_hist_prevalence.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$mtct_risk$model_value <-
  median_model_values[grepl("^mtct_risk.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$progression_rates$model_value <-
  median_model_values[grepl("^progression_rates.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$mortality_curves$model_value <-
  median_model_values[grepl("^mortality_curves.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$odds_ratios$model_value <-
  median_model_values[grepl("^odds_ratios.*",names(median_model_values))]
out_mat_median[[1]]$mapped_output$mapped_liver_disease_demography$model_value <-
  median_model_values[grepl("^mapped_liver_disease_demography.*",names(median_model_values))]
# out_mat_ci_lower
out_mat_ci_lower[[1]]$mapped_output$globocan_hcc_incidence$model_value <-
  ci_lower_model_values[grepl("^globocan_hcc_incidence.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$gbd_cirrhosis_mortality$model_value <-
  ci_lower_model_values[grepl("^gbd_cirrhosis_mortality.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$risk_of_chronic_carriage$model_value <-
  ci_lower_model_values[grepl("^risk_of_chronic_carriage.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$seromarker_prevalence$model_value <-
  ci_lower_model_values[grepl("^seromarker_prevalence.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$nat_hist_prevalence$model_value <-
  ci_lower_model_values[grepl("^nat_hist_prevalence.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$mtct_risk$model_value <-
  ci_lower_model_values[grepl("^mtct_risk.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$progression_rates$model_value <-
  ci_lower_model_values[grepl("^progression_rates.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$mortality_curves$model_value <-
  ci_lower_model_values[grepl("^mortality_curves.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$odds_ratios$model_value <-
  ci_lower_model_values[grepl("^odds_ratios.*",names(ci_lower_model_values))]
out_mat_ci_lower[[1]]$mapped_output$mapped_liver_disease_demography$model_value <-
  ci_lower_model_values[grepl("^mapped_liver_disease_demography.*",names(ci_lower_model_values))]
# out_mat_ci_upper
out_mat_ci_upper[[1]]$mapped_output$globocan_hcc_incidence$model_value <-
  ci_upper_model_values[grepl("^globocan_hcc_incidence.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$gbd_cirrhosis_mortality$model_value <-
  ci_upper_model_values[grepl("^gbd_cirrhosis_mortality.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$risk_of_chronic_carriage$model_value <-
  ci_upper_model_values[grepl("^risk_of_chronic_carriage.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$seromarker_prevalence$model_value <-
  ci_upper_model_values[grepl("^seromarker_prevalence.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$nat_hist_prevalence$model_value <-
  ci_upper_model_values[grepl("^nat_hist_prevalence.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$mtct_risk$model_value <-
  ci_upper_model_values[grepl("^mtct_risk.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$progression_rates$model_value <-
  ci_upper_model_values[grepl("^progression_rates.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$mortality_curves$model_value <-
  ci_upper_model_values[grepl("^mortality_curves.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$odds_ratios$model_value <-
  ci_upper_model_values[grepl("^odds_ratios.*",names(ci_upper_model_values))]
out_mat_ci_upper[[1]]$mapped_output$mapped_liver_disease_demography$model_value <-
  ci_upper_model_values[grepl("^mapped_liver_disease_demography.*",names(ci_upper_model_values))]


# Plot the median and 95% CI for each summary statistic
out_mat <- out_mat_median

library(ggplot2)
library(gridExtra)
pdf(file = here("calibration", "output", "median_95ci_fit_123_180520.pdf"), paper="a4r")
plot_list = list()
for (i in 1:length(out_mat)) {

  # OUTPUTS

  # Seromarker prevalence ----

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
                      geom_line(data =  subset(out_mat_ci_lower[[i]]$mapped_output$seromarker_prevalence,
                                               outcome == "HBsAg_prevalence"),
                                aes(x = age, y = model_value,
                                    colour = sex), linetype = "dashed") +
                      geom_line(data =  subset(out_mat_ci_upper[[i]]$mapped_output$seromarker_prevalence,
                                               outcome == "HBsAg_prevalence"),
                                aes(x = age, y = model_value,
                                    colour = sex), linetype = "dashed") +
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
      geom_line(data =  subset(out_mat_ci_lower[[i]]$mapped_output$seromarker_prevalence,
                               outcome == "Anti_HBc_prevalence"),
                aes(x = age, y = model_value,
                    colour = sex), linetype = "dashed") +
      geom_line(data =  subset(out_mat_ci_upper[[i]]$mapped_output$seromarker_prevalence,
                               outcome == "Anti_HBc_prevalence"),
                aes(x = age, y = model_value,
                    colour = sex), linetype = "dashed") +
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
      geom_line(data =  subset(out_mat_ci_lower[[i]]$mapped_output$seromarker_prevalence,
                               outcome == "HBeAg_prevalence"),
                aes(x = age, y = model_value,
                    colour = sex), linetype = "dashed") +
      geom_line(data =  subset(out_mat_ci_upper[[i]]$mapped_output$seromarker_prevalence,
                               outcome == "HBeAg_prevalence"),
                aes(x = age, y = model_value,
                    colour = sex), linetype = "dashed") +
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

  # Liver disease rates ----

  ## GLOBOCAN PAF-adjusted cancer incidence and mortality in 2018
  globocan_outcome_facet_labels <- c("HCC case incidence", "HCC mortality")
  names(globocan_outcome_facet_labels) <- c("hcc_incidence", "hcc_mortality")

  p_globocan1 <- print(ggplot(data = subset(out_mat[[i]]$mapped_output$globocan_hcc_incidence,
                                            time == 2018)) +
                         geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                         geom_errorbar(aes(x = paste(age_min,"-",age_max),
                                           ymin = subset(out_mat_ci_lower[[i]]$mapped_output$globocan_hcc_incidence,
                                                         time == 2018)$model_value*100000,
                                           ymax = subset(out_mat_ci_upper[[i]]$mapped_output$globocan_hcc_incidence,
                                                         time == 2018)$model_value*100000),
                                       col = "black", width = 0.2) +
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
                               legend.margin=margin(t = 0, unit="cm")))

  ## GLOBOCAN PAF-adjusted cancer incidence in 1988 and 1998
  p_globocan2 <- print(ggplot(data = subset(out_mat[[i]]$mapped_output$globocan_hcc_incidence,
                                            time != 2018)) +
                         geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                         geom_errorbar(aes(x = paste(age_min,"-",age_max),
                                           ymin = subset(out_mat_ci_lower[[1]]$mapped_output$globocan_hcc_incidence,
                                                         time != 2018)$model_value*100000,
                                           ymax = subset(out_mat_ci_upper[[1]]$mapped_output$globocan_hcc_incidence,
                                                         time != 2018)$model_value*100000),
                                       col = "black", width = 0.2) +
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
                               legend.margin = margin(t = 0, unit="cm")))

  ## GBD HBV-related cirrhosis mortality rate
  p_gbd <- print(ggplot(data = out_mat[[i]]$mapped_output$gbd_cirrhosis_mortality) +
                   geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000, fill = "Model")) +
                   geom_errorbar(aes(x = paste(age_min,"-",age_max),
                                     ymin = out_mat_ci_lower[[i]]$mapped_output$gbd_cirrhosis_mortality$model_value*100000,
                                     ymax = out_mat_ci_upper[[i]]$mapped_output$gbd_cirrhosis_mortality$model_value*100000),
                                 col = "black", width = 0.2) +
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
                         legend.margin = margin(t = 0, unit="cm")))

  # Demographic characteristics of HBV-related liver disease patients ----
  # Proportion male
  plot_ld_prop_male <-
    ggplot(data = subset(out_mat[[i]]$mapped_output$mapped_liver_disease_demography,
                         grepl("prop_male$",outcome))) +
    geom_col(aes(x = gsub("_.*$","",outcome), y = model_value)) +
    geom_errorbar(aes(x = gsub("_.*$","",outcome),
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$mapped_liver_disease_demography,
                                    grepl("prop_male$",outcome))$model_value,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$mapped_liver_disease_demography,
                                    grepl("prop_male$",outcome))$model_value),
                  col = "black", width = 0.2) +
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
    geom_errorbar(aes(x = gsub("_.*$","",outcome),
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$mapped_liver_disease_demography,
                                    grepl("mean_age$",outcome))$model_value,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$mapped_liver_disease_demography,
                                    grepl("mean_age$",outcome))$model_value),
                  col = "black", width = 0.2) +
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

  # Risk of chronic carriage ----

  ## Risk of chronic carriage: change to author and year and add caption which ones are from Edmunds
  p_p_chronic <- print(ggplot(data = out_mat[[i]]$mapped_output$risk_of_chronic_carriage) +
                         geom_line(aes(x = age, y = model_value, group = "Model", linetype = "Model")) +
                         geom_line(data = out_mat_ci_lower[[i]]$mapped_output$risk_of_chronic_carriage,
                                   aes(x = age, y = model_value), linetype = "dashed") +
                         geom_line(data = out_mat_ci_upper[[i]]$mapped_output$risk_of_chronic_carriage,
                                   aes(x = age, y = model_value), linetype = "dashed") +
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

  ## Mortality curves ----
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
                           geom_errorbar(aes(x = time_interval_years,
                                             ymin = rbind(out_mat_ci_lower[[i]]$mapped_output$mortality_curves,
                                                          mortality_curves_zeros)$model_value,
                                             ymax = rbind(out_mat_ci_upper[[i]]$mapped_output$mortality_curves,
                                                          mortality_curves_zeros)$model_value),
                                         col = "grey", width = 0.2) +
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

  ## ORs ----
  p_or <- print(ggplot(data = out_mat[[i]]$mapped_output$odds_ratios) +
                  geom_col(aes(x = gsub("odds_ratio_", "", outcome), y = log(model_value),
                               fill = "Model")) +
                  geom_errorbar(aes(x = gsub("odds_ratio_", "", outcome),
                                    ymin = log(out_mat_ci_lower[[i]]$mapped_output$odds_ratios$model_value),
                                    ymax = log(out_mat_ci_upper[[i]]$mapped_output$odds_ratios$model_value)),
                                col = "black", width = 0.2) +
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

  # Natural history prevalence plots ----
  out_mat[[i]]$mapped_output$nat_hist_prevalence$model_num <-
    gsub(".*[[:digit:]]{4}_", "",out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique)
  out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence$model_num <-
    gsub(".*[[:digit:]]{4}_", "",out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence$id_unique)
  out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence$model_num <-
    gsub(".*[[:digit:]]{4}_", "",out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence$id_unique)

  # GMB1 PROLIFICA plots: infection phase in chronic carriers
  gmb1_facet_labels <- c("Male blood donors", "Community screening pop.")
  names(gmb1_facet_labels) <- c("Male", "Mixed")

  plot1_gmb1 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                     id_paper == "GMB1" &
                                       model_num != "cc_dcc" & model_num != "hcc"),
                       aes(x = model_num)) +
    geom_col(aes(y = model_value))+
    geom_errorbar(aes(x = model_num,
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "GMB1" &
                                      model_num != "cc_dcc" & model_num != "hcc")$model_value,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "GMB1" &
                                      model_num != "cc_dcc" & model_num != "hcc")$model_value),
                  col = "black", width = 0.2) +
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
    geom_errorbar(aes(ymin = subset(out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "GMB1" &
                                      (model_num == "cc_dcc" | model_num == "hcc"))$model_value,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "GMB1" &
                                      (model_num == "cc_dcc" | model_num == "hcc"))$model_value),
                  col = "black", width = 0.2) +
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
          axis.title.y = element_text(size = 8))

  # 1-1 plots: infection phase in chronic carriers without liver disease
  study_1_facet_labels <- c("1986: median age 11 years", "2013: median age 38 years")
  names(study_1_facet_labels) <- c(1986, 2013)

  plot1_1 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                  grepl(".*it,_ir,_ic_and_enchb$", out_mat[[i]]$mapped_output$nat_hist_prevalence$outcome)),
                    aes(x = toupper(model_num))) +
    geom_col(aes(y = model_value))+
    geom_errorbar(aes(ymin = subset(out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence,
                                    grepl(".*it,_ir,_ic_and_enchb$",
                                          out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence$outcome))$model_value,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence,
                                    grepl(".*it,_ir,_ic_and_enchb$",
                                          out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence$outcome))$model_value),
                  col = "black", width = 0.2) +
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
    geom_errorbar(aes(ymin = subset(out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "1" &
                                      (outcome == "hcc_prevalence_in_chronic_carriers" |
                                         outcome == "cc_and_dcc_prevalence_in_chronic_carriers"))$model_value,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "1" &
                                      (outcome == "hcc_prevalence_in_chronic_carriers" |
                                         outcome == "cc_and_dcc_prevalence_in_chronic_carriers"))$model_value),
                  col = "black", width = 0.2) +
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
          axis.title.y = element_text(size = 8))

  # 1-1 plots: chronic carriers by age
  plot3_1 <- ggplot(data = out_mat[[i]]$mapped_output$nat_hist_prevalence[
    out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_2013_ir_enchb_cc_dcc",]) +
    geom_col(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = model_value, fill = "Model"))+
    geom_errorbar(aes(x = reorder(paste(age_min,"-",age_max), age_min),
                      ymin = out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence[
                        out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_2013_ir_enchb_cc_dcc",]$model_value,
                      ymax = out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence[
                        out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_2013_ir_enchb_cc_dcc",]$model_value),
                  col = "black", width = 0.2) +
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
          plot.caption = element_text(hjust = 0, size = 8))

  ## Natural history prevalence PLOT 1
  p_nat_hist_prev1 <- grid.arrange(plot1_gmb1, plot1_1, plot2_gmb1, plot2_1,
                                   plot3_1, nrow = 3,
                                   layout_matrix = rbind(c(1,2), c(3,4), c(5,5)),
                                   top = "Prevalence measures in chronic carriers")

  # A4 Proportion of deaths from DCC and HCC in comp. cirrhosis cohort
  plot_nat_hist_a4 <- ggplot(data = subset(out_mat[[i]]$mapped_output$nat_hist_prevalence,
                                           id_paper == "A4")) +
    geom_col(aes(x = model_num, y = model_value, fill = "Model"))+
    geom_errorbar(aes(x = model_num,
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "A4")$model_value,
                      ymax =  subset(out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence,
                                     id_paper == "A4")$model_value),
                  col = "black", width = 0.2) +
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
    geom_errorbar(aes(x = model_num,
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence,
                                    id_paper == "GMB2")$model_value,
                      ymax =  subset(out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence,
                                     id_paper == "GMB2")$model_value),
                  col = "black", width = 0.2) +
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
    geom_errorbar(aes(x = reorder(paste(age_min,"-",age_max), age_min),
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence,
                                    outcome == "hbeag_prevalence_in_hcc" |
                                      outcome == "hbeag_prevalence_in_cirrhosis")$model_value,
                      ymax =subset(out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence,
                                   outcome == "hbeag_prevalence_in_hcc" |
                                     outcome == "hbeag_prevalence_in_cirrhosis")$model_value),
                  col = "black", width = 0.2) +
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

  # Vertical transmission plots ----
  # 1-1 plots: chronic infections due to vertical transmission
  plot_nat_hist_1 <- ggplot(data = out_mat[[i]]$mapped_output$nat_hist_prevalence[
    out_mat[[i]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_1986_incident_chronic_births",]) +
    geom_col(aes(x = model_num, y = model_value))+
    geom_errorbar(aes(x = model_num,
                      ymin = out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence[
                        out_mat_ci_lower[[i]]$mapped_output$nat_hist_prevalence$id_unique ==
                          "id_1_1_1986_incident_chronic_births",]$model_value,
                      ymax =out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence[
                        out_mat_ci_upper[[i]]$mapped_output$nat_hist_prevalence$id_unique ==
                          "id_1_1_1986_incident_chronic_births",]$model_value),
                  col = "black", width = 0.2) +
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
    geom_errorbar(aes(ymin = out_mat_ci_lower[[i]]$mapped_output$mtct_risk$model_value,
                      ymax = out_mat_ci_upper[[i]]$mapped_output$mtct_risk$model_value),
                  col = "black", width = 0.2) +
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

  ## Progression rates ----
  out_mat[[i]]$mapped_output$progression_rates$type <-
    gsub("shadow[[:digit:]].{0,1}_", "", out_mat[[i]]$mapped_output$progression_rates$outcome)
  out_mat[[i]]$mapped_output$progression_rates$type <- gsub("_.$", "", out_mat[[i]]$mapped_output$progression_rates$type)

  out_mat_ci_lower[[i]]$mapped_output$progression_rates$type <-
    gsub("shadow[[:digit:]].{0,1}_", "", out_mat_ci_lower[[i]]$mapped_output$progression_rates$outcome)
  out_mat_ci_lower[[i]]$mapped_output$progression_rates$type <- gsub("_.$", "", out_mat_ci_lower[[i]]$mapped_output$progression_rates$type)

  out_mat_ci_upper[[i]]$mapped_output$progression_rates$type <-
    gsub("shadow[[:digit:]].{0,1}_", "", out_mat_ci_upper[[i]]$mapped_output$progression_rates$outcome)
  out_mat_ci_upper[[i]]$mapped_output$progression_rates$type <- gsub("_.$", "", out_mat_ci_upper[[i]]$mapped_output$progression_rates$type)


  # Study 1: HCC incidence
  plot_1_hcc_incidence <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                               id_paper == "1" & type == "hcc_incidence")) +
    geom_col(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = model_value*100000))+
    geom_errorbar(aes(x = paste(bl_age_min_years,"-",bl_age_max_years),
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$progression_rates,
                                    id_paper == "1" & type == "hcc_incidence")$model_value*100000,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$progression_rates,
                                    id_paper == "1" & type == "hcc_incidence")$model_value*100000),
                  col = "black", width = 0.2) +
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
    geom_errorbar(aes(x = outcome,
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$progression_rates,
                                    id_paper == "1" & type == "dcc_incidence")$model_value*100000,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$progression_rates,
                                    id_paper == "1" & type == "dcc_incidence")$model_value*100000),
                  col = "black", width = 0.2) +
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
          axis.text.x = element_blank())


  # Study 1: Mortality
  plot_1_mortality <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                           id_paper == "1" & type == "mortality")) +
    geom_col(aes(x = sex, y = model_value*100000))+
    geom_errorbar(aes(x = sex,
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$progression_rates,
                                    id_paper == "1" & type == "mortality")$model_value*100000,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$progression_rates,
                                    id_paper == "1" & type == "mortality")$model_value*100000),
                  col = "black", width = 0.2) +
    geom_point(aes(x = sex, y = data_value*100000),
               col = "red", shape = 4, size = 3, stroke = 2) +
    geom_errorbar(aes(x = sex, ymax = ci_upper*100000, ymin = ci_lower*100000),
                  col = "red", width = 0.1)  +
    labs(title = "All-cause mortality rate in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Deaths per 100000 PY", x = "Sex") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 8))

  # Study A6: Mortality in cirrhosis cohort
  plot_a6_mortality <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                            id_paper == "A6")) +
    geom_col(aes(x = outcome, y = model_value*100, fill = "Model"))+
    geom_errorbar(aes(x = outcome,
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$progression_rates,
                                    id_paper == "A6")$model_value*100,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$progression_rates,
                                    id_paper == "A6")$model_value*100),
                  col = "black", width = 0.2) +
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
          legend.position = "bottom")

  ## Combined liver disease incidence and mortality rates
  p_prog_rates1 <- grid.arrange(plot_1_hcc_incidence, plot_1_dcc_incidence,
                                plot_1_mortality, plot_a6_mortality, nrow = 2, widths = 4:3,
                                top = "Liver disease-related rates")

  # Study 1: HBeAg loss
  plot_1_eag_loss <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                          id_paper == "1" & type == "eag_loss")) +
    geom_col(aes(x = sex, y = model_value*100))+
    geom_errorbar(aes(x = sex,
                      ymin =subset(out_mat_ci_lower[[i]]$mapped_output$progression_rates,
                                   id_paper == "1" & type == "eag_loss")$model_value*100,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$progression_rates,
                                    id_paper == "1" & type == "eag_loss")$model_value*100),
                  col = "black", width = 0.2) +
    geom_point(aes(x = sex, y = data_value*100),
               shape = 4, size = 3, stroke = 2, col = "red") +
    geom_errorbar(aes(x = sex, ymax = ci_upper*100, ymin = ci_lower*100),
                  col = "red", width = 0.1)  +
    labs(title = "Rate of HBeAg loss in chronic carriers",
         subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
         y = "Cases per 100 PY", x = "Sex") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8))

  # Study 6: HBsAg loss
  plot_6_sag_loss <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                          id_paper == "6")) +
    geom_col(aes(x = outcome, y = model_value*100, fill = "Model"))+
    geom_errorbar(aes(x = outcome,
                      ymin =subset(out_mat_ci_lower[[i]]$mapped_output$progression_rates,
                                   id_paper == "6")$model_value*100,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$progression_rates,
                                    id_paper == "6")$model_value*100),
                  col = "black", width = 0.2) +
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
          axis.text.x = element_blank())

  ## Combined seromarker loss rates
  p_prog_rates2 <- grid.arrange(plot_1_eag_loss, plot_6_sag_loss, nrow = 1, widths = 2:1,
                                top = "Seromarker clearance rates")

  # Transmission-related data from GMB6 and GMB7 ----
  plot_horizontal_transmission <- ggplot(data = subset(out_mat[[i]]$mapped_output$progression_rates,
                                                       id_paper == "GMB6" | id_paper == "GMB7")) +
    geom_col(aes(x = gsub(".*_","",outcome), y = model_value))+
    geom_errorbar(aes(x = gsub(".*_","",outcome),
                      ymin = subset(out_mat_ci_lower[[i]]$mapped_output$progression_rates,
                                    id_paper == "GMB6" | id_paper == "GMB7")$model_value,
                      ymax = subset(out_mat_ci_upper[[i]]$mapped_output$progression_rates,
                                    id_paper == "GMB6" | id_paper == "GMB7")$model_value),
                  col = "black", width = 0.2) +
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
    theme(plot.title = element_text(hjust = 0.5))


  ## Combined transmission-related plot
  p_transmission_rates <- grid.arrange(plot_mtct, plot_nat_hist_1, plot_horizontal_transmission,
                                       layout_matrix = rbind(c(1,1),
                                                             c(2,3)),
                                       top = "Transmission-related measures")

  # List of all plots
  plot_list[[i]] <- list(p_hbsag1, p_antihbc, p_hbeag,
                         p_globocan1, p_globocan2,
                         p_gbd, p_ld_demog,
                         p_p_chronic, p_mort_curves, p_or,
                         p_nat_hist_prev1, p_nat_hist_prev2,
                         p_prog_rates1, p_prog_rates2, p_transmission_rates)

}
dev.off()

### Simulate validation outcomes with best parmsets ----
# Calculate validation outcomes
validation_out_mat <- apply(params_mat_targets5,1,
                            function(x)
                              simulate_validation_outcomes(default_parameter_list = parameter_list,
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
#save(validation_out_mat, file = here("calibration", "output", "validation_output_123_180520.Rdata"))


# HCC mortality rate in 2018
quantile(sapply(validation_out_mat, "[[", "total_hcc_deaths_2018_f")/
           sapply(validation_out_mat, "[[", "pop_female_2018point5"), probs = c(0.05, 0.5, 0.95))*100000
quantile(sapply(validation_out_mat, "[[", "total_hcc_deaths_2018_m")/
           sapply(validation_out_mat, "[[", "pop_male_2018point5"), probs = c(0.05, 0.5, 0.95))*100000

quantile((sapply(validation_out_mat, "[[", "total_hcc_deaths_2018_f")+sapply(validation_out_mat, "[[", "total_hcc_deaths_2018_m"))/
           (sapply(validation_out_mat, "[[", "pop_female_2018point5")+sapply(validation_out_mat, "[[", "pop_male_2018point5")), probs = c(0.05, 0.5, 0.95))*100000

# GBD cirrhosis mortality rate in 2018
quantile(sapply(validation_out_mat, "[[", "total_cirrhosis_deaths_2018_f")/
           sapply(validation_out_mat, "[[", "pop_female_2018point5"), probs = c(0.05, 0.5, 0.95))*100000
quantile(sapply(validation_out_mat, "[[", "total_cirrhosis_deaths_2018_m")/
           sapply(validation_out_mat, "[[", "pop_male_2018point5"), probs = c(0.05, 0.5, 0.95))*100000

quantile((sapply(validation_out_mat, "[[", "total_cirrhosis_deaths_2018_f")+sapply(validation_out_mat, "[[", "total_cirrhosis_deaths_2018_m"))/
           (sapply(validation_out_mat, "[[", "pop_female_2018point5")+sapply(validation_out_mat, "[[", "pop_male_2018point5")), probs = c(0.05, 0.5, 0.95))*100000

plot(x = ages,
     y = (apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_chb_2019_m")),2,median)+
            apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_cc_2019_m")),2,median)+
            apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_dcc_2019_m")),2,median))/
       apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m")),2,median),
     ylab = "")
lines(x = ages,
      y = (apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_chb_2019_f")),2,median)+
             apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_cc_2019_f")),2,median)+
             apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_dcc_2019_f")),2,median))/
        apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_f")),2,median))


prop_eligible_m <- (do.call("rbind", lapply(validation_out_mat, "[[", "number_chb_2019_m"))+
                      do.call("rbind", lapply(validation_out_mat, "[[", "number_cc_2019_m"))+
                      do.call("rbind", lapply(validation_out_mat, "[[", "number_dcc_2019_m")))/
  do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m"))

prop_eligible_f <- (do.call("rbind", lapply(validation_out_mat, "[[", "number_chb_2019_f"))+
                      do.call("rbind", lapply(validation_out_mat, "[[", "number_cc_2019_f"))+
                      do.call("rbind", lapply(validation_out_mat, "[[", "number_dcc_2019_f")))/
  do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_f"))

num_eligible_m <- (do.call("rbind", lapply(validation_out_mat, "[[", "number_chb_2019_m"))+
                     do.call("rbind", lapply(validation_out_mat, "[[", "number_cc_2019_m"))+
                     do.call("rbind", lapply(validation_out_mat, "[[", "number_dcc_2019_m")))

num_eligible_f <- (do.call("rbind", lapply(validation_out_mat, "[[", "number_chb_2019_f"))+
                     do.call("rbind", lapply(validation_out_mat, "[[", "number_cc_2019_f"))+
                     do.call("rbind", lapply(validation_out_mat, "[[", "number_dcc_2019_f")))

plot(ages, y=apply(prop_eligible_m,2,median), ylim = c(0,1), ylab = "", type = "l")
lines(ages, y=apply(prop_eligible_m,2,quantile, probs = 0.05), col = "blue")
lines(ages, y=apply(prop_eligible_m,2,quantile, probs = 0.95), col = "blue")
lines(ages, y=apply(prop_eligible_f,2,median), col = "red")
lines(ages, y=apply(prop_eligible_f,2,quantile, probs = 0.05), col = "pink")
lines(ages, y=apply(prop_eligible_f,2,quantile, probs = 0.95), col = "pink")

# Adult carriers treatment eligibility (18+)
quantile((apply(num_eligible_f[,which(ages == 18):200],1,sum)+
            apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_it_2019_f"))[,which(ages == 30):200],1,sum))/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_f"))[,which(ages == 18):200],1,sum), probs = c(0.05,0.5,0.95))

quantile(apply(num_eligible_f[,which(ages == 18):200],1,sum)/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_f"))[,which(ages == 18):200],1,sum), probs = c(0.05,0.5,0.95))

quantile((apply(num_eligible_m[,which(ages == 18):200],1,sum)+
            apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_it_2019_m"))[,which(ages == 30):200],1,sum))/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m"))[,which(ages == 18):200],1,sum), probs = c(0.05,0.5,0.95))

quantile(apply(num_eligible_m[,which(ages == 18):200],1,sum)/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m"))[,which(ages == 18):200],1,sum), probs = c(0.05,0.5,0.95))


# Older carriers treatment eligibility (30+)
quantile((apply(num_eligible_f[,which(ages == 30):200],1,sum)+
            apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_it_2019_f"))[,which(ages == 30):200],1,sum))/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_f"))[,which(ages == 30):200],1,sum), probs = c(0.05,0.5,0.95))

quantile((apply(num_eligible_m[,which(ages == 30):200],1,sum)+
            apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_it_2019_m"))[,which(ages == 30):200],1,sum))/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m"))[,which(ages == 30):200],1,sum), probs = c(0.05,0.5,0.95))

# Without IT:
quantile(apply(num_eligible_f[,which(ages == 30):200],1,sum)/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_f"))[,which(ages == 30):200],1,sum), probs = c(0.05,0.5,0.95))
quantile(apply(num_eligible_m[,which(ages == 30):200],1,sum)/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m"))[,which(ages == 30):200],1,sum), probs = c(0.05,0.5,0.95))


# Which compartments mainly contribute?
quantile((apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_hcc_2019_m"))[,which(ages == 30):200],1,sum))/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m"))[,which(ages == 30):200],1,sum), probs = c(0.05,0.5,0.95))

quantile((apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_hcc_2019_f"))[,which(ages == 30):200],1,sum))/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_f"))[,which(ages == 30):200],1,sum), probs = c(0.05,0.5,0.95))

# Men: 11% CC, 0.3% DCC, 3.4% IT, 0.5% IR, 6% ENCHB,
# Women: 7% CC, 0.1% DCC, 2.8% IT, 0.5% IR, 22.5% ENCHB!


# Women of childbearing age
quantile(apply(num_eligible_f[,which(ages == 15):which(ages == 45)],1,sum)/
           apply(do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_f"))[,which(ages == 15):which(ages == 45)],1,sum), probs = c(0.05,0.5,0.95))


median(do.call("rbind", lapply(validation_out_mat, "[[", "number_cc_2019_m"))[,which(ages == 40)]/
         do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m"))[,which(ages == 40)])
median(do.call("rbind", lapply(validation_out_mat, "[[", "number_dcc_2019_m"))[,which(ages == 40)]/
         do.call("rbind", lapply(validation_out_mat, "[[", "number_carriers_2019_m"))[,which(ages == 40)])


### Project into the future ----
tic()
sim <- apply(params_mat_targets5, 1,
                    function(x)
                      run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                                parms_to_change = list(b1 = as.list(x)$b1,
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
                                                       vacc_eff = as.list(x)$vacc_eff),
                                scenario = "vacc"))
toc()

out <- lapply(sim, code_model_output)

# Extract HBsAg prevalence

hbsag_prev <- data.frame(time = out[[1]]$time)
hbsag_prev <- cbind(hbsag_prev,
                    sapply(lapply(out, "[[", "infectioncat_total"), "[[", "carriers")/
                      sapply(lapply(out, "[[", "pop_total"), "[[", "pop_total"))
hbsag_prev$median <- apply(hbsag_prev, 1, median)
hbsag_prev$lower <- apply(hbsag_prev, 1, quantile, prob = 0.025)
hbsag_prev$upper <- apply(hbsag_prev, 1,quantile, prob = 0.975)
hbsag_prev_long <-gather(hbsag_prev, key = "sim", value = "prev", -time)

ggplot(hbsag_prev_long) +
  geom_line(aes(x=time, y = prev*100, group = sim), colour = "grey80") +
  geom_line(data = hbsag_prev, aes(x=time, y = median*100), col = "red") +
  geom_line(data = hbsag_prev, aes(x=time, y = lower*100), col = "blue") +
  geom_line(data = hbsag_prev, aes(x=time, y = upper*100), col = "blue") +
  theme_bw() +
  xlab("Year") +
  ylab("HBsAg prevalence (%)") +
  labs(title = "Status quo projection (maintaining 93% infant vaccine coverage)",
       caption = "Red line = median projection, blue line = 95% credible interval") +
  xlim(c(1980,2080))

# Current prev (2019): 10.8 (6.7-15.7)
# Pre-vaccination prev (1990):  15.6 (11.1-21.0)
# 31% reduction
# Prev in 2030: 9.3 (5.4-13.9)

# WHO prev in 2016: 5.8 (4.7-7.1)
# WHO pre-vacc: 11 (9-13.4)
# 47% reduction

# HBsAg prev in under 5 year olds
hbsag_prev_u5 <- data.frame(time = out[[1]]$time)
hbsag_prev <- cbind(hbsag_prev,
                    sapply(lapply(out, "[[", "infectioncat_total"), "[[", "carriers")/
                      sapply(lapply(out, "[[", "pop_total"), "[[", "pop_total"))
hbsag_prev$median <- apply(hbsag_prev, 1, median)
hbsag_prev$lower <- apply(hbsag_prev, 1, quantile, prob = 0.025)
hbsag_prev$upper <- apply(hbsag_prev, 1,quantile, prob = 0.975)
hbsag_prev_long <-gather(hbsag_prev, key = "sim", value = "prev", -time)
