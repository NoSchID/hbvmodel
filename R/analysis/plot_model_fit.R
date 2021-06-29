# Code for plotting fit to the calibration targets

### Load packages and accepted simulations ----
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
#load(file = here("output", "fits", "best_fits_50_of_100000_wed_domain_weights_210819.Rdata"))  # this was from LSR
#load(here("calibration", "output", "model_fit_output_123_180520.Rdata")) # out_mat
load(here("calibration", "output", "model_fit_output_kmeans_221220.Rdata")) # out_mat
input_out_mat <- out_mat

### Function to calculate 95% confidence intervals for all datasets ----
# Except on mortality curves
# Calculate 95% CIs for proportions using Wilson method (better for small samples)
# For data points that don't have 95% CIs
# For rate datasets, exclude datasets where rate = 0
# Function: not working right now
library(binom)
calculate_95_ci <- function(input_dataset, data_type) {

  if(data_type == "proportion") {


    input_dataset[is.na(input_dataset$data_value)==FALSE,] <-
      filter(input_dataset,
             !(is.na(data_value))) %>%
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

### Calculate average model values ----

#best_fits_out_mat <- out_mat_wed_domain_50
best_fits_out_mat <- out_mat

# Extract mapped output
best_fits_mapped_output <- lapply(best_fits_out_mat, "[[", "mapped_output")

#best_fits_out_mat <- out_mat_wed_domain_50
best_fits_out_mat <- out_mat

# Extract mapped output
best_fits_mapped_output <- lapply(best_fits_out_mat, "[[", "mapped_output")
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

# Prepare a mock out_mat to work on
out_mat <- best_fits_out_mat[[1]]$mapped_output
# Remove model_value column to avoid confusion
#out_mat <- lapply(out_mat, function(x) {x$model_value <- NULL ; x})

# Also average Keneba-Manduar force of infection:
#mapped_output_for_error$progression_rates$data_value[mapped_output_for_error$progression_rates$outcome ==
#                                                       "gmb6_1_a_foi" |
#                                                       mapped_output_for_error$progression_rates$outcome ==
#                                                       "gmb6_1_b_foi"] <- (0.115500*83+0.309500*31)/(83+31)
# NOTE: weights on these are halved later, but left sample size as it is


### Split into different datasets for plotting and recalculate confidence intervals where necessary ----

# Globocan rates
# Remove data ci_lower and ci_upper to avoid confusion
globocan_out_mat <- select(out_mat$globocan_hcc_incidence, -ci_lower, -ci_upper)
# Turn into long format then remove the mock value column
globocan_out_mat_long <- gather(globocan_out_mat, key = "type", value = "value",
                                c("data_value", "model_value"))
#globocan_out_mat_long$ci_lower[globocan_out_mat_long$type == "model_value"] <- NA
#globocan_out_mat_long$ci_upper[globocan_out_mat_long$type == "model_value"] <- NA
globocan_out_mat_long$value <- NULL

globocan_out_mat_long$median <- c(out_mat$globocan_hcc_incidence$data_value,
                                  median_model_values[grepl("^globocan_hcc_incidence.*",names(median_model_values))])
globocan_out_mat_long$ci_lower <- c(out_mat$globocan_hcc_incidence$ci_lower,
                                    ci_lower_model_values[grepl("^globocan_hcc_incidence.*",names(ci_lower_model_values))])
globocan_out_mat_long$ci_upper<- c(out_mat$globocan_hcc_incidence$ci_upper,
                                    ci_upper_model_values[grepl("^globocan_hcc_incidence.*",names(ci_upper_model_values))])


gbd_cirrhosis_out_mat <- select(out_mat$gbd_cirrhosis_mortality, -ci_lower, -ci_upper)
# Turn into long format then remove the mock value column
gbd_cirrhosis_out_mat_long <- gather(gbd_cirrhosis_out_mat, key = "type", value = "value",
                                c("data_value", "model_value"))
#globocan_out_mat_long$ci_lower[globocan_out_mat_long$type == "model_value"] <- NA
#globocan_out_mat_long$ci_upper[globocan_out_mat_long$type == "model_value"] <- NA
gbd_cirrhosis_out_mat_long$value <- NULL

gbd_cirrhosis_out_mat_long$median <- c(out_mat$gbd_cirrhosis_mortality$data_value,
                                  median_model_values[grepl("^gbd_cirrhosis_mortality.*",names(median_model_values))])
gbd_cirrhosis_out_mat_long$ci_lower <- c(out_mat$gbd_cirrhosis_mortality$ci_lower,
                                    ci_lower_model_values[grepl("^gbd_cirrhosis_mortality.*",names(ci_lower_model_values))])
gbd_cirrhosis_out_mat_long$ci_upper<- c(out_mat$gbd_cirrhosis_mortality$ci_upper,
                                   ci_upper_model_values[grepl("^gbd_cirrhosis_mortality.*",names(ci_upper_model_values))])

# Seromarker prevalence
seromarker_out_mat <- out_mat$seromarker_prevalence
seromarker_out_mat$model_median <- median_model_values[grepl("^seromarker_prevalence.*",names(median_model_values))]
seromarker_out_mat$model_ci_lower <- ci_lower_model_values[grepl("^seromarker_prevalence.*",names(ci_lower_model_values))]
seromarker_out_mat$model_ci_upper <- ci_upper_model_values[grepl("^seromarker_prevalence.*",names(ci_upper_model_values))]

# Average Keneba-Manduar values:
seromarker_out_mat <- seromarker_out_mat %>%
  group_by(outcome, id_paper, id_group, time, sex, age) %>%
  mutate(weighted_mean = ifelse(test = (study_link == "KM vaccine cohort" & id_proc != "x"),
                                yes = weighted.mean(data_value, sample_size),
                                no = NA),
         sample_size_sum = ifelse(test = (study_link == "KM vaccine cohort" & id_proc != "x"),
                                  yes = sum(sample_size),
                                  no = NA)) %>%
  mutate(weighted_mean = coalesce(weighted_mean, data_value),
         sample_size_sum = coalesce(sample_size_sum, sample_size))

seromarker_out_mat$data_value <- seromarker_out_mat$weighted_mean
seromarker_out_mat$weighted_mean <- NULL
seromarker_out_mat$sample_size <- seromarker_out_mat$sample_size_sum
seromarker_out_mat$sample_size_sum <- NULL
seromarker_out_mat <- data.frame(seromarker_out_mat)
# Calculate 95% CI
seromarker_out_mat <- calculate_95_ci(seromarker_out_mat, "proportion")

# Risk of chronic carriage
p_chronic_out_mat <- out_mat$risk_of_chronic_carriage
p_chronic_out_mat$model_median <- median_model_values[grepl("^risk_of_chronic_carriage.*",names(median_model_values))]
p_chronic_out_mat$model_ci_lower <- ci_lower_model_values[grepl("^risk_of_chronic_carriage.*",names(ci_lower_model_values))]
p_chronic_out_mat$model_ci_upper <- ci_upper_model_values[grepl("^risk_of_chronic_carriage.*",names(ci_upper_model_values))]

# Natural history
nat_hist_prev_out_mat <- select(out_mat$nat_hist_prevalence, -ci_lower, -ci_upper)
# Turn into long format then remove the mock value column
nat_hist_prev_out_mat  <- gather(nat_hist_prev_out_mat, key = "type", value = "value",
                                c("data_value", "model_value"))
nat_hist_prev_out_mat$value <- NULL
nat_hist_prev_out_mat$median <- c(out_mat$nat_hist_prevalence$data_value,
                                  median_model_values[grepl("^nat_hist_prevalence.*",names(median_model_values))])
nat_hist_prev_out_mat$ci_lower <- c(out_mat$nat_hist_prevalence$ci_lower,
                                    ci_lower_model_values[grepl("^nat_hist_prevalence.*",names(ci_lower_model_values))])
nat_hist_prev_out_mat$ci_upper<- c(out_mat$nat_hist_prevalence$ci_upper,
                                   ci_upper_model_values[grepl("^nat_hist_prevalence.*",names(ci_upper_model_values))])

# Liver disease demography
ld_demog_out_mat <- select(out_mat$mapped_liver_disease_demography, -ci_lower, -ci_upper)
# Turn into long format then remove the mock value column
ld_demog_out_mat<- gather(ld_demog_out_mat, key = "type", value = "value",
                                 c("data_value", "model_value"))
ld_demog_out_mat$value <- NULL
ld_demog_out_mat$median <- c(out_mat$mapped_liver_disease_demography$data_value,
                                  median_model_values[grepl("^mapped_liver_disease_demography.*",names(median_model_values))])
ld_demog_out_mat$ci_lower <- c(out_mat$mapped_liver_disease_demography$ci_lower,
                                    ci_lower_model_values[grepl("^mapped_liver_disease_demography.*",names(ci_lower_model_values))])
ld_demog_out_mat$ci_upper<- c(out_mat$mapped_liver_disease_demography$ci_upper,
                                   ci_upper_model_values[grepl("^mapped_liver_disease_demography.*",names(ci_upper_model_values))])

# Add labels
ld_demog_out_mat$outcome2 <- "Proportion male"
ld_demog_out_mat$outcome2[ld_demog_out_mat$outcome %in%
                            c("hcc_mean_age", "cirrhosis_mean_age")] <- "Mean age"

# Mortality curves
mort_curves_out_mat <- select(out_mat$mortality_curves, -ci_lower, -ci_upper)
# Turn into long format then remove the mock value column
mort_curves_out_mat  <- gather(mort_curves_out_mat, key = "type", value = "value",
                                 c("data_value", "model_value"))
mort_curves_out_mat$value <- NULL
mort_curves_out_mat$median <- c(out_mat$mortality_curves$data_value,
                                  median_model_values[grepl("^mortality_curves.*",names(median_model_values))])
mort_curves_out_mat$ci_lower <- c(out_mat$mortality_curves$ci_lower,
                                ci_lower_model_values[grepl("^mortality_curves.*",names(ci_lower_model_values))])
mort_curves_out_mat$ci_upper <- c(out_mat$mortality_curves$ci_upper,
                                  ci_upper_model_values[grepl("^mortality_curves.*",names(ci_upper_model_values))])


# Progression rates
# Remove data ci_lower and ci_upper to avoid confusion
prog_rates_out_mat <- select(out_mat$progression_rates, -ci_lower, -ci_upper)
# Turn into long format then remove the mock value column
prog_rates_out_mat_long <- gather(prog_rates_out_mat, key = "type", value = "value",
                                c("data_value", "model_value"))
prog_rates_out_mat_long$value <- NULL

prog_rates_out_mat_long$median <- c(out_mat$progression_rates$data_value,
                                  median_model_values[grepl("^progression_rates.*",names(median_model_values))])
prog_rates_out_mat_long$ci_lower <- c(out_mat$progression_rates$ci_lower,
                                    ci_lower_model_values[grepl("^progression_rates.*",names(ci_lower_model_values))])
prog_rates_out_mat_long$ci_upper<- c(out_mat$progression_rates$ci_upper,
                                   ci_upper_model_values[grepl("^progression_rates.*",names(ci_upper_model_values))])


# Plots
# HBsAg prevalence by time, age and sex* ----

# Define study labels

hbsag_studies <- unique(data.frame(time = subset(seromarker_out_mat,
                                                 outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$time,
                                   paper_first_author = subset(seromarker_out_mat,
                                                               outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                   paper_year = subset(seromarker_out_mat,
                                                       outcome == "HBsAg_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                   study_link = subset(seromarker_out_mat,
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

# Plots with ribbons for CIs
#tiff(here("output", "fit_plots", "hbsag_plot.tiff"), height = 6.4, width =10, units = 'in', res=300)
ggplot(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence")) +
  geom_line(aes(x = age, y = model_median, group = sex, linetype = "Model"), size = 1) +
  geom_point(aes(x = age, y = data_value, fill = "Data"), col = "red",
             shape = 4, stroke = 2) +
  geom_ribbon(aes(x=age, ymin=model_ci_lower, ymax=model_ci_upper, group = sex), alpha = 0.1) +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower), col = "red") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "red"), labels = "Observed data") +
  facet_wrap(~ time, ncol = 4) +
  geom_text(size = 3.5, data = hbsag_study_labels,
            mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
  labs(y = "HBsAg prevalence (proportion)", x = "Age (years)",
       colour = "Sex:") +
  theme_classic() +
  xlim(0,80) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3)) +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15))
#dev.off()


# LSR plot:
#tiff(here("output", "fits", "lsr_plots", "hbsag_plot2.tiff"), height = 15, width =10, units = 'in', res=300)
ggplot(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence")) +
  geom_line(aes(x = age, y = model_median, linetype = "Model", colour = sex), size = 1) +
  geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
             shape = 4, stroke = 2) +
  geom_line(aes(x = age, y = model_ci_lower, colour = sex), linetype = "dashed", size = 1) +
  geom_line(aes(x = age, y = model_ci_upper, colour = sex), linetype = "dashed", size = 1) +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, group = sex, colour = sex)) +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "black"), labels = "Observed data") +
  scale_colour_manual(values = c("Mixed" = "gray35", "Male" = "navyblue", "Female" = "steelblue"),
                      labels = c("Mixed" = "Both sexes", "Male" = "Male", "Female" = "Female")) +
  facet_wrap(~ time, ncol = 2) +
  geom_text(size = 3.5, data = hbsag_study_labels,
            mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
  labs(title = "HBsAg prevalence in The Gambia",
       y = "HBsAg prevalence (proportion)", x = "Age (years)",
       colour = "Sex:") +
  theme_classic() +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3)) +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15))
#dev.off()

# THESIS PLOT
hbsag_study_labels_floor <- hbsag_study_labels
hbsag_study_labels_floor$label[hbsag_study_labels_floor$label=="Thursz 1995, Bellamy 1998"] <-
  "Chotard 1992, Thursz 1995,\nBellamy 1998"
hbsag_study_labels_floor <- hbsag_study_labels_floor[
  hbsag_study_labels_floor$label!="Chotard 1992",]

#tiff(here("output", "fit_plots", "hbsag_plot.tiff"), height = 6.4, width =10, units = 'in', res=300)
sero_p1 <- ggplot(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence")) +
  geom_line(aes(x = age, y = model_median*100, group = sex, linetype = "Model"), size = 1) +
  geom_point(aes(x = age, y = data_value*100, fill = "Data"), col = "red",
             shape = 4, stroke = 2) +
  geom_ribbon(aes(x=age, ymin=model_ci_lower*100, ymax=model_ci_upper*100, group = sex), alpha = 0.1) +
  geom_errorbar(aes(x = age, ymax = ci_upper*100, ymin = ci_lower*100), col = "red") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "red"), labels = "Observed data") +
  facet_wrap(~ floor(time), ncol = 3) +
  geom_text(size = 3.5, data = hbsag_study_labels_floor,
            mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
  labs(y = "HBsAg prevalence (%)", x = "Age (years)") +
  theme_classic() +
  xlim(0,80) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3)) +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15))
#dev.off()


# Anti-HBc prevalence by time and age* ----
anti_hbc_studies <- unique(data.frame(time = subset(seromarker_out_mat,
                                                    outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$time,
                                      paper_first_author = subset(seromarker_out_mat,
                                                                  outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                      paper_year = subset(seromarker_out_mat,
                                                          outcome == "Anti_HBc_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                      study_link = subset(seromarker_out_mat,
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
antihbc_study_labels$label[antihbc_study_labels$label=="Thursz 1995, Bellamy 1998"] <-
  "Thursz 1995,\nBellamy 1998"
antihbc_study_labels$label[1:5] <- paste0("\n", antihbc_study_labels$label[1:5])

# Datapoints here are mixed, male or female populations and model projection is mixed.
#tiff(here("output", "fits", "lsr_plots", "antihbc_plot.tiff"), height = 9.4, width = 12, units = 'in', res=300)
ggplot(data = subset(seromarker_out_mat, outcome == "Anti_HBc_prevalence")) +
  geom_line(aes(x = age, y = model_median, linetype = "Model"), size = 1) +
  geom_point(aes(x = age, y = data_value, fill = "Data"), col = "red",
             shape = 4, stroke = 2) +
  geom_ribbon(aes(x=age, ymin=model_ci_lower, ymax=model_ci_upper), alpha = 0.1) +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower), col = "red") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"),
                        labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "red"), labels = "Observed data") +
#  scale_colour_manual(values = c("Mixed" = "gray35", "Male" = "navyblue", "Female" = "steelblue")) +
  facet_wrap(~ time, ncol = 3) +
  geom_text(size = 3.5, data = antihbc_study_labels,
            mapping = aes(x = Inf, y = 0.4, label = label), hjust=1.05, vjust=1.5) +
  labs(title = "Anti-HBc prevalence in The Gambia",
       y = "Anti-HBc prevalence (proportion)", x = "Age (years)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15)) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3))
#dev.off()

# THESIS PLOT (all in one)
#tiff(here("output", "fits", "lsr_plots", "antihbc_plot.tiff"), height = 9.4, width = 12, units = 'in', res=300)
sero_p2 <- ggplot(data = subset(seromarker_out_mat, outcome == "Anti_HBc_prevalence")) +
  geom_line(aes(x = age, y = model_median*100, linetype = "Model"), size = 1) +
  geom_point(aes(x = age, y = data_value*100, fill = "Data"), col = "red",
             shape = 4, stroke = 2) +
  geom_ribbon(aes(x=age, ymin=model_ci_lower*100, ymax=model_ci_upper*100), alpha = 0.1) +
  geom_errorbar(aes(x = age, ymax = ci_upper*100, ymin = ci_lower*100), col = "red") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"),
                        labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "red"), labels = "Observed data") +
  #  scale_colour_manual(values = c("Mixed" = "gray35", "Male" = "navyblue", "Female" = "steelblue")) +
 # facet_wrap(~ time, ncol = 3) +
 # geom_text(size = 3.5, data = antihbc_study_labels,
 #           mapping = aes(x = Inf, y = 0.4, label = label), hjust=1.05, vjust=1.5) +
  labs(y = "Anti-HBc prevalence (%)", x = "Age (years)") +
  theme_classic() +
  xlim(0,80) +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0, size = 6),
        #legend.margin=margin(t = 0, unit="cm"),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15)) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3))
#dev.off()

# Risk of chronic carriage* ----
#tiff(here("output", "fits", "lsr_plots", "p_chronic_plot.tiff"), height = 9, width = 14, units = 'in', res=300)
ggplot(p_chronic_out_mat) +
  geom_line(aes(x = age, y = model_median, group = "Model",
                linetype = "Model"), size = 1) +
  geom_ribbon(aes(x=age, ymin=model_ci_lower, ymax=model_ci_upper), alpha = 0.1) +
  geom_point(aes(x = age, y = data_value, fill = "Data"),
             shape = 4, stroke = 2, colour = "red") +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower), col = "red") +
  scale_linetype_manual(name = "", values = c("Model" = "solid"), labels = c("Model" = "Model projection")) +
  scale_fill_manual(name = NULL, values = c("Data" = "red"), labels = "Observed data") +
  labs(title = "Risk of chronic carriage by age at infection",
       y = "Risk of chronic carriage (proportion)", x = "Age at infection (years)",
       fill = "") +
  theme_classic() +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3)) +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15)) +
  ylim(0,1) +
  xlim(0,30)
#dev.off()

# THESIS PLOT
#tiff(here("output", "fits", "lsr_plots", "p_chronic_plot.tiff"), height = 9, width = 14, units = 'in', res=300)
sero_p3 <- ggplot(p_chronic_out_mat) +
  geom_line(aes(x = age, y = model_median*100, group = "Model",
                linetype = "Model"), size = 1) +
  geom_ribbon(aes(x=age, ymin=model_ci_lower*100, ymax=model_ci_upper*100), alpha = 0.1) +
  geom_point(aes(x = age, y = data_value*100, fill = "Data"),
             shape = 4, stroke = 2, colour = "red") +
  geom_errorbar(aes(x = age, ymax = ci_upper*100, ymin = ci_lower*100), col = "red") +
  scale_linetype_manual(name = "", values = c("Model" = "solid"), labels = c("Model" = "Model projection")) +
  scale_fill_manual(name = NULL, values = c("Data" = "red"), labels = "Observed data") +
  labs(y = "Risk of chronic carriage (%)", x = "Age at infection (years)",
       fill = "") +
  theme_classic() +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3)) +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15)) +
  xlim(0,30)
#dev.off()

# HBeAg prevalence by time and age and overall* ----

# Define study labels
hbeag_studies <- unique(data.frame(time = subset(seromarker_out_mat,
                                                 outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$time,
                                   paper_first_author = subset(seromarker_out_mat,
                                                               outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$paper_first_author,
                                   paper_year = subset(seromarker_out_mat,
                                                       outcome == "HBeAg_prevalence" & is.na(data_value) == FALSE)$paper_year,
                                   study_link = subset(seromarker_out_mat,
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

# Seminar plot: HBeAg prevalence at one timepoint (no distinction between studies)
#tiff(here("output", "fit_plots", "hbeag_plot.tiff"), height = 8, width = 8, units = 'in', res=300)
ggplot(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence")) +
  geom_line(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
            aes(x = age,  y = model_median, linetype = "Model"), size = 1) +
  geom_point(aes(x = age, y = data_value, fill = "Data"), col = "red2",
             shape = 4, stroke = 2) +
  geom_ribbon(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
              aes(x=age, ymin=model_ci_lower, ymax=model_ci_upper, group = sex), alpha = 0.1) +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower), col = "red2") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "red"), labels = "Observed data") +
  labs(title = "HBeAg prevalence in chronic carriers",
       y = "HBeAg prevalence (proportion)", x = "Age (years)",
       fill = "") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3))
#dev.off()

# LSR plots:
ggplot(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence")) +
  geom_line(aes(x = age, y = model_median, linetype = "Model", colour = sex)) +
  geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
             shape = 4, stroke = 1.5) +
  geom_line(aes(x = age, y = model_ci_lower, colour = sex), linetype = "dashed") +
  geom_line(aes(x = age, y = model_ci_upper, colour = sex), linetype = "dashed") +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "black"), labels = "Observed prevalence") +
  facet_wrap(~ time, ncol = 3) +
  geom_text(size = 3, data = hbeag_study_labels,
            mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
  labs(title = "HBeAg prevalence over time and by age",
       y = "HBeAg prevalence (proportion)", x = "Age (years)",
       colour = "Sex",
       caption = "Keneba Manduar cohort: Whittle studies, Mendy 1999 & 2008, Van der Sande 2006, Shimakawa 2016 |\nGHIS: Chotard 1992, Fortuin 1993, Whittle 1995, Mendy 1999, Peto 2014 | GLCS: Mendy 2010 | PROLIFICA: Lemoine 2016") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm")) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3))

# HBeAg prevalence at one timepoint
ggplot(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence")) +
  geom_line(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
            aes(x = age,  y = model_median, linetype = "Model")) +
  geom_line(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
            aes(x = age, y = model_ci_lower, linetype = "Percentile")) +
  geom_line(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
            aes(x = age, y = model_ci_upper), linetype = "dashed") +
  geom_point(aes(x = age, y = data_value, colour = study_link),
             shape = 4, stroke = 1.5) +
 # geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = study_link)) +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid", "Percentile" = "dashed"),
                        labels = c("Median model projection", "Model 5th and 95th percentile")) +
  scale_colour_discrete(name = "Observed prevalence (studies):", labels = c("GHIS", "GLCS (subset)",
                                                 "Keneba Manduar chronic carrier cohort",
                                                 "Keneba Manduar chronic carrier cohort+GHIS",
                                                 "Keneba Manduar (others)",
                                                 "PROLIFICA", "Other studies")) +
  labs(title = "HBeAg prevalence over time and by age",
       y = "HBeAg prevalence (proportion)", x = "Age (years)",
       caption = "Keneba Manduar studies: Whittle studies, Mendy 1999 & 2008, Van der Sande 2006, Shimakawa 2016 | GHIS: Chotard 1992, Fortuin 1993, Whittle 1995, Mendy 1999, Peto 2014 |\nGLCS: Mendy 2010 | PROLIFICA: Lemoine 2016 | Others: Ryder 1984, Ryder 1992, Mayans 1990") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm")) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3))

# THESIS PLOT
nat_hist_p1 <- ggplot(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence")) +
  geom_line(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
            aes(x = age,  y = model_median*100, linetype = "Model"), size = 1) +
  geom_point(aes(x = age, y = data_value*100, fill = "Data"), col = "red2",
             shape = 4, stroke = 2) +
  geom_ribbon(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
              aes(x=age, ymin=model_ci_lower*100, ymax=model_ci_upper*100, group = sex), alpha = 0.1) +
  geom_errorbar(aes(x = age, ymax = ci_upper*100, ymin = ci_lower*100), col = "red2") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "red"), labels = "Observed data") +
  labs(y = "HBeAg prevalence (%)", x = "Age (years)",
       fill = "") +
  theme_classic() +
  xlim(0,80) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 15),
        legend.position = "bottom") +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3))

# Nat hist prev* ----
nat_hist_prev_out_mat$model_num <-
  gsub(".*[[:digit:]]{4}_", "",nat_hist_prev_out_mat$id_unique)

# GMB1 PROLIFICA plots: infection phase in chronic carriers
gmb1_facet_labels <- c("Male blood donors", "Community screening")
names(gmb1_facet_labels) <- c("Male", "Mixed")

# Remove data 95% CI for plot
#nat_hist_prev_out_mat$ci_lower[nat_hist_prev_out_mat$type == "data_value"] <- NA
#nat_hist_prev_out_mat$ci_upper[nat_hist_prev_out_mat$type == "data_value"] <- NA

#tiff(here("output", "fit_plots", "infection_phase_prolifica_plot.tiff"), height = 6.7, width = 12, units = 'in', res=300)
ggplot(data = subset(nat_hist_prev_out_mat,
                     id_paper == "GMB1" &
                       model_num != "cc_dcc" & model_num != "hcc"),
       aes(x = model_num)) +
  geom_col(aes(y = median, group = type, fill= type), position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(~sex, scales = "free", labeller = labeller(sex = gmb1_facet_labels)) +
  scale_x_discrete(breaks=c("ic", "ir_enchb"),
                   labels=c("HBeAg-negative\ninfection", "Chronic\nhepatitis B")) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  labs(title = "Infection phase in chronic carriers",
       subtitle = "PROLIFICA (Lemoine 2016)",
       y = "Prevalence (proportion)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))+
  ylim(0,1)
#dev.off()

# 1-1 plots: infection phase in chronic carriers without liver disease
study_1_facet_labels <- c("1986: median age 11 years", "2013: median age 38 years")
names(study_1_facet_labels) <- c(1986, 2013)
study_1_x_labels <- c("HBeAg+\ninfection", "HBeAg+\nCHB", "HBeAg-\ninfection", "HBeAg-\nCHB")

#tiff(here("output", "fit_plots", "infection_phase_shimakawa_plot.tiff"), height = 6.7, width = 12, units = 'in', res=300)
ggplot(data = subset(nat_hist_prev_out_mat,
                     grepl(".*it,_ir,_ic_and_enchb$", nat_hist_prev_out_mat$outcome)),
       aes(x = toupper(model_num))) +
  geom_col(aes(y = median, group = type, fill= type), position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(~time, scales = "free", labeller = labeller(time = study_1_facet_labels)) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  scale_x_discrete(labels = c("IT" = "HBeAg+\ninfection", "IR" = "HBeAg+\nCHB",
                              "IC" = "HBeAg-\ninfection", "ENCHB" = "HBeAg-\nCHB")) +
  labs(title = "Infection phase in chronic carriers\nwithout liver disease",
       subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
       y = "Prevalence (proportion)", x = "") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))
#dev.off()

# THESIS PLOTS
nat_hist_p2 <- ggplot(data = subset(nat_hist_prev_out_mat,
                                    id_paper == "GMB1" &
                                      model_num != "cc_dcc" & model_num != "hcc"),
                      aes(x = model_num)) +
  geom_col(aes(y = median*100, group = type, fill= type), position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100, ymax = ci_upper*100,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(~sex, scales = "free", labeller = labeller(sex = gmb1_facet_labels)) +
  scale_x_discrete(breaks=c("ic", "ir_enchb"),
                   labels=c("HBeAg-\ninfection", "CHB")) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  labs(y = "Prevalence (%)", x = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  ylim(0,100)

study_1_facet_labels <- c("Chronic carrier cohort\n(median age 11 years)", "Chronic carrier cohort\n(median age 38 years)")
names(study_1_facet_labels) <- c(1986, 2013)

shimakawa_dist <- subset(nat_hist_prev_out_mat,
                         grepl(".*it,_ir,_ic_and_enchb$", nat_hist_prev_out_mat$outcome))
shimakawa_dist$model_num <- factor(shimakawa_dist$model_num)
levels(shimakawa_dist$model_num) <- list("HBeAg+\ninfection"="it",
                                         "HBeAg+\nCHB" = "ir",
                                         "HBeAg-\ninfection" = "ic",
                                         "HBeAg-\nCHB" = "enchb")

nat_hist_p3 <- ggplot(data =shimakawa_dist,
                      aes(x =model_num)) +
  geom_col(aes(y = median*100, group = type, fill= type), position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100, ymax = ci_upper*100,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(~time, scales = "free", labeller = labeller(time = study_1_facet_labels)) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  scale_x_discrete(labels = c("IT" = "HBeAg+\ninfection", "IR" = "HBeAg+\nCHB",
                              "IC" = "HBeAg-\ninfection", "ENCHB" = "HBeAg-\nCHB")) +
  labs(y = "Prevalence (%)", x = "") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))+
  ylim(0,100)

# GLOBOCAN rates* ----
globocan_outcome_facet_labels <- c("HCC case incidence", "HCC mortality")
names(globocan_outcome_facet_labels) <- c("hcc_incidence", "hcc_mortality")
# GLOBOCAN PAF-adjusted cancer incidence and mortality in 2018
# Remove 95%CI for data
#globocan_out_mat_long$ci_lower[globocan_out_mat_long$type == "data_value"] <- NA
#globocan_out_mat_long$ci_upper[globocan_out_mat_long$type == "data_value"] <- NA

# New style:
ggplot(data = globocan_out_mat_long[globocan_out_mat_long$time == 2018 &
                                      globocan_out_mat_long$outcome == "hcc_mortality",],
       aes(x = paste(age_min,"-",age_max))) +
  geom_col(data= subset(globocan_out_mat_long, time == 2018 &
                                        outcome == "hcc_mortality" & type == "data_value"),
           aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_line(data= subset(globocan_out_mat_long, time == 2018 &
                          outcome == "hcc_mortality" & type == "model_value"),
           aes(y = median*100000, group = type)) +
  geom_errorbar(data= subset(globocan_out_mat_long, time == 2018 &
                               outcome == "hcc_mortality" & type == "data_value"),
                aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  geom_ribbon(data= subset(globocan_out_mat_long, time == 2018 &
                             outcome == "hcc_mortality" & type == "model_value"),
              aes(ymin=ci_lower*100000, ymax=ci_upper*100000, group = type,
                  colour = type), alpha = 0.2) +
  facet_grid(~ sex) +
  #  facet_grid(outcome ~ sex,
  #             labeller = labeller(outcome =globocan_outcome_facet_labels)) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  labs(title = "HBV-related HCC mortality rates in 2018",
       y = "Deaths per 100,000", x = "Age group (years)",
       subtitle = "GLOBOCAN estimated rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))

# Seminar plot: HCC mortality
# Could change last x axis label on these
#tiff(here("output", "fit_plots", "globocan_hcc_mort_2018_plot.tiff"), height = 7, width = 10, units = 'in', res=300)
ggplot(data = globocan_out_mat_long[globocan_out_mat_long$time == 2018 & globocan_out_mat_long$outcome == "hcc_mortality",],
       aes(x = paste(age_min,"-",age_max))) +
  geom_col(aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(~ sex) +
#  facet_grid(outcome ~ sex,
#             labeller = labeller(outcome =globocan_outcome_facet_labels)) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  labs(title = "HBV-related HCC mortality rates in 2018",
       y = "Deaths per 100000", x = "Age group (years)",
       subtitle = "GLOBOCAN estimated rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))
#dev.off()

# Seminar plot: HCC incidence
#tiff(here("output", "fit_plots", "globocan_hcc_inc_2018_plot.tiff"), height = 7, width = 10, units = 'in', res=300)
ggplot(data = globocan_out_mat_long[globocan_out_mat_long$time == 2018 & globocan_out_mat_long$outcome == "hcc_incidence",],
       aes(x = paste(age_min,"-",age_max))) +
  geom_col(aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(~ sex) +
  #  facet_grid(outcome ~ sex,
  #             labeller = labeller(outcome =globocan_outcome_facet_labels)) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  labs(title = "HBV-related HCC incidence rates in 2018",
       y = "HCC cases per 100000", x = "Age group (years)",
       subtitle = "GLOBOCAN estimated rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))
#dev.off()

# Seminar plot: Older HCC incidence rates
#tiff(here("output", "fit_plots", "globocan_hcc_inc_1988_1998_plot.tiff"), height = 10, width = 10, units = 'in', res=300)
ggplot(data = globocan_out_mat_long[globocan_out_mat_long$time != 2018,],
       aes(x = paste(age_min,"-",age_max))) +
  geom_col(aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(time ~ sex) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  labs(title = "HBV-related HCC incidence rates in 1988 and 1998",
       y = "HCC cases per 100000", x = "Age group (years)",
       subtitle = "GLOBOCAN estimated rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 12, angle = 45, vjust=0.7),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))
#dev.off()


# LSR plot
#tiff(here("output", "fits", "lsr_plots", "globocan_hcc_2018_plot.tiff"), height = 9, width = 15, units = 'in', res=300)
ggplot(data = globocan_out_mat_long[globocan_out_mat_long$time == 2018,],
       aes(x = paste(age_min,"-",age_max))) +
  geom_col(aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(outcome ~ sex,
             labeller = labeller(outcome =globocan_outcome_facet_labels)) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "steelblue"),
                    labels = c("model_value" = "Median model projection\nwith 5th and 95th percentile\nerror bars",
                               "data_value" = "GLOBOCAN observed rate")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "white")) +
  theme_bw() +
  labs(title = "GLOBOCAN HBV-related HCC incidence and mortality rates in 2018",
       y = "Cases/deaths per 100000 PY", x = "Age (years)",
       subtitle = "GLOBOCAN rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.margin=margin(t = 0, unit="cm"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))
#dev.off()

## GLOBOCAN PAF-adjusted cancer incidence in 1988 and 1998
ggplot(data = globocan_out_mat_long[globocan_out_mat_long$time != 2018,],
       aes(x = paste(age_min,"-",age_max))) +
  geom_col(aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3)+
  facet_grid(time ~ sex) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "coral"),
                    labels = c("Median model projection", "GLOBOCAN observed rate")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "red"),
                      labels = c("Model 5th and 95th percentile", "Observed rate 95% CI")) +
  labs(title = "GLOBOCAN HBV-related HCC incidence rates in 1988 and 1998",
       y = "Cases per 100000 PY", x = "Age (years)",
       subtitle = "GLOBOCAN rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.margin=margin(t = 0, unit="cm"))

# THESIS PLOT: HCC INCIDENCE IN 2018
hcc_inc_2018 <- globocan_out_mat_long[globocan_out_mat_long$time == 2018 &
                                        globocan_out_mat_long$outcome == "hcc_incidence",]
hcc_inc_2018$sex_type <- paste(hcc_inc_2018$sex, hcc_inc_2018$type)
hcc_inc_2018$age_group <- paste0(hcc_inc_2018$age_min,"-",hcc_inc_2018$age_max)
hcc_inc_2018$age_group[hcc_inc_2018$age_group=="70-99.5"] <- "70+"

hcc_cirr_p1 <- ggplot(data = hcc_inc_2018,
                       aes(x = age_group)) +
  geom_line(data= subset(hcc_inc_2018, type == "model_value"),
            aes(y = median*100000, group = sex_type, colour= sex_type),
            size = 1) +
  geom_point(data= subset(hcc_inc_2018, type == "data_value"),
             aes(y = median*100000, group = sex_type, colour= sex_type),
             shape=4, stroke = 2) +
  geom_errorbar(data= subset(hcc_inc_2018, type == "data_value"),
                aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = sex_type, colour= sex_type), width = 0.15,
                show.legend = FALSE)+
  geom_ribbon(data= subset(hcc_inc_2018, type == "model_value"),
              aes(ymin=ci_lower*100000, ymax=ci_upper*100000,
                  group = sex_type, fill = sex_type),
              alpha = 0.1) +
  guides(color = guide_legend(override.aes = list(shape = c(4, NA, 4,NA),
                                                  linetype =c(NA,1,NA,1))),
         fill=FALSE) +
  scale_colour_manual(labels=c("Female data_value" = "Observed data, women",
                               "Female model_value" = "Model output, women",
                               "Male data_value" = "Observed data, men",
                               "Male model_value" = "Model output, men"),
                      values=c("Female data_value" = "red4",
                               "Female model_value" = "black",
                               "Male data_value" = "red",
                               "Male model_value" = "grey50")) +
  scale_fill_manual(values=c("Female data_value" = "red4",
                             "Female model_value" = "black",
                             "Male data_value" = "red",
                             "Male model_value" = "grey60")) +
  #scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  #scale_shape_manual(name = NULL, values = c("Data" = 4), labels = "Observed data") +
  theme_classic() +
  labs(y = "HBV-related HCC cases\nper 100,000", x = "Age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = c(0.3,0.85),
        axis.title = element_text(size = 20),  #14
        axis.text = element_text(size = 15),  #14
        legend.text = element_text(size = 15),  #13
        legend.title = element_blank())

# GBD rates* ----
ggplot(data = gbd_cirrhosis_out_mat_long[gbd_cirrhosis_out_mat_long$time == 2017,],
       aes(x = paste(age_min,"-",age_max))) +
  geom_col(data= subset(gbd_cirrhosis_out_mat_long, time == 2017 &
                          type == "data_value"),
           aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_line(data= subset(gbd_cirrhosis_out_mat_long, time == 2017 &
                            type == "model_value"),
            aes(y = median*100000, group = type)) +
  geom_errorbar(data= subset(gbd_cirrhosis_out_mat_long, time == 2017 &
                               type == "data_value"),
                aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  geom_ribbon(data= subset(gbd_cirrhosis_out_mat_long, time == 2017 &
                              type == "model_value"),
              aes(ymin=ci_lower*100000, ymax=ci_upper*100000, group = type,
                  colour = type), alpha = 0.2) +
  facet_grid(~ sex) +
  #  facet_grid(outcome ~ sex,
  #             labeller = labeller(outcome =globocan_outcome_facet_labels)) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  labs(title = "HBV-related cirrhosis mortality rates in 2018",
       y = "Deaths per 100,000", x = "Age group (years)"
       #subtitle = "GLOBOCAN estimated rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)"
       ) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))

# THESIS PLOT: Cirrhosis mortality 2017
cirr_mort_2017 <- gbd_cirrhosis_out_mat_long[gbd_cirrhosis_out_mat_long$time == 2017,]
cirr_mort_2017$sex_type <- paste(cirr_mort_2017$sex, cirr_mort_2017$type)
cirr_mort_2017$age_group <- paste0(cirr_mort_2017$age_min,"-",floor(cirr_mort_2017$age_max))
cirr_mort_2017$age_group[cirr_mort_2017$age_group=="70-99"] <- "70+"

hcc_cirr_p2 <-ggplot(data = cirr_mort_2017,
       aes(x = age_group)) +
  geom_line(data= subset(cirr_mort_2017, type == "model_value"),
            aes(y = median*100000, group = sex_type, colour= sex_type),
            size = 1) +
  geom_point(data= subset(cirr_mort_2017, type == "data_value"),
             aes(y = median*100000, group = sex_type, colour= sex_type),
             shape=4, stroke = 2) +
  geom_errorbar(data= subset(cirr_mort_2017, type == "data_value"),
                aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = sex_type, colour= sex_type), width = 0.15,
                show.legend = FALSE)+
  geom_ribbon(data= subset(cirr_mort_2017, type == "model_value"),
              aes(ymin=ci_lower*100000, ymax=ci_upper*100000,
                  group = sex_type, fill = sex_type),
              alpha = 0.1) +
  guides(color = guide_legend(override.aes = list(shape = c(4, NA, 4,NA),
                                                  linetype =c(NA,1,NA,1))),
         fill=FALSE) +
  scale_colour_manual(labels=c("Female data_value" = "Observed data, women",
                               "Female model_value" = "Model output, women",
                               "Male data_value" = "Observed data, men",
                               "Male model_value" = "Model output, men"),
                      values=c("Female data_value" = "red4",
                               "Female model_value" = "black",
                               "Male data_value" = "red",
                               "Male model_value" = "grey50")) +
  scale_fill_manual(values=c("Female data_value" = "red4",
                             "Female model_value" = "black",
                             "Male data_value" = "red",
                             "Male model_value" = "grey60")) +
  #scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  #scale_shape_manual(name = NULL, values = c("Data" = 4), labels = "Observed data") +
  theme_classic() +
  labs(y = "HBV-related cirrhosis deaths\nper 100,000", x = "Age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "none",
        axis.title = element_text(size = 20),  #14
        axis.text = element_text(size = 15),  #14
        legend.text = element_text(size = 15),  #13
        legend.title = element_blank())


# Survival curves* ----
## Mortality curves

# Add articifial zeros at first timestep to allow plotting of step curves
mortality_curves_zeros <- mort_curves_out_mat
mortality_curves_zeros$time_interval_years <- 0
mortality_curves_zeros$median <- 0
mortality_curves_zeros$ci_lower <- NA
mortality_curves_zeros$ci_upper <- NA
mortality_curves_zeros$number_at_risk <- mortality_curves_zeros$sample_size
mortality_curves_zeros <- unique(mortality_curves_zeros)

mort_curves_out_mat <- rbind(mort_curves_out_mat, mortality_curves_zeros)

# Add labels for panels with reference and study population
mort_curves_labels <- c("Mortality in compensated\ncirrhosis patients\n(Shimakawa, 2016)",
                        "Mortality in HCC patients\n(Yang, 2017)",
                        "HCC incidence in\ncirrhosis patients\n(Diarra, 2010)",
                        "Mortality in\ncirrhosis patients\n(Diarra, 2010)",
                        "Mortality in HCC patients\n(Bah, 2011)")
names(mort_curves_labels) <- c("shadow4_cum_mortality", "shadow5_cum_mortality",
                               "shadow6_cum_hcc_incidence", "shadow6_cum_mortality",
                               "shadow7_cum_mortality")

mort_curve_cirrhosis <- ggplot() +
  geom_step(data=subset(mort_curves_out_mat, outcome %in% c("shadow4_cum_mortality",
                    "shadow6_cum_hcc_incidence", "shadow6_cum_mortality") &
                    type == "model_value"),
            aes(x = time_interval_years, y = median, linetype = "Model projection")) +
  geom_step(data=subset(mort_curves_out_mat, outcome %in% c("shadow4_cum_mortality",
                                                            "shadow6_cum_hcc_incidence", "shadow6_cum_mortality") &
                          type == "model_value"),
            aes(x = time_interval_years, y = ci_lower), colour="grey", linetype="dashed") +
  geom_step(data=subset(mort_curves_out_mat, outcome %in% c("shadow4_cum_mortality",
                                                            "shadow6_cum_hcc_incidence", "shadow6_cum_mortality") &
                          type == "model_value"),
            aes(x = time_interval_years, y = ci_upper), colour="grey", linetype="dashed") +
  geom_point(data=subset(mort_curves_out_mat, outcome %in% c("shadow4_cum_mortality",
                                                             "shadow6_cum_hcc_incidence", "shadow6_cum_mortality") &
                           type == "data_value" &
                           time_interval_years != 0),
               aes(x = time_interval_years, y = median, colour = "Observed data"),
             shape = 4, stroke = 2) +
  facet_grid(~ outcome,
             labeller = labeller(outcome = mort_curves_labels)) +
  scale_colour_manual(name = "", values = c("Observed data" = "red")) +
  scale_linetype_manual(name = "", values = c("Model projection" = "solid")) +
  labs(y = "Cumulative probability", x = "Follow-up time (years)") +
  theme_classic() +
  ylim(0,1) +
  theme(legend.margin=margin(t = 0, unit="cm"),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 17))

mort_curve_hcc <- ggplot() +
  geom_step(data=subset(mort_curves_out_mat, outcome %in% c("shadow5_cum_mortality",
                                                            "shadow7_cum_mortality") &
                          type == "model_value"),
            aes(x = time_interval_years, y = median, linetype = "Model projection")) +
  geom_step(data=subset(mort_curves_out_mat, outcome %in% c("shadow5_cum_mortality",
                                                            "shadow7_cum_mortality") &
                          type == "model_value"),
            aes(x = time_interval_years, y = ci_lower), colour="grey", linetype="dashed") +
  geom_step(data=subset(mort_curves_out_mat, outcome %in% c("shadow5_cum_mortality",
                                                            "shadow7_cum_mortality") &
                          type == "model_value"),
            aes(x = time_interval_years, y = ci_upper), colour="grey", linetype="dashed") +
  geom_point(data=subset(mort_curves_out_mat, outcome %in% c("shadow5_cum_mortality",
                                                            "shadow7_cum_mortality") &
                           type == "data_value" &
                           time_interval_years != 0),
             aes(x = time_interval_years, y = median, colour = "Observed data"),
             shape = 4, stroke = 2) +
  facet_grid(~ outcome, scales="free_x",
             labeller = labeller(outcome = mort_curves_labels)) +
  scale_colour_manual(name = "", values = c("Observed data" = "red")) +
  scale_linetype_manual(name = "", values = c("Model projection" = "solid")) +
  labs(y = "Cumulative probability", x = "Follow-up time (years)") +
  theme_classic() +
  theme(legend.margin=margin(t = 0, unit="cm"),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "bottom",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 17))

# Proportion of births due to MTCT ----

# 1-1 plots: chronic infections due to vertical transmission
#tiff(here("output", "fit_plots", "proportion_from_mtct_shimakawa_plot.tiff"), height = 7, width = 10, units = 'in', res=300)
ggplot(data = subset(nat_hist_prev_out_mat, id_unique=="id_1_1_1986_incident_chronic_births"),
       aes(x = model_num)) +
  geom_col(aes(y = median, group = type, fill= type), position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  labs(title = "Proportion of chronic infection due to vertical transmission\nin unvaccinated population",
       subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
       y = "Proportion", x = "") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))
#dev.off()

# Shimakawa rates ----

# Subset HCC rates
shimakawa_hcc_rates <- prog_rates_out_mat_long[c(3:6, 17:20),]
# Replace values of 0 with very small value so they are visible on the plot
shimakawa_hcc_rates$median[c(1,3,4)] <- 0.00001

# Seminar plot: HCC incidence
#tiff(here("output", "fit_plots", "shimakawa_hcc_rate_plot.tiff"), height = 6.4, width = 8, units = 'in', res=300)
ggplot(data = shimakawa_hcc_rates,
       aes(x = paste(bl_age_min_years,"-",bl_age_max_years))) +
  geom_col(aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  facet_grid(~sex, scales = "free") +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  labs(title = "HCC incidence rate in chronic carriers",
       subtitle = "Keneba Manduar chronic carrier cohort (Shimakawa 2016)",
       y = "Cases per 100,000 person-years", x = "Baseline age group (years)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))
#dev.off()





# Liver disease demography ----
# Proportion male

plot_ld_mean_age <- ggplot(data = subset(ld_demog_out_mat, outcome2 == "Mean age"),
       aes(x = pop_group_clinical)) +
  geom_col(aes(y = median, group = type, fill= type), position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, group = type, colour= type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
#  facet_wrap(~outcome2, scales = "free_y") +
  scale_x_discrete(breaks=c("HBsAg-positive cirrhosis patients", "HBsAg-positive HCC patients"),
                   labels=c("Cirrhosis", "HCC")) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  guides(fill=FALSE) +
  ylim(0,70) +
  labs(#title = "Characteristics of HBV-related liver disease patients",
       #subtitle = "Gambia Liver Cancer Study (Mendy 2010)",
       y = "Mean age at presentation (years)", x = "Disease outcome") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))

plot_ld_prop_male <-
  ggplot(data = subset(ld_demog_out_mat, outcome2 == "Proportion male"),
         aes(x = pop_group_clinical)) +
  geom_col(aes(y = median, group = type, fill= type), position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, group = type, colour= type), position = position_dodge(width=0.9), width = 0.3,
                show.legend = FALSE)+
  scale_x_discrete(breaks=c("HBsAg-positive cirrhosis patients", "HBsAg-positive HCC patients"),
                   labels=c("Cirrhosis", "HCC")) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "red2"),
                    labels = c("model_value" = "Model projection",
                               "data_value" = "Observed data")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "black")) +
  theme_classic() +
  #guides(fill=FALSE) +
  ylim(0,1) +
  labs(#title = "Characteristics of HBV-related liver disease patients",
       #subtitle = "Gambia Liver Cancer Study (Mendy 2010)",
       y = "Proportion male", x = "Disease outcome") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))


## Combined liver disease demography
p_ld_demog <- grid.arrange(plot_ld_mean_age, plot_ld_prop_male, nrow = 1, widths= 1:2,
                           top = "Characteristics of HBV-related liver disease patients\nGambia Liver Cancer Study (Mendy 2010)")

# PAPER PLOT (sAg prevalence and HCC incidence) ----

# For paper: combine pre-vaccination prevalence plot with post-vaccination ones

# Pre-vaccination prevalence by age
# Plotting only male model projection here but they are very similar
hbsag_p1_facet <- "Pre-vaccination period"
names(hbsag_p1_facet) <- c("HBsAg_prevalence")
# Random choice of outcome to be able to add a facet label

hbsag_p1 <- ggplot(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence" & time < 1990)) +
  geom_line(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence" & time == 1989 & sex == "Male"),
            aes(x = age, y = model_median*100, linetype = "Model"), size = 1) +
  geom_ribbon(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence" & time == 1989 & sex == "Male"),
              aes(x=age, ymin=model_ci_lower*100, ymax=model_ci_upper*100),  fill= "grey50",
              alpha = 0.1) +
  geom_point(aes(x = age, y = data_value*100, fill = "Data"), col = "#458EBF", # green=#55C667FF
             shape = 4, stroke = 2) +   # blue = #31688EFF
  geom_errorbar(aes(x = age, ymax = ci_upper*100, ymin = ci_lower*100), col = "#458EBF") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model output") +
  scale_fill_manual(name = NULL, values = c("Data" = "#458EBF"), labels = "Observed data") +
  labs(y = "HBsAg prevalence (%)", x = "Age (years)") +
  scale_y_continuous(breaks = c(0,20,40)) +
  theme_classic() +
  facet_wrap(~outcome,
             labeller = labeller(outcome =hbsag_p1_facet)) +
  xlim(0,65) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3)) +
  guides(linetype=FALSE, fill = FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.text = element_text(size = 14),
        #axis.title.x = element_text(size = 14),
        axis.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14))


# Post-vaccination age groups affected by vaccination
vacc_age_groups <- subset(seromarker_out_mat, outcome == "HBsAg_prevalence" & time >= 1990 &
                            ((paper_first_author == "Peto") |
                               (paper_first_author == "Lemoine" & time==2013) |
                               (paper_first_author == "Bittaye")))

hbsag_p2_facet <- c("17 years post-vaccination", "23 years post-vaccination", "25 years post-vaccination")
names(hbsag_p2_facet) <- c(2007.5, 2013,2015)

hbsag_p2 <- ggplot(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence" & time %in%
                                   c(2007.5, 2013,2015))) +
  geom_line(aes(x = age, y = model_median*100, group = sex, linetype = "Model"), size = 1) +
  geom_point(data = vacc_age_groups, aes(x = age, y = data_value*100, fill = "Data"),
             col = "#458EBF",shape = 4, stroke = 2) +
  geom_ribbon(aes(x=age, ymin=model_ci_lower*100, ymax=model_ci_upper*100, group = sex),
              fill= "grey50",
              alpha = 0.1) +
  geom_errorbar(data = vacc_age_groups, aes(x = age, ymax = ci_upper*100, ymin = ci_lower*100),
                col = "#458EBF") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model output") +
  scale_fill_manual(name = NULL, values = c("Data" = "#458EBF"), labels = "Observed data") +
  facet_wrap(~ time, ncol = 1, labeller = labeller(time =hbsag_p2_facet)) +
  labs(y = "", x = "Age (years)") +
  scale_y_continuous(breaks = c(0,20,40)) +
  theme_classic() +
  xlim(0,65) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3)) +
  #guides(linetype=FALSE, fill = FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = "bottom",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

#grid.arrange(hbsag_p1, hbsag_p2, ncol=2)

# For other plot
library(grid)
hbsag_p <- grid.arrange(hbsag_p1, hbsag_p2, nrow = 2, heights = c(1,3),
                        left = textGrob("Chronic infection prevalence (%)",
                                        rot = 90,
                                        gp=gpar(fontsize=14)))

# Line chart of 2018 HCC incidence for paper
hcc_inc_2018 <- globocan_out_mat_long[globocan_out_mat_long$time == 2018 &
                                        globocan_out_mat_long$outcome == "hcc_incidence",]
hcc_inc_2018$sex_type <- paste(hcc_inc_2018$sex, hcc_inc_2018$type)
hcc_inc_2018$age_group <- paste0(hcc_inc_2018$age_min,"-",hcc_inc_2018$age_max)
hcc_inc_2018$age_group[hcc_inc_2018$age_group=="70-99.5"] <- "70+"

hcc_inc_plot <- ggplot(data = hcc_inc_2018,
       aes(x = age_group)) +
  geom_line(data= subset(hcc_inc_2018, type == "model_value"),
            aes(y = median*100000, group = sex_type, colour= sex_type),
            size = 1) +
  geom_point(data= subset(hcc_inc_2018, type == "data_value"),
             aes(y = median*100000, group = sex_type, colour= sex_type),
             shape=4, stroke = 2) +
  geom_errorbar(data= subset(hcc_inc_2018, type == "data_value"),
                aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = sex_type, colour= sex_type), width = 0.15,
                show.legend = FALSE)+
  geom_ribbon(data= subset(hcc_inc_2018, type == "model_value"),
              aes(ymin=ci_lower*100000, ymax=ci_upper*100000,
                  group = sex_type, fill = sex_type),
              alpha = 0.1) +
  guides(color = guide_legend(override.aes = list(shape = c(4, NA, 4,NA),
                                                  linetype =c(NA,1,NA,1))),
         fill=FALSE) +
  scale_colour_manual(labels=c("Female data_value" = "Observed data, women",
                               "Female model_value" = "Model output, women",
                               "Male data_value" = "Observed data, men",
                               "Male model_value" = "Model output, men"),
                      values=c("Female data_value" = "#9ecae1",  # purple=#440154
                               "Female model_value" = "#9ecae1",
                               "Male data_value" = "#1E4158",  #blue original = #31688EFF
                               "Male model_value" = "#1E4158")) +
  scale_fill_manual(values=c("Female data_value" = "#9ecae1",  # purple=#440154
                               "Female model_value" = "#9ecae1",  # blue #458EBF"
                               "Male data_value" = "#1E4158",  #blue = #31688EFF
                               "Male model_value" = "#1E4158")) +
  #scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  #scale_shape_manual(name = NULL, values = c("Data" = 4), labels = "Observed data") +
  theme_classic() +
  labs(y = "HBV-related HCC cases\nper 100,000 in 2018", x = "Age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = c(0.25,0.82),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_blank())

# HBeAg plot
hbeag_p <- ggplot(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence")) +
  geom_line(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
            aes(x = age,  y = model_median*100, linetype = "Model"), size = 1) +
  geom_point(aes(x = age, y = data_value*100, fill = "Data"), col = "#458EBF",
             shape = 4, stroke = 2) +
  geom_ribbon(data = subset(seromarker_out_mat, outcome == "HBeAg_prevalence" & time == 1992),
              aes(x=age, ymin=model_ci_lower*100, ymax=model_ci_upper*100, group = sex),
              fill= "grey50", alpha = 0.1) +
  geom_errorbar(aes(x = age, ymax = ci_upper*100, ymin = ci_lower*100),col = "#458EBF") +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model output") +
  scale_fill_manual(name = NULL, values = c("Data" = "#458EBF"), labels = "Observed data") +
  labs(y = "HBeAg prevalence (%)", x = "Age (years)",
       fill = "") +
  theme_classic() +
  xlim(0,80) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3)) +
  #guides(linetype=FALSE, fill = FALSE) +
 theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      legend.margin=margin(t = 0, unit="cm"),
      legend.position = c(0.8,0.8),
      axis.title.x = element_text(size = 14),
      axis.title.y =element_text(size = 14),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      strip.text.x = element_text(size = 14))


# COMBINED PLOT
plot_a <- arrangeGrob(hbsag_p, top = textGrob("A", x = unit(0.01, "npc")
                                              , y   = unit(1, "npc"), just=c("left","top"),
                                              gp=gpar(col="black", fontsize=20)))

plot_b <- arrangeGrob(hbeag_p, top = textGrob("B", x = unit(0.01, "npc")
                                              , y   = unit(1, "npc"), just=c("left","top"),
                                              gp=gpar(col="black", fontsize=20)))
plot_c <- arrangeGrob(hcc_inc_plot, top = textGrob("C", x = unit(0.01, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=20)))

hbeag_hcc_p <- grid.arrange(plot_b, plot_c,nrow=2)

model_fits <- grid.arrange(plot_a,
             hbeag_hcc_p, ncol=2)
#ggsave("model_fits.pdf", plot = model_fits,
#       width= 30, height=20, units = "cm")
#ggsave("model_fits.png", plot = model_fits,
#       width= 30, height=20, units = "cm")

# Previously included corrplot instead of HBeAg (Tim suggested to remove)
# Correlation plot of data and accepted model values
model_values_mat_accepted <- readRDS(file = here("calibration","output", "accepted_model_output_for_calibration_correlation_plot.rds"))

# All values between 0 and 1 (all but 4)
corrplot <- ggplot(model_values_mat_accepted[model_values_mat_accepted$data_value<1,])+
  stat_summary(aes(x=data_value, y = value,
                   colour = reorder(as.factor(weight), desc(as.factor(weight)))),
               fun="median", geom="point") +
  stat_summary(aes(x=data_value, y = value, colour = reorder(as.factor(weight), desc(as.factor(weight)))),
               fun.min= function(x) quantile(x,0.025),
               fun.max= function(x) quantile(x,0.975),
               geom="errorbar") +
  geom_abline(slope=1,intercept=0) +
  scale_colour_viridis_d(end=0.6,
                         labels=c("0.1" = "Lower data quality", "1" = "Higher data quality")) +
  scale_x_continuous(breaks=c(0,0.5,1)) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  ylab("\nModel output") +
  xlab("Observed data") +
  labs(colour="Weight") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.margin=margin(t = 0, unit="cm"),
        legend.position = c(0.25,0.9),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 13.5),
        legend.title = element_blank())


# CALIBRATION PLOTS FOR THESIS ----
library(grid)

# HBsAg, anti-HBC, risk of chronic carriage
sero_p1a <- arrangeGrob(sero_p1, top = textGrob("A", x = unit(0.01, "npc"),
                                                  y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=22)))
sero_p2b <- arrangeGrob(sero_p2, top = textGrob("B", x = unit(0.01, "npc"),
                                                y   = unit(1, "npc"), just=c("left","top"),
                                                gp=gpar(col="black", fontsize=22)))
sero_p3c <- arrangeGrob(sero_p3, top = textGrob("C", x = unit(0.01, "npc"),
                                                y   = unit(1, "npc"), just=c("left","top"),
                                                gp=gpar(col="black", fontsize=22)))


sero_p2_3 <- grid.arrange(sero_p2b, sero_p3c, nrow=2)

#tiff(file = "seromarker_transmission_fit.tiff", width=300, height=230, units = "mm", res=300, pointsize = 0.99)
grid.arrange(sero_p1a, sero_p2_3, ncol = 2, widths=c(3,2))
#dev.off()

# HCC incidence and cirrhosis mortality
hcc_cirr_p1a <- arrangeGrob(hcc_cirr_p1, top = textGrob("A", x = unit(0.01, "npc"),
                                                y   = unit(1, "npc"), just=c("left","top"),
                                                gp=gpar(col="black", fontsize=22)))
hcc_cirr_p2b <- arrangeGrob(hcc_cirr_p2, top = textGrob("B", x = unit(0.01, "npc"),
                                                y   = unit(1, "npc"), just=c("left","top"),
                                                gp=gpar(col="black", fontsize=22)))

#tiff(file = "disease_outcomes_fit.tiff", width=300, height=130, units = "mm", res=300, pointsize = 0.99)
grid.arrange(hcc_cirr_p1a, hcc_cirr_p2b, ncol = 2)
#dev.off()


# HBeAg prevalence and distribution of carrier states
# HBsAg, anti-HBC, risk of chronic carriage
nat_hist_p1a <- arrangeGrob(nat_hist_p1, top = textGrob("A", x = unit(0.01, "npc"),
                                                y   = unit(1, "npc"), just=c("left","top"),
                                                gp=gpar(col="black", fontsize=22)))

nat_hist_p2_3 <- grid.arrange(nat_hist_p2, nat_hist_p3, nrow=2)

nat_hist_p2_3b <- arrangeGrob(nat_hist_p2_3, top = textGrob("B", x = unit(0.01, "npc"),
                                                y   = unit(1, "npc"), just=c("left","top"),
                                                gp=gpar(col="black", fontsize=22)))


#tiff(file = "natural_history_fit.tiff", width=300, height=130, units = "mm", res=300, pointsize = 0.99)
grid.arrange(nat_hist_p1a, nat_hist_p2_3b, ncol = 2, widths=c(2,3))
#dev.off()

# Survival curves
#tiff(file = "mortality_curves_fit.tiff", width=300, height=190, units = "mm", res=200, pointsize = 0.99)
grid.arrange(mort_curve_cirrhosis, mort_curve_hcc, nrow=2)
#dev.off()

# Correlation plot for thesis
#tiff(file = "correlation_plot_fit.tiff", width=200, height=150, units = "mm", res=200, pointsize = 0.99)
corrplot
#dev.off()



