# Code for plotting fit to the calibration targets

### Load packages and accepted simulations ----
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
#load(file = here("output", "fits", "best_fits_50_of_100000_wed_domain_weights_210819.Rdata"))  # this was from LSR
load(here("calibration", "output", "model_fit_output_119_060120.Rdata")) # out_mat
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
# HBsAg prevalence by time, age and sex ----

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
  labs(title = "HBsAg prevalence in The Gambia",
       y = "HBsAg prevalence (proportion)", x = "Age (years)",
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

ggplot(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence" & time < 1990)) +
  geom_line(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence" & time == 1989),
            aes(x = age, y = model_median, group = sex, colour = sex)) +
  geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
             shape = 4, stroke = 1.5)

# Anti-HBc prevalence by time and age ----
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

tiff(here("output", "fits", "lsr_plots", "antihbc_plot.tiff"), height = 9.4, width = 12, units = 'in', res=300)
ggplot(data = subset(seromarker_out_mat, outcome == "Anti_HBc_prevalence")) +
  geom_line(aes(x = age, y = model_median, linetype = "Model", colour = sex), size = 1) +
  geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
             shape = 4, stroke = 2) +
  geom_line(aes(x = age, y = model_ci_lower, colour = sex), linetype = "dashed", size = 1) +
  geom_line(aes(x = age, y = model_ci_upper, colour = sex), linetype = "dashed", size = 1) +
  #geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "black"), labels = "Observed prevalence") +
  scale_colour_manual(values = c("Mixed" = "gray35", "Male" = "navyblue", "Female" = "steelblue")) +
  facet_wrap(~ time, ncol = 3) +
  geom_text(size = 3.5, data = antihbc_study_labels,
            mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
  labs(title = "Anti-HBc prevalence in The Gambia",
       y = "Anti-HBc prevalence (proportion)", x = "Age (years)",
       colour = "Sex") +
  theme_bw() +
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
dev.off()

# HBeAg prevalence by time and age and overall ----

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
tiff(here("output", "fit_plots", "hbeag_plot.tiff"), height = 8, width = 8, units = 'in', res=300)
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
dev.off()

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




# GLOBOCAN rates ----
globocan_outcome_facet_labels <- c("HCC case incidence", "HCC mortality")
names(globocan_outcome_facet_labels) <- c("hcc_incidence", "hcc_mortality")
# GLOBOCAN PAF-adjusted cancer incidence and mortality in 2018
# Remove 95%CI for data
#globocan_out_mat_long$ci_lower[globocan_out_mat_long$type == "data_value"] <- NA
#globocan_out_mat_long$ci_upper[globocan_out_mat_long$type == "data_value"] <- NA

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


# Risk of chronic carriage ----
tiff(here("output", "fits", "lsr_plots", "p_chronic_plot.tiff"), height = 9, width = 14, units = 'in', res=300)
ggplot(data = p_chronic_out_mat) +
  geom_line(aes(x = age, y = model_median, group = "Model", linetype = "Model"), size = 1, colour = "gray35") +
  geom_line(aes(x = age, y = model_ci_lower), linetype = "dashed", size = 1, colour = "gray35") +
  geom_line(aes(x = age, y = model_ci_upper), linetype = "dashed", size = 1, colour = "gray35") +
  geom_point(aes(x = age, y = data_value, fill = "Observed risk"),
             shape = 4, stroke = 2, colour = "gray35") +
  scale_linetype_manual(name = "", values = c("Model" = "solid"), labels = c("Model" = "Model projection")) +
  labs(title = "Risk of chronic carriage by age at infection",
       y = "Risk of chronic carriage (proportion)", x = "Age at infection (years)",
       fill = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.title = element_text(size = 9),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  ylim(0,1) +
  xlim(0,30)
dev.off()

# Nat hist prev ----
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
tiff(here("output", "fit_plots", "shimakawa_hcc_rate_plot.tiff"), height = 6.4, width = 8, units = 'in', res=300)
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
dev.off()




