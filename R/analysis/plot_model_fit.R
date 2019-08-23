# Code for plotting fit to the calibration targets

# Need to load an out_mat with all the accepted simulations
library(here)
load(file = here("output", "fits", "best_fits_50_of_100000_wed_domain_weights_210819.Rdata"))

best_fits_out_mat <- out_mat_wed_domain_50
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
ci_lower_model_values <- apply(best_fits_model_values_mat,1,quantile, probs = 0.05)
ci_upper_model_values <- apply(best_fits_model_values_mat,1,quantile, probs = 0.95)

# Prepare a mock out_mat to work on
out_mat <- best_fits_out_mat[[1]]$mapped_output
# Remove model_value column to avoid confusion
#out_mat <- lapply(out_mat, function(x) {x$model_value <- NULL ; x})

# Split into different datasets for plotting

# Globocan rates
# Remove data ci_lower and ci_upper to avoid confusion
globocan_out_mat <- select(out_mat$globocan_hcc_incidence, -ci_lower, -ci_upper)
# Turn into long format then remove the mock value column
globocan_out_mat_long <- gather(globocan_out_mat, key = "type", value = "value",
                                c("data_value", "model_value"))
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
ggplot(data = subset(seromarker_out_mat, outcome == "HBsAg_prevalence")) +
  geom_line(aes(x = age, y = model_median, linetype = "Model", colour = sex)) +
  geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
             shape = 4, stroke = 1.5) +
  geom_line(aes(x = age, y = model_ci_lower, colour = sex), linetype = "dashed") +
  geom_line(aes(x = age, y = model_ci_upper, colour = sex), linetype = "dashed") +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "black"), labels = "Observed prevalence") +
  facet_wrap(~ time, ncol = 3) +
  geom_text(size = 3, data = hbsag_study_labels,
            mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
  labs(title = "HBsAg prevalence over time and by age",
       y = "HBsAg prevalence (proportion)", x = "Age (years)",
       colour = "Sex",
       caption = "Keneba Manduar cohort: Whittle studies, Van der Sande 2005 | GHIS: Chotard 1992, Fortuin 1993, Wild 1993 | GLCS: Kirk 2004 | PROLIFICA: Lemoine 2016") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm")) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3))

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

ggplot(data = subset(seromarker_out_mat, outcome == "Anti_HBc_prevalence")) +
  geom_line(aes(x = age, y = model_median, linetype = "Model", colour = sex)) +
  geom_point(aes(x = age, y = data_value, fill = "Data", colour = sex),
             shape = 4, stroke = 1.5) +
  geom_line(aes(x = age, y = model_ci_lower, colour = sex), linetype = "dashed") +
  geom_line(aes(x = age, y = model_ci_upper, colour = sex), linetype = "dashed") +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, colour = sex)) +
  scale_linetype_manual(name = NULL, values = c("Model" = "solid"), labels = "Model projection") +
  scale_fill_manual(name = NULL, values = c("Data" = "black"), labels = "Observed prevalence") +
  facet_wrap(~ time, ncol = 3) +
  geom_text(size = 3, data = antihbc_study_labels,
            mapping = aes(x = Inf, y = Inf, label = label), hjust=1.05, vjust=1.5) +
  labs(title = "Anti-HBc prevalence over time and by age",
       y = "Anti-HBc prevalence (proportion)", x = "Age (years)",
       colour = "Sex",
       caption = "Keneba Manduar cohort: Whittle studies") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.margin=margin(t = 0, unit="cm")) +
  guides(linetype = guide_legend(order = 1),
         fill = guide_legend(order = 2),
         colour = guide_legend(order = 3))

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
ggplot(data = globocan_out_mat_long[globocan_out_mat_long$time == 2018,],
       aes(x = paste(age_min,"-",age_max))) +
  geom_col(aes(y = median*100000, group = type, fill= type),
           position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower*100000, ymax = ci_upper*100000,
                    group = type, colour = type), position = position_dodge(width=0.9), width = 0.3)+
  facet_grid(outcome ~ sex,
             labeller = labeller(outcome =globocan_outcome_facet_labels)) +
  scale_fill_manual("", values = c("model_value" = "gray35", "data_value" = "coral"),
                    labels = c("Median model projection", "GLOBOCAN observed rate")) +
  scale_colour_manual("", values = c("model_value" = "black", "data_value" = "red"),
                      labels = c("Model 5th and 95th percentile", "Observed rate 95% CI")) +
  theme_bw() +
  labs(title = "GLOBOCAN HBV-related HCC incidence and mortality rates in 2018",
       y = "Cases/deaths per 100000 PY", x = "Age (years)",
       subtitle = "GLOBOCAN rates were multiplied by PAF from Ryder 1992 and Kirk 2004 (GLCS)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.margin=margin(t = 0, unit="cm"))

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

