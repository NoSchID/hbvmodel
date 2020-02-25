# Impact of HBsAg screening frequency 27/01/20
require(here)  # for setting working directory
library(tidyr)
library(dplyr)
library(ggplot2)

# Function to extract total incident HBV-related deaths at each timestep for each simulation ----
extract_hbv_deaths <- function(output_file) {

  # Extract absolute incident HBV-related deaths per timestep
  proj_deaths <- cbind(output_file[[1]]$time,
                       (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                          sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                          sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total")))
  proj_deaths <- as.data.frame(proj_deaths)
  colnames(proj_deaths)[1] <- "time"

  return(proj_deaths)

}

extract_hcc_cases <- function(output_file) {

  # Extract absolute incident HBV-related deaths per timestep
  proj_hcc <- cbind(output_file[[1]]$time,
                       (sapply(lapply(output_file,"[[", "incident_hcc"), "[[", "incident_number_total")+
                          sapply(lapply(output_file,"[[", "screened_incident_hcc"), "[[", "incident_number_total")+
                          sapply(lapply(output_file,"[[", "treated_incident_hcc"), "[[", "incident_number_total")))
  proj_hcc <- as.data.frame(proj_hcc)
  colnames(proj_hcc)[1] <- "time"

  return(proj_hcc)

}

# Function to extract outcomes of interest (e.g. infection incidence, HBV-related deaths) from output of
# simulations with multiple parameter sets (requires running code_model_output first)
extract_outcomes <- function(output_file, scenario_label) {

  # HBsAg prevalence

  proj_prev <- cbind(output_file[[1]]$time,
                     (sapply(lapply(output_file,"[[", "infectioncat_total"), "[[", "carriers")/
                        sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total")))
  colnames(proj_prev)[1] <- "time"
  proj_prev_summary <- data.frame(time = output_file[[1]]$time)
  proj_prev_summary$median <- apply(proj_prev[,-1],1,median)
  proj_prev_summary$lower <- apply(proj_prev[,-1],1,quantile, prob = 0.025)
  proj_prev_summary$upper <- apply(proj_prev[,-1],1,quantile, prob = 0.975)
  proj_prev_long <- gather(as.data.frame(proj_prev), key = "iteration", value = "hbsag_prev", -time)

  # Absolute number of new cases of chronic HBV carriage per timestep
  proj_inc <- cbind(output_file[[1]]$time,
                    (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                       sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births")+
                       sapply(lapply(output_file,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections")))
  colnames(proj_inc)[1] <- "time"

  proj_inc_summary <- data.frame(time = output_file[[1]]$time)
  proj_inc_summary$median <- apply(proj_inc[,-1],1,median)
  proj_inc_summary$lower <- apply(proj_inc[,-1],1,quantile, prob = 0.025)
  proj_inc_summary$upper <- apply(proj_inc[,-1],1,quantile, prob = 0.975)
  proj_inc_long <- gather(as.data.frame(proj_inc), key = "iteration", value = "chronic_cases", -time)

  # Incidence rate of chronic HBV carriage per timestep
  proj_inc_rate <- cbind(output_file[[1]]$time,
                         (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                            sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births")+
                            sapply(lapply(output_file,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections"))/
                           sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total"))
  colnames(proj_inc_rate)[1] <- "time"
  proj_inc_rate_summary <- data.frame(time = output_file[[1]]$time)
  proj_inc_rate_summary$median <- apply(proj_inc_rate[,-1],1,median)
  proj_inc_rate_summary$lower <- apply(proj_inc_rate[,-1],1,quantile, prob = 0.025)
  proj_inc_rate_summary$upper <- apply(proj_inc_rate[,-1],1,quantile, prob = 0.975)
  proj_inc_rate_long <- gather(as.data.frame(proj_inc_rate), key = "iteration", value = "chronic_cases_rate", -time)

  # Absolute incident HBV-related deaths per timestep
  proj_deaths <- cbind(output_file[[1]]$time,
                       (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                          sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                          sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total")))
  colnames(proj_deaths)[1] <- "time"
  proj_deaths_summary <- data.frame(time = output_file[[1]]$time)
  proj_deaths_summary$median <- apply(proj_deaths[,-1],1,median)
  proj_deaths_summary$lower <- apply(proj_deaths[,-1],1,quantile, prob = 0.025)
  proj_deaths_summary$upper <- apply(proj_deaths[,-1],1,quantile, prob = 0.975)
  proj_deaths_long <- gather(as.data.frame(proj_deaths), key = "iteration", value = "deaths", -time)

  # Incidence rate of  HBV-related deaths per timestep
  proj_deaths_rate <- cbind(output_file[[1]]$time,
                            (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                               sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                               sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total"))/
                              sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total"))
  colnames(proj_deaths_rate)[1] <- "time"
  proj_deaths_rate_summary <- data.frame(time = output_file[[1]]$time)
  proj_deaths_rate_summary$median <- apply(proj_deaths_rate[,-1],1,median)
  proj_deaths_rate_summary$lower <- apply(proj_deaths_rate[,-1],1,quantile, prob = 0.025)
  proj_deaths_rate_summary$upper <- apply(proj_deaths_rate[,-1],1,quantile, prob = 0.975)
  proj_deaths_rate_long <- gather(as.data.frame(proj_deaths_rate), key = "iteration", value = "deaths_rate", -time)

  # Age-standardised chronic infection incidence per 100000 per timestep!!

  # Age-standardised HBV-related deaths incidence per 100000 per timestep
  #proj_deaths_standardised <- sapply(output_file, calculate_age_standardised_hbvdeaths_rate) # per person per timestep
  #proj_deaths_standardised_summary <- data.frame(time = output_file[[1]]$time)
  #proj_deaths_standardised_summary$median <- apply(proj_deaths_standardised,1,median)
  #proj_deaths_standardised_summary$lower <- apply(proj_deaths_standardised,1,quantile, prob = 0.025)
  #proj_deaths_standardised_summary$upper <- apply(proj_deaths_standardised,1,quantile, prob = 0.975)
  #  proj_deaths_standardised_long <- gather(as.data.frame(proj_deaths_standardised), key = "iteration", value = "hbv_death_rate", -time)

  # Label dataframes with scenario:
  proj_prev_summary$scenario <- scenario_label
  proj_inc_summary$scenario <- scenario_label
  proj_inc_rate_summary$scenario <- scenario_label
  proj_deaths_summary$scenario <- scenario_label
  proj_deaths_rate_summary$scenario <- scenario_label
  #proj_deaths_standardised_summary$scenario <- scenario_label

  return(list(hbsag_prev_summary = proj_prev_summary,
              chronic_incidence_summary = proj_inc_summary,
              chronic_incidence_rate_summary = proj_inc_rate_summary,
              hbv_deaths_summary = proj_deaths_summary,
              #proj_deaths_standardised_summary = proj_deaths_standardised_summary,
              hbv_deaths_rate_summary = proj_deaths_rate_summary))
}

load(here("output", "HBsAg screening intervals", "output_scenario_status_quo_220120.RData"))
deaths_sq <- extract_hbv_deaths(out_sq[[1]])
hcc_sq <- extract_hcc_cases(out_sq[[1]])
out_plots_sq <- extract_outcomes(output_file = out_sq[[1]], scenario_label = "status_quo")
rm(out_sq)
gc()

load(here("output", "HBsAg screening intervals", "screening_freq_once_270120.RData"))
deaths_screen_once <- extract_hbv_deaths(screen_once)
hcc_screen_once <- extract_hcc_cases(screen_once)
out_plots_screen_once <- extract_outcomes(output_file = screen_once, scenario_label = "screen_once")
rm(screen_once)
gc()

load(here("output", "HBsAg screening intervals", "screening_freq_0point5_230120.RData"))
deaths_screen_0point5 <- extract_hbv_deaths(screen_0point5)
hcc_screen_0point5 <- extract_hcc_cases(screen_0point5)
rm(screen_0point5)
gc()

load(here("output", "HBsAg screening intervals", "screening_freq_1_250120.RData"))
deaths_screen_1 <- extract_hbv_deaths(screen_1)
hcc_screen_1 <- extract_hcc_cases(screen_1)
rm(screen_1)
gc()

load(here("output", "HBsAg screening intervals", "screening_freq_5_250120.RData"))
deaths_screen_5 <- extract_hbv_deaths(screen_5)
hcc_screen_5 <- extract_hcc_cases(screen_5)
out_plots_screen_5 <- extract_outcomes(output_file = screen_5, scenario_label = "screen_5")
rm(screen_5)
gc()

load(here("output", "HBsAg screening intervals", "screening_freq_10_250120.RData"))
deaths_screen_10 <- extract_hbv_deaths(screen_10)
hcc_screen_10 <- extract_hcc_cases(screen_10)
rm(screen_10)
gc()

load(here("output", "HBsAg screening intervals", "screening_freq_20_260120.RData"))
deaths_screen_20 <- extract_hbv_deaths(screen_20)
hcc_screen_20 <- extract_hcc_cases(screen_20)
rm(screen_20)
gc()



# Calculate cumulative number of HBV-related deaths averted compared to status quo scenario ----
# Since 2020.0 (or any other time)

# Function to calculate total deaths averted by scenario compared to counterfactual
# from_year is the timestep-0.5 from which to count deaths (e.g. 2020.5 = deaths that have occured since 2020)
# by_year_vector is the vector of timesteps up until which to count cumulative deaths (i.e. 2031 includes all deaths in 2030)
# Returns the median, 2.5h and 97.5h percentile of all parameter sets
calculate_cumulative_hbv_deaths_averted <- function(scenario_deaths, counterfactual_deaths,
                                                    from_year, by_year_vector, scenario_label) {

  cum_deaths_counterfactual <- list()
  cum_deaths_scenario <- list()
  cum_deaths_averted <- list()
  cum_deaths_averted_percent <- list()
  cum_deaths_averted_summary <- data.frame(V1 = rep(0, 3))
  cum_deaths_averted_percent_summary <- data.frame(V1 = rep(0, 3))

  for (i in 1:length(by_year_vector)) {
  cum_deaths_counterfactual[[i]] <- apply(counterfactual_deaths[which(counterfactual_deaths$time==from_year):
                                                 which(counterfactual_deaths$time==by_year_vector[i]),-1],2,sum)
  cum_deaths_scenario[[i]] <- apply(scenario_deaths[which(scenario_deaths$time==from_year):
                                                            which(scenario_deaths$time==by_year_vector[i]),-1],2,sum)
  cum_deaths_averted[[i]] <- cum_deaths_counterfactual[[i]]-cum_deaths_scenario[[i]]
  cum_deaths_averted_percent[[i]] <- (cum_deaths_counterfactual[[i]]-cum_deaths_scenario[[i]])/cum_deaths_counterfactual[[i]]

  # Summarise into median, 2.5th and 97.5th percentile
  cum_deaths_averted_summary[,i] <- quantile(cum_deaths_averted[[i]], prob = c(0.025,0.5,0.975))
  # Summarise percent reduction into median, 2.5th and 97.5th percentile
  cum_deaths_averted_percent_summary[,i] <- quantile(cum_deaths_averted_percent[[i]], prob = c(0.025,0.5,0.975))
  }

  cum_deaths_averted_summary <- as.data.frame(t(cum_deaths_averted_summary))
  colnames(cum_deaths_averted_summary) <- c("ui_lower", "median", "ui_upper")
  cum_deaths_averted_summary$by_year <- by_year_vector
  cum_deaths_averted_summary$scenario <- scenario_label
  rownames(cum_deaths_averted_summary) <- NULL

  cum_deaths_averted_prop_summary <- as.data.frame(t(cum_deaths_averted_percent_summary))
  colnames(cum_deaths_averted_prop_summary) <- c("ui_lower", "median", "ui_upper")
  cum_deaths_averted_prop_summary$by_year <- by_year_vector
  cum_deaths_averted_prop_summary$scenario <- scenario_label
  rownames(cum_deaths_averted_prop_summary) <- NULL

  cum_deaths_averted <- as.data.frame(matrix(unlist(cum_deaths_averted),
                                             nrow=length(by_year_vector), byrow = TRUE))
  cum_deaths_averted$by_year <- by_year_vector
  cum_deaths_averted$scenario <- scenario_label


  return(list("cum_deaths_averted_summary" = cum_deaths_averted_summary,
              "cum_deaths_averted" = cum_deaths_averted,
              "cum_deaths_averted_prop_summary" = cum_deaths_averted_prop_summary))
}


# Percent reduction from status quo
cumulative_deaths_averted_prop <- rbind(calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_5,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_5")[[3]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_once,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_once")[[3]])

# Percent reduction from status quo
cumulative_deaths_averted <- rbind(calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_5,
                                                                                counterfactual_deaths = deaths_sq,
                                                                                from_year = 2020.5,
                                                                                by_year_vector = c(2030, 2050, 2099),
                                                                                scenario_label = "screen_5")[[1]],
                                        calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_once,
                                                                                counterfactual_deaths = deaths_sq,
                                                                                from_year = 2020.5,
                                                                                by_year_vector = c(2030, 2050, 2099),
                                                                                scenario_label = "screen_once")[[1]])


# Counterfactual: status quo
cumulative_deaths_averted <- rbind(calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_0point5,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_0point5")[[1]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_1,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_1")[[1]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_5,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_5")[[1]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_10,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_10")[[1]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_20,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_20")[[1]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_once,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_once")[[1]])

level_order <- c("screen_0point5", "screen_1", "screen_5",
                 "screen_10","screen_20","screen_once")

cumulative_deaths_averted$by_year_label <- paste("By", cumulative_deaths_averted$by_year)
cumulative_deaths_averted$by_year_label[cumulative_deaths_averted$by_year_label == "By 2099"] <- "By 2100"
cumulative_deaths_averted$value[cumulative_deaths_averted$by_year == 2030 &
                                      cumulative_deaths_averted$scenario %in% c("screen_10", "screen_20")] <- NA
cumulative_deaths_averted$scenario[cumulative_deaths_averted$by_year == 2030 &
                                  cumulative_deaths_averted$scenario %in% c("screen_10", "screen_20")] <- NA
# Remove continuous scenario
cumulative_deaths_averted$scenario[cumulative_deaths_averted$scenario == "screen_0point5"] <- NA
cumulative_deaths_averted$value[cumulative_deaths_averted$scenario == "screen_0point5"] <- NA

ggplot(data = cumulative_deaths_averted[is.na(cumulative_deaths_averted$scenario)==FALSE &
                                          cumulative_deaths_averted$by_year != 2050,]) +
  geom_point(aes(x = factor(scenario, level = level_order), y = median)) +
  geom_errorbar(aes(x = factor(scenario, level = level_order), ymin=ui_lower, ymax=ui_upper)) +
  scale_x_discrete(labels = c("screen_0point5" = "Continuous", "screen_1" = "1 year", "screen_5" = "5 years",
                              "screen_10" = "10 years", "screen_20" = "20 years", "screen_once" = "One-off")) +
  facet_wrap(~by_year_label, scales = "free") +
  xlab("Screening frequency") +
  ylab("Cumulative number of HBV-related deaths averted") +
  labs(title = "Effect of screening frequency on HBV-related deaths averted compared to status quo scenario") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = cumulative_deaths_averted[is.na(cumulative_deaths_averted$scenario)==FALSE &
                                          cumulative_deaths_averted$by_year == 2030,]) +
  geom_point(aes(x = factor(scenario, level = level_order), y = median), size = 3) +
  geom_errorbar(aes(x = factor(scenario, level = level_order), ymin=ui_lower, ymax=ui_upper), width = 0.3) +
  scale_x_discrete(labels = c("screen_0point5" = "Continuous", "screen_1" = "1 year", "screen_5" = "5 years",
                              "screen_10" = "10 years", "screen_20" = "20 years", "screen_once" = "One-off")) +
  xlab("Screening frequency") +
  ylab("Cumulative number of HBV-related deaths averted") +
  labs(title = "By 2030") +
  theme_classic() +
  ylim(0,6500) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5))


ggplot(data = cumulative_deaths_averted[is.na(cumulative_deaths_averted$scenario)==FALSE &
                                          cumulative_deaths_averted$by_year == 2099,]) +
  geom_point(aes(x = factor(scenario, level = level_order), y = median), size = 3) +
  geom_errorbar(aes(x = factor(scenario, level = level_order), ymin=ui_lower, ymax=ui_upper), width = 0.3) +
  scale_x_discrete(labels = c("screen_0point5" = "Continuous", "screen_1" = "1 year", "screen_5" = "5 years",
                              "screen_10" = "10 years", "screen_20" = "20 years", "screen_once" = "One-off")) +
  xlab("Screening frequency") +
  ylab("Cumulative number of HBV-related deaths averted") +
  labs(title = "By 2100") +
  theme_classic() +
  ylim(0,65000) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5))

# For by 2031, a screening frequency of 20 years represents a one off screen
# For by 2030, a screening frequency of 10 year represents one-off screen
# For by 2050, a screening frequency of 10 year represents 2 screens

# By 2030: freq 1 = 10 screens, freq 5 = 2 screens, freq 10 = 1 screen , freq 20 = 1 screen
# By 2050: freq 1 = 30 screens, freq 5 = 6 screens, freq 10 = 3 screens, freq 20 = 2 screens
# By 2099: freq 1 = 79 screens , freq 5 = 16 screens, freq 10 = 8 screens, freq 20 = 4 screens

seq(2020,2100, by = 20)

ggplot(data = cumulative_deaths_averted) +
  geom_point(aes(x = factor(scenario, level = level_order), y = median)) +
  geom_errorbar(aes(x = factor(scenario, level = level_order), ymin=ui_lower, ymax=ui_upper)) +
  facet_wrap(~by_year,scales = "free")


cumulative_deaths_averted_bp <- rbind(calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_0point5,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_0point5")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_1,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_1")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_5,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_5")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_10,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_10")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_20,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_20")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_once,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_once")[[2]])


cumulative_deaths_averted_bp2 <- gather(cumulative_deaths_averted_bp, key = "sim", value = "value", -scenario, - by_year)

cumulative_deaths_averted_bp2$value[cumulative_deaths_averted_bp2$by_year == 2030 &
                               cumulative_deaths_averted_bp2$scenario %in% c("screen_10", "screen_20")] <- NA

cumulative_deaths_averted_bp2$by_year_label <- paste("By", cumulative_deaths_averted_bp2$by_year)
cumulative_deaths_averted_bp2$by_year_label[cumulative_deaths_averted_bp2$by_year_label == "By 2099"] <- "By 2100"

ggplot(data = cumulative_deaths_averted_bp2) +
  geom_boxplot(aes(x = factor(scenario, level = level_order), y = value)) +
  scale_x_discrete(labels = c("screen_0point5" = "Continuous", "screen_1" = "1 year", "screen_5" = "5 years",
                              "screen_10" = "10 years", "screen_20" = "20 years", "screen_once" = "One-off")) +
  facet_wrap(~by_year_label,scales = "free") +
  xlab("Screening frequency") +
  ylab("Cumulative number of HBV-related deaths averted") +
  labs(title = "Effect of screening frequency on HBV-related deaths averted compared to status quo scenario") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = cumulative_deaths_averted_bp[,c(1,120,121)]) +
  geom_point(aes(x = factor(scenario, level = level_order), y = V1)) +
  facet_wrap(~by_year, scales = "free")


###
x_6m <- out
x_6m_deaths <- x_6m$hbv_deaths$incident_number_total+x_6m$screened_hbv_deaths$incident_number_total+
  x_6m$treated_hbv_deaths$incident_number_total

x_sq <- out
x_sq_deaths <- x_sq$hbv_deaths$incident_number_total+x_sq$screened_hbv_deaths$incident_number_total+
  x_sq$treated_hbv_deaths$incident_number_total

x_1 <- out
x_1_deaths <- x_1$hbv_deaths$incident_number_total+x_1$screened_hbv_deaths$incident_number_total+
  x_1$treated_hbv_deaths$incident_number_total

sum(x_sq_deaths[which(x_sq$time==2020.5):which(x_sq$time==2030)])-
  sum(x_6m_deaths[which(x_sq$time==2020.5):which(x_sq$time==2030)])

sum(x_sq_deaths[which(x_sq$time==2020.5):which(x_sq$time==2099)])-
  sum(x_6m_deaths[which(x_sq$time==2020.5):which(x_sq$time==2099)])

sum(x_sq_deaths[which(x_sq$time==2020.5):which(x_sq$time==2030)])-
  sum(x_1_deaths[which(x_sq$time==2020.5):which(x_sq$time==2030)])

sum(x_sq_deaths[which(x_sq$time==2020.5):which(x_sq$time==2099)])-
  sum(x_1_deaths[which(x_sq$time==2020.5):which(x_sq$time==2099)])


cum_deaths_counterfactual[[i]] <- apply(counterfactual_deaths[which(counterfactual_deaths$time==from_year):
                                                                which(counterfactual_deaths$time==by_year_vector[i]),-1],2,sum)
cum_deaths_scenario[[i]] <- apply(scenario_deaths[which(scenario_deaths$time==from_year):
                                                    which(scenario_deaths$time==by_year_vector[i]),-1],2,sum)
cum_deaths_averted[[i]] <- cum_deaths_counterfactual[[i]]-cum_deaths_scenario[[i]]


# HCC cases

cumulative_hcc_cases_averted_bp <- rbind(calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_0point5,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_0point5")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_1,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_1")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_5,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_5")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_10,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_10")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_20,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_20")[[2]],
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_once,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_once")[[2]])


cumulative_hcc_cases_averted_bp <- gather(cumulative_hcc_cases_averted_bp, key = "sim", value = "value", -scenario, - by_year)

ggplot(data = cumulative_hcc_cases_averted_bp) +
  geom_boxplot(aes(x = factor(scenario, level = level_order), y = value)) +
  facet_wrap(~by_year,scales = "free")


# Counterfactual: one-off screening in 2020
cumulative_deaths_averted2 <- rbind(calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_0point5,
                                                                            counterfactual_deaths = deaths_screen_once,
                                                                            from_year = 2020.5,
                                                                            by_year_vector = c(2030, 2050, 2099),
                                                                            scenario_label = "screen_0point5"),
                                    calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_1,
                                                                            counterfactual_deaths = deaths_screen_once,
                                                                            from_year = 2020.5,
                                                                            by_year_vector = c(2030, 2050, 2099),
                                                                            scenario_label = "screen_1"),
                                    calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_5,
                                                                            counterfactual_deaths = deaths_screen_once,
                                                                            from_year = 2020.5,
                                                                            by_year_vector = c(2030, 2050, 2099),
                                                                            scenario_label = "screen_5"),
                                    calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_10,
                                                                            counterfactual_deaths = deaths_screen_once,
                                                                            from_year = 2020.5,
                                                                            by_year_vector = c(2030, 2050, 2099),
                                                                            scenario_label = "screen_10"),
                                    calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_20,
                                                                            counterfactual_deaths = deaths_screen_once,
                                                                            from_year = 2020.5,
                                                                            by_year_vector = c(2030, 2050, 2099),
                                                                            scenario_label = "screen_20"))




# Plot status quo vs screen-and-treat scenarios : HBV related deaths ----
# Combine projection summaries from different scenarios
proj_prev_total <- rbind(out_plots_sq$hbsag_prev_summary, out_plots_screen_once$hbsag_prev_summary,
                         out_plots_screen_5$hbsag_prev_summary)
proj_inc_total <- rbind(out_plots_sq$chronic_incidence_summary, out_plots_screen_once$chronic_incidence_summary,
                        out_plots_screen_5$chronic_incidence_summary)
proj_deaths_total <- rbind(out_plots_sq$hbv_deaths_summary, out_plots_screen_once$hbv_deaths_summary,
                           out_plots_screen_5$hbv_deaths_summary)

# HBV-related deaths
ggplot(proj_deaths_total) +
  geom_ribbon(data = proj_deaths_total[proj_deaths_total$scenario=="status_quo",],
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, group = scenario, colour = scenario,fill = scenario),
              colour = NA, alpha = 0.1)+
  geom_ribbon(data = proj_deaths_total[proj_deaths_total$scenario!="status_quo",],
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, group = scenario, colour = scenario, fill = scenario),
              linetype = "dashed", alpha = 0.15)+
  geom_line(aes(x=time, y = median/0.5, group = scenario, colour = scenario), size =1)+
  labs(title = "HBV-related deaths",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "Status quo: historical infant vaccine coverage since 1990 and maintaining 93% coverage after 2018\n
       One-off screen+treat: infant vaccine and screening once in 2020\n
       Repeat screen+treat: Infant vaccine and screening every 5 years starting in 2020") +
  scale_color_manual(limits = c("status_quo", "screen_once", "screen_5"),
                     labels = c("status_quo"="Status quo",
                                "screen_once"="One-off screen+treat",
                                "screen_5"="Repeat screen+treat"),
  values = c("screen_once"="steelblue", "status_quo"="orange", "screen_5"="deeppink")) +
  scale_fill_manual(limits = c("status_quo", "screen_once", "screen_5"),
                    labels = c("status_quo"="Status quo",
                               "screen_once"="One-off screen+treat",
                               "screen_5"="Repeat screen+treat"),
                    values = c("screen_once"="steelblue", "status_quo"="orange", "screen_5"="deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(2015,2100)) +
  ylab("Annual number of HBV-related deaths")+
  xlab("Year")+
  ylim(0,2000) +
  theme_classic()

# HBV-related deaths
ggplot(proj_deaths_total[proj_deaths_total$scenario!="screen_5",]) +
  geom_ribbon(data = proj_deaths_total[proj_deaths_total$scenario=="status_quo",],
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, group = scenario, colour = scenario,fill = scenario),
              colour = NA, alpha = 0.1)+
  geom_ribbon(data = proj_deaths_total[proj_deaths_total$scenario=="screen_once",],
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, group = scenario, colour = scenario, fill = scenario),
              linetype = "dashed", alpha = 0.15)+
  geom_line(aes(x=time, y = median/0.5, group = scenario, colour = scenario), size =1)+
  labs(title = "HBV-related deaths",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "Status quo: historical infant vaccine coverage since 1990 and maintaining 93% coverage after 2018\n
       One-off screen+treat: infant vaccine and screening once in 2020\n
       Repeat screen+treat: Infant vaccine and screening every 5 years starting in 2020") +
  scale_color_manual(limits = c("status_quo", "screen_once"),
                     labels = c("status_quo"="Status quo",
                                "screen_once"="One-off screen+treat"),
                     values = c("screen_once"="steelblue", "status_quo"="orange")) +
  scale_fill_manual(limits = c("status_quo", "screen_once"),
                    labels = c("status_quo"="Status quo",
                               "screen_once"="One-off screen+treat"),
                    values = c("screen_once"="steelblue", "status_quo"="orange")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(2015,2030)) +
  ylab("Annual number of HBV-related deaths")+
  xlab("Year")+
  ylim(0,2000) +
  theme_classic()

# black, steelblue, deeppink, orange


# Plot status quo vs screen-and-treat scenarios : HCC incidence ----

# Summarise HCC values
hcc_sq_summary <- data.frame(time = hcc_sq$time, scenario = "status_quo")
hcc_sq_summary$median <- apply(hcc_sq[,-1],1,median)
hcc_sq_summary$ui_lower <- apply(hcc_sq[,-1],1,quantile, prob = 0.025)
hcc_sq_summary$ui_upper <- apply(hcc_sq[,-1],1,quantile, prob = 0.975)

hcc_screen_once_summary <- data.frame(time = hcc_screen_once$time, scenario = "screen_once")
hcc_screen_once_summary$median <- apply(hcc_screen_once[,-1],1,median)
hcc_screen_once_summary$ui_lower <- apply(hcc_screen_once[,-1],1,quantile, prob = 0.025)
hcc_screen_once_summary$ui_upper <- apply(hcc_screen_once[,-1],1,quantile, prob = 0.975)

hcc_screen_5_summary <- data.frame(time = hcc_screen_5$time, scenario = "screen_5")
hcc_screen_5_summary$median <- apply(hcc_screen_5[,-1],1,median)
hcc_screen_5_summary$ui_lower <- apply(hcc_screen_5[,-1],1,quantile, prob = 0.025)
hcc_screen_5_summary$ui_upper <- apply(hcc_screen_5[,-1],1,quantile, prob = 0.975)

# Combine into 1 dataframe
proj_hcc_total <- rbind(hcc_sq_summary, hcc_screen_once_summary, hcc_screen_5_summary)

ggplot(proj_hcc_total) +
  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
  geom_ribbon(data = proj_hcc_total[proj_hcc_total$scenario == "status_quo",],
              aes(x=time, ymin=ui_lower, ymax=ui_upper, group = scenario, fill = scenario),
              colour = NA, alpha = 0.1)+
  geom_ribbon(data = proj_hcc_total[proj_hcc_total$scenario != "status_quo",],
              aes(x=time, ymin=ui_lower, ymax=ui_upper, group = scenario, colour = scenario), fill = NA, linetype = "dashed", alpha = 0.3)+
  # geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "HCC incidence",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical infant vaccine coverage since 1990 and maintains 93% coverage after 2018\n
       **Birth dose vaccine introduced in 2020 with 80% coverage\n
       ***Birth dose as in ** and one-off screening+continuous treatment of 80% of adults in 2020") +
  #  scale_color_manual(labels = c("vacc"="Status quo infant vaccine*",
  #                                "bdvacc"="Infant+birth dose vaccine**",
  #                                "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
  #                     values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  #  scale_fill_manual(labels = c("vacc"="Status quo infant vaccine*",
  #                               "bdvacc"="Infant+birth dose vaccine**",
  #                               "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
  #                    values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Deaths per 6 months")+
  xlab("Year")+
  ylim(0,400) +
  theme_classic()

ggplot(proj_hcc_total) +
  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
  geom_ribbon(data = proj_hcc_total,
              aes(x=time, ymin=ui_lower, ymax=ui_upper, group = scenario, fill = scenario),
              alpha = 0.1)+
  # geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "HCC incidence",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical infant vaccine coverage since 1990 and maintains 93% coverage after 2018\n
       **Birth dose vaccine introduced in 2020 with 80% coverage\n
       ***Birth dose as in ** and one-off screening+continuous treatment of 80% of adults in 2020") +
  #  scale_color_manual(labels = c("vacc"="Status quo infant vaccine*",
  #                                "bdvacc"="Infant+birth dose vaccine**",
  #                                "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
  #                     values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  #  scale_fill_manual(labels = c("vacc"="Status quo infant vaccine*",
  #                               "bdvacc"="Infant+birth dose vaccine**",
  #                               "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
  #                    values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(2015,2100)) +
  ylab("Deaths per 6 months")+
  xlab("Year")+
  ylim(0,400) +
  theme_classic()
