# Impact of HBsAg screening frequency 27/01/20

require(here)  # for setting working directory

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
  cum_deaths_averted_summary <- data.frame(V1 = rep(0, 3))

  for (i in 1:length(by_year_vector)) {
  cum_deaths_counterfactual[[i]] <- apply(counterfactual_deaths[which(counterfactual_deaths$time==from_year):
                                                 which(counterfactual_deaths$time==by_year_vector[i]),-1],2,sum)
  cum_deaths_scenario[[i]] <- apply(scenario_deaths[which(scenario_deaths$time==from_year):
                                                            which(scenario_deaths$time==by_year_vector[i]),-1],2,sum)
  cum_deaths_averted[[i]] <- cum_deaths_counterfactual[[i]]-cum_deaths_scenario[[i]]

  # Summarise into median, 2.5th and 97.5th percentile
  cum_deaths_averted_summary[,i] <- quantile(cum_deaths_averted[[i]], prob = c(0.025,0.5,0.975))
  }

  cum_deaths_averted_summary <- as.data.frame(t(cum_deaths_averted_summary))
  colnames(cum_deaths_averted_summary) <- c("ui_lower", "median", "ui_upper")
  cum_deaths_averted_summary$by_year <- by_year_vector
  cum_deaths_averted_summary$scenario <- scenario_label
  rownames(cum_deaths_averted_summary) <- NULL

  return(cum_deaths_averted_summary)
}

# Counterfactual: status quo
cumulative_deaths_averted <- rbind(calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_0point5,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_0point5"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_1,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_1"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_5,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_5"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_10,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_10"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_20,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_20"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = deaths_screen_once,
                                                                           counterfactual_deaths = deaths_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_once"))

reps <- 3
cumulative_deaths_averted$screening_frequency <- c(rep(0.5,reps), rep(1,reps), rep(5,reps), rep(10,reps), rep(20,reps), rep(0,reps))

# For by 2031, a screening frequency of 20 years represents a one off screen
# For by 2030, a screening frequency of 10 year represents one-off screen
# For by 2050, a screening frequency of 10 year represents 2 screens

# By 2030: freq 1 = 10 screens, freq 5 = 2 screens, freq 10 = 1 screen , freq 20 = 1 screen
# By 2050: freq 1 = 30 screens, freq 5 = 6 screens, freq 10 = 3 screens, freq 20 = 2 screens
# By 2099: freq 1 = 79 screens , freq 5 = 16 screens, freq 10 = 8 screens, freq 20 = 4 screens

seq(2020,2100, by = 20)

ggplot(data = cumulative_deaths_averted) +
  geom_point(aes(x = screening_frequency, y = median)) +
  geom_errorbar(aes(x = screening_frequency, ymin=ui_lower, ymax=ui_upper)) +
  facet_wrap(~by_year, scales = "free")

cumulative_hcc_cases_averted <- rbind(calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_0point5,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_0point5"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_1,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_1"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_5,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_5"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_10,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_10"),
                                   calculate_cumulative_hbv_deaths_averted(scenario_deaths = hcc_screen_20,
                                                                           counterfactual_deaths = hcc_sq,
                                                                           from_year = 2020.5,
                                                                           by_year_vector = c(2030, 2050, 2099),
                                                                           scenario_label = "screen_20"))

reps <- 3
cumulative_hcc_cases_averted$screening_frequency <- c(rep(0.5,reps), rep(1,reps), rep(5,reps), rep(10,reps), rep(20,reps))

ggplot(data = cumulative_hcc_cases_averted) +
  geom_point(aes(x = screening_frequency, y = median)) +
  geom_errorbar(aes(x =screening_frequency, ymin=ui_lower, ymax=ui_upper)) +
  facet_wrap(~by_year, scales = "free")

ggplot(data = cumulative_deaths_averted[cumulative_deaths_averted$by_year == 2099,]) +
  geom_point(aes(x = screening_frequency, y = median)) +
  geom_errorbar(aes(x =screening_frequency, ymin=ui_lower, ymax=ui_upper)) +
  ylim(0,100000)
ggplot(data = cumulative_hcc_cases_averted[cumulative_hcc_cases_averted$by_year == 2099,]) +
  geom_point(aes(x = screening_frequency, y = median)) +
  geom_errorbar(aes(x =screening_frequency, ymin=ui_lower, ymax=ui_upper)) +
  ylim(0,35000)

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




# Plot status quo vs screen-and-treat scenarios ----
# Combine projection summaries from different scenarios
proj_prev_total <- rbind(out_plots_sq$hbsag_prev_summary, out_plots_screen_once$hbsag_prev_summary,
                         out_plots_screen_5$hbsag_prev_summary)
proj_inc_total <- rbind(out_plots_sq$chronic_incidence_summary, out_plots_screen_once$chronic_incidence_summary,
                        out_plots_screen_5$chronic_incidence_summary)
proj_deaths_total <- rbind(out_plots_sq$hbv_deaths_summary, out_plots_screen_once$hbv_deaths_summary,
                           out_plots_screen_5$hbv_deaths_summary)

# HBV-related deaths
ggplot(proj_deaths_total) +
  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, group = scenario, colour = scenario), fill = NA, linetype = "dashed", alpha = 0.3)+
  # geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "HBV-related deaths",
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
  ylim(0,1000) +
  theme_classic()
