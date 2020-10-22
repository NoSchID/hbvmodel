# Incremental analysis of different monitoring strategies
# Varying frequency and age group

# What is the optimal monitoring frequency for different age groups?
# Exploring this from a cohort perspective, assuming 15-60 year olds have been screened (A1)

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
library(BCEA)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Function to plot boxplot whiskers as 95% percentile
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
# Function to find dominated and extended dominated strategies on cost-effectiveness frontier
find_dominated_strategies <- function(df, exposure, outcome) {
  plot <- ceef.plot(bcea(e=cbind(matrix(rep(0, 183)),
                                 as.matrix(select(df, scenario, sim, outcome) %>%
                                             filter(sim != "X") %>%
                                             spread(key="scenario", value = outcome)%>%
                                             select(-sim))),
                         c=cbind(matrix(rep(0, 183)),
                                 as.matrix(select(df, scenario, sim, exposure) %>%
                                             filter(sim != "X") %>%
                                             spread(key="scenario", value = exposure)%>%
                                             select(-sim))),
                         ref=1,
                         interventions=c("No monitoring", colnames(select(df, scenario, sim, outcome) %>%
                                                                     filter(sim != "X") %>%
                                                                     spread(key="scenario", value = outcome)%>%
                                                                     select(-sim))),
                         Kmax=50000000,
                         plot=FALSE),
                    graph="base", relative = FALSE)
  return(plot)
}
# Function to create dataframe for incremental monitoring strategies plots
create_incremental_plot_df <- function(interactions_df, py_on_treatment_df,
                                       deaths_averted_df, ly_saved_df,
                                       monitoring_assessment_cost = 15.77,
                                       treatment_py_cost = 84.88,
                                       scenario_labels_obj = NULL,
                                       ref_label = "No monitoring") {


  # interactions_df has columns: scenario, sim, monitoring_assessments, treatment_initiations
  # py_on_treatment_df has columns:scenario, sim, py_on_treatment
  # deaths_averted_df and ly_saved_df has columns: scenario, sim, value (sim already in correct format)
  # Ensure all inputs have the same scenarios!
  interactions <- interactions_df
  interactions$total_interactions <- interactions$monitoring_assessments+
    interactions$treatment_initiations
  interactions$sim <- gsub("[^0-9]", "", interactions$sim)

  # Combine interactions and assign costs based on Shevanthi's Gambia paper:
  # $15.77 per monitoring assessment and $84.88 per person-year of treatment
  interactions <- left_join(interactions, py_on_treatment_df,
                            by=c("scenario", "sim")) %>%
    mutate(monitoring_cost = monitoring_assessments*monitoring_assessment_cost,
           treatment_cost = py_on_treatment*treatment_py_cost,
           total_cost = monitoring_cost+treatment_cost)

  # Combine with deaths averted
  df <- deaths_averted_df

  df <- df %>%
    left_join(interactions, by = c("scenario", "sim")) %>%
    drop_na()
  colnames(df)[colnames(df)=="value"] <- "deaths_averted"

  # Combine with LY saved
  df_ly <- ly_saved_df

  df <- df_ly %>%
    left_join(df, by = c("scenario", "sim")) %>%
    drop_na()
  colnames(df)[colnames(df)=="value"] <- "ly_saved"

  # Scenario labels
  if(length(scenario_labels_obj) > 0L) {
    levels(df$scenario) <- scenario_labels_obj
  }

  # Separate columns for frequency and age group
  df$scenario2 <- df$scenario
  df <- extract(df, col="scenario2", into = c("frequency","age_group"), regex="^(\\S+)\\s+(.*)")

  # Add reference point at 0
  ref_df <- data.frame(matrix(c(ref_label, "X", rep(0,(ncol(df)-4))), nrow = 1))
  ref_df <- ref_df[rep(seq_len(nrow(ref_df)), each = nrow(crossing(df$frequency, df$age_group))),]
  ref_df <- cbind(ref_df, data.frame(crossing(df$frequency, df$age_group)))
  colnames(ref_df) <- colnames(df)

  df$scenario <- as.character(df$scenario)
  df <- rbind(df, ref_df)
  df[,-c(1:2,ncol(df)-1,ncol(df))] <- apply(df[,-c(1:2,ncol(df)-1,ncol(df))], 2, as.numeric)
  df$scenario <- factor(df$scenario)

  return(df)


}

# Load files ----

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

# Status quo
out1 <- readRDS(paste0(out_path, "a1_out1_status_quo_cohort_240920.rds"))
out1 <- out1[[1]]
out2 <- readRDS(paste0(out_path, "out2_status_quo_180820.rds"))
out2 <- out2[[1]]

# No monitoring
out3 <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_201020.rds"))
out3 <- out3[[1]]

# Monitoring all age groups
out4 <- readRDS(paste0(out_path, "a1_out4_screen_2020_monit_10_290920.rds"))
out4 <- out4[[1]]
out5 <- readRDS(paste0(out_path, "a1_out5_screen_2020_monit_5_201020.rds"))
out5 <- out5[[1]]
out5a <- readRDS(paste0(out_path, "a1_out5a_screen_2020_monit_6_290920.rds"))
out5a <- out5a[[1]]
out5b <- readRDS(paste0(out_path, "a1_out5b_screen_2020_monit_7_290920.rds"))
out5b <- out5b[[1]]
out5c <- readRDS(paste0(out_path, "a1_out5c_screen_2020_monit_8_290920.rds"))
out5c <- out5c[[1]]
out5d <- readRDS(paste0(out_path, "a1_out5d_screen_2020_monit_9_290920.rds"))
out5d <- out5d[[1]]
out6 <- readRDS(paste0(out_path, "a1_out6_screen_2020_monit_1_201020.rds"))
out6 <- out6[[1]]
out6a <- readRDS(paste0(out_path, "a1_out6a_screen_2020_monit_2_280920.rds"))
out6a <- out6a[[1]]
out6b <- readRDS(paste0(out_path, "a1_out6b_screen_2020_monit_3_280920.rds"))
out6b <- out6b[[1]]
out6c <- readRDS(paste0(out_path, "a1_out6c_screen_2020_monit_4_280920.rds"))
out6c <- out6c[[1]]

# Monitoring different age groups until they age out (yearly)
monit_out1 <- readRDS(paste0(out_path, "a1_monit_out1_201020.rds"))
monit_out1 <- monit_out1[[1]]
monit_out2 <- readRDS(paste0(out_path, "a1_monit_out2_201020.rds"))
monit_out2 <- monit_out2[[1]]
monit_out3 <- readRDS(paste0(out_path, "a1_monit_out3_240920.rds"))
monit_out3 <- monit_out3[[1]]
monit_out4 <- readRDS(paste0(out_path, "a1_monit_out4_240920.rds"))
monit_out4 <- monit_out4[[1]]
monit_out5 <- readRDS(paste0(out_path, "a1_monit_out5_240920.rds"))
monit_out5 <- monit_out5[[1]]

# Monitoring different age groups until they age out (every 5 years)
monit_out6 <- readRDS(paste0(out_path, "a1_monit_out6_201020.rds"))
monit_out6 <- monit_out6[[1]]
monit_out7 <- readRDS(paste0(out_path, "a1_monit_out7_201020.rds"))
monit_out7 <- monit_out7[[1]]
monit_out8 <- readRDS(paste0(out_path, "a1_monit_out8_240920.rds"))
monit_out8 <- monit_out8[[1]]
monit_out9 <- readRDS(paste0(out_path, "a1_monit_out9_240920.rds"))
monit_out9 <- monit_out9[[1]]
monit_out10 <- readRDS(paste0(out_path, "a1_monit_out10_240920.rds"))
monit_out10 <- monit_out10[[1]]

# Labels of age groups being monitored
# sim1, sim6 = 15-30
# sim2, sim7 = 15-45
# sim3, sim8 = 30-45
# sim4, sim 9 = 30+
# sim5, sim10 = 45+

# Monitoring different age cohorts until they die (every year)
# A2 = screen and treat 45-60 year olds. A3 = screen and treat 30-60 year olds
a2_out3 <- readRDS(paste0(out_path, "a2_out3_screen_2020_monit_0_161020.rds"))
a2_out3 <- a2_out3[[1]]
a2_out5 <- readRDS(paste0(out_path, "a2_out5_screen_2020_monit_5_161020.rds"))
a2_out5 <- a2_out5[[1]]
a2_out6 <- readRDS(paste0(out_path, "a2_out6_screen_2020_monit_1_161020.rds"))
a2_out6 <- a2_out6[[1]]
a3_out3 <- readRDS(paste0(out_path, "a3_out3_screen_2020_monit_0_161020.rds"))
a3_out3 <- a3_out3[[1]]
a3_out5 <- readRDS(paste0(out_path, "a3_out5_screen_2020_monit_5_161020.rds"))
a3_out5 <- a3_out5[[1]]
a3_out6 <- readRDS(paste0(out_path, "a3_out6_screen_2020_monit_1_161020.rds"))
a3_out6 <- a3_out6[[1]]

scenario_labels <- list("No treatment" = "status_quo_cohort",
                        "No monitoring" = "screen_2020_monit_0",
                        "5-yearly all ages"="screen_2020_monit_5",
                        "Yearly all ages"="screen_2020_monit_1",
                        "Yearly 15-30"="screen_2020_monit_sim1",
                        "Yearly 15-45"="screen_2020_monit_sim2",
                        "Yearly 30-45"="screen_2020_monit_sim3",
                        "Yearly 30+"="screen_2020_monit_sim4",
                        "Yearly 45+"="screen_2020_monit_sim5",
                        "5-yearly 15-30"="screen_2020_monit_sim6",
                        "5-yearly 15-45"= "screen_2020_monit_sim7",
                        "5-yearly 30-45"="screen_2020_monit_sim8",
                        "5-yearly 30+"="screen_2020_monit_sim9",
                        "5-yearly 45+"="screen_2020_monit_sim10")

# Subset yearly frequencies
sub_yearly <- names(scenario_labels)[4:9]
# Subset 5-yearly frequencies
sub_5yearly <- names(scenario_labels)[c(3,10:14)]
# Selection for plot
sub_mixed <- names(scenario_labels)[c(3,4,5,8,10,13)]
# Since the best are 5-yearly 45+, 5-yearly 30+ and 5-yearly all ages, chose here <30, >30 and all

# Age group comparisons
sub_age_groups1 <- names(scenario_labels)[c(5,7,9,10,12,14)]
# 15-30, 30-45, 45+
sub_age_groups2 <- names(scenario_labels)[c(5,8,10,13)]
# 15-30, 30+
sub_age_groups3 <- names(scenario_labels)[c(6,9,11,14)]
# 15-45, 45+
sub_age_groups_all <- names(scenario_labels)[c(3,4)]

# 1) Optimal monitoring frequency across all ages ----
# Including status quo
#out3 = 0, out4 = 10, out5d = 9, out5c = 8, out5b = 7, out5a = 6, out5 = 5,
#out6c = 4, out6b = 3, out6a = 2, out6 =1
freq_interactions <- rbind(
  cbind(scenario = "No monitoring",
        out3$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 10 years",
        out4$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 9 years",
        out5d$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 8 years",
        out5c$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 7 years",
        out5b$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 6 years",
        out5a$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 5 years",
        out5$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 4 years",
        out6c$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 3 years",
        out6b$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 2 years",
        out6a$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Every 1 year",
        out6$interactions[[16]]$total_interactions[-c(1:3)]))

freq_interactions_long <- gather(freq_interactions, key = "sim",
                                 value = "interactions", - scenario)

freq_cohort_deaths_averted <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out1,
                                 scenario_objects = list(out3,
                                                         out4,
                                                         out5,
                                                         out5a,
                                                         out5b,
                                                         out5c,
                                                         out5d,
                                                         out6,
                                                         out6a,
                                                         out6b,
                                                         out6c),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")
monitoring_freq_label <- list("No monitoring"="screen_2020_monit_0",
                              "Every 10 years"="screen_2020_monit_10",
                              "Every 9 years"="screen_2020_monit_9",
                              "Every 8 years"="screen_2020_monit_8",
                              "Every 7 years"="screen_2020_monit_7",
                              "Every 6 years"="screen_2020_monit_6",
                              "Every 5 years"="screen_2020_monit_5",
                              "Every 4 years"="screen_2020_monit_4",
                              "Every 3 years"="screen_2020_monit_3",
                              "Every 2 years"="screen_2020_monit_2",
                              "Every 1 year"="screen_2020_monit_1")
levels(freq_cohort_deaths_averted$scenario) <- monitoring_freq_label

freq_cohort_deaths_averted <- filter(freq_cohort_deaths_averted, type == "number_averted") %>%
  select(scenario, sim, value)
freq_cohort_deaths_averted$sim <- paste0("X", as.numeric(gsub("[^0-9]", "", freq_cohort_deaths_averted$sim)))
# Check simulations match:
unique(freq_cohort_deaths_averted$sim) == colnames(freq_interactions[-1])

incremental_df_freq <- left_join(freq_interactions_long,
                                 freq_cohort_deaths_averted, by = c("scenario", "sim"))
colnames(incremental_df_freq)[4] <- "deaths_averted"

# Add status quo point at 0 (reference)
incremental_df_freq <- rbind(data.frame(scenario = rep("Status quo", 183),
                                        sim = unique(incremental_df_freq$sim),
                                        interactions = rep(0,183),
                                        deaths_averted = rep(0,183)),
                             incremental_df_freq)

# Find extended dominated strategies:
find_dominated_strategies(incremental_df_freq, "interactions", "deaths_averted")
# No monitoring, Every 10 years, Every 9 years
incremental_df_freq$frontier <- "Include"
incremental_df_freq$frontier[incremental_df_freq$scenario %in% c("No monitoring",
                                                                 "Every 10 years", "Every 9 years")] <-
  "Dominated"

incremental_df_freq_summary <- group_by(incremental_df_freq, scenario, frontier) %>%
  summarise(median_interactions = median(interactions),
            median_deaths_averted = median(deaths_averted))

# CORRECT PLOT showing all cumulative deaths averted and incremental interactions
# compared to status quo of no treatment
# The x and y value of a point shows the cumulative number of deaths averted and interactions
# invested compared to no treatment. The slope of the lines show the increment in this
# between consecutive strategies.
# However in interpreting this need to keep in mind the type of interactions varies between
# different strategies
ggplot(data = incremental_df_freq) +
  geom_line(data= subset(incremental_df_freq, frontier == "Include"),
            aes(y = deaths_averted, x = interactions, group = sim),
            colour = "grey", alpha = 0.3) +
  geom_point(aes(y = deaths_averted, x= interactions,
                 group = scenario, colour = scenario), alpha = 0.15) +
  stat_ellipse(data = incremental_df_freq[incremental_df_freq$scenario != "Status quo",],
               aes(y=deaths_averted,x=interactions,
                   group = scenario, fill= scenario),
               geom = "polygon",
               alpha = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_line(data = subset(incremental_df_freq_summary, frontier == "Include"),
            aes(x = median_interactions,
                                                    y = median_deaths_averted), size = 1) +
  geom_point(data = incremental_df_freq_summary, aes(x = median_interactions,
                                                     y = median_deaths_averted,
                                                     group = scenario, colour = scenario),
             size = 5) +
  geom_point(data = incremental_df_freq_summary, aes(x = median_interactions,
                                                     y = median_deaths_averted,
                                                     group = scenario),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual(values = brewer.pal(11,"RdYlBu")) +
  scale_colour_manual("Monitoring frequency",
                      values = c("black", brewer.pal(11,"RdYlBu")),
                      labels = c("Status quo" = "No treatment")) +
  guides(fill=FALSE) +
  xlab("Incremental total clinical interactions") +
  ylab("Incremental HBV-related deaths averted") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  coord_flip()
# Need to use coord_flip to get scenario lines connected in the right order because of some negative values

# Same for life-years
freq_cohort_ly_saved <-
  plot_ly_gained_cohort(counterfactual_object = out1,
                        scenario_objects = list(out3,
                                                out4,
                                                out5,
                                                out5a,
                                                out5b,
                                                out5c,
                                                out5d,
                                                out6,
                                                out6a,
                                                out6b,
                                                out6c),
                        outcome_to_plot = "number_averted",
                        counterfactual_label = "treatment programme without monitoring")
colnames(freq_cohort_ly_saved)[1:2] <- c("scenario", "counterfactual")
levels(freq_cohort_ly_saved$scenario) <- monitoring_freq_label

freq_cohort_ly_saved <- filter(freq_cohort_ly_saved, type == "number_averted") %>%
  select(scenario, sim, value)
freq_cohort_ly_saved$sim <- paste0("X", as.numeric(gsub("[^0-9]", "", freq_cohort_ly_saved$sim)))
# Check simulations match:
unique(freq_cohort_ly_saved$sim) == colnames(freq_interactions[-1])

incremental_df_ly_freq <- left_join(freq_interactions_long,
                                    freq_cohort_ly_saved, by = c("scenario", "sim"))
colnames(incremental_df_ly_freq)[4] <- "ly_saved"

# Add status quo point at 0 (reference)
incremental_df_ly_freq <- rbind(data.frame(scenario = rep("Status quo", 183),
                                           sim = unique(incremental_df_ly_freq$sim),
                                           interactions = rep(0,183),
                                           ly_saved = rep(0,183)),
                                incremental_df_ly_freq)

# Find extended dominated strategies:
find_dominated_strategies(incremental_df_ly_freq, "interactions", "ly_saved")
# No monitoring, Every 10 years, Every 9 years
incremental_df_ly_freq$frontier <- "Include"
incremental_df_ly_freq$frontier[incremental_df_freq$scenario %in% c("No monitoring")] <-
  "Dominated"


incremental_df_ly_freq_summary <- group_by(incremental_df_ly_freq, scenario, frontier) %>%
  summarise(median_interactions = median(interactions),
            median_ly_saved = median(ly_saved))

# CORRECT PLOT showing all cumulative LY gained and incremental interactions
# compared to status quo of no treatment
# However in interpreting this need to keep in mind the type of interactions varies between
# different strategies
ggplot(data = incremental_df_ly_freq) +
  geom_line(data=subset(incremental_df_ly_freq,frontier == "Include"),
            aes(y = ly_saved, x = interactions, group = sim),
            colour = "grey", alpha = 0.3) +
  geom_point(aes(y = ly_saved, x= interactions,
                 group = scenario, colour = scenario), alpha = 0.15) +
  stat_ellipse(data = incremental_df_ly_freq[incremental_df_ly_freq$scenario != "Status quo",],
               aes(y=ly_saved,x=interactions,
                   group = scenario, fill= scenario),
               geom = "polygon",
               alpha = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_line(data = subset(incremental_df_ly_freq_summary, frontier == "Include"),
            aes(x = median_interactions,
                y = median_ly_saved), size = 1) +
  geom_point(data = incremental_df_ly_freq_summary, aes(x = median_interactions,
                                                        y = median_ly_saved,
                                                        group = scenario, colour = scenario),
             size = 5) +
  geom_point(data = incremental_df_ly_freq_summary, aes(x = median_interactions,
                                                        y = median_ly_saved,
                                                        group = scenario),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual(values = brewer.pal(11,"RdYlBu")) +
  scale_colour_manual("Monitoring frequency",
                      values = c("black", brewer.pal(11,"RdYlBu")),
                      labels = c("Status quo" = "No treatment")) +
  guides(fill=FALSE) +
  xlab("Incremental total clinical interactions") +
  ylab("Incremental life-years saved") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  coord_flip()

# 2) Incremental benefits and interactions between a set of strategies (monitoring by age) ----

## Approach 1: monitoring given age group from entry until aging out ----
# Incrementally intensive strategies are: 15-30 years, 15-45 years, 15+ (all ages)

# INCREMENTAL INTERACTIONS IN EACH SCENARIO BY 2100 COMPARED TO NO MONITORING
# out3, monit_out6, monit_out7, out5 (5 yearly)

# Order could be:
# 5-yearly in 15-30 => yearly in 15-30 (better or worse?) => if better:
# yearly in 15-30 + 5-yearly in 30-45 => yearly in 15-30 + yearly in 30-45 etc
# Right now, going from Yearly 15-30 to 5-yearly 15-45, which makes no sense
# Can I derive the increment for this from my current simulations?
# However, interesting observation that going from yearly monitoring in 15-45
# to monitoring 5-yearly across all ages leads to more deaths averted at a lower cost!

# Create data frame with all interactions and outcomes of interest
interactions1 <- rbind(
  cbind(scenario = "screen_2020_monit_5",
        left_join(gather(out5$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(out5$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_1",
        left_join(gather(out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(out6$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim1",
        left_join(gather(monit_out1$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out1$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim2",
        left_join(gather(monit_out2$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out2$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim6",
        left_join(gather(monit_out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out6$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim7",
        left_join(gather(monit_out7$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out7$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim"))
)

# Extract person-years on treatment
interactions1_py_on_treatment <- rbind(
  data.frame(scenario = "screen_2020_monit_5",
             sim = names(out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_1",
             sim = names(out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment= out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim1",
             sim = names(monit_out1$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out1$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim2",
             sim = names(monit_out2$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out2$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim6",
             sim = names(monit_out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim7",
             sim = names(monit_out7$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = monit_out7$py_on_treatment[[16]]-out3$py_on_treatment[[16]]))

# Outcome 1: HBV related deaths averted
cohort_deaths_averted_long <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
                                 scenario_objects = list(out5,
                                                         out6,
                                                         monit_out1,
                                                         monit_out2,
                                                         monit_out6,
                                                         monit_out7),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "treatment programme without monitoring")
#levels(cohort_deaths_averted_long$scenario) <- scenario_labels
cohort_deaths_averted_long <- subset(cohort_deaths_averted_long, type == "number_averted") %>%
  select(scenario, sim, value)
cohort_deaths_averted_long$sim <- gsub("[^0-9]", "", cohort_deaths_averted_long$sim)

# Outcome 2: Life-years saved
cohort_ly_gained_long <-
  plot_ly_gained_cohort(counterfactual_object = out3,
                        scenario_objects = list(out5,
                                                out6,
                                                monit_out1,
                                                monit_out2,
                                                monit_out6,
                                                monit_out7),
                        outcome_to_plot = "number_averted",
                        counterfactual_label = "treatment programme without monitoring")
colnames(cohort_ly_gained_long)[1:2] <- c("scenario", "counterfactual")

cohort_ly_gained_long <- subset(cohort_ly_gained_long, type == "number_averted") %>%
  select(scenario, sim, value)
cohort_ly_gained_long$sim <- gsub("[^0-9]", "", cohort_ly_gained_long$sim)

# Combine into full dataframe
df1 <- create_incremental_plot_df(interactions_df=interactions1,
                                          py_on_treatment_df=interactions1_py_on_treatment,
                                          deaths_averted_df=cohort_deaths_averted_long,
                                          ly_saved_df = cohort_ly_gained_long,
                                  scenario_labels_obj = scenario_labels)

# Use BCEA package to find dominated and extended dominated strategies
find_dominated_strategies(df1, "total_interactions", "deaths_averted")
find_dominated_strategies(df1,"total_cost", "deaths_averted")
find_dominated_strategies(df1,"total_interactions", "ly_saved")
find_dominated_strategies(df1, "total_cost", "ly_saved")
# Deaths averted and interactions: Yearly 15-45 (A), Yearly 15-30 (E)
# Deaths averted and cost: Yearly 15-45 (A), Yearly 15-30 (E), 5-yearly 15-45 (E), 5-yearly 15-30 (E),
# LY saved and interactions: Yearly 15-45 (A), Yearly 15-30 (E)
# LY saved and cost: Yearly 15-45 (A), Yearly 15-30 (E), 5-yearly 15-30 (E)

# Add labels for this to dataframe
df1$frontier_interactions <- "Include"
df1$frontier_interactions[df1$scenario %in% c("Yearly 15-45", "Yearly 15-30")] <- "Dominated"
df1$frontier_deaths_averted_cost <- "Include"
df1$frontier_deaths_averted_cost[df1$scenario %in% c("Yearly 15-45", "Yearly 15-30",
                                                     "5-yearly 15-45", "5-yearly 15-30")] <- "Dominated"
df1$frontier_ly_saved_cost <- "Include"
df1$frontier_ly_saved_cost[df1$scenario %in% c("Yearly 15-45", "Yearly 15-30", "5-yearly 15-30")] <-
  "Dominated"

df1_summary <- df1 %>%
  group_by(scenario, frequency, age_group, frontier_interactions, frontier_deaths_averted_cost, frontier_ly_saved_cost) %>%
  summarise(median_deaths_averted = median(deaths_averted),
            median_ly_saved = median(ly_saved),
            median_interactions = median(total_interactions),
            median_cost = median(total_cost))

# Plots ----
# Plot: interactions vs deaths averted
ggplot(df1) +
  geom_line(data= subset(df1, frontier_interactions== "Include"),
            aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  stat_ellipse(data=subset(df1, scenario != "No monitoring"),
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted),
                   fill= reorder(scenario, deaths_averted)),
               geom = "polygon", na.rm = FALSE, alpha = 0.3) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  # Overlay median
  geom_line(data = subset(df1_summary, frontier_interactions == "Include"),
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual("Monitoring strategies\nstarting at entry\ninto cohort",
                    values=brewer.pal(6,"RdYlBu")) +
  scale_colour_manual("Monitoring strategies\nstarting at entry\ninto cohort",
                      values=c("black", brewer.pal(6,"RdYlBu"))) +
  guides(fill=FALSE) +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs deaths averted (separated by frequency)
ggplot(df1) +
  geom_line(aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = df1_summary,
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~frequency) +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs deaths averted
ggplot(df1) +
  geom_line(data= subset(df1, frontier_deaths_averted_cost == "Include"),
            aes(x = deaths_averted, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_cost,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_cost,
                   group = reorder(scenario, deaths_averted),
                   fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df1_summary, frontier_deaths_averted_cost == "Include"),
            aes(y = median_cost,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_cost,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_cost,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental cost (USD)") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs LY saved
ggplot(df1) +
  geom_line(data= subset(df1, frontier_interactions == "Include"),
            aes(x = ly_saved, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_interactions,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_interactions,
                   group = reorder(scenario, ly_saved), fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df1_summary, frontier_interactions == "Include"),
            aes(y = median_interactions,
                x = median_ly_saved), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental number of clinical interactions") +
  #  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs LY saved
ggplot(df1) +
  geom_line(data= subset(df1, frontier_ly_saved_cost == "Include"),
            aes(x = ly_saved, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_cost,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_cost,
                   group = reorder(scenario, ly_saved),
                   fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df1_summary, frontier_ly_saved_cost == "Include"),
            aes(y = median_cost,
                x = median_ly_saved), size = 1) +
  geom_point(data = df1_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df1_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  #  geom_abline(intercept =0, slope = 240, linetype = "dashed") +   # WTP: need to compare to status quo?
  #  geom_abline(intercept =0, slope = 1460, linetype = "dashed") +  # WTP: need to compare to status quo?
  #  geom_abline(intercept =0, slope = 487, linetype = "dashed") +   # WTP: need to compare to status quo?
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental cost (USD)") +
  # xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

## Approach 2: monitoring given age cohort from entry until death ----
# Incrementally intensive strategies are: 45+ years, 30+ years, 15+ (all ages)
# Start with 5-yearly frequency
# Increment between no monitoring and status quo: out3-out1
# Increment between monitoring 45+ year olds and no monitoring in any age group: a2_out5-a2_out3
# Increment between adding 30-45 year old to monitoring vs only 45+ year olds:
# (a3_out5-a3_out3)-(a2_out5-a2_out3)
# Increment between adding 30-45 year old to monitoring vs only 30+ year olds:
#(a1_out5-a1_out3)-(a3_out5-a3_out3)

# Create data frame with all interactions and outcomes of interest
interactions2 <- rbind(
  cbind(scenario = "5-yearly 15+ (all ages)",
        left_join(gather(out5$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(out5$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Yearly 15+ (all ages)",
        left_join(gather(out6$interactions[[16]]$total_assessed[-c(1:3)]-out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(out6$interactions[[16]]$total_treated[-c(1:3)]-out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "5-yearly 45+",
        left_join(gather(a2_out5$interactions[[16]]$total_assessed[-c(1:3)]-a2_out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(a2_out5$interactions[[16]]$total_treated[-c(1:3)]-a2_out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "5-yearly 30+",
        left_join(gather(a3_out5$interactions[[16]]$total_assessed[-c(1:3)]-a3_out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(a3_out5$interactions[[16]]$total_treated[-c(1:3)]-a3_out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Yearly 45+",
        left_join(gather(a2_out6$interactions[[16]]$total_assessed[-c(1:3)]-a2_out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(a2_out6$interactions[[16]]$total_treated[-c(1:3)]-a2_out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "Yearly 30+",
        left_join(gather(a3_out6$interactions[[16]]$total_assessed[-c(1:3)]-a3_out3$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(a3_out6$interactions[[16]]$total_treated[-c(1:3)]-a3_out3$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim"))
)

# Extract person-years on treatment
interactions2_py_on_treatment <- rbind(
  data.frame(scenario = "5-yearly 15+ (all ages)",
             sim = names(out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment = out5$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "Yearly 15+ (all ages)",
             sim = names(out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
             py_on_treatment= out6$py_on_treatment[[16]]-out3$py_on_treatment[[16]]),
  data.frame(scenario = "5-yearly 45+",
             sim = names(a2_out5$py_on_treatment[[16]]-a2_out3$py_on_treatment[[16]]),
             py_on_treatment = a2_out5$py_on_treatment[[16]]-a2_out3$py_on_treatment[[16]]),
  data.frame(scenario = "5-yearly 30+",
             sim = names(a3_out5$py_on_treatment[[16]]-a3_out3$py_on_treatment[[16]]),
             py_on_treatment = a3_out5$py_on_treatment[[16]]-a3_out3$py_on_treatment[[16]]),
  data.frame(scenario = "Yearly 45+",
             sim = names(a2_out6$py_on_treatment[[16]]-a2_out3$py_on_treatment[[16]]),
             py_on_treatment = a2_out6$py_on_treatment[[16]]-a2_out3$py_on_treatment[[16]]),
  data.frame(scenario = "Yearly 30+",
             sim = names(a3_out6$py_on_treatment[[16]]-a3_out3$py_on_treatment[[16]]),
             py_on_treatment = a3_out6$py_on_treatment[[16]]-a3_out3$py_on_treatment[[16]]))

# Outcome 1: HBV related deaths averted
cohort_deaths_averted2 <- rbind(
  cbind(scenario = "5-yearly 45+",
        a2_out3$cum_hbv_deaths[[16]][-c(1:3)]-a2_out5$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "5-yearly 30+",
        a3_out3$cum_hbv_deaths[[16]][-c(1:3)]-a3_out5$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "5-yearly 15+ (all ages)",
        out3$cum_hbv_deaths[[16]][-c(1:3)]-out5$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 45+",
        a2_out3$cum_hbv_deaths[[16]][-c(1:3)]-a2_out6$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 30+",
        a3_out3$cum_hbv_deaths[[16]][-c(1:3)]-a3_out6$cum_hbv_deaths[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 15+ (all ages)",
        out3$cum_hbv_deaths[[16]][-c(1:3)]-out6$cum_hbv_deaths[[16]][-c(1:3)]))
colnames(cohort_deaths_averted2)[-1] <- as.numeric(gsub("[^0-9]", "",
                                                        colnames(cohort_deaths_averted2)[-1]))
cohort_deaths_averted2_long <- gather(cohort_deaths_averted2, key = "sim",
                                      value = "deaths_averted", - scenario)

# Outcome 2: Life-years saved

# Life-years saved
cohort_ly_gained2 <- rbind(
  cbind(scenario = "5-yearly 45+",
        a2_out5$ly[[16]][-c(1:3)]-a2_out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "5-yearly 30+",
        a3_out5$ly[[16]][-c(1:3)]-a3_out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "5-yearly 15+ (all ages)",
        out5$ly[[16]][-c(1:3)]-out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 45+",
        a2_out6$ly[[16]][-c(1:3)]-a2_out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 30+",
        a3_out6$ly[[16]][-c(1:3)]-a3_out3$ly[[16]][-c(1:3)]),
  cbind(scenario = "Yearly 15+ (all ages)",
        out6$ly[[16]][-c(1:3)]-out3$ly[[16]][-c(1:3)]))
colnames(cohort_ly_gained2)[-1] <- as.numeric(gsub("[^0-9]", "",
                                                   colnames(cohort_ly_gained2)[-1]))
cohort_ly_gained2_long <- gather(cohort_ly_gained2, key = "sim",
                                 value = "ly_saved", - scenario)

# Combine into full dataframe
df2 <- create_incremental_plot_df(interactions_df=interactions2,
                                  py_on_treatment_df=interactions2_py_on_treatment,
                                  deaths_averted_df=cohort_deaths_averted2_long,
                                  ly_saved_df = cohort_ly_gained2_long)

# Find dominated strategies for each combination of exposure and outcome using BCEA package
find_dominated_strategies(df2, "total_interactions", "deaths_averted")
find_dominated_strategies(df2,"total_cost", "deaths_averted")
find_dominated_strategies(df2,"total_interactions", "ly_saved")
find_dominated_strategies(df2, "total_cost", "ly_saved")
# Deaths averted & interactions, LY saved & any exp.:
# Yearly 30+, Yearly 45+, 5-yearly 30+, 5-yearly 45+
# Deaths averted & cost: Yearly 30+, 5-yearly 45+, Yearly 45+

# Add labels to this dataframe
df2$frontier_deaths_averted_interactions <- "Include"
df2$frontier_deaths_averted_interactions[df2$scenario %in% c("Yearly 30+", "5-yearly 30+",
                                                             "5-yearly 45+", "Yearly 45+")] <- "Dominated"
df2$frontier_deaths_averted_cost <- "Include"
df2$frontier_deaths_averted_cost[df2$scenario %in% c("Yearly 30+", "5-yearly 45+", "Yearly 45+")] <-
  "Dominated"
df2$frontier_ly_saved <- "Include"
df2$frontier_ly_saved[df2$scenario %in% c("Yearly 30+", "5-yearly 30+",
                                          "5-yearly 45+", "Yearly 45+")] <- "Dominated"

df2_summary <- df2 %>%
  group_by(scenario, frequency, age_group, frontier_deaths_averted_interactions,
           frontier_deaths_averted_cost, frontier_ly_saved) %>%
  summarise(median_deaths_averted = median(deaths_averted),
            median_ly_saved = median(ly_saved),
            median_interactions = median(total_interactions),
            median_cost = median(total_cost))

# Plots ----
# Plot: interactions vs deaths averted (combined)
ggplot(df2) +
  geom_line(data= subset(df2, frontier_deaths_averted_interactions == "Include"),
            aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)),
               alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df2_summary, frontier_deaths_averted_interactions == "Include"),
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  #xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs deaths averted (separated by birth cohort)
ggplot(df2) +
  geom_line(data= df2,
            aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = df2_summary,
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~age_group) +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  #xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs deaths averted (separate by frequency)
ggplot(df2) +
  geom_line(data= df2,
            aes(x = deaths_averted, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_interactions,
                   group = reorder(scenario, deaths_averted), fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = df2_summary,
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~frequency) +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs deaths averted
ggplot(df2) +
  geom_line(data= subset(df2, frontier_cost == "Include"),
            aes(x = deaths_averted, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = deaths_averted, y = total_cost,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=deaths_averted,y=total_cost,
                   group = reorder(scenario, deaths_averted),
                   fill= reorder(scenario, deaths_averted)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df2_summary, frontier_cost == "Include"),
            aes(y = median_cost,
                x = median_deaths_averted), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental cost (USD)") +
  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: interactions vs LY saved
ggplot(df2) +
  geom_line(data= subset(df2, frontier_interactions == "Include"),
            aes(x = ly_saved, y= total_interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_interactions,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)), alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_interactions,
                   group = reorder(scenario, ly_saved), fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df2_summary, frontier_interactions == "Include"),
            aes(y = median_interactions,
                x = median_ly_saved), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental number of clinical interactions") +
  #  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs LY saved
ggplot(df2) +
  geom_line(data= subset(df2, frontier_cost == "Include"),
            aes(x = ly_saved, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_cost,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_cost,
                   group = reorder(scenario, ly_saved),
                   fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = subset(df2_summary, frontier_cost == "Include"),
            aes(y = median_cost,
                x = median_ly_saved), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  geom_abline(intercept =0, slope = 240, linetype = "dashed") +
  geom_abline(intercept =0, slope = 1460, linetype = "dashed") +
  geom_abline(intercept =0, slope = 487, linetype = "dashed") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental cost (USD)") +
  # xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# Plot: cost vs LY saved (separated by frequency - no dominated strategies)
ggplot(df2) +
  geom_line(data= df2,
            aes(x = ly_saved, y= total_cost,
                group = sim), colour = "grey", alpha = 0.5) +
  geom_point(aes(x = ly_saved, y = total_cost,
                 group =reorder(scenario, ly_saved), colour = reorder(scenario, ly_saved)),
             alpha = 0.4) +
  stat_ellipse(geom = "polygon",
               aes(x=ly_saved,y=total_cost,
                   group = reorder(scenario, ly_saved),
                   fill= reorder(scenario, ly_saved)), alpha = 0.3) +
  # Overlay median
  geom_line(data = df2_summary,
            aes(y = median_cost,
                x = median_ly_saved), size = 1) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved),
                 colour = reorder(scenario, median_ly_saved)), size = 5) +
  geom_point(data = df2_summary,
             aes(y = median_cost,
                 x = median_ly_saved,
                 group = reorder(scenario, median_ly_saved)),
             size = 5, shape = 1, colour = "black") +
  facet_wrap(~frequency) +
  geom_abline(intercept =0, slope = 240, linetype = "dashed") +
  geom_abline(intercept =0, slope = 1460, linetype = "dashed") +
  geom_abline(intercept =0, slope = 487, linetype = "dashed") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental cost (USD)") +
  # xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# This plot shows the increase in deaths averted and # total interactions between each 2 strategies,
# with "treatment without monitoring" being at 0 (status quo of no treatment is shown as negative to that).
# The assumption here is that everyone is monitored from the start, but it stops at some age (e.g. 30, 45).
# Would ideally exclude the lines/points outside the 95% interval

## Combine approaches 1 and 2 ----
df1 <- select(df, scenario, sim, interactions, value) %>%
  filter(scenario %in% c("No treatment", "No monitoring",
                         "5-yearly 15-30", "Yearly 15-30",
                         "5-yearly 15-45", "5-yearly all ages",
                         "Yearly all ages",  "Yearly 15-45"))
colnames(df1)[4] <- "deaths_averted"
df1$scenario <- factor(paste("By age:", df1$scenario))
# Remove duplicate scenarios (all ages) from approach 2:
df2 <- filter(age_cohort_df, scenario != "5-yearly 15+ (all ages)" & scenario != "Yearly 15+ (all ages)")
df2$scenario <- factor(paste("By birth cohort:", df2$scenario))
combined_df <- rbind(df1, df2)
# Label dominated strategies (based on median):
combined_df$scenario <- as.character(combined_df$scenario)
combined_df$scenario <- factor(combined_df$scenario,
                               levels = c("By birth cohort: No monitoring",
                                          "By age: 5-yearly 15-30",
                                          "By birth cohort: 5-yearly 45+",
                                          "By age: Yearly 15-30",
                                          "By age: 5-yearly 15-45",
                                          "By birth cohort: 5-yearly 30+",
                                          "By age: 5-yearly all ages",
                                          "By age: Yearly all ages",
                                          "By birth cohort: Yearly 45+",
                                          "By age: Yearly 15-45",
                                          "By birth cohort: Yearly 30+"))
levels(combined_df$scenario) <- list("No monitoring" = "By birth cohort: No monitoring",
                                     "By age: 15-30,\nevery 5 years" = "By age: 5-yearly 15-30",
                                     "By birth cohort: 45+,\nevery 5 years" = "By birth cohort: 5-yearly 45+",
                                     "By age: 15-30,\nevery 1 year" = "By age: Yearly 15-30",
                                     "By age: 15-45,\nevery 5 years" = "By age: 5-yearly 15-45",
                                     "By birth cohort:\n30+, every 5 years" = "By birth cohort: 5-yearly 30+",
                                     "All ages,\nevery 5 years" = "By age: 5-yearly all ages",
                                     "All ages,\nevery 1 year" = "By age: Yearly all ages",
                                     "By birth cohort: 45+,\nevery 1 year" = "By birth cohort: Yearly 45+",
                                     "By age: 15-45,\nevery 1 year" = "By age: Yearly 15-45",
                                     "By birth cohort: 30+,\nevery 1 year" = "By birth cohort: Yearly 30+")

combined_df$frontier <- "Include"
combined_df$frontier[combined_df$scenario %in% c("By birth cohort: 45+,\nevery 1 year",
                                                 "By age: 15-45,\nevery 1 year",
                                                 "By birth cohort: 30+,\nevery 1 year")] <- "Dominated"


# Summary:
combined_df_summary <- group_by(combined_df, scenario, frontier) %>%
  summarise(median_interactions = median(interactions),
            median_deaths_averted = median(deaths_averted))

# Plot:
ggplot(combined_df) +
  geom_point(aes(x = deaths_averted,interactions,
                 group = scenario, colour = scenario), alpha = 0.4) +
  geom_line(data=subset(combined_df, frontier == "Include"),
            aes(x = deaths_averted,interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  stat_ellipse(geom = "polygon",
               aes(deaths_averted,interactions,
                   group = scenario, fill= scenario), alpha = 0.3) +
  geom_line(data = subset(combined_df_summary, frontier == "Include"),
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = combined_df_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = scenario, colour = scenario), size = 5) +
  geom_point(data = combined_df_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = scenario),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  #  xlim(-100,7100) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

# This one might look quite different in life-years or if using costs (treatment time)

# Life-years saved
df_ly1 <- select(df_ly, scenario, sim, interactions, value) %>%
  filter(scenario %in% c("No treatment", "No monitoring",
                         "5-yearly 15-30", "Yearly 15-30",
                         "5-yearly 15-45", "5-yearly all ages",
                         "Yearly all ages",  "Yearly 15-45"))
colnames(df_ly1)[4] <- "ly_saved"
df_ly1$scenario <- factor(paste("By age:", df_ly1$scenario))
# Remove duplicate scenarios (all ages) from approach 2:
df_ly2 <- filter(age_cohort_ly_df, scenario != "5-yearly 15+ (all ages)" &
                   scenario != "Yearly 15+ (all ages)")
df_ly2$scenario <- factor(paste("By birth cohort:", df_ly2$scenario))
combined_df_ly <- rbind(df_ly1, df_ly2)
# Label dominated strategies (based on median):
combined_df_ly$scenario <- as.character(combined_df_ly$scenario)
combined_df_ly$scenario <- factor(combined_df_ly$scenario,
                                  levels = c("By birth cohort: No monitoring",
                                             "By age: 5-yearly 15-30",
                                             "By birth cohort: 5-yearly 45+",
                                             "By age: Yearly 15-30",
                                             "By age: 5-yearly 15-45",
                                             "By birth cohort: 5-yearly 30+",
                                             "By age: 5-yearly all ages",
                                             "By age: Yearly all ages",
                                             "By birth cohort: Yearly 45+",
                                             "By age: Yearly 15-45",
                                             "By birth cohort: Yearly 30+"))
levels(combined_df_ly$scenario) <- list("No monitoring" = "By birth cohort: No monitoring",
                                        "By age: 15-30,\nevery 5 years" = "By age: 5-yearly 15-30",
                                        "By birth cohort: 45+,\nevery 5 years" = "By birth cohort: 5-yearly 45+",
                                        "By age: 15-30,\nevery 1 year" = "By age: Yearly 15-30",
                                        "By age: 15-45,\nevery 5 years" = "By age: 5-yearly 15-45",
                                        "By birth cohort:\n30+, every 5 years" = "By birth cohort: 5-yearly 30+",
                                        "All ages,\nevery 5 years" = "By age: 5-yearly all ages",
                                        "All ages,\nevery 1 year" = "By age: Yearly all ages",
                                        "By birth cohort: 45+,\nevery 1 year" = "By birth cohort: Yearly 45+",
                                        "By age: 15-45,\nevery 1 year" = "By age: Yearly 15-45",
                                        "By birth cohort: 30+,\nevery 1 year" = "By birth cohort: Yearly 30+")

combined_df_ly$frontier <- "Include"
combined_df_ly$frontier[combined_df_ly$scenario %in% c("By birth cohort: 45+,\nevery 1 year",
                                                       "By age: 15-45,\nevery 1 year",
                                                       "By birth cohort: 30+,\nevery 1 year",
                                                       "By birth cohort: 45+,\nevery 5 years",
                                                       "By birth cohort:\n30+, every 5 years")] <- "Dominated"


# Summary:
combined_df_ly_summary <- group_by(combined_df_ly, scenario, frontier) %>%
  summarise(median_interactions = median(interactions),
            median_ly_saved = median(ly_saved))

# Plot:
ggplot(combined_df_ly) +
  geom_point(aes(x = ly_saved,interactions,
                 group = scenario, colour = scenario), alpha = 0.4) +
  geom_line(data=subset(combined_df_ly, frontier == "Include"),
            aes(x = ly_saved,interactions,
                group = sim), colour = "grey", alpha = 0.5) +
  stat_ellipse(geom = "polygon",
               aes(ly_saved,interactions,
                   group = scenario, fill= scenario), alpha = 0.3) +
  geom_line(data = subset(combined_df_ly_summary, frontier == "Include"),
            aes(y = median_interactions,
                x = median_ly_saved), size = 1) +
  geom_point(data = combined_df_ly_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = scenario, colour = scenario), size = 5) +
  geom_point(data = combined_df_ly_summary,
             aes(y = median_interactions,
                 x = median_ly_saved,
                 group = scenario),
             size = 5, shape = 1, colour = "black") +
  scale_fill_brewer("Monitoring strategies\nstarting at entry\ninto cohort", palette="RdYlBu") +
  scale_colour_brewer("Monitoring strategies\nstarting at entry\ninto cohort",
                      palette="RdYlBu") +
  xlab("Incremental life-years saved") +
  ylab("Incremental number of clinical interactions") +
  #  xlim(-100,7100) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))


# BCEA package test ----
library(BCEA)
?bcea

interactions_test <- t(interactions_by_2100[,-1])
interactions_test <- cbind(interactions_test[,2]-interactions_test[,1],
                           interactions_test[,3]-interactions_test[,1],
                           interactions_test[,4]-interactions_test[,1])
interactions_test <- cbind(interactions_test, rep(0,183))
outcome_test <- t(deaths_averted_by_2100[,-1])
outcome_test <- cbind(outcome_test[,2]-outcome_test[,1],
                      outcome_test[,3]-outcome_test[,1],
                      outcome_test[,4]-outcome_test[,1])
outcome_test <- cbind(outcome_test, rep(0,183))


test <- bcea(e=outcome_test,c=interactions_test,
             ref=4,
             interventions=c(as.character(interactions_by_2100$scenario[-1]), "No monitoring"),
             Kmax=50000,
             plot=TRUE)

ceef.plot(test, graph="base", relative = FALSE)

#
ceef_interactions_by_2100 <- rbind(
  cbind(scenario = "No monitoring",
        out3$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "5-Yearly 15-30",
        monit_out6$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "5-Yearly 15-45",
        monit_out7$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "5-Yearly 15+ (all ages)",
        out5$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Yearly 15-30",
        monit_out1$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Yearly 15-45",
        monit_out2$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]),
  cbind(scenario = "Yearly 15+ (all ages)",
        out6$interactions[[16]]$total_interactions[-c(1:3)]-out3$interactions[[16]]$total_interactions[-c(1:3)]))
ceef_interactions <- t(ceef_interactions_by_2100[,-1])

ceef_deaths_averted_by_2100 <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3,
                                 scenario_objects = list(out3, monit_out6, monit_out7, out5,
                                                         monit_out1, monit_out2, out6),
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "no monitoring")

ceef_deaths_averted_by_2100<- filter(ceef_deaths_averted_by_2100, type == "number_averted") %>%
  select(scenario, sim, value) %>%
  spread(key = "scenario", value = "value") %>%
  select(-sim)

test2 <- bcea(e=as.matrix(ceef_deaths_averted_by_2100),c=ceef_interactions,
              ref=1,
              interventions=c(as.character(ceef_interactions_by_2100$scenario)),
              Kmax=50000,
              plot=TRUE)

ceef.plot(test2, graph="base", relative = FALSE)




