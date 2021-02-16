# Analyse the effect of screening/monitoring/treatment in different age groups

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
#source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# For all these plots, check values line up with total!!!

# Load data ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/"

out_path_monit <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

# Screening (no monitoring) in separate age groups
# A4 = 15-30
# A5 = 30-45
# A2 = 45-65
# All simulations involve IT treatment
a4_out3 <- readRDS(paste0(out_path, "a4_it_out3_screen_2020_monit_0_treatment_effect_270121.rds"))
a4_out3 <- a4_out3[[1]]
a5_out3 <- readRDS(paste0(out_path, "a5_it_out3_screen_2020_monit_0_treatment_effect_270121.rds"))
a5_out3 <- a5_out3[[1]]
a2_out3 <- readRDS(paste0(out_path, "a2_it_out3_screen_2020_monit_0_treatment_effect_270121.rds"))
a2_out3 <- a2_out3[[1]]

# Same cohort simulations without treatment:
a4_out1 <- readRDS(paste0(out_path, "a4_it_out1_screen_2020_monit_0_cohort_treatment_effect_270121.rds"))
a4_out1 <- a4_out1[[1]]
a5_out1 <- readRDS(paste0(out_path, "a5_it_out1_screen_2020_monit_0_cohort_treatment_effect_270121.rds"))
a5_out1 <- a5_out1[[1]]
a2_out1 <- readRDS(paste0(out_path, "a2_it_out1_screen_2020_monit_0_cohort_treatment_effect_270121.rds"))
a2_out1 <- a2_out1[[1]]

# Normal population simulations of effect of screening by age
out2 <- readRDS(paste0(out_path_monit, "out2_status_quo_301120.rds"))
out2 <- out2[[1]]
# Carriers by age in sq scenario:
out2_carriers <- readRDS(paste0(out_path, "out_sq_carriers.rds"))
out2_comps_by_age <- readRDS(paste0(out_path, "out_sq_compartment_prevalence_by_age.rds"))
da <- 0.5
ages <- seq(0,99.5,0.5)

a2_out3_pop <- readRDS(paste0(out_path_monit, "a2_it_out3_screen_2020_monit_0_180121.rds"))
a2_out3_pop <- a2_out3_pop[[1]]
a4_out3_pop <- readRDS(paste0(out_path_monit, "a4_it_out3_screen_2020_monit_0_190121.rds"))
a4_out3_pop <- a4_out3_pop[[1]]
a5_out3_pop <- readRDS(paste0(out_path_monit, "a5_it_out3_screen_2020_monit_0_190121.rds"))
a5_out3_pop <- a5_out3_pop[[1]]

# With/without monitoring
a1_out3_pop <- readRDS(paste0(out_path_monit, "a1_it_out3_screen_2020_monit_0_180121.rds"))
a1_out3_pop <- a1_out3_pop[[1]]   # No monitoring
a1_out6_pop <- readRDS(paste0(out_path_monit, "a1_it_out6_screen_2020_monit_1_130121.rds"))
a1_out6_pop <- a1_out6_pop[[1]]  # 1 year


# Simulate the effect of treatment in the screened+treated cohort ----
# Total incident HBV/HCC deaths averted over time ----
# Prepare time vector used in simulations
timevec <- seq(1960, 2120.5,0.5)

# Create dataframe of deaths averted
averted_outcomes <- rbind(
  data.frame(outcome="hbv_deaths_averted", age="15-30", sex="male", time=timevec,
           a4_out1$inc_cohort_hbv_deaths_male_over_time-
             a4_out3$inc_cohort_hbv_deaths_male_over_time),
  data.frame(outcome="hbv_deaths_averted", age="15-30", sex="female", time=timevec,
             a4_out1$inc_cohort_hbv_deaths_female_over_time-
               a4_out3$inc_cohort_hbv_deaths_female_over_time),
  data.frame(outcome="hcc_deaths_averted", age="15-30", sex="male", time=timevec,
             a4_out1$inc_cohort_hcc_deaths_male_over_time-
               a4_out3$inc_cohort_hcc_deaths_male_over_time),
  data.frame(outcome="hcc_deaths_averted", age="15-30", sex="female", time=timevec,
             a4_out1$inc_cohort_hcc_deaths_female_over_time-
               a4_out3$inc_cohort_hcc_deaths_female_over_time),
  # Cirrhosis deaths = total HBV deaths - HCC deaths
  data.frame(outcome="cirrhosis_deaths_averted", age="15-30", sex="male", time=timevec,
             (a4_out1$inc_cohort_hbv_deaths_male_over_time-a4_out1$inc_cohort_hcc_deaths_male_over_time)-
               (a4_out3$inc_cohort_hbv_deaths_male_over_time-a4_out3$inc_cohort_hcc_deaths_male_over_time)),
  data.frame(outcome="cirrhosis_deaths_averted", age="15-30", sex="female", time=timevec,
             (a4_out1$inc_cohort_hbv_deaths_female_over_time-a4_out1$inc_cohort_hcc_deaths_female_over_time)-
               (a4_out3$inc_cohort_hbv_deaths_female_over_time-a4_out3$inc_cohort_hcc_deaths_female_over_time)),
  data.frame(outcome="hcc_cases_averted", age="15-30", sex="male", time=timevec,
             a4_out1$inc_cohort_hcc_cases_male_over_time-
               a4_out3$inc_cohort_hcc_cases_male_over_time),
  data.frame(outcome="hcc_cases_averted", age="15-30", sex="female", time=timevec,
             a4_out1$inc_cohort_hcc_cases_female_over_time-
               a4_out3$inc_cohort_hcc_cases_female_over_time),
  data.frame(outcome="hbv_deaths_averted", age="30-45", sex="male", time=timevec,
             a5_out1$inc_cohort_hbv_deaths_male_over_time-
               a5_out3$inc_cohort_hbv_deaths_male_over_time),
  data.frame(outcome="hbv_deaths_averted", age="30-45", sex="female", time=timevec,
             a5_out1$inc_cohort_hbv_deaths_female_over_time-
               a5_out3$inc_cohort_hbv_deaths_female_over_time),
  data.frame(outcome="hcc_deaths_averted", age="30-45", sex="male", time=timevec,
             a5_out1$inc_cohort_hcc_deaths_male_over_time-
               a5_out3$inc_cohort_hcc_deaths_male_over_time),
  data.frame(outcome="hcc_deaths_averted", age="30-45", sex="female", time=timevec,
             a5_out1$inc_cohort_hcc_deaths_female_over_time-
               a5_out3$inc_cohort_hcc_deaths_female_over_time),
  data.frame(outcome="cirrhosis_deaths_averted", age="30-45", sex="male", time=timevec,
             (a5_out1$inc_cohort_hbv_deaths_male_over_time-a5_out1$inc_cohort_hcc_deaths_male_over_time)-
               (a5_out3$inc_cohort_hbv_deaths_male_over_time-a5_out3$inc_cohort_hcc_deaths_male_over_time)),
  data.frame(outcome="cirrhosis_deaths_averted", age="30-45", sex="female", time=timevec,
             (a5_out1$inc_cohort_hbv_deaths_female_over_time-a5_out1$inc_cohort_hcc_deaths_female_over_time)-
               (a5_out3$inc_cohort_hbv_deaths_female_over_time-a5_out3$inc_cohort_hcc_deaths_female_over_time)),
  data.frame(outcome="hcc_cases_averted", age="30-45", sex="male", time=timevec,
             a5_out1$inc_cohort_hcc_cases_male_over_time-
               a5_out3$inc_cohort_hcc_cases_male_over_time),
  data.frame(outcome="hcc_cases_averted", age="30-45", sex="female", time=timevec,
             a5_out1$inc_cohort_hcc_cases_female_over_time-
               a5_out3$inc_cohort_hcc_cases_female_over_time),
  data.frame(outcome="hbv_deaths_averted", age="45-65", sex="male", time=timevec,
             a2_out1$inc_cohort_hbv_deaths_male_over_time-
               a2_out3$inc_cohort_hbv_deaths_male_over_time),
  data.frame(outcome="hbv_deaths_averted", age="45-65", sex="female", time=timevec,
             a2_out1$inc_cohort_hbv_deaths_female_over_time-
               a2_out3$inc_cohort_hbv_deaths_female_over_time),
  data.frame(outcome="hcc_deaths_averted", age="45-65", sex="male", time=timevec,
             a2_out1$inc_cohort_hcc_deaths_male_over_time-
               a2_out3$inc_cohort_hcc_deaths_male_over_time),
  data.frame(outcome="hcc_deaths_averted", age="45-65", sex="female", time=timevec,
             a2_out1$inc_cohort_hcc_deaths_female_over_time-
               a2_out3$inc_cohort_hcc_deaths_female_over_time),
  data.frame(outcome="cirrhosis_deaths_averted", age="45-65", sex="male", time=timevec,
             (a2_out1$inc_cohort_hbv_deaths_male_over_time-a2_out1$inc_cohort_hcc_deaths_male_over_time)-
               (a2_out3$inc_cohort_hbv_deaths_male_over_time-a2_out3$inc_cohort_hcc_deaths_male_over_time)),
  data.frame(outcome="cirrhosis_deaths_averted", age="45-65", sex="female", time=timevec,
             (a2_out1$inc_cohort_hbv_deaths_female_over_time-a2_out1$inc_cohort_hcc_deaths_female_over_time)-
               (a2_out3$inc_cohort_hbv_deaths_female_over_time-a2_out3$inc_cohort_hcc_deaths_female_over_time)),
  data.frame(outcome="hcc_cases_averted", age="45-65", sex="male", time=timevec,
             a2_out1$inc_cohort_hcc_cases_male_over_time-
               a2_out3$inc_cohort_hcc_cases_male_over_time),
  data.frame(outcome="hcc_cases_averted", age="45-65", sex="female", time=timevec,
             a2_out1$inc_cohort_hcc_cases_female_over_time-
               a2_out3$inc_cohort_hcc_cases_female_over_time)
)
# Remove irrelevant timesteps
averted_outcomes <-filter(averted_outcomes, time>=2015 & time <=2100)

averted_outcomes_long <- gather(averted_outcomes, key = "sim", value = "value",
                                -outcome, -age, -sex,-time)


# New facet label names
facet_labels_deaths <- c("Averted cirrhosis deaths", "Averted HCC deaths")
names(facet_labels_deaths) <- c("cirrhosis_deaths_averted", "hcc_deaths_averted")

#averted_outcomes_long$y_min[averted_outcomes_long$sex=="male" &
#                        averted_outcomes_long$outcome == "cirrhosis_deaths_averted"] <- -14
#averted_outcomes_long$y_max[averted_outcomes_long$sex=="male" &
#                              averted_outcomes_long$outcome == "cirrhosis_deaths_averted"] <- 150
#averted_outcomes_long$y_min[averted_outcomes_long$sex=="male" &
#                              averted_outcomes_long$outcome == "hcc_deaths_averted"] <- -14
#averted_outcomes_long$y_max[averted_outcomes_long$sex=="male" &
#                              averted_outcomes_long$outcome == "hcc_deaths_averted"] <- 51

# Men
ggplot(subset(averted_outcomes_long,sex == "male" & outcome %in%
                c("cirrhosis_deaths_averted", "hcc_deaths_averted"))) +
  geom_line(aes(x=time, y = value/0.5, group = sim), col = "grey") +
  stat_summary(aes(x=time, y = value), fun="median", geom = "line", col = "red")+
  facet_wrap(outcome~age,ncol = 3, scales = "free_y",
             labeller = labeller(outcome = facet_labels_deaths)) +
  ylab("HBV-related deaths averted per year") +
  geom_hline(yintercept=0) +
  labs(title="HBV-related deaths averted by treatment in different age groups in men") +
#  geom_blank(aes(y = y_min)) +
#  geom_blank(aes(y = y_max)) +
  theme_classic()


# Women
ggplot(subset(averted_outcomes_long,sex == "female" & outcome %in%
                c("cirrhosis_deaths_averted", "hcc_deaths_averted"))) +
  geom_line(aes(x=time, y = value/0.5, group = sim), col = "grey") +
  stat_summary(aes(x=time, y = value), fun="median", geom = "line", col = "red")+
  facet_wrap(outcome~age,ncol = 3, scales = "free_y",
             labeller = labeller(outcome = facet_labels_deaths)) +
  ylab("HBV-related deaths averted per year") +
  geom_hline(yintercept=0) +
  labs(title="HBV-related deaths averted by treatment in different age groups in women") +
  theme_classic()

# Cumulative HBV/HCC deaths by age by 2100 to show shift in pattern of HCC deaths by age in the treated cohort ----
# What would happen in those who need treatment but would not get it?
# Following the screened and treated cohort shows this, but gives the total in both
# Need to substract from this the deaths that occur in those screened, to get the deaths in
# treated untreated only.
ages <- seq(0,99.5,0.5)

# X1 = screened+treated HCC deaths in out1
# X2 = screened HCC deaths in out3
# X3 = X1-X2
# X4 = treated HCC deaths in out3
# x1 = deaths that would occur in the screened and treated people in total
# if those who need treatment don't get it
# x2 = deaths in those screened but not treated if those who need treatment receive it
# x3 = deaths that would occur in those who need treatment but don't receive it
# x4 = deaths that occur with treatment in those treated

# HCC CASES

# Treat ages 15-30:
a4_x1 <- a4_out1$cum_cohort_hcc_cases_male_by_age_2100+
  a4_out1$cum_cohort_hcc_cases_female_by_age_2100
a4_x2 <- a4_out3$cum_screened_hcc_cases_male_by_age_2100 +
  a4_out3$cum_screened_hcc_cases_female_by_age_2100
a4_x3 <- a4_x1-a4_x2
a4_x4 <- a4_out3$cum_treated_hcc_cases_male_by_age_2100+
  a4_out3$cum_treated_hcc_cases_female_by_age_2100
a4_x3$sim <- rownames(a4_x3)
a4_x4$sim <- rownames(a4_x4)
# Treat ages 30-45
a5_x1 <- a5_out1$cum_cohort_hcc_cases_male_by_age_2100+
  a5_out1$cum_cohort_hcc_cases_female_by_age_2100
a5_x2 <- a5_out3$cum_screened_hcc_cases_male_by_age_2100 +
  a5_out3$cum_screened_hcc_cases_female_by_age_2100
a5_x3 <- a5_x1-a5_x2
a5_x4 <- a5_out3$cum_treated_hcc_cases_male_by_age_2100+
  a5_out3$cum_treated_hcc_cases_female_by_age_2100
a5_x3$sim <- rownames(a5_x3)
a5_x4$sim <- rownames(a5_x4)
# Treat ages 45-65
a2_x1 <- a2_out1$cum_cohort_hcc_cases_male_by_age_2100+
  a2_out1$cum_cohort_hcc_cases_female_by_age_2100
a2_x2 <- a2_out3$cum_screened_hcc_cases_male_by_age_2100 +
  a2_out3$cum_screened_hcc_cases_female_by_age_2100
a2_x3 <- a2_x1-a2_x2
a2_x4 <- a2_out3$cum_treated_hcc_cases_male_by_age_2100+
  a2_out3$cum_treated_hcc_cases_female_by_age_2100
a2_x3$sim <- rownames(a2_x3)
a2_x4$sim <- rownames(a2_x4)


df <- rbind(
  data.frame(type="No treatment", treated_age="15-30",
             gather(a4_x3,key="age", value = "value",-sim)),
  data.frame(type="With treatment", treated_age="15-30",
             gather(a4_x4,key="age", value = "value",-sim)),
  data.frame(type="No treatment", treated_age="30-45",
             gather(a5_x3,key="age", value = "value",-sim)),
  data.frame(type="With treatment", treated_age="30-45",
             gather(a5_x4,key="age", value = "value",-sim)),
  data.frame(type="No treatment", treated_age="45-65",
             gather(a2_x3,key="age", value = "value",-sim)),
  data.frame(type="With treatment", treated_age="45-65",
             gather(a2_x4,key="age", value = "value",-sim))
)
df$age <- rep(ages, each = 183)

# Divide deaths by 0.5 to get value in each 1-year age group
ggplot(df) +
  facet_wrap(~treated_age, ncol = 3) +
  geom_line(data=subset(df, type == "No treatment"),
            aes(x=age, y = value/0.5, group = sim), col = "grey70") +
  geom_line(data=subset(df, type == "With treatment"),
            aes(x=age, y = value/0.5, group = sim), col = "grey90") +
  stat_summary(aes(x=age, y = value/0.5, colour=type), fun="median", geom = "line", size = 2)+
  #facet_wrap(outcome~age,ncol = 3, scales = "free_y") +
  ylab("Cumulative HBV-related HCC cases 2020-2100") +
  theme_bw()

# Ideally make separate plots and start them at age at enrollment into cohort (15/30/45)
# or have no data at 0-entry and a vertical line
age_p1 <-ggplot(subset(df, treated_age == "15-30")) +
  stat_summary(aes(x=age, y = value/0.5, fill=type, colour = type),
               fun.min= function(x) quantile(x,0.025),
               fun.max= function(x) quantile(x,0.975),
               geom = "ribbon", alpha = 0.05, lty ="dashed")+
  stat_summary(aes(x=age, y = value/0.5, colour=type), fun="median", geom = "line", size = 1)+
  ylab("Cumulative HBV-related HCC cases 2020-2100") +
  xlab("Age at HCC onset (years)") +
  labs(title = "15-30 year old cohort") +
  geom_rect(aes(xmin=0, xmax=15, ymin=-Inf, ymax=Inf), fill = "grey") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction="vertical",
        legend.text = element_text(color = "white"),   # use these lines to make legend disappear
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white")) +
  scale_color_discrete(guide = guide_legend(override.aes = list(color = "white")))

age_p2 <-ggplot(subset(df, treated_age == "30-45")) +
  stat_summary(aes(x=age, y = value/0.5, fill=type, colour = type),
               fun.min= function(x) quantile(x,0.025),
               fun.max= function(x) quantile(x,0.975),
               geom = "ribbon", alpha = 0.05, lty ="dashed")+
  stat_summary(aes(x=age, y = value/0.5, colour=type), fun="median", geom = "line", size = 1)+
  ylab("Cumulative HBV-related HCC cases 2020-2100") +
  xlab("Age at HCC onset (years)") +
  labs(fill = "Scenario:", colour = "Scenario:", title = "30-45 year old cohort") +
  geom_rect(aes(xmin=0, xmax=30, ymin=-Inf, ymax=Inf), fill = "grey") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.direction="vertical")

age_p3 <-ggplot(subset(df, treated_age == "45-65")) +
  stat_summary(aes(x=age, y = value/0.5, fill=type, colour = type),
               fun.min= function(x) quantile(x,0.025),
               fun.max= function(x) quantile(x,0.975),
               geom = "ribbon", alpha = 0.05, lty ="dashed")+
  stat_summary(aes(x=age, y = value/0.5, colour=type), fun="median", geom = "line", size = 1)+
  ylab("Cumulative HBV-related HCC cases 2020-2100") +
  xlab("Age at HCC onset (years)") +
  labs(title = "45-65 year old cohort") +
  geom_rect(aes(xmin=0, xmax=45, ymin=-Inf, ymax=Inf), fill = "grey") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction="vertical",
        legend.text = element_text(color = "white"),   # use these lines to make legend disappear
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white")) +
  scale_color_discrete(guide = guide_legend(override.aes = list(color = "white")))

grid.arrange(age_p1,age_p2,age_p3,ncol = 3,
             top ="Effect of antiviral therapy on HCC incidence in treatment-eligible HBV carriers")




# Population effects: plot with resource utilisation ----
# Number of deaths averted in each age group
deaths_averted_by_age <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a4_out3_pop,
                                                  a5_out3_pop,
                                                  a2_out3_pop),
                          counterfactual_label = "treatment programme without monitoring")

deaths_averted_by_age_proportion <- rbind(
  data.frame(age_group = "15-30",
             filter(deaths_averted_by_age, scenario == "a4_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age, scenario == "a5_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(deaths_averted_by_age, scenario == "a2_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value))
)

deaths_averted_by_age <- rbind(
  data.frame(age_group = "15-30",
             filter(deaths_averted_by_age, scenario == "a4_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age, scenario == "a5_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(deaths_averted_by_age, scenario == "a2_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value))
)
deaths_averted_by_age$sim <- gsub("[^0-9]", "", deaths_averted_by_age$sim)

# Number of DALYS averted in each age group
dalys_averted_by_age <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a4_out3_pop,
                                                  a5_out3_pop,
                                                  a2_out3_pop),
                          outcome_to_avert = "dalys",
                          counterfactual_label = "treatment programme without monitoring")

dalys_averted_by_age <- rbind(
  data.frame(age_group = "15-30",
             filter(dalys_averted_by_age, scenario == "a4_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(dalys_averted_by_age, scenario == "a5_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(dalys_averted_by_age, scenario == "a2_screen_2020_monit_0" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value))
)
dalys_averted_by_age$sim <- gsub("[^0-9]", "", dalys_averted_by_age$sim)
colnames(dalys_averted_by_age)[colnames(dalys_averted_by_age)=="value"] <- "dalys_averted"

# Add carriers by age in 2020
# Calculate number of carriers by age group
carriers_female <- lapply(out2_carriers, "[[", "carriers_female")
carriers_male <- lapply(out2_carriers, "[[", "carriers_male")
total_carriers_by_age_2020 <- do.call(rbind.data.frame, (lapply(carriers_female, function(x) x[which(out2_carriers[[1]]$time==2020),])))+
  do.call(rbind.data.frame, (lapply(carriers_male, function(x) x[which(out2_carriers[[1]]$time==2020),])))

# Add treatment eligible carriers by age in 2020
carriers_eligible_by_age_2020 <-
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$ir_female,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$ir_male,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$enchb_female,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$enchb_male,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$cc_female,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$cc_male,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$dcc_female,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$dcc_male,
                 function(x) x[which(out2_comps_by_age$time==2020),]))

carriers_it_by_age_2020 <-
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$it_female,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$it_male,
                 function(x) x[which(out2_comps_by_age$time==2020),]))

# Combine into dataframes
carriers_by_age_group_2020 <-
  rbind(data.frame(age_group = "15-30",
                   sim = rownames(total_carriers_by_age_2020[,which(ages==15): which(ages==30-da)]),
                   carriers = rowSums(total_carriers_by_age_2020[,which(ages==15): which(ages==30-da)])),
        data.frame(age_group = "30-45",
                   sim = rownames(total_carriers_by_age_2020[,which(ages==30): which(ages==45-da)]),
                   carriers = rowSums(total_carriers_by_age_2020[,which(ages==30): which(ages==45-da)])),
        data.frame(age_group = "45-65",
                   sim = rownames(total_carriers_by_age_2020[,which(ages==45): which(ages==65-da)]),
                   carriers = rowSums(total_carriers_by_age_2020[,which(ages==45): which(ages==65-da)]))
  )

carriers_eligible_by_age_group_2020 <-
  rbind(data.frame(age_group = "15-30",
                   sim = rownames(carriers_eligible_by_age_2020[,which(ages==15): which(ages==30-da)]),
                   carriers_eligible = rowSums(carriers_eligible_by_age_2020[,which(ages==15): which(ages==30-da)])),
        data.frame(age_group = "30-45",
                   sim = rownames(carriers_eligible_by_age_2020[,which(ages==30): which(ages==45-da)]),
                   carriers_eligible = rowSums(carriers_eligible_by_age_2020[,which(ages==30): which(ages==45-da)])+
                     rowSums(carriers_it_by_age_2020[,which(ages==30): which(ages==45-da)])),
        data.frame(age_group = "45-65",
                   sim = rownames(carriers_eligible_by_age_2020[,which(ages==45): which(ages==65-da)]),
                   carriers_eligible = rowSums(carriers_eligible_by_age_2020[,which(ages==45): which(ages==65-da)])+
                     rowSums(carriers_it_by_age_2020[,which(ages==45): which(ages==65-da)]))
  )

# Total interactions by age and type (counting treatment initiation only)
interactions_by_age <- rbind(
  data.frame(age_group = "15-30",
             sim = names(a4_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a4_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a4_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a4_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a4_out3_pop$py_on_treatment[[16]])),
  data.frame(age_group = "30-45",
             sim = names(a5_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a5_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a5_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a5_out3_pop$py_on_treatment[[16]])),
  data.frame(age_group = "45-65",
             sim = names(a2_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a2_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a2_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a2_out3_pop$py_on_treatment[[16]]))
)


interactions_by_age <- interactions_by_age %>%
  mutate(total_interactions = hbsag_tests + clinical_assessments + treatment_initiations)
interactions_by_age <- gather(interactions_by_age, key = "interaction_type", value = "value",
                              -age_group, - sim)

# Merge deaths averted by age with carrier population by age and with DALYs averted
outcomes_by_age <- left_join(left_join(left_join(deaths_averted_by_age, carriers_by_age_group_2020,
                                       by = c("age_group", "sim")),
                             dalys_averted_by_age, by = c("age_group","scenario", "sim")),
                             carriers_eligible_by_age_group_2020, by = c("age_group", "sim")) %>%
  mutate(deaths_averted_per_1000_carrier = value/carriers*1000,
         dalys_averted_per_1000_carrier = dalys_averted/carriers*1000,
         deaths_averted_per_eligible_carrier = value/carriers_eligible,
         dalys_averted_per_eligible_carrier = dalys_averted/carriers_eligible)
colnames(outcomes_by_age)[colnames(outcomes_by_age)=="value"] <- "deaths_averted"
outcomes_by_age <- gather(outcomes_by_age, key = "outcome", value = "value",
                          -age_group, -scenario, -sim)

# New facet label names
outcome_labels <- c("DALYs", "DALYs\nper 1,000 carriers", "HBV deaths",
                    "HBV deaths\nper 1,000 carriers")
names(outcome_labels) <- c("dalys_averted", "dalys_averted_per_1000_carrier",
                      "deaths_averted", "deaths_averted_per_1000_carrier")

outcomes_by_age$outcome <- factor(outcomes_by_age$outcome, levels =
                                    c("deaths_averted","deaths_averted_per_1000_carrier",
                                      "dalys_averted", "dalys_averted_per_1000_carrier",
                                      "deaths_averted_per_eligible_carrier", "dalys_averted_per_eligible_carrier",
                                      "carriers", "carriers_eligible"))

# Differences in deaths averted by age group are mostly driven by differences in the
# carrier population, which is highest in the 30-45 group.
# Deaths averted per carrier is very similar across age groups.

# Plot of outcomes absolute and per carrier
ggplot(outcomes_by_age[outcomes_by_age$outcome %in% c("deaths_averted",
                                                      "deaths_averted_per_1000_carrier",
                                                      "dalys_averted",
                                                      "dalys_averted_per_1000_carrier"),]) +
  stat_summary(aes(x=age_group, y= value, group=outcome, fill = outcome),
               fun="median", geom="bar", position = "dodge2", width = 0.95, colour="black")+
  stat_summary(aes(x=age_group, y= value, group=outcome),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", position = "dodge2", width = 0.15)+
  scale_fill_manual(values=c("dalys_averted" = "#A180A9",
                             "dalys_averted_per_1000_carrier" = "#A180A9",
                             "deaths_averted" = "#1F968BFF",
                             "deaths_averted_per_1000_carrier" = "#1F968BFF")) +
  guides(fill=FALSE) +
  xlab("Screened age group (years)")+
  ylab("Number averted") +
  facet_wrap(~outcome, scales="free_y",
             labeller = labeller(outcome = outcome_labels)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))
# Note that outcomes here are for the AGE GROUP TARGETED IN SCREEN (NOT the outcome in that age group)
# DALYS averted per carrier are more similar across age groups than per eligible carrier.


total_interactions_errorbar <- subset(interactions_by_age, interaction_type == "total_interactions") %>%
  group_by(age_group) %>%
  summarise(lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

interactions_median <- subset(interactions_by_age) %>%
  group_by(age_group, interaction_type) %>%
  summarise(value= median(value))

library(viridis)
options(scipen=1000000)
#median_carriers <- subset(outcomes_by_age, outcome=="carriers") %>% group_by(age_group) %>%
#  summarise(median = median(value))

ggplot() +
  geom_col(data=subset(interactions_median, !(interaction_type %in% c("total_interactions", "py_on_treatment"))),
           aes(x=age_group, y= value/1000, fill = reorder(interaction_type, -value)),
           width = 0.95, colour = "black") +
  geom_errorbar(data=total_interactions_errorbar,
                aes(x = age_group, ymin = lower/1000, ymax = upper/1000), width = 0.15) +
  labs(fill = "Resource utilisation") +
  scale_fill_viridis_d(direction=-1,
                       labels = c("hbsag_tests" = "Serological tests",
                                  "clinical_assessments" = "Assessments",
                                  "treatment_initiations" = "Treatment")) +
  #scale_linetype_manual(values = c("HBV carriers" = "dashed")) +
  guides(linetype=guide_legend(title=NULL),
         fill=guide_legend(title=NULL)) +
  theme_classic() +
  theme(legend.position=c(.78,.8)) +
  ylab("Resources utilised\n(thousands)") +
  xlab("Screened age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Same plot with costs:
interactions_by_age$cost <- NA
interactions_by_age$cost[interactions_by_age$interaction_type=="hbsag_tests"] <-
  interactions_by_age$value[interactions_by_age$interaction_type=="hbsag_tests"]*8.3
interactions_by_age$cost[interactions_by_age$interaction_type=="clinical_assessments"] <-
  interactions_by_age$value[interactions_by_age$interaction_type=="clinical_assessments"]*33
interactions_by_age$cost[interactions_by_age$interaction_type=="py_on_treatment"] <-
  interactions_by_age$value[interactions_by_age$interaction_type=="py_on_treatment"]*66.5

interactions_cost_median <- subset(interactions_by_age) %>%
  group_by(age_group, interaction_type) %>%
  summarise(cost= median(cost))

# Need to change order of colour here to match:
ggplot(subset(interactions_cost_median, !(interaction_type %in% c("total_interactions", "treatment_initiations")))) +
  geom_col(aes(x=age_group, y= cost/1000000, fill = interaction_type),
           width = 0.95, colour = "black") +
#  geom_errorbar(data=total_interactions_errorbar,
#                aes(x = age_group, ymin = lower/1000, ymax = upper/1000),width = 0.25) +
  labs(fill = "Resource utilisation") +
  scale_fill_viridis_d(direction=-1,
                       labels = c("hbsag_tests" = "Serological tests",
                                  "clinical_assessments" = "Assessments",
                                  "py_on_treatment" = "Treatment")) +
  guides(linetype=guide_legend(title=NULL),
         fill=guide_legend(title=NULL)) +
  theme_classic() +
  theme(legend.position=c(.78,.8)) +
  ylab("Cost (millions)") +
  xlab("Screened age group")+
  theme(axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13))

# Population effects: Treatment effect over time ----

hbv_deaths_over_time <-
  rbind(gather(out2$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario),
    gather(a1_out3_pop$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario),
        gather(a1_out6_pop$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario))

hbv_deaths_over_time <- hbv_deaths_over_time %>%
  group_by(time, scenario) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value,0.975))

# Plot options:
ggplot() +
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "status_quo"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5),
              linetype = "dashed",fill = "grey80",colour="black",  alpha = 0.5) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "status_quo"),
            aes(x=time, y = median/0.5), colour="black", fill = "grey",
            linetype = "dashed") +
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "screen_2020_monit_0" &
                            time >=2019.5),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5),
              linetype = "solid",fill = "grey80",colour="black",  alpha = 0.5) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "screen_2020_monit_0" &
                          time >=2019.5),
            aes(x=time, y = median/0.5), colour="black", fill = "grey",
            linetype = "solid") +
  geom_vline(xintercept=2019.5) +
  xlim(1990,2090)+
  theme_classic()

ggplot() +
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "status_quo"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5),
              linetype = "dashed",fill = "#92C5DE",colour="#2166AC",  alpha = 0.5) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "status_quo"),
            aes(x=time, y = median/0.5), colour="#2166AC", fill = "grey",
            linetype = "dashed") +
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "screen_2020_monit_0" &
                            time >=2019.5),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5),
              linetype = "solid",fill = "#92C5DE",colour="#2166AC",  alpha = 0.5) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "screen_2020_monit_0" &
                          time >=2019.5),
            aes(x=time, y = median/0.5), colour="#2166AC", fill = "grey",
            linetype = "solid") +
  geom_vline(xintercept=2019.5) +
  xlim(1990,2090)+
  theme_classic()

ggplot() +
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "status_quo"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5,fill = scenario, colour=scenario),
              linetype = "solid", alpha = 0.5) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "status_quo"),
            aes(x=time, y = median/0.5), colour="#B2182B",
            linetype = "solid", size=0.75) +
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "screen_2020_monit_0"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, fill = scenario, colour=scenario),
              linetype = "solid",alpha = 0.8) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "screen_2020_monit_0"),
            aes(x=time, y = median/0.5), colour="#2166AC",
            linetype = "solid", size=0.75) +
  scale_fill_manual("Scenario",
                    labels = c("status_quo" = "Infant vaccination only",
                               "screen_2020_monit_0" = "Screening & treatment in 2020\n(no monitoring)"),
                    values=c("status_quo" = "#D6604D",
                             "screen_2020_monit_0" = "#92C5DE")) +
  scale_colour_manual("Scenario",
                      labels = c("status_quo" = "Infant vaccination only",
                                 "screen_2020_monit_0" = "Screening & treatment in 2020\n(no monitoring)"),
                      values=c("status_quo" = "#B2182B",
                             "screen_2020_monit_0" = "#2166AC")) +
 # guides(fill=FALSE, colour = FALSE) +
  geom_vline(xintercept=2019.5, linetype="dashed") +
  ylab("HBV-related deaths per year") +
  scale_x_continuous("Year",
                     breaks=c(1990, 2010, 2019.5, 2040, 2070, 2090),
                     labels=c(1990, 2010, 2020, 2040, 2070, 2090),
                     limits=c(1990,2090)) +
  theme_classic() +
  theme(legend.position=c(.78,.8),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13))

