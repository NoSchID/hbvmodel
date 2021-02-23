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

## Different output:

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

## Standard output:

# Normal population simulations of effect of screening by age
out2 <- readRDS(paste0(out_path_monit, "out2_status_quo_301120.rds"))
out2 <- out2[[1]]
# Carriers by age in sq scenario:
out2_carriers <- readRDS(paste0(out_path, "out_sq_carriers.rds"))
out2_comps_by_age <- readRDS(paste0(out_path, "out_sq_compartment_prevalence_by_age.rds"))
out2_disease_outcomes <- readRDS(paste0(out_path, "out_sq_disease_outcomes_by_age.rds"))

da <- 0.5
ages <- seq(0,99.5,0.5)

# Simulations by age group screened, no monitoring:
a2_out3_pop <- readRDS(paste0(out_path_monit, "a2_it_out3_screen_2020_monit_0_180121.rds"))
a2_out3_pop <- a2_out3_pop[[1]]
a4_out3_pop <- readRDS(paste0(out_path_monit, "a4_it_out3_screen_2020_monit_0_190121.rds"))
a4_out3_pop <- a4_out3_pop[[1]]
a5_out3_pop <- readRDS(paste0(out_path_monit, "a5_it_out3_screen_2020_monit_0_190121.rds"))
a5_out3_pop <- a5_out3_pop[[1]]
# Simulations by age group screened, with optimal monitoring:
a2_monit_sim7_pop <- readRDS(paste0(out_path_monit, "a2_it_screen_2020_monit_out7_190221.rds"))
a2_monit_sim7_pop <- a2_monit_sim7_pop[[1]]
a4_monit_sim7_pop <- readRDS(paste0(out_path_monit, "a4_it_screen_2020_monit_out7_190221.rds"))
a4_monit_sim7_pop <- a4_monit_sim7_pop[[1]]
a5_monit_sim7_pop <- readRDS(paste0(out_path_monit, "a5_it_screen_2020_monit_out7_190221.rds"))
a5_monit_sim7_pop <- a5_monit_sim7_pop[[1]]

# With/without monitoring
a1_out3_pop <- readRDS(paste0(out_path_monit, "a1_it_out3_screen_2020_monit_0_180121.rds"))
a1_out3_pop <- a1_out3_pop[[1]]   # No monitoring

a1_out6_pop <- readRDS(paste0(out_path_monit, "a1_it_out6_screen_2020_monit_1_130121.rds"))
a1_out6_pop <- a1_out6_pop[[1]]  # 1 year

# With monitoring in separate age groups
# With/without monitoring
out3_it <- readRDS(paste0(out_path_monit, "a1_it_out3_screen_2020_monit_0_180121.rds"))
out3_it <- out3_it[[1]]   # No monitoring

monit_out10 <- readRDS(paste0(out_path_monit, "a1_it_monit_out10_200121.rds"))
monit_out10 <- monit_out10[[1]]

monit_out6 <- readRDS(paste0(out_path_monit, "a1_it_monit_out6_200121.rds"))
monit_out6 <- monit_out6[[1]]

monit_out8 <- readRDS(paste0(out_path_monit, "a1_it_monit_out8_210121.rds"))
monit_out8 <- monit_out8[[1]]

# Different monitoring frequencies across all ages
# out3_it above
# a1_out6_pop
out4b_it <- readRDS(paste0(out_path_monit, "a1_it_out4b_screen_2020_monit_20_230221.rds"))
out4b_it <- out4b_it[[1]]   # 20 years
out4_it <- readRDS(paste0(out_path_monit, "a1_it_out4_screen_2020_monit_10_140121.rds"))
out4_it <- out4_it[[1]]   # 10 years
out5_it <- readRDS(paste0(out_path_monit, "a1_it_out5_screen_2020_monit_5_161220.rds"))
out5_it <- out5_it[[1]]   # 5 years
out6a_it <- readRDS(paste0(out_path_monit, "a1_it_out6a_screen_2020_monit_2_161220.rds"))
out6a_it <- out6a_it[[1]]  # 2 years

## 1) Simulate the effect of treatment in the screened+treated cohort ----
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
age_p1 <- ggplot(subset(df, treated_age == "15-30")) +
  stat_summary(aes(x=age, y = value/0.5, fill=type, colour = type),
               fun.min= function(x) quantile(x,0.025),
               fun.max= function(x) quantile(x,0.975),
               geom = "ribbon", alpha = 0.05, lty ="dashed")+
  stat_summary(aes(x=age, y = value/0.5, colour=type), fun="median", geom = "line", size = 1)+
  ylab("Cumulative HBV-related HCC cases 2020-2100") +
  xlab("Age at HCC onset (years)") +
  labs(title = "15-30 year old cohort") +
  geom_rect(aes(xmin=0, xmax=15, ymin=-Inf, ymax=Inf), fill = "grey") +
  scale_fill_manual(values= c("No treatment"= "#A180A9",
                              "With treatment" ="#1F968BFF")) +
  scale_color_manual(values= c("No treatment"= "#A180A9",
                               "With treatment" ="#1F968BFF"),
                     guide = guide_legend(override.aes = list(color = "white"))) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction="vertical",
        legend.text = element_text(color = "white"),   # use these lines to make legend disappear
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"))

age_p2 <-ggplot(subset(df, treated_age == "30-45")) +
  stat_summary(aes(x=age, y = value/0.5, fill=type, colour = type),
               fun.min= function(x) quantile(x,0.025),
               fun.max= function(x) quantile(x,0.975),
               geom = "ribbon", alpha = 0.05, lty ="dashed")+
  stat_summary(aes(x=age, y = value/0.5, colour=type), fun="median", geom = "line", size = 1)+
  ylab("Cumulative HBV-related HCC cases 2020-2100") +
  xlab("Age at HCC onset (years)") +
  scale_fill_manual(values= c("No treatment"= "#A180A9",
                              "With treatment" ="#1F968BFF")) +
  scale_color_manual(values= c("No treatment"= "#A180A9",
                               "With treatment" ="#1F968BFF")) +
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
  scale_fill_manual(values= c("No treatment"= "#A180A9",
                              "With treatment" ="#1F968BFF")) +
  scale_color_manual(values= c("No treatment"= "#A180A9",
                               "With treatment" ="#1F968BFF"),
                     guide = guide_legend(override.aes = list(color = "white"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction="vertical",
        legend.text = element_text(color = "white"),   # use these lines to make legend disappear
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"))

grid.arrange(age_p1,age_p2,age_p3,ncol = 3,
             top ="Effect of antiviral therapy on HCC incidence in treatment-eligible HBV carriers")




# HBV deaths absolute number and rate in SQ scenario ----

hbv_deaths_2020 <-
  (do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
       function(x) x[which(out2_disease_outcomes$time==2021),]))-
  do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                          function(x) x[which(out2_disease_outcomes$time==2020),])))+
  (do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_female,
                          function(x) x[which(out2_disease_outcomes$time==2021),]))-
  do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_female,
                          function(x) x[which(out2_disease_outcomes$time==2020),])))

hbv_deaths_2020_group <- rbind(
  data.frame(age = "<15",
        value = apply(hbv_deaths_2020[,which(ages==0):which(ages==14.5)],1,sum)),
  data.frame(age = "15-30",
        value = apply(hbv_deaths_2020[,which(ages==15):which(ages==29.5)],1,sum)),
  data.frame(age = "30-45",
        value = apply(hbv_deaths_2020[,which(ages==30):which(ages==44.5)],1,sum)),
  data.frame(age = "45-65",
        value = apply(hbv_deaths_2020[,which(ages==45):which(ages==64.5)],1,sum)),
  data.frame(age = "65+",
        value = apply(hbv_deaths_2020[,which(ages==65):which(ages==99.5)],1,sum))
)

hcc_2020 <-
  (do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_cases_male,
                           function(x) x[which(out2_disease_outcomes$time==2021),]))-
     do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_cases_male,
                             function(x) x[which(out2_disease_outcomes$time==2020),])))+
  (do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_cases_female,
                           function(x) x[which(out2_disease_outcomes$time==2021),]))-
     do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_cases_female,
                             function(x) x[which(out2_disease_outcomes$time==2020),])))

hcc_2020_group <- rbind(
  data.frame(age = "<15",
             value = apply(hcc_2020[,which(ages==0):which(ages==14.5)],1,sum)),
  data.frame(age = "15-30",
             value = apply(hcc_2020[,which(ages==15):which(ages==29.5)],1,sum)),
  data.frame(age = "30-45",
             value = apply(hcc_2020[,which(ages==30):which(ages==44.5)],1,sum)),
  data.frame(age = "45-65",
             value = apply(hcc_2020[,which(ages==45):which(ages==64.5)],1,sum)),
  data.frame(age = "65+",
             value = apply(hcc_2020[,which(ages==65):which(ages==99.5)],1,sum))
)


hbv_deaths_2020 <-
  gather(hbv_deaths_2020,
         key="age", value = "value")
hbv_deaths_2020$sim <- rep(1:183, 200)
hbv_deaths_2020$age <- rep(seq(0,99.5,0.5), each=183)

carriers_male <- lapply(out2_carriers, "[[", "carriers_male")
carriers_female <- lapply(out2_carriers, "[[", "carriers_female")

carriers_2020 <-
  do.call(rbind.data.frame, (lapply(carriers_male, function(x)
    x[which(out2_carriers[[1]]$time==2020.5),])))+
    do.call(rbind.data.frame, (lapply(carriers_female, function(x)
      x[which(out2_carriers[[1]]$time==2020.5),])))

carriers_2020_group <- rbind(
  data.frame(age = "<15",
             value = apply(carriers_2020[,which(ages==0):which(ages==14.5)],1,sum)),
  data.frame(age = "15-30",
             value = apply(carriers_2020[,which(ages==15):which(ages==29.5)],1,sum)),
  data.frame(age = "30-45",
             value = apply(carriers_2020[,which(ages==30):which(ages==44.5)],1,sum)),
  data.frame(age = "45-65",
             value = apply(carriers_2020[,which(ages==45):which(ages==64.5)],1,sum)),
  data.frame(age = "65+",
             value = apply(carriers_2020[,which(ages==65):which(ages==99.5)],1,sum))
)

carriers_2020 <- gather(carriers_2020, key="age", value = "value")
hbv_deaths_2020$carriers <- carriers_2020$value
hbv_deaths_2020$rate <- hbv_deaths_2020$value/hbv_deaths_2020$carriers

hbv_deaths_2020_group$carriers <- carriers_2020_group$value
hbv_deaths_2020_group$rate <- hbv_deaths_2020_group$value/hbv_deaths_2020_group$carriers

p1 <- ggplot(hbv_deaths_2020_group) +
  stat_summary(aes(x=age, y = value), fun = "median", geom = "col") +
  stat_summary(aes(x=age, y = value), fun.min = {function(x) quantile(x, 0.025)},
               fun.max = {function(x) quantile(x, 0.975)},
               geom = "errorbar", width = 0.15) +
  ylab("HBV-related deaths in 2020") +
  xlab("Age (years)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

p2 <- ggplot(hbv_deaths_2020) +
  stat_summary(aes(x=age, y = rate*10000), fun = "median", geom = "line") +
  stat_summary(aes(x=age, y = rate*10000), fun.min = {function(x) quantile(x, 0.025)},
               fun.max = {function(x) quantile(x, 0.975)},
               geom = "ribbon", alpha=0.1) +
  ylab("HBV-related mortality rate\nper 10,000 carriers in 2020") +
  xlab("Age (years)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

grid.arrange(p1, p2, ncol = 2)

ggplot(hcc_2020_group) +
  stat_summary(aes(x=age, y = value), fun = "median", geom = "col") +
  stat_summary(aes(x=age, y = value), fun.min = {function(x) quantile(x, 0.025)},
               fun.max = {function(x) quantile(x, 0.975)},
               geom = "errorbar", width = 0.15) +
  ylab("Incident HCC cases in 2020") +
  xlab("Age (years)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

## Population effects of screening and monitoring by age and over time ----
# Population effects: plot with resource utilisation ----
# Number of deaths averted in each age group
deaths_averted_by_age <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a4_out3_pop,
                                                  a5_out3_pop,
                                                  a2_out3_pop,
                                                  a4_monit_sim7_pop,
                                                  a5_monit_sim7_pop,
                                                  a2_monit_sim7_pop),
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
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "15-30",
             filter(deaths_averted_by_age, scenario == "a4_screen_2020_monit_sim7" &
                      by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age, scenario == "a5_screen_2020_monit_sim7" &
                      by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(deaths_averted_by_age, scenario == "a2_screen_2020_monit_sim7" &
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
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "15-30",
             filter(deaths_averted_by_age, scenario == "a4_screen_2020_monit_sim7" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age, scenario == "a5_screen_2020_monit_sim7" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(deaths_averted_by_age, scenario == "a2_screen_2020_monit_sim7" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value))
)
deaths_averted_by_age$sim <- gsub("[^0-9]", "", deaths_averted_by_age$sim)

# Number of DALYS averted in each age group
dalys_averted_by_age <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a4_out3_pop,
                                                  a5_out3_pop,
                                                  a2_out3_pop,
                                                  a4_monit_sim7_pop,
                                                  a5_monit_sim7_pop,
                                                  a2_monit_sim7_pop),
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
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "15-30",
             filter(dalys_averted_by_age, scenario == "a4_screen_2020_monit_sim7" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(dalys_averted_by_age, scenario == "a5_screen_2020_monit_sim7" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(dalys_averted_by_age, scenario == "a2_screen_2020_monit_sim7" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value))
)
dalys_averted_by_age$sim <- gsub("[^0-9]", "", dalys_averted_by_age$sim)
colnames(dalys_averted_by_age)[colnames(dalys_averted_by_age)=="value"] <- "dalys_averted"

# Add carriers by age in 2020
# Calculate number of carriers by age group
carriers_female <- lapply(out2_carriers, "[[", "carriers_female")
carriers_male <- lapply(out2_carriers, "[[", "carriers_male")
pop_female <- lapply(out2_carriers, "[[", "pop_female")
pop_male <- lapply(out2_carriers, "[[", "pop_male")
total_carriers_by_age_2020 <- do.call(rbind.data.frame, (lapply(carriers_female, function(x) x[which(out2_carriers[[1]]$time==2020),])))+
  do.call(rbind.data.frame, (lapply(carriers_male, function(x) x[which(out2_carriers[[1]]$time==2020),])))
prevalence_by_age_2020 <- (do.call(rbind.data.frame, (lapply(carriers_female, function(x) x[which(out2_carriers[[1]]$time==2020),])))+
  do.call(rbind.data.frame, (lapply(carriers_male, function(x) x[which(out2_carriers[[1]]$time==2020),]))))/
  (do.call(rbind.data.frame, (lapply(pop_female, function(x) x[which(out2_carriers[[1]]$time==2020),])))+
  do.call(rbind.data.frame, (lapply(pop_male, function(x) x[which(out2_carriers[[1]]$time==2020),]))))

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

prevalence_by_age_group_2020 <-
  rbind(data.frame(age_group = "15-30",
                   sim = rownames(prevalence_by_age_2020[,which(ages==15): which(ages==30-da)]),
                   carriers = rowSums(prevalence_by_age_2020[,which(ages==15): which(ages==30-da)])),
        data.frame(age_group = "30-45",
                   sim = rownames(prevalence_by_age_2020[,which(ages==30): which(ages==45-da)]),
                   carriers = rowSums(prevalence_by_age_2020[,which(ages==30): which(ages==45-da)])),
        data.frame(age_group = "45-65",
                   sim = rownames(prevalence_by_age_2020[,which(ages==45): which(ages==65-da)]),
                   carriers = rowSums(prevalence_by_age_2020[,which(ages==45): which(ages==65-da)]))
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
             scenario = "a4_screen_2020_monit_0",
             sim = names(a4_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a4_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a4_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a4_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a4_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a4_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a4_out3_pop$py_on_treatment[[16]])),
  data.frame(age_group = "30-45",
             scenario = "a5_screen_2020_monit_0",
             sim = names(a5_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a5_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a5_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a5_out3_pop$py_on_treatment[[16]])),
  data.frame(age_group = "45-65",
             scenario = "a2_screen_2020_monit_0",
             sim = names(a2_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a2_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a2_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a2_out3_pop$py_on_treatment[[16]])),
  data.frame(age_group = "15-30",
             scenario = "a4_screen_2020_monit_sim7",
             sim = names(a4_monit_sim7_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a4_monit_sim7_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a4_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a4_monit_sim7_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a4_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a4_monit_sim7_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a4_monit_sim7_pop$py_on_treatment[[16]])),
  data.frame(age_group = "30-45",
             scenario = "a5_screen_2020_monit_sim7",
             sim = names(a5_monit_sim7_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a5_monit_sim7_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a5_monit_sim7_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a5_monit_sim7_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a5_monit_sim7_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a5_monit_sim7_pop$py_on_treatment[[16]])),
  data.frame(age_group = "45-65",
             scenario = "a2_screen_2020_monit_sim7",
             sim = names(a2_monit_sim7_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a2_monit_sim7_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a2_monit_sim7_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a2_monit_sim7_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a2_monit_sim7_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a2_monit_sim7_pop$py_on_treatment[[16]]))
)

interactions_by_age <- interactions_by_age %>%
  mutate(total_interactions = hbsag_tests + clinical_assessments + monitoring_assessments + treatment_initiations,
         treatment_initiations_per_assessment = treatment_initiations/(clinical_assessments+monitoring_assessments),
         py_on_treatment_per_initiation = py_on_treatment/treatment_initiations)
interactions_by_age <- gather(interactions_by_age, key = "interaction_type", value = "value",
                              -age_group, -scenario, - sim)

interactions_by_age_rel <- subset(interactions_by_age, interaction_type %in%
                                    c("treatment_initiations_per_assessment",
                                      "py_on_treatment_per_initiation"))
interactions_by_age <- subset(interactions_by_age, !(interaction_type %in%
                                    c("treatment_initiations_per_assessment",
                                      "py_on_treatment_per_initiation")))

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
ggplot(subset(outcomes_by_age, scenario %in% c("a2_screen_2020_monit_sim7",
                                               "a4_screen_2020_monit_sim7",
                                               "a5_screen_2020_monit_sim7") &
                                               outcome %in% c("deaths_averted",
                                                      "deaths_averted_per_1000_carrier",
                                                      "dalys_averted",
                                                      "dalys_averted_per_1000_carrier"))) +
  stat_summary(aes(x=age_group, y= value, group=outcome, fill = outcome),
               fun="median", geom="bar", position = "dodge2", width = 0.95, colour="black")+
  # Add median without monitoring:
  stat_summary(data= subset(outcomes_by_age, scenario %in% c("a2_screen_2020_monit_0",
                                                             "a4_screen_2020_monit_0",
                                                             "a5_screen_2020_monit_0") &
                              outcome %in% c("deaths_averted",
                                             "deaths_averted_per_1000_carrier",
                                             "dalys_averted",
                                             "dalys_averted_per_1000_carrier")),
               aes(x=age_group, y= value, group=outcome, fill = outcome),
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

# Only absolutes:
x3 <- ggplot(subset(outcomes_by_age, scenario %in% c("a2_screen_2020_monit_sim7",
                                                     "a4_screen_2020_monit_sim7",
                                                     "a5_screen_2020_monit_sim7") &
                      outcome %in% c("deaths_averted",
                                     "dalys_averted"))) +
  stat_summary(aes(x=age_group, y= value, group=outcome, fill = outcome),
               fun="median", geom="bar", position = "dodge2", width = 0.95, colour="black",
               alpha=0.5)+
  stat_summary(data= subset(outcomes_by_age, scenario %in% c("a2_screen_2020_monit_0",
                                                             "a4_screen_2020_monit_0",
                                                             "a5_screen_2020_monit_0") &
                              outcome %in% c("deaths_averted",
                                             "dalys_averted")),
               aes(x=age_group, y= value, group=outcome, fill = outcome),
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
  facet_wrap(~outcome, ncol=1, scales="free_y",
             labeller = labeller(outcome = outcome_labels)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))


# Plot distribution of carriers
x1<- ggplot(carriers_by_age_group_2020) +
  stat_summary(aes(x=age_group, y= carriers/1000),
               fun="median", geom="bar", position = "dodge2", width = 0.95,
               colour="black", fill = "grey70") +
  stat_summary(aes(x=age_group, y= carriers/1000),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", position = "dodge2", width = 0.15)+
#  scale_fill_manual(values=c("dalys_averted" = "#A180A9",
#                             "dalys_averted_per_1000_carrier" = "#A180A9",
#                             "deaths_averted" = "#1F968BFF",
#                             "deaths_averted_per_1000_carrier" = "#1F968BFF")) +
#  guides(fill=FALSE) +
  xlab("Age group (years)")+
  ylab("HBV carriers in 2020 (thousands)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

# Plot of % prevalence
x2 <- ggplot(prevalence_by_age_group_2020) +
  stat_summary(aes(x=age_group, y= carriers),
               fun="median", geom="bar", position = "dodge2", width = 0.95, colour="black",
               fill = "grey70") +
  stat_summary(aes(x=age_group, y= carriers),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", position = "dodge2", width = 0.15)+
  #  scale_fill_manual(values=c("dalys_averted" = "#A180A9",
  #                             "dalys_averted_per_1000_carrier" = "#A180A9",
  #                             "deaths_averted" = "#1F968BFF",
  #                             "deaths_averted_per_1000_carrier" = "#1F968BFF")) +
  #  guides(fill=FALSE) +
  xlab("Age group (years)")+
  ylab("HBV prevalence in 2020 (%)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

ggplot(gather(prevalence_by_age_2020, key = "age", value = "value") %>%
  mutate(age=rep(seq(0,99.5,0.5), each = 183))) +
  stat_summary(aes(x=age, y= value*100),
               fun="median", geom="line", colour="black") +
  stat_summary(aes(x=age, y= value*100),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="ribbon", alpha = 0.1)  +
  xlab("Age (years)")+
  ylab("HBV prevalence (%)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

# Plot of treatment initiations per assessment
interactions_by_age_rel$with_monitoring <- "No"
interactions_by_age_rel$with_monitoring[interactions_by_age_rel$scenario %in% c("a4_screen_2020_monit_sim7",
                                                                "a5_screen_2020_monit_sim7",
                                                                "a2_screen_2020_monit_sim7")] <- "Yes"

x4 <- ggplot(subset(interactions_by_age_rel, with_monitoring == "No" & interaction_type==
                "treatment_initiations_per_assessment")) +
  geom_boxplot(aes(x=age_group, y = value), fill = "grey70") +
  ylim(0,0.35) +
  xlab("Age group (years)")+
  ylab("Treatment initiations\nper initial assessment") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

x5 <- ggplot(subset(interactions_by_age_rel, with_monitoring == "No" & interaction_type==
                "py_on_treatment_per_initiation")) +
  geom_boxplot(aes(x=age_group, y = value), fill = "grey70") +
  ylim(0,40) +
  xlab("Age group (years)")+
  ylab("Person-years on treatment\nper treatment initiation") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

# Combined plot
grid.arrange(x3, grid.arrange(x1,x2, ncol = 1), grid.arrange(x4,x5, ncol = 1), ncol = 3)

# Distribution of resource utilisation ----
total_interactions_errorbar <- subset(interactions_by_age, interaction_type == "total_interactions") %>%
  group_by(scenario, age_group) %>%
  summarise(lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

interactions_median <- subset(interactions_by_age) %>%
  group_by(scenario, age_group, interaction_type) %>%
  summarise(value= median(value))
interactions_median$interaction_type <- factor(interactions_median$interaction_type,
                                               levels= c("total_interactions","monitoring_assessments",
                                                         "py_on_treatment", "treatment_initiations",
                                                  "clinical_assessments","hbsag_tests"))

library(viridis)
options(scipen=1000000)
#median_carriers <- subset(outcomes_by_age, outcome=="carriers") %>% group_by(age_group) %>%
#  summarise(median = median(value))

x2_1 <- ggplot() +
  geom_col(data=subset(interactions_median, scenario %in% c("a2_screen_2020_monit_sim7",
                                                            "a4_screen_2020_monit_sim7",
                                                            "a5_screen_2020_monit_sim7") &
                         !(interaction_type %in% c("total_interactions", "py_on_treatment"))),
           aes(x=age_group, y= value/1000, fill = interaction_type),
           width = 0.95, colour = "black") +
#geom_errorbar(data=subset(total_interactions_errorbar,scenario %in% c("a2_screen_2020_monit_sim7",
#                                                                        "a4_screen_2020_monit_sim7",
#                                                                        "a5_screen_2020_monit_sim7")),
#                aes(x = age_group, ymin = lower/1000, ymax = upper/1000), width = 0.15) +
  labs(fill = "Resource utilisation") +
  scale_fill_manual(labels = c("hbsag_tests" = "Serological tests",
                                  "clinical_assessments" = "Assessments",
                                  "monitoring_assessments" = "Monitoring",
                                  "treatment_initiations" = "Treatment"),
                    values = c("hbsag_tests" = "#FDE725FF",
                               "clinical_assessments" = "#35B779FF",
                               "monitoring_assessments" = "#31688EFF",
                               "treatment_initiations" = "#440154FF")) +
  #scale_linetype_manual(values = c("HBV carriers" = "dashed")) +
  guides(linetype=guide_legend(title=NULL),
         fill=guide_legend(title=NULL)) +
  theme_classic() +
  theme(legend.position=c(.78,.8)) +
  ylab("Resources utilised (thousands)") +
  xlab("Screened age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Same plot with DISCOUNTED costs:

object_list1 <- list(a4_out3_pop,
                    a4_monit_sim7_pop)
object_list2 <- list(a5_out3_pop,
                     a5_monit_sim7_pop)
object_list3 <- list(a2_out3_pop,
                     a2_monit_sim7_pop)

# Calculate discounted interactions
age_interactions1 <- list()
age_interactions_py_on_treatment1 <- list()
age_interactions2 <- list()
age_interactions_py_on_treatment2 <- list()
age_interactions3 <- list()
age_interactions_py_on_treatment3 <- list()
annual_discounting_rate <- 0.03

for (i in 1:length(object_list1)) {
  age_interactions1[[i]] <-
    cbind(scenario = object_list1[[i]]$cohort_age_at_death$scenario,
          assemble_discounted_interactions_for_monitoring_frequencies(object_list1[[i]],
                                                                      no_monitoring_object = a4_out3_pop))
  age_interactions_py_on_treatment1[[i]] <-
    data.frame(scenario = object_list1[[i]]$cohort_age_at_death$scenario,
               discount_outcome_2020_to_2100(scenario_object=object_list1[[i]],
                                             object_to_subtract=NULL,
                                             outcome="py_on_treatment",
                                             yearly_discount_rate=annual_discounting_rate))
  age_interactions2[[i]] <-
    cbind(scenario = object_list2[[i]]$cohort_age_at_death$scenario,
          assemble_discounted_interactions_for_monitoring_frequencies(object_list2[[i]],
                                                                      no_monitoring_object = a5_out3_pop))
  age_interactions_py_on_treatment2[[i]] <-
    data.frame(scenario = object_list2[[i]]$cohort_age_at_death$scenario,
               discount_outcome_2020_to_2100(scenario_object=object_list2[[i]],
                                             object_to_subtract=NULL,
                                             outcome="py_on_treatment",
                                             yearly_discount_rate=annual_discounting_rate))
  age_interactions3[[i]] <-
    cbind(scenario = object_list3[[i]]$cohort_age_at_death$scenario,
          assemble_discounted_interactions_for_monitoring_frequencies(object_list3[[i]],
                                                                      no_monitoring_object = a2_out3_pop))
  age_interactions_py_on_treatment3[[i]] <-
    data.frame(scenario = object_list3[[i]]$cohort_age_at_death$scenario,
               discount_outcome_2020_to_2100(scenario_object=object_list3[[i]],
                                             object_to_subtract=NULL,
                                             outcome="py_on_treatment",
                                             yearly_discount_rate=annual_discounting_rate))

}
age_interactions <-rbind(do.call("rbind", age_interactions1),
                           do.call("rbind", age_interactions2),
                           do.call("rbind", age_interactions3))
age_interactions_py_on_treatment <- rbind(do.call("rbind", age_interactions_py_on_treatment1),
                                          do.call("rbind", age_interactions_py_on_treatment2),
                                          do.call("rbind", age_interactions_py_on_treatment3))

discounted_interactions <- left_join(age_interactions, age_interactions_py_on_treatment,
                                     by=c("scenario", "sim"))
discounted_interactions$age_group <- "15-30"
discounted_interactions$age_group[discounted_interactions$scenario %in%
                                    c("a5_screen_2020_monit_0","a5_screen_2020_monit_sim7")] <- "30-45"
discounted_interactions$age_group[discounted_interactions$scenario %in%
                                    c("a2_screen_2020_monit_0","a2_screen_2020_monit_sim7")] <- "45-65"
discounted_interactions <- gather(discounted_interactions, key="interaction_type",
                                  value="value", -scenario, -age_group,-sim)

discounted_interactions$cost <- NA
discounted_interactions$cost[discounted_interactions$interaction_type=="hbsag_tests"] <-
  discounted_interactions$value[discounted_interactions$interaction_type=="hbsag_tests"]*8.3
discounted_interactions$cost[discounted_interactions$interaction_type=="clinical_assessments"] <-
  discounted_interactions$value[discounted_interactions$interaction_type=="clinical_assessments"]*33
discounted_interactions$cost[discounted_interactions$interaction_type=="py_on_treatment"] <-
  discounted_interactions$value[discounted_interactions$interaction_type=="py_on_treatment"]*66.5
discounted_interactions$cost[discounted_interactions$interaction_type=="monitoring_assessments"] <-
  discounted_interactions$value[discounted_interactions$interaction_type=="monitoring_assessments"]*25.5

discounted_interactions_cost_median <- subset(discounted_interactions) %>%
  group_by(scenario, age_group, interaction_type) %>%
  summarise(cost= median(cost))
discounted_interactions_cost_median$interaction_type <- factor(discounted_interactions_cost_median$interaction_type,
                                               levels= c("total_interactions","monitoring_assessments",
                                                         "py_on_treatment", "treatment_initiations",
                                                         "clinical_assessments","hbsag_tests"))

# Need to change order of colour here to match:
x2_2 <- ggplot(subset(discounted_interactions_cost_median, scenario %in% c("a2_screen_2020_monit_sim7",
                                                              "a4_screen_2020_monit_sim7",
                                                              "a5_screen_2020_monit_sim7") &
                      !(interaction_type %in% c("total_interactions", "treatment_initiations")))) +
  geom_col(aes(x=age_group, y= cost/1000000, fill = interaction_type),
           width = 0.95, colour = "black") +
#  geom_errorbar(data=total_interactions_errorbar,
#                aes(x = age_group, ymin = lower/1000, ymax = upper/1000),width = 0.25) +
  labs(fill = "Resource utilisation") +
  scale_fill_manual(labels = c("hbsag_tests" = "Serological tests",
                                  "clinical_assessments" = "Assessments",
                                  "monitoring_assessments" = "Monitoring",
                                  "py_on_treatment" = "Treatment"),
                       values = c("hbsag_tests" = "#FDE725FF",
                                  "clinical_assessments" = "#35B779FF",
                                  "monitoring_assessments" = "#31688EFF",
                                  "py_on_treatment" = "#440154FF")) +
  guides(linetype=guide_legend(title=NULL),
         #fill=guide_legend(title=NULL),
         fill = FALSE) +
  theme_classic() +
  theme(legend.position=c(.78,.8)) +
  ylab("Discounted cost (millions)") +
  xlab("Screened age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Resource plots
grid.arrange(x2_1,x2_2, ncol = 2)


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


# Diagnostic plots for screening and monitoring by age ----
## Plots of breakdown of impact of 5-yearly monitoring and screening by separate age group (previously in monitoring frequency file) ----

# EFFECT OF MONITORING BY AGE (monitoring compared to no monitoring)

# DALYs averted by monitoring
dalys_averted_by_age_group <-
  plot_hbv_deaths_averted(counterfactual_object = out3_it,
                          scenario_objects = list(monit_out6,
                                                  monit_out8,
                                                  monit_out10),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "treatment programme without monitoring")
dalys_averted_by_age_group <- subset(dalys_averted_by_age_group, type == "number_averted" &
                                       by_year==2100) %>%
  select(scenario, sim, value)
dalys_averted_by_age_group$sim <- gsub("[^0-9]", "", dalys_averted_by_age_group$sim)
colnames(dalys_averted_by_age_group)[3] <- "dalys_averted"

deaths_averted_by_age_group <-
  plot_hbv_deaths_averted(counterfactual_object = out3_it,
                          scenario_objects = list(monit_out6,
                                                  monit_out8,
                                                  monit_out10),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "treatment programme without monitoring")
deaths_averted_by_age_group <- subset(deaths_averted_by_age_group, type == "number_averted" &
                                        by_year==2100) %>%
  select(scenario, sim, value)
# ISSUE WITH SIM NAMES HERE, USE THOSE FROM DALYS
#deaths_averted_by_age_group$sim <- gsub("[^0-9]", "", deaths_averted_by_age_group$sim)
deaths_averted_by_age_group$sim <-dalys_averted_by_age_group$sim
colnames(deaths_averted_by_age_group)[3] <- "deaths_averted"


# In terms of interactions, want the extra monitoring assessments, treatment initiations and
# person years of treatment due to monitoring

interactions_by_age_group <- rbind(
  cbind(scenario = "screen_2020_monit_sim6",
        left_join(gather(monit_out6$interactions[[16]]$total_assessed[-c(1:3)]-
                           out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out6$interactions[[16]]$total_treated[-c(1:3)]-
                           out3_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim8",
        left_join(gather(monit_out8$interactions[[16]]$total_assessed[-c(1:3)]-
                           out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out8$interactions[[16]]$total_treated[-c(1:3)]-
                           out3_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_sim10",
        left_join(gather(monit_out10$interactions[[16]]$total_assessed[-c(1:3)]-
                           out3_it$interactions[[16]]$total_assessed[-c(1:3)],
                         key = "sim", value = "monitoring_assessments"),
                  gather(monit_out10$interactions[[16]]$total_treated[-c(1:3)]-
                           out3_it$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")))
interactions_by_age_group$sim <- gsub("[^0-9]", "", interactions_by_age_group$sim)

py_on_treatment_by_age_group <-rbind(
  data.frame(scenario = "screen_2020_monit_sim6",
             sim = names(monit_out6$py_on_treatment[[16]]),
             py_on_treatment = monit_out6$py_on_treatment[[16]]-
               out3_it$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim8",
             sim = names(monit_out8$py_on_treatment[[16]]),
             py_on_treatment = monit_out8$py_on_treatment[[16]]-
               out3_it$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_sim10",
             sim = names(monit_out10$py_on_treatment[[16]]),
             py_on_treatment = monit_out10$py_on_treatment[[16]]-
               out3_it$py_on_treatment[[16]]))

df_by_age_group <- left_join(left_join(
  left_join(interactions_by_age_group, py_on_treatment_by_age_group,
            by = c("scenario", "sim")),
  dalys_averted_by_age_group, by = c("scenario", "sim")),
  deaths_averted_by_age_group, by = c("scenario", "sim"))
df_by_age_group$scenario <- factor(df_by_age_group$scenario)
levels(df_by_age_group$scenario) <- list("15-30" = "screen_2020_monit_sim6",
                                         "30-45" = "screen_2020_monit_sim8",
                                         "45+" = "screen_2020_monit_sim10")


# DALYs averted
p <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted)) +
  ylab("Inc. DALYs averted") +
  xlab("Monitored age group") +
  ylim(0,61000)

px <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted)) +
  ylab("Inc. HBV deaths averted") +
  xlab("Monitored age group") +
  ylim(0,4100)

p0 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = monitoring_assessments)) +
  ylab("Monitoring assessments") +
  xlab("Monitored age group") +
  ylim(0,400000)

# Incremental DALYs averted per monitoring assessment
p1 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/monitoring_assessments)) +
  ylab("Inc. DALYs averted\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,1)

# Incremental DALYs averted per treatment initiations
p2 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/treatment_initiations)) +
  ylab("Inc. DALYs averted\nper inc. treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,30)

# Incremental DALYs averted per PY on treatment
p3 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/py_on_treatment)) +
  ylab("Inc. DALYs averted\nper inc. PY on treatment") +
  xlab("Monitored age group") +
  ylim(0,0.75)

# Incremental deaths averted per monitoring assessment
p4 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/monitoring_assessments)) +
  ylab("Inc. HBV-related deaths averted\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,0.04)

# Incremental deaths  averted per treatment initiations
p5 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/treatment_initiations)) +
  ylab("Inc. HBV-related deaths averted\nper inc. treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,0.7)

# Incremental deaths  averted per PY on treatment
p6 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/py_on_treatment)) +
  ylab("Inc. HBV-related deaths averted\nper inc. PY on treatment") +
  xlab("Monitored age group") +
  ylim(0,0.05)

# Incremental treatment initiations per monitoring assessment
p7 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = treatment_initiations/monitoring_assessments)) +
  ylab("New treatment initiations\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,0.07)

# Incremental PY on treatment per monitoring assessment
p8 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/monitoring_assessments)) +
  ylab("Inc. PY on treatment\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,1.75)

# Incremental PY on treatment per treatment initiation
p9 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/treatment_initiations)) +
  ylab("Inc. PY on treatment\nper inc. treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,35)

grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 3)
grid.arrange(p7,p8, ncol = 1)

# Compare this with previous result!! Totally different because of change in
# initial treatment initiations (out3_it) - this would be case in all scenarios including
#

# Additional plots per carrier (especially monitoring assessments per carrier)
# would be helpful (this would be basically the cohort size over time in the different age groups)
# But monitoring assessments is kind of a proxy for carriers over time (except
# that age groups are not equally broad)
# Carriers at entry would also be interesting.
# Also per remaining life expectancy/average age at death? But this would
# need to be calculated in group of given age only (not entire cohort)

# Average age at death in screened+treated cohort
#ggplot(gather(rbind(out3_it$cohort_age_at_death, out6_it$cohort_age_at_death,
#                    out1_it$cohort_age_at_death),
#              key = "sim", value = "age", -scenario))+
#  geom_boxplot(aes(x=reorder(scenario, age),y=age))
#
#quantile(out6_it$cohort_age_at_death[,-1]-out3_it$cohort_age_at_death[,-1],
#         c(0.5,0.025,0.975))*12
# Average age at death postponed by 8 (3-20) months by monitoring every year
#quantile(out3_it$cohort_age_at_death[,-1]-out1_it$cohort_age_at_death[,-1],
#         c(0.5,0.025,0.975))
# Average age at death postponed by 1.5 (0.7-3) years by treating without monitoring
#quantile(((out3_it$cohort_cum_hbv_deaths[,-1]/out3_it$cohort_size[,-1])-
#            (out6_it$cohort_cum_hbv_deaths[,-1]/out6_it$cohort_size[,-1]))/
#           (out3_it$cohort_cum_hbv_deaths[,-1]/out3_it$cohort_size[,-1]),
#         c(0.5,0.025,0.975))
# Lifetime risk of dying from HBV reduced by 4% (1-8%) by monitoring every year
# in absolute terms or 52% (29-74%) relative to lifetime risk with no monitoring!

#quantile(((out1_it$cohort_cum_hbv_deaths[,-1]/out1_it$cohort_size[,-1])-
#            (out3_it$cohort_cum_hbv_deaths[,-1]/out3_it$cohort_size[,-1]))/
#           (out1_it$cohort_cum_hbv_deaths[,-1]/out1_it$cohort_size[,-1]),
#         c(0.5,0.025,0.975))
# Lifetime risk of dying from HBV reduced by 6% (3-12) by treating without monitoring
# in absolute terms or 45% (31-59%) relative to lifetime risk with no treatment!


## EFFECT OF SCREENING WITHOUT MONITORING BY AGE

# DALYs averted
dalys_averted_by_screened_age_group <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a2_out3_pop,
                                                  a4_out3_pop,
                                                  a5_out3_pop),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "No treatment")
dalys_averted_by_screened_age_group <- subset(dalys_averted_by_screened_age_group,
                                              type == "number_averted" &
                                                by_year==2100) %>%
  select(scenario, sim, value)
dalys_averted_by_screened_age_group$sim <- gsub("[^0-9]", "",
                                                dalys_averted_by_screened_age_group$sim)
colnames(dalys_averted_by_screened_age_group)[3] <- "dalys_averted"

deaths_averted_by_screened_age_group <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a2_out3_pop,
                                                  a4_out3_pop,
                                                  a5_out3_pop),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "No treatment")
deaths_averted_by_screened_age_group <- subset(deaths_averted_by_screened_age_group,
                                               type == "number_averted" &
                                                 by_year==2100) %>%
  select(scenario, sim, value)
deaths_averted_by_screened_age_group$sim <- gsub("[^0-9]", "", deaths_averted_by_screened_age_group$sim)
colnames(deaths_averted_by_screened_age_group)[3] <- "deaths_averted"

# In terms of interactions, want the HBsAg tests, initial clinical assessments,
# initial treatment initiations and person years of treatment due to initial treatment

interactions_by_screened_age_group <- rbind(
  cbind(scenario = "a2_screen_2020_monit_0",
        left_join(left_join(
          gather(a2_out3_pop$interactions[[16]]$total_screened[-c(1:3)],
                 key = "sim", value = "hbsag_tests"),
          gather(a2_out3_pop$interactions[[16]]$total_assessed[-c(1:3)],
                 key = "sim", value = "clinical_assessments"), by = "sim"),
          gather(a2_out3_pop$interactions[[16]]$total_treated[-c(1:3)],
                 key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "a4_screen_2020_monit_0",
        left_join(left_join(
          gather(a4_out3_pop$interactions[[16]]$total_screened[-c(1:3)],
                 key = "sim", value = "hbsag_tests"),
          gather(a4_out3_pop$interactions[[16]]$total_assessed[-c(1:3)],
                 key = "sim", value = "clinical_assessments"), by = "sim"),
          gather(a4_out3_pop$interactions[[16]]$total_treated[-c(1:3)],
                 key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "a5_screen_2020_monit_0",
        left_join(left_join(
          gather(a5_out3_pop$interactions[[16]]$total_screened[-c(1:3)],
                 key = "sim", value = "hbsag_tests"),
          gather(a5_out3_pop$interactions[[16]]$total_assessed[-c(1:3)],
                 key = "sim", value = "clinical_assessments"), by = "sim"),
          gather(a5_out3_pop$interactions[[16]]$total_treated[-c(1:3)],
                 key = "sim", value = "treatment_initiations"), by = "sim")))
interactions_by_screened_age_group$sim <- gsub("[^0-9]", "", interactions_by_screened_age_group$sim)

py_on_treatment_by_screened_age_group <-rbind(
  data.frame(scenario = "a2_screen_2020_monit_0",
             sim = names(a2_out3_pop$py_on_treatment[[16]]),
             py_on_treatment = a2_out3_pop$py_on_treatment[[16]]),
  data.frame(scenario = "a4_screen_2020_monit_0",
             sim = names(a4_out3_pop$py_on_treatment[[16]]),
             py_on_treatment = a4_out3_pop$py_on_treatment[[16]]),
  data.frame(scenario = "a5_screen_2020_monit_0",
             sim = names(a5_out3_pop$py_on_treatment[[16]]),
             py_on_treatment = a5_out3_pop$py_on_treatment[[16]]))

df_by_screened_age_group <- left_join(left_join(
  left_join(interactions_by_screened_age_group, py_on_treatment_by_screened_age_group,
            by = c("scenario", "sim")),
  dalys_averted_by_screened_age_group, by = c("scenario", "sim")),
  deaths_averted_by_screened_age_group, by = c("scenario", "sim"))
df_by_screened_age_group$scenario <- factor(df_by_screened_age_group$scenario)
levels(df_by_screened_age_group$scenario) <- list("15-30" = "a4_screen_2020_monit_0",
                                                  "30-45" = "a5_screen_2020_monit_0",
                                                  "45-65" = "a2_screen_2020_monit_0")

# Add per sAg test

# DALYs averted
ps <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted)) +
  ylab("DALYs averted") +
  xlab("Screened age group") +
  ylim(0,120000)

psx <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted)) +
  ylab("HBV deaths averted") +
  xlab("Screened age group") +
  ylim(0,4600)


ps0 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = clinical_assessments)) +
  ylab("Clinical assessments") +
  xlab("Screened age group") +
  ylim(0,50000)

# Incremental DALYs averted per clinical assessment
ps1 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/clinical_assessments)) +
  ylab("DALYs averted\nper clinical assessment") +
  xlab("Screened age group") +
  ylim(0,4.5)

# Incremental DALYs averted per treatment initiations
ps2 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/treatment_initiations)) +
  ylab("DALYs averted\nper treatment initiation") +
  xlab("Screened age group") +
  ylim(0,30)

# Incremental DALYs averted per PY on treatment
ps3 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/py_on_treatment)) +
  ylab("DALYs averted\nper PY on treatment") +
  xlab("Screened age group") +
  ylim(0,1)

# Incremental deaths averted per clinical assessment
ps4 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/clinical_assessments)) +
  ylab("HBV-related deaths averted\nper clinical assessment") +
  xlab("Screened age group") +
  ylim(0,0.16)

# Incremental deaths  averted per treatment initiations
ps5 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/treatment_initiations)) +
  ylab("HBV-related deaths averted\nper treatment initiation") +
  xlab("Screened age group") +
  ylim(0,0.8)

# Incremental deaths  averted per PY on treatment
ps6 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/py_on_treatment)) +
  ylab("HBV-related deaths averted\nper PY on treatment") +
  xlab("Screened age group") +
  ylim(0,0.04)

# Incremental treatment initiations per clinical assessment
ps7 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = treatment_initiations/clinical_assessments)) +
  ylab("Treatment initiations\nper clinical assessment") +
  xlab("Screened age group") +
  ylim(0,0.31)

# Incremental PY on treatment per clinical assessment
ps8 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/clinical_assessments)) +
  ylab("PY on treatment\nper clinical assessment") +
  xlab("Screened age group") +
  ylim(0,7.5)

# Incremental PY on treatment per treatment initiation
ps9 <- ggplot(df_by_screened_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/treatment_initiations)) +
  ylab("PY on treatment\nper treatment initiation") +
  xlab("Screened age group") +
  ylim(0,45)

grid.arrange(ps1,ps2,ps3,ps4,ps5,ps6, ncol = 3)

grid.arrange(ps7,ps9,ps8,p7,p9, p8, ncol = 3)

grid.arrange(ps1,ps2,ps3,
             p1,p2,p3,
             ps4,ps5,ps6,
             p4,p5,p6, ncol = 3)

grid.arrange(ps7,ps9,ps1,ps2,ps3,ps4,ps5,ps6,
             p7,p9,p1,p2,p3,p4,p5,p6, ncol = 8)

# Interactions
grid.arrange(ps7,ps9,p7,p9, ncol = 2)
# DALYS
grid.arrange(ps7,ps9,ps1,ps2,ps3,p7,p9,p1,p2,p3, ncol = 5)
# Deaths
grid.arrange(ps7,ps9,ps4,ps5,ps6,p7,p9,p4,p5,p6, ncol = 5)

# Screen:
grid.arrange(ps,psx,ps0,ps7,ps9,
             ps1,ps2,ps3,
             ps4,ps5,ps6,
             layout_matrix = rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5),
                                   c(6,6,6,6,6,7,7,7,7,7,8,8,8,8,8),
                                   c(9,9,9,9,9,10,10,10,10,10,11,11,11,11,11)))
# Monitor
grid.arrange(p,px,p0,p7,p9,
             p1,p2,p3,
             p4,p5,p6,
             layout_matrix = rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5),
                                   c(6,6,6,6,6,7,7,7,7,7,8,8,8,8,8),
                                   c(9,9,9,9,9,10,10,10,10,10,11,11,11,11,11)))


# Clinical assessments are proportional to carriers in the population in 2020.
# Monitoring assessments are proportional to carriers of that age over time.

# Age patterns in DALYS averted due to monitoring seem to be mostly explained
# by PY on treatment (young people have higher DALYS averted but also longer lifespan)
# Nevertheless DALYS averted per PY on treatment due to initial screening are still
# higher in the younger age group than among others. So far the only
# potential explanation I have found in sim 5 for this is that the proportion of treatment-eligible
# carriers with cirrhosis varies in that age group between the initial screen vs the
# monitoring. In 2020 the proportion of cirrhosis among treatment eligible carriers by age is:
# 27% (15-30), 16% (30-45), 17% (45-65). Of the monitoring events, the proportion in cirrhotic
# compartments of all treatment eligible compartments by age is:
# 9% (15-30), 8% (30-45), 12% (45-65). So possibly both the PY on treatment
# and DALYS averted depend on cirrhotic status.
# The higher proportion of treatment eligible carriers being cirrhotic in the younger
# age group is because IT is not treatment eligible. Yet also still depends on
# absolute number. Only way to find this out would be to simulate the
# treatment effect in those cirrhotic at first assessment vs non-cirrhotic.
# Nevertheless would expect a higher effect and longer life among those non-cirrhotic
# at treatment initiation. Pattern for YLL averted is the same.
# Clinical assessments corresponds to number of carriers and treatment initiations
# to number of treatment eligible carriers.

# Might want to show something like: % of total treatment initiations
# after first vs monitoring assessments by age

# Why does increasing monitoring frequency have diminishing returns? ----

# Also look at sensitivity analysis for this

# Hypothesis 1: incremental new treatment initiations per incremental assessments
# decline with increasing frequency
# (= monitoring more often does not identify more people to treat)

# Compare: every 20 years vs no monitoring, every 10 vs every 20, every 5 vs every 10,
# every 2 vs every 5, every 1 vs every 2

# out3_it, a1_out6_pop, out4b_it, out4_it, out5_it, out6a_it

# UPDATE THESE WITH RESIMULATED YEARLY AND 10-YEARLY MONITORING #

incremental_df <-rbind(
  data.frame(scenario = "screen_2020_monit_20",
             comparator = "screen_2020_monit_0",
             cohort_dalys = unlist(out3_it$cohort_dalys[,-1]-out4b_it$cohort_dalys[,-1]),
             py_on_treatment = unlist(out4b_it$py_on_treatment[[16]])-
                                  unlist(out3_it$py_on_treatment[[16]]),
             treatment_initiations = (unlist(out4b_it$interactions[[16]]$total_treated[,-c(1:3)])-
                                        unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)])),
             monitoring_assessments = unlist(out4b_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_10",
             comparator = "screen_2020_monit_20",
             cohort_dalys = unlist(out4b_it$cohort_dalys[,-1]-out4_it$cohort_dalys[,-1]),
             py_on_treatment = unlist(out4_it$py_on_treatment[[16]])-
               unlist(out4b_it$py_on_treatment[[16]]),
             treatment_initiations = (unlist(out4_it$interactions[[16]]$total_treated[,-c(1:3)])-
                                        unlist(out4b_it$interactions[[16]]$total_treated[,-c(1:3)])),
             monitoring_assessments = unlist(out4_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out4b_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             comparator = "screen_2020_monit_10",
             cohort_dalys = unlist(out4_it$cohort_dalys[,-1]-out5_it$cohort_dalys[,-1]),
             py_on_treatment = unlist(out5_it$py_on_treatment[[16]])-
               unlist(out4_it$py_on_treatment[[16]]),
             treatment_initiations = (unlist(out5_it$interactions[[16]]$total_treated[,-c(1:3)])-
                                        unlist(out4_it$interactions[[16]]$total_treated[,-c(1:3)])),
             monitoring_assessments = unlist(out5_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out4_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_2",
             comparator = "screen_2020_monit_5",
             cohort_dalys = unlist(out5_it$cohort_dalys[,-1]-out6a_it$cohort_dalys[,-1]),
             py_on_treatment = unlist(out6a_it$py_on_treatment[[16]])-
               unlist(out5_it$py_on_treatment[[16]]),
             treatment_initiations = (unlist(out6a_it$interactions[[16]]$total_treated[,-c(1:3)])-
                                        unlist(out5_it$interactions[[16]]$total_treated[,-c(1:3)])),
             monitoring_assessments = unlist(out6a_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out5_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             comparator = "screen_2020_monit_2",
             cohort_dalys = unlist(out6a_it$cohort_dalys[,-1]-a1_out6_pop$cohort_dalys[,-1]),
             py_on_treatment = unlist(a1_out6_pop$py_on_treatment[[16]])-
               unlist(out6a_it$py_on_treatment[[16]]),
             treatment_initiations = (unlist(a1_out6_pop$interactions[[16]]$total_treated[,-c(1:3)])-
                                        unlist(out6a_it$interactions[[16]]$total_treated[,-c(1:3)])),
             monitoring_assessments = unlist(a1_out6_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out6a_it$interactions[[16]]$total_assessed[,-c(1:3)]))
)

incremental_df$treatment_initiations_per_assessment <- incremental_df$treatment_initiations/
  incremental_df$monitoring_assessments
incremental_df$py_on_treatment_per_assessment <- incremental_df$py_on_treatment/
  incremental_df$monitoring_assessments
incremental_df$dalys_averted_per_treatment_initiations <- incremental_df$cohort_dalys/
  incremental_df$treatment_initiations
incremental_df$dalys_averted_per_py_on_treatment <- incremental_df$cohort_dalys/
  incremental_df$py_on_treatment
incremental_df$scenario <- factor(incremental_df$scenario,
                                               levels = c("screen_2020_monit_20",
                                                          "screen_2020_monit_10",
                                                          "screen_2020_monit_5",
                                                          "screen_2020_monit_2",
                                                          "screen_2020_monit_1"))

ggplot(incremental_df) +
  geom_boxplot(aes(x=scenario, y = treatment_initiations_per_assessment))


ggplot(incremental_df) +
  geom_boxplot(aes(x=scenario, y = treatment_initiations_per_assessment))

# From 10 yearly monitoring onwards, the ratio of new treatment initiations per
# incremental monitoring assessment decreases
# Though it increases between monitoring every 20 vs every 10 years
# The incremental person-time on treatment per assessment decreases from 20-yearly

# other way to visualise this:
# Here everything is compared to no monitoring
interactions_df <-rbind(
  data.frame(scenario = "screen_2020_monit_20",
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-out4b_it$cohort_dalys[,-1]),
             treatment_initiations = unlist(out4b_it$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(out4b_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_10",
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-out4_it$cohort_dalys[,-1]),
             treatment_initiations = unlist(out4_it$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(out4_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-out5_it$cohort_dalys[,-1]),
             treatment_initiations = unlist(out5_it$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(out5_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_2",
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-out6a_it$cohort_dalys[,-1]),
             treatment_initiations = unlist(out6a_it$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(out6a_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-a1_out6_pop$cohort_dalys[,-1]),
             treatment_initiations = unlist(a1_out6_pop$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(a1_out6_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)]))
)

ggplot(interactions_df) +
  geom_point(aes(x=monitoring_assessments, y = treatment_initiations, colour = scenario))

# Hypothesis 2: The incremental impact per additional treatment initiations declines with
# increasing frequency
# (= the disease progression identified by more monitoring is less affected by treatment,
# e.g. because cirrhotic cases are treated immediately).

ggplot(incremental_df) +
  geom_boxplot(aes(x=scenario, y = dalys_averted_per_treatment_initiations))

ggplot(incremental_df) +
  geom_boxplot(aes(x=scenario, y = dalys_averted_per_py_on_treatment))

# Other visualisation:
ggplot(interactions_df) +
  geom_point(aes(x=treatment_initiations, y = cohort_dalys_averted, colour = scenario))


boxplot(unlist(out3_it$cohort_dalys[,-1]-out4_it$cohort_dalys[,-1])/
          (unlist(out4_it$interactions[[16]]$total_treated[,-c(1:3)])-
             unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)])),
        ylim=c(0,13))

boxplot(unlist(out4_it$cohort_dalys[,-1]-out5_it$cohort_dalys[,-1])/
          (unlist(out5_it$interactions[[16]]$total_treated[,-c(1:3)])-
             unlist(out4_it$interactions[[16]]$total_treated[,-c(1:3)])),
        ylim=c(0,13))

boxplot(unlist(out5_it$cohort_dalys[,-1]-out6_it$cohort_dalys[,-1])/
          (unlist(out6_it$interactions[[16]]$total_treated[,-c(1:3)])-
             unlist(out5_it$interactions[[16]]$total_treated[,-c(1:3)])),
        ylim=c(0,13))

t <-data.frame(scenario = rep(c("10", "5", "2", "1"), each=183),
               dalys=c(unlist(out3_it$cohort_dalys[,-1]-out4_it$cohort_dalys[,-1]),
                       unlist(out3_it$cohort_dalys[,-1]-out5_it$cohort_dalys[,-1]),
                       unlist(out3_it$cohort_dalys[,-1]-out6a_it$cohort_dalys[,-1]),
                       unlist(out3_it$cohort_dalys[,-1]-out6_it$cohort_dalys[,-1])),
               treatment = c((unlist(out4_it$py_on_treatment[[16]])-
                                unlist(out3_it$py_on_treatment[[16]])),
                             (unlist(out5_it$py_on_treatment[[16]])-
                                unlist(out3_it$py_on_treatment[[16]])),
                             (unlist(out6a_it$py_on_treatment[[16]])-
                                unlist(out3_it$py_on_treatment[[16]])),
                             (unlist(out6_it$py_on_treatment[[16]])-
                                unlist(out3_it$py_on_treatment[[16]]))))

ggplot(t) +
  stat_ellipse(aes(x=treatment, y = dalys, colour=scenario)) +
  theme_classic()

t2 <-data.frame(dalys_10=unlist(out3_it$cohort_dalys[,-1]-out4_it$cohort_dalys[,-1]),
                dalys_5=unlist(out4_it$cohort_dalys[,-1]-out5_it$cohort_dalys[,-1]),
                dalys_2=unlist(out5_it$cohort_dalys[,-1]-out6a_it$cohort_dalys[,-1]),
                dalys_1=unlist(out6a_it$cohort_dalys[,-1]-out6_it$cohort_dalys[,-1]),
                py_10 = (unlist(out4_it$py_on_treatment[[16]])-
                           unlist(out3_it$py_on_treatment[[16]])),
                py_5=(unlist(out5_it$py_on_treatment[[16]])-
                        unlist(out4_it$py_on_treatment[[16]])),
                py_2=(unlist(out6a_it$py_on_treatment[[16]])-
                        unlist(out5_it$py_on_treatment[[16]])),
                py_1=(unlist(out6_it$py_on_treatment[[16]])-
                        unlist(out6a_it$py_on_treatment[[16]])))
t2$ratio_10 <- t2$dalys_10/t2$py_10
t2$ratio_5 <- t2$dalys_5/t2$py_5
t2$ratio_2 <- t2$dalys_2/t2$py_2
t2$ratio_1 <- t2$dalys_1/t2$py_1

boxplot(t2$ratio_10, ylim=c(0,1))
boxplot(t2$ratio_5, ylim=c(0,1))
boxplot(t2$ratio_2, ylim=c(0,1))
boxplot(t2$ratio_1, ylim=c(0,1))

# The ratio of incremental DALYS averted per incremental time on treatment  declines
# with increasing monitoring frequencies.


