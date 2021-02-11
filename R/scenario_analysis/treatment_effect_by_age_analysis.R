# Analyse the effect of screening/monitoring/treatment in different age groups

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)

# For all these plots, check values line up with total!!!

# Load data ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/"

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



