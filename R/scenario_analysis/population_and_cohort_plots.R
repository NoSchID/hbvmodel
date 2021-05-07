# Analyse the effect of screening/monitoring/treatment in different age groups

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(grid)
#source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Load data ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/"

out_path_monit <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

out_path_repeat_screen <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/repeat_screening_analysis/"


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
out1_it <- readRDS(paste0(out_path_monit, "a1_it_out1_status_quo_cohort_200121.rds"))
out1_it <- out1_it[[1]]

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

a1_out6_pop <- readRDS(paste0(out_path_monit, "a1_it_out6_screen_2020_monit_1_240221.rds"))
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

# Every 5 years in <45 year olds
monit_out7 <- readRDS(paste0(out_path_monit, "a1_it_screen_2020_monit_out7_050321.rds"))
monit_out7 <- monit_out7[[1]]

# Different monitoring frequencies across all ages
# out3_it above
# a1_out6_pop
out4b_it <- readRDS(paste0(out_path_monit, "a1_it_out4b_screen_2020_monit_20_230221.rds"))
out4b_it <- out4b_it[[1]]   # 20 years
out4_it <- readRDS(paste0(out_path_monit, "a1_it_out4_screen_2020_monit_10_240221.rds"))
out4_it <- out4_it[[1]]   # 10 years
out5_it <- readRDS(paste0(out_path_monit, "a1_it_out5_screen_2020_monit_5_250221.rds"))
out5_it <- out5_it[[1]]   # 5 years
out6a_it <- readRDS(paste0(out_path_monit, "a1_it_out6a_screen_2020_monit_2_250221.rds"))
out6a_it <- out6a_it[[1]]  # 2 years

# Repeat screening strategies
# monit_sim7 repeated in 2020+2030
out8b_2030_monit_sim7 <- readRDS(paste0(out_path_repeat_screen, "a1_it_out8b_monit_sim7_screen_10b_2030_290321.rds"))
out8b_2030_monit_sim7 <- out8b_2030_monit_sim7[[1]]

### Treatment effect in the screened+treated cohort (HCC delayed vs prevented) ----
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

# HBV deaths combined by age
a1_x1 <- a4_out1$cum_cohort_hbv_deaths_male_by_age_2100+
  a4_out1$cum_cohort_hbv_deaths_female_by_age_2100+
  a5_out1$cum_cohort_hbv_deaths_male_by_age_2100+
  a5_out1$cum_cohort_hbv_deaths_female_by_age_2100+
  a2_out1$cum_cohort_hbv_deaths_male_by_age_2100+
  a2_out1$cum_cohort_hbv_deaths_female_by_age_2100
a1_x2 <- a4_out3$cum_screened_hbv_deaths_male_by_age_2100 +
  a4_out3$cum_screened_hbv_deaths_female_by_age_2100+
  a5_out3$cum_screened_hbv_deaths_male_by_age_2100 +
  a5_out3$cum_screened_hbv_deaths_female_by_age_2100+
  a2_out3$cum_screened_hbv_deaths_male_by_age_2100 +
  a2_out3$cum_screened_hbv_deaths_female_by_age_2100
a1_x3 <- a1_x1-a1_x2
a1_x4 <- a4_out3$cum_treated_hbv_deaths_male_by_age_2100+
  a4_out3$cum_treated_hbv_deaths_female_by_age_2100+
  a5_out3$cum_treated_hbv_deaths_male_by_age_2100+
  a5_out3$cum_treated_hbv_deaths_female_by_age_2100+
  a2_out3$cum_treated_hbv_deaths_male_by_age_2100+
  a2_out3$cum_treated_hbv_deaths_female_by_age_2100
a1_x3$sim <- rownames(a1_x3)
a1_x4$sim <- rownames(a1_x4)

df_a1 <- rbind(
  data.frame(type="No treatment",
             gather(a1_x3,key="age", value = "value",-sim)),
  data.frame(type="With treatment",
             gather(a1_x4,key="age", value = "value",-sim)))
df_a1$age <- rep(ages, each = 183)

ggplot(df_a1) +
  stat_summary(aes(x=age, y = value/0.5, colour=type), fun="median",
               geom = "line", size = 1)+
  stat_summary(aes(x=age, y = value/0.5, fill=type, colour = type),
               fun.min= function(x) quantile(x,0.025),
               fun.max= function(x) quantile(x,0.975),
               geom = "ribbon", alpha = 0.1, lty ="dashed")+
  ylab("Cumulative HBV-related deaths by 2100") +
  xlab("Age at HBV death (years)") +
  #labs(title="In the cohort") +
  geom_vline(xintercept=30) +
  geom_vline(xintercept=45) +
  geom_rect(aes(xmin=0, xmax=15, ymin=-Inf, ymax=Inf), fill = "grey") +
  scale_fill_manual("Scenario", values= c("No treatment"= "#A180A9",
                                          "With treatment" ="#1F968BFF")) +
  scale_color_manual("Scenario", values= c("No treatment"= "#A180A9",
                                           "With treatment" ="#1F968BFF")) +
  theme_bw()+
  theme(legend.position = c(0.85, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

group_by(df_a1,type, age) %>%
  summarise(median(value/0.5)) %>%
  filter(age==45)

ggplot(df_a1) +
  stat_summary(data=subset(df_a1, type=="No treatment"),
      aes(x=age, y = value/0.5, fill="Averted deaths"), fun="median",
               geom = "area", size = 1, alpha = 1)+
  stat_summary(data=subset(df_a1, type=="With treatment"),
               aes(x=age, y = value/0.5), fun="median", fill="white",
               geom = "area", size = 1, alpha = 1)+
  stat_summary(aes(x=age, y = value/0.5, colour=type), fun="median",
               geom = "line", size = 2)+
#  stat_summary(aes(x=age, y = value/0.5, fill=type, colour = type),
 #              fun.min= function(x) quantile(x,0.025),
#               fun.max= function(x) quantile(x,0.975),
#               geom = "ribbon", alpha = 0, lty ="dashed")+
  ylab("Cumulative HBV-related deaths by 2100") +
  #labs(title="In the cohort") +
  geom_segment(aes(x = 30, y = 3.93, xend = 30, yend = 48.6)) +
  geom_segment(aes(x = 45, y = 16.6, xend = 45, yend = 131)) +
  geom_hline(yintercept=0) +
  scale_fill_manual("", values= c("Averted deaths" = "grey70")) +
  scale_color_manual("Scenario", values= c("No treatment"= "#A180A9",
                                           "With treatment" ="#1F968BFF")) +
  scale_x_continuous("Age at HBV death (years)", breaks=c(15,30,45,90),
                     limits=c(15,90)) +
  theme_bw()+
  theme(legend.position = c(0.85, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))


# What proportion of averted deaths would have occurred before the age of 30/45?
deaths_averted_distribution_by_age <- rbind(
  data.frame(age_group = "15-30",
             value = apply(a1_x3[,which(ages==15):which(ages==29.5)]-
                            a1_x4[,which(ages==15):which(ages==29.5)],1,sum)), # total deaths averted by age
  data.frame(age_group = "30-45",
             value = apply(a1_x3[,which(ages==30):which(ages==44.5)]-
                             a1_x4[,which(ages==30):which(ages==44.5)],1,sum)),
  data.frame(age_group = "45+",
             value = apply(a1_x3[,which(ages==45):which(ages==99.5)]-
                             a1_x4[,which(ages==45):which(ages==99.5)],1,sum))
  )


quantile(apply(a1_x3[,which(ages==15):which(ages==29.5)]-a1_x4[,which(ages==15):which(ages==29.5)],1,sum)/
  apply(a1_x3[,1:200]-a1_x4[,1:200],1,sum), c(0.5,0.025,0.975))
quantile(apply(a1_x3[,which(ages==15):which(ages==44.5)]-a1_x4[,which(ages==15):which(ages==44.5)],1,sum)/
           apply(a1_x3[,1:200]-a1_x4[,1:200],1,sum), c(0.5,0.025,0.975))

deaths_averted_distribution_by_age2 <- gather(a1_x3[,1:200]-
                             a1_x4[,1:200], key = "age", value = "value")
deaths_averted_distribution_by_age2$age <- rep(ages, each = 183)
deaths_averted_distribution_by_age2$age_group <- "15-30"
deaths_averted_distribution_by_age2$age_group[deaths_averted_distribution_by_age2$age>=30 &
                                                deaths_averted_distribution_by_age2$age<45] <- "30-45"
deaths_averted_distribution_by_age2$age_group[deaths_averted_distribution_by_age2$age>=45] <- "45+"


pp7 <- ggplot(subset(deaths_averted_distribution_by_age2, age %in% 0:99)) +
  stat_summary(aes(x=age, y = value*2, fill = age_group), fun="median",
               geom = "col", width = 1) +
  ylab("Averted HBV-related deaths") +
  xlab("Age at death (years)") +
  scale_x_continuous(breaks=c(15,30,45,60,75),
                     limits=c(14,90)) +
  scale_fill_viridis_d("", begin=0.6, end = 0.8) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.70),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        plot.margin = unit(c(5.5,35,5.5,5.5), "pt"))
# COuld change this plot to show deaths averted BY screening among the different groups

ggplot(deaths_averted_distribution_by_age) +
  stat_summary(aes(x="15-65", y = value, fill = reorder(age_group,-value)), fun="median",
               geom = "bar", position = "fill") +
  ylab("Percentage of averted HBV deaths") +
  xlab("Screened age group (years)") +
  scale_fill_viridis_d("Age") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

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




# Distribution of HBV deaths by age in 2020: absolute number and rate in SQ scenario ----

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

hcc_deaths_2020 <-
  (do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                           function(x) x[which(out2_disease_outcomes$time==2021),]))-
     do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                             function(x) x[which(out2_disease_outcomes$time==2020),])))+
  (do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                           function(x) x[which(out2_disease_outcomes$time==2021),]))-
     do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                             function(x) x[which(out2_disease_outcomes$time==2020),])))


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

# Proportion of HBV deaths due to HCC:
total_hcc_deaths_2020 <-
  # Incident HCC deaths in men
    apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                           function(x) x[which(out2_disease_outcomes$time==2021),]))-
     do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                             function(x) x[which(out2_disease_outcomes$time==2020),]))),1,sum) +
  # Incident HCC deaths in women
  apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                           function(x) x[which(out2_disease_outcomes$time==2021),]))-
     do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                             function(x) x[which(out2_disease_outcomes$time==2020),]))),1,sum)

total_hbv_deaths_2020 <-
  # Incident HBV deaths in men
  apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                                 function(x) x[which(out2_disease_outcomes$time==2021),]))-
           do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                                   function(x) x[which(out2_disease_outcomes$time==2020),]))),1,sum) +
  # Incident HBV deaths in women
  apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_female,
                                 function(x) x[which(out2_disease_outcomes$time==2021),]))-
           do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_female,
                                   function(x) x[which(out2_disease_outcomes$time==2020),]))),1,sum)

prop_deaths_from_hcc <- total_hcc_deaths_2020/total_hbv_deaths_2020
quantile(prop_deaths_from_hcc, c(0.50,0.025,0.975))

# Proportion of deaths in men:
quantile(apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                               function(x) x[which(out2_disease_outcomes$time==2021),]))-
         do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                                 function(x) x[which(out2_disease_outcomes$time==2020),]))),1,sum)/
  total_hbv_deaths_2020, c(0.5,0.025,0.975))

# Cirrhosis mortality compared to GBD (2018)
hcc_deaths_2018 <- (apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                                                   function(x) x[which(out2_disease_outcomes$time==2019),]))-
                             do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                                                     function(x) x[which(out2_disease_outcomes$time==2018),]))),1,sum) +
                      apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                                                     function(x) x[which(out2_disease_outcomes$time==2019),]))-
                               do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                                                       function(x) x[which(out2_disease_outcomes$time==2018),]))),1,sum))

cirrhosis_deaths_2018 <- (apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                                                        function(x) x[which(out2_disease_outcomes$time==2019),]))-
                                  do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                                                          function(x) x[which(out2_disease_outcomes$time==2018),]))),1,sum) +
  apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_female,
                                 function(x) x[which(out2_disease_outcomes$time==2019),]))-
           do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_female,
                                   function(x) x[which(out2_disease_outcomes$time==2018),]))),1,sum))-
  hcc_deaths_2018

pop_male <- lapply(out2_carriers, "[[", "pop_male")
pop_female <- lapply(out2_carriers, "[[", "pop_female")

pop_2018 <-
  apply(do.call(rbind.data.frame, (lapply(pop_male, function(x)
    x[which(out2_carriers[[1]]$time==2018.5),]))),1,sum)+
  apply(do.call(rbind.data.frame, (lapply(pop_female, function(x)
    x[which(out2_carriers[[1]]$time==2018.5),]))),1,sum)

quantile(cirrhosis_deaths_2018/pop_2018*100000, c(0.5,0.025,0.975))
quantile(hcc_deaths_2018/pop_2018*100000, c(0.5,0.025,0.975))

# Male HCC mort rate:
quantile(apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                                function(x) x[which(out2_disease_outcomes$time==2019),]))-
          do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                                  function(x) x[which(out2_disease_outcomes$time==2018),]))),1,sum)/
    apply(do.call(rbind.data.frame, (lapply(pop_male, function(x)
      x[which(out2_carriers[[1]]$time==2018.5),]))),1,sum)*100000, c(0.5,0.025,0.975))
# Female HCC mort rate:
quantile(apply((do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                                        function(x) x[which(out2_disease_outcomes$time==2019),]))-
                  do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                                          function(x) x[which(out2_disease_outcomes$time==2018),]))),1,sum)/
           apply(do.call(rbind.data.frame, (lapply(pop_female, function(x)
             x[which(out2_carriers[[1]]$time==2018.5),]))),1,sum)*100000, c(0.5,0.025,0.975))


# Cirrhosis deaths by age
cirrhosis_deaths_by_age_2018 <- ((do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                                                         function(x) x[which(out2_disease_outcomes$time==2019),]))-
                                   do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_male,
                                                           function(x) x[which(out2_disease_outcomes$time==2018),]))) +
                            (do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_female,
                                                           function(x) x[which(out2_disease_outcomes$time==2019),]))-
                                     do.call("rbind", lapply(out2_disease_outcomes$cum_hbv_deaths_female,
                                                             function(x) x[which(out2_disease_outcomes$time==2018),]))))-
  ((do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                                  function(x) x[which(out2_disease_outcomes$time==2019),]))-
            do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_male,
                                    function(x) x[which(out2_disease_outcomes$time==2018),]))) +
     (do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                                    function(x) x[which(out2_disease_outcomes$time==2019),]))-
              do.call("rbind", lapply(out2_disease_outcomes$cum_hcc_deaths_female,
                                      function(x) x[which(out2_disease_outcomes$time==2018),]))))

carriers_male <- lapply(out2_carriers, "[[", "carriers_male")
carriers_female <- lapply(out2_carriers, "[[", "carriers_female")

carriers_by_age_2018 <-
  do.call(rbind.data.frame, (lapply(carriers_male, function(x)
    x[which(out2_carriers[[1]]$time==2018.5),])))+
  do.call(rbind.data.frame, (lapply(carriers_female, function(x)
    x[which(out2_carriers[[1]]$time==2018.5),])))

cirrhosis_mort_carriers <- cirrhosis_deaths_by_age_2018/carriers_by_age_2018
cirrhosis_mort_carriers$sim <- rownames(cirrhosis_mort_carriers)
cirrhosis_mort_carriers <- gather(cirrhosis_mort_carriers, key = "age", value = "value", -sim)
cirrhosis_mort_carriers$age <- rep(ages, each = 183)

ggplot(cirrhosis_mort_carriers) +
  geom_line(aes(x=age, y = value, group = sim), col = "grey") +
  stat_summary(aes(x=age, y = value), fun = "median", geom="line", col = "red") +
  theme_classic()

# Carriers:

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

### Population effects of screening and monitoring by age and over time ----
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

# Proportion eligible among adults:
quantile((apply(carriers_eligible_by_age_2020[,which(ages==15):which(ages==99.5)],1,sum)+
  apply(carriers_it_by_age_2020[,which(ages==30):which(ages==99.5)],1,sum))/
  apply(total_carriers_by_age_2020[,which(ages==15):which(ages==99.5)],1,sum),
  c(0.5,0.025,0.975))

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

# Test for DALYS averted per PY on treatment

test <- data.frame(subset(interactions_by_age, scenario %in% c("a4_screen_2020_monit_sim7",
                                            "a5_screen_2020_monit_sim7",
                                            "a2_screen_2020_monit_sim7") &
         interaction_type=="py_on_treatment"),
      dalys_averted = subset(outcomes_by_age, scenario %in% c("a4_screen_2020_monit_sim7",
                                              "a5_screen_2020_monit_sim7",
                                              "a2_screen_2020_monit_sim7") &
               outcome=="dalys_averted")$value)
test$dalys_averted_per_py_on_treatment <- test$dalys_averted/test$value

ggplot(test) +
  geom_boxplot(aes(x=age_group, y = dalys_averted_per_py_on_treatment)) +
  ylim(0,0.8)

# Test for distribution of carriers by treatment eligibility and cirrhosis status

# Add treatment eligible carriers by age in 2020
carriers_cirrhotic_male_by_age_2020 <-
#  do.call(rbind.data.frame,
#          lapply(out2_comps_by_age$cc_female,
#                 function(x) x[which(out2_comps_by_age$time==2020),]))+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$cc_male,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
#  do.call(rbind.data.frame,
 #         lapply(out2_comps_by_age$dcc_female,
 #                function(x) x[which(out2_comps_by_age$time==2020),]))#+
  do.call(rbind.data.frame,
          lapply(out2_comps_by_age$dcc_male,
                 function(x) x[which(out2_comps_by_age$time==2020),]))

carriers_cirrhotic_female_by_age_2020 <-
   do.call(rbind.data.frame,
          lapply(out2_comps_by_age$cc_female,
                 function(x) x[which(out2_comps_by_age$time==2020),]))+
   do.call(rbind.data.frame,
         lapply(out2_comps_by_age$dcc_female,
                function(x) x[which(out2_comps_by_age$time==2020),]))

# Combine into dataframes
carriers_cirrhotic_by_age_group_2020 <-
  rbind(data.frame(age_group = "15-30",
                   sim = rownames(carriers_cirrhotic_by_age_2020[,which(ages==15): which(ages==30-da)]),
                   carriers_cirrhotic_male = rowSums(carriers_cirrhotic_male_by_age_2020[,which(ages==15): which(ages==30-da)]),
        carriers_cirrhotic_female = rowSums(carriers_cirrhotic_female_by_age_2020[,which(ages==15): which(ages==30-da)])),
        data.frame(age_group = "30-45",
                   sim = rownames(carriers_cirrhotic_by_age_2020[,which(ages==30): which(ages==45-da)]),
                   carriers_cirrhotic_male  = rowSums(carriers_cirrhotic_male_by_age_2020[,which(ages==30): which(ages==45-da)]),
                   carriers_cirrhotic_female = rowSums(carriers_cirrhotic_female_by_age_2020[,which(ages==30): which(ages==45-da)])),
        data.frame(age_group = "45-65",
                   sim = rownames(carriers_cirrhotic_by_age_2020[,which(ages==45): which(ages==65-da)]),
                   carriers_cirrhotic_male  = rowSums(carriers_cirrhotic_male_by_age_2020[,which(ages==45): which(ages==65-da)]),
                   carriers_cirrhotic_female = rowSums(carriers_cirrhotic_female_by_age_2020[,which(ages==45): which(ages==65-da)]))
  )

test2 <- left_join(left_join(carriers_by_age_group_2020, carriers_eligible_by_age_group_2020,
          by = c("age_group", "sim")),
          carriers_cirrhotic_by_age_group_2020, by = c("age_group", "sim"))
test2$carriers_eligible_noncirrhotic <- test2$carriers_eligible-test2$carriers_cirrhotic_male-
  test2$carriers_cirrhotic_female
#test2$prop_eligible_cirrhotic <- test2$carriers_cirrhotic/test2$carriers
#test2$prop_eligible_noncirrhotic <- test2$carriers_eligible_noncirrhotic/test2$carriers
test2p <- gather(test2, key="outcome", value = "value", -age_group, -sim)

ggplot(subset(test2p, outcome  %in% c("carriers_eligible_noncirrhotic",
                                      "carriers_cirrhotic_male","carriers_cirrhotic_female",
                                      "carriers"))) +
  stat_summary(aes(x=age_group, y = value, fill = reorder(outcome, -value)),
               fun= "median", geom="bar",
               position="fill") +
  geom_hline(yintercept=0.04)

ggplot(subset(test2p, outcome %in% c("carriers_eligible_noncirrhotic",
                                     "carriers_cirrhotic"))) +
  stat_summary(aes(x=age_group, y = value, fill = reorder(outcome, -value)), fun= "median", geom="bar",
               position="fill")

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


# Paper panel plot of no monitoring programme ----

# PAPER PLOT / THESIS PLOT

# Need to run: functions in moniotirng freq script, previous 2 sections,
# pp7 in Cumulative HBV/HCC deaths by age by 2100 to show shift in pattern section

pp1 <- ggplot(data= subset(outcomes_by_age, scenario %in% c("a2_screen_2020_monit_0",
                                                     "a4_screen_2020_monit_0",
                                                     "a5_screen_2020_monit_0") &
                      outcome == "deaths_averted")) +
  stat_summary(aes(x=age_group, y= value/1000, fill = outcome),
               fun="median", geom="bar", position = "dodge2", width = 0.95, colour="black")+
  stat_summary(aes(x=age_group, y= value/1000),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", position = "dodge2", width = 0.15)+
  scale_fill_manual(values=c("deaths_averted" = "#31688EFF")) +
  guides(fill=FALSE) +
  xlab("Screened age group (years)")+
  ylab("HBV-related deaths averted\n(thousands)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),  #13.5
        axis.title = element_text(size = 15), #13.5
        strip.text = element_text(size = 14))

pp2 <- ggplot(data= subset(outcomes_by_age, scenario %in% c("a2_screen_2020_monit_0",
                                                     "a4_screen_2020_monit_0",
                                                     "a5_screen_2020_monit_0") &
                      outcome == "dalys_averted")) +
  stat_summary(aes(x=age_group, y= value/1000, fill = outcome),
               fun="median", geom="bar", position = "dodge2", width = 0.95, colour="black")+
  stat_summary(aes(x=age_group, y= value/1000),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", position = "dodge2", width = 0.15)+
  scale_fill_manual(values=c("dalys_averted" = "#31688EFF")) +
  guides(fill=FALSE) +
  xlab("Screened age group (years)")+
  ylab("DALYs averted \n(thousands)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14))

# Affected population by stage of care
pp3 <- ggplot() +
  geom_col(data=subset(interactions_median, scenario %in% c("a2_screen_2020_monit_0",
                                                            "a4_screen_2020_monit_0",
                                                            "a5_screen_2020_monit_0") &
                         !(interaction_type %in% c("total_interactions", "py_on_treatment",
                                                   "monitoring_assessments"))),
           aes(x=age_group, y= value/1000, fill = interaction_type),
           width = 0.95, colour = "black") +
  #geom_errorbar(data=subset(total_interactions_errorbar,scenario %in% c("a2_screen_2020_monit_sim7",
  #                                                                        "a4_screen_2020_monit_sim7",
  #                                                                        "a5_screen_2020_monit_sim7")),
  #                aes(x = age_group, ymin = lower/1000, ymax = upper/1000), width = 0.15) +
  labs(fill = "Resource utilisation") +
  scale_fill_viridis_d(labels = c("hbsag_tests" = "Screening",
                               "clinical_assessments" = "Assessment",
                               "treatment_initiations" = "Treatment")) +
  #scale_linetype_manual(values = c("HBV carriers" = "dashed")) +
  guides(linetype=guide_legend(title=NULL),
         fill=guide_legend(title=NULL)) +
  theme_classic() +
  theme(legend.position=c(.75,.8)) +
  ylab("Population affected\n(thousands)") +
  xlab("Screened age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12.5),  #12
        legend.title = element_text(size = 14))

# Cost by stage of care
pp4 <- ggplot(subset(discounted_interactions_cost_median, scenario %in% c("a2_screen_2020_monit_0",
                                                                   "a4_screen_2020_monit_0",
                                                                   "a5_screen_2020_monit_0") &
                !(interaction_type %in% c("total_interactions", "treatment_initiations",
                                          "monitoring_assessments")))) +
  geom_col(aes(x=age_group, y= cost/1000000, fill = interaction_type),
           width = 0.95, colour = "black") +
  #  geom_errorbar(data=total_interactions_errorbar,
  #                aes(x = age_group, ymin = lower/1000, ymax = upper/1000),width = 0.25) +
  labs(fill = "Resource utilisation") +
  scale_fill_viridis_d(labels = c("hbsag_tests" = "Serological tests",
                                  "clinical_assessments" = "Clinical assessments",
                                  "py_on_treatment" = "Treatment")) +
#  scale_fill_manual(labels = c("hbsag_tests" = "Serological tests",
#                               "clinical_assessments" = "Assessments",
#                               "monitoring_assessments" = "Monitoring",
#                               "py_on_treatment" = "Treatment"),
#                    values = c("hbsag_tests" = "#FDE725FF",
#                               "clinical_assessments" = "#35B779FF",
#                               "monitoring_assessments" = "#31688EFF",
#                               "py_on_treatment" = "#440154FF")) +
  scale_y_continuous(breaks=c(0,2,4,6,8)) +
  guides(linetype=guide_legend(title=NULL),
         #fill=guide_legend(title=NULL),
         fill = FALSE) +
  theme_classic() +
  theme(legend.position=c(.78,.8)) +
  ylab("Discounted cost\n(million US$)") +
  xlab("Screened age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12.5),
        legend.title = element_text(size = 14))

prevalence_by_age_group_2020_summary <- prevalence_by_age_group_2020 %>%
  group_by(age_group) %>%
  summarise(median=median(carriers),
            lower = quantile(carriers, 0.025),
            upper = quantile(carriers, 0.975))

# Plot distribution of carriers
pp5 <- ggplot(carriers_by_age_group_2020) +
  stat_summary(aes(x=age_group, y= carriers/1000),
               fun="median", geom="bar", position = "dodge2", #width = 0.95,
               colour="black", fill = "grey95") +
  stat_summary(aes(x=age_group, y= carriers/1000),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", position = "dodge2", width = 0.15)+
 geom_pointrange(data = prevalence_by_age_group_2020_summary,
                 aes(x=c(1.15,2.15,3.15), y= median*10, ymin=lower *10,
                     ymax=upper*10, colour = age_group),
                size = 1)+    # "#55C667FF", colour="#31688EFF",
  geom_errorbar(data = prevalence_by_age_group_2020_summary,
                  aes(x=c(1.15,2.15,3.15), y= median*10, ymin=lower *10,
                      ymax=upper*10, colour = age_group), width=0.15,
                  size = 1)+
  scale_colour_viridis_d("", begin=0.6, end = 0.8) +
  guides(colour=FALSE) +
  xlab("Age group (years)")+
  scale_y_continuous(name = "HBV carriers (thousands)",
    sec.axis = sec_axis(~./10, name="HBV prevalence (%)")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.y.right = element_text(color = "#1C8C6F"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

# Plot of treatment initiations per assessment
interactions_by_age_rel$with_monitoring <- "No"
interactions_by_age_rel$with_monitoring[interactions_by_age_rel$scenario %in% c("a4_screen_2020_monit_sim7",
                                                                                "a5_screen_2020_monit_sim7",
                                                                                "a2_screen_2020_monit_sim7")] <- "Yes"

pp6 <- ggplot(subset(interactions_by_age_rel, with_monitoring == "No" & interaction_type==
                      "treatment_initiations_per_assessment")) +
  stat_summary(aes(x=age_group, y= value*100, colour=age_group),
               fun="median",
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="pointrange", size = 1)+
  stat_summary(aes(x=age_group, y= value*100, colour=age_group),
               fun="median",
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", width = 0.1, size = 1)+
  scale_colour_viridis_d("", begin=0.6, end = 0.8) +
  guides(colour=FALSE) +
  coord_cartesian(ylim=c(0, 25)) +
  xlab("Age group (years)")+
  ylab("Treatment eligibility (%)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        plot.margin = unit(c(5.5,35,5.5,5.5), "pt"))

# Run pp7 in Cumulative HBV/HCC deaths by age by 2100 to show shift in pattern section


# PANEL PLOT:

programme_plot <- grid.arrange(pp1,pp2,pp3,pp4,ncol=2)

ppa <- arrangeGrob(programme_plot, top = textGrob("A", x = unit(0.01, "npc"),
                                       y   = unit(1, "npc"), just=c("left","top"),
                                       gp=gpar(col="black", fontsize=20)))

ppb <- arrangeGrob(pp5, top = textGrob("B", x = unit(0.1, "npc"),
                                                  y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20)))


ppc <- arrangeGrob(pp6, top = textGrob("C", x = unit(0.1, "npc"),
                                                  y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20)))


ppd <- arrangeGrob(pp7, top = textGrob("D", x = unit(0.1, "npc"),
                                                  y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20)))

expl_plots <- grid.arrange(ppb,ppc,ppd, ncol=1)

#tiff(file = "basic_programme_explanatory_plot.tiff", width=320, height=213, units = "mm", res=300, pointsize = 0.99)
grid.arrange(ppa, expl_plots,
             ncol =2, widths=2:1)
#dev.off()


### Population effects: Treatment effect over time ----

hbv_deaths_over_time <-
  rbind(gather(out2$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario),
    gather(a1_out3_pop$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario),
        gather(a1_out6_pop$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario),
    gather(monit_out7$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario),
    gather(out8b_2030_monit_sim7$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario))

hbv_deaths_over_time <- hbv_deaths_over_time %>%
  group_by(time, scenario) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value,0.975))

hbv_inc_over_time <-
  rbind(gather(out2$timeseries$total_chronic_infections, key="sim", value = "value", -time,-scenario),
        gather(a1_out3_pop$timeseries$total_chronic_infections, key="sim", value = "value", -time,-scenario),
        gather(a1_out6_pop$timeseries$total_chronic_infections, key="sim", value = "value", -time,-scenario),
        gather(monit_out7$timeseries$total_chronic_infections, key="sim", value = "value", -time,-scenario))

hbv_inc_over_time <- hbv_inc_over_time %>%
  group_by(time, scenario) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value,0.975))

# WITH UNCERTAINTY: #

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
  geom_ribbon(data=subset(hbv_deaths_over_time, scenario == "screen_2020_monit_sim7"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, fill = scenario, colour=scenario),
              linetype = "solid",alpha = 0.8) +
  geom_line(data=subset(hbv_deaths_over_time, scenario == "screen_2020_monit_sim7"),
            aes(x=time, y = median/0.5), colour="#2166AC",
            linetype = "solid", size=0.75) +
  scale_fill_manual("Scenario",
                    labels = c("status_quo" = "Infant vaccination only",
                               "screen_2020_monit_sim7" = "Screening & treatment in 2020\n(monitor every 5 years\nin <45 year olds)"),
                    values=c("status_quo" = "#D6604D",
                             "screen_2020_monit_sim7" = "#92C5DE")) +
  scale_colour_manual("Scenario",
                      labels = c("status_quo" = "Infant vaccination only",
                                 "screen_2020_monit_sim7" = "Screening & treatment in 2020\n(monitor every 5 years\nin <45 year olds)"),
                      values=c("status_quo" = "#B2182B",
                             "screen_2020_monit_sim7" = "#2166AC")) +
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

# SQ plots of deaths and incidence (comparison to no historical intervention in sensitivity analysis script)
outcomes_over_time <- rbind(
  data.frame(outcome="HBV-related deaths",
            gather(out2$timeseries$total_hbv_deaths, key="sim", value = "value", -time,-scenario)),
  data.frame(outcome="New chronic HBV infections",
             gather(out2$timeseries$total_chronic_infections, key="sim", value = "value", -time,-scenario))        )

outcomes_over_time <- outcomes_over_time %>%
  group_by(time, outcome) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value,0.975))

ggplot() +
  geom_ribbon(data=outcomes_over_time,
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5),
              linetype = "solid", alpha = 0.2) +
  geom_line(data=outcomes_over_time,
            aes(x=time, y = median/0.5), colour="#B2182B",
            linetype = "solid", size=0.75) +
  facet_wrap(~outcome, scales="free_y", nrow=2) +
  geom_vline(xintercept=1990-0.5, linetype="dashed") +
  ylab("Incident number per year") +
  scale_x_continuous("Year",
                     breaks=c(1990, 2019.5, 2040, 2070, 2090),
                     labels=c(1990, 2020, 2040, 2070, 2090),
                     limits=c(1960,2090)) +
  theme_classic() +
  theme(legend.position=c(.78,.8),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15))

# Thesis plot: Incidence over time and sources of infection ----

# Run previous section for incidence

# Incidence plot:
#png(file = "incidence_timeplot.tiff", width=300, height=125, units = "mm", res=300, pointsize =0.99)
ggplot() +
  geom_ribbon(data=subset(hbv_inc_over_time, scenario == "status_quo"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5,fill = scenario, colour=scenario),
              linetype = "solid", size= 0.75, alpha = 0.25) +
  geom_line(data=subset(hbv_inc_over_time, scenario == "status_quo"),
            aes(x=time, y = median/0.5), colour="#B2182B",
            linetype = "solid", size=1.1) +
  geom_ribbon(data=subset(hbv_inc_over_time, scenario == "screen_2020_monit_0"),
              aes(x=time, ymin=lower/0.5, ymax=upper/0.5, fill = scenario, colour=scenario),
              linetype = "solid", size= 0.75, alpha = 0.2) +
  geom_line(data=subset(hbv_inc_over_time, scenario == "screen_2020_monit_0"),
            aes(x=time, y = median/0.5), colour="#2166AC",
            linetype = "solid", size=1.1) +
  scale_fill_manual("Scenario",
                    labels = c("status_quo" = "Base case",
                               "screen_2020_monit_0" = "Screening and treatment programme"),
                    values=c("status_quo" = "#D6604D",
                             "screen_2020_monit_0" = "#285686")) +
  scale_colour_manual("Scenario",
                      labels = c("status_quo" = "Base case",
                                 "screen_2020_monit_0" = "Screening and treatment programme"),
                      values=c("status_quo" = "#B2182B",
                               "screen_2020_monit_0" = "#285686")) +
  # guides(fill=FALSE, colour = FALSE) +
  #geom_vline(xintercept=2019.5, linetype="dashed") +
  ylab("New chronic HBV infections per year") +
  guides(fill=guide_legend(rev=T),
         colour=guide_legend(rev=T)) +
  scale_x_continuous("Year",
                     breaks=c(2019.5, 2050, 2080),
                     labels=c(2020, 2050, 2080),
                     limits=c(2010,2100)) +
  ylim(0,3250) +
  theme_classic() +
  theme(legend.position=c(.75,.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_blank())
#dev.off()

## SOURCES OF INFECTION
inf_source <- out2 <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/a1_it_screen_2020_monit_0_infection_sources_160321.rds")
inf_source <- inf_source [[1]]

# Rows = time, columns = simulations
timevec <- inf_source[[1]]$time
new_inf_i2 <- as.data.frame(sapply(inf_source, "[[", "new_inf_i2"))
new_inf_i2$time <- timevec
new_inf_i2 <- gather(new_inf_i2, key="sim", value = "i2", -time)

new_inf_i3 <- as.data.frame(sapply(inf_source, "[[", "new_inf_i3"))
new_inf_i3$time <- timevec
new_inf_i3 <- gather(new_inf_i3, key="sim", value = "i3", -time)

new_inf_i4 <- as.data.frame(sapply(inf_source, "[[", "new_inf_i4"))
new_inf_i4$time <- timevec
new_inf_i4 <- gather(new_inf_i4, key="sim", value = "i4", -time)

new_inf_mtct <- as.data.frame(sapply(inf_source, "[[", "mtct"))
new_inf_mtct$time <- timevec
new_inf_mtct <- gather(new_inf_mtct, key="sim", value = "mtct", -time)

new_inf_source <- left_join(left_join(
  left_join(new_inf_i2, new_inf_i3, by = c("time", "sim")),
  new_inf_i4, by = c("time", "sim")),
  new_inf_mtct, by =  c("time", "sim"))
new_inf_source$mtct_i2 <- new_inf_source$mtct+new_inf_source$i2
new_inf_source$mtct_i2_i3 <- new_inf_source$mtct+new_inf_source$i2+new_inf_source$i3
new_inf_source$total <- new_inf_source$mtct+new_inf_source$i2+new_inf_source$i3+
  new_inf_source$i4

new_inf_source_summary <- group_by(new_inf_source, time) %>%
  summarise(median_mtct = median(mtct),
            median_mtct_i2 = median(mtct_i2),
            median_mtct_i2_i3 = median(mtct_i2_i3),
            median_total = median(total))

quantile(new_inf_source$mtct[new_inf_source$time==2020]/
  new_inf_source$total[new_inf_source$time==2020], c(0.5,0.025,0.975))
quantile(new_inf_source$i2[new_inf_source$time==2020]/
           new_inf_source$total[new_inf_source$time==2020], c(0.5,0.025,0.975))
quantile(new_inf_source$i3[new_inf_source$time==2020]/
           new_inf_source$total[new_inf_source$time==2020], c(0.5,0.025,0.975))
quantile(new_inf_source$i4[new_inf_source$time==2020]/
           new_inf_source$total[new_inf_source$time==2020], c(0.5,0.025,0.975))

ggplot(new_inf_source_summary) +
  geom_line(aes(x=time, y = median_mtct)) +
  geom_line(aes(x=time, y = median_mtct_i2)) +
  geom_line(aes(x=time, y = median_mtct_i2_i3)) +
  geom_line(aes(x=time, y = median_total)) +
  theme_classic() +
  xlim(1985,2050)

# Paper plot: HBV deaths and treatment need over time (with repeat screening) ----

# Scenarios are status quo and treatment programme with optimal monitoring (or maybe different options)

# Rerun previous time section for deaths!

# Left side: MEDIAN HBV deaths #

hbv_deaths_over_time$scenario <- factor(hbv_deaths_over_time$scenario)

# V1 (no repeat screening)
timeplot1 <- ggplot(data=subset(hbv_deaths_over_time, scenario != "screen_2020_monit_1" &
                                  scenario != "screen_2020_monit_sim7" &
                                  scenario != "monit_sim7_screen_10b_2030")) +
  geom_line(aes(x=time, y = median/0.5, colour = reorder(scenario,-median)), size = 1) +
  ylab("HBV-related deaths per year") +
  scale_x_continuous("Year",
                     breaks=c(1990, 2019.5, 2050, 2080),
                     labels=c(1990, 2020, 2050, 2080),
                     limits=c(1990,2100)) +
  scale_colour_viridis_d(end=0.6, "",
                         labels =  c("status_quo" = "Base case (no treatment)",
                                     "screen_2020_monit_sim7" = "Treatment programme,\nmonitor 5-yearly in <45 year olds",
                                     "screen_2020_monit_0" = "Treatment programme,\nno monitoring")) +
#  guides(colour=FALSE) +
  theme_classic() +
  theme(legend.position=c(.7,.9),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 13))

# V2 (with repeat screening)
# Get viridis colours for 4 scenarios:
scales::viridis_pal(end=0.9)(length(1:4))
# "#440154FF" "#35608DFF" "#22A884FF" "#BBDF27FF"

timeplot1_v2 <- ggplot(data=subset(hbv_deaths_over_time, scenario != "screen_2020_monit_1")) +
  geom_line(aes(x=time, y = median/0.5, colour = reorder(scenario,-median)), size = 1) +
  ylab("HBV-related deaths per year") +
  scale_x_continuous("Year",
                     breaks=c(1990, 2019.5, 2050, 2080),
                     labels=c(1990, 2020, 2050, 2080),
                     limits=c(1990,2100)) +
  scale_colour_viridis_d(end=0.9, "",
                         labels =  c("status_quo" = "Base case",
                                     "screen_2020_monit_sim7" = "Screening 2020,\nmonitor 5-yearly in <45 year olds",
                                     "screen_2020_monit_0" = "Screening 2020, no monitoring",
                                     "monit_sim7_screen_10b_2030" = "Screening 2020+2030,\nmonitor 5-yearly in <45 year olds")) +
  #  guides(colour=FALSE) +
  theme_classic() +
  theme(legend.position=c(.75,.9),
        axis.text = element_text(size = 14.5),
        axis.title = element_text(size = 14.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))


# Right side: % of unmet treatment need over time
# (line or maybe bar chart).
# % of unmet treatment need among the whole Gambian population indicates there
# will be less treatment need in the future.


##

# TREATMENT NEED PLOT WITH OUT2 BUT NO REPEAT SCREENING AND NO COVERAGE
out_path_full_output <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/"
out2_it_tn <- readRDS(paste0(out_path_full_output, "out2_status_quo_repeat_screens_comparison_280121.rds"))
out2_it_tn <- out2_it_tn[[1]]
out3_it_tn <- readRDS(paste0(out_path_full_output, "a1_it_out3_screen_2020_monit_0_repeat_screens_comparison_280121.rds"))
out3_it_tn <- out3_it_tn[[1]]
monit_out7_it_tn <- readRDS(paste0(out_path_full_output, "a1_it_out3_screen_2020_monit_sim7_repeat_screens_comparison_290121.rds"))
monit_out7_it_tn <- monit_out7_it_tn[[1]]
out8b_it_2030_monit_sim7_tn <- readRDS(paste0(out_path_full_output, "a1_it_out8b_monit_sim7_screen_10b_2030_repeat_screens_comparison_290121.rds"))
out8b_it_2030_monit_sim7_tn  <- out8b_it_2030_monit_sim7_tn[[1]]


unmet_need_df2 <- rbind(
  data.frame(scenario = "screen_2020_monit_0", at_time = "2020",type ="proportion",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2020]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2020])/
               out3_it_tn$total_pop_15_to_65_over_time[out3_it_tn$time==2020]),
  data.frame(scenario = "screen_2020_monit_0", at_time = "2030",type ="proportion",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2030]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2030])/
               out3_it_tn$total_pop_15_to_65_over_time[out3_it_tn$time==2030]),
  data.frame(scenario = "screen_2020_monit_0", at_time = "2040",type ="proportion",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2040]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2040])/
               out3_it_tn$total_pop_15_to_65_over_time[out3_it_tn$time==2040]),
  data.frame(scenario = "screen_2020_monit_0", at_time = "2050",type ="proportion",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2050]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2050])/
               out3_it_tn$total_pop_15_to_65_over_time[out3_it_tn$time==2050]),
  data.frame(scenario = "screen_2020_monit_0", at_time = "2060",type ="proportion",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2060]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2060])/
               out3_it_tn$total_pop_15_to_65_over_time[out3_it_tn$time==2060]),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2020",type ="proportion",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2020]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2020])/
               monit_out7_it_tn$total_pop_15_to_65_over_time[monit_out7_it_tn$time==2020]),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2030",type ="proportion",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2030]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2030])/
               monit_out7_it_tn$total_pop_15_to_65_over_time[monit_out7_it_tn$time==2030]),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2040",type ="proportion",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2040]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2040])/
               monit_out7_it_tn$total_pop_15_to_65_over_time[monit_out7_it_tn$time==2040]),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2050",type ="proportion",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2050]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2050])/
               monit_out7_it_tn$total_pop_15_to_65_over_time[monit_out7_it_tn$time==2050]),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2060",type ="proportion",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2060]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2060])/
               monit_out7_it_tn$total_pop_15_to_65_over_time[monit_out7_it_tn$time==2060]),

    data.frame(scenario = "screen_2020_2030_monit_sim7", at_time = "2020",type ="proportion",
             remaining_treatment_need=(out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2020]+
                                         out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2020])/
               out8b_it_2030_monit_sim7_tn$total_pop_15_to_65_over_time[out8b_it_2030_monit_sim7_tn$time==2020]),
  data.frame(scenario = "screen_2020_2030_monit_sim7", at_time = "2030",type ="proportion",
             remaining_treatment_need=(out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2030]+
                                         out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2030])/
               out8b_it_2030_monit_sim7_tn$total_pop_15_to_65_over_time[out8b_it_2030_monit_sim7_tn$time==2030]),
  data.frame(scenario = "screen_2020_2030_monit_sim7", at_time = "2040",type ="proportion",
             remaining_treatment_need=(out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2040]+
                                         out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2040])/
               out8b_it_2030_monit_sim7_tn$total_pop_15_to_65_over_time[out8b_it_2030_monit_sim7_tn$time==2040]),
  data.frame(scenario = "screen_2020_2030_monit_sim7", at_time = "2050",type ="proportion",
             remaining_treatment_need=(out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2050]+
                                         out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2050])/
               out8b_it_2030_monit_sim7_tn$total_pop_15_to_65_over_time[out8b_it_2030_monit_sim7_tn$time==2050]),
  data.frame(scenario = "screen_2020_2030_monit_sim7", at_time = "2060",type ="proportion",
             remaining_treatment_need=(out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2060]+
                                         out8b_it_2030_monit_sim7_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out8b_it_2030_monit_sim7_tn$time==2060])/
               out8b_it_2030_monit_sim7_tn$total_pop_15_to_65_over_time[out8b_it_2030_monit_sim7_tn$time==2060]),

  data.frame(scenario = "sq", at_time = "2020",type ="proportion",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2020]/
               out2_it_tn$total_pop_15_to_65_over_time[out2_it_tn$time==2020]),
  data.frame(scenario = "sq", at_time = "2030",type ="proportion",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2030]/
               out2_it_tn$total_pop_15_to_65_over_time[out2_it_tn$time==2030]),
  data.frame(scenario = "sq", at_time = "2040",type ="proportion",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2040]/
               out2_it_tn$total_pop_15_to_65_over_time[out2_it_tn$time==2040]),
  data.frame(scenario = "sq", at_time = "2050",type ="proportion",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2050]/
               out2_it_tn$total_pop_15_to_65_over_time[out2_it_tn$time==2050]),
  data.frame(scenario = "sq", at_time = "2060",type ="proportion",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2060]/
               out2_it_tn$total_pop_15_to_65_over_time[out2_it_tn$time==2060])
)

unmet_need_df2$scenario <- factor(unmet_need_df2$scenario,
                                  levels = c("sq", "screen_2020_monit_0", "screen_2020_monit_sim7",
                                             "screen_2020_2030_monit_sim7"))

# V1: no repeat screening
timeplot2 <- ggplot(data=subset(unmet_need_df2,scenario != "screen_2020_monit_sim7" &
                                  scenario != "screen_2020_2030_monit_sim7")) +
  stat_summary(aes(x=reorder(scenario, -remaining_treatment_need),
                   y = remaining_treatment_need*100,
                   fill = scenario), colour="black",
               fun="median", geom="bar") +
  facet_wrap(~at_time, ncol = 5, strip.position="bottom") +
  scale_fill_viridis_d(end=0.6, "Scenario",
                       labels =  c("sq" = "Base case (no treatment)",
                                   "screen_2020_monit_sim7" = "Treatment programme,\nmonitor 5-yearly in <45 year olds",
                                   "screen_2020_monit_0" = "Treatment programme,\nno monitoring")) +
  guides(fill=FALSE) +
  ylab("Unmet treatment need in total population (%)") +
  xlab("Year") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = c(0.75, 0.75),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

# V2: with repeat screening
timeplot2_v2 <- ggplot(data=subset(unmet_need_df2)) +
  stat_summary(aes(x=reorder(scenario, -remaining_treatment_need),
                   y = remaining_treatment_need*100,
                   fill = scenario), colour="black",
               fun="median", geom="bar") +
  facet_wrap(~at_time, ncol = 5, strip.position="bottom") +
  scale_fill_viridis_d(end=0.8, "Scenario",
                       labels =  c("sq" = "Base case",
                                   "screen_2020_monit_sim7" = "Screening 2020,\nmonitor 5-yearly in <45 year olds",
                                   "screen_2020_2030_monit_sim7" = "Screening 2020+2030,\nmonitor 5-yearly in <45 year olds",
                                   "screen_2020_monit_0" = "Screening 2020,\nno monitoring")) +
  guides(fill=FALSE) +
  ylab("Unmet treatment need in total population (%)") +
  xlab("Year") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = c(0.75, 0.75),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 14.5),
        axis.title = element_text(size = 14.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))


# Combined plot:
grid.arrange(timeplot1, timeplot2, ncol = 2, widths = 3:2)

# PAPER PLOT:
#tiff(file = "timeplot.tiff", width=300, height=125, units = "mm", res=300, pointsize = 0.99)
grid.arrange(timeplot1_v2, timeplot2_v2, ncol = 2, widths = 3:2)
#dev.off()

# THESIS PLOT (combined with plot from repeat_screening_analysis.R):
#png(file = "scenarios_timeplot.tiff", width=300/1.2, height=125/1.2, units = "mm", res=300, pointsize =0.99)
print(timeplot1_v2)
#dev.off()

# Number instead of %
unmet_need_df_number <- rbind(
  data.frame(scenario = "screen_2020_monit_0", at_time = "2020",type ="number",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2020]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2020])),
  data.frame(scenario = "screen_2020_monit_0", at_time = "2030",type ="number",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2030]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2030])),
  data.frame(scenario = "screen_2020_monit_0", at_time = "2040",type ="number",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2040]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2040])),
  data.frame(scenario = "screen_2020_monit_0", at_time = "2050",type ="number",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2050]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2050])),
  data.frame(scenario = "screen_2020_monit_0", at_time = "2060",type ="number",
             remaining_treatment_need=(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out3_it_tn$time==2060]+
                                         out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[out3_it_tn$time==2060])),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2020",type ="number",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2020]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2020])),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2030",type ="number",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2030]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2030])),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2040",type ="number",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2040]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2040])),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2050",type ="number",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2050]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2050])),
  data.frame(scenario = "screen_2020_monit_sim7", at_time = "2060",type ="number",
             remaining_treatment_need=(monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[monit_out7_it_tn$time==2060]+
                                         monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65[monit_out7_it_tn$time==2060])),
  data.frame(scenario = "sq", at_time = "2020",type ="number",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2020]),
  data.frame(scenario = "sq", at_time = "2030",type ="number",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2030]),
  data.frame(scenario = "sq", at_time = "2040",type ="number",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2040]),
  data.frame(scenario = "sq", at_time = "2050",type ="number",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2050]),
  data.frame(scenario = "sq", at_time = "2060",type ="number",
             remaining_treatment_need=out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65[out2_it_tn$time==2060])
)

unmet_need_df_number$scenario <- factor(unmet_need_df_number$scenario,
                                  levels = c("sq", "screen_2020_monit_0", "screen_2020_monit_sim7"))

ggplot(data=subset(unmet_need_df_number,scenario != "screen_2020_monit_sim7")) +
  stat_summary(aes(x=reorder(scenario, -remaining_treatment_need),
                   y = remaining_treatment_need*100,
                   fill = scenario), colour="black",
               fun="median", geom="bar") +
  facet_wrap(~at_time, ncol = 5, strip.position="bottom") +
  scale_fill_viridis_d(end=0.6, "Scenario",
                       labels =  c("sq" = "Base case (no treatment)",
                                   "screen_2020_monit_sim7" = "Treatment programme,\nmonitor 5-yearly in <45 year olds",
                                   "screen_2020_monit_0" = "Treatment programme,\nno monitoring")) +
  guides(fill=FALSE) +
  ylab("Unmet treatment need in total population (%)") +
  xlab("Year") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = c(0.75, 0.75),
        strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))


test <- cbind(subset(unmet_need_df_number, scenario == "sq"),
              with_treatment = subset(unmet_need_df_number, scenario == "screen_2020_monit_0")$remaining_treatment_need)

# Number on treatment over time (those removed from unmet treatment need by the treatment)
# divided by all treatment eligible in out2
plot(x=out2_it_tn$time,
       y=apply((out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65-
(out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65+
  out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65))/
  out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65,1,median)*100)

# Percent of all treatment eligible carriers by 2100 identified by initial screen
quantile(apply((out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65-
         (out3_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65+
            out3_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65))/
        out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65,2,sum),
        c(0.5,0.025,0.975))

quantile(apply((out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65-
                  (monit_out7_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65+
                     monit_out7_it_tn$treatment_eligible_carriers_screened_over_time_15_to_65))/
                 out2_it_tn$treatment_eligible_carriers_undiagnosed_over_time_15_to_65,2,sum),
         c(0.5,0.025,0.975))



### Diagnostic plots for screening and monitoring by age, or monitoring frequency ----
# Plots of breakdown of impact of 5-yearly monitoring and screening by separate age group (previously in monitoring frequency file) ----


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

# All outcomes and interactions are incremental to the no monitoring scenario.

# DALYs averted
p <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/1000), fill = "grey", lwd=1) +
  ylab("DALYs averted\n(thousands)") +
  xlab("Monitored age group") +
  ylim(0,61) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

px <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/1000), fill = "grey", lwd=1) +
  ylab("HBV deaths averted\n(thousands)") +
  xlab("Monitored age group") +
  ylim(0,4.1) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

p0 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = monitoring_assessments/1000), fill = "grey", lwd=1) +
  ylab("Monitoring assessments\n(thousands)") +
  xlab("Monitored age group") +
  ylim(0,400)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental DALYs averted per monitoring assessment
p1 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/monitoring_assessments), fill = "grey", lwd=1) +
  ylab("DALYs averted\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,1)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental DALYs averted per treatment initiations
p2 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/treatment_initiations), fill = "grey", lwd=1) +
  ylab("DALYs averted\nper treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,30)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental DALYs averted per PY on treatment
p3 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = dalys_averted/py_on_treatment), fill = "grey", lwd=1) +
  ylab("DALYs averted\nper PY on treatment") +
  xlab("Monitored age group") +
  ylim(0,0.75)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental deaths averted per monitoring assessment
p4 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/monitoring_assessments), fill = "grey", lwd=1) +
  ylab("HBV deaths averted\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,0.04)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental deaths  averted per treatment initiations
p5 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/treatment_initiations), fill = "grey", lwd=1) +
  ylab("HBV deaths averted\nper treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,0.7)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental deaths  averted per PY on treatment
p6 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = deaths_averted/py_on_treatment), fill = "grey", lwd=1) +
  ylab("HBV deaths averted\nper PY on treatment") +
  xlab("Monitored age group") +
  ylim(0,0.05)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental treatment initiations per monitoring assessment
p7 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = treatment_initiations/monitoring_assessments), fill = "grey", lwd=1) +
  ylab("Treatment initiations\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,0.07)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental PY on treatment per monitoring assessment
p8 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/monitoring_assessments), fill = "grey", lwd=1) +
  ylab("PY on treatment\nper monitoring assessment") +
  xlab("Monitored age group") +
  ylim(0,1.75)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

# Incremental PY on treatment per treatment initiation
p9 <- ggplot(df_by_age_group) +
  geom_boxplot(aes(x=scenario, y = py_on_treatment/treatment_initiations), fill = "grey", lwd=1) +
  ylab("PY on treatment\nper treatment initiation") +
  xlab("Monitored age group") +
  ylim(0,35)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 15))

grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 3)
grid.arrange(p7,p8, ncol = 1)

# THESIS PLOT
# Monitor
library(grid)
#png(file = "monitoring_by_age_plot.png", width=360, height=220, units = "mm", res=300, pointsize = 0.99)
grid.arrange(p,px,p0,p7,p9,
             p1,p2,p3,
             p4,p5,p6,
             layout_matrix = rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5),
                                   c(6,6,6,6,6,7,7,7,7,7,8,8,8,8,8),
                                   c(9,9,9,9,9,10,10,10,10,10,11,11,11,11,11)),
             bottom = textGrob("Monitored age group", gp=gpar(fontsize=15)))
#dev.off()

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

# Theoretical monitoring plots (cohort): monitoring frequency, diminishing returns and reasons why ----

# Use simulations: every 20, 10, 5, 2 and 1 years
# Sims = out4b_it, out4_it, out5_it, out6a_it, a1_out6_pop
# Comparator is NO MONITORING: out3_it
# Maybe change this to doubling frequencies of 16, 8, 4, 2 and 1 years

# Also look at sensitivity analysis for this!

# Plot of diminishing returns: % of DALYS averted in the cohort
monitoring_prop_dalys_averted <-
  plot_hbv_deaths_averted_cohort(counterfactual_object = out3_it,
                                 scenario_objects = list(out4b_it, out4_it, out5_it,
                                                          out6a_it, a1_out6_pop),
                                 outcome_to_avert = "cohort_dalys",
                                 outcome_to_plot = "number_averted",
                                 counterfactual_label = "no monitoring")

monitoring_prop_dalys_averted_summary <- subset(monitoring_prop_dalys_averted,
                                                type == "proportion_averted") %>%
  group_by(scenario) %>%
  summarise(median = median(value*100),
            lower  = quantile(value*100, 0.025),
            upper = quantile(value*100, 0.975))
monitoring_prop_dalys_averted_summary$outcome <- "percentage_dalys_averted"


# Interactions compared to no monitoring
interactions_df <-rbind(
  data.frame(scenario = "screen_2020_monit_20",
             sim = colnames(out3_it$cohort_dalys[,-1]),
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-out4b_it$cohort_dalys[,-1]),
             treatment_initiations = unlist(out4b_it$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(out4b_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_10",
             sim = colnames(out3_it$cohort_dalys[,-1]),
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-out4_it$cohort_dalys[,-1]),
             treatment_initiations = unlist(out4_it$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(out4_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             sim = colnames(out3_it$cohort_dalys[,-1]),
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-out5_it$cohort_dalys[,-1]),
             treatment_initiations = unlist(out5_it$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(out5_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_2",
             sim = colnames(out3_it$cohort_dalys[,-1]),
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-out6a_it$cohort_dalys[,-1]),
             treatment_initiations = unlist(out6a_it$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(out6a_it$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             sim = colnames(out3_it$cohort_dalys[,-1]),
             cohort_dalys_averted = unlist(out3_it$cohort_dalys[,-1]-a1_out6_pop$cohort_dalys[,-1]),
             treatment_initiations = unlist(a1_out6_pop$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_treated[,-c(1:3)]),
             monitoring_assessments = unlist(a1_out6_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(out3_it$interactions[[16]]$total_assessed[,-c(1:3)]))
)

interactions_df_summary <- group_by(interactions_df, scenario) %>%
  summarise(median = median(treatment_initiations/monitoring_assessments*1000),
            lower = quantile(treatment_initiations/monitoring_assessments*1000, 0.025),
            upper = quantile(treatment_initiations/monitoring_assessments*1000, 0.975))
interactions_df_summary$outcome <- "treatment_initiations_per_1000_monitoring_assessments"

monitoring_prop_dalys_averted_summary <- rbind(monitoring_prop_dalys_averted_summary,
                                               interactions_df_summary, interactions_df_summary2)
monitoring_prop_dalys_averted_summary$scenario <- factor(monitoring_prop_dalys_averted_summary$scenario,
                                                         levels = c("screen_2020_monit_20",
                                                                    "screen_2020_monit_10",
                                                                    "screen_2020_monit_5",
                                                                    "screen_2020_monit_2",
                                                                    "screen_2020_monit_1"))

# Plot doesn't really work with both cause the outcome units are not the same (e.g. percentage)
ggplot(subset(monitoring_prop_dalys_averted_summary, outcome == "percentage_dalys_averted")) +
  geom_point(aes(x=reorder(scenario, desc(scenario)), y = median), size = 3,
             position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=reorder(scenario, desc(scenario)), ymin = lower, ymax = upper),
                width = 0.15,  position=position_dodge(width=0.9)) +
  theme_classic() +
  guides(colour=FALSE) +
  scale_x_discrete(labels=c("screen_2020_monit_20" = "20 years",
                            "screen_2020_monit_10" = "10 years",
                            "screen_2020_monit_5" = "5 years",
                            "screen_2020_monit_2" = "2 years",
                            "screen_2020_monit_1" = "1 year")) +
  xlab("Monitoring\ninterval") +
  ylim(0,80) +
  ylab("DALYs averted (%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  coord_flip()

dalys_against_assessments <- subset(monitoring_prop_dalys_averted_summary, outcome == "percentage_dalys_averted")
dalys_against_assessments$median_monitoring_assessments <- data.frame(group_by(interactions_df, scenario) %>%
  summarise(median = median(monitoring_assessments)))[,"median"]

# THESIS PLOT

#png(file = "monitoring_freq_plot.png", width=300, height=80, units = "mm", res=300, pointsize = 0.99)
ggplot(dalys_against_assessments) +
  geom_point(aes(x= median_monitoring_assessments/1000,
                 y = median, colour = scenario), size = 4) +
  geom_errorbar(aes(x=median_monitoring_assessments/1000,
                    ymin = lower, ymax = upper, colour = scenario), width = 20) +
  theme_classic() +
  #guides(colour=FALSE) +
  scale_colour_viridis_d("Monitoring\ninterval",
                        labels=c("screen_2020_monit_20" = "20 years",
                            "screen_2020_monit_10" = "10 years",
                            "screen_2020_monit_5" = "5 years",
                            "screen_2020_monit_2" = "2 years",
                            "screen_2020_monit_1" = "1 year")) +
  xlab("Monitoring assessments (thousands)") +
  ylim(0,80) +
  xlim(0,1500)+
  ylab("DALYs averted (%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
#dev.off()


# Hypothesis 1: incremental new treatment initiations per incremental assessments
# decline with increasing frequency
# (= monitoring more often does not identify more people to treat)

# Compare: every 20 years vs no monitoring, every 10 vs every 20, every 5 vs every 10,
# every 2 vs every 5, every 1 vs every 2

# out3_it, a1_out6_pop, out4b_it, out4_it, out5_it, out6a_it

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

# From 20 yearly monitoring onwards, the ratio of new treatment initiations per
# incremental monitoring assessment decreases
# The incremental person-time on treatment per assessment decreases from 20-yearly

median_values <- group_by(interactions_df, scenario) %>%
  summarise(monitoring_assessments=median(monitoring_assessments),
            treatment_initiations=median(treatment_initiations),
            cohort_dalys_averted = median(cohort_dalys_averted))

monit_exp_p1 <- ggplot(interactions_df) +
  geom_point(aes(x=monitoring_assessments/1000, y = treatment_initiations/1000,
                 colour = reorder(scenario, monitoring_assessments))) +
  geom_point(data=median_values, aes(x=monitoring_assessments/1000,
                                     y = treatment_initiations/1000,
                                     group = reorder(scenario, monitoring_assessments),
                 colour = reorder(scenario, monitoring_assessments)), size = 5) +
  geom_point(data=median_values, aes(x=monitoring_assessments/1000,
                                     y = treatment_initiations/1000,
                                     group = reorder(scenario, monitoring_assessments)),
             size = 5, shape = 21, colour="black") +
 scale_colour_viridis_d("Average monitoring interval",
                       labels=c("screen_2020_monit_20" = "20 years",
                                "screen_2020_monit_10" = "10 years",
                                "screen_2020_monit_5" = "5 years",
                                "screen_2020_monit_2" = "2 years",
                                "screen_2020_monit_1" = "1 year")) +
  ylab("Treatment initiations (thousands)") +
  xlab("Monitoring assessments (thousands)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
                axis.title = element_text(size = 15),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 15))




# Hypothesis 2: The incremental impact per additional treatment initiations declines with
# increasing frequency
# (= the disease progression identified by more monitoring is less affected by treatment,
# e.g. because cirrhotic cases are treated immediately).

ggplot(incremental_df) +
  geom_boxplot(aes(x=scenario, y = dalys_averted_per_treatment_initiations))

ggplot(incremental_df) +
  geom_boxplot(aes(x=scenario, y = dalys_averted_per_py_on_treatment))

# Other visualisation:
monit_exp_p2 <- ggplot(interactions_df) +
  geom_point(aes(x=treatment_initiations/1000, y = cohort_dalys_averted/1000,
                 colour = reorder(scenario, treatment_initiations))) +
  geom_point(data=median_values,
             aes(x = treatment_initiations/1000,
                 y = cohort_dalys_averted/1000,
                 group = reorder(scenario, treatment_initiations),
                 colour = reorder(scenario, treatment_initiations)), size = 5) +
  geom_point(data=median_values,
             aes(x = treatment_initiations/1000,
                 y = cohort_dalys_averted/1000,
                 group = reorder(scenario, treatment_initiations),
                 colour = reorder(scenario, treatment_initiations)), size = 5,
             shape=21, col = "black") +
  scale_colour_viridis_d("Average monitoring interval",
                        labels=c("screen_2020_monit_20" = "20 years",
                                 "screen_2020_monit_10" = "10 years",
                                 "screen_2020_monit_5" = "5 years",
                                 "screen_2020_monit_2" = "2 years",
                                 "screen_2020_monit_1" = "1 year")) +
  xlab("Treatment initiations (thousands)") +
  ylab("DALYs averted (thousands)") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

# The ratio of incremental DALYS averted per incremental time on treatment  declines
# with increasing monitoring frequencies.

# THESIS PLOT:
library(ggpubr)
#png(file = "monitoring_freq_explanatory_plot.png", width=300, height=150, units = "mm", res=300, pointsize = 0.99)
ggarrange(monit_exp_p1, monit_exp_p2,nrow=1, common.legend = TRUE, legend="bottom")
#dev.off()

# Base case stats ----
quantile(out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2020,-c(1,2)]+
  out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2020.5,-c(1,2)],
  c(0.5,0.025,0.975))

quantile(out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2020,-c(1,2)]+
           out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2020.5,-c(1,2)],
         c(0.5,0.025,0.975))

quantile(out2$timeseries$prev[out2$timeseries$prev$time==2020,-c(1,2)],
         c(0.5,0.025,0.975))

# Deaths and DALYS occurring between 2020 and 2100
quantile(out2$cum_hbv_deaths[[16]][,-c(1:3)],c(0.5,0.025,0.975))
quantile(out2$dalys[[16]][,-c(1:3)],c(0.5,0.025,0.975))

# Deaths and DALYS averted between 2020 and 2100 without monitoring
quantile((out2$cum_hbv_deaths[[16]][,-c(1:3)]-out3_it$cum_hbv_deaths[[16]][,-c(1:3)])/
           out2$cum_hbv_deaths[[16]][,-c(1:3)],c(0.5,0.025,0.975))
quantile((out2$dalys[[16]][,-c(1:3)]-out3_it$dalys[[16]][,-c(1:3)])/
           out2$dalys[[16]][,-c(1:3)],c(0.5,0.025,0.975))

quantile(out3_it$interactions[[16]]$total_screened[,-c(1:3)],c(0.5,0.025,0.975))
quantile(out3_it$interactions[[16]]$total_treated[,-c(1:3)],c(0.5,0.025,0.975))

quantile((out2$dalys[[16]][,-c(1:3)]-out3_it$dalys[[16]][,-c(1:3)])/
  (out2$dalys[[16]][,-c(1:3)]-a1_out6_pop$dalys[[16]][,-c(1:3)]), c(0.5,0.025,0.975))

quantile((out1_it$cohort_dalys[,-1]-out3_it$cohort_dalys[,-1])/
           (out1_it$cohort_dalys[,-1]-a1_out6_pop$cohort_dalys[,-1]), c(0.5,0.025,0.975))

# Cohort outcomes table ----

cohort_dalys <- rbind(out1_it$cohort_dalys,
                      out3_it$cohort_dalys,
                      out4_it$cohort_dalys,
                      out5_it$cohort_dalys,
                      out6a_it$cohort_dalys,
                      a1_out6_pop$cohort_dalys)
cohort_dalys <- gather(cohort_dalys, key = "sim", value="dalys",-scenario)

cohort_size <- rbind(out1_it$cohort_size,
                      out3_it$cohort_size,
                      out4_it$cohort_size,
                      out5_it$cohort_size,
                      out6a_it$cohort_size,
                      a1_out6_pop$cohort_size)
cohort_size <- gather(cohort_size, key = "sim", value="cohort_size",-scenario)

cohort_hbv_deaths <- rbind(out1_it$cohort_cum_hbv_deaths,
                      out3_it$cohort_cum_hbv_deaths,
                      out4_it$cohort_cum_hbv_deaths,
                      out5_it$cohort_cum_hbv_deaths,
                      out6a_it$cohort_cum_hbv_deaths,
                      a1_out6_pop$cohort_cum_hbv_deaths)
cohort_hbv_deaths <- gather(cohort_hbv_deaths, key = "sim", value="deaths",-scenario)

cohort_age_at_death <-rbind(out1_it$cohort_age_at_death,
                            out3_it$cohort_age_at_death,
                            out4_it$cohort_age_at_death,
                            out5_it$cohort_age_at_death,
                            out6a_it$cohort_age_at_death,
                            a1_out6_pop$cohort_age_at_death)
cohort_age_at_death <- gather(cohort_age_at_death, key = "sim", value="age_at_death",-scenario)

cohort_outcomes <- cohort_dalys
cohort_outcomes$dalys_per_person <- cohort_outcomes$dalys/cohort_size$cohort_size
cohort_outcomes$hbv_deaths <- cohort_hbv_deaths$deaths
cohort_outcomes$hbv_death_risk <- cohort_hbv_deaths$deaths/cohort_size$cohort_size
cohort_outcomes$average_age_at_death <- cohort_age_at_death$age_at_death

quantile(out3_it$cohort_age_at_death[,-1]-out1_it$cohort_age_at_death[,-1],
         c(0.5,0.025,0.975))
quantile(a1_out6_pop$cohort_age_at_death[,-1]-out1_it$cohort_age_at_death[,-1],
         c(0.5,0.025,0.975))


View(cohort_outcomes %>% group_by(scenario) %>%
  summarise(median_dalys = round(median(dalys),0),
            lower_dalys = round(quantile(dalys, 0.025),0),
            upper_dalys = round(quantile(dalys, 0.975),0),
            median_dalys_pp = round(median(dalys_per_person),1),
            lower_dalys_pp = round(quantile(dalys_per_person, 0.025),1),
            upper_dalys_pp = round(quantile(dalys_per_person, 0.975),1),
            median_hbv_deaths = round(median(hbv_deaths),0),
            lower_hbv_deaths = round(quantile(hbv_deaths , 0.025),0),
            upper_hbv_deaths = round(quantile(hbv_deaths , 0.975),0),
            median_hbv_death_risk = round(median(hbv_death_risk*100),1),
            lower_hbv_death_risk = round(quantile(hbv_death_risk*100, 0.025),1),
            upper_hbv_death_risk = round(quantile(hbv_death_risk*100, 0.975),1),
            median_average_age_at_death = round(median(average_age_at_death),1),
            lower_average_age_at_death = round(quantile(average_age_at_death, 0.025),1),
            upper_average_age_at_death = round(quantile(average_age_at_death, 0.975),1)))




# Testplot : plot total DALYs with base case
cohort_dalys_summary <- cohort_dalys %>%
  group_by(scenario) %>%
  summarise(median = median(value),
            lower  = quantile(value, 0.025),
            upper = quantile(value, 0.975))

ggplot(subset(cohort_dalys_summary, scenario != "status_quo_cohort")) +
  geom_point(aes(x=reorder(scenario,median), y = median), size = 3,
             position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=reorder(scenario, median), ymin = lower, ymax = upper),
                width = 0.15,  position=position_dodge(width=0.9)) +
  theme_classic() +
  guides(colour=FALSE) +
  scale_x_discrete(labels=c("screen_2020_monit_20" = "20 years",
                            "screen_2020_monit_10" = "10 years",
                            "screen_2020_monit_5" = "5 years",
                            "screen_2020_monit_2" = "2 years",
                            "screen_2020_monit_1" = "1 year")) +
  xlab("Monitoring\ninterval") +
  ylim(0,175000) +
  #ylab("DALYs averted (%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  coord_flip()
