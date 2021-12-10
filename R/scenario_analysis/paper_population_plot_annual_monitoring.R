# Paper panel plot of effect for annual monitoring scenario
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

# Screening (no monitoring) in separate age groups
# A4 = 15-30
# A5 = 30-45
# A2 = 45-65

## Different output:

# All simulations involve IT treatment
a4_out6 <- readRDS(paste0(out_path, "a4_it_out6_screen_2020_monit_1_treatment_effect_091221.rds"))
a4_out6 <- a4_out6[[1]]
a5_out6 <- readRDS(paste0(out_path, "a5_it_out6_screen_2020_monit_1_treatment_effect_091221.rds"))
a5_out6 <- a5_out6[[1]]
a2_out6 <- readRDS(paste0(out_path, "a2_it_out6_screen_2020_monit_1_treatment_effect_101221.rds"))
a2_out6 <- a2_out6[[1]]

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
# Simulations by age group screened, with annual monitoring:
a2_out6_pop <- readRDS(paste0(out_path_monit, "a2_it_out6_screen_2020_monit_1_091221.rds"))
a2_out6_pop <- a2_out6_pop[[1]]
a4_out6_pop <- readRDS(paste0(out_path_monit, "a4_it_out6_screen_2020_monit_1_091221.rds"))
a4_out6_pop <- a4_out6_pop[[1]]
a5_out6_pop <- readRDS(paste0(out_path_monit, "a5_it_out6_screen_2020_monit_1_091221.rds"))
a5_out6_pop <- a5_out6_pop[[1]]

# Need to run: functions in monitoring_frequency_by_age_incremental_analysis.r script ----
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
a1_x2 <- a4_out6$cum_screened_hbv_deaths_male_by_age_2100 +
  a4_out6$cum_screened_hbv_deaths_female_by_age_2100+
  a5_out6$cum_screened_hbv_deaths_male_by_age_2100 +
  a5_out6$cum_screened_hbv_deaths_female_by_age_2100+
  a2_out6$cum_screened_hbv_deaths_male_by_age_2100 +
  a2_out6$cum_screened_hbv_deaths_female_by_age_2100
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
        plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"))
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


# Population effects: plot with resource utilisation ----
# Number of deaths averted in each age group
deaths_averted_by_age <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a4_out3_pop,
                                                  a5_out3_pop,
                                                  a2_out3_pop,
                                                  a4_out6_pop,
                                                  a5_out6_pop,
                                                  a2_out6_pop),
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
             filter(deaths_averted_by_age, scenario == "a4_screen_2020_monit_1" &
                      by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age, scenario == "a5_screen_2020_monit_1" &
                      by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(deaths_averted_by_age, scenario == "a2_screen_2020_monit_1" &
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
             filter(deaths_averted_by_age, scenario == "a4_screen_2020_monit_1" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age, scenario == "a5_screen_2020_monit_1" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(deaths_averted_by_age, scenario == "a2_screen_2020_monit_1" &
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
                                                  a4_out6_pop,
                                                  a5_out6_pop,
                                                  a2_out6_pop),
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
             filter(dalys_averted_by_age, scenario == "a4_screen_2020_monit_1" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(dalys_averted_by_age, scenario == "a5_screen_2020_monit_1" &
                      by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-65",
             filter(dalys_averted_by_age, scenario == "a2_screen_2020_monit_1" &
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
             treatment_initiations_due_to_assessment = unlist(a4_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_monitoring = 0,
             py_on_treatment = unlist(a4_out3_pop$py_on_treatment[[16]])),
  data.frame(age_group = "30-45",
             scenario = "a5_screen_2020_monit_0",
             sim = names(a5_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a5_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a5_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_assessment = unlist(a5_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_monitoring = 0,
             py_on_treatment = unlist(a5_out3_pop$py_on_treatment[[16]])),
  data.frame(age_group = "45-65",
             scenario = "a2_screen_2020_monit_0",
             sim = names(a2_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a2_out3_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a2_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_assessment = unlist(a2_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_monitoring = 0,
             py_on_treatment = unlist(a2_out3_pop$py_on_treatment[[16]])),
  data.frame(age_group = "15-30",
             scenario = "a4_screen_2020_monit_1",
             sim = names(a4_out6_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a4_out6_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a4_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a4_out6_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a4_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a4_out6_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_assessment = unlist(a4_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_monitoring = unlist(a4_out6_pop$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(a4_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a4_out6_pop$py_on_treatment[[16]])),
  data.frame(age_group = "30-45",
             scenario = "a5_screen_2020_monit_1",
             sim = names(a5_out6_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a5_out6_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a5_out6_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a5_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a5_out6_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_assessment = unlist(a5_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_monitoring = unlist(a5_out6_pop$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(a5_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a5_out6_pop$py_on_treatment[[16]])),
  data.frame(age_group = "45-65",
             scenario = "a2_screen_2020_monit_1",
             sim = names(a2_out6_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a2_out6_pop$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             monitoring_assessments =  unlist(a2_out6_pop$interactions[[16]]$total_assessed[,-c(1:3)])-
               unlist(a2_out3_pop$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a2_out6_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_assessment = unlist(a2_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             treatment_initiations_due_to_monitoring = unlist(a2_out6_pop$interactions[[16]]$total_treated[,-c(1:3)])-
               unlist(a2_out3_pop$interactions[[16]]$total_treated[,-c(1:3)]),
             py_on_treatment = unlist(a2_out6_pop$py_on_treatment[[16]]))
)



interactions_by_age <- interactions_by_age %>%
  mutate(total_interactions = hbsag_tests + clinical_assessments + monitoring_assessments + treatment_initiations,
         treatment_initiations_per_assessment = treatment_initiations/(clinical_assessments+monitoring_assessments),
         treatment_initiations_per_clinical_assessment = treatment_initiations_due_to_assessment/clinical_assessments,
         treatment_initiations_per_monitoring_assessment = treatment_initiations_due_to_monitoring/monitoring_assessments,
         py_on_treatment_per_initiation = py_on_treatment/treatment_initiations)
interactions_by_age <- gather(interactions_by_age, key = "interaction_type", value = "value",
                              -age_group, -scenario, - sim)

interactions_by_age_rel <- subset(interactions_by_age, interaction_type %in%
                                    c("treatment_initiations_per_assessment",
                                      "treatment_initiations_per_clinical_assessment",
                                      "treatment_initiations_per_monitoring_assessment",
                                      "py_on_treatment_per_initiation"))
interactions_by_age <- subset(interactions_by_age, !(interaction_type %in%
                                                       c("treatment_initiations_per_assessment",
                                                         "treatment_initiations_per_clinical_assessment",
                                                         "treatment_initiations_per_monitoring_assessment",
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


# Same plot with DISCOUNTED costs:

object_list1 <- list(a4_out3_pop,
                     a4_out6_pop)
object_list2 <- list(a5_out3_pop,
                     a5_out6_pop)
object_list3 <- list(a2_out3_pop,
                     a2_out6_pop)

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
                                    c("a5_screen_2020_monit_0","a5_screen_2020_monit_1")] <- "30-45"
discounted_interactions$age_group[discounted_interactions$scenario %in%
                                    c("a2_screen_2020_monit_0","a2_screen_2020_monit_1")] <- "45-65"
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

# Paper panel plot of no monitoring programme ----

# PAPER PLOT / THESIS PLOT

# Need to run: functions in monitoring_frequency_by_age_incremental_analysis.r script,
# previous 2 sections (Population effects, Distribution of resource utilisation),
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

# Alternative: distinguish effect of monitoring
# Errorbars relate to total (with monitoring)
outcomes_by_age$with_monitoring <- "Yes"
outcomes_by_age$with_monitoring[outcomes_by_age$scenario %in% c("a2_screen_2020_monit_0",
                                                                "a4_screen_2020_monit_0",
                                                                "a5_screen_2020_monit_0")] <- "No"
outcomes_by_age$with_monitoring <- factor(outcomes_by_age$with_monitoring, levels=c("Yes", "No"))

pp1 <- ggplot(data= subset(outcomes_by_age,
                      outcome == "deaths_averted")) +
  stat_summary(aes(x=age_group, y= value/1000, fill = with_monitoring),
               fun="median", geom="bar", position = "identity", width = 0.95, colour="black")+
  stat_summary(data=subset(outcomes_by_age, with_monitoring=="Yes" &
                      outcome == "deaths_averted"),
               aes(x=age_group, y= value/1000),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", position = "dodge2", width = 0.15)+
  scale_fill_manual(values=c("No" = "#35B779FF",
                             "Yes" = "#31688EFF")) +
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

#scales::viridis_pal()(4)
#green = "#35B779FF",
#purple = "#440154FF"
#blue="#31688EFF"
#yellow="#FDE725FF"

pp2 <- ggplot(data= subset(outcomes_by_age, scenario %in% c("a2_screen_2020_monit_1",
                                                            "a4_screen_2020_monit_1",
                                                            "a5_screen_2020_monit_1") &
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

# New plot:
pp2 <- ggplot(data= subset(outcomes_by_age,
                      outcome == "dalys_averted")) +
  stat_summary(aes(x=age_group, y= value/1000, fill = with_monitoring),
               fun="median", geom="bar", position = "identity", width = 0.95, colour="black")+
  stat_summary(data=subset(outcomes_by_age, with_monitoring=="Yes" &
                             outcome == "dalys_averted"),
               aes(x=age_group, y= value/1000),
               fun.min = function(z) {quantile(z,0.025)},
               fun.max = function(z) {quantile(z,0.975)},
               geom="errorbar", position = "dodge2", width = 0.15)+
  scale_fill_manual(values=c("No" = "#35B779FF",
                             "Yes" = "#31688EFF")) +
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
  geom_col(data=subset(interactions_median, scenario %in% c("a2_screen_2020_monit_1",
                                                            "a4_screen_2020_monit_1",
                                                            "a5_screen_2020_monit_1") &
                         !(interaction_type %in% c("total_interactions", "py_on_treatment"))),
           aes(x=age_group, y= value/1000000, fill = interaction_type),
           width = 0.95, colour = "black") +
  #geom_errorbar(data=subset(total_interactions_errorbar,scenario %in% c("a2_screen_2020_monit_sim7",
  #                                                                        "a4_screen_2020_monit_sim7",
  #                                                                        "a5_screen_2020_monit_sim7")),
  #                aes(x = age_group, ymin = lower/1000, ymax = upper/1000), width = 0.15) +
  labs(fill = "Resource utilisation") +
 # scale_fill_viridis_d(labels = c("hbsag_tests" = "Screening",
 #                                 "clinical_assessments" = "Assessment",
 #                                 "treatment_initiations" = "Treatment",
 #                                 "monitoring_assessments" = "Monitoring")) +
   scale_fill_manual(labels = c("hbsag_tests" = "Screening",
                                   "clinical_assessments" = "Assessment",
                                   "treatment_initiations" = "Treatment",
                                   "monitoring_assessments" = "Monitoring"),
                     values= c("hbsag_tests" = "#FDE725FF",
                               "clinical_assessments" = "#35B779FF",
                               "treatment_initiations" = "#440154FF",
                               "monitoring_assessments" = "#31688EFF")) +
   #scale_linetype_manual(values = c("HBV carriers" = "dashed")) +
  scale_y_continuous(breaks=c(0,1,2), limits=c(0,2)) +
  guides(linetype=guide_legend(title=NULL),
         fill=guide_legend(title=NULL)) +
  theme_classic() +
  theme(legend.position=c(.8225,.82)) +
  ylab("Population affected\n(millions)") +
  xlab("Screened age group (years)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),  #12
        legend.title = element_text(size = 14))

# Cost by stage of care
pp4 <- ggplot(subset(discounted_interactions_cost_median, scenario %in% c("a2_screen_2020_monit_1",
                                                                          "a4_screen_2020_monit_1",
                                                                          "a5_screen_2020_monit_1") &
                       !(interaction_type %in% c("total_interactions", "treatment_initiations")))) +
  geom_col(aes(x=age_group, y= cost/1000000, fill = interaction_type),
           width = 0.95, colour = "black") +
  #  geom_errorbar(data=total_interactions_errorbar,
  #                aes(x = age_group, ymin = lower/1000, ymax = upper/1000),width = 0.25) +
  labs(fill = "Resource utilisation") +
  # scale_fill_viridis_d(labels = c("hbsag_tests" = "Serological tests",
  #                                 "clinical_assessments" = "Clinical assessments",
  #                                 "py_on_treatment" = "Treatment",
  #                                 "monitoring_assessments" = "Monitoring")) +
  scale_fill_manual(labels = c("hbsag_tests" = "Screening",
                               "clinical_assessments" = "Assessment",
                               "py_on_treatment" = "Treatment",
                               "monitoring_assessments" = "Monitoring"),
                    values= c("hbsag_tests" = "#FDE725FF",
                              "clinical_assessments" = "#35B779FF",
                              "py_on_treatment" = "#440154FF",
                              "monitoring_assessments" = "#31688EFF")) +
 # scale_y_continuous(breaks=c(0,2,4,6,8)) +
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
interactions_by_age_rel$with_monitoring[interactions_by_age_rel$scenario %in% c("a4_screen_2020_monit_1",
                                                                                "a5_screen_2020_monit_1",
                                                                                "a2_screen_2020_monit_1")] <- "Yes"
## New
# Want  treatment initiations break up into initial+monitoring

pp6 <- ggplot(subset(interactions_by_age_rel, with_monitoring == "No" & interaction_type==
                       "treatment_initiations_per_clinical_assessment")) +
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
        plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"))

# NEW:
pp6a <-ggplot(subset(interactions_by_age_rel, with_monitoring == "Yes" & interaction_type==
                "treatment_initiations_per_clinical_assessment")) +
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
  labs(title="Assessment") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14),
        strip.text = element_text(size = 14),
        plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"))


pp6b <- ggplot(subset(interactions_by_age_rel, with_monitoring == "Yes" & interaction_type==
                  "treatment_initiations_per_monitoring_assessment")) +
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
    scale_y_continuous(breaks=c(0,1,2)) +
    guides(colour=FALSE) +
    coord_cartesian(ylim=c(0,2))+
    xlab("Age group (years)")+
    ylab("Treatment eligibility (%)") +
    labs(title="Monitoring") +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 14),
          #axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size=14),
          axis.title=element_blank(),  #remove y axis labels
          strip.text = element_text(size = 14),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"))

pp6 <- grid.arrange(pp6a,pp6b,ncol=2,
             left=textGrob("Treatment eligibility (%)",rot=90,gp=gpar(fontsize=14)),
             bottom=textGrob("Age groups (years)",gp=gpar(fontsize=14)))

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
treatment_panel_plot <- grid.arrange(ppa, expl_plots,
                                     ncol=2, widths=2:1)
#dev.off()

# NEED TO ADD IN LEGEND THAT Population affected for monitoring
# actually means monitoring assessments

# ggsave("basic_programme_explanatory_plot_annual_monitoring_to_update.pdf", plot = treatment_panel_plot ,
#        width= 36.5, height=21.3, units = "cm")
# ggsave("basic_programme_explanatory_plot_annual_monitoring_to_update.png", plot = treatment_panel_plot ,
#        width= 36.5, height=21.3, units = "cm")







