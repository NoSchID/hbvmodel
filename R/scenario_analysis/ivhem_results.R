# Results for IVHEM poster
require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
library(BCEA)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# May need to add monitoring on treatment here!! #

# Load data ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

out_path2 <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/ivhem_simulations/"

# Status quo
out1 <- readRDS(paste0(out_path, "a1_out1_status_quo_cohort_240920.rds"))
out1 <- out1[[1]]
out2 <- readRDS(paste0(out_path, "out2_status_quo_180820.rds"))
out2 <- out2[[1]]

# Age-specific cohorts without treatment
a5_out1 <- readRDS(paste0(out_path2, "a5_out1_status_quo_cohort_181120.rds"))
a5_out1 <- a5_out1[[1]]
a4_out1 <- readRDS(paste0(out_path2, "a4_out1_status_quo_cohort_191120.rds"))
a4_out1 <- a4_out1[[1]]
a2_out1 <- readRDS(paste0(out_path2, "a2_out1_status_quo_cohort_191120.rds"))
a2_out1 <- a2_out1[[1]]
a6_out1 <- readRDS(paste0(out_path2, "a6_out1_status_quo_cohort_191120.rds"))
a6_out1 <- a6_out1[[1]]

# Status quo cohorts by age group

# Status quo full code_model_output
out2_carriers <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_carriers.rds")

# No monitoring by age
# A4 (15-30) A5 (30-45) A2 (45-60) A6 (60-65) A7 (65-70)
# A1 = 15-60
a4_out3 <- readRDS(paste0(out_path, "a4_out3_screen_2020_monit_0_231020.rds"))
a4_out3 <- a4_out3[[1]]
a5_out3 <- readRDS(paste0(out_path, "a5_out3_screen_2020_monit_0_231020.rds"))
a5_out3 <- a5_out3[[1]]
a2_out3 <- readRDS(paste0(out_path, "a2_out3_screen_2020_monit_0_161020.rds"))
a2_out3 <- a2_out3[[1]]
a6_out3 <- readRDS(paste0(out_path2, "a6_out3_screen_2020_monit_0_041120.rds"))
a6_out3 <- a6_out3[[1]]
a7_out3 <- readRDS(paste0(out_path2, "a7_out3_screen_2020_monit_0_041120.rds"))
a7_out3 <- a7_out3[[1]]
a1_out3 <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_201020.rds"))
a1_out3 <- a1_out3[[1]]

# 5-yearly monitoring by age
a2_out5 <- readRDS(paste0(out_path, "a2_out5_screen_2020_monit_5_161020.rds"))
a2_out5 <- a2_out5[[1]]
a5_out5 <- readRDS(paste0(out_path, "a5_out5_screen_2020_monit_5_231020.rds"))
a5_out5 <- a5_out5[[1]]
a4_out5 <- readRDS(paste0(out_path, "a4_out5_screen_2020_monit_5_231020.rds"))
a4_out5 <- a4_out5[[1]]
a6_out5 <- readRDS(paste0(out_path2, "a6_out5_screen_2020_monit_5_041120.rds"))
a6_out5 <- a6_out5[[1]]
a7_out5 <- readRDS(paste0(out_path2, "a7_out5_screen_2020_monit_5_041120.rds"))
a7_out5 <- a7_out5[[1]]
a1_out5 <- readRDS(paste0(out_path, "a1_out5_screen_2020_monit_5_201020.rds"))
a1_out5 <- a1_out5[[1]]

# Yearly monitoring by age
a2_out6 <- readRDS(paste0(out_path, "a2_out6_screen_2020_monit_1_161020.rds"))
a2_out6 <- a2_out6[[1]]
a5_out6 <- readRDS(paste0(out_path, "a5_out6_screen_2020_monit_1_271020.rds"))
a5_out6 <- a5_out6[[1]]
a4_out6 <- readRDS(paste0(out_path2, "a4_out6_screen_2020_monit_1_041120.rds"))
a4_out6 <- a4_out6[[1]]
a6_out6 <- readRDS(paste0(out_path2, "a6_out6_screen_2020_monit_1_041120.rds"))
a6_out6 <- a6_out6[[1]]
a7_out6 <- readRDS(paste0(out_path2, "a7_out6_screen_2020_monit_1_041120.rds"))
a7_out6 <- a7_out6[[1]]
a1_out6 <- readRDS(paste0(out_path, "a1_out6_screen_2020_monit_1_201020.rds"))
a1_out6 <- a1_out6[[1]]

# REALISTIC UPTAKE SIMS (15-65)
cx_out3 <- readRDS(paste0(out_path2, "cx_out3_screen_2020_monit_0_041120.rds"))
cx_out3 <- cx_out3[[1]]
cx_out5 <- readRDS(paste0(out_path2, "cx_out5_screen_2020_monit_5_041120.rds"))
cx_out5 <- cx_out5[[1]]
cx_out6 <- readRDS(paste0(out_path2, "cx_out6_screen_2020_monit_1_041120.rds"))
cx_out6 <- cx_out6[[1]]

# Function to plot boxplot whiskers as 95% percentile
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# Age group analysis ----

# Number of deaths averted in each age group
deaths_averted_by_age1 <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                                 scenario_objects = list(a4_out3,
                                                         a4_out5,
                                                         a4_out6),
                                 counterfactual_label = "treatment programme without monitoring")
deaths_averted_by_age2 <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a5_out3,
                                                  a5_out5,
                                                  a5_out6),
                          counterfactual_label = "treatment programme without monitoring")
deaths_averted_by_age3 <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a2_out3,
                                                  a2_out5,
                                                  a2_out6),
                          counterfactual_label = "treatment programme without monitoring")

deaths_averted_by_age4 <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a6_out3,
                                                  a6_out5,
                                                  a6_out6),
                          counterfactual_label = "treatment programme without monitoring")

deaths_averted_by_age5 <-
  plot_hbv_deaths_averted(counterfactual_object = out2,
                          scenario_objects = list(a7_out3,
                                                  a7_out5,
                                                  a7_out6),
                          counterfactual_label = "treatment programme without monitoring")

# Add 4 and 5 together:
#deaths_averted_by_age4 <- left_join((filter(deaths_averted_by_age4, by_year == 2100 &
#         type == "number_averted") %>%
#           select(scenario, sim, value)),
#         (filter(deaths_averted_by_age5, by_year == 2100 &
#                   type == "number_averted") %>%
#            select(scenario, sim, value)), by = c("scenario", "sim")) %>%
#  mutate(value = value.x+value.y) %>%
#  select(-value.x, -value.y)


deaths_averted_by_age_all_scenarios <- rbind(
  data.frame(age_group = "15-30",
                    filter(deaths_averted_by_age1, by_year == 2100 &
                       type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age2, by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-60",
             filter(deaths_averted_by_age3, by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "60-65",
             filter(deaths_averted_by_age4, by_year == 2100 &
                      type == "number_averted") %>% select(scenario, sim, value))
  )

deaths_averted_by_age_all_scenarios_prop <- rbind(
  data.frame(age_group = "15-30",
             filter(deaths_averted_by_age1, by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_averted_by_age2, by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "45-60",
             filter(deaths_averted_by_age3, by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value)),
  data.frame(age_group = "60-65",
             filter(deaths_averted_by_age4, by_year == 2100 &
                      type == "proportion_averted") %>% select(scenario, sim, value))
)

deaths_averted_by_age <- subset(deaths_averted_by_age_all_scenarios, scenario == "screen_2020_monit_0")
deaths_averted_by_age$sim <- gsub("[^0-9]", "", deaths_averted_by_age$sim)

# Number of life-years saved in each age group
# CALCULATE LY IN THE COHORT to avoid inclusion of extra births
ly_saved_by_age1 <-
  plot_ly_gained_cohort(counterfactual_object = a4_out1,
                          scenario_objects = list(a4_out3,
                                                  a4_out5,
                                                  a4_out6),
                          counterfactual_label = "treatment programme without monitoring")
ly_saved_by_age2 <-
  plot_ly_gained_cohort(counterfactual_object = a5_out1,
                          scenario_objects = list(a5_out3,
                                                  a5_out5,
                                                  a5_out6),
                          counterfactual_label = "treatment programme without monitoring")
ly_saved_by_age3 <-
  plot_ly_gained_cohort(counterfactual_object = a2_out1,
                          scenario_objects = list(a2_out3,
                                                  a2_out5,
                                                  a2_out6),
                          counterfactual_label = "treatment programme without monitoring")

ly_saved_by_age4 <-
  plot_ly_gained_cohort(counterfactual_object = a6_out1,
                          scenario_objects = list(a6_out3,
                                                  a6_out5,
                                                  a6_out6),
                          counterfactual_label = "treatment programme without monitoring")
#ly_saved_by_age5 <-
#  plot_ly_gained_cohort(counterfactual_object = out2,
#                          scenario_objects = list(a7_out3,
#                                                  a7_out5,
#                                                  a7_out6),
#                          counterfactual_label = "treatment programme without monitoring")

# Add 4 and 5 together:
#ly_saved_by_age4 <- left_join((filter(ly_saved_by_age4, by_year == 2100 &
#                                              type == "number_averted") %>%
#                                       select(counterfactual, sim, value)),
#                                    (filter(ly_saved_by_age5, by_year == 2100 &
#                                              type == "number_averted") %>%
#                                       select(counterfactual, sim, value)), by = c("counterfactual", "sim")) %>%
#  mutate(value = value.x+value.y) %>%
#  select(-value.x, -value.y)


ly_saved_by_age_all_scenarios <- rbind(
  data.frame(age_group = "15-30",
             filter(ly_saved_by_age1,
                      type == "number_averted") %>% select(counterfactual, sim, value)),
  data.frame(age_group = "30-45",
             filter(ly_saved_by_age2,
                      type == "number_averted") %>% select(counterfactual, sim, value)),
  data.frame(age_group = "45-60",
             filter(ly_saved_by_age3,
                      type == "number_averted") %>% select(counterfactual, sim, value)),
  data.frame(age_group = "60-65",
             filter(ly_saved_by_age4,
                      type == "number_averted") %>% select(counterfactual, sim, value))
)
colnames(ly_saved_by_age_all_scenarios)[colnames(ly_saved_by_age_all_scenarios)=="counterfactual"] <- "scenario"
ly_saved_by_age <- subset(ly_saved_by_age_all_scenarios, scenario == "screen_2020_monit_0")

ly_saved_by_age$sim <- gsub("[^0-9]", "", ly_saved_by_age$sim)
colnames(ly_saved_by_age)[colnames(ly_saved_by_age)=="value"] <- "ly_saved"

# Add carriers by age in 2020
# Calculate number of carriers by age group
carriers_female <- lapply(out2_carriers, "[[", "carriers_female")
carriers_male <- lapply(out2_carriers, "[[", "carriers_male")
total_carriers_by_age_2020 <- do.call(rbind.data.frame, (lapply(carriers_female, function(x) x[which(out2_carriers[[1]]$time==2020),])))+
  do.call(rbind.data.frame, (lapply(carriers_male, function(x) x[which(out2_carriers[[1]]$time==2020),])))

carriers_by_age_group_2020 <-
  rbind(data.frame(age_group = "15-30",
                   sim = rownames(total_carriers_by_age_2020[,which(ages==15): which(ages==30-da)]),
                   carriers = rowSums(total_carriers_by_age_2020[,which(ages==15): which(ages==30-da)])),
        data.frame(age_group = "30-45",
                   sim = rownames(total_carriers_by_age_2020[,which(ages==30): which(ages==45-da)]),
                   carriers = rowSums(total_carriers_by_age_2020[,which(ages==30): which(ages==45-da)])),
        data.frame(age_group = "45-60",
                   sim = rownames(total_carriers_by_age_2020[,which(ages==45): which(ages==60)]),
                   carriers = rowSums(total_carriers_by_age_2020[,which(ages==45): which(ages==60)])),
        data.frame(age_group = "60-65",
                   sim = rownames(total_carriers_by_age_2020[,which(ages==60): which(ages==65-da)]),
              carriers = rowSums(total_carriers_by_age_2020[,which(ages==60): which(ages==65-da)]))
  )

# Total interactions by age and type (counting treatment initiation only)
interactions_by_age <- rbind(
  data.frame(age_group = "15-30",
           sim = names(a4_out3$interactions[[16]]$total_screened[,-c(1:3)]),
           hbsag_tests = unlist(a4_out3$interactions[[16]]$total_screened[,-c(1:3)]),
           clinical_assessments = unlist(a4_out3$interactions[[16]]$total_assessed[,-c(1:3)]),
           treatment_initiations = unlist(a4_out3$interactions[[16]]$total_treated[,-c(1:3)])),
  data.frame(age_group = "30-45",
             sim = names(a5_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a5_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a5_out3$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a5_out3$interactions[[16]]$total_treated[,-c(1:3)])),
  data.frame(age_group = "45-60",
             sim = names(a2_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a2_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a2_out3$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a2_out3$interactions[[16]]$total_treated[,-c(1:3)])),
  data.frame(age_group = "60-65",
             sim = names(a6_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             hbsag_tests = unlist(a6_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             clinical_assessments = unlist(a6_out3$interactions[[16]]$total_assessed[,-c(1:3)]),
             treatment_initiations = unlist(a6_out3$interactions[[16]]$total_treated[,-c(1:3)]))
)
interactions_by_age <- interactions_by_age %>%
  mutate(total_interactions = hbsag_tests + clinical_assessments + treatment_initiations)
interactions_by_age <- gather(interactions_by_age, key = "interaction_type", value = "value",
                              -age_group, - sim)

# Merge deaths averted by age with carrier population by age and with LY saved
outcomes_by_age <- left_join(left_join(deaths_averted_by_age, carriers_by_age_group_2020,
                                   by = c("age_group", "sim")),
                             ly_saved_by_age, by = c("age_group","scenario", "sim")) %>%
  mutate(deaths_averted_per_carrier = value/carriers,
         ly_saved_per_carrier = ly_saved/carriers)
colnames(outcomes_by_age)[colnames(outcomes_by_age)=="value"] <- "deaths_averted"
outcomes_by_age <- gather(outcomes_by_age, key = "outcome", value = "value",
                                -age_group, -scenario, -sim, -deaths_averted_per_carrier, -ly_saved_per_carrier)


# Deaths averted (absolute)
ggplot(subset(outcomes_by_age, outcome=="deaths_averted")) +
  geom_boxplot(aes(x=age_group, y=value)) +
  ylim(0,4500)
ggplot(subset(outcomes_by_age, outcome=="carriers")) +
  geom_boxplot(aes(x=age_group, y=value))
# Deaths averted per 10000 carriers
ggplot(outcomes_by_age) +
  geom_boxplot(aes(x=age_group, y=deaths_averted_per_carrier*10000))+
  ylim(0,1100)
ggplot(outcomes_by_age) +
  geom_boxplot(aes(x=age_group, y=ly_saved_per_carrier*10000))
# Differences in deaths averted by age group are mostly driven by differences in the
# carrier population, which is highest in the 30-45 group.
# Deaths averted per carrier is very similar across age groups.
ggplot(outcomes_by_age) +
  geom_bar(aes(x=age_group, y= value),
         position = "dodge", stat = "summary", fun = "median") +
  geom_errorbar(mapping = aes(x=age_group, y= value, group = outcome),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                width = 0.25) +
  facet_wrap(~outcome, scales="free_y") +
  theme_classic()
# Note that outcomes here are for the AGE GROUP TARGETED IN SCREEN (NOT the outcome in that age group)
# Maybe make a stacked bar chart of outcome by age group? And then same for interactions?

p1 <- ggplot(subset(outcomes_by_age, outcome == "deaths_averted")) +
  geom_bar(aes(x=age_group, y= value/1000, group = outcome), fill = "grey40",
           stat = "summary", fun = "median", width = 0.9) +
  geom_errorbar(mapping = aes(x=age_group, y= value/1000, group = outcome),
                stat = "summary",
                fun.min = function(z) {quantile(z,0.025)},
                fun.max = function(z) {quantile(z,0.975)},
                width = 0.25) +
  ylab("HBV deaths\naverted\n(thousands)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))

p2<- ggplot(subset(outcomes_by_age, outcome == "ly_saved")) +
  geom_bar(aes(x=age_group, y= value/1000), fill = "grey40",
           position = "dodge", stat = "summary", fun = "median", width = 0.9) +
  geom_errorbar(mapping = aes(x=age_group, y= value/1000, group = outcome),
                stat = "summary",
                fun.min = function(z) {quantile(z,0.025)},
                fun.max = function(z) {quantile(z,0.975)},
                width = 0.25) +
  ylab("Life-years\nsaved\n(thousands)") +
  theme_classic()+
  theme(axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))

total_interactions_errorbar <- subset(interactions_by_age, interaction_type == "total_interactions") %>%
  group_by(age_group) %>%
  summarise(lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

library(viridis)
options(scipen=1000000)
median_carriers <- subset(outcomes_by_age, outcome=="carriers") %>% group_by(age_group) %>%
  summarise(median = median(value))

p3 <- ggplot(subset(interactions_by_age, interaction_type != "total_interactions")) +
  geom_bar(aes(x=age_group, y= value/1000, fill = reorder(interaction_type, -value)),
           stat = "summary", fun = "median", width = 0.9) +
  geom_errorbar(data=median_carriers,
           aes(x=age_group, y= median/1000, ymin= median/1000, ymax= median/1000, linetype = "HBV carriers"),
           width = 0.9, size = 1) +
  geom_errorbar(data=total_interactions_errorbar,
                aes(x = age_group, ymin = lower/1000, ymax = upper/1000),width = 0.25) +
  labs(fill = "Resource utilisation") +
  scale_fill_viridis_d(direction=-1,
                       labels = c("hbsag_tests" = "HBsAg tests",
                                  "clinical_assessments" = "Clinical\nassessments",
                                  "treatment_initiations" = "Treatment\ninitiations")) +
  scale_linetype_manual(values = c("HBV carriers" = "dashed")) +
  guides(linetype=guide_legend(title=NULL),
         fill=guide_legend(title=NULL)) +
  theme_classic() +
  theme(legend.position=c(.78,.8)) +
  ylab("Resources utilised\n(thousands)") +
  xlab("Screened age group")+
  theme(axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13))

#grid.arrange(grid.arrange(p1,p2,ncol=2), p3, nrow=2)

grid.arrange(p1,p2,p3,ncol=1, heights = c(1,1,3))
# Rethink carrier line

deaths_per_int1 <- plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                     scenario_objects = list(a4_out3),
                                                     interaction_type = "total_interactions",
                                                     counterfactual_label = "no treatment")
deaths_per_int2 <- plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                                      scenario_objects = list(a5_out3),
                                                                      interaction_type = "total_interactions",
                                                                      counterfactual_label = "no treatment")
deaths_per_int3 <- plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                                      scenario_objects = list(a2_out3),
                                                                      interaction_type = "total_interactions",
                                                                      counterfactual_label = "no treatment")


deaths_per_int <- rbind(
  data.frame(age_group = "15-30",
             filter(deaths_per_int1, by_year == 2100) %>% select(scenario, sim, value)),
  data.frame(age_group = "30-45",
             filter(deaths_per_int2, by_year == 2100) %>% select(scenario, sim, value)),
  data.frame(age_group = "45-60",
             filter(deaths_per_int3, by_year == 2100) %>% select(scenario, sim, value))
)

ggplot(deaths_per_int) +
  geom_boxplot(aes(x=age_group, y=value*10000)) +
  ylab("HBV-related deaths averted per 10,000 interactions")

# Deaths averted per interaction by far the lowest in the 15-30 age group,
# reflecting the large number of people of that age group yet low prevalence (many people to test
# but proportionally fewer identified)
# Age groups in abstract ----
# 15-65 = A1 + A6
# 30-65 = A5+A2+A6
# 15-45 = A4+A5
# 45-65 = A2+A6

# Proportion of deaths averted
deaths_averted_by_age_all_scenarios_prop$sim <- paste0("X", gsub("[^0-9]", "", deaths_averted_by_age_all_scenarios_prop$sim))
deaths_averted_by_age_prop_wide <- spread(deaths_averted_by_age_all_scenarios_prop, key = "age_group",
                                     value = "value")
deaths_averted_by_age_prop_wide$age15to65 <- deaths_averted_by_age_prop_wide$`15-30`+
  deaths_averted_by_age_prop_wide$`30-45`+deaths_averted_by_age_prop_wide$`45-60`+
  deaths_averted_by_age_prop_wide$`60-65`
deaths_averted_by_age_prop_wide$age30to65 <- deaths_averted_by_age_prop_wide$`30-45`+
  deaths_averted_by_age_prop_wide$`45-60`+
  deaths_averted_by_age_prop_wide$`60-65`
deaths_averted_by_age_prop_wide$age15to45 <- deaths_averted_by_age_prop_wide$`15-30`+
  deaths_averted_by_age_prop_wide$`30-45`
deaths_averted_by_age_prop_wide$age45to65 <- deaths_averted_by_age_prop_wide$`45-60`+
  deaths_averted_by_age_prop_wide$`60-65`
deaths_averted_by_age_prop_wide <- deaths_averted_by_age_prop_wide[, -c(3:6)]
deaths_averted_by_age_prop_long <- gather(deaths_averted_by_age_prop_wide, key = "age_group", value = "deaths_averted",
                                     -scenario, -sim)

deaths_averted_by_age_prop_long  %>%
  group_by(scenario, age_group) %>%
  summarise(median = median(deaths_averted),
            lower = quantile(deaths_averted, prob = 0.025),
            upper = quantile(deaths_averted, prob = 0.975))

# Other outcomes

deaths_averted_by_age_all_scenarios$sim <- paste0("X", gsub("[^0-9]", "", deaths_averted_by_age_all_scenarios$sim))
deaths_averted_by_age_wide <- spread(deaths_averted_by_age_all_scenarios, key = "age_group",
                                     value = "value")
deaths_averted_by_age_wide$age15to65 <- deaths_averted_by_age_wide$`15-30`+
  deaths_averted_by_age_wide$`30-45`+deaths_averted_by_age_wide$`45-60`+
  deaths_averted_by_age_wide$`60-65`
deaths_averted_by_age_wide$age30to65 <- deaths_averted_by_age_wide$`30-45`+
  deaths_averted_by_age_wide$`45-60`+
  deaths_averted_by_age_wide$`60-65`
deaths_averted_by_age_wide$age15to45 <- deaths_averted_by_age_wide$`15-30`+
  deaths_averted_by_age_wide$`30-45`
deaths_averted_by_age_wide$age45to65 <- deaths_averted_by_age_wide$`45-60`+
  deaths_averted_by_age_wide$`60-65`
deaths_averted_by_age_wide <- deaths_averted_by_age_wide[, -c(3:6)]
deaths_averted_by_age_long <- gather(deaths_averted_by_age_wide, key = "age_group", value = "deaths_averted",
                                     -scenario, -sim)

# LY saved
ly_saved_by_age_all_scenarios$sim <- paste0("X", gsub("[^0-9]", "", ly_saved_by_age_all_scenarios$sim))
ly_saved_by_age_wide <- spread(ly_saved_by_age_all_scenarios, key = "age_group",
                                     value = "value")
ly_saved_by_age_wide$age15to65 <- ly_saved_by_age_wide$`15-30`+
  ly_saved_by_age_wide$`30-45`+ly_saved_by_age_wide$`45-60`+
  ly_saved_by_age_wide$`60-65`
ly_saved_by_age_wide$age30to65 <- ly_saved_by_age_wide$`30-45`+
  ly_saved_by_age_wide$`45-60`+
  ly_saved_by_age_wide$`60-65`
ly_saved_by_age_wide$age15to45 <- ly_saved_by_age_wide$`15-30`+
  ly_saved_by_age_wide$`30-45`
ly_saved_by_age_wide$age45to65 <- ly_saved_by_age_wide$`45-60`+
  ly_saved_by_age_wide$`60-65`
ly_saved_by_age_wide <- ly_saved_by_age_wide[, -c(3:6)]
ly_saved_by_age_long <- gather(ly_saved_by_age_wide, key = "age_group", value = "ly_saved",
                                     -scenario, -sim)

interactions_by_age_long <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age15to65",
             sim = names(a1_out3$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a1_out3$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out3$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age15to65",
             sim = names(a1_out5$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a1_out5$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age15to65",
             sim = names(a1_out6$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a1_out6$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age30to65",
             sim = names(a5_out3$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a5_out3$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a2_out3$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out3$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age30to65",
             sim = names(a5_out5$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a5_out5$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a2_out5$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age30to65",
             sim = names(a5_out6$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a5_out6$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a2_out6$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age15to45",
             sim = names(a4_out3$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a4_out3$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a5_out3$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age15to45",
             sim = names(a4_out5$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a4_out5$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a5_out5$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age15to45",
             sim = names(a4_out6$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a4_out6$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a5_out6$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age45to65",
             sim = names(a2_out3$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a2_out3$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out3$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age45to65",
             sim = names(a2_out5$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a2_out5$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_interactions[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age45to65",
             sim = names(a2_out6$interactions[[16]]$total_interactions[,-c(1:3)]),
             interactions = unlist(a2_out6$interactions[[16]]$total_interactions[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_interactions[,-c(1:3)]))
)

screening_by_age_long <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age15to65",
             sim = names(a1_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a1_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out3$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age15to65",
             sim = names(a1_out5$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a1_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age15to65",
             sim = names(a1_out6$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a1_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age30to65",
             sim = names(a5_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a5_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a2_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out3$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age30to65",
             sim = names(a5_out5$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a5_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a2_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age30to65",
             sim = names(a5_out6$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a5_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a2_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age15to45",
             sim = names(a4_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a4_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a5_out3$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age15to45",
             sim = names(a4_out5$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a4_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a5_out5$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age15to45",
             sim = names(a4_out6$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a4_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a5_out6$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age45to65",
             sim = names(a2_out3$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a2_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out3$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age45to65",
             sim = names(a2_out5$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a2_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age45to65",
             sim = names(a2_out6$interactions[[16]]$total_screened[,-c(1:3)]),
             screening = unlist(a2_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_screened[,-c(1:3)]))
)

assessment_and_treatment_by_age_long <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age15to65",
             sim = names(a1_out3$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a1_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a1_out3$interactions[[16]]$total_screened[,-c(1:3)]) +
               unlist(a6_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out3$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age15to65",
             sim = names(a1_out5$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a1_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a1_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out5$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age15to65",
             sim = names(a1_out6$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a1_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a1_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out6$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age30to65",
             sim = names(a5_out3$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a5_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a5_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a2_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a2_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out3$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age30to65",
             sim = names(a5_out5$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a5_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a5_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a2_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a2_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out5$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age30to65",
             sim = names(a5_out6$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a5_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a5_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a2_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a2_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out6$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age15to45",
             sim = names(a4_out3$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a4_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a4_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a5_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a5_out3$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age15to45",
             sim = names(a4_out5$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a4_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a4_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a5_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a5_out5$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age15to45",
             sim = names(a4_out6$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a4_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a4_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a5_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a5_out6$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_0",
             age_group = "age45to65",
             sim = names(a2_out3$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a2_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a2_out3$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out3$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out3$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_5",
             age_group = "age45to65",
             sim = names(a2_out5$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a2_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a2_out5$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out5$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out5$interactions[[16]]$total_screened[,-c(1:3)])),
  data.frame(scenario = "screen_2020_monit_1",
             age_group = "age45to65",
             sim = names(a2_out6$interactions[[16]]$total_interactions[,-c(1:3)]),
             assessment_and_treatment = unlist(a2_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a2_out6$interactions[[16]]$total_screened[,-c(1:3)])+
               unlist(a6_out6$interactions[[16]]$total_interactions[,-c(1:3)])-
               unlist(a6_out6$interactions[[16]]$total_screened[,-c(1:3)]))
)

outcome_per_int_by_original_groups <- left_join(left_join(left_join(
    left_join(deaths_averted_by_age_long, ly_saved_by_age_long,
              by =  c("scenario", "age_group", "sim")),
    interactions_by_age_long, by = c("scenario", "age_group", "sim")),
    assessment_and_treatment_by_age_long,  by = c("scenario", "age_group", "sim")),
    screening_by_age_long, by = c("scenario", "age_group", "sim")) %>%
    mutate(deaths_averted_per_interaction = deaths_averted/interactions,
           ly_saved_per_interaction = ly_saved/interactions,
           deaths_averted_per_assessment_and_treatment = deaths_averted/assessment_and_treatment,
           ly_saved_per_assessment_and_treatment = ly_saved/assessment_and_treatment,
           deaths_averted_per_screening = deaths_averted/screening,
           ly_saved_per_screening = ly_saved/screening)

# Deaths averted per 10,000 total interactions
outcome_per_int_by_original_groups  %>%
  group_by(scenario, age_group) %>%
  summarise(median_per_10000 = median(deaths_averted_per_interaction)*10000,
            lower_per_10000 = quantile(deaths_averted_per_interaction, prob = 0.025)*10000,
            upper_per_10000 = quantile(deaths_averted_per_interaction, prob = 0.975)*10000)

# LY saved per 10,000 total interactions
outcome_per_int_by_original_groups  %>%
  group_by(scenario, age_group) %>%
  summarise(median_per_10000 = median(ly_saved_per_interaction)*10000,
            lower_per_10000 = quantile(ly_saved_per_interaction, prob = 0.025)*10000,
            upper_per_10000 = quantile(ly_saved_per_interaction, prob = 0.975)*10000)

# Deaths averted per 1000 screens
outcome_per_int_by_original_groups  %>%
  group_by(scenario, age_group) %>%
  summarise(median_per_10000 = median(deaths_averted_per_screening)*1000,
            lower_per_10000 = quantile(deaths_averted_per_screening, prob = 0.025)*1000,
            upper_per_10000 = quantile(deaths_averted_per_screening, prob = 0.975)*1000)

# LY saved per 1000 screens
outcome_per_int_by_original_groups  %>%
  group_by(scenario, age_group) %>%
  summarise(median_per_10000 = median(ly_saved_per_screening)*1000,
            lower_per_10000 = quantile(ly_saved_per_screening, prob = 0.025)*1000,
            upper_per_10000 = quantile(ly_saved_per_screening, prob = 0.975)*1000)

# Deaths averted per 1000 clinical assessments and treatments
outcome_per_int_by_original_groups  %>%
  group_by(scenario, age_group) %>%
  summarise(median_per_10000 = median(deaths_averted_per_assessment_and_treatment)*1000,
            lower_per_10000 = quantile(deaths_averted_per_assessment_and_treatment, prob = 0.025)*1000,
            upper_per_10000 = quantile(deaths_averted_per_assessment_and_treatment, prob = 0.975)*1000)

# LY saved per 1000 clinical assessments and treatments
outcome_per_int_by_original_groups  %>%
  group_by(scenario, age_group) %>%
  summarise(median_per_10000 = median(ly_saved_per_assessment_and_treatment)*1000,
            lower_per_10000 = quantile(ly_saved_per_assessment_and_treatment, prob = 0.025)*1000,
            upper_per_10000 = quantile(ly_saved_per_assessment_and_treatment, prob = 0.975)*1000)

# Absolute deaths averted
outcome_per_int_by_original_groups  %>%
  group_by(scenario, age_group) %>%
  summarise(median = median(deaths_averted),
            lower = quantile(deaths_averted, prob = 0.025),
            upper = quantile(deaths_averted, prob = 0.975))

# Absolute LY saved
outcome_per_int_by_original_groups  %>%
  group_by(scenario, age_group) %>%
  summarise(median = median(ly_saved),
            lower = quantile(ly_saved, prob = 0.025),
            upper = quantile(ly_saved, prob = 0.975))

median(subset(outcome_per_int_by_original_groups, scenario == "screen_2020_monit_0" &
         age_group == "age15to65")$deaths_averted/
  subset(outcome_per_int_by_original_groups, scenario == "screen_2020_monit_0" &
                                                  age_group == "age30to65")$deaths_averted)

quantile(subset(outcome_per_int_by_original_groups, scenario == "screen_2020_monit_0" &
                age_group == "age15to65")$ly_saved/
         subset(outcome_per_int_by_original_groups, scenario == "screen_2020_monit_0" &
                  age_group == "age30to65")$ly_saved, prob = c(0.5,0.025,0.975))

# Population impact of no monitoring vs no uptake ----
# Show 15-65, low uptake, no monitoring (cx_out3)
# 15-65, low uptake, with 5 monitoring (cx_out5)
# 15-65, low uptake, with yearly monitoring  (cx_out6)
# 15-65, high uptake, no monitoring/5 yearly monitoring (A1+A6) - monitoring to add!

# Absolute deaths averted
cx_deaths_averted <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                             scenario_objects = list(cx_out3,
                                                                     cx_out5,
                                                                     cx_out6),
                                             counterfactual_label = "no treatment")
a1_deaths_averted <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                             scenario_objects = list(a1_out3,
                                                                     a1_out5,
                                                                     a1_out6),
                                             counterfactual_label = "no treatment")
a6_deaths_averted <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                             scenario_objects = list(a6_out3,
                                                                     a6_out5,
                                                                     a6_out6),
                                             counterfactual_label = "no treatment")

deaths_averted_15to65 <- left_join(a1_deaths_averted, a6_deaths_averted,
                                   by = c("from_year", "by_year", "counterfactual", "scenario",
                                          "type", "sim")) %>%
  mutate(value = value.x+value.y) %>%
  filter(type == "number_averted", by_year == 2100) %>%
  select(scenario, sim, value)

deaths_averted_by_uptake <- rbind(
  cbind(deaths_averted_15to65, uptake = "ambitious"),
  cbind(filter(cx_deaths_averted, type == "number_averted", by_year == 2100) %>%
          select(scenario, sim, value), uptake = "realistic")
)

# Deaths averted per total interactions
cx_deaths_averted_per_int <- plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                   scenario_objects = list(cx_out3,
                                                                           cx_out5,
                                                                           cx_out6),
                                                   interaction_type = "total_interactions",
                                                   counterfactual_label = "no treatment") %>%
  filter(by_year == 2100) %>%
  select(scenario, sim, value)

# For combination of A1 and A6, need to calculate this manually
# For no monitoring:
deaths_averted_per_interaction_15to65a <-
  data.frame(scenario = "Ambitious no monitoring",
             sim = colnames(out2$cum_hbv_deaths[[16]][,-c(1:3)]),
             deaths_averted = unlist((out2$cum_hbv_deaths[[16]][,-c(1:3)]-a1_out3$cum_hbv_deaths[[16]][,-c(1:3)])+
             (out2$cum_hbv_deaths[[16]][,-c(1:3)]-a6_out3$cum_hbv_deaths[[16]][,-c(1:3)])),
             total_interactions = unlist(a1_out3$interactions[[16]]$total_interactions[,-c(1:3)]+
               a6_out3$interactions[[16]]$total_interactions[,-c(1:3)]))   %>%
  mutate(value = deaths_averted/total_interactions) %>%
  select(-deaths_averted, - total_interactions)

deaths_averted_per_interaction_15to65b <-
  data.frame(scenario = "Ambitious 5-yearly",
             sim = colnames(out2$cum_hbv_deaths[[16]][,-c(1:3)]),
             deaths_averted = unlist((out2$cum_hbv_deaths[[16]][,-c(1:3)]-a1_out5$cum_hbv_deaths[[16]][,-c(1:3)])+
                                       (out2$cum_hbv_deaths[[16]][,-c(1:3)]-a6_out5$cum_hbv_deaths[[16]][,-c(1:3)])),
             total_interactions = unlist(a1_out5$interactions[[16]]$total_interactions[,-c(1:3)]+
                                           a6_out5$interactions[[16]]$total_interactions[,-c(1:3)]))   %>%
  mutate(value = deaths_averted/total_interactions) %>%
  select(-deaths_averted, - total_interactions)

deaths_averted_per_interaction_15to65c <-
  data.frame(scenario = "Ambitious yearly",
             sim = colnames(out2$cum_hbv_deaths[[16]][,-c(1:3)]),
             deaths_averted = unlist((out2$cum_hbv_deaths[[16]][,-c(1:3)]-a1_out6$cum_hbv_deaths[[16]][,-c(1:3)])+
                                       (out2$cum_hbv_deaths[[16]][,-c(1:3)]-a6_out6$cum_hbv_deaths[[16]][,-c(1:3)])),
             total_interactions = unlist(a1_out6$interactions[[16]]$total_interactions[,-c(1:3)]+
                                           a6_out6$interactions[[16]]$total_interactions[,-c(1:3)]))   %>%
  mutate(value = deaths_averted/total_interactions) %>%
  select(-deaths_averted, - total_interactions)


deaths_averted_per_interaction_by_uptake <- rbind(cx_deaths_averted_per_int,
                                                  deaths_averted_per_interaction_15to65a,
                                                  deaths_averted_per_interaction_15to65b,
                                                  deaths_averted_per_interaction_15to65c)
deaths_averted_per_interaction_by_uptake$uptake <- "Low"
deaths_averted_per_interaction_by_uptake$uptake[
  deaths_averted_per_interaction_by_uptake$scenario %in% c("Ambitious no monitoring",
                                                           "Ambitious 5-yearly",
                                                           "Ambitious yearly")] <- "Ambitious"
deaths_averted_per_interaction_by_uptake$scenario[deaths_averted_per_interaction_by_uptake$scenario ==
                                                    "Ambitious no monitoring"] <- "screen_2020_monit_0"
deaths_averted_per_interaction_by_uptake$scenario[deaths_averted_per_interaction_by_uptake$scenario ==
                                                    "Ambitious 5-yearly"] <- "screen_2020_monit_5"
deaths_averted_per_interaction_by_uptake$scenario[deaths_averted_per_interaction_by_uptake$scenario ==
                                                    "Ambitious yearly"] <- "screen_2020_monit_1"
deaths_averted_per_interaction_by_uptake$scenario <- factor(deaths_averted_per_interaction_by_uptake$scenario,
                                                            levels = c("screen_2020_monit_0",
                                                                       "screen_2020_monit_5",
                                                                       "screen_2020_monit_1"))


px1 <- ggplot(deaths_averted_by_uptake,
       aes(x = scenario, y = value, colour = reorder(uptake, value), fill = reorder(uptake, value))) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge2", width = 0.5, alpha = 0.4, lwd = 1) +
  scale_x_discrete(labels = c("screen_2020_monit_0" = "No\nmonitoring",
                              "screen_2020_monit_5" = "Every\n5 years",
                              "screen_2020_monit_1" = "Every\n1 year")) +
  scale_fill_manual(labels = c("realistic" = "Feasible",
                                 "ambitious" = "Ambitious"),
                      values = c("realistic" = "grey40",
                                 "ambitious" = "#440154")) +
  scale_colour_manual(labels = c("realistic" = "Feasible",
                               "ambitious" = "Ambitious"),
                    values = c("realistic" = "grey40",
                               "ambitious" = "#440154")) +
  labs(fill = "Programme coverage") +
  ylab("HBV-related deaths averted") +
  xlab("Monitoring frequency") +
  guides(colour = "none") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  ylim(0,12000)

px2 <- ggplot(deaths_averted_per_interaction_by_uptake, aes(x = scenario, y = value*10000,
                                                     colour = reorder(uptake, value),
                                                     fill = reorder(uptake, value))) +
  stat_summary(fun.data=f, geom="boxplot", position = "dodge2", width = 0.5, alpha = 0.4, lwd = 1) +
  scale_x_discrete(labels = c("screen_2020_monit_0" = "No\nmonitoring",
                              "screen_2020_monit_5" = "Every\n5 years",
                              "screen_2020_monit_1" = "Every\n1 year")) +
  scale_fill_manual(labels = c("Low" = "Feasible"),
                               values = c("Low" = "grey40",
                               "Ambitious" = "#440154")) +
  scale_colour_manual(labels = c("Low" = "Feasible"),
                      values = c("Low" = "grey40",
                                 "Ambitious" = "#440154")) +
  labs(fill = "Programme coverage") +
  ylab("HBV-related deaths averted\nper 10,000 resources utilised") +
  xlab("Monitoring frequency") +
  guides(colour = "none") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  ylim(0,80)

library(ggpubr)
ggarrange(px1, px2, ncol=2, common.legend = TRUE, legend="bottom")

grid.arrange(px1, px2, ncol = 2)



