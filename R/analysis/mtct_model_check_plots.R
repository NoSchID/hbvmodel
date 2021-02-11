# Plots for check of proportion of new chronic infections due to MTCT ----

library(here)
library(ggplot2)

# Read in chronic_infection_incidence_status_quo_250920.rds
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/"
chronic_infection_incidence <- readRDS(paste0(out_path, "sq_chronic_infection_incidence_status_quo_250920.rds"))

prop_mtct_long <- chronic_infection_incidence$prop_mtct
# Extract the other 2...

prop_mtct_summary <- group_by(prop_mtct_long, time) %>%
  summarise(median = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

ggplot(prop_mtct_summary) +
  geom_line(aes(x = time, y = median*100), col = "red") +
  geom_ribbon(aes(x=time, ymin = lower*100, ymax = upper*100), alpha = 0.2) +
  geom_vline(xintercept = 1990) +
  scale_x_continuous("Time", breaks = seq(1990,2100,  by = 10), limits = c(1986,2100)) +
  scale_y_continuous("New chronic infections due to MTCT (%)", breaks = seq(0,100,  by = 10)) +
  theme_classic()

inf_horizontal_summary <- group_by(inf_horizontal_long, time) %>%
  summarise(median = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

ggplot(inf_horizontal_summary) +
  geom_line(aes(x = time, y = median*2), col = "red") +
  geom_ribbon(aes(x=time, ymin = lower*2, ymax = upper*2), alpha = 0.2) +
  #  xlim(1989,2100) +
  xlim(1985,2100) +
  ylim(0,8000) +
  theme_classic()

inf_mtct_summary <- group_by(inf_mtct_long, time) %>%
  summarise(median = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

ggplot(inf_mtct_summary) +
  geom_line(aes(x = time, y = median*2), col = "red") +
  geom_ribbon(aes(x=time, ymin = lower*2, ymax = upper*2), alpha = 0.2) +
  geom_vline(xintercept = 2005) +
  geom_vline(xintercept = 2035) +
  xlim(1985,2100) +
  ylim(0,3000) +
  theme_classic()

inf_all <- inf_horizontal_long + inf_mtct_long
inf_all$time <- rep(inf[[1]]$time, each = 183)

inf_all_summary <- group_by(inf_all, time) %>%
  summarise(median = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

chronic_inf_combi <- rbind(cbind(inf_horizontal_summary, route = "horizontal"),
                           cbind(inf_mtct_summary, route = "mtct"))

ggplot(chronic_inf_combi) +
  geom_line(aes(x = time, y = median*2, group = route, colour = route)) +
  geom_ribbon(aes(x=time, ymin = lower*2, ymax = upper*2, group = route,fill = route), alpha = 0.2) +
  geom_vline(xintercept = 1990) +
  labs(colour = "Route", fill = "Route") +
  scale_x_continuous("Time", breaks = seq(1990,2100,  by = 10), limits = c(1986,2100)) +
  scale_y_continuous("Annual new chronic HBV infections", breaks = seq(0,8000,  by = 1000),
                     limits = c(0,3000)) +
  theme_classic()

# Model check for effective vaccine coverage ----
load(file = here("calibration", "output", "model_fit_output_kmeans_221220.Rdata"))
load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans

# Check effective vaccine coverage in 1 year olds (i.e. the % of 1 year olds in immune compartment)
# - the proportion that would have been there historically due to recovery

new_vacc_meth <- out_mat
#output1 <- new_vacc_meth2
#output2 <- new_vacc_meth
#new_vacc_meth <- output2

# For vaccine coverage, have to substract earlier recovery % (ca 4 %)
# Numerator = number in immune compartment = ever_infected - carriers
# Denominator = total population aged 1 year
vacc_cov_2020 <- 183L
vacc_cov_1997 <- 183L
vacc_cov_2010 <- 183L
rec_1990 <- 183L
immune_prop_over_time <- matrix(0, ncol = 183, nrow = 280)

for (i in 1:183) {
  # In 1989.5
  rec_1990[i] <- (new_vacc_meth[[i]]$full_output$ever_infected[220,which(ages==1)]-
                    new_vacc_meth[[i]]$full_output$carriers[220,which(ages==1)])/
    new_vacc_meth[[i]]$full_output$pop[220,which(ages==1)]

  immune_prop_over_time[,i] <- (new_vacc_meth[[i]]$full_output$ever_infected[,which(ages==1)]-
                                  new_vacc_meth[[i]]$full_output$carriers[,which(ages==1)])/
    new_vacc_meth[[i]]$full_output$pop[,which(ages==1)]

  vacc_cov_1997[i] <- ((new_vacc_meth[[i]]$full_output$ever_infected[235,which(ages==1)]-
                          new_vacc_meth[[i]]$full_output$carriers[235,which(ages==1)])/
                         new_vacc_meth[[i]]$full_output$pop[235,which(ages==1)])-rec_1990[i]

  vacc_cov_2010[i] <- ((new_vacc_meth[[i]]$full_output$ever_infected[261,which(ages==1)]-
                          new_vacc_meth[[i]]$full_output$carriers[261,which(ages==1)])/
                         new_vacc_meth[[i]]$full_output$pop[261,which(ages==1)])-rec_1990[i]

  vacc_cov_2020[i] <- ((new_vacc_meth[[i]]$full_output$ever_infected[280,which(ages==1)]-
                          new_vacc_meth[[i]]$full_output$carriers[280,which(ages==1)])/
                         new_vacc_meth[[i]]$full_output$pop[280,which(ages==1)])-rec_1990[i]
}

test <- as.data.frame(cbind(vacc_cov_2020,unlist(params_mat_accepted_kmeans$vacc_eff)))
test$high_eff <- 1
test$high_eff[test$V2<0.97] <- 0

plot(x=1:183, y = vacc_cov_1997, ylim = c(0,1))
points(x=1:183, y = vacc_cov_2010, col = "blue")
points(x=1:183, y = vacc_cov_2020, col = "red")
points(x=1:183, y = params_mat_accepted_kmeans$vacc_eff, col = "green")

boxplot(vacc_cov_2020~high_eff, test, ylim = c(0.7,1))

# Effective coverage over time (subtract recovered proportion)
# Drop timesteps before 1990
immune_prop_over_time <- immune_prop_over_time[which(out_mat[[1]]$full_output$time>=1990),]

# Subtract rec_1990 from each row
eff_vacc_cov <- sweep(immune_prop_over_time, 2, rec_1990, `-`)
# Could adjust by vaccine efficacy:
#eff_vacc_cov  <- sweep(immune_prop_over_time, 2, params_mat_accepted_kmeans$vacc_eff, `/`)
eff_vacc_cov <- data.frame(year = seq(1990, max(out_mat[[1]]$full_output$time), by =0.5),
                          eff_vacc_cov)
eff_vacc_cov <- gather(eff_vacc_cov, key = "sim", value = "value", -year)

p1 <- ggplot(eff_vacc_cov) +
  geom_line(aes(x=year-1, y = value*100, group = sim, col = "Model")) +
  stat_summary(aes(x=year-1, y = value*100),
               fun="median", geom="line", col = "black", size = 1) +
  geom_point(data=input_who_vaccine_coverage,
             aes(x=year, y = coverage_proportion*100, col = "Data"), size=3) +
  scale_colour_manual("",labels = c("Model" = "Simulated effective vaccine coverage",
                                 "Data" = "Input vaccine coverage"),
                      values = c("Model" = "grey",
                                 "Data" = "red")) +
  scale_x_continuous(breaks=c(1990,2000,2010,2018)) +
  theme_classic() +
  ylab("Percentage immunised among 1-year olds (%)") +
  xlab("Year") +
  ylim(0,100) +
  theme(legend.position=c(0.7,0.2),
        axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14))

p2 <- ggplot(params_mat_accepted_kmeans) +
  geom_boxplot(aes(x="", y=vacc_eff*100), fill = "grey") +
  theme_classic() +
  ylab("Vaccine efficacy (%)") +
  xlab("Posterior") +
  ylim(0,100) +
  theme(legend.position=c(0.7,0.2),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

grid.arrange(p1,p2,widths=c(3,1))
