# Plots for check of proportion of new chronic infections due to MTCT

# Read in chronic_infection_incidence_status_quo_250920.rds

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
