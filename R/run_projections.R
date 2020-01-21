# Run projections under different scenarios

### Load packages and source file ----
require(here)  # for setting working directory
source(here("R/imperial_model_interventions.R"))

# Simulate intervention model ----

load(here("calibration", "input", "accepted_parmsets_119_060120.Rdata")) # params_mat_targets5

sim <- apply(params_mat_targets5[1,],1,
             function(x)
               run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                         parms_to_change =
                           list(b1 = as.list(x)$b1,
                                b2 = as.list(x)$b2,
                                b3 = as.list(x)$b3,
                                mtct_prob_s = as.list(x)$mtct_prob_s,
                                mtct_prob_e = as.list(x)$mtct_prob_e,
                                alpha = as.list(x)$alpha,
                                p_chronic_in_mtct = as.list(x)$p_chronic_in_mtct,
                                p_chronic_function_r = as.list(x)$p_chronic_function_r,
                                p_chronic_function_s = as.list(x)$p_chronic_function_s,
                                pr_it_ir = as.list(x)$pr_it_ir,
                                pr_ir_ic = as.list(x)$pr_ir_ic,
                                eag_prog_function_rate = as.list(x)$eag_prog_function_rate,
                                pr_ir_enchb = as.list(x)$pr_ir_enchb,
                                pr_ir_cc_female = as.list(x)$pr_ir_cc_female,
                                pr_ir_cc_age_threshold = as.list(x)$pr_ir_cc_age_threshold,
                                pr_ic_enchb = as.list(x)$pr_ic_enchb,
                                sag_loss_slope = as.list(x)$sag_loss_slope,
                                pr_enchb_cc_female = as.list(x)$pr_enchb_cc_female,
                                cirrhosis_male_cofactor = as.list(x)$cirrhosis_male_cofactor,
                                pr_cc_dcc = as.list(x)$pr_cc_dcc,
                                cancer_prog_coefficient_female = as.list(x)$cancer_prog_coefficient_female,
                                cancer_age_threshold = as.list(x)$cancer_age_threshold,
                                cancer_male_cofactor = as.list(x)$cancer_male_cofactor,
                                hccr_it = as.list(x)$hccr_it,
                                hccr_ir = as.list(x)$hccr_ir,
                                hccr_enchb = as.list(x)$hccr_enchb,
                                hccr_cc = as.list(x)$hccr_cc,
                                hccr_dcc = as.list(x)$hccr_dcc,
                                mu_cc = as.list(x)$mu_cc,
                                mu_dcc = as.list(x)$mu_dcc,
                                mu_hcc = as.list(x)$mu_hcc,
                                vacc_eff = as.list(x)$vacc_eff),
                         scenario = "vacc"))


out_0point1_2 <- out

out <- code_model_output(sim[[1]])
outpath <- out

# Proportion in each infection compartment per timestep
plot(outpath$time,outpath$infectioncat_total$carriers/outpath$pop_total$pop_total,type = "l",
     ylim = c(0,0.2), xlim = c(1960,2100))
lines(outpath$time,outpath$infectioncat_total$carriers/outpath$pop_total$pop_total, col = "pink")

outpath <- out_0point5
plot(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births+
        outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections),
     type = "l", ylim = c(0, 4000),
     xlab = "Time", ylab = "New cases of chronic HBV carriage per timestep")
outpath <- out_0point1
lines(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births+
        outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections)*5, col = "red")
outpath <- out_0point1_2
lines(outpath$time,
      (outpath$incident_chronic_infections$horizontal_chronic_infections+
         outpath$incident_chronic_infections$chronic_births+
         outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections)*5, col = "pink")

outpath <- out_0point5
plot(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births+
        outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections)*100000/
       outpath$pop_total$pop_total,
     type = "l", ylim = c(0, 400),
     xlab = "Time", ylab = "New cases of chronic HBV carriage per timestep")
outpath <- out_0point1
lines(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births+
        outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections)*100000*5/
       outpath$pop_total$pop_total, col = "red")
outpath <- out_0point1_2
lines(outpath$time,
      (outpath$incident_chronic_infections$horizontal_chronic_infections+
         outpath$incident_chronic_infections$chronic_births+
         outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections)*100000*5/
        outpath$pop_total$pop_total, col = "pink")

outpath <- out_0point5
plot(outpath$time, outpath$pop_total$pop_total,
     xlab = "Year", ylab = "Population size", type = "l", xlim = c(2000, 2020))
points(popsize_total$time, popsize_total$pop, col = "red")
outpath <- out_0point1
lines(outpath$time, outpath$pop_total$pop_total, col = "red")
outpath <- out_0point_2
lines(outpath$time, outpath$pop_total$pop_total, col = "pink")


#sim_no_vacc <- sim
#sim_vacc_bdvacc_screen <- sim
#save(sim_vacc_bdvacc_screen, file = here("output", "sims_scenario_vacc_bdvacc80_screenadults80_140120.RData"))
#out_vacc_bdvacc_screen <- lapply(sim_vacc_bdvacc_screen,code_model_output)
#save(out_vacc_bdvacc_screen, file = here("output", "sims_output_scenario_vacc_bdvacc_screen_140120.RData"))

# Check for negative numbers => need to switch solver to lsoda!!
o <- sim[[1]]
any(o[,1:(2*n_agecat*n_infectioncat)]<0)

# Test: Calculate all cancer rates

cancer_prog_function <- matrix(0, ncol = 119, nrow = 200)
cancer_prog_female_ic <- matrix(0, ncol = 119, nrow = 200)
cancer_prog_male_ic <- matrix(0, ncol = 119, nrow = 200)
cancer_prog_female_enchb <- matrix(0, ncol = 119, nrow = 200)
cancer_prog_male_enchb <- matrix(0, ncol = 119, nrow = 200)

for (i in c(1:119)) {
  cancer_prog_function[,i] <- (parameters$cancer_prog_coefficient_female[i] * (ages - parameters$cancer_age_threshold[i]))^2  # Rate in females
  cancer_prog_function[,i] <- cancer_prog_function[,i] *
    c(rep(0, times = parameters$cancer_age_threshold[i]/da),
      rep(1, times = n_agecat - parameters$cancer_age_threshold[i]/da))  # Set transition to 0 in <10 year olds
  cancer_prog_female_ic[,i] <- sapply(cancer_prog_function[,i], function(x) min(x,1)) # Set maximum annual rate is 1
  cancer_prog_male_ic[,i] <- sapply(parameters$cancer_male_cofactor[i]*cancer_prog_female_ic[,i], function(x) min(x,1))  # Rate in males, cannot exceed 1
  cancer_prog_female_enchb[,i] <- cancer_prog_female_ic[,i] * parameters$hccr_enchb[i]
  cancer_prog_male_enchb[,i] <- cancer_prog_male_ic[,i] * parameters$hccr_enchb[i]
}

# Rate in 50 year olds
quantile(cancer_prog_female_ic[which(ages == 50),], prob = c(0.025,0.5,0.975))
quantile(cancer_prog_male_ic[which(ages == 50),], prob = c(0.025,0.5,0.975))
quantile(parameter_list$thccr_chb * cancer_prog_female_enchb[which(ages == 50),], prob = c(0.025,0.5,0.975))
quantile(parameter_list$thccr_chb * cancer_prog_male_enchb[which(ages == 50),], prob = c(0.025,0.5,0.975))
quantile(params_mat_targets5$hccr_dcc, prob = c(0.025,0.5,0.975))


# 1) Plots for single iteration  ----
out <- code_model_output(sim[[1]])
outpath <- out

out3 <- code_model_output(sim[[1]])
outpath <- out3

lines(outpath$time,
      (outpath$incident_chronic_infections$horizontal_chronic_infections+
         outpath$incident_chronic_infections$chronic_births+
         outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections)/5,
      col = "blue")

# Total number in each infection compartment per timestep
plot(outpath$time,outpath$infectioncat_total$carriers,type = "l", ylim = c(0,4000000))
lines(outpath$time,outpath$infectioncat_total$sus,col= "red")
lines(outpath$time,outpath$infectioncat_total$immune,col= "blue")

# Proportion in each infection compartment per timestep
plot(outpath$time,outpath$infectioncat_total$carriers/outpath$pop_total$pop_total,type = "l",
     ylim = c(0,0.7), xlim = c(1960,2100))
lines(outpath$time,outpath$infectioncat_total$sus/outpath$pop_total$pop_total,col= "red")
lines(outpath$time,outpath$infectioncat_total$immune/outpath$pop_total$pop_total,col= "blue")

# Number of new cases of chronic HBV carriage at each timestep
plot(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births+
        outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections),
     type = "l", xlim = c(1960,2100), ylim = c(0, 3000),
     xlab = "Time", ylab = "New cases of chronic HBV carriage per timestep")
lines(outpath$time,
      outpath$incident_chronic_infections$chronic_births,
      lty = "dashed")
lines(outpath$time,
      outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections,
      lty = "dashed", col = "red")

# Number of new cases of chronic HBV carriage at each timestep per 100000 people
plot(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births+
        outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections)*100000/
       outpath$pop_total$pop_total,
     type = "l", xlim = c(1960,2100), ylim = c(0, 300),
     xlab = "Time", ylab = "New cases of chronic HBV carriage per timestep")
lines(outpath$time,
      outpath$incident_chronic_infections$chronic_births,
      lty = "dashed")
lines(outpath$time,
      outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections,
      lty = "dashed", col = "red")

# Number of HBV-related deaths at each timestep
plot(outpath$time,
     outpath$hbv_deaths$incident_number_total+outpath$screened_hbv_deaths$incident_number_total+
       outpath$treated_hbv_deaths$incident_number_total,
     type = "l", xlim = c(1960,2100), ylim = c(0, 600),
     xlab = "Time", ylab = "HBV-related deaths per timestep")
lines(outpath$time, outpath$screened_hbv_deaths$incident_number_total,
      lty = "dashed", col = "red")
lines(outpath$time, outpath$treated_hbv_deaths$incident_number_total,
      lty = "dashed", col = "pink")

# Number of HBV-related deaths at each timestep per 100000 people
plot(x=outpath$time,
     y=(outpath$hbv_deaths$incident_number_total+outpath$screened_hbv_deaths$incident_number_total+
          outpath$treated_hbv_deaths$incident_number_total)*100000/outpath$pop_total$pop_total,
     type = "l", xlim = c(1960,2100), ylim = c(0, 20),
     xlab = "Time", ylab = "HBV-related deaths per timestep per 100000 people")
lines(outpath$time, outpath$screened_hbv_deaths$incident_number_total,
      lty = "dashed", col = "red")
lines(outpath$time, outpath$treated_hbv_deaths$incident_number_total,
      lty = "dashed", col = "pink")

lines(x=outpath$time,
      y=(outpath$hbv_deaths$incident_number_total+outpath$screened_hbv_deaths$incident_number_total+
           outpath$treated_hbv_deaths$incident_number_total)*100000/outpath$pop_total$pop_total, col = "blue")

unscreen_pop <- rowSums(sim[[1]][,c(2:1801, 4602:6401)])
screen_pop <- rowSums(select(sim[[1]], starts_with("S_")))
treat_pop <- rowSums(select(sim[[1]], starts_with("T_")))
plot(x=outpath$time, y = unscreen_pop)
plot(x=outpath$time, y = screen_pop)
plot(x=outpath$time, y = treat_pop)

deaths1 <- rowSums(outpath$out_cum_hbv_deathsf)+rowSums(outpath$out_cum_hbv_deathsm)
deaths2 <- rowSums(outpath$out_cum_screened_hbv_deathsf)+rowSums(outpath$out_cum_screened_hbv_deathsm)
deaths3 <- rowSums(outpath$out_cum_treated_hbv_deathsf)+rowSums(outpath$out_cum_treated_hbv_deathsm)
plot(out$time,deaths1+deaths2+deaths3, type = "l", xlim=c(2015,2100))
lines(out$time,deaths1+deaths2+deaths3, col = "blue")
(deaths1+deaths2+deaths3)[which(out$time==2020.5)]
(deaths1+deaths2+deaths3)[which(out$time==2030)]
(deaths1+deaths2+deaths3)[which(out$time==2050)]

plot(out$time[1:499], c(0, diff(deaths1[1:which(out$time == 2020)], lag=1),
                        diff(deaths1[which(out$time ==2020.5):which(out$time == 2099.5)], lag=1)),
     ylim = c(0,1000))
plot(out$time[1:499], c(0, diff(deaths2[1:which(out$time == 2020)], lag=1),
                        diff(deaths2[which(out$time ==2020.5):which(out$time == 2099.5)], lag=1)),
     ylim = c(0,100))
diff(deaths2[1:which(out$time == 2020)], lag=1)
diff(deaths2[which(out$time ==2020.5):which(out$time == 2099.5)], lag=1)
round(diff(deaths2, lag=1),2)
head(deaths1)

plot(x=out$time, y=rowSums(outpath$out_cum_hbv_deathsf)+rowSums(outpath$out_cum_hbv_deathsm), type = "l")
lines(x=out$time, y=rowSums(outpath$out_cum_screened_hbv_deathsf)+rowSums(outpath$out_cum_screened_hbv_deathsm), type = "l", col = "blue")
plot(x=out$time, y=rowSums(outpath$out_cum_screened_hbv_deathsf)+rowSums(outpath$out_cum_screened_hbv_deathsm), type = "l")

plot(x=outpath$time, y = outpath$hbv_deaths$incident_number_total*100000/
       unscreen_pop,ylim=c(0,30))
plot(x=outpath$time, y = outpath$screened_hbv_deaths$incident_number_total*100000/
       screen_pop,ylim=c(0,12))
plot(x=outpath$time[outpath$time>2015], y = outpath$treated_hbv_deaths$incident_number_total[outpath$time>2015]*100000/
       treat_pop_2015plus,ylim=c(0,12))
plot(x=outpath$time[outpath$time>2015],
     y = (outpath$hbv_deaths$incident_number_total+outpath$screened_hbv_deaths$incident_number_total+
            outpath$treated_hbv_deaths$incident_number_total)[outpath$time>2015]*100000/
       (unscreen_pop_2015plus+screen_pop_2015plus+treat_pop_2015plus),ylim=c(0,12))


# Plot total population size over timesteps
plot(outpath$time, outpath$pop_total$pop_total,
     xlab = "Year", ylab = "Population size", type = "l", ylim = c(0,7000000))
points(popsize_total$time, popsize_total$pop, col = "red")

# Plot total number of births over time periods
plot(outpath$births_group5$timegroup, outpath$births_group5$incident_number,
     xlab = "5-year time periods", ylab = "Total number of births",
     type = "l", ylim = c(0,600000))
points(x = as.numeric(strtrim(input_births_clean$time, width = 4)),
       y = input_births_clean$births, col = "red")

# Plot total number of deaths over time periods
plot(x = outpath$deaths_total_group5$timegroup,
     y = outpath$deaths_total_group5$total,
     ylim = c(0,400000), type = "l", xlab = "5-year time periods", ylab = "Total number of deaths")
points(x = as.numeric(strtrim(deaths_total$time, width = 4)),
       y = deaths_total$deaths, col = "red")


# 2) Plots for multiple parameter sets ----

#load(here("output", "sims_scenario_vacc_120120.RData"))   # sim
#load(here("output", "sims_scenario_no_vacc_120120.RData")) # sim_no_vacc
#out_no_vacc <- lapply(sim_no_vacc,code_model_output)
#out_vacc <- out
#save(out_vacc, file = here("output", "sims_output_scenario_vacc_130120.RData"))
#save(out_no_vacc, file = here("output", "sims_output_scenario_no_vacc_130120.RData"))
#out_vacc_bdvacc <- lapply(sim_vacc_bdvacc,code_model_output)
#save(out_vacc_bdvacc, file = here("output", "sims_output_scenario_vacc_bdvacc_140120.RData"))
# Load the outputs from a scenario with and without vaccine
# THESE HAVE MOVED - change path
#load(here("output", "sims_output_scenario_vacc_130120.RData"))
#load(here("output", "sims_output_scenario_no_vacc_130120.RData"))

# No historical intervention scenario: out_no_vacc

# HBsAg prevalence
proj_prev_no_vacc <- cbind(out_no_vacc[[1]]$time,
                           (sapply(lapply(out_no_vacc,"[[", "infectioncat_total"), "[[", "carriers")/
                              sapply(lapply(out_no_vacc,"[[", "pop_total"), "[[", "pop_total")))
colnames(proj_prev_no_vacc)[1] <- "time"
proj_prev_summary_no_vacc <- data.frame(time = out_no_vacc[[1]]$time)
proj_prev_summary_no_vacc$median <- apply(proj_prev_no_vacc[,-1],1,median)
proj_prev_summary_no_vacc$lower <- apply(proj_prev_no_vacc[,-1],1,quantile, prob = 0.025)
proj_prev_summary_no_vacc$upper <- apply(proj_prev_no_vacc[,-1],1,quantile, prob = 0.975)
proj_prev_no_vacc_long <- gather(as.data.frame(proj_prev_no_vacc), key = "iteration", value = "hbsag_prev", -time)

# Number of new cases of chronic HBV carriage per 100000 per timestep
proj_inc_no_vacc <- cbind(out_no_vacc[[1]]$time,
                          (sapply(lapply(out_no_vacc,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                             sapply(lapply(out_no_vacc, "[[", "incident_chronic_infections"), "[[", "chronic_births")+
                             sapply(lapply(out_no_vacc,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections"))*100000/
                            sapply(lapply(out_no_vacc,"[[", "pop_total"), "[[", "pop_total"))
colnames(proj_inc_no_vacc)[1] <- "time"
proj_inc_summary_no_vacc <- data.frame(time = out_no_vacc[[1]]$time)
proj_inc_summary_no_vacc$median <- apply(proj_inc_no_vacc[,-1],1,median)
proj_inc_summary_no_vacc$lower <- apply(proj_inc_no_vacc[,-1],1,quantile, prob = 0.025)
proj_inc_summary_no_vacc$upper <- apply(proj_inc_no_vacc[,-1],1,quantile, prob = 0.975)
proj_inc_no_vacc_long <- gather(as.data.frame(proj_inc_no_vacc), key = "iteration", value = "chronic_cases_per_100000", -time)

# Incident HBV-related deaths per 100000 per timestep
proj_deaths_no_vacc <- cbind(out_no_vacc[[1]]$time,
                             (sapply(lapply(out_no_vacc,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                                sapply(lapply(out_no_vacc,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                                sapply(lapply(out_no_vacc,"[[", "treated_hbv_deaths"), "[[", "incident_number_total"))*100000/
                               sapply(lapply(out_no_vacc,"[[", "pop_total"), "[[", "pop_total"))
colnames(proj_deaths_no_vacc)[1] <- "time"
proj_deaths_summary_no_vacc <- data.frame(time = out_no_vacc[[1]]$time)
proj_deaths_summary_no_vacc$median <- apply(proj_deaths_no_vacc[,-1],1,median)
proj_deaths_summary_no_vacc$lower <- apply(proj_deaths_no_vacc[,-1],1,quantile, prob = 0.025)
proj_deaths_summary_no_vacc$upper <- apply(proj_deaths_no_vacc[,-1],1,quantile, prob = 0.975)
proj_deaths_no_vacc_long <- gather(as.data.frame(proj_deaths_no_vacc), key = "iteration", value = "deaths_per_100000", -time)

# HCC incidence per 100000 per timestep
proj_hcc_no_vacc <- cbind(out_no_vacc[[1]]$time,
                          (sapply(lapply(out_no_vacc,"[[", "incident_hcc"), "[[", "incident_number_total")*100000/
                             sapply(lapply(out_no_vacc,"[[", "pop_total"), "[[", "pop_total")))
colnames(proj_hcc_no_vacc)[1] <- "time"
proj_hcc_summary_no_vacc <- data.frame(time = out_no_vacc[[1]]$time)
proj_hcc_summary_no_vacc$median <- apply(proj_hcc_no_vacc[,-1],1,median)
proj_hcc_summary_no_vacc$lower <- apply(proj_hcc_no_vacc[,-1],1,quantile, prob = 0.025)
proj_hcc_summary_no_vacc$upper <- apply(proj_hcc_no_vacc[,-1],1,quantile, prob = 0.975)

# Label dataframes with scenario:
proj_prev_summary_no_vacc$scenario <- "no_vacc"
proj_inc_summary_no_vacc$scenario <- "no_vacc"
proj_deaths_summary_no_vacc$scenario <- "no_vacc"
proj_hcc_summary_no_vacc$scenario <- "mo_vacc"

# Vaccine scenario (status quo): out_vacc
out_vacc <- out_vacc_bdvacc

proj_prev_vacc <- cbind(out_vacc[[1]]$time,
                        (sapply(lapply(out_vacc,"[[", "infectioncat_total"), "[[", "carriers")/
                           sapply(lapply(out_vacc,"[[", "pop_total"), "[[", "pop_total")))
colnames(proj_prev_vacc)[1] <- "time"
proj_prev_summary_vacc <- data.frame(time = out_vacc[[1]]$time)
proj_prev_summary_vacc$median <- apply(proj_prev_vacc[,-1],1,median)
proj_prev_summary_vacc$lower <- apply(proj_prev_vacc[,-1],1,quantile, prob = 0.025)
proj_prev_summary_vacc$upper <- apply(proj_prev_vacc[,-1],1,quantile, prob = 0.975)

proj_deaths_vacc <- cbind(out_vacc[[1]]$time,
                          (sapply(lapply(out_vacc,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                             sapply(lapply(out_vacc,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                             sapply(lapply(out_vacc,"[[", "treated_hbv_deaths"), "[[", "incident_number_total"))*100000/
                            sapply(lapply(out_vacc,"[[", "pop_total"), "[[", "pop_total"))
colnames(proj_deaths_vacc)[1] <- "time"
proj_deaths_summary_vacc <- data.frame(time = out_vacc[[1]]$time)
proj_deaths_summary_vacc$median <- apply(proj_deaths_vacc[,-1],1,median)
proj_deaths_summary_vacc$lower <- apply(proj_deaths_vacc[,-1],1,quantile, prob = 0.025)
proj_deaths_summary_vacc$upper <- apply(proj_deaths_vacc[,-1],1,quantile, prob = 0.975)

proj_inc_vacc <- cbind(out_vacc[[1]]$time,
                       (sapply(lapply(out_vacc,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                          sapply(lapply(out_vacc,"[[", "incident_chronic_infections"), "[[", "chronic_births")+
                          sapply(lapply(out_vacc,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections"))*100000/
                         sapply(lapply(out_vacc,"[[", "pop_total"), "[[", "pop_total"))
colnames(proj_inc_vacc)[1] <- "time"
proj_inc_summary_vacc <- data.frame(time = out_vacc[[1]]$time)
proj_inc_summary_vacc$median <- apply(proj_inc_vacc[,-1],1,median)
proj_inc_summary_vacc$lower <- apply(proj_inc_vacc[,-1],1,quantile, prob = 0.025)
proj_inc_summary_vacc$upper <- apply(proj_inc_vacc[,-1],1,quantile, prob = 0.975)

proj_hcc_vacc <- cbind(out_vacc[[1]]$time,
                       (sapply(lapply(out_vacc,"[[", "incident_hcc"), "[[", "incident_number_total")*100000/
                          sapply(lapply(out_vacc,"[[", "pop_total"), "[[", "pop_total")))
colnames(proj_hcc_vacc)[1] <- "time"
proj_hcc_summary_vacc <- data.frame(time = out_vacc[[1]]$time)
proj_hcc_summary_vacc$median <- apply(proj_hcc_vacc[,-1],1,median)
proj_hcc_summary_vacc$lower <- apply(proj_hcc_vacc[,-1],1,quantile, prob = 0.025)
proj_hcc_summary_vacc$upper <- apply(proj_hcc_vacc[,-1],1,quantile, prob = 0.975)

# Label dataframes with scenario:
proj_prev_summary_vacc$scenario <- "vacc"
proj_inc_summary_vacc$scenario <- "vacc"
proj_deaths_summary_vacc$scenario <- "vacc"
proj_hcc_summary_vacc$scenario <- "vacc"

proj_prev_summary_bdvacc <- proj_prev_summary_vacc
proj_inc_summary_bdvacc <- proj_inc_summary_vacc
proj_deaths_summary_bdvacc <- proj_deaths_summary_vacc
proj_hcc_summary_bdvacc <- proj_hcc_summary_vacc

proj_prev_summary_bdvacc$scenario <- "vacc_bdvacc"
proj_inc_summary_bdvacc$scenario <- "vacc_bdvacc"
proj_deaths_summary_bdvacc$scenario <- "vacc_bdvacc"
proj_hcc_summary_bdvacc$scenario <- "vacc_bdvacc"

# 2a) PLOT SINGLE SCENARIO ----

out <- out_vacc_bdvacc

# HBsAg prevalence
proj_prev <- cbind(out[[1]]$time,
                   (sapply(lapply(out,"[[", "infectioncat_total"), "[[", "carriers")/
                      sapply(lapply(out,"[[", "pop_total"), "[[", "pop_total")))
colnames(proj_prev)[1] <- "time"
proj_prev_summary <- data.frame(time = out[[1]]$time)
proj_prev_summary$median <- apply(proj_prev[,-1],1,median)
proj_prev_summary$lower <- apply(proj_prev[,-1],1,quantile, prob = 0.025)
proj_prev_summary$upper <- apply(proj_prev[,-1],1,quantile, prob = 0.975)
proj_prev_long <- gather(as.data.frame(proj_prev), key = "iteration", value = "hbsag_prev", -time)

ggplot(proj_prev_long) +
  geom_line(aes(x=time, y = hbsag_prev*100, group = iteration), col = "grey90")+
  geom_line(data = proj_prev_summary, aes(x = time, y = median*100), col = "red") +
  geom_line(data = proj_prev_summary, aes(x = time, y = lower*100), col = "blue") +
  geom_line(data = proj_prev_summary, aes(x = time, y = upper*100), col = "blue") +
  xlim(1960,2100) +
  labs(title = "HBsAg prevalence (%)") +
  ylim(0,25) +
  theme_bw()

# HBV-related mortality per 100000 per timestep (0.5 years)
proj_deaths <- cbind(out[[1]]$time,
                     (sapply(lapply(out,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                        sapply(lapply(out,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                        sapply(lapply(out,"[[", "treated_hbv_deaths"), "[[", "incident_number_total"))*100000/
                       sapply(lapply(out,"[[", "pop_total"), "[[", "pop_total"))
colnames(proj_deaths)[1] <- "time"
proj_deaths_summary <- data.frame(time = out[[1]]$time)
proj_deaths_summary$median <- apply(proj_deaths[,-1],1,median)
proj_deaths_summary$lower <- apply(proj_deaths[,-1],1,quantile, prob = 0.025)
proj_deaths_summary$upper <- apply(proj_deaths[,-1],1,quantile, prob = 0.975)
proj_deaths_long <- gather(as.data.frame(proj_deaths), key = "iteration", value = "deaths_per_100000", -time)

ggplot(proj_deaths_long) +
  geom_line(aes(x=time, y = deaths_per_100000, group = iteration), col = "grey90")+
  geom_line(data = proj_deaths_summary, aes(x = time, y = median), col = "red") +
  geom_line(data = proj_deaths_summary, aes(x = time, y = lower), col = "blue") +
  geom_line(data = proj_deaths_summary, aes(x = time, y = upper), col = "blue") +
  labs(title = "HBV-related deaths per 100,000 people per 6-month timestep") +
  xlim(1960,2100) +
  ylim(0,50) +
  theme_bw()

# Age-standardised HBV-related mortality rate per 100000 per timestep (0.5 years)
# a) Calculate crude age-specific rates per person-year: need age-specific number of deaths and age-specific population size
# b) Multiply by reference population (Gambian pop in 2020)
# c) Calculate total expected deaths (sum of age-specific numbers) and divide by total Gambian popsize in 2020

# Number of new cases of chronic HBV carriage per 100000 per timestep
proj_inc <- cbind(out[[1]]$time,
                  (sapply(lapply(out,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                     sapply(lapply(out,"[[", "incident_chronic_infections"), "[[", "chronic_births")+
                     sapply(lapply(out,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections"))*100000/
                    sapply(lapply(out,"[[", "pop_total"), "[[", "pop_total"))
colnames(proj_inc)[1] <- "time"
proj_inc_summary <- data.frame(time = out[[1]]$time)
proj_inc_summary$median <- apply(proj_inc[,-1],1,median)
proj_inc_summary$lower <- apply(proj_inc[,-1],1,quantile, prob = 0.025)
proj_inc_summary$upper <- apply(proj_inc[,-1],1,quantile, prob = 0.975)
proj_inc_long <- gather(as.data.frame(proj_inc), key = "iteration", value = "chronic_cases_per_100000", -time)

ggplot(proj_inc_long) +
  geom_line(aes(x=time, y = chronic_cases_per_100000, group = iteration), col = "grey90")+
  geom_line(data = proj_inc_summary, aes(x = time, y = median), col = "red") +
  geom_line(data = proj_inc_summary, aes(x = time, y = lower), col = "blue") +
  geom_line(data = proj_inc_summary, aes(x = time, y = upper), col = "blue") +
  xlim(1960,2100) +
  labs(title = "New cases of chronic HBV carriage per 100,000 people per 6-month timestep") +
  theme_bw()

# Absolute number of new cases of chronic HBV carriage per timestep
proj_inc2 <- cbind(out[[1]]$time,
                   (sapply(lapply(out,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                      sapply(lapply(out,"[[", "incident_chronic_infections"), "[[", "chronic_births")+
                      sapply(lapply(out,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections")))
colnames(proj_inc2)[1] <- "time"
proj_inc2_summary <- data.frame(time = out[[1]]$time)
proj_inc2_summary$median <- apply(proj_inc2[,-1],1,median)
proj_inc2_summary$lower <- apply(proj_inc2[,-1],1,quantile, prob = 0.025)
proj_inc2_summary$upper <- apply(proj_inc2[,-1],1,quantile, prob = 0.975)
proj_inc2_long <- gather(as.data.frame(proj_inc2), key = "iteration", value = "chronic_cases", -time)

ggplot(proj_inc2_long) +
  geom_line(aes(x=time, y = chronic_cases, group = iteration), col = "grey90")+
  geom_line(data = proj_inc2_summary, aes(x = time, y = median), col = "red") +
  geom_line(data = proj_inc2_summary, aes(x = time, y = lower), col = "blue") +
  geom_line(data = proj_inc2_summary, aes(x = time, y = upper), col = "blue") +
  xlim(1960,2100) +
  labs(title = "New cases of chronic HBV carriage per 6-month timestep") +
  theme_bw()

# HCC incidence per 100000 per timestep (0.5 years)
# Note this needs to be changed to include with treatment
proj_hcc <- cbind(out[[1]]$time,
                  (sapply(lapply(out,"[[", "incident_hcc"), "[[", "incident_number_total")*100000/
                     sapply(lapply(out,"[[", "pop_total"), "[[", "pop_total")))
colnames(proj_hcc)[1] <- "time"
proj_hcc_summary <- data.frame(time = out[[1]]$time)
proj_hcc_summary$median <- apply(proj_hcc[,-1],1,median)
proj_hcc_summary$lower <- apply(proj_hcc[,-1],1,quantile, prob = 0.025)
proj_hcc_summary$upper <- apply(proj_hcc[,-1],1,quantile, prob = 0.975)
proj_hcc_long <- gather(as.data.frame(proj_hcc), key = "iteration", value = "hcc_cases_per_100000", -time)

ggplot(proj_hcc_long) +
  geom_line(aes(x=time, y = hcc_cases_per_100000, group = iteration), col = "grey90")+
  geom_line(data = proj_hcc_summary, aes(x = time, y = median), col = "red") +
  geom_line(data = proj_hcc_summary, aes(x = time, y = lower), col = "blue") +
  geom_line(data = proj_hcc_summary, aes(x = time, y = upper), col = "blue") +
  labs(title = "HBV-related HCC cases per 100,000 people per 6-month timestep") +
  xlim(1960,2100) +
  ylim(0,20) +
  theme_bw()

# 2b) PLOTS FOR NO HISTORICAL INTERVENTION/STATUS QUO ----

# Prepare output data to plot
#load(here("output", "sims_output_scenario_vacc_130120.RData"))  # out_vacc
out_plot_vacc <- extract_outcomes(output_file = out_vacc, scenario_label = "vacc")
rm(out_vacc)
gc()
#load(here("output", "sims_output_scenario_no_vacc_130120.RData")) # out_no_vacc
out_plot_no_vacc <- extract_outcomes(output_file = out_no_vacc, scenario_label = "no_vacc")
rm(out_no_vacc)
gc()

# Combine projection summaries from different scenarios
proj_prev_total <- rbind(out_plot_vacc$hbsag_prev_summary, out_plot_no_vacc$hbsag_prev_summary)
proj_inc_total <- rbind(out_plot_vacc$chronic_incidence_summary, out_plot_no_vacc$chronic_incidence_summary)
proj_deaths_total <- rbind(out_plot_vacc$hbv_deaths_summary, out_plot_no_vacc$hbv_deaths_summary)
proj_deaths_standard_rate_total <- rbind(out_plot_vacc$proj_deaths_standardised_summary,
                                         out_plot_no_vacc$proj_deaths_standardised_summary)

#proj_hcc_total <- rbind(proj_hcc_summary, proj_hcc_summary_no_vacc_no_vacc)

# Prevalence
ggplot(proj_prev_total) +
  geom_line(aes(x=time, y = median*100, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower*100, ymax=upper*100, group = scenario, fill = scenario), alpha = 0.1)+
  #  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Effect of infant vaccination on HBsAg prevalence",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Prevalence (%)")+
  xlab("Year")+
  ylim(0,25) +
  theme_classic()

# Chronic HBV incidence (absolute)
ggplot(proj_inc_total) +
  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, group = scenario, fill = scenario), alpha = 0.1)+
  #  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Effect of infant vaccination on number of new chronic HBV infections",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Incident cases of chronic HBV infection per 6 months")+
  xlab("Year")+
  #  ylim(0,25) +
  theme_classic()

# HBV-related deaths (absolute)
ggplot(proj_deaths_total) +
  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, group = scenario, fill = scenario), alpha = 0.1)+
  # geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Effect of infant vaccination on number of HBV-related deaths",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Deaths per 6 months")+
  xlab("Year")+
  ylim(0,3000) +
  theme_classic()

# HBV-related age-standardised death rate per 100000
ggplot(proj_deaths_standard_rate_total) +
  geom_line(aes(x=time, y = median*100000, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower*100000, ymax=upper*100000, group = scenario, fill = scenario), alpha = 0.1)+
  # geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Effect of infant vaccination on age-standardised HBV-related death rate",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Age-standardised death rate per 100,000 people per 6 months")+
  xlab("Year")+
  ylim(0,45) +
  theme_classic()

# HCC incidence
#ggplot(proj_hcc_total) +
#  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
#  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, group = scenario, fill = scenario), alpha = 0.1)+
#  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
#  labs(title = "Effect of infant vaccination on HBV-related HCC incidence",
#       colour = "Modelled scenario", fill = "Modelled scenario",
#       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
#  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
#  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
#  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
#  ylab("Incident HCC cases per 100,000 people per 6 months")+
#  xlab("Year")+
#  ylim(0,15) +
#  theme_classic()

#save(proj_prev_total, proj_inc_total, proj_deaths_total, proj_hcc_total, file = here("output", "elimination_proj_2_scenarios_130120.RData"))

# 2c) PLOTS FOR STATUS QUO/OTHER INTERVENTIONS ----

# Prepare output data to plot
# Status quo scenario: continue infant vaccine at 93% coverage
#load(here("output", "sims_output_scenario_vacc_130120.RData"))  # out_vacc
out_plot_vacc <- extract_outcomes(output_file = out_vacc, scenario_label = "vacc")
rm(out_vacc)
gc()
# Infant + birth dose vaccine: introduce in 2020 at 80% coverage
#load(here("output", "sims_output_scenario_vacc_bdvacc_140120.RData")) # out_vacc_bdvacc
out_plot_bdvacc <- extract_outcomes(output_file = out_vacc_bdvacc, scenario_label = "bdvacc")
rm(out_vacc_bdvacc)
gc()
# Infant + birth dose vaccine (introduce in 2020 at 80% coverage) + screean-and-treat (one-off screen of 80% of 80+ adults)
#load(here("output", "sims_output_scenario_vacc_bdvacc_screen_140120.RData")) # out_vacc_bdvacc_screen
out_plot_bdvacc_screen <- extract_outcomes(output_file = out_vacc_bdvacc_screen, scenario_label = "bdvacc_screen")
rm(out_vacc_bdvacc_screen)
gc()

# Combine projection summaries from different scenarios
proj_prev_total <- rbind(out_plot_vacc$hbsag_prev_summary, out_plot_bdvacc$hbsag_prev_summary,
                         out_plot_bdvacc_screen$hbsag_prev_summary)
proj_inc_total <- rbind(out_plot_vacc$chronic_incidence_summary, out_plot_bdvacc$chronic_incidence_summary,
                        out_plot_bdvacc_screen$chronic_incidence_summary)
proj_deaths_total <- rbind(out_plot_vacc$hbv_deaths_summary, out_plot_bdvacc$hbv_deaths_summary,
                           out_plot_bdvacc_screen$hbv_deaths_summary)
proj_deaths_standard_rate_total <- rbind(out_plot_vacc$proj_deaths_standardised_summary,
                                         out_plot_bdvacc$proj_deaths_standardised_summary,
                                         out_plot_bdvacc_screen$proj_deaths_standardised_summary)

# Create plots

# HBsAg prevalence
ggplot(proj_prev_total) +
  geom_line(aes(x=time, y = median*100, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower*100, ymax=upper*100, group = scenario, fill = scenario), alpha = 0.1)+
  #  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "HBsAg prevalence",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical infant vaccine coverage since 1990 and maintains 93% coverage after 2018\n
       **Birth dose vaccine introduced in 2020 with 80% coverage\n
       ***Birth dose as in ** and one-off screening+continuous treatment of 80% of adults in 2020") +
  scale_color_manual(labels = c("vacc"="Status quo infant vaccine*",
                                "bdvacc"="Infant+birth dose vaccine**",
                                "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
                     values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_fill_manual(labels = c("vacc"="Status quo infant vaccine*",
                               "bdvacc"="Infant+birth dose vaccine**",
                               "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
                    values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Prevalence (%)")+
  xlab("Year")+
  ylim(0,25) +
  theme_classic()

# Chronic HBV incidence (absolute number)
ggplot(proj_inc_total) +
  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, group = scenario, fill = scenario), alpha = 0.1)+
  #  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Chronic HBV infection incidence",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical infant vaccine coverage since 1990 and maintains 93% coverage after 2018\n
       **Birth dose vaccine introduced in 2020 with 80% coverage\n
       ***Birth dose as in ** and one-off screening+continuous treatment of 80% of adults in 2020") +
  scale_color_manual(labels = c("vacc"="Status quo infant vaccine*",
                                "bdvacc"="Infant+birth dose vaccine**",
                                "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
                     values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_fill_manual(labels = c("vacc"="Status quo infant vaccine*",
                               "bdvacc"="Infant+birth dose vaccine**",
                               "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
                    values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Incident cases of chronic HBV infection per 6 months")+
  xlab("Year")+
  ylim(0,7000) +
  theme_classic()

# HBV-related deaths
ggplot(proj_deaths_total) +
  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, group = scenario, fill = scenario), alpha = 0.1)+
  # geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "HBV-related deaths",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical infant vaccine coverage since 1990 and maintains 93% coverage after 2018\n
       **Birth dose vaccine introduced in 2020 with 80% coverage\n
       ***Birth dose as in ** and one-off screening+continuous treatment of 80% of adults in 2020") +
  scale_color_manual(labels = c("vacc"="Status quo infant vaccine*",
                                "bdvacc"="Infant+birth dose vaccine**",
                                "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
                     values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_fill_manual(labels = c("vacc"="Status quo infant vaccine*",
                               "bdvacc"="Infant+birth dose vaccine**",
                               "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
                    values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Deaths per 6 months")+
  xlab("Year")+
  ylim(0,1000) +
  theme_classic()


# HBV-related age-standardised death rate per 100000
ggplot(proj_deaths_standard_rate_total) +
  geom_line(aes(x=time, y = median*100000, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower*100000, ymax=upper*100000, group = scenario, fill = scenario), alpha = 0.1)+
  # geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Age-standardised HBV-related death rate",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical infant vaccine coverage since 1990 and maintains 93% coverage after 2018\n
       **Birth dose vaccine introduced in 2020 with 80% coverage\n
       ***Birth dose as in ** and one-off screening+continuous treatment of 80% of adults in 2020") +
  scale_color_manual(labels = c("vacc"="Status quo infant vaccine*",
                                "bdvacc"="Infant+birth dose vaccine**",
                                "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
                     values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_fill_manual(labels = c("vacc"="Status quo infant vaccine*",
                               "bdvacc"="Infant+birth dose vaccine**",
                               "bdvacc_screen"="Infant+birth dose vaccine\nand screening***"),
                    values = c("bdvacc"="turquoise", "vacc"="deeppink", "bdvacc_screen"="orange")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Age-standardised death rate per 100,000 people per 6 months")+
  xlab("Year")+
  ylim(0,45) +
  theme_classic()


