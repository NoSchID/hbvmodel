# Run projections under different scenarios

### Load packages and source file ----
require(here)  # for setting working directory
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Simulate intervention model ----

#load(here("calibration", "input", "accepted_parmsets_123_180520.Rdata")) # params_mat_targets5
load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans
load(here("analysis_input", "scenario_a1_it_parms.Rdata"))
load(here("analysis_input", "scenario_a1_parms.Rdata"))
load(here("analysis_input", "scenario_anc1_it_parms.Rdata"))
load(here("analysis_input", "scenario_wpl1_parms.Rdata"))

sim <- apply(params_mat_accepted_kmeans[2,],1,
             function(x)
               run_model(sim_duration = runtime, default_parameter_list =scenario_a1_it_parms,
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
                                vacc_eff = as.list(x)$vacc_eff,
                                screening_years = c(2020),
                                #screening_coverage = 0.9,
                                #apply_treat_it = 1,
                                #prop_negative_to_remove_from_rescreening = 0,
                                apply_screen_not_treat = 0,
                                #monitoring_rate = monit_rate_vec,
                                monitoring_rate = c(rep(0, length(which(ages == 0):which(ages == 15-da))),
                                                    rep(1/10, length(which(ages == 15):which(ages == 45-da))),
                                                    rep(0, length(which(ages == 45):which(ages == 100-da)))),
                                apply_repeat_screen = 0,
                                #apply_lifetime_monitoring = 1,
                                #monitoring_prob = 1,
                                #treatment_initiation_prob = 1,
                                #lifetime_monitoring_event_rate = monit_rate_vec,
                                #min_age_to_screen = 15,
                                #max_age_to_screen = 65,
                                #min_age_to_repeat_screen = 15,
                                #max_age_to_repeat_screen = 49.5,
                                repeat_screening_years = seq(2020.5,2040,0.5)),
                                drop_timesteps_before = 1960,
                         scenario = "vacc_screen"))

out3 <- code_model_output(sim[[1]])
outpath <- out


load(here("output", "sims_output_scenario_vacc_130120.RData"))
out <- out_vacc
rm(out_vacc)
gc()

out <- lapply(out, function(x) x[!(names(x) %in% c("births_group5", "deaths_total_group5"))])

format(object.size(out), units = "Gb")
# First turn time into dataframe object
out_sub <- lapply(out, function(x) {x$time <- as.data.frame(x$time) ; x})
out_sub <- lapply(out_sub, function(x) x[-c(length(x))])

tic()
for (i in 1:length(out_sub)) {
  out_sub[[i]] <- lapply(out_sub[[i]], function(x) x[-c(1:221),])
}
toc()

out_sub <- lapply(out_sub, function(x) lapply(x, function(y) y[-c(1:221),]))
out_sub <- lapply(out_sub, function(x) {x$time <- as.numeric(x$time) ; x})
final_list <- Map(c, out_sub, inparms)
format(object.size(final_list), units = "Mb")


# Try to save population output in 1955
sim <- sim[[1]]$out
model_pop1955 <- sim[which(sim$time==1955),1:(2*n_infectioncat*n_agecat)+1]
#save(here("output/simulated_inits_1955.RData"))
# Load initial population saved from previous model run
# Check numbers in screened and treatment comps are really 0
# Check if output storage needs to be taken over too
#load(here("data/simulated_inits_1880.RData"))  # this is saved from previous model run
init_pop_sim <- c("Sf" = select(model_pop1955, starts_with("Sf")),
                  "ITf" = select(model_pop1955, starts_with("ITf")),
                  "IRf" = select(model_pop1955, starts_with("IRf")),
                  "ICf" = select(model_pop1955, starts_with("ICf")),
                  "ENCHBf" = select(model_pop1955, starts_with("ENCHBf")),
                  "CCf" = select(model_pop1955, starts_with("CCf")),
                  "DCCf" = select(model_pop1955, starts_with("DCCf")),
                  "HCCf" = select(model_pop1955, starts_with("HCCf")),
                  "Rf" = select(model_pop1955, starts_with("Rf")),
                  "S_Sf" = rep(0,n_agecat),
                  "S_ITf" = rep(0,n_agecat),
                  "S_IRf" = rep(0,n_agecat),
                  "S_ICf" = rep(0,n_agecat),
                  "S_ENCHBf" = rep(0,n_agecat),
                  "S_CCf" = rep(0,n_agecat),
                  "S_DCCf" = rep(0,n_agecat),
                  "S_HCCf" = rep(0,n_agecat),
                  "S_Rf" = rep(0,n_agecat),
                  "T_CHBf" = rep(0,n_agecat),
                  "T_CCf" = rep(0,n_agecat),
                  "T_DCCf" = rep(0,n_agecat),
                  "T_HCCf" = rep(0,n_agecat),
                  "T_Rf" = rep(0,n_agecat),
                  "Sm" = select(model_pop1955, starts_with("Sm")),
                  "ITm" = select(model_pop1955, starts_with("ITm")),
                  "IRm" = select(model_pop1955, starts_with("IRm")),
                  "ICm" = select(model_pop1955, starts_with("ICm")),
                  "ENCHBm" = select(model_pop1955, starts_with("ENCHBm")),
                  "CCm" = select(model_pop1955, starts_with("CCm")),
                  "DCCm" = select(model_pop1955, starts_with("DCCm")),
                  "HCCm" = select(model_pop1955, starts_with("HCCm")),
                  "Rm" = select(model_pop1955, starts_with("Rm")),
                  "S_Sm" = rep(0,n_agecat),
                  "S_ITm" = rep(0,n_agecat),
                  "S_IRm" = rep(0,n_agecat),
                  "S_ICm" = rep(0,n_agecat),
                  "S_ENCHBm" = rep(0,n_agecat),
                  "S_CCm" = rep(0,n_agecat),
                  "S_DCCm" = rep(0,n_agecat),
                  "S_HCCm" = rep(0,n_agecat),
                  "S_Rm" = rep(0,n_agecat),
                  "T_CHBm" = rep(0,n_agecat),
                  "T_CCm" = rep(0,n_agecat),
                  "T_DCCm" = rep(0,n_agecat),
                  "T_HCCm" = rep(0,n_agecat),
                  "T_Rm" =rep(0,n_agecat),
                  output_storage)
init_pop_sim <- unlist(init_pop_sim)

from_1850 <- out
from_1955 <- out
sim2_1850 <- out
sim2_1955 <- out


plot(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births+
        outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections),
     type = "l", xlim = c(1955,2100), ylim = c(0, 4000),
     xlab = "Time", ylab = "New cases of chronic HBV carriage per timestep")
lines(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births+
        outpath$screened_incident_chronic_infections$screened_horizontal_chronic_infections), col = "pink")

outpath <- sim2_1850
outpath <- sim2_1955

plot(outpath$time,outpath$infectioncat_total$carriers,type = "l", ylim = c(0,400000), xlim = c(1955,2100))
lines(outpath$time,outpath$infectioncat_total$carriers,col = "pink")

# Plot total population size over timesteps
plot(outpath$time, outpath$pop_total$pop_total,
     xlab = "Year", ylab = "Population size", type = "l", ylim = c(0,7000000))
points(popsize_total$time, popsize_total$pop, col = "red")
lines(outpath$time, outpath$pop_total$pop_total, col = "pink")



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
o <- sim[[1]]$out
any(o[,1:(2*n_agecat*n_infectioncat)]<0)

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
     ylim = c(0,1), xlim = c(1960,2100))
lines(outpath$time,outpath$infectioncat_total$sus/outpath$pop_total$pop_total,col= "red")
lines(outpath$time,outpath$infectioncat_total$immune/outpath$pop_total$pop_total,col= "blue")

# Number of new cases of chronic HBV carriage at each timestep
plot(outpath$time,
     (outpath$incident_chronic_infections$horizontal_chronic_infections+
        outpath$incident_chronic_infections$chronic_births),
     type = "l", xlim = c(1955,2100), ylim = c(0, 3000),
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
load(here("output", "sims_output_scenario_vacc_130120.RData"))  # out_vacc
out_plot_vacc <- extract_outcomes(output_file = out_vacc, scenario_label = "vacc")
rm(out_vacc)
gc()
load(here("output", "sims_output_scenario_no_vacc_130120.RData")) # out_no_vacc
out_plot_no_vacc <- extract_outcomes(output_file = out_no_vacc, scenario_label = "no_vacc")
rm(out_no_vacc)
gc()

# Combine projection summaries from different scenarios
proj_prev_total <- rbind(out_plot_vacc$hbsag_prev_summary, out_plot_no_vacc$hbsag_prev_summary)
proj_number_inf_total <- rbind(out_plot_vacc$number_infected_summary, out_plot_no_vacc$number_infected_summary)
proj_inc_total <- rbind(out_plot_vacc$chronic_incidence_summary, out_plot_no_vacc$chronic_incidence_summary)
proj_inc_rate_total <- rbind(out_plot_vacc$chronic_incidence_rate_summary,
                             out_plot_no_vacc$chronic_incidence_rate_summary)
proj_deaths_total <- rbind(out_plot_vacc$hbv_deaths_summary, out_plot_no_vacc$hbv_deaths_summary)
proj_deaths_rate_total <- rbind(out_plot_vacc$hbv_deaths_rate_summary, out_plot_no_vacc$hbv_deaths_rate_summary)
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

# Number of people living with chronic HBV
ggplot(proj_number_inf_total) +
  geom_line(aes(x=time, y = median, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, group = scenario, fill = scenario), alpha = 0.1)+
  #  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Effect of infant vaccination on the number of people living with chronic HBV infection",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Number living with chronic HBV infection")+
  xlab("Year")+
#  ylim(0,25) +
  theme_classic()


# Chronic HBV incidence (absolute)
ggplot(proj_inc_total) +
  geom_line(aes(x=time, y = median/da, group = scenario, colour = scenario))+
  #  geom_line(aes(x=time, y = median/da, group = scenario, colour = scenario))+
    geom_ribbon(aes(x=time, ymin=lower/da, ymax=upper/da, group = scenario, fill = scenario), alpha = 0.1)+
  #  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Effect of infant vaccination on number of new chronic HBV infections",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Annual incident cases of chronic HBV infection")+
  xlab("Year")+
  #  ylim(0,25) +
  theme_classic()

# Annual chronic HBV incidence per 100,000 people
ggplot(proj_inc_rate_total) +
  geom_line(aes(x=time, y = (median*100000)/da, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=(lower*100000)/da, ymax=(upper*100000)/da, group = scenario, fill = scenario), alpha = 0.1)+
  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  geom_vline(aes(xintercept = 2015), col = "grey80", linetype = "dashed") +
  labs(title = "Effect of infant vaccination on incidence rate of chronic HBV infections",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Annual incidence of chronic HBV infections per 100,000 persons")+
  xlab("Year")+
  #  ylim(0,25) +
  theme_classic()

# HBV-related deaths (absolute)
ggplot(proj_deaths_total) +
  geom_line(aes(x=time, y = median/da, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=lower/da, ymax=upper/da, group = scenario, fill = scenario), alpha = 0.1)+
  # geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  labs(title = "Effect of infant vaccination on number of HBV-related deaths",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Annual number of HBV-related deaths")+
  xlab("Year")+
  ylim(0,6000) +
  theme_classic()

# Annual HBV-related deaths per 100000 persons
ggplot(proj_deaths_rate_total) +
  geom_line(aes(x=time, y = (median*100000)/da, group = scenario, colour = scenario))+
  geom_ribbon(aes(x=time, ymin=(lower*100000)/da, ymax=(upper*100000)/da, group = scenario, fill = scenario), alpha = 0.1)+
  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  geom_vline(aes(xintercept = 2015), col = "grey80", linetype = "dashed") +
  labs(title = "Effect of infant vaccination on HBV-related mortality rate",
       colour = "Modelled scenario", fill = "Modelled scenario",
       caption = "*Status quo scenario reflects historical vaccine coverage since 1990 and maintains 93% coverage after 2018") +
  scale_color_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_fill_manual(labels = c("No vaccine", "Status quo with vaccine*"), values = c("steelblue", "deeppink")) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Annual deaths per 100,000 persons")+
  xlab("Year")+
  ylim(0,100) +
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



# Analyse impact of infant vaccination compared to no historical intervention ----
load(here("output", "sims_output_scenario_vacc_130120.RData"))  # out_vacc
output_file <- out_vacc

# Extract absolute number of new cases of chronic HBV carriage per timestep for each iteration
proj_inc_vacc <- cbind(output_file[[1]]$time,
                  (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                     sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births")+
                     sapply(lapply(output_file,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections")))
colnames(proj_inc_vacc)[1] <- "time"

# Extract absolute incident HBV-related deaths per timestep
proj_deaths_vacc <- cbind(output_file[[1]]$time,
                     (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                        sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                        sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total")))
colnames(proj_deaths_vacc)[1] <- "time"

# Extract population size per timestep
proj_pop_vacc <- cbind(output_file[[1]]$time,
                       (sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total")))
colnames(proj_pop_vacc)[1] <- "time"
proj_pop_vacc <- as.data.frame(proj_pop_vacc)

# Extract number living with HBV per timestep
proj_prev_vacc <- cbind(output_file[[1]]$time,
                   (sapply(lapply(output_file,"[[", "infectioncat_total"), "[[", "carriers")))
colnames(proj_prev_vacc)[1] <- "time"
proj_prev_vacc <- as.data.frame(proj_prev_vacc)

# Extract number in need of treatment in 2019 (of 30+ year olds)
treat_comp_index <- c(which(grepl("^IRf.", names(init_pop))),
                      which(grepl("^IRm.", names(init_pop))),
                      which(grepl("^ENCHBf.", names(init_pop))),
                      which(grepl("^ENCHBm.", names(init_pop))),
                      which(grepl("^CCf.", names(init_pop))),
                      which(grepl("^CCm.", names(init_pop))),
                      which(grepl("^DCCf.", names(init_pop))),
                      which(grepl("^DCCm.", names(init_pop))))+1

full_output <- lapply(output_file,"[[", "full_output")
treat_comps <- lapply(full_output, "[", treat_comp_index)
treat_comps <- lapply(treat_comps, sum_pop_by_age, time = full_output[[1]]$time)
treat_comps_30plus <- lapply(treat_comps, "[", which(ages ==30):which(ages==100-da))
treatment_need_number_30plus <- cbind(output_file[[1]]$time, sapply(treat_comps_30plus,rowSums))
treatment_need_number_30plus  <- as.data.frame(treatment_need_number_30plus)
colnames(treatment_need_number_30plus)[1] <- "time"
# Summarise
quantile(treatment_need_number_30plus[treatment_need_number_30plus$time == 2019,-1], prob = c(0.025,0.5,0.975))

# Extract carriers aged 30+ years
carriers_30plus <- lapply(lapply(output_file,"[[", "carriers"), "[", which(ages==30):which(ages==100-da))
carriers_30plus <- cbind(output_file[[1]]$time, sapply(carriers_30plus,rowSums))
carriers_30plus <- as.data.frame(carriers_30plus)
colnames(carriers_30plus)[1] <- "time"

# Calculate treatment eligible %
quantile(treatment_need_number_30plus[treatment_need_number_30plus$time == 2019,-1]/
           carriers_30plus[carriers_30plus$time == 2019,-1], prob = c(0.025,0.5,0.975))


rm(out_vacc)
rm(output_file)
gc()

#load(here("output", "sims_output_scenario_no_vacc_130120.RData")) # out_no_vacc
output_file <- out_no_vacc
# Extract absolute number of new cases of chronic HBV carriage per timestep for each iteration
proj_inc_no_vacc <- cbind(output_file[[1]]$time,
                       (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                          sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births")+
                          sapply(lapply(output_file,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections")))
colnames(proj_inc_no_vacc)[1] <- "time"

# Extract absolute incident HBV-related deaths per timestep
proj_deaths_no_vacc <- cbind(output_file[[1]]$time,
                          (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                             sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                             sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total")))
colnames(proj_deaths_no_vacc)[1] <- "time"

# Extract population size per timestep
proj_pop_no_vacc <- cbind(output_file[[1]]$time,
                       (sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total")))
colnames(proj_pop_no_vacc)[1] <- "time"
proj_pop_no_vacc <- as.data.frame(proj_pop_no_vacc)

# Extract number living with HBV per timestep
proj_prev_no_vacc <- cbind(output_file[[1]]$time,
                        (sapply(lapply(output_file,"[[", "infectioncat_total"), "[[", "carriers")))
colnames(proj_prev_no_vacc)[1] <- "time"
proj_prev_no_vacc <- as.data.frame(proj_prev_no_vacc)

rm(out_no_vacc)
rm(output_file)
gc()

proj_inc_no_vacc <- as.data.frame(proj_inc_no_vacc)
proj_inc_vacc <- as.data.frame(proj_inc_vacc)
proj_deaths_no_vacc <- as.data.frame(proj_deaths_no_vacc)
proj_deaths_vacc <- as.data.frame(proj_deaths_vacc)

# Compare status quo scenario and counterfactual for each parameter set separately

# Cumulative number of new chronic infections averted since 1990:
# Cumulative number of new chronic infections between 1990 and 2020 without vaccine:

cum_inc_no_vacc <- apply(proj_inc_no_vacc[which(proj_inc_no_vacc$time==1990.5):which(proj_inc_no_vacc$time==2020),-1],2,sum)
cum_inc_vacc <- apply(proj_inc_vacc[which(proj_inc_vacc$time==1990.5):which(proj_inc_vacc$time==2020),-1],2,sum)
cum_inc_averted <- cum_inc_no_vacc-cum_inc_vacc
# Summarise the number of new chronic infections averted
quantile(cum_inc_averted, prob = c(0.025,0.5,0.975))

# Cumulative number of HBV-related deaths averted since 1990:

cum_deaths_no_vacc <- apply(proj_deaths_no_vacc[which(proj_deaths_no_vacc$time==1990.5):which(proj_deaths_no_vacc$time==2020),-1],2,sum)
cum_deaths_vacc <- apply(proj_deaths_vacc[which(proj_deaths_vacc$time==1990.5):which(proj_deaths_vacc$time==2020),-1],2,sum)
cum_deaths_averted <- cum_deaths_no_vacc-cum_deaths_vacc
# Summarise the number of HBV-related deaths averted
quantile(cum_deaths_averted, prob = c(0.025,0.5,0.975))

# How many incident cases of chronic infection and HBV-related deaths in 2019? Status quo scenario
# Incidence was calculated as cases from t-1 to t
inc_in_2019 <- apply(proj_inc_vacc[which(proj_inc_vacc$time==2019.5):which(proj_inc_vacc$time==2020),-1],2,sum)
quantile(inc_in_2019, prob = c(0.025,0.5,0.975))
# Deaths in 2019
deaths_in_2019 <- apply(proj_deaths_vacc[which(proj_deaths_vacc$time==2019.5):which(proj_deaths_vacc$time==2020),-1],2,sum)
quantile(deaths_in_2019, prob = c(0.025,0.5,0.975))
# Number living with chronic HBV in 2019
n_hbsag_in_2019 <- apply(proj_prev_vacc[which(proj_prev_vacc$time==2019):which(proj_prev_vacc$time==2019.5),-1],2,mean)
quantile(n_hbsag_in_2019, prob = c(0.025,0.5,0.975))
# % prevalence in 2019
prev_hbsag_in_2019 <- apply(proj_prev_vacc[which(proj_prev_vacc$time==2019):which(proj_prev_vacc$time==2019.5),-1],2,mean)/
  apply(proj_pop_vacc[which(proj_pop_vacc$time==2019):which(proj_pop_vacc$time==2019.5),-1],2,mean)
quantile(prev_hbsag_in_2019, prob = c(0.025,0.5,0.975))

# Reduction in new cases of chronic infection per 100,000 population in 2030 compared to 2015:
# Cases in 2030: calculate cases from 2030-2030.5+cases from 2030.5-2031, divide by average of population in 2030 and 2030.5
inc_per_100000_2030 <- apply(proj_inc_vacc[which(proj_inc_vacc$time==2030.5):which(proj_inc_vacc$time==2031),-1],2,sum)*100000/
  apply(proj_pop_vacc[which(proj_pop_vacc$time==2030):which(proj_pop_vacc$time==2030.5),-1],2,mean)
inc_per_100000_2015 <- apply(proj_inc_vacc[which(proj_inc_vacc$time==2015.5):which(proj_inc_vacc$time==2016),-1],2,sum)*100000/
  apply(proj_pop_vacc[which(proj_pop_vacc$time==2015):which(proj_pop_vacc$time==2015.5),-1],2,mean)
# reduction percentage:
quantile((inc_per_100000_2015-inc_per_100000_2030)/inc_per_100000_2015, prob =c(0.025,0.5,0.975))

# Reduction in HBV-related deaths per 100,000 population in 2030 compared to 2015:
deaths_per_100000_2030 <- apply(proj_deaths_vacc[which(proj_deaths_vacc$time==2030.5):which(proj_deaths_vacc$time==2031),-1],2,sum)*100000/
  apply(proj_pop_vacc[which(proj_pop_vacc$time==2030):which(proj_pop_vacc$time==2030.5),-1],2,mean)
deaths_per_100000_2015 <- apply(proj_deaths_vacc[which(proj_deaths_vacc$time==2015.5):which(proj_deaths_vacc$time==2016),-1],2,sum)*100000/
  apply(proj_pop_vacc[which(proj_pop_vacc$time==2015):which(proj_pop_vacc$time==2015.5),-1],2,mean)
# reduction percentage:
quantile((deaths_per_100000_2015-deaths_per_100000_2030)/deaths_per_100000_2015, prob =c(0.025,0.5,0.975))

# Reduction in HBsAg prevalence % in 2030 compared to 2015:
prev_2030 <- apply(proj_prev_vacc[which(proj_prev_vacc$time==2030):which(proj_prev_vacc$time==2030.5),-1],2,mean)/
  apply(proj_pop_vacc[which(proj_pop_vacc$time==2030):which(proj_pop_vacc$time==2030.5),-1],2,mean)
prev_2015 <- apply(proj_prev_vacc[which(proj_prev_vacc$time==2015):which(proj_prev_vacc$time==2015.5),-1],2,mean)/
  apply(proj_pop_vacc[which(proj_pop_vacc$time==2015):which(proj_pop_vacc$time==2015.5),-1],2,mean)
quantile((prev_2015-prev_2030)/prev_2015, prob = c(0.025,0.5,0.975))

# Time to elimination:
calculate_reduction_compared_to_2015 <- function(inc_data, pop_data, year) {
  inc_per_100000_pred <- apply(inc_data[which(inc_data$time==year+0.5):which(inc_data$time==year+1),-1],2,sum)*100000/
    apply(pop_data[which(pop_data$time==year):which(pop_data$time==year+0.5),-1],2,mean)
  inc_per_100000_baseline <- apply(inc_data[which(inc_data$time==2015.5):which(inc_data$time==2016),-1],2,sum)*100000/
    apply(pop_data[which(pop_data$time==2015):which(pop_data$time==2015.5),-1],2,mean)

  # reduction percentage:
  red <- quantile((inc_per_100000_baseline-inc_per_100000_pred)/inc_per_100000_baseline, prob =c(0.025,0.5,0.975))

  return(red)
}

calculate_year_of_elimination <- function(inc_data, pop_data, years) {

  inc_per_100000_pred <- list()

  inc_per_100000_baseline <- apply(inc_data[which(inc_data$time==2015.5):which(inc_data$time==2016),-1],2,sum)*100000/
    apply(pop_data[which(pop_data$time==2015):which(pop_data$time==2015.5),-1],2,mean)

  for (i in years) {
  inc_per_100000_pred[[i]] <- apply(inc_data[which(inc_data$time==i+0.5):which(inc_data$time==i+1),-1],2,sum)*100000/
    apply(pop_data[which(pop_data$time==i):which(pop_data$time==i+0.5),-1],2,mean)
  }

  #names(inc_per_100000_pred) <- seq_along(inc_per_100000_pred)
  inc_per_100000_pred[sapply(inc_per_100000_pred, is.null)] <- NULL

  red <- list()
  for (i in 1:length(years)) {
    red[[i]] <- (inc_per_100000_baseline-inc_per_100000_pred[[i]])/inc_per_100000_baseline
  }

  red <- matrix(unlist(red), ncol = 119, byrow = TRUE)
  rownames(red) <- years

  # reduction percentage:
  #red <- quantile((inc_per_100000_baseline-inc_per_100000_pred)/inc_per_100000_baseline, prob =c(0.025,0.5,0.975))

  return(red)
}
test <- calculate_year_of_elimination(proj_deaths_vacc, proj_pop_vacc, c(2030,2050, 2060, 2090, 2098))
apply(test,2,function(x){min(which(x>0.65))})

# Very few parameter sets achieve elimination before 2100! Would have to run model for longer to calculate median
# time to elimination!



# Test ----
#Calculate vaccine coverage
# Extract for 1 year olds:
ever_infected <- as.data.frame(sapply(lapply(out, "[[", "ever_infected"), "[", 3))
carriers <- as.data.frame(sapply(lapply(out, "[[", "carriers"), "[", 3))
pop <- as.data.frame(sapply(lapply(out, "[[", "pop"), "[", 3))

ever_infected_1990 <- ever_infected[which(out$`9035`$time==1990):nrow(ever_infected),]
carriers_1990 <- carriers[which(out$`9035`$time==1990):nrow(carriers),]
pop_1990 <- pop[which(out$`9035`$time==1990):nrow(pop),]

# Calculate proportion recovered:
prop_rec_1989 <- unlist(ever_infected[which(out$`9035`$time==1990),]/pop[which(out$`9035`$time==1990),])
prop_in_r <- ((ever_infected_1990-carriers_1990)/pop_1990)
prop_vacc <- sweep(prop_in_r, 2, prop_rec_1989, "-")
prop_vacc$time <- seq(1990, 2100.5, by = 0.5)
prop_vacc <- as.data.frame(prop_vacc)
prop_vacc <- gather(prop_vacc,key = "sim", value = "prop", -time)

ggplot(data = prop_vacc) +
  geom_line(aes(x = time, y = prop, group = sim, col = "Simulated effective coverage")) +
  geom_point(data = vaccine_coverage, aes(x = year, y = coverage, col = "Input coverage")) +
  scale_colour_manual(values = c("Simulated effective coverage" = "grey", "Input coverage" = "red")) +
  labs(colour = "") +
  ylab("Proportion immunised among 1-year olds") +
  theme_classic()



# Prepare output data to plot
load(here("output", "sims_output_scenario_vacc_new_method_110220.RData"))
out_vacc <- out_vacc_new_method
out_plot_vacc <- extract_outcomes(output_file = out_vacc, scenario_label = "vacc")
rm(out_vacc)
gc()

# New method:
# proj_inc_summary, proj_inc_rate_summary, proj_deaths_summary, proj_deaths_rate
output_file <- out_vacc

proj_inc <- cbind(output_file[[1]]$time,
                  (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                     sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births")))
colnames(proj_inc)[1] <- "time"

proj_inc_summary <- data.frame(time = output_file[[1]]$time)
proj_inc_summary$median <- apply(proj_inc[,-1],1,median)
proj_inc_summary$lower <- apply(proj_inc[,-1],1,quantile, prob = 0.025)
proj_inc_summary$upper <- apply(proj_inc[,-1],1,quantile, prob = 0.975)
proj_inc_long <- gather(as.data.frame(proj_inc), key = "iteration", value = "chronic_cases", -time)


# Absolute incident HBV-related deaths per timestep
proj_deaths <- cbind(output_file[[1]]$time,
                     (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                        sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                        sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total")))
colnames(proj_deaths)[1] <- "time"
proj_deaths_summary <- data.frame(time = output_file[[1]]$time)
proj_deaths_summary$median <- apply(proj_deaths[,-1],1,median)
proj_deaths_summary$lower <- apply(proj_deaths[,-1],1,quantile, prob = 0.025)
proj_deaths_summary$upper <- apply(proj_deaths[,-1],1,quantile, prob = 0.975)
proj_deaths_long <- gather(as.data.frame(proj_deaths), key = "iteration", value = "deaths", -time)

ggplot(proj_inc_summary) +
  geom_line(aes(x=time, y = median), col = "red")+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper), alpha = 0, col = "blue")+
  #  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Annual incident cases of chronic HBV infection")+
  xlab("Year")+
  ylim(0,5000) +
  theme_classic()

ggplot() +
  geom_line(data= proj_inc_long, aes(x=time, y = chronic_cases, group = iteration), col = "grey")+
  geom_line(data= proj_inc_summary, aes(x=time, y = median), col = "red", size = 2) +
  geom_ribbon(data= proj_inc_summary, aes(x=time, ymin=lower, ymax=upper), alpha = 0, col = "blue")+
  #  geom_vline(aes(xintercept = 2030), col = "grey80", linetype = "dashed") +
  #  xlim(1960,2100) +
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Annual incident cases of chronic HBV infection")+
  xlab("Year")+
    ylim(0,5000) +
  theme_classic()

# Plots for mock parms (to test) ----
# Absolute number of new cases of chronic HBV carriage per timestep
output_file <- out
()
proj_inc <- cbind(output_file[[1]]$time,
                  (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                     sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births")))
colnames(proj_inc) <- c("time", "median", "lower", "upper")

# Incidence rate of chronic HBV carriage per timestep
proj_inc_rate <- cbind(output_file[[1]]$time,
                       (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                          sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births"))/
                         sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total"))
colnames(proj_inc_rate) <- c("time", "median", "lower", "upper")

# Absolute incident HBV-related deaths per timestep
proj_deaths <- cbind(output_file[[1]]$time,
                     (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                        sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                        sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total")))
colnames(proj_deaths) <- c("time", "median", "lower", "upper")

# Incidence rate of  HBV-related deaths per timestep
proj_deaths_rate <- cbind(output_file[[1]]$time,
                          (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                             sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                             sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total"))/
                            sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total"))
colnames(proj_deaths_rate) <- c("time", "median", "lower", "upper")

# Chronic HBV incidence (absolute number)
ggplot(as.data.frame(proj_inc)) +
  geom_line(aes(x=time, y = median))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper), alpha = 0.1)+
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Incident cases of chronic HBV infection per 6 months")+
  xlab("Year")+
  ylim(0,5000) +
  theme_classic()

# Chronic HBV incidence rate
ggplot(as.data.frame(proj_inc_rate)) +
  geom_line(aes(x=time, y = median*100000))+
  geom_ribbon(aes(x=time, ymin=lower*100000, ymax=upper*100000), alpha = 0.1)+
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("Chronic infection incidence per 100,000 per 6 months")+
  xlab("Year")+
  ylim(0,500) +
  theme_classic()

# HBV deaths (absolute number)
ggplot(as.data.frame(proj_deaths)) +
  geom_line(aes(x=time, y = median))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper), alpha = 0.1)+
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("HBV deaths per 6 months")+
  xlab("Year")+
  ylim(0,500) +
  theme_classic()

# HBV death rate
ggplot(as.data.frame(proj_deaths_rate)) +
  geom_line(aes(x=time, y = median*100000))+
  geom_ribbon(aes(x=time, ymin=lower*100000, ymax=upper*100000), alpha = 0.1)+
  scale_x_continuous(breaks=seq(1960, 2100, by = 10), limits = c(1960,2100)) +
  ylab("HBV deaths per 100,000 per 6 months")+
  xlab("Year")+
  ylim(0,30) +
  theme_classic()
