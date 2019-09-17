# Future predictions with calibrated model (LSR stage - weighted Euclidian distance, 50 of 100000)

# Need to load an out_mat with all the accepted simulations
library(here)
library(tidyr)
library(dplyr)
load(file = here("output", "fits", "best_fits_50_of_100000_wed_domain_weights_210819.Rdata"))

source(here("R", "imperial_model_main.R"))

# Save parameter sets
parmsets <- as.data.frame(do.call(rbind, lapply(out_mat_wed_domain_50, "[[", "parameter_set")))

sim <- apply(parmsets,1,
                 function(x)
                   run_model(sim_duration = 250,
                             default_parameter_list = parameter_list,
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
                                    vacc_eff = as.list(x)$vacc_eff
                               )))

#save(sim, file = here("output", "fits","lsr_plots", "sims_with_accepted_parms.Rdata"))
load(file = here("output", "fits","lsr_plots", "sims_with_accepted_parms.Rdata"))


#sim2 <- do.call(rbind.data.frame, sim)
out <- lapply(sim, code_model_output)

# Calculate HBsAg prevalence over time
prev <- list()
for(i in 1:50) {
  prev[[i]] <- lapply(lapply(out, "[[", "infectioncat_total"), "[[", "carriers")[[i]]/
    lapply(lapply(out, "[[", "pop_total"), "[[", "pop_total")[[i]]
}
prev <- as.data.frame(do.call(cbind, prev))
prev$time <- out[[1]]$time

prev_long <- gather(prev, key = "sim", "value" = prev, -time)

# Calculate median and 5th and 95th percentile at each timepoint
plot(x = out[[1]]$time, y = apply(prev,1,median), type = "l", xlim = c(1980,2050), ylim = c(0,0.35))
lines(x = out[[1]]$time, y = apply(prev,1,quantile, probs = 0.05))
lines(x = out[[1]]$time, y = apply(prev,1,quantile, probs = 0.95))

library(ggplot2)

ggplot() +
  geom_ribbon(data = prev, aes(x = time, ymin = apply(prev,1,quantile, probs = 0.05),
                               ymax = apply(prev,1,quantile, probs = 0.95)), fill = "steelblue",
              alpha = 0.5)+
  geom_line(data = prev, aes(x = time, y = apply(prev,1,median)), col = "black")+
  geom_line(data = prev, aes(x = time, y = apply(prev,1,quantile, probs = 0.05)), col = "black",
            linetype = "dashed")+
  geom_line(data = prev, aes(x = time, y = apply(prev,1,quantile, probs = 0.95)), col = "black",
            linetype = "dashed")+
  theme_bw() +
  xlim(c(1980,2100)) +
  ylim(c(0,0.3)) +
  xlab("Time") + ylab("HBsAg prevalence (proportion)")



# Plots

outpath <- out[[1]]

## INFECTION DYNAMICS

# Total number in each infection compartment per timestep
plot(outpath$time,outpath$infectioncat_total$carriers,type = "l", ylim = c(0,4000000))
lines(outpath$time,outpath$infectioncat_total$sus,col= "red")
lines(outpath$time,outpath$infectioncat_total$immune,col= "blue")

# Proportion in each infection compartment per timestep
plot(outpath$time,outpath$infectioncat_total$carriers/outpath$pop_total$pop_total,type = "l",
     ylim = c(0,0.7), xlim = c(1880,2100))
lines(outpath$time,outpath$infectioncat_total$sus/outpath$pop_total$pop_total,col= "red")
lines(outpath$time,outpath$infectioncat_total$immune/outpath$pop_total$pop_total,col= "blue")

plot(outpath$time,outpath$infectioncat_total$carriers/outpath$pop_total$pop_total,type = "l",
     xlim = c(1980,2100))

# Carrier prevalence in 1980
outpath$infectioncat_total$carriers[which(outpath$infectioncat_total$time == 1980)]/
  outpath$pop_total$pop_total[which(outpath$pop_total$time == 1980)]
# Carrier prevalence in 1990
outpath$infectioncat_total$carriers[which(outpath$infectioncat_total$time == 1990)]/
  outpath$pop_total$pop_total[which(outpath$pop_total$time == 1990)]
# Carrier prevalence in 2015
outpath$infectioncat_total$carriers[which(outpath$infectioncat_total$time == 2015)]/
  outpath$pop_total$pop_total[which(outpath$pop_total$time == 2015)]
# Carrier prevalence in 2020
outpath$infectioncat_total$carriers[which(outpath$infectioncat_total$time == 2020)]/
  outpath$pop_total$pop_total[which(outpath$pop_total$time == 2020)]
# Reduction:
(outpath$infectioncat_total$carriers[which(outpath$infectioncat_total$time == 2015)]/
    outpath$pop_total$pop_total[which(outpath$pop_total$time == 2015)])/(outpath$infectioncat_total$carriers[which(outpath$infectioncat_total$time == 1990)]/
                                                                           outpath$pop_total$pop_total[which(outpath$pop_total$time == 1990)])

# Carrier prevalence over time in different age groups
plot(x = outpath$time, y = as.numeric(unlist(outpath$carriers[,2]/outpath$pop[,2]))) # age 0.5
abline(v = 1991)
plot(x = outpath$time, y = as.numeric(unlist(outpath$carriers[,3]/outpath$pop[,3]))) # age 1
abline(v = 1991)
plot(x = outpath$time, y = as.numeric(unlist(outpath$carriers[,4]/outpath$pop[,4]))) # age 1.5
abline(v = 1991)
plot(x = outpath$time, y = as.numeric(unlist(outpath$carriers[,5]/outpath$pop[,5]))) # age 2
abline(v = 1991)
plot(x = outpath$time, y = as.numeric(unlist(outpath$carriers[,22]/outpath$pop[,22])))
abline(v = 1991)

# Carrier prevalence by age in 1980
plot(ages, outpath$carriers[which(outpath$time == 1980),]/
       outpath$pop[which(outpath$time == 1980),], type = "l", ylim = c(0,0.3))
#points(gambia_prevdata$age, gambia_prevdata$edmunds_prev, col = "red")

# anti-HBc prevalence by age in 1980
plot(ages, outpath$ever_infected[which(outpath$time == 1980),]/
       outpath$pop[which(outpath$time == 1980),], type = "l", ylim = c(0,1))

# HBeAg prevalence in chronic carriers by age in 1980
plot(ages, outpath$eag_positive[which(outpath$time == 1980),]/
       outpath$carriers[which(outpath$time == 1980),], type = "l", ylim = c(0,1))
points(input_hbeag_dataset$age, input_hbeag_dataset$data_value)

# Carrier prevalence by age in 2015
plot(ages, outpath$carriers[which(outpath$time == 2015),]/
       outpath$pop[which(outpath$time == 2015),], type = "l", ylim = c(0,0.3))
plot(ages, outpath$carriers[which(outpath$time == 2020),]/
       outpath$pop[which(outpath$time == 2020),], type = "l", ylim = c(0,0.3))

# HBeAg prevalence in chronic carriers by age in 2015
plot(ages, outpath$eag_positive[which(outpath$time == 2015),]/
       outpath$carriers[which(outpath$time == 2015),], type = "l", ylim = c(0,1))

