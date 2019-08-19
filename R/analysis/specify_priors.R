# First prior specification
# Some of these are quite informative and not always in line with the african data

library(truncnorm)

# 1) Draw parameter sets randomly from prior distribution
n_sims <- 10000  # number of simulations/parameter sets
# First sample parameters where prior distributions depend on each other
b1 <- runif(n_sims, 0.03, 0.7)
b2 <- runif(n_sims, 0, b1)
b3 <- runif(n_sims, 0, b1)
mtct_prob_s <- rbeta(n_sims, 1.5,13.5)
mtct_prob_e <- runif(n_sims, mtct_prob_s, 0.9)
hccr_it <- rtruncnorm(n_sims, a=1, b=Inf, mean=6, sd=3) # normal truncated at 1
hccr_cc <- runif(n_sims, hccr_it, 100)
hccr_enchb <- runif(n_sims, hccr_it, hccr_cc)
hccr_ir <- runif(n_sims, hccr_enchb, hccr_cc)
# Combine in dataframe and sample remaining parameters
params_mat <- data.frame(b1 = b1,
                         b2 = b2,
                         b3 = b3,
                         mtct_prob_s = mtct_prob_s,
                         mtct_prob_e = mtct_prob_e,
                         alpha = runif(n_sims, 1.5,10),
                         p_chronic_in_mtct = rbeta(n_sims, 10.49,1.3),
                         p_chronic_function_r = rnorm(n_sims,0.65,0.1),
                         p_chronic_function_s = rnorm(n_sims,0.46,0.1),
                         pr_it_ir = rgamma(n_sims,3.63,26.27),
                         pr_ir_ic = runif(n_sims, 0,1),
                         eag_prog_function_rate = runif(n_sims,0,0.01),
                         pr_ir_enchb = rgamma(n_sims, 1.49, 97.58),
                         pr_ir_cc_female = runif(n_sims, 0.005, 0.05),
                         pr_ir_cc_age_threshold = sample(0:15,n_sims,replace=TRUE),
                         pr_ic_enchb = rgamma(n_sims, 3.12, 141.30),
                         sag_loss_slope = rnorm(n_sims, 0.0004106, 0.00005),
                         pr_enchb_cc_female = rgamma(n_sims, 2.3, 123.8),
                         cirrhosis_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                         pr_cc_dcc = rgamma(n_sims,17.94,423.61),
                         cancer_prog_coefficient_female = runif(n_sims, 0.0001, 0.0003),
                         cancer_age_threshold = sample(0:15,n_sims,replace=TRUE),
                         cancer_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                         hccr_it = hccr_it,
                         hccr_ir = hccr_ir,
                         hccr_enchb = hccr_enchb,
                         hccr_cc = hccr_cc,
                         hccr_dcc = rgamma(n_sims, 8.09, 101.33),
                         mu_cc = rgamma(n_sims, 4.25, 124.91),
                         mu_dcc = rgamma(n_sims, 1.49, 0.98),
                         mu_hcc = rgamma(n_sims, 1.49, 0.98),
                         vacc_eff = rbeta(n_sims, 7.07, 0.37))


plot_prior <- function(parm) {
  plot(density(prior[,parm]), main = parm, lwd=3, lty=2, col="blue")
}

prior <- params_mat

par(mfrow=c(1,1))
plot_prior("b1")
plot_prior("b2")
plot_prior("b3")
plot_prior("alpha")
plot_prior("mtct_prob_s")
plot_prior("mtct_prob_e")
plot_prior("p_chronic_in_mtct")
plot_prior("p_chronic_function_r")
plot_prior("p_chronic_function_s")
# up until here they are good
plot_prior("pr_it_ir")
plot_prior("pr_ir_ic")
plot_prior("eag_prog_function_rate")
# these are also OK
plot_prior("pr_ir_enchb")  # this one is quite specific
plot_prior("pr_ir_cc_female")
plot_prior("pr_ir_cc_age_threshold")
plot_prior("pr_ic_enchb")  # need to change this one to reflect African data more
plot_prior("sag_loss_slope")  # OK based on AFrican data
plot_prior("pr_enchb_cc_female")  # quite strong
plot_prior("cirrhosis_male_cofactor")
plot_prior("pr_cc_dcc")   # OK cause wouldn't expect this one to be different
# other 2 are ok
plot_prior("cancer_prog_coefficient_female")
plot_prior("cancer_age_threshold")
plot_prior("cancer_male_cofactor")
# ok
plot_prior("hccr_it")
plot_prior("hccr_ir")
plot_prior("hccr_enchb")
plot_prior("hccr_cc")
# OK
plot_prior("hccr_dcc")  # can think about this one
plot_prior("mu_cc")  # ok
plot_prior("mu_dcc")  # adapt this one from African data
plot_prior("mu_hcc")  # adapt this one from African data
plot_prior("vacc_eff")

# Update some priors

# Changed priors to reflect African data better
# pr_ic_enchb
# before:
plot(density(rgamma(n_sims, 3.12, 141.30)), ylim = c(0,35), xlim = c(0,0.15))
abline(v = 0.045)
abline(v = 0.015)
# BETTER: center on 0.03, range 0.01 to 0.1
plot(density(rgamma(n_sims, 3.49, 83.05)), ylim = c(0,35), xlim = c(0,0.15))
abline(v = 0.045)
abline(v = 0.03)
abline(v = 0.015)

# mu_dcc: center on 1 because of African data, range unchanged
# before:
plot(density(rgamma(n_sims, 1.49, 0.98)))
abline(v = 1)
abline(v = 0.5)
# better:
plot(density(rgamma(n_sims, 2.18, 1.18)))
abline(v = 1)
abline(v = 0.5)
# chose mu_hcc same as mu_dcc


# To check

# pr_ir_enchb
# this one I took from Shevanthi's because I wanted to reflect it is rare
# and we have no data
plot(density(rgamma(n_sims, 1.49, 97.58)), ylim=c(0,45), xlim =c(0,0.2))
# double the range (0-0.05) now (0-0.1), keep center at 0.005
# BETTER
plot(density(rgamma(n_sims, 1.22, 44.20)), ylim=c(0,45), xlim =c(0,0.2))

# pr_enchb_cc_female could go up higher to 0.2 but keep the center
# before
plot(density(rgamma(n_sims, 2.3, 123.8)), ylim = c(0,40), xlim = c(0,0.3))
# BETTER
plot(density(rgamma(n_sims, 1.23, 22.33)), ylim = c(0,40), xlim = c(0,0.3))

# hccr_dcc
plot(density(rgamma(n_sims, 8.09, 101.33)), ylim = c(0,15), xlim = c(0,0.5))
# this one is OK to be specific because the source is from several different locations
# and the data should be able to give information on this
# but just to ensure consistency with hccr_cc, double the range (0.004 to 0.24)
# keep center
# BETTER
plot(density(rgamma(n_sims, 3.08, 29.76)), ylim = c(0,15), xlim = c(0,0.5))


# Plot all priors with the new ones
params_mat <- data.frame(b1 = b1,
                         b2 = b2,
                         b3 = b3,
                         mtct_prob_s = mtct_prob_s,
                         mtct_prob_e = mtct_prob_e,
                         alpha = runif(n_sims, 1.5,10),
                         p_chronic_in_mtct = rbeta(n_sims, 10.49,1.3),
                         p_chronic_function_r = rnorm(n_sims,0.65,0.1),
                         p_chronic_function_s = rnorm(n_sims,0.46,0.1),
                         pr_it_ir = rgamma(n_sims,3.63,26.27),
                         pr_ir_ic = runif(n_sims, 0,1),
                         eag_prog_function_rate = runif(n_sims,0,0.01),
                         pr_ir_enchb = rgamma(n_sims, 1.22, 44.20),
                         pr_ir_cc_female = runif(n_sims, 0.005, 0.05),
                         pr_ir_cc_age_threshold = sample(0:15,n_sims,replace=TRUE),
                         pr_ic_enchb = rgamma(n_sims, 3.49, 83.05),
                         sag_loss_slope = rnorm(n_sims, 0.0004106, 0.00005),
                         pr_enchb_cc_female = rgamma(n_sims, 1.23, 22.33),
                         cirrhosis_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                         pr_cc_dcc = rgamma(n_sims,17.94,423.61),
                         cancer_prog_coefficient_female = runif(n_sims, 0.0001, 0.0003),
                         cancer_age_threshold = sample(0:15,n_sims,replace=TRUE),
                         cancer_male_cofactor = rtruncnorm(n_sims, a = 1, mean = 3.5, sd = 4),
                         hccr_it = hccr_it,
                         hccr_ir = hccr_ir,
                         hccr_enchb = hccr_enchb,
                         hccr_cc = hccr_cc,
                         hccr_dcc = rgamma(n_sims, 3.08, 29.76),
                         mu_cc = rgamma(n_sims, 4.25, 124.91),
                         mu_dcc = rgamma(n_sims, 2.18, 1.18),
                         mu_hcc = rgamma(n_sims, 2.18, 1.18),
                         vacc_eff = rbeta(n_sims, 7.07, 0.37))

prior <- params_mat

plot_prior("b1")
plot_prior("b2")
plot_prior("b3")
plot_prior("alpha")
plot_prior("mtct_prob_s")
plot_prior("mtct_prob_e")
plot_prior("p_chronic_in_mtct")
plot_prior("p_chronic_function_r")
plot_prior("p_chronic_function_s")
plot_prior("pr_it_ir")
plot_prior("pr_ir_ic")
plot_prior("eag_prog_function_rate")
plot_prior("pr_ir_enchb")
plot_prior("pr_ir_cc_female")
plot_prior("pr_ir_cc_age_threshold")
plot_prior("pr_ic_enchb")
plot_prior("sag_loss_slope")
plot_prior("pr_enchb_cc_female")
plot_prior("cirrhosis_male_cofactor")
plot_prior("pr_cc_dcc")
plot_prior("cancer_prog_coefficient_female")
plot_prior("cancer_age_threshold")
plot_prior("cancer_male_cofactor")
plot_prior("hccr_it")
plot_prior("hccr_ir")
plot_prior("hccr_enchb")
plot_prior("hccr_cc")
plot_prior("hccr_dcc")
plot_prior("mu_cc")
plot_prior("mu_dcc")
plot_prior("mu_hcc")
plot_prior("vacc_eff")
