# Natural history model uncertainty analysis
library(here)
library(tidyr)
library(dplyr)
library(corrplot)
library(ggplot2)

# Comparison of prior and posterios (initially run in run_calibration_local.R) ----

# Load parameter sets chosen based on kmeans clustering
load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans
# Load priors
load(here("calibration", "input", "lhs_samples_1000000.Rdata"))

# Table of median, 2.5th and 97.5th quantile of priors and posteriors
posterior_summary <- gather(params_mat_accepted_kmeans, key = "parameter", value = "sim") %>%
  group_by(parameter) %>%
  summarise(post_median = round(median(sim),4),
            post_cri_lower = round(quantile(sim, prob = 0.025),4),
            post_cri_upper = round(quantile(sim, prob = 0.975),4))

prior_summary <- gather(params_mat, key = "parameter", value = "sim") %>%
  group_by(parameter) %>%
  summarise(prior_median = round(median(sim),4),
            prior_cri_lower = round(quantile(sim, prob = 0.025),4),
            prior_cri_upper = round(quantile(sim, prob = 0.975),4))

prior_posterior_summary <- left_join(prior_summary, posterior_summary, by ="parameter")
#write.csv(prior_posterior_summary, file=here("calibration", "output", "prior_posterior_summary_table_161120.csv"), row.names = FALSE)

# Histograms of prior and posterior distribution
plot_prior_posterior <- function(parm) {
  plot(density(posterior[,parm]), xlim = c(min(min(prior[,parm]),min((posterior[,parm]))), max(max(prior[,parm]),max((posterior[,parm])))),
       ylim = c(0, max(max(density(prior[,parm])$y),max((density(posterior[,parm])$y)))), main= parm,
       lwd=3, col="red")
  if (parm %in% c("pr_ir_cc_age_threshold", "cancer_age_threshold")) {
    lines(density(prior[,parm], bw = 1), lwd=3, lty=2, col="blue")
  } else {
    lines(density(prior[,parm]), lwd=3, lty=2, col="blue")
  }
  #  legend("bottomleft", legend=c("prior density","posterior density"),
  #         col=c("blue","red"), lty=c(3,1), lwd=c(3,3), cex = 1)
}

prior <- params_mat
posterior <- params_mat_accepted_kmeans

plot(density(prior[,"pr_ir_cc_age_threshold"], bw = 1))
hist(prior[,"pr_ir_cc_age_threshold"], breaks=seq(min(prior[,"pr_ir_cc_age_threshold"])-0.5,
                                                  max(prior[,"pr_ir_cc_age_threshold"])+0.5, by=1))


#pdf(file = here("calibration", "output", "prior_posterior_density_plots_240720.pdf"), title="Prior and posterior density")
par(mfrow=c(3,3))
plot(x=0,y=0, col = "white", xlab = "", ylab = "")
legend("center", legend=c("prior density","posterior density"),
       col=c("blue","red"), lty=c(3,1), lwd=c(3,3), cex = 1)
plot_prior_posterior("b1")
plot_prior_posterior("b2")
plot_prior_posterior("b3")
plot_prior_posterior("alpha")
plot_prior_posterior("mtct_prob_s")
plot_prior_posterior("mtct_prob_e")
plot_prior_posterior("p_chronic_in_mtct")
plot_prior_posterior("p_chronic_function_r")
plot_prior_posterior("p_chronic_function_s")
plot_prior_posterior("pr_it_ir")
plot_prior_posterior("pr_ir_ic")
plot_prior_posterior("eag_prog_function_rate")
plot_prior_posterior("pr_ir_enchb")
plot_prior_posterior("pr_ir_cc_female")
plot_prior_posterior("pr_ir_cc_age_threshold")
plot_prior_posterior("pr_ic_enchb")
plot_prior_posterior("sag_loss_slope")
plot_prior_posterior("pr_enchb_cc_female")
plot_prior_posterior("cirrhosis_male_cofactor")
plot_prior_posterior("pr_cc_dcc")
plot_prior_posterior("cancer_prog_coefficient_female")
plot_prior_posterior("cancer_age_threshold")
plot_prior_posterior("cancer_male_cofactor")
plot_prior_posterior("hccr_it")
plot_prior_posterior("hccr_ir")
plot_prior_posterior("hccr_enchb")
plot_prior_posterior("hccr_cc")
plot_prior_posterior("hccr_dcc")
plot_prior_posterior("mu_cc")
plot_prior_posterior("mu_dcc")
plot_prior_posterior("mu_hcc")
plot_prior_posterior("vacc_eff")
par(mfrow=c(1,1))
#dev.off()



# Analysis of the prior and posterior interquartile range ----
iqr_prior_post <- rbind(t(apply(prior, 2, quantile, prob = c(0.25, 0.5, 0.75))),
                         t(apply(posterior, 2, quantile, prob = c(0.25, 0.5, 0.75))))
colnames(iqr_prior_post) <- c("lower", "median", "upper")
iqr_prior_post <- as.data.frame(iqr_prior_post)
iqr_prior_post$type <- c(rep("prior", 32), rep("post", 32))
iqr_prior_post$parm <- c(colnames(prior), colnames(posterior))
iqr_prior_post$iqr_width <- iqr_prior_post$upper-iqr_prior_post$lower

# Calculate reduction in IQR width between prior and posterior for each parameter
iqr_width_reduction <- select(iqr_prior_post, parm, type, iqr_width) %>%
  spread(key = "type", value = "iqr_width")
iqr_width_reduction$abs_red <- round(iqr_width_reduction$prior-iqr_width_reduction$post,5)
iqr_width_reduction$rel_red <- round(iqr_width_reduction$abs_red/iqr_width_reduction$prior,4)
  # relative to prior width
iqr_width_reduction$rel_red_rank <- rank(-iqr_width_reduction$rel_red)

# Parameters with a more thsn 10% reduction in the IQR range compared to prior width
iqr_width_reduction$parm[iqr_width_reduction$rel_red>0.1]

# Posterior width compared to posterior median
posterior_iqr_width <- filter(iqr_prior_post, type=="post") %>%
  select(parm, median, iqr_width) %>%
  mutate(relative_posterior_width = iqr_width/median)

posterior_iqr_width$parm[rank(posterior_iqr_width$relative_posterior_width)>18]

# Plots
# Median and IQR prior and posterior for all parameters
ggplot(iqr_prior_post) +
  geom_point(aes(x=type, y = median)) +
  geom_errorbar(aes(x=type, ymin = lower, ymax = upper)) +
  facet_wrap(~parm, scales = "free")

# Parameters with at least a 10% reduction in IQR width and ranked by reduction
iqr_prior_post <- left_join(iqr_prior_post, select(iqr_width_reduction, parm, rel_red_rank), by = "parm")
ggplot(iqr_prior_post[iqr_prior_post$rel_red_rank <=20,]) +
  geom_point(aes(x=reorder(type, desc(type)), y = median)) +
  geom_errorbar(aes(x=type, ymin = lower, ymax = upper))+
  facet_wrap(~reorder(parm, rel_red_rank), scales = "free") +
  ylab("Parameter value") +
  xlab("") +
  scale_x_discrete(labels = c("prior" = "Prior", "post" = "Posterior")) +
  theme_classic()
# Order this somehow by largest reduction in width

# Plot posteriors with the largest IQR
posterior_long <- gather(posterior, key = "parm", value = "posterior_value")

ggplot(posterior_long) +
  geom_boxplot(aes(x=parm, y = posterior_value)) +
  facet_wrap(~parm, scales= "free")

(((0.000238 * (ages - 9))^2)  * c(rep(0, times = 9/da),rep(1, times = n_agecat - 9/da)))[41:121]
# 7.341062e-05 female 45 year old IC, 0.003517926 IR, 0.002378231 ENCHB
# 0.0013 male 45 year old IC, 0.06229759 IR, 0.04211516 ENCHB
# In mixed female cohort assuming 90% IC, 5% IR and 5% ENCHB: 0.00036
# In mixed male cohort: 0.0064


# Posterior correlation plot ----
corrmatrix_post <- cor(posterior) # methods Spearman and Pearson look similar, only Kendall focuses on the highest correlations
# Correlation plot
corrplot(corrmatrix_post, method = "circle")
# Clustered correlation plot (https://stackoverflow.com/questions/45896231/r-corrplot-with-clustering-default-dissimilarity-measure-for-correlation-matrix)
corrplot(corrmatrix_post, method = "circle", order="hclust")
# Look into different hclust.method

corrmatrix_prior <- cor(prior)
corrplot(corrmatrix_prior, method = "circle")

# Maybe here only show the subset of correlated parameter values?

# Plot posterior composite parameters ----
ages <- seq(0,100-0.5,0.5)
n_agecat <- length(ages)
da <- 0.5

# Risk of chronic carriage by age is a calibration target: Omit here
#p_chronic_function <- matrix(NA, ncol = 183, nrow = n_agecat)
#for (i in 1:183) {
#  p_chronic_function[,i] <- c(rep(posterior$p_chronic_in_mtct[i],0.5/da),
#                              exp(-posterior$p_chronic_function_r[i] *
#                                    ages[which(ages == 0.5):n_agecat]^posterior$p_chronic_function_s[i]))
#}

#p_chronic_function <- gather(data.frame(p_chronic_function), key = "sim", value = "value")
#p_chronic_function$age <- rep(ages)
# Check:
#nrow(p_chronic_function[p_chronic_function$age == 50,]) == 183
#ggplot(p_chronic_function) +
#  geom_line(aes(x=age, y = value, group = sim), col = "grey") +
#  theme_classic()

# Age-specific function of progression through IT and IR (IT=>IR and IR=>IC)
# Rate from IT to IR
eag_prog_function_it_ir <- matrix(NA, ncol = 183, nrow = n_agecat)
for (i in 1:183) {
  eag_prog_function_it_ir[,i] <- posterior$pr_it_ir[i] * exp(posterior$eag_prog_function_rate[i] * ages)
}
eag_prog_function_it_ir <- gather(data.frame(eag_prog_function_it_ir), key = "sim", value = "value")
eag_prog_function_it_ir$age <- rep(ages)
# Check:
nrow(eag_prog_function_it_ir[eag_prog_function_it_ir$age == 50,]) == 183

# Rate from IR to IC
eag_prog_function_ir_ic <- matrix(NA, ncol = 183, nrow = n_agecat)
for (i in 1:183) {
  eag_prog_function_ir_ic[,i] <- posterior$pr_ir_ic[i] * exp(posterior$eag_prog_function_rate[i] * ages)
}
eag_prog_function_ir_ic <- gather(data.frame(eag_prog_function_ir_ic), key = "sim", value = "value")
eag_prog_function_ir_ic$age <- rep(ages)
# Check:
nrow(eag_prog_function_ir_ic[eag_prog_function_ir_ic$age == 50,]) == 183

ggplot(eag_prog_function_it_ir) +
  stat_summary(aes(x=age, y = value),
               fun.min = function(x) quantile(x, 0.025),
               fun.max = function(x) quantile(x, 0.975),
               geom="ribbon", fill = "red", alpha = 0.1) +
  stat_summary(aes(x=age, y = value), fun = "median", geom="line", col = "red") +
  theme_classic() +
  xlab("Age") +
  ylab("Progression rate from HBeAg-positive infection\nto HBeAg-positive CHB\n(per person-year)") +
  ylim(0,0.3)

ggplot(eag_prog_function_ir_ic) +
  stat_summary(aes(x=age, y = value),
               fun.min = function(x) quantile(x, 0.025),
               fun.max = function(x) quantile(x, 0.975),
               geom="ribbon", fill = "red", alpha = 0.1) +
  stat_summary(aes(x=age, y = value), fun = "median", geom="line", col = "red") +
  theme_classic() +
  xlab("Age") +
  ylab("Progression rate from HBeAg-positive CHB\nto HBeAg-negative infection\n(per person-year)") +
  ylim(0,2.7)

# Observation: age effect is very small for progression through eAg-positive compartments
# Though depends also what was possible in prior

# Age-specific progression to HCC from all carrier compartments other than DCC
# Represented by a shifted quadratic function that increases with age and
# prevents people younger than 10 years to progress to HCC
# ADAPTATION 18/06/19: add an intercept to allow switching off of age dependence
#cancer_prog_function <- (parameters$cancer_prog_coefficient_female * (ages - parameters$cancer_age_threshold))^2  # Rate in females
#cancer_prog_function <- cancer_prog_function *
#  c(rep(0, times = parameters$cancer_age_threshold/da),
#    rep(1, times = n_agecat - parameters$cancer_age_threshold/da))  # Set transition to 0 in <10 year olds
#cancer_prog_female <- sapply(cancer_prog_function, function(x) min(x,1)) # Set maximum annual rate is 1
#cancer_prog_male <- sapply(parameters$cancer_male_cofactor*cancer_prog_female, function(x) min(x,1))  # Rate in males, cannot exceed 1
#cancer_prog_rates <- matrix(data = c(cancer_prog_female, cancer_prog_male),
#                            nrow = n_agecat, ncol = 2)  # store in a matrix to apply to compartment
#parameters$cancer_prog_rates <- cancer_prog_rates

# Age-specific progression from IR to CC (HBeAg-positive cirrhosis)
# ADAPTATION 18/06/19: addition of the cirrhosis male cofactor
#pr_ir_cc_function <- c(rep(0, times = parameters$pr_ir_cc_age_threshold/da),
#                       rep(parameters$pr_ir_cc_female, times = n_agecat - parameters$pr_ir_cc_age_threshold/da))
#pr_ir_cc_function_female <- sapply(pr_ir_cc_function, function(x) min(x,5))  # annual rate cannot exceed 5
#pr_ir_cc_function_male <- sapply(pr_ir_cc_function_female*parameters$cirrhosis_male_cofactor,
#                                 function(x) min(x,5))
#parameters$pr_ir_cc_function <- matrix(data = c(pr_ir_cc_function_female,
#                                                pr_ir_cc_function_male),
#                                       nrow = n_agecat, ncol = 2)  # store in a matrix to apply to compartment

# Sex-specific progression from ENCHB to CC (no age effect)
# Shimakawa 2016 found no association between current age and development of significant liver fibrosis,
# so I remove this age dependence
#pr_enchb_cc_function <- c(rep(parameters$pr_enchb_cc_female, n_agecat))  # same rate at each age
#pr_enchb_cc_function_female <- sapply(pr_enchb_cc_function, function(x) min(x,5)) # Set maximum annual rate to 5
#pr_enchb_cc_function_male <- sapply(parameters$cirrhosis_male_cofactor *
#                                      pr_enchb_cc_function_female, function(x) min(x,5))  # Rate in males, cannot exceed 5
#pr_enchb_cc_rates <- matrix(data = c(pr_enchb_cc_function_female,
#                                     pr_enchb_cc_function_male),
#                            nrow = n_agecat, ncol = 2)
#parameters$pr_enchb_cc_rates <- pr_enchb_cc_rates

# Age-specific HBsAg loss (addition by me)
# ADAPTATION 26/06/19: express this as a linear function with age based on analysis of Yusuke's data
#parameters$sag_loss <- parameters$sag_loss_slope * ages


# Explore key uncertainties in African natural history using PRCC ----
# Proportion of chronic infections due to MTCT in 2020 ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/"
chronic_infection_incidence <-
  readRDS(paste0(out_path, "sq_chronic_infection_incidence_status_quo_250920.rds"))
library(sensitivity)

# Check order of sims is the samr:
gsub("\\..*", "",
     chronic_infection_incidence$prop_mtct$sim[chronic_infection_incidence$prop_mtct$time == 2020])==
rownames(params_mat_accepted_kmeans)

# Run PRCC
prcc_prop_mtct <- pcc(params_mat_accepted_kmeans,
                      chronic_infection_incidence$prop_mtct$value[
                        chronic_infection_incidence$prop_mtct$time == 2020],
                 rank = TRUE, nboot = 100)
plot(prcc_prop_mtct)
abline(h=0)
rownames(prcc_prop_mtct$PRCC)[which(abs(prcc_prop_mtct$PRCC$original)>=quantile(abs(prcc_prop_mtct$PRCC$original),
                                                                      prob = 0.9))]
# Most influential parms in 2020 (and 2030) are:
# "b3" ; "mtct_prob_e" ; "mtct_prob_s" ; "vacc_eff"
# In 1990: "mtct_prob_e" ; "mtct_prob_s"; "p_chronic_in_mtct" ; "pr_it_ir"

# Maybe here look at parms that don't overlap with 0


### TEST: Sensitivity analysis on out2 (no treatment) ----
load(here("analysis_input/accepted_parmsets_123_180520.Rdata"))
library(sensitivity)
# Maybe try epi.prcc in epiR package which is based on the HIV paper and calculates p-value.

# PRCC
# Check order is the same
rownames(params_mat_accepted)==
  colnames(out2$timeseries$total_hbv_deaths_rate[out2$timeseries$total_hbv_deaths_rate$time == 2020,-c(1:2)])

out_vec <- c(unlist(out2$timeseries$total_hbv_deaths_rate[
  out2$timeseries$total_hbv_deaths_rate$time == 2050,-c(1:2)]))


# Run PRCC
test_prcc <- pcc(params_mat_accepted,
                 out_vec,
                 rank = TRUE, nboot = 100)
#View(test_prcc$PRCC)
rownames(test_prcc$PRCC)[which(abs(test_prcc$PRCC$original)>=quantile(abs(test_prcc$PRCC$original), prob = 0.9))]
# 1990: "pr_ir_cc_female", "pr_ir_cc_age_threshold", "pr_ic_enchb", "cirrhosis_male_cofactor"
# 2010: "pr_ir_ic", "pr_ir_cc_female", "pr_ic_enchb", "pr_enchb_cc_female"
# 2020: "pr_ir_ic", "pr_ir_cc_female", "pr_ic_enchb", "pr_enchb_cc_female"
# 2030: "mtct_prob_e", "p_chronic_function_r", "pr_ic_enchb", "pr_enchb_cc_female"
# 2050: "mtct_prob_e", "mtct_prob_s", "p_chronic_in_mtct", "pr_ic_enchb"

plot(test_prcc)
abline(h=0)
plot(x= params_mat_accepted$pr_ic_enchb, y = out_vec)  # Example of a correlated parameter
plot(x= params_mat_accepted$hccr_dcc, y = out_vec)     # Example of an uncorrelated parameter

# Assess monotonic relationships: not always clear
for(i in 1:ncol(params_mat_accepted)) {
  plot(x= params_mat_accepted[,i], y = out_vec, xlab = colnames(params_mat_accepted)[i])
}

# pr_it_ir, pr_enchb_cc_female, pr_cc_dcc, hccr_dcc are all fine
x <- params_mat_accepted$hccr_dcc
y <- out_vec
cor(rank(y), rank(x))^2 # almost 0
summary(lm(rank(y) ~ poly(rank(x), 2)))
# A low squared Spearman's rank correlation but high R-squared from such regression indicates a
# strong non-monotonic relationship.
# Need to check this is the correct way of testing this




