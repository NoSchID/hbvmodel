# Natural history model uncertainty analysis
library(here)
library(tidyr)
library(dplyr)
library(corrplot)
library(ggplot2)
library(sensitivity)
library(ggrepel)

# Load parameter sets chosen based on kmeans clustering
load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans
# Load priors
load(here("calibration", "input", "lhs_samples_1000000.Rdata"))

# Load output files
carriers_by_age <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_carriers.rds")
compartments_by_age <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_compartment_prevalence_by_age.rds")
disease_outcomes_by_age <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_disease_outcomes_by_age.rds")


prior <- params_mat
posterior <- params_mat_accepted_kmeans


# Comparison of prior and posterios (initially run in run_calibration_local.R) ----
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
# Proportion of chronic infections due to MTCT PRCC analysis ----

# Note that bootstrapping calculates a 95% confidence interval, which is a
# measure of uncertainty/accuracy of computed sensitivity parameter. It
# computes the sample distribution via resampling andthen pick inner 95 % quantile.
# Interpret as: significant input parameters: 0 not contained inconfidence interval
# Resampling:generate bootstrap sample (M≈1000 times) bydrawing from existing simulated sample (sizeN) new samples ofsizeN  with replacement
# I think it picks a subset of parameter set/output combinations and calculates the PRCC. Then does this many times for different subsets.
# Maybe try epi.prcc in epiR package which is based on the HIV paper and calculates p-value.

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/"
chronic_infection_incidence <-
  readRDS(paste0(out_path, "sq_chronic_infection_incidence_status_quo_250920.rds"))

# Check order of sims is the samr:
gsub("\\..*", "",
     chronic_infection_incidence$prop_mtct$sim[chronic_infection_incidence$prop_mtct$time == 2020])==
rownames(params_mat_accepted_kmeans)

# Run PRCC for 2020
prcc_prop_mtct_2020 <- pcc(params_mat_accepted_kmeans,
                      chronic_infection_incidence$prop_mtct$value[
                        chronic_infection_incidence$prop_mtct$time == 2020],
                 rank = TRUE, nboot = 100)   # checked and there was no difference between 100 and 1000 nboot
plot(prcc_prop_mtct_2020)
abline(h=0)

# Significant parameters (95% CI not including 0):
prop_mtct_2020_significant_parms <- data.frame(parm = rownames(prcc_prop_mtct_2020$PRCC)[
  which(sign(prcc_prop_mtct_2020$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2020$PRCC$`min. c.i.`))],
  prcc_mean = prcc_prop_mtct_2020$PRCC$original[
    which(sign(prcc_prop_mtct_2020$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2020$PRCC$`min. c.i.`))],
  prcc_ci_lower = prcc_prop_mtct_2020$PRCC$`min. c.i.`[
    which(sign(prcc_prop_mtct_2020$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2020$PRCC$`min. c.i.`))],
  prcc_ci_upper = prcc_prop_mtct_2020$PRCC$`max. c.i.`[
    which(sign(prcc_prop_mtct_2020$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2020$PRCC$`min. c.i.`))])
prop_mtct_2020_significant_parms <- arrange(prop_mtct_2020_significant_parms, -abs(prcc_mean))

# Run PRCC for 2040
prcc_prop_mtct_2040 <- pcc(params_mat_accepted_kmeans,
                           chronic_infection_incidence$prop_mtct$value[
                             chronic_infection_incidence$prop_mtct$time == 2040],
                           rank = TRUE, nboot = 100)   # checked and there was no difference between 100 and 1000 nboot
plot(prcc_prop_mtct_2040)
abline(h=0)

# Significant parameters (95% CI not including 0):
prop_mtct_2040_significant_parms <- data.frame(parm = rownames(prcc_prop_mtct_2040$PRCC)[
  which(sign(prcc_prop_mtct_2040$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2040$PRCC$`min. c.i.`))],
  prcc_mean = prcc_prop_mtct_2040$PRCC$original[
    which(sign(prcc_prop_mtct_2040$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2040$PRCC$`min. c.i.`))],
  prcc_ci_lower = prcc_prop_mtct_2040$PRCC$`min. c.i.`[
    which(sign(prcc_prop_mtct_2040$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2040$PRCC$`min. c.i.`))],
  prcc_ci_upper = prcc_prop_mtct_2040$PRCC$`max. c.i.`[
    which(sign(prcc_prop_mtct_2040$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2040$PRCC$`min. c.i.`))])
prop_mtct_2040_significant_parms <- arrange(prop_mtct_2040_significant_parms, -abs(prcc_mean))

# Run PRCC for 1989
prcc_prop_mtct_1989 <- pcc(params_mat_accepted_kmeans,
                           chronic_infection_incidence$prop_mtct$value[
                             chronic_infection_incidence$prop_mtct$time == 1989],
                           rank = TRUE, nboot = 100)   # checked and there was no difference between 100 and 1000 nboot
plot(prcc_prop_mtct_1989)
abline(h=0)

# Significant parameters (95% CI not including 0):
prop_mtct_1989_significant_parms <- data.frame(parm = rownames(prcc_prop_mtct_1989$PRCC)[
  which(sign(prcc_prop_mtct_1989$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_1989$PRCC$`min. c.i.`))],
  prcc_mean = prcc_prop_mtct_1989$PRCC$original[
    which(sign(prcc_prop_mtct_1989$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_1989$PRCC$`min. c.i.`))],
  prcc_ci_lower = prcc_prop_mtct_1989$PRCC$`min. c.i.`[
    which(sign(prcc_prop_mtct_1989$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_1989$PRCC$`min. c.i.`))],
  prcc_ci_upper = prcc_prop_mtct_1989$PRCC$`max. c.i.`[
    which(sign(prcc_prop_mtct_1989$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_1989$PRCC$`min. c.i.`))])
prop_mtct_1989_significant_parms <- arrange(prop_mtct_1989_significant_parms, -abs(prcc_mean))

prop_mtct_2020_significant_parms
prop_mtct_2040_significant_parms
prop_mtct_1989_significant_parms
# Interesting that vacc_eff actually affects prop_mtct in the pre-vaccination period
# This is because value of vacc_eff and post-vaccination data influences calibration of other
# parameters. Though it's not correlated with other parms in the correlation matrix.

# Parameters affecting prop_mtct at all times:
int_vec <- intersect(intersect(prop_mtct_2020_significant_parms$parm,
          prop_mtct_2040_significant_parms$parm),
          prop_mtct_1989_significant_parms$parm)
data.frame(parm = int_vec,
           prop_1989 = subset(prop_mtct_1989_significant_parms, parm %in% int_vec)$prcc_mean,
           prop_2020 = subset(prop_mtct_2020_significant_parms, parm %in% int_vec)$prcc_mean,
           prop_2040 = subset(prop_mtct_2040_significant_parms, parm %in% int_vec)$prcc_mean)
# "mtct_prob_s", "vacc_eff", "mtct_prob_e", "p_chronic_in_mtct", "pr_it_ir"
# Note that for "p_chronic_in_mtct", the direction of association varies pre- and post-vaccination
# "p_chronic_in_mtct" and "pr_it_ir"  have a weaker association (<0.5)

# Parameters unique to the post-vaccination prop
unique_postvacc <- unique(c(setdiff(prop_mtct_2020_significant_parms$parm,
          prop_mtct_1989_significant_parms$parm),
         setdiff(prop_mtct_2040_significant_parms$parm,
        prop_mtct_1989_significant_parms$parm)))
arrange(subset(prop_mtct_2020_significant_parms, parm %in% unique_postvacc), parm)
arrange(subset(prop_mtct_2040_significant_parms, parm %in% unique_postvacc), parm)
# of these only b3 is >0.5

# Parameters unique to the pre-vaccination prop:
unique_prevacc <- unique(c(setdiff(prop_mtct_1989_significant_parms$parm,
        prop_mtct_2020_significant_parms$parm),
        setdiff(prop_mtct_1989_significant_parms$parm,
        prop_mtct_2040_significant_parms$parm)))
arrange(subset(prop_mtct_1989_significant_parms, parm %in% unique_prevacc), parm)
# "b1", "alpha", "b2", "pr_ir_cc_female", "pr_ir_ic",
# but none of these is >0.5 (including CI)

# Change in prop_mtct between 1989 and 2020
prcc_prop_mtct_2020_increase <- pcc(params_mat_accepted_kmeans,
                           chronic_infection_incidence$prop_mtct$value[
                             chronic_infection_incidence$prop_mtct$time == 2020]-chronic_infection_incidence$prop_mtct$value[
                               chronic_infection_incidence$prop_mtct$time == 1989],
                           rank = TRUE, nboot = 100)   # checked and there was no difference between 100 and 1000 nboot
plot(prcc_prop_mtct_2020_increase)
abline(h=0)

# Significant parameters (95% CI not including 0):
prop_mtct_2020_increase_significant_parms <- data.frame(parm = rownames(prcc_prop_mtct_2020_increase$PRCC)[
  which(sign(prcc_prop_mtct_2020_increase$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2020_increase$PRCC$`min. c.i.`))],
  prcc_mean = prcc_prop_mtct_2020_increase$PRCC$original[
    which(sign(prcc_prop_mtct_2020_increase$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2020_increase$PRCC$`min. c.i.`))],
  prcc_ci_lower = prcc_prop_mtct_2020_increase$PRCC$`min. c.i.`[
    which(sign(prcc_prop_mtct_2020_increase$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2020_increase$PRCC$`min. c.i.`))],
  prcc_ci_upper = prcc_prop_mtct_2020_increase$PRCC$`max. c.i.`[
    which(sign(prcc_prop_mtct_2020_increase$PRCC$`max. c.i.`)==sign(prcc_prop_mtct_2020_increase$PRCC$`min. c.i.`))])
prop_mtct_2020_increase_significant_parms <- arrange(prop_mtct_2020_increase_significant_parms, -abs(prcc_mean))

# Compare PRCC between looking at prop in 2020 vs change between 1989 and 2020
prop_mtct_2020_increase_significant_parms
prop_mtct_2020_significant_parms
df <- full_join(select(prop_mtct_2020_significant_parms, parm, prcc_mean),
                select(prop_mtct_2020_increase_significant_parms, parm, prcc_mean),
                by = "parm")
colnames(df) <- c("parm", "mean_2020", "mean_2020_change")
df$class <- ifelse(abs(df$mean_2020 - df$mean_2020_change) >0.05, "change", "no_change")
df$class[is.na(df$class)] <- "not_significant"

# Slope chart
ggplot(df) +
  geom_point(aes(x=1, y=mean_2020, col = class), size=2, show.legend=F) +
  geom_point(aes(x=2, y=mean_2020_change, col = class), size=2, show.legend=F) +
  geom_segment(aes(x=1, xend=2, y=mean_2020, yend=mean_2020_change, col = class), size=.75, show.legend=F) +
  geom_vline(xintercept=1, linetype="dashed", size=.1) +
  geom_vline(xintercept=2, linetype="dashed", size=.1) +
  labs(x="", y="PRCC") +  # Axis labels
  ylim(-0.6,1) +
  xlim(.5, 2.5) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("change"="red", "no_change"="darkblue", "not_significant" = "red"))  +
  geom_text_repel(label=df$parm, y=df$mean_2020, x=rep(1, NROW(df)), hjust=1.1, size=3.5, direction = "y") +
  geom_text_repel(label=df$parm, y=df$mean_2020_change, x=rep(2, NROW(df)), hjust=-0.1, size=3.5, direction = "y") +
  geom_text(label="2020", x=1, y=0.9, hjust=1.2, size=5) +  # title
  geom_text(label="Change 1989-2020", x=2, y=0.9, hjust=-0.1, size=5) +
  theme_bw()
# Main difference is that mtct_prob_s, mtct_prob_e and to a lesser degree p_chronic_function_s and_r
# have more influence on prop in 2020 than on the change in prop 1989-2020
# Parameters influencing 2020 but not change: p_chronic_in_mtct, pr_ir_enchb, pr_it_ir
# Parameters influencing change but not 2020: b1, b2, alpha, sag_loss_slope, mu_cc
# Most influential for both (>0.5): vacc_eff, b3

# Correlated parameter example
plot(x = params_mat_accepted_kmeans$mtct_prob_s,
      y = chronic_infection_incidence$prop_mtct$value[
        chronic_infection_incidence$prop_mtct$time == 1989])

# Assess that all parameters have monotonic relationship to outcome
# Most commonly there is no obvious relationship rather than monotonic/non-monotonic one
out_vec <- chronic_infection_incidence$prop_mtct$value[
  chronic_infection_incidence$prop_mtct$time == 2020]
for(i in 1:ncol(params_mat_accepted_kmeans)) {
  plot(x= params_mat_accepted_kmeans[,i], y = out_vec, xlab = colnames(params_mat_accepted_kmeans)[i])
}

#x <- params_mat_accepted_kmeans$b1
#y <- out_vec
#cor(rank(y), rank(x))^2 # almost 0
#summary(lm(rank(y) ~ poly(rank(x), 2)))
# A low squared Spearman's rank correlation but high R-squared from such regression indicates a
# strong non-monotonic relationship.
# Need to check this is the correct way of testing this

# Proportion of chronic infections due to MTCT: Qualitative comparison of posterior distributions ----
par(mfrow=c(1,2))
hist(params_mat_accepted_kmeans$vacc_eff[
  which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]>=0.75)],
  xlim = c(0.5,1))
hist(params_mat_accepted_kmeans$vacc_eff[
  which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]<0.75)],
  xlim = c(0.5,1))

# Kolmogorov-Smirnov test
ks.test(params_mat_accepted_kmeans$vacc_eff[
  which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]>=0.75)],
  params_mat_accepted_kmeans$vacc_eff[
    which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]<0.75)])
# p-value = 0.007189
# There is some evidence to reject the null hypothesis that the 2 samples are drawn from
# the same continuous distribution

hist(params_mat_accepted_kmeans$b3[
  which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]>=0.75)],
  xlim = c(0,0.5))
hist(params_mat_accepted_kmeans$b3[
  which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]<0.75)],
  xlim = c(0,0.5))

ks.test(params_mat_accepted_kmeans$b3[
  which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]>=0.75)],
  params_mat_accepted_kmeans$b3[
    which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]<0.75)])
# p-value = 0.001452

# Perform a Kolmogorov-Smirnov test on all posterior distributions, distinguishing between
# 2020 MTCT prob contribution of >/<75%
p_values <- data.frame(parm = 0, p_value = 0)
for(i in 1:ncol(params_mat_accepted_kmeans)) {
  x[i] <- ks.test(params_mat_accepted_kmeans[,i][
    which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]>=0.75)],
    params_mat_accepted_kmeans[,i][
      which(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]<0.75)])$p.value
  p_values[i,] <- c(colnames(params_mat_accepted_kmeans)[i], x[i])
}
p_values[p_values$p_value<=0.05,]




# Correlation of proportion MTCT with HBeAg prevalence in pregnant women and proportion ever infected ----

ever_inf_1990 <- do.call("rbind",
    lapply(compartments_by_age$ever_inf_female, function(x) x[which(compartments_by_age$time==1990),]))+
  do.call("rbind",
          lapply(compartments_by_age$ever_inf_male, function(x) x[which(compartments_by_age$time==1990),]))

pop_female <- lapply(carriers_by_age, "[[", "pop_female")
pop_male <- lapply(carriers_by_age, "[[", "pop_male")

pop_1990 <- do.call("rbind",
                         lapply(pop_female, function(x) x[which(compartments_by_age$time==1990),]))+
  do.call("rbind",
          lapply(pop_male, function(x) x[which(compartments_by_age$time==1990),]))

prop_ever_infected <- ever_inf_1990/pop_1990
prop_ever_infected$sim <- rownames(prop_ever_infected)
prop_ever_infected <- gather(prop_ever_infected, key = "age", value = "value", -sim)
prop_ever_infected$age <- rep(seq(0,99.5,0.5), each = 183)

prop_ever_infected$mtct_prob_2020 <- "Low"
prop_ever_infected$mtct_prob_2020[prop_ever_infected$sim %in%
                                    gsub("\\..*", "", chronic_infection_incidence$prop_mtct$sim[
                                      chronic_infection_incidence$prop_mtct$time == 2020 &
                                        chronic_infection_incidence$prop_mtct$value >=0.65])] <- "High"

prop_ever_infected$mtct_prob_2040 <- "Low"
prop_ever_infected$mtct_prob_2040[prop_ever_infected$sim %in%
                                    gsub("\\..*", "", chronic_infection_incidence$prop_mtct$sim[
                                      chronic_infection_incidence$prop_mtct$time == 2040 &
                                        chronic_infection_incidence$prop_mtct$value >=0.5])] <- "High"

ggplot(prop_ever_infected) +
  geom_line(aes(x=age, y = value, group = sim, colour = mtct_prob_2040))

quantile(prop_ever_infected$value[prop_ever_infected$age == 50 & prop_ever_infected$mtct_prob_2020 == "High"])
quantile(prop_ever_infected$value[prop_ever_infected$age == 50 & prop_ever_infected$mtct_prob_2020 == "Low"])
quantile(prop_ever_infected$value[prop_ever_infected$age == 30 & prop_ever_infected$mtct_prob_2040 == "High"])
quantile(prop_ever_infected$value[prop_ever_infected$age == 30 & prop_ever_infected$mtct_prob_2040 == "Low"])


par(mfrow=c(1,2))
# Prop MTCT in 2020
boxplot(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2020 & gsub("\\..*", "", chronic_infection_incidence$prop_mtct$sim) %in%
                                              prop_ever_infected$sim[which(prop_ever_infected$age==50 & prop_ever_infected$value>=0.85)]],
        ylim= c(0.25,1), ylab = "Prop. chronic infections due to MTCT in 2020",
        main = "Prop. ever infected >=85% at age 50\nin 1990")
boxplot(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2020 & gsub("\\..*", "", chronic_infection_incidence$prop_mtct$sim) %in%
                                                        prop_ever_infected$sim[which(prop_ever_infected$age==50 & prop_ever_infected$value<0.85)]],
        ylim= c(0.25,1), main = "Prop. ever infected <85% at age 50\nin 1990")
# Proportion of chronic infections due to MTCT in 2020 is slightly lower in simulations where
# prop. ever infected in 50 year olds in 1990 was calibrated as over 85%
# 85% also corresponds to about the lower quartile

# Prop MTCT in 2040
boxplot(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040 & gsub("\\..*", "", chronic_infection_incidence$prop_mtct$sim) %in%
                                                      prop_ever_infected$sim[which(prop_ever_infected$age==50 & prop_ever_infected$value>=0.85)]],
        ylim= c(0,1), ylab = "Prop. chronic infections due to MTCT in 2040",
        main = "Prop. ever infected >=85% at age 50\nin 1990")
boxplot(chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040 & gsub("\\..*", "", chronic_infection_incidence$prop_mtct$sim) %in%
                                                      prop_ever_infected$sim[which(prop_ever_infected$age==50 & prop_ever_infected$value<0.85)]],
        ylim= c(0,1), main = "Prop. ever infected <85% at age 50\nin 1990")

# Look at direct correlation
par(mfrow=c(1,1))
plot(x = prop_ever_infected$value[prop_ever_infected$age == 50],
     y = chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040],
     xlim = c(0,1), ylim = c(0,1),
     ylab = "Proportion of new chronic infections due to MTCT in 2040",
     xlab = "Anti-HBc prevalence in 50 year olds in 1990")
# But note b3 also correlates with alpha and b1 so there might be a further relationship here
# Could do PCA or machine learning approach to find relationship between anti-HBc prev, prop MTCT
# b3, b1 and alpha

# PCA
pca_obj <- data.frame(
  prop_ever_inf = prop_ever_infected$value[prop_ever_infected$age == 50],
  prop_mtct = chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040],
  b1 = params_mat_accepted_kmeans$b1,
  b3 = params_mat_accepted_kmeans$b3,
  alpha = params_mat_accepted_kmeans$alpha)

pca <- prcomp(pca_obj, center = TRUE, scale = TRUE)
summary(pca)

library(factoextra)
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

mtct_group <- chronic_infection_incidence$prop_mtct$value[chronic_infection_incidence$prop_mtct$time == 2040]
mtct_group[mtct_group <0.5] <- "<50%"
mtct_group[mtct_group >=0.5] <- ">50%"


fviz_pca_ind(pca, geom.ind = "point", pointshape = 21,
             pointsize = 2,
             fill.ind = mtct_group,
             col.ind = "black",
             palette = "jco",
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Proportion of new chronic infections\ndue to MTCT in 2040") +
  ggtitle("2D PCA-plot from 30 feature dataset") +
  theme(plot.title = element_text(hjust = 0.5))
