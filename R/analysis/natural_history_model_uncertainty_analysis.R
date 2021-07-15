# Natural history model uncertainty analysis
library(here)
library(tidyr)
library(dplyr)
library(corrplot)
library(ggplot2)
library(sensitivity)
library(ggrepel)
library(gridExtra)

# Load parameter sets chosen based on kmeans clustering
load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans
# Load priors
load(here("calibration", "input", "lhs_samples_1000000.Rdata"))

# Load output files
carriers_by_age <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_carriers.rds")
compartments_by_age <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_compartment_prevalence_by_age.rds")
disease_outcomes_by_age <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_disease_outcomes_by_age.rds")

out2 <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/out2_status_quo_301120.rds")
out2 <- out2[[1]]

prior <- params_mat
posterior <- params_mat_accepted_kmeans


# Test for cost averted: Incident HCC cases in 1999 (ignore) ----
# Cum cases in 1992-cum cases in 1991
quantile((sapply(disease_outcomes_by_age$cum_hcc_cases_male,
       function(x) rowSums(x[which(disease_outcomes_by_age$time==1992),]))-
  sapply(disease_outcomes_by_age$cum_hcc_cases_male,
         function(x) rowSums(x[which(disease_outcomes_by_age$time==1991),])))+
  (sapply(disease_outcomes_by_age$cum_hcc_cases_female,
          function(x) rowSums(x[which(disease_outcomes_by_age$time==1992),]))-
     sapply(disease_outcomes_by_age$cum_hcc_cases_female,
            function(x) rowSums(x[which(disease_outcomes_by_age$time==1991),]))),c(0.5,0.025,0.975))

# In 2012:
quantile((sapply(disease_outcomes_by_age$cum_hcc_cases_male,
                 function(x) rowSums(x[which(disease_outcomes_by_age$time==2013),]))-
            sapply(disease_outcomes_by_age$cum_hcc_cases_male,
                   function(x) rowSums(x[which(disease_outcomes_by_age$time==2012),])))+
           (sapply(disease_outcomes_by_age$cum_hcc_cases_female,
                   function(x) rowSums(x[which(disease_outcomes_by_age$time==2013),]))-
              sapply(disease_outcomes_by_age$cum_hcc_cases_female,
                     function(x) rowSums(x[which(disease_outcomes_by_age$time==2012),]))),c(0.5,0.025,0.975))


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


# Assign new names for density plots (load parameter_names below)
prior2 <- prior
posterior2 <- posterior

parm_names_table <- data.frame(parm = names(prior2),
                               new_name = names(prior2))
parm_names_table$new_name <- factor(parm_names_table$new_name)
levels(parm_names_table$new_name) <- parameter_names

# Histograms of prior and posterior distribution
plot_prior_posterior <- function(parm) {
  plot(density(posterior2[,parm]), xlim = c(min(min(prior2[,parm]),min((posterior2[,parm]))), max(max(prior2[,parm]),max((posterior2[,parm])))),
       ylim = c(0, max(max(density(prior2[,parm])$y),max((density(posterior2[,parm])$y)))),
       main= parm_names_table$new_name[parm_names_table$parm==parm],
       lwd=3, col="red")
  if (parm %in% c("pr_ir_cc_age_threshold", "cancer_age_threshold")) {
    lines(density(prior2[,parm], bw = 1), lwd=3, lty=5, col="black")
  } else {
    lines(density(prior2[,parm]), lwd=3, lty=5, col="black")
  }
  #  legend("bottomleft", legend=c("prior density","posterior density"),
  #         col=c("blue","red"), lty=c(3,1), lwd=c(3,3), cex = 1)
}

plot(density(prior[,"pr_ir_cc_age_threshold"], bw = 1))
hist(prior[,"pr_ir_cc_age_threshold"], breaks=seq(min(prior[,"pr_ir_cc_age_threshold"])-0.5,
                                                  max(prior[,"pr_ir_cc_age_threshold"])+0.5, by=1))


#pdf(file = here("calibration", "output", "prior_posterior_density_plots_010721.pdf"), title="Model parameters: prior and posterior density")
par(mfrow=c(3,3))
plot(x=0,y=0, col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
legend("center", legend=c("Prior density","Posterior density"),
       col=c("black","red"), lty=c(5,1), lwd=c(3,3), cex=1.3)

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

# Save as png:
#png(file = here("calibration", "output", "prior_posterior_density_plots1_020721.png"),
#    width=2100, height=2800, res=300)
par(mfrow=c(5,3))
plot(x=0,y=0, col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
legend("center", legend=c("Prior density","Posterior density"),
       col=c("black","red"), lty=c(5,1), lwd=c(3,3), cex=1.3)
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
#dev.off()

#png(file = here("calibration", "output", "prior_posterior_density_plots2_020721.png"),
#    width=2100, height=3500, res=300)
par(mfrow=c(6,3))
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
#dev.off()


# Plot of prior and posterior median and 95% CrI
# Name for plots
transmission_parameters <- c("alpha", "b1", "b2", "b3", "mtct_prob_s",
                             "mtct_prob_e", "vacc_eff", "p_chronic_function_r",
                             "p_chronic_function_s", "p_chronic_in_mtct")
parameter_names <- list(
                    # Transmission parameters
                     "beta1" = "b1",
                     "beta2" = "b2",
                     "beta3" = "b3",
                     "Relative infectiousness\nin HBeAg+" = "alpha",
                     "MTCT risk from\nHBeAg+ mother" = "mtct_prob_e",
                     "MTCT risk from\nHBeAg- mother" = "mtct_prob_s",
                     "Risk of chronic\ncarriage at birth" = "p_chronic_in_mtct",
                     "Coefficient for risk\nof chronic carriage (cr)" = "p_chronic_function_r",
                     "Coefficient for risk\nof chronic carriage (cs)" = "p_chronic_function_s",
                     "Vaccine efficacy" = "vacc_eff",
                     # Natural history
                     "Rate from HBeAg+\ninfection to CHB at age 0" = "pr_it_ir",
                     "Rate from HBeAg+ CHB to\nHBeAg- infection at age 0" = "pr_ir_ic",
                     "Coefficient for progression\nthrough HBeAg+ compartments" = "eag_prog_function_rate",
                     "Rate from HBeAg+\nto HBeAg- CHB" = "pr_ir_enchb",
                     "Parameter for\nHBsAg loss" = "sag_loss_slope",
                     "Rate from HBeAg-\ninfection to CHB" = "pr_ic_enchb",
                     # Liver disease
                     "Rate from HBeAg+ CHB\nto CC in women" = "pr_ir_cc_female",
                     "Minimum age for\ncirrhosis (HBeAg+)" = "pr_ir_cc_age_threshold",
                     "Rate from HBeAg- CHB\nto CC in women" = "pr_enchb_cc_female",
                     "Rate ratio for\ncirrhosis in men" = "cirrhosis_male_cofactor",
                     "Rate of decompensation" = "pr_cc_dcc",
                     "Coefficient for progression\nto HCC in women" = "cancer_prog_coefficient_female",
                     "Minimum age for HCC" = "cancer_age_threshold",
                     "Rate ratio for\nHCC in men" = "cancer_male_cofactor",
                     "Rate ratio for HCC\nin HBeAg+ infection" = "hccr_it",
                     "Rate ratio for HCC\nin HBeAg+ CHB" = "hccr_ir",
                     "Rate ratio for HCC\nin HBeAg- CHB" = "hccr_enchb",
                     "Rate ratio for\nHCC in CC" = "hccr_cc",
                     "Rate from\nDCC to HCC" = "hccr_dcc",
                     "Mortality rate\nfrom CC" = "mu_cc",
                     "Mortality rate\nfrom DCC" = "mu_dcc",
                     "Mortality rate\nfrom HCC" = "mu_hcc")


# Rename all levels, by name
prior_posterior_summary2 <- prior_posterior_summary
prior_posterior_summary2$parameter_category <- "natural_history"
prior_posterior_summary2$parameter_category[prior_posterior_summary2$parameter %in%
                                              transmission_parameters] <- "transmission"
prior_posterior_summary2$parameter <- factor(prior_posterior_summary2$parameter)
levels(prior_posterior_summary2$parameter) <- parameter_names


pp_plot1 <- ggplot(subset(prior_posterior_summary2, parameter_category=="transmission")) +
  geom_point(aes(x=parameter, y = prior_median), col = "grey70") +
  geom_errorbar(aes(x=parameter, ymin=prior_cri_lower, ymax=prior_cri_upper), width= 0.15,
                col = "grey70", linetype="dashed") +
  geom_point(aes(x=1.1, y = post_median)) +
  geom_errorbar(aes(x=1.1, ymin=post_cri_lower, ymax=post_cri_upper), width= 0.15) +
  facet_wrap(~parameter, scales = "free", ncol =5) +
  labs(title="Transmission parameters") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.title = element_blank())

pp_plot2 <- ggplot(subset(prior_posterior_summary2, parameter_category=="natural_history")) +
  geom_point(aes(x=parameter, y = prior_median, col = "Prior"), size=2) +
  geom_errorbar(aes(x=parameter, ymin=prior_cri_lower, ymax=prior_cri_upper, col = "Prior"),
                width= 0.15,
                linetype="dashed") +
  geom_point(aes(x=1.1, y = post_median, col = "Posterior"), size = 2) +
  geom_errorbar(aes(x=1.1, ymin=post_cri_lower, ymax=post_cri_upper, col = "Posterior"),
                width= 0.15) +
  scale_colour_manual(values=c("Prior"="grey70",
                               "Posterior"="black")) +
  guides(color = guide_legend(reverse = TRUE)) +
  facet_wrap(~parameter, scales = "free", ncol =5) +
  labs(title="Natural history parameters") +
  theme_classic() +
  theme(legend.position=c(0.65,0.05),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.title = element_blank())

#tiff(file = "prior_posterior_plot.tiff", width=260, height=300, units = "mm", res=200, pointsize = 0.99)
grid.arrange(pp_plot1, pp_plot2, nrow=2, heights= c(2,5))
#dev.off()

prior_posterior_summary$prior_range <- prior_posterior_summary$prior_cri_upper-
  prior_posterior_summary$prior_cri_lower
prior_posterior_summary$posterior_range <- prior_posterior_summary$post_cri_upper-
  prior_posterior_summary$post_cri_lower
prior_posterior_summary$abs_red <- prior_posterior_summary$prior_range-
  prior_posterior_summary$posterior_range
prior_posterior_summary$rel_red <- prior_posterior_summary$abs_red/
  prior_posterior_summary$prior_range

# Analysis of the prior and posterior standard deviation/coefficient of variation ----
# Coefficient of variation = relative standard deviation

posterior_sd <- gather(params_mat_accepted_kmeans, key = "parameter", value = "sim") %>%
  group_by(parameter) %>%
  summarise(post_sd = sd(sim),
            post_cov = sd(sim)/mean(sim))

prior_sd <- gather(params_mat, key = "parameter", value = "sim") %>%
  group_by(parameter) %>%
  summarise(prior_sd = sd(sim),
            prior_cov = sd(sim)/mean(sim))

prior_posterior_sd <- left_join(prior_sd, posterior_sd, by ="parameter")
prior_posterior_sd$abs_red <- prior_posterior_sd$prior_sd-
  prior_posterior_sd$post_sd
prior_posterior_sd$rel_red <- prior_posterior_sd$abs_red/
  prior_posterior_sd$prior_sd
# Note: reduction in dispersion is calculated on reduction in standard deviation compared to
# prior SD. Overall variance in prior and posterior is shown by COV.


# Analysis of the prior and posterior interquartile range ----

posterior_iqr <- gather(params_mat_accepted_kmeans, key = "parameter", value = "sim") %>%
  group_by(parameter) %>%
  summarise(post_iqr = quantile(sim, prob = 0.75)-quantile(sim, prob = 0.25),
            post_rel_iqr = (quantile(sim, prob = 0.75)-quantile(sim, prob = 0.25))/
              median(sim))

prior_iqr <- gather(params_mat, key = "parameter", value = "sim") %>%
  group_by(parameter) %>%
  summarise(prior_iqr = quantile(sim, prob = 0.75)-quantile(sim, prob = 0.25),
            prior_rel_iqr = (quantile(sim, prob = 0.75)-quantile(sim, prob = 0.25))/
              median(sim))

prior_posterior_iqr <- left_join(prior_iqr, posterior_iqr, by ="parameter")
prior_posterior_iqr$abs_red <- prior_posterior_iqr$prior_iqr-
  prior_posterior_iqr$post_iqr
prior_posterior_iqr$rel_red <- prior_posterior_iqr$abs_red/
  prior_posterior_iqr$prior_iqr

# Kolmogorov-Smirnov Test
p <- rep(0,32)
for (i in 1:32) {
  p[i] <- ks.test(prior[,colnames(prior)==unique(colnames(prior))[i]],
                  posterior[,colnames(posterior)==unique(colnames(prior))[i]])$p.value
  names(p)[i] <- unique(colnames(prior))[i]
}
p_adjusted <- p.adjust(p, method = "holm", n = length(p))
p_adjusted[which(p_adjusted==0)] # these should be <0.001

# Combine in table
prior_posterior_summary$parameter==prior_posterior_iqr$parameter

prior_posterior_summary$parameter==colnames(prior) # FALSE

# Table should be: relative reduction in CRI, relative reduction in IQR,
# prior relative IQR, posterior relative IQR, prior median, posterior median
result_df <- data.frame(parameter=prior_posterior_summary$parameter,
                        rel_red_cri = prior_posterior_summary$rel_red, # reduction in credible interval between prior and posterior
                        rel_red_iqr = prior_posterior_iqr$rel_red,# reduction in IQR between prior and posterior
                        prior_rel_iqr = prior_posterior_iqr$prior_rel_iqr,  # relative IQR of priors
                        posterior_rel_iqr  = prior_posterior_iqr$post_rel_iqr, # relative IQR of posteriors
                        prior_median=prior_posterior_summary$prior_median,
                        posterior_median=prior_posterior_summary$post_median)

# Add KS test with Holm's correction
p_adjusted <- data.frame(parameter=names(p_adjusted), adjusted_p_value = p_adjusted)
result_df <- left_join(result_df, p_adjusted, by = "parameter")
result_df$parameter <- factor(result_df$parameter)
levels(result_df$parameter) <- parameter_names
result_df$parameter <- gsub("\n", " ", result_df$parameter)
result_df <- arrange(result_df, -rel_red_cri)
#write.csv(result_df, file=here("calibration", "output", "prior_posterior_stats_iqr_150421.csv"), row.names = FALSE)

result_df$parameter[result_df$rel_red_iqr>=0.5]
result_df$parameter[result_df$posterior_rel_iqr>=0.7]
# Is reduction associated with variability of prior? Only somewhat - some with
# higher relative IQR
# also have larger reduction.
plot(y=result_df$rel_red_iqr, x = result_df$prior_rel_iqr)
abline(h=0)
points(x=result_df$prior_rel_iqr[result_df$rel_red_iqr>=0.5],
       y=result_df$rel_red_iqr[result_df$rel_red_iqr>=0.5], col = "red")

plot(y=result_df$posterior_rel_iqr, x = result_df$prior_rel_iqr)

cor.test(result_df$prior_rel_iqr, result_df$rel_red_iqr)

View(result_df[result_df$prior_rel_iqr>1,])



## Previous calc:
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

# (ignore) Fligner Killeen test for equality of variance corrected for multiple testing and summary table ----

# Mood's median test: drop this one
median.test <- function(x, y){
  z <- c(x, y)
  g <- rep(1:2, c(length(x), length(y)))
  m <- median(z)
  fisher.test(z < m, g)$p.value
}

p <- rep(0,32)
for (i in 1:32) {
  p[i] <- median.test(prior[,i], posterior[,i])
  names(p)[i] <- colnames(prior[,i])
}
p_adjusted <- p.adjust(p, method = "holm", n = length(p))

View(cbind(colnames(posterior), p_adjusted))

# Fligner Killeen test for homogeneity of variance
# Combine prior and posterior into 1 df
prior_posterior_full_df <- rbind(cbind(type="prior", gather(prior, key = "parameter", value = "sim")),
  cbind(type="posterior", gather(posterior, key = "parameter", value = "sim")))

p <- rep(0,32)
for (i in 1:32) {
  p[i] <- fligner.test(sim ~ type, data = subset(prior_posterior_full_df,
                                                 parameter==unique(prior_posterior_full_df$parameter)[i]))$p.value
  names(p)[i] <- unique(prior_posterior_full_df$parameter)[i]
}
p_adjusted <- p.adjust(p, method = "holm", n = length(p))

# Add this manually in Excel table
p_adjusted
test <- data.frame(parameter=names(p_adjusted), adjusted_p_value = p_adjusted)
test$parameter <- factor(test$parameter)
levels(test$parameter) <- parameter_names

# Combine in table
prior_posterior_summary$parameter==prior_posterior_sd$parameter
prior_posterior_summary$parameter==colnames(prior)

result_df <- data.frame(parameter=prior_posterior_summary$parameter,
                        prior_median=prior_posterior_summary$prior_median,
                        posterior_median=prior_posterior_summary$post_median,
                        rel_red_cri = prior_posterior_summary$rel_red, # reduction in credible interval between prior and posterior
                        rel_red_sd = prior_posterior_sd$rel_red,# reduction in SD between prior and posterior
                        prior_cov = prior_posterior_sd$prior_cov,  # coefficient of variation of priors
                        posterior_cov  = prior_posterior_sd$post_cov) # coefficient of variation of posteriors

# Add Mood test with Holm's correction
p_adjusted <- data.frame(parameter=colnames(prior), adjusted_p_value = p_adjusted)
result_df <- left_join(result_df, p_adjusted, by = "parameter")
result_df$parameter <- factor(result_df$parameter)
levels(result_df$parameter) <- parameter_names
result_df <- arrange(result_df, -rel_red_cri)
#write.csv(result_df, file=here("calibration", "output", "prior_posterior_stats_300321.csv"), row.names = FALSE)

result_df$parameter[result_df$rel_red_sd>=0.5]
result_df$parameter[result_df$posterior_cov>=0.7]
# Is reduction associated with COV of prior? Yes somewhat - higher COV usually leads to larger reduction.
plot(y=result_df$rel_red_sd, x = result_df$prior_cov)
abline(h=0)
points(x=result_df$prior_cov[result_df$rel_red_sd>=0.4],
       y=result_df$rel_red_sd[result_df$rel_red_sd>=0.4], col = "red")

plot(y=result_df$posterior_cov, x = result_df$prior_cov)

cor.test(result_df$prior_cov, result_df$rel_red_sd)
cor.test(result_df$prior_cov, result_df$posterior_cov)

# Posterior correlation plot ----
corrmatrix_post <- cor(posterior) # methods Spearman and Pearson look similar, only Kendall focuses on the highest correlations

which(apply(corrmatrix_post, 2, function(r) any(r > 0.4 & r <1)))
which(apply(corrmatrix_post, 2, function(r) any(r < (-0.4))))
# Over 25% correlation:
# b1, b2, b3, mtct_prob_e, p_chronic_function_s, pr_ic_enchb, cancer_prog_coefficient_female,
# hccr_ir, hccr_enchb, hccr_cc, mu_cc
# alpha, p_chronic_function_r, pr_it_ir, eag_prog_function_rate, cirrhosis_male_cofactor
# Over 40% correlation:
# b1, b2, b3, alpha, hccr_ir, hccr_enchb, hccr_cc

# Correlation plot with over 0.4 correlation coefficient
corrmatrix_post_40 <- corrmatrix_post[rownames(corrmatrix_post) %in%
                                        c("b1", "b2", "b3", "alpha", "hccr_ir",
                                          "hccr_enchb", "hccr_cc"),colnames(corrmatrix_post) %in%
                                        c("b1", "b2", "b3", "alpha", "hccr_ir",
                                          "hccr_enchb", "hccr_cc")]
rownames(corrmatrix_post_40) <- c("beta1", "beta2", "beta3", "Relative infectiousness\nin HBeAg+",
                                  "Rate ratio for HCC\nin HBeAg+ CHB",
                                  "Rate ratio for HCC\nin HBeAg- CHB",
                                  "Rate ratio for\nHCC in CC")
colnames(corrmatrix_post_40) <- c("beta1", "beta2", "beta3", "Relative infectiousness\nin HBeAg+",
                                  "Rate ratio for HCC\nin HBeAg+ CHB",
                                  "Rate ratio for HCC\nin HBeAg- CHB",
                                  "Rate ratio for\nHCC in CC")
#tiff(file = "posterior_correlation_plot.tiff", width=300, height=300, units = "mm", res=300, pointsize = 0.99)
corrplot(corrmatrix_post_40, method = "circle", tl.col="black",
         tl.cex = 2, cl.cex = 1.5) # no clustering
#dev.off()

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

# Same with random prior sample:
set.seed(123)
random_sample <- sample(c(1:1000000), 183)
prior_sample <- prior[random_sample,]

eag_prog_function_it_ir2 <- matrix(NA, ncol = 183, nrow = n_agecat)
for (i in 1:183) {
  eag_prog_function_it_ir2[,i] <- prior_sample$pr_it_ir[i] * exp(prior_sample$eag_prog_function_rate[i] * ages)
}
eag_prog_function_it_ir2 <- gather(data.frame(eag_prog_function_it_ir2), key = "sim", value = "value")
eag_prog_function_it_ir2$age <- rep(ages)

# Rate from IR to IC
eag_prog_function_ir_ic2 <- matrix(NA, ncol = 183, nrow = n_agecat)
for (i in 1:183) {
  eag_prog_function_ir_ic2[,i] <- prior_sample$pr_ir_ic[i] * exp(prior_sample$eag_prog_function_rate[i] * ages)
}
eag_prog_function_ir_ic2 <- gather(data.frame(eag_prog_function_ir_ic2), key = "sim", value = "value")
eag_prog_function_ir_ic2$age <- rep(ages)

# Age-specific progression from IR to CC (HBeAg-positive cirrhosis)

pr_ir_cc_female_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
pr_ir_cc_female_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
pr_ir_cc_male_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
pr_ir_cc_male_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)

for (i in 1:183) {
  pr_ir_cc_female_prior[,i] <- c(rep(0, times = prior_sample$pr_ir_cc_age_threshold[i]/da),
                                                        rep(prior_sample$pr_ir_cc_female[i], times = n_agecat - prior_sample$pr_ir_cc_age_threshold[i]/da))
  pr_ir_cc_female_prior[,i] <- sapply(pr_ir_cc_female_prior[,i], function(x) min(x,5))  # annual rate cannot exceed 5
  pr_ir_cc_male_prior[,i] <-  sapply(pr_ir_cc_female_prior[,i]*
                                       prior_sample$cirrhosis_male_cofactor[i],
                                       function(x) min(x,5))

  pr_ir_cc_female_posterior[,i] <- c(rep(0, times = posterior$pr_ir_cc_age_threshold[i]/da),
                                 rep(posterior$pr_ir_cc_female[i], times = n_agecat - posterior$pr_ir_cc_age_threshold[i]/da))
  pr_ir_cc_female_posterior[,i] <- sapply(pr_ir_cc_female_posterior[,i], function(x) min(x,5))  # annual rate cannot exceed 5
  pr_ir_cc_male_posterior[,i] <-  sapply(pr_ir_cc_female_posterior[,i]*
                                           posterior$cirrhosis_male_cofactor[i],
                                     function(x) min(x,5))
}
pr_ir_cc_female_prior <- gather(data.frame(pr_ir_cc_female_prior), key = "sim", value = "value")
pr_ir_cc_female_prior$age <- rep(ages)
pr_ir_cc_male_prior <- gather(data.frame(pr_ir_cc_male_prior), key = "sim", value = "value")
pr_ir_cc_male_prior$age <- rep(ages)
pr_ir_cc_female_posterior <- gather(data.frame(pr_ir_cc_female_posterior), key = "sim", value = "value")
pr_ir_cc_female_posterior$age <- rep(ages)
pr_ir_cc_male_posterior <- gather(data.frame(pr_ir_cc_male_posterior), key = "sim", value = "value")
pr_ir_cc_male_posterior$age <- rep(ages)

# Progression from ENCHB to CC varies by sex but not by age in the model

# Combine all:
progression_functions_combined <- rbind(
  data.frame(type="Posterior", transition = "HBeAg+ infection to\nHBeAg+ CHB",
             eag_prog_function_it_ir),
  data.frame(type="Prior", transition = "HBeAg+ infection to\nHBeAg+ CHB",
             eag_prog_function_it_ir2),
  data.frame(type="Posterior", transition= "HBeAg+ CHB to\nHBeAg- infection",
             eag_prog_function_ir_ic),
  data.frame(type="Prior", transition="HBeAg+ CHB to\nHBeAg- infection",
             eag_prog_function_ir_ic2),
  data.frame(type="Posterior", transition= "HBeAg+ CHB to CC\n(female)",
             pr_ir_cc_female_posterior),
  data.frame(type="Prior", transition="HBeAg+ CHB to CC\n(female)",
             pr_ir_cc_female_prior),
  data.frame(type="Posterior", transition= "HBeAg+ CHB to CC\n(male)",
             pr_ir_cc_male_posterior),
  data.frame(type="Prior", transition="HBeAg+ CHB to CC\n(male)",
             pr_ir_cc_male_prior)
)
progression_functions_combined$transition <- factor(progression_functions_combined$transition,
                                                levels=c("HBeAg+ infection to\nHBeAg+ CHB",
                                                         "HBeAg+ CHB to\nHBeAg- infection",
                                                         "HBeAg+ CHB to CC\n(female)",
                                                         "HBeAg+ CHB to CC\n(male)"))


# Plot for thesis:
#tiff(file = "prior_posterior_progression_rates.tiff", width=300, height=180, units = "mm", res=300, pointsize = 0.99)
ggplot(progression_functions_combined) +
  stat_summary(aes(x=age, y = value, fill = type, colour=type),
               fun.min = function(x) quantile(x, 0.025),
               fun.max = function(x) quantile(x, 0.975),
               geom="ribbon", alpha = 0.1, linetype="dashed") +
  stat_summary(aes(x=age, y = value, colour=type), fun = "median", geom="line") +
  scale_fill_manual(values=c("Prior" = "grey30",
                             "Posterior" = "red")) +
  scale_colour_manual(values=c("Prior" = "black",
                               "Posterior" = "red")) +
  guides(color = guide_legend(reverse = TRUE),
         fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~transition, ncol=2, scales="free") +
  theme_classic() +
  xlab("Age (years)") +
  xlim(0,75) +
  ylab("Progression rate (per person-year)") +
  theme(legend.title=element_blank(),
      legend.text=element_text(size=13),
      strip.text = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13))
#dev.off()

# Just posterior:
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

cancer_rates_female_ic_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_ic_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_ic_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_ic_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_it_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_it_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_it_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_it_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_ir_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_ir_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_ir_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_ir_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_enchb_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_enchb_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_enchb_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_enchb_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_cc_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_female_cc_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_cc_prior <- matrix(NA, ncol = 183, nrow = n_agecat)
cancer_rates_male_cc_posterior <- matrix(NA, ncol = 183, nrow = n_agecat)

for (i in 1:183) {
  cancer_rates_female_ic_prior[,i] <- (prior_sample$cancer_prog_coefficient_female[i] * (ages - prior_sample$cancer_age_threshold[i]))^2  # Rate in females
  cancer_rates_female_ic_prior[,i]  <- cancer_rates_female_ic_prior[,i] *
    c(rep(0, times = prior_sample$cancer_age_threshold[i]/da),
      rep(1, times = n_agecat - prior_sample$cancer_age_threshold[i]/da))  # Set transition to 0 in <10 year olds
  cancer_rates_female_ic_prior[,i]  <- sapply(cancer_rates_female_ic_prior[,i], function(x) min(x,1)) # Set maximum annual rate is 1

  cancer_rates_male_ic_prior[,i] <- sapply(prior_sample$cancer_male_cofactor[i]*
                                             cancer_rates_female_ic_prior[,i], function(x) min(x,1))  # Rate in males, cannot exceed 1

  cancer_rates_female_it_prior[,i] <- prior_sample$hccr_it[i]*
                                             cancer_rates_female_ic_prior[,i]
  cancer_rates_male_it_prior[,i] <- prior_sample$hccr_it[i]*
                                               cancer_rates_male_ic_prior[,i]

  cancer_rates_female_ir_prior[,i] <- prior_sample$hccr_ir[i]*
    cancer_rates_female_ic_prior[,i]
  cancer_rates_male_ir_prior[,i] <- prior_sample$hccr_ir[i]*
    cancer_rates_male_ic_prior[,i]

  cancer_rates_female_enchb_prior[,i] <- prior_sample$hccr_enchb[i]*
    cancer_rates_female_ic_prior[,i]
  cancer_rates_male_enchb_prior[,i] <- prior_sample$hccr_enchb[i]*
    cancer_rates_male_ic_prior[,i]

  cancer_rates_female_cc_prior[,i] <- prior_sample$hccr_cc[i]*
    cancer_rates_female_ic_prior[,i]
  cancer_rates_male_cc_prior[,i] <- prior_sample$hccr_cc[i]*
    cancer_rates_male_ic_prior[,i]

  cancer_rates_female_ic_posterior[,i] <- (posterior$cancer_prog_coefficient_female[i] * (ages - posterior$cancer_age_threshold[i]))^2  # Rate in females
  cancer_rates_female_ic_posterior[,i]  <- cancer_rates_female_ic_posterior[,i] *
    c(rep(0, times = posterior$cancer_age_threshold[i]/da),
      rep(1, times = n_agecat - posterior$cancer_age_threshold[i]/da))  # Set transition to 0 in <10 year olds
  cancer_rates_female_ic_posterior[,i]  <- sapply(cancer_rates_female_ic_posterior[,i], function(x) min(x,1)) # Set maximum annual rate is 1

  cancer_rates_male_ic_posterior[,i] <- sapply(posterior$cancer_male_cofactor[i]*
                                             cancer_rates_female_ic_posterior[,i], function(x) min(x,1))  # Rate in males, cannot exceed 1

  cancer_rates_female_it_posterior[,i] <- posterior$hccr_it[i]*
    cancer_rates_female_ic_posterior[,i]
  cancer_rates_male_it_posterior[,i] <- posterior$hccr_it[i]*
    cancer_rates_male_ic_posterior[,i]

  cancer_rates_female_ir_posterior[,i] <- posterior$hccr_ir[i]*
    cancer_rates_female_ic_posterior[,i]
  cancer_rates_male_ir_posterior[,i] <- posterior$hccr_ir[i]*
    cancer_rates_male_ic_posterior[,i]

  cancer_rates_female_enchb_posterior[,i] <- posterior$hccr_enchb[i]*
    cancer_rates_female_ic_posterior[,i]
  cancer_rates_male_enchb_posterior[,i] <- posterior$hccr_enchb[i]*
    cancer_rates_male_ic_posterior[,i]

  cancer_rates_female_cc_posterior[,i] <- posterior$hccr_cc[i]*
    cancer_rates_female_ic_posterior[,i]
  cancer_rates_male_cc_posterior[,i] <- posterior$hccr_cc[i]*
    cancer_rates_male_ic_posterior[,i]

}
# Long format:
cancer_rates_female_ic_prior <- gather(data.frame(cancer_rates_female_ic_prior), key = "sim", value = "value")
cancer_rates_female_ic_prior$age <- rep(ages)
cancer_rates_male_ic_prior <- gather(data.frame(cancer_rates_male_ic_prior), key = "sim", value = "value")
cancer_rates_male_ic_prior$age <- rep(ages)
cancer_rates_female_it_prior <- gather(data.frame(cancer_rates_female_it_prior), key = "sim", value = "value")
cancer_rates_female_it_prior$age <- rep(ages)
cancer_rates_male_it_prior <- gather(data.frame(cancer_rates_male_it_prior), key = "sim", value = "value")
cancer_rates_male_it_prior$age <- rep(ages)
cancer_rates_female_ir_prior <- gather(data.frame(cancer_rates_female_ir_prior), key = "sim", value = "value")
cancer_rates_female_ir_prior$age <- rep(ages)
cancer_rates_male_ir_prior <- gather(data.frame(cancer_rates_male_ir_prior), key = "sim", value = "value")
cancer_rates_male_ir_prior$age <- rep(ages)
cancer_rates_female_enchb_prior <- gather(data.frame(cancer_rates_female_enchb_prior), key = "sim", value = "value")
cancer_rates_female_enchb_prior$age <- rep(ages)
cancer_rates_male_enchb_prior <- gather(data.frame(cancer_rates_male_enchb_prior), key = "sim", value = "value")
cancer_rates_male_enchb_prior$age <- rep(ages)
cancer_rates_female_cc_prior <- gather(data.frame(cancer_rates_female_cc_prior), key = "sim", value = "value")
cancer_rates_female_cc_prior$age <- rep(ages)
cancer_rates_male_cc_prior <- gather(data.frame(cancer_rates_male_cc_prior), key = "sim", value = "value")
cancer_rates_male_cc_prior$age <- rep(ages)

cancer_rates_female_ic_posterior <- gather(data.frame(cancer_rates_female_ic_posterior), key = "sim", value = "value")
cancer_rates_female_ic_posterior$age <- rep(ages)
cancer_rates_male_ic_posterior <- gather(data.frame(cancer_rates_male_ic_posterior), key = "sim", value = "value")
cancer_rates_male_ic_posterior$age <- rep(ages)
cancer_rates_female_it_posterior <- gather(data.frame(cancer_rates_female_it_posterior), key = "sim", value = "value")
cancer_rates_female_it_posterior$age <- rep(ages)
cancer_rates_male_it_posterior <- gather(data.frame(cancer_rates_male_it_posterior), key = "sim", value = "value")
cancer_rates_male_it_posterior$age <- rep(ages)
cancer_rates_female_ir_posterior <- gather(data.frame(cancer_rates_female_ir_posterior), key = "sim", value = "value")
cancer_rates_female_ir_posterior$age <- rep(ages)
cancer_rates_male_ir_posterior <- gather(data.frame(cancer_rates_male_ir_posterior), key = "sim", value = "value")
cancer_rates_male_ir_posterior$age <- rep(ages)
cancer_rates_female_enchb_posterior <- gather(data.frame(cancer_rates_female_enchb_posterior), key = "sim", value = "value")
cancer_rates_female_enchb_posterior$age <- rep(ages)
cancer_rates_male_enchb_posterior <- gather(data.frame(cancer_rates_male_enchb_posterior), key = "sim", value = "value")
cancer_rates_male_enchb_posterior$age <- rep(ages)
cancer_rates_female_cc_posterior <- gather(data.frame(cancer_rates_female_cc_posterior), key = "sim", value = "value")
cancer_rates_female_cc_posterior$age <- rep(ages)
cancer_rates_male_cc_posterior <- gather(data.frame(cancer_rates_male_cc_posterior), key = "sim", value = "value")
cancer_rates_male_cc_posterior$age <- rep(ages)

# Combine into dataframe: (type, sex, disease_state)
cancer_rates <-rbind(
  data.frame(type="Prior", sex="Female", disease_state="ic",
             cancer_rates_female_ic_prior),
  data.frame(type="Prior", sex="Female", disease_state="it",
             cancer_rates_female_it_prior),
  data.frame(type="Prior", sex="Female", disease_state="ir",
             cancer_rates_female_ir_prior),
  data.frame(type="Prior", sex="Female", disease_state="enchb",
             cancer_rates_female_enchb_prior),
  data.frame(type="Prior", sex="Female", disease_state="cc",
             cancer_rates_female_cc_prior),
  data.frame(type="Prior", sex="Male", disease_state="ic",
             cancer_rates_male_ic_prior),
  data.frame(type="Prior", sex="Male", disease_state="it",
             cancer_rates_male_it_prior),
  data.frame(type="Prior", sex="Male", disease_state="ir",
             cancer_rates_male_ir_prior),
  data.frame(type="Prior", sex="Male", disease_state="enchb",
             cancer_rates_male_enchb_prior),
  data.frame(type="Prior", sex="Male", disease_state="cc",
             cancer_rates_male_cc_prior),
  data.frame(type="Posterior", sex="Female", disease_state="ic",
             cancer_rates_female_ic_posterior),
  data.frame(type="Posterior", sex="Female", disease_state="it",
             cancer_rates_female_it_posterior),
  data.frame(type="Posterior", sex="Female", disease_state="ir",
             cancer_rates_female_ir_posterior),
  data.frame(type="Posterior", sex="Female", disease_state="enchb",
             cancer_rates_female_enchb_posterior),
  data.frame(type="Posterior", sex="Female", disease_state="cc",
             cancer_rates_female_cc_posterior),
  data.frame(type="Posterior", sex="Male", disease_state="ic",
             cancer_rates_male_ic_posterior),
  data.frame(type="Posterior", sex="Male", disease_state="it",
             cancer_rates_male_it_posterior),
  data.frame(type="Posterior", sex="Male", disease_state="ir",
             cancer_rates_male_ir_posterior),
  data.frame(type="Posterior", sex="Male", disease_state="enchb",
             cancer_rates_male_enchb_posterior),
  data.frame(type="Posterior", sex="Male", disease_state="cc",
             cancer_rates_male_cc_posterior)
)

cancer_rates$disease_state <- factor(cancer_rates$disease_state)
levels(cancer_rates$disease_state) <- list("HBeAg- infection"="ic",
                        "HBeAg+ infection"="it",
                        "HBeAg- CHB"="enchb",
                        "HBeAg+ CHB"="ir",
                        "CC"="cc")

#tiff(file = "prior_posterior_cancer_rates.tiff", width=300, height=130, units = "mm", res=300, pointsize = 0.99)
ggplot(cancer_rates) +
  stat_summary(aes(x=age, y = value*100, fill = type, colour=type),
               fun.min = function(x) quantile(x, 0.025),
               fun.max = function(x) quantile(x, 0.975),
               geom="ribbon", alpha = 0.1, linetype="dashed") +
  stat_summary(aes(x=age, y = value*100, colour=type), fun = "median", geom="line") +
  scale_fill_manual(values=c("Prior" = "grey30",
                             "Posterior" = "red")) +
  scale_colour_manual(values=c("Prior" = "black",
                               "Posterior" = "red")) +
  guides(color = guide_legend(reverse = TRUE),
         fill = guide_legend(reverse = TRUE)) +
  facet_wrap(sex~disease_state, ncol=5, scales="free_y") +
  theme_classic() +
  xlab("Age (years)") +
  xlim(0,75) +
  ylab("Progression rates to hepatocellular carcinoma\n(per 100 person-years)") +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=13),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))
#dev.off()

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

eag_pos_2020 <- sapply(compartments_by_age$eag_positive, function(x) x[which(compartments_by_age$time==2020),
                                                      which(seq(0,99.5,0.5)==25),])/
(sapply(lapply(carriers_by_age, "[[", "carriers_female"), function(x) x[which(carriers_by_age[[1]]$time==2020),
                                                       which(seq(0,99.5,0.5)==25),])+
  sapply(lapply(carriers_by_age, "[[", "carriers_male"), function(x) x[which(carriers_by_age[[1]]$time==2020),
                                                                    which(seq(0,99.5,0.5)==25),]))

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
  eag_pos = eag_pos_2020,
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

# k-means cluster
library(factoextra)
cluster_obj <- as.matrix(scale(pca_obj))

fviz_nbclust(cluster_obj, kmeans, method = "silhouette")

km.res <- kmeans(cluster_obj, 2, nstart = 100)
print(km.res)
km.res$size
km.res$center
km.res$cluster

plot_obj <- data.frame(prop_mtct = pca_obj$prop_mtct,
                       prop_ever_inf = pca_obj$prop_ever_inf,
                       prop_eag = pca_obj$eag_pos,
                       b1 = pca_obj$b1,
                       b3= pca_obj$b3,
                       alpha = pca_obj$alpha,
                         cluster = km.res$cluster)

plot_obj_long <- gather(plot_obj, key = "variable", value = "value", -cluster)

ggplot(plot_obj) +
  geom_boxplot(aes(x = cluster, y = prop_mtct, group = cluster))

ggplot(plot_obj) +
  geom_point(aes(x = prop_mtct, y = prop_ever_inf, group = as.factor(cluster), colour = as.factor(cluster)))

ggplot(plot_obj) +
  geom_point(aes(x = prop_mtct, y = b3, group = as.factor(cluster), colour = as.factor(cluster)))

ggplot(plot_obj) +
  geom_point(aes(x = prop_mtct, y = alpha, group = as.factor(cluster), colour = as.factor(cluster)))

ggplot(plot_obj) +
  geom_point(aes(x = prop_mtct, y = prop_eag, group = as.factor(cluster), colour = as.factor(cluster)))


ggplot(plot_obj) +
  geom_boxplot(aes(x = cluster, y = prop_eag, group = cluster))

ggplot(plot_obj_long[plot_obj_long$variable != "alpha",]) +
  geom_path(aes(x=variable, y = value, group = as.factor(cluster), colour = as.factor(cluster)),
            lineend = "round", linejoin = "round")


# Scatterplot matrix/correlation matrix
pairs(pca_obj, pch=19)

# Proportion of treatment eligibility and of cirrhosis by age ----

# % of treatment eligible carriers with cirrhosis by age in 2020
timeind <- which(compartments_by_age$time==2020)
prop_cirr_by_age <- (do.call("rbind", lapply(compartments_by_age$cc_female, function(x) x[timeind,])) +
                       do.call("rbind", lapply(compartments_by_age$cc_male, function(x) x[timeind,])) +
                       do.call("rbind", lapply(compartments_by_age$dcc_female, function(x) x[timeind,])) +
                       do.call("rbind", lapply(compartments_by_age$dcc_male, function(x) x[timeind,])))/
  (do.call("rbind", lapply(compartments_by_age$cc_female, function(x) x[timeind,])) +
     do.call("rbind", lapply(compartments_by_age$cc_male, function(x) x[timeind,])) +
     do.call("rbind", lapply(compartments_by_age$dcc_female, function(x) x[timeind,])) +
     do.call("rbind", lapply(compartments_by_age$dcc_male, function(x) x[timeind,]))+
     do.call("rbind", lapply(compartments_by_age$ir_female, function(x) x[timeind,]))+
     do.call("rbind", lapply(compartments_by_age$ir_male, function(x) x[timeind,]))+
     do.call("rbind", lapply(compartments_by_age$enchb_female, function(x) x[timeind,]))+
     do.call("rbind", lapply(compartments_by_age$enchb_male, function(x) x[timeind,])))
prop_cirr_by_age$sim <- rownames(prop_cirr_by_age)
prop_cirr_by_age <- gather(prop_cirr_by_age, key = "age", value = "prop",-sim)
prop_cirr_by_age$age <- rep(ages, each = 183)

n_cirr <- do.call("rbind", lapply(compartments_by_age$cc_female, function(x) x[timeind,])) +
  do.call("rbind", lapply(compartments_by_age$cc_male, function(x) x[timeind,])) +
  do.call("rbind", lapply(compartments_by_age$dcc_female, function(x) x[timeind,])) +
  do.call("rbind", lapply(compartments_by_age$dcc_male, function(x) x[timeind,]))
n_eligible <- do.call("rbind", lapply(compartments_by_age$cc_female, function(x) x[timeind,])) +
  do.call("rbind", lapply(compartments_by_age$cc_male, function(x) x[timeind,])) +
  do.call("rbind", lapply(compartments_by_age$dcc_female, function(x) x[timeind,])) +
  do.call("rbind", lapply(compartments_by_age$dcc_male, function(x) x[timeind,]))+
  do.call("rbind", lapply(compartments_by_age$ir_female, function(x) x[timeind,]))+
  do.call("rbind", lapply(compartments_by_age$ir_male, function(x) x[timeind,]))+
  do.call("rbind", lapply(compartments_by_age$enchb_female, function(x) x[timeind,]))+
  do.call("rbind", lapply(compartments_by_age$enchb_male, function(x) x[timeind,]))
n_it <- do.call("rbind", lapply(compartments_by_age$it_female, function(x) x[timeind,])) +
  do.call("rbind", lapply(compartments_by_age$it_male, function(x) x[timeind,]))

n_dcc <-   do.call("rbind", lapply(compartments_by_age$dcc_female, function(x) x[timeind,])) +
  do.call("rbind", lapply(compartments_by_age$dcc_male, function(x) x[timeind,]))

n_carriers <- (do.call("rbind", lapply(lapply(carriers_by_age, "[[", "carriers_female"),function(x) x[timeind,]))+
                 do.call("rbind", lapply(lapply(carriers_by_age, "[[", "carriers_male"),function(x) x[timeind,])))

quantile(apply(n_cirr[,which(ages==15):which(ages==30-da)], 1,sum)/
           apply(n_eligible[,which(ages==15):which(ages==30-da)], 1,sum), prob=c(0.5,0.025,0.975))

quantile(apply(n_cirr[,which(ages==45):which(ages==65-da)], 1,sum)/
           apply(n_eligible[,which(ages==45):which(ages==65-da)], 1,sum), prob=c(0.5,0.025,0.975))

quantile((apply(n_cirr[,which(ages==45):which(ages==65-da)], 1,sum)/
            apply(n_eligible[,which(ages==45):which(ages==65-da)], 1,sum))/
           (apply(n_cirr[,which(ages==15):which(ages==30-da)], 1,sum)/
              apply(n_eligible[,which(ages==15):which(ages==30-da)], 1,sum)), prob=c(0.5,0.025,0.975))
# Prop cirrhosis among treatment eligible is not significantly higher in 45-65 year olds than in
# 15-30 year olds

# with IT counting as eligible:
quantile(apply(n_cirr[,which(ages==15):which(ages==30-da)], 1,sum)/
           apply(n_eligible[,which(ages==15):which(ages==30-da)], 1,sum), prob=c(0.5,0.025,0.975))

quantile(apply(n_cirr[,which(ages==45):which(ages==65-da)], 1,sum)/
           (apply(n_eligible[,which(ages==45):which(ages==65-da)], 1,sum)+apply(n_it[,which(ages==45):which(ages==65-da)], 1,sum)),
         prob=c(0.5,0.025,0.975))

# DCC only
quantile(apply(n_dcc[,which(ages==15):which(ages==30-da)], 1,sum)/
           apply(n_eligible[,which(ages==15):which(ages==30-da)], 1,sum), prob=c(0.5,0.025,0.975))

quantile(apply(n_dcc[,which(ages==45):which(ages==65-da)], 1,sum)/
           apply(n_eligible[,which(ages==45):which(ages==65-da)], 1,sum), prob=c(0.5,0.025,0.975))

quantile((apply(n_dcc[,which(ages==45):which(ages==65-da)], 1,sum)/
            apply(n_eligible[,which(ages==45):which(ages==65-da)], 1,sum))/
           (apply(n_dcc[,which(ages==15):which(ages==30-da)], 1,sum)/
              apply(n_eligible[,which(ages==15):which(ages==30-da)], 1,sum)), prob=c(0.5,0.025,0.975))
# Neither is prop DCC

quantile(apply(n_cirr[,which(ages==15):which(ages==30-da)], 1,sum)/
           apply(n_carriers[,which(ages==15):which(ages==30-da)], 1,sum), prob=c(0.5,0.025,0.975))

quantile(apply(n_cirr[,which(ages==45):which(ages==65-da)], 1,sum)/
           apply(n_carriers[,which(ages==45):which(ages==65-da)], 1,sum), prob=c(0.5,0.025,0.975))

prop_cirr_carrier_by_age <- (do.call("rbind", lapply(compartments_by_age$cc_female, function(x) x[timeind,])) +
                               do.call("rbind", lapply(compartments_by_age$cc_male, function(x) x[timeind,])) +
                               do.call("rbind", lapply(compartments_by_age$dcc_female, function(x) x[timeind,])) +
                               do.call("rbind", lapply(compartments_by_age$dcc_male, function(x) x[timeind,])))/
  (do.call("rbind", lapply(lapply(carriers_by_age, "[[", "carriers_female"),function(x) x[timeind,]))+
     do.call("rbind", lapply(lapply(carriers_by_age, "[[", "carriers_male"),function(x) x[timeind,])))
prop_cirr_carrier_by_age$sim <- rownames(prop_cirr_carrier_by_age)
prop_cirr_carrier_by_age <- gather(prop_cirr_carrier_by_age, key = "age", value = "prop",-sim)
prop_cirr_carrier_by_age$age <- rep(ages, each = 183)

ggplot(prop_cirr_by_age, aes(x=age, y = prop)) +
  geom_line(aes(group = sim), col = "grey") +
  stat_summary(fun=median, geom="line", col = "red") +
  geom_vline(xintercept=15) +
  geom_vline(xintercept=30) +
  geom_vline(xintercept=45) +
  geom_vline(xintercept=65)

ggplot(prop_cirr_carrier_by_age, aes(x=age, y = prop)) +
  geom_line(aes(group = sim), col = "grey") +
  stat_summary(fun=median, geom="line", col = "red") +
  geom_vline(xintercept=15) +
  geom_vline(xintercept=30) +
  geom_vline(xintercept=45) +
  geom_vline(xintercept=65)

# PRCC for current HBV-related mortality and infection incidence ----

# Check names are the same
names(out2$timeseries$total_chronic_infections[
  out2$timeseries$total_chronic_infections$time==2020,-c(1,2)])==rownames(params_mat_accepted_kmeans)
names(out2$timeseries$total_hbv_deaths[
  out2$timeseries$total_hbv_deaths$time==2020,-c(1,2)])==rownames(params_mat_accepted_kmeans)


library(epiR)
prcc_inc <- epi.prcc(cbind(params_mat_accepted_kmeans,
                           t(out2$timeseries$total_chronic_infections[
  out2$timeseries$total_chronic_infections$time==2020,-c(1,2)])),
                              sided.test = 2, conf.level = 0.95)
prcc_inc_rate <- epi.prcc(cbind(params_mat_accepted_kmeans,
                           t(out2$timeseries$total_chronic_infections_rate[
                             out2$timeseries$total_chronic_infections_rate$time==2020,-c(1,2)])),
                     sided.test = 2, conf.level = 0.95)
# Same for incidence rate and absolute number of new cases in 2020!

prcc_mort <- epi.prcc(cbind(params_mat_accepted_kmeans,
                           t(out2$timeseries$total_hbv_deaths[
                             out2$timeseries$total_hbv_deaths$time==2020,-c(1,2)])),
                     sided.test = 2, conf.level = 0.95)

prcc_df <- data.frame(parameter= colnames(params_mat_accepted_kmeans),
                      inc_prcc= prcc_inc$est,
                      inc_p_value = prcc_inc$p.value,
                      inc_rate_prcc= prcc_inc_rate$est,
                      inc_rate_p_value = prcc_inc_rate$p.value,
                      mort_prcc= prcc_mort$est,
                      mort_p_value = prcc_mort$p.value)
prcc_df <- arrange(prcc_df, -abs(inc_prcc))
