# Elimination analysis

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
library(epiR)
source(here("R/imperial_model_interventions.R"))
#source(here("R/scenario_analysis/calculate_outcomes.R"))

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"
out2 <- readRDS(paste0(out_path, "out2_status_quo_301120.rds"))
out2 <- out2[[1]]

# When would elimination be achieved with only continued infant vaccination?
carriers_by_age <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_carriers.rds")

# When would elimination be achieved with linear scale up of BD (in relative terms)?
#out_bd <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/elimination_analysis/testsim_bd_vacc_linear_231220.rds")
#out_bd <- out_bd[[1]]

# 80% birth dose coverage from 2020 (Shevanthi comparator scenario)
out_bd <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/elimination_analysis/elimination_bd_cov80_220421.rds")
out_bd <- out_bd[[1]]

# 80% birth dose+PPT coverage from 2020 (Shevanthi comparator scenario)
out_bd_ppt <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/elimination_analysis/elimination_bd_ppt_cov80_220421.rds")
out_bd_ppt <- out_bd_ppt[[1]]

# 80% birth dose+PPT+treatment coverage from 2020 (Shevanthi comparator scenario)
out_bd_ppt_treat <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/elimination_analysis/elimination_bd_ppt_cov80_treatment_230421.rds")
out_bd_ppt_treat <- out_bd_ppt_treat[[1]]

# When would elimination be achieved with the one-off treatment programme?
# Monitoring every 5 years til age 45
#out_treat <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/a1_it_screen_2020_monit_out7_050321.rds")
#out_treat <- out_treat[[1]]

load(here("calibration", "input", "accepted_parmsets_kmeans_170820.Rdata")) # params_mat_accepted_kmeans

# Comparison of status quo to no historical intervention is in Sensitivity analysis script ----
# Target A1: Reduce chronic infection incidence by 90% from 2015 ----

# Value in 2015
round(quantile(apply(out2$timeseries$total_chronic_infections[
  out2$timeseries$total_chronic_infections$time %in% c(2015,2015.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(apply(out2$timeseries$total_hbv_deaths[
  out2$timeseries$total_hbv_deaths$time %in% c(2015,2015.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)

# VALUES IN 2030

# New chronic infections in 2030
round(quantile(apply(out2$timeseries$total_chronic_infections[
  out2$timeseries$total_chronic_infections$time %in% c(2030,2030.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(out2$timeseries$total_chronic_infections_rate[
  out2$timeseries$total_chronic_infections$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*1000000,0)
round(quantile(apply(out_bd$timeseries$total_chronic_infections[
  out_bd$timeseries$total_chronic_infections$time %in% c(2030,2030.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(out_bd$timeseries$total_chronic_infections_rate[
  out_bd$timeseries$total_chronic_infections$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*1000000,0)
round(quantile(apply(out_bd_ppt$timeseries$total_chronic_infections[
  out_bd_ppt$timeseries$total_chronic_infections$time %in% c(2030,2030.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(out_bd_ppt$timeseries$total_chronic_infections_rate[
  out_bd_ppt$timeseries$total_chronic_infections$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*1000000,0)

round(quantile(apply(out_bd_ppt_treat$timeseries$total_chronic_infections[
  out_bd_ppt_treat$timeseries$total_chronic_infections$time %in% c(2030,2030.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(out_bd_ppt_treat$timeseries$total_chronic_infections_rate[
  out_bd_ppt_treat$timeseries$total_chronic_infections$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*1000000,0)


round(quantile(out_bd$timeseries$prev_under5[
  out_bd$timeseries$prev_under5$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*100,2)
round(quantile(out_bd_ppt$timeseries$prev_under5[
  out_bd_ppt$timeseries$prev_under5$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*100,2)
round(quantile(out_bd_ppt_treat$timeseries$prev_under5[
  out_bd_ppt_treat$timeseries$prev_under5$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*100,2)

# %Reduction by given year:
round(quantile((out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,
                                         -c(1:2)]-
  out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2030,
                                           -c(1:2)])/
  out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,
                                           -c(1:2)],
  c(0.5,0.025,0.975)),2)
# In 2030: 60% (41-74%)
# Median >=90% in 2054
# Lower percentile >=90% in year 2077

# Could also look at rate (<10 per million people)
round(quantile(out2$timeseries$total_chronic_infections_rate[
  out2$timeseries$total_chronic_infections_rate$time==2056,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)

# With BD:
round(quantile((out_bd$timeseries$total_chronic_infections[
  out_bd$timeseries$total_chronic_infections$time==2015,-c(1:2)]-
  out_bd$timeseries$total_chronic_infections[
    out_bd$timeseries$total_chronic_infections$time==2030,-c(1:2)])/
  out_bd$timeseries$total_chronic_infections[
    out_bd$timeseries$total_chronic_infections$time==2015,-c(1:2)],
  c(0.5,0.025,0.975)),2)
# In 2030: 84% (69-90%)
# Median >=90% in 2038
# Lower percentile >=90% in year 2055

# Rate (<10 per million people)
round(quantile(out_bd$timeseries$total_chronic_infections_rate[
  out_bd$timeseries$total_chronic_infections_rate$time==2067,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)

# With BD+PPT:
round(quantile((out_bd_ppt$timeseries$total_chronic_infections[
  out_bd_ppt$timeseries$total_chronic_infections$time==2015,-c(1:2)]-
    out_bd_ppt$timeseries$total_chronic_infections[
      out_bd_ppt$timeseries$total_chronic_infections$time==2030,-c(1:2)])/
    out_bd_ppt$timeseries$total_chronic_infections[
      out_bd_ppt$timeseries$total_chronic_infections$time==2015,-c(1:2)],
  c(0.5,0.025,0.975)),2)
# In 2030: 84% (69-90%)
# Median >=90% in 2038
# Lower percentile >=90% in year 2055

# Rate (<10 per million people)
round(quantile(out_bd_ppt$timeseries$total_chronic_infections_rate[
  out_bd_ppt$timeseries$total_chronic_infections_rate$time==2064,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)

# With BD+PPT+treatment:
round(quantile((out_bd_ppt_treat$timeseries$total_chronic_infections[
  out_bd_ppt_treat$timeseries$total_chronic_infections$time==2015,-c(1:2)]-
    out_bd_ppt_treat$timeseries$total_chronic_infections[
      out_bd_ppt_treat$timeseries$total_chronic_infections$time==2053,-c(1:2)])/
    out_bd_ppt_treat$timeseries$total_chronic_infections[
      out_bd_ppt_treat$timeseries$total_chronic_infections$time==2015,-c(1:2)],
  c(0.5,0.025,0.975)),2)

# Rate (<10 per million people)
round(quantile(out_bd_ppt_treat$timeseries$total_chronic_infections_rate[
  out_bd_ppt_treat$timeseries$total_chronic_infections_rate$time==2065,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)


# With treatment programme:
quantile((out_treat$timeseries$total_chronic_infections[out_treat$timeseries$total_chronic_infections$time==2015,
                                                     -c(1:2)]-
            out_treat$timeseries$total_chronic_infections[out_treat$timeseries$total_chronic_infections$time==2076,
                                                       -c(1:2)])/
           out_treat$timeseries$total_chronic_infections[out_treat$timeseries$total_chronic_infections$time==2015,
                                                      -c(1:2)], c(0.5,0.025,0.975))
# In 2030: 61% (42-76)
# Median >=90% in 2054
# Lower percentile >=90% in year 2077
# => unchanged from vacc only

# Target A2: Achieve <=0.1% HBsAg prevalence among <5 year olds ----
prev_2020 <- rowSums((do.call("rbind", lapply(lapply(carriers_by_age, "[[", "carriers_female"), function(x)
  x[which(carriers_by_age[[1]]$time==2020),which(ages==0):which(ages==4.5)])) +
    do.call("rbind", sapply(lapply(carriers_by_age, "[[", "carriers_male"), function(x)
    x[which(carriers_by_age[[1]]$time==2020),which(ages==0):which(ages==4.5)]))))/
  rowSums((do.call("rbind", lapply(lapply(carriers_by_age, "[[", "pop_female"), function(x)
    x[which(carriers_by_age[[1]]$time==2020),which(ages==0):which(ages==4.5)])) +
     do.call("rbind", sapply(lapply(carriers_by_age, "[[", "pop_male"), function(x)
       x[which(carriers_by_age[[1]]$time==2020),which(ages==0):which(ages==4.5)]))))
quantile(prev_2020*100, c(0.5,0.025,0.975))
# In 2020 the prevalence in <5 year olds is 1.09% (0.60-1.95)

year <- 2030
prev <- rowSums((do.call("rbind", lapply(lapply(carriers_by_age, "[[", "carriers_female"), function(x)
  x[which(carriers_by_age[[1]]$time==year),which(ages==0):which(ages==4.5)])) +
    do.call("rbind", sapply(lapply(carriers_by_age, "[[", "carriers_male"), function(x)
      x[which(carriers_by_age[[1]]$time==year),which(ages==0):which(ages==4.5)]))))/
  rowSums((do.call("rbind", lapply(lapply(carriers_by_age, "[[", "pop_female"), function(x)
    x[which(carriers_by_age[[1]]$time==year),which(ages==0):which(ages==4.5)])) +
      do.call("rbind", sapply(lapply(carriers_by_age, "[[", "pop_male"), function(x)
        x[which(carriers_by_age[[1]]$time==year),which(ages==0):which(ages==4.5)]))))
round(quantile(prev*100, c(0.5,0.025,0.975)),2)
# Median <=0.1% in 2055
# Upper percentile <0.1% in year 2081

# In 2000: median of 2.77
# In 2010: median of 1.99
# In 2020: median of 1.09
# In 2030: 0.49
# In 2040: 0.24
# In 2050: 0.14
# In 2060: 0.07

plot(x=c(2000,2010,2020,2030,2040,2050, 2060),
     y=c(2.77,1.99,1.09,0.49,0.24,0.14, 0.07), ylim = c(0,3))
# Prevalence declines by about 50% every 10 years but this levels
# off at very low prevalences (in absolute terms)

# With BD:
round(quantile(out_bd$timeseries$prev_under5[
  out_bd$timeseries$prev_under$time==2050,
  -c(1:2)]*100, c(0.5,0.025,0.975)),1)

# With BD+PPT:
round(quantile(out_bd_ppt$timeseries$prev_under5[
  out_bd_ppt$timeseries$prev_under$time==2046,
  -c(1:2)]*100, c(0.5,0.025,0.975)),1)

# With BD+PPT+treatment:
round(quantile(out_bd_ppt$timeseries$prev_under5[
  out_bd_ppt$timeseries$prev_under$time==2046,
  -c(1:2)]*100, c(0.5,0.025,0.975)),1)

# CDA suggests <=0.1% among 1 year olds ----
# Target B1: Reduce HBV deaths by 65% from 2015 ----

# VALUES IN 2030
round(quantile(apply(out2$timeseries$total_hbv_deaths[
  out2$timeseries$total_hbv_deaths$time %in% c(2030,2030.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(out2$timeseries$total_hbv_deaths_rate[
  out2$timeseries$total_hbv_deaths_rate$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*1000000,0)
round(quantile(apply(out_bd$timeseries$total_hbv_deaths[
  out_bd$timeseries$total_hbv_deaths$time %in% c(2030,2030.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(out_bd$timeseries$total_hbv_deaths_rate[
  out_bd$timeseries$total_hbv_deaths_rate$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*1000000,0)
round(quantile(apply(out_bd_ppt$timeseries$total_hbv_deaths[
  out_bd_ppt$timeseries$total_hbv_deaths$time %in% c(2030,2030.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(out_bd_ppt$timeseries$total_hbv_deaths_rate[
  out_bd_ppt$timeseries$total_hbv_deaths_rate$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*1000000,0)

round(quantile(apply(out_bd_ppt_treat$timeseries$total_hbv_deaths[
  out_bd_ppt_treat$timeseries$total_hbv_deaths$time %in% c(2030,2030.5),-c(1:2)],2,sum),
  c(0.5,0.025,0.975)),0)
round(quantile(out_bd_ppt_treat$timeseries$total_hbv_deaths_rate[
  out_bd_ppt_treat$timeseries$total_hbv_deaths_rate$time==2030,-c(1:2)],
  c(0.5,0.025,0.975))*1000000,0)

# HBV deaths in 2015
out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,
                                 -c(1:2)]

# %Reduction by given year:
round(quantile((out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,
                                           -c(1:2)]-
            out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2030,
                                             -c(1:2)])/
           out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,
                                            -c(1:2)],
         c(0.5,0.025,0.975)),2)
# In 2030: 10% (-7-34%)

# Median >=65% in 2062
# Lower percentile >=65% in year 2090

# <50 per million people
round(quantile(out2$timeseries$total_hbv_deaths_rate[
  out2$timeseries$total_hbv_deaths_rate$time==2054,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)

# Does BD bring forward death reduction?
round(quantile((out_bd$timeseries$total_hbv_deaths[out_bd$timeseries$total_hbv_deaths$time==2015,
                                           -c(1:2)]-
            out_bd$timeseries$total_hbv_deaths[out_bd$timeseries$total_hbv_deaths$time==2030,
                                             -c(1:2)])/
           out_bd$timeseries$total_hbv_deaths[out_bd$timeseries$total_hbv_deaths$time==2015,
                                            -c(1:2)],
         c(0.5,0.025,0.975)),2)
# In 2030: 10% (-7-34%)
# Median >=65% in 2059
# Lower percentile >=65% in year  2077
# => a little

# <50 per million people
round(quantile(out_bd$timeseries$total_hbv_deaths_rate[
  out_bd$timeseries$total_hbv_deaths_rate$time==2052,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)

# With BD+PPT
round(quantile((out_bd_ppt$timeseries$total_hbv_deaths[
  out_bd_ppt$timeseries$total_hbv_deaths$time==2015,-c(1:2)]-
                  out_bd_ppt$timeseries$total_hbv_deaths[
                    out_bd_ppt$timeseries$total_hbv_deaths$time==2030, -c(1:2)])/
                 out_bd_ppt$timeseries$total_hbv_deaths[
                   out_bd_ppt$timeseries$total_hbv_deaths$time==2015,-c(1:2)],
               c(0.5,0.025,0.975)),2)
# In 2030: 10% (-7-34%)
# Median >=65% in 2059
# Lower percentile >=65% in year  2077
# => a little

# <50 per million people
round(quantile(out_bd_ppt$timeseries$total_hbv_deaths_rate[
  out_bd_ppt$timeseries$total_hbv_deaths_rate$time==2052,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)

# With BD+PPT+treatment
round(quantile((out_bd_ppt_treat$timeseries$total_hbv_deaths[
  out_bd_ppt_treat$timeseries$total_hbv_deaths$time==2015,-c(1:2)]-
    out_bd_ppt_treat$timeseries$total_hbv_deaths[
      out_bd_ppt_treat$timeseries$total_hbv_deaths$time==2030, -c(1:2)])/
    out_bd_ppt_treat$timeseries$total_hbv_deaths[
      out_bd_ppt_treat$timeseries$total_hbv_deaths$time==2015,-c(1:2)],
  c(0.5,0.025,0.975)),2)
# <50 per million people
round(quantile(out_bd_ppt_treat$timeseries$total_hbv_deaths_rate[
  out_bd_ppt_treat$timeseries$total_hbv_deaths_rate$time==2020,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)

plot(x=out2$timeseries$total_hbv_deaths$time,
     y=apply(out2$timeseries$total_hbv_deaths[,-c(1,2)],1,median),
     xlim = c(2020,2080), type="l")
lines(x=out_bd_ppt_treat$timeseries$total_hbv_deaths$time,
     y=apply(out_bd_ppt_treat$timeseries$total_hbv_deaths[,-c(1,2)],1,median), col = "red")

# <50 per million people
round(quantile(out_bd_ppt_treat$timeseries$total_hbv_deaths_rate[
  out_bd_ppt_treat$timeseries$total_hbv_deaths_rate$time==2052,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)


# With treatment:
quantile((out_treat$timeseries$total_hbv_deaths[out_treat$timeseries$total_hbv_deaths$time==2015,
                                             -c(1:2)]-
            out_treat$timeseries$total_hbv_deaths[out_treat$timeseries$total_hbv_deaths$time==2090,
                                               -c(1:2)])/
           out_treat$timeseries$total_hbv_deaths[out_treat$timeseries$total_hbv_deaths$time==2015,
                                              -c(1:2)],
         c(0.5,0.025,0.975))
# In 2030: 49% (36-64%)

# Median >=65% in 2060
# Lower percentile >=65% in year2090
# THIS TARGET IS NOT AFFECTED BY ONE-OFF TREATMENT

# Target B2: Achieve <=3.5 HBV deaths per 100,000 ----
quantile(out2$timeseries$total_hbv_deaths_rate[
  out2$timeseries$total_hbv_deaths_rate$time==2062,-c(1:2)]*
  100000, c(0.5,0.025,0.975))
# In 2030: 7.4 (3.7-13.0)

# Median <= 3.5 in 2048
# Upper percentile <= 3.5 in 2062

# With BD
quantile(out_bd$timeseries$total_hbv_deaths_rate[
  out_bd$timeseries$total_hbv_deaths_rate$time==2059,-c(1:2)]*
    100000, c(0.5,0.025,0.975))
# In 2030: 7.4 (3.7-13.0)

# Median <= 3.5 in 2048
# Upper percentile <= 3.5 in 2059
# not as much change

# With treatment
quantile(out_treat$timeseries$total_hbv_deaths_rate[
  out_treat$timeseries$total_hbv_deaths_rate$time==2058,-c(1:2)]*
    100000, c(0.5,0.025,0.975))
# In 2030: 4 (2-7.7)  (nearly halved)

# Median <= 3.5 in 2039 => quite a bit earlier
# Upper percentile <= 3.5 in 2059 => not changed

# With this metric, treatment on average brings elimination forward
# But the problem is that one-time treatment programme even with monitoring
# goes back to same levels as vaccination eventually.
# So there would need to be some regular alternative in addition.
# Maybe in combination with BD this effect would be sustained.


# CDA suggests <= 5 per 100,000 ----

# Incidence plots ----

# Check: correlation between prop_mtct and relative impact of BD
effect <- (out2$timeseries$total_chronic_infections[
  out2$timeseries$total_chronic_infections$time==2030,-c(1:2)]-
  out_bd$timeseries$total_chronic_infections[
    out_bd$timeseries$total_chronic_infections$time==2030,-c(1:2)])/
  out2$timeseries$total_chronic_infections[
    out2$timeseries$total_chronic_infections$time==2030,-c(1:2)]

median_prop <- apply(out2$timeseries$chronic_births[
  out2$timeseries$chronic_births$time %in% seq(2020,2030,0.5),-c(1:2)]/
        out2$timeseries$total_chronic_infections[
          out2$timeseries$total_chronic_infections$time%in% seq(2020,2030,0.5),-c(1:2)],2,median)
median_prop <- out2$timeseries$chronic_births[
  out2$timeseries$chronic_births$time==2020,-c(1:2)]/
    out2$timeseries$total_chronic_infections[
      out2$timeseries$total_chronic_infections$time==2020,-c(1:2)]
plot(x=c(median_prop),
  y = effect, ylim =c(0,1), xlim=c(0,1))

cor.test(as.numeric(median_prop), as.numeric(effect))

# What is the reduction in infection incidence in 2030 and 2040 achieved with BD compared to SQ?
quantile((out2$timeseries$total_chronic_infections[
  out2$timeseries$total_chronic_infections$time==2040,-c(1:2)]-
    out_bd$timeseries$total_chronic_infections[
      out_bd$timeseries$total_chronic_infections$time==2040,-c(1:2)])/
    out2$timeseries$total_chronic_infections[
      out2$timeseries$total_chronic_infections$time==2040,-c(1:2)],
  c(0.5,0.025,0.975))

# What is the reduction in deaths in 2050 and 2070 achieved with BD compared to SQ?
quantile((out2$timeseries$total_hbv_deaths[
  out2$timeseries$total_hbv_deaths$time==2070,-c(1:2)]-
    out_bd$timeseries$total_hbv_deaths[
      out_bd$timeseries$total_hbv_deaths$time==2070,-c(1:2)])/
    out2$timeseries$total_hbv_deaths[
      out2$timeseries$total_hbv_deaths$time==2070,-c(1:2)],
  c(0.5,0.025,0.975))

# What is the reduction in deaths in 2030 and 2050 achieved with treatment scenario compared to SQ?
quantile((out2$timeseries$total_hbv_deaths[
  out2$timeseries$total_hbv_deaths$time==2050,-c(1:2)]-
    out_bd_ppt_treat$timeseries$total_hbv_deaths[
      out_bd_ppt_treat$timeseries$total_hbv_deaths$time==2050,-c(1:2)])/
    out2$timeseries$total_hbv_deaths[
      out2$timeseries$total_hbv_deaths$time==2050,-c(1:2)],
  c(0.5,0.025,0.975))

inc_summary <- gather(rbind(out2$timeseries$total_chronic_infections,
      out_bd$timeseries$total_chronic_infections,
      out_bd_ppt$timeseries$total_chronic_infections,
      out_bd_ppt_treat$timeseries$total_chronic_infections), key="sim", value = "value",
      -time,-scenario) %>%
  group_by(scenario, time) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value, 0.975))

deaths_summary <- gather(rbind(out2$timeseries$total_hbv_deaths,
                            out_bd$timeseries$total_hbv_deaths,
                            out_bd_ppt$timeseries$total_hbv_deaths,
                            out_bd_ppt_treat$timeseries$total_hbv_deaths), key="sim", value = "value",
                      -time,-scenario) %>%
  group_by(scenario, time) %>%
  summarise(median=median(value),
            lower=quantile(value, 0.025),
            upper=quantile(value, 0.975))

inc_summary$scenario <- factor(inc_summary$scenario, levels =
                                 c("status_quo", "elimination_birth_dose_80",
                                   "elimination_birth_dose_ppt_80",
                                   "elimination_birth_dose_ppt_treatment"))
deaths_summary$scenario <- factor(deaths_summary$scenario, levels =
                                 c("status_quo", "elimination_birth_dose_80",
                                   "elimination_birth_dose_ppt_80",
                                   "elimination_birth_dose_ppt_treatment"))

inc_plot <- ggplot(inc_summary) +
  geom_line(aes(x=time, y = median/0.5, colour=scenario), size=1) +
  #geom_ribbon(aes(x=time, ymin = lower, ymax=upper, fill=scenario),alpha=0.5) +
  scale_colour_manual("", labels=c("status_quo" = "Base case (infant vaccination)",
                               "elimination_birth_dose_80" = "Infant+birth dose vaccination",
                               "elimination_birth_dose_ppt_80" = "Infant+birth dose vaccination+PAP",
                               "elimination_birth_dose_ppt_treatment" = "Infant+birth dose vaccination+PAP+treatment"),
                      values=c("status_quo" = "#DC663A",
                               "elimination_birth_dose_80" = "#EBAA0D",
                               "elimination_birth_dose_ppt_80" = "#731D84",
                               "elimination_birth_dose_ppt_treatment" = "#6BA51D")) +
  ylab("Annual incident chronic HBV infections") +
  scale_x_continuous("Year", breaks=seq(2020,2080,10),
                     limits=c(2015,2080)) +
  theme_classic() +
  theme(legend.position = c(0.62, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17)) +
  ylim(0,1200)

deaths_plot <- ggplot(deaths_summary) +
  geom_line(aes(x=time, y = median/0.5, colour=scenario), size=1) +
  #geom_ribbon(aes(x=time, ymin = lower, ymax=upper, fill=scenario),alpha=0.5) +
  scale_colour_manual("", labels=c("status_quo" = "Base case (infant vaccination)",
                                   "elimination_birth_dose_80" = "Infant+birth dose vaccination",
                                   "elimination_birth_dose_ppt_80" = "Infant+birth dose vaccination+PAP",
                                   "elimination_birth_dose_ppt_treatment" = "Infant+birth dose vaccination+PAP+treatment"),
                      values=c("status_quo" = "#DC663A",
                               "elimination_birth_dose_80" = "#EBAA0D",
                               "elimination_birth_dose_ppt_80" = "#731D84",
                               "elimination_birth_dose_ppt_treatment" = "#6BA51D")) +
  ylab("Annual HBV-related deaths") +
  scale_x_continuous("Year", breaks=seq(2020,2080,10),
                     limits=c(2015,2080)) +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17)) +
  ylim(0,500)

library(grid)
inc_plot_a <- arrangeGrob(inc_plot, top = textGrob("A", x = unit(0.01, "npc"),
                                                  y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20)))
deaths_plot_b <- arrangeGrob(deaths_plot, top = textGrob("B", x = unit(0.01, "npc"),
                                                   y   = unit(1, "npc"), just=c("left","top"),
                                                   gp=gpar(col="black", fontsize=20)))

#png(file = "elimination_scenarios_plot.png", width=300, height=260, units = "mm", res=300, pointsize = 0.99)
grid.arrange(inc_plot_a, deaths_plot_b, nrow=2)
#dev.off()

## SENSITIVITY ANALYSIS -----
# PRCC for baseline projections (epiR package)----
# Possible outcomes is: reduction between 2015 and 2030 or time to elimination according to
# this definition.

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

names(baseline_inf_red_2030)==rownames(params_mat_accepted_kmeans)

baseline_inf_red_2030 <-
 (out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,-c(1:2)]-
 out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2030,-c(1:2)])/
 out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,-c(1:2)]

all_inf_red <- data.frame(time=out2$timeseries$total_chronic_infections$time[
  out2$timeseries$total_chronic_infections$time>=2015],
  sweep(sweep(out2$timeseries$total_chronic_infections[
  out2$timeseries$total_chronic_infections$time>=2015,-c(1:2)],2,
  as.numeric(out2$timeseries$total_chronic_infections[
    out2$timeseries$total_chronic_infections$time==2015,-c(1:2)]), `-`)*-1,2,
  as.numeric(out2$timeseries$total_chronic_infections[
    out2$timeseries$total_chronic_infections$time==2015,-c(1:2)]), `/`))

inf_elimination_year <- 0
for (i in 2:183) {
  inf_elimination_year[i] <- min(all_inf_red$time[which(all_inf_red[,i]>=0.9)])
}

baseline_deaths_red_2030 <-
  (out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,-c(1:2)]-
     out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2030,-c(1:2)])/
  out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,-c(1:2)]


# Calculate PRCCs
prcc_baseline_inf_red_2030 <- epi.prcc(cbind(params_mat_accepted_kmeans, t(baseline_inf_red_2030)),
                              sided.test = 2, conf.level = 0.95)
prcc_baseline_inf_elimination_year <- epi.prcc(cbind(params_mat_accepted_kmeans,
                                                     inf_elimination_year),
                                       sided.test = 2, conf.level = 0.95)
# Actually the PRCC does not seem to work/be associated with any parameters for the year of elimination!
prcc_baseline_deaths_red_2030 <- epi.prcc(cbind(params_mat_accepted_kmeans, t(baseline_deaths_red_2030)),
                                       sided.test = 2, conf.level = 0.95)


prcc_baseline_inf_red_2030 <- data.frame(parameter= colnames(params_mat_accepted_kmeans),
                     baseline_inf_red_prcc = prcc_baseline_inf_red_2030$est,
                     baseline_inf_red_p_value = prcc_baseline_inf_red_2030$p.value)
prcc_baseline_inf_red_2030 <- arrange(prcc_baseline_inf_red_2030, -abs(baseline_inf_red_prcc))
prcc_baseline_deaths_red_2030 <- data.frame(parameter= colnames(params_mat_accepted_kmeans),
                                         baseline_deaths_red_prcc = prcc_baseline_deaths_red_2030$est,
                                         baseline_deaths_red_p_value = prcc_baseline_deaths_red_2030$p.value)
prcc_baseline_deaths_red_2030 <- arrange(prcc_baseline_deaths_red_2030, -abs(baseline_deaths_red_prcc))

prcc_baseline_inf_red_2030$parameter <- factor(prcc_baseline_inf_red_2030$parameter)
levels(prcc_baseline_inf_red_2030$parameter) <- parameter_names
prcc_baseline_deaths_red_2030$parameter <- factor(prcc_baseline_deaths_red_2030$parameter)
levels(prcc_baseline_deaths_red_2030$parameter) <- parameter_names

prcc_baseline_inf_red_2030$sign <- "Positive"
prcc_baseline_inf_red_2030$sign[
  prcc_baseline_inf_red_2030$baseline_inf_red_prcc<0] <- "Negative"
prcc_baseline_deaths_red_2030$sign <- "Positive"
prcc_baseline_deaths_red_2030$sign[
  prcc_baseline_deaths_red_2030$baseline_deaths_red_prcc<0] <- "Negative"

plot(x=params_mat_accepted_kmeans$mtct_prob_s, y = baseline_inf_red_2030)
plot(x=params_mat_accepted_kmeans$b3, y = baseline_inf_red_2030)


# Plot 10 most influential parameters:
ep1 <- ggplot(subset(prcc_baseline_inf_red_2030, rank(abs(baseline_inf_red_prcc))>22)) +
  geom_col(aes(x=reorder(parameter, abs(baseline_inf_red_prcc)),
               y = abs(baseline_inf_red_prcc), fill = sign), colour="black", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual("Correlation", values = c("Negative" = "#E41A1C",
                                              "Positive" = "#377EB8")) +
  ylab("Partial rank correlation coefficient") +
  xlab("Parameter") +
  ylim(0,1) +
  #scale_x_discrete(labels = c()) +
  labs(fill="Correlation", title = "Reduction in chronic infection incidence") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = c(0.8, 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 14.5),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

ep2 <- ggplot(subset(prcc_baseline_deaths_red_2030, rank(abs(baseline_deaths_red_prcc))>22)) +
  geom_col(aes(x=reorder(parameter, abs(baseline_deaths_red_prcc)),
               y = abs(baseline_deaths_red_prcc), fill = sign), colour="black", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual("Correlation", values = c("Negative" = "#E41A1C",
                                              "Positive" = "#377EB8")) +
  ylab("Partial rank correlation coefficient") +
  xlab("Parameter") +
  ylim(0,1) +
  #scale_x_discrete(labels = c()) +
  labs(title = "Reduction in HBV-related deaths") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 14.5),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

#png(file = "elimination_baseline_prcc.png", width=310, height=130, units = "mm", res=300, pointsize = 0.99)
grid.arrange(ep1, ep2, nrow = 1)
#dev.off()

#MTCT risk from HBeAg-negative and -positive mothers, the vaccine efficacy,
#and the transmission coefficient in adults, β3
# the rate ratio for cirrhosis in men, and the progression rates to compensated cirrhosis in
# women, MTCT risk from HBeAg-negative and -positive mothers,
#the progression rate from HBeAg-positive CHB to HBeAg-negative infection and the
#rate from HBeAg-negative infection to HBeAg-negative CHB.


# TESTS
effect <- (out2$timeseries$total_chronic_infections[
  out2$timeseries$total_chronic_infections$time==2030,-c(1:2)]-
    out_bd$timeseries$total_chronic_infections[
      out_bd$timeseries$total_chronic_infections$time==2030,-c(1:2)])/
  out2$timeseries$total_chronic_infections[
    out2$timeseries$total_chronic_infections$time==2030,-c(1:2)]

effect <- (out2$timeseries$total_hbv_deaths[
  out2$timeseries$total_hbv_deaths$time==2050,-c(1:2)]-
    out_bd$timeseries$total_hbv_deaths[
      out_bd$timeseries$total_hbv_deaths$time==2050,-c(1:2)])/
  out2$timeseries$total_hbv_deaths[
    out2$timeseries$total_hbv_deaths$time==2050,-c(1:2)]

prcc <- epi.prcc(cbind(params_mat_accepted_kmeans, t(effect)),
                                       sided.test = 2, conf.level = 0.95)

prcc <- data.frame(parameter= colnames(params_mat_accepted_kmeans),
                                         prcc = prcc$est,
                                         p_value = prcc$p.value)
prcc <- arrange(prcc, -abs(prcc))


plot(x=as.numeric(median_prop), y = as.numeric(effect), xlim =c(0,1), ylim = c(0,1))


# WOULD ELIMINATION BE ACHIEVED WITH TREATMENT PROGRAMME ALONE? ----
out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/repeat_screening_anc_analysis/"

monit_out7 <- readRDS(paste0(out_path, "a1_it_screen_2020_monit_out7_050321.rds"))
monit_out7 <- monit_out7[[1]]

pop_2020_anc_2050_sim7 <- readRDS(paste0(out_path, "pop_2020_anc_2050_no_rescreen_monit_sim7_300421.rds"))
pop_2020_anc_2050_sim7 <- pop_2020_anc_2050_sim7[[1]]

# When is a 90% reduction in incidence achieved?
round(quantile((monit_out7$timeseries$total_chronic_infections[monit_out7$timeseries$total_chronic_infections$time==2015,
                                                         -c(1:2)]-
                  monit_out7$timeseries$total_chronic_infections[monit_out7$timeseries$total_chronic_infections$time==2053,
                                                           -c(1:2)])/
                 monit_out7$timeseries$total_chronic_infections[monit_out7$timeseries$total_chronic_infections$time==2015,
                                                          -c(1:2)],
               c(0.5,0.025,0.975)),2)
round(quantile((pop_2020_anc_2050_sim7$timeseries$total_chronic_infections[pop_2020_anc_2050_sim7$timeseries$total_chronic_infections$time==2015,
                                                               -c(1:2)]-
                  pop_2020_anc_2050_sim7$timeseries$total_chronic_infections[pop_2020_anc_2050_sim7$timeseries$total_chronic_infections$time==2052,
                                                                 -c(1:2)])/
                 pop_2020_anc_2050_sim7$timeseries$total_chronic_infections[pop_2020_anc_2050_sim7$timeseries$total_chronic_infections$time==2015,
                                                                -c(1:2)],
               c(0.5,0.025,0.975)),2)


round(quantile((monit_out7$timeseries$total_hbv_deaths[monit_out7$timeseries$total_hbv_deaths$time==2015,
                                                 -c(1:2)]-
                  monit_out7$timeseries$total_hbv_deaths[monit_out7$timeseries$total_hbv_deaths$time==2022,
                                                   -c(1:2)])/
                 monit_out7$timeseries$total_hbv_deaths[monit_out7$timeseries$total_hbv_deaths$time==2015,
                                                  -c(1:2)],
               c(0.5,0.025,0.975)),2)

round(quantile(monit_out7$timeseries$total_hbv_deaths_rate[
  monit_out7$timeseries$total_hbv_deaths_rate$time==2022,
  -c(1:2)]*1000000, c(0.5,0.025,0.975)),0)

plot(x=monit_out7$timeseries$total_hbv_deaths$time,
     y=apply(monit_out7$timeseries$total_hbv_deaths[,-c(1,2)],1,median),
     xlim =c(2010,2070), type="l")
lines(x=out_bd_ppt_treat$timeseries$total_hbv_deaths$time,
     y=apply(out_bd_ppt_treat$timeseries$total_hbv_deaths[,-c(1,2)],1,median),
     xlim =c(2010,2070), type="l", col = "red")

plot(x=monit_out7$timeseries$total_hbv_deaths_rate$time,
     y=apply(monit_out7$timeseries$total_hbv_deaths_rate[,-c(1,2)],1,median)*1000000,
     xlim =c(2010,2070), type="l")
lines(x=out_bd_ppt_treat$timeseries$total_hbv_deaths_rate$time,
      y=apply(out_bd_ppt_treat$timeseries$total_hbv_deaths_rate[,-c(1,2)],1,median)*1000000,
      xlim =c(2010,2070), type="l", col = "red")
