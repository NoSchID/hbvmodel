# Elimination analysis

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
source(here("R/imperial_model_interventions.R"))
#source(here("R/scenario_analysis/calculate_outcomes.R"))

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"
out2 <- readRDS(paste0(out_path, "out2_status_quo_301120.rds"))
out2 <- out2[[1]]

# When would elimination be achieved with only continued infant vaccination?
carriers_by_age <- readRDS("C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/kmeans_full_output/out_sq_carriers.rds")

# Target A1: Reduce chronic infection incidence by 90% from 2015 ----

# New chronic infections in 2015
out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,
                                         -c(1:2)]
# %Reduction by given year:
quantile((out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,
                                         -c(1:2)]-
  out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2054,
                                           -c(1:2)])/
  out2$timeseries$total_chronic_infections[out2$timeseries$total_chronic_infections$time==2015,
                                           -c(1:2)],
  c(0.5,0.025,0.975))
# Median >=90% in 2054
# Lower percentile >=90% in year 2077

# Could also look at rate

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

year <- 2060
prev <- rowSums((do.call("rbind", lapply(lapply(carriers_by_age, "[[", "carriers_female"), function(x)
  x[which(carriers_by_age[[1]]$time==year),which(ages==0):which(ages==4.5)])) +
    do.call("rbind", sapply(lapply(carriers_by_age, "[[", "carriers_male"), function(x)
      x[which(carriers_by_age[[1]]$time==year),which(ages==0):which(ages==4.5)]))))/
  rowSums((do.call("rbind", lapply(lapply(carriers_by_age, "[[", "pop_female"), function(x)
    x[which(carriers_by_age[[1]]$time==year),which(ages==0):which(ages==4.5)])) +
      do.call("rbind", sapply(lapply(carriers_by_age, "[[", "pop_male"), function(x)
        x[which(carriers_by_age[[1]]$time==year),which(ages==0):which(ages==4.5)]))))
quantile(prev*100, c(0.5,0.025,0.975))
# Median <0.1% in 2056
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

# CDA suggests <=0.1% among 1 year olds ----
# Target B1: Reduce HBV deaths by 65% from 2015 ----
# HBV deaths in 2015
out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,
                                 -c(1:2)]

# %Reduction by given year:
quantile((out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,
                                           -c(1:2)]-
            out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2090,
                                             -c(1:2)])/
           out2$timeseries$total_hbv_deaths[out2$timeseries$total_hbv_deaths$time==2015,
                                            -c(1:2)],
         c(0.5,0.025,0.975))
# In 2030: 10% (-7-34%)

# Median >=65% in 2062
# Lower percentile >=65% in year 2090

# Target B2: Achieve <=3.5 HBV deaths per 100,000 ----
quantile(out2$timeseries$total_hbv_deaths_rate[
  out2$timeseries$total_hbv_deaths_rate$time==2062,-c(1:2)]*
  100000, c(0.5,0.025,0.975))
# In 2030: 7.4 (3.7-13.0)

# Median <= 3.5 in 20148
# Upper percentile <= 3.5 in 2062





# CDA suggests <= 5 per 100,000 ----
