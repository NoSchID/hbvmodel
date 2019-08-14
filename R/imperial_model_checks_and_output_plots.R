#########################################################
### Imperial HBV model: model checks and output plots ###
#########################################################

### Load packages ----
require(here)
### Load main script with model, data and other functions ----
source(here("R/imperial_model_main.R"))
### Run the simulation: 1 SCENARIO (vacc/no_vacc) ----
# Default scenario: infant vaccine (apply_vacc = 1)

# Set all infection parms to zero:
#parameter_list <- lapply(parameter_list, FUN= function(x) x*0)
tic()
sim <- run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                 parms_to_change = list(b1 = 0.04, b2 = 0.04, mtct_prob_s = 0.05,
                                        mtct_prob_e = 0.6,  # decrease
                                        alpha = 10,
                                        b3 = 0.01,
                                        eag_prog_function_rate = 0,
                                        pr_it_ir = 0.1,
                                        pr_ir_ic = 0.8,
                                        pr_ir_cc_female = 0.1,
                                        pr_ir_cc_age_threshold = 30,
                                        pr_ir_enchb = 0.005,
                                        pr_ic_enchb = 0.01,
                                        pr_enchb_cc_female = 0.005, # 0.005, 0.016
                                        sag_loss_slope = 0.0004106,
                                        hccr_dcc = 0.2,  # 5 times increase
                                        hccr_it = 1,
                                        hccr_ir = 8,  # doubled
                                        hccr_enchb = 4,
                                        hccr_cc = 20,
                                        cirrhosis_male_cofactor = 5,  # increase, 20
                                        cancer_prog_coefficient_female = 0,  # doubled 0.0002
                                        cancer_age_threshold = 10,
                                        cancer_male_cofactor = 5,
                                        mu_cc = 0.005,
                                        mu_hcc = 0.8,
                                        mu_dcc = 1.5),
                 scenario = "vacc")
out <- code_model_output(sim)
toc()

outpath <- out

# Save numbers in each compartment in given year
#model_pop1880 <- out$full_output[out$time == 1880,1:(2*n_infectioncat*n_agecat)+1]
#save(model_pop1880, file = here("data/simulated_inits_1880.RData"))
#model_pop1960 <- out$full_output[221,1:(2*n_infectioncat*n_agecat)+1]
#save(model_pop1960, file = here("data/simulated_inits_1960.RData"))

### Run the simulation: 2 SCENARIOS (vacc and no_vacc) ----
tic()
out <- run_scenarios(default_parameter_list = parameter_list,
                     parms_to_change = list(b1 = 0.07, sim_starttime = starttime))
toc()


### Output plots for 2 scenarios with and without vaccine impact ----
# Number of people living with chronic HBV at each timestep with/without infant vaccine
plot(out$scenario_no_vacc$time,out$scenario_no_vacc$infectioncat_total$carriers,
     type = "l", xlim = c(1960,2100),
     xlab = "Time", ylab = "Number of people living with chronic HBV", main = "No vaccine (black) vs infant vaccine (red)")
lines(out$scenario_vacc$time, out$scenario_vacc$infectioncat_total$carriers, col = "red")

# Carrier prevalence over time with/without infant vaccine
plot(out$scenario_vacc$time, out$scenario_vacc$infectioncat_total$carriers/
       out$scenario_vacc$pop_total$pop_total,
     type = "l", xlim = c(1850,2100), ylim = c(0,0.4), col = "red",
     xlab = "Time", ylab = "Chronic carrier prevalence", main = "No vaccine (black) vs infant vaccine (red)")
lines(out$scenario_no_vacc$time,out$scenario_no_vacc$infectioncat_total$carriers/
        out$scenario_no_vacc$pop_total$pop_total, col = "black")

# Number of new cases of chronic HBV carriage at each timestep with/without infant vaccine
plot(out$scenario_no_vacc$time,
     (out$scenario_no_vacc$incident_chronic_infections$horizontal_chronic_infections+
        out$scenario_no_vacc$incident_chronic_infections$chronic_births),
     type = "l", xlim = c(1960,2100), ylim = c(0, 12000),
     xlab = "Time", ylab = "New cases of chronic HBV carriage per timestep",
     main = "No vaccine (black) vs infant vaccine (red), from MTCT = dashed")
lines(out$scenario_vacc$time,
      (out$scenario_vacc$incident_chronic_infections$horizontal_chronic_infections+
         out$scenario_vacc$incident_chronic_infections$chronic_births),
      col = "red")
lines(out$scenario_no_vacc$time,
      out$scenario_no_vacc$incident_chronic_infections$chronic_births,
      lty = "dashed")
lines(out$scenario_vacc$time,
      out$scenario_vacc$incident_chronic_infections$chronic_births,
      lty = "dashed", col = "red")

# Incidence risk of chronic HBV carriage per timestep with/without infant vaccine
plot(out$scenario_no_vacc$time,
     (out$scenario_no_vacc$incident_chronic_infections$horizontal_chronic_infections+
        out$scenario_no_vacc$incident_chronic_infections$chronic_births)*100000/
       apply(out$scenario_no_vacc$sus,1,sum),
     type = "l", xlim = c(1960,2100), ylim = c(0, 700),
     xlab = "Time", ylab = "Incidence of chronic HBV carriage per 100000 per timestep",
     main = "No vaccine (black) vs infant vaccine (red), from MTCT = dashed")
lines(out$scenario_no_vacc$time,
      out$scenario_no_vacc$incident_chronic_infections$chronic_births*100000/
        apply(out$scenario_no_vacc$sus,1,sum),
      lty = "dashed")
lines(out$scenario_vacc$time,
      (out$scenario_vacc$incident_chronic_infections$horizontal_chronic_infections+
         out$scenario_vacc$incident_chronic_infections$chronic_births)*100000/
        apply(out$scenario_vacc$sus,1,sum),
      col = "red")
lines(out$scenario_vacc$time,
      out$scenario_vacc$incident_chronic_infections$chronic_births*100000/
        apply(out$scenario_vacc$sus,1,sum),
      lty = "dashed", col = "red")

# Number of HBV-related deaths at each timestep with/without infant vaccine
plot(out$scenario_no_vacc$time, out$scenario_no_vacc$hbv_deaths$incident_number_total,
     type = "l", xlim = c(1960,2100), ylim = c(0, 1600),
     xlab = "Time", ylab = "HBV-related deaths per timestep",
     main = "No vaccine (black) vs infant vaccine (red), deaths among males = dashed")
lines(out$scenario_vacc$time, out$scenario_vacc$hbv_deaths$incident_number_total, col = "red")
lines(out$scenario_no_vacc$time, out$scenario_no_vacc$hbv_deaths$incident_number_male,
      lty = "dashed")
lines(out$scenario_no_vacc$time, out$scenario_vacc$hbv_deaths$incident_number_male,
      col = "red", lty = "dashed")

# Susceptible prevalence over time with/without infant vaccine
plot(out$scenario_vacc$time, out$scenario_vacc$infectioncat_total$sus/
       out$scenario_vacc$pop_total$pop_total,
     type = "l", xlim = c(1850,2100), ylim = c(0,1), col = "red",
     xlab = "Time", ylab = "Prevalence in susceptible compartment", main = "No vaccine (black) vs infant vaccine (red)")
lines(out$scenario_no_vacc$time,out$scenario_no_vacc$infectioncat_total$sus/
        out$scenario_no_vacc$pop_total$pop_total, col = "black")

# Immune prevalence over time with/without infant vaccine
plot(out$scenario_vacc$time, out$scenario_vacc$infectioncat_total$immune/
       out$scenario_vacc$pop_total$pop_total,
     type = "l", xlim = c(1850,2100), ylim = c(0,1), col = "red",
     xlab = "Time", ylab = "Prevalence in immune compartment", main = "No vaccine (black) vs infant vaccine (red)")
lines(out$scenario_no_vacc$time,out$scenario_no_vacc$infectioncat_total$immune/
        out$scenario_no_vacc$pop_total$pop_total, col = "black")

# Contribution of MTCT to incidence
# No vaccine
View(out$scenario_no_vacc$incident_infections %>%
       filter(time > 1950.5) %>%
       mutate(total_infections = horizontal_infections + infected_births,
              prop_mtct = infected_births/total_infections))
# Vaccine
View(out$scenario_vacc$incident_infections %>%
       filter(time > 1950.5) %>%
       mutate(total_infections = horizontal_infections + infected_births,
              prop_mtct = infected_births/total_infections))

### Output plots for 1 scenario ----
outpath <- out$scenario_vacc

## INFECTION DYNAMICS

# Total number in each infection compartment per timestep
plot(outpath$time,outpath$infectioncat_total$carriers,type = "l", ylim = c(0,4000000))
lines(outpath$time,outpath$infectioncat_total$sus,col= "red")
lines(outpath$time,outpath$infectioncat_total$immune,col= "blue")

# Proportion in each infection compartment per timestep
plot(outpath$time,outpath$infectioncat_total$carriers/outpath$pop_total$pop_total,type = "l",
     ylim = c(0,0.7), xlim = c(1880,2020))
lines(outpath$time,outpath$infectioncat_total$sus/outpath$pop_total$pop_total,col= "red")
lines(outpath$time,outpath$infectioncat_total$immune/outpath$pop_total$pop_total,col= "blue")

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

### Model checks: basics ----
outpath <- out

# Are there any negative numbers in the output?
any(unlist(outpath$full_output) < 0)

# Do susceptibles + carriers + immunes = total population (age-specific numbers)?
all.equal(outpath$sus + outpath$carriers + outpath$immune,
          outpath$pop, check.names = FALSE)
# Do susceptibles + carriers + immunes = total population (total numbers)?
all.equal(outpath$infectioncat_total$sus + outpath$infectioncat_total$carriers +
            outpath$infectioncat_total$immune,
          outpath$pop_total$pop_total, check.names = FALSE)

### Model checks: demography plots 1 ----
outpath <- out

## Total population, births, deaths

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

### Model checks: demography plots 2 - age structure ----

# Plots 1950-1990
par(mfrow=c(4,2))
# Female age structure in 1950
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 1950.5),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1950 - women")
points(x = seq(2,82,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1950"]/5,
       col = "red")
# Male age structure in 1950
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 1950.5),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1950 - men")
points(x = seq(2,82,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "1950"]/5,
       col = "red")
# Female age structure in 1970
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 1970),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1970 - women")
points(x = seq(2,82,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1970"]/5,
       col = "red")
# Male age structure in 1970
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 1970),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1970 - men")
points(x = seq(2,82,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "1970"]/5,
       col = "red")
# Female age structure in 1980
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 1980),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1980 - women")
points(x = seq(2,82,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1980"]/5,
       col = "red")
# Male age structure in 1980
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 1980),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1980 - men")
points(x = seq(2,82,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "1980"]/5,
       col = "red")
# Female age structure in 1990
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 1990),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1990 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1990"]/5,
       col = "red")
# Male age structure in 1990
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 1990),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1990 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "1990"]/5,
       col = "red")

# Plots 2000-2050
par(mfrow=c(4,2))
# Female age structure in 2000
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 2000),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2000 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2000"]/5,
       col = "red")
# Male age structure in 2000
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 2000),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2000 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2000"]/5,
       col = "red")
# Female age structure in 2010
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 2010),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2010 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2010"]/5,
       col = "red")
# Male age structure in 2010
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 2010),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2010 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2010"]/5,
       col = "red")
# Female age structure in 2020
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 2020),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2020 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2020"]/5,
       col = "red")
# Male age structure in 2020
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 2020),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2020 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2020"]/5,
       col = "red")
# Female age structure in 2050
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 2050),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2050 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2050"]/5,
       col = "red")
# Male age structure in 2050
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 2050),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2050 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2050"]/5,
       col = "red")

# Plots 2080-2100
par(mfrow=c(1,2))
# Female age structure in 2080
plot(x = ages,
     y = outpath$pop_female[which(outpath$time == 2080),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2080 - women", ylim = c(0, 55000))
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2080"]/5,
       col = "red")
# Male age structure in 2080
plot(x = ages,
     y = outpath$pop_male[which(outpath$time == 2080),index$ages_all]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2080 - men", ylim = c(0, 55000))
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2080"]/5,
       col = "red")
par(mfrow=c(1,1))


### To do: model checks ----
# devtools::test()

# DEBUGGING CODE

# Where do negative numbers appear at first timestep?
out2 <- out[2,]
min(out2)
colnames(out2)[as.numeric(out2) == min(out2)]
colnames(out2)[as.numeric(out2) <0]

colnames(out[4,])[as.numeric(out[4,]) < 0]
# ICf1, ICm1, nearly all HCC

# DEBUG: Check which compartments are negative
has.neg.col <- apply(out, 2, function(col) any(col < 0))

View(has.neg.col[has.neg.col == TRUE])
which(has.neg.col)

apply(out, 1, function(row) any(row < 0))


which(apply(out[4,], 2, function(col) any(col < 0)))

## @knitr part3
