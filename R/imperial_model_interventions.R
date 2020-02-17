#####################################
### Imperial HBV projection model ###
#####################################
# Model based on imperial_model_main 07/10/19, copied on 07/01/2020
# Copied over the natural history model used in calibration for addition of PMTCT and treatment

### Load packages ----
require(here)  # for setting working directory
require(tidyr)  # for data processing
require(dplyr)  # for data processing
require(deSolve)  # ODE solver
#require(lhs)  # Latin Hypercube sampling
#require(parallel)  # for running code in parallel

# Packages for making plots
#require(ggplot2)  # for calibration plots
#require(gridExtra)  # for calibration plots
#require(grid)  # for calibration plots

# Packages for timing, profiling and making code more efficient
#library(tictoc)  # for timing code
#library(profvis)  # for profiling
#library(efficient)  # for profiling


### Define simulation parameters ----
## Country
countryname <- "gambia"

## Times
dt <- 0.5                              # timestep (years)
# dt/da can only be 0.5 to match the WAIFW matrix at the moment
# For demography, it can be up to 1, or multiple of 5 thereafter (because of women of childbearing age)
starttime <- 1850
runtime <- 251                     # number of years to run the model for
#times <- round((0:(runtime/dt))*dt,2) # vector of timesteps
#times_labels <- times+starttime       # year labels for timestep vector

## Age groups
da <- dt                              # time spent in each age group (years)
ages <- round((0:((100-da)/da))*da,2)  # vector of all age groups
n_agecat <- length(ages)               # number of age groups

ages_wocba <- round(((15/da):((50-da)/da))*da,2)        # age groups 15-49 years (women of childbearing age)

## Infection compartments
n_nathistcat <- 9
n_screencat <- 9
n_treatcat <- 5
n_infectioncat <- n_nathistcat+n_screencat+n_treatcat
# Natural history (calibration) model has 9 infection compartments (S, IT, IR, IC, ENCHB, CC, DCC, HCC, R)
# Now added 9 screened compartments + 5 ever treated compartments

## Definition of indices
index <- list("infcat_all" = 1:n_infectioncat,            # index for all infection status compartments
              "ages_all" = 1:n_agecat,                    # index for all age groups
              "ages_wocba" = which(ages == min(ages_wocba)):which(ages == max(ages_wocba)), # index for age group 15-49 years (women of childbearing age)
              "ages_0to1" = which(ages == 0):(1/da),
              "ages_1to5" = which(ages == 1):which(ages == 6-da),       # index for age groups 1-5 years not in use
              "ages_6to15" = which(ages == 6):which(ages == 16-da),     # index for age groups 6-15 years not in use
              "ages_16to100" = which(ages == 16):n_agecat)              # index for age groups 16-100 years not in use

### Load and clean demographic data ----
# Load inputs
#source(here("R/imperial_model_load.R"))
# Clean demographic data
#source(here("R/imperial_model_clean_demography.R"))

# Load preformatted data for dt = da = 0.5
if (da == 0.5) {
  load(here("data/demogdata_0point5.RData"))
} else {
  print("Error: can only preload data for 0.5 time/agestep")
}

### Load and clean vaccine coverage data ----
# Data downloaded from WHO reported vaccine coverage timeseries (HepB3)
input_who_vaccine_coverage <- read.csv(here("data-raw", "infant_vaccine_coverage.csv"),
                                       stringsAsFactors = FALSE)

# Linearly interpolate coverage for missing years
input_who_vaccine_coverage$coverage_interp <- approx(input_who_vaccine_coverage$year, input_who_vaccine_coverage$coverage_proportion,
                                                     xout = input_who_vaccine_coverage$year, method = "linear", rule = 2)$y

# Fill in a dataframe with 0.01 timesteps with constant yearly coverage values until 2100
# Assuming last coverage from 2017 stays the same
vaccine_coverage <- data.frame(year = seq(min(input_who_vaccine_coverage$year), 2100, 0.1),
                               coverage = approx(input_who_vaccine_coverage$year, input_who_vaccine_coverage$coverage_interp,
                                                 xout = seq(min(input_who_vaccine_coverage$year), 2100, 0.1),
                                                 method = "constant", rule = 2)$y)
vaccine_coverage_list <- as.data.frame(list(times = vaccine_coverage$year,
                                            coverage = vaccine_coverage$coverage))
timevary_vaccine_coverage <- approxfun(vaccine_coverage_list, method = "linear", rule = 2)

### Infection data preparation ----
gambia_prevdata <- read.csv(here("testdata", "edmunds_gambia_prev.csv"), stringsAsFactors = FALSE)

# Interpolate prevalence and prop. ever infected
gambia_prev <- approx(x = gambia_prevdata$age, y = gambia_prevdata$edmunds_prev,
                      xout = ages, method = "linear", rule = 2)
gambia_prev <- data.frame(age = gambia_prev$x, prev = gambia_prev$y)
gambia_ever_inf <- approx(x = gambia_prevdata$age, y = gambia_prevdata$edmunds_prop_ever_infected,
                          xout = ages, method = "linear", rule = 2)
gambia_ever_inf <- data.frame(age = gambia_ever_inf$x, ever_inf = gambia_ever_inf$y)

# Calculate the number of susceptibles, acutely infected, chronic infected and recovered
gambia_immune <- gambia_ever_inf$ever_inf - gambia_prev$prev
gambia_infected <- gambia_prev$prev
gambia_sus <- 1-gambia_ever_inf$ever_inf

# Fill HBeAg prevalence (in HBsAg-positives) in with data from Shimakawa paper
gambia_eag <- rep(c(0.95, 0.95, 0.9, 0.9, 0.65, 0.65, 0.6, 0.6, 0.6, 0.6,
                    0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.2, 0.2,
                    0.1, 0.1, rep(0.1, 20), rep (0.05, 60)), each = 1/da) # made up last value

### THE MODEL FUNCTION ----

# NOTE: I have changed the infant vaccine application from what it was in the calibrated model
imperial_model <- function(timestep, pop, parameters, sim_starttime) {

  with(as.list(parameters), {

    # PREPARATION

    # Set up population array with infection compartments
    # Matrix 1 = females, matrix 2 = males, rows = agesteps, columns = infection compartments
    # Example: pop[ages,infectionstatus,sex]
    pop <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    # Demography: define time-varying parameters (calls the row corresponding to current timestep)

    # This approach for varying timesteps calls the rate for the closest timestep:
    # Births, migration and mortality can be switched off through parameter list
    if (births_on == 0) {
      fertility_rate <- rep(0, times = ncol(fert_rates[,-1]))
    } else {
      fertility_rate <- fert_rates[min(which(abs(fert_rates[,1]-(timestep+sim_starttime))==
                                               min(abs(fert_rates[,1]-(timestep+sim_starttime))))),-1]
    }

    if (migration_on == 0) {
      migration_rate <- matrix(rep(0, times = 2*n_agecat), ncol = 2)
    } else {
      migration_rate <- matrix(c(migration_rates_female[min(which(abs(migration_rates_female[,1]-(timestep+sim_starttime))==
                                                                    min(abs(migration_rates_female[,1]-(timestep+sim_starttime))))),-1],
                                 migration_rates_male[min(which(abs(migration_rates_male[,1]-(timestep+sim_starttime))==
                                                                  min(abs(migration_rates_male[,1]-(timestep+sim_starttime))))),-1]),
                               ncol = 2)    # 2 columns for sex-specific rates
    }

    if (mortality_on == 0) {
      mortality_rate <- matrix(rep(0, times = 2*n_agecat), ncol = 2)
    } else {
      mortality_rate <- matrix(c(mort_rates_female[min(which(abs(mort_rates_female[,1]-(timestep+sim_starttime))==
                                                               min(abs(mort_rates_female[,1]-(timestep+sim_starttime))))),-1],
                                 mort_rates_male[min(which(abs(mort_rates_male[,1]-(timestep+sim_starttime))==
                                                             min(abs(mort_rates_male[,1]-(timestep+sim_starttime))))),-1]),
                               ncol = 2)    # 2 columns for sex-specific rates
    }

    # Notation:
    # Indices for infection compartments
    # Natural history compartments
    S <- 1                           # Susceptible
    IT <- 2                           # Chronic infection: immune tolerant
    IR <- 3                           # Chronic infection: immune reactive
    IC <- 4                           # Chronic infection: inactive carrier
    ENCHB <- 5                        # Chronic infection: HBeAg-negative CHB
    CC <- 6                           # Chronic disease: compensated cirrhosis
    DCC <- 7                          # Chronic disease: decompensated cirrhosis
    HCC <- 8                          # Chronic disease: hepatocellular carcinoma
    R <- 9                           # Immune
    # Screened compartments
    S_S <- 10                        # Susceptible
    IT_S <- 11                        # Chronic infection: immune tolerant
    IR_S <- 12                        # Chronic infection: immune reactive
    IC_S <- 13                        # Chronic infection: inactive carrier
    ENCHB_S <- 14                     # Chronic infection: HBeAg-negative CHB
    CC_S <- 15                        # Chronic disease: compensated cirrhosis
    DCC_S <- 16                       # Chronic disease: decompensated cirrhosis
    HCC_S <- 17                       # Chronic disease: hepatocellular carcinoma
    R_S <- 18                        # Immune
    # Ever treated compartments
    CHB_T <- 19                       # Treated IR and ENCHB
    CC_T <- 20                        # Treated CC
    DCC_T <- 21                       # Treated DCC
    HCC_T <- 22                       # HCC developing in treated
    R_T <- 23                         # Recovered after treatment

    # Grouping of infected compartments for transmission purposes
    HBeAg_neg <- c(IC:HCC, IC_S:HCC_S)   # HBeAg-negative infected compartments
    HBeAg_pos <- c(IT,IR,IT_S, IR_S)     # HBeAg-positive infected compartments
    Treated_carriers <- CHB_T:HCC_T      # treated carrier compartments

    # Infant vaccination: set date for introduction
    # Vaccination coverage is 0 until the specified starttime and only if vaccine switch is on
    # Vaccine coverage varies over time
    if (apply_vacc == 1 & timestep >= (vacc_introtime-sim_starttime)) {
      vacc_cov = timevary_vaccine_coverage(timestep + sim_starttime)
    } else {
      vacc_cov = 0
    }

    # Birth dose vaccination: set date for introduction
    # Vaccination coverage is 0 until the specified starttime and only if vaccine switch is on
    # Need to adapt vaccine coverage to vary over time
    if (apply_bdvacc == 1 & timestep >= (bdvacc_introtime-sim_starttime)) {
      bdvacc_cov = bdvacc_cov
    } else {
      bdvacc_cov = 0
    }

    # Initialise arrays for storage of outputs
    dpop <- array(rep(0,2 * n_infectioncat * n_agecat),
                  dim=c(n_agecat,n_infectioncat,2))      # female and male population in each infection comp
    deaths <- array(rep(0,2 * n_infectioncat * n_agecat),
                    dim=c(n_agecat,n_infectioncat,2))    # female and male incident deaths in each infection comp
    migrants <- array(rep(0,2 * n_infectioncat * n_agecat),
                      dim=c(n_agecat,n_infectioncat,2))  # female and male incident migrants in each infection comp
    dcum_infections <- matrix(rep(0, 2* n_agecat),
                              ncol = 2, nrow = n_agecat)      # female and male incident infections
    dcum_chronic_infections <- matrix(rep(0, 2* n_agecat),
                                      ncol = 2, nrow = n_agecat)      # female and male incident chronic infections
    dcum_hbv_deaths <- matrix(rep(0, 2* n_agecat),
                              ncol = 2, nrow = n_agecat)  # female and male incident HBV-related deaths (CC+DCC+HCC deaths)
    dcum_dcc_deaths <- matrix(rep(0, 2* n_agecat),
                              ncol = 2, nrow = n_agecat)  # female and male DCC-related deaths
    dcum_hcc_deaths <- matrix(rep(0, 2* n_agecat),
                              ncol = 2, nrow = n_agecat)  # female and male incident HCC-related deaths
    dcum_cc_deaths <- matrix(rep(0, 2* n_agecat),
                              ncol = 2, nrow = n_agecat)  # female and male incident CC-related deaths
    dcum_eag_loss <- matrix(rep(0, 2* n_agecat),
                            ncol = 2, nrow = n_agecat)  # female and male incident HBeAg loss
    dcum_sag_loss <- matrix(rep(0, 2* n_agecat),
                            ncol = 2, nrow = n_agecat)  # female and male incident HBsAg loss
    dcum_dcc <- matrix(rep(0, 2* n_agecat),
                       ncol = 2, nrow = n_agecat)  # female and male incident DCC cases
    dcum_hcc <- matrix(rep(0, 2* n_agecat),
                       ncol = 2, nrow = n_agecat)  # female and male incident HCC cases

    # Female and male transitions to HCC
    it_to_hcc_transitions <- matrix(rep(0, 2* n_agecat),
                                    ncol = 2, nrow = n_agecat)  # from immune tolerant
    ir_to_hcc_transitions <- matrix(rep(0, 2* n_agecat),
                                    ncol = 2, nrow = n_agecat)  # from immune reactive
    ic_to_hcc_transitions <- matrix(rep(0, 2* n_agecat),
                                    ncol = 2, nrow = n_agecat)  # from inactive carrier
    enchb_to_hcc_transitions <- matrix(rep(0, 2* n_agecat),
                                       ncol = 2, nrow = n_agecat)  # from HBeAg-negative CHB
    cc_to_hcc_transitions <- matrix(rep(0, 2* n_agecat),
                                    ncol = 2, nrow = n_agecat)  # from compensated cirrhosis
    dcc_to_hcc_transitions <- matrix(rep(0, 2* n_agecat),
                                     ncol = 2, nrow = n_agecat)  # from decompensated cirrhosis

    # Female and male transition from IR (HBeAg-positive) to CC
    ir_to_cc_transitions <- matrix(rep(0, 2* n_agecat),
                                   ncol = 2, nrow = n_agecat)

    # Female and male transition from ENCHB (HBeAg-negative) to CC
    enchb_to_cc_transitions <- matrix(rep(0, 2* n_agecat),
                                      ncol = 2, nrow = n_agecat)

    # Outcomes in the screened compartments
    dcum_screened_infections <- matrix(rep(0, 2* n_agecat),
                                ncol = 2, nrow = n_agecat)      # female and male incident infections after screening
    dcum_screened_chronic_infections <- matrix(rep(0, 2* n_agecat),
                                        ncol = 2, nrow = n_agecat)      # female and male incident chronic infections after screening
    dcum_screened_hbv_deaths <- matrix(rep(0, 2* n_agecat),
                              ncol = 2, nrow = n_agecat)  # female and male incident HBV-related deaths (CC+DCC+HCC deaths) after screening
    dcum_screened_cc_deaths <- matrix(rep(0, 2* n_agecat),
                                ncol = 2, nrow = n_agecat)  # female and male incident CC-related deaths after screening
    dcum_screened_dcc_deaths <- matrix(rep(0, 2* n_agecat),
                                      ncol = 2, nrow = n_agecat)  # female and male incident DCC-related deaths after screening
    dcum_screened_hcc_deaths <- matrix(rep(0, 2* n_agecat),
                                      ncol = 2, nrow = n_agecat)  # female and male incident HCC-related deaths after screening
    dcum_screened_dcc <- matrix(rep(0, 2* n_agecat),
                       ncol = 2, nrow = n_agecat)  # female and male incident DCC cases after screening
    dcum_screened_hcc <- matrix(rep(0, 2* n_agecat),
                       ncol = 2, nrow = n_agecat)  # female and male incident HCC cases after screening

    # Outcomes in the treated compartments
    dcum_treated_hbv_deaths <- matrix(rep(0, 2* n_agecat),
                                       ncol = 2, nrow = n_agecat)  # female and male incident HBV-related deaths (CC+DCC+HCC deaths) after treatment
    # No new infections in those treated
    dcum_treated_dcc_deaths <- matrix(rep(0, 2* n_agecat),
                                      ncol = 2, nrow = n_agecat)  # female and male incident DCC-related deaths after treatment
    dcum_treated_hcc_deaths <- matrix(rep(0, 2* n_agecat),
                                      ncol = 2, nrow = n_agecat)  # female and male incident HCC-related deaths after treatment
    dcum_treated_hcc <- matrix(rep(0, 2* n_agecat),
                                ncol = 2, nrow = n_agecat)  # female and male incident HCC cases after treatment
    dcum_treated_sag_loss <- matrix(rep(0, 2* n_agecat),
                            ncol = 2, nrow = n_agecat)  # female and male incident HBsAg loss after treatment

    # TRANSMISSION

    # NEED TO ADD TREATED COMPARTMENTS TO CONTRIBUTE TO MTCT?

    # Mother-to-child transmission and births
    dcum_infected_births <- sum(fertility_rate *
                                  (rowSums((1-bdvacc_cov) * mtct_prob_e * pop[index$ages_wocba,HBeAg_pos,1]) +
                                     rowSums(bdvacc_cov * mtct_prob_ebd * pop[index$ages_wocba,HBeAg_pos,1])+
                                     rowSums((1-bdvacc_cov) * mtct_prob_s * pop[index$ages_wocba,HBeAg_neg,1])+
                                     rowSums(bdvacc_cov * mtct_prob_sbd * pop[index$ages_wocba,HBeAg_neg,1])+
                                     rowSums((1-bdvacc_cov) * mtct_prob_treat_cofactor * mtct_prob_s * pop[index$ages_wocba,Treated_carriers,1])+
                                     rowSums(bdvacc_cov * mtct_prob_treatbd * pop[index$ages_wocba,Treated_carriers,1])))

    # Previous version without birth dose:
    #dcum_infected_births <- sum(fertility_rate *
    #                              (rowSums(mtct_prob_e * pop[index$ages_wocba,HBeAg_pos,1]) +
    #                                 rowSums(mtct_prob_s * pop[index$ages_wocba,HBeAg_neg,1])+
    #                                   rowSums(mtct_prob_treat * pop[index$ages_wocba,Treated_carriers,1])))

    # infected births come from acute and chronic women of childbearing age
    dcum_chronic_births <- p_chronic_function[1] * dcum_infected_births # infected babies becoming chronic carriers
    dcum_uninfected_births <- sum(fertility_rate * pop[index$ages_wocba,index$infcat_all,1]) -
      dcum_infected_births  # uninfected births = all births - infected babies
    dcum_nonchronic_births <- sum(fertility_rate * pop[index$ages_wocba,index$infcat_all,1]) -
      dcum_chronic_births
    dcum_births <- dcum_infected_births + dcum_uninfected_births

    # Alternative force of infection definition with WAIFW matrix

    # Define WAIFW matrix
    # Age-dependent mixing between 4 age groups: 0 year olds, 0.5-5 years, 6-15 years, 16-100 years
    # Assuming no effective contact between children (0.5-5 years) and adults (>15 years)
    # Assuming no horizontal transmission from and to infants (0 years old)
    beta <- matrix(0, nrow = 4, ncol = 4)  # matrix of transmission parameters
    beta[2,2] <- b1                        # transmission among children 0.5-5 years
    beta[3,3] <- b2                        # transmission among juveniles 6-15 years
    beta[4,4] <- b3                        # transmission among adults 16-100 years
    beta[2,3] <- b2                        # transmission from juveniles to children
    beta[3,2] <- b2                        # = transmission from children to juveniles
    beta[3,4] <- b3                        # transmission from adults to juveniles
    beta[4,3] <- b3                        # = transmission from juveniles to adults

    # 0.5 year olds can get infected:
    i_1to4 <- which(ages == 0.5):which(ages == (5-da))
    i_5to14 <- which(ages == 5):which(ages == (15-da))
    i_15to100 <- which(ages == 15):which(ages == (100-da))

    # Calculate infectious proportion of the population (HBeAg-negatives+HBeAg-positives+Treated carriers) in age groups
    infectious_vector <- c(sum(pop[1,HBeAg_neg,1:2])/sum(pop[1,index$infcat_all,1:2]) +
                             (alpha * sum(pop[1,HBeAg_pos,1:2])/sum(pop[1,index$infcat_all,1:2]))+
                             (alpha2 * sum(pop[1,Treated_carriers,1:2])/sum(pop[1,index$infcat_all,1:2])), # 0 year olds
                           sum(pop[i_1to4,HBeAg_neg,1:2])/sum(pop[i_1to4,index$infcat_all,1:2]) +
                             (alpha * sum(pop[i_1to4,HBeAg_pos,1:2])/sum(pop[i_1to4,index$infcat_all,1:2]))+
                             (alpha2 * sum(pop[i_1to4,Treated_carriers,1:2])/sum(pop[i_1to4,index$infcat_all,1:2])), # 1-5 year olds
                           sum(pop[i_5to14,HBeAg_neg,1:2])/sum(pop[i_5to14,index$infcat_all,1:2]) +
                             (alpha * sum(pop[i_5to14,HBeAg_pos,1:2])/sum(pop[i_5to14,index$infcat_all,1:2]))+
                             (alpha2 * sum(pop[i_5to14,Treated_carriers,1:2])/sum(pop[i_5to14,index$infcat_all,1:2])), # 6-15 year olds
                           sum(pop[i_15to100,HBeAg_neg,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]) +
                             (alpha * sum(pop[i_15to100,HBeAg_pos,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]))+
                             (alpha2 * sum(pop[i_15to100,Treated_carriers,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]))) # 16-100 year olds

    # Previous version: without treated compartments (treated compartments are not infectious)
    #infectious_vector <- c(sum(pop[1,HBeAg_neg,1:2])/sum(pop[1,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[1,HBeAg_pos,1:2])/sum(pop[1,index$infcat_all,1:2])), # 0 year olds
    #                       sum(pop[i_1to4,HBeAg_neg,1:2])/sum(pop[i_1to4,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[i_1to4,HBeAg_pos,1:2])/sum(pop[i_1to4,index$infcat_all,1:2])), # 1-5 year olds
    #                       sum(pop[i_5to14,HBeAg_neg,1:2])/sum(pop[i_5to14,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[i_5to14,HBeAg_pos,1:2])/sum(pop[i_5to14,index$infcat_all,1:2])), # 6-15 year olds
    #                       sum(pop[i_15to100,HBeAg_neg,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]) +
    #                        (alpha * sum(pop[i_15to100,HBeAg_pos,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]))) # 16-100 year olds


    # Multiply WAIFW matrix by the age-specific proportion of infectious individuals
    # Returns a vector with force of infection for every age - 4 different values:
    # 0 in 0-year olds, different values for 0.5-5, 5-15 and 15-100 year olds
    foi_unique <- beta %*% infectious_vector
    # Repeat these values for every 1 year age group (assuming 0.5 year olds can't get horizontally infected)
    foi <- c(rep(foi_unique[1], times = length(which(ages == 0):which(ages == 0.5-da))),
             rep(foi_unique[2], times = length(i_1to4)),
             rep(foi_unique[3], times = length(i_5to14)),
             rep(foi_unique[4], times = length(i_15to100)))

    # SIMULATE PROGRESSION: DIFFERENTIAL EQUATIONS (solving for each sex separately)
    # Note: age-specific progression functions are defined within run_model function

    for (i in 1:2) {        # i = sex [1 = female, 2 = male]

      # Define transitions:

      # Demography: Incident deaths and migrants
      deaths[index$ages_all,index$infcat_all,i] <- mortality_rate[index$ages_all,i] *
        pop[index$ages_all,index$infcat_all,i]

      migrants[index$ages_all,index$infcat_all,i] <- migration_rate[index$ages_all,i] *
        pop[index$ages_all,index$infcat_all,i]
      # Applying same age-specific mortality rate to every infection compartment
      # Returns an array with indicent deaths and net migrants for every age (rows), infection state (columns) and sex (arrays)

      # Natural history transitions (unscreened)

      # Infection: Incident infections (all) and incident chronic infections
      dcum_infections[index$ages_all,i] <- foi * pop[index$ages_all,S,i]
      dcum_chronic_infections[index$ages_all,i] <-
        p_chronic_function * dcum_infections[index$ages_all,i]
      # Returns a matrix with incident infections for every age (rows) and every sex (columns)

      # Incident deaths due to HBV (from cirrhosis and HCC)
      dcum_hbv_deaths[index$ages_all,i] <-
        mu_cc * pop[index$ages_all,CC,i] +
        mu_dcc * pop[index$ages_all,DCC,i] +
        mu_hcc * pop[index$ages_all,HCC,i]
      # Returns a matrix with incident HBV deaths for every age (rows) and every sex (columns)

      # Incident deaths due to HCC only
      dcum_hcc_deaths[index$ages_all,i] <- mu_hcc * pop[index$ages_all,HCC,i]
      # Returns a matrix with incident HCC deaths for every age (rows) and every sex (columns)

      # Incident deaths due to DCC only
      dcum_dcc_deaths[index$ages_all,i] <- mu_dcc * pop[index$ages_all,DCC,i]
      # Returns a matrix with incident DCC deaths for every age (rows) and every sex (columns)

      # Incident deaths due to CC only
      dcum_cc_deaths[index$ages_all,i] <- mu_cc * pop[index$ages_all,CC,i]
      # Returns a matrix with incident DCC deaths for every age (rows) and every sex (columns)

      # Incidence of HBeAg loss: transition from IR to IC and IR to ENCHB
      dcum_eag_loss[index$ages_all,i] <-
        pr_ir_ic * eag_prog_function * pop[index$ages_all,IR,i] +
        pr_ir_enchb * pop[index$ages_all,IR,i]

      # Transition from IC to R (sAg loss)
      dcum_sag_loss[index$ages_all,i] <- sag_loss * pop[index$ages_all,IC,i]

      # Transitions to HCC
      ic_to_hcc_transitions[index$ages_all,i] <-
        cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IC,i]

      it_to_hcc_transitions[index$ages_all,i] <-
        hccr_it * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IT,i]

      ir_to_hcc_transitions[index$ages_all,i] <-
        hccr_ir * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IR,i]

      enchb_to_hcc_transitions[index$ages_all,i] <-
        hccr_enchb * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB,i]

      # CC to HCC transitions = HCC incidence in compensated cirrhotics
      cc_to_hcc_transitions[index$ages_all,i] <-
        hccr_cc * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CC,i]

      # DCC to HCC transitions = HCC incidence in decompensated cirrhotics
      dcc_to_hcc_transitions[index$ages_all,i] <-
        hccr_dcc * pop[index$ages_all,DCC,i]

      # DCC incidence
      dcum_dcc[index$ages_all,i] <- pr_cc_dcc * pop[index$ages_all,CC,i]

      # Total HCC incidence
      dcum_hcc[index$ages_all,i] <-
        it_to_hcc_transitions[index$ages_all,i] +
        ir_to_hcc_transitions[index$ages_all,i] +
        ic_to_hcc_transitions[index$ages_all,i] +
        enchb_to_hcc_transitions[index$ages_all,i] +
        cc_to_hcc_transitions[index$ages_all,i] +
        dcc_to_hcc_transitions[index$ages_all,i]

      # IR to CC transitions = cirrhosis incidence in HBeAg-positives
      ir_to_cc_transitions[index$ages_all,i] <-
        pr_ir_cc_function[index$ages_all,i] * pop[index$ages_all,IR,i]

      # ENCHB to CC transitions = cirrhosis incidence in HBeAg-negatives
      enchb_to_cc_transitions[index$ages_all,i] <-
        pr_enchb_cc_rates[index$ages_all,i] * pop[index$ages_all,ENCHB,i]

      # Post-screening transitions (untreated compartments)

      # Infection: Incident infections (all) and incident chronic infections after screening (untreated)
      dcum_screened_infections[index$ages_all,i] <- foi * pop[index$ages_all,S_S,i]
      dcum_screened_chronic_infections[index$ages_all,i] <-
        p_chronic_function * dcum_screened_infections[index$ages_all,i]
      # Returns a matrix with incident infections for every age (rows) and every sex (columns)

      # Incident deaths due to HBV (from cirrhosis and HCC) after screening (untreated)
      dcum_screened_hbv_deaths[index$ages_all,i] <-
        mu_cc * pop[index$ages_all,CC_S,i] +
        mu_dcc * pop[index$ages_all,DCC_S,i] +
        mu_hcc * pop[index$ages_all,HCC_S,i]
      # Returns a matrix with incident HBV deaths for every age (rows) and every sex (columns)

      # Incident deaths due to HCC only after screening (untreated)
      dcum_screened_hcc_deaths[index$ages_all,i] <- mu_hcc * pop[index$ages_all,HCC_S,i]
      # Returns a matrix with incident HCC deaths for every age (rows) and every sex (columns)

      # Incident deaths due to DCC only after screening (untreated)
      dcum_screened_dcc_deaths[index$ages_all,i] <- mu_dcc * pop[index$ages_all,DCC_S,i]
      # Returns a matrix with incident DCC deaths for every age (rows) and every sex (columns)

      # Incident deaths due to CC only after screening (untreated)
      dcum_screened_cc_deaths[index$ages_all,i] <- mu_cc * pop[index$ages_all,CC_S,i]
      # Returns a matrix with incident DCC deaths for every age (rows) and every sex (columns)

      # DCC incidence after screening (untreated)
      dcum_screened_dcc[index$ages_all,i] <- pr_cc_dcc * pop[index$ages_all,CC_S,i]

      # Total HCC incidence after screening (untreated)
      dcum_screened_hcc[index$ages_all,i] <-
        hccr_it * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IT_S,i] +
        hccr_ir * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IR_S,i] +
        cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IC_S,i] +
        hccr_enchb * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB_S,i] +
        hccr_cc * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CC_S,i] +
        hccr_dcc * pop[index$ages_all,DCC_S,i]

      # Post-treatment transitions

      # Incident deaths due to HBV (from cirrhosis and HCC) after treatment
      dcum_treated_hbv_deaths[index$ages_all,i] <-
        tmu_dcc * pop[index$ages_all,DCC_T,i] +
        tmu_hcc_cofactor * mu_hcc * pop[index$ages_all,HCC_T,i]
      # Returns a matrix with incident HBV deaths for every age (rows) and every sex (columns)

      # Incident deaths due to HCC only after treatment
      dcum_treated_hcc_deaths[index$ages_all,i] <- tmu_hcc_cofactor * mu_hcc * pop[index$ages_all,HCC_T,i]
      # Returns a matrix with incident HCC deaths for every age (rows) and every sex (columns)

      # Incident deaths due to DCC only after treatment
      dcum_treated_dcc_deaths[index$ages_all,i] <- tmu_dcc * pop[index$ages_all,DCC_T,i]
      # Returns a matrix with incident DCC deaths for every age (rows) and every sex (columns)

      # Total HCC incidence after treatment
      dcum_treated_hcc[index$ages_all,i] <-
        thccr_chb * hccr_enchb * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CHB_T,i] +
        thccr_cc * hccr_cc * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CC_T,i] +
        thccr_dcc * hccr_dcc * pop[index$ages_all,DCC_T,i]

      # Incident HBsAg loss after treatment
      dcum_treated_sag_loss[index$ages_all,i] <- sag_loss * pop[index$ages_all,CHB_T,i]

      # DIFFERENTIAL EQUATIONS
      # Transitions between natural history compartments (no screening):

      # Susceptibles
      dpop[index$ages_all,S,i] <- -(diff(c(0,pop[index$ages_all,S,i]))/da) -
        dcum_chronic_infections[index$ages_all,i] -
        (1-p_chronic_function) * dcum_infections[index$ages_all,i] -
        deaths[index$ages_all,S,i] +
        migrants[index$ages_all,S,i]

      # Immune tolerant
      dpop[index$ages_all,IT,i] <- -(diff(c(0,pop[index$ages_all,IT,i]))/da) +
        dcum_chronic_infections[index$ages_all,i] -
        pr_it_ir * eag_prog_function * pop[index$ages_all,IT,i] -
        it_to_hcc_transitions[index$ages_all,i] -
        deaths[index$ages_all,IT,i] + migrants[index$ages_all,IT,i]

      # Immune reactive
      dpop[index$ages_all,IR,i] <- -(diff(c(0,pop[index$ages_all,IR,i]))/da) +
        pr_it_ir * eag_prog_function * pop[index$ages_all,IT,i] -
        pr_ir_ic * eag_prog_function * pop[index$ages_all,IR,i] -
        pr_ir_enchb * pop[index$ages_all,IR,i] -
        ir_to_cc_transitions[index$ages_all,i] -
        ir_to_hcc_transitions[index$ages_all,i] -
        deaths[index$ages_all,IR,i] + migrants[index$ages_all,IR,i]

      # Inactive carrier
      dpop[index$ages_all,IC,i] <- -(diff(c(0,pop[index$ages_all,IC,i]))/da) +
        pr_ir_ic * eag_prog_function * pop[index$ages_all,IR,i] -
        pr_ic_enchb * pop[index$ages_all,IC,i] -
        dcum_sag_loss[index$ages_all,i] -
        ic_to_hcc_transitions[index$ages_all,i] -
        deaths[index$ages_all,IC,i] + migrants[index$ages_all,IC,i]

      # HBeAg-negative CHB
      dpop[index$ages_all,ENCHB,i] <- -(diff(c(0,pop[index$ages_all,ENCHB,i]))/da) +
        pr_ir_enchb * pop[index$ages_all,IR,i] +
        pr_ic_enchb * pop[index$ages_all,IC,i] -
        enchb_to_cc_transitions[index$ages_all,i] -
        enchb_to_hcc_transitions[index$ages_all,i] -
        deaths[index$ages_all,ENCHB,i] + migrants[index$ages_all,ENCHB,i]

      # Compensated cirrhosis
      dpop[index$ages_all,CC,i] <- -(diff(c(0,pop[index$ages_all,CC,i]))/da) +
        enchb_to_cc_transitions[index$ages_all,i] +
        ir_to_cc_transitions[index$ages_all,i] -
        dcum_dcc[index$ages_all,i] -
        cc_to_hcc_transitions[index$ages_all,i] -
        dcum_cc_deaths[index$ages_all,i] -
        deaths[index$ages_all,CC,i] + migrants[index$ages_all,CC,i]

      # Decompensated cirrhosis
      dpop[index$ages_all,DCC,i] <- -(diff(c(0,pop[index$ages_all,DCC,i]))/da) +
        dcum_dcc[index$ages_all,i] -
        dcc_to_hcc_transitions[index$ages_all,i] -
        dcum_dcc_deaths[index$ages_all,i] -
        deaths[index$ages_all,DCC,i] + migrants[index$ages_all,DCC,i]

      # HCC
      dpop[index$ages_all,HCC,i] <- -(diff(c(0,pop[index$ages_all,HCC,i]))/da) +
        dcum_hcc[index$ages_all,i] -      # includes transitions from IT, IR, IC, ENCHB, CC + DCC
        dcum_hcc_deaths[index$ages_all,i] -
        deaths[index$ages_all,HCC,i] + migrants[index$ages_all,HCC,i]

      # Immunes
      dpop[index$ages_all,R,i] <- -(diff(c(0,pop[index$ages_all,R,i]))/da) +
        (1-p_chronic_function) * dcum_infections[index$ages_all,i] +
        dcum_sag_loss[index$ages_all,i] -
        deaths[index$ages_all,R,i] + migrants[index$ages_all,R,i]

      # Babies are born susceptible or infected (age group 1 = age 0)
      dpop[1,S,i] <- dpop[1,S,i] + sex_ratio[i] * dcum_nonchronic_births
      dpop[1,IT,i] <- dpop[1,IT,i] + sex_ratio[i] * dcum_chronic_births

      # New vaccination approach:
       dpop[which(ages == 0):which(ages == 1-da),S,i] <- dpop[which(ages == 0):which(ages == 1-da),S,i] -
        (-log(1-(vacc_cov * vacc_eff))/0.5 * pop[which(ages == 0):which(ages == 1-da),S,i])
       dpop[which(ages == 0):which(ages == 1-da),R,i] <- dpop[which(ages == 0):which(ages == 1-da),R,i] +
        (-log(1-(vacc_cov * vacc_eff))/0.5 * pop[which(ages == 0):which(ages == 1-da),S,i])

       # Previous vaccination approach: applied at 0.5 years of age
      # dpop[which(ages == 0.5),S,i] <- dpop[which(ages == 0.5),S,i] -
      #  (vacc_cov * vacc_eff * pop[which(ages == 0.5),S,i])
      # dpop[which(ages == 0.5),R,i] <- dpop[which(ages == 0.5),R,i] +
      #  (vacc_cov * vacc_eff * pop[which(ages == 0.5),S,i])

      # Transitions between screened compartments (same as natural history but without births)

      # Screened Susceptibles
      dpop[index$ages_all,S_S,i] <- -(diff(c(0,pop[index$ages_all,S_S,i]))/da) -
        dcum_screened_chronic_infections[index$ages_all,i] -
        (1-p_chronic_function) * dcum_screened_infections[index$ages_all,i] -
        deaths[index$ages_all,S_S,i] +
        migrants[index$ages_all,S_S,i]

      # Screened Immune tolerant
      dpop[index$ages_all,IT_S,i] <- -(diff(c(0,pop[index$ages_all,IT_S,i]))/da) +
        dcum_screened_chronic_infections[index$ages_all,i] -
        pr_it_ir * eag_prog_function * pop[index$ages_all,IT_S,i] -
        hccr_it * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IT_S,i] -
        deaths[index$ages_all,IT_S,i] + migrants[index$ages_all,IT_S,i]

      # Screened Immune reactive
      dpop[index$ages_all,IR_S,i] <- -(diff(c(0,pop[index$ages_all,IR_S,i]))/da) +
        pr_it_ir * eag_prog_function * pop[index$ages_all,IT_S,i] -
        pr_ir_ic * eag_prog_function * pop[index$ages_all,IR_S,i] -
        pr_ir_enchb * pop[index$ages_all,IR_S,i] -
        pr_ir_cc_function[index$ages_all,i] * pop[index$ages_all,IR_S,i] -
        hccr_ir * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IR_S,i] -
        monitoring_rate * monitoring_prob * treatment_initiation_prob * pop[index$ages_all,IR_S,i] -       # transition into treated CHB
        deaths[index$ages_all,IR_S,i] + migrants[index$ages_all,IR_S,i]

      # Screened Inactive carrier
      dpop[index$ages_all,IC_S,i] <- -(diff(c(0,pop[index$ages_all,IC_S,i]))/da) +
        pr_ir_ic * eag_prog_function * pop[index$ages_all,IR_S,i] -
        pr_ic_enchb * pop[index$ages_all,IC_S,i] -
        sag_loss * pop[index$ages_all,IC_S,i] -
        cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IC_S,i] -
        deaths[index$ages_all,IC_S,i] + migrants[index$ages_all,IC_S,i]

      # Screened HBeAg-negative CHB
      dpop[index$ages_all,ENCHB_S,i] <- -(diff(c(0,pop[index$ages_all,ENCHB_S,i]))/da) +
        pr_ir_enchb * pop[index$ages_all,IR_S,i] +
        pr_ic_enchb * pop[index$ages_all,IC_S,i] -
        pr_enchb_cc_rates[index$ages_all,i] * pop[index$ages_all,ENCHB_S,i] -
        hccr_enchb * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB_S,i] -
        monitoring_rate * monitoring_prob * treatment_initiation_prob * pop[index$ages_all,ENCHB_S,i] -     # transition into treated CHB
        deaths[index$ages_all,ENCHB_S,i] + migrants[index$ages_all,ENCHB_S,i]

      # Screened Compensated cirrhosis
      dpop[index$ages_all,CC_S,i] <- -(diff(c(0,pop[index$ages_all,CC_S,i]))/da) +
        pr_enchb_cc_rates[index$ages_all,i] * pop[index$ages_all,ENCHB_S,i] +
        pr_ir_cc_function[index$ages_all,i] * pop[index$ages_all,IR_S,i] -
        dcum_screened_dcc[index$ages_all,i] -
        hccr_cc * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CC_S,i] -
        dcum_screened_cc_deaths[index$ages_all,i] -
        monitoring_rate * monitoring_prob * treatment_initiation_prob * pop[index$ages_all,CC_S,i] -    # transition into treated CC
        deaths[index$ages_all,CC_S,i] + migrants[index$ages_all,CC_S,i]

      # Screened Decompensated cirrhosis
      dpop[index$ages_all,DCC_S,i] <- -(diff(c(0,pop[index$ages_all,DCC_S,i]))/da) +
        dcum_screened_dcc[index$ages_all,i] -
        hccr_dcc * pop[index$ages_all,DCC_S,i] -
        dcum_screened_dcc_deaths[index$ages_all,i] -
        monitoring_rate * monitoring_prob * treatment_initiation_prob  * pop[index$ages_all,DCC_S,i] -   # transition into treated DCC
        deaths[index$ages_all,DCC_S,i] + migrants[index$ages_all,DCC_S,i]

      # Screened HCC
      dpop[index$ages_all,HCC_S,i] <- -(diff(c(0,pop[index$ages_all,HCC_S,i]))/da) +
        dcum_screened_hcc[index$ages_all,i] -
        dcum_screened_hcc_deaths[index$ages_all,i] -
        deaths[index$ages_all,HCC_S,i] + migrants[index$ages_all,HCC_S,i]

      # Screened Immunes
      dpop[index$ages_all,R_S,i] <- -(diff(c(0,pop[index$ages_all,R_S,i]))/da) +
        (1-p_chronic_function) * dcum_screened_infections[index$ages_all,i] +
        sag_loss * pop[index$ages_all,IC_S,i] -
        deaths[index$ages_all,R_S,i] + migrants[index$ages_all,R_S,i]

      # Transitions between (ever) treated compartments

      # Ever treated CHB
      dpop[index$ages_all,CHB_T,i] <- -(diff(c(0,pop[index$ages_all,CHB_T,i]))/da) +
        monitoring_rate * monitoring_prob * treatment_initiation_prob  * pop[index$ages_all,IR_S,i] +
        monitoring_rate * monitoring_prob * treatment_initiation_prob * pop[index$ages_all,ENCHB_S,i] -
        sag_loss * pop[index$ages_all,CHB_T,i] -
        thccr_chb * hccr_enchb * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CHB_T,i] -
        deaths[index$ages_all,CHB_T,i] + migrants[index$ages_all,CHB_T,i]

      # Ever treated CC
      dpop[index$ages_all,CC_T,i] <- -(diff(c(0,pop[index$ages_all,CC_T,i]))/da) +
        monitoring_rate * monitoring_prob * treatment_initiation_prob * pop[index$ages_all,CC_S,i] -
        thccr_cc * hccr_cc * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CC_T,i] -
        deaths[index$ages_all,CC_T,i] + migrants[index$ages_all,CC_T,i]

      # Ever treated DCC
      dpop[index$ages_all,DCC_T,i] <- -(diff(c(0,pop[index$ages_all,DCC_T,i]))/da) +
        monitoring_rate * monitoring_prob * treatment_initiation_prob * pop[index$ages_all,DCC_S,i] -
        thccr_dcc * hccr_dcc * pop[index$ages_all,DCC_T,i] -
        dcum_treated_dcc_deaths[index$ages_all,i]-
        deaths[index$ages_all,DCC_T,i] + migrants[index$ages_all,DCC_T,i]

      # HCC after treatment
      dpop[index$ages_all,HCC_T,i] <- -(diff(c(0,pop[index$ages_all,HCC_T,i]))/da) +
        dcum_treated_hcc[index$ages_all,i] -
        dcum_treated_hcc_deaths[index$ages_all,i] -
        deaths[index$ages_all,HCC_T,i] + migrants[index$ages_all,HCC_T,i]

     # Recovery after treatment
      dpop[index$ages_all,R_T,i] <- -(diff(c(0,pop[index$ages_all,R_T,i]))/da) +
        dcum_treated_sag_loss[index$ages_all,i] -
        deaths[index$ages_all,R_T,i] + migrants[index$ages_all,R_T,i]

    }

    # OUTPUT


    # Sum age-specific number of incident background deaths across infection compartments for output
    dcum_deaths <- cbind(rowSums(deaths[,,1]), rowSums(deaths[,,2]))
    # Age-specific number of incident background deaths among liver disease patients (CC, DCC and HCC)
    dcum_background_deaths_ld <- cbind(rowSums(deaths[index$ages_all,CC:HCC,1]),
                                       rowSums(deaths[index$ages_all,CC:HCC,2]))

    # Return results
    res <- c(dpop,
             dcum_deaths,
             dcum_infections,
             dcum_chronic_infections,
             dcum_births,
             dcum_infected_births,
             dcum_chronic_births,
             dcum_hbv_deaths,
             dcum_hcc_deaths,
             dcum_eag_loss,
             dcum_sag_loss,
             dcum_dcc,
             dcum_hcc,
             dcc_to_hcc_transitions,
             dcum_background_deaths_ld,
             dcum_dcc_deaths,
             dcum_cc_deaths,
             ir_to_cc_transitions,
             enchb_to_cc_transitions,
             cc_to_hcc_transitions,
             it_to_hcc_transitions,
             ir_to_hcc_transitions,
             ic_to_hcc_transitions,
             enchb_to_hcc_transitions,
             dcum_screened_infections,
             dcum_screened_chronic_infections,
             dcum_screened_hbv_deaths,
             dcum_screened_cc_deaths,
             dcum_screened_dcc_deaths,
             dcum_screened_hcc_deaths,
             dcum_screened_dcc,
             dcum_screened_hcc,
             dcum_treated_hbv_deaths,
             dcum_treated_dcc_deaths,
             dcum_treated_hcc_deaths,
             dcum_treated_hcc,
             dcum_treated_sag_loss,
             total_screened_uninfected = 0,
             total_screened_chb = 0,
             total_screened_cirrhosis = 0,
             total_screened_ineligible = 0)

    list(res, p_chronic_function = p_chronic_function)

  })

}

### Model-related functions ----

## Event function: reset population size to initial (1850) size in 1950
reset_pop_1950 <- function(timestep, pop, parameters){
  with (as.list(pop, parameters),{

    pop_to_reset <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))
    initialpop <- array(unlist(init_pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    current_pop <- cbind(rowSums(pop_to_reset[,,1]), rowSums(pop_to_reset[,,2]))
    pop_increase <- cbind(rowSums(initialpop[,,1]), rowSums(initialpop[,,2]))/current_pop
    scaler <- array(c(rep(pop_increase[,1], n_infectioncat),
                      rep(pop_increase[,2], n_infectioncat)), dim = c(n_agecat,n_infectioncat,2))
    pop_to_reset <- scaler * pop_to_reset
    return(c(pop_to_reset, unlist(init_pop[(2*n_infectioncat*n_agecat+1):length(init_pop)])))

  })
}

screen_pop <- function(timestep, pop, parameters){
  with (as.list(pop),{

    # Define indices
    #nat_hist_f <- 1:(n_nathistcat*n_agecat)
    #screened_f <- (max(nat_hist_f)+1):((n_nathistcat+n_screencat)*n_agecat)
    #treated_f <- (max(screened_f)+1):(n_infectioncat*n_agecat)
    #nat_hist_m <- (max(treated_f)+1):((n_infectioncat+n_nathistcat)*n_agecat)
    #screened_m <- (max(nat_hist_m)+1):((n_infectioncat+n_nathistcat+n_screencat)*n_agecat)
    #treated_m <- (max(screened_m)+1):(2*n_infectioncat*n_agecat)

    # Select the compartments to screen: infection compartments 1-9 (S-R), women and men
    #comps_to_screen <- array(unlist(pop[c(nat_hist_f, nat_hist_m)]),
    #                         dim=c(n_agecat,n_nathistcat,2))
    # Select the screened compartments: infection compartments 10-18 (S_S-R_S)
    #comps_screened <- array(unlist(pop[c(screened_f, screened_m)]),
    #                      dim=c(n_agecat,n_screencat,2))
    #comps_treated <- array(unlist(pop[c(treated_f, treated_m)]),
    #                       dim=c(n_agecat,n_treatcat,2))

    # Define indices for age groupts to screen
    age_groups_to_screen <- which(ages==parameters$min_age_to_screen):which(ages==parameters$max_age_to_screen)

    # Select state variables array
    total_pop <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    # Record the total number of people screened at each screening event
    # Note this is not a cumulative output in the model and stays the same for years where there
    # is no screening programme (sum unique values to calculate total HBsAg tests)
    pop_to_screen_uninfected <- parameters$screening_coverage * total_pop[age_groups_to_screen,c(1,9),1:2]
    pop_to_screen_chb <- parameters$screening_coverage * total_pop[age_groups_to_screen,c(3,5),1:2]
    pop_to_screen_cirrhosis <- parameters$screening_coverage * total_pop[age_groups_to_screen,c(6,7),1:2]
    pop_to_screen_ineligible <- parameters$screening_coverage * total_pop[age_groups_to_screen,c(2,4,8),1:2]

    total_screened_uninfected <- sum(pop_to_screen_uninfected)
    total_screened_chb <- sum(pop_to_screen_chb)
    total_screened_cirrhosis <- sum(pop_to_screen_cirrhosis)
    total_screened_ineligible <- sum(pop_to_screen_ineligible)

    # Calculate the population to screen then treat in the treatment eligible compartments:
    # IR (3), ENCHB (5), CC (6), DCC (7)
    pop_to_treat_ir <- parameters$screening_coverage * parameters$link_to_care_prob *
      parameters$treatment_initiation_prob * total_pop[age_groups_to_screen,3,1:2]

    pop_to_treat_enchb <- parameters$screening_coverage * parameters$link_to_care_prob *
      parameters$treatment_initiation_prob * total_pop[age_groups_to_screen,5,1:2]

    pop_to_treat_cc <- parameters$screening_coverage * parameters$link_to_care_prob *
      parameters$treatment_initiation_prob * total_pop[age_groups_to_screen,6,1:2]

    pop_to_treat_dcc <- parameters$screening_coverage * parameters$link_to_care_prob *
      parameters$treatment_initiation_prob * total_pop[age_groups_to_screen,7,1:2]

    # Calculate the population identified as treatment-ineligible (eligible for monitoring):
    # IT (2), IC (4)
    pop_to_monitor <- parameters$screening_coverage * parameters$link_to_care_prob *
      total_pop[age_groups_to_screen,c(2,4),1:2]

    # Explore treating IT>30 years
    # Currently susceptibles are screened but remain in susceptible compartment: explore vaccinating

    # Move the screened population
#    total_pop[age_groups_to_screen,1:n_nathistcat,1:2] <- total_pop[age_groups_to_screen,1:n_nathistcat,1:2]-pop_to_screen
#    total_pop[age_groups_to_screen,(n_nathistcat+1):(n_nathistcat+n_screencat),1:2] <-
#    total_pop[age_groups_to_screen,(n_nathistcat+1):(n_nathistcat+n_screencat),1:2]+
#    pop_to_screen

    # Remove the population to treat from undiagnosed compartments
    total_pop[age_groups_to_screen,3,1:2] <- total_pop[age_groups_to_screen,3,1:2] -
      pop_to_treat_ir
    total_pop[age_groups_to_screen,5,1:2] <- total_pop[age_groups_to_screen,5,1:2] -
      pop_to_treat_enchb
    total_pop[age_groups_to_screen,6,1:2] <- total_pop[age_groups_to_screen,6,1:2] -
      pop_to_treat_cc
    total_pop[age_groups_to_screen,7,1:2] <- total_pop[age_groups_to_screen,7,1:2] -
      pop_to_treat_dcc

    # Add the population to treat to treated compartments
    # CHB_T (19), CC_T (20), DCC_T (21)
    total_pop[age_groups_to_screen,19,1:2] <- total_pop[age_groups_to_screen,19,1:2] + pop_to_treat_ir +
      pop_to_treat_enchb
    total_pop[age_groups_to_screen,20,1:2] <- total_pop[age_groups_to_screen,20,1:2] + pop_to_treat_cc
    total_pop[age_groups_to_screen,21,1:2] <- total_pop[age_groups_to_screen,21,1:2] + pop_to_treat_dcc

    # Move the population from undiagnosed to treatment-ineligible (to monitor) compartments
    # IT_S (11) and IC_S (13)
    total_pop[age_groups_to_screen,c(2,4),1:2] <- total_pop[age_groups_to_screen,c(2,4),1:2] - pop_to_monitor
    total_pop[age_groups_to_screen,c(11,13),1:2] <- total_pop[age_groups_to_screen,c(11,13),1:2] + pop_to_monitor

    return(c(total_pop, unlist(init_pop[(2*n_infectioncat*n_agecat+1):(length(init_pop)-4)]),
             total_screened_uninfected, total_screened_chb,
             total_screened_cirrhosis, total_screened_ineligible))
  })
}

event_func <- function(timestep, pop, parameters){

#    if(timestep >= 1950-parameters$sim_starttime & timestep <= 1950.5-parameters$sim_starttime) {
  if(timestep == 1950-parameters$sim_starttime) {
    reset_pop_1950(timestep, pop, parameters)
  } else if(timestep %in% (parameters$screening_years-parameters$sim_starttime)) {
    screen_pop(timestep, pop, parameters)
  } else {
    return(c(unlist(pop)))
  }

}

## Functions to run the model

# Function to define/update parameter values for model run
generate_parameters <- function(..., default_parameter_list, parms_to_change = list(...)) {

  # Default situation: using input parameters
  defaults <- default_parameter_list

  if (length(parms_to_change) == 0L) {
    print("Using default parameter values")
    return(defaults)
  }

  # Alternative: parameter values to change are specified in function call
  # Check they are named correctly
  if (is.null(names(parms_to_change)) || !all(nzchar(names(parms_to_change)))) {
    stop("All arguments must be named")
  }

  # Check if specified parameters exist in default list
  new_parms_added <- setdiff(names(parms_to_change), names(defaults))

  # If not in list, print warning message for newly added parameters
  if (length(new_parms_added) > 0L) {
    print(paste0("WARNING! New parameters added: ", new_parms_added))
  } else {
    print("Using parameter values specified in function call")
  }
  # NOTE: could change this to check if the value has been updated

  # Final parameter set to use in model run: updated default parameter list
  final_parms <- modifyList(defaults, parms_to_change)

  # Turn age thresholds into integers:
  final_parms$cancer_age_threshold <- round(final_parms$cancer_age_threshold,0)
  final_parms$pr_ir_cc_age_threshold <- round(final_parms$pr_ir_cc_age_threshold,0)

  return(final_parms)

}

# Function to run the model once for a given scenario
# This includes calculation of age-specific progression functions
run_model <- function(..., sim_duration = runtime,
                      init_pop_vector = init_pop,
                      default_parameter_list, parms_to_change = list(...),
                      scenario = "vacc") {

  ## Define parameter values for model run:
  # Using default input parameter list or with updated values specified in parms_to_change
  parameters <- generate_parameters(default_parameter_list = default_parameter_list,
                                    parms_to_change = parms_to_change)

  input_parms <- parameters

  ## DEFINE FUNCTIONS FOR AGE-SPECIFIC NATURAL HISTORY PROGESSION
  # Using given parameters, and save in parameter list for input

  # Age-specific risk of becoming a chronic carrier after infection (Edmunds function)
  # The risk is fixed at 0.89 for under 0.5 year olds (perinatal infection)
  parameters$p_chronic_function <- c(rep(parameters$p_chronic_in_mtct,0.5/da),
                                     exp(-parameters$p_chronic_function_r *
                                           ages[which(ages == 0.5):n_agecat]^parameters$p_chronic_function_s))

  # Age-specific function of progression through IT and IR (IT=>IR and IR=>IC)
  # Represented by an exponential function that decreases with age
  parameters$eag_prog_function <- exp(parameters$eag_prog_function_rate * ages)

  # Age-specific progression to HCC from all carrier compartments other than DCC
  # Represented by a shifted quadratic function that increases with age and
  # prevents people younger than 10 years to progress to HCC
  # Margaret's version:
  #  cancer_prog_function <- (parameters$cancer_prog_coefficient * (ages - parameters$cancer_age_threshold))^2  # Rate in females
  #  cancer_prog_function <- cancer_prog_function *
  #    c(rep(0, times = which(ages == parameters$cancer_age_threshold-da)),
  #      rep(1, times = n_agecat - which(ages == parameters$cancer_age_threshold-da)))  # Set transition to 0 in <10 year olds
  #  cancer_prog_female <- sapply(cancer_prog_function, function(x) min(x,1)) # Set maximum annual rate is 1
  #  cancer_prog_male <- sapply(parameters$cancer_prog_male_cofactor*cancer_prog_female, function(x) min(x,1))  # Rate in males, cannot exceed 1
  #  cancer_prog_rates <- matrix(data = c(cancer_prog_female, cancer_prog_male),
  #                              nrow = n_agecat, ncol = 2)  # store in a matrix to apply to compartment
  #  parameters$cancer_prog_rates <- cancer_prog_rates

  # ADAPTATION 18/06/19: add an intercept to allow switching off of age dependence
  cancer_prog_function <- (parameters$cancer_prog_coefficient_female * (ages - parameters$cancer_age_threshold))^2  # Rate in females
  cancer_prog_function <- cancer_prog_function *
    c(rep(0, times = parameters$cancer_age_threshold/da),
      rep(1, times = n_agecat - parameters$cancer_age_threshold/da))  # Set transition to 0 in <10 year olds
  cancer_prog_female <- sapply(cancer_prog_function, function(x) min(x,1)) # Set maximum annual rate is 1
  cancer_prog_male <- sapply(parameters$cancer_male_cofactor*cancer_prog_female, function(x) min(x,1))  # Rate in males, cannot exceed 1
  cancer_prog_rates <- matrix(data = c(cancer_prog_female, cancer_prog_male),
                              nrow = n_agecat, ncol = 2)  # store in a matrix to apply to compartment
  parameters$cancer_prog_rates <- cancer_prog_rates

  # Age-specific progression from IR to ENCHB: REMOVED
  # Is set to 0 in under 20 year olds and a constant in over 20 year olds
  #  pr_ir_enchb_function <- c(rep(0, times = which(ages == 0.5-da)),
  #                            rep(parameters$pr_ir_enchb, times = n_agecat - which(ages == 0.5-da)))

  # Age-specific progression from IR to CC (HBeAg-positive cirrhosis)

  # ADAPTATION 18/06/19: addition of the cirrhosis male cofactor
  pr_ir_cc_function <- c(rep(0, times = parameters$pr_ir_cc_age_threshold/da),
                         rep(parameters$pr_ir_cc_female, times = n_agecat - parameters$pr_ir_cc_age_threshold/da))
  pr_ir_cc_function_female <- sapply(pr_ir_cc_function, function(x) min(x,5))  # annual rate cannot exceed 5
  pr_ir_cc_function_male <- sapply(pr_ir_cc_function_female*parameters$cirrhosis_male_cofactor,
                                   function(x) min(x,5))
  parameters$pr_ir_cc_function <- matrix(data = c(pr_ir_cc_function_female,
                                                  pr_ir_cc_function_male),
                                         nrow = n_agecat, ncol = 2)  # store in a matrix to apply to compartment

  # Age-specific progression from ENCHB to CC: REMOVED
  # Margaret represents this using a shifted quadratic function that increases with age and
  # prevents people younger than 25 years to progress to HCC:
  #   cirrhosis_prog_function <- parameters$pr_enchb_cc * (parameters$cirrhosis_prog_coefficient * (ages - parameters$pr_enchb_cc_age_threshold))^2  # Rate in females
  #  cirrhosis_prog_function <- parameters$pr_enchb_cc *  # rate in females
  #    c(rep(0, times = parameters$pr_enchb_cc_age_threshold/da),
  #      rep(1, times = n_agecat - parameters$pr_enchb_cc_age_threshold/da))  # Set rate to 0 in those aged < age threshold

  # Sex-specific progression from ENCHB to CC (no age effect)
  # Shimakawa 2016 found no association between current age and development of significant liver fibrosis,
  # so I remove this age dependence
  pr_enchb_cc_function <- c(rep(parameters$pr_enchb_cc_female, n_agecat))  # same rate at each age
  pr_enchb_cc_function_female <- sapply(pr_enchb_cc_function, function(x) min(x,5)) # Set maximum annual rate to 5
  pr_enchb_cc_function_male <- sapply(parameters$cirrhosis_male_cofactor *
                                        pr_enchb_cc_function_female, function(x) min(x,5))  # Rate in males, cannot exceed 5
  pr_enchb_cc_rates <- matrix(data = c(pr_enchb_cc_function_female,
                                       pr_enchb_cc_function_male),
                              nrow = n_agecat, ncol = 2)
  parameters$pr_enchb_cc_rates <- pr_enchb_cc_rates

  # Age-specific HBsAg loss (addition by me)
  # ADAPTATION 26/06/19: express this as a linear function with age based on analysis of Yusuke's data
  #sag_loss <- rep(0.01, n_agecat)
  parameters$sag_loss <- parameters$sag_loss_slope * ages

  # Update parameters for intervention scenario: vaccine (= default) or no vaccine (counterfactual)
  # vacc_screen = status quo vaccine + screening programme
  # vacc_bdvacc = status quo infant vaccine + birth dose vaccine
  # vacc_bdvacc_screen = status quo infant vaccine + birth dose vaccine + screening programme
  if (scenario == "vacc" | scenario == "vacc_screen") {
    parameters$apply_vacc <- 1
  } else if (scenario == "vacc_bdvacc" | scenario == "vacc_bdvacc_screen") {
    parameters$apply_vacc <- 1
    parameters$apply_bdvacc <- 1
  } else if (scenario == "no_vacc") {
    parameters$apply_vacc <- 0
    parameters$apply_bdvacc <- 0
  } else {
    print("Not a valid scenario. Options: vacc, no_vacc, vacc_bdvacc, vacc_screen, vacc_bdvacc_screen")
  }

  ## Run model simulation
  timestep_vector <- round((0:((sim_duration-dt)/dt))*dt,2)
  timestep_labels <- timestep_vector + parameters$sim_starttime

#  if (1950 %in% timestep_labels) {
#    timestep_1950 <- timestep_vector[which(timestep_labels == 1950)]
#    out <- as.data.frame(ode.1D(y = init_pop_vector, times = timestep_vector, func = imperial_model,
#                                parms = parameters, nspec = 1, method = "lsoda",
#                                events = list(func = reset_pop_1950, time = timestep_1950)))
#  } else {
#    out <- as.data.frame(ode.1D(y = init_pop_vector, times = timestep_vector, func = imperial_model,
#                                parms = parameters, nspec = 1, method = "lsoda"))
#  }

  timestep_1950 <- timestep_vector[which(timestep_labels == 1950)]
  timesteps_for_screening <- timestep_vector[which(timestep_labels %in% parameters$screening_years)]

  if (1950 %in% timestep_labels & scenario %in% c("vacc_screen", "vacc_bdvacc_screen")) {

    out <- as.data.frame(ode.1D(y = init_pop_vector, times = timestep_vector, func = imperial_model,
                                parms = parameters, nspec = 1, method = "ode45",
                                events = list(func = event_func,
                                              time = c(timestep_1950, timesteps_for_screening))))

  } else if (1950 %in% timestep_labels & !(scenario %in% c("vacc_screen", "vacc_bdvacc_screen"))) {

    out <- as.data.frame(ode.1D(y = init_pop_vector, times = timestep_vector, func = imperial_model,
                                parms = parameters, nspec = 1, method = "ode45",
                                events = list(func = reset_pop_1950, time = timestep_1950)))

  } else if (!(1950 %in% timestep_labels) & scenario %in% c("vacc_screen", "vacc_bdvacc_screen")) {

    out <- as.data.frame(ode.1D(y = init_pop_vector, times = timestep_vector, func = imperial_model,
                                parms = parameters, nspec = 1, method = "ode45",
                                events = list(func = screen_pop, time = timesteps_for_screening)))

  } else if (!(1950 %in% timestep_labels) & !(scenario %in% c("vacc_screen", "vacc_bdvacc_screen"))) {

    out <- as.data.frame(ode.1D(y = init_pop_vector, times = timestep_vector, func = imperial_model,
                                parms = parameters, nspec = 1, method = "ode45"))

  }

  # Add year label to timestep
  out$time   <-  out$time + parameters$sim_starttime

  return(list(out=out, input_parameters = input_parms))

  #list(func = positive_fun, time = times)
  #events = list(func = reset_pop_1950, time = 100)

}

# Function to the model twice under different scenarios
run_scenarios <- function(..., default_parameter_list, parms_to_change = list(...)) {

  sim_vacc <- run_model(sim_duration = runtime,
                     default_parameter_list = default_parameter_list,
                     parms_to_change = parms_to_change,
                     scenario = "vacc")
  out_vacc <- code_model_output(sim_vacc)

  sim_no_vacc <- run_model(sim_duration = runtime,
                           default_parameter_list = default_parameter_list,
                           parms_to_change = parms_to_change,
                           scenario = "no_vacc")
  out_no_vacc <- code_model_output(sim_no_vacc)

  outlist <- list("scenario_vacc" = out_vacc,
                  "scenario_no_vacc" = out_no_vacc)

  return(outlist)
}


### Output-related functions ----
# Function to sum numbers from different compartments for each age and time step
sum_pop_by_age <- function(time, pop_output_file) {
  pop_output <- data.frame(time = time, pop_output_file) %>%
    gather(key = "agegroup", value = "pop", -time) %>%       # turn into wide format
    arrange(time) %>%                                        # order by timestep
    mutate(agegroup = as.numeric(replace(agegroup,           # remove infection information
                                         values = ages))) %>%
    group_by(time, agegroup) %>%
    summarise(pop = sum(pop)) %>%                            # sum numbers for each age group at each timestep
    spread(key = "agegroup", value = "pop")                  # return to wide format
  # Note: this code is quicker than using do.call(cbind, by(t(output), rep(0:99,9), FUN = colSums))
  # slowest step is ordering by time

  return(as.data.frame(pop_output[,-1]))
}

# Function to calculate incidence per timestep from cumulative number output
calculate_incident_numbers <- function(cumulative_output) {
  # Takes as input cumulative number (transition) output from the model

  # First check if the input is a vector
  if (is.null(dim(cumulative_output)) == TRUE) {

    incident_numbers <- c(cumulative_output[1],  # number at first timestep
                          diff(cumulative_output, lag = 1))
    # number at current timestep - number at previous timestep
    # Replace negative numbers by the cumulative value at that timestep (since events reset count to 0)
    incident_numbers[incident_numbers<0] <- cumulative_output[which(incident_numbers<0)]

  } else {  # if not use operation on whole data frame

    incident_numbers <- rbind(cumulative_output[1,],  # number at first timestep
                              apply(cumulative_output, 2, diff, lag = 1))
    # number at current timestep - number at previous timestep
    # Replace negative numbers by the cumulative value at that timestep (since events reset count to 0)
    incident_numbers[incident_numbers<0] <- cumulative_output[incident_numbers<0]

  }

  # Returns a vector/dataframe of new cases since the last timestep
  return(incident_numbers)
}


# Function to code relevant model output (stored in list)
code_model_output <- function(output) {

  input_parms <- output$input_parameters
  output <- output$out

  ## Extract separate outputs: state variables (number at every timestep)
  out <- output[,2:(n_agecat*n_infectioncat*2+1)]

  # Infection compartments
  out_sf <- out[,grepl("^Sf.",names(out))]
  out_sm <- out[,grepl("^Sm.",names(out))]
  out_itf <- out[,grepl("^ITf.",names(out))]
  out_itm <- out[,grepl("^ITm.",names(out))]
  out_irf <- out[,grepl("^IRf.",names(out))]
  out_irm <- out[,grepl("^IRm.",names(out))]
  out_icf <- out[,grepl("^ICf.",names(out))]
  out_icm <- out[,grepl("^ICm.",names(out))]
  out_enchbf <- out[,grepl("^ENCHBf.",names(out))]
  out_enchbm <- out[,grepl("^ENCHBm.",names(out))]
  out_ccf <- out[,grepl("^CCf.",names(out))]
  out_ccm <- out[,grepl("^CCm.",names(out))]
  out_dccf <- out[,grepl("^DCCf.",names(out))]
  out_dccm <- out[,grepl("^DCCm.",names(out))]
  out_hccf <- out[,grepl("^HCCf.",names(out))]
  out_hccm <- out[,grepl("^HCCm.",names(out))]
  out_rf <- out[,grepl("^Rf.",names(out))]
  out_rm <- out[,grepl("^Rm.",names(out))]
  # Screened compartments
  out_screen_sf <- out[,grepl("^S_Sf.",names(out))]
  out_screen_sm <- out[,grepl("^S_Sm.",names(out))]
  out_screen_itf <- out[,grepl("^S_ITf.",names(out))]
  out_screen_itm <- out[,grepl("^S_ITm.",names(out))]
  out_screen_irf <- out[,grepl("^S_IRf.",names(out))]
  out_screen_irm <- out[,grepl("^S_IRm.",names(out))]
  out_screen_icf <- out[,grepl("^S_ICf.",names(out))]
  out_screen_icm <- out[,grepl("^S_ICm.",names(out))]
  out_screen_enchbf <- out[,grepl("^S_ENCHBf.",names(out))]
  out_screen_enchbm <- out[,grepl("^S_ENCHBm.",names(out))]
  out_screen_ccf <- out[,grepl("^S_CCf.",names(out))]
  out_screen_ccm <- out[,grepl("^S_CCm.",names(out))]
  out_screen_dccf <- out[,grepl("^S_DCCf.",names(out))]
  out_screen_dccm <- out[,grepl("^S_DCCm.",names(out))]
  out_screen_hccf <- out[,grepl("^S_HCCf.",names(out))]
  out_screen_hccm <- out[,grepl("^S_HCCm.",names(out))]
  out_screen_rf <- out[,grepl("^S_Rf.",names(out))]
  out_screen_rm <- out[,grepl("^S_Rm.",names(out))]
  # Treated compartments
  out_treat_chbf <- out[,grepl("^T_CHBf.",names(out))]
  out_treat_chbm <- out[,grepl("^T_CHBm.",names(out))]
  out_treat_ccf <- out[,grepl("^T_CCf.",names(out))]
  out_treat_ccm <- out[,grepl("^T_CCm.",names(out))]
  out_treat_dccf <- out[,grepl("^T_DCCf.",names(out))]
  out_treat_dccm <- out[,grepl("^T_DCCm.",names(out))]
  out_treat_hccf <- out[,grepl("^T_HCCf.",names(out))]
  out_treat_hccm <- out[,grepl("^T_HCCm.",names(out))]
  out_treat_rf <- out[,grepl("^T_Rf.",names(out))]
  out_treat_rm <- out[,grepl("^T_Rm.",names(out))]

  # Total population
  out_popf <- select(out, contains("f"))
  out_popm <- select(out, contains("m"))
  out_pop <- cbind(out_popf, out_popm)

  ## Extract separate outputs: incident variables (transitions between states)

  # Demographic transitions per timestep (cumulative number of births and deaths)
  out_cum_deathsf <- output[,grepl("^cum_deathsf.",names(output))]
  out_cum_deathsm <- output[,grepl("^cum_deathsm.",names(output))]
  out_cum_births <- unlist(select(output, contains("cum_births")))

  # Cumulative HBV incidence from horizontal transmission
  out_cum_infectionsf <- output[,grepl("^cum_infectionsf.",names(output))]
  out_cum_infectionsm <- output[,grepl("^cum_infectionsm.",names(output))]

  # Cumulative HBV incidence from MTCT (number of infected births)
  out_cum_infected_births <- unlist(output[,grepl("^cum_infected_births",names(output))])

  # Cumulative chronic infection incidence from horizontal transmission
  out_cum_chronic_infectionsf <- output[,grepl("^cum_chronic_infectionsf.",names(output))]
  out_cum_chronic_infectionsm <- output[,grepl("^cum_chronic_infectionsm.",names(output))]

  # Cumulative chronic infection incidence from MTCT (number of chronically infected births)
  out_cum_chronic_births <- unlist(output[,grepl("^cum_chronic_births",names(output))])

  # Cumulative number of HBV-related deaths (from cirrhosis and HCC)
  out_cum_hbv_deathsf <- output[,grepl("^cum_hbv_deathsf.",names(output))]
  out_cum_hbv_deathsm <- output[,grepl("^cum_hbv_deathsm.",names(output))]

  # Cumulative number of HCC cases (from all possible compartments)
  out_cum_hccf <- output[,grepl("^cum_incident_hccf.",names(output))]
  out_cum_hccm <- output[,grepl("^cum_incident_hccm.",names(output))]

  # Cumulative HBV incidence from horizontal transmission in the screened group
  out_cum_screened_infectionsf <- output[,grepl("^cum_screened_infectionsf.",names(output))]
  out_cum_screened_infectionsm <- output[,grepl("^cum_screened_infectionsm.",names(output))]

  # Cumulative chronic infection incidence from horizontal transmission in the screened group
  out_cum_screened_chronic_infectionsf <- output[,grepl("^cum_screened_chronic_infectionsf.",names(output))]
  out_cum_screened_chronic_infectionsm <- output[,grepl("^cum_screened_chronic_infectionsm.",names(output))]

  # Cumulative number of HBV-related deaths (from cirrhosis and HCC) in the screened group
  out_cum_screened_hbv_deathsf <- output[,grepl("^cum_screened_hbv_deathsf.",names(output))]
  out_cum_screened_hbv_deathsm <- output[,grepl("^cum_screened_hbv_deathsm.",names(output))]

  # Cumulative number of HCC cases (from all possible compartments) in the screened group
  out_cum_screened_hccf <- output[,grepl("^cum_screened_incident_hccf.",names(output))]
  out_cum_screened_hccm <- output[,grepl("^cum_screened_incident_hccm.",names(output))]

  # Cumulative number of HBV-related deaths (from cirrhosis and HCC) in the treated group
  out_cum_treated_hbv_deathsf <- output[,grepl("^cum_treated_hbv_deathsf.",names(output))]
  out_cum_treated_hbv_deathsm <- output[,grepl("^cum_treated_hbv_deathsm.",names(output))]

  # Cumulative number of HCC cases (from all possible compartments) in the treated group
  out_cum_treated_hccf <- output[,grepl("^cum_treated_incident_hccf.",names(output))]
  out_cum_treated_hccm <- output[,grepl("^cum_treated_incident_hccm.",names(output))]

  ## Process infection outputs
  # Combine into data frames with outputs of interest for further analysis

  # Age-specific number in each infection compartment at each time step
  sus <- data.frame(pop = out_sf + out_sm + out_screen_sf + out_screen_sm)   # need to change the column names
  carriers <- data.frame(pop = (out_itf + out_itm +
                                  out_irf + out_irm +
                                  out_icf+out_icm+
                                  out_enchbf+out_enchbm+
                                  out_ccf+out_ccm+
                                  out_dccf+out_dccm+
                                  out_hccf+out_hccm+
                                  out_screen_itf + out_screen_itm +
                                  out_screen_irf + out_screen_irm +
                                  out_screen_icf+out_screen_icm+
                                  out_screen_enchbf+out_screen_enchbm+
                                  out_screen_ccf+out_screen_ccm+
                                  out_screen_dccf+out_screen_dccm+
                                  out_screen_hccf+out_screen_hccm+
                                  out_treat_chbf+out_treat_chbm+
                                  out_treat_ccf+out_treat_ccm+
                                  out_treat_dccf+out_treat_dccm+
                                  out_treat_hccf+out_treat_hccm))
  carriers_female <- data.frame(pop = (out_itf+
                                         out_irf+
                                         out_icf+
                                         out_enchbf+
                                         out_ccf+out_dccf+out_hccf+
                                         out_screen_itf +
                                         out_screen_irf +
                                         out_screen_icf+
                                         out_screen_enchbf+
                                         out_screen_ccf+
                                         out_screen_dccf+
                                         out_screen_hccf+
                                         out_treat_chbf+
                                         out_treat_ccf+
                                         out_treat_dccf+
                                         out_treat_hccf))
  carriers_male <- data.frame(pop = (out_itm+
                                       out_irm+
                                       out_icm+
                                       out_enchbm+
                                       out_ccm+out_dccm+out_hccm+
                                       out_screen_itm +
                                       out_screen_irm +
                                       out_screen_icm+
                                       out_screen_enchbm+
                                       out_screen_ccm+
                                       out_screen_dccm+
                                       out_screen_hccm+
                                       out_treat_chbm+
                                       out_treat_ccm+
                                       out_treat_dccm+
                                       out_treat_hccm))
  immune <- data.frame(pop = out_rf + out_rm + out_screen_rf + out_screen_rm + out_treat_rf + out_treat_rm)
  ever_infected <- data.frame(pop = carriers + immune)
  ever_infected_female <- data.frame(pop = carriers_female + out_rf + out_screen_rf + out_treat_rf)
  ever_infected_male <- data.frame(pop = carriers_male + out_rm + out_screen_rm + out_treat_rm)
  eag_positive <- data.frame(pop = (out_itf + out_itm +
                                      out_irf + out_irm +
                                       out_screen_itf + out_screen_itm +
                                         out_screen_irf + out_screen_irm))
  eag_positive_female <- data.frame(pop = (out_itf + out_irf + out_screen_itf + out_screen_irf))
  eag_positive_male <- data.frame(pop = (out_itm + out_irm + out_screen_itm + out_screen_irm))

  # Total number in each infection compartment per time step
  infectioncat_total <- data.frame(time = output$time,
                                   sus = rowSums(sus),
                                   carriers = rowSums(carriers),
                                   immune = rowSums(immune),
                                   ever_infected = rowSums(ever_infected))


  # Calculate number of new cases per timestep from cumulative number output

  # Age-specific HBV incidence from horizontal transmission - for women, men and both (total)
  horizontal_infections_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_infectionsf))
  names(horizontal_infections_female) <- sprintf("incident_number%g",ages)

  horizontal_infections_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_infectionsm))
  names(horizontal_infections_male) <- sprintf("incident_number%g",ages)

  # Total number of incident infections from horizontal transmission and MTCT per time step
  incident_infections <- data.frame(time = output$time,
                                    horizontal_infections = rowSums(horizontal_infections_female) +
                                      rowSums(horizontal_infections_male),
                                    infected_births = calculate_incident_numbers(out_cum_infected_births))

  # Age-specific chronic infection incidence from horizontal transmission - for women, men and both (total)
  horizontal_chronic_infections_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_chronic_infectionsf))
  names(horizontal_chronic_infections_female) <- sprintf("incident_number%g",ages)

  horizontal_chronic_infections_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_chronic_infectionsm))
  names(horizontal_chronic_infections_male) <- sprintf("incident_number%g",ages)


  # Total number of incident chronic infections from horizontal transmission and MTCT per time step
  incident_chronic_infections <- data.frame(time = output$time,
                                            horizontal_chronic_infections = rowSums(horizontal_chronic_infections_female) +
                                              rowSums(horizontal_chronic_infections_male),
                                            chronic_births = calculate_incident_numbers(out_cum_chronic_births))

  # Age-specific number of HBV-related deaths - for women and men
  hbv_deaths_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_hbv_deathsf))
  names(hbv_deaths_female) <- sprintf("incident_number%g",ages)

  hbv_deaths_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_hbv_deathsm))
  names(hbv_deaths_male) <- sprintf("incident_number%g",ages)

  # Total number of HBV deaths per time step
  hbv_deaths <- data.frame(time = output$time,
                           incident_number_female = rowSums(hbv_deaths_female),
                           incident_number_male = rowSums(hbv_deaths_male))
  hbv_deaths$incident_number_total <- hbv_deaths$incident_number_female + hbv_deaths$incident_number_male

  # Age-specific number of total HCC cases - for women and men
  incident_hcc_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_hccf))
  names(incident_hcc_female) <- sprintf("incident_number%g",ages)

  incident_hcc_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_hccm))
  names(incident_hcc_male) <- sprintf("incident_number%g",ages)

  # Total number of total HCC cases per time step
  incident_hcc <- data.frame(time = output$time,
                             incident_number_female = rowSums(incident_hcc_female),
                             incident_number_male = rowSums(incident_hcc_male))
  incident_hcc$incident_number_total <- incident_hcc$incident_number_female +
    incident_hcc$incident_number_male

  # Screened compartents:

  # Age-specific HBV incidence from horizontal transmission after screening - for women, men and both (total)
  screened_horizontal_infections_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_screened_infectionsf))
  names(screened_horizontal_infections_female) <- sprintf("incident_number%g",ages)

  screened_horizontal_infections_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_screened_infectionsm))
  names(screened_horizontal_infections_male) <- sprintf("incident_number%g",ages)

  # Total number of incident infections after screening (only from horizontal transmission) per time step
  screened_incident_infections <- data.frame(time = output$time,
                                    screened_horizontal_infections = rowSums(screened_horizontal_infections_female) +
                                      rowSums(screened_horizontal_infections_male))

  # Age-specific chronic infection incidence from horizontal transmission after screening - for women, men and both (total)
  screened_horizontal_chronic_infections_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_screened_chronic_infectionsf))
  names(screened_horizontal_chronic_infections_female) <- sprintf("incident_number%g",ages)

  screened_horizontal_chronic_infections_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_screened_chronic_infectionsm))
  names(screened_horizontal_chronic_infections_male) <- sprintf("incident_number%g",ages)


  # Total number of incident chronic infections after screening (from horizontal transmission only) per time step
  screened_incident_chronic_infections <- data.frame(time = output$time,
                                              screened_horizontal_chronic_infections = rowSums(screened_horizontal_chronic_infections_female) +
                                              rowSums(screened_horizontal_chronic_infections_male))

  # Age-specific number of HBV-related deaths after screening - for women and men
  screened_hbv_deaths_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_screened_hbv_deathsf))
  names(screened_hbv_deaths_female) <- sprintf("incident_number%g",ages)

  screened_hbv_deaths_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_screened_hbv_deathsm))
  names(screened_hbv_deaths_male) <- sprintf("incident_number%g",ages)

  # Total number of HBV deaths per time step after screening
  screened_hbv_deaths <- data.frame(time = output$time,
                                    incident_number_female = rowSums(screened_hbv_deaths_female),
                                    incident_number_male = rowSums(screened_hbv_deaths_male))
  screened_hbv_deaths$incident_number_total <- screened_hbv_deaths$incident_number_female + screened_hbv_deaths$incident_number_male

  # Age-specific number of total HCC cases after screening - for women and men
  screened_incident_hcc_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_screened_hccf))
  names(screened_incident_hcc_female) <- sprintf("incident_number%g",ages)

  screened_incident_hcc_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_screened_hccm))
  names(screened_incident_hcc_male) <- sprintf("incident_number%g",ages)

  # Total number of total HCC cases per time step
  screened_incident_hcc <- data.frame(time = output$time,
                             incident_number_female = rowSums(screened_incident_hcc_female),
                             incident_number_male = rowSums(screened_incident_hcc_male))
  screened_incident_hcc$incident_number_total <- screened_incident_hcc$incident_number_female +
    screened_incident_hcc$incident_number_male


  # Treated compartments:

  # Age-specific number of HBV-related deaths after treatment - for women and men
  treated_hbv_deaths_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_treated_hbv_deathsf))
  names(treated_hbv_deaths_female) <- sprintf("incident_number%g",ages)

  treated_hbv_deaths_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_treated_hbv_deathsm))
  names(treated_hbv_deaths_male) <- sprintf("incident_number%g",ages)

  # Total number of HBV deaths per time step after treatment
  treated_hbv_deaths <- data.frame(time = output$time,
                                    incident_number_female = rowSums(treated_hbv_deaths_female),
                                    incident_number_male = rowSums(treated_hbv_deaths_male))
  treated_hbv_deaths$incident_number_total <- treated_hbv_deaths$incident_number_female + treated_hbv_deaths$incident_number_male

  # Age-specific number of total HCC cases after treatment - for women and men
  treated_incident_hcc_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_treated_hccf))
  names(treated_incident_hcc_female) <- sprintf("incident_number%g",ages)

  treated_incident_hcc_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_treated_hccm))
  names(treated_incident_hcc_male) <- sprintf("incident_number%g",ages)

  # Total number of total HCC cases per time step
  treated_incident_hcc <- data.frame(time = output$time,
                                      incident_number_female = rowSums(treated_incident_hcc_female),
                                      incident_number_male = rowSums(treated_incident_hcc_male))
  treated_incident_hcc$incident_number_total <- treated_incident_hcc$incident_number_female +
    treated_incident_hcc$incident_number_male


  ## Code demography outputs

  # Population:

  # Age-specific and total (last column) population per time step
  pop_female <- sum_pop_by_age(time = output$time, pop_output_file = out_popf)
  pop_male <- sum_pop_by_age(time = output$time, pop_output_file = out_popm)

  pop <- data.frame(pop = (pop_female + pop_male))

  # Total female, male and both population per time step
  pop_total <- data.frame(time = output$time,
                          pop_female = rowSums(pop_female),
                          pop_male = rowSums(pop_male)) %>%
    mutate(pop_total = pop_female + pop_male)

  # Births:

  # Total number of births at each timestep
  births <- data.frame(time = output$time,
                       incident_number = calculate_incident_numbers(out_cum_births))
  names(births) <- c("time", "incident_number")

  # Total number of births grouped in 5-year time periods
  births_group5 <- births %>%
    mutate(timegroup = floor(time / 5) * 5) %>%
    group_by(timegroup) %>%
    summarise_all(sum) %>%
    select(-time)

  # Deaths:

  # Age-specific and total deaths per time step - for women, men and both (total)
  deaths_female <- data.frame(incident_number = calculate_incident_numbers(out_cum_deathsf))
  names(deaths_female) <- sprintf("incident_number%g",ages)
  deaths_female$incident_number_total <- rowSums(deaths_female)

  deaths_male <- data.frame(incident_number = calculate_incident_numbers(out_cum_deathsm))
  names(deaths_male) <- sprintf("incident_number%g",ages)
  deaths_male$incident_number_total <- rowSums(deaths_male)

  deaths_total <- data.frame(time = output$time,
                             deaths = deaths_female + deaths_male)
  names(deaths_total)[-1] <- c(sprintf("incident_number%g",ages), "total")

  # Total number of deaths grouped in 5-year time periods
  deaths_total_group5 <- deaths_total %>%
    mutate(timegroup = floor(time / 5) * 5) %>%
    group_by(timegroup) %>%
    summarise_all(sum) %>%
    select(-time)

  toreturn <- list("time" = output$time,
                   "sus" = sus,
                   "carriers_female" = carriers_female,
                   "carriers_male" = carriers_male,
                   "carriers" = carriers,
                   "eag_positive_female" = eag_positive_female,
                   "eag_positive_male" = eag_positive_male,
                   "eag_positive" = eag_positive,
                   "immune" = immune,
                   "ever_infected" = ever_infected,
                   "ever_infected_female" = ever_infected_female,
                   "ever_infected_male" = ever_infected_male,
                   "infectioncat_total" = infectioncat_total,
                   "pop_female" = pop_female,
                   "pop_male" = pop_male,
                   "pop" = pop,
                   "pop_total" = pop_total,
                   "deaths_total_group5" = deaths_total_group5,
                   "births_group5" =  births_group5,
                   "incident_infections" = incident_infections,                   # only without screening/treatment
                   "incident_chronic_infections" = incident_chronic_infections,   # only without screening/treatment
                   "hbv_deaths" = hbv_deaths,                                     # only without screening/treatment
                   "incident_hcc" = incident_hcc,                                 # only without screening/treatment
                   "screened_incident_infections" = screened_incident_infections,
                   "screened_incident_chronic_infections" = screened_incident_chronic_infections,
                   "screened_hbv_deaths" = screened_hbv_deaths,
                   "screened_incident_hcc" = screened_incident_hcc,
                   "treated_hbv_deaths" = treated_hbv_deaths,
                   "treated_incident_hcc" = treated_incident_hcc,
                   "full_output" = output,
                   "input_parameters" = input_parms
  )
  return(toreturn)

}

# Shorter/quicker function to process model outputs:
# calculates carriers, ever infected, eAg positives and total pop
code_model_output_summary <- function(output) {

  input_parms <- output$input_parameters
  output <- output$out

  ## Extract separate outputs: state variables (number at every timestep)
  out <- output[,2:(n_agecat*n_infectioncat*2+1)]

  # Infection compartments
  out_sf <- out[,grepl("^Sf.",names(out))]
  out_sm <- out[,grepl("^Sm.",names(out))]
  out_itf <- out[,grepl("^ITf.",names(out))]
  out_itm <- out[,grepl("^ITm.",names(out))]
  out_irf <- out[,grepl("^IRf.",names(out))]
  out_irm <- out[,grepl("^IRm.",names(out))]
  out_icf <- out[,grepl("^ICf.",names(out))]
  out_icm <- out[,grepl("^ICm.",names(out))]
  out_enchbf <- out[,grepl("^ENCHBf.",names(out))]
  out_enchbm <- out[,grepl("^ENCHBm.",names(out))]
  out_ccf <- out[,grepl("^CCf.",names(out))]
  out_ccm <- out[,grepl("^CCm.",names(out))]
  out_dccf <- out[,grepl("^DCCf.",names(out))]
  out_dccm <- out[,grepl("^DCCm.",names(out))]
  out_hccf <- out[,grepl("^HCCf.",names(out))]
  out_hccm <- out[,grepl("^HCCm.",names(out))]
  out_rf <- out[,grepl("^Rf.",names(out))]
  out_rm <- out[,grepl("^Rm.",names(out))]
  # Screened compartments
  out_screen_sf <- out[,grepl("^S_Sf.",names(out))]
  out_screen_sm <- out[,grepl("^S_Sm.",names(out))]
  out_screen_itf <- out[,grepl("^S_ITf.",names(out))]
  out_screen_itm <- out[,grepl("^S_ITm.",names(out))]
  out_screen_irf <- out[,grepl("^S_IRf.",names(out))]
  out_screen_irm <- out[,grepl("^S_IRm.",names(out))]
  out_screen_icf <- out[,grepl("^S_ICf.",names(out))]
  out_screen_icm <- out[,grepl("^S_ICm.",names(out))]
  out_screen_enchbf <- out[,grepl("^S_ENCHBf.",names(out))]
  out_screen_enchbm <- out[,grepl("^S_ENCHBm.",names(out))]
  out_screen_ccf <- out[,grepl("^S_CCf.",names(out))]
  out_screen_ccm <- out[,grepl("^S_CCm.",names(out))]
  out_screen_dccf <- out[,grepl("^S_DCCf.",names(out))]
  out_screen_dccm <- out[,grepl("^S_DCCm.",names(out))]
  out_screen_hccf <- out[,grepl("^S_HCCf.",names(out))]
  out_screen_hccm <- out[,grepl("^S_HCCm.",names(out))]
  out_screen_rf <- out[,grepl("^S_Rf.",names(out))]
  out_screen_rm <- out[,grepl("^S_Rm.",names(out))]
  # Treated compartments
  out_treat_chbf <- out[,grepl("^T_CHBf.",names(out))]
  out_treat_chbm <- out[,grepl("^T_CHBm.",names(out))]
  out_treat_ccf <- out[,grepl("^T_CCf.",names(out))]
  out_treat_ccm <- out[,grepl("^T_CCm.",names(out))]
  out_treat_dccf <- out[,grepl("^T_DCCf.",names(out))]
  out_treat_dccm <- out[,grepl("^T_DCCm.",names(out))]
  out_treat_hccf <- out[,grepl("^T_HCCf.",names(out))]
  out_treat_hccm <- out[,grepl("^T_HCCm.",names(out))]
  out_treat_rf <- out[,grepl("^T_Rf.",names(out))]
  out_treat_rm <- out[,grepl("^T_Rm.",names(out))]

  # Total population
  out_popf <- select(out, contains("f"))
  out_popm <- select(out, contains("m"))
  out_pop <- cbind(out_popf, out_popm)

  ## Process infection outputs
  # Combine into data frames with outputs of interest for further analysis

  # Age-specific number in each infection compartment at each time step
  carriers <- data.frame(pop = (out_itf + out_itm +
                                  out_irf + out_irm +
                                  out_icf+out_icm+
                                  out_enchbf+out_enchbm+
                                  out_ccf+out_ccm+
                                  out_dccf+out_dccm+
                                  out_hccf+out_hccm+
                                  out_screen_itf + out_screen_itm +
                                  out_screen_irf + out_screen_irm +
                                  out_screen_icf+out_screen_icm+
                                  out_screen_enchbf+out_screen_enchbm+
                                  out_screen_ccf+out_screen_ccm+
                                  out_screen_dccf+out_screen_dccm+
                                  out_screen_hccf+out_screen_hccm+
                                  out_treat_chbf+out_treat_chbm+
                                  out_treat_ccf+out_treat_ccm+
                                  out_treat_dccf+out_treat_dccm+
                                  out_treat_hccf+out_treat_hccm))
  carriers_female <- data.frame(pop = (out_itf+
                                         out_irf+
                                         out_icf+
                                         out_enchbf+
                                         out_ccf+out_dccf+out_hccf+
                                         out_screen_itf +
                                         out_screen_irf +
                                         out_screen_icf+
                                         out_screen_enchbf+
                                         out_screen_ccf+
                                         out_screen_dccf+
                                         out_screen_hccf+
                                         out_treat_chbf+
                                         out_treat_ccf+
                                         out_treat_dccf+
                                         out_treat_hccf))
  carriers_male <- data.frame(pop = (out_itm+
                                       out_irm+
                                       out_icm+
                                       out_enchbm+
                                       out_ccm+out_dccm+out_hccm+
                                       out_screen_itm +
                                       out_screen_irm +
                                       out_screen_icm+
                                       out_screen_enchbm+
                                       out_screen_ccm+
                                       out_screen_dccm+
                                       out_screen_hccm+
                                       out_treat_chbm+
                                       out_treat_ccm+
                                       out_treat_dccm+
                                       out_treat_hccm))
  immune <- data.frame(pop = out_rf + out_rm + out_screen_rf + out_screen_rm + out_treat_rf + out_treat_rm)
  ever_infected <- data.frame(pop = carriers + immune)
  ever_infected_female <- data.frame(pop = carriers_female + out_rf + out_screen_rf + out_treat_rf)
  ever_infected_male <- data.frame(pop = carriers_male + out_rm + out_screen_rm + out_treat_rm)
  eag_positive <- data.frame(pop = (out_itf + out_itm +
                                      out_irf + out_irm +
                                      out_screen_itf + out_screen_itm +
                                      out_screen_irf + out_screen_irm))
  eag_positive_female <- data.frame(pop = (out_itf + out_irf + out_screen_itf + out_screen_irf))
  eag_positive_male <- data.frame(pop = (out_itm + out_irm + out_screen_itm + out_screen_irm))

  ## Code demography outputs

  # Population:

  # Age-specific and total (last column) population per time step
  pop_female <- sum_pop_by_age(time = output$time, pop_output_file = out_popf)
  pop_male <- sum_pop_by_age(time = output$time, pop_output_file = out_popm)

  pop <- data.frame(pop = (pop_female + pop_male))

  toreturn <- list("time" = output$time,
                   "carriers_female" = carriers_female,
                   "carriers_male" = carriers_male,
                   "carriers" = carriers,
                   "eag_positive_female" = eag_positive_female,
                   "eag_positive_male" = eag_positive_male,
                   "eag_positive" = eag_positive,
                   "ever_infected" = ever_infected,
                   "ever_infected_female" = ever_infected_female,
                   "ever_infected_male" = ever_infected_male,
                   "pop_female" = pop_female,
                   "pop_male" = pop_male,
                   "pop" = pop,
                   "input_parameters" = input_parms)
  return(toreturn)

}

# Function to extract outcomes of interest (e.g. infection incidence, HBV-related deaths) from output of
# simulations with multiple parameter sets (requires running code_model_output first)
extract_outcomes <- function(output_file, scenario_label) {

  # HBsAg prevalence

  proj_prev <- cbind(output_file[[1]]$time,
                     (sapply(lapply(output_file,"[[", "infectioncat_total"), "[[", "carriers")/
                        sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total")))
  colnames(proj_prev)[1] <- "time"
  proj_prev_summary <- data.frame(time = output_file[[1]]$time)
  proj_prev_summary$median <- apply(proj_prev[,-1],1,median)
  proj_prev_summary$lower <- apply(proj_prev[,-1],1,quantile, prob = 0.025)
  proj_prev_summary$upper <- apply(proj_prev[,-1],1,quantile, prob = 0.975)
  proj_prev_long <- gather(as.data.frame(proj_prev), key = "iteration", value = "hbsag_prev", -time)

  # Absolute number of new cases of chronic HBV carriage per timestep
  proj_inc <- cbind(output_file[[1]]$time,
                    (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                       sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births")+
                       sapply(lapply(output_file,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections")))
  colnames(proj_inc)[1] <- "time"

  proj_inc_summary <- data.frame(time = output_file[[1]]$time)
  proj_inc_summary$median <- apply(proj_inc[,-1],1,median)
  proj_inc_summary$lower <- apply(proj_inc[,-1],1,quantile, prob = 0.025)
  proj_inc_summary$upper <- apply(proj_inc[,-1],1,quantile, prob = 0.975)
  proj_inc_long <- gather(as.data.frame(proj_inc), key = "iteration", value = "chronic_cases", -time)

  # Incidence rate of chronic HBV carriage per timestep
  proj_inc_rate <- cbind(output_file[[1]]$time,
                    (sapply(lapply(output_file,"[[", "incident_chronic_infections"), "[[", "horizontal_chronic_infections")+
                       sapply(lapply(output_file, "[[", "incident_chronic_infections"), "[[", "chronic_births")+
                       sapply(lapply(output_file,"[[", "screened_incident_chronic_infections"), "[[", "screened_horizontal_chronic_infections"))/
                       sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total"))
  colnames(proj_inc_rate)[1] <- "time"
  proj_inc_rate_summary <- data.frame(time = output_file[[1]]$time)
  proj_inc_rate_summary$median <- apply(proj_inc_rate[,-1],1,median)
  proj_inc_rate_summary$lower <- apply(proj_inc_rate[,-1],1,quantile, prob = 0.025)
  proj_inc_rate_summary$upper <- apply(proj_inc_rate[,-1],1,quantile, prob = 0.975)
  proj_inc_rate_long <- gather(as.data.frame(proj_inc_rate), key = "iteration", value = "chronic_cases_rate", -time)

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

  # Incidence rate of  HBV-related deaths per timestep
  proj_deaths_rate <- cbind(output_file[[1]]$time,
                       (sapply(lapply(output_file,"[[", "hbv_deaths"), "[[", "incident_number_total")+
                          sapply(lapply(output_file,"[[", "screened_hbv_deaths"), "[[", "incident_number_total")+
                          sapply(lapply(output_file,"[[", "treated_hbv_deaths"), "[[", "incident_number_total"))/
                         sapply(lapply(output_file,"[[", "pop_total"), "[[", "pop_total"))
  colnames(proj_deaths_rate)[1] <- "time"
  proj_deaths_rate_summary <- data.frame(time = output_file[[1]]$time)
  proj_deaths_rate_summary$median <- apply(proj_deaths_rate[,-1],1,median)
  proj_deaths_rate_summary$lower <- apply(proj_deaths_rate[,-1],1,quantile, prob = 0.025)
  proj_deaths_rate_summary$upper <- apply(proj_deaths_rate[,-1],1,quantile, prob = 0.975)
  proj_deaths_rate_long <- gather(as.data.frame(proj_deaths_rate), key = "iteration", value = "deaths_rate", -time)

  # Age-standardised chronic infection incidence per 100000 per timestep!!

  # Age-standardised HBV-related deaths incidence per 100000 per timestep
  proj_deaths_standardised <- sapply(output_file, calculate_age_standardised_hbvdeaths_rate) # per person per timestep
  proj_deaths_standardised_summary <- data.frame(time = output_file[[1]]$time)
  proj_deaths_standardised_summary$median <- apply(proj_deaths_standardised,1,median)
  proj_deaths_standardised_summary$lower <- apply(proj_deaths_standardised,1,quantile, prob = 0.025)
  proj_deaths_standardised_summary$upper <- apply(proj_deaths_standardised,1,quantile, prob = 0.975)
  #  proj_deaths_standardised_long <- gather(as.data.frame(proj_deaths_standardised), key = "iteration", value = "hbv_death_rate", -time)

  # Label dataframes with scenario:
  proj_prev_summary$scenario <- scenario_label
  proj_inc_summary$scenario <- scenario_label
  proj_inc_rate_summary$scenario <- scenario_label
  proj_deaths_summary$scenario <- scenario_label
  proj_deaths_rate_summary$scenario <- scenario_label
  proj_deaths_standardised_summary$scenario <- scenario_label

  return(list(hbsag_prev_summary = proj_prev_summary,
              chronic_incidence_summary = proj_inc_summary,
              chronic_incidence_rate_summary = proj_inc_rate_summary,
              hbv_deaths_summary = proj_deaths_summary,
              hbv_deaths_rate_summary = proj_deaths_rate_summary,
              proj_deaths_standardised_summary = proj_deaths_standardised_summary))
}

# Function to calculate age-standardised rate of HBV-related deaths (called by extract_outcomes)
calculate_age_standardised_hbvdeaths_rate <- function(output_file) {
  # Age-standardised incidence of HBV-related deaths per 100000 per timestep
  # a) Calculate crude age-specific rates per person-year: need age-specific number of deaths and age-specific population size
  deaths_by_age <- output_file$out_cum_hbv_deathsf+
    output_file$out_cum_hbv_deathsm+
    output_file$out_cum_screened_hbv_deathsf+
    output_file$out_cum_screened_hbv_deathsm+
    output_file$out_cum_treated_hbv_deathsf+
    output_file$out_cum_treated_hbv_deathsm
  deaths_by_age <- calculate_incident_numbers(deaths_by_age)

  # Crude rate
  deaths_rate <- deaths_by_age/output_file$pop
  deaths_rate[is.na(deaths_rate)==TRUE] <- 0

  # b) Multiply by reference population (Gambian pop in 2020)
  ref_pop <- output_file$pop[output_file$time==2020,]
  deaths_rate_ref_pop <- sweep(as.matrix(deaths_rate), MARGIN=2, as.matrix(ref_pop), "*")

  # c) Calculate total expected deaths (sum of age-specific numbers) and divide by total Gambian popsize in 2020
  age_stand_rate <- rowSums(deaths_rate_ref_pop)/rowSums(ref_pop)

  return(age_stand_rate)
}

### Functions to run projection scenarios ----
### These functions apply run_model and code_model_output to a collection of parameter sets
run_hbsag_screening_scenarios <- function(..., default_parameter_list, calibrated_parameter_sets,
                                          parms_to_change = list(...)) {

  # Screen every 6 months
  sim_0point5 <- apply(calibrated_parameter_sets, 1,
                          function(x) run_model(sim_duration = runtime,
                           default_parameter_list = default_parameter_list,
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
                                  screening_years =  seq(2020,2100, by = 0.5)),
                           scenario = "vacc_screen"))
  out_0point5 <- lapply(sim_0point5,code_model_output)

  # Screen every year
  sim_1 <- apply(calibrated_parameter_sets, 1,
                       function(x) run_model(sim_duration = runtime,
                                             default_parameter_list = default_parameter_list,
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
                                                    screening_years =  seq(2020,2100, by = 1)),
                                             scenario = "vacc_screen"))
  out_1 <- lapply(sim_1,code_model_output)

  # Screen every 2 years
  sim_2 <- apply(calibrated_parameter_sets, 1,
                    function(x) run_model(sim_duration = runtime,
                                          default_parameter_list = default_parameter_list,
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
                                                 screening_years =  seq(2020,2100, by = 2)),
                                          scenario = "vacc_screen"))
  out_2 <- lapply(sim_2,code_model_output)

  # Screen every 5 years
  sim_5 <- apply(calibrated_parameter_sets, 1,
                    function(x) run_model(sim_duration = runtime,
                                          default_parameter_list = default_parameter_list,
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
                                                 screening_years =  seq(2020,2100, by = 5)),
                                          scenario = "vacc_screen"))
  out_5 <- lapply(sim_5,code_model_output)

  # Screen every 10 years
  sim_10 <- apply(calibrated_parameter_sets, 1,
                    function(x) run_model(sim_duration = runtime,
                                          default_parameter_list = default_parameter_list,
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
                                                 screening_years =  seq(2020,2100, by = 10)),
                                          scenario = "vacc_screen"))
  out_10 <- lapply(sim_10,code_model_output)

  # Screen every 15 years
  sim_15 <- apply(calibrated_parameter_sets, 1,
                    function(x) run_model(sim_duration = runtime,
                                          default_parameter_list = default_parameter_list,
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
                                                 screening_years =  seq(2020,2100, by = 15)),
                                          scenario = "vacc_screen"))
  out_15 <- lapply(sim_15,code_model_output)


  # Screen every 20 years
  sim_20 <- apply(calibrated_parameter_sets, 1,
                    function(x) run_model(sim_duration = runtime,
                                          default_parameter_list = default_parameter_list,
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
                                                 screening_years =  seq(2020,2100, by = 20)),
                                          scenario = "vacc_screen"))
  out_20 <- lapply(sim_20,code_model_output)

  # Screen every 25 years
  sim_25 <- apply(calibrated_parameter_sets, 1,
                    function(x) run_model(sim_duration = runtime,
                                          default_parameter_list = default_parameter_list,
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
                                                 screening_years =  seq(2020,2100, by = 25)),
                                          scenario = "vacc_screen"))

  out_25 <- lapply(sim_25,code_model_output)

  # Screen every 30 years
  sim_30 <- apply(calibrated_parameter_sets, 1,
                    function(x) run_model(sim_duration = runtime,
                                          default_parameter_list = default_parameter_list,
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
                                                 screening_years =  seq(2020,2100, by = 30)),
                                          scenario = "vacc_screen"))
  out_30 <- lapply(sim_30,code_model_output)

  # One off screen in 2020
  sim_once <- apply(calibrated_parameter_sets, 1,
                     function(x) run_model(sim_duration = runtime,
                                           default_parameter_list = default_parameter_list,
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
                                                  screening_years =  c(2020)),
                                           scenario = "vacc_screen"))
  out_once <- lapply(sim_once,code_model_output)


  # Status quo scenario: no screening
  sim_sq <- apply(calibrated_parameter_sets, 1,
                     function(x) run_model(sim_duration = runtime,
                                           default_parameter_list = default_parameter_list,
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
  out_sq <- lapply(sim_sq,code_model_output)


  outlist <- list("screen_0point5" = out_0point5,
                  "screen_1" = out_1,
                  "screen_2" = out_2,
                  "screen_5" = out_5,
                  "screen_10" = out_10,
                  "screen_15" = out_15,
                  "screen_20" = out_20,
                  "screen_25" = out_25,
                  "screen_30" = out_30,
                  "screen_once" = out_once,
                  "status_quo" = out_sq)

  return(outlist)
}

# Scenario run is a screening and treatment streategy - year of screening can be specified
run_one_hbsag_screening_scenario <- function(..., default_parameter_list, calibrated_parameter_sets,
                                          parms_to_change = list(...), years_of_test, label) {

  sim <- apply(calibrated_parameter_sets, 1,
                       function(x) run_model(sim_duration = runtime,
                                             default_parameter_list = default_parameter_list,
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
                                                    screening_years = years_of_test),
                                             scenario = "vacc_screen"))

  out <- lapply(sim, code_model_output)

  outlist <- list("screen" = out)
  names(outlist) <- label

  return(outlist)
}

# Scenario can be specified
run_one_scenario <- function(..., default_parameter_list, calibrated_parameter_sets,
                                             parms_to_change = list(...), scenario) {

  # Status quo scenario: no screening
  sim <- apply(calibrated_parameter_sets, 1,
                  function(x) run_model(sim_duration = runtime,
                                        default_parameter_list = default_parameter_list,
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
                                        scenario = scenario))
  out <- lapply(sim, code_model_output)

  label <- paste0("out_", scenario)

  outlist <- list("out" = out)
  names(outlist) <- label

  return(outlist)
}

### MODEL INPUT ----

# DEMOGRAPHY: SEE ABOVE

# INITIAL POPULATION
# Set up initial population: age- and sex-specific population size in 1950
# Note: names in initial population vector is reproduced in output
# guessed proportion in each compartment

output_storage <- c("cum_deathsf" = rep(0,n_agecat), "cum_deathsm" = rep(0,n_agecat),
                    "cum_infectionsf" = rep(0,n_agecat), "cum_infectionsm" = rep(0,n_agecat),
                    "cum_chronic_infectionsf" = rep(0,n_agecat), "cum_chronic_infectionsm" = rep(0,n_agecat),
                    "cum_births" = 0, "cum_infected_births" = 0, "cum_chronic_births" = 0,
                    "cum_hbv_deathsf" = rep(0,n_agecat), "cum_hbv_deathsm" = rep(0,n_agecat),
                    "cum_hcc_deathsf" = rep(0,n_agecat), "cum_hcc_deathsm" = rep(0,n_agecat),
                    "cum_eag_lossf" = rep(0,n_agecat), "cum_eag_lossm" = rep(0,n_agecat),
                    "cum_sag_lossf" = rep(0,n_agecat), "cum_sag_lossm" = rep(0,n_agecat),
                    "cum_incident_dccf" = rep(0,n_agecat), "cum_incident_dccm" = rep(0,n_agecat),
                    "cum_incident_hccf" = rep(0,n_agecat), "cum_incident_hccm" = rep(0,n_agecat),
                    "cum_dcc_to_hccf" = rep(0,n_agecat), "cum_dcc_to_hccm" = rep(0,n_agecat),
                    "cum_background_deaths_ldf" =  rep(0,n_agecat),
                    "cum_background_deaths_ldm" =  rep(0,n_agecat),
                    "cum_dcc_deathsf" = rep(0,n_agecat), "cum_dcc_deathsm" = rep(0,n_agecat),
                    "cum_cc_deathsf" = rep(0,n_agecat), "cum_cc_deathsm" = rep(0,n_agecat),
                    "cum_ir_to_ccf" = rep(0,n_agecat), "cum_ir_to_ccm" = rep(0,n_agecat),
                    "cum_enchb_to_ccf" = rep(0,n_agecat), "cum_enchb_to_ccm" = rep(0,n_agecat),
                    "cum_cc_to_hccf" = rep(0,n_agecat), "cum_cc_to_hccm" = rep(0,n_agecat),
                    "cum_it_to_hccf" = rep(0,n_agecat), "cum_it_to_hccm" = rep(0,n_agecat),
                    "cum_ir_to_hccf" = rep(0,n_agecat), "cum_ir_to_hccm" = rep(0,n_agecat),
                    "cum_ic_to_hccf" = rep(0,n_agecat), "cum_ic_to_hccm" = rep(0,n_agecat),
                    "cum_enchb_to_hccf" = rep(0,n_agecat), "cum_enchb_to_hccm" = rep(0,n_agecat),
                    "cum_screened_infectionsf" = rep(0,n_agecat), "cum_screened_infectionsm" = rep(0,n_agecat),
                    "cum_screened_chronic_infectionsf" = rep(0,n_agecat), "cum_screened_chronic_infectionsm" = rep(0,n_agecat),
                    "cum_screened_hbv_deathsf" = rep(0,n_agecat), "cum_screened_hbv_deathsm" = rep(0,n_agecat),
                    "cum_screened_cc_deathsf" = rep(0,n_agecat), "cum_screened_cc_deathsm" = rep(0,n_agecat),
                    "cum_screened_dcc_deathsf" = rep(0,n_agecat), "cum_screened_dcc_deathsm" = rep(0,n_agecat),
                    "cum_screened_hcc_deathsf" = rep(0,n_agecat), "cum_screened_hcc_deathsm" = rep(0,n_agecat),
                    "cum_screened_incident_dccf" = rep(0,n_agecat), "cum_screened_incident_dccm" = rep(0,n_agecat),
                    "cum_screened_incident_hccf" = rep(0,n_agecat), "cum_screened_incident_hccm" = rep(0,n_agecat),
                    "cum_treated_hbv_deathsf" = rep(0,n_agecat), "cum_treated_hbv_deathsm" = rep(0,n_agecat),
                    "cum_treated_dcc_deathsf" = rep(0,n_agecat), "cum_treated_dcc_deathsm" = rep(0,n_agecat),
                    "cum_treated_hcc_deathsf" = rep(0,n_agecat), "cum_treated_hcc_deathsm" = rep(0,n_agecat),
                    "cum_treated_incident_hccf" = rep(0,n_agecat), "cum_treated_incident_hccm" = rep(0,n_agecat),
                    "cum_treated_sag_lossf" = rep(0,n_agecat), "cum_treated_sag_lossm" = rep(0,n_agecat),
                    "total_screened_uninfected" = 0, "total_screened_chb" = 0,
                    "total_screened_cirrhosis" = 0, "total_screened_ineligible" = 0)

init_pop <- c("Sf" = popsize_1950$pop_female*(1-gambia_infected),
              "ITf" = popsize_1950$pop_female*gambia_infected*gambia_eag*0.5,
              "IRf" = popsize_1950$pop_female*gambia_infected*gambia_eag*0.5,
              "ICf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
              "ENCHBf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
              "CCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
              "DCCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
              "HCCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
              "Rf" = rep(0,n_agecat),
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
              "T_Rf" =rep(0,n_agecat),
              "Sm" = popsize_1950$pop_male*(1-gambia_infected),
              "ITm" = popsize_1950$pop_male*gambia_infected*gambia_eag*0.5,
              "IRm" = popsize_1950$pop_male*gambia_infected*gambia_eag*0.5,
              "ICm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
              "ENCHBm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
              "CCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
              "DCCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
              "HCCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
              "Rm" = rep(0,n_agecat),
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


# Total population in 1950:
N0 <- sum(init_pop[1:(n_infectioncat * n_agecat * 2)])

## TRANSMISSION, NATURAL HISTORY AND INTERVENTION PARAMETERS

parameter_list <- list(
  # TRANSMISSION PARAMETERS
  b1 = 0.1,     # Parameter not the same as in original model. Margaret value in Ethiopia is 0.2027
  b2 = 0.009,   # Parameter not the same as in original model. Margaret value in Ethiopia is 0.001
  b3 = 0.001,   # Parameter not the same as in original model. Margaret value in Ethiopia
  alpha = 15,         # Relative infectiousness of eAg-positives (Shevanthi value)
  mtct_prob_e = 0.9,  # Shevanthi value, probability of perinatal transmission from HBeAg-positive mother
  mtct_prob_s = 0.2,  # S value, Margaret value in Ethiopia is 0.3681, probability of perinatal transmission from HBeAg-negative infected mother
  # NATURAL HISTORY PROGRESSION RATES AND PARAMETERS
  p_chronic_in_mtct = 0.89,  # Edmunds value for risk of chronic carriage in <0.5 year olds (infected perinatally)
  p_chronic_function_r = 0.65,  # Edmunds decay rate parameter in exponential function of age-specific chronic carriage risk
  p_chronic_function_s = 0.46,  # Edmunds s parameter in exponential function of age-specific chronic carriage risk
  pr_it_ir = 0.20374,  # S
  pr_ir_ic = 0.10187,  # S
  eag_prog_function_rate = 0.001,  # S Ethiopia, 0.1281 in Margaret's Ethiopia fit
  pr_ir_enchb = 0.005,  # S
  # pr_ir_enchb_age_threshold = 20,  # M, not present in S model
  pr_ir_cc_female = 0.028,  # S but not present in published model
  pr_ir_cc_age_threshold = 20,  # M, not present in S model
  pr_ic_enchb = 0.01,   # S
  #  sag_loss_intercept = -0.001697,  # intercept of linear function with age
  sag_loss_slope = 0.0004106,  # slope of linear function with age
  pr_enchb_cc_female = 0.04,  # Progression to CC (from ENCHB) S value
  # cirrhosis_prog_coefficient = 0.0341,  # Margaret value
  # pr_enchb_cc_age_threshold = 20,  # M value, not present in S model
  cirrhosis_male_cofactor = 12.32,  # M value
  pr_cc_dcc = 0.04,  # Progression from CC to DCC, S value
  # PROGRESSION RATES TO HEPATOCELLULAR CARCINOMA
  cancer_prog_coefficient_female = 4.0452e-05,  # value from Margaret
  cancer_age_threshold = 10,  # value from Margaret
  cancer_male_cofactor = 5.2075,  # value from Margaret
  hccr_it = 1,  # S
  hccr_ir = 2, # S
  hccr_enchb = 2, # S
  hccr_cc = 13, # S
  hccr_dcc = 0.04, # S
  # HBV-RELATED MORTALITY RATES (MORTALITY FROM LIVER DISEASE)
  mu_cc = 0.039,  # S
  mu_dcc = 0.314,  # S
  mu_hcc = 0.5,  # S
  # INFANT VACCINATION PARAMETERS
  # Vaccine coverage is read in through a function
  vacc_eff = 0.95,                           # vaccine efficacy, S
  vacc_introtime = 1990,                     # year of vaccine introduction
  # BIRTH DOSE VACCINATION PARAMETERS
  bdvacc_introtime = 2020,                   # year of birth dose vaccine introduction
  bdvacc_cov = 0.05,                          # birth dose vaccine coverage
  mtct_prob_ebd = 0.32,                      # MTCT risk from eAg-positive mother with birth dose (95% CI from Keane = 0.1-0.58)
  mtct_prob_sbd = 0,                         # MTCT risk from eAg-negative mother with birth dose (Keane)
  mtct_prob_treatbd = 0,                     # MTCT risk from treated carrier mother with birth dose (assumption, not peripartum therapy)
  # TREATMENT PARAMETERS
  screening_years = c(2020),                 # vectors for years to implement the screening programme
  screening_coverage = 0.7,                  # proportion of population in given age group to screen
  min_age_to_screen = 30,                    # Minimum age group to screen
  max_age_to_screen = 70,                    # Maximum age group to screen
  link_to_care_prob = 0.81,                  # probability of linkage to care (liver disease assessment) after HBsAg test
  treatment_initiation_prob = 1,             # probability of initiating treatment after diagnosis of treatment eligibility
  monitoring_rate = 1/5,                     # annual rate of monitoring for treatment eligibility
  monitoring_prob = 1,                       # probability of being monitored (1-proportion lost to follow-up)
  alpha2 = 1,                                # relative infectiousness with treatment compared to HBeAg-negatives
  mtct_prob_treat_cofactor = 1,              # relative infectiousness of mother-to-child transmission risk from treated mother (NOT peripartum therapy) compared to HBeAg-negative mother
#  tr_vir_supp = 2,                           # rate of achieveing virological suppression after treatment initiation
  thccr_chb = 0.27,                          # hazard ratio for progression to HCC from CHB on treatment
  thccr_cc = 0.23,                           # hazard ratio for progression to HCC from CC on treatment
  thccr_dcc = 0.23,                          # hazard ratio for progression to HCC from DCC on treatment
  tmu_dcc = 0.18,                            # mortality rate from DCC after treatment
  tmu_hcc_cofactor = 1,                      # hazard ratio for reduction of HCC mortality after treatment compared to no treatment
  # SIMULATION PARAMETERS
  sim_starttime = starttime,
  # INTERVENTION ON/OFF SWITCH (1/0)
  apply_vacc = 1,
  apply_bdvacc = 0,
  # DEMOGRAPHY ON/OFF SWITCH (1/0)
  births_on = 1,
  migration_on = 1,
  mortality_on = 1)

# Store names of all parameters
parameter_names <- names(parameter_list)



