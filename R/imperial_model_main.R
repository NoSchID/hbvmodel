###################################
### Imperial HBV model 25/06/19 ###
###################################
# Model described in Shevanthi's thesis and adapted by Margaret
# New adaptations in this version:
# - male cofactor on IR=>CC progression
# - 0.5 year olds participating in horizontal transmission
# Currently only infant vaccination in second age group, no birth dose or treatment
# Solving using ODE45

### Load packages ----
require(here)  # for setting working directory
require(tidyr)  # for data processing
require(dplyr)  # for data processing
require(deSolve)  # ODE solver
require(lhs)  # Latin Hypercube sampling
require(parallel)  # for running code in parallel

# Packages for making plots
require(ggplot2)  # for calibration plots
require(gridExtra)  # for calibration plots
require(grid)  # for calibration plots

# Packages for timing, profiling and making code more efficient
library(tictoc)  # for timing code
library(profvis)  # for profiling
library(efficient)  # for profiling


### Define simulation parameters ----
## Country
countryname <- "gambia"

## Times
dt <- 0.5                              # timestep (years)
# dt/da can only be 0.5 to match the WAIFW matrix at the moment
# For demography, it can be up to 1, or multiple of 5 thereafter (because of women of childbearing age)
starttime <- 1850
runtime <- 250                     # number of years to run the model for
#times <- round((0:(runtime/dt))*dt,2) # vector of timesteps
#times_labels <- times+starttime       # year labels for timestep vector

## Age groups
da <- dt                              # time spent in each age group (years)
ages <- round((0:((100-da)/da))*da,2)  # vector of all age groups
n_agecat <- length(ages)               # number of age groups

ages_wocba <- round(((15/da):((50-da)/da))*da,2)        # age groups 15-49 years (women of childbearing age)

## Infection compartments
n_infectioncat <- 9                      # Number of infection compartments

## Definition of indices
index <- list("infcat_all" = 1:n_infectioncat,            # index for all infection status compartments
              "ages_all" = 1:n_agecat,                    # index for all age groups
              "ages_wocba" = which(ages == min(ages_wocba)):which(ages == max(ages_wocba)), # index for age group 15-49 years (women of childbearing age)
              "ages_0to1" = which(ages == 0):(1/da),
              "ages_1to5" = which(ages == 1):which(ages == 6-da),       # index for age groups 1-5 years
              "ages_6to15" = which(ages == 6):which(ages == 16-da),     # index for age groups 6-15 years
              "ages_16to100" = which(ages == 16):n_agecat)          # index for age groups 16-100 years

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
       xout = input_who_vaccine_coverage$year, method = "linear")$y

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
    S <- 1                           # Susceptible
    IT <- 2                           # Chronic infection: immune tolerant
    IR <- 3                           # Chronic infection: immune reactive
    IC <- 4                           # Chronic infection: inactive carrier
    ENCHB <- 5                        # Chronic infection: HBeAg-negative CHB
    CC <- 6                           # Chronic disease: compensated cirrhosis
    DCC <- 7                          # Chronic disease: decompensated cirrhosis
    HCC <- 8                          # Chronic disease: hepatocellular carcinoma
    R <- 9                           # Immune
    # Grouping of infected compartments
    HBeAg_neg <- IC:HCC               # HBeAg-negative infected compartments
    HBeAg_pos <- c(IT,IR)             # HBeAg-positive infected compartments

    # Infant vaccination: set date for introduction
    # Vaccination coverage is 0 until the specified starttime and only if vaccine switch is on
    # Vaccine coverage varies over time
    if (apply_vacc == 1 & timestep >= (vacc_introtime-sim_starttime)) {
      vacc_cov = timevary_vaccine_coverage(timestep + sim_starttime)
    } else {
      vacc_cov = 0
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



    # TRANSMISSION

    # Mother-to-child transmission and births
    dcum_infected_births <- sum(fertility_rate *
                                  (rowSums(mtct_prob_e * pop[index$ages_wocba,c(IT,IR),1]) +
                                     rowSums(mtct_prob_s * pop[index$ages_wocba,IC:HCC,1])))
    # infected births come from acute and chronic women of childbearing age
    dcum_chronic_births <- p_chronic_function[1] * dcum_infected_births # infected babies becoming chronic carriers
    dcum_uninfected_births <- sum(fertility_rate * pop[index$ages_wocba,index$infcat_all,1]) -
      dcum_infected_births  # uninfected births = all births - infected babies
    dcum_nonchronic_births <- sum(fertility_rate * pop[index$ages_wocba,index$infcat_all,1]) -
      dcum_chronic_births
    dcum_births <- dcum_infected_births + dcum_uninfected_births

    #births <- sum(fertility_rate * pop[index$ages_wocba,index$infcat_all,1])
    # applying the same age-specific fertility rate to every infection compartment

    # Horizontal transmission: Age-specific force of infection (same for men and women)

    # Define indices for contact groups
    #i_1to4 <- which(ages == 1):which(ages == (5-da))
    #i_1to14 <- which(ages == 1):which(ages == (15-da))
    #i_5plus <- which(ages == 5):which(ages == (100-da))

    # Imperial model FOI

    # Set up vector to store the age-specific FOI
    #foi <- rep(0,n_agecat)
    # FOI experienced by 1-4 year olds
    #foi[i_1to4] <- b1 * sum(apply(pop[i_1to4,HBeAg_neg,1:2],1,sum))/sum(pop[i_1to4,index$infcat_all,1:2]) +
    #                        min(1,b1 * alpha) * sum(apply(pop[i_1to4,HBeAg_pos,1:2],1,sum))/sum(pop[i_1to4,index$infcat_all,1:2])
    # FOI experienced by 1-14 year olds
    #foi[i_1to14] <- foi[i_1to14] + b2 * sum(apply(pop[i_1to14,HBeAg_neg,1:2],1,sum))/sum(pop[i_1to14,index$infcat_all,1:2]) +
    #                          (b2 * alpha) * sum(apply(pop[i_1to14,HBeAg_pos,1:2],1,sum))/sum(pop[i_1to14,index$infcat_all,1:2])
    # FOI experienced by 5+ year olds
    #foi[i_5plus] <- foi[i_5plus] + b3 * sum(apply(pop[i_5plus,HBeAg_neg,1:2],1,sum))/sum(pop[i_5plus,index$infcat_all,1:2]) +
    #                            (b3 * alpha) * sum(apply(pop[i_5plus,HBeAg_pos,1:2],1,sum))/sum(pop[i_5plus,index$infcat_all,1:2])

    # Alternative force of infection definition with WAIFW matrix

    # Define WAIFW matrix
    # Age-dependent mixing between 4 age groups: 0 year olds, 1-5 years, 6-15 years, 16-100 years
    # Assuming no effective contact between children (1-5 years) and adults (>15 years)
    # Assuming no horizontal transmission from and to infants (0 years old)
    beta <- matrix(0, nrow = 4, ncol = 4)  # matrix of transmission parameters
    beta[2,2] <- b1                        # transmission among children 1-5 years
    beta[3,3] <- b2                        # transmission among juveniles 6-15 years
    beta[4,4] <- b3                        # transmission among adults 16-100 years
    beta[2,3] <- b2                        # transmission from juveniles to children
    beta[3,2] <- b2                        # = transmission from children to juveniles
    beta[3,4] <- b3                        # transmission from adults to juveniles
    beta[4,3] <- b3                        # = transmission from juveniles to adults

    # Define a vector of the age-specific prevalence of infectious individuals:
    # Infectious compartments are IT, IR, IC, ENCHB, CC, DCC, HCC
    # HBeAg-positive individuals (IT, IR) are more infectious than HBeAg-negatives (multiply by alpha)
    # Sum prevalence in HBeAg-negatives and HBeAg-positives multiplied by alpha
    # Returns 1 number per transmission age group (4 total)
    i_1to4 <- which(ages == 1):which(ages == (5-da))
    i_5to14 <- which(ages == 5):which(ages == (15-da))
    i_15to100 <- which(ages == 15):which(ages == (100-da))

    infectious_vector <- c(sum(pop[index$ages_0to1,HBeAg_neg,1:2])/sum(pop[index$ages_0to1,index$infcat_all,1:2]) +
                             (alpha * sum(pop[index$ages_0to1,HBeAg_pos,1:2])/sum(pop[index$ages_0to1,index$infcat_all,1:2])), # 0 year olds
                           sum(pop[i_1to4,HBeAg_neg,1:2])/sum(pop[i_1to4,index$infcat_all,1:2]) +
                             (alpha * sum(pop[i_1to4,HBeAg_pos,1:2])/sum(pop[i_1to4,index$infcat_all,1:2])), # 1-4 year olds
                           sum(pop[i_5to14,HBeAg_neg,1:2])/sum(pop[i_5to14,index$infcat_all,1:2]) +
                             (alpha * sum(pop[i_5to14,HBeAg_pos,1:2])/sum(pop[i_5to14,index$infcat_all,1:2])), # 5-14 year olds
                           sum(pop[i_15to100,HBeAg_neg,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]) +
                             (alpha * sum(pop[i_15to100,HBeAg_pos,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]))) # 15-100 year olds

    # 0.5 year olds can get infected:
    #i_1to4 <- which(ages == 0.5):which(ages == (5-da))
    #i_5to14 <- which(ages == 5):which(ages == (15-da))
    #i_15to100 <- which(ages == 15):which(ages == (100-da))

    #infectious_vector <- c(sum(pop[1,HBeAg_neg,1:2])/sum(pop[1,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[1,HBeAg_pos,1:2])/sum(pop[1,index$infcat_all,1:2])), # 0 year olds
    #                       sum(pop[i_1to4,HBeAg_neg,1:2])/sum(pop[i_1to4,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[i_1to4,HBeAg_pos,1:2])/sum(pop[i_1to4,index$infcat_all,1:2])), # 1-5 year olds
    #                       sum(pop[i_5to14,HBeAg_neg,1:2])/sum(pop[i_5to14,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[i_5to14,HBeAg_pos,1:2])/sum(pop[i_5to14,index$infcat_all,1:2])), # 6-15 year olds
    #                       sum(pop[i_15to100,HBeAg_neg,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[i_15to100,HBeAg_pos,1:2])/sum(pop[i_15to100,index$infcat_all,1:2]))) # 16-100 year olds


    # Multiply WAIFW matrix by the age-specific proportion of infectious individuals
    # Returns a vector with force of infection for every age - 4 different values:
    # 0 in 0-year olds, different values for 1-5, 6-15 and 16-100 year olds
    foi_unique <- beta %*% infectious_vector
    # Repeat these values for every 1 year age group (assuming 0.5 year olds can't get horizontally infected)
    foi <- c(rep(foi_unique[1], times = 2),
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

      # Infection: Incident infections (all) and incident chronic infections
      dcum_infections[index$ages_all,i] <- foi * pop[index$ages_all,S,i]
      dcum_chronic_infections[index$ages_all,i] <-
        p_chronic_function * dcum_infections[index$ages_all,i]
      # Returns a matrix with incident infections for every age (rows) and every sex (columns)

      # Natural history transitions
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

      # Incidence of HBeAg loss: transition from IR to IC and IR to ENCHB
      dcum_eag_loss[index$ages_all,i] <-
        pr_ir_ic * eag_prog_function * pop[index$ages_all,IR,i] +
        pr_ir_enchb * pop[index$ages_all,IR,i]

      # Transition from IC to R (sAg loss)
      dcum_sag_loss[index$ages_all,i] <- sag_loss * pop[index$ages_all,IC,i]

      # DCC incidence
      dcum_dcc[index$ages_all,i] <- pr_cc_dcc * pop[index$ages_all,CC,i]

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

      # Transitions between compartments:

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
        mu_cc * pop[index$ages_all,CC,i] -
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
        # vacc_cov * vacc_eff * pop[index$ages_all,S,i] +
        (1-p_chronic_function) * dcum_infections[index$ages_all,i] +
        dcum_sag_loss[index$ages_all,i] -
        deaths[index$ages_all,R,i] + migrants[index$ages_all,R,i]

      # Babies are born susceptible or infected (age group 1)
      dpop[1,S,i] <- dpop[1,S,i] + sex_ratio[i] * dcum_nonchronic_births
      dpop[1,IT,i] <- dpop[1,IT,i] + sex_ratio[i] * dcum_chronic_births
      #dpop[1,R,i] <- dpop[1,R,i] + sex_ratio[i] * (1-p_chronic_function[1]) * infected_births

      # Vaccination: applied at 0.5 years of age (this only makes sense if
      # age step is 0.5!)
      dpop[2,S,i] <- dpop[2,S,i] - (vacc_cov * vacc_eff * pop[2,S,i])
      dpop[2,R,i] <- dpop[2,R,i] + (vacc_cov * vacc_eff * pop[2,S,i])

    }

    # OUTPUT


    # Sum age-specific number of incident background deaths across infection compartments for output
    dcum_deaths <- cbind(rowSums(deaths[,,1]), rowSums(deaths[,,2]))
    # Age-specific number of incident background deaths among liver disease patients (CC, DCC and HCC)
    dcum_background_deaths_ld <- cbind(rowSums(deaths[index$ages_all,CC:HCC,1]),
                                       rowSums(deaths[index$ages_all,CC:HCC,2]))

    # Return results
    res <- c(dpop, dcum_deaths, dcum_infections, dcum_chronic_infections,
             dcum_births, dcum_infected_births, dcum_chronic_births,
             dcum_hbv_deaths, dcum_hcc_deaths, dcum_eag_loss,
             dcum_sag_loss, dcum_dcc, dcum_hcc, dcc_to_hcc_transitions,
             dcum_background_deaths_ld, dcum_dcc_deaths, ir_to_cc_transitions,
             enchb_to_cc_transitions, cc_to_hcc_transitions,
             it_to_hcc_transitions, ir_to_hcc_transitions, ic_to_hcc_transitions,
             enchb_to_hcc_transitions)
    list(res, p_chronic_function = p_chronic_function)

  })

}

### Model-related functions ----

## Function to interpolate demographic parameters over time - specific to these datasets
timevary_parameters_old <- function(timestep, dataset) {
  # Input datasets are matrices of age-specific mortality rates, birth rate and migration rate for every 5-year period
  res <- rep(0,ncol(dataset))
  for (i in 2:ncol(dataset)) {
    res[i] <- spline(x = dataset[,1], y = dataset[,i], xout = timestep)[[2]]
  }
  return(res[-1])
} # for loop instead of apply: 69.36s minimally quicker

## Event function: reset population size to initial (1850) size in 1950
reset_pop_1950 <- function(timestep, pop, parameters){
  with (as.list(pop),{
    pop_to_reset <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))
    initialpop <- array(unlist(init_pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    current_pop <- cbind(rowSums(pop_to_reset[,,1]), rowSums(pop_to_reset[,,2]))
    pop_increase <- cbind(rowSums(initialpop[,,1]), rowSums(initialpop[,,2]))/current_pop
    scaler <- array(c(rep(pop_increase[,1], n_infectioncat),
                      rep(pop_increase[,2], n_infectioncat)), dim = c(n_agecat,n_infectioncat,2))
    pop_to_reset <- scaler * pop_to_reset
    return(c(pop_to_reset, unlist(init_pop[(2*n_infectioncat*n_agecat+1):length(init_pop)])))
    #   return(c(pop_to_reset, rep(0,n_agecat), rep(0,n_agecat), rep(0,n_agecat),
    #            rep(0,n_agecat), rep(0,n_agecat), rep(0,n_agecat),
    #            0, 0, 0,
    #            rep(0,n_agecat), rep(0,n_agecat)))
  })
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
  if (scenario == "vacc") {
    parameters$apply_vacc <- 1
  } else if (scenario == "no_vacc") {
    parameters$apply_vacc <- 0
  } else {
    print("Not a valid scenario. Options: vacc, no_vacc")
  }

  ## Run model simulation
  timestep_vector <- round((0:((sim_duration-dt)/dt))*dt,2)
  timestep_labels <- timestep_vector + parameters$sim_starttime

  if (1950 %in% timestep_labels) {
    timestep_1950 <- timestep_vector[which(timestep_labels == 1950)]
    out <- as.data.frame(ode.1D(y = init_pop_vector, times = timestep_vector, func = imperial_model,
                                parms = parameters, nspec = 1, method = "ode45",
                                events = list(func = reset_pop_1950, time = timestep_1950)))
  } else {
  out <- as.data.frame(ode.1D(y = init_pop_vector, times = timestep_vector, func = imperial_model,
                              parms = parameters, nspec = 1, method = "ode45"))
  }

  # Add year label to timestep
  out$time   <-  out$time + parameters$sim_starttime

  return(out)

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

  } else {  # if not use operation on whole data frame

    incident_numbers <- rbind(cumulative_output[1,],  # number at first timestep
                              apply(cumulative_output, 2, diff, lag = 1))
    # number at current timestep - number at previous timestep

  }

  # Returns a vector/dataframe of new cases since the last timestep
  return(incident_numbers)
}

# Function to code relevant model output (stored in list)
code_model_output <- function(output) {

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

  ## Process infection outputs
  # Combine into data frames with outputs of interest for further analysis

  # Age-specific number in each infection compartment at each time step
  sus <- data.frame(pop = out_sf + out_sm)   # need to change the column names
  carriers <- data.frame(pop = (out_itf + out_itm +
                                  out_irf + out_irm +
                                  out_icf+out_icm+
                                  out_enchbf+out_enchbm+
                                  out_ccf+out_ccm+
                                  out_dccf+out_dccm+
                                  out_hccf+out_hccm))
  carriers_female <- data.frame(pop = (out_itf+
                                         out_irf+
                                         out_icf+
                                         out_enchbf+
                                         out_ccf+out_dccf+out_hccf))
  carriers_male <- data.frame(pop = (out_itm+
                                       out_irm+
                                       out_icm+
                                       out_enchbm+
                                       out_ccm+out_dccm+out_hccm))
  immune <- data.frame(pop = out_rf + out_rm)
  ever_infected <- data.frame(pop = carriers + immune)
  ever_infected_female <- data.frame(pop = carriers_female + out_rf)
  ever_infected_male <- data.frame(pop = carriers_male + out_rm)
  eag_positive <- data.frame(pop = (out_itf + out_itm +
                                      out_irf + out_irm))
  eag_positive_female <- data.frame(pop = (out_itf + out_irf))
  eag_positive_male <- data.frame(pop = (out_itm + out_irm))

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
                   "incident_infections" = incident_infections,
                   "incident_chronic_infections" = incident_chronic_infections,
                   "hbv_deaths" = hbv_deaths,
                   "incident_hcc" = incident_hcc,
                   "full_output" = output)
  return(toreturn)

}

# Shorter/quicker function to process model outputs:
# calculates carriers, ever infected, eAg positives and total pop
code_model_output_summary <- function(output) {

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
                                  out_hccf+out_hccm))
  carriers_female <- data.frame(pop = (out_itf+
                                         out_irf+
                                         out_icf+
                                         out_enchbf+
                                         out_ccf+out_dccf+out_hccf))
  carriers_male <- data.frame(pop = (out_itm+
                                       out_irm+
                                       out_icm+
                                       out_enchbm+
                                       out_ccm+out_dccm+out_hccm))

  immune <- data.frame(pop = out_rf + out_rm)

  ever_infected <- data.frame(pop = carriers + immune)
  ever_infected_female <- data.frame(pop = carriers_female + out_rf)
  ever_infected_male <- data.frame(pop = carriers_male + out_rm)

  eag_positive <- data.frame(pop = (out_itf + out_itm +
                                      out_irf + out_irm))
  eag_positive_female <- data.frame(pop = (out_itf + out_irf))
  eag_positive_male <- data.frame(pop = (out_itm + out_irm))

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
                   "pop" = pop)
  return(toreturn)

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
                    "cum_ir_to_ccf" = rep(0,n_agecat), "cum_ir_to_ccm" = rep(0,n_agecat),
                    "cum_enchb_to_ccf" = rep(0,n_agecat), "cum_enchb_to_ccm" = rep(0,n_agecat),
                    "cum_cc_to_hccf" = rep(0,n_agecat), "cum_cc_to_hccm" = rep(0,n_agecat),
                    "cum_it_to_hccf" = rep(0,n_agecat), "cum_it_to_hccm" = rep(0,n_agecat),
                    "cum_ir_to_hccf" = rep(0,n_agecat), "cum_ir_to_hccm" = rep(0,n_agecat),
                    "cum_ic_to_hccf" = rep(0,n_agecat), "cum_ic_to_hccm" = rep(0,n_agecat),
                    "cum_enchb_to_hccf" = rep(0,n_agecat), "cum_enchb_to_hccm" = rep(0,n_agecat))
init_pop <- c("Sf" = popsize_1950$pop_female*(1-gambia_infected),
               "ITf" = popsize_1950$pop_female*gambia_infected*gambia_eag*0.5,
               "IRf" = popsize_1950$pop_female*gambia_infected*gambia_eag*0.5,
               "ICf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "ENCHBf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "CCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "DCCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "HCCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "Rf" = rep(0,n_agecat),
               "Sm" = popsize_1950$pop_male*(1-gambia_infected),
               "ITm" = popsize_1950$pop_male*gambia_infected*gambia_eag*0.5,
               "IRm" = popsize_1950$pop_male*gambia_infected*gambia_eag*0.5,
               "ICm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "ENCHBm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "CCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "DCCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "HCCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "Rm" = rep(0,n_agecat),
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
  # SIMULATION PARAMETERS
  sim_starttime = starttime,
  # INTERVENTION ON/OFF SWITCH (1/0)
  apply_vacc = 1,
  # DEMOGRAPHY ON/OFF SWITCH (1/0)
  births_on = 1,
  migration_on = 1,
  mortality_on = 1)

# Store names of all parameters
parameter_names <- names(parameter_list)
