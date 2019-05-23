###################################
### Imperial HBV model 22/05/19 ###
###################################
# Model described in Shevanthi's thesis and adapted by Margaret
# Currently only infant vaccination in second age group, no birth dose or treatment
# Solving using ODE45

### Load packages ----
require(tidyr)  # for data processing
require(dplyr)  # for data processing
require(deSolve)  # ODE solver
library(tictoc)  # for timing code
require(here)  # for setting working directory
library(profvis)  # for code profiling
require(ggplot2)

### Simulation parameters ----
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

# Need to do this in prepared datasets!!
fert_rates[,1] <- seq(1850, (2100-0.5), 0.5)
mort_rates_female[,1] <- seq(1850, (2100-0.5), 0.5)
mort_rates_male[,1] <- seq(1850, (2100-0.5), 0.5)
migration_rates_female[,1] <- seq(1850, (2100-0.5), 0.5)
mogration_rates_male[,1] <- seq(1850, (2100-0.5), 0.5)

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

### THE MODEL ----
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
    i_1to4 <- which(ages == 1):which(ages == (5-da))
    i_1to14 <- which(ages == 1):which(ages == (15-da))
    i_5plus <- which(ages == 5):which(ages == (100-da))

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
    #infectious_vector <- c(sum(pop[1,HBeAg_neg,1:2])/sum(pop[1,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[1,HBeAg_pos,1:2])/sum(pop[1,index$infcat_all,1:2])), # 0 year olds
    #                       sum(pop[c(2,index$ages_1to5),HBeAg_neg,1:2])/sum(pop[c(2,index$ages_1to5),index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[c(2,index$ages_1to5),HBeAg_pos,1:2])/sum(pop[c(2,index$ages_1to5),index$infcat_all,1:2])), # 1-5 year olds
    #                       sum(pop[index$ages_6to15,HBeAg_neg,1:2])/sum(pop[index$ages_6to15,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[index$ages_6to15,HBeAg_pos,1:2])/sum(pop[index$ages_6to15,index$infcat_all,1:2])), # 6-15 year olds
    #                       sum(pop[index$ages_16to100,HBeAg_neg,1:2])/sum(pop[index$ages_16to100,index$infcat_all,1:2]) +
    #                         (alpha * sum(pop[index$ages_16to100,HBeAg_pos,1:2])/sum(pop[index$ages_16to100,index$infcat_all,1:2]))) # 16-100 year olds


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
      dcum_dcc[index$ages_all,i] <- dccrate * pop[index$ages_all,CC,i]

      # Transitions to HCC
      it_to_hcc_transitions[index$ages_all,i] <-
        hccr_it * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IT,i]

      ir_to_hcc_transitions[index$ages_all,i] <-
        hccr_ir * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IR,i]

      ic_to_hcc_transitions[index$ages_all,i] <-
        hccr_ic * cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IC,i]

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
        pr_ir_cc_function * pop[index$ages_all,IR,i]

      # ENCHB to CC transitions = cirrhosis incidence in HBeAg-negatives
      enchb_to_cc_transitions[index$ages_all,i] <-
        cirrhosis_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB,i]

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
        pr_ir_enchb_function * pop[index$ages_all,IR,i] -
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
        pr_ir_enchb_function * pop[index$ages_all,IR,i] +
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

# Run the model
# Old function, not used anymore but has likelihood
run_model_old <- function(b1 = b1, b2 = b2, b3 = b3, alpha = alpha,
                          mtct_prob_e = mtct_prob_e, mtct_prob_s = mtct_prob_s,
                          p_chronic = p_chronic, pr_it_ir = pr_it_ir, pr_ir_ic = pr_ir_ic,
                          pr_ir_enchb = pr_ir_enchb, pr_ic_enchb = pr_ic_enchb, sag_loss = sag_loss,
                          ccrate = ccrate, dccrate = dccrate,
                          hccr_it = hccr_it, hccr_ir = hccr_ir, hccr_ic = hccr_ic,
                          hccr_enchb = hccr_enchb, hccr_cc = hccr_cc, hccr_dcc = hccr_dcc,
                          mu_cc = mu_cc, mu_dcc = mu_dcc, mu_hcc = mu_hcc,
                          vacc_cov = vacc_cov, vacc_eff = vacc_eff, vacc_introtime = vacc_introtime) {

  # Add parameters into list
  parameters <- list(b1 = b1, b2 = b2, b3 = b3, alpha = alpha,
                     mtct_prob_e = mtct_prob_e, mtct_prob_s = mtct_prob_s,
                     p_chronic = p_chronic, pr_it_ir = pr_it_ir, pr_ir_ic = pr_ir_ic,
                     pr_ir_enchb = pr_ir_enchb, pr_ic_enchb = pr_ic_enchb, sag_loss = sag_loss,
                     ccrate = ccrate, dccrate = dccrate,
                     hccr_it, hccr_ir, hccr_ic, hccr_enchb, hccr_cc, hccr_dcc,
                     mu_cc = mu_cc, mu_dcc = mu_dcc, mu_hcc = mu_hcc,
                     vacc_cov = vacc_cov, vacc_eff = vacc_eff, vacc_introtime = vacc_introtime)

  # Run simulation
  out <- as.data.frame(ode.1D(y = init_pop, times = times, func = imperial_model,
                              parms = parameters, nspec = 1, method = "lsoda"))
  out$time   <-  out$time + starttime

  # Code carrier prevalence as output
  #pop_by_age <- out[,1+sindex] + out[,1+aindex] + out[,1+iindex] + out[,1+rindex]
  #prev_by_age <- out[,1+iindex]/pop_by_age
  #prop_everinf_by_age <- (out[,1+aindex] + out[,1+iindex] + out[,1+rindex])/pop_by_age

  # Data of number infected to fit to
  #data_prev <- as.numeric(edmunds_prev_by_age$y*pop_by_age[2000,])
  #data_everinf <- as.numeric(edmunds_everinf_by_age$y*pop_by_age[2000,])

  # Log likelihood
  #  LL <- sum(dbinom(x = round(data_prev), size = round(as.numeric(pop_by_age[2000,])),
  #                   prob = as.numeric(prev_by_age[2000,]), log = TRUE))

  toreturn <- out
  #toreturn <- list(out = out, prev_by_age = prev_by_age, loglikelihood = LL)
  #toReturn <- c(modelprev = as.numeric(prev_by_age[2000,]*100), LL = LL)

  return(toreturn)
}

# Try to prevent negative numbers using event (not working so far)
positive_fun <- function(timestep, pop, parameters){
  with(as.list(pop), {
    pop <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))
    pop[pop<1e-6] <- 0   # also tried with < 0
    return(c(pop, rep(0,n_agecat), rep(0,n_agecat), rep(0,n_agecat), rep(0,n_agecat)))
  })
}

# Function to run the model once for a given scenario
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
  # The risk is fixed at 0.89 for under-0.5 year olds (perinatal infection)
  parameters$p_chronic_function <- c(rep(0.89, 0.5/da),
                            exp(-parameters$p_chronic_function_r *
                                  ages[which(ages == 0.5):n_agecat]^parameters$p_chronic_function_s))

  # Age-specific function of progression through IT and IR (IT=>IR and IR=>IC)
  # Represented by an exponential function that decreases with age
  parameters$eag_prog_function <- parameters$eag_prog_function_intercept *
                                  exp(-parameters$eag_prog_function_rate * ages)

  # Age-specific progression to HCC from all carrier compartments other than DCC
  # Represented by a shifted quadratic function that increases with age and
  # prevents people younger than 10 years to progress to HCC
  cancer_prog_function <- (parameters$cancer_prog_coefficient * (ages - parameters$cancer_age_treshold))^2  # Rate in females
  cancer_prog_function <- cancer_prog_function *
    c(rep(0, times = which(ages == parameters$cancer_age_treshold-da)),
      rep(1, times = n_agecat - which(ages == parameters$cancer_age_treshold-da)))  # Set transition to 0 in <10 year olds
  cancer_prog_female <- sapply(cancer_prog_function, function(x) min(x,1)) # Set maximum annual rate is 1
  cancer_prog_male <- sapply(parameters$cancer_prog_male_cofactor*cancer_prog_female, function(x) min(x,1))  # Rate in males, cannot exceed 1
  cancer_prog_rates <- matrix(data = c(cancer_prog_female, cancer_prog_male),
                              nrow = n_agecat, ncol = 2)  # store in a matrix to apply to compartment
  parameters$cancer_prog_rates <- cancer_prog_rates

  # Age-specific progression from IR to ENCHB and IR to CC
  # Is set to 0 in under 20 year olds and a constant in over 20 year olds
  pr_ir_enchb_function <- c(rep(0, times = which(ages == 20-da)),
                            rep(parameters$pr_ir_enchb, times = n_agecat - which(ages == 20-da)))
  pr_ir_enchb_function <- sapply(pr_ir_enchb_function, function(x) min(x,5))  # annual rate cannot exceed 5
  parameters$pr_ir_enchb_function <- pr_ir_enchb_function

  pr_ir_cc_function <- c(rep(0, times = which(ages == 20-da)),
                         rep(parameters$pr_ir_cc, times = n_agecat - which(ages == 20-da)))
  pr_ir_cc_function <- sapply(pr_ir_cc_function, function(x) min(x,5))  # annual rate cannot exceed 5
  parameters$pr_ir_cc_function <- pr_ir_cc_function

  # Age-specific progression from ENCHB to CC
  # Represented by a shifted quadratic function that increases with age and
  # prevents people younger than 25 years to progress to HCC
  cirrhosis_prog_function <- parameters$pr_enchb_cc * (parameters$cirrhosis_prog_coefficient * (ages - 25))^2  # Rate in females
  cirrhosis_prog_function <- cirrhosis_prog_function *
    c(rep(0, times = which(ages == 25-da)),
      rep(1, times = n_agecat - which(ages == 25-da)))  # Set rate to 0 in <25 year olds
  cirrhosis_prog_function_female <- sapply(cirrhosis_prog_function, function(x) min(x,5)) # Set maximum annual rate to 5
  cirrhosis_prog_function_male <- sapply(parameters$cirrhosis_prog_male_cofactor *
                                           cirrhosis_prog_function_female, function(x) min(x,5))  # Rate in males, cannot exceed 1
  cirrhosis_prog_rates <- matrix(data = c(cirrhosis_prog_function_female,
                                          cirrhosis_prog_function_male),
                                 nrow = n_agecat, ncol = 2)  # store in a matrix to apply to compartment
  parameters$cirrhosis_prog_rates <- cirrhosis_prog_rates

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
  out_sf <- select(out, starts_with("Sf"))
  out_sm <- select(out, starts_with("Sm"))
  out_itf <- select(out, starts_with("ITf"))
  out_itm <- select(out, starts_with("ITm"))
  out_irf <- select(out, starts_with("IRf"))
  out_irm <- select(out, starts_with("IRm"))
  out_icf <- select(out, starts_with("ICf"))
  out_icm <- select(out, starts_with("ICm"))
  out_enchbf <- select(out, starts_with("ENCHBf"))
  out_enchbm <- select(out, starts_with("ENCHBm"))
  out_ccf <- select(out, starts_with("CCf"))
  out_ccm <- select(out, starts_with("CCm"))
  out_dccf <- select(out, starts_with("DCCf"))
  out_dccm <- select(out, starts_with("DCCm"))
  out_hccf <- select(out, starts_with("HCCf"))
  out_hccm <- select(out, starts_with("HCCm"))
  out_rf <- select(out, starts_with("Rf"))
  out_rm <- select(out, starts_with("Rm"))

  # Total population
  out_popf <- select(out, contains("f"))
  out_popm <- select(out, contains("m"))
  out_pop <- cbind(out_popf, out_popm)

  ## Extract separate outputs: incident variables (transitions between states)

  # Demographic transitions per timestep (cumulative number of births and deaths)
  out_cum_deathsf <- select(output, starts_with("cum_deathsf"))
  out_cum_deathsm <- select(output, starts_with("cum_deathsm"))
  out_cum_births <- unlist(select(output, contains("cum_births")))

  # Cumulative HBV incidence from horizontal transmission
  out_cum_infectionsf <- select(output, starts_with("cum_infectionsf"))
  out_cum_infectionsm <- select(output, starts_with("cum_infectionsm"))

  # Cumulative HBV incidence from MTCT (number of infected births)
  out_cum_infected_births <- unlist(select(output, starts_with("cum_infected_births")))

  # Cumulative chronic infection incidence from horizontal transmission
  out_cum_chronic_infectionsf <- select(output, starts_with("cum_chronic_infectionsf"))
  out_cum_chronic_infectionsm <- select(output, starts_with("cum_chronic_infectionsm"))

  # Cumulative chronic infection incidence from MTCT (number of chronically infected births)
  out_cum_chronic_births <- unlist(select(output, starts_with("cum_chronic_births")))

  # Cumulative number of HBV-related deaths (from cirrhosis and HCC)
  out_cum_hbv_deathsf <- select(output, starts_with("cum_hbv_deathsf"))
  out_cum_hbv_deathsm <- select(output, starts_with("cum_hbv_deathsm"))

  # Cumulative number of HCC cases (from all possible compartments)
  out_cum_hccf <- select(output, starts_with("cum_incident_hccf"))
  out_cum_hccm <- select(output, starts_with("cum_incident_hccm"))

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
  b1 = 0.2027,  # Margaret value in Ethiopia, for my definition of FOI use 0.07
  b2 = 0.001,   # Margaret value in Ethiopia
  b3 = 0.001,   # Margaret value in Ethiopia
  alpha = 15,         # Shevanthi value, relative infectiousness of eAg-positives
  mtct_prob_e = 0.9,  # Shevanthi value, probability of perinatal transmission from HBeAg-positive mother
  mtct_prob_s = 0.3681, # Margaret value in Ethiopia, probability of perinatal transmission from HBeAg-negative infected mother
  # NATURAL HISTORY PROGRESSION RATES AND PARAMETERS
  p_chronic_function_r = 0.65,  # Edmunds decay rate parameter in exponential function of age-specific chronic carriage risk
  p_chronic_function_s = 0.46,  # Edmunds s parameter in exponential function of age-specific chronic carriage risk
  pr_it_ir = 0.1,
  pr_ir_ic = 0.05,
  eag_prog_function_intercept = 9.5,  # 9.5 in Margaret's Ethiopia fit, 19.8873 in Shevanthi's
  eag_prog_function_rate = 0.1281,  # 0.1281 in Margaret's Ethiopia fit, 0.977 in Shevanthi's
  pr_ir_enchb = 0.005,
  pr_ir_cc = 0.028,
  pr_ic_enchb = 0.01,
  sag_loss = 0.01,  # Shevanthi value, inactive carrier to recovered transition
  #sag_loss = sagloss_rates,
  pr_enchb_cc = 0.04,  # Progression to CC (from ENCHB)
  cirrhosis_prog_coefficient = 0.0341,
  cirrhosis_prog_male_cofactor = 12.32,
  dccrate = 0.04,  # Progression to DCC (from CC)
  # PROGRESSION RATES TO HEPATOCELLULAR CARCINOMA
  cancer_prog_coefficient = 4.0452e-05,  # value from Margaret
  cancer_age_treshold = 10,  # value from Margaret
  cancer_prog_male_cofactor = 5.2075,  # value from Margaret
  hccr_it = 1,
  hccr_ir = 2,
  hccr_ic = 0.5,
  hccr_enchb = 2,
  hccr_cc = 13,
  hccr_dcc = 0.04,
  # HBV-RELATED MORTALITY RATES (MORTALITY FROM LIVER DISEASE)
  mu_cc = 0.039,
  mu_dcc = 0.314,
  mu_hcc = 0.5,
  # INFANT VACCINATION PARAMETERS
  # Vaccine coverage is read in through a function
  vacc_eff = 0.95,                           # vaccine efficacy
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

### Run the simulation: 1 SCENARIO (vacc/no_vacc) ----
# Default scenario: infant vaccine (apply_vacc = 1)

# Set all infection parms to zero:
#parameter_list <- lapply(parameter_list, FUN= function(x) x*0)
tic()
sim <- run_model(sim_duration = runtime, default_parameter_list = parameter_list,
                 parms_to_change = list(b1 = 0.1, b2 = 0.009, mtct_prob_s = 0.14),
                 scenario = "vacc")
out <- code_model_output(sim)
toc()


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
plot(outpath$time,outpath$infectioncat_total$carriers/outpath$pop_total$pop_total,type = "l", ylim = c(0,0.7))
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

# Carrier prevalence by age in 2015
plot(ages, outpath$carriers[which(outpath$time == 2015),]/
       outpath$pop[which(outpath$time == 2015),], type = "l", ylim = c(0,0.3))

# HBeAg prevalence in chronic carriers by age in 2015
plot(ages, outpath$eag_positive[which(outpath$time == 2015),]/
       outpath$carriers[which(outpath$time == 2015),], type = "l", ylim = c(0,1))

### MODEL CHECK: DEMOGRAPHY PLOTS ----
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

## Age structure plots ----

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

## Save model output
model_pop1960 <- out$full_output[221,1:(2*n_infectioncat*n_agecat)+1]
save(model_pop1960, file = here("data/simulated_inits_1960.RData"))

model_pop1880 <- out$full_output[out$time == 1880,1:(2*n_infectioncat*n_agecat)+1]
save(model_pop1880, file = here("data/simulated_inits_1880.RData"))

### Model calibration ----
require(lhs)
require("parallel")

# Function to simulate shadow models within fitting function
run_shadow_model <- function(init_age_from, init_age_to, init_sex,
                             init_compartment_from, init_compartment_to,
                             shadow_init, shadowsim_duration,
                             shadow_default_parameter_list,
                             shadow_parms_to_change) {

  # Input arguments explanation
  # The shadow model allows to follow a cohort of individuals in:
  # given compartments: init_compartment_from to init_compartment_to
  # > Compartment numbers are: 1 = Susceptible, 2 = IT, 3 = IR, 4 = IC, 5 = ENCHB,
  # > 6 = CC, 7 = DCC, 8 = HCC, 9 = Recovered
  # given age groups: init_age_from to init_age_to
  # of a given sex: init_sex takes "both", "female" or "male" as input
  # To mimic an empirical study, the cohort has a given starting year (defined in shadow_parms_to_change)
  # Define shadow_init is the population distribution simulated by the overall fitting algorithm in the starting year
  # shadowsim_duration is the length of follow-up defined in the study + 1
  # shadow_default_parameter_list should be the parameters used by overall fitting algorithm
  # shadow_parms_to_change needs to be adapted to turn off any influx into the cohort to observe
  # That involves turning off at least births and migration
  # Vaccination is automatically switched off

  # Set up cohort to follow:
  # Define index for groups of interest (infection compartments and age groups)
  if (init_sex == "both") {
    shadow_index <- c(t(mapply(seq, from = (which(ages %in% seq(init_age_from,init_age_to,da))+(init_compartment_from-1)*n_agecat),
                               to = n_agecat*(init_compartment_to), by =n_agecat)),  # for women
                      t(mapply(seq, from = (which(ages %in% seq(init_age_from,init_age_to,da))+(9+init_compartment_from-1)*n_agecat),
                               to= 2*n_agecat*(n_infectioncat-(n_infectioncat-init_compartment_to)/2), by =n_agecat)))  # for men
  } else if (init_sex == "female") {
    shadow_index <- c(t(mapply(seq, from = (which(ages %in% seq(init_age_from,init_age_to,da))+(init_compartment_from-1)*n_agecat),
                               to = n_agecat*(init_compartment_to), by =n_agecat)))
  } else if (init_sex == "male") {
    shadow_index <- c(t(mapply(seq, from = (which(ages %in% seq(init_age_from,init_age_to,da))+(9+init_compartment_from-1)*n_agecat),
                               to= 2*n_agecat*(n_infectioncat-(n_infectioncat-init_compartment_to)/2), by =n_agecat)))
  } else  {
    return(print("init_sex can be both, female or male"))
  }

  # Get init pop from full model output in a given year (defined in function call)
  shadow_init_pop <- c(shadow_init, output_storage)
  shadow_init_pop <- unlist(shadow_init_pop)

  # Set all age groups and compartments other than those to follow to nearly 0
  # (model does not seem to run if it is exactly 0)
  shadow_init_pop[-c(shadow_index,((2*n_agecat*n_infectioncat+1):length(shadow_init_pop)))] <-
    0.000000001

  # Special cases: representing a cohort of liver disease patients
  # In studies of cirrhosis, decompensated patients are overrepresented
  # To simulate these cohorts, increase DCC patients so the cohort is
  # 50% CC and 50% DCC patients
  # If the cohort includes CC, DCC and HCC, distribute as 20, 20 and 60% (Olubuyide study A6)
  if (init_compartment_from == 6 & init_compartment_to == 7) {

    cc_index <- grep("^CC", names(shadow_init_pop))  # index for CC compartments
    dcc_index <- grep("^DCC", names(shadow_init_pop))  # index for DCC compartments
    cc_to_dcc_ratio <- sum(shadow_init_pop[cc_index])/sum(shadow_init_pop[dcc_index])  # current ratio of CC to DCC patients
    shadow_init_pop[dcc_index] <- shadow_init_pop[dcc_index] * cc_to_dcc_ratio  # Increase DCC to obtain 50-50 ratio

  } else if (init_compartment_from == 6 & init_compartment_to == 8) {

    cc_index <- grep("^CC", names(shadow_init_pop))  # index for CC compartments
    dcc_index <- grep("^DCC", names(shadow_init_pop))  # index for DCC compartments
    hcc_index <- grep("^HCC", names(shadow_init_pop))  # index for HCC compartments
    cc_to_dcc_ratio <- sum(shadow_init_pop[cc_index])/sum(shadow_init_pop[dcc_index])  # current ratio of CC to DCC patients
    shadow_init_pop[dcc_index] <- shadow_init_pop[dcc_index] * cc_to_dcc_ratio  # Increase DCC to obtain 50-50 ratio
    dcc_to_hcc_ratio <- sum(shadow_init_pop[dcc_index])/sum(shadow_init_pop[hcc_index])  # current ratio of HCC to DCC patients
    shadow_init_pop[hcc_index] <- shadow_init_pop[hcc_index] * dcc_to_hcc_ratio * 3  # Increase HCC to obtain 60-40 ratio HCC to cirrhosis
  }

  # Set possible transitions within cohort to follow (no new additions) in function call

  # Simulate the cohort
  shadow_sim <- run_model(sim_duration = shadowsim_duration, init_pop_vector = shadow_init_pop,
                          default_parameter_list = shadow_default_parameter_list,
                          parms_to_change = shadow_parms_to_change,
                          scenario = "no_vacc")

  return(shadow_sim)

}

# Function to fit to age- and sex-specific seroprevalence data by mapping
# matching model output to year and age
# Use to fit to HBsAg, HBeAg and anti-HBc seroprevalence data
map_seromarker_prev <- function(seromarker_num, seromarker_denom, prev_dataset, model_output) {

  # For HBsAg, seromarker_num = carriers and seromarker_denom = pop

  seromarker_num_female <- paste0(seromarker_num, "_female")
  seromarker_num_male <- paste0(seromarker_num, "_male")
  seromarker_denom_female <- paste0(seromarker_denom, "_female")
  seromarker_denom_male <- paste0(seromarker_denom, "_male")

  # Filter the output dataset by the year of interest and calculate prevalence for all ages

  # For data from both sexes:
  model_prev_subset_both <- data.frame(time = model_output$time[model_output$time %in%
                                                                      prev_dataset$time[prev_dataset$sex == "Mixed"]],
                                       sex = "Mixed",
                                       prev = model_output[[seromarker_num]][model_output$time %in%
                                                             prev_dataset$time[prev_dataset$sex == "Mixed"],]/
                                         model_output[[seromarker_denom]][model_output$time %in% prev_dataset$time[prev_dataset$sex == "Mixed"],])


  # For women:
  model_prev_subset_female <- data.frame(time = model_output$time[model_output$time %in%
                                                                               prev_dataset$time[prev_dataset$sex == "Female"]],
                                       sex = "Female",
                                       prev = model_output[[seromarker_num_female]][model_output$time %in%
                                                                               prev_dataset$time[prev_dataset$sex == "Female"],]/
                                         model_output[[seromarker_denom_female]][model_output$time %in% prev_dataset$time[prev_dataset$sex == "Female"],])


  # For men:
  model_prev_subset_male <- data.frame(time = model_output$time[model_output$time %in%
                                               prev_dataset$time[prev_dataset$sex == "Male"]],
                                         sex = "Male",
                                         prev = model_output[[seromarker_num_male]][model_output$time %in%
                                                                                        prev_dataset$time[prev_dataset$sex == "Male"],]/
                                           model_output[[seromarker_denom_male]][model_output$time %in% prev_dataset$time[prev_dataset$sex == "Male"],])


  # Assign all columns the same names to combine into 1 dataframe
  names(model_prev_subset_both) <- c("time", "sex", paste0("prev", index$ages_all))
  names(model_prev_subset_female) <- c("time", "sex", paste0("prev", index$ages_all))
  names(model_prev_subset_male) <- c("time", "sex", paste0("prev", index$ages_all))

  # Combine sex-specific dataframes and turn into long format
  model_prev_subset <- rbind(model_prev_subset_female,
                             model_prev_subset_male,
                             model_prev_subset_both)

  model_prev_subset <- gather(model_prev_subset,
                              key = "age", value = "model_value", -time, -sex)
  model_prev_subset$age <- ages[as.numeric(gsub("\\D", "", model_prev_subset$age))]  # Assign ages as column

  # Merge with the dataset to fit to (transform factor to character vector)
  model_prev_subset$sex <- as.character(model_prev_subset$sex)

  mapped_output_seromarker <- full_join(prev_dataset,
                                        model_prev_subset,
                                        by = c("sex", "time", "age"))

  return(mapped_output_seromarker)

}

# Trial function to calculate sum of least squares for overall prevalence at given time point
fit_model_sse <- function(..., default_parameter_list, parms_to_change = list(...),
                          scenario = "vacc", data_to_fit) {

  # Simulation with given parameter set ----

  # Simulation parameters for fitting procedure:
  # Simulation starts in 1880, runs for 140 years
  # This is because data starts in 1980, so the population can reach equilibrium with
  # given parameter set
  parameters_for_fit <- generate_parameters(default_parameter_list = default_parameter_list,
                                            parms_to_change = c(parms_to_change,
                                                                sim_starttime = 1880))

  sim <- run_model(sim_duration = 140, init_pop_vector = init_pop_sim,
                   default_parameter_list = parameters_for_fit,
                   parms_to_change = NULL,
                   scenario = "vacc")
  out <- code_model_output(sim)

  # Save population distributions for shadow models
  model_pop1978 <- out$full_output[which(out$time==1978),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1983 <- out$full_output[which(out$time==1983),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1985 <- out$full_output[which(out$time==1985),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1995 <- out$full_output[which(out$time==1995),1:(2*n_infectioncat*n_agecat)+1]
  model_pop2005 <- out$full_output[which(out$time==2005),1:(2*n_infectioncat*n_agecat)+1]
  model_pop2012 <- out$full_output[which(out$time==2012),1:(2*n_infectioncat*n_agecat)+1]

  # Matching datasets to model ouput ----

  # Define my datapoints to fit to:
  data <- data_to_fit

  ## 1) GLOBOCAN age-specific HCC incidence in The Gambia in 2018 ----

  # For women:
  globocan_hcc_incidencef <- data.frame(outcome = "hcc_incidence",
                                        time = 2018,
                                        sex = "Female",
                                        age_min = unique(data_to_fit$globocan_incidence_data$age_min[data_to_fit$globocan_incidence_data$sex == "Female"]),
                                        age_max = unique(data_to_fit$globocan_incidence_data$age_max[data_to_fit$globocan_incidence_data$sex == "Female"]),
                                        model_value = 0)

  for (i in 1:nrow(globocan_hcc_incidencef)){
    globocan_hcc_incidencef$model_value[i] <-
      (sum(select(sim, starts_with("cum_incident_hccf"))[sim$time == 2019,which(ages == globocan_hcc_incidencef$age_min[i]):which(ages == globocan_hcc_incidencef$age_max[i])])-
        sum(select(sim, starts_with("cum_incident_hccf"))[sim$time == 2018,which(ages == globocan_hcc_incidencef$age_min[i]):which(ages == globocan_hcc_incidencef$age_max[i])]))/
         sum(out$pop_female[out$time == 2018.5,which(ages == globocan_hcc_incidencef$age_min[i]):which(ages == globocan_hcc_incidencef$age_max[i])])
  }


  # For men:
  globocan_hcc_incidencem <- data.frame(outcome = "hcc_incidence",
                                        time = 2018,
                                        sex = "Male",
                                        age_min = unique(data_to_fit$globocan_incidence_data$age_min[data_to_fit$globocan_incidence_data$sex == "Male"]),
                                        age_max = unique(data_to_fit$globocan_incidence_data$age_max[data_to_fit$globocan_incidence_data$sex == "Male"]),
                                        model_value = 0)

  for (i in 1 : nrow(globocan_hcc_incidencem)){
    globocan_hcc_incidencem$model_value[i] <-
      (sum(select(sim, starts_with("cum_incident_hccm"))[sim$time == 2019,which(ages == globocan_hcc_incidencem$age_min[i]):which(ages == globocan_hcc_incidencem$age_max[i])])-
        sum(select(sim, starts_with("cum_incident_hccm"))[sim$time == 2018,which(ages == globocan_hcc_incidencem$age_min[i]):which(ages == globocan_hcc_incidencem$age_max[i])]))/
          sum(out$pop_male[out$time == 2018.5,which(ages == globocan_hcc_incidencem$age_min[i]):which(ages == globocan_hcc_incidencem$age_max[i])])
  }

  # GLOBOCAN age-specific HCC mortality rate in The Gambia in 2018:

  # For women:
  globocan_hcc_mortalityf <- data.frame(outcome = "hcc_mortality",
                                        time = 2018,
                                        sex = "Female",
                                        age_min = unique(data_to_fit$globocan_incidence_data$age_min[data_to_fit$globocan_incidence_data$sex == "Female"]),
                                        age_max = unique(data_to_fit$globocan_incidence_data$age_max[data_to_fit$globocan_incidence_data$sex == "Female"]),
                                        model_value = 0)

  for (i in 1:nrow(globocan_hcc_mortalityf)){
    globocan_hcc_mortalityf$model_value[i] <-
      (sum(select(sim, starts_with("cum_hcc_deathsf"))[sim$time == 2019,which(ages == globocan_hcc_mortalityf$age_min[i]):which(ages == globocan_hcc_mortalityf$age_max[i])])-
        sum(select(sim, starts_with("cum_hcc_deathsf"))[sim$time == 2018,which(ages == globocan_hcc_mortalityf$age_min[i]):which(ages == globocan_hcc_mortalityf$age_max[i])]))/
          sum(out$pop_female[out$time == 2018.5,which(ages == globocan_hcc_mortalityf$age_min[i]):which(ages == globocan_hcc_mortalityf$age_max[i])])
  }

  # For men:
  globocan_hcc_mortalitym <- data.frame(outcome = "hcc_mortality",
                                        time = 2018,
                                        sex = "Male",
                                        age_min = unique(data_to_fit$globocan_incidence_data$age_min[data_to_fit$globocan_incidence_data$sex == "Male"]),
                                        age_max = unique(data_to_fit$globocan_incidence_data$age_max[data_to_fit$globocan_incidence_data$sex == "Male"]),
                                        model_value = 0)

  for (i in 1:nrow(globocan_hcc_mortalitym)){
    globocan_hcc_mortalitym$model_value[i] <-
      (sum(select(sim, starts_with("cum_hcc_deathsm"))[sim$time == 2019,which(ages == globocan_hcc_mortalitym$age_min[i]):which(ages == globocan_hcc_mortalitym$age_max[i])])-
         sum(select(sim, starts_with("cum_hcc_deathsm"))[sim$time == 2018,which(ages == globocan_hcc_mortalitym$age_min[i]):which(ages == globocan_hcc_mortalitym$age_max[i])]))/
      sum(out$pop_male[out$time == 2018.5,which(ages == globocan_hcc_mortalitym$age_min[i]):which(ages == globocan_hcc_mortalitym$age_max[i])])
  }

  # Combine HCC incidence and mortality sets and map to GLOBOCAN input data
  globocan_incidence_output <- rbind(globocan_hcc_incidencef, globocan_hcc_incidencem,
                                     globocan_hcc_mortalityf, globocan_hcc_mortalitym)

  # For merging, transform factors to character objects
  globocan_incidence_output$outcome <- as.character(globocan_incidence_output$outcome)
  globocan_incidence_output$sex <- as.character(globocan_incidence_output$sex)

  mapped_globocan_incidence <- left_join(data_to_fit$globocan_incidence_data,
                                         globocan_incidence_output,
                                         by = c("outcome", "time", "sex", "age_min", "age_max"))

  # MULTIPLY DATA VALUE BY PAF to get HBV-related HCC only
  paf_under50 <- 0.83
  paf_over50 <- 0.32
  mapped_globocan_incidence$data_value[mapped_globocan_incidence$age_max < 55] <-
    mapped_globocan_incidence$data_value[mapped_globocan_incidence$age_max < 55]*paf_under50
  mapped_globocan_incidence$data_value[mapped_globocan_incidence$age_min >= 55] <-
    mapped_globocan_incidence$data_value[mapped_globocan_incidence$age_min >= 55]*paf_over50

  ## 2) RISK OF CHRONIC CARRIAGE BY AGE AT INFECTION ----
  model_p_chronic <- data.frame(outcome = "p_chronic",
                                age = ages,
                                model_value = unlist(head(select(sim, starts_with("p_chronic_function")),1)))
  model_p_chronic$outcome <- as.character(model_p_chronic$outcome)

  mapped_p_chronic <- full_join(data_to_fit$p_chronic,
                                model_p_chronic, by = c("outcome", "age"))

  ## 3) AGE- AND SEX- SPECIFIC SEROMARKER PREVALENCE ----
  # Mapping matching model output to specific year and age

  # HBsAg prevalence in the population:
  mapped_output_hbsag <- map_seromarker_prev(seromarker_num = "carriers",
                                             seromarker_denom = "pop",
                                             prev_dataset = data_to_fit$hbsag_prevalence,
                                             model_output = out)
  mapped_output_hbsag$outcome <- "HBsAg_prevalence"

  # HBeAg prevalence in chronic carriers:
  mapped_output_hbeag <- map_seromarker_prev(seromarker_num = "eag_positive",
                                              seromarker_denom = "carriers",
                                              prev_dataset = data_to_fit$hbeag_prevalence,
                                              model_output = out)
  mapped_output_hbeag$outcome <- "HBeAg_prevalence"

  # anti-HBc prevalence in the population:
  mapped_output_antihbc <- map_seromarker_prev(seromarker_num = "ever_infected",
                                             seromarker_denom = "pop",
                                             prev_dataset = data_to_fit$antihbc_prevalence,
                                             model_output = out)
  mapped_output_antihbc$outcome <- "Anti_HBc_prevalence"

  # Combine into a single dataframe
  mapped_seromarker_prevalence <- rbind(mapped_output_hbsag, mapped_output_hbeag, mapped_output_antihbc)

  ## 4) NATURAL HISTORY-RELATED PREVALENCE ----
  # Various estimates in chronic carriers/liver disease subgroups
  # Outputs calculated manually and matched to specific agegroup and year

  # Prepare output dataframe:
  model_output_nat_hist <- data.frame(id_unique = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                       age_min = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                       age_max = rep(0, nrow(data_to_fit$natural_history_prevalence)),
                                       model_value = rep(0, nrow(data_to_fit$natural_history_prevalence)))

  # Prepare denominators/numerators in common to several studies:

  # Male carriers, aged 27-35.5 years, in 2013
  denom_gmb1_2_2013 <- (sum(out$carriers_male[which(out$time == 2013),(which(ages ==27):which(ages ==35.5))]))

  # Carriers without liver disease, aged 4.5 to 21.5 years, in 1986
  denom_1_1_1986 <- sum(select(sim, starts_with("ITf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          select(sim, starts_with("IRf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          select(sim, starts_with("ICf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          select(sim, starts_with("ENCHBf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
                          select(sim, starts_with("ITm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          select(sim, starts_with("IRm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          select(sim, starts_with("ICm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
                          select(sim, starts_with("ENCHBm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])

  # Carriers, aged 33 to 47 years, in 2012
  denom_gmb1_1_2012 <- (sum(out$carriers[which(out$time == 2012),(which(ages ==33):which(ages ==47))]))

  # Carriers without liver disease, aged 8 to 95.5 years, in 2013
  denom_1_1_2013 <- sum(select(sim, starts_with("ITf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("IRf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("ICf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("ENCHBf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
                          select(sim, starts_with("ITm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("IRm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("ICm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("ENCHBm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])

  # Carriers with significant liver fibrosis or cirrhosis in 2013
  num_1_1_2013 <- select(sim, starts_with("IRf"))[which(sim$time == 2013),]+
    select(sim, starts_with("ENCHBf"))[which(sim$time == 2013),]+
    select(sim, starts_with("CCf"))[which(sim$time == 2013),]+
    select(sim, starts_with("DCCf"))[which(sim$time == 2013),]+
    select(sim, starts_with("IRm"))[which(sim$time == 2013),]+
    select(sim, starts_with("ENCHBm"))[which(sim$time == 2013),]+
    select(sim, starts_with("CCm"))[which(sim$time == 2013),]+
    select(sim, starts_with("DCCm"))[which(sim$time == 2013),]

  # Incident CC cases in 1999
  denom_gmb15_2_1999 <-
    select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1999),]+
    select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1999),]+
    select(sim, starts_with("cum_enchb_to_ccf"))[which(sim$time == 1999),]+
    select(sim, starts_with("cum_enchb_to_ccm"))[which(sim$time == 1999),]-
    select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1998),]-
    select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1998),]-
    select(sim, starts_with("cum_enchb_to_ccf"))[which(sim$time == 1998),]-
    select(sim, starts_with("cum_enchb_to_ccm"))[which(sim$time == 1998),]

  # Incident HBeAg-positive non-cirrhotic HCC cases in 1990
  # (same intermediate timepoint used for 2 studies)
  num_gmb12_gmb15_1990 <-
    select(sim, starts_with("cum_it_to_hccf"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_it_to_hccm"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_ir_to_hccf"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_ir_to_hccm"))[which(sim$time == 1991),]-
    select(sim, starts_with("cum_it_to_hccf"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_it_to_hccm"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_ir_to_hccf"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_ir_to_hccm"))[which(sim$time == 1990),]

  # Incident non-cirrhotic HCC cases in 1990
  # (same intermediate timepoint used for 2 studies)
  denom_gmb12_gmb15_1990 <-
    select(sim, starts_with("cum_it_to_hccf"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_it_to_hccm"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_ir_to_hccf"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_ir_to_hccm"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_ic_to_hccf"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_ic_to_hccm"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_enchb_to_hccf"))[which(sim$time == 1991),]+
    select(sim, starts_with("cum_enchb_to_hccm"))[which(sim$time == 1991),] -
    select(sim, starts_with("cum_it_to_hccf"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_it_to_hccm"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_ir_to_hccf"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_ir_to_hccm"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_ic_to_hccf"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_ic_to_hccm"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_enchb_to_hccf"))[which(sim$time == 1990),]-
    select(sim, starts_with("cum_enchb_to_hccm"))[which(sim$time == 1990),]

  # Calculate each datapoint manually and distinguish by minimum and maximum age and
  # a unique ID as follows:
  # id_*study ID*_*group ID*_*datapoint time*_*numerator*

  model_output_nat_hist[1,] <-
    c("id_gmb1_2_2013_ic",
      age_min = 27,
      age_max = 35,
      (sum(select(sim, starts_with("ICm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[2,] <-
    c("id_gmb1_2_2013_it_ic",
      age_min = 27,
      age_max = 35,
      (sum(select(sim, starts_with("ITm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
         sum(select(sim, starts_with("ICm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[3,] <-
    c("id_gmb1_2_2013_ir_enchb",
      age_min = 27,
      age_max = 35,
      (sum(select(sim, starts_with("IRm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
         sum(select(sim, starts_with("ENCHBm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[4,] <-
    c("id_gmb1_2_2013_cc_dcc",
      age_min = 27,
      age_max = 35,
      (sum(select(sim, starts_with("CCm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
         sum(select(sim, starts_with("DCCm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[5,] <-
    c("id_gmb1_2_2013_hcc",
      age_min = 27,
      age_max = 35,
      (sum(select(sim, starts_with("HCCm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
        denom_gmb1_2_2013)

  model_output_nat_hist[6,] <-
    c("id_1_1_1986_it",
      age_min = 4.5,
      age_max = 21.5,
      sum(select(sim, starts_with("ITf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            select(sim, starts_with("ITm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[7,] <-
    c("id_1_1_1986_ir",
      age_min = 4.5,
      age_max = 21.5,
      sum(select(sim, starts_with("IRf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            select(sim, starts_with("IRm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[8,] <-
    c("id_1_1_1986_enchb",
      age_min = 4.5,
      age_max = 21.5,
      sum(select(sim, starts_with("ENCHBf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            select(sim, starts_with("ENCHBm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[9,] <-
    c("id_1_1_1986_ic",
      age_min = 4.5,
      age_max = 21.5,
      sum(select(sim, starts_with("ICf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            select(sim, starts_with("ICm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        denom_1_1_1986)

  model_output_nat_hist[10,] <-
    c("id_1_1_1986_hcc",
      age_min = 4.5,
      age_max = 21.5,
      sum(select(sim, starts_with("HCCf"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]+
            select(sim, starts_with("HCCm"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
        (sum(out$carriers[which(out$time == 1986),(which(ages ==4.5):which(ages ==21.5))])))

  model_output_nat_hist[11,] <-
    c("id_gmb1_1_2012_ic",
      age_min = 33,
      age_max = 47,
      (sum(select(sim, starts_with("ICf"))[which(sim$time == 2012),which(ages ==33):which(ages ==47)]+
             select(sim, starts_with("ICm"))[which(sim$time == 2012),which(ages ==33):which(ages ==47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[12,] <-
    c("id_gmb1_1_2012_it_ic",
      age_min = 33,
      age_max = 47,
      (sum(select(sim, starts_with("ICf"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
         select(sim, starts_with("ITf"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
         select(sim, starts_with("ICm"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
         select(sim, starts_with("ITm"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[13,] <-
    c("id_gmb1_1_2012_ir_enchb",
      age_min = 33,
      age_max = 47,
      (sum(select(sim, starts_with("IRf"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
         select(sim, starts_with("ENCHBf"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
           select(sim, starts_with("IRm"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
           select(sim, starts_with("ENCHBm"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[14,] <-
    c("id_gmb1_1_2012_cc_dcc",
      age_min = 33,
      age_max = 47,
      (sum(select(sim, starts_with("CCf"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             select(sim, starts_with("DCCf"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             select(sim, starts_with("CCm"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             select(sim, starts_with("DCCm"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[15,] <-
    c("id_gmb1_1_2012_hcc",
      age_min = 33,
      age_max = 47,
      (sum(select(sim, starts_with("HCCf"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]+
             select(sim, starts_with("HCCm"))[which(sim$time == 2012),which(ages == 33):which(ages == 47)]))/
        denom_gmb1_1_2012)

  model_output_nat_hist[16,] <-
    c("id_1_1_2013_it",
      age_min = 8,
      age_max = 95,
      sum(select(sim, starts_with("ITf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            select(sim, starts_with("ITm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[17,] <-
    c("id_1_1_2013_ir",
      age_min = 8,
      age_max = 95,
      sum(select(sim, starts_with("IRf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            select(sim, starts_with("IRm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[18,] <-
    c("id_1_1_2013_enchb",
      age_min = 8,
      age_max = 95,
      sum(select(sim, starts_with("ENCHBf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            select(sim, starts_with("ENCHBm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[19,] <-
    c("id_1_1_2013_ic",
      age_min = 8,
      age_max = 95,
      sum(select(sim, starts_with("ICf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
            select(sim, starts_with("ICm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
        denom_1_1_2013)

  model_output_nat_hist[20,] <-
    c("id_1_1_2013_cc_dcc",
      age_min = 8,
      age_max = 95,
      (sum(select(sim, starts_with("CCf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
         select(sim, starts_with("DCCf"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
         select(sim, starts_with("CCm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]+
         select(sim, starts_with("DCCm"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==8):which(ages ==95.5))])))

  model_output_nat_hist[21,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 8,
      age_max = 29,
      (sum(num_1_1_2013[,which(ages ==8):which(ages == 29.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==8):which(ages ==29.5))])))

  model_output_nat_hist[22,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 30,
      age_max = 39,
      (sum(num_1_1_2013[,which(ages ==30):which(ages == 39.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==30):which(ages ==39.5))])))

  model_output_nat_hist[23,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 40,
      age_max = 49,
      (sum(num_1_1_2013[,which(ages ==40):which(ages == 49.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==40):which(ages ==49.5))])))

  model_output_nat_hist[24,] <-
    c("id_1_1_2013_ir_enchb_cc_dcc",
      age_min = 50,
      age_max = 95,
      (sum(num_1_1_2013[,which(ages ==50):which(ages == 95.5)]))/
        (sum(out$carriers[which(out$time == 2013),(which(ages ==50):which(ages ==95.5))])))

  # CC prevalence in HCC
  # Approximate as proportion of incident HCC cases in 1999 originating from CC
  model_output_nat_hist[25,] <-
    c("id_gmb2_1_1999_incident_hcc_cases_from_cc",
      age_min = 15,
      age_max = 67,
      (sum(select(sim, starts_with("cum_cc_to_hcc"))[which(sim$time == 1999),])-
         sum(select(sim, starts_with("cum_cc_to_hcc"))[which(sim$time == 1998),]))/
        (sum(select(sim, starts_with("cum_incident_hcc"))[which(sim$time == 1999),])-
           sum(select(sim, starts_with("cum_incident_hcc"))[which(sim$time == 1998),])))

  # DCC prevalence in HCC, 1999
  # Approximate as proportion of incident HCC cases in 1999 originating from DCC
  model_output_nat_hist[26,] <-
    c("id_gmb2_1_1999_incident_hcc_cases_from_dcc",
      age_min = 15,
      age_max = 67,
      (sum(select(sim, starts_with("cum_dcc_to_hcc"))[which(sim$time == 1999),])-
         sum(select(sim, starts_with("cum_dcc_to_hcc"))[which(sim$time == 1998),]))/
        (sum(select(sim, starts_with("cum_incident_hcc"))[which(sim$time == 1999),])-
           sum(select(sim, starts_with("cum_incident_hcc"))[which(sim$time == 1998),])))

  # HBeAg prevalence in cirrhosis patients, 1999
  # Approximate as proportion of incident CC cases in 1999 originating from IR
  model_output_nat_hist[27,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 17,
      age_max = 34,
      (sum(select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1999),which(ages == 17):which(ages == 34.5)])+
         sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1999),which(ages == 17):which(ages == 34.5)])-
            sum(select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1998),which(ages == 17):which(ages == 34.5)])-
               sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1998),which(ages == 17):which(ages == 34.5)]))/
               sum(denom_gmb15_2_1999[which(ages == 17):which(ages == 34.5)]))

  model_output_nat_hist[28,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 35,
      age_max = 44,
      (sum(select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1999),which(ages == 35):which(ages == 44.5)])+
         sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1999),which(ages == 35):which(ages == 44.5)])-
         sum(select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1998),which(ages == 35):which(ages == 44.5)])-
         sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1998),which(ages == 35):which(ages == 44.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 35):which(ages == 44.5)]))

  model_output_nat_hist[29,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 45,
      age_max = 54,
      (sum(select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1999),which(ages == 45):which(ages == 54.5)])+
         sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1999),which(ages == 45):which(ages == 54.5)])-
         sum(select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1998),which(ages == 45):which(ages == 54.5)])-
         sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1998),which(ages == 45):which(ages == 54.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 45):which(ages == 54.5)]))

  model_output_nat_hist[30,] <-
    c("id_gmb15_2_1999_hbeag_cirrhosis",
      age_min = 55,
      age_max = 64,
      (sum(select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1999),which(ages == 55):which(ages == 64.5)])+
         sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1999),which(ages == 55):which(ages == 64.5)])-
         sum(select(sim, starts_with("cum_ir_to_ccf"))[which(sim$time == 1998),which(ages == 55):which(ages == 64.5)])-
         sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1998),which(ages == 55):which(ages == 64.5)]))/
        sum(denom_gmb15_2_1999[which(ages == 55):which(ages == 64.5)]))

  # HBeAg prevalence in HCC patients, 1982
  # Approximate as proportion of non-cirrhotic incident HCC cases in 1990
  # originating from IT and IR
  model_output_nat_hist[31,] <-
    c("id_gmb12_1_1982_hbeag_hcc",
      age_min = 15,
      age_max = 49,
      sum(num_gmb12_gmb15_1990[which(ages == 15):which(ages==49.5)])/
       sum(denom_gmb12_gmb15_1990[which(ages == 15):which(ages == 49.5)]))

  model_output_nat_hist[32,] <-
    c("id_gmb12_1_1982_hbeag_hcc",
      age_min = 50,
      age_max = 72,
      sum(num_gmb12_gmb15_1990[which(ages == 50):which(ages==72.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 50):which(ages == 72.5)]))

  model_output_nat_hist[33,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 17,
      age_max = 34,
      sum(num_gmb12_gmb15_1990[which(ages == 17):which(ages == 34.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 17):which(ages == 34.5)]))

  model_output_nat_hist[34,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 35,
      age_max = 44,
      sum(num_gmb12_gmb15_1990[which(ages == 35):which(ages == 44.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 35):which(ages == 44.5)]))

  model_output_nat_hist[35,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 45,
      age_max = 54,
      sum(num_gmb12_gmb15_1990[which(ages == 45):which(ages == 54.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 45):which(ages == 54.5)]))

  model_output_nat_hist[36,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 55,
      age_max = 64,
      sum(num_gmb12_gmb15_1990[which(ages == 55):which(ages == 64.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 55):which(ages == 64.5)]))

  model_output_nat_hist[37,] <-
    c("id_gmb15_1_1999_hbeag_hcc",
      age_min = 65,
      age_max = 87,
      sum(num_gmb12_gmb15_1990[which(ages == 65):which(ages == 87.5)])/
        sum(denom_gmb12_gmb15_1990[which(ages == 65):which(ages == 87.5)]))

  # Proportion of chronic carriers attributable to mother-to-child transmission
  model_output_nat_hist[38,] <-
    c("id_1_1_1986_incident_chronic_births",
      age_min = 0,
      age_max = 99.5,
      sum(select(sim, starts_with("cum_chronic_births"))[which(sim$time == 1985),])/
        (sum(select(sim, starts_with("cum_chronic_births"))[which(sim$time == 1985),])+
           sum(select(sim, starts_with("cum_chronic_infections"))[which(sim$time == 1985),])))

  # Turn age into numeric format:
  model_output_nat_hist$age_min <- as.numeric(model_output_nat_hist$age_min)
  model_output_nat_hist$age_max <- as.numeric(model_output_nat_hist$age_max)
  model_output_nat_hist$model_value <- as.numeric(model_output_nat_hist$model_value)

  # Merge with the dataset to fit to:
  mapped_nat_hist_prevalence <- left_join(data_to_fit$natural_history_prevalence,
                                      model_output_nat_hist,
                                      by = c("id_unique", "age_min", "age_max"))

  ## 5) MOTHER-TO-CHILD TRANSMISSION RISK ----
  # Overall mother-to-child transmission risk in given year (irrespective of maternal HBeAg)
  mtct_risk <- data.frame(time = out$time[out$time %in% data_to_fit$mtct_risk$time],
                          model_value = rowSums(out$eag_positive_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba] * parameters_for_fit$mtct_prob_e +
                                                (out$carriers_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]-
                                                   out$eag_positive_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]) * parameters_for_fit$mtct_prob_s)/
                            rowSums(out$carriers_female[out$time %in% data_to_fit$mtct_risk$time,index$ages_wocba]))

  mapped_mtct_risk <- left_join(data_to_fit$mtct_risk, mtct_risk, by = "time")

  ## 6) NATURAL HISTORY PROGRESSION RATES AND MORTALITY CURVES ----
  # Some are calculated from full model output and others require a shadow model

  # Prepare output storage for progression rates and mortality curves:
  progression_rates <- data.frame(outcome = c("shadow1a_eag_loss_m",
                                              "shadow1b_eag_loss_m",
                                              "shadow1a_eag_loss_f",
                                              "shadow1b_eag_loss_f",
                                              "shadow1a_hcc_incidence_m",
                                              "shadow1b_hcc_incidence_m",
                                              "shadow1a_hcc_incidence_f",
                                              "shadow1b_hcc_incidence_f",
                                              "shadow1_dcc_incidence",
                                              "shadow1_mortality_m",
                                              "shadow1_mortality_f",
                                              "shadow3_mortality",
                                              "shadow2_sag_loss",
                                              "gmb6_1_a_foi",
                                              "gmb6_1_b_foi",
                                              "gmb7_1_chronic_infection_incidence"),
                                  model_value = 0,
                                  stringsAsFactors = FALSE)

  model_mort_curves <- data.frame(outcome = data_to_fit$mortality_curves$outcome,
                                  model_value = 0,
                                  stringsAsFactors = FALSE)

  ## 6a) TRANSMISSION-REALTED RATES

  # Force of infection (any infection) in 0.5-8.5 year olds (GMB6)
  # Numerator = cumulative incidence of chronic infections and transitions to immune compartment (= any horizontal infection) over follow-up
  # Denominator = person-time in susceptible compartment
  progression_rates[progression_rates$outcome=="gmb6_1_a_foi","model_value"] <-
    (sum(select(sim, starts_with("cum_infectionsf"))[which(sim$time == 1984),(which(ages == 0.5):which(ages == 8.5))]+
                     select(sim, starts_with("cum_infectionsm"))[which(sim$time == 1984),(which(ages == 0.5):which(ages == 8.5))]) -
                  sum(select(sim, starts_with("cum_infectionsf"))[which(sim$time == 1980),(which(ages == 0.5):which(ages == 8.5))]+
                        select(sim, starts_with("cum_infectionsm"))[which(sim$time == 1980),(which(ages == 0.5):which(ages == 8.5))]))/
    ((sum(select(sim, starts_with("Sf"))[(which(sim$time == 1980):which(sim$time == 1983.5)),(which(ages == 0.5):which(ages == 8.5))]+
            select(sim, starts_with("Sm"))[(which(sim$time == 1980):which(sim$time == 1983.5)),(which(ages == 0.5):which(ages == 8.5))]))*dt)
  progression_rates[progression_rates$outcome=="gmb6_1_b_foi","model_value"] <-
    progression_rates[progression_rates$outcome=="gmb6_1_a_foi","model_value"]

  # Incidence rate of chronic infections in 0.5-7.5 year olds (GMB7)
  # Numerator = cumulative incidence of chronic infections over follow-up
  # Denominator = person-time in susceptible compartment
  progression_rates[progression_rates$outcome=="gmb7_1_chronic_infection_incidence","model_value"] <-
    (sum(select(sim, starts_with("cum_chronic_infectionsf"))[which(sim$time == 1982),(which(ages == 0.5):which(ages == 7.5))]+
                                   select(sim, starts_with("cum_chronic_infectionsm"))[which(sim$time == 1982),(which(ages == 0.5):which(ages == 7.5))]) -
                               sum(select(sim, starts_with("cum_chronic_infectionsf"))[which(sim$time == 1981),(which(ages == 0.5):which(ages == 7.5))]+
                                     select(sim, starts_with("cum_chronic_infectionsm"))[which(sim$time == 1981),(which(ages == 0.5):which(ages == 7.5))]))/
    ((sum(select(sim, starts_with("Sf"))[(which(sim$time == 1981):which(sim$time == 1981.5)),(which(ages == 0.5):which(ages == 7.5))]+
            select(sim, starts_with("Sm"))[(which(sim$time == 1981):which(sim$time == 1981.5)),(which(ages == 0.5):which(ages == 7.5))]))*dt)

  ## 6b) SHADOW MODELS 1 a and b: SHIMAKAWA NATURAL HISTORY COHORT

  # Follow 2 cohorts of chronic carriers - 1 of 0-19 year olds (1a) and
  # one of 20-29 year olds (1b) for 28 years, starting in 1985
  # no one had HCC at baseline and we assume no one had DCC at baseline
  # Compartments: IT, IR, IC, ENCHB, CC
  # Switch off births, migation, betas and vaccination to prevent influx of new carriers
  shadow1a_sim <- run_shadow_model(init_age_from = 0, init_age_to = 19.5, init_sex = "both",
                                   init_compartment_from = 2, init_compartment_to = 6,
                                   shadow_default_parameter_list = parameters_for_fit,
                                   shadow_init = model_pop1985, shadowsim_duration = 29,
                                   shadow_parms_to_change = list(sim_starttime = 1985,
                                                                 births_on = 0,
                                                                 migration_on = 0,
                                                                 b1 = 0,
                                                                 b2 = 0,
                                                                 b3 = 0))
  # Total initial population size in shadow model 1a
  shadow1a_init_pop <- sum(shadow1a_sim[1,1:(2*n_infectioncat*n_agecat)])

  shadow1b_sim <- run_shadow_model(init_age_from = 20, init_age_to = 29.5, init_sex = "both",
                                   init_compartment_from = 2, init_compartment_to = 6,
                                   shadow_default_parameter_list = parameters_for_fit,
                                   shadow_init = model_pop1985, shadowsim_duration = 29,
                                   shadow_parms_to_change = list(sim_starttime = 1985,
                                                                 births_on = 0,
                                                                 migration_on = 0,
                                                                 b1 = 0,
                                                                 b2 = 0,
                                                                 b3 = 0))
  # Total initial population size in shadow model 1b
  shadow1b_init_pop <- sum(shadow1b_sim[1,1:(2*n_infectioncat*n_agecat)])

  # Calculate outputs to fit to:
  shadow1a_out <- code_model_output(shadow1a_sim)
  shadow1b_out <- code_model_output(shadow1b_sim)

  # Total HCC incidence rate per person-year over follow-up
  # Numerator = cumulative number of incident HCC cases over follow-up
  # Denominator = person-timestep at risk * dt (carrier compartments other than HCC)
  # MODEL 1a
  # In women:
  progression_rates[progression_rates$outcome=="shadow1a_hcc_incidence_f","model_value"] <-
    sum(select(tail(shadow1a_sim,1), starts_with("cum_incident_hccf")))/
    ((sum(head(shadow1a_out$carriers_female,-1)) - sum(head(select(shadow1a_out$full_output,
                                                                        starts_with("HCCf")),-1)))*dt)

  # In men:
  progression_rates[progression_rates$outcome=="shadow1a_hcc_incidence_m","model_value"] <-
    sum(select(tail(shadow1a_sim,1), starts_with("cum_incident_hccm")))/
    ((sum(head(shadow1a_out$carriers_male,-1)) - sum(head(select(shadow1a_out$full_output,
                                                                        starts_with("HCCm")),-1)))*dt)

  # MODEL 1b
  # In women:
  progression_rates[progression_rates$outcome=="shadow1b_hcc_incidence_f","model_value"] <-
    sum(select(tail(shadow1b_sim,1), starts_with("cum_incident_hccf")))/
    ((sum(head(shadow1b_out$carriers_female,-1)) - sum(head(select(shadow1b_out$full_output,
                                                                        starts_with("HCCf")),-1)))*dt)

  # In men:
  progression_rates[progression_rates$outcome=="shadow1b_hcc_incidence_m","model_value"] <-
    sum(select(tail(shadow1b_sim,1), starts_with("cum_incident_hccm")))/
    ((sum(head(shadow1b_out$carriers_male,-1)) - sum(head(select(shadow1b_out$full_output,
                                                                      starts_with("HCCm")),-1)))*dt)


  # Total incidence of non-malignant ESLD (DCC) per person-year (both sexes)
  # Numerator = cumulative number of incident DCC cases over follow-up (at last timestep)
  # - cumulative number of transitions from DCC to HCC
  # Denominator = person-timestep at risk * dt (carrier compartments other than DCC and HCC)
  # MODEL 1a:
  shadow1a_dcc_rate <- (sum(select(tail(shadow1a_sim,1), starts_with("cum_incident_dcc")))-
                          sum(select(tail(shadow1a_sim,1), starts_with("cum_dcc_to_hcc"))))/
    (sum(head(shadow1a_out$carriers,-1)) -
       (sum(select(head(shadow1a_sim,-1), starts_with("HCC")))) -
       (sum(select(head(shadow1a_sim,-1),starts_with("DCC"))))*dt)

  # MODEL 1b:
  shadow1b_dcc_rate <- (sum(select(tail(shadow1b_sim,1), starts_with("cum_incident_dcc")))-
                          sum(select(tail(shadow1b_sim,1), starts_with("cum_dcc_to_hcc"))))/
    (sum(head(shadow1b_out$carriers,-1)) -
       (sum(select(head(shadow1b_sim,-1), starts_with("HCC")))) -
       (sum(select(head(shadow1b_sim,-1),starts_with("DCC"))))*dt)
  # USE AVERAGE ACROSS AGE GROUPS:
  progression_rates[progression_rates$outcome=="shadow1_dcc_incidence","model_value"] <-
    weighted.mean(x = c(shadow1a_dcc_rate, shadow1b_dcc_rate),
                      w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  # Incidence rate of eAg loss per person-year
  # Numerator = cumulative number of cases of eAg loss at last timestep
  # Denominator = person-time spent in IT and IR compartments
  # MODEL 1a
  # In women:
  progression_rates[progression_rates$outcome=="shadow1a_eag_loss_f","model_value"] <-
    (sum(select(tail(shadow1a_sim,1), starts_with("cum_eag_lossf"))))/
    ((sum(select(head(shadow1a_sim,-1),starts_with("ITf")))+
        sum(select(head(shadow1a_sim,-1),starts_with("IRf"))))*dt)
  # In men:
  progression_rates[progression_rates$outcome=="shadow1a_eag_loss_m","model_value"] <-
    (sum(select(tail(shadow1a_sim,1), starts_with("cum_eag_lossm"))))/
    ((sum(select(head(shadow1a_sim,-1),starts_with("ITm")))+
        sum(select(head(shadow1a_sim,-1),starts_with("IRm"))))*dt)

  # MODEL 1b
  # In women:
  progression_rates[progression_rates$outcome=="shadow1b_eag_loss_f","model_value"] <-
    (sum(select(tail(shadow1b_sim,1), starts_with("cum_eag_lossf"))))/
    ((sum(select(head(shadow1b_sim,-1),starts_with("ITf")))+
        sum(select(head(shadow1b_sim,-1),starts_with("IRf"))))*dt)
  # In men:
  progression_rates[progression_rates$outcome=="shadow1b_eag_loss_m","model_value"] <-
    (sum(select(tail(shadow1b_sim,1), starts_with("cum_eag_lossm"))))/
    ((sum(select(head(shadow1b_sim,-1),starts_with("ITm")))+
        sum(select(head(shadow1b_sim,-1),starts_with("IRm"))))*dt)
  # Mortality rate from any cause (HBV-related deaths + background mortality)
  # This was measured in the whole cohort (no matter where they progressed to)
  # MODEL 1a
  # In women:
  shadow1a_mortality_ratef <- (sum(select(tail(shadow1a_sim,1), starts_with("cum_hbv_deathsf"))) +
                                 sum(select(tail(shadow1a_sim,1), starts_with("cum_deathsf"))))/
    (sum(head(shadow1a_out$pop_female,-1))*dt)
  # In men:
  shadow1a_mortality_ratem <- (sum(select(tail(shadow1a_sim,1), starts_with("cum_hbv_deathsm"))) +
                                 sum(select(tail(shadow1a_sim,1), starts_with("cum_deathsm"))))/
    (sum(head(shadow1a_out$pop_male,-1))*dt)
  # MODEL 1b
  # In women:
  shadow1b_mortality_ratef <- (sum(select(tail(shadow1b_sim,1), starts_with("cum_hbv_deathsf"))) +
                                 sum(select(tail(shadow1b_sim,1), starts_with("cum_deathsf"))))/
    (sum(head(shadow1b_out$pop_female,-1))*dt)
  # In men:
  shadow1b_mortality_ratem <- (sum(select(tail(shadow1b_sim,1), starts_with("cum_hbv_deathsm"))) +
                                 sum(select(tail(shadow1b_sim,1), starts_with("cum_deathsm"))))/
    (sum(head(shadow1b_out$pop_male,-1))*dt)

  # USE AVERAGE ACROSS AGE GROUPS
  progression_rates[progression_rates$outcome=="shadow1_mortality_f","model_value"] <-
    weighted.mean(x = c(shadow1a_mortality_ratef, shadow1b_mortality_ratef),
                                           w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))
  progression_rates[progression_rates$outcome=="shadow1_mortality_m","model_value"] <-
    weighted.mean(x = c(shadow1a_mortality_ratem, shadow1b_mortality_ratem),
                                           w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  ## 6c) SHADOW MODEL 2: COURSAGET CHRONIC CARRIER COHORT

  # Follow a cohort of chronic carriers aged 0-2 years
  # starting in 1978, for 7 years
  # Assume no one had HCC or DCC at baseline based on their age
  # compartments: IT, IR, IC, ENCHB, CC
  # Only represent progression within chronic carriers, so switch off births,
  # migation, betas and vaccination
  shadow2_sim <- run_shadow_model(init_age_from = 0, init_age_to = 2.5, init_sex = "both",
                           init_compartment_from = 2, init_compartment_to = 6,
                           shadow_default_parameter_list = parameters_for_fit,
                           shadow_init = model_pop1978, shadowsim_duration = 8,
                           shadow_parms_to_change = list(sim_starttime = 1978,
                                                         births_on = 0,
                                                         migration_on = 0,
                                                         b1 = 0,
                                                         b2 = 0,
                                                         b3 = 0))

  # Calculate output for fitting:
  # Overall rate of sAg loss
  # Numerator = cumulative number of incident transitions from IC to immune
  # Denominator = person-time in IC compartment at risk
  progression_rates[progression_rates$outcome=="shadow2_sag_loss","model_value"] <-
    (sum(select(tail(shadow2_sim,1), starts_with("cum_sag_loss"))))/
    (sum(select(head(shadow2_sim,-1),starts_with("IC")))*dt)

  ## 6d) SHADOW MODEL 3: OLUBUYIDE

  # Olubuyide (A6) cohort of CC, DCC and HCC patients, followed from 1983 for 6 years
  # Scaling initial population to achieve distribution: 60% HCC, 20% CC, 20% DCC
  shadow3_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1983, shadowsim_duration = 7,
                                  shadow_parms_to_change = list(sim_starttime = 1983,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Output: Mortality rate in CC, DCC and HCC from any cause (HBV-related or background mortality)
  # Numerator = sum of cumulative number of incident deaths between 1989 and 1983
  # Denominator = person-time in compartments at risk (compensated cirrhosis, decompensated cirrhosis, HCC) * dt
  progression_rates[progression_rates$outcome=="shadow3_mortality","model_value"] <-
    (sum(select(tail(shadow3_sim,1), starts_with("cum_hbv_deaths"))) +
                               sum(select(tail(shadow3_sim,1), starts_with("cum_background_deaths_ld"))))/
    ((sum(select(head(shadow3_sim,-1), starts_with("CC"))) + sum(select(head(shadow3_sim,-1), starts_with("DCC"))) +
       sum(select(head(shadow3_sim,-1), starts_with("HCC"))))*dt)

  ## 6e) SHADOW MODEL 4: SHIMAKAWA COMPENSATED CIRRHOSIS COHORT

  # To fit survival curve at time interval 0.5 years
  # Follow a cohort of compensated cirrhosis patients
  # starting in 2012, for 0.5 years
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow4_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 6,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop2012, shadowsim_duration = 1,
                                  shadow_parms_to_change = list(sim_starttime = 2012,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from HBV-related cause or background) after 0.5 years
  # Numerator = sum of incident deaths at next timestep
  # Can use background deaths from all LD patients because DCC and HCC compartments are empty
  # Denominator = number in compensated cirrhosis compartment at first timestep
  model_mort_curves[model_mort_curves$outcome=="shadow4_cum_mortality","model_value"] <-
    (sum(select(shadow4_sim, starts_with("cum_hbv_deaths"))[which(shadow4_sim$time == 2012.5),]) +
      sum(select(shadow4_sim, starts_with("cum_background_deaths_ld"))[which(shadow4_sim$time == 2012.5),]))/
    (sum(select(shadow4_sim, starts_with("CC"))[which(shadow4_sim$time == 2012),]))

  # Proportion of deaths due to DCC and HCC

  mapped_nat_hist_prevalence$model_value[
    mapped_nat_hist_prevalence$id_unique == "id_a4_1_2014_shadow_incident_deaths"] <-
    (sum(select(shadow4_sim, starts_with("cum_hcc_deaths"))[which(shadow4_sim$time == 2012.5),]) +
      sum(select(shadow4_sim, starts_with("cum_dcc_deaths"))[which(shadow4_sim$time == 2012.5),]))/
      (sum(select(shadow4_sim, starts_with("cum_hbv_deaths"))[which(shadow4_sim$time == 2012.5),]) +
        sum(select(shadow4_sim, starts_with("cum_background_deaths_ld"))[which(shadow4_sim$time == 2012.5),]))

  ## 6f) SHADOW MODEL 5: YANG HCC COHORT

  # To fit survival curve at time intervals 0.5, 1 and 1.5 years
  # Follow a cohort of HCC patients, starting in 2012, for 1.5 years
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow5_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                           init_compartment_from = 8, init_compartment_to = 8,
                           shadow_default_parameter_list = parameters_for_fit,
                           shadow_init = model_pop2012, shadowsim_duration = 2,
                           shadow_parms_to_change = list(sim_starttime = 2012,
                                                         births_on = 0,
                                                         migration_on = 0,
                                                         b1 = 0,
                                                         b2 = 0,
                                                         b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from HCC or background) after 0.5, 1 and 1.5 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because DCC and HCC compartments are empty
  # Denominator = number in compensated cirrhosis compartment at t0

  for (i in 1:3) {  # i = timestep
    model_mort_curves[model_mort_curves$outcome=="shadow5_cum_mortality","model_value"][i] <-
      ((sum(select(shadow5_sim, starts_with("cum_hcc_deaths"))[(i+1),])) +
    (sum(select(shadow5_sim, starts_with("cum_background_deaths_ld"))[(i+1),])))/
    sum(select(shadow5_sim, starts_with("HCC"))[which(shadow5_sim$time == 2012),])
  }

  ## 6g) SHADOW MODEL 6: DIARRA CIRRHOSIS COHORT

  # To fit survival curve at time intervals 0.5 and 1
  # Follow a cohort of CC and DCC patients, starting in 2005, for 1 year
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow6_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 6, init_compartment_to = 7,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop2005, shadowsim_duration = 1.5,
                                  shadow_parms_to_change = list(sim_starttime = 2005,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))

  # Calculate output for fitting:
  # Cumulative mortality probability (from cirrhosis, HCC or background) after 0.5 and 1 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because HCC compartment is empty at baseline
  # Denominator = number in CC and DCC compartments at t0

  for (i in 1:2) {  # i = timestep
    model_mort_curves[model_mort_curves$outcome=="shadow6_cum_mortality","model_value"][i] <-
      ((sum(select(shadow6_sim, starts_with("cum_hbv_deaths"))[(i+1),])) +
                                   (sum(select(shadow6_sim, starts_with("cum_background_deaths_ld"))[(i+1),])))/
      (sum(select(shadow6_sim, starts_with("CC"))[which(shadow6_sim$time == 2005),]) +
         sum(select(shadow6_sim, starts_with("DCC"))[which(shadow6_sim$time == 2005),]))
  }

  # Cumulative HCC incidence after 0.5 and 1 years
  # Numerator = sum of icnident HCC cases at given timestep
  # Can use HCC cases coming from any compartment because there are no people in the other
  # compartments to transition
  # Denominator = number in CC and DCC compartments at t0

  for (i in 1:2) {  # i = timestep index, so 1 = t0
    model_mort_curves[model_mort_curves$outcome=="shadow6_cum_hcc_incidence","model_value"][i] <-
      (sum(select(shadow6_sim, starts_with("cum_incident_hcc"))[(i+1),]))/
      (sum(select(shadow6_sim, starts_with("CC"))[which(shadow6_sim$time == 2005),]) +
         sum(select(shadow6_sim, starts_with("DCC"))[which(shadow6_sim$time == 2005),]))
  }

  ## 6h) SHADOW MODEL 7: GLOBOCAN SURVIVAL CURVE

  # To fit survival curve at time intervals 1, 3 and 5
  # Follow a cohort of HCC patients, starting in 1993, until 1997.5
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow7_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 8, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1995, shadowsim_duration = 5.5,
                                  shadow_parms_to_change = list(sim_starttime = 1995,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))


  # Calculate output for fitting:
  # Cumulative mortality probability (from HCC or background) after 1, 3 and 5 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because CC and DCC compartments are empty at baseline
  # Denominator = number in HCC compartment at t01
  mapped_globocan_mortality_curve <- data_to_fit$globocan_mortality_curve
  mapped_globocan_mortality_curve$model_value[1] <- ((sum(select(shadow7_sim, starts_with("cum_hcc_deaths"))[which(shadow7_sim$time == 1996),])) +
                                   (sum(select(shadow7_sim, starts_with("cum_background_deaths_ld"))[which(shadow7_sim$time == 1996),])))/
      (sum(select(shadow7_sim, starts_with("HCC"))[which(shadow7_sim$time == 1995),]))

  mapped_globocan_mortality_curve$model_value[2] <- ((sum(select(shadow7_sim, starts_with("cum_hcc_deaths"))[which(shadow7_sim$time == 1998),])) +
        (sum(select(shadow7_sim, starts_with("cum_background_deaths_ld"))[which(shadow7_sim$time == 1998),])))/
      (sum(select(shadow7_sim, starts_with("HCC"))[which(shadow7_sim$time == 1995),]))

  mapped_globocan_mortality_curve$model_value[3] <- ((sum(select(shadow7_sim, starts_with("cum_hcc_deaths"))[which(shadow7_sim$time == 2000),])) +
                                 (sum(select(shadow7_sim, starts_with("cum_background_deaths_ld"))[which(shadow7_sim$time == 2000),])))/
      (sum(select(shadow7_sim, starts_with("HCC"))[which(shadow7_sim$time == 1995),]))


  ## Combine model predictions with input datasets
  mapped_progression_rates <- left_join(data_to_fit$progression_rates, progression_rates, by = "outcome")
  mapped_mortality_curves <- cbind(data_to_fit$mortality_curves, model_value = model_mort_curves$model_value)
  mapped_mortality_curves <- bind_rows(mapped_mortality_curves,
                                       mapped_globocan_mortality_curve) # add GLOBOCAN shadow model


  ## 7) ASSOCIATION DATA ----

  # Prepare dataframe
  association_output <- data.frame(outcome = c("odds_ratio_current_hbeag_positivity_and_hcc",
                                               "odds_ratio_current_hbeag_positivity_and_cirrhosis",
                                               "odds_ratio_male_sex_and_significant_liver_fibrosis_or_cirrhosis"))

  # GMB15 GLCS: Association of current HBeAg status and HCC
  # Outcome = HCC, exposure = concurrent HBeAg status
  # Cases = HCC compartment, controls = all other carrier compartments minus CC and DCC
  # because their HBeAg status is not known
  # Cases and controls were recruited in 1999 (midpoint) and were aged 15-83.5 years
  # 83% of participants were male so assume this represents the OR in males
  # Age range was derived from mean age +/- 2*SD
  # Use the proportion of incident non-cirrhotic HCC cases originating from IT and IR
  # as a proxy for the proportion HBeAg-positive among cases and multiply with the total cases

  assoc_hbeag_hcc <- data.frame(
    total_cases = sum(select(sim, starts_with("HCCm"))[which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]),
    # Exposed controls = IT and IR compartment, which are HCC-free by definition
    exposed_controls = sum(select(sim, starts_with("ITm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                             select(sim, starts_with("IRm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    # Unexposed controls = IC and ENCHB compartment
    unexposed_controls = sum(select(sim, starts_with("ICm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                               select(sim, starts_with("ENCHBm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    prop_exposed_cases = (sum(select(sim, starts_with("cum_it_to_hccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
                            sum(select(sim, starts_with("cum_ir_to_hccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]))/
      (sum(select(sim, starts_with("cum_it_to_hccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(select(sim, starts_with("cum_ir_to_hccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(select(sim, starts_with("cum_ic_to_hccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])+
         sum(select(sim, starts_with("cum_enchb_to_hccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]))
  )

  # Exposed cases = HBeAg-positive HCC patients (approximated)
  assoc_hbeag_hcc$exposed_cases <- assoc_hbeag_hcc$total_cases *
    assoc_hbeag_hcc$prop_exposed_cases
  # Unexposed cases = HBeAg-negative HCC patients
  assoc_hbeag_hcc$unexposed_cases <- assoc_hbeag_hcc$total_cases -
    assoc_hbeag_hcc$exposed_cases

  hbeag_hcc_odds_ratio <- (
    assoc_hbeag_hcc$exposed_cases * assoc_hbeag_hcc$unexposed_controls)/
    (assoc_hbeag_hcc$unexposed_cases * assoc_hbeag_hcc$exposed_controls)

  # GMB15 GLCS: Association of current HBeAg status and cirrhosis
  # Outcome = CC or DCC, exposure = concurrent HBeAg status
  # Cases = CC+DCC compartment, controls = all other carrier compartments minus HCC
  # because their HBeAg status is not known
  # Cases and controls were recruited in 1999 (midpoint) and were aged 15-83.5 years
  # 80% of participants were male so assume this represents the OR in males
  # Age range was derived from mean age +/- 2*SD
  # Use the proportion of incident CC cases originating from IR as a proxy for
  # the proportion HBeAg-positive among cases and multiply this with the total cases

  assoc_hbeag_cirrhosis <- data.frame(
    total_cases = sum(select(sim, starts_with("CCm"))[which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]+
                        select(sim, starts_with("DCCm"))[which(sim$time == 1999), which(ages == 15):which(ages == 83.5)]),
    total_pop = sum(out$carriers_male[which(out$time == 1999),which(ages == 15):which(ages == 83.5)])-
      sum(select(sim, starts_with("HCCm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    # Exposed controls = IT and IR compartment, which are cirrhosis-free by definition
    exposed_controls = sum(select(sim, starts_with("ITm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
                             select(sim, starts_with("IRm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]),
    prop_exposed_cases = sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])/
      sum(select(sim, starts_with("cum_ir_to_ccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)]+
            select(sim, starts_with("cum_enchb_to_ccm"))[which(sim$time == 1999),which(ages == 15):which(ages == 83.5)])
  )

  # Exposed cases = HBeAg-positive CC or DCC patients (approximated)
  assoc_hbeag_cirrhosis$exposed_cases <- assoc_hbeag_cirrhosis$total_cases *
    assoc_hbeag_cirrhosis$prop_exposed_cases
  # Unexposed cases = HBeAg-negative CC or DCC patients
  assoc_hbeag_cirrhosis$unexposed_cases <- assoc_hbeag_cirrhosis$total_cases -
    assoc_hbeag_cirrhosis$exposed_cases
  assoc_hbeag_cirrhosis$total_controls <- assoc_hbeag_cirrhosis$total_pop -
    assoc_hbeag_cirrhosis$total_cases
  # Unexposed controls = HBeAg-negative non-cirrhotic carriers
  assoc_hbeag_cirrhosis$unexposed_controls <- assoc_hbeag_cirrhosis$total_controls -
    assoc_hbeag_cirrhosis$exposed_controls

  hbeag_cirrhosis_odds_ratio <- (
    assoc_hbeag_cirrhosis$exposed_cases * assoc_hbeag_cirrhosis$unexposed_controls)/
    (assoc_hbeag_cirrhosis$unexposed_cases * assoc_hbeag_cirrhosis$exposed_controls)

  # 1-1 Shimakawa: Association of sex and significant liver fibrosis or cirrhosis
  # Use outcome as a proxy for CC/DCC, exposure = male sex
  # Cases = CC+DCC compartment, controls = all other carrier compartments minus HCC
  # Cases and controls were recruited in 2013 and were aged 8-95.5 years
  assoc_sex_cirrhosis <- data.frame(
    # Exposed cases = males with CC or DCC
    exposed_cases = sum(select(sim, starts_with("CCm"))[which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]+
                          select(sim, starts_with("DCCm"))[which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    # Unexposed cases = females with CC or DCC
    unexposed_cases = sum(select(sim, starts_with("CCf"))[which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]+
                            select(sim, starts_with("DCCf"))[which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    total_exposed = sum(out$carriers_male[which(out$time == 2013),which(ages == 8):which(ages == 95.5)])-
      sum(select(sim, starts_with("HCCm"))[which(sim$time == 2013),which(ages == 8):which(ages == 95.5)]),
    total_unexposed = sum(out$carriers_female[which(out$time == 2013),which(ages == 8):which(ages == 95.5)])-
      sum(select(sim, starts_with("HCCf"))[which(sim$time == 2013),which(ages == 8):which(ages == 95.5)])
  )

  # Exposed controls = males in IT, IR, IC or ENCHB compartment
  assoc_sex_cirrhosis$exposed_controls <- assoc_sex_cirrhosis$total_exposed -
    assoc_sex_cirrhosis$exposed_cases
  # Unexposed controls = females in IT, IR, IC or ENCHB compartment
  assoc_sex_cirrhosis$unexposed_controls <- assoc_sex_cirrhosis$total_unexposed -
    assoc_sex_cirrhosis$unexposed_cases

  sex_cirrhosis_odds_ratio <-
    (assoc_sex_cirrhosis$exposed_cases * assoc_sex_cirrhosis$unexposed_controls)/
    (assoc_sex_cirrhosis$unexposed_cases * assoc_sex_cirrhosis$exposed_controls)


  # Combine odds ratios in dataframe
  association_output$model_value <- c(hbeag_hcc_odds_ratio,
                                      hbeag_cirrhosis_odds_ratio,
                                      sex_cirrhosis_odds_ratio)

  # Merge with data to fit to (first transform factor to character vector)
  association_output$outcome <- as.character(association_output$outcome)
  mapped_odds_ratios <- left_join(data_to_fit$odds_ratios,
                                  association_output,
                                  by = "outcome")

  # Combine all mapped outputs ----
  mapped_output_complete <- list(globocan_hcc_incidence = mapped_globocan_incidence,
                                 risk_of_chronic_carriage = mapped_p_chronic,
                                 seromarker_prevalence = mapped_seromarker_prevalence,
                                 nat_hist_prevalence = mapped_nat_hist_prevalence,
                                 mtct_risk = mapped_mtct_risk,
                                 progression_rates = mapped_progression_rates,
                                 mortality_curves = mapped_mortality_curves,
                                 odds_ratios = mapped_odds_ratios)


  # Calculate sum of least squares
  data_model_diff <- as.numeric(unlist(lapply(mapped_output_complete, function(x) x$data_value)))-
    as.numeric(unlist(lapply(mapped_output_complete, function(x) x$model_value)))

  sse <- sum(data_model_diff^2)

  # Return relevant info (SSE and the matched datapoints and outputs)
  res <- list(sse = sse, mapped_output = mapped_output_complete, full_output = out)

  return(res)

}



load(here("data/simulated_inits_1880.RData"))  # this is saved from previous model run
init_pop_sim <- c("Sf" = select(model_pop1880, starts_with("Sf")),
                  "ITf" = select(model_pop1880, starts_with("ITf")),
                  "IRf" = select(model_pop1880, starts_with("IRf")),
                  "ICf" = select(model_pop1880, starts_with("ICf")),
                  "ENCHBf" = select(model_pop1880, starts_with("ENCHBf")),
                  "CCf" = select(model_pop1880, starts_with("CCf")),
                  "DCCf" = select(model_pop1880, starts_with("DCCf")),
                  "HCCf" = select(model_pop1880, starts_with("HCCf")),
                  "Rf" = select(model_pop1880, starts_with("Rf")),
                  "Sm" = select(model_pop1880, starts_with("Sm")),
                  "ITm" = select(model_pop1880, starts_with("ITm")),
                  "IRm" = select(model_pop1880, starts_with("IRm")),
                  "ICm" = select(model_pop1880, starts_with("ICm")),
                  "ENCHBm" = select(model_pop1880, starts_with("ENCHBm")),
                  "CCm" = select(model_pop1880, starts_with("CCm")),
                  "DCCm" = select(model_pop1880, starts_with("DCCm")),
                  "HCCm" = select(model_pop1880, starts_with("HCCm")),
                  "Rm" = select(model_pop1880, starts_with("Rm")),
                  output_storage)
init_pop_sim <- unlist(init_pop_sim)


# Define my datapoints to fit to: HBsAg prevalence dataset
calibration_data_path <- "data/calibration_data"
input_hbsag_dataset <- read.csv(here(calibration_data_path,
                                     "hbsag_prevalence.csv"),
                             header = TRUE, check.names = FALSE,
                             stringsAsFactors = FALSE)
input_antihbc_dataset <- read.csv(here(calibration_data_path,
                                     "antihbc_prevalence.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)
input_hbeag_dataset <- read.csv(here(calibration_data_path,
                                     "hbeag_prevalence.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)
input_natural_history_prev_dataset <- read.csv(here(calibration_data_path,
                                     "natural_history_prevalence.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)

input_mtct_risk_dataset <- read.csv(here(calibration_data_path,
                                         "mtct_risk.csv"),
                                         header = TRUE, check.names = FALSE,
                                         stringsAsFactors = FALSE)

input_progression_rates <- read.csv(here(calibration_data_path,
                                         "progression_rates.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

input_mortality_curves <- read.csv(here(calibration_data_path,
                                         "mortality_curves.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

input_risk_of_chronic_carriage <- read.csv(here(calibration_data_path,
                                        "risk_of_chronic_carriage.csv"),
                                   header = TRUE, check.names = FALSE,
                                   stringsAsFactors = FALSE)

input_globocan_incidence_data <- read.csv(here(calibration_data_path,
                                               "globocan_incidence_data.csv"),
                                          header = TRUE, check.names = FALSE,
                                          stringsAsFactors = FALSE)

input_globocan_mortality_curve <- read.csv(here(calibration_data_path,
                                               "globocan_mortality_curve.csv"),
                                          header = TRUE, check.names = FALSE,
                                          stringsAsFactors = FALSE,
                                          fileEncoding="UTF-8-BOM")

input_odds_ratios <- read.csv(here(calibration_data_path,
                                   "odds_ratios.csv"),
                              header = TRUE, check.names = FALSE,
                              stringsAsFactors = FALSE)


# Need to change name of this list
calibration_datasets_list <- list(hbsag_prevalence = input_hbsag_dataset,
                                 antihbc_prevalence = input_antihbc_dataset,
                                 hbeag_prevalence = input_hbeag_dataset,
                                 natural_history_prevalence = input_natural_history_prev_dataset,
                                 mtct_risk = input_mtct_risk_dataset,
                                 progression_rates = input_progression_rates,
                                 mortality_curves = input_mortality_curves,
                                 globocan_incidence_data = input_globocan_incidence_data,
                                 globocan_mortality_curve = input_globocan_mortality_curve,
                                 p_chronic = input_risk_of_chronic_carriage,
                                 odds_ratios = input_odds_ratios)


# Input datasets
prior_vaccine_efficacy <- read.csv(here("data-raw", "input_infant_vaccine_efficacy.csv"),
                                   header = TRUE, check.names = FALSE,
                                   stringsAsFactors = FALSE)
prior_mtct_risk <- read.csv(here("data-raw", "input_mtct_risk.csv"),
                            header = TRUE, check.names = FALSE,
                            stringsAsFactors = FALSE)
prior_paf_liver_disease <- read.csv(here("data-raw", "input_paf_liver_disease.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)
prior_progression_rates <- read.csv(here("data-raw", "input_progression_rates.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

# Using LHS
n_sims <- 1  # number of simulations
n_parms_to_vary <- 3  # number of parameters to infer - this requires manual adaptations below
lhs_samples <- randomLHS(n_sims, n_parms_to_vary) # draw 100 samples from uniform distribution U(0,1) using a Latin Hypercube design
params_mat <- data.frame(b1 = lhs_samples[,1],
                         b2 = lhs_samples[,2],
                         mtct_prob_s = lhs_samples[,3])
params_mat$b1 <- 0 + (0.2-0) * params_mat$b1 # rescale U(0,1) to be U(0,0.2)
params_mat$b2 <- 0 + (0.01-0) * params_mat$b2 # rescale U(0,1) to be U(0,0.01)
params_mat$mtct_prob_s <- 0 + (0.5-0) * params_mat$mtct_prob_s # rescale U(0,1) to be U(0,0.5)

#params_mat <- data.frame(b1 = 0.1, b2 = 0.009, mtct_prob_s = 0.14)

# Run without parallelising
time1 <- proc.time()
out_mat <- apply(params_mat,1,
                    function(x)
                      fit_model_sse(default_parameter_list = parameter_list,
                                    data_to_fit = calibration_datasets_list,
                                    parms_to_change =
                                      list(b1 = as.list(x)$b1,
                                           b2 = as.list(x)$b2,
                                           mtct_prob_s = as.list(x)$mtct_prob_s)))
sim_duration = proc.time() - time1
sim_duration["elapsed"]/60

# Matrix of parameter values, model estimates for prevalence in 1980 and 2015, and SSE
out_mat_subset <- sapply(out_mat, "[[", "sse")
res_mat <- cbind(params_mat, sse = out_mat_subset)
res_mat[res_mat$sse == min(res_mat$sse),]


# Output plots
library(gridExtra)

# Legend for plots: black lines = model, red cross is data

# HBsAg prevalence by time and age
plot_hbsag <- ggplot(data = out_mat[[1]]$mapped_output$seromarker_prevalence[
  out_mat[[1]]$mapped_output$seromarker_prevalence$outcome == "HBsAg_prevalence",]) +
  geom_line(aes(x = age, y = model_value, group  = sex, colour = sex)) +
  geom_point(aes(x = age, y = data_value, group = sex, colour = sex), shape = 4, stroke = 1.5) +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, group = sex, colour = sex)) +
  facet_wrap(~ time, ncol = 3) +
  labs(title = "HBsAg", y = "Prevalence (proportion)", x = "Age (years)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,0.6)

# Anti-HBc prevalence by time and age
plot_antihbc <- ggplot(data = out_mat[[1]]$mapped_output$seromarker_prevalence[
  out_mat[[1]]$mapped_output$seromarker_prevalence$outcome == "Anti_HBc_prevalence",]) +
  geom_line(aes(x = age, y = model_value, group  = sex, colour = sex)) +
  geom_point(aes(x = age, y = data_value, group = sex, colour = sex), shape = 4, stroke = 1.5) +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, group = sex, colour = sex)) +
  facet_wrap(~ time, ncol = 3) +
  labs(title = "Anti-HBc", y = "Prevalence (proportion)", x = "Age (years)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

## Combined HBsAg and anti-HBc prevalence by time and age
grid.arrange(plot_hbsag, plot_antihbc, nrow = 1)

## HBeAg prevalence by time and age
ggplot(data = out_mat[[1]]$mapped_output$seromarker_prevalence[
  out_mat[[1]]$mapped_output$seromarker_prevalence$outcome == "HBeAg_prevalence",]) +
  geom_line(aes(x = age, y = model_value, group  = sex, colour = sex)) +
  geom_point(aes(x = age, y = data_value, group = sex, colour = sex), shape = 4, stroke = 1.5) +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower, group = sex, colour = sex)) +
  facet_wrap(~ time, ncol = 3) +
  labs(title = "HBeAg", y = "Prevalence (proportion)", x = "Age (years)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

## GLOBOCAN PAF-adjusted cancer incidence and mortality
ggplot(data = out_mat[[1]]$mapped_output$globocan_hcc_incidence) +
  geom_col(aes(x = paste(age_min,"-",age_max), y = model_value*100000)) +
#  geom_errorbar(aes(x = paste(age_min,"-",age_max), ymax = ci_upper, ymin = ci_lower)) +
  geom_point(aes(x = paste(age_min,"-",age_max), y = data_value*100000), col = "red",
             shape = 4, stroke = 1.5) +
  facet_grid(outcome ~ sex) +
  labs(title = "GLOBOCAN HBV-related HCC incidence and mortality rates 2018",
       y = "Cases/deaths per 100000", x = "Age (years)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,100)

## Risk of chronic carriage
ggplot(data = out_mat[[1]]$mapped_output$risk_of_chronic_carriage) +
  geom_point(aes(x = age, y = data_value), col = "red", shape = 4, stroke = 1.5) +
  geom_errorbar(aes(x = age, ymax = ci_upper, ymin = ci_lower), col = "red") +
  geom_line(aes(x = age, y = model_value)) +
  labs(title = "Risk of chronic carriage by age at infection",
       y = "Risk (proportion)", x = "Age (years)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1) + xlim(0,30)

## Mortality curves
# Add articifial zeros at first timestep to allow plotting of step curves
mortality_curves_zeros <- out_mat[[1]]$mapped_output$mortality_curves
mortality_curves_zeros$time_interval_years <- 0
mortality_curves_zeros$data_value <- 0
mortality_curves_zeros$model_value <- 0
mortality_curves_zeros$number_at_risk <- mortality_curves_zeros$sample_size
mortality_curves_zeros <- unique(mortality_curves_zeros)

ggplot(data = rbind(out_mat[[1]]$mapped_output$mortality_curves, mortality_curves_zeros)) +
  geom_step(aes(x = time_interval_years, y = model_value)) +
  geom_point(aes(x = time_interval_years, y = data_value), col = "red",
             shape = 4, stroke = 1.5) +
  facet_grid(~ outcome, scales = "free") +
  labs(title = "Cumulative probability of death/HCC over time",
       y = "Cumulative probability", x = "Follow-up time (years)") +
  theme(plot.title = element_text(hjust = 0.5))

## ORs
ggplot(data = out_mat[[1]]$mapped_output$odds_ratios) +
  geom_col(aes(x = gsub("odds_ratio_", "", outcome), y = model_value)) +
  geom_point(aes(x = gsub("odds_ratio_", "", outcome), y = data_value),
             col = "red", shape = 4, size = 3, stroke = 2) +
  geom_errorbar(aes(x = gsub("odds_ratio_", "", outcome),
                    ymax = ci_upper, ymin = ci_lower), col= "red", width = 0.2) +
  labs(title = "Odds ratios",
       y = "OR", x = "Exposure and outcome") +
  theme(plot.title = element_text(hjust = 0.5))

# Natural history prevalence plots
out_mat[[1]]$mapped_output$nat_hist_prevalence$model_num <-
  gsub(".*[[:digit:]]{4}_", "",out_mat[[1]]$mapped_output$nat_hist_prevalence$id_unique)

# GMB1 plots: infection phase in chronic carriers
plot1_gmb1 <- ggplot(data = subset(out_mat[[1]]$mapped_output$nat_hist_prevalence,
                     id_paper == "GMB1" &
                     model_num != "cc_dcc" & model_num != "hcc")) +
  geom_col(aes(x = model_num, y = model_value))+
  geom_point(aes(x = model_num, y = data_value), shape = 4, size = 3, stroke = 2, col = "red") +
  facet_grid(~sex, scales = "free") +
  labs(title = "GMB1 Infection phase in chronic carriers",
       y = "Prevalence (proportion)", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

# GMB1 plots: liver disease in chronic carriers
plot2_gmb1 <-ggplot(data = subset(out_mat[[1]]$mapped_output$nat_hist_prevalence,
                     id_paper == "GMB1" &
                       (model_num == "cc_dcc" | model_num == "hcc"))) +
  geom_col(aes(x = model_num, y = model_value))+
  geom_point(aes(x = model_num, y = data_value), shape = 4, size = 3, stroke = 2, col = "red") +
  facet_grid(~sex, scales = "free") +
  labs(title = "GMB1 Liver disease in chronic carriers",
       y = "Prevalence (proportion)", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,0.1)

# 1-1 plots: infection phase in chronic carriers without liver disease
plot1_1 <- ggplot(data = subset(out_mat[[1]]$mapped_output$nat_hist_prevalence,
                     grepl(".*it,_ir,_ic_and_enchb$", out_mat[[1]]$mapped_output$nat_hist_prevalence$outcome))) +
  geom_col(aes(x = model_num, y = model_value))+
  geom_point(aes(x = model_num, y = data_value), shape = 4, size = 3, stroke = 2, col = "red") +
  facet_grid(~time, scales = "free") +
  labs(title = "1 Infection phase in chronic carriers without liver disease",
       y = "Prevalence (proportion)", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

# 1-1 plots: liver disease in chronic carriers
plot2_1 <-ggplot(data = subset(out_mat[[1]]$mapped_output$nat_hist_prevalence,
                     id_paper == "1" &
                     (outcome == "hcc_prevalence_in_chronic_carriers" |
                      outcome == "cc_and_dcc_prevalence_in_chronic_carriers"))) +
  geom_col(aes(x = model_num, y = model_value))+
  geom_point(aes(x = model_num, y = data_value), shape = 4, size = 3, stroke = 2, col = "red") +
  facet_grid(~time, scales = "free_x") +
  labs(title = "1 Liver disease in chronic carriers",
       y = "Prevalence (proportion)", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,0.1)

# 1-1 plots: chronic carriers by age
plot3_1 <- ggplot(data = out_mat[[1]]$mapped_output$nat_hist_prevalence[
  out_mat[[1]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_2013_ir_enchb_cc_dcc",]) +
  geom_col(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = model_value))+
  geom_point(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = data_value),
             shape = 4, size = 3, stroke = 2, col = "red") +
  labs(title = "1 Significant liver fibrosis or cirrhosis in chronic carriers",
       y = "Prevalence (proportion)", x = "Age") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,0.75)

## Natural history prevalence PLOT 1
nat_hist_prev_plot1 <- grid.arrange(plot1_gmb1, plot1_1, plot2_gmb1, plot2_1, plot3_1, nrow = 3)

# A4
plot_nat_hist_a4 <- ggplot(data = subset(out_mat[[1]]$mapped_output$nat_hist_prevalence,
                     id_paper == "A4")) +
  geom_col(aes(x = model_num, y = model_value))+
  geom_point(aes(x = model_num, y = data_value),
             shape = 4, size = 3, stroke = 2, col = "red") +
  labs(title = "Proportion of deaths from DCC and HCC\n in A4 compensated cirrhosis cohort",
       y = "Proportion", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

# GMB2
plot_nat_hist_gmb2 <- ggplot(data = subset(out_mat[[1]]$mapped_output$nat_hist_prevalence,
                     id_paper == "GMB2")) +
  geom_col(aes(x = model_num, y = model_value))+
  geom_point(aes(x = model_num, y = data_value),
             shape = 4, size = 3, stroke = 2, col = "red") +
  labs(title = "Cirrhosis in HBsAg-positive HCC patients",
       y = "Proportion", x = "Age") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

# HBeAg prevalence in liver disease patients: GMB12 and GMB15
plot_nat_hist_glcs <- ggplot(data = subset(out_mat[[1]]$mapped_output$nat_hist_prevalence,
                     outcome == "hbeag_prevalence_in_hcc" |
                       outcome == "hbeag_prevalence_in_cirrhosis")) +
  geom_col(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = model_value))+
  geom_point(aes(x = reorder(paste(age_min,"-",age_max), age_min), y = data_value),
             shape = 4, size = 3, stroke = 2, col = "red") +
  facet_grid(~id_unique, scales = "free_x") +
  labs(title = "HBeAg prevalence in HBsAg-positive HCC/cirrhosis patients",
       y = "Prevalence (proportion)", x = "Age") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

## Natural history prevalence PLOT 2
nat_hist_prev_plot2 <- grid.arrange(plot_nat_hist_a4,  plot_nat_hist_gmb2, plot_nat_hist_glcs, nrow = 2,
             layout_matrix = rbind(c(1,2),
                                   c(3,3)))

# Vertical transmission plots
# 1-1 plots: chronic infections due to vertical transmission
plot_nat_hist_1 <- ggplot(data = out_mat[[1]]$mapped_output$nat_hist_prevalence[
  out_mat[[1]]$mapped_output$nat_hist_prevalence$id_unique == "id_1_1_1986_incident_chronic_births",]) +
  geom_col(aes(x = model_num, y = model_value))+
  geom_point(aes(x = model_num, y = data_value),
             shape = 4, size = 3, stroke = 2, col = "red") +
  geom_errorbar(aes(x = model_num, ymax = ci_upper, ymin = ci_lower), col = "red", width = 0.1) +
  labs(title = "Chronic infections due to vertical transmission",
       y = "Proportion", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

# MTCT risk
plot_mtct <- ggplot(data = out_mat[[1]]$mapped_output$mtct_risk) +
  geom_col(aes(x = id_paper, y = model_value))+
  geom_point(aes(x = id_paper, y = data_value),
             shape = 4, size = 3, stroke = 2, col = "red") +
  facet_grid(~time, scales = "free_x") +
  labs(title = "Mother-to-child transmission risk",
       y = "Risk (proportion)", x = "Study ID") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1)

## Progression rates
out_mat[[1]]$mapped_output$progression_rates$type <-
  gsub("shadow[[:digit:]].{0,1}_", "", out_mat[[1]]$mapped_output$progression_rates$outcome)
out_mat[[1]]$mapped_output$progression_rates$type <- gsub("_.$", "", out_mat[[1]]$mapped_output$progression_rates$type)

# Study 1: HCC incidence
plot_1_hcc_incidence <- ggplot(data = subset(out_mat[[1]]$mapped_output$progression_rates,
                     id_paper == "1" & type == "hcc_incidence")) +
  geom_col(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = model_value*100000))+
  geom_point(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = data_value*100000),
             shape = 4, size = 3, stroke = 2, col = "red") +
  geom_errorbar(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), ymax = ci_upper*100000, ymin = ci_lower*100000),
                col = "red", width = 0.1)  +
  facet_grid(~sex, scales = "free") +
  labs(title = "HCC incidence in chronic carrier cohort",
       y = "Cases per 100000 PY", x = "Baseline age group (years)") +
  theme(plot.title = element_text(hjust = 0.5))

# Study 1: DCC incidence
plot_1_dcc_incidence <- ggplot(data = subset(out_mat[[1]]$mapped_output$progression_rates,
                     id_paper == "1" & type == "dcc_incidence")) +
  geom_col(aes(x = outcome, y = model_value*100000))+
  geom_point(aes(x = outcome, y = data_value*100000),
             shape = 4, size = 3, stroke = 2, col = "red") +
  geom_errorbar(aes(x = outcome, ymax = ci_upper*100000, ymin = ci_lower*100000),
                col = "red", width = 0.1)  +
  labs(title = "DCC incidence in chronic carrier cohort",
       y = "Cases per 100000 PY", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,100)


# Study 1: Mortality
plot_1_mortality <- ggplot(data = subset(out_mat[[1]]$mapped_output$progression_rates,
                     id_paper == "1" & type == "mortality")) +
  geom_col(aes(x = sex, y = model_value*100000))+
  geom_point(aes(x = sex, y = data_value*100000),
             shape = 4, size = 3, stroke = 2, col = "red") +
  geom_errorbar(aes(x = sex, ymax = ci_upper*100000, ymin = ci_lower*100000),
                col = "red", width = 0.1)  +
  labs(title = "All-cause mortality in chronic carrier cohort",
       y = "Deaths per 100000 PY", x = "Sex") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,1000)

# Study A6: Mortality in cirrhosis cohort
plot_a6_mortality <- ggplot(data = subset(out_mat[[1]]$mapped_output$progression_rates,
                     id_paper == "A6")) +
  geom_col(aes(x = outcome, y = model_value*100))+
  geom_point(aes(x = outcome, y = data_value*100),
             shape = 4, size = 3, stroke = 2, col = "red") +
  geom_errorbar(aes(x = outcome, ymax = ci_upper*100, ymin = ci_lower*100),
                col = "red", width = 0.1)  +
  labs(title = "Mortality in cirrhosis cohort",
       y = "Deaths per 100 PY", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,60)

## Combined liver disease incidence and mortality rates
grid.arrange(plot_1_hcc_incidence, plot_1_dcc_incidence,
             plot_1_mortality, plot_a6_mortality, nrow = 2, widths = 4:3)

# Study 1: HBeAg loss
plot_1_eag_loss <- ggplot(data = subset(out_mat[[1]]$mapped_output$progression_rates,
                                        id_paper == "1" & type == "eag_loss")) +
  geom_col(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = model_value*100))+
  geom_point(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), y = data_value*100),
             shape = 4, size = 3, stroke = 2, col = "red") +
  geom_errorbar(aes(x = paste(bl_age_min_years,"-",bl_age_max_years), ymax = ci_upper*100, ymin = ci_lower*100),
                col = "red", width = 0.1)  +
  facet_grid(~sex, scales = "free") +
  labs(title = "HBeAg loss in chronic carrier cohort",
       y = "Cases per 100 PY", x = "Baseline age group (years)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,12)

# Study 6: HBsAg loss
plot_6_sag_loss <- ggplot(data = subset(out_mat[[1]]$mapped_output$progression_rates,
                     id_paper == "6")) +
  geom_col(aes(x = outcome, y = model_value*100))+
  geom_point(aes(x = outcome, y = data_value*100),
             shape = 4, size = 3, stroke = 2, col = "red") +
  geom_errorbar(aes(x = outcome, ymax = ci_upper*100, ymin = ci_lower*100),
                col = "red", width = 0.1)  +
  labs(title = "HBsAg loss in chronic carrier children",
       y = "Cases per 100 PY", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,5)

## Combined seromarker loss rates
grid.arrange(plot_1_eag_loss, plot_6_sag_loss, nrow = 1, widths = 2:1)

# Transmission-related data from GMB6 and GMB7
plot_horizontal_transmission <- ggplot(data = subset(out_mat[[1]]$mapped_output$progression_rates,
                     id_paper == "GMB6" | id_paper == "GMB7")) +
  geom_col(aes(x = gsub(".*_","",outcome), y = model_value))+
  geom_point(aes(x = gsub(".*_","",outcome), y = data_value),
             shape = 4, size = 3, stroke = 2, col = "red") +
  geom_errorbar(aes(x = gsub(".*_","",outcome), ymax = ci_upper, ymin = ci_lower),
                col = "red", width = 0.1)  +
  labs(title = "Horizontal transmission-related rates",
       y = "Rate (per PY)", x = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Force of infection", "Chronic infection incidence")) +
  ylim(0,2)


## Combined transmission-related plot
grid.arrange(plot_mtct, plot_nat_hist_1, plot_horizontal_transmission,
             layout_matrix = rbind(c(1,1),
                                 c(2,3)))

# HBsAg over time?
plot(x = out_mat[[1]]$full_output$time,
  y = apply(out_mat[[1]]$full_output$carriers,1,sum)/apply(out_mat[[1]]$full_output$pop,1,sum),
  ylim = c(0,0.5))

# Parallelised code ----
# Set up cluster
cl <- makeCluster(4)
clusterEvalQ(cl, {library(dplyr); library(tidyr); library(deSolve)})
clusterExport(cl, ls())

time1 <- proc.time()
out_mat <- parApply(cl = cl, params_mat,1,
                    function(x) fit_model_sse(default_parameter_list = parameter_list,
                                              data_to_fit = calibration_datasets_list,
                                              parms_to_change = list(b1 = as.list(x)$b1,
                                                                     b2 = as.list(x)$b2,
                                                                     mtct_prob_s = as.list(x)$mtct_prob_s)))
sim_duration = proc.time() - time1
sim_duration["elapsed"]/60
# Timing: 4.9 min for 20 sim unparallelised, 2.5 min when parallelised,
# 12 min for 100 sims in parallel

# Important: stop cluster!!
stopCluster(cl)

#res_mat <- cbind(params_mat, do.call(rbind.data.frame, out_mat_subset)) # this would work for a list
#out_mat_subset <- as.data.frame(t(sapply(out_mat, "[", c("prev_est_1980", "prev_est_2015", "sse"))))
#res_mat <- cbind(params_mat, unnest(out_mat_subset))

# Fit to overall prevalence for minimum SSE
plot(x = seq(1960,2019.5, by = 0.5),
     y = out_mat[[which(res_mat$sse == min(res_mat$sse))]]$carrier_prev_total,
     type = "l", ylim = c(0,0.5))
points(x = hbsag_dataset$time,
       y = hbsag_dataset$prev,
       col = "red")

# Fit to age-specific prevalence for minimum SSE
#plot(x = ages,
#     y = out_mat[[which(res_mat$sse == min(res_mat$sse))]]$prev_by_age_1980,
#     type = "l", ylim = c(0,0.5))
#points(x = hbsag_by_age_dataset_1980$age,
#       y = hbsag_by_age_dataset_1980$prev,
#       col = "red")


# Target fitting approach / filtration
res_mat$fit <- 0

res_mat[(res_mat$prev_est_1980 >= hbsag_dataset$ci_lower[1]) &
          (res_mat$prev_est_1980 <= hbsag_dataset$ci_upper[1]) &
          (res_mat$prev_est_2015 >= hbsag_dataset$ci_lower[2]) &
          (res_mat$prev_est_2015 <= hbsag_dataset$ci_upper[2]),]$fit <- 1
table(res_mat$fit)

# Box plot of parameter estimates
boxplot(subset(res_mat,fit==1)$b1,subset(res_mat,fit==1)$b2,subset(res_mat,fit==1)$mtct_prob_s,
        names = c("b1", "b2", "mtct_prob_s"),ylim=c(0,0.5))

# Box plot of priors and posteriors
par(mfrow=c(1,3))
boxplot(res_mat$b1, subset(res_mat,fit==1)$b1, col=c(grey(0.6),2),ylim=c(0,0.6),
        names=c("Prior","Posterior"), main="LHS for 'b1':\nprior/posteriors distributions")

boxplot(res_mat$b2, subset(res_mat,fit==1)$b2, col=c(grey(0.6),2),ylim=c(0,0.6),
        names=c("Prior","Posterior"), main="LHS for 'b2':\nprior/posteriors distributions")

boxplot(res_mat$mtct_prob_s, subset(res_mat,fit==1)$mtct_prob_s, col=c(grey(0.6),2),
        ylim=c(0,0.6), names=c("Prior","Posterior"),
        main="LHS for 'mtct_prob_s':\nprior/posteriors distributions")
par(mfrow=c(1,1))





