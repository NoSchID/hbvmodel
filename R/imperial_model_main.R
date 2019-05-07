###################################
### Imperial HBV model 07/05/19 ###
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

### Load and clean data ----
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

# Annual rate (proportion) of HBsAg loss by age from Shimakawa paper
shimakawa_sagloss <- data.frame(age = c("0-9", "10-19", "20-29", "30-39", "40-49",
                                        "50-70"),
                                sagloss_rate = c(0.001, 0.0046, 0.0101,
                                                 0.0105, 0.0232, 0.0239))

# Assume everyone in the 10/20 year age groups has the same rate and 70-80 year
# olds have the same rate as 50-70 year olds
sagloss_rates <- c(rep(shimakawa_sagloss$sagloss_rate,each = 10/da),
                   rep(shimakawa_sagloss$sagloss_rate[which(shimakawa_sagloss == "50-70")],
                       40/da))

## Calculate age-specific cancer rate progression function (Shevanthi)
#cancer_prog_female <- 1e-07 * ((ages * 100* 0.2) + 2 * exp(0.0953 * ages))
#cancer_prog_male <- 5 * cancer_prog_female


### THE MODEL ----
imperial_model <- function(timestep, pop, parameters, sim_starttime) {

  with(as.list(parameters), {

    # PREPARATION

    # Set up population array with infection compartments
    # Matrix 1 = females, matrix 2 = males, rows = agesteps, columns = infection compartments
    # Example: pop[ages,infectionstatus,sex]
    pop <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    # Demography: define time-varying parameters (calls the row corresponding to current timestep)

    # Ignore: this approach only works for Euler (solving at fixed timestep):
    #fertility_rate <- fert_rates[which(times == timestep),-1]
    #migration_rate <- matrix(c(migration_rates_female[which(times == timestep),-1],
    #                           migration_rates_male[which(times == timestep),-1]),
    #                         ncol = 2)    # 2 columns for sex-specific rates
    #mortality_rate <- matrix(c(mort_rates_female[which(times == timestep),-1],
    #                             mort_rates_male[which(times == timestep),-1]),
    #                         ncol = 2)    # 2 columns for sex-specific rates

    # This approach for varying timesteps calls the rate for the closest timestep:
    # Births, migration and mortality can be switched off through parameter list
    if (births_on == 0) {
      fertility_rate <- rep(0, times = ncol(fert_rates[,-1]))
    } else {
      fertility_rate <- fert_rates[min(which(abs(fert_rates[,1]-timestep)==min(abs(fert_rates[,1]-timestep)))),-1]
    }

    if (migration_on == 0) {
      migration_rate <- matrix(rep(0, times = 2*n_agecat), ncol = 2)
    } else {
    migration_rate <- matrix(c(migration_rates_female[min(which(abs(migration_rates_female[,1]-timestep)==min(abs(migration_rates_female[,1]-timestep)))),-1],
                               migration_rates_male[min(which(abs(migration_rates_male[,1]-timestep)==min(abs(migration_rates_male[,1]-timestep)))),-1]),
                             ncol = 2)    # 2 columns for sex-specific rates
    }

    if (mortality_on == 0) {
      mortality_rate <- matrix(rep(0, times = 2*n_agecat), ncol = 2)
    } else {
    mortality_rate <- matrix(c(mort_rates_female[min(which(abs(mort_rates_female[,1]-timestep)==min(abs(mort_rates_female[,1]-timestep)))),-1],
                               mort_rates_male[min(which(abs(mort_rates_male[,1]-timestep)==min(abs(mort_rates_male[,1]-timestep)))),-1]),
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
    if (apply_vacc == 1 & timestep >= (vacc_introtime-sim_starttime)) {
      vacc_cov = vacc_cov
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
                              ncol = 2, nrow = n_agecat)  # female and male incident HBV-related deaths
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
    dflow_dcc_to_hcc <- matrix(rep(0, 2* n_agecat),
                            ncol = 2, nrow = n_agecat)  # female and male transition from DCC to HCC


    # TRANSMISSION

    # Mother-to-child transmission and births
    dcum_infected_births <- sum(fertility_rate *
                                  (apply(mtct_prob_e * pop[index$ages_wocba,c(IT,IR),1],1,sum) +
                                     apply(mtct_prob_s * pop[index$ages_wocba,IC:HCC,1],1,sum)))
    # infected births come from acute and chronic women of childbearing age
    dcum_chronic_births <- p_chronic[1] * dcum_infected_births # infected babies becoming chronic carriers
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

      # Define some transitions:

      # Demography: Incident deaths and migrants
      deaths[index$ages_all,index$infcat_all,i] <- mortality_rate[index$ages_all,i] *
        pop[index$ages_all,index$infcat_all,i]

      migrants[index$ages_all,index$infcat_all,i] <- migration_rate[index$ages_all,i] *
        pop[index$ages_all,index$infcat_all,i]
      # Applying same age-specific mortality rate to every infection compartment
      # Returns an array with indicent deaths and net migrants for every age (rows), infection state (columns) and sex (arrays)

      # Infection: Incident infections (all) and incident chronic infections
      dcum_infections[index$ages_all,i] <- foi * pop[index$ages_all,S,i]
      dcum_chronic_infections[index$ages_all,i] <-  p_chronic * dcum_infections[index$ages_all,i]
      # Returns a matrix with incident infections for every age (rows) and every sex (columns)

      # Natural history transitions
      # Incident deaths due to HBV (from cirrhosis and HCC)
      dcum_hbv_deaths[index$ages_all,i] <- mu_cc * pop[index$ages_all,CC,i] +
        mu_dcc * pop[index$ages_all,DCC,i] + mu_hcc * pop[index$ages_all,HCC,i]
      # Returns a matrix with incident HBV deaths for every age (rows) and every sex (columns)

      # Incident deaths due to HCC only
      dcum_hcc_deaths[index$ages_all,i] <- mu_hcc * pop[index$ages_all,HCC,i]
      # Returns a matrix with incident HBV deaths for every age (rows) and every sex (columns)

      # Incidence of HBeAg loss: transition from IR to IC and IR to ENCHB
      dcum_eag_loss[index$ages_all,i] <- pr_ir_ic*eag_loss_function * pop[index$ages_all,IR,i] +
        pr_ir_enchb * pop[index$ages_all,IR,i]

      # Transition from IC to R (sAg loss)
      dcum_sag_loss[index$ages_all,i] <- sag_loss * pop[index$ages_all,IC,i]

      # DCC incidence
      dcum_dcc[index$ages_all,i] <- dccrate * pop[index$ages_all,CC,i]

      # Total HCC incidence
      dcum_hcc[index$ages_all,i] <- hccr_it*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IT,i] +
        hccr_ir*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IR,i] +
        hccr_ic*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IC,i] +
        hccr_enchb*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB,i] +
        hccr_cc*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CC,i] +
        hccr_dcc * pop[index$ages_all,DCC,i]

      # HCC incidence in decompensated cirrhotics
      dflow_dcc_to_hcc[index$ages_all,i] <- hccr_dcc * pop[index$ages_all,DCC,i]

      # Transitions between compartments:

      # Susceptibles
      dpop[index$ages_all,S,i] <- -(diff(c(0,pop[index$ages_all,S,i]))/da) -
        # (vacc_cov * vacc_eff * pop[index$ages_all,S,i]) -
        dcum_chronic_infections[index$ages_all,i] -
        (1-p_chronic) * dcum_infections[index$ages_all,i] -
        deaths[index$ages_all,S,i] + migrants[index$ages_all,S,i]

      # Immune tolerant
      dpop[index$ages_all,IT,i] <- -(diff(c(0,pop[index$ages_all,IT,i]))/da) +
        dcum_chronic_infections[index$ages_all,i] -
        pr_it_ir*eag_loss_function * pop[index$ages_all,IT,i] -
        hccr_it*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IT,i] -
        deaths[index$ages_all,IT,i] + migrants[index$ages_all,IT,i]

      # Immune reactive
      dpop[index$ages_all,IR,i] <- -(diff(c(0,pop[index$ages_all,IR,i]))/da) +
        pr_it_ir*eag_loss_function * pop[index$ages_all,IT,i] -
        pr_ir_ic*eag_loss_function * pop[index$ages_all,IR,i] -
        pr_ir_enchb_function * pop[index$ages_all,IR,i] -
        pr_ir_cc_function * pop[index$ages_all,IR,i] -
        hccr_ir*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IR,i] -
        deaths[index$ages_all,IR,i] + migrants[index$ages_all,IR,i]

      # Inactive carrier
      dpop[index$ages_all,IC,i] <- -(diff(c(0,pop[index$ages_all,IC,i]))/da) +
        pr_ir_ic*eag_loss_function * pop[index$ages_all,IR,i] -
        pr_ic_enchb * pop[index$ages_all,IC,i] -
        dcum_sag_loss[index$ages_all,i] -
        hccr_ic*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IC,i] -
        deaths[index$ages_all,IC,i] + migrants[index$ages_all,IC,i]

      # HBeAg-negative CHB
      dpop[index$ages_all,ENCHB,i] <- -(diff(c(0,pop[index$ages_all,ENCHB,i]))/da) +
        pr_ir_enchb_function * pop[index$ages_all,IR,i] +
        pr_ic_enchb * pop[index$ages_all,IC,i] -
        cirrhosis_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB,i] -
        hccr_enchb*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB,i] -
        deaths[index$ages_all,ENCHB,i] + migrants[index$ages_all,ENCHB,i]

      # Compensated cirrhosis
      dpop[index$ages_all,CC,i] <- -(diff(c(0,pop[index$ages_all,CC,i]))/da) +
        cirrhosis_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB,i] +
        pr_ir_cc_function * pop[index$ages_all,IR,i] -
        dcum_dcc[index$ages_all,i] -
        hccr_cc*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CC,i] -
        mu_cc * pop[index$ages_all,CC,i] -
        deaths[index$ages_all,CC,i] + migrants[index$ages_all,CC,i]

      # Decompensated cirrhosis
      dpop[index$ages_all,DCC,i] <- -(diff(c(0,pop[index$ages_all,DCC,i]))/da) +
        dcum_dcc[index$ages_all,i] -
        dflow_dcc_to_hcc[index$ages_all,i] -
        mu_dcc * pop[index$ages_all,DCC,i] -
        deaths[index$ages_all,DCC,i] + migrants[index$ages_all,DCC,i]

      # HCC
      dpop[index$ages_all,HCC,i] <- -(diff(c(0,pop[index$ages_all,HCC,i]))/da) +
        hccr_it*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IT,i] +
        hccr_ir*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IR,i] +
        hccr_ic*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,IC,i] +
        hccr_enchb*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,ENCHB,i] +
        hccr_cc*cancer_prog_rates[index$ages_all,i] * pop[index$ages_all,CC,i] +
        dflow_dcc_to_hcc[index$ages_all,i] -
        dcum_hcc_deaths[index$ages_all,i] -
        deaths[index$ages_all,HCC,i] + migrants[index$ages_all,HCC,i]

      # Immunes
      dpop[index$ages_all,R,i] <- -(diff(c(0,pop[index$ages_all,R,i]))/da) +
        # vacc_cov * vacc_eff * pop[index$ages_all,S,i] +
        (1-p_chronic) * dcum_infections[index$ages_all,i] +
        sag_loss * pop[index$ages_all,IC,i] -
        deaths[index$ages_all,R,i] + migrants[index$ages_all,R,i]

      # Babies are born susceptible or infected (age group 1)
      dpop[1,S,i] <- dpop[1,S,i] + sex_ratio[i] * dcum_nonchronic_births
      dpop[1,IT,i] <- dpop[1,IT,i] + sex_ratio[i] * dcum_chronic_births
      #dpop[1,R,i] <- dpop[1,R,i] + sex_ratio[i] * (1-p_chronic[1]) * infected_births

      # Vaccination: applied at 0.5 years of age (this only makes sense if
      # age step is 0.5!)
      dpop[2,S,i] <- dpop[2,S,i] - (vacc_cov * vacc_eff * pop[2,S,i])
      dpop[2,R,i] <- dpop[2,R,i] + (vacc_cov * vacc_eff * pop[2,S,i])

    }

    # OUTPUT

    # Sum age-specific number of incident background deaths across infection compartments for output
    dcum_deaths <- apply(deaths,c(1,3),sum)
    # Age-specific number of incident background deaths among liver disease patients (CC, DCC and HCC)
    dcum_background_deaths_ld <- apply(deaths[index$ages_all,CC:HCC,1:2],c(1,3),sum)

    # Return results
    res <- c(dpop, dcum_deaths, dcum_infections, dcum_chronic_infections,
             dcum_births, dcum_infected_births, dcum_chronic_births,
             dcum_hbv_deaths, dcum_hcc_deaths, dcum_eag_loss,
             dcum_sag_loss, dcum_dcc, dcum_hcc, dflow_dcc_to_hcc,
             dcum_background_deaths_ld = dcum_background_deaths_ld)
    list(res)

  })

}

### Model-related functions ----

## Function to interpolate demographic parameters over time - specific to these datasets
timevary_parameters_old <- function(timestep, dataset) {
  # Input datasets are matrices of age-specific mortality rates, birth rate and migration rate for every 5-year period
  res <- rep(0,100)
  for (i in 2:101) {
    res[i] <- spline(x = dataset[,1], y = dataset[,i], xout= timestep)[[2]]
  }
  return(res[-1])
} # for loop instead of apply: 69.36s minimally quicker

## Event function: reset population size to initial (1850) size in 1950
reset_pop_1950 <- function(timestep, pop, parameters){
  with (as.list(pop),{
    pop_to_reset <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))
    initialpop <- array(unlist(init_pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    current_pop <- apply(pop_to_reset,c(1,3),sum)
    pop_increase <- apply(initialpop, c(1,3), sum)/current_pop
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

  # Age-specific function of progression through IT and IR (IT=>IR and IR=>IC)
  # Represented by an exponential function that decreases with age
  parameters$eag_loss_function <- parameters$eag_loss_accelerator * exp(-parameters$eag_loss * ages)

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
sum_pop_by_age <- function(time = out$time, pop_output_file) {
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

  return(as.data.frame(pop_output))
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

  # Infection-related transitions per timestep
  # Delete this comments:
  #cum_infections <- select(output, contains("cum_infections"))
  #cum_chronic_infections <- select(output, contains("cum_chronic_infections"))
  #cum_infected_births <- select(output, contains("cum_infected_births"))
  #cum_chronic_births <- select(output, contains("cum_chronic_births"))
  #cum_hbv_deaths <- select(output, contains("cum_hbv_deaths"))

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
  sus <- data.frame(time = output$time,
                    pop = out_sf + out_sm)   # need to change the column names
  carriers <- data.frame(time = output$time,
                         pop = (out_itf + out_itm +
                                  out_irf + out_irm +
                                  out_icf+out_icm+
                                  out_enchbf+out_enchbm+
                                  out_ccf+out_ccm+
                                  out_dccf+out_dccm+
                                  out_hccf+out_hccm))
  carriers_female <- data.frame(time = output$time,
                                pop = (out_itf+
                                         out_irf+
                                         out_icf+
                                         out_enchbf+
                                         out_ccf+out_dccf+out_hccf))
  carriers_male <- data.frame(time = output$time,
                              pop = (out_itm+
                                       out_irm+
                                       out_icm+
                                       out_enchbm+
                                       out_ccm+out_dccm+out_hccm))
  immune <- data.frame(time = output$time, pop = out_rf + out_rm)
  ever_infected <- data.frame(time = output$time, pop = carriers[,-1] + immune[,-1])
  ever_infected_female <- data.frame(time = output$time, pop = carriers_female[,-1] + out_rf)
  ever_infected_male <- data.frame(time = output$time, pop = carriers_male[,-1] + out_rm)
  eag_positive <- data.frame(time = output$time,
                             pop = (out_itf + out_itm +
                                      out_irf + out_irm))
  eag_positive_female <- data.frame(time = output$time,
                                    pop = (out_itf + out_irf))
  eag_positive_male <- data.frame(time = output$time,
                                    pop = (out_itm + out_irm))

  # Total number in each infection compartment per time step
  infectioncat_total <- data.frame(time = output$time,
                                   sus = apply(sus[,-1], 1, sum),
                                   carriers = apply(carriers[,-1], 1, sum),
                                   immune = apply(immune[,-1], 1, sum),
                                   ever_infected = apply(ever_infected[,-1],1,sum))


  # Calculate number of new cases per timestep from cumulative number output

  # Age-specific HBV incidence from horizontal transmission - for women, men and both (total)
  horizontal_infections_female <- data.frame(time = output$time,
                                             incident_number = calculate_incident_numbers(out_cum_infectionsf))
  names(horizontal_infections_female)[-1] <- sprintf("incident_number%g",ages)

  horizontal_infections_male <- data.frame(time = output$time,
                                           incident_number = calculate_incident_numbers(out_cum_infectionsm))
  names(horizontal_infections_male)[-1] <- sprintf("incident_number%g",ages)

  # Total number of incident infections from horizontal transmission and MTCT per time step
  incident_infections <- data.frame(time = output$time,
                                    horizontal_infections = apply(horizontal_infections_female[,-1], 1, sum) +
                                      apply(horizontal_infections_male[,-1], 1, sum),
                                    infected_births = calculate_incident_numbers(out_cum_infected_births))

  # Age-specific chronic infection incidence from horizontal transmission - for women, men and both (total)
  horizontal_chronic_infections_female <- data.frame(time = output$time,
                                                     incident_number = calculate_incident_numbers(out_cum_chronic_infectionsf))
  names(horizontal_chronic_infections_female)[-1] <- sprintf("incident_number%g",ages)

  horizontal_chronic_infections_male <- data.frame(time = output$time,
                                                   incident_number = calculate_incident_numbers(out_cum_chronic_infectionsm))
  names(horizontal_chronic_infections_male)[-1] <- sprintf("incident_number%g",ages)


  # Total number of incident chronic infections from horizontal transmission and MTCT per time step
  incident_chronic_infections <- data.frame(time = output$time,
                                            horizontal_chronic_infections = apply(horizontal_chronic_infections_female[,-1], 1, sum) +
                                              apply(horizontal_chronic_infections_male[,-1], 1, sum),
                                            chronic_births = calculate_incident_numbers(out_cum_chronic_births))

  # Age-specific number of HBV-related deaths - for women and men
  hbv_deaths_female <- data.frame(time = output$time,
                                  incident_number = calculate_incident_numbers(out_cum_hbv_deathsf))
  names(hbv_deaths_female)[-1] <- sprintf("incident_number%g",ages)

  hbv_deaths_male <- data.frame(time = output$time,
                                incident_number = calculate_incident_numbers(out_cum_hbv_deathsm))
  names(hbv_deaths_male)[-1] <- sprintf("incident_number%g",ages)

  # Total number of HBV deaths per time step
  hbv_deaths <- data.frame(time = output$time,
                           incident_number_female = apply(hbv_deaths_female[,-1], 1, sum),
                           incident_number_male = apply(hbv_deaths_male[,-1], 1, sum))
  hbv_deaths$incident_number_total <- hbv_deaths$incident_number_female + hbv_deaths$incident_number_male

  # Age-specific number of total HCC cases - for women and men
  incident_hcc_female <- data.frame(time = output$time,
                                  incident_number = calculate_incident_numbers(out_cum_hccf))
  names(incident_hcc_female)[-1] <- sprintf("incident_number%g",ages)

  incident_hcc_male <- data.frame(time = output$time,
                                incident_number = calculate_incident_numbers(out_cum_hccm))
  names(incident_hcc_male)[-1] <- sprintf("incident_number%g",ages)

  # Total number of total HCC cases per time step
  incident_hcc <- data.frame(time = output$time,
                           incident_number_female = apply(incident_hcc_female[,-1], 1, sum),
                           incident_number_male = apply(incident_hcc_male[,-1], 1, sum))
  incident_hcc$incident_number_total <- incident_hcc$incident_number_female +
    incident_hcc$incident_number_male


  ## Code demography outputs

  # Population:

  # Age-specific and total (last column) population per time step
  pop_female <- sum_pop_by_age(time = output$time, pop_output_file = out_popf)
  pop_male <- sum_pop_by_age(time = output$time, pop_output_file = out_popm)

  pop <- cbind(time = output$time, pop_female[,-1] + pop_male[,-1])

  # Total female, male and both population per time step
  pop_total <- data.frame(time = output$time,
                          pop_female = apply(pop_female[,-1], 1, sum),
                          pop_male = apply(pop_male[,-1], 1, sum)) %>%
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
  deaths_female <- data.frame(time = output$time,
                              incident_number = calculate_incident_numbers(out_cum_deathsf))
  names(deaths_female)[-1] <- sprintf("incident_number%g",ages)
  deaths_female$incident_number_total <- apply(deaths_female[,-1], 1, sum)

  deaths_male <- data.frame(time = output$time,
                            incident_number = calculate_incident_numbers(out_cum_deathsm))
  names(deaths_male)[-1] <- sprintf("incident_number%g",ages)
  deaths_male$incident_number_total <- apply(deaths_male[,-1], 1, sum)

  deaths_total <- data.frame(time = output$time,
                             deaths = deaths_female[,-1] + deaths_male[,-1])
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
                    "cum_hcc_from_dccf" = rep(0,n_agecat), "cum_hcc_from_dccm" = rep(0,n_agecat),
                    "cum_background_deaths_ldf" =  rep(0,n_agecat),
                    "cum_background_deaths_ldm" =  rep(0,n_agecat))
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


# Check initial population vector is the right size
length(init_pop) == n_infectioncat*n_agecat*2+4*n_agecat+2

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
  # NATURAL HISTORY PROGRESSION RATES
  p_chronic = c(0.89, exp(-0.65*ages[-1]^0.46)),     # Age-dependent probability of chronic carriage. Adapted Edmunds
  pr_it_ir = 0.1,
  pr_ir_ic = 0.05,
  eag_loss_accelerator = 9.5, # 9.5 in Margaret's Ethiopia fit, 19.8873 in Shevanthi's
  eag_loss = 0.1281, # 0.1281 in Margaret's Ethiopia fit, 0.977 in Shevanthi's
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
  #vacc_cov = c(0, 0.92, rep(0,n_agecat-2)),  # vaccine is only applied in 1-year olds, need to get time-varying data
  vacc_cov = 0.92,
  vacc_eff = 0.95,                           # vaccine efficacy
  vacc_introtime = 1991,                     # year of vaccine introduction
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


#b1 = 0.07

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
plot(x = outpath$carriers[,1], y = as.numeric(unlist(outpath$carriers[,2]/outpath$pop[,2]))) # age 0.5
abline(v = 1991)
plot(x = outpath$carriers[,1], y = as.numeric(unlist(outpath$carriers[,3]/outpath$pop[,3]))) # age 1
abline(v = 1991)
plot(x = outpath$carriers[,1], y = as.numeric(unlist(outpath$carriers[,4]/outpath$pop[,4]))) # age 1.5
abline(v = 1991)
plot(x = outpath$carriers[,1], y = as.numeric(unlist(outpath$carriers[,5]/outpath$pop[,5]))) # age 2
abline(v = 1991)
plot(x = outpath$carriers[,1], y = as.numeric(unlist(outpath$carriers[,22]/outpath$pop[,22]))) #
abline(v = 1991)

# Carrier prevalence by age in 1980
plot(ages, outpath$carriers[which(outpath$carriers$time == 1980),-1]/
       outpath$pop[which(outpath$pop$time == 1980),-1], type = "l", ylim = c(0,0.3))
#points(gambia_prevdata$age, gambia_prevdata$edmunds_prev, col = "red")

# anti-HBc prevalence by age in 1980
plot(ages, outpath$ever_infected[which(outpath$ever_infected$time == 1980),-1]/
       outpath$pop[which(outpath$pop$time == 1980),-1], type = "l", ylim = c(0,1))

# HBeAg prevalence in chronic carriers by age in 1980
plot(ages, outpath$eag_positive[which(outpath$eag_positive$time == 1980),-1]/
       outpath$carriers[which(outpath$carriers$time == 1980),-1], type = "l", ylim = c(0,1))

# Carrier prevalence by age in 2015
plot(ages, outpath$carriers[which(outpath$carriers$time == 2015),-1]/
       outpath$pop[which(outpath$pop$time == 2015),-1], type = "l", ylim = c(0,0.3))

# anti-HBc prevalence by age in 1980
plot(ages, outpath$ever_infected[which(outpath$ever_infected$time == 2015),-1]/
       outpath$pop[which(outpath$pop$time == 2015),-1], type = "l", ylim = c(0,1))

# HBeAg prevalence in chronic carriers by age in 2015
plot(ages, outpath$eag_positive[which(outpath$eag_positive$time == 2015),-1]/
       outpath$carriers[which(outpath$carriers$time == 2015),-1], type = "l", ylim = c(0,1))

### MODEL CHECK: DEMOGRAPHY PLOTS ----
outpath <- out

# Are there any negative numbers in the output?
any(unlist(outpath$full_output) < 0)

# Do susceptibles + carriers + immunes = total population?
all.equal(sus[,-1] + carriers[,-1] + immune[,-1], pop[,-1], check.names = FALSE)
all.equal(infectioncat_total$sus + infectioncat_total$carriers + infectioncat_total$immune,
          pop_total$pop_total, check.names = FALSE)

# Do susceptibles + carriers + immunes = total population (age-specific numbers)?
all.equal(outpath$sus[,-1] + outpath$carriers[,-1] + outpath$immune[,-1],
          outpath$pop[,-1], check.names = FALSE)
# Do susceptibles + carriers + immunes = total population (total numbers)?
all.equal(outpath$infectioncat_total$sus + outpath$infectioncat_total$carriers +
            outpath$infectioncat_total$immune,
          outpath$pop_total$pop_total, check.names = FALSE)

## Total population, births, deaths

# Plot total population size over timesteps
plot(outpath$pop_total$time, outpath$pop_total$pop_total,
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
     y = outpath$pop_female[which(outpath$pop_female$time == 1950.5),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1950 - women")
points(x = seq(2,82,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1950"]/5,
       col = "red")
# Male age structure in 1950
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 1950.5),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1950 - men")
points(x = seq(2,82,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "1950"]/5,
       col = "red")
# Female age structure in 1970
plot(x = ages,
     y = outpath$pop_female[which(outpath$pop_female$time == 1970),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1970 - women")
points(x = seq(2,82,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1970"]/5,
       col = "red")
# Male age structure in 1970
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 1970),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1970 - men")
points(x = seq(2,82,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "1970"]/5,
       col = "red")
# Female age structure in 1980
plot(x = ages,
     y = outpath$pop_female[which(outpath$pop_female$time == 1980),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1980 - women")
points(x = seq(2,82,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1980"]/5,
       col = "red")
# Male age structure in 1980
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 1980),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1980 - men")
points(x = seq(2,82,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "1980"]/5,
       col = "red")
# Female age structure in 1990
plot(x = ages,
     y = outpath$pop_female[which(outpath$pop_female$time == 1990),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1990 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1990"]/5,
       col = "red")
# Male age structure in 1990
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 1990),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "1990 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "1990"]/5,
       col = "red")

# Plots 2000-2050
par(mfrow=c(4,2))
# Female age structure in 2000
plot(x = ages,
     y = outpath$pop_female[which(outpath$pop_female$time == 2000),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2000 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2000"]/5,
       col = "red")
# Male age structure in 2000
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 2000),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2000 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2000"]/5,
       col = "red")
# Female age structure in 2010
plot(x = ages,
     y = outpath$pop_female[which(outpath$pop_female$time == 2010),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2010 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2010"]/5,
       col = "red")
# Male age structure in 2010
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 2010),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2010 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2010"]/5,
       col = "red")
# Female age structure in 2020
plot(x = ages,
     y = outpath$pop_female[which(outpath$pop_female$time == 2020),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2020 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2020"]/5,
       col = "red")
# Male age structure in 2020
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 2020),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2020 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2020"]/5,
       col = "red")
# Female age structure in 2050
plot(x = ages,
     y = outpath$pop_female[which(outpath$pop_female$time == 2050),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2050 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2050"]/5,
       col = "red")
# Male age structure in 2050
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 2050),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2050 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2050"]/5,
       col = "red")

# Plots 2080-2100
par(mfrow=c(1,2))
# Female age structure in 2080
plot(x = ages,
     y = outpath$pop_female[which(outpath$pop_female$time == 2080),index$ages_all+1]/da,
     type = "l", xlab = "Age", ylab = "Population", main = "2080 - women", ylim = c(0, 55000))
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2080"]/5,
       col = "red")
# Male age structure in 2080
plot(x = ages,
     y = outpath$pop_male[which(outpath$pop_male$time == 2080),index$ages_all+1]/da,
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

### First attempts at fitting using least squares ----
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
  model_prev_subset_both <- data.frame(time = model_output[[seromarker_num]][model_output[[seromarker_num]]$time %in%
                                                                      prev_dataset$time[prev_dataset$sex == "Mixed"],1],
                                       sex = "Mixed",
                                       prev = model_output[[seromarker_num]][model_output[[seromarker_num]]$time %in%
                                                             prev_dataset$time[prev_dataset$sex == "Mixed"],-1]/
                                         model_output[[seromarker_denom]][model_output[[seromarker_denom]]$time %in% prev_dataset$time[prev_dataset$sex == "Mixed"],-1])


  # For women:
  model_prev_subset_female <- data.frame(time = model_output[[seromarker_num_female]][model_output[[seromarker_num_female]]$time %in%
                                                                               prev_dataset$time[prev_dataset$sex == "Female"],1],
                                       sex = "Female",
                                       prev = model_output[[seromarker_num_female]][model_output[[seromarker_num_female]]$time %in%
                                                                               prev_dataset$time[prev_dataset$sex == "Female"],-1]/
                                         model_output[[seromarker_denom_female]][model_output[[seromarker_denom_female]]$time %in% prev_dataset$time[prev_dataset$sex == "Female"],-1])


  # For men:
  model_prev_subset_male <- data.frame(time = model_output[[seromarker_num_male]][model_output[[seromarker_num_male]]$time %in%
                                                                                        prev_dataset$time[prev_dataset$sex == "Male"],1],
                                         sex = "Male",
                                         prev = model_output[[seromarker_num_male]][model_output[[seromarker_num_male]]$time %in%
                                                                                        prev_dataset$time[prev_dataset$sex == "Male"],-1]/
                                           model_output[[seromarker_denom_male]][model_output[[seromarker_denom_male]]$time %in% prev_dataset$time[prev_dataset$sex == "Male"],-1])


  # Assign all columns the same names to combine into 1 dataframe
  names(model_prev_subset_both) <- c("time", "sex", paste0("prev", index$ages_all))
  names(model_prev_subset_female) <- c("time", "sex", paste0("prev", index$ages_all))
  names(model_prev_subset_male) <- c("time", "sex", paste0("prev", index$ages_all))

  # Combine sex-specific dataframes and turn into long format
  model_prev_subset <- rbind(model_prev_subset_female, model_prev_subset_male, model_prev_subset_both)

  model_prev_subset <- gather(model_prev_subset, key = "age", value = "model_prev", -time, -sex)
  model_prev_subset$age <- ages[as.numeric(gsub("\\D", "", model_prev_subset$age))]  # Assign ages as column

  # Merge with the dataset to fit to
  mapped_output_seromarker <- left_join(prev_dataset, model_prev_subset, by = c("sex", "time", "age"))

  return(mapped_output_seromarker)

}


# Trial function to calculate sum of least squares for overall prevalence at given time point
fit_model_sse <- function(..., default_parameter_list, parms_to_change = list(...),
                          scenario = "vacc", data_to_fit) {

  # Simulation parameters for fitting procedure:
  # Simulation starts in 1960, runs for 60 years
  parameters_for_fit <- generate_parameters(default_parameter_list = default_parameter_list,
                                            parms_to_change = c(parms_to_change,
                                                                sim_starttime = 1960))

  sim <- run_model(sim_duration = 60, init_pop_vector = init_pop_sim,
                   default_parameter_list = parameters_for_fit,
                   parms_to_change = NULL,
                   scenario = "vacc")
  out <- code_model_output(sim)

  # Save population distributions for shadow models
  model_pop1974 <- out$full_output[which(out$time==1974),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1978 <- out$full_output[which(out$time==1978),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1983 <- out$full_output[which(out$time==1983),1:(2*n_infectioncat*n_agecat)+1]
  model_pop1993 <- out$full_output[which(out$time==1993),1:(2*n_infectioncat*n_agecat)+1]
  model_pop2005 <- out$full_output[which(out$time==2005),1:(2*n_infectioncat*n_agecat)+1]
  model_pop2012 <- out$full_output[which(out$time==2012),1:(2*n_infectioncat*n_agecat)+1]

  # Define my prevalence datapoints:
  data <- data_to_fit
  data_antihbc_prev <- data_to_fit$antihbc_prevalence
  data_nat_hist_prev <- data_to_fit$natural_history_prevalence

  # Define my model prediction matching the data to fit to:

  # Fitting to age-specific prevalence data: map matching model output to year and age

  ## Age- and sex-specific HBsAg prevalence:
  mapped_output_hbsag <- map_seromarker_prev(seromarker_num = "carriers",
                                             seromarker_denom = "pop",
                                             prev_dataset = data_to_fit$hbsag_prevalence,
                                             model_output = out)

  ## Age- and sex-specific HBeAg prevalence:
  mapped_output_hbeag <- map_seromarker_prev(seromarker_num = "eag_positive",
                                              seromarker_denom = "carriers",
                                              prev_dataset = data_to_fit$hbeag_prevalence,
                                              model_output = out)

  ## Age- and sex-specific anti-HBc prevalence:
  # THIS IS NOT WORKING PROPERLY: is it because it doesn't stabilise??
  mapped_output_antihbc <- map_seromarker_prev(seromarker_num = "ever_infected",
                                             seromarker_denom = "pop",
                                             prev_dataset = data_to_fit$antihbc_prevalence,
                                             model_output = out)

  return(list(mapped_output_antihbc = mapped_output_antihbc, out = out))

  ## Various natural history prevalence estimates

  # Prepare output vectors
  id_gmb1_2 <- c(ic = 0, ic_it = 0, ir_enchb = 0, cc_dcc = 0, hcc = 0)
  id_1_1_1986 <- c(it = 0, ir = 0, enchb = 0, ic = 0, hcc = 0)
  id_gmb1_1 <- c(ic = 0, ic_it = 0, ir_enchb = 0, cc_dcc = 0, hcc = 0)
  id_1_1_2013 <- c(it = 0, ir = 0, enchb = 0, ic = 0, cc_dcc = 0,
                   ir_enchb_cc_dcc_age1 = 0, ir_enchb_cc_dcc_age2 = 0,
                   ir_enchb_cc_dcc_age3 = 0, ir_enchb_cc_dcc_age4 = 0)

  # GMB1-2
  # Row 1 IC
  id_gmb1_2["ic"] <- (sum(select(sim, starts_with("ICm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
    (sum(out$carriers_male[which(out$time == 2013),(which(ages ==27):which(ages ==35.5))+1]))
  # Row 2 IC and IT
  id_gmb1_2["ic_it"] <- (sum(select(sim, starts_with("ICm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
      sum(select(sim, starts_with("ITm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
    (sum(out$carriers_male[which(out$time == 2013),(which(ages ==27):which(ages ==35.5))+1]))
  # Row 3 IR and ENCHB
  id_gmb1_2["ir_enchb"] <- (sum(select(sim, starts_with("IRm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
      sum(select(sim, starts_with("ENCHBm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
    (sum(out$carriers_male[which(out$time == 2013),(which(ages ==27):which(ages ==35.5))+1]))
  # Row 4 CC and DCC
  id_gmb1_2["cc_dcc"] <-(sum(select(sim, starts_with("CCm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)])+
      sum(select(sim, starts_with("DCCm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
    (sum(out$carriers_male[which(out$time == 2013),(which(ages ==27):which(ages ==35.5))+1]))
  # Row 5 HCC
  id_gmb1_2["hcc"] <-(sum(select(sim, starts_with("HCCm"))[which(sim$time == 2013),which(ages ==27):which(ages ==35.5)]))/
    (sum(out$carriers_male[which(out$time == 2013),(which(ages ==27):which(ages ==35.5))+1]))

  # Study 1-1 in 1986
  denom_1_1_1986 <- sum(select(sim, starts_with("IT"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
    select(sim, starts_with("IR"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
    select(sim, starts_with("IC"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)] +
    select(sim, starts_with("ENCHB"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])

  # Row 6
  id_1_1_1986["it"] <- sum(select(sim, starts_with("IT"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
    denom_1_1_1986

  # Row 7
  id_1_1_1986["ir"] <- sum(select(sim, starts_with("IR"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
    denom_1_1_1986

  # Row 8
  id_1_1_1986["enchb"] <- sum(select(sim, starts_with("ENCHB"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
    denom_1_1_1986

  id_1_1_1986["ic"] <- sum(select(sim, starts_with("IC"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)])/
    denom_1_1_1986

  # Row 9 HCC prevalence in chronic carriers
  id_1_1_1986["hcc"] <-(sum(select(sim, starts_with("HCC"))[which(sim$time == 1986),which(ages ==4.5):which(ages ==21.5)]))/
    (sum(out$carriers[which(out$time == 1986),(which(ages ==4.5):which(ages ==21.5))+1]))

  # GMB1-1
  # Row 10 IC
  id_gmb1_1["ic"] <- (sum(select(sim, starts_with("IC"))[which(sim$time == 2012.5),which(ages ==33):which(ages ==47)]))/
    (sum(out$carriers[which(out$time == 2012.5),(which(ages ==33):which(ages ==47))+1]))
  # Row 11 IC and IT
  id_gmb1_1["ic_it"] <- (sum(select(sim, starts_with("IC"))[which(sim$time == 2012.5),which(ages == 33):which(ages == 47)])+
                           sum(select(sim, starts_with("IT"))[which(sim$time == 2012.5),which(ages == 33):which(ages == 47)]))/
    (sum(out$carriers[which(out$time == 2012.5),(which(ages == 33):which(ages ==47))+1]))
  # Row 12 IR and ENCHB
  id_gmb1_1["ir_enchb"] <- (sum(select(sim, starts_with("IR"))[which(sim$time == 2012.5),which(ages ==33):which(ages ==47)])+
                              sum(select(sim, starts_with("ENCHB"))[which(sim$time == 2012.5),which(ages ==33):which(ages ==47)]))/
    (sum(out$carriers[which(out$time == 2012.5),(which(ages ==33):which(ages ==47))+1]))
  # Row 13 CC and DCC
  id_gmb1_1["cc_dcc"] <-(sum(select(sim, starts_with("CC"))[which(sim$time == 2012.5),which(ages ==33):which(ages ==47)])+
                           sum(select(sim, starts_with("DCC"))[which(sim$time == 2012.5),which(ages ==33):which(ages ==47)]))/
    (sum(out$carriers[which(out$time == 2012.5),(which(ages ==33):which(ages ==47))+1]))
  # Row 14 HCC
  id_gmb1_1["hcc"] <-(sum(select(sim, starts_with("HCC"))[which(sim$time == 2012.5),which(ages ==33):which(ages ==47)]))/
    (sum(out$carriers[which(out$time == 2012.5),(which(ages ==33):which(ages ==47))+1]))

  # Study 1-1 in 2013
  denom_1_1_2013 <- sum(select(sim, starts_with("IT"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("IR"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("IC"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)] +
                          select(sim, starts_with("ENCHB"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])

  # Row 16
  id_1_1_2013["it"] <- sum(select(sim, starts_with("IT"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
    denom_1_1_2013

  # Row 17
  id_1_1_2013["ir"] <- sum(select(sim, starts_with("IR"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
    denom_1_1_2013

  # Row 18
  id_1_1_2013["enchb"] <- sum(select(sim, starts_with("ENCHB"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
    denom_1_1_2013

  # Row 19
  id_1_1_2013["ic"] <- sum(select(sim, starts_with("IC"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])/
    denom_1_1_2013

  # Row 20 CC and DCC prevalence in chronic carriers
  id_1_1_2013["cc_dcc"] <-(sum(select(sim, starts_with("CC"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)])+
                             sum(select(sim, starts_with("DCC"))[which(sim$time == 2013),which(ages ==8):which(ages ==95.5)]))/
    (sum(out$carriers[which(out$time == 2013),(which(ages ==8):which(ages ==95.5))+1]))

  # Row 21 IR+ENCHB+CC+DCC
  num_1_1_2013 <- select(sim, starts_with("IR"))[which(sim$time == 2013),]+select(sim, starts_with("ENCHB"))[which(sim$time == 2013),]+
    select(sim, starts_with("CC"))[which(sim$time == 2013),]+select(sim, starts_with("DCC"))[which(sim$time == 2013),]
  id_1_1_2013["ir_enchb_cc_dcc_age1"] <- sum(num_1_1_2013[,which(ages ==8):which(ages == 29.5)])/
    (sum(out$carriers[which(out$time == 2013),(which(ages ==8):which(ages ==29.5))+1]))
  # Row 22
  id_1_1_2013["ir_enchb_cc_dcc_age2"] <- sum(num_1_1_2013[,which(ages ==30):which(ages == 39.5)])/
    (sum(out$carriers[which(out$time == 2013),(which(ages ==30):which(ages ==39.5))+1]))
  # Row 23
  id_1_1_2013["ir_enchb_cc_dcc_age3"] <- sum(num_1_1_2013[,which(ages ==40):which(ages == 49.5)])/
    (sum(out$carriers[which(out$time == 2013),(which(ages ==40):which(ages ==49.5))+1]))
  # Row 24
  id_1_1_2013["ir_enchb_cc_dcc_age4"] <- sum(num_1_1_2013[,which(ages ==50):which(ages == 95.5)])/
    (sum(out$carriers[which(out$time == 2013),(which(ages ==50):which(ages ==95.5))+1]))

  # Combine all data points with model predictions:
  mapped_output_complete <- rbind(mapped_output_hbsag, mapped_output_hbeag)

  # Calculate sum of least squares
  sse <- sum((mapped_output_complete$value - mapped_output_complete$prev_model)^2)

  # Return relevant info (SSE and the matched datapoints and outputs)
  #res <- list(sse = sse, mapped_output = mapped_output,
  #            carrier_prev_total = sim$infectioncat_total$carriers/sim$pop_total$pop_total)

  # Fitting to infection incidence rate:

  # Force of infection (any infection) in 0.5-8.5 year olds (GMB6)
  # Numerator = cumulative incidence of chronic infections and transitions to immune compartment (= any horizontal infection) over follow-up
  # Denominator = person-time in susceptible compartment
  foi_rate <- (sum(select(sim, starts_with("cum_infections"))[which(sim$time == 1984),(which(ages == 0.5):which(ages == 8.5))]) -
                  sum(select(sim, starts_with("cum_infections"))[which(sim$time == 1980),(which(ages == 0.5):which(ages == 8.5))]))/
    ((sum(select(sim, starts_with("S"))[(which(sim$time == 1980):which(sim$time == 1983.5)),(which(ages == 0.5):which(ages == 8.5))]))*dt)

  # Incidence rate of chronic infections in 0.5-7.5 year olds (GMB7)
  # Numerator = cumulative incidence of chronic infections over follow-up
  # Denominator = person-time in susceptible compartment
  chronic_infection_rate <- (sum(select(sim, starts_with("cum_chronic_infections"))[which(sim$time == 1982),(which(ages == 0.5):which(ages == 7.5))]) -
                               sum(select(sim, starts_with("cum_chronic_infections"))[which(sim$time == 1981),(which(ages == 0.5):which(ages == 7.5))]))/
    ((sum(select(sim, starts_with("S"))[(which(sim$time == 1981):which(sim$time == 1981.5)),(which(ages == 0.5):which(ages == 7.5))]))*dt)

  ## SHADOW MODELS 1 a and b: SHIMAKAWA NATURAL HISTORY COHORT
  # Follow 2 cohorts of chronic carriers - 1 of 0-19 year olds (1a) and
  # one of 20-29 year olds (1b) for 28 years, starting in 1974
  # no one had HCC at baseline and we assume no one had DCC at baseline
  # Compartments: IT, IR, IC, ENCHB, CC
  # Switch off births, migation, betas and vaccination to prevent influx of new carriers
  shadow1a_sim <- run_shadow_model(init_age_from = 0, init_age_to = 19.5, init_sex = "both",
                                   init_compartment_from = 2, init_compartment_to = 6,
                                   shadow_default_parameter_list = parameters_for_fit,
                                   shadow_init = model_pop1974, shadowsim_duration = 29,
                                   shadow_parms_to_change = list(sim_starttime = 1974,
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
                                   shadow_init = model_pop1974, shadowsim_duration = 29,
                                   shadow_parms_to_change = list(sim_starttime = 1974,
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
  shadow1a_hcc_ratef <- sum(select(tail(shadow1a_sim,1), starts_with("cum_incident_hccf")))/
    ((sum(head(shadow1a_out$carriers_female[,-1],-1)) - sum(head(select(shadow1a_out$full_output,
                                                                        starts_with("HCCf")),-1)))*dt)

  # In men:
  shadow1a_hcc_ratem <- sum(select(tail(shadow1a_sim,1), starts_with("cum_incident_hccm")))/
    ((sum(head(shadow1a_out$carriers_male[,-1],-1)) - sum(head(select(shadow1a_out$full_output,
                                                                        starts_with("HCCm")),-1)))*dt)

  # MODEL 1b
  # In women:
  shadow1b_hcc_ratef <- sum(select(tail(shadow1b_sim,1), starts_with("cum_incident_hccf")))/
    ((sum(head(shadow1b_out$carriers_female[,-1],-1)) - sum(head(select(shadow1b_out$full_output,
                                                                        starts_with("HCCf")),-1)))*dt)

  # In men:
  shadow1b_hcc_ratem <- sum(select(tail(shadow1b_sim,1), starts_with("cum_incident_hccm")))/
    ((sum(head(shadow1b_out$carriers_male[,-1],-1)) - sum(head(select(shadow1b_out$full_output,
                                                                      starts_with("HCCm")),-1)))*dt)



  # Total incidence of non-malignant ESLD (DCC) per person-year (both sexes)
  # Numerator = cumulative number of incident DCC cases over follow-up (at last timestep)
  # - cumulative number of transitions from DCC to HCC
  # Denominator = person-timestep at risk * dt (carrier compartments other than DCC and HCC)
  # MODEL 1a:
  shadow1a_dcc_rate <- (sum(select(tail(shadow1a_sim,1), starts_with("cum_incident_dcc")))-
                          sum(select(tail(shadow1a_sim,1), starts_with("cum_hcc_from_dcc"))))/
    (sum(head(shadow1a_out$carriers[,-1],-1)) -
       (sum(select(head(shadow1a_sim,-1), starts_with("HCC")))) -
       (sum(select(head(shadow1a_sim,-1),starts_with("DCC"))))*dt)

  # MODEL 1b:
  shadow1b_dcc_rate <- (sum(select(tail(shadow1b_sim,1), starts_with("cum_incident_dcc")))-
                          sum(select(tail(shadow1b_sim,1), starts_with("cum_hcc_from_dcc"))))/
    (sum(head(shadow1b_out$carriers[,-1],-1)) -
       (sum(select(head(shadow1b_sim,-1), starts_with("HCC")))) -
       (sum(select(head(shadow1b_sim,-1),starts_with("DCC"))))*dt)
  # AVERAGE ACROSS AGE GROUPS:
  shadow1_dcc_rate <- weighted.mean(x = c(shadow1a_dcc_rate, shadow1b_dcc_rate),
                                    w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))

  # Incidence rate of eAg loss per person-year
  # Numerator = cumulative number of cases of eAg loss at last timestep
  # Denominator = person-time spent in IT and IR compartments
  # MODEL 1a
  # In women:
  shadow1a_eag_loss_ratef <- (sum(select(tail(shadow1a_sim,1), starts_with("cum_eag_lossf"))))/
    ((sum(select(head(shadow1a_sim,-1),starts_with("ITf")))+
        sum(select(head(shadow1a_sim,-1),starts_with("IRf"))))*dt)
  # In men:
  shadow1a_eag_loss_ratem <- (sum(select(tail(shadow1a_sim,1), starts_with("cum_eag_lossm"))))/
    ((sum(select(head(shadow1a_sim,-1),starts_with("ITm")))+
        sum(select(head(shadow1a_sim,-1),starts_with("IRm"))))*dt)

  # MODEL 1b
  # In women:
  shadow1b_eag_loss_ratef <- (sum(select(tail(shadow1b_sim,1), starts_with("cum_eag_lossf"))))/
    ((sum(select(head(shadow1b_sim,-1),starts_with("ITf")))+
        sum(select(head(shadow1b_sim,-1),starts_with("IRf"))))*dt)
  # In men:
  shadow1b_eag_loss_ratem <- (sum(select(tail(shadow1b_sim,1), starts_with("cum_eag_lossm"))))/
    ((sum(select(head(shadow1b_sim,-1),starts_with("ITm")))+
        sum(select(head(shadow1b_sim,-1),starts_with("IRm"))))*dt)
  # Mortality rate from any cause (HBV-related deaths + background mortality)
  # This was measured in the whole cohort (no matter where they progressed to)
  # MODEL 1a
  # In women:
  shadow1a_mortality_ratef <- (sum(select(tail(shadow1a_sim,1), starts_with("cum_hbv_deathsf"))) +
                                 sum(select(tail(shadow1a_sim,1), starts_with("cum_deathsf"))))/
    (sum(head(shadow1a_out$pop_female[,-1],-1))*dt)
  # In men:
  shadow1a_mortality_ratem <- (sum(select(tail(shadow1a_sim,1), starts_with("cum_hbv_deathsm"))) +
                                 sum(select(tail(shadow1a_sim,1), starts_with("cum_deathsm"))))/
    (sum(head(shadow1a_out$pop_male[,-1],-1))*dt)
  # MODEL 1b
  # In women:
  shadow1b_mortality_ratef <- (sum(select(tail(shadow1b_sim,1), starts_with("cum_hbv_deathsf"))) +
                                 sum(select(tail(shadow1b_sim,1), starts_with("cum_deathsf"))))/
    (sum(head(shadow1b_out$pop_female[,-1],-1))*dt)
  # In men:
  shadow1b_mortality_ratem <- (sum(select(tail(shadow1b_sim,1), starts_with("cum_hbv_deathsm"))) +
                                 sum(select(tail(shadow1b_sim,1), starts_with("cum_deathsm"))))/
    (sum(head(shadow1b_out$pop_male[,-1],-1))*dt)
  # AVERAGE ACROSS AGE GROUPS
  shadow1_mortality_ratef <- weighted.mean(x = c(shadow1a_mortality_ratef, shadow1b_mortality_ratef),
                                           w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))
  shadow1_mortality_ratem <- weighted.mean(x = c(shadow1a_mortality_ratem, shadow1b_mortality_ratem),
                                           w = c(sum(shadow1a_init_pop), sum(shadow1b_init_pop)))


  # List of outputs to fit to
  shadow1_outputs <- list(shadow1a_eag_loss_ratef = shadow1a_eag_loss_ratef,
                          shadow1a_eag_loss_ratem = shadow1a_eag_loss_ratem,
                          shadow1b_eag_loss_ratef = shadow1b_eag_loss_ratef,
                          shadow1b_eag_loss_ratem  = shadow1b_eag_loss_ratem,
                          shadow1a_hcc_ratef = shadow1a_hcc_ratef,
                          shadow1a_hcc_ratem = shadow1a_hcc_ratem,
                          shadow1b_hcc_ratef = shadow1b_hcc_ratef,
                          shadow1b_hcc_ratem = shadow1b_hcc_ratem,
                          shadow1_dcc_rate = shadow1_dcc_rate,
                          shadow1_mortality_ratef = shadow1_mortality_ratef,
                          shadow1_mortality_ratem = shadow1_mortality_ratem)

  ## SHADOW MODEL 2: COURSAGET CHRONIC CARRIER COHORT
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
  shadow2_sag_loss_rate <- (sum(select(tail(shadow2_sim,1), starts_with("cum_sag_loss"))))/
    (sum(select(head(shadow2_sim,-1),starts_with("IC")))*dt)

  # SHADOW MODEL 3: OLUBUYIDE
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
  shadow3_mortality_rate <- (sum(select(tail(shadow3_sim,1), starts_with("cum_hbv_deaths"))) +
                               sum(select(tail(shadow3_sim,1), starts_with("cum_background_deaths_ld"))))/
    ((sum(select(head(shadow3_sim,-1), starts_with("CC"))) + sum(select(head(shadow3_sim,-1), starts_with("DCC"))) +
       sum(select(head(shadow3_sim,-1), starts_with("HCC"))))*dt)

  # SHADOW MODEL 4: SHIMAKAWA COMPENSATED CIRRHOSIS COHORT
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
  shadow4_cum_mortality <- (sum(select(shadow4_sim, starts_with("cum_hbv_deaths"))[which(shadow4_sim$time == 2012.5),]) +
      sum(select(shadow4_sim, starts_with("cum_background_deaths_ld"))[which(shadow4_sim$time == 2012.5),]))/
    (sum(select(shadow4_sim, starts_with("CC"))[which(shadow4_sim$time == 2012),]))

  # SHADOW MODEL 5: YANG HCC COHORT
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

  shadow5_cum_mortality <- c(0,0,0)  # define vector for 3 timesteps
  for (i in 1:3) {  # i = timestep
  shadow5_cum_mortality[i] <- ((sum(select(shadow5_sim, starts_with("cum_hcc_deaths"))[(i+1),])) +
    (sum(select(shadow5_sim, starts_with("cum_background_deaths_ld"))[(i+1),])))/
    sum(select(shadow5_sim, starts_with("HCC"))[which(shadow5_sim$time == 2012),])
  }

  # SHADOW MODEL 6: DIARRA CIRRHOSIS COHORT
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
  shadow6_cum_mortality <- c(0,0)  # define vector for 2 timesteps
  for (i in 1:2) {  # i = timestep
    shadow6_cum_mortality[i] <- ((sum(select(shadow6_sim, starts_with("cum_hbv_deaths"))[(i+1),])) +
                                   (sum(select(shadow6_sim, starts_with("cum_background_deaths_ld"))[(i+1),])))/
      (sum(select(shadow6_sim, starts_with("CC"))[which(shadow6_sim$time == 2005),]) +
         sum(select(shadow6_sim, starts_with("DCC"))[which(shadow6_sim$time == 2005),]))
  }

  # Cumulative HCC incidence after 0.5 and 1 years
  # Numerator = sum of icnident HCC cases at given timestep
  # Can use HCC cases coming from any compartment because there are no people in the other
  # compartments to transition
  # Denominator = number in CC and DCC compartments at t0
  shadow6_cum_hcc <- c(0,0)  # define vector for 2 timesteps
  for (i in 1:2) {  # i = timestep index, so 1 = t0
    shadow6_cum_hcc[i] <- (sum(select(shadow6_sim, starts_with("cum_incident_hcc"))[(i+1),]))/
      (sum(select(shadow6_sim, starts_with("CC"))[which(shadow6_sim$time == 2005),]) +
         sum(select(shadow6_sim, starts_with("DCC"))[which(shadow6_sim$time == 2005),]))
  }

  # SHADOW MODEL 7: GLOBOCAN SURVIVAL CURVE
  # To fit survival curve at time intervals 1, 3 and 5
  # Follow a cohort of HCC patients, starting in 1993, until 1997.5
  # Switch off births, migration, betas and vaccination
  # Since there is no one in other chronic carrier compartments, no need to switch of transitions
  shadow7_sim <- run_shadow_model(init_age_from = 0, init_age_to = 99.5, init_sex = "both",
                                  init_compartment_from = 8, init_compartment_to = 8,
                                  shadow_default_parameter_list = parameters_for_fit,
                                  shadow_init = model_pop1993, shadowsim_duration = 5.5,
                                  shadow_parms_to_change = list(sim_starttime = 1993,
                                                                births_on = 0,
                                                                migration_on = 0,
                                                                b1 = 0,
                                                                b2 = 0,
                                                                b3 = 0))


  # Calculate output for fitting:
  # Cumulative mortality probability (from HCC or background) after 1, 3 and 5 years
  # Numerator = sum of incident deaths at given timestep
  # Can use background deaths from all LD patients because CC and DCC compartments are empty at baseline
  # Denominator = number in HCC compartment at t0
    shadow7_cum_mortality1 <- ((sum(select(shadow7_sim, starts_with("cum_hcc_deaths"))[which(shadow7_sim$time == 1994),])) +
                                   (sum(select(shadow7_sim, starts_with("cum_background_deaths_ld"))[which(shadow7_sim$time == 1994),])))/
      (sum(select(shadow7_sim, starts_with("HCC"))[which(shadow7_sim$time == 1993),]))

    shadow7_cum_mortality3 <- ((sum(select(shadow7_sim, starts_with("cum_hcc_deaths"))[which(shadow7_sim$time == 1996),])) +
        (sum(select(shadow7_sim, starts_with("cum_background_deaths_ld"))[which(shadow7_sim$time == 1996),])))/
      (sum(select(shadow7_sim, starts_with("HCC"))[which(shadow7_sim$time == 1993),]))

    shadow7_cum_mortality5 <- ((sum(select(shadow7_sim, starts_with("cum_hcc_deaths"))[which(shadow7_sim$time == 1998),])) +
                                 (sum(select(shadow7_sim, starts_with("cum_background_deaths_ld"))[which(shadow7_sim$time == 1998),])))/
      (sum(select(shadow7_sim, starts_with("HCC"))[which(shadow7_sim$time == 1993),]))

    shadow7_cum_mortality <- c(shadow7_cum_mortality1, shadow7_cum_mortality3, shadow7_cum_mortality5)


  # Incorporate test: incident infections and births need to be 0!
  return(shadow7_cum_mortality)

}



load(here("data/simulated_inits_1960.RData"))  # this is saved from previous model run
init_pop_sim <- c("Sf" = select(model_pop1960, starts_with("Sf")),
                  "ITf" = select(model_pop1960, starts_with("ITf")),
                  "IRf" = select(model_pop1960, starts_with("IRf")),
                  "ICf" = select(model_pop1960, starts_with("ICf")),
                  "ENCHBf" = select(model_pop1960, starts_with("ENCHBf")),
                  "CCf" = select(model_pop1960, starts_with("CCf")),
                  "DCCf" = select(model_pop1960, starts_with("DCCf")),
                  "HCCf" = select(model_pop1960, starts_with("HCCf")),
                  "Rf" = select(model_pop1960, starts_with("Rf")),
                  "Sm" = select(model_pop1960, starts_with("Sm")),
                  "ITm" = select(model_pop1960, starts_with("ITm")),
                  "IRm" = select(model_pop1960, starts_with("IRm")),
                  "ICm" = select(model_pop1960, starts_with("ICm")),
                  "ENCHBm" = select(model_pop1960, starts_with("ENCHBm")),
                  "CCm" = select(model_pop1960, starts_with("CCm")),
                  "DCCm" = select(model_pop1960, starts_with("DCCm")),
                  "HCCm" = select(model_pop1960, starts_with("HCCm")),
                  "Rm" = select(model_pop1960, starts_with("Rm")),
                  output_storage)
init_pop_sim <- unlist(init_pop_sim)


# Define my datapoints to fit to: HBsAg prevalence dataset
input_hbsag_dataset <- read.csv(here("data",
                                     "hbsag_prevalence.csv"),
                             header = TRUE, check.names = FALSE,
                             stringsAsFactors = FALSE)
input_antihbc_dataset <- read.csv(here("data",
                                     "antihbc_prevalence.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)
input_hbeag_dataset <- read.csv(here("data",
                                     "hbeag_prevalence.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)
input_natural_history_prev_dataset <- read.csv(here("data",
                                     "natural_history_prevalence.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)


prevalence_datasets_list <- list(hbsag_prevalence = input_hbsag_dataset,
                                 antihbc_prevalence = input_antihbc_dataset,
                                 hbeag_prevalence = input_hbeag_dataset,
                                 natural_history_prevalence = input_natural_history_prev_dataset)


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
# get no fits if I do 100 simulations from U(0,1) for b1 and b2

# Run without parallelising
time1 <- proc.time()
out_mat <- apply(params_mat,1,
                    function(x) fit_model_sse(default_parameter_list = parameter_list,
                                              data_to_fit = prevalence_datasets_list,
                                              parms_to_change = list(b1 = as.list(x)$b1,
                                                                     b2 = as.list(x)$b2,
                                                                     mtct_prob_s = as.list(x)$mtct_prob_s)))
sim_duration = proc.time() - time1
sim_duration["elapsed"]/60

# Parallelised code ----
# Set up cluster
cl <- makeCluster(4)
clusterEvalQ(cl, {library(dplyr); library(tidyr); library(deSolve)})
clusterExport(cl, ls())

time1 <- proc.time()
out_mat <- parApply(cl = cl, params_mat,1,
                    function(x) fit_model_sse(default_parameter_list = parameter_list,
                                              data_to_fit = hbsag_dataset_list,
                                              parms_to_change = list(b1 = as.list(x)$b1,
                                                                     b2 = as.list(x)$b2,
                                                                     mtct_prob_s = as.list(x)$mtct_prob_s)))
sim_duration = proc.time() - time1
sim_duration["elapsed"]/60
# Timing: 4.9 min for 20 sim unparallelised, 2.5 min when parallelised,
# 12 min for 100 sims in parallel

# Important: stop cluster!!
stopCluster(cl)

# Matrix of parameter values, model estimates for prevalence in 1980 and 2015, and SSE
#res_mat <- cbind(params_mat, do.call(rbind.data.frame, out_mat_subset)) # this would work for a list
#out_mat_subset <- as.data.frame(t(sapply(out_mat, "[", c("prev_est_1980", "prev_est_2015", "sse"))))
#res_mat <- cbind(params_mat, unnest(out_mat_subset))
out_mat_subset <- sapply(out_mat, "[[", "sse")
res_mat <- cbind(params_mat, sse = out_mat_subset)
res_mat[res_mat$sse == min(res_mat$sse),]

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





