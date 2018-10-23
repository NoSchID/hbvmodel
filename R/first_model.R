#################################
### Simple HBV model 22/10/18 ###
#################################
# Version 2: SIR model with demography and age structure
# 4 compartments: Susceptible, Acute infection, Chronic infection, and Immune

#### Load packages ----
require(deSolve)
require(here)
library(tictoc)

#### Data preparation ----
gambia_prevdata <- read.csv(here("testdata", "edmunds_gambia_prev.csv"))
gambia_lifeexpectancy <- read.csv(here("testdata", "gambia_1980-1985_lifeexpectancy.csv"))

# Calculate age-specific mortality rates
# Interpolate life expectancy for missing yearly age groups (only available in 5-year groups)
missing_age_exp <- 0:80
life_expectancy <- as.data.frame(approx(x = gambia_lifeexpectancy$age, y = gambia_lifeexpectancy$life_expectancy,
                   xout = missing_age_exp))
mortality_rates_0to80 <- 1/life_expectancy$y
mortality_rates_0to80 <- mortality_rates_0to80[-length(mortality_rates_0to80)]

# Fertility rates of women of childbearing age from UN WPP
gambia_fertility <- data.frame(age = c("15-19", "20-24", "25-29", "30-34", "35-39",
                                             "40-44", "45-49"),
                              fertility_rate = c(0.2013, 0.2822, 0.2760,
                                                 0.2255, 0.1588, 0.0903, 0.0243))
# Assume all women in the 5 year age groups have the same fertility rate
fertility_rates_0to80 <- c(rep(0,15),
                           rep(gambia_fertility$fertility_rate,each = 5),
                           rep(0,30))

# Annual rate (proportion) of HBsAg loss by age from Shimakawa paper
shimakawa_sagloss <- data.frame(age = c("0-9", "10-19", "20-29", "30-39", "40-49",
                                        "50-70"),
                                sagloss_rate = c(0.001, 0.0046, 0.0101,
                                                0.0105, 0.0232, 0.0239))
# Assume everyone in the 10/20 year age groups has the same rate and 70-80 year
# olds have the same rate as 50-70 year olds
sagloss_rates_0to80 <- c(rep(shimakawa_sagloss$sagloss_rate,each = 10),
                         rep(shimakawa_sagloss$sagloss_rate[which(shimakawa_sagloss == "50-70")], 20))

# Prepare data to fit to:
# age-specific prevalence for each year age group (0-79)
edmunds_prev_by_age <- as.data.frame(approx(x = gambia_prevdata$age, y = gambia_prevdata$edmunds_prev,
                                            xout = 0:79))
edmunds_prev_by_age$y <- round(edmunds_prev_by_age$y, 2)
edmunds_prev_by_age$y[c(1,77:80)] <- edmunds_prev_by_age$y[76]
# age-specific proportion ever infected (0-79)
edmunds_everinf_by_age <- as.data.frame(approx(x = gambia_prevdata$age, y = gambia_prevdata$edmunds_prop_ever_infected,
                                               xout = 0:79))
edmunds_everinf_by_age$y <- round(edmunds_everinf_by_age$y, 2)
edmunds_everinf_by_age$y[1] <- 0.10
edmunds_everinf_by_age$y[77:80] <- edmunds_everinf_by_age$y[76]

#### Preparation of simulation parameters and notation ----

## TIMES
stepsize <- 0.1
runtime <- 210 # years
times <- seq(0, 0+runtime, stepsize) # Model simulates HBV epidemic from 1890 to 2100

## AGE GROUPS
da <- 1                     # time spent in each age group = 1 year
ages <- seq(0,80-da,da)     # Age at the beginning of each age group - max age reached is 80
n_agecat <- length(ages)    # 100 age groups of 1 year

## DEFINITION OF INDICES
# Infection compartments
sindex <- 1:n_agecat                      # Susceptible
aindex <- (n_agecat+1):(2*n_agecat)       # Acute infection
iindex <- (2*n_agecat+1):(3*n_agecat)     # Chronic infection
rindex <- (3*n_agecat+1):(4*n_agecat)     # Immune
# Age groupings
childindex <- 1:which(ages == 5)                 # Age groups 0-5 years
juvindex <- which(ages == 6):which(ages == 15)   # Age groups 6-15 years
adultindex <- which(ages == 16):n_agecat         # Age groups 16-80 years

#### Functions ----
### The model: specify ODEs
hbv_model <- function(times, pop, parameters){

  # Set up vector with the 4 compartments (all containing n_agecat age groups)
  S  <-  pop[sindex]          # Susceptible compartment
  A  <-  pop[aindex]          # Acute infection compartment
  I  <-  pop[iindex]          # Chronic infection compartment
  R  <-  pop[rindex]          # Immune compartment
  pop_by_age  <- S + A + I + R # Total population in each age group
  N <- sum(pop)                     # Total population size

  with(as.list(parameters), {

    # Horizontal transmission: WAIFW matrix
    # Assuming no effective contact between children and adults
    beta <- matrix(0, nrow = n_agecat, ncol = n_agecat)# matrix of transmission rates
    beta[childindex,childindex] <- b1                  # transmission among children (1-5 years)
    beta[juvindex,juvindex] <- b2                      # transmission among juveniles (6-15 years)
    beta[adultindex,adultindex] <- b3                  # transmission among adults (16-100 years)
    beta[childindex,juvindex] <- b2                    # transmission from children to juveniles
    beta[juvindex, childindex] <- b2                   # transmission from juveniles to children
    beta[juvindex, adultindex] <- b3                   # transmission from juveniles to adults
    beta[adultindex, juvindex] <- b3                   # transmission from adults to juveniles

    # Force of infection
    foi <- beta %*% ((A + alpha*I)/pop_by_age)

    # Partial differential equations
    dS <- - (diff(c(0,S))/da) - (foi * S) - (mu * S)
    dA <- - (diff(c(0,A))/da) + (foi * S) - (p_chronic * gamma_acute * A) - ((1-p_chronic) * gamma_acute * A) - (mu * A)
    dI <- - (diff(c(0,I))/da) + (p_chronic * gamma_acute * A) - (sag_loss * I) - (mu * I) - (mu_hbv * I)
    dR <- - (diff(c(0,R))/da) + ((1-p_chronic) * gamma_acute * A) + (sag_loss * I) - (mu * R)

    # Demography
    # Putting mortality instead of fertility rates and not dividing by 2 for now to keep pop constant
    infected_births <- ((mtct_prob_a * A + mtct_prob_i * I)) * mortality_rates_0to80
    uninfected_births <- ((pop_by_age) * mortality_rates_0to80) - infected_births
    #  births <- fertility_rates_0to80 * total_pop_byage/2
    #  births <- b * total_pop_byage

    # Births
    # Restore additional deaths from last age group as births for constant population size
    dS[1] <- dS[1] + sum(uninfected_births) + S[n_agecat]/da + A[n_agecat]/da + I[n_agecat]/da + R[n_agecat]/da
    dA[1] <- dA[1] + sum(infected_births)

    # Return results
    res <-  cbind(dS, dA, dI, dR)
    list(res)
  })
}

### Function to run the model
run_model <- function(b1, b2, b3,
                      alpha, gamma_acute, p_chronic,
                      sag_loss, mu_hbv,
                      mtct_prob_a, mtct_prob_i, mu, b) {

  # Add parameters into list
  parameters <- list(b1 = b1, b2 = b2, b3 = b3,
                     alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
                     sag_loss = sag_loss, mu_hbv = mu_hbv,
                     mtct_prob_a = mtct_prob_a, mtct_prob_i = mtct_prob_i, mu = mu, b = b)

  # Run simulation
  out <- as.data.frame(ode.1D(y = init_pop, times = times, func = hbv_model,
                              parms = parameters, nspec = 4, names = c("S", "A", "I", "R")))

  # Code carrier prevalence as output
  pop_by_age <- out[,1+sindex] + out[,1+aindex] + out[,1+iindex] + out[,1+rindex]
  prev_by_age <- out[,1+iindex]/pop_by_age
  prop_everinf_by_age <- (out[,1+aindex] + out[,1+iindex] + out[,1+rindex])/pop_by_age

  # Data of number infected to fit to
  data_prev <- as.numeric(edmunds_prev_by_age$y*pop_by_age[2000,])
  data_everinf <- as.numeric(edmunds_everinf_by_age$y*pop_by_age[2000,])

  # Log likelihood
  LL <- sum(dbinom(x = round(data_prev), size = round(as.numeric(pop_by_age[2000,])),
                   prob = as.numeric(prev_by_age[2000,]), log = TRUE))

  toreturn <- list(out = out, prev_by_age = prev_by_age, loglikelihood = LL)
  #toReturn <- c(modelprev = as.numeric(prev_by_age[2000,]*100), LL = LL)

  return(toreturn)
}

### Functions to return output of interest
return_compartment_output <- function(b1, b2, b3,
                                      alpha, gamma_acute, p_chronic,
                                      sag_loss, mu_hbv,
                                      mtct_prob_a, mtct_prob_i, mu, b) {
  temp <- run_model(b1 = b1, b2 = b2, b3 = b3,
                    alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
                    sag_loss = sag_loss, mu_hbv = mu_hbv,
                    mtct_prob_a = mtct_prob_a, mtct_prob_i = mtct_prob_i, mu = mu, b = b)
  out <- temp$out
  return(out) }

loglikelihood_function <- function(parms_to_estimate) {
  temp <- run_model(b1 = parms_to_estimate[1], b2 = parms_to_estimate[2], b3 = parms_to_estimate[3],
                    alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
                    sag_loss = sag_loss, mu_hbv = mu_hbv,
                    mtct_prob_a = mtct_prob_a, mtct_prob_i = mtct_prob_i, mu = mu, b = b)
  LL <- temp$loglikelihood
  return(LL) }


#### Input ----

## DEMOGRAPHY
# Set up initial population
init_pop <- c(
  S = c(rep(17000,6), rep(5050,10), rep(1350,30), rep(186,34)),
  A = c(rep(500,6), rep(200,10), rep(50,30), rep(14,34)),
  I = c(rep(4959,6), rep(4375,10), rep(1400,30), rep(160,34)),
  R = c(rep(6708,6), rep(7875,10), rep(6533,30), rep(1640,34))
)
N0 <- sum(init_pop)

# Mortality and birth input parameters
mu <- mortality_rates_0to80   # background mortality rate
b <- mu                       # birth rate is assumed to equal the mortality rate

## TRANSMISSION PARAMETERS
b1 <- 0.45     # beta-child (up to 5-year olds)
b2 <- 0.01     # beta-young (up to 15-year olds)
b3 <- 0.005    # beta-all (over 5-year olds)

## NATURAL HISTORY PARAMETERS (annual rates parameterised from Edmunds and Shimakawa)
# foi = force of infection, p_chronic = probability of becoming a chronic carrier,
# gamma_acute = rate of recovery from acute infection, sag_loss = rate of HBsAg loss (recovery),
# mu_hbv = rate of HBV-specific deaths, alpha = relative infectiousness of carriers,
# mtct_prob_a = probability of perinatal transmission from acute mother,
# mtct_prob_i = probability of perinatal transmission from carrier mother.
# Age-dependent probability of becoming a chronic carrier: Edmunds approach
# except 0.89 for whole first year instead of just 0.5 years).
alpha <- 0.16
gamma_acute <- 8 # gamma_acute changed manually to fit prevalence, was 4 in Edmunds
p_chronic <- c(0.89, exp(-0.65*ages[-1]^0.46))
sag_loss <- sagloss_rates_0to80 # 0.01, 0.025
mu_hbv <- 0.0003
mtct_prob_a <- 0.711
mtct_prob_i <- 0.109

#### Try fitting transmission parameters to prevalence ----
#beta_guess <- c(0.1, 0.01, 0.01)
#tic()
#optim(fn = loglikelihood_function, par = beta_guess, control = list(fnscale=-1))
#toc()
# estimates were b1 = 0.45, b2 = 0.01, b3 = 0.005, takes 217.79 sec


#### Output ----
out <- return_compartment_output(b1 = b1, b2 = b2, b3 = b3,
                 alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
                 sag_loss = sag_loss, mu_hbv = mu_hbv,
                 mtct_prob_a = mtct_prob_a, mtct_prob_i = mtct_prob_i, mu = mu, b = b)

## Tables/vectors with output
time <- out[,1]
sus <- out[,(1+sindex)]
acute <- out[,(1+aindex)]
carriers <- out[,(1+iindex)]
immune <- out[,(1+rindex)]
ever_infected <- acute + carriers + immune
totalpop <- apply(out[,-1], 1, sum)
agespec_pop <- sus + acute + carriers + immune

## Visualisation of population dynamics and aging
# Evolution of aging in 3 broad groups (children, juveniles, adults)
all_children <- (apply(carriers[,childindex],1,sum)
                 + apply(sus[,childindex],1,sum)
                 + apply(acute[,childindex],1,sum)
                 + apply(immune[,childindex],1,sum))
all_juveniles <- (apply(carriers[,juvindex],1,sum)
                  + apply(sus[,juvindex],1,sum)
                  + apply(acute[,juvindex],1,sum)
                  + apply(immune[,juvindex],1,sum))
all_adults <- (apply(carriers[,adultindex],1,sum)
               + apply(sus[,adultindex],1,sum)
               + apply(acute[,adultindex],1,sum)
               + apply(immune[,adultindex],1,sum))
plot(time, all_children, ylim = c(0,N0))
lines(time, all_juveniles, col = "red")
lines(time, all_adults, col = "blue")

# Age distribution after equilibrium is reached (takes long to stabilise)
plot(1:n_agecat, agespec_pop[1000,])

## Visualisation of infection dynamics
# Susceptibles, infectious and immune (all age groups added)
plot(time,apply(carriers,1,sum),type = "l", ylim = c(0,N0))
lines(time,apply(acute,1,sum),col= "green")
lines(time,apply(sus,1,sum),col= "red")
lines(time,apply(immune,1,sum),col= "blue")

# Chronic carrier numbers by age group over time
plot(time, apply(carriers[,childindex],1,sum), type = "l", ylim = c(0,N0))
lines(time, apply(carriers[,juvindex],1,sum), col = "red")
lines(time, apply(carriers[,adultindex],1,sum), col = "blue")

# Prevalence by age at 1 time point (after equilibrium is reached)
plot(1:n_agecat, carriers[1000,]/agespec_pop[1000,], type = "l", ylim = c(0,0.3))
# Add Gambia data
points(gambia_prevdata$age, gambia_prevdata$edmunds_prev)

# Proportion ever infected at 1 time point (after equilibrium is reached)
plot(1:n_agecat, ever_infected[1000,]/agespec_pop[1000,], type = "l", ylim = c(0,1))
points(gambia_prevdata$age, gambia_prevdata$edmunds_prop_ever_infected)

### Run model checks
# devtools::test()

