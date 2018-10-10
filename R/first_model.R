#################################
### Simple HBV model 10/10/18 ###
#################################
# Version 2: SIR model with demography and age structure
# 4 compartments: Susceptible, Acute infection, Chronic infection, and Immune

#### Load packages ----
require(deSolve)
library(here)

#### Data ----
pop_size <- read.csv(here("testdata", "gambia_1985_pop.csv"))


#### Preparation of simulation parameters and notation ----

## TIMES
stepsize <- 0.1
runtime <- 210 # years
times <- seq(1890, 1890+runtime, stepsize) # Model simulates HBV epidemic from 1890 to 2100

## AGE GROUPS
da <- 1                     # time spent in each age group = 1 year
ages <- seq(1,80,da)       # Ages 1-80
n_agecat <- length(ages)    # 100 age groups of 1 year

## DEFINITION OF INDICES
# Infection compartments
sindex <- 1:n_agecat                      # Susceptible
aindex <- (n_agecat+1):(2*n_agecat)       # Acute infection
iindex <- (2*n_agecat+1):(3*n_agecat)     # Chronic infection
rindex <- (3*n_agecat+1):(4*n_agecat)     # Immune
# Age groupings
childindex <- 1:5                            # Age groups 1-5 years
juvindex <- 6:15                             # Age groups 6-15 years
adultindex <- 16:80                         # Age groups 16-100 years


#### Functions ----
### The model: specify ODEs
hbv_model <- function(times, pop, parameters){

  # Vector with the 4 compartments (all containing n_agecat age groups)
  S  <-  pop[sindex]          # Susceptible compartment
  A  <-  pop[aindex]          # Acute infection compartment
  I  <-  pop[iindex]          # Chronic infection compartment
  R  <-  pop[rindex]          # Immune compartment
  total_pop_byage  <- S + A + I + R # Total population in each age group
  N <- sum(pop)                     # Total population size

  with(as.list(parameters), {

    foi <- beta %*% ((A + alpha*I)/total_pop_byage)

    # Differential equations
    dS <- - (diff(c(0,S))/da) - (foi * S) - (mu * S)
    dA <- - (diff(c(0,A))/da) + (foi * S) - (p_chronic * A) - (gamma_acute * A) - (mu * A)
    dI <- - (diff(c(0,I))/da) + (p_chronic * A) - (sag_loss * I) - (mu * I) - (mu_hbv * I)
    dR <- - (diff(c(0,R))/da) + (gamma_acute * A) + (sag_loss * I) - (mu * R)

    # Demography: restore additional deaths from last age group as births for constant pop
    dS[1] <- dS[1] + (b * N) + S[n_agecat]/da + A[n_agecat]/da + I[n_agecat]/da + R[n_agecat]/da
    # all babies are born susceptible

    # Output of interest
 #   dcum_incid  <-  foi * S   # check how to add this to output

    # Return results
    res <-  cbind(dS, dA, dI, dR)  #dcum_incid
    list(res)
  })
}

### Function to run the model (solver r4k, with initial conditions and timesteps)
run_model <- function(init_pop, times, parameters) {

  out <- as.data.frame(ode.1D(y = init_pop, times = times, func = hbv_model,
                              parms = parameters, nspec = 4, names = c("S", "A", "I", "R")))

  # Label output table columns   !!!REDO THIS!!!
#  names(out) <- c("time", "susceptible", "acute", "chronic", "immune", "cum_incid")
#  out$total_pop <- (out$susceptible + out$acute + out$chronic + out$immune)
#  out$chronic_prevalence <- (out$chronic/out$total_pop)

  return(out)
}


#### Input ----

## DEMOGRAPHY
# Initial population
init_pop <- c(
  S = c(rep(20400,5), rep(5050,10), rep(1350,30), rep(186,35)),
  A = c(rep(600,5), rep(200,10), rep(50,30), rep(14,35)),
  I = c(rep(5950,5), rep(4375,10), rep(1400,30), rep(160,35)),
  R = c(rep(8050,5), rep(7875,10), rep(6533,30), rep(1640,35))
  )

#init_pop <- c(
#  S = c(rep(20400,5), rep(5050,10), rep(1350,30), rep(186,35)),
#  A = c(rep(600,5), rep(200,10), rep(50,30), rep(14,35)),
#  I = c(rep(5950,5), rep(4375,10), rep(1400,30), rep(160,35)),
#  R = c(rep(8050,5), rep(7875,10), rep(6533,30), rep(1640,35)))
#  cum_incid = 0)

N0 <- sum(init_pop)

mu <- 0.014   # background mortality rate, assumed based on 70 years life expectancy
b <- mu       # birth rate is assumed to equal the mortality rate

## TRANSMISSION
b1 <- 2     # beta-child (up to 5-year olds)
b2 <- 0.5   # beta-young (up to 15-year olds)
b3 <- 0.5   # beta-all (over 5-year olds)
# leads to very quick depletion of susceptibles but high prevalence (nearly 10%)
# alternative: 1,0.1,0.1 maintains susceptible pool but only 5% prevalence

# WAIFW matrix
# Assuming no effective contact between children and adults
beta <- matrix(0, nrow=n_agecat,ncol=n_agecat)     # matrix of transmission rates
beta[childindex,childindex] <- b1                  # transmission among children (1-5 years)
beta[juvindex,juvindex] <- b2                      # transmission among juveniles (6-15 years)
beta[adultindex,adultindex] <- b3                  # transmission among adults (16-100 years)
beta[childindex,juvindex] <- b2                    # transmission from children to juveniles
beta[juvindex, childindex] <- b2                   # transmission from juveniles to children
beta[juvindex, adultindex] <- b3                   # transmission from juveniles to adults
beta[adultindex, juvindex] <- b3                   # transmission from adults to juveniles

## NATURAL HISTORY (parameterised from Edmunds)
# foi = force of infection, p_chronic = probability of becoming a chronic carrier,
# gamma_acute = rate of recovery from acute infection, sag_loss = rate of HBsAg loss (recovery),
# mu_hbv = rate of HBV-specific deaths,
# rate unit is annual so each time in the model is 1 year
parameters <- c(gamma_acute = 4, sag_loss = 0.01,
                alpha = 0.16, mu = mu, mu_hbv = 0, b = b) # mu_hbv =0.0003

# Age-dependent probability of becoming a chronic carrier:
# 0.9 in first group (corresponding to 0 years?)
# 0.4 until age 5 years
# 0.05 for everyone older
p_chronic <- c(0.9, rep(0.4, max(childindex)-1), rep(0.05, n_agecat-max(childindex)))

#### Output ----
out <- run_model(init_pop, times, parameters)

## Tables/vectors with output
time <- out[,1]
sus <- out[,(1+sindex)]
acute <- out[,(1+aindex)]
carriers <- out[,(1+iindex)]
immune <- out[,(1+rindex)]
totalpop <- apply(out[,-1], 1, sum)
agespec_pop <- sus + acute + carriers + immune

## Visualisation of population dynamics/aging
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
plot(1:80, agespec_pop[2000,])

## Visualisation of infection dynamics
# Susceptibles, infectious and immune (all age groups added)
plot(time,apply(carriers,1,sum),type = "l", ylim = c(0,N0))
lines(time,apply(acute,1,sum),col= "green")
lines(time,apply(sus,1,sum),col= "red")
lines(time,apply(immune,1,sum),col= "blue")

# Chronic prevalence by age group over time
plot(time, apply(carriers[,childindex],1,sum), ylim = c(0,N0))
lines(time, apply(carriers[,juvindex],1,sum), col = "red")
lines(time, apply(carriers[,adultindex],1,sum), col = "blue")

# Prevalence by age at 1 time point (after equilibrium is reached)
plot(1:n_agecat, carriers[2000,]/agespec_pop[2000,], ylim = c(0,0.2))

# Proportion immune at 1 time point (after equilibrium is reached)
plot(1:n_agecat, immune[2000,]/agespec_pop[2000,], ylim = c(0,1))

### Run model checks
# devtools::test()

