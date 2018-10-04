#################################
### Simple HBV model 04/10/18 ###
#################################
# Version 2: SIR model with demography and age structure
# 4 compartments: Susceptible, Acute infection, Chronic infection, and Immune

#### Load packages ----
require(deSolve)

#### Preparation of simulation parameters and notation ----

## TIMES
stepsize <- 0.1
runtime <- 210 # years
times <- seq(1890, 1890+runtime, stepsize) # Model simulates HBV epidemic from 1890 to 2100

## AGE GROUPS
da <- 1                     # time spent in each age group = 1 year
ages <- seq(1,100,da)       # Ages 1-100
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
adultindex <- 16:100                         # Age groups 16-100 years


#### Functions ----
### The model: specify ODEs
hbv_model <- function(times, pop, parameters){

  # Size of the compartments (creates the vector pop of length = number of compartments)

  S  <-  pop[sindex]  # Susceptible compartment
  A  <-  pop[aindex]  # Acute infection compartment
  I  <-  pop[iindex]  # Chronic infection compartment
  R  <-  pop[rindex]  # Immune compartment
  N  <- S + A + I + R  # Total population  CHECK IF THIS WORKS

  with(as.list(parameters), {  # find values within parameters vector

   # foi <- beta * ((A + alpha*I)/total_pop)         # need to adapt this
   # foi <- beta%*%I

    fluxS <- -diff(c(S[1], S, S[num_agecat]))/da
    fluxA <- -diff(c(A[1], A, A[num_agecat]))/da
    fluxI <- -diff(c(I[1], I, I[num_agecat]))/da
    fluxR <- -diff(c(R[1], R, R[num_agecat]))/da

    # Differential equations
    dS <- - (diff(fluxS)/da) - (foi * S) - (mu * S)
    dA <- - (diff(fluxA)/da) + (foi * S) - (p_chronic * A) - (gamma_acute * A) - (mu * A)
    dI <- - (diff(fluxI)/da) + (p_chronic * A) - (sag_loss * I) - (mu * I) - (mu_hbv * I)
    dR <- - (diff(fluxR)/da) + (gamma_acute * A) + (sag_loss * I) - (mu * R)
    dS[1] <- dS[1] + (b * N)                    # all babies are born susceptible
    dcum_incid  <-  foi * S

    # Return results
    res <-  cbind(dS, dA, dI, dR, dcum_incid)
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


# NEED TO ADAPT EVERYTHING FROM HERE AND ADD WAIFW MATRIX
#### Input ----

## DEMOGRAPHY
# Initial population
init_pop <- c(
  S0 = 180000,
  A0 = 2000,
  I0 = 68000,
  R0 = 250000,
  cum_incid = 0)
# Initial number in each compartment was estimated by output obtained with foi = 0.111
# corresponding to a prevalence of 4%
mu <- 0.014
b <- mu
# mu = background mortality rate is assumed based on 70 years life expectancy
# b = birth rate is assumed to equal the mortality rate

## NATURAL HISTORY (parameterised from Edmunds)
# foi = force of infection, p_chronic = probability of becoming a chronic carrier,
# gamma_acute = rate of recovery from acute infection, sag_loss = rate of HBsAg loss (recovery),
# mu_hbv = rate of HBV-specific deaths,
# rate unit is annual so each time in the model is 1 year
parameters <- c(beta = 9, p_chronic = 0.3, gamma_acute = 4, sag_loss = 0.01,
                alpha = 0.16, mu = mu, mu_hbv = 0.0003, b = b) # foi = 0.111

#### Run the model ----
out <- run_model(init_pop, times, parameters)

### Present output
plot(out[, "time"], out[, "immune"], type = "l", ylim = c(0,500000))
lines(out[, "time"], out[, "susceptible"], col = "green")
lines(out[, "time"], out[, "acute"], col = "red")
lines(out[, "time"], out[, "chronic"], col = "blue")

plot(out[, "time"], out[, "chronic_prevalence"], type = "l", ylim = c(0,0.15))
plot(out[, "time"], out[, "cum_incid"], type = "l")

### Run model checks
# devtools::test()

