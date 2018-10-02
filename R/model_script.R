#################################
### Simple HBV model 01/10/18 ###
#################################
# Version 1: SIR model with demography
# 4 compartments: Susceptible, Acute infection, Chronic infection, and Immune

#### Load packages ----
require(deSolve)

#### Functions ----
### The model: specify ODEs
hbv_model <- function(times, pop, parameters){

  # Size of the compartments (creates the vector pop of length = number of compartments)

  S  <-  pop[1] # Susceptible compartment
  A  <-  pop[2] # Acute infection compartment
  I  <-  pop[3] # Chronic infection compartment
  R  <-  pop[4] # Immune compartment
  total_pop <- S + A + I + R  # total population size

  with(as.list(parameters), {  # this tells R to look into parameters vector to find values

    foi <- beta * ((A + alpha*I)/total_pop)

    # Differential equations
    dS <- (-foi * S) - (mu * S) + (b * total_pop)
    dA <- (foi * S) - (p_chronic * A) - (gamma_acute * A) - (mu * A)
    dI <- (p_chronic * A) - (sag_loss * I) - (mu * I) - (mu_hbv * I)
    dR <- (gamma_acute * A) + (sag_loss * I) - (mu * R)
    dcum_incid  <-  foi * S

    # Return results
    res <-  cbind(dS, dA, dI, dR, dcum_incid)
    list(res)
  })
}

### Function to run the model (solver r4k, with initial conditions and timesteps)
run_model <- function(init_pop, parameters) {

  # Running the Model For Defined Time Period & Stepsize
  stepsize <- 0.25
  runtime <- 210 # years
  times <- seq(1890, 1890+runtime, stepsize) # Model simulates HBV epidemic from 1890 to 2100

  out <- as.data.frame(rk4(y = init_pop, times = times, func = hbv_model, parms = parameters))

  # Label output table columns
  names(out) <- c("time", "susceptible", "acute", "chronic", "immune", "cum_incid")
  out$total_pop <- (out$susceptible + out$acute + out$chronic + out$immune)
  out$chronic_prevalence <- (out$chronic/out$total_pop)

  return(out)
}

#### Specify paramaters ----

## DEMOGRAPHY
# Initial population
init_pop <- c(S0 = 180000, A0 = 2000, I0 = 68000, R0 = 250000, cum_incid = 0)
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
out <- run_model(init_pop, parameters)

### Present output
plot(out[, "time"], out[, "immune"], type = "l", ylim = c(0,500000))
lines(out[, "time"], out[, "susceptible"], col = "green")
lines(out[, "time"], out[, "acute"], col = "red")
lines(out[, "time"], out[, "chronic"], col = "blue")

plot(out[, "time"], out[, "chronic_prevalence"], type = "l", ylim = c(0,0.15))
plot(out[, "time"], out[, "cum_incid"], type = "l")

### Run model checks
# devtools::test()

