#################################
### Simple HBV model 09/10/18 ###
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

  # Size of the compartments (creates the vector pop of length = number of compartments)

  S  <-  pop[sindex]  # Susceptible compartment
  A  <-  pop[aindex]  # Acute infection compartment
  I  <-  pop[iindex]  # Chronic infection compartment
  R  <-  pop[rindex]  # Immune compartment
  total_pop  <- S + A + I + R  # Total population  CHECK IF THIS WORKS

  with(as.list(parameters), {  # find values within parameters vector

   # foi <- beta * ((A + alpha*I)/total_pop)         # need to adapt this
    foi <- beta %*% ((A + alpha*I)/total_pop) # check first that pop size is constant

    # Differential equations
    dS <- - (diff(c(0,S))/da) - (foi * S) - (mu * S)
    dA <- - (diff(c(0,A))/da) + (foi * S) - (p_chronic * A) - (gamma_acute * A) - (mu * A)
    dI <- - (diff(c(0,I))/da) + (p_chronic * A) - (sag_loss * I) - (mu * I) - (mu_hbv * I)
    dR <- - (diff(c(0,R))/da) + (gamma_acute * A) + (sag_loss * I) - (mu * R)

    # Demography: restore additional deaths from last age group and add births
    dS[n_agecat] <- dS[n_agecat] + S[n_agecat]/da
    dA[n_agecat] <- dA[n_agecat] + A[n_agecat]/da
    dI[n_agecat] <- dI[n_agecat] + I[n_agecat]/da
    dR[n_agecat] <- dR[n_agecat] + R[n_agecat]/da
    dS[1] <- dS[1] + (b * N)                         # all babies are born susceptible

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
  R = c(rep(8050,5), rep(7875,10), rep(6533,30), rep(1640,35)))
#  cum_incid = 0)

N <- sum(init_pop)
# Initial number in each compartment was estimated by output obtained with foi = 0.111
# corresponding to a prevalence of 4%

mu <- 0.014
b <- mu
# mu = background mortality rate is assumed based on 70 years life expectancy
# b = birth rate is assumed to equal the mortality rate

## TRANSMISSION
#b1 <- 0.0985   # beta-child
#b2 <- 0.001   # beta-young
#b3 <- 0.001   # beta-all
b1 <- 1   # beta-child
b2 <- 0.1   # beta-young
b3 <- 0.1   # beta-all

# WAIFW matrix
beta <- matrix(0, nrow=n_agecat,ncol=n_agecat) # matrix of transmission rates
beta[childindex,childindex] <- b1                  # transmission among children (1-5 years)
beta[juvindex,juvindex] <- b2                      # transmission among juveniles (6-15 years)
beta[adultindex,adultindex] <- b3                  # transmission among adults (16-100 years)
beta[childindex,juvindex] <- b2                    # transmission from children to juveniles
beta[juvindex, childindex] <- b2                   # transmission from juveniles to children
beta[juvindex, adultindex] <- b3                   # transmission from juveniles to adults
beta[adultindex, juvindex] <- b3                   # transmission from adults to juveniles
# Assuming no effective contact between children and adults

## NATURAL HISTORY (parameterised from Edmunds)
# foi = force of infection, p_chronic = probability of becoming a chronic carrier,
# gamma_acute = rate of recovery from acute infection, sag_loss = rate of HBsAg loss (recovery),
# mu_hbv = rate of HBV-specific deaths,
# rate unit is annual so each time in the model is 1 year
parameters <- c(gamma_acute = 4, sag_loss = 0.01,
                alpha = 0.16, mu = mu, mu_hbv = 0.0003, b = b) # foi = 0.111, beta = 9

# Age-dependent probability of becoming a chronic carrier
# 0.9 in first group (corresponding to 0 years?)
# 0.4 until age 5 years
# 0.05 for everyone older
p_chronic <- c(0.9, rep(0.4, max(childindex)-1), rep(0.05, n_agecat-max(childindex)))
## is this working? calculate prev proportion

#### Run the model ----
out <- run_model(init_pop, times, parameters)

time <- out[,1]
sus <- out[,(1+sindex)]
acute <- out[,(1+aindex)]
carriers <- out[,(1+iindex)]
immune <- out[,(1+rindex)]
total_pop <- apply(out[,-1], 1, sum)

# Susceptibles, infectious and immune (all age groups added)
plot(time,apply(carriers,1,sum),type = "l", ylim = c(0,N))
lines(time,apply(acute,1,sum),col= "green")
lines(time,apply(sus,1,sum),col= "red")
lines(time,apply(immune,1,sum),col= "blue")

# about 13% prevalence with beta: 1,0.1,0.1

# Evolution of aging
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
plot(time, all_children, ylim = c(0,N))
lines(time, all_juveniles, col = "red")
lines(time, all_adults, col = "blue")

# Prevalence by age group over time
plot(time, apply(carriers[,childindex],1,sum), ylim = c(0,N))
lines(time, apply(carriers[,juvindex],1,sum), col = "red")
lines(time, apply(carriers[,adultindex],1,sum), col = "blue")

# Prevalence by age at one time point
agespec_pop <- c(T= c(init_pop[sindex] + init_pop[aindex] + init_pop[iindex] + init_pop[rindex]))
plot(1:80, carriers[1300,]/agespec_pop, ylim = c(0,1))
# there is a problem with the buildup of carriers in the last age category!


### Run model checks
# devtools::test()

