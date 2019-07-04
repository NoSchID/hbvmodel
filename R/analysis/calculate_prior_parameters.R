# Truncated normal distribution
library(truncnorm)
hist(rtruncnorm(100000, a=1, b=Inf, mean = 6, sd = 3))

# Gamma distribution
# check at https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html
mode = 0.5
min = 0
max = 5
sd = (max - min)/4
# Here are the corresponding rate and shape parameter values:
ra = (mode + sqrt( mode^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
sh = 1 + mode * ra
show(sh)
show(ra)

hist(rgamma(10000,sh,ra))
pgamma(0.1, sh, ra, lower.tail = TRUE) # prob lower than 0.005
pgamma(0.2, sh, ra, lower.tail= FALSE) # prob higher than 0.08

# Beta distribution
# check at https://homepage.divms.uiowa.edu/~mbognar/applets/beta.html
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

mean = 0.95
min = 0.7
max = 1
sd = (max - min)/4
variance = sd^2
variance < (mean*(1-mean))  # function only works if var < mean*(1-mean)

alpha <- estBetaParams(mean, variance)$alpha
beta <- estBetaParams(mean, variance)$beta
alpha
beta

pbeta(0.8, alpha, beta, lower.tail = TRUE) # prob of getting a value lower than x
pbeta(0.95, alpha, beta, lower.tail= FALSE) # prob of getting a value higher than x
