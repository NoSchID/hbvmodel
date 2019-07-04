# Shimakawa 2016 sAg loss data by age
# Age groups are: 0-9 years, 10-19, 20-29, 30-39, 40-49, 50-70
# and are represented by their midpoint
sag_loss_data <- data.frame(age = c(5,15,25,35,45,60),
                            rate = c(0.001,0.0046,0.0101,0.0105,0.0232,0.0239),
                            events = c(1,10,24,16,19,15))

sag_loss_data$ci_lower <- sag_loss_data$rate/(exp(1.96/sqrt(sag_loss_data$events)))
sag_loss_data$ci_upper <- sag_loss_data$rate*(exp(1.96/sqrt(sag_loss_data$events)))

# Plot data
plot(sag_loss_data$age, sag_loss_data$ci_upper, col = "red", type = "l", ylim = c(0,0.04),
     lty = "dotted", xlab = "Age", ylab = "Rate of sAg loss (per person-year)")
points(sag_loss_data$age, sag_loss_data$rate, col = "red")
lines(sag_loss_data$age, sag_loss_data$ci_lower, col = "red", lty = "dotted")


# BEST RESULT
# Linear model with intercept fixed at 0
plot(ages, 4.106e-04*ages, type = "l", ylim = c(0,0.05))
lines(sag_loss_data$age, sag_loss_data$ci_upper, col = "red", lty = "dotted")
lines(sag_loss_data$age, sag_loss_data$ci_lower, col = "red", lty = "dotted")
points(sag_loss_data$age, sag_loss_data$rate, col = "red")
# Chose priors empirically - which coefficients fit within 95% CIs while not exceed 0.05?
lines(ages, 5e-04*ages, col = "blue")
lines(ages, 3e-04*ages, col = "blue")

# EXPLORE AND ANALYSE

plot(ages, -1.697e-03+4.513e-04*ages, type = "l")
lines(ages, 4.106e-04*ages, col = "green")
lines(ages, 4.2e-04*ages, col = "blue")
lines(sag_loss_data$age, sag_loss_data$ci_upper, col = "red",
     lty = "dotted", xlab = "Age", ylab = "Rate of sAg loss (per person-year)")
points(sag_loss_data$age, sag_loss_data$rate, col = "red")
lines(sag_loss_data$age, sag_loss_data$ci_lower, col = "red", lty = "dotted")


# Linear model
linear_model <- lm(rate ~ age, data = sag_loss_data)
summary(linear_model)
linear_pred <- predict(linear_model)
# Add straight line to plot
abline(lm(rate ~ age, data = sag_loss_data), col = "blue")
abline(2*-1.697e-03, 2*4.513e-04)
abline(0.5*-1.697e-03, 0.5*4.513e-04)

lines(sag_loss_data$age, sag_loss_data$ci_upper, col = "red", type = "l", lty = "dotted")
points(sag_loss_data$age, sag_loss_data$rate)
lines(sag_loss_data$age, sag_loss_data$ci_lower, col = "red", lty = "dotted")

# Linear model without intercept
linear_model2 <- lm(rate ~ 0 + age, data = sag_loss_data)
summary(linear_model2)
linear_pred <- predict(linear_model2)
# Add straight line to plot
abline(lm(rate ~ 0 + age, data = sag_loss_data), col = "green")
abline(2*-1.697e-03, 2*4.513e-04)
abline(0.5*-1.697e-03, 0.5*4.513e-04)



# Quadratic model
age2 <- sag_loss_data$age^2
quadratic_model <-lm(rate ~ age + age2, data = sag_loss_data)
quad_pred <- predict(quadratic_model)
points(sag_loss_data$age, quad_pred, col = "green")
summary(quadratic_model)

# Exponential model
exponential_model <- lm(log(rate)~ age), data = sag_loss_data)
summary(exponential_model)
exp_pred <- exp(predict(exponential_model))

lines(sag_loss_data$age, exp_pred, col = "brown")
lines(sag_loss_data$age, 0.001518549*exp(0.05407*sag_loss_data$age))
lines(sag_loss_data$age, 0.001518549*exp(0.1*sag_loss_data$age))

plot(sag_loss_data$age, sag_loss_data$ci_upper, col = "red", type = "l", ylim = c(0,0.04), lty = "dotted")
points(sag_loss_data$age, sag_loss_data$rate)
lines(sag_loss_data$age, sag_loss_data$ci_lower, col = "red", lty = "dotted")
lines(sag_loss_data$age, 0.001602709* exp(0.05407*sag_loss_data$age))  # estimated exponential model
lines(sag_loss_data$age, 0.002 * exp(0.09*sag_loss_data$age^0.84), col = "blue")

# Poisson model
pois <- glm(rate ~ age, data = sag_loss_data, family = poisson)
summary(pois)
points(sag_loss_data$age, exp(predict(pois)))

plot(sag_loss_data$age, 0.002 * 1/exp(0.09*sag_loss_data$age^0.84))


log_model <- lm(rate~log(age), data = sag_loss_data)
log_pred <- predict(log_model)
summary(log_model)

plot(sag_loss_data$age, log_pred, ylim = c(0,0.04), type = "l")
points(sag_loss_data$age, sag_loss_data$rate, col = "red")
# The equation is
-0.017455+0.009374*log(sag_loss_data$age)

plot(ages,-0.017455+0.009374*log(ages), ylim = c(0,0.04))
points(sag_loss_data$age, sag_loss_data$rate, col = "red")
lines(sag_loss_data$age, sag_loss_data$ci_lower, col = "red", lty = "dotted")
lines(sag_loss_data$age, sag_loss_data$ci_upper, col = "red", lty = "dotted")

myeq <- sapply(-0.017455+0.009*log(ages), function(x) max(x,0))
plot(ages,myeq, ylim = c(0,0.04))
points(sag_loss_data$age, sag_loss_data$rate, col = "red")
lines(sag_loss_data$age, sag_loss_data$ci_lower, col = "red", lty = "dotted")
lines(sag_loss_data$age, sag_loss_data$ci_upper, col = "red", lty = "dotted")
lines(ages,sapply(-1.697e-03+0.0004513*ages, function(x) max(x,0)))

