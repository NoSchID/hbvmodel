# Process dataset of remaining life expectancy by age, time and sex for DALY calculation

# Read in dataset from Margaret (processed from UN WPP)
library(here)
life_ex <- read.csv(file=here("data-raw/unwpp_2019_remaining_life_expectancy.csv"))
ages <- seq(0,100-0.5,by=0.5)

# Ages from 0 to 100, time from 1950 to 2100
# Interpolate linearly between age groups
life_ex_female <- life_ex[life_ex$gender=="female",] %>%
  select(country, year, age_from, value)

life_ex_male <- life_ex[life_ex$gender=="male",] %>%
  select(country, year, age_from, value)

# Female dataset

# Linear interpolation for missing ages
life_ex_interp_f <- list()
for(i in 1:length(unique(life_ex_female$year))) {
  life_ex_interp_f[[i]] <- approx(x = life_ex_female$age_from[life_ex_female$year==
                                                                   unique(life_ex_female$year)[i]],
                                      y = life_ex_female$value[life_ex_female$year==
                                                                 unique(life_ex_female$year)[i]],
                                      xout = c(ages),
                                      method = "linear",
                                      rule = 2)
}

life_ex_interp_f <- matrix(unlist(lapply(life_ex_interp_f, "[", 2)), ncol = 200, byrow = TRUE)
life_ex_interp_f <- cbind(as.character(unique(life_ex_female$year)), as.data.frame(life_ex_interp_f))
names(life_ex_interp_f) <- c("time", ages)

# Linear interpolation over time between 1950 and 2120
demog_times <- seq(1950, 2120, by = 0.5)
time_interp_f <- apply(life_ex_interp_f, 2, FUN = approx,
                       x = as.numeric(life_ex_interp_f[,1]),
                       xout = demog_times, method = "linear", rule = 2)

time_interp_f  <- do.call("cbind", lapply(time_interp_f, "[[", "y"))
time_interp_f[,1] <- demog_times

# Male dataset

# Linear interpolation for missing ages
life_ex_interp_m <- list()
for(i in 1:length(unique(life_ex_male$year))) {
  life_ex_interp_m[[i]] <- approx(x = life_ex_male$age_from[life_ex_male$year==
                                                                unique(life_ex_male$year)[i]],
                                  y = life_ex_male$value[life_ex_male$year==
                                                             unique(life_ex_male$year)[i]],
                                  xout = c(ages),
                                  method = "linear",
                                  rule = 2)
}

life_ex_interp_m <- matrix(unlist(lapply(life_ex_interp_m, "[", 2)), ncol = 200, byrow = TRUE)
life_ex_interp_m <- cbind(as.character(unique(life_ex_male$year)), as.data.frame(life_ex_interp_m))
names(life_ex_interp_m) <- c("time", ages)

# Linear interpolation over time between 1950 and 2120
demog_times <- seq(1950, 2120, by = 0.5)
time_interp_m <- apply(life_ex_interp_m, 2, FUN = approx,
                       x = as.numeric(life_ex_interp_m[,1]),
                       xout = demog_times, method = "linear", rule = 2)

time_interp_m  <- do.call("cbind", lapply(time_interp_m, "[[", "y"))
time_interp_m[,1] <- demog_times

#write.csv(time_interp_f, file = here("data", "remaining_life_expectancy_female.csv"),
#          row.names=FALSE)
#write.csv(time_interp_m, file = here("data", "remaining_life_expectancy_male.csv"),
#          row.names=FALSE)
