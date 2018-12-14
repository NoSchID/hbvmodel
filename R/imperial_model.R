###################################
### Imperial HBV model 12/12/18 ###
###################################
# Model described in Shevanthi's thesis with some adaptations
# Currently only infant vaccination at 1 year of age, no birth dose or treatment


# 4 compartments: Susceptible, Acute infection, Chronic infection, and Immune
# With horizontal transmission (FOI) and MTCT
# With infant vaccine

### Load packages ----
require(tidyr)
require(dplyr)
require(deSolve)
require(tictoc)
require(here)

### Simulation parameters ----
## Country (gambia or senegal)
countryname <- "gambia"

## Times
dt <- 0.1                             # timesteps
starttime <- 1950
runtime <- 65-dt                      # number of years to run the model for
times <- seq(0, runtime, dt)          # vector of all timesteps
times_labels <- times+starttime

## Age groups
da <- 1                               # time spent in each age group
ages <- seq(0,100-da,da)              # vector of all age groups
n_agecat <- length(ages)              # number of age groups

ages_wocba <- seq(15,50-da,da)        # age groups 15-49 years (women of childbearing age)

## Infection compartments
n_infectioncat <- 9                   # Number of infection compartments

## Definition of indices
index <- list("infection" = 1:n_infectioncat,         # index for all infection compartments
              "ages" = 1:n_agecat,                    # index for all age groups
              "ages_wocba" = ages_wocba+1,            # index for age group 15-49 years (women of childbearing age)
              "ages_1to5" = which(ages == 1):which(ages == 5),       # index for age groups 1-5 years
              "ages_6to15" = which(ages == 6):which(ages == 15),     # index for age groups 6-15 years
              "ages_16to100" = which(ages == 16):n_agecat)          # index for age groups 16-100 years

# Age groupings
childindex <- 1:which(ages == 5)                 # Age groups 0-5 years
juvindex <- which(ages == 6):which(ages == 15)   # Age groups 6-15 years
adultindex <- which(ages == 16):n_agecat         # Age groups 16-80 years

### Infection data preparation ----
gambia_prevdata <- read.csv(here("testdata", "edmunds_gambia_prev.csv"), stringsAsFactors = FALSE)
gambia_lifeexpectancy <- read.csv(here("testdata", "gambia_1980-1985_lifeexpectancy.csv"))

# Interpolate prevalence and prop. ever infected
gambia_prevdata$age[gambia_prevdata$age == 0.5] <- 0
gambia_prevdata <- rbind(gambia_prevdata, data.frame(age = 99, edmunds_prev = 0.05, edmunds_prop_ever_infected = 0.95))
gambia_prev <- approx(x = gambia_prevdata$age, y = gambia_prevdata$edmunds_prev, xout = ages)
gambia_prev <- data.frame(age = gambia_prev$x, prev = gambia_prev$y)
gambia_ever_inf <- approx(x = gambia_prevdata$age, y = gambia_prevdata$edmunds_prop_ever_infected, xout = ages)
gambia_ever_inf <- data.frame(age = gambia_ever_inf$x, ever_inf = gambia_ever_inf$y)

# Calculate the number of susceptibles, acutely infected, chronic infected and recovered
gambia_immune <- gambia_ever_inf$ever_inf - gambia_prev$prev
gambia_infected <- gambia_prev$prev
gambia_sus <- 1-gambia_ever_inf$ever_inf
# Fill HBeAg prevalence (in HBsAg-positives) in with data from Shimakawa paper
gambia_eag <- c(0.95, 0.95, 0.9, 0.9, 0.65, 0.65, 0.6, 0.6, 0.6, 0.6,
                0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.2, 0.2,
                0.1, 0.1, rep(0.1, 20), rep (0.05, 60)) # made up last value

# Annual rate (proportion) of HBsAg loss by age from Shimakawa paper
shimakawa_sagloss <- data.frame(age = c("0-9", "10-19", "20-29", "30-39", "40-49",
                                        "50-70"),
                                sagloss_rate = c(0.001, 0.0046, 0.0101,
                                                0.0105, 0.0232, 0.0239))
# Assume everyone in the 10/20 year age groups has the same rate and 70-80 year
# olds have the same rate as 50-70 year olds
sagloss_rates_0to80 <- c(rep(shimakawa_sagloss$sagloss_rate,each = 10),
                         rep(shimakawa_sagloss$sagloss_rate[which(shimakawa_sagloss == "50-70")], 40))

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

## Calculate age-specific HBeAg-loss function (Shevanthi)
eag_loss <- 19.9996 * exp(-1.1078 * ages)

## Calculate age-specific cancer rate progression function (Shevanthi)
cancer_prog_female <- exp(-7) * (ages * 100 * 0.2 + 2 * exp(0.1077 * ages))
cancer_prog_male <- 5 * cancer_prog_female

### Load demographic datasets ----
# Manual data preparation in Excel:
# copy relevant data in a new file
# all datasets need the first column being labelled "time"
# set cells with estimates to "number" format
input_birthrate_data <- read.csv(here("data", countryname, "birthrate.csv"),
                                 stringsAsFactors = FALSE)  # already divided this by 1000

input_migration_data <- read.csv(here("data", countryname, "migration_rates.csv"),
                                 stringsAsFactors = FALSE)

# Age- and sex-specific datasets
input_popsize_female <- read.csv(here("data", countryname, "popsize_female.csv"),
                                 header = TRUE, check.names = FALSE,
                                 stringsAsFactors = FALSE)
input_popsize_male <- read.csv(here("data", countryname, "popsize_male.csv"),
                               header = TRUE, check.names = FALSE,
                               stringsAsFactors = FALSE)
input_mortality_female <- read.csv(here("data", countryname, "mortality_female.csv"),
                                   header = TRUE, check.names = FALSE,
                                   stringsAsFactors = FALSE)
input_mortality_male <- read.csv(here("data", countryname, "mortality_male.csv"),
                                 header = TRUE, check.names = FALSE,
                                 stringsAsFactors = FALSE)

# 5-year survival rates (survival ratio) from abridged life tables for migration rates
input_survivalratio_female <- read.csv(here("data", countryname, "survivalratio_female.csv"),
                                       header = TRUE, check.names = FALSE,
                                       stringsAsFactors = FALSE)
input_survivalratio_male <- read.csv(here("data", countryname, "survivalratio_male.csv"),
                                     header = TRUE, check.names = FALSE,
                                     stringsAsFactors = FALSE)

# Age-specific fertility rates 1950-2015
input_fertility_data <- read.csv(here("data", countryname, "fertility_rates.csv"),
                                 header = TRUE, check.names = FALSE,
                                 stringsAsFactors = FALSE)

# For validation
input_deaths_female <- read.csv(here("data", countryname, "deaths_female.csv"),
                                header = TRUE, check.names = FALSE,
                                stringsAsFactors = FALSE)
input_deaths_male <- read.csv(here("data", countryname, "deaths_male.csv"),
                              header = TRUE, check.names = FALSE,
                              stringsAsFactors = FALSE)
input_births_total <- read.csv(here("data", countryname, "births_total.csv"),
                               header = TRUE, check.names = FALSE,
                               stringsAsFactors = FALSE)


### Functions for demographic data cleaning and preparation ----
# Function to prepare age-specific (columns) mortality rates by broad time period (rows)
prepare_mort_rates <- function(mortality_dataset) {
  # Input dataset is a mortality rate for ages 0, 1-5, 5-10, ..., 80-85, 85-100
  # Add a row for the final age group (99) and assign same value as for final data point (85)
  lastrow <- group_by(mortality_dataset, time) %>%
    filter(age == max(age)) %>%
    mutate(age = replace(age, values = max(ages)))

  mort_rates <- bind_rows(mortality_dataset, lastrow) %>%
    arrange(time, age)

  # Constant interpolation for missing ages
  mort_rates_interp <- 0
  for(i in mort_rates$time) {
    mort_rates_interp[i] <- list(approx(x = mort_rates$age[mort_rates$time == i],
                                        y = mort_rates$mortality_rate[mort_rates$time == i],
                                        xout = c(ages),
                                        method = "constant"))
  }

  out <- matrix(unlist(lapply(mort_rates_interp[-1], "[", 2)), ncol = n_agecat, byrow = TRUE)
  out <- cbind(time = as.character(unique(mort_rates$time)), as.data.frame(out))
  names(out) <- c("time", ages)
  # Output is in wide format (rows = time period, columns = 1-year age groups)

  return(out)
}

# Function to prepare age-specific (columns) fertility rates by broad time period (rows)
prepare_fert_rates  <- function(fertility_dataset){
  # Input dataset is fertility rates for age groups 15-49 years in 5 years over 5-year time periods
  fert_rates <- fertility_dataset %>%
    gather(key = "age", value = "fertility_rate", -time) %>%    # turn into long format
    arrange(time, age) %>%
    mutate(fertility_rate = replace(fertility_rate, values = fertility_rate/1000)) %>% # rates were imported as per 1000
    mutate(age = replace(age, values = as.numeric(strtrim(age, width = 2)))) # change age to the first age in each group

  # Add a row for the final age group (49) and assign same value as for final data point (45)
  lastrow <- group_by(fert_rates, time) %>%
    filter(age == max(age)) %>%
    mutate(age = replace(age, values = max(ages_wocba)))

  fert_rates <- bind_rows(fert_rates, lastrow) %>%
    arrange(time, age)

  # Constant interpolation for missing ages
  fert_rates_interp <- 0
  for(i in fert_rates$time) {
    fert_rates_interp[i] <- list(approx(x = fert_rates$age[fert_rates$time == i],
                                        y = fert_rates$fertility_rate[fert_rates$time == i],
                                        xout = c(ages_wocba),
                                        method = "constant"))
  }

  out <- matrix(unlist(lapply(fert_rates_interp[-1], "[", 2)), ncol = length(ages_wocba), byrow = TRUE)
  out <- cbind(time = as.character(unique(fert_rates$time)), as.data.frame(out))
  names(out) <- c("time", ages_wocba)
  # Output is in wide format (rows = time period, columns = 1-year age groups)

  return(out)
}

## Function to calculate and prepare age- and sex-specific migration rates
# Calculates age-specific number of net migrants using the cohort component forward method
# based on survival ratio from abridged life tables
calculate_migration_rates <- function(survival_dataset, pop_dataset) {
  survival_dataset$survival_ratio <- as.numeric(survival_dataset$survival_ratio)  # turn survival rates into numeric format

  # Prepare a dataset with population size for every 5 years
  popsize_5years <- clean_number_dataset(pop_dataset, type_label = "pop") %>%
    filter(time %in% seq(1950,2015,5))
  popsize_5years$age[popsize_5years$age == "80+"] <- "80-84"
  popsize_5years <- popsize_5years[popsize_5years$age != "85-89" &
                                     popsize_5years$age != "90-94" &
                                     popsize_5years$age != "95-99" &
                                     popsize_5years$age != "100+",]     # remove age groups >84 years

  # Define labels corresponding to age groups in the population size dataset
  age_labels <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
                  "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
                  "70-74", "75-79", "80-84")

  # Drop survival ratio for 0-1 years (assume 0-1 year olds have
  # same survival ratio as 1-5 year olds) and over 85 year olds (for simplicity):
  survivalratio <- survival_dataset[survival_dataset$age != 0 &
                                      survival_dataset$age != 85,]
  # Label age groups and time period for merging
  survivalratio$age <- age_labels
  survivalratio$time <- as.numeric(substr(survivalratio$time, 1,4)) # restrict to first year in period

  # Merge population size and survival rate datasets and
  # calculate number of net migrants for each age group
  migrants <- left_join(popsize_5years, survivalratio, by = c("age", "time")) %>%
    select(-age_interval) %>%
    mutate(survivors = survival_ratio * pop) %>%
    mutate(lagged_survivors = lag(survivors,18)) %>%
    mutate(migrants = pop-lagged_survivors)
  # survivors calculates the number of people in given age group surviving to the next time period
  # = population of age a at time 0 * survival ratio
  # population in age group a+t at time t is population in the next age group at the next time period
  # Then calculate the net number of migrants as the  the difference between actual
  # population at time t and the population at time 0 survived to time t (forward estimation)

  # Further assumptions
  migrants$migrants[migrants$age == "0-4"] <- NA
  migrants$migrants[is.na(migrants$migrants) == TRUE] <- 0 # assume no migrants aged <5 years
  # Reallocate correct time periods to lagged migration risks
  migrants <- migrants[migrants$time != "1950",]
  time_labels <- c("1950-1955", "1955-1960", "1960-1965", "1965-1970", "1970-1975",
                   "1975-1980", "1980-1985", "1985-1990", "1990-1995", "1995-2000",
                   "2000-2005", "2005-2010", "2010-2015")
  migrants$time <- rep(time_labels, each = length(age_labels))

  # Calculate migration rates from risks
  migrants <- mutate(migrants, migration_rate = (-log(1-migrants/pop))/5) %>%
    select(time, age, migration_rate)
  # Migration rates were calculated from data for each 5-year age group and time period

  return(migrants)
}

calculate_migration_rates_under5 <- function(survival_dataset_female, pop_dataset_female) {
  survival_dataset_female$survival_ratio <- as.numeric(survival_dataset_female$survival_ratio)  # turn survival rates into numeric format

  # Prepare a dataset with population size for every 5 years
  popsize_5years <- clean_number_dataset(pop_dataset_female, type_label = "pop") %>%
    filter(time %in% seq(1950,2015,5))
  popsize_5years$age[popsize_5years$age == "80+"] <- "80-84"
  popsize_5years <- popsize_5years[popsize_5years$age != "85-89" &
                                     popsize_5years$age != "90-94" &
                                     popsize_5years$age != "95-99" &
                                     popsize_5years$age != "100+",]     # remove age groups >84 years

  # Define labels corresponding to age groups in the population size dataset
  age_labels <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
                  "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
                  "70-74", "75-79", "80-84")

  # Drop survival ratio for 0-1 years (assume 0-1 year olds have
  # same survival ratio as 1-5 year olds) and over 85 year olds (for simplicity):
  survivalratio <- survival_dataset_female[survival_dataset_female$age != 0 &
                                             survival_dataset_female$age != 85,]
  # Label age groups and time period for merging
  survivalratio$age <- age_labels
  survivalratio$time <- as.numeric(substr(survivalratio$time, 1,4)) # restrict to first year in period

  # Merge population size and survival rate datasets and
  # calculate number of net migrants for each age group
  migrants <- left_join(popsize_5years, survivalratio, by = c("age", "time")) %>%
    select(-age_interval) %>%
    mutate(survivors = survival_ratio * pop) %>%
    mutate(lagged_survivors = lag(survivors,18)) %>%
    mutate(migrants = pop-lagged_survivors)

  # Calculate the child woman ratio for every 5 years
  cwr <- group_by(popsize_5years, time) %>%
    filter(age == "0-4" | age %in% c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44")) %>%
    summarise(cwr = pop[age == "0-4"]/sum(pop[age %in% c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")]))

  popsize_under5 <- filter(popsize_5years, age == "0-4") %>%
    select(pop)

  # Calculate number of migrants for 0-4 year olds using the child woman ratio
  migrants_under5 <- filter(migrants, age %in% c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44")) %>%
    group_by(time) %>%
    summarise(migrants_wocba = sum(migrants)) %>%
    drop_na() %>%
    mutate(time = replace(time, values = seq(1950,2010,5)))  %>%
    mutate(migrants_under5 = as.numeric(migrants_wocba) * cwr$cwr[-length(cwr$cwr)] * 0.25) %>%
    mutate(migration_risk_under5 = migrants_under5/popsize_under5$pop[-length(popsize_under5$pop)]) %>%
    mutate(migration_rate_under5 = -log(1-migration_risk_under5)/5) %>%
    select(time, migration_rate_under5)

  return(as.data.frame(migrants_under5))
}

prepare_migration_rates <- function(survival_dataset, pop_dataset, survival_dataset_female, pop_dataset_female) {
  migrants <- calculate_migration_rates(survival_dataset, pop_dataset) %>%
    mutate(age = replace(age, values = seq(0,84,5)))   # assign first age of age group
  # Migration rates were calculated from data for each 5-year age group and time period

  migrants_under5 <- calculate_migration_rates_under5(survival_dataset_female, pop_dataset_female)
  migrants$migration_rate[migrants$age == 0] <- migrants_under5$migration_rate_under5

  # Prepare age-specific rates for model
  # Add a row for the final age group (99) and assign same value as for final data point (80)
  lastrow <- group_by(migrants, time) %>%
    filter(age == max(age)) %>%
    mutate(age = replace(age, values = max(ages)))

  migration_rates <- bind_rows(migrants, lastrow) %>%
    arrange(time, age)

  # Constant interpolation for missing ages
  migration_rates_interp <- 0
  for(i in migration_rates$time) {
    migration_rates_interp[i] <- list(approx(x = migration_rates$age[migration_rates$time == i],
                                             y = migration_rates$migration_rate[migration_rates$time == i],
                                             xout = c(ages),
                                             method = "constant"))
  }

  out <- matrix(unlist(lapply(migration_rates_interp[-1], "[", 2)), ncol = n_agecat, byrow = TRUE)
  out <- cbind(time = as.character(unique(migration_rates$time)), as.data.frame(out))
  names(out) <- c("time", ages)
  # Output is in wide format (rows = time period, columns = 1-year age groups)

  return(out)

}

clean_number_dataset <- function(dataset, type_label) {
  ## Clean input datasets (delete age groups with missing data and display integers)
  size <- dataset %>%
    gather(key = "age", value = "thousands", -time) %>%     # turn into long format
    arrange(time)
  size$thousands <- as.numeric(size$thousands)    # change numbers into numeric format
  size <- size %>%
    drop_na %>%                                                 # delete age groups with missing values
    mutate(number = as.numeric(thousands) * 1000) %>%          # numbers were imported as 1000s
    select(-thousands)
  names(size) <- c("time", "age", type_label)
  return(size)
}

prepare_popsize <- function(pop_dataset) {
  ## Clean input datasets (delete age groups with missing data and display integers)
  dataset_clean <- clean_number_dataset(pop_dataset, type_label = "pop")

  ## Prepare numbers for each age step
  # Need to distinguish between pre-1990 data (final age group 80+)
  # and from 1990 onwards (final age group 100+)
  # Delete age group 100+
  popsize <- dataset_clean[dataset_clean$age != "100+",]

  # Repeat each population size by the age group interval multiplied by 1/da to get all age steps
  popsize <- popsize[rep(seq_len(nrow(popsize)), each = 5/da),]

  # Divide numbers in 5-year age groups by time spent in that group, assuming uniform distribution
  popsize$pop <- popsize$pop/(5/da)

  # Label ages:
  # For pre-1990 data, this assumes the maximum age of the population is 84 years
  # For data from 1990 onwards, this assumes the maximum age is 99 years
  popsize$age[popsize$time < 1990] <- rep(seq(0,85-da,da),
                                          times = length(unique(popsize$time[popsize$time < 1990])))
  popsize$age[popsize$time >= 1990] <- rep(seq(0,100-da,da),
                                           times = length(unique(popsize$time[popsize$time >= 1990])))
  popsize$age <- as.numeric(popsize$age)   # age needs to be numeric for ordering
  # Extend dataset for every year to the full age range
  popsize <- merge(popsize,
                   data.frame(time = rep(unique(popsize$time), each = n_agecat),
                              age = as.numeric(rep(ages, times = length(unique(popsize$time))))),
                   by = c("time", "age"), all = TRUE)
  # Assign population size of 0 to all newly added age groups (>84 years) with no data
  popsize[is.na(popsize == TRUE)] <- 0

  return(popsize)
}

### Prepare demographic data ----

## Prepare age- and sex-specific mortality rates by broad time period
# Input: central death rate in abridged life tables
mort_rates_female <- prepare_mort_rates(input_mortality_female)
mort_rates_male <- prepare_mort_rates(input_mortality_male)

## Prepare age-specific fertility rates by broad time period
# Input: age-specific fertility rates for age groups 15-50
fert_rates <- prepare_fert_rates(input_fertility_data)

## Prepare age- and sex-specific migration rates by broad time period
# Input: age-specific survival ratio from abridged life tables and annual population size
# by 5-year age group and sex
migration_rates_female <- prepare_migration_rates(input_survivalratio_female, input_popsize_female,
                                                  input_survivalratio_female, input_popsize_female)
migration_rates_male <- prepare_migration_rates(input_survivalratio_male, input_popsize_male,
                                                input_survivalratio_female, input_popsize_female)

## Clean popsize dataset for output comparison
input_popsize_female_clean <- clean_number_dataset(input_popsize_female, type_label = "pop")
input_popsize_male_clean <- clean_number_dataset(input_popsize_male, type_label = "pop")

## Prepare annual age- and sex-specific population size data
# Input: annual population size by 5-year age group and sex
popsize_female <- prepare_popsize(input_popsize_female)
popsize_male <- prepare_popsize(input_popsize_male)

# Extract 1950 population data
popsize_1950 <- left_join(popsize_female[popsize_female$time == "1950",],
                          popsize_male[popsize_female$time == "1950",],
                          by = c("time", "age")) %>%
  select(-time) %>%
  rename(pop_female = pop.x, pop_male = pop.y)

## Datasets for model validations

# Total population size (both sexes) per year
popsize_total <- popsize_female %>%
  mutate(pop = pop + popsize_male$pop) %>%
  group_by(time) %>%
  summarise(pop = sum(pop))

## Prepare total number of deaths for each year 1950-2015 (for model check)
# Input: age- and sex-specific number of deaths by 5-year time period (female and male dataset)
input_deaths_female_clean <- clean_number_dataset(input_deaths_female, type_label = "deaths")
input_deaths_male_clean <- clean_number_dataset(input_deaths_male, type_label = "deaths")

# Total number of deaths (both sexes) per 5 years
deaths_total <- input_deaths_female_clean %>%
  mutate(deaths = deaths + input_deaths_male_clean$deaths) %>%
  group_by(time) %>%
  summarise(deaths = sum(deaths))

# Extract deaths in 1950
deaths_1950 <- left_join(input_deaths_female_clean[input_deaths_female_clean$time == "1950-1955",],
                         input_deaths_male_clean[input_deaths_male_clean$time == "1950-1955",],
                         by = c("time", "age")) %>%
  rename(deaths_female = deaths.x, deaths_male = deaths.y) %>%
  mutate(deaths_female = replace(deaths_female, values = deaths_female/5),  # get deaths for a single year (1950) in 5-year time period
         deaths_male = replace(deaths_male, values = deaths_male/5)) %>%
  select(-time)
deaths_1950 <- deaths_1950[rep(seq_len(nrow(deaths_1950)), each = 5/da),] %>%
  mutate(deaths_female = replace(deaths_female, values = deaths_female/(5/da)),  # get deaths for a single age in 5-year age groups
         deaths_male = replace(deaths_male, values = deaths_male/(5/da))) %>%
  mutate(age = replace(age, values = ages))

# Total number of births per 5 years
births_total <- select(clean_number_dataset(input_births_total, type_label = "births"), -age)

### Model-related functions ----

## Function to interpolate demographic parameters over time
timevary_parameters <- function(timestep, dataset) {
  # Input datasets are age-specific mortality rates, birth rate and migration rate for every 5-year period
  dataset$time <- seq(0+2, runtime+2, 5)    # convert time to number starting from 0, this is specific to the datasets used in the model
  # dataset$time <- strtrim(dataset$time, width = 4)
  lastrow <- dataset[nrow(dataset),]
  lastrow$time <- times[length(times)]
  dataset <- rbind(dataset, lastrow)
  res <- apply(dataset, 2, FUN = spline, x = dataset$time, xout = timestep)
  res <- unlist(lapply(res, "[", "y"))
  res <- res[-1]
  return(res)
} # spline

## The model
imperial_model <- function(timestep, pop, parameters){

  with(as.list(parameters), {

    # Define time-varying parameters
    mortality_rate <- matrix(c(timevary_parameters(timestep, dataset = mort_rates_female),
                               timevary_parameters(timestep, dataset = mort_rates_male)),
                             ncol = 2)    # 2 columns for sex-specific rates
    #migration_rate <- timevary_parameters(timestep, dataset = input_migration_data)
    #birth_rate <- timevary_parameters(timestep, dataset = input_birthrate_data)
    fertility_rate <- timevary_parameters(timestep, dataset = fert_rates)

    migration_rate <- matrix(c(timevary_parameters(timestep, dataset = migration_rates_female),
                               timevary_parameters(timestep, dataset = migration_rates_male)),
                             ncol = 2)    # 2 columns for sex-specific rates

    # Set up population array with infection compartments
    # matrix 1 = females, matrix 2 = males, rows = agesteps, columns = infection compartments
    pop <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    # Label infection compartments
    S <- 1                           # Susceptible
    IT <- 2                           # Chronic infection: immune tolerant
    IR <- 3                           # Chronic infection: immune reactive
    IC <- 4                           # Chronic infection: inactive carrier
    ENCHB <- 5                        # Chronic infection: HBeAg-negative CHB
    CC <- 6                           # Chronic disease: compensated cirrhosis
    DCC <- 7                          # Chronic disease: decompensated cirrhosis
    HCC <- 8                          # Chronic disease: hepatocellular carcinoma
    R <- 9                           # Immune

    # Horizontal transmission: WAIFW matrix
    # Assuming no effective contact between children and adults
    beta <- matrix(0, nrow = n_agecat, ncol = n_agecat) # matrix of transmission rates
    beta[index$ages_1to5, index$ages_1to5] <- b1        # transmission among children (0-5 years)
    beta[index$ages_6to15, index$ages_6to15] <- b2      # transmission among juveniles (6-15 years)
    beta[index$ages_16to100, index$ages_16to100] <- b3  # transmission among adults (16-100 years)
    beta[index$ages_1to5, index$ages_6to15] <- b2       # transmission from children to juveniles
    beta[index$ages_6to15, index$ages_1to5] <- b2       # transmission from juveniles to children
    beta[index$ages_6to15, index$ages_16to100] <- b3    # transmission from juveniles to adults
    beta[index$ages_16to100, index$ages_6to15] <- b3    # transmission from adults to juveniles

    # Vaccination:
    if (timestep < (vacc_introtime-starttime)) {
      vacc_cov = 0
    } else {
      vacc_cov = vacc_cov
    }

    # Initialise output arrays
    dpop <- array(rep(0,2 * n_infectioncat * n_agecat),dim=c(n_agecat,n_infectioncat,2))      # female and male population in each infection comp
    deaths <- array(rep(0,2 * n_infectioncat * n_agecat),dim=c(n_agecat,n_infectioncat,2))    # female and male incident deaths in each infection comp
    migrants <- array(rep(0,2 * n_infectioncat * n_agecat),dim=c(n_agecat,n_infectioncat,2))  # female and male incident migrants in each infection comp
    infections <- matrix(rep(0, 2* n_agecat), ncol = 2, nrow = n_agecat)                      # female and male incident infections

    # Model equations

    #births <- sum(fertility_rate * pop[index$ages_wocba,index$infection,1])   # applying the same age-specific fertility rate to every infection compartment
    infected_births <- sum(p_chronic[1] * fertility_rate * (apply(mtct_prob_e * pop[index$ages_wocba,c(IT,IR),1],1,sum) + apply(mtct_prob_s * pop[index$ages_wocba,IC:HCC,1],1,sum)))    # infected births come from acute and chronic women of childbearing age
    uninfected_births <- sum(fertility_rate * pop[index$ages_wocba,index$infection,1]) - infected_births  # uninfected births = all births - infected babies
    births <- infected_births + uninfected_births

    # Force of infection - same for women and men TO TEST
    foi <- (beta %*% (apply(pop[index$ages,IC:HCC,1:2],1,sum) + apply(alpha * pop[index$ages,c(IT,IR),1:2],1,sum)))/sum(pop[index$ages,index$infection,1:2])
    # gives a vector with force of infection for every age

    for (i in 1:2) {        # i = sex [1 = female, 2 = male]

      # Incident deaths and migrants
      deaths[index$ages,index$infection,i] <- mortality_rate[index$ages,i] * pop[index$ages,index$infection,i] # works: applying same age-specific mortality rate to every infection compartment
      migrants[index$ages,index$infection,i] <-  migration_rate[index$ages,i] * pop[index$ages,index$infection,i]

      # Incident infections: for every sex gives a vector with incident infections for every age, store as 2 columns
      infections[index$ages,i] <- foi * pop[index$ages,S,i]

      # Partial differential equations

      # Susceptibles
      dpop[index$ages,S,i] <- -(diff(c(0,pop[-length(index$ages),S,i],0))/da) -
                            deaths[index$ages,S,i] + migrants[index$ages,S,i] -
                            p_chronic * infections[index$ages,i] -
                            (1-p_chronic) * infections[index$ages,i] -
                            (vacc_cov * vacc_eff * pop[index$ages,S,i])
      # Immune tolerant
      dpop[index$ages,IT,i] <- -(diff(c(0,pop[-length(index$ages),IT,i],0))/da) -
                            deaths[index$ages,IT,i] + migrants[index$ages,IT,i] +
                            p_chronic * infections[index$ages,i] -
                            pr1 * dpop[index$ages,IT,i] -
                            hccr1[index$ages,i] * dpop[index$ages,IT,i]

      # Immune reactive
      dpop[index$ages,IR,i] <- -(diff(c(0,pop[-length(index$ages),IR,i],0))/da) -
        deaths[index$ages,IR,i] + migrants[index$ages,IR,i] +
        pr1 * dpop[index$ages,IT,i] -
        pr2 * dpop[index$ages,IR,i] -
        pr3 * dpop[index$ages,IR,i] -
        hccr2[index$ages,i] * dpop[index$ages,IR,i]

      # Inactive carrier
      dpop[index$ages,IC,i] <- -(diff(c(0,pop[-length(index$ages),IC,i],0))/da) -
        deaths[index$ages,IC,i] + migrants[index$ages,IC,i] +
        pr2 * dpop[index$ages,IR,i] -
        pr4 * dpop[index$ages,IC,i] -
        sag_loss * dpop[index$ages,IC,i] -
        hccr3[index$ages,i] * dpop[index$ages,IC,i]

      # HBeAg-negative CHB
      dpop[index$ages,ENCHB,i] <- -(diff(c(0,pop[-length(index$ages),ENCHB,i],0))/da) -
        deaths[index$ages,ENCHB,i] + migrants[index$ages,ENCHB,i] +
        pr3 * dpop[index$ages,IR,i] +
        pr4 * dpop[index$ages,IC,i] -
        ccrate * dpop[index$ages,ENCHB,i] -
        hccr4[index$ages,i] * dpop[index$ages,ENCHB,i]

      # Compensated cirrhosis
      dpop[index$ages,CC,i] <- -(diff(c(0,pop[-length(index$ages),CC,i],0))/da) -
        deaths[index$ages,CC,i] + migrants[index$ages,CC,i] +
        ccrate * dpop[index$ages,ENCHB,i] -
        dccrate * dpop[index$ages,CC,i] -
        hccr5[index$ages,i] * dpop[index$ages,CC,i] -
        mu_cc * dpop[index$ages,CC,i]

      # Decompensated cirrhosis
      dpop[index$ages,DCC,i] <- -(diff(c(0,pop[-length(index$ages),DCC,i],0))/da) -
        deaths[index$ages,DCC,i] + migrants[index$ages,DCC,i] +
        dccrate * dpop[index$ages,CC,i] -
        hccr6 * dpop[index$ages,DCC,i] -
        mu_dcc * dpop[index$ages,DCC,i]

      # HCC
      dpop[index$ages,HCC,i] <- -(diff(c(0,pop[-length(index$ages),HCC,i],0))/da) -
        deaths[index$ages,HCC,i] + migrants[index$ages,HCC,i] +
        hccr1[index$ages,i] * dpop[index$ages,IT,i] +
        hccr2[index$ages,i] * dpop[index$ages,IR,i] +
        hccr3[index$ages,i] * dpop[index$ages,IC,i] +
        hccr4[index$ages,i] * dpop[index$ages,ENCHB,i] +
        hccr5[index$ages,i] * dpop[index$ages,CC,i] +
        hccr6 * dpop[index$ages,DCC,i] -
        mu_hcc * dpop[index$ages,HCC,i]

      # Immunes
      dpop[index$ages,R,i] <- -(diff(c(0,pop[-length(index$ages),R,i],0))/da) -
        deaths[index$ages,R,i] + migrants[index$ages,R,i] +
        (1-p_chronic) * infections[index$ages,i] +
        sag_loss * dpop[index$ages,IC,i] +
        (vacc_cov * vacc_eff * pop[index$ages,S,i])

      # Babies are born susceptible or acutely infected (age group 1)
      dpop[1,S,i] <- dpop[1,S,i] + sex_ratio[i] * uninfected_births
      dpop[1,IT,i] <- dpop[1,IT,i] + sex_ratio[i] * infected_births

    }

    # Return results
    res <- c(dpop, births)
    #res <- c(dpop, deaths, migrants, births)
    list(res)
  })
}


## Function to run the model
run_model <- function(b1 = b1, b2 = b2, b3 = b3,
                      alpha, gamma_acute, p_chronic,
                      sag_loss, mu_hbv,
                      mtct_prob_e, mtct_prob_s,
                      vacc_cov = vacc_cov, vacc_eff = vacc_eff, vacc_introtime = vacc_introtime) {

  # Add parameters into list
  parameters <- list(b1 = b1, b2 = b2, b3 = b3,
                     alpha = alpha, p_chronic = p_chronic, sag_loss = sag_loss,
                     mtct_prob_e = mtct_prob_e, mtct_prob_s = mtct_prob_s,
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

## Functions to return output of interest
#return_compartment_output <- function(b1, b2, b3,
#                                      alpha, gamma_acute, p_chronic,
#                                      sag_loss, mu_hbv,
#                                      mtct_prob_a, mtct_prob_i, mu, b) {
#  temp <- run_model(b1 = b1, b2 = b2, b3 = b3,
#                    alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
#                    sag_loss = sag_loss, mu_hbv = mu_hbv,
#                    mtct_prob_a = mtct_prob_a, mtct_prob_i = mtct_prob_i, mu = mu, b = b)
#  out <- temp$out
#  return(out) }

#loglikelihood_function <- function(parms_to_estimate) {
#  temp <- run_model(b1 = parms_to_estimate[1], b2 = parms_to_estimate[2], b3 = parms_to_estimate[3],
#                    alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
#                    sag_loss = sag_loss, mu_hbv = mu_hbv,
#                    mtct_prob_a = mtct_prob_a, mtct_prob_i = mtct_prob_i, mu = mu, b = b)
#  LL <- temp$loglikelihood
#  return(LL) }


### Model input ----

## DEMOGRAPHY
# Set up initial population
# Initial population = age- and sex-specific population size in 1950
#init_pop <- c("Sf" = popsize_1950$pop_female, "Sm" = popsize_1950$pop_male,
#              "cum_deathsf" = deaths_1950$deaths_female*dt, "cum_deathsm"t = deaths_1950$deaths_female*dt,
#              "cum_migrantsf" = rep(0, n_agecat), "cum_migrantsm" = rep(0, n_agecat),
#              "cum_births" = 0)

# Note: names in initial population vector is reproduced in output
init_pop <- c("Sf" = popsize_1950$pop_female*gambia_sus,
              "ITf" = popsize_1950$pop_female*gambia_infected*gambia_eag*0.6,
              "IRf" = popsize_1950$pop_female*gambia_infected*gambia_eag*0.4,
              "ICf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.6,
              "ENCHBf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.25,
              "CCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.1,
              "DCCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.04,
              "HCCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.01,
              "Rf" = popsize_1950$pop_female*gambia_immune,
              "Sm" = popsize_1950$pop_male*gambia_sus,
              "ITm" = popsize_1950$pop_male*gambia_infected*gambia_eag*0.6,
              "IRm" = popsize_1950$pop_male*gambia_infected*gambia_eag*0.4,
              "ICm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.6,
              "ENCHBm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.25,
              "CCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.1,
              "DCCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.04,
              "HCCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.01,
              "Rm" = popsize_1950$pop_male*gambia_immune,
              "cum_births" = 0)
# made up percentages in each compartment
N0 <- sum(init_pop[1:(n_infectioncat * n_agecat * 2)])

## TRANSMISSION PARAMETERS
b1 <- 15     # beta-child (up to 5-year olds)
b2 <- 5      # beta-young (up to 15-year olds)
b3 <- 0.95      # beta-all (over 5-year olds)

## NATURAL HISTORY PARAMETERS (annual rates parameterised from Edmunds and Shimakawa)
# foi = force of infection, p_chronic = probability of becoming a chronic carrier,
# gamma_acute = rate of recovery from acute infection, sag_loss = rate of HBsAg loss (recovery),
# mu_hbv = rate of HBV-specific deaths, alpha = relative infectiousness of carriers,
# mtct_prob_a = probability of perinatal transmission from acute mother,
# mtct_prob_i = probability of perinatal transmission from carrier mother.
# Age-dependent probability of becoming a chronic carrier: Edmunds approach
# except 0.89 for whole first year instead of just 0.5 years).
alpha <- 15 # Shevanthi value
p_chronic <- c(0.89, exp(-0.65*ages[-1]^0.46))
#sag_loss <- 0.01 # Shevanthi value
sag_loss <- sagloss_rates_0to80
mtct_prob_e <- 0.9  # Shevanthi value
mtct_prob_s <- 0.05   # Shevanthi model value
sex_ratio <- c(0.4926, 0.5074)
pr1 <- 0.1 * eag_loss
pr2 <- 0.05 * eag_loss
pr3 <- 0.005
pr4 <- 0.01
ccrate <- 0.04
dccrate <- 0.04
hccr1 <- matrix(data = c(cancer_prog_female, cancer_prog_male), nrow = n_agecat, ncol = 2)
hccr2 <- matrix(data = c(2*cancer_prog_female, 2*cancer_prog_male), nrow = n_agecat, ncol = 2)
hccr3 <- matrix(data = c(0.5*cancer_prog_female, 0.5*cancer_prog_male), nrow = n_agecat, ncol = 2)
hccr4 <- matrix(data = c(2*cancer_prog_female, 2*cancer_prog_male), nrow = n_agecat, ncol = 2)
hccr5 <- matrix(data = c(13*cancer_prog_female, 13*cancer_prog_male), nrow = n_agecat, ncol = 2)
hccr6 <- 0.04
mu_cc <- 0.039
mu_dcc <- 0.314
mu_hcc <- mu_dcc
#mu_hcc <- 0.5
# Add these new parameters into run model function

## INTERVENTION PARAMETERS
vacc_cov <- c(0, 0.92, rep(0,n_agecat-2)) # vaccine is only applied in 1-year olds, need to get time-varying data
vacc_eff <- 0.95
vacc_introtime <- 1991 # year of vaccine introduction

### Try fitting transmission parameters to prevalence ----
#beta_guess <- c(0.1, 0.01, 0.01)
#tic()
#optim(fn = loglikelihood_function, par = beta_guess, control = list(fnscale=-1))
#toc()
# estimates were b1 = 0.45, b2 = 0.01, b3 = 0.005, takes 217.79 sec


### Run the model ----
tic()
out <- run_model(b1 = b1, b2 = b2, b3 = b3,
                 alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
                 sag_loss = sag_loss, mu_hbv = mu_hbv,
                 mtct_prob_e = mtct_prob_e, mtct_prob_s = mtct_prob_s,
                 vacc_cov = vacc_cov, vacc_eff = vacc_eff, vacc_introtime = vacc_introtime)
toc()

# DEBUG: Check which compartments are negative
has.neg.col <- apply(out, 2, function(col) any(col < 0))

apply(out, 2, function(x) if(!all(x >= 0)) print(names(x)))



View(has.neg.col[has.neg.col == TRUE])
which(has.neg.col)
# Where do negative numbers appear at first timestep?
out2 <- out[2,]
colnames(out2)[as.numeric(out2) < 0]
# ICf1, ICm1, nearly all HCC


#out <- return_compartment_output(b1 = b1, b2 = b2, b3 = b3,
#                 alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
#                 sag_loss = sag_loss, mu_hbv = mu_hbv,
#                 mtct_prob_a = mtct_prob_a, mtct_prob_i = mtct_prob_i, mu = mu, b = b)

### Output-related functions ----
# Function to sum numbers from different compartments for each age and time step
sum_pop_by_age <- function(time = out$time, output) {
  out <- data.frame(time = time, output) %>%
    gather(key = "agegroup", value = "pop", -time) %>%       # turn into wide format
    arrange(time) %>%                                        # order by timestep
    mutate(agegroup = as.numeric(replace(agegroup,           # remove infection information
                                         values = as.numeric(gsub("\\D", "", agegroup))))) %>%
    group_by(time, agegroup) %>%
    summarise(pop = sum(pop)) %>%                            # sum numbers for each age group at each timestep
    spread(key = "agegroup", value = "pop")                  # return to wide format

  return(out)
}

### Code output ----

## Extract separate outputs

# Infection compartments
out_sf <- select(out, starts_with("Sf"))
out_sm <- select(out, starts_with("Sm"))
out_af <- select(out, starts_with("Af"))
out_am <- select(out, starts_with("Am"))
out_if <- select(out, starts_with("If"))
out_im <- select(out, starts_with("Im"))
out_rf <- select(out, starts_with("Rf"))
out_rm <- select(out, starts_with("Rm"))

# Population
out_popf <- select(out[,2:(n_agecat*n_infectioncat*2+1)], contains("f"))
out_popm <- select(out[,2:(n_agecat*n_infectioncat*2+1)], contains("m"))
out_pop <- cbind(out_popf, out_popm)

## Code infection outputs
# Age-specific number in each infection compartment per time step
model_sus <- data.frame(time = out$time, pop = out_sf + out_sm)   # need to change the column names
model_acute <- data.frame(time = out$time, pop = out_af + out_am)
model_carriers <- data.frame(time = out$time, pop = out_if + out_im)
model_immune <- data.frame(time = out$time, pop = out_rf + out_rm)

# Total number in each infection compartment per time step
model_infectioncat_total <- data.frame(time = out$time,
                                sus = apply(model_sus[,-1], 1, sum),
                                acute = apply(model_acute[,-1], 1, sum),
                                carriers = apply(model_carriers[,-1], 1, sum),
                                immune = apply(model_immune[,-1], 1, sum))
#ever_infected <- acute + carriers + immune

## Code demography outputs

# Population:

# Age-specific and total (last column) population per time step
model_pop_female <- sum_pop_by_age(output = out_popf)
model_pop_male <- sum_pop_by_age(output = out_popm)
model_pop <- sum_pop_by_age(output = out_pop)

# Total female, male and both population per time step
model_pop_total <- data.frame(time = out$time,
                              pop_female = apply(model_pop_female[,-1], 1, sum),
                              pop_male = apply(model_pop_male[,-1], 1, sum)) %>%
  mutate(pop_total = pop_female + pop_male)

# Births:

# Total number of births at each timestep
model_births <- data.frame(time = out$time,
                           births = c(out$cum_births[1], diff(out$cum_births, lag = 1)))
names(model_births) <- c("time", "births")

# Total number of births grouped in 5-year time periods
model_births_group5 <- model_births %>%
  mutate(timegroup = floor(time / 5) * 5) %>%
  group_by(timegroup) %>%
  summarise_all(sum) %>%
  select(-time)

## Plots

### Plots ----

## DEMOGRAPHY

## Total population, births, deaths

# Plot total population size over timesteps
plot(popsize_total$time, popsize_total$pop, col = "red",
     xlab = "Year", ylab = "Population size")
lines(model_pop_total$time, model_pop_total$pop_total)

# Plot total number of births over time periods
plot(x = as.numeric(strtrim(births_total$time, width = 4)),
     y = births_total$births, col = "red",
     xlab = "5-year time periods", ylab = "Total number of births")
lines(model_births_group5$timegroup, model_births_group5$births)

## Age structure

# Plot female age structure in 1970
plot(x = ages,
     y = model_pop_female[201,index$ages+1],
     type = "l", xlab = "Age", ylab = "Population", main = "1970 - women")
points(x = seq(2,82,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1970"]/5,
       col = "red")

# Plot male age structure in 2005
plot(x = ages,
     y = model_pop_male[551,index$ages+1],
     type = "l", xlab = "Age", ylab = "Population", main = "2005 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2005"]/5,
       col = "red")

## INFECTION DYNAMICS

# Total number in each infection compartment per timestep
plot(out$time,model_infectioncat_total$carriers,type = "l", ylim = c(0,2000000))
lines(out$time,model_infectioncat_total$acute,col= "green")
lines(out$time,model_infectioncat_total$sus,col= "red")
lines(out$time,model_infectioncat_total$immune,col= "blue")

plot(out$time,model_infectioncat_total$acute/model_pop_total$pop_total,col= "green")
plot(out$time,model_infectioncat_total$carriers/model_pop_total$pop_total)

plot(x = out$time, y = as.numeric(unlist(model_acute[,2]/model_pop[,2])))
plot(x = out$time, y = as.numeric(unlist(model_acute[,3]/model_pop[,3])))

plot(x = model_carriers[,1], y = as.numeric(unlist(model_carriers[,2]/model_pop[,2])))
plot(x = model_carriers[,1], y = as.numeric(unlist(model_carriers[,3]/model_pop[,3])))

# Running model with vaccine and no MTCT: acute and carrier prevalence go to 0
# With vaccine AND MTCT: no visual difference to without MTCT
# No vaccine: carrier and acute prevalence stabilise at about 0.15/0.008

# Carrier prevalence by age in 1980
plot(ages, model_carriers[which(model_carriers$time == 1980),-1]/model_pop[which(model_pop$time == 1980),-1], type = "l", ylim = c(0,0.3))
# Add Gambia data
points(gambia_prevdata$age, gambia_prevdata$edmunds_prev, col = "red")

# Proportion ever infected at 1 time point (after equilibrium is reached)
plot(1:n_agecat, ever_infected[1000,]/agespec_pop[1000,], type = "l", ylim = c(0,1))
points(gambia_prevdata$age, gambia_prevdata$edmunds_prop_ever_infected)

### Run model checks
# devtools::test()

