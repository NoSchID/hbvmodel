###################################
### Imperial HBV model 15/01/19 ###
###################################
# Model described in Shevanthi's thesis with some adaptations
# Currently only infant vaccination at 1 year of age, no birth dose or treatment


### Load packages ----
require(tidyr)
require(dplyr)
require(deSolve)
library(tictoc)
require(here)
library(profvis)

### Simulation parameters ----
## Country
countryname <- "gambia"

## Times
dt <- 0.1                             # timesteps
starttime <- 1850
runtime <- 250-dt                     # number of years to run the model for
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
index <- list("infcat_all" = 1:n_infectioncat,            # index for all infection status compartments
              "ages_all" = 1:n_agecat,                    # index for all age groups
              "ages_wocba" = ages_wocba+1,                # index for age group 15-49 years (women of childbearing age)
              "ages_1to5" = which(ages == 1):which(ages == 5),       # index for age groups 1-5 years
              "ages_6to15" = which(ages == 6):which(ages == 15),     # index for age groups 6-15 years
              "ages_16to100" = which(ages == 16):n_agecat)          # index for age groups 16-100 years

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
eag_loss <- 19.8873 * exp(-0.977 * ages)

## Calculate age-specific cancer rate progression function (Shevanthi)
cancer_prog_female <- 1e-07 * ((ages * 100* 0.2) + 2 * exp(0.0953 * ages))
cancer_prog_male <- 5 * cancer_prog_female

### Load demographic datasets ----
# Manual data preparation in Excel for 1950-2015 data:
# copy relevant data in a new file
# all datasets need the first column being labelled "time"
# set cells with estimates to "number" format and display full decimals
# save as CSV
# Manual data preparation in Excel for 2015-2100 data:
# deleted top rows and logo up until the header
# set cells with estimates to "number" format and display:
# 4 decimals for rates, all decimals for popsize
input_birthrate_data <- read.csv(here("data", countryname, "birthrate.csv"),
                                 stringsAsFactors = FALSE)  # already divided this by 1000

input_migration_data <- read.csv(here("data", countryname, "migration_rates.csv"),
                                 stringsAsFactors = FALSE)

# Age- and sex-specific datasets
input_popsize_female_1950to2015 <- read.csv(here("data",
                                                 countryname,
                                                 "WPP2017_POP_F15_3_ANNUAL_POPULATION_BY_AGE_FEMALE_1950_2015.csv"),
                                            header = TRUE, check.names = FALSE,
                                            stringsAsFactors = FALSE)
input_popsize_female_2015to2100 <- read.csv(here("data",
                                                 countryname,
                                                 "WPP2017_POP_F15_3_ANNUAL_POPULATION_BY_AGE_FEMALE_MEDIUM_2015_2100.csv"),
                                            header = TRUE, check.names = FALSE,
                                            stringsAsFactors = FALSE)

input_popsize_male_1950to2015 <- read.csv(here("data",
                                                 countryname,
                                                 "WPP2017_POP_F15_2_ANNUAL_POPULATION_BY_AGE_MALE_1950_2015.csv"),
                                            header = TRUE, check.names = FALSE,
                                            stringsAsFactors = FALSE)
input_popsize_male_2015to2100 <- read.csv(here("data",
                                                 countryname,
                                                 "WPP2017_POP_F15_2_ANNUAL_POPULATION_BY_AGE_MALE_MEDIUM_2015_2100.csv"),
                                            header = TRUE, check.names = FALSE,
                                            stringsAsFactors = FALSE)


# Abridged life tables for mortality rates and survival ratio (to calculate migration rates)
input_lifetables_male_1950to2015 <- read.csv(here("data",
                                                  countryname,
                                                  "WPP2017_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE_1950_2015.csv"),
                                             header = TRUE, check.names = FALSE,
                                             stringsAsFactors = FALSE)
input_lifetables_male_2015to2050 <- read.csv(here("data",
                                                  countryname,
                                                  "WPP2017_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE_MEDIUM_2015_2050.csv"),
                                             header = TRUE, check.names = FALSE,
                                             stringsAsFactors = FALSE)
input_lifetables_male_2050to2100 <- read.csv(here("data",
                                                  countryname,
                                                  "WPP2017_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE_MEDIUM_2050_2100.csv"),
                                             header = TRUE, check.names = FALSE,
                                             stringsAsFactors = FALSE)

input_lifetables_female_1950to2015 <- read.csv(here("data",
                                                    countryname,
                                                    "WPP2017_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE_1950_2015.csv"),
                                               header = TRUE, check.names = FALSE,
                                               stringsAsFactors = FALSE)
input_lifetables_female_2015to2050 <- read.csv(here("data",
                                                    countryname,
                                                    "WPP2017_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE_MEDIUM_2015_2050.csv"),
                                               header = TRUE, check.names = FALSE,
                                               stringsAsFactors = FALSE)
input_lifetables_female_2050to2100 <- read.csv(here("data",
                                                    countryname,
                                                    "WPP2017_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE_MEDIUM_2050_2100.csv"),
                                               header = TRUE, check.names = FALSE,
                                               stringsAsFactors = FALSE)

# 5-year survival rates (survival ratio) from abridged life tables for migration rates

# Age-specific fertility rates 1950-2015
input_fertility_1950to2015 <- read.csv(here("data",
                                            countryname,
                                            "WPP2017_FERT_F07_AGE_SPECIFIC_FERTILITY_1950_2015.csv"),
                                 header = TRUE, check.names = FALSE,
                                 stringsAsFactors = FALSE)
input_fertility_2015to2100 <- read.csv(here("data",
                                            countryname,
                                            "WPP2017_FERT_F07_AGE_SPECIFIC_FERTILITY_MEDIUM_2015_2100.csv"),
                                       header = TRUE, check.names = FALSE,
                                       stringsAsFactors = FALSE)

# For validation
# Deaths by age
input_deaths_female_1950to2015 <- read.csv(here("data",
                                                    countryname,
                                                    "WPP2017_MORT_F04_3_DEATHS_BY_AGE_FEMALE_1950_2015.csv"),
                                               header = TRUE, check.names = FALSE,
                                               stringsAsFactors = FALSE)
input_deaths_female_2015to2100 <- read.csv(here("data",
                                                countryname,
                                                "WPP2017_MORT_F04_3_DEATHS_BY_AGE_FEMALE_MEDIUM_2015_2100.csv"),
                                           header = TRUE, check.names = FALSE,
                                           stringsAsFactors = FALSE)

input_deaths_male_1950to2015 <- read.csv(here("data",
                                                countryname,
                                                "WPP2017_MORT_F04_2_DEATHS_BY_AGE_MALE_1950_2015.csv"),
                                           header = TRUE, check.names = FALSE,
                                           stringsAsFactors = FALSE)
input_deaths_male_2015to2100 <- read.csv(here("data",
                                                countryname,
                                                "WPP2017_MORT_F04_2_DEATHS_BY_AGE_MALE_MEDIUM_2015_2100.csv"),
                                           header = TRUE, check.names = FALSE,
                                           stringsAsFactors = FALSE)

# Total births
input_births_1950to2015 <- read.csv(here("data",
                                         countryname,
                                         "WPP2017_FERT_F01_BIRTHS_BOTH_SEXES_1950_2015.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)
input_births_2015to2100 <- read.csv(here("data",
                                              countryname,
                                              "WPP2017_FERT_F01_BIRTHS_BOTH_SEXES_MEDIUM_2015_2100.csv"),
                                         header = TRUE, check.names = FALSE,
                                         stringsAsFactors = FALSE)

### Functions for demographic data cleaning and preparation ----

## Functions for CLEANING loaded datasets from UN WPP

# Function to clean fertility rates datasets
clean_fertility_data <- function(fertility_dataset, country) {
  colnames(fertility_dataset) <- c("index", "variant", "location",
                                   "notes", "country_code", "time",
                                   fertility_dataset[1,7:13])  # set correct headers
  fertility_dataset_clean <- fertility_dataset[-1,] %>%  # delete header row from data
    filter(location == country) %>%  # subset data from desired country. This can be made more generic by looking for country code
    select(-index, -variant, -location, -notes, -country_code)

  fertility_dataset_clean[,-1] <- apply(fertility_dataset_clean[,-1], 2, function(x) as.numeric(x))

  # Now the dataset is the same as input_fertility_data
  return(fertility_dataset_clean)
}
# NOTE THAT FERTILITY RATES ARE STILL IN THOUSANDS AT THIS POINT

# Function to clean abridged life tables
clean_lifetables <- function(lifetable_dataset, country) {
  colnames(lifetable_dataset) <- c("index", "variant", "location",
                                   "notes", "country_code", "time",
                                   "age", "age_interval", "mortality_rate",
                                   lifetable_dataset[1,10:14],
                                   "survival_ratio", lifetable_dataset[1,16:ncol(lifetable_dataset)])  # set correct headers
  lifetable_dataset_clean <- lifetable_dataset[-1,] %>% # delete header row from data
    filter(location == country) %>% # subset desired country
    select(time, age, age_interval, mortality_rate, survival_ratio)  # remove unused data columns

  lifetable_dataset_clean[,-1] <- apply(lifetable_dataset_clean[,-1], 2, function(x) as.numeric(x))

  return(lifetable_dataset_clean)
}

# Function to clean population size/number datasets (e.g. population size, deaths, births -
# specifiy type (column name) in type_label)
clean_number_dataset <- function(dataset, country, type_label) {
  ## Clean input datasets (delete age groups with missing data and display integers)
  ## type_label = type of number in the dataset ("pop", "deaths", "births")

  if(type_label == "pop" | type_label == "deaths") {

  colnames(dataset) <- c("index", "variant", "location",
                                   "notes", "country_code", "time",
                                   dataset[1,7:ncol(dataset)])  # set correct headers
  dataset_clean <- dataset[-1,] %>% # delete header row from data
    filter(location == country) %>% # subset desired country
    select(-index, -variant, -location, - notes, -country_code) %>% # remove unused data columns
    gather(key = "age", value = "thousands", -time) %>%     # turn into long format
    arrange(time)

  dataset_clean$thousands <- as.numeric(dataset_clean$thousands)    # change numbers into numeric format

  # When switching to numeric format, cells with "..." turn to NA
  dataset_clean <- dataset_clean %>%
    drop_na %>%                                               # delete age groups with missing values
    mutate(number = as.numeric(thousands) * 1000) %>%         # numbers were imported as 1000s
    select(-thousands)
  names(dataset_clean) <- c("time", "age", type_label)

  return(dataset_clean)

  } else if(type_label == "births") {

    colnames(dataset) <- c("index", "variant", "location",
                           "notes", "country_code",
                           dataset[1,6:ncol(dataset)])  # set correct headers

    dataset_clean <- dataset[-1,] %>% # delete header row from data
      filter(location == country) %>% # subset desired country
      select(-index, -variant, -location, - notes, -country_code) %>% # remove unused data columns
      gather(key = "time", value = "thousands") %>%                   # turn into long format
      arrange(time)

    dataset_clean <- dataset_clean %>%
      mutate(number = as.numeric(thousands) * 1000) %>%         # numbers were imported as 1000s
      select(-thousands)

    names(dataset_clean) <- c("time", type_label)

    return(dataset_clean)

  } else {

    print("Not a valid data type (pop, deaths or births)")

  }

}

## Functions for PREPARING clean datasets for model input

# Function to prepare age-specific (columns) mortality rates by broad time period (rows)
prepare_mort_rates <- function(mortality_dataset) {
  mortality_dataset <- select(mortality_dataset, time, age, age_interval, mortality_rate)
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

  # Turn into matrix with time mid-points for more efficient interpolation over time
  out <- timevary_demog_rates(out)
#  out <- prepare_demogrates_for_timevary(out)

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

  # Turn into matrix with time mid-points for more efficient interpolation over time
  out <- timevary_demog_rates(out)

  return(out)
}

# Function to calculate and prepare age- and sex-specific migration rates
# Calculates age-specific number of net migrants using the cohort component forward method
# based on survival ratio from abridged life tables
calculate_migration_rates <- function(survival_dataset, pop_dataset) {
  survival_dataset <- select(survival_dataset, time, age, age_interval, survival_ratio) # remove mortality rates from life table

  # Prepare a dataset with population size for every 5 years
  popsize_5years <- pop_dataset %>%
    filter(time %in% seq(1950,2100,5))
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
  popsize_5years$time <- as.numeric(popsize_5years$time)    # turn into numeric format for merging

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
                   "2000-2005", "2005-2010", "2010-2015", "2015-2020", "2020-2025",
                   "2025-2030", "2030-2035", "2035-2040", "2040-2045", "2045-2050",
                   "2050-2055", "2055-2060", "2060-2065", "2065-2070", "2070-2075",
                   "2075-2080", "2080-2085", "2085-2090", "2090-2095", "2095-2100")
  migrants$time <- rep(time_labels, each = length(age_labels))


  # Calculate migration rates from risks
  migrants <- mutate(migrants, migration_rate = (-log(1-migrants/pop))/5) %>%
    select(time, age, migration_rate)
  # Migration rates were calculated from data for each 5-year age group and time period

  return(migrants)
}

calculate_migration_rates_under5 <- function(survival_dataset_female, pop_dataset_female) {
  survival_dataset_female <- select(survival_dataset_female, time, age, age_interval, survival_ratio) # remove mortality rates from life table

  # Prepare a dataset with population size for every 5 years
  popsize_5years <- pop_dataset_female %>%
    filter(time %in% seq(1950,2100,5))
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
  popsize_5years$time <- as.numeric(popsize_5years$time)    # turn into numeric format for merging

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
    summarise(cwr = pop[age == "0-4"]/sum(pop[age %in% c("15-19", "20-24", "25-29",
                                                         "30-34", "35-39", "40-44", "45-49")]))

  popsize_under5 <- filter(popsize_5years, age == "0-4") %>%
    select(pop)

  # Calculate number of migrants for 0-4 year olds using the child woman ratio
  migrants_under5 <- filter(migrants, age %in% c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44")) %>%
    group_by(time) %>%
    summarise(migrants_wocba = sum(migrants)) %>%
    drop_na() %>%
    mutate(time = replace(time, values = seq(1950,2095,5)))  %>%
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
  # Output is in wide format (rows = time period, columns = age groups)

  # Turn into matrix with time mid-points for more efficient interpolation over time
  out <- timevary_demog_rates(out)

  return(out)

}

# Prepare a dataset with population size for every age group
prepare_popsize <- function(pop_dataset) {
  ## Clean input datasets (delete age groups with missing data and display integers)
  pop_dataset$time <- as.numeric(pop_dataset$time)    # change time into numeric format

  ## Prepare numbers for each age step
  ## Need to distinguish between pre-1990 data (final age group 80+)
  ## and from 1990 onwards (final age group 100+)

  # Delete age group 100+
  popsize <- pop_dataset[pop_dataset$age != "100+",]

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

## Sub-functions called in above functions

# Interpolate demographic rates for every timestep (more efficient than calling within model)
timevary_demog_rates <- function(dataset) {
  # Input data frames are age-specific mortality rates, birth rate and migration rate for every 5-year period
  dataset <- data.matrix(dataset)   # convert data frame to matrix for faster manipulation

  # Create a complete dataset for model run time:
  # to run model for a whole aging cycle (100 years) for stabilisation,
  # need to duplicate 1950-2050 rates for the 1850-1950 period
  rates1950to2015 <- dataset[1:20,]
  dataset_complete <- rbind(rates1950to2015, dataset)
  dataset_complete[,1] <- seq(1852-1850, 2097-1850, 5)   # convert time to number starting from 0 and midpoint of time period
  lastrow <- dataset_complete[nrow(dataset_complete),]
  lastrow[1] <- times[length(times)]                     # repeat the last time point for the final timestep in model
  dataset_complete <- rbind(dataset_complete, lastrow)

  # Interpolate rates for every timestep in model:
  # Linear interpolation over time, but extrapolation is constant to nearest data point (for first/last timesteps)
  rates_interp <- apply(dataset_complete, 2, FUN = approx, x = dataset_complete[,1],
                        xout = times, method = "linear", rule = 2)
  rates_interp <- unlist(lapply(rates_interp, "[", "y"))
  rates_interp <- matrix(rates_interp, nrow = length(times)) # store in matrix (columns = age, rows = timesteps)
  rates_interp[,1] <- times  # relabel timesteps to start at 0
  rownames(rates_interp) <- rates_interp[,1]+starttime  # label timesteps according to year

  return(rates_interp)
}

### Prepare demographic data ----

## Clean population size datasets (also for output comparison)
# Females
input_popsize_female_1950to2015_clean <- clean_number_dataset(input_popsize_female_1950to2015,
                                                              "Gambia", "pop")
input_popsize_female_2015to2100_clean <- clean_number_dataset(input_popsize_female_2015to2100,
                                                              "Gambia", "pop") %>%
  filter(time != 2015)  # remove duplicated data for 2015 (present in both datasets)
input_popsize_female_clean <- rbind(input_popsize_female_1950to2015_clean,
                                    input_popsize_female_2015to2100_clean)
# Males
input_popsize_male_1950to2015_clean <- clean_number_dataset(input_popsize_male_1950to2015,
                                                            "Gambia", "pop")
input_popsize_male_2015to2100_clean <- clean_number_dataset(input_popsize_male_2015to2100,
                                                            "Gambia", "pop") %>%
  filter(time != 2015)  # remove duplicated data for 2015 (present in both datasets)
input_popsize_male_clean <- rbind(input_popsize_male_1950to2015_clean,
                                  input_popsize_male_2015to2100_clean)

## Prepare annual age- and sex-specific population size data
# Input: annual population size by 5-year age group and sex
popsize_female <- prepare_popsize(input_popsize_female_clean)
popsize_male <- prepare_popsize(input_popsize_male_clean)

# Extract 1950 population data
popsize_1950 <- left_join(popsize_female[popsize_female$time == "1950",],
                          popsize_male[popsize_female$time == "1950",],
                          by = c("time", "age")) %>%
  select(-time) %>%
  rename(pop_female = pop.x, pop_male = pop.y)

## Clean and prepare age-specific fertility rates by broad time period
input_fertility_clean <- rbind(clean_fertility_data(fertility_dataset =
                                                      input_fertility_1950to2015,
                                                    country = "Gambia"),
                               clean_fertility_data(fertility_dataset =
                                                      input_fertility_2015to2100,
                                                    country = "Gambia"))

# Input: age-specific fertility rates for age groups 15-50
fert_rates <- prepare_fert_rates(input_fertility_clean)
colnames(fert_rates) <- c("time", 15:49)

## Sex ratio at birth (Gambia-specific, constant over time)
sex_ratio <- c(0.4926, 0.5074)

## Clean lifetables datasets for mortality rates and survival ratio (calculation of migration rates)
# Males
input_lifetables_male_clean <- rbind(clean_lifetables(input_lifetables_male_1950to2015, "Gambia"),
                                     clean_lifetables(input_lifetables_male_2015to2050, "Gambia"),
                                     clean_lifetables(input_lifetables_male_2050to2100, "Gambia"))
# Females
input_lifetables_female_clean <- rbind(clean_lifetables(input_lifetables_female_1950to2015, "Gambia"),
                                       clean_lifetables(input_lifetables_female_2015to2050, "Gambia"),
                                       clean_lifetables(input_lifetables_female_2050to2100, "Gambia"))

## Prepare age- and sex-specific mortality rates by broad time period
# Input: central death rate in abridged life tables
mort_rates_female <- prepare_mort_rates(input_lifetables_female_clean)
colnames(mort_rates_female) <- c("time", ages)
mort_rates_male <- prepare_mort_rates(input_lifetables_male_clean)
colnames(mort_rates_male) <- c("time", ages)

## Prepare age- and sex-specific migration rates by broad time period
# Input: age-specific survival ratio from abridged life tables and annual population size
# by 5-year age group and sex
migration_rates_female <- prepare_migration_rates(input_lifetables_female_clean, input_popsize_female_clean,
                                                  input_lifetables_female_clean, input_popsize_female_clean)
colnames(migration_rates_female) <- c("time", ages)
migration_rates_male <- prepare_migration_rates(input_lifetables_male_clean, input_popsize_male_clean,
                                                input_lifetables_female_clean, input_popsize_female_clean)
colnames(migration_rates_male) <- c("time", ages)

## Datasets for model validations

# Total population size (both sexes) per year
popsize_total <- popsize_female %>%
  mutate(pop = pop + popsize_male$pop) %>%
  group_by(time) %>%
  summarise(pop = sum(pop))
# Note the numbers are not exactly the same as total UN WPP in later years because
# I deleted the 100+ age group

## Prepare total number of deaths for each year 1950-2015 (for model check)
# Input: age- and sex-specific number of deaths by 5-year time period (female and male dataset)
# Females
input_deaths_female_clean <- rbind(clean_number_dataset(input_deaths_female_1950to2015,
                                                        country = "Gambia",
                                                        type_label = "deaths"),
                                   clean_number_dataset(input_deaths_female_2015to2100,
                                                        country = "Gambia",
                                                        type_label = "deaths"))
# Males
input_deaths_male_clean <- rbind(clean_number_dataset(input_deaths_male_1950to2015,
                                                      country = "Gambia",
                                                      type_label = "deaths"),
                                 clean_number_dataset(input_deaths_male_2015to2100,
                                                      country = "Gambia",
                                                      type_label = "deaths"))


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

# Total number of births (both sexes) per 5 years
input_births_clean <- rbind(clean_number_dataset(input_births_1950to2015, "Gambia", "births"),
                             clean_number_dataset(input_births_2015to2100, "Gambia", "births"))

### Output-related functions ----
# Function to sum numbers from different compartments for each age and time step
sum_pop_by_age <- function(time = out$time, output) {
  out <- data.frame(time = time, output) %>%
    gather(key = "agegroup", value = "pop", -time) %>%       # turn into wide format
    arrange(time) %>%                                        # order by timestep
    mutate(agegroup = as.numeric(replace(agegroup,           # remove infection information
                                         values = ages))) %>%
    group_by(time, agegroup) %>%
    summarise(pop = sum(pop)) %>%                            # sum numbers for each age group at each timestep
    spread(key = "agegroup", value = "pop")                  # return to wide format

  # Note: this code is quicker than using do.call(cbind, by(t(output), rep(0:99,9), FUN = colSums))
  # slowest step is ordering by time

  return(as.data.frame(out))
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
    pop <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))
    initialpop <- array(unlist(init_pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    current_pop <- apply(pop,c(1,3),sum)
    pop_increase <- apply(initialpop, c(1,3), sum)/current_pop
    scaler <- array(c(rep(pop_increase[,1], n_infectioncat),
                      rep(pop_increase[,2], n_infectioncat)), dim = c(n_agecat,n_infectioncat,2))
    pop <- scaler * pop
    return(c(pop, 0, 0, rep(0,100), rep(0,100), rep(0,100), rep(0,100)))
  })
}

## THE MODEL
imperial_model <- function(timestep, pop, parameters){

  with(as.list(parameters), {

    # PREPARATION

    # Set up population array with infection compartments
    # Matrix 1 = females, matrix 2 = males, rows = agesteps, columns = infection compartments
    # Example: pop[ages,infectionstatus,sex]
    pop <- array(unlist(pop[1:(2 * n_infectioncat * n_agecat)]),dim=c(n_agecat,n_infectioncat,2))

    # Demography: define time-varying parameters (calls the row corresponding to current timestep)
    mortality_rate <- matrix(c(mort_rates_female[timestep*10+1,-1],
                               mort_rates_male[timestep*10+1,-1]),
                             ncol = 2)    # 2 columns for sex-specific rates
    fertility_rate <- fert_rates[timestep*10+1,-1]
    migration_rate <- matrix(c(migration_rates_female[timestep*10+1,-1],
                               migration_rates_male[timestep*10+1,-1]),
                             ncol = 2)    # 2 columns for sex-specific rates

    # Notation: label infection compartments
    S <- 1                           # Susceptible
    IT <- 2                           # Chronic infection: immune tolerant
    IR <- 3                           # Chronic infection: immune reactive
    IC <- 4                           # Chronic infection: inactive carrier
    ENCHB <- 5                        # Chronic infection: HBeAg-negative CHB
    CC <- 6                           # Chronic disease: compensated cirrhosis
    DCC <- 7                          # Chronic disease: decompensated cirrhosis
    HCC <- 8                          # Chronic disease: hepatocellular carcinoma
    R <- 9                           # Immune
    HBeAg_neg <- IC:HCC               # HBeAg-negative infected compartments
    HBeAg_pos <- c(IT,IR)             # HBeAg-positive infected compartments

    # Horizontal transmission: define WAIFW matrix
    # Assuming no effective contact between children and adults
    #beta <- matrix(0, nrow = n_agecat, ncol = n_agecat) # matrix of transmission rates
    #beta[index$ages_1to5, index$ages_1to5] <- b1        # transmission among children (1-5 years)
    #beta[index$ages_6to15, index$ages_6to15] <- b2      # transmission among juveniles (6-15 years)
    #beta[index$ages_16to100, index$ages_16to100] <- b3  # transmission among adults (16-100 years)
    #beta[index$ages_1to5, index$ages_6to15] <- b2       # transmission from children to juveniles
    #beta[index$ages_6to15, index$ages_1to5] <- b2       # transmission from juveniles to children
    #beta[index$ages_6to15, index$ages_16to100] <- b3    # transmission from juveniles to adults
    #beta[index$ages_16to100, index$ages_6to15] <- b3    # transmission from adults to juveniles

    # Infant vaccination: set date for introduction
    # Vaccination coverage is 0 until the specified starttime and only if vaccine switch is on
    if (apply_vacc == 1 & timestep >= (vacc_introtime-starttime)) {
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
    infections <- matrix(rep(0, 2* n_agecat),
                         ncol = 2, nrow = n_agecat)      # female and male incident infections


    # MODEL EQUATIONS

    # Mother-to-child transmission and births
    infected_births <- sum(p_chronic[1] * fertility_rate * (apply(mtct_prob_e * pop[index$ages_wocba,c(IT,IR),1],1,sum) + apply(mtct_prob_s * pop[index$ages_wocba,IC:HCC,1],1,sum)))    # infected births come from acute and chronic women of childbearing age
    uninfected_births <- sum(fertility_rate * pop[index$ages_wocba,index$infcat_all,1]) - infected_births  # uninfected births = all births - infected babies
    births <- infected_births + uninfected_births
    # applying the same age-specific fertility rate to every infection compartment

    # Age-specific force of infection (same for men and women)

 #   foi <- (beta %*%
#              (apply(pop[index$ages_all,HBeAg_neg,1:2],1,sum) +
#               apply(alpha * pop[index$ages_all,HBeAg_pos,1:2],1,sum)))/
#               sum(pop[index$ages_all,index$infcat_all,1:2])

    # Multiply WAIFW matrix by proportion in infectious compartments (IT, IR, IC, ENCHB, CC, DCC, HCC)
    # (= proportion of total population size (1 number) at each time step)
    # HBeAg-positive individuals (IT, IR) are more infectious than HBeAg-negatives (multiply by alpha)
    # Returns a vector with force of infection for every age - 4 different values:
    # 0 in 0-year olds, different values for 1-5, 6-15 and 16-100 year olds

    # Imperial model FOI

    foi <- rep(0,100)
    foi[index$ages_1to5] <- b1 * sum(apply(pop[index$ages_1to5,HBeAg_neg,1:2],1,sum))/sum(pop[index$ages_1to5,index$infcat_all,1:2]) +
                            min(1,b1 * alpha) * sum(apply(pop[index$ages_1to5,HBeAg_pos,1:2],1,sum))/sum(pop[index$ages_1to5,index$infcat_all,1:2])

    foi[2:16] <- foi[2:16] + b2 * sum(apply(pop[2:16,HBeAg_neg,1:2],1,sum))/sum(pop[2:16,index$infcat_all,1:2]) +
                            (b2 * alpha) * sum(apply(pop[2:16,HBeAg_pos,1:2],1,sum))/sum(pop[2:16,index$infcat_all,1:2])

    foi[6:100] <- foi[6:100] + b3 * sum(apply(pop[6:100,HBeAg_neg,1:2],1,sum))/sum(pop[6:100,index$infcat_all,1:2]) +
                              (b3 * alpha) * sum(apply(pop[6:100,HBeAg_pos,1:2],1,sum))/sum(pop[6:100,index$infcat_all,1:2])

    # Partial differential equations (solving for each sex separately)

      for (i in 1:2) {        # i = sex [1 = female, 2 = male]

        # Demography: Incident deaths and migrants
        deaths[index$ages_all,index$infcat_all,i] <- mortality_rate[index$ages_all,i] *
          pop[index$ages_all,index$infcat_all,i]
        migrants[index$ages_all,index$infcat_all,i] <-  migration_rate[index$ages_all,i] * pop[index$ages_all,index$infcat_all,i]
        # Applying same age-specific mortality rate to every infection compartment
        # Returns an array with indicent deaths and net migrants for every age (rows), infection state (columns) and sex (arrays)

        # Infection: Incident infections
        infections[index$ages_all,i] <- foi * pop[index$ages_all,S,i]
        # Returns a matrix with incident infections for every age (rows) and every sex (columns)

        # Transitions between compartments:

        # Susceptibles
        dpop[index$ages_all,S,i] <- -(diff(c(0,pop[-length(index$ages_all),S,i],0))/da) -
          p_chronic * infections[index$ages_all,i] -
          (1-p_chronic) * infections[index$ages_all,i] -
          deaths[index$ages_all,S,i] + migrants[index$ages_all,S,i]
        #  - (vacc_cov * vacc_eff * pop[index$ages_all,S,i])

        # Immune tolerant
        dpop[index$ages_all,IT,i] <- -(diff(c(0,pop[-length(index$ages_all),IT,i],0))/da) +
          p_chronic * infections[index$ages_all,i] -
          pr_it_ir * pop[index$ages_all,IT,i] -
          hccr_it[index$ages_all,i] * pop[index$ages_all,IT,i] -
          deaths[index$ages_all,IT,i] + migrants[index$ages_all,IT,i]

        # Immune reactive
        dpop[index$ages_all,IR,i] <- -(diff(c(0,pop[-length(index$ages_all),IR,i],0))/da) +
          pr_it_ir * pop[index$ages_all,IT,i] -
          pr_ir_ic * pop[index$ages_all,IR,i] -
          pr_ir_enchb * pop[index$ages_all,IR,i] -
          hccr_ir[index$ages_all,i] * pop[index$ages_all,IR,i] -
          deaths[index$ages_all,IR,i] + migrants[index$ages_all,IR,i]

        # Inactive carrier
        dpop[index$ages_all,IC,i] <- -(diff(c(0,pop[-length(index$ages_all),IC,i],0))/da) +
          pr_ir_ic * pop[index$ages_all,IR,i] -
          pr_ic_enchb * pop[index$ages_all,IC,i] -
          sag_loss * pop[index$ages_all,IC,i] -
          hccr_ic[index$ages_all,i] * pop[index$ages_all,IC,i] -
          deaths[index$ages_all,IC,i] + migrants[index$ages_all,IC,i]

        # HBeAg-negative CHB
        dpop[index$ages_all,ENCHB,i] <- -(diff(c(0,pop[-length(index$ages_all),ENCHB,i],0))/da) +
          pr_ir_enchb * pop[index$ages_all,IR,i] +
          pr_ic_enchb * pop[index$ages_all,IC,i] -
          ccrate * pop[index$ages_all,ENCHB,i] -
          hccr_enchb[index$ages_all,i] * pop[index$ages_all,ENCHB,i] -
          deaths[index$ages_all,ENCHB,i] + migrants[index$ages_all,ENCHB,i]

        # Compensated cirrhosis
        dpop[index$ages_all,CC,i] <- -(diff(c(0,pop[-length(index$ages_all),CC,i],0))/da) +
          ccrate * pop[index$ages_all,ENCHB,i] -
          dccrate * pop[index$ages_all,CC,i] -
          hccr_cc[index$ages_all,i] * pop[index$ages_all,CC,i] -
          mu_cc * pop[index$ages_all,CC,i] -
          deaths[index$ages_all,CC,i] + migrants[index$ages_all,CC,i]

        # Decompensated cirrhosis
        dpop[index$ages_all,DCC,i] <- -(diff(c(0,pop[-length(index$ages_all),DCC,i],0))/da) +
          dccrate * pop[index$ages_all,CC,i] -
          hccr_dcc * pop[index$ages_all,DCC,i] -
          mu_dcc * pop[index$ages_all,DCC,i] -
          deaths[index$ages_all,DCC,i] + migrants[index$ages_all,DCC,i]

        # HCC
        dpop[index$ages_all,HCC,i] <- -(diff(c(0,pop[-length(index$ages_all),HCC,i],0))/da) +
          hccr_it[index$ages_all,i] * pop[index$ages_all,IT,i] +
          hccr_ir[index$ages_all,i] * pop[index$ages_all,IR,i] +
          hccr_ic[index$ages_all,i] * pop[index$ages_all,IC,i] +
          hccr_enchb[index$ages_all,i] * pop[index$ages_all,ENCHB,i] +
          hccr_cc[index$ages_all,i] * pop[index$ages_all,CC,i] +
          hccr_dcc * pop[index$ages_all,DCC,i] -
          mu_hcc * pop[index$ages_all,HCC,i] -
          deaths[index$ages_all,HCC,i] + migrants[index$ages_all,HCC,i]

        # Immunes
        dpop[index$ages_all,R,i] <- -(diff(c(0,pop[-length(index$ages_all),R,i],0))/da) +
          (1-p_chronic) * infections[index$ages_all,i] +
          sag_loss * pop[index$ages_all,IC,i] -
          deaths[index$ages_all,R,i] + migrants[index$ages_all,R,i]
        #+ vacc_cov * vacc_eff * pop[index$ages_all,S,i]

        # Babies are born susceptible or infected (age group 1)
        dpop[1,S,i] <- dpop[1,S,i] + sex_ratio[i] * uninfected_births
        dpop[1,IT,i] <- dpop[1,IT,i] + sex_ratio[i] * infected_births
        #dpop[1,R,i] <- dpop[1,R,i] + sex_ratio[i] * (1-p_chronic[1]) * infected_births

        # Vaccination: applied at 1 year of age
        dpop[2,S,i] <- dpop[2,S,i] - (vacc_cov * vacc_eff * pop[2,S,i])
        dpop[2,R,i] <- dpop[2,R,i] + (vacc_cov * vacc_eff * pop[2,S,i])

      }

    # OUTPUT

    # Sum age-specific number of incident deaths across infection compartments
    deaths_out <- apply(deaths,c(1,3),sum)

    # Return results
    res <- c(dpop, births, infected_births, deaths_out, infections)
    #res <- c(dpop, deaths, migrants, births)
    list(res)
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
# Old function, not used anymore
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

run_model <- function(..., default_parameter_list, parms_to_change = list(...)) {

  ## Define parameter values for model run:
  # Using default input parameter list or with updated values specified in parms_to_change
  parameters <- generate_parameters(default_parameter_list = default_parameter_list,
                                    parms_to_change = parms_to_change)

  ## Run model simulation
  out <- as.data.frame(ode.1D(y = init_pop, times = times, func = imperial_model,
                              parms = parameters, nspec = 1, method = "euler",
                              events = list(func = reset_pop_1950, time = 100)))
  # events = list(func = eventfun, time = 100)
  out$time   <-  out$time + starttime

  ## Store different types of outputs as a list
  pop <- data.frame(time = out$time, out[,2:(n_agecat*n_infectioncat*2+1)])
  births <- data.frame(time = out$time, select(out, contains("births")))
  deaths <- data.frame(time = out$time, select(out, contains("deaths")))
  infections <- data.frame(time = out$time, select(out, contains("incidence")))

  toreturn <- list(out = pop, births = births, deaths = deaths, incidence = infections)

  return(toreturn)
}


### Model input ----

# DEMOGRAPHY: SEE ABOVE

# INITIAL POPULATION
# Set up initial population: age- and sex-specific population size in 1950
# Note: names in initial population vector is reproduced in output
# guessed proportion in each compartment
init_pop <- c("Sf" = popsize_1950$pop_female*(1-gambia_infected),
               "ITf" = popsize_1950$pop_female*gambia_infected*gambia_eag*0.5,
               "IRf" = popsize_1950$pop_female*gambia_infected*gambia_eag*0.5,
               "ICf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "ENCHBf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "CCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "DCCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "HCCf" = popsize_1950$pop_female*gambia_infected*(1-gambia_eag)*0.2,
               "Rf" = rep(0,100),
               "Sm" = popsize_1950$pop_male*(1-gambia_infected),
               "ITm" = popsize_1950$pop_male*gambia_infected*gambia_eag*0.5,
               "IRm" = popsize_1950$pop_male*gambia_infected*gambia_eag*0.5,
               "ICm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "ENCHBm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "CCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "DCCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "HCCm" = popsize_1950$pop_male*gambia_infected*(1-gambia_eag)*0.2,
               "Rm" = rep(0,100),
               "cum_births" = 0, "cum_infected_births" = 0,
               "cum_deathsf" = rep(0,100), "cum_deathsm" = rep(0,100),
               "cum_incidencef" = rep(0,100), "cum_incidencem" = rep(0,100))
# Total population in 1950:
N0 <- sum(init_pop[1:(n_infectioncat * n_agecat * 2)])

## TRANSMISSION, NATURAL HISTORY AND INTERVENTION PARAMETERS
# p_chronic: Edmunds approach adapted by using 0.89 for whole first year instead of just 0.5 years

parameter_list <- list(
  # TRANSMISSION PARAMETERS
  b1 = 0.04,
  b2 = 0.001,
  b3 = 0.001,
  alpha = 15,         # Shevanthi value, relative infectiousness of eAg-positives
  mtct_prob_e = 0.9,  # Shevanthi value, probability of perinatal transmission from HBeAg-positive mother
  mtct_prob_s = 0.05, # Shevanthi model value, probability of perinatal transmission from HBeAg-negative infected mother
  # NATURAL HISTORY PROGRESSION RATES
  p_chronic = c(0.89, exp(-0.65*ages[-1]^0.46)),     # Age-dependent probability of chronic carriage. Adapted Edmunds
  pr_it_ir = 0.1 * eag_loss,
  pr_ir_ic = 0.05 * eag_loss,
  pr_ir_enchb = 0.005,
  pr_ic_enchb = 0.01,
  sag_loss = 0.01,  # Shevanthi value, inactive carrier to recovered transition
  #sag_loss = sagloss_rates_0to80,
  ccrate = 0.04,  # Progression to CC (from ENCHB)
  dccrate = 0.04,  # Progression to DCC (from CC)
  # PROGRESSION RATES TO HEPATOCELLULAR CARCINOMA
  hccr_it = matrix(data = c(cancer_prog_female, cancer_prog_male), nrow = n_agecat, ncol = 2),
  hccr_ir = matrix(data = c(2*cancer_prog_female, 2*cancer_prog_male), nrow = n_agecat, ncol = 2),
  hccr_ic = matrix(data = c(0.5*cancer_prog_female, 0.5*cancer_prog_male), nrow = n_agecat, ncol = 2),
  hccr_enchb = matrix(data = c(2*cancer_prog_female, 2*cancer_prog_male), nrow = n_agecat, ncol = 2),
  hccr_cc = matrix(data = c(13*cancer_prog_female, 13*cancer_prog_male), nrow = n_agecat, ncol = 2),
  hccr_dcc = 0.04,
  # HBV-RELATED MORTALITY RATES (MORTALITY FROM LIVER DISEASE)
  mu_cc = 0.039,
  mu_dcc = 0.314,
  mu_hcc = 0.5,
  # INFANT VACCINATION PARAMETERS
#  vacc_cov = c(0, 0.92, rep(0,n_agecat-2)),  # vaccine is only applied in 1-year olds, need to get time-varying data
  vacc_cov = 0.92,
  vacc_eff = 0.95,                           # vaccine efficacy
  vacc_introtime = 1991,                     # year of vaccine introduction
  # INTERVENTION ON/OFF SWITCH (1/0)
  apply_vacc = 1)

# Store names of all parameters
parameter_names <- names(parameter_list)

### Run the model ----

# Default intervention: infant vaccine (apply_vacc = 1)
# Switch off by setting to 0 in function call
tic()
output <- run_model(default_parameter_list = parameter_list,
                 parms_to_change = list(apply_vacc = 0))
toc()


#out <- return_compartment_output(b1 = b1, b2 = b2, b3 = b3,
#                 alpha = alpha, gamma_acute = gamma_acute, p_chronic = p_chronic,
#                 sag_loss = sag_loss, mu_hbv = mu_hbv,
#                 mtct_prob_a = mtct_prob_a, mtct_prob_i = mtct_prob_i, mu = mu, b = b)

### Code output ----

## Extract separate outputs

# Infection compartments
out <- output$out
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

# Population
out_popf <- select(out[,2:(n_agecat*n_infectioncat*2+1)], contains("f"))
out_popm <- select(out[,2:(n_agecat*n_infectioncat*2+1)], contains("m"))
out_pop <- cbind(out_popf, out_popm)

# Cumulative deaths
out_cum_deathsf <- select(output$deaths, starts_with("cum_deathsf"))
out_cum_deathsm <- select(output$deaths, starts_with("cum_deathsm"))

# Cumulative HBV incidence from horizontal transmission only
out_cum_incidencef <- select(output$incidence, starts_with("cum_incidencef"))
out_cum_incidencem <- select(output$incidence, starts_with("cum_incidencem"))


## Code infection outputs

# Age-specific number in each infection compartment per time step
model_sus <- data.frame(time = out$time, pop = out_sf + out_sm)   # need to change the column names
model_carriers <- data.frame(time = out$time,
                             pop = (out_itf + out_itm +
                             out_irf + out_irm +
                             out_icf+out_icm+
                             out_enchbf+out_enchbm+
                             out_ccf+out_ccm+
                             out_dccf+out_dccm+
                             out_hccf+out_hccm))
model_carriers_female <- data.frame(time = out$time,
                                    pop = (out_itf+
                                             out_irf+
                                             out_icf+
                                             out_enchbf+
                                             out_ccf+out_dccf+out_hccf))
model_carriers_male <- data.frame(time = out$time,
                                    pop = (out_itm+
                                             out_irm+
                                             out_icm+
                                             out_enchbm+
                                             out_ccm+out_dccm+out_hccm))
model_immune <- data.frame(time = out$time, pop = out_rf + out_rm)
model_ever_inf <- data.frame(time = out$time,
                             pop = model_carriers[,-1] + model_immune[,-1])

# Total number in each infection compartment per time step
model_infectioncat_total <- data.frame(time = out$time,
                                sus = apply(model_sus[,-1], 1, sum),
                                carriers = apply(model_carriers[,-1], 1, sum),
                                immune = apply(model_immune[,-1], 1, sum),
                                ever_infected = apply(model_ever_inf[,-1],1,sum))

# Age-specific HBV incidence from horizontal transmission only - for women, men and both (total)
model_horizontal_infections_female <- data.frame(time = out$time,
                                                 incidence = rbind(out_cum_incidencef[1,],
                                                                   apply(out_cum_incidencef,
                                                                         2, diff, lag = 1)))
names(model_horizontal_infections_female)[-1] <- sprintf("incidence%d",ages)

model_horizontal_infections_male <- data.frame(time = out$time,
                                               incidence = rbind(out_cum_incidencem[1,],
                                                                 apply(out_cum_incidencem,
                                                                       2, diff, lag = 1)))
names(model_horizontal_infections_male)[-1] <- sprintf("incidence%d",ages)


#ever_infected <- acute + carriers + immune

## Code demography outputs

# Population:

# Age-specific and total (last column) population per time step
model_pop_female <- sum_pop_by_age(output = out_popf)
model_pop_male <- sum_pop_by_age(output = out_popm)
model_pop <- cbind(time = model_pop_female$time, model_pop_female[,-1] + model_pop_male[,-1])

# Total female, male and both population per time step
model_pop_total <- data.frame(time = out$time,
                              pop_female = apply(model_pop_female[,-1], 1, sum),
                              pop_male = apply(model_pop_male[,-1], 1, sum)) %>%
  mutate(pop_total = pop_female + pop_male)

# Births:

# Total number of births at each timestep
model_births <- data.frame(time = out$time,
                           births = c(output$births$cum_births[1],
                                      diff(output$births$cum_births, lag = 1)))
names(model_births) <- c("time", "births")

# Total number of births grouped in 5-year time periods
model_births_group5 <- model_births %>%
  mutate(timegroup = floor(time / 5) * 5) %>%
  group_by(timegroup) %>%
  summarise_all(sum) %>%
  select(-time)

# Deaths:

# Age-specific and total deaths per time step - for women, men and both (total)
model_deaths_female <- data.frame(time = out$time,
                                  deaths = rbind(out_cum_deathsf[1,],
                                                 apply(out_cum_deathsf, 2, diff, lag = 1)))
names(model_deaths_female)[-1] <- sprintf("deaths%d",ages)
model_deaths_female$total <- apply(model_deaths_female[,-1], 1, sum)
model_deaths_male <- data.frame(time = out$time,
                                deaths = rbind(out_cum_deathsm[1,],
                                               apply(out_cum_deathsm, 2, diff, lag = 1)))
names(model_deaths_male)[-1] <- sprintf("deaths%d",ages)
model_deaths_male$total <- apply(model_deaths_male[,-1], 1, sum)
model_deaths_total <- data.frame(time = out$time,
                                 deaths = model_deaths_female[,-1] + model_deaths_male[,-1])
names(model_deaths_total)[-1] <- c(sprintf("deaths%d",ages), "total")

# Total number of deaths grouped in 5-year time periods
model_deaths_total_group5 <- model_deaths_total %>%
  mutate(timegroup = floor(time / 5) * 5) %>%
  group_by(timegroup) %>%
  summarise_all(sum) %>%
  select(-time)


# Check:
all.equal(model_sus[,-1] + model_carriers[,-1] + model_immune[,-1], model_pop[,-1], check.names = FALSE)
all.equal(model_infectioncat_total$sus + model_infectioncat_total$carriers + model_infectioncat_total$immune,
          model_pop_total$pop_total, check.names = FALSE)



### Plots ----

## DEMOGRAPHY

## Total population, births, deaths

# Plot total population size over timesteps
plot(model_pop_total$time, model_pop_total$pop_total,
     xlab = "Year", ylab = "Population size")
points(popsize_total$time, popsize_total$pop, col = "red")

# Plot total number of births over time periods
plot(x = as.numeric(strtrim(input_births_clean$time, width = 4)),
     y = input_births_clean$births, col = "red",
     xlab = "5-year time periods", ylab = "Total number of births")
lines(model_births_group5$timegroup, model_births_group5$births)

# Plot total number of deaths over time periods
plot(x = as.numeric(strtrim(deaths_total$time, width = 4)),
     y = deaths_total$deaths, col = "red",
     xlab = "5-year time periods", ylab = "Total number of deaths")
lines(model_deaths_total_group5$timegroup, model_deaths_total_group5$total)

## Age structure

# Plot female age structure in 1970
plot(x = ages,
     y = model_pop_female[which(model_pop_female$time == 1970),index$ages_all+1],
     type = "l", xlab = "Age", ylab = "Population", main = "1970 - women")
points(x = seq(2,82,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "1970"]/5,
       col = "red")

# Plot male age structure in 2005
plot(x = ages,
     y = model_pop_male[which(model_pop_male$time == 2005),index$ages_all+1],
     type = "l", xlab = "Age", ylab = "Population", main = "2005 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2005"]/5,
       col = "red")

# Plot male age structure in 2070
plot(x = ages,
     y = model_pop_male[which(model_pop_male$time == 2070),index$ages_all+1],
     type = "l", xlab = "Age", ylab = "Population", main = "2070 - men")
points(x = seq(2,102,5),
       y = input_popsize_male_clean$pop[input_popsize_male_clean$time == "2070"]/5,
       col = "red")

# Plot female age structure in 2050
plot(x = ages,
     y = model_pop_female[which(model_pop_female$time == 2050),index$ages_all+1],
     type = "l", xlab = "Age", ylab = "Population", main = "2050 - women")
points(x = seq(2,102,5),
       y = input_popsize_female_clean$pop[input_popsize_female_clean$time == "2050"]/5,
       col = "red")

## INFECTION DYNAMICS

# Total number in each infection compartment per timestep
plot(out$time,model_infectioncat_total$carriers,type = "l", ylim = c(0,4000000))
lines(out$time,model_infectioncat_total$sus,col= "red")
lines(out$time,model_infectioncat_total$immune,col= "blue")

# Proportion in each infection compartment per timestep
plot(out$time,model_infectioncat_total$carriers/model_pop_total$pop_total,type = "l", ylim = c(0,0.7))
lines(out$time,model_infectioncat_total$sus/model_pop_total$pop_total,col= "red")
lines(out$time,model_infectioncat_total$immune/model_pop_total$pop_total,col= "blue")

# Carrier prevalence over time
plot(out$time,model_infectioncat_total$carriers/model_pop_total$pop_total, ylim = c(0,0.5))
lines(out$time,apply(model_carriers_female[,-1],1,sum)/model_pop_total$pop_female, col = "pink")
lines(out$time,apply(model_carriers_male[,-1],1,sum)/model_pop_total$pop_male, col = "blue")

model_infectioncat_total$carriers[which(model_infectioncat_total$time == 1990)]/
  model_pop_total$pop_total[which(model_pop_total$time == 1990)]

model_infectioncat_total$carriers[which(model_infectioncat_total$time == 2015)]/
  model_pop_total$pop_total[which(model_pop_total$time == 2015)]

# Carrier prevalence over time in different age groups
plot(x = model_carriers[,1], y = as.numeric(unlist(model_carriers[,2]/model_pop[,2]))) # age 0
plot(x = model_carriers[,1], y = as.numeric(unlist(model_carriers[,3]/model_pop[,3]))) # age 1
plot(x = model_carriers[,1], y = as.numeric(unlist(model_carriers[,4]/model_pop[,4]))) # age 2
plot(x = model_carriers[,1], y = as.numeric(unlist(model_carriers[,22]/model_pop[,22]))) # age 20

# Carrier prevalence by age in 1980
plot(ages, model_carriers[which(model_carriers$time == 1980),-1]/
       model_pop[which(model_pop$time == 1980),-1], type = "l", ylim = c(0,0.3))
points(gambia_prevdata$age, gambia_prevdata$edmunds_prev, col = "red")

# Carrier prevalence by age in 2015
plot(ages, model_carriers[which(model_carriers$time == 2015),-1]/
       model_pop[which(model_pop$time == 2015),-1], type = "l", ylim = c(0,0.3))

# Proportion ever infected in 1980
plot(ages, model_ever_inf[which(model_ever_inf$time == 1980),-1]/
       model_pop[which(model_pop$time == 1980),-1], type = "l", ylim = c(0,1))
points(gambia_prevdata$age, gambia_prevdata$edmunds_prop_ever_infected, col = "red")


plot(x = out$time, apply(out_sf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.5))
plot(x = out$time, apply(out_sm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.5))
plot(x = out$time, apply(out_itf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.5))
plot(x = out$time, apply(out_itm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.5))
plot(x = out$time, apply(out_irf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.2))
plot(x = out$time, apply(out_irm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.2))
plot(x = out$time, apply(out_icf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.2))
plot(x = out$time, apply(out_icm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.2))
plot(x = out$time, apply(out_enchbf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.1))
plot(x = out$time, apply(out_enchbm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.1))
plot(x = out$time, apply(out_ccf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.05))
plot(x = out$time, apply(out_ccm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.05))
plot(x = out$time, apply(out_dccf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.05))
plot(x = out$time, apply(out_dccm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.05))
plot(x = out$time, apply(out_hccf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.05))
plot(x = out$time, apply(out_hccm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.05))
plot(x = out$time, apply(out_rf,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.5))
plot(x = out$time, apply(out_rm,1,sum)/apply(out[,-1],1,sum), ylim = c(0,0.5))

### Run model checks
# devtools::test()

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

