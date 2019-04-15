########################################
### Imperial HBV model               ###
### Clean input demographic data     ###
### Source: UN WPP                   ###
########################################

### Prepare datasets for all years between 1850 and 2100
# UN WPP 2017 provides data from 1950 to 2015, and projections from 2015 to 2100
# Copy 1950-2050 data for 1850-1950 period
demog_times <- round((0:(250/dt))*dt,2)
demog_times_labels <- demog_times + 1850

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

## Functions for PREPARING clean datasets for model input (e.g. interpolation over age and time)

# Function to prepare age-specific (columns) mortality rates by broad time period (rows)
prepare_mort_rates <- function(mortality_dataset) {
  mortality_dataset <- select(mortality_dataset, time, age, age_interval, mortality_rate)
  # Input dataset is a mortality rate for ages 0, 1-5, 5-10, ..., 80-85, 85-100

  # Constant interpolation for missing ages
  mort_rates_interp <- 0
  for(i in mortality_dataset$time) {
    mort_rates_interp[i] <- list(approx(x = mortality_dataset$age[mortality_dataset$time == i],
                                        y = mortality_dataset$mortality_rate[mortality_dataset$time == i],
                                        xout = c(ages),
                                        method = "constant",
                                        rule = 2))
  }

  out <- matrix(unlist(lapply(mort_rates_interp[-1], "[", 2)), ncol = n_agecat, byrow = TRUE)
  out <- cbind(time = as.character(unique(mortality_dataset$time)), as.data.frame(out))
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

  # Constant interpolation for missing ages
  fert_rates_interp <- 0
  for(i in fert_rates$time) {
    fert_rates_interp[i] <- list(approx(x = fert_rates$age[fert_rates$time == i],
                                        y = fert_rates$fertility_rate[fert_rates$time == i],
                                        xout = c(ages_wocba),
                                        method = "constant", rule = 2))
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
  migration_rates <- migrants

  # Constant interpolation for missing ages
  migration_rates_interp <- 0
  for(i in migration_rates$time) {
    migration_rates_interp[i] <- list(approx(x = migration_rates$age[migration_rates$time == i],
                                             y = migration_rates$migration_rate[migration_rates$time == i],
                                             xout = c(ages),
                                             method = "constant", rule = 2))
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
  lastrow[1] <- demog_times[length(demog_times)]                     # repeat the last time point for the final timestep in model
  dataset_complete <- rbind(dataset_complete, lastrow)

  # Interpolate rates for every timestep in model:
  # Linear interpolation over time, but extrapolation is constant to nearest data point (for first/last timesteps)
  rates_interp <- apply(dataset_complete, 2, FUN = approx, x = dataset_complete[,1],
                        xout = demog_times, method = "linear", rule = 2)
  rates_interp <- unlist(lapply(rates_interp, "[", "y"))
  rates_interp <- matrix(rates_interp, nrow = length(demog_times)) # store in matrix (columns = age, rows = timesteps)
  rates_interp[,1] <- demog_times  # relabel timesteps to start at 0
  rownames(rates_interp) <- rates_interp[,1]+1850  # label timesteps according to year

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
colnames(fert_rates) <- c("time", ages_wocba)

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


# Save the cleaned datasets (need to specify time/agestep in filename)
#save(input_popsize_female_clean = input_popsize_female_clean,
#     input_popsize_male_clean = input_popsize_male_clean,
#     popsize_female = popsize_female,
#     popsize_male = popsize_male,
#     popsize_1950 = popsize_1950,
#     input_fertility_clean = input_fertility_clean,
#     fert_rates = fert_rates,
#     sex_ratio = sex_ratio,
#     input_lifetables_male_clean = input_lifetables_male_clean,
#     input_lifetables_female_clean = input_lifetables_female_clean,
#     mort_rates_female = mort_rates_female,
#     mort_rates_male = mort_rates_male,
#     migration_rates_female = migration_rates_female,
#     migration_rates_male = migration_rates_male,
#     popsize_total = popsize_total,
#     input_deaths_female_clean = input_deaths_female_clean,
#     input_deaths_male_clean = input_deaths_male_clean,
#     deaths_total = deaths_total,
#     deaths_1950 = deaths_1950,
#     input_births_clean = input_births_clean,
#     file = here("data/demogdata_0point5.RData"))

# Combine all datasets into a list
#input_demographic_data_clean <- list(
#  input_popsize_female_clean = input_popsize_female_clean,
#  input_popsize_male_clean = input_popsize_male_clean,
#  popsize_female = popsize_female,
#  popsize_male = popsize_male,
#  popsize_1950 = popsize_1950,
#  input_fertility_clean = input_fertility_clean,
#  fert_rates = fert_rates,
#  sex_ratio = sex_ratio,
#  input_lifetables_male_clean = input_lifetables_male_clean,
#  input_lifetables_female_clean = input_lifetables_female_clean,
#  mort_rates_female = mort_rates_female,
#  mort_rates_male = mort_rates_male,
#  migration_rates_female = migration_rates_female,
#  migration_rates_male = migration_rates_male,
#  popsize_total = popsize_total,
#  input_deaths_female_clean = input_deaths_female_clean,
#  input_deaths_male_clean = input_deaths_male_clean,
#  deaths_total = deaths_total,
#  deaths_1950 = deaths_1950,
#  input_births_clean = input_births_clean
#)

