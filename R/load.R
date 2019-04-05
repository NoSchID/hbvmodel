########################################
### Imperial HBV model               ###
### Load required input              ###
########################################

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

inpath_demography <- "data-raw/demography"

#input_birthrate_data <- read.csv(here(inpath_demography, countryname, "birthrate.csv"),
#                                 stringsAsFactors = FALSE)  # already divided this by 1000

#input_migration_data <- read.csv(here(inpath_demography, countryname, "migration_rates.csv"),
#                                 stringsAsFactors = FALSE)

# Age- and sex-specific datasets
input_popsize_female_1950to2015 <- read.csv(here(inpath_demography,
                                                 countryname,
                                                 "WPP2017_POP_F15_3_ANNUAL_POPULATION_BY_AGE_FEMALE_1950_2015.csv"),
                                            header = TRUE, check.names = FALSE,
                                            stringsAsFactors = FALSE)
input_popsize_female_2015to2100 <- read.csv(here(inpath_demography,
                                                 countryname,
                                                 "WPP2017_POP_F15_3_ANNUAL_POPULATION_BY_AGE_FEMALE_MEDIUM_2015_2100.csv"),
                                            header = TRUE, check.names = FALSE,
                                            stringsAsFactors = FALSE)

input_popsize_male_1950to2015 <- read.csv(here(inpath_demography,
                                               countryname,
                                               "WPP2017_POP_F15_2_ANNUAL_POPULATION_BY_AGE_MALE_1950_2015.csv"),
                                          header = TRUE, check.names = FALSE,
                                          stringsAsFactors = FALSE)
input_popsize_male_2015to2100 <- read.csv(here(inpath_demography,
                                               countryname,
                                               "WPP2017_POP_F15_2_ANNUAL_POPULATION_BY_AGE_MALE_MEDIUM_2015_2100.csv"),
                                          header = TRUE, check.names = FALSE,
                                          stringsAsFactors = FALSE)


# Abridged life tables for mortality rates and survival ratio (to calculate migration rates)
input_lifetables_male_1950to2015 <- read.csv(here(inpath_demography,
                                                  countryname,
                                                  "WPP2017_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE_1950_2015.csv"),
                                             header = TRUE, check.names = FALSE,
                                             stringsAsFactors = FALSE)
input_lifetables_male_2015to2050 <- read.csv(here(inpath_demography,
                                                  countryname,
                                                  "WPP2017_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE_MEDIUM_2015_2050.csv"),
                                             header = TRUE, check.names = FALSE,
                                             stringsAsFactors = FALSE)
input_lifetables_male_2050to2100 <- read.csv(here(inpath_demography,
                                                  countryname,
                                                  "WPP2017_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE_MEDIUM_2050_2100.csv"),
                                             header = TRUE, check.names = FALSE,
                                             stringsAsFactors = FALSE)

input_lifetables_female_1950to2015 <- read.csv(here(inpath_demography,
                                                    countryname,
                                                    "WPP2017_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE_1950_2015.csv"),
                                               header = TRUE, check.names = FALSE,
                                               stringsAsFactors = FALSE)
input_lifetables_female_2015to2050 <- read.csv(here(inpath_demography,
                                                    countryname,
                                                    "WPP2017_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE_MEDIUM_2015_2050.csv"),
                                               header = TRUE, check.names = FALSE,
                                               stringsAsFactors = FALSE)
input_lifetables_female_2050to2100 <- read.csv(here(inpath_demography,
                                                    countryname,
                                                    "WPP2017_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE_MEDIUM_2050_2100.csv"),
                                               header = TRUE, check.names = FALSE,
                                               stringsAsFactors = FALSE)

# 5-year survival rates (survival ratio) from abridged life tables for migration rates

# Age-specific fertility rates 1950-2015
input_fertility_1950to2015 <- read.csv(here(inpath_demography,
                                            countryname,
                                            "WPP2017_FERT_F07_AGE_SPECIFIC_FERTILITY_1950_2015.csv"),
                                       header = TRUE, check.names = FALSE,
                                       stringsAsFactors = FALSE)
input_fertility_2015to2100 <- read.csv(here(inpath_demography,
                                            countryname,
                                            "WPP2017_FERT_F07_AGE_SPECIFIC_FERTILITY_MEDIUM_2015_2100.csv"),
                                       header = TRUE, check.names = FALSE,
                                       stringsAsFactors = FALSE)

# For validation
# Deaths by age
input_deaths_female_1950to2015 <- read.csv(here(inpath_demography,
                                                countryname,
                                                "WPP2017_MORT_F04_3_DEATHS_BY_AGE_FEMALE_1950_2015.csv"),
                                           header = TRUE, check.names = FALSE,
                                           stringsAsFactors = FALSE)
input_deaths_female_2015to2100 <- read.csv(here(inpath_demography,
                                                countryname,
                                                "WPP2017_MORT_F04_3_DEATHS_BY_AGE_FEMALE_MEDIUM_2015_2100.csv"),
                                           header = TRUE, check.names = FALSE,
                                           stringsAsFactors = FALSE)

input_deaths_male_1950to2015 <- read.csv(here(inpath_demography,
                                              countryname,
                                              "WPP2017_MORT_F04_2_DEATHS_BY_AGE_MALE_1950_2015.csv"),
                                         header = TRUE, check.names = FALSE,
                                         stringsAsFactors = FALSE)
input_deaths_male_2015to2100 <- read.csv(here(inpath_demography,
                                              countryname,
                                              "WPP2017_MORT_F04_2_DEATHS_BY_AGE_MALE_MEDIUM_2015_2100.csv"),
                                         header = TRUE, check.names = FALSE,
                                         stringsAsFactors = FALSE)

# Total births
input_births_1950to2015 <- read.csv(here(inpath_demography,
                                         countryname,
                                         "WPP2017_FERT_F01_BIRTHS_BOTH_SEXES_1950_2015.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)
input_births_2015to2100 <- read.csv(here(inpath_demography,
                                         countryname,
                                         "WPP2017_FERT_F01_BIRTHS_BOTH_SEXES_MEDIUM_2015_2100.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

