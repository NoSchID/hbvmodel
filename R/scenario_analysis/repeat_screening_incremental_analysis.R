# Incremental repeat screening comparison

# Population-level impact of repeat screening

require(here)  # for setting working directory
require(ggplot2)
require(tidyr)
require(dplyr)
require(gridExtra)
require(BCEA)
require(RColorBrewer)
source(here("R/imperial_model_interventions.R"))
source(here("R/scenario_analysis/calculate_outcomes.R"))

# Function to plot boxplot whiskers as 95% percentile
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# Load files (A1/E1 (10% coverage)/EA1 (50% screening coverage)) ----

out_path <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/repeat_screening_analysis/"
out_path_monit <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

# NEW FILE WITH DALYS (replace older ones later when everything has been re-simulated)

# Status quo
out2n <- readRDS(paste0(out_path, "out2_status_quo_301120.rds"))
out2n <- out2n[[1]]
out1n <- readRDS(paste0(out_path_monit, "a1_out1_status_quo_cohort_301120.rds"))
out1n <- out1n[[1]]
# This is comparison object for population-level screen in 2020
# No repeat screening (2020 only) 90% coverage
out3n <- readRDS(paste0(out_path_monit, "a1_out3_screen_2020_monit_0_011220.rds"))
out3n <- out3n[[1]]

# 2020 screen with 10-yearly monitoring, 90% coverage
out4 <- readRDS(paste0(out_path_monit, "a1_out4_screen_2020_monit_10_301120.rds"))
out4 <- out4[[1]]
out5 <- readRDS(paste0(out_path_monit, "a1_out5_screen_2020_monit_5_301120.rds"))
out5 <- out5[[1]]

# 2020 screen with no/5-yearly monitoring, 90% coverage, IT treated
out3_it <- readRDS(paste0(out_path_monit, "a1_it_out3_screen_2020_monit_0_180121.rds"))
out3_it <- out3_it[[1]]
out5_it <- readRDS(paste0(out_path_monit, "a1_it_out5_screen_2020_monit_5_161220.rds"))
out5_it <- out5_it[[1]]
monit_out7 <- readRDS(paste0(out_path_monit, "a1_it_monit_out7_161220.rds"))
monit_out7 <- monit_out7[[1]]

# Population-based repeat screen in 2030 with/without rescreening - with IT treated
out8a_it_2030_monit_0 <- readRDS(paste0(out_path, "a1_it_out8a_monit_0_screen_10a_2030_080121.rds"))
out8a_it_2030_monit_0 <- out8a_it_2030_monit_0[[1]]
out8a_it_2030_monit_5 <- readRDS(paste0(out_path, "a1_it_out8a_monit_5_screen_10a_2030_080121.rds"))
out8a_it_2030_monit_5 <- out8a_it_2030_monit_5[[1]]
out8b_it_2030_monit_0 <- readRDS(paste0(out_path, "a1_it_out8b_monit_0_screen_10b_2030_080121.rds"))
out8b_it_2030_monit_0 <- out8b_it_2030_monit_0[[1]]
out8b_it_2030_monit_5 <- readRDS(paste0(out_path, "a1_it_out8b_monit_5_screen_10b_2030_080121.rds"))
out8b_it_2030_monit_5 <- out8b_it_2030_monit_5[[1]]

# With repeat screening & 90% cov
out8b_2030n <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2030_041220.rds"))
out8b_2030n <- out8b_2030n[[1]]
out8b_2040n <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2040_031220.rds"))
out8b_2040n <- out8b_2040n[[1]]
out8b_2050n <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2050_031220.rds"))
out8b_2050n <- out8b_2050n[[1]]
out8b_2060n <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2060_031220.rds"))
out8b_2060n <- out8b_2060n[[1]]
out8b_2030_monit10n <- readRDS(paste0(out_path, "a1_out8b_monit_10_screen_10b_2030_071220.rds"))
out8b_2030_monit10n <- out8b_2030_monit10n[[1]]
out8b_2040_monit10n <- readRDS(paste0(out_path, "a1_out8b_monit_10_screen_10b_2040_071220.rds"))
out8b_2040_monit10n <- out8b_2040_monit10n[[1]]

# No repeat screening (2020 only), 50% coverage
out3_cov50n <- readRDS(paste0(out_path, "ea1_out3_screen_2020_monit_0_091220.rds"))
out3_cov50n <- out3_cov50n[[1]]
# 2020 screen with 10-yearly monitoring, 50% coverage
out4_cov50n <- readRDS(paste0(out_path, "ea1_out4_screen_2020_monit_10_091220.rds"))
out4_cov50n <- out4_cov50n[[1]]
# With repeat screening & 50% coverage
out8b_2030_cov50n <- readRDS(paste0(out_path, "ea1_out8b_monit_0_screen_10b_2030_101220.rds"))
out8b_2030_cov50n <- out8b_2030_cov50n[[1]]
out8b_2040_cov50n <- readRDS(paste0(out_path, "ea1_out8b_monit_0_screen_10b_2040_091220.rds"))
out8b_2040_cov50n <- out8b_2040_cov50n[[1]]
out8b_2030_monit10_cov50n <- readRDS(paste0(out_path, "ea1_out8b_monit_10_screen_10b_2030_101220.rds"))
out8b_2030_monit10_cov50n <- out8b_2030_monit10_cov50n[[1]]
out8b_2040_monit10_cov50n <- readRDS(paste0(out_path, "ea1_out8b_monit_10_screen_10b_2040_091220.rds"))
out8b_2040_monit10_cov50n <- out8b_2040_monit10_cov50n[[1]]

# Screening in antenatal care WITH IT TREATED (76% of pregnant women)
# One-time screen in 2020, no monitoring
anc_it_2020_monit_0 <- readRDS(paste0(out_path, "anc1_out3_screen_2020_monit_0_150121.rds"))
anc_it_2020_monit_0 <- anc_it_2020_monit_0[[1]]
# One-time screen in 2020, 5-yearly monitoring
#anc_it_2020_monit_5 <- readRDS(paste0(out_path, "anc1_out5_screen_2020_monit_5_050121.rds"))
#anc_it_2020_monit_5 <- anc_it_2020_monit_5[[1]]

# Regular screening at 76% coverage until specific year (2030/2040)
# With or without re-testing of those previously confirmed negative, no monitoring or every 5 years
# a=with rescreen, b=no rescreen
anc_it_2030_monit_0b <- readRDS(paste0(out_path, "anc1_2030_no_rescreen_monit_0_150121.rds"))
anc_it_2030_monit_0b <- anc_it_2030_monit_0b[[1]]
anc_it_2030_monit_5b <- readRDS(paste0(out_path, "anc1_2030_no_rescreen_monit_5_140121.rds"))
anc_it_2030_monit_5b <- anc_it_2030_monit_5b[[1]]
anc_it_2030_monit_0a <- readRDS(paste0(out_path, "anc1_2030_with_rescreen_monit_0_140121.rds"))
anc_it_2030_monit_0a <- anc_it_2030_monit_0a[[1]]
anc_it_2040_monit_0b <- readRDS(paste0(out_path, "anc1_2040_no_rescreen_monit_0_150121.rds"))
anc_it_2040_monit_0b <- anc_it_2040_monit_0b[[1]]
anc_it_2040_monit_5b <- readRDS(paste0(out_path, "anc1_2040_no_rescreen_monit_5_140121.rds"))
anc_it_2040_monit_5b <- anc_it_2040_monit_5b[[1]]
anc_it_2040_monit_0a <- readRDS(paste0(out_path, "anc1_2040_with_rescreen_monit_0_140121.rds"))
anc_it_2040_monit_0a <- anc_it_2040_monit_0a[[1]]

anc_it_2030_monit_0b_cohort <- readRDS(paste0(out_path, "anc1_2030_no_rescreen_monit_0_cohort_140121.rds"))
anc_it_2030_monit_0b_cohort <- anc_it_2030_monit_0b_cohort[[1]]
anc_it_2030_monit_0b_bdvacc <- readRDS(paste0(out_path, "anc1_2030_no_rescreen_monit_0_bdvacc_140121.rds"))
anc_it_2030_monit_0b_bdvacc <- anc_it_2030_monit_0b_bdvacc[[1]]

# Combination of population- and ANC screening (one-time pop followed by continuous ANC)
# IT is treated
pop_2020_anc_2030_monit_0b <- readRDS(paste0(out_path, "pop_2020_anc_2030_no_rescreen_monit_0_210121.rds"))
pop_2020_anc_2030_monit_0b <- pop_2020_anc_2030_monit_0b[[1]]
pop_2020_anc_2030_monit_5b <- readRDS(paste0(out_path, "pop_2020_anc_2030_no_rescreen_monit_5_210121.rds"))
pop_2020_anc_2030_monit_5b <- pop_2020_anc_2030_monit_5b[[1]]
pop_2020_anc_2030_monit_out7b <- readRDS(paste0(out_path, "pop_2020_anc_2030_no_rescreen_monit_sim7_210121.rds"))
pop_2020_anc_2030_monit_out7b <- pop_2020_anc_2030_monit_out7b[[1]]
pop_2020_anc_2030_monit_0a <- readRDS(paste0(out_path, "pop_2020_anc_2030_with_rescreen_monit_0_210121.rds"))
pop_2020_anc_2030_monit_0a <- pop_2020_anc_2030_monit_0a[[1]]
pop_2020_anc_2040_monit_0b <- readRDS(paste0(out_path, "pop_2020_anc_2040_no_rescreen_monit_0_210121.rds"))
pop_2020_anc_2040_monit_0b <- pop_2020_anc_2040_monit_0b[[1]]
pop_2020_anc_2040_monit_0a <- readRDS(paste0(out_path, "pop_2020_anc_2040_with_rescreen_monit_0_210121.rds"))
pop_2020_anc_2040_monit_0a <- pop_2020_anc_2040_monit_0a[[1]]
pop_2020_anc_2040_monit_5b <- readRDS(paste0(out_path, "pop_2020_anc_2040_no_rescreen_monit_5_220121.rds"))
pop_2020_anc_2040_monit_5b <- pop_2020_anc_2040_monit_5b[[1]]
pop_2020_anc_2040_monit_out7b <- readRDS(paste0(out_path, "pop_2020_anc_2040_no_rescreen_monit_sim7_220121.rds"))
pop_2020_anc_2040_monit_out7b <- pop_2020_anc_2040_monit_out7b[[1]]

# Screening in ANC until 2050 removed because lack of effect past 2030 already visible
#anc_it_2050_monit_5a <- readRDS(paste0(out_path, "anc1_2050_with_rescreen_monit_5_090121.rds"))
#anc_it_2050_monit_5a <- anc_it_2050_monit_5a[[1]]
#anc_it_2050_monit_0a <- readRDS(paste0(out_path, "anc1_2050_with_rescreen_monit_0_100121.rds"))
#anc_it_2050_monit_0a <- anc_it_2050_monit_0a[[1]]
#anc_it_2050_monit_0b <- readRDS(paste0(out_path, "anc1_2050_no_rescreen_monit_0_100121.rds"))
#anc_it_2050_monit_0b <- anc_it_2050_monit_0b[[1]]

# Yearly screening instead of continuous - don't use anymore as these were just temporary:
#anc_it_2030_monit_0b <- readRDS(paste0(out_path, "anc1_2030x_no_rescreen_monit_0_070121.rds"))
#anc_it_2030_monit_0b <- anc_it_2030_monit_0b[[1]]
#anc_it_2030_monit_5b <- readRDS(paste0(out_path, "anc1_2030x_no_rescreen_monit_5_070121.rds"))
#anc_it_2030_monit_5b <- anc_it_2030_monit_5b[[1]]
#anc_it_2030_monit_0a <- readRDS(paste0(out_path, "anc1_2030x_with_rescreen_monit_0_070121.rds"))
#anc_it_2030_monit_0a <- anc_it_2030_monit_0a[[1]]
#anc_it_2040_monit_0b <- readRDS(paste0(out_path, "anc1_2040x_no_rescreen_monit_0_070121.rds"))
#anc_it_2040_monit_0b <- anc_it_2040_monit_0b[[1]]

# More realistic implementation of the population-based programme (with IT treated):
# Gradual screening modelled as screening 10% of the population every year between 2020
# and 2028.
# No monitoring or every 5 years
e1_it_out3_gradual <- readRDS(paste0(out_path, "e1_it_out3_screen_2020_gradual_monit_0_060121.rds"))
e1_it_out3_gradual <- e1_it_out3_gradual[[1]]
e1_it_out5_gradual <- readRDS(paste0(out_path, "e1_it_out5_screen_2020_gradual_monit_5_050121.rds"))
e1_it_out5_gradual <- e1_it_out5_gradual[[1]]

# OLD FILES

# No monitoring:

# Status quo
#out2 <- readRDS(paste0(out_path, "out2_status_quo_180820.rds"))
#out2 <- out2[[1]]

# Screened untreated different cohorts
# No repeat screening
# 90% cov
#out1 <- readRDS(paste0(out_path_monit, "a1_out1_status_quo_cohort_240920.rds"))
#out1 <- out1[[1]]
# 10% cov
out8b_2020_cov10_sq <- readRDS(paste0(out_path, "e1_out1_screen_2020_status_quo_cohort_cov10_231120.rds"))
out8b_2020_cov10_sq <- out8b_2020_cov10_sq[[1]]
# With repeat screening
# 90% cov
out8b_2030_sq <- readRDS(paste0(out_path, "a1_out1_screen_2030_status_quo_cohort_231120.rds"))
out8b_2030_sq <- out8b_2030_sq[[1]]
out8b_2040_sq <- readRDS(paste0(out_path, "a1_out1_screen_2040_status_quo_cohort_231120.rds"))
out8b_2040_sq <- out8b_2040_sq[[1]]
# 10% cov
out8b_2030_cov10_sq <- readRDS(paste0(out_path, "e1_out1_screen_2030_status_quo_cohort_cov10_231120.rds"))
out8b_2030_cov10_sq <- out8b_2030_cov10_sq[[1]]
out8b_2040_cov10_sq <- readRDS(paste0(out_path, "e1_out1_screen_2040_status_quo_cohort_cov10_231120.rds"))
out8b_2040_cov10_sq <- out8b_2040_cov10_sq[[1]]

# No repeat screening (2020 only)
#out3 <- readRDS(paste0(out_path, "a1_out3_screen_2020_monit_0_240920.rds"))
#out3 <- out3[[1]]

# With repeat screening in the same age group (15-60) and rescreening of HBsAg-negatives
out8 <- readRDS(paste0(out_path, "a1_out8_monit_0_screen_10_290920.rds"))
out8 <- out8[[1]]
out9 <- readRDS(paste0(out_path, "a1_out9_monit_0_screen_5_290920.rds"))
out9 <- out9[[1]]
out10 <- readRDS(paste0(out_path, "a1_out10_monit_0_screen_1_290920.rds"))
out10 <- out10[[1]]

# With repeat screening in the unscreened age group (depends on repeat screening frequency)
out8a <- readRDS(paste0(out_path, "a1_out8a_monit_0_screen_10a_021020.rds"))
out8a <- out8a[[1]]
out9a <- readRDS(paste0(out_path, "a1_out9a_monit_0_screen_5a_021020.rds"))
out9a <- out9a[[1]]
out10a <- readRDS(paste0(out_path, "a1_out10a_monit_0_screen_1a_021020.rds"))
out10a <- out10a[[1]]

# With repeat screening in the same age group (15-60) but without rescreening of HBsAg-negatives
out8b <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_061020.rds"))
out8b <- out8b[[1]]
out9b <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_061020.rds"))
out9b <- out9b[[1]]
out10b <- readRDS(paste0(out_path, "a1_out10b_monit_0_screen_1b_061020.rds"))
out10b <- out10b[[1]]

### VARY NUMBER OF REPEAT SCREENING EVENTS:

# Vary number of repeat screening events for a 10-year frequency and without rescreening
out8b_2030 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2030_071020.rds"))
out8b_2030 <- out8b_2030[[1]]
out8b_2040 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2040_071020.rds"))
out8b_2040 <- out8b_2040[[1]]
out8b_2050 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2050_071020.rds"))
out8b_2050 <- out8b_2050[[1]]
out8b_2060 <- readRDS(paste0(out_path, "a1_out8b_monit_0_screen_10b_2060_071020.rds"))
out8b_2060 <- out8b_2060[[1]]
# 10-yearly monitoring
#out4 <- readRDS(paste0(out_path_monit, "a1_out4_screen_2020_monit_10_290920.rds"))  # One-off screen
#out4 <- out4[[1]]
out8b_2030_monit10 <- readRDS(paste0(out_path, "a1_out8b_monit_10_screen_10b_2030_201120.rds"))
out8b_2030_monit10 <- out8b_2030_monit10[[1]]
out8b_2040_monit10 <- readRDS(paste0(out_path, "a1_out8b_monit_10_screen_10b_2040_201120.rds"))
out8b_2040_monit10 <- out8b_2040_monit10[[1]]
# 5 yearly monitoring
#out5 <- readRDS(paste0(out_path_monit, "a1_out5_screen_2020_monit_5_201020.rds"))  # One-off screen
#out5 <- out5[[1]]
out8b_2030_monit5 <- readRDS(paste0(out_path, "a1_out8b_monit_5_screen_10b_2030_201120.rds"))
out8b_2030_monit5 <- out8b_2030_monit5[[1]]
out8b_2040_monit5 <- readRDS(paste0(out_path, "a1_out8b_monit_5_screen_10b_2040_231120.rds"))
out8b_2040_monit5 <- out8b_2040_monit5[[1]]

## Vary number of repeat screening events for a 10-year frequency and without rescreening
# 10% screening coverage!
out8b_2020_cov10 <- readRDS(paste0(out_path, "e1_out3_screen_2020_monit_0_191120.rds"))
out8b_2020_cov10 <- out8b_2020_cov10[[1]]
out8b_2030_cov10 <- readRDS(paste0(out_path, "e1_out8b_monit_0_screen_10b_2030_191120.rds"))
out8b_2030_cov10 <- out8b_2030_cov10[[1]]
out8b_2040_cov10 <- readRDS(paste0(out_path, "e1_out8b_monit_0_screen_10b_2040_191120.rds"))
out8b_2040_cov10 <- out8b_2040_cov10[[1]]
out8b_2050_cov10 <- readRDS(paste0(out_path, "e1_out8b_monit_0_screen_10b_2050_191120.rds"))
out8b_2050_cov10 <- out8b_2050_cov10[[1]]
out8b_2060_cov10 <- readRDS(paste0(out_path, "e1_out8b_monit_0_screen_10b_2060_191120.rds"))
out8b_2060_cov10 <- out8b_2060_cov10[[1]]
# 10-yearly monitoring
out8b_2020_cov10_monit10 <- readRDS(paste0(out_path, "e1_out4_screen_2020_monit_10_201120.rds"))
out8b_2020_cov10_monit10 <- out8b_2020_cov10_monit10[[1]]
out8b_2030_cov10_monit10 <- readRDS(paste0(out_path, "e1_out8b_monit_10_screen_10b_2030_201120.rds"))
out8b_2030_cov10_monit10 <- out8b_2030_cov10_monit10[[1]]
out8b_2040_cov10_monit10 <- readRDS(paste0(out_path, "e1_out8b_monit_10_screen_10b_2040_201120.rds"))
out8b_2040_cov10_monit10 <- out8b_2040_cov10_monit10[[1]]
# 5-yearly monitoring
out8b_2020_cov10_monit5 <- readRDS(paste0(out_path, "e1_out5_screen_2020_monit_5_201120.rds"))
out8b_2020_cov10_monit5 <- out8b_2020_cov10_monit5[[1]]
out8b_2030_cov10_monit5 <- readRDS(paste0(out_path, "e1_out8b_monit_5_screen_10b_2030_201120.rds"))
out8b_2030_cov10_monit5 <- out8b_2030_cov10_monit5[[1]]
out8b_2040_cov10_monit5 <- readRDS(paste0(out_path, "e1_out8b_monit_5_screen_10b_2040_201120.rds"))
out8b_2040_cov10_monit5 <- out8b_2040_cov10_monit5[[1]]

# 50% screening coverage!
out8b_2020_cov50 <- readRDS(paste0(out_path, "ea1_out3_screen_2020_monit_0_251120.rds"))
out8b_2020_cov50 <- out8b_2020_cov50[[1]]
out8b_2030_cov50 <- readRDS(paste0(out_path, "ea1_out8b_monit_0_screen_10b_2030_251120.rds"))
out8b_2030_cov50 <- out8b_2030_cov50[[1]]
out8b_2040_cov50 <- readRDS(paste0(out_path, "ea1_out8b_monit_0_screen_10b_2040_251120.rds"))
out8b_2040_cov50 <- out8b_2040_cov50[[1]]
# 10-yearly monitoring
out8b_2020_cov50_monit10 <- readRDS(paste0(out_path, "ea1_out4_screen_2020_monit_10_251120.rds"))
out8b_2020_cov50_monit10 <- out8b_2020_cov50_monit10[[1]]
out8b_2030_cov50_monit10 <- readRDS(paste0(out_path, "ea1_out8b_monit_10_screen_10b_2030_251120.rds"))
out8b_2030_cov50_monit10 <- out8b_2030_cov50_monit10[[1]]
out8b_2040_cov50_monit10 <- readRDS(paste0(out_path, "ea1_out8b_monit_10_screen_10b_2040_251120.rds"))
out8b_2040_cov50_monit10 <- out8b_2040_cov50_monit10[[1]]
# 5-yearly monitoring
out8b_2020_cov50_monit5 <- readRDS(paste0(out_path, "ea1_out5_screen_2020_monit_5_251120.rds"))
out8b_2020_cov50_monit5 <- out8b_2020_cov50_monit5[[1]]
out8b_2030_cov50_monit5 <- readRDS(paste0(out_path, "ea1_out8b_monit_5_screen_10b_2030_251120.rds"))
out8b_2030_cov50_monit5 <- out8b_2030_cov50_monit5[[1]]
out8b_2040_cov50_monit5 <- readRDS(paste0(out_path, "ea1_out8b_monit_5_screen_10b_2040_251120.rds"))
out8b_2040_cov50_monit5 <- out8b_2040_cov50_monit5[[1]]

# Vary number of repeat screening events for a 5-year frequency and without rescreening
out9b_2030 <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_2030_121020.rds"))
out9b_2030 <- out9b_2030[[1]]
out9b_2040 <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_2040_121020.rds"))
out9b_2040 <- out9b_2040[[1]]
out9b_2050 <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_2050_121020.rds"))
out9b_2050 <- out9b_2050[[1]]
out9b_2060 <- readRDS(paste0(out_path, "a1_out9b_monit_0_screen_5b_2060_121020.rds"))
out9b_2060 <- out9b_2060[[1]]

# Access channels analysis
# Workplace screening
#out3_wpl <-  readRDS(paste0(out_path, "wpl1_out3_screen_2020_monit_0_301020.rds"))
#out3_wpl <- out3_wpl[[1]]
# ANC screening
#out3_anc <-  readRDS(paste0(out_path, "anc1_out3_screen_2020_monit_0_021120.rds"))
#out3_anc <- out3_anc[[1]]


# Labels
timepoints <- c(2050,2100)
period_labs <- c(paste0("2020-",timepoints[1]), paste0("2020-",timepoints[2]))
names(period_labs) <- c(as.character(timepoints[1]), as.character(timepoints[2]))






# For 90% coverage and no monitoring, is a repeat screen cost-effective? ----
# Comparing: out3, out8b_2030, out8b_2040, out8b_2050
interactions_no_monit <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        left_join(left_join(left_join(gather(out3n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out3n$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out3n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_0_screen_10b_2030",
        left_join(left_join(left_join(gather(out8b_2030n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2030n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2030n$interactions[[16]]$total_assessed[-c(1:3)]-out8b_2030n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2030n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_0_screen_10b_2040",
        left_join(left_join(left_join(gather(out8b_2040n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2040n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2040n$interactions[[16]]$total_assessed[-c(1:3)]-out8b_2040n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2040n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_0_screen_10b_2050",
        left_join(left_join(left_join(gather(out8b_2050n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2050n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2050n$interactions[[16]]$total_assessed[-c(1:3)]-out8b_2050n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2050n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_10",
        left_join(left_join(left_join(gather(out4$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out4$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out4$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_5",
        left_join(left_join(left_join(gather(out5$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out5$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out5$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim"))
)

py_on_treatment_no_monit <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             sim = names(out3n$py_on_treatment[[16]]),
             py_on_treatment = out3n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_0_screen_10b_2030",
             sim = names(out8b_2030n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2030n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_0_screen_10b_2040",
             sim = names(out8b_2040n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2040n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_0_screen_10b_2050",
             sim = names(out8b_2050n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2050n$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_10",
             sim = names(out4$py_on_treatment[[16]]),
             py_on_treatment = out4$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_5",
             sim = names(out5$py_on_treatment[[16]]),
             py_on_treatment = out5$py_on_treatment[[16]])
)

# Outcome 1: DALYs
dalys_averted_no_monit <-
  plot_hbv_deaths_averted(counterfactual_object = out2n,
                          scenario_objects = list(out3n,out8b_2030n,out8b_2040n,out8b_2050n, out4, out5),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
dalys_averted_no_monit <- subset(dalys_averted_no_monit, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
dalys_averted_no_monit$sim <- gsub("[^0-9]", "", dalys_averted_no_monit$sim)

# Outcome 2: HBV-related deaths
deaths_averted_no_monit <-
  plot_hbv_deaths_averted(counterfactual_object = out2n,
                          scenario_objects = list(out3n,out8b_2030n,out8b_2040n,out8b_2050n, out4, out5),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
deaths_averted_no_monit <- subset(deaths_averted_no_monit, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
deaths_averted_no_monit$sim <- gsub("[^0-9]", "", deaths_averted_no_monit$sim)

# Combine into dataframe (No discounting)
incremental_df_no_monit <- create_incremental_plot_df(interactions_df=interactions_no_monit,
                                      py_on_treatment_df=py_on_treatment_no_monit,
                                      deaths_averted_df=deaths_averted_no_monit,
                                      ly_saved_df = dalys_averted_no_monit, # replace LY by DALYs
                                      hbsag_test_cost = 8.3,
                                      clinical_assessment_cost = 33,
                                      monitoring_assessment_cost = 25.5,
                                      treatment_py_cost = 66.5,
                                      ref_label = "No treatment")
colnames(incremental_df_no_monit)[colnames(incremental_df_no_monit)=="ly_saved"] <- "dalys_averted"

incremental_df_no_monit2 <- incremental_df_no_monit
# Remove monitoring in the 2020 screen
incremental_df_no_monit <- subset(incremental_df_no_monit, scenario != "screen_2020_monit_10" &
                                    scenario != "screen_2020_monit_5")

# Fid dominated strategies and ICERs
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- incremental_df_no_monit[which(incremental_df_no_monit$sim==
                                              unique(incremental_df_no_monit$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))

# All strategies 100% non-dominated

icer_list <- list()
for(i in 1:183) {
  print(i)
  icer_list[[i]] <- incremental_df_no_monit[which(incremental_df_no_monit$sim==
                                    unique(incremental_df_no_monit$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result

ggplot(subset(incremental_df_no_monit, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_point(data= group_by(incremental_df_no_monit, scenario) %>% summarise(dalys_averted=median(dalys_averted),
                                                                    total_cost = median(total_cost)),
             aes(x=dalys_averted, y = total_cost, group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), size = 5) +
  scale_fill_manual(values = rev(brewer.pal(4,"RdYlBu"))) +
  scale_colour_manual("Screening periods\n(no monitoring)",
                      labels = c("screen_2020_monit_0" = "2020",
                                 "monit_0_screen_10b_2030" = "2020-2030",
                                 "monit_0_screen_10b_2040" = "2020-2040",
                                 "monit_0_screen_10b_2050" = "2020-2050"),
                      values = c("black", rev(brewer.pal(4,"RdYlBu")))) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  geom_abline(slope=391, linetype = "dashed") +
  geom_abline(slope=518, linetype = "dashed") +
  theme_bw()

# 2030 repeat screen is cost-effective if no monitoring is available.
# Does this remain the case if monitoring becomes an alternative option?
# Fid dominated strategies and ICERs
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- incremental_df_no_monit2[which(incremental_df_no_monit2$sim==
                                                              unique(incremental_df_no_monit2$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# 10- and 5-yearly monitoring is dominated! So repeat screen remains cost-effective.

ggplot(subset(incremental_df_no_monit2, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_point(data= group_by(incremental_df_no_monit2, scenario) %>% summarise(dalys_averted=median(dalys_averted),
                                                                             total_cost = median(total_cost)),
             aes(x=dalys_averted, y = total_cost, group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), size = 5) +
  scale_fill_manual(values = rev(brewer.pal(6,"RdYlBu"))) +
  scale_colour_manual("Screening periods\n(no monitoring)",
                      labels = c("screen_2020_monit_0" = "2020",
                                 "monit_0_screen_10b_2030" = "2020-2030",
                                 "monit_0_screen_10b_2040" = "2020-2040",
                                 "monit_0_screen_10b_2050" = "2020-2050"),
                      values = c("black", rev(brewer.pal(6,"RdYlBu")))) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  geom_abline(slope=391, linetype = "dashed") +
  geom_abline(slope=518, linetype = "dashed") +
  theme_bw()

# This shows that 2030 repeat screen is slightly more cost-effective (though almost the same)
# than one-time screen with 10-yearly monitoring. 2040 repeat screen was below WTP but
# still dominates one-time screen with 5-yearly monitoring.

# Nevertheless, will likely find that combination of repeat screen and monitoring is the most cost-effective (see next section)

# For 90% coverage and no monitoring, compare monitoring vs repeat screen (no discounting) ----
interactions <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        left_join(left_join(left_join(gather(out3n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out3n$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out3n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_0_screen_10b_2030",
        left_join(left_join(left_join(gather(out8b_2030n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2030n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2030n$interactions[[16]]$total_assessed[-c(1:3)]-out8b_2030n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2030n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_0_screen_10b_2040",
        left_join(left_join(left_join(gather(out8b_2040n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2040n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2040n$interactions[[16]]$total_assessed[-c(1:3)]-out8b_2040n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2040n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_10",
        left_join(left_join(left_join(gather(out4$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out4$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out4$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_10_screen_10b_2030",
        left_join(left_join(left_join(gather(out8b_2030_monit10n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2030n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2030_monit10n$interactions[[16]]$total_assessed[-c(1:3)]-
                                     out8b_2030n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2030_monit10n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_10_screen_10b_2040",
        left_join(left_join(left_join(gather(out8b_2040_monit10n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2040n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2040_monit10n$interactions[[16]]$total_assessed[-c(1:3)]-
                                     out8b_2040n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2040_monit10n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim"))
)



py_on_treatment <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             sim = names(out3n$py_on_treatment[[16]]),
             py_on_treatment = out3n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_0_screen_10b_2030",
             sim = names(out8b_2030n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2030n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_0_screen_10b_2040",
             sim = names(out8b_2040n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2040n$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_10",
             sim = names(out4$py_on_treatment[[16]]),
             py_on_treatment = out4$py_on_treatment[[16]]),
  data.frame(scenario = "monit_10_screen_10b_2030",
             sim = names(out8b_2030_monit10n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2030_monit10n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_10_screen_10b_2040",
             sim = names(out8b_2040_monit10n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2040_monit10n$py_on_treatment[[16]])
)

# Outcome 1: DALYs
dalys_averted <-
  plot_hbv_deaths_averted(counterfactual_object = out2n,
                          scenario_objects = list(out3n,out8b_2030n,out8b_2040n,
                                                  out4, out8b_2030_monit10n,
                                                  out8b_2040_monit10n),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
dalys_averted <- subset(dalys_averted, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
dalys_averted$sim <- gsub("[^0-9]", "", dalys_averted$sim)

# For cohort DALYS, would need to compare repeat screens to a different out1 (with repeat screen)

# Outcome 2: HBV-related deaths
deaths_averted <-
  plot_hbv_deaths_averted(counterfactual_object = out2n,
                          scenario_objects = list(out3n,out8b_2030n,out8b_2040n,
                                                  out4, out8b_2030_monit10n,
                                                  out8b_2040_monit10n),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
deaths_averted <- subset(deaths_averted, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
deaths_averted$sim <- gsub("[^0-9]", "", deaths_averted$sim)

# Combine into dataframe (No discounting)
incremental_df <- create_incremental_plot_df(interactions_df=interactions,
                                                      py_on_treatment_df=py_on_treatment,
                                                      deaths_averted_df=deaths_averted,
                                                      ly_saved_df = dalys_averted, # replace LY by DALYs
                                                      hbsag_test_cost = 8.3,
                                                      clinical_assessment_cost = 84.4,
                                                      monitoring_assessment_cost = 40.1,
                                                      treatment_py_cost = 60,
                                                      ref_label = "No treatment")
colnames(incremental_df)[colnames(incremental_df)=="ly_saved"] <- "dalys_averted"


# Find dominated strategies and ICERs
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- incremental_df[which(incremental_df$sim==
                                                              unique(incremental_df$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Only 2020 screen with no monitoring and 2030 and 2040 screens with 10-yearly monitoring
# has >50% probability of being non-dominated.

incremental_df_non_dom <- subset(incremental_df, scenario %in% c("screen_2020_monit_0",
                                                                 "monit_10_screen_10b_2030",
                                                                 "monit_10_screen_10b_2040"))

icer_list <- list()
for(i in 1:183) {
  print(i)
  icer_list[[i]] <- incremental_df_non_dom[which(incremental_df_non_dom$sim==
                                                    unique(incremental_df_non_dom$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result
# This suggests that monitoring becomes more cost-effective if there are also repeat screens?

ggplot(subset(incremental_df, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_point(data= group_by(incremental_df, scenario) %>% summarise(dalys_averted=median(dalys_averted),
                                                                    total_cost = median(total_cost)),
             aes(x=dalys_averted, y = total_cost, group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), size = 5) +
 scale_fill_manual(values = rev(brewer.pal(6,"RdYlBu"))) +
 scale_colour_manual("Scenarios",
                     labels = c("screen_2020_monit_0" = "2020 screen, no monitoring",
                                "screen_2020_monit_10" = "2020 screen, monitor every 10 years",
                                "monit_0_screen_10b_2030" = "2020+2030 screen, no monitoring",
                                "monit_0_screen_10b_2040" = "2020+2030+2040 screen, no monitoring",
                                "monit_10_screen_10b_2030" = "2020+2030 screen, monitor every 10 years",
                                "monit_10_screen_10b_2040" = "2020+2030+2040 screen, monitor every 10 years"),
                      values = c("black", rev(brewer.pal(6,"RdYlBu")))) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  geom_abline(slope=391, linetype = "dashed") +
  geom_abline(slope=518, linetype = "dashed") +
  theme_bw()

# Most cost-effective startegy <WTP was 2040 repeat screen+ 10-yearly monitoring.
# Plot shows that repeat screen + monitoring dominates repeat screen without monitoring as well as
# one-time screen with monitoring. Need to confirm this with discounting.

# Why does monitoring avert more DALYs if there is more than 1 screen?
# In the cohort
x1 <- plot_hbv_deaths_averted_cohort(counterfactual_object = out3n,
                        scenario_objects = list(out4),
                        outcome_to_avert = "cohort_dalys",
                        outcome_to_plot = "proportion_averted",
                        counterfactual_label = "no treatment")
x2 <- plot_hbv_deaths_averted_cohort(counterfactual_object = out8b_2030n,
                        scenario_objects = list(out8b_2030_monit10n),
                        outcome_to_avert = "cohort_dalys",
                        outcome_to_plot = "proportion_averted",
                        counterfactual_label = "no treatment")
# In the population
x3 <- plot_hbv_deaths_averted(counterfactual_object = out3n,
                                     scenario_objects = list(out4),
                                     outcome_to_avert = "dalys",
                                     outcome_to_plot = "proportion_averted",
                                     counterfactual_label = "no treatment")

x4 <- plot_hbv_deaths_averted(counterfactual_object = out8b_2030n,
                                     scenario_objects = list(out8b_2030_monit10n),
                                     outcome_to_avert = "dalys",
                                     outcome_to_plot = "proportion_averted",
                                     counterfactual_label = "no treatment")

quantile(x1$value[x1$type=="proportion_averted"], c(0.5,0.025,0.975))
quantile(x2$value[x2$type=="proportion_averted"], c(0.5,0.025,0.975))
quantile(x1$value[x1$type=="number_averted"], c(0.5,0.025,0.975))
quantile(x2$value[x2$type=="number_averted"], c(0.5,0.025,0.975))
# The proportion of DALYs averted in the diagnosed cohort by monitoring is pretty much the same
# whether there is 1 or 2 screens.
# However the number averted is quite a bit higher because there are more people in the cohort
# to be monitored.
quantile(x3$value[x3$type=="proportion_averted" & x3$by_year == 2100], c(0.5,0.025,0.975))
quantile(x4$value[x4$type=="proportion_averted"& x4$by_year == 2100], c(0.5,0.025,0.975))
quantile(x3$value[x3$type=="number_averted" & x3$by_year == 2100], c(0.5,0.025,0.975))
quantile(x4$value[x4$type=="number_averted"& x4$by_year == 2100], c(0.5,0.025,0.975))
# Therefore in the population, the proportion  of DALYs averted by monitoring is a little higher if there
# is a repeat screen.

# What about costs?
quantile(incremental_df$total_cost[incremental_df$scenario == "monit_10_screen_10b_2030"]-
  incremental_df$total_cost[incremental_df$scenario == "monit_0_screen_10b_2030"],
  c(0.5,0.025,0.975))
quantile(incremental_df$total_cost[incremental_df$scenario == "screen_2020_monit_10"]-
           incremental_df$total_cost[incremental_df$scenario == "screen_2020_monit_0"],
         c(0.5,0.025,0.975))
# Absolute cost is slightly higher with repeat screen because more people to monitor.
quantile((incremental_df$total_cost[incremental_df$scenario == "monit_10_screen_10b_2030"]-
           incremental_df$total_cost[incremental_df$scenario == "monit_0_screen_10b_2030"])/
           incremental_df$total_cost[incremental_df$scenario == "monit_0_screen_10b_2030"],
         c(0.5,0.025,0.975))
quantile((incremental_df$total_cost[incremental_df$scenario == "screen_2020_monit_10"]-
           incremental_df$total_cost[incremental_df$scenario == "screen_2020_monit_0"])/
           incremental_df$total_cost[incremental_df$scenario == "screen_2020_monit_0"],
         c(0.5,0.025,0.975))
# Proportionally the cost of monitoring is higher if there is only the 2020 screen.
quantile((incremental_df$treatment_initiations[incremental_df$scenario == "monit_10_screen_10b_2030"]-
            incremental_df$treatment_initiations[incremental_df$scenario == "monit_0_screen_10b_2030"])/
           incremental_df$treatment_initiations[incremental_df$scenario == "monit_0_screen_10b_2030"],
         c(0.5,0.025,0.975))
quantile((incremental_df$treatment_initiations[incremental_df$scenario == "screen_2020_monit_10"]-
            incremental_df$treatment_initiations[incremental_df$scenario == "screen_2020_monit_0"])/
           incremental_df$treatment_initiations[incremental_df$scenario == "screen_2020_monit_0"],
         c(0.5,0.025,0.975))
# However the person-time in treatment as a result of monitoring is proportionally almost the same.
# Though the treatment initiations are proportionally higher in the one-time screen.
# Monitoring initiates a larger proportion of the cohort on treatment in the first screen than
# if there are 2 screens.

# For 90% coverage and no monitoring, compare monitoring vs repeat screen (with discounting) ----
annual_discounting_rate <- 0.03

assemble_discounted_interactions_for_screening_strategies <- function(scenario_object,
                                                                        discount_rate = annual_discounting_rate,
                                                                      assessment_object) {
  # Compares monitoring to out3 and everything else to status quo
  # at discount rate of 3%
  out <- left_join(
    left_join(
      left_join(discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_screened",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "hbsag_tests"),
                discount_outcome_2020_to_2100(scenario_object=assessment_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_assessed",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "clinical_assessments")),
      discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                    object_to_subtract=assessment_object,
                                    outcome="interactions",
                                    interaction_outcome="total_assessed",
                                    yearly_discount_rate=discount_rate,
                                    interaction_colname = "monitoring_assessments")),
    discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                  object_to_subtract=NULL,
                                  outcome="interactions",
                                  interaction_outcome="total_treated",
                                  yearly_discount_rate=discount_rate,
                                  interaction_colname = "treatment_initiations"))

  return(out)
}

interactions_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        assemble_discounted_interactions_for_screening_strategies(out3n, assessment_object = out3n)),
  cbind(scenario = "monit_0_screen_10b_2030",
        assemble_discounted_interactions_for_screening_strategies(out8b_2030n, assessment_object = out8b_2030n)),
  cbind(scenario = "monit_0_screen_10b_2040",
        assemble_discounted_interactions_for_screening_strategies(out8b_2040n, assessment_object = out8b_2040n)),
  cbind(scenario = "screen_2020_monit_10",
        assemble_discounted_interactions_for_screening_strategies(out4, assessment_object = out3n)),
  cbind(scenario = "monit_10_screen_10b_2030",
        assemble_discounted_interactions_for_screening_strategies(out8b_2030_monit10n, assessment_object = out8b_2030n)),
  cbind(scenario = "monit_10_screen_10b_2040",
        assemble_discounted_interactions_for_screening_strategies(out8b_2040_monit10n, assessment_object = out8b_2040n))
)

py_on_treatment_disc <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out3n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out8b_2030n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2040",
             discount_outcome_2020_to_2100(scenario_object=out8b_2040n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_10",
             discount_outcome_2020_to_2100(scenario_object=out4,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out8b_2030_monit10n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2040",
             discount_outcome_2020_to_2100(scenario_object=out8b_2040_monit10n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate))
)
py_on_treatment_disc$sim <- gsub("[^0-9]", "", py_on_treatment_disc$sim)

dalys_averted_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3n,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2030n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2040",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2040n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_10",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out4,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2030_monit10n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2040",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2040_monit10n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate))
)
dalys_averted_disc$sim <- gsub("[^0-9]", "", dalys_averted_disc$sim)
colnames(dalys_averted_disc)[colnames(dalys_averted_disc) == "dalys"] <-
  "value"

deaths_averted_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3n,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2030n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2040",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2040n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_10",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out4,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2030_monit10n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2040",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2040_monit10n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate))
)
deaths_averted_disc$sim <- gsub("[^0-9]", "", deaths_averted_disc$sim)
colnames(deaths_averted_disc)[colnames(deaths_averted_disc) == "cum_hbv_deaths"] <-
  "value"

# Combine into dataframe (No discounting)
incremental_df_disc <- create_incremental_plot_df(interactions_df=interactions_disc,
                                             py_on_treatment_df=py_on_treatment_disc,
                                             deaths_averted_df=deaths_averted_disc,
                                             ly_saved_df = dalys_averted_disc, # replace LY by DALYs
                                             hbsag_test_cost = 8.3,
                                             clinical_assessment_cost = 84.4,
                                             monitoring_assessment_cost = 40.1,
                                             treatment_py_cost = 60,
                                             ref_label = "No treatment")
colnames(incremental_df_disc)[colnames(incremental_df_disc)=="ly_saved"] <- "dalys_averted"

# Find dominated strategies and ICERs
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- incremental_df_disc[which(incremental_df_disc$sim==
                                                     unique(incremental_df_disc$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Only 2020 screen with no monitoring and 2030 and 2040 screens with 10-yearly monitoring
# has >50% probability of being non-dominated.

incremental_df_non_dom_disc <- subset(incremental_df_disc, scenario %in% c("screen_2020_monit_0",
                                                                 "monit_10_screen_10b_2030",
                                                                 "monit_10_screen_10b_2040"))

icer_list <- list()
for(i in 1:183) {
  print(i)
  icer_list[[i]] <- incremental_df_non_dom_disc[which(incremental_df_non_dom_disc$sim==
                                                   unique(incremental_df_non_dom_disc$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result

# CER plot
ggplot(subset(incremental_df_disc, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_point(data= group_by(incremental_df_disc, scenario) %>% summarise(dalys_averted=median(dalys_averted),
                                                                         total_cost = median(total_cost)),
             aes(x=dalys_averted, y = total_cost, group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), size = 5) +
  scale_fill_manual(values = rev(brewer.pal(7,"RdYlBu"))) +
  scale_colour_manual("Scenarios",
                      labels = c("screen_2020_monit_0" = "2020 screen, no monitoring",
                                 "screen_2020_monit_10" = "2020 screen, monitor every 10 years",
                                 "monit_0_screen_10b_2030" = "2020+2030 screen, no monitoring",
                                 "monit_0_screen_10b_2040" = "2020+2030+2040 screen, no monitoring",
                                 "monit_10_screen_10b_2030" = "2020+2030 screen, monitor every 10 years",
                                 "monit_10_screen_10b_2040" = "2020+2030+2040 screen, monitor every 10 years"),
                      values = c("black", rev(brewer.pal(7,"RdYlBu")))) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  geom_abline(slope=391, linetype = "dashed") +
  geom_abline(slope=518, linetype = "dashed") +
  theme_bw()


# For 50% coverage and no monitoring, compare monitoring vs repeat screen (no discounting) ----
interactions <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        left_join(left_join(left_join(gather(out3n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out3n$interactions[[16]]$total_assessed[-c(1:3)]-out3n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out3n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_0_cov50",
        left_join(left_join(left_join(gather(out3_cov50n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out3_cov50n$interactions[[16]]$total_assessed[-c(1:3)]-out3_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out3_cov50n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_0_screen_10b_2030_cov50",
        left_join(left_join(left_join(gather(out8b_2030_cov50n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2030_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2030_cov50n$interactions[[16]]$total_assessed[-c(1:3)]-out8b_2030_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2030_cov50n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_0_screen_10b_2040_cov50",
        left_join(left_join(left_join(gather(out8b_2040_cov50n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2040_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2040_cov50n$interactions[[16]]$total_assessed[-c(1:3)]-out8b_2040_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2040_cov50n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "screen_2020_monit_10_cov50",
        left_join(left_join(left_join(gather(out4_cov50n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out3_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out4_cov50n$interactions[[16]]$total_assessed[-c(1:3)]-out3_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out4_cov50n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_10_screen_10b_2030_cov50",
        left_join(left_join(left_join(gather(out8b_2030_monit10_cov50n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2030_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2030_monit10_cov50n$interactions[[16]]$total_assessed[-c(1:3)]-
                                     out8b_2030_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2030_monit10_cov50n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim")),
  cbind(scenario = "monit_10_screen_10b_2040_cov50",
        left_join(left_join(left_join(gather(out8b_2040_monit10_cov50n$interactions[[16]]$total_screened[-c(1:3)],
                                             key = "sim", value = "hbsag_tests"),
                                      gather(out8b_2040_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                             key = "sim", value = "clinical_assessments"), by = "sim"),
                            gather(out8b_2040_monit10_cov50n$interactions[[16]]$total_assessed[-c(1:3)]-
                                     out8b_2040_cov50n$interactions[[16]]$total_assessed[-c(1:3)],
                                   key = "sim", value = "monitoring_assessments"), by = "sim"),
                  gather(out8b_2040_monit10_cov50n$interactions[[16]]$total_treated[-c(1:3)],
                         key = "sim", value = "treatment_initiations"), by = "sim"))
)


py_on_treatment <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             sim = names(out3n$py_on_treatment[[16]]),
             py_on_treatment = out3n$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_0_cov50",
             sim = names(out3_cov50n$py_on_treatment[[16]]),
             py_on_treatment = out3_cov50n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_0_screen_10b_2030_cov50",
             sim = names(out8b_2030_cov50n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2030_cov50n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_0_screen_10b_2040_cov50",
             sim = names(out8b_2040_cov50n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2040_cov50n$py_on_treatment[[16]]),
  data.frame(scenario = "screen_2020_monit_10_cov50",
             sim = names(out4_cov50n$py_on_treatment[[16]]),
             py_on_treatment = out4_cov50n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_10_screen_10b_2030_cov50",
             sim = names(out8b_2030_monit10_cov50n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2030_monit10_cov50n$py_on_treatment[[16]]),
  data.frame(scenario = "monit_10_screen_10b_2040_cov50",
             sim = names(out8b_2040_monit10_cov50n$py_on_treatment[[16]]),
             py_on_treatment = out8b_2040_monit10_cov50n$py_on_treatment[[16]])
)

# Outcome 1: DALYs
dalys_averted <-
  plot_hbv_deaths_averted(counterfactual_object = out2n,
                          scenario_objects = list(out3n, out3_cov50n,out8b_2030_cov50n,
                                                  out8b_2040_cov50n,
                                                  out4_cov50n, out8b_2030_monit10_cov50n,
                                                  out8b_2040_monit10_cov50n),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
dalys_averted <- subset(dalys_averted, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
dalys_averted$sim <- gsub("[^0-9]", "", dalys_averted$sim)

# For cohort DALYS, would need to compare repeat screens to a different out1 (with repeat screen)

# Outcome 2: HBV-related deaths
deaths_averted <-
  plot_hbv_deaths_averted(counterfactual_object = out2n,
                          scenario_objects = list(out3n, out3_cov50n,out8b_2030_cov50n,
                                                  out8b_2040_cov50n,
                                                  out4_cov50n, out8b_2030_monit10_cov50n,
                                                  out8b_2040_monit10_cov50n),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
deaths_averted <- subset(deaths_averted, type == "number_averted" & by_year == 2100) %>%
  select(scenario, sim, value)
deaths_averted$sim <- gsub("[^0-9]", "", deaths_averted$sim)

# Combine into dataframe (No discounting)
incremental_df <- create_incremental_plot_df(interactions_df=interactions,
                                             py_on_treatment_df=py_on_treatment,
                                             deaths_averted_df=deaths_averted,
                                             ly_saved_df = dalys_averted, # replace LY by DALYs
                                             hbsag_test_cost = 8.3,
                                             clinical_assessment_cost = 84.4,
                                             monitoring_assessment_cost = 40.1,
                                             treatment_py_cost = 60,
                                             ref_label = "No treatment")
colnames(incremental_df)[colnames(incremental_df)=="ly_saved"] <- "dalys_averted"


# Find dominated strategies and ICERs
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- incremental_df[which(incremental_df$sim==
                                                     unique(incremental_df$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()

# Only 2020 screen with no monitoring and 2030 and 2040 screens with 10-yearly monitoring
# has >50% probability of being non-dominated. (same as for 90% cov)
# What if additionally adding the one-time screen without monitoring at 90% coverage?
# Then the only non-dominated strategies are the one-time 90% coverage and the 2040
# repeat screen with 10-yearly monitoring at 50% coverage!

#incremental_df_non_dom <- subset(incremental_df, scenario %in% c("screen_2020_monit_0_cov50",
#                                                                 "monit_10_screen_10b_2030_cov50",
#                                                                 "monit_10_screen_10b_2040_cov50"))

incremental_df_non_dom <- subset(incremental_df, scenario %in% c("screen_2020_monit_0",
                                                                 "monit_10_screen_10b_2040_cov50"))

icer_list <- list()
for(i in 1:183) {
  print(i)
  icer_list[[i]] <- incremental_df_non_dom[which(incremental_df_non_dom$sim==
                                                   unique(incremental_df_non_dom$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result
# This suggests that monitoring becomes more cost-effective if there are also repeat screens?
# If 90% coverage one-time is included, the 2040 50% cov with 10-yearly monitoring
# is less likely to still be cost-effective.

ggplot(subset(incremental_df, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_point(data= group_by(incremental_df, scenario) %>% summarise(dalys_averted=median(dalys_averted),
                                                                    total_cost = median(total_cost)),
             aes(x=dalys_averted, y = total_cost, group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), size = 5) +
  scale_fill_manual(values = rev(brewer.pal(7,"RdYlBu"))) +
  scale_colour_manual("Scenarios",
                      labels = c("screen_2020_monit_0" = "2020 screen, no monitoring",
                                 "screen_2020_monit_10" = "2020 screen, monitor every 10 years",
                                 "monit_0_screen_10b_2030" = "2020+2030 screen, no monitoring",
                                 "monit_0_screen_10b_2040" = "2020+2030+2040 screen, no monitoring",
                                 "monit_10_screen_10b_2030" = "2020+2030 screen, monitor every 10 years",
                                 "monit_10_screen_10b_2040" = "2020+2030+2040 screen, monitor every 10 years"),
                      values = c("black", rev(brewer.pal(7,"RdYlBu")))) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  geom_abline(slope=391, linetype = "dashed") +
  geom_abline(slope=518, linetype = "dashed") +
  theme_bw()

# Most cost-effective startegy <WTP again was 2040 repeat screen+10-yearly monitoring.
# With 90% coverage, 2030 repeat screen and 10-yearly monitoring were almost identical in
# cost and effect, but with 50% coverage the repeat screen has more cost and effect
# because there are not as many people to monitor in the cohort. It is at the second stage (2040)
# that they are pretty much the same again.

# With one-time 90% coverage included, this strategy dominates 2020 screens
# with 50% coverage (with/without monitoring) and 2030 repeat screen without monitoring.
# 2040 with monitoring dominates 2040 without monitoring and 2030 with monitoring.

# For 50% coverage and no monitoring, compare monitoring vs repeat screen (with discounting) ----
annual_discounting_rate <- 0.03

assemble_discounted_interactions_for_screening_strategies <- function(scenario_object,
                                                                      discount_rate = annual_discounting_rate,
                                                                      assessment_object) {
  # Compares monitoring to out3 and everything else to status quo
  # at discount rate of 3%
  out <- left_join(
    left_join(
      left_join(discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_screened",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "hbsag_tests"),
                discount_outcome_2020_to_2100(scenario_object=assessment_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_assessed",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "clinical_assessments")),
      discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                    object_to_subtract=assessment_object,
                                    outcome="interactions",
                                    interaction_outcome="total_assessed",
                                    yearly_discount_rate=discount_rate,
                                    interaction_colname = "monitoring_assessments")),
    discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                  object_to_subtract=NULL,
                                  outcome="interactions",
                                  interaction_outcome="total_treated",
                                  yearly_discount_rate=discount_rate,
                                  interaction_colname = "treatment_initiations"))

  return(out)
}

interactions_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        assemble_discounted_interactions_for_screening_strategies(out3n, assessment_object = out3n)),
  cbind(scenario = "screen_2020_monit_0_cov50",
        assemble_discounted_interactions_for_screening_strategies(out3_cov50n,
                                                                  assessment_object = out3_cov50n)),
  cbind(scenario = "monit_0_screen_10b_2030_cov50",
        assemble_discounted_interactions_for_screening_strategies(out8b_2030_cov50n,
                                                                  assessment_object = out8b_2030_cov50n)),
  cbind(scenario = "monit_0_screen_10b_2040_cov50",
        assemble_discounted_interactions_for_screening_strategies(out8b_2040_cov50n,
                                                                  assessment_object = out8b_2040_cov50n)),
  cbind(scenario = "screen_2020_monit_10_cov50",
        assemble_discounted_interactions_for_screening_strategies(out4_cov50n,
                                                                  assessment_object = out3_cov50n)),
  cbind(scenario = "monit_10_screen_10b_2030_cov50",
        assemble_discounted_interactions_for_screening_strategies(out8b_2030_monit10_cov50n,
                                                                  assessment_object = out8b_2030_cov50n)),
  cbind(scenario = "monit_10_screen_10b_2040_cov50",
        assemble_discounted_interactions_for_screening_strategies(out8b_2040_monit10_cov50n,
                                                                  assessment_object = out8b_2040_cov50n))
)

py_on_treatment_disc <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out3n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_0_cov50",
             discount_outcome_2020_to_2100(scenario_object=out3_cov50n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030_cov50",
             discount_outcome_2020_to_2100(scenario_object=out8b_2030_cov50n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2040_cov50",
             discount_outcome_2020_to_2100(scenario_object=out8b_2040_cov50n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_10_cov50",
             discount_outcome_2020_to_2100(scenario_object=out4_cov50n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2030_cov50",
             discount_outcome_2020_to_2100(scenario_object=out8b_2030_monit10_cov50n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2040_cov50",
             discount_outcome_2020_to_2100(scenario_object=out8b_2040_monit10_cov50n,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate))
)
py_on_treatment_disc$sim <- gsub("[^0-9]", "", py_on_treatment_disc$sim)

dalys_averted_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3n,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "screen_2020_monit_0_cov50",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3_cov50n,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2030_cov50n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2040_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2040_cov50n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_10_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out4_cov50n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2030_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2030_monit10_cov50n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2040_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2040_monit10_cov50n,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate))
)
dalys_averted_disc$sim <- gsub("[^0-9]", "", dalys_averted_disc$sim)
colnames(dalys_averted_disc)[colnames(dalys_averted_disc) == "dalys"] <-
  "value"

deaths_averted_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3n,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  cbind(scenario = "screen_2020_monit_0_cov50",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3_cov50n,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2030_cov50n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2040_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2040_cov50n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_10_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out4_cov50n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2030_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2030_monit10_cov50n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_10_screen_10b_2040_cov50",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_2040_monit10_cov50n,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate))
)
deaths_averted_disc$sim <- gsub("[^0-9]", "", deaths_averted_disc$sim)
colnames(deaths_averted_disc)[colnames(deaths_averted_disc) == "cum_hbv_deaths"] <-
  "value"

# Combine into dataframe (No discounting)
incremental_df_disc <- create_incremental_plot_df(interactions_df=interactions_disc,
                                                  py_on_treatment_df=py_on_treatment_disc,
                                                  deaths_averted_df=deaths_averted_disc,
                                                  ly_saved_df = dalys_averted_disc, # replace LY by DALYs
                                                  hbsag_test_cost = 8.3,
                                                  clinical_assessment_cost = 84.4,
                                                  monitoring_assessment_cost = 40.1,
                                                  treatment_py_cost = 60,
                                                  ref_label = "No treatment")
colnames(incremental_df_disc)[colnames(incremental_df_disc)=="ly_saved"] <- "dalys_averted"

# Find dominated strategies and ICERs
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- incremental_df_disc[which(incremental_df_disc$sim==
                                                          unique(incremental_df_disc$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw()


incremental_df_non_dom_disc <- subset(incremental_df_disc, scenario %in% c("screen_2020_monit_0",
                                                                           "monit_10_screen_10b_2040_cov50"))
# Even the 2040 repeat screen with monitoring has nearly 50% probability of being dominated!

icer_list <- list()
for(i in 1:183) {
  print(i)
  icer_list[[i]] <- incremental_df_non_dom_disc[which(incremental_df_non_dom_disc$sim==
                                                        unique(incremental_df_non_dom_disc$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result

# CER plot
ggplot(subset(incremental_df_disc, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_point(data= group_by(incremental_df_disc, scenario) %>% summarise(dalys_averted=median(dalys_averted),
                                                                    total_cost = median(total_cost)),
             aes(x=dalys_averted, y = total_cost, group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), size = 5) +
  scale_fill_manual(values = rev(brewer.pal(7,"RdYlBu"))) +
  scale_colour_manual("Scenarios",
                      labels = c("screen_2020_monit_0" = "90% 2020 screen,\nno monitoring",
                                 "screen_2020_monit_0_cov50" = "50% 2020 screen,\nno monitoring",
                                 "screen_2020_monit_10_cov50" = "50% 2020 screen,\nmonitor every 10 years",
                                 "monit_0_screen_10b_2030_cov50" = "50% 2020+2030 screen,\nno monitoring",
                                 "monit_0_screen_10b_2040_cov50" = "50% 2020+2030+2040 screen,\nno monitoring",
                                 "monit_10_screen_10b_2030_cov50" = "50% 2020+2030 screen,\nmonitor every 10 years",
                                 "monit_10_screen_10b_2040_cov50" = "50% 2020+2030+2040 screen,\nmonitor every 10 years"),
                      values = c("black", rev(brewer.pal(7,"RdYlBu")))) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  geom_abline(slope=391, linetype = "dashed") +
  geom_abline(slope=518, linetype = "dashed") +
  theme_bw()


# Compare repeat screens, screening coverage and monitoring ----

# Ambitious coverage scenarios =
# 2020 no monit = out3, 10-yearly = out4, 5-yearly = out5
# 2020+2030 no monit = out8b_2030, 10-yearly = out8b_2030_monit10, 5-yearly = out8b_2030_monit5
# 2020+2030+2-4- no monit = out8b_2040, 10-yearly = out8b_2040_monit10, 5-yearly = out8b_2040_monit5

# Deaths averted
equity_deaths_averted <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                 scenario_objects = list(out3,
                                                                         out4,
                                                                         out5,
                                                                         out8b_2030,
                                                                         out8b_2030_monit10,
                                                                         out8b_2030_monit5,
                                                                         out8b_2040,
                                                                         out8b_2040_monit10,
                                                                         out8b_2040_monit5,
                                                                         out8b_2020_cov10,
                                                                         out8b_2020_cov10_monit10,
                                                                         out8b_2020_cov10_monit5,
                                                                         out8b_2030_cov10,
                                                                         out8b_2030_cov10_monit10,
                                                                         out8b_2030_cov10_monit5,
                                                                         out8b_2040_cov10,
                                                                         out8b_2040_cov10_monit10,
                                                                         out8b_2040_cov10_monit5,
                                                                         out8b_2020_cov50,
                                                                         out8b_2020_cov50_monit10,
                                                                         out8b_2020_cov50_monit5,
                                                                         out8b_2030_cov50,
                                                                         out8b_2030_cov50_monit10,
                                                                         out8b_2030_cov50_monit5,
                                                                         out8b_2040_cov50,
                                                                         out8b_2040_cov50_monit10,
                                                                         out8b_2040_cov50_monit5),
                                                 outcome_to_plot = "number_averted",
                                                 counterfactual_label = "no treatment")
equity_deaths_averted$sim <- gsub("[^0-9]", "", equity_deaths_averted$sim)
equity_deaths_averted <- filter(equity_deaths_averted, by_year == 2100 &
                                  type == "number_averted")

equity_deaths_averted$monitoring_frequency <- "No monitoring"
equity_deaths_averted$monitoring_frequency[equity_deaths_averted$scenario %in%
                                             c("screen_2020_monit_10",
                                               "monit_10_screen_10b_2030",
                                               "monit_10_screen_10b_2040",
                                               "screen_2020_monit_10_cov10",
                                               "monit_10_screen_10b_2030_cov10",
                                               "monit_10_screen_10b_2040_cov10",
                                               "screen_2020_monit_10_cov50",
                                               "monit_10_screen_10b_2030_cov50",
                                               "monit_10_screen_10b_2040_cov50")] <- "Every 10 years"
equity_deaths_averted$monitoring_frequency[equity_deaths_averted$scenario %in%
                                             c("screen_2020_monit_5",
                                               "monit_5_screen_10b_2030",
                                               "monit_5_screen_10b_2040",
                                               "screen_2020_monit_5_cov10",
                                               "monit_5_screen_10b_2030_cov10",
                                               "monit_5_screen_10b_2040_cov10",
                                               "screen_2020_monit_5_cov50",
                                               "monit_5_screen_10b_2030_cov50",
                                               "monit_5_screen_10b_2040_cov50")] <-
  "Every 5 years"

equity_deaths_averted$end_repeat_screen <- "10% 2020"
equity_deaths_averted$end_repeat_screen[equity_deaths_averted$scenario %in%
                                          c("monit_0_screen_10b_2030_cov10",
                                            "monit_10_screen_10b_2030_cov10",
                                            "monit_5_screen_10b_2030_cov10")] <- "10% 2030"
equity_deaths_averted$end_repeat_screen[equity_deaths_averted$scenario %in%
                                          c("monit_0_screen_10b_2040_cov10",
                                            "monit_10_screen_10b_2040_cov10",
                                            "monit_5_screen_10b_2040_cov10")] <- "10% 2040"
equity_deaths_averted$end_repeat_screen[equity_deaths_averted$scenario %in%
                                          c("screen_2020_monit_0",
                                            "screen_2020_monit_10",
                                            "screen_2020_monit_5")] <- "90% 2020"
equity_deaths_averted$end_repeat_screen[equity_deaths_averted$scenario %in%
                                          c("monit_0_screen_10b_2030",
                                            "monit_10_screen_10b_2030",
                                            "monit_5_screen_10b_2030")] <- "90% 2030"
equity_deaths_averted$end_repeat_screen[equity_deaths_averted$scenario %in%
                                          c("monit_0_screen_10b_2040",
                                            "monit_10_screen_10b_2040",
                                            "monit_5_screen_10b_2040")] <-
  "90% 2040"
equity_deaths_averted$end_repeat_screen[equity_deaths_averted$scenario %in%
                                          c("screen_2020_monit_0_cov50",
                                            "screen_2020_monit_10_cov50",
                                            "screen_2020_monit_5_cov50")] <- "50% 2020"
equity_deaths_averted$end_repeat_screen[equity_deaths_averted$scenario %in%
                                          c("monit_0_screen_10b_2030_cov50",
                                            "monit_10_screen_10b_2030_cov50",
                                            "monit_5_screen_10b_2030_cov50")] <- "50% 2030"
equity_deaths_averted$end_repeat_screen[equity_deaths_averted$scenario %in%
                                          c("monit_0_screen_10b_2040_cov50",
                                            "monit_10_screen_10b_2040_cov50",
                                            "monit_5_screen_10b_2040_cov50")] <-
  "50% 2040"

# LY saved (cohort) only for 90% coverage so far
equity_ly_saved1 <- plot_ly_gained_cohort(counterfactual_object = out1,
                                          scenario_objects = list(out3,
                                                                  out4,
                                                                  out5),
                                          outcome_to_plot = "number_averted",
                                          counterfactual_label = "no treatment")

equity_ly_saved2 <- plot_ly_gained_cohort(counterfactual_object = out8b_2030_sq,
                                          scenario_objects = list(out8b_2030,
                                                                  out8b_2030_monit10,
                                                                  out8b_2030_monit5),
                                          outcome_to_plot = "number_averted",
                                          counterfactual_label = "no treatment")

equity_ly_saved3 <- plot_ly_gained_cohort(counterfactual_object = out8b_2040_sq,
                                          scenario_objects = list(out8b_2040,
                                                                  out8b_2040_monit10,
                                                                  out8b_2040_monit5),
                                          outcome_to_plot = "number_averted",
                                          counterfactual_label = "no treatment")
equity_ly_saved <- rbind(equity_ly_saved1, equity_ly_saved2, equity_ly_saved3)
colnames(equity_ly_saved)[colnames(equity_ly_saved) %in% c("counterfactual", "scenario")] <-
  c("scenario", "counterfactual")
equity_ly_saved$sim <- gsub("[^0-9]", "", equity_ly_saved$sim)
equity_ly_saved <- filter(equity_ly_saved, type == "number_averted")

equity_ly_saved$monitoring_frequency <- "No monitoring"
equity_ly_saved$monitoring_frequency[equity_ly_saved$scenario %in%
                                       c("screen_2020_monit_10",
                                         "monit_10_screen_10b_2030",
                                         "monit_10_screen_10b_2040"#,
                                         #"screen_2020_monit_10_cov10",
                                         #"monit_10_screen_10b_2030_cov10",
                                         #"monit_10_screen_10b_2040_cov10"
                                       )] <- "Every 10 years"
equity_ly_saved$monitoring_frequency[equity_ly_saved$scenario %in%
                                       c("screen_2020_monit_5",
                                         "monit_5_screen_10b_2030",
                                         "monit_5_screen_10b_2040"#,
                                         #"screen_2020_monit_5_cov10",
                                         #"monit_5_screen_10b_2030_cov10",
                                         #"monit_5_screen_10b_2040_cov10"
                                       )] <- "Every 5 years"

equity_ly_saved$end_repeat_screen <- NA
equity_ly_saved$end_repeat_screen[equity_ly_saved$scenario %in%
                                    c("screen_2020_monit_0",
                                      "screen_2020_monit_10",
                                      "screen_2020_monit_5")] <- "90% 2020"
equity_ly_saved$end_repeat_screen[equity_ly_saved$scenario %in%
                                    c("monit_0_screen_10b_2030",
                                      "monit_10_screen_10b_2030",
                                      "monit_5_screen_10b_2030")] <- "90% 2030"
equity_ly_saved$end_repeat_screen[equity_ly_saved$scenario %in%
                                    c("monit_0_screen_10b_2040",
                                      "monit_10_screen_10b_2040",
                                      "monit_5_screen_10b_2040")] <-
  "90% 2040"

# Barchart of deaths averted - all very similar
# Increase going from 2020-no monitoring to 2030-10 yearly monitoring is larger than
# to 2040-no monitoring

View(equity_deaths_averted %>%
       group_by(end_repeat_screen, monitoring_frequency) %>%
       summarise(median(value)))

ggplot(equity_deaths_averted,
       aes(x=end_repeat_screen, y = value,
           fill = reorder(monitoring_frequency, -value))) +
  stat_summary(fun= median,position=position_dodge(width=0.7),
               geom="bar", width = 0.7)+
  geom_segment(aes(x=2.7, y=1368, xend=10, yend=1368), lty="dashed") +
  geom_segment(aes(x=3.7, y=3067, xend=10, yend=3067), lty="dashed") +
  geom_segment(aes(x=6.1, y=3803, xend=10, yend=3803), lty="dashed") +
  geom_segment(aes(x=6.9, y=5027, xend=10, yend=5027), lty="dashed") +
  geom_vline(xintercept = 3.5) +
  geom_vline(xintercept = 6.5) +
  ylab("HBV deaths averted") +
  xlab("Screening coverage & end year of repeat screen") +
  labs(fill="Monitoring frequency") +
  #  stat_summary(fun.min = function(x)quantile(x,0.025),
  #               fun.max = function(x) quantile(x,0.975),
  #               position=position_dodge(width=0.95),
  #               geom="errorbar", width = 0.5)+
  #  ylim(0,15000) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust = 1))
#  scale_colour_manual(values = c("black", "blue", "red"))
# On the median, shows that all simulations with 10% coverage get nowhere near the 90% effect of any scenario.
# For 90% coverage: 2020 screen with 10/5 yearly monitoring has about the same effect as
# 2040 repeat screen without monitoring.
# 2030 repeat screen with 5 yearly monitoring has nearly the same effect as 2040 repeat screen with 5-yearly
# monitoring. Note errorbars are completely overlapping for all
# 50% coverage 2020-2040 is about the same as 90% coverage in 2020! Is this due to declining
# carrier numbers over time (but then if infections go down so do deaths)? Or rather more HBV deaths to
# be averted in absolute terms? Check deaths averted per (untreated) carrier!

# Impact of adding 10-yearly monitoring vs repeating screen in 2030 is the same (actually slightly better).
# With repeat screen in 2030, impact of adding 10 yearly monitoring is far more than
# adding a second repeat screen. This is the case for 90% coverage and also for 50% coverage
# if monitoring in 2020 is 5-yearly (for 50% coverage, better to add 1 repeat screen than
# 10 yearly monitoring)

# From lower coverage also need to consider increasing coverage as an example of expanding access.
# Both increasing to 90% coverage and to 2 repeat screens have higher impact than one-off screen
# with 5-yearly monitoring.

# Next question is: what is the maximum cost-effective strategy here?

# SLOPE CHART TEST
#equity_deaths_averted2 <- separate(data = equity_deaths_averted,
#                                   col = end_repeat_screen, into = c("cov", "end_year"), sep = " ")
#equity_deaths_averted2$cov_monit <- paste(equity_deaths_averted2$cov, equity_deaths_averted2$monitoring_frequency)
#equity_deaths_averted2 <- select(equity_deaths_averted2, sim, end_year, cov_monit, value) #%>%
#  pivot_wider(names_from = end_year,
#              values_from = value) %>%
#  group_by(cov_monit) %>%
#  summarise(median_2020 = median(`2020`),
#            median_2030 = median(`2030`),
#            median_2040 = median(`2040`))
#equity_deaths_averted2_long <- pivot_longer(equity_deaths_averted2,
#                                            !cov_monit, names_to = "end_year", values_to = "median")#
#
#equity_deaths_averted2_long$end_year[equity_deaths_averted2_long$end_year=="median_2020"] <- 1
#equity_deaths_averted2_long$end_year[equity_deaths_averted2_long$end_year=="median_2030"] <- 2
#equity_deaths_averted2_long$end_year[equity_deaths_averted2_long$end_year=="median_2040"] <- 3##

#test <- filter(equity_deaths_averted2_long, cov_monit %in%
#                 c("90% Every 10 years", "90% Every 5 years", "90% No monitoring")) %>%
#  select(end_year, median, cov_monit)
#comb <- combn(nrow(test), 2)
#connections <- data.frame(
#  from = test[comb[1, ], 1:2],
#  to   = test[comb[2, ], 1:3]
#)
#names(connections) <- c("x1", "y1", "x2", "y2", "label")
#ggplot(test) +
#  geom_segment(data=connections, aes(x=x1, y=y1, xend=x2, yend=y2), col="grey") +
#  geom_point(aes(end_year, median, col = cov_monit), size=5) +
#  scale_x_discrete(labels = c("1" = "2020",
#                              "2" = "2020-2030",
#                              "3" = "2020-2040")) +
#  ylim(0,8200)


# Barchart of LY saved for 90% coverage
ggplot(equity_ly_saved,
       aes(x=end_repeat_screen, y = value,
           fill = reorder(monitoring_frequency, -value))) +
  stat_summary(fun= median,position=position_dodge(width=0.95),
               geom="bar")+
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               position=position_dodge(width=0.95),
               geom="errorbar", width = 0.5, alpha = 0.25)+
  #  ylim(0,17500) +
  theme_bw()

# With cost-effectiveness (test)
equity_deaths_df <- equity_deaths_averted %>%
  filter(end_repeat_screen %in% c("90% 2020", "90% 2030", "90% 2040",
                                  "50% 2020", "50% 2030", "50% 2040")) %>%
  select(scenario, sim, value) %>%
  pivot_wider(names_from = "scenario", values_from = "value") %>%
  select(-sim)
intervention_labels <- colnames(equity_deaths_df)
equity_deaths_df <- cbind(as.matrix(equity_deaths_df), rep(0,183))

# ADD 50% COVERAGE HERE

equity_int_df <- cbind(unlist(out3$interactions[[16]]$total_interactions[,-c(1:3)]),
                       unlist(out4$interactions[[16]]$total_interactions[,-c(1:3)]),
                       unlist(out5$interactions[[16]]$total_interactions[,-c(1:3)]),
                       unlist(out8b_2030$interactions[[16]]$total_interactions[,-c(1:3)]),
                       unlist(out8b_2030_monit10$interactions[[16]]$total_interactions[,-c(1:3)]),
                       unlist(out8b_2030_monit5$interactions[[16]]$total_interactions[,-c(1:3)]),
                       unlist(out8b_2040$interactions[[16]]$total_interactions[,-c(1:3)]),
                       unlist(out8b_2040_monit10$interactions[[16]]$total_interactions[,-c(1:3)]),
                       unlist(out8b_2040_monit5$interactions[[16]]$total_interactions[,-c(1:3)]),
                       rep(0,183))
equity_int_df <- as.matrix(equity_int_df)

ceef.plot.median(bcea(e=equity_deaths_df,c=equity_int_df,ref=ncol(equity_deaths_df),
                      interventions=c(intervention_labels, "Status quo")),
                 graph="base")

# Assuming 22 years on treatment on average
cs <- 8.3 # cost of screen per person
ca <- 72.3 # cost of assessment per person
cm <- 43.1 # cost of monitoring per person
ct <- 60*22 # cost per person for 22 years on treatment

equity_cost_df <- cbind((unlist(out3$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out3$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           unlist(out3$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        (unlist(out4$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out3$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           (unlist(out4$interactions[[16]]$total_assessed[,-c(1:3)])-
                              unlist(out3$interactions[[16]]$total_assessed[,-c(1:3)]))*cm+
                           unlist(out4$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        (unlist(out5$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out3$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           (unlist(out5$interactions[[16]]$total_assessed[,-c(1:3)])-
                              unlist(out3$interactions[[16]]$total_assessed[,-c(1:3)]))*cm+
                           unlist(out5$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        (unlist(out8b_2030$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out8b_2030$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           unlist(out8b_2030$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        (unlist(out8b_2030_monit10$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out8b_2030$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           (unlist(out8b_2030_monit10$interactions[[16]]$total_assessed[,-c(1:3)])-
                              unlist(out8b_2030$interactions[[16]]$total_assessed[,-c(1:3)]))*cm+
                           unlist(out8b_2030_monit10$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        (unlist(out8b_2030_monit5$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out8b_2030$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           (unlist(out8b_2030_monit5$interactions[[16]]$total_assessed[,-c(1:3)])-
                              unlist(out8b_2030$interactions[[16]]$total_assessed[,-c(1:3)]))*cm+
                           unlist(out8b_2030_monit5$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        (unlist(out8b_2040$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out8b_2040$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           unlist(out8b_2040$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        (unlist(out8b_2040_monit10$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out8b_2040$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           (unlist(out8b_2040_monit10$interactions[[16]]$total_assessed[,-c(1:3)])-
                              unlist(out8b_2040$interactions[[16]]$total_assessed[,-c(1:3)]))*cm+
                           unlist(out8b_2040_monit10$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        (unlist(out8b_2040_monit5$interactions[[16]]$total_screened[,-c(1:3)])*cs+
                           unlist(out8b_2040$interactions[[16]]$total_assessed[,-c(1:3)])*ca+
                           (unlist(out8b_2040_monit5$interactions[[16]]$total_assessed[,-c(1:3)])-
                              unlist(out8b_2040$interactions[[16]]$total_assessed[,-c(1:3)]))*cm+
                           unlist(out8b_2040_monit5$interactions[[16]]$total_treated[,-c(1:3)])*ct),
                        rep(0,183))
equity_cost_df <- as.matrix(equity_cost_df)

ceef.plot.median(bcea(e=equity_deaths_df[,c(1:9,ncol(equity_deaths_df))],
                      c=equity_cost_df,ref=10,
                      interventions=c(intervention_labels[1:9], "Status quo")),
                 graph="base")

equity_ly_df <- equity_ly_saved %>%
  select(scenario, sim, value) %>%
  pivot_wider(names_from = "scenario", values_from = "value") %>%
  select(-sim)
intervention_labels <- colnames(equity_ly_df)
equity_ly_df <- cbind(as.matrix(equity_ly_df), rep(0,183))

ceef.plot.median(bcea(e=equity_ly_df,c=equity_cost_df,ref=ncol(equity_ly_df),
                      interventions=c(intervention_labels, "Status quo")),
                 graph="base")
# For the Fibroscan cost, on frontier are (for 90% only): 2020 no monit, 2040 monit 10,
# 2040 monit 5, yet others are very close so need to do cost-effectiveness on each sim separately.
# In terms of the median ICER, even 2040 5 years would be considered cost-effective!
# But similar average value to looking at individual sims.
# Checked in the monitoring freq. analysis and even there 5-yearly monitoring has a
# very low ICER compared to SQ! it is just the incremental to every 6 years is bad.

# Testplot: Which strategy has highest deaths averted per cost?
deaths_averted_per_cost <- data.frame(equity_deaths_df[,c(1:9)]/equity_cost_df[,1:9])
deaths_averted_per_cost <- pivot_longer(deaths_averted_per_cost, everything(),
                                        names_to = "scenario", values_to = "value")

deaths_averted_per_cost$monitoring_frequency <- "No monitoring"
deaths_averted_per_cost$monitoring_frequency[deaths_averted_per_cost$scenario %in%
                                               c("screen_2020_monit_10",
                                                 "monit_10_screen_10b_2030",
                                                 "monit_10_screen_10b_2040"#,
                                                 #"screen_2020_monit_10_cov10",
                                                 #"monit_10_screen_10b_2030_cov10",
                                                 #"monit_10_screen_10b_2040_cov10"
                                               )] <- "Every 10 years"
deaths_averted_per_cost$monitoring_frequency[deaths_averted_per_cost$scenario %in%
                                               c("screen_2020_monit_5",
                                                 "monit_5_screen_10b_2030",
                                                 "monit_5_screen_10b_2040"#,
                                                 #"screen_2020_monit_5_cov10",
                                                 #"monit_5_screen_10b_2030_cov10",
                                                 #"monit_5_screen_10b_2040_cov10"
                                               )] <- "Every 5 years"

deaths_averted_per_cost$end_repeat_screen <- NA
deaths_averted_per_cost$end_repeat_screen[deaths_averted_per_cost$scenario %in%
                                            c("screen_2020_monit_0",
                                              "screen_2020_monit_10",
                                              "screen_2020_monit_5")] <- "90% 2020"
deaths_averted_per_cost$end_repeat_screen[deaths_averted_per_cost$scenario %in%
                                            c("monit_0_screen_10b_2030",
                                              "monit_10_screen_10b_2030",
                                              "monit_5_screen_10b_2030")] <- "90% 2030"
deaths_averted_per_cost$end_repeat_screen[deaths_averted_per_cost$scenario %in%
                                            c("monit_0_screen_10b_2040",
                                              "monit_10_screen_10b_2040",
                                              "monit_5_screen_10b_2040")] <-
  "90% 2040"

ggplot(deaths_averted_per_cost,
       aes(x=end_repeat_screen, y = value,
           fill = reorder(monitoring_frequency, -value))) +
  stat_summary(fun= median,position=position_dodge(width=0.95),
               geom="bar")+
  stat_summary(fun.min = function(x)quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               position=position_dodge(width=0.95),
               geom="errorbar", width = 0.5)+
  #   ylim(0,15000) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust = 1))

## Compare regular antenatal care screening (76% antenatal coverage) with the population-based screen ----
# Note that screening in ANC only requires the cost of the test.
# Repeat screening simulations here don't include IT treatment yet so need to compare this
# to the one-off screen from monitoring_frequency:

# Files:
# out3_it: Standard pop-based screen in 2020, no monitoring
# out5_it: Standard pop-based screen in 2020, 5-yearly monitoring
# e1_it_out3_gradual: same pop-based screen but with 10% coverage each year spread over
# 2020-2028, no monitoring
# e1_it_out5_gradual: same pop-based screen but with 10% coverage each year spread over
# 2020-2028, 5-yearly monitoring
# anc_it_2020_monit_0: screen 90% of pregnant women in 2020, no monitoring
# anc_it_2030_monit_0b: screen 90% of pregnant women all the time until 2030, no monitoring,
# no retesting of those previosuly found to be HBsAg-negative
# anc_it_2030_monit_5a: : screen 90% of pregnant women all the time until 2030, 5-yearly monitoring,
# with retesting of those previosuly found to be HBsAg-negative
# anc_it_2030_monit_0a: screen 90% of pregnant women all the time until 2030, no monitoring,
# with retesting of those previosuly found to be HBsAg-negative
# This is basically testing all women presenting to antenatal care any time until 2030
# anc_it_2040_monit_5a: screen 90% of pregnant women all the time until 2040, 5-yearly monitoring,
# with retesting of those previosuly found to be HBsAg-negative

# Scenario names + object
# screen_2020_monit_0 - out3_it
# screen_2020_monit_5 - out5_it
# screen_2020_anc_monit_0 - anc_it_2020_monit_0
# anc_2030_no_rescreen_monit_0 - anc_it_2030_monit_0b
# anc_2030_with_rescreen_monit_5 - anc_it_2030_monit_5a
# anc_2030_with_rescreen_monit_0 - anc_it_2030_monit_0a
# anc_2040_with_rescreen_monit_5 - anc_it_2040_monit_5a

annual_discounting_rate <- 0.03

assemble_discounted_interactions_for_screening_strategies <- function(scenario_object,
                                                                      discount_rate = annual_discounting_rate,
                                                                      assessment_object) {
  # Compares monitoring to out3 and everything else to status quo
  # at discount rate of 3%
  out <- left_join(
    left_join(
      left_join(discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_screened",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "hbsag_tests"),
                discount_outcome_2020_to_2100(scenario_object=assessment_object,
                                              object_to_subtract=NULL,
                                              outcome="interactions",
                                              interaction_outcome="total_assessed",
                                              yearly_discount_rate=discount_rate,
                                              interaction_colname = "clinical_assessments")),
      discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                    object_to_subtract=assessment_object,
                                    outcome="interactions",
                                    interaction_outcome="total_assessed",
                                    yearly_discount_rate=discount_rate,
                                    interaction_colname = "monitoring_assessments")),
    discount_outcome_2020_to_2100(scenario_object=scenario_object,
                                  object_to_subtract=NULL,
                                  outcome="interactions",
                                  interaction_outcome="total_treated",
                                  yearly_discount_rate=discount_rate,
                                  interaction_colname = "treatment_initiations"))

  return(out)
}

# Here the assessment object only has to change for simulations that involve monitoring
anc_interactions_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        assemble_discounted_interactions_for_screening_strategies(out3_it,
                                                                  assessment_object = out3_it)),
  cbind(scenario = "screen_2020_monit_5",
        assemble_discounted_interactions_for_screening_strategies(out5_it,
                                                                  assessment_object = out3_it)),
  cbind(scenario = "screen_2020_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(monit_out7,
                                                                  assessment_object = out3_it)),
  cbind(scenario = "screen_2020_anc_monit_0",
        assemble_discounted_interactions_for_screening_strategies(anc_it_2020_monit_0,
                                                                  assessment_object = anc_it_2020_monit_0)),
  cbind(scenario = "anc_2030_no_rescreen_monit_0",
        assemble_discounted_interactions_for_screening_strategies(anc_it_2030_monit_0b,
                                                                  assessment_object = anc_it_2030_monit_0b)),
  cbind(scenario = "anc_2030_no_rescreen_monit_5",
        assemble_discounted_interactions_for_screening_strategies(anc_it_2030_monit_5b,
                                                                  assessment_object = anc_it_2030_monit_0b)),
  cbind(scenario = "anc_2030_with_rescreen_monit_0",
        assemble_discounted_interactions_for_screening_strategies(anc_it_2030_monit_0a,
                                                                  assessment_object = anc_it_2030_monit_0a)),
  cbind(scenario = "anc_2040_no_rescreen_monit_5",
        assemble_discounted_interactions_for_screening_strategies(anc_it_2040_monit_5b,
                                                                  assessment_object = anc_it_2040_monit_5b)),
  cbind(scenario = "anc_2040_no_rescreen_monit_0",
        assemble_discounted_interactions_for_screening_strategies(anc_it_2040_monit_0b,
                                                                  assessment_object = anc_it_2040_monit_0b)),
  cbind(scenario = "anc_2040_with_rescreen_monit_0",
        assemble_discounted_interactions_for_screening_strategies(anc_it_2040_monit_0a,
                                                                  assessment_object = anc_it_2040_monit_0a)),
# Population-based repeat screening:
  cbind(scenario = "monit_0_screen_10a_2030",
        assemble_discounted_interactions_for_screening_strategies(out8a_it_2030_monit_0,
                                                                  assessment_object = out8a_it_2030_monit_0)),
  cbind(scenario = "monit_5_screen_10a_2030",
        assemble_discounted_interactions_for_screening_strategies(out8a_it_2030_monit_5,
                                                                  assessment_object = out8a_it_2030_monit_0)),
  cbind(scenario = "monit_0_screen_10b_2030",
        assemble_discounted_interactions_for_screening_strategies(out8b_it_2030_monit_0,
                                                                  assessment_object = out8b_it_2030_monit_0)),
  cbind(scenario = "monit_5_screen_10b_2030",
        assemble_discounted_interactions_for_screening_strategies(out8b_it_2030_monit_5,
                                                                  assessment_object = out8b_it_2030_monit_0)),
  # Combination of pop and ANC testing:
  cbind(scenario = "pop_2020_anc_2030_no_rescreen_monit_0",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2030_monit_0b,
                                                                  assessment_object = pop_2020_anc_2030_monit_0b)),
  cbind(scenario = "pop_2020_anc_2030_no_rescreen_monit_5",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2030_monit_5b,
                                                                  assessment_object = pop_2020_anc_2030_monit_0b)),
  cbind(scenario = "pop_2020_anc_2030_no_rescreen_monit_sim7",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2030_monit_out7b,
                                                                  assessment_object = pop_2020_anc_2030_monit_0b)),
  cbind(scenario = "pop_2020_anc_2030_with_rescreen_monit_0",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2030_monit_0a,
                                                                  assessment_object = pop_2020_anc_2030_monit_0a)),
  cbind(scenario = "pop_2020_anc_2040_no_rescreen_monit_0",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2040_monit_0b,
                                                                  assessment_object = pop_2020_anc_2040_monit_0b)),
cbind(scenario = "pop_2020_anc_2040_no_rescreen_monit_5",
      assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2040_monit_5b,
                                                                assessment_object = pop_2020_anc_2040_monit_0b)),
cbind(scenario = "pop_2020_anc_2040_no_rescreen_monit_sim7",
      assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2040_monit_out7b,
                                                                assessment_object = pop_2020_anc_2040_monit_0b)),
  cbind(scenario = "pop_2020_anc_2040_with_rescreen_monit_0",
        assemble_discounted_interactions_for_screening_strategies(pop_2020_anc_2040_monit_0a,
                                                                  assessment_object = pop_2020_anc_2040_monit_0a))
  )


anc_py_on_treatment_disc <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out3_it,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out5_it,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=monit_out7,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_anc_monit_0",
             discount_outcome_2020_to_2100(scenario_object=anc_it_2020_monit_0,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=anc_it_2030_monit_0b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=anc_it_2030_monit_5b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=anc_it_2030_monit_0a,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=anc_it_2040_monit_5b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=anc_it_2040_monit_0b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=anc_it_2040_monit_0a,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  # Population-based repeat screening:
  data.frame(scenario = "monit_0_screen_10a_2030",
             discount_outcome_2020_to_2100(scenario_object=out8a_it_2030_monit_0,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_5_screen_10a_2030",
             discount_outcome_2020_to_2100(scenario_object=out8a_it_2030_monit_5,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out8b_it_2030_monit_0,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_5_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out8b_it_2030_monit_5,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  # Combination of pop and ANC testing:
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=pop_2020_anc_2030_monit_0b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=pop_2020_anc_2030_monit_5b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=pop_2020_anc_2030_monit_out7b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=pop_2020_anc_2030_monit_0a,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=pop_2020_anc_2040_monit_0b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=pop_2020_anc_2040_monit_5b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=pop_2020_anc_2040_monit_out7b,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=pop_2020_anc_2040_monit_0a,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate))
  )
anc_py_on_treatment_disc$sim <- gsub("[^0-9]", "", anc_py_on_treatment_disc$sim)


# DALYS averted
anc_dalys_averted_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3_it,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out5_it,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=monit_out7,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_anc_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2020_monit_0,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2030_monit_0b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2030_monit_5b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2030_monit_0a,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2040_monit_5b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2040_monit_0b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2040_monit_0a,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  # Population-based repeat screening:
  data.frame(scenario = "monit_0_screen_10a_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8a_it_2030_monit_0,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_5_screen_10a_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8a_it_2030_monit_5,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_it_2030_monit_0,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_5_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_it_2030_monit_5,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  # Combination of pop and ANC testing:
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2030_monit_0b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2030_monit_5b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2030_monit_out7b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2030_monit_0a,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2040_monit_0b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2040_monit_5b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2040_monit_out7b,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2040_monit_0a,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate))
)
anc_dalys_averted_disc$sim <- gsub("[^0-9]", "", anc_dalys_averted_disc$sim)
colnames(anc_dalys_averted_disc)[colnames(anc_dalys_averted_disc) == "dalys"] <-
  "value"

# Deaths averted
anc_deaths_averted_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3_it,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out5_it,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract= monit_out7,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_anc_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2020_monit_0,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2030_monit_0b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2030_monit_5b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2030_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2030_monit_0a,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2040_monit_5b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2040_monit_0b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "anc_2040_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=anc_it_2040_monit_0a,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  # Population-based repeat screening:
  data.frame(scenario = "monit_0_screen_10a_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8a_it_2030_monit_0,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_5_screen_10a_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8a_it_2030_monit_5,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_0_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_it_2030_monit_0,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "monit_5_screen_10b_2030",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out8b_it_2030_monit_5,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  # Combination of pop and ANC testing:
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2030_monit_0b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2030_monit_5b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_no_rescreen_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2030_monit_out7b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2030_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2030_monit_0a,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2040_monit_0b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2040_monit_5b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_no_rescreen_monit_sim7",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2040_monit_out7b,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "pop_2020_anc_2040_with_rescreen_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=pop_2020_anc_2040_monit_0a,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate))
)
anc_deaths_averted_disc$sim <- gsub("[^0-9]", "", anc_deaths_averted_disc$sim)
colnames(anc_deaths_averted_disc)[colnames(anc_deaths_averted_disc) == "cum_hbv_deaths"] <-
  "value"

# Combine into dataframe with updated costs!
hbsag_test_cost <- 8.3  ## ensure hbsag_test_cost is the same as below
anc_incremental_df_disc <- create_incremental_plot_df(interactions_df=anc_interactions_disc,
                                                  py_on_treatment_df=anc_py_on_treatment_disc,
                                                  deaths_averted_df=anc_deaths_averted_disc,
                                                  ly_saved_df = anc_dalys_averted_disc, # replace LY by DALYs
                                                  hbsag_test_cost = 8.3,
                                                  clinical_assessment_cost = 33, # 84.4,
                                                  monitoring_assessment_cost = 25.5, # 40.1,
                                                  treatment_py_cost = 66.5, #60,
                                                  ref_label = "No treatment")
colnames(anc_incremental_df_disc)[colnames(anc_incremental_df_disc)=="ly_saved"] <- "dalys_averted"

# Correct the cost of HBsAg test in ANC to 1.7 per HBsAg test:
# Strategies with only ANC testing:
anc_only_strategies <- c("screen_2020_anc_monit_0", "anc_2030_no_rescreen_monit_0",
                         "anc_2030_no_rescreen_monit_5", "anc_2030_with_rescreen_monit_0",
                         "anc_2040_no_rescreen_monit_5", "anc_2040_no_rescreen_monit_0",
                         "anc_2040_with_rescreen_monit_0")
anc_incremental_df_disc[
  anc_incremental_df_disc$scenario %in% anc_only_strategies,]$screening_cost <-
  anc_incremental_df_disc[
    anc_incremental_df_disc$scenario %in% anc_only_strategies,]$hbsag_tests*1.7
anc_incremental_df_disc[
  anc_incremental_df_disc$scenario %in% anc_only_strategies,]$total_cost <-
  anc_incremental_df_disc[
    anc_incremental_df_disc$scenario %in% anc_only_strategies,]$screening_cost +
  anc_incremental_df_disc[
    anc_incremental_df_disc$scenario %in% anc_only_strategies,]$assessment_cost +
  anc_incremental_df_disc[
    anc_incremental_df_disc$scenario %in% anc_only_strategies,]$monitoring_cost+
  anc_incremental_df_disc[
    anc_incremental_df_disc$scenario %in% anc_only_strategies,]$treatment_cost
# Strategies combining population and ANC testing
# Need to distinguish between the screens in the population test and those conducted in ANC
# (compare to screening cost in one-time pop test)
# Normal cost from the pop-based test + reduced cost for subsequent ANC screens
pop_anc_combi_strategies <- c("pop_2020_anc_2030_no_rescreen_monit_0",
                              "pop_2020_anc_2030_no_rescreen_monit_5",
                              "pop_2020_anc_2030_no_rescreen_monit_sim7",
                              "pop_2020_anc_2030_with_rescreen_monit_0",
                              "pop_2020_anc_2040_no_rescreen_monit_0",
                              "pop_2020_anc_2040_with_rescreen_monit_0")

for (i in 1:length(pop_anc_combi_strategies)) {
  anc_incremental_df_disc[
    anc_incremental_df_disc$scenario == pop_anc_combi_strategies[i],]$screening_cost <-
    anc_incremental_df_disc[anc_incremental_df_disc$scenario ==
                              "screen_2020_monit_0",]$hbsag_tests*hbsag_test_cost +
    (anc_incremental_df_disc[
      anc_incremental_df_disc$scenario == pop_anc_combi_strategies[i],]$hbsag_tests-
       anc_incremental_df_disc[
         anc_incremental_df_disc$scenario == "screen_2020_monit_0",]$hbsag_tests)*1.7
  anc_incremental_df_disc[
    anc_incremental_df_disc$scenario == pop_anc_combi_strategies[i],]$total_cost <-
    anc_incremental_df_disc[
      anc_incremental_df_disc$scenario == pop_anc_combi_strategies[i],]$screening_cost +
    anc_incremental_df_disc[
      anc_incremental_df_disc$scenario == pop_anc_combi_strategies[i],]$assessment_cost +
    anc_incremental_df_disc[
      anc_incremental_df_disc$scenario == pop_anc_combi_strategies[i],]$monitoring_cost+
    anc_incremental_df_disc[
      anc_incremental_df_disc$scenario == pop_anc_combi_strategies[i],]$treatment_cost
}


# Cost-effectiveness analysis ----

# Find dominated strategies and ICERs
dominance_prob_list <- list()

for(i in 1:183) {
  print(i)
  dominance_prob_list[[i]] <- anc_incremental_df_disc[which(anc_incremental_df_disc$sim==
                                                          unique(anc_incremental_df_disc$sim)[i]),]
  dominance_prob_list[[i]] <- assign_dominated_strategies(dominance_prob_list[[i]],
                                                          exposure="total_cost",
                                                          outcome="dalys_averted")
}
dominance_prob_df <- do.call("rbind", dominance_prob_list)
dominance_prob_result <- group_by(dominance_prob_df, scenario, dominated) %>%
  tally() %>%
  spread(key = "dominated", value = "n") %>%
  replace_na(list(No = 0, Yes = 0))
dominance_prob_result$prob_non_dominated <- dominance_prob_result$No/
  (dominance_prob_result$Yes+dominance_prob_result$No)
# Any less than 50% chance of being non-dominated?
any(dominance_prob_result$prob_non_dominated<0.5)

ggplot(dominance_prob_result) +
  geom_col(aes(x=reorder(scenario, desc(prob_non_dominated)), y = prob_non_dominated)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1))

# With 3% discounting rate and 76% antenatal coverage (same as with 90% cov),
# and with pop/ANC combination strategies added in + updated costs from Liem:
# The non-dominated strategies are:
# monit_5_screen_10a_2030, monit_5_screen_10b_2030, screen_2020_anc_monit_0,
# anc_2030_no_rescreen_monit_0, anc_2040_no_rescreen_monit_0,
# pop_2020_anc_2040_no_rescreen_monit_0,
# pop_2020_anc_2030_no_rescreen_monit_sim7

# But some of these are very close to 50%!
# Still need to add in optimal monitoring frequencies for ANC and pop repeat screening
# strategies.

anc_incremental_df_non_dom_disc <- subset(anc_incremental_df_disc,
                                          scenario %in% c("monit_5_screen_10a_2030",
                                                          "monit_5_screen_10b_2030",
                                                          "screen_2020_anc_monit_0",
                                                          "anc_2030_no_rescreen_monit_0",
                                                          "anc_2040_no_rescreen_monit_0",
                                                          "pop_2020_anc_2030_no_rescreen_monit_sim7",
                                                          "pop_2020_anc_2040_no_rescreen_monit_0"))

icer_list <- list()
for(i in 1:183) {
  print(i)
  icer_list[[i]] <- anc_incremental_df_non_dom_disc[which(anc_incremental_df_non_dom_disc$sim==
                                                        unique(anc_incremental_df_non_dom_disc$sim)[i]),]
  icer_list[[i]] <- calculate_icer_per_sim(icer_list[[i]],
                                           exposure="total_cost",
                                           outcome="dalys_averted")
}
icer_df <- do.call("rbind", icer_list)
icer_result <- group_by(icer_df, scenario, comparator) %>%
  arrange(sim,total_cost) %>%
  summarise(icer_median = median(icer),
            icer_lower = quantile(icer, 0.025),
            icer_upper = quantile(icer, 0.975)) %>%
  arrange(icer_median)
icer_result

anc_incremental_df_disc$frontier <- "Dominated"
anc_incremental_df_disc$frontier[anc_incremental_df_disc$scenario %in%
                                   c("monit_5_screen_10a_2030",
                                     "monit_5_screen_10b_2030",
                                     "screen_2020_anc_monit_0",
                                     "anc_2030_no_rescreen_monit_0",
                                     "anc_2040_no_rescreen_monit_0",
                                     "pop_2020_anc_2030_no_rescreen_monit_sim7",
                                     "pop_2020_anc_2040_no_rescreen_monit_0",
                                     "No treatment")] <- "Non-dominated"

anc_incremental_df_disc_median <- group_by(anc_incremental_df_disc, scenario, frontier) %>%
  summarise(deaths_averted = median(deaths_averted),
            total_interactions = median(total_interactions),
            dalys_averted = median(dalys_averted),
            total_cost = median(total_cost))

# CER plot
ggplot(subset(anc_incremental_df_disc, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = subset(anc_incremental_df_disc_median, frontier == "Non-dominated"),
            aes(x = dalys_averted, y = total_cost), size = 1) +
  geom_point(data = anc_incremental_df_disc_median,
             aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost),
                 shape = frontier),
             size = 5) +
  geom_point(data = subset(anc_incremental_df_disc_median, frontier == "Dominated"),
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5, shape = 1,
             colour = "black") +
  geom_point(data = subset(anc_incremental_df_disc_median, frontier == "Non-dominated"),
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5, shape = 2,
             colour = "black") +
  scale_fill_manual(values = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 11))) +
  scale_colour_manual("Screening strategy\n(monitoring frequency,re-testing)",
                      values = c("black",  brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 11)),
                      labels = c("anc_2030_no_rescreen_monit_0"="ANC 2020-2030 (0,no)",
                                 "anc_2030_with_rescreen_monit_0"="ANC 2020-2030 (0,yes)",
                                 "anc_2030_no_rescreen_monit_5"="ANC 2020-2030 (5,no)",
                                 "anc_2040_no_rescreen_monit_0"="ANC 2020-2040 (0,no)",
                                 "anc_2040_with_rescreen_monit_0"="ANC 2020-2040 (0,yes)",
                                 "anc_2040_no_rescreen_monit_5"="ANC 2020-2040 (5,no)",
                                 "monit_0_screen_10a_2030"="Pop 2020+2030 (0,yes)",
                                 "monit_0_screen_10b_2030" ="Pop 2020+2030 (0,no)",
                                 "monit_5_screen_10a_2030"="Pop 2020+2030 (5,yes)",
                                 "monit_5_screen_10b_2030"="Pop 2020+2030 (5,no)",
                                 "screen_2020_anc_monit_0"="ANC 2020 (0)",
                                 "screen_2020_anc_monit_5"="ANC 2020 (5)",
                                 "screen_2020_monit_0"="Pop 2020 (0)",
                                 "screen_2020_monit_5"="Pop 2020 (5)",
                                 "anc_2030_no_rescreen_monit_0_bdvacc"="ANC 2020-2030 (0,no)\n+ birth dose",
                                 "No treatment"="No treatment")) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  guides(shape=FALSE) +
#  geom_abline(slope=391, linetype = "dashed") +
#  geom_abline(slope=518, linetype = "dashed") +
  theme_bw()

# CER plot of non-dominated strategies
ggplot(subset(anc_incremental_df_disc, (frontier=="Non-dominated" | scenario %in%
                c("screen_2020_monit_sim7", "screen_2020_monit_0")) & scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_line(data = subset(anc_incremental_df_disc_median, frontier == "Non-dominated"),
            aes(x = dalys_averted, y = total_cost), size = 1) +
  geom_point(data =  subset(anc_incremental_df_disc_median, frontier == "Non-dominated" |
                              scenario %in%
                              c("screen_2020_monit_sim7", "screen_2020_monit_0")),
             aes(x = dalys_averted, y = total_cost,
                 group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost),
                 shape = frontier),
             size = 5) +
  geom_point(data = subset(anc_incremental_df_disc_median,
                           scenario %in%
                             c("screen_2020_monit_sim7", "screen_2020_monit_0")),
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5, shape = 1,
             colour = "black") +
  geom_point(data = subset(anc_incremental_df_disc_median, frontier == "Non-dominated"),
             aes(x = dalys_averted, y = total_cost,
                 group =reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)),
             size = 5, shape = 2,
             colour = "black") +
  scale_fill_manual(values = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 11))) +
  scale_colour_manual("Screening strategy\n(monitoring frequency,re-testing)",
                      values = c("black",  brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 11)),
                      labels = c("anc_2030_no_rescreen_monit_0"="ANC 2020-2030 (0,no)",
                                 "anc_2040_no_rescreen_monit_0"="ANC 2020-2040 (0,no)",
                                 "monit_5_screen_10a_2030"="Pop 2020+2030 (5,yes)",
                                 "monit_5_screen_10b_2030"="Pop 2020+2030 (5,no)",
                                 "screen_2020_anc_monit_0"="ANC 2020 (0)",
                                 "screen_2020_monit_0"="Pop 2020 (0)",
                                 "screen_2020_monit_sim7"="Pop 2020\n(5 in <45 year olds)",
                                 "pop_2020_anc_2030_no_rescreen_monit_sim7"="Pop 2020+\nANC 2020-2030\n(5 in <45 year olds,no)",
                                 "pop_2020_anc_2040_no_rescreen_monit_0"="Pop 2020+\nANC 2020-2040 (0,no)",
                                 "No treatment"="No treatment")) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  guides(shape=FALSE) +
  #geom_abline(slope=391, linetype = "dashed") +
  #geom_abline(slope=518, linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))


# Impact analysis in cohort and with added birth dose ----

# WHAT PROPORTION OF DALYS AVERTED BY ANC SCREENING IS IN COHORT?

# DALYS averted in cohort
da_in_cohort <- anc_it_2030_monit_0b_cohort$cohort_dalys[,-1]-anc_it_2030_monit_0b$cohort_dalys[,-1]

# DALYS averted in the population (could also use out2n as reference here, same)
da_in_pop <- anc_it_2030_monit_0b_cohort$dalys[[16]][,-c(1:3)]-anc_it_2030_monit_0b$dalys[[16]][,-c(1:3)]

# Ratio:
quantile(da_in_cohort/da_in_pop, c(0.5,0.025,0.975))
# 97% (90-100%)

# WITH added birth dose:
da_in_cohort_bd <- anc_it_2030_monit_0b_cohort$cohort_dalys[,-1]-anc_it_2030_monit_0b_bdvacc$cohort_dalys[,-1]
da_in_pop_bd <- anc_it_2030_monit_0b_cohort$dalys[[16]][,-c(1:3)]-anc_it_2030_monit_0b_bdvacc$dalys[[16]][,-c(1:3)]
# Ratio:
quantile(da_in_cohort_bd/da_in_pop_bd, c(0.5,0.025,0.975))
# 61% (25-89%)
# Makes sense!

# How many extra DALYS and deaths are averted through addition of the BD?
dalys_averted_by_pmtct <-
  plot_hbv_deaths_averted(counterfactual_object = out2n,
                          scenario_objects = list(anc_it_2030_monit_0b,
                                                  anc_it_2030_monit_0b_bdvacc,
                                                  out3_it,
                                                  out8b_it_2030_monit_0,
                                                  out8b_it_2030_monit_5),
                          outcome_to_avert = "dalys",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
dalys_averted_by_pmtct <- subset(dalys_averted_by_pmtct) %>%
  select(scenario, sim,type,by_year, value)
dalys_averted_by_pmtct$sim <- gsub("[^0-9]", "", dalys_averted_by_pmtct$sim)
colnames(dalys_averted_by_pmtct)[5] <- "dalys_averted"

ggplot(dalys_averted_by_pmtct[dalys_averted_by_pmtct$type == "number_averted" &
                                dalys_averted_by_pmtct$by_year==2100,],
       aes(x=reorder(scenario, dalys_averted), y = dalys_averted)) +
  stat_summary(fun.data=f, geom="boxplot", width = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_x_discrete("Screening strategy (monitoring frequency)",
                   labels = c("anc_2030_no_rescreen_monit_0" = "ANC 2020-2030 (0)",
                              "anc_2030_no_rescreen_monit_0_bdvacc" = "ANC 2020-2030 (0) +\nbirth dose scale-up",
                              "screen_2020_monit_0" = "Population-based\n2020 (0)",
                              "monit_0_screen_10b_2030" = "Population-based\n2020+2030 (0)",
                              "monit_5_screen_10b_2030" = "Population-based\n2020+2030 (5)"))+
  ylab("DALYs averted compared to status quo by 2100")

ggplot(dalys_averted_by_pmtct[dalys_averted_by_pmtct$type == "proportion_averted" &
                                dalys_averted_by_pmtct$by_year==2100,],
       aes(x=reorder(scenario, dalys_averted), y = dalys_averted)) +
  stat_summary(fun.data=f, geom="boxplot", width = 0.5) +
  ylim(0,0.75)
# So ANC with birth dose scale-up is quite close in effect to the
# one-time screen without monitoring (add in here the optimal monitoring frequency in young ages).
# Note though that timescale of impact is different:
# By 2030:
ggplot(dalys_averted_by_pmtct[dalys_averted_by_pmtct$type == "proportion_averted" &
                                dalys_averted_by_pmtct$by_year==2030,],
       aes(x=reorder(scenario, dalys_averted), y = dalys_averted)) +
  stat_summary(fun.data=f, geom="boxplot", width = 0.5) +
  ylim(0,0.75)
# By 2050:
ggplot(dalys_averted_by_pmtct[dalys_averted_by_pmtct$type == "proportion_averted" &
                                dalys_averted_by_pmtct$by_year==2050,],
       aes(x=reorder(scenario, dalys_averted), y = dalys_averted)) +
  stat_summary(fun.data=f, geom="boxplot", width = 0.5) +
  ylim(0,0.75)

# Deaths averted
deaths_averted_by_pmtct <-
  plot_hbv_deaths_averted(counterfactual_object = out2n,
                          scenario_objects = list(anc_it_2030_monit_0b,
                                                  anc_it_2030_monit_0b_bdvacc,
                                                  out3_it,
                                                  out8b_it_2030_monit_0,
                                                  out8b_it_2030_monit_5),
                          outcome_to_avert = "cum_hbv_deaths",
                          outcome_to_plot = "number_averted",
                          counterfactual_label = "no treatment")
deaths_averted_by_pmtct <- subset(deaths_averted_by_pmtct,
                                 by_year==2100) %>%
  select(scenario, sim,type, value)
deaths_averted_by_pmtct$sim <- gsub("[^0-9]", "",deaths_averted_by_pmtct$sim)
colnames(deaths_averted_by_pmtct)[4] <- "deaths_averted"

ggplot(deaths_averted_by_pmtct[deaths_averted_by_pmtct$type == "number_averted",],
       aes(x=scenario, y = deaths_averted)) +
  stat_summary(fun.data=f, geom="boxplot", width = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

ggplot(deaths_averted_by_pmtct[deaths_averted_by_pmtct$type == "proportion_averted",],
       aes(x=scenario, y = deaths_averted)) +
  stat_summary(fun.data=f, geom="boxplot", width = 0.5) +
  ylim(0,0.75)
# Pattern for deaths averted is basically the same as for DALYS


# Also look at infections averted and compare this to pop-based




## Compare one-off vs gradual programme (check) ----

# Note that screening in ANC only requires the cost of the test.
# Repeat screening simulations here don't include IT treatment yet so need to compare this
# to the one-off screen from monitoring_frequency:
out3_it <- readRDS(paste0(out_path_monit, "a1_it_out3_screen_2020_monit_0_161220.rds"))
out3_it <- out3_it[[1]]
out5_it <- readRDS(paste0(out_path_monit, "a1_it_out5_screen_2020_monit_5_161220.rds"))
out5_it <- out5_it[[1]]

# Files:
# out3_it: Standard pop-based screen in 2020, no monitoring
# out5_it: Standard pop-based screen in 2020, 5-yearly monitoring
# e1_it_out3_gradual: same pop-based screen but with 10% coverage each year spread over
# 2020-2028, no monitoring
# e1_it_out5_gradual: same pop-based screen but with 10% coverage each year spread over
# 2020-2028, 5-yearly monitoring

annual_discounting_rate <- 0

gradual_interactions_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        assemble_discounted_interactions_for_screening_strategies(out3_it,
                                                                  assessment_object = out3_it)),
  cbind(scenario = "screen_2020_monit_5",
        assemble_discounted_interactions_for_screening_strategies(out5_it,
                                                                  assessment_object = out3_it)),
  cbind(scenario = "screen_2020_gradual_monit_0",
        assemble_discounted_interactions_for_screening_strategies(e1_it_out3_gradual,
                                                                  assessment_object = anc_it_2020_monit_0)),
  cbind(scenario = "screen_2020_gradual_monit_5",
        assemble_discounted_interactions_for_screening_strategies(e1_it_out5_gradual,
                                                                  assessment_object = e1_it_out3_gradual))
)

gradual_py_on_treatment_disc <- rbind(
  data.frame(scenario = "screen_2020_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out3_it,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out5_it,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_gradual_monit_0",
             discount_outcome_2020_to_2100(scenario_object=e1_it_out3_gradual,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_gradual_monit_5",
             discount_outcome_2020_to_2100(scenario_object=e1_it_out5_gradual,
                                           object_to_subtract=NULL,
                                           outcome="py_on_treatment",
                                           yearly_discount_rate=annual_discounting_rate))
)
gradual_py_on_treatment_disc$sim <- gsub("[^0-9]", "", gradual_py_on_treatment_disc$sim)

gradual_dalys_averted_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3_it,
                                      outcome="dalys",
                                      yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out5_it,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_gradual_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=e1_it_out3_gradual,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_gradual_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=e1_it_out5_gradual,
                                           outcome="dalys",
                                           yearly_discount_rate=annual_discounting_rate))
)
gradual_dalys_averted_disc$sim <- gsub("[^0-9]", "", gradual_dalys_averted_disc$sim)
colnames(gradual_dalys_averted_disc)[colnames(gradual_dalys_averted_disc) == "dalys"] <-
  "value"

# Deaths averted
gradual_deaths_averted_disc <- rbind(
  cbind(scenario = "screen_2020_monit_0",
        discount_outcome_2020_to_2100(scenario_object=out2n,
                                      object_to_subtract=out3_it,
                                      outcome="cum_hbv_deaths",
                                      yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=out5_it,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_gradual_monit_0",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=e1_it_out3_gradual,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate)),
  data.frame(scenario = "screen_2020_gradual_monit_5",
             discount_outcome_2020_to_2100(scenario_object=out2n,
                                           object_to_subtract=e1_it_out5_gradual,
                                           outcome="cum_hbv_deaths",
                                           yearly_discount_rate=annual_discounting_rate))
)
gradual_deaths_averted_disc$sim <- gsub("[^0-9]", "", gradual_deaths_averted_disc$sim)
colnames(gradual_deaths_averted_disc)[colnames(gradual_deaths_averted_disc) == "cum_hbv_deaths"] <-
  "value"

# Combine into dataframe
gradual_incremental_df_disc <- create_incremental_plot_df(interactions_df=gradual_interactions_disc,
                                                      py_on_treatment_df=gradual_py_on_treatment_disc,
                                                      deaths_averted_df=gradual_deaths_averted_disc,
                                                      ly_saved_df = gradual_dalys_averted_disc, # replace LY by DALYs
                                                      hbsag_test_cost = 8.3,
                                                      clinical_assessment_cost = 84.4,
                                                      monitoring_assessment_cost = 40.1,
                                                      treatment_py_cost = 60,
                                                      ref_label = "No treatment")
colnames(gradual_incremental_df_disc)[colnames(gradual_incremental_df_disc)=="ly_saved"] <- "dalys_averted"

# Plot ----

# Plot below shows that 10% screening for 8 years only achieves just over half
# the impact and cost as the immediate 90% for a given monitoring frequency
# (with or without discounting) but having a delayed introduction does not appear to
# make the programme less cost-effective, it just requires screening for a longer time.

# Check with no discounting: is this because fewer people have been screened? YES
ggplot(gradual_incremental_df_disc[gradual_incremental_df_disc$scenario != "No treatment",]) +
  geom_boxplot(aes(x=scenario, y = hbsag_tests))
# Is the effect directly proportional to the number of people screened?
gradual_incremental_df_disc$dalys_averted_per_test <-
  gradual_incremental_df_disc$dalys_averted/gradual_incremental_df_disc$hbsag_tests
ggplot(gradual_incremental_df_disc[gradual_incremental_df_disc$scenario != "No treatment",]) +
  geom_boxplot(aes(x=scenario, y = dalys_averted_per_test))
# Tiny bit higher with the immediate programme but overall very similar
# So 10% of the population every year until 2028 is a smaller number than
# 90% of the 2020 population (about 81% of it)

# CER plot
ggplot(subset(gradual_incremental_df_disc, scenario != "No treatment")) +
  stat_ellipse(aes(x = dalys_averted, y = total_cost,
                   group = reorder(scenario, total_cost),
                   fill= reorder(scenario, total_cost)),
               geom = "polygon",
               alpha = 0.2) +
  geom_point(data= group_by(gradual_incremental_df_disc, scenario) %>%
               summarise(dalys_averted=median(dalys_averted),
                         total_cost = median(total_cost)),
             aes(x=dalys_averted, y = total_cost, group = reorder(scenario, total_cost),
                 colour = reorder(scenario, total_cost)), size = 5) +
  scale_fill_manual(values = rev(brewer.pal(4,"RdYlBu"))) +
  scale_colour_manual("Scenarios",
                      #                      labels = c("screen_2020_monit_0" = "2020 screen, no monitoring",
                      #                                 "screen_2020_monit_10" = "2020 screen, monitor every 10 years",
                      #                                 "monit_0_screen_10b_2030" = "2020+2030 screen, no monitoring",
                      #                                 "monit_0_screen_10b_2040" = "2020+2030+2040 screen, no monitoring",
                      #                                 "monit_10_screen_10b_2030" = "2020+2030 screen, monitor every 10 years",
                      #                                 "monit_10_screen_10b_2040" = "2020+2030+2040 screen, monitor every 10 years"),
                      values = c("black", rev(brewer.pal(4,"RdYlBu")))) +
  ylab("Total cost (USD 2019)") +
  xlab("DALYs averted") +
  guides(fill=FALSE) +
  geom_abline(slope=391, linetype = "dashed") +
  geom_abline(slope=518, linetype = "dashed") +
  theme_bw()

# Test: One-off population based vs workplace screening (ignore) ----
plot_hbv_deaths_averted(counterfactual_object = out2,
                        scenario_objects = list(out3,
                                                out3_wpl),
                        counterfactual_label = "no treatment programme",
                        x_axis = "screening")

plot_hbv_deaths_averted_per_healthcare_interaction(counterfactual_object = out2,
                                                   scenario_objects = list(out3,
                                                                           out3_wpl),
                                                   interaction_type = "total_interactions",
                                                   counterfactual_label = "no treatment programme",
                                                   x_axis = "screening")

plot_ly_gained(counterfactual_object = out2,
               scenario_objects = list(out3,
                                       out3_wpl),
               counterfactual_label = "no treatment programme",
               x_axis = "screening")

plot_ly_gained_per_healthcare_interaction(counterfactual_object = out2,
                                          scenario_objects = list(out3,
                                                                  out3_wpl),
                                          interaction_type = "total_interactions",
                                          counterfactual_label = "no treatment programme",
                                          x_axis = "screening")



# Compare this to out5
d1 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                              scenario_objects = list(out3),
                              outcome_to_plot = "number_averted",
                              counterfactual_label = "no treatment")
d2 <- plot_hbv_deaths_averted(counterfactual_object = out2,
                              scenario_objects = list(out3_wpl),
                              outcome_to_plot = "number_averted",
                              counterfactual_label = "no treatment")

i1 <- out3$interactions[[16]]$total_interactions[,-c(1:3)]
i2<- out3_wpl$interactions[[16]]$total_interactions[,-c(1:3)]

deaths_df <- cbind(d1$value[d1$by_year==2100 &
                              d1$type == "number_averted"],
                   d2$value[d2$by_year==2100 &
                              d2$type == "number_averted"],
                   rep(0,183))
int_df <- cbind(unlist(i1), unlist(i2),rep(0,183))


ceef.plot(bcea(e=deaths_df,
               c=int_df,ref=3,interventions=c("General population", "Workplace", "No treatment")),graph="base")


# Test: 1 REPEAT SCREEN VS ALL AGE MONITORING AND ACCESS CHANNELS (ignore) ----
out_path2 <-
  "C:/Users/Nora Schmit/Documents/Model development/hbvmodel - analysis output/monitoring_frequency/"

# No monitoring and 5-yearly monitoring
out3 <- readRDS(paste0(out_path2, "a1_out3_screen_2020_monit_0_201020.rds"))
out3 <- out3[[1]]
out5 <- readRDS(paste0(out_path2, "a1_out5_screen_2020_monit_5_201020.rds"))
out5 <- out5[[1]]

# Compare this to out5
monit_scenario_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                 scenario_objects = list(out5),
                                                 outcome_to_plot = "number_averted",
                                                 counterfactual_label = "no treatment")
screen_scenario_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                  scenario_objects = list(out3),
                                                  outcome_to_plot = "number_averted",
                                                  counterfactual_label = "no treatment")
repscreen_scenario_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                                     scenario_objects = list(out8b_2030),
                                                     outcome_to_plot = "number_averted",
                                                     counterfactual_label = "no treatment")

wpl_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                      scenario_objects = list(out3_wpl),
                                      outcome_to_plot = "number_averted",
                                      counterfactual_label = "no treatment")

anc_deaths <- plot_hbv_deaths_averted(counterfactual_object = out2,
                                      scenario_objects = list(out3_anc),
                                      outcome_to_plot = "number_averted",
                                      counterfactual_label = "no treatment")


monit_scenario_int <- out5$interactions[[16]]$total_interactions[,-c(1:3)]
screen_scenario_int <- out3$interactions[[16]]$total_interactions[,-c(1:3)]
repscreen_scenario_int <- out8b_2030$interactions[[16]]$total_interactions[,-c(1:3)]
wpl_int <- out3_wpl$interactions[[16]]$total_interactions[,-c(1:3)]
anc_int <- out3_anc$interactions[[16]]$total_interactions[,-c(1:3)]

# NOTE monit_scenario is only 1 with monitoring

deaths_df <- cbind(monit_scenario_deaths$value[monit_scenario_deaths$by_year==2100 &
                                                 monit_scenario_deaths$type == "number_averted"],
                   screen_scenario_deaths$value[screen_scenario_deaths$by_year==2100 &
                                                  screen_scenario_deaths$type == "number_averted"],
                   repscreen_scenario_deaths$value[repscreen_scenario_deaths$by_year==2100 &
                                                     repscreen_scenario_deaths$type == "number_averted"],
                   wpl_deaths$value[wpl_deaths$by_year==2100 &
                                      wpl_deaths$type == "number_averted"],
                   anc_deaths$value[anc_deaths$by_year==2100 &
                                      anc_deaths$type == "number_averted"],
                   rep(0,183))
int_df <- cbind(unlist(monit_scenario_int), unlist(screen_scenario_int), unlist(repscreen_scenario_int),
                unlist(wpl_int), unlist(anc_int),
                rep(0,183))

scenario_labels2 <- c("2020 population screen,\n5-yearly monitoring",
                      "2020 population screen,\nno monitoring",
                      "2020+2030 population screen,\nno monitoring",
                      "2020 workplace screen,\nno monitoring",
                      "2020 antenatal screen,\nno monitoring",
                      "No treatment")

ceef.plot.median(bcea(e=deaths_df,c=int_df,ref=ncol(deaths_df),
                      interventions=scenario_labels2),
                 graph="base")

deaths_df <- data.frame(deaths_df)
colnames(deaths_df) <- scenario_labels2
deaths_df$sim <- rownames(deaths_df)
deaths_df <- gather(deaths_df, key = "scenario", value = "deaths_averted", -sim)

int_df <- data.frame(int_df)
colnames(int_df) <- scenario_labels2
int_df$sim <- as.character(seq(1:183))
int_df <- gather(int_df, key = "scenario", value = "interactions", -sim)

combi <- left_join(deaths_df, int_df, by = c("scenario", "sim"))
# Checked that simulation order is the same

combi$frontier <- "Include"
combi$frontier[combi$scenario %in% c("2020 population screen,\nno monitoring",
                                     "2020+2030 population screen,\nno monitoring",
                                     "2020 antenatal screen,\nno monitoring")] <- "Dominated"

combi_summary <- combi %>%
  group_by(scenario, frontier) %>%
  summarise(median_deaths_averted = median(deaths_averted),
            median_interactions = median(interactions))

# Plot
ggplot(combi) +
  geom_line(data= subset(combi, frontier== "Include"),
            aes(x = deaths_averted, y= interactions,
                group = sim), colour = "grey", alpha = 0.3) +
  stat_ellipse(data=subset(combi, scenario != "No treatment"),
               aes(x=deaths_averted,y=interactions,
                   group = reorder(scenario, deaths_averted),
                   fill= reorder(scenario, deaths_averted)),
               geom = "polygon", na.rm = FALSE, alpha = 0.2) +
  geom_point(aes(x = deaths_averted, y = interactions,
                 group =reorder(scenario, deaths_averted), colour = reorder(scenario, deaths_averted)),
             alpha = 0.15) +
  # Overlay median
  geom_line(data = subset(combi_summary, frontier == "Include"),
            aes(y = median_interactions,
                x = median_deaths_averted), size = 1) +
  geom_point(data = combi_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted),
                 colour = reorder(scenario, median_deaths_averted)), size = 5) +
  geom_point(data = combi_summary,
             aes(y = median_interactions,
                 x = median_deaths_averted,
                 group = reorder(scenario, median_deaths_averted)),
             size = 5, shape = 1, colour = "black") +
  scale_fill_manual("Screening scenarios",
                    values=rev(brewer.pal(5,"RdYlBu"))) +
  scale_colour_manual("Screening scenarios",
                      values=c("black", rev(brewer.pal(5,"RdYlBu")))) +
  guides(fill=FALSE,
         colour=guide_legend(keywidth=0.1,keyheight=0.5,default.unit="inch")) +
  xlab("Incremental HBV-related deaths averted") +
  ylab("Incremental number of clinical interactions") +
  #  xlim(-150,7200) +
  #  ylim(0,2500000) +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        title = element_text(size = 15))

###

# In out5, average time on treatment is 22.3 years - assume the same for repeat screen

monit_scenario_cost <- out5$interactions[[16]]$total_screened[,-c(1:3)]*10.38+
  out3$interactions[[16]]$total_assessed[,-c(1:3)]*120+
  (out5$interactions[[16]]$total_assessed[,-c(1:3)]-out3$interactions[[16]]$total_assessed[,-c(1:3)])*15.77+
  out3$interactions[[16]]$total_treated[,-c(1:3)]*22.3*84.88

repscreen_scenario_cost <- out8b_2030$interactions[[16]]$total_screened[,-c(1:3)]*10.38+
  out8b_2030$interactions[[16]]$total_assessed[,-c(1:3)]*120+
  out8b_2030$interactions[[16]]$total_treated[,-c(1:3)]*22.3*84.88

cost_df <- cbind(unlist(monit_scenario_cost), unlist(repscreen_scenario_cost),rep(0,183))

ceef.plot(bcea(e=deaths_df,c=cost_df,ref=3,interventions=c("Monitoring", "Repeat screen", "No treatment")),graph="base")



# Test: Plot HBV deaths rate with re-screening ----
hbv_mort_rate <- rbind(
  gather(out2$timeseries$total_hbv_deaths_rate, key = "sim", value = "rate", -time,-scenario),
  gather(out3$timeseries$total_hbv_deaths_rate, key = "sim", value = "rate", -time,-scenario),
  gather(out8b_2040$timeseries$total_hbv_deaths_rate, key = "sim", value = "rate", -time,-scenario),
  gather(out8b$timeseries$total_hbv_deaths_rate, key = "sim", value = "rate", -time,-scenario))

ggplot(hbv_mort_rate, aes(x=time, y = rate, colour = scenario, fill = scenario)) +
  stat_summary(fun = median, geom = "line") +
  # stat_summary(fun.min = function(x) quantile(x,0.025),
  #               fun.max = function(x) quantile(x,0.975),
  #               geom = "ribbon", alpha = 0.2) +
  xlim(1985,2080) +
  theme_bw()

# Reduction in mortality rate between 2015-2030 (target = 65%)
hbv_mort_rate_wide <- hbv_mort_rate %>%
  pivot_wider(names_from = sim, values_from = rate)

quantile((hbv_mort_rate_wide[hbv_mort_rate_wide$scenario=="status_quo" &
                               hbv_mort_rate_wide$time==2015,-c(1:2)]-
            hbv_mort_rate_wide[hbv_mort_rate_wide$scenario=="status_quo" &
                                 hbv_mort_rate_wide$time==2030,-c(1:2)])/
           hbv_mort_rate_wide[hbv_mort_rate_wide$scenario=="status_quo" &
                                hbv_mort_rate_wide$time==2015,-c(1:2)],
         prob=c(0.5,0.025,0.975))

# Absolute target (CDA suggestion): < ≤5 per 100 000
quantile((hbv_mort_rate_wide[hbv_mort_rate_wide$scenario=="status_quo" &
                               hbv_mort_rate_wide$time==2030,-c(1:2)]),
         prob=c(0.5,0.025,0.975))*100000

# With screen and treat
quantile((hbv_mort_rate_wide[hbv_mort_rate_wide$scenario=="screen_2020_monit_0" &
                               hbv_mort_rate_wide$time==2015,-c(1:2)]-
            hbv_mort_rate_wide[hbv_mort_rate_wide$scenario=="screen_2020_monit_0" &
                                 hbv_mort_rate_wide$time==2030,-c(1:2)])/
           hbv_mort_rate_wide[hbv_mort_rate_wide$scenario=="screen_2020_monit_0" &
                                hbv_mort_rate_wide$time==2015,-c(1:2)],
         prob=c(0.5,0.025,0.975))
# Almost reaches reduction elimination target

# Absolute target (CDA suggestion): < ≤5 per 100 000
quantile((hbv_mort_rate_wide[hbv_mort_rate_wide$scenario=="screen_2020_monit_0" &
                               hbv_mort_rate_wide$time==2030,-c(1:2)]),
         prob=c(0.5,0.025,0.975))*100000
# Absolute target is reached with the one-time screening programme (reduced from status quo)

quantile((out2$timeseries$total_chronic_infections_rate[
  out2$timeseries$total_chronic_infections_rate$time==2015,-c(1:2)]-
    out2$timeseries$total_chronic_infections_rate[
      out2$timeseries$total_chronic_infections_rate$time==2030,-c(1:2)])/
    out2$timeseries$total_chronic_infections_rate[
      out2$timeseries$total_chronic_infections_rate$time==2015,-c(1:2)],
  prob = c(0.5,0.025,0.975))

# out8b
# out8b_2040


