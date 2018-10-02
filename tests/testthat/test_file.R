context("Basic transmission dynamics")

# Transmission dynamics tests

# 1) The number of people in each compartment,
# the cumulative incidence and the prevalence are always positive
test_that("all compartments positive", {
  test_output <- run_model(init_pop, parameters)
  expect_true(all(test_output >= 0))
})

# 2) If there is no transmission, the incidence is 0 throughout
test_that("no transmission = no infections", {
  parameters["beta"] <- 0
  test_output <- run_model(init_pop, parameters)
  expect_true(all(test_output$cum_incid == 0))
})

# 3) If there are no infectious individuals at t0, the incidence is 0 throughout
test_that("no infectious individuals = no infections", {
  init_pop["A0"] <- 0
  init_pop["I0"] <- 0
  test_output <- run_model(init_pop, parameters)
  expect_true(all(test_output$cum_incid == 0))
})

# 4) The population size remains constant if there is no HBV-specific mortality
test_that("no HBV-specific mortality = constant population size", {
  parameters["mu_hbv"] <- 0
  test_output <- run_model(init_pop, parameters)
  expect_equal(test_output[length(test_output$total_pop),]$total_pop, test_output[1,]$total_pop)
})

