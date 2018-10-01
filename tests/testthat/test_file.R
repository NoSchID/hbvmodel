context("Model works?")


# Test 1: Check that the number of people in each compartment is always positive
test_that("all compartments positive", {
  output <- run_model(init_pop, edmunds_parameters)
  expect_true(all(output >= 0))
})
