context("Results on example dataset")

test_that("Current output matches expected results", {
  fit <- mlm(biomarkers ~ .^2, data = patients, fit = TRUE)
  expect_equal_to_reference(fit, "example_ref.rds")
})
