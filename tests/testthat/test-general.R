context("I/O format")

test_that("Response of class 'data.frame' raises an error", {
  expect_error(manta(as.data.frame(biomarkers) ~ ., data = patients), 
               "invalid type \\(list\\)")
})

test_that("A single response raises an error", {
  expect_error(manta(biomarkers[,1] ~ ., data = patients), 
               "number of response variables")
})

context("Models")

test_that("Null and intercept-only models raise an error", {
  e <- "at least one predictor"
  expect_error(manta(biomarkers ~ -1, data = patients), e)
  expect_error(manta(biomarkers ~ 0, data = patients), e)
  expect_error(manta(biomarkers ~ 1, data = patients), e)
})

test_that("na.action in manta is 'na.omit'", {
  na.default <- options("na.action")$na.action
  options(na.action = "na.fail")
  expect_error(manta(biomarkers ~ ., data = patients), NA)
  options(na.action = na.default)
})

context("Other arguments")

test_that("Transformation generating NaNs raises an error", {
  expect_error(manta(biomarkers ~ ., data = patients, transform = "log"),
               "transformation requires")
  expect_error(manta(biomarkers ~ ., data = patients, transform = "sqrt"),
               "transformation requires")
})

test_that("Summary stats are identical for subsets of terms", {
  res.full <- manta(biomarkers ~ .^2, data = patients)
  subset <- c("age", "gender:status")
  res.subset <- manta(biomarkers ~ .^2, data = patients, subset = subset)
  expect_identical(res.full$aov.tab[subset, ], res.subset$aov.tab[subset, ])
})

test_that("Subset with unknown terms raises an error", {
  subset <- c("age", "UNKOWN")
  expect_error(manta(biomarkers ~ .^2, data = patients, subset = subset), 
               "Unknown terms in subset")
})

context("Results on example dataset")

test_that("Current output matches expected results", {
  res <- manta(biomarkers ~ .^2, data = patients, fit = TRUE)
  expect_equal_to_reference(res, "example_ref.rds")
  
  res <- manta(biomarkers ~ .^2, data = patients, type = "I", fit = TRUE)
  expect_equal_to_reference(res, "example_ref_I.rds")
  
  res <- manta(biomarkers ~ .^2, data = patients, type = "III", fit = TRUE)
  expect_equal_to_reference(res, "example_ref_III.rds")
})

context("Asymptotic P-value calculation")

lambda <- c(0.507478, 0.4269077, 0.3067685)
eps0 <- 1e-16
pv <- p.asympt(ss = 100, df = 1, lambda = lambda, eps = eps0)
test_that("Precision is updated and AS204 retried", expect_gt(pv[2], eps0))
test_that("P-values below precision become precision", expect_identical(pv[1], pv[2]))

test_that("Precision above a threshold raises an error", {
  expect_error(p.asympt(ss = 100, df = 1, lambda = lambda, eps = eps0, eps.stop = eps0), 
               "Precision of asymptotic P-value")
})

test_that("AS 204 fault indicators other than 0, 4, 5 and 9 raise an error", {
  expect_error(AS204(1, c(-1, 2)), "fault indicator: -1")
  expect_error(AS204(0, 1:3), "fault indicator: 2")
  expect_error(AS204(20, c(30, 1), mult = c(1, 30), mode = 0), "fault indicator: 3")
  expect_error(AS204(100, c(0.241, 0.103, 0.099), mode = 0), "fault indicator: 6")
})

test_that("AS 204 fault indicators 4, 5 and 9 return NULL",{
  expect_null(AS204(100, 1:3, maxit = 10)) # fault indicator 4
  expect_null(AS204(100, c(0.507478, 0.4269077, 0.3067685), eps = 1e-16)) # fault indicator 5
  expect_null(AS204(100, c(0.7334573, 0.7150620, 0.7078593), eps = 1e-16, maxit = 10)) # fault indicator 9
})
