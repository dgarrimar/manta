context("eigenvalue computation")

test_that("eigenvalues are positive", {
  dmat <- dist(cbind(rnorm(100),rnorm(100)), method = "canberra")
  expect_error(mlm(dmat ~ runif(100)), "all eigenvalues of G should be > 0")
})

