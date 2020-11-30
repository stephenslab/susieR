context("test_small_resid.R")

# This test comes from Issue #7.
test_that("susie works when residual variance is small",{
  data <- readRDS("full_data_1_sim_gaussian_null_1.rds")
  fit <- susie(data$data$X,data$data$Y,L = 10,null_weight = 0,
               residual_variance_lowerbound = 1e-15)
  expect_true(all(is.finite(fit$elbo)))
})
