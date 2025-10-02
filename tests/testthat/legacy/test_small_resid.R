# context("test_small_resid.R")

# This test comes from Issue #7.
# test_that("susie works when residual variance is small",{
#   dat <- readRDS("full_data_1_sim_gaussian_null_1.rds")$data
#   fit <- susie(dat$X,dat$Y,L = 10,null_weight = 0)
#   expect_equal(fit$sigma2,1e-4)
#   expect_true(all(is.finite(fit$elbo)))
# })
