context("test_init.R")

test_that("susie init", with(simulate(sparse=T), {
  original.res  = readRDS('susiefit_original_res.rds')
  expect_warning(susie(X, y, L=2, tol=1E-2, s_init = original.res, estimate_prior_variance = FALSE))
  expect_error(susie(X, y, L=20, tol=1E-2, s_init = original.res, estimate_prior_variance = FALSE))
}))