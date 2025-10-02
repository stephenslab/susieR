context("test_init.R")

test_that("susie init", with(simulate(sparse=T), {
  original.res  = readRDS('susiefit_original_res.rds')
  original.res$lbf_variable = matrix(0, nrow(original.res$alpha), ncol(original.res$alpha))
  expect_message(capture_output(susie(X, y, L=2, tol=1E-2, s_init = original.res, estimate_prior_variance = FALSE,verbose=T)))
  expect_message(capture_output(susie(X, y, L=20, tol=1E-2, s_init = original.res, estimate_prior_variance = FALSE, verbose=T)))
}))
