context("test_susie_small.R")

test_that(paste("check that ELBO is monotonically increasing for ",
                "estimate_residual_method = 'Servin_Stephens', ",
                "with L = 1"),{
  set.seed(1)
  data(data_small)
  y <- data_small$y
  X <- data_small$X
  fit <- susie(X,y,L = 1,estimate_residual_method = "Servin_Stephens",
               alpha0 = 0.1,beta0 = 0.1,tol = 1e-6,verbose = TRUE)
  expect_true(all(diff(fit$elbo) >= 0))
})
