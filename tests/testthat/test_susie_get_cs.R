context("test_susie_get_cs.R")

test_that("susie_get_cs purity calculations are correct",{

  # Simulate a data set with correlated variables.
  dat <- simulate(40,100)
  X   <- dat$X
  y   <- dat$y
  X   <- cbind(X,
               X + matrix(rnorm(4000),40,100)/20,
               X + matrix(rnorm(4000),40,100)/20)

  # Fit a susie model.
  fit <- susie(X,y,estimate_prior_variance = FALSE)

  # The purity calculations should be the same whether or not the
  # Rfast package functions are used, and all the purity statistics
  # should be positive.
  set.seed(1)
  purity1 <- susie_get_cs(fit,X,min_abs_corr = 0,use_rfast = FALSE)$purity
  expect_gt(min(purity1),0)
  skip_if_not_installed("Rfast")
  set.seed(1)
  purity2 <- susie_get_cs(fit,X,min_abs_corr = 0,use_rfast = TRUE)$purity
  expect_equal(purity1,purity2,scale = 1,tolerance = 1e-15)
})
