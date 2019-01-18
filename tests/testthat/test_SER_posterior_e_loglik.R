context("test_SER_posterior_e_loglik.R")

test_that("SER_posterior_e_loglik agrees with version 0.3", with(simulate(sparse=T), {
  original.res = readRDS('SER_original_res.rds')
  Eb = rep(1, p)
  Eb2 = rep(1, p)
  s2 = s$sigma2
  
  scaledX = set_X_attributes(X)
  scaledX.sparse = set_X_attributes(X.sparse)
  
  dense.res = SER_posterior_e_loglik(scaledX,y,s2,Eb,Eb2)
  sparse.res = SER_posterior_e_loglik(scaledX.sparse,y,s2,Eb,Eb2)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
  
}))
