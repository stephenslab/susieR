context("test_Eloglik.R")

test_that("Eloglik agrees with version 0.3", with(simulate(sparse = TRUE), {
  original.res = readRDS('Eloglik_original_res.rds')
  
  
  scaledX = set_X_attributes(X)
  scaledX.sparse = set_X_attributes(X.sparse)
  
  dense.res  = Eloglik(scaledX, y, s)
  sparse.res = Eloglik(scaledX.sparse, y, s)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
}))
