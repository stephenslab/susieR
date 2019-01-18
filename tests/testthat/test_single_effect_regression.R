context("test_single_effect_regression.R")

test_that("single_effect_regression agrees with version 0.3", with(simulate(sparse=T), {
  original.res = readRDS('singleReg_original_res.rds')
  V = 0.2
  
  scaledX = set_X_attributes(X)
  scaledX.sparse = set_X_attributes(X.sparse)
  
  dense.res = single_effect_regression(y,scaledX,V)
  sparse.res = single_effect_regression(y,scaledX.sparse,V)
  
  expect_equal_SER(sparse.res, original.res)
  expect_equal_SER(dense.res, original.res)
}))

