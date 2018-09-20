context("test_get_ER2.R")

test_that("get_ER2 agrees with version 0.3",{
  original.res = load_data('ER2_original_res.rds')
  simulate(sparse = TRUE)
  
  scaledX = safe_colScale(X)
  scaledX.sparse = safe_colScale(X.sparse)
  
  dense.res = get_ER2(scaledX, y, s)
  sparse.res = get_ER2(scaledX.sparse, y, s)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
})
