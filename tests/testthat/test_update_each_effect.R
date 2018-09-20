context("test_update_each_effect.R")

test_that("update_each_effect agrees with version 0.3",{
  original.res = load_data('vbupdate_original_res.rds')
  simulate(sparse=T)
  
  scaledX = safe_colScale(X)
  scaledX.sparse = safe_colScale(X.sparse)
  
  dense.res = update_each_effect(scaledX,y,s)
  sparse.res = update_each_effect(scaledX.sparse,y,s)
  
  expect_equal_susie_update(sparse.res, original.res)
  expect_equal_susie_update(dense.res, original.res)
})
