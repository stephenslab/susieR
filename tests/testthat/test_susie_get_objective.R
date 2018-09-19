test_that("susie_get_objective agrees with version 0.3",{
  original.res = load_data('objective_original_res.rds')
  simulate(sparse = TRUE)
  
  scaledX = safe_colScale(X)
  scaledX.sparse = safe_colScale(X.sparse)
  
  dense.res = susie_get_objective(scaledX, y, s)
  sparse.res = susie_get_objective(scaledX.sparse, y, s)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
 
})
