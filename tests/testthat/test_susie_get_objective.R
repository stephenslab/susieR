test_that("sparse version susie_get_objective",{
  original.res = load_data('objective_original_res.rds')
  simulate(sparse=T)
  
  scaledX = susieR:::safe_colScale(X)
  scaledX.sparse = susieR:::safe_colScale(X.sparse)
  
  dense.res = susieR:::susie_get_objective(scaledX, y, s)
  sparse.res = susieR:::susie_get_objective(scaledX.sparse, y, s)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
  
})
