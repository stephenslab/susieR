test_that("sparse version get_ER2",{
  original.res = load_data('ER2_original_res.rds')
  simulate(sparse=T)
  
  scaledX = susieR:::safe_colScale(X)
  scaledX.sparse = susieR:::safe_colScale(X.sparse)
  
  dense.res = susieR:::get_ER2(scaledX, y, s)
  sparse.res = susieR:::get_ER2(scaledX.sparse, y, s)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
  
})
