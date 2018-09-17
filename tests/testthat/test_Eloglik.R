test_that("sparse version Eloglik",{
  original.res = load_data('Eloglik_original_res.rds')
  simulate(sparse=T)
  
  scaledX = susieR:::safe_colScale(X)
  scaledX.sparse = susieR:::safe_colScale(X.sparse)
  
  dense.res = susieR:::Eloglik(scaledX, y, s)
  sparse.res = susieR:::Eloglik(scaledX.sparse, y, s)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
  
})
