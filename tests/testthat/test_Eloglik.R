test_that("sparse version Eloglik",{
  original.res = load_data('Eloglik_original_res.rds')
  simulate(sparse=T)
  
  scaledX = safe_colScale(X)
  scaledX.sparse = safe_colScale(X.sparse)
  
  dense.res = Eloglik(scaledX, y, s)
  sparse.res = Eloglik(scaledX.sparse, y, s)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
  
})
