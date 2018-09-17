test_that("sparse version update_each_effect",{
  original.res = load_data('vbupdate_original_res.rds')
  simulate(sparse=T)
  
  scaledX = susieR:::safe_colScale(X)
  scaledX.sparse = susieR:::safe_colScale(X.sparse)
  
  dense.res = susieR:::update_each_effect(scaledX,y,s)
  sparse.res = susieR:::update_each_effect(scaledX.sparse,y,s)
  
  is_equal_s(sparse.res, original.res)
  is_equal_s(dense.res, original.res)
})
