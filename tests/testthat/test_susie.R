test_that("sparse version susie",{
  original.res = load_data('susiefit_original_res.rds')
  original.res2 = load_data('susiefit_original_res2.rds')
  original.res3 = load_data('susiefit_original_res3.rds')
  original.res4 = load_data('susiefit_original_res4.rds')
  simulate(sparse=T)
  
  dense.res = susie(X, y)
  sparse.res = susie(X.sparse, y)
  
  dense.res2 = susie(X, y, standardize=TRUE, intercept = FALSE)
  sparse.res2 = susie(X.sparse, y, standardize=TRUE, intercept = FALSE)
  
  dense.res3 = susie(X, y, standardize=FALSE, intercept = TRUE)
  sparse.res3 = susie(X.sparse, y, standardize=FALSE, intercept = TRUE)
  
  dense.res4 = susie(X, y, standardize=FALSE, intercept = FALSE)
  sparse.res4 = susie(X.sparse, y, standardize=FALSE, intercept = FALSE)
  
  is_equal_susie(sparse.res, original.res)
  is_equal_susie(dense.res, original.res)
  is_equal_susie(sparse.res2, original.res2)
  is_equal_susie(dense.res2, original.res2)
  is_equal_susie(sparse.res3, original.res3)
  is_equal_susie(dense.res3, original.res3)
  is_equal_susie(sparse.res4, original.res4)
  is_equal_susie(dense.res4, original.res4)
})
