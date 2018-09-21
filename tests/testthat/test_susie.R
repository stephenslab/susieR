context("test_susie.R")

test_that("susie agrees with version 0.3",{
  original.res = load_data('susiefit_original_res.rds')
  original.res2 = load_data('susiefit_original_res2.rds')
  original.res3 = load_data('susiefit_original_res3.rds')
  original.res4 = load_data('susiefit_original_res4.rds')
  simulate(sparse=T)
  
  dense.res = susie(X, y, tol=1E-2)
  sparse.res = susie(X.sparse, y, tol=1E-2)
  
  dense.res2 = susie(X, y, standardize=TRUE, intercept = FALSE, tol=1E-2)
  sparse.res2 = susie(X.sparse, y, standardize=TRUE, intercept = FALSE,
      tol=1E-2)
  
  dense.res3 = susie(X, y, standardize=FALSE, intercept = TRUE, tol=1E-2)
  sparse.res3 = susie(X.sparse, y, standardize=FALSE, intercept = TRUE,
      tol=1E-2)
  
  dense.res4 = susie(X, y, standardize=FALSE, intercept = FALSE, tol=1E-2)
  sparse.res4 = susie(X.sparse, y, standardize=FALSE, intercept = FALSE,
      tol=1E-2)
  
  expect_equal_susie(sparse.res, original.res)
  expect_equal_susie(dense.res, original.res)
  expect_equal_susie(sparse.res2, original.res2)
  expect_equal_susie(dense.res2, original.res2)
  expect_equal_susie(sparse.res3, original.res3)
  expect_equal_susie(dense.res3, original.res3)
  expect_equal_susie(sparse.res4, original.res4)
  expect_equal_susie(dense.res4, original.res4)
})
