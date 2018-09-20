context("test_sparse_safe_colScale.R")

test_that("sparse version safe_colScale",{
  simulate(sparse=T)
  
  dense.res = safe_colScale(X)
  sparse.res = safe_colScale(X.sparse)
  
  dense.faceX = dense.res
  attributes(dense.faceX) = NULL
  attributes(X) = NULL
  
  sparse.faceX = sparse.res
  attributes(sparse.faceX) = NULL
  attributes(X.sparse) = NULL
  
  expect_equal(sparse.faceX, X.sparse)
  expect_equal(attr(dense.res, 'scaled:center'), attr(sparse.res, 'scaled:center'))
  expect_equal(attr(dense.res, 'scaled:scale'), attr(sparse.res, 'scaled:scale'))
  expect_equal(attr(dense.res, 'd'), attr(sparse.res, 'd'))
})
