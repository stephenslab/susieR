context("test_sparse_safe_colScale.R")

test_that("sparse version safe_colScale", with(simulate(sparse=T), {
  
  dense.res = safe_colScale(X)
  sparse.res = safe_colScale(X.sparse)
  
  expect_equal(sparse.res, X.sparse, check.attributes=FALSE)
  expect_equal(attr(dense.res, 'scaled:center'), attr(sparse.res, 'scaled:center'))
  expect_equal(attr(dense.res, 'scaled:scale'), attr(sparse.res, 'scaled:scale'))
  expect_equal(attr(dense.res, 'd'), attr(sparse.res, 'd'))
}))
