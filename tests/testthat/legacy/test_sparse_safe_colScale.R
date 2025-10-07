context("test_sparse_set_X_attributes.R")

test_that("sparse version set_X_attributes", with(simulate(sparse=T), {
  
  dense.res = set_X_attributes(X)
  sparse.res = set_X_attributes(X.sparse)
  
  expect_equal(sparse.res, X.sparse, check.attributes=FALSE)
  expect_equal(attr(dense.res, 'scaled:center'), attr(sparse.res, 'scaled:center'))
  expect_equal(attr(dense.res, 'scaled:scale'), attr(sparse.res, 'scaled:scale'))
  expect_equal(attr(dense.res, 'd'), attr(sparse.res, 'd'))
}))
