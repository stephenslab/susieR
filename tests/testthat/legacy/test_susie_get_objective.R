context("test_get_objective.R")

test_that("get_objective agrees with version 0.3", with(simulate(sparse = TRUE), {
  original.res = readRDS('objective_original_res.rds')

  scaledX = set_X_attributes(X)
  scaledX.sparse = set_X_attributes(X.sparse)
  s$Xr = colSums(compute_MXt(s$alpha*s$mu, scaledX))

  dense.res = get_objective(scaledX, y, s)
  sparse.res = get_objective(scaledX.sparse, y, s)

  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)

}))