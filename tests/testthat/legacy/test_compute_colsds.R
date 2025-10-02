context("test_compute_colsds.R")

test_that(paste("Check that compute_colSds gives same result as",
                "matrixStats::colSds for a dense matrix"),{
  set.seed(1)
  X  <- matrix(rnorm(24),4,6)
  y1 <- compute_colSds(X)
  y2 <- matrixStats::colSds(X)
  expect_equal(y1,y2,scale = 1,tolerance = 1e-15)
})

test_that(paste("Check that compute_colSds gives same result as",
                "matrixStats::colSds for a sparse matrix"),{
  set.seed(1)
  X  <- matrix(rnorm(24),4,6)
  X[runif(24) < 0.3] <- 0
  Y  <- as(X,"CsparseMatrix")
  y1 <- compute_colSds(Y)
  y2 <- matrixStats::colSds(X)
  expect_equal(y1,y2,scale = 1,tolerance = 1e-15)
})
