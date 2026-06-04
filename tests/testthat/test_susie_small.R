context("test_susie_small")

# The NIG + L=1 ELBO monotonicity check on data_small has been merged into
# test_susie.R ("susie ELBO is monotonically increasing for NIG residuals
# with L = 1 (data_small)").  This file is retained as required but adds
# nothing new.
test_that("data_small is loadable and has expected structure", {
  data(data_small)
  expect_true(!is.null(data_small$X))
  expect_true(!is.null(data_small$y))
  expect_equal(length(data_small$y), nrow(data_small$X))
})
