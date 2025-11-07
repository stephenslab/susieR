context("test_univariate_regression.R")

test_that(paste("Check that variants of univaraite_regression produce",
                "the same result"),{

  dat <- simulate(200,500)
  X <- dat$X
  y <- dat$y
  res1 <- univariate_regression(X,y,method = "lmfit")
  res2 <- univariate_regression(X,y,method = "sumstats")
  expect_equal(res1,res2,tolerance = 1e-8)
})
