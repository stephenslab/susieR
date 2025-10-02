context("test_intercept_standardize.R")

test_that("scaling and intercept works as expected", with(simulate(200,1000), {
  s1 = susie(X,y,intercept= TRUE, standardize=TRUE)
  s2 = susie(X,y,intercept = FALSE, standardize = FALSE)
  s3 = susie(X,y,intercept =TRUE, standardize = FALSE)
  s4 = susie(X,y,intercept = FALSE,standardize = TRUE)
  
  expect_equal(predict(s2),predict(s2,X))
  expect_equal(predict(s4),predict(s4,X))
  expect_equal(predict(s1),predict(s1,X))
  expect_equal(predict(s3),predict(s3,X))
  
  expect_equal(s2$intercept, 0)
  expect_equal(s4$intercept, 0)
}))
