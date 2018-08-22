test_that("scaling and intercept works as expected",{
  set.seed(1)
  x = matrix(rnorm(2000,3,4),ncol=10)
  b = rnorm(10)
  y = x %*% b + rnorm(200,0,0.1)
  s1 = susie(x,y,intercept= TRUE, standardize=TRUE)
  s2 = susie(x,y,intercept = FALSE, standardize = FALSE)
  s3 = susie(x,y,intercept =TRUE, standardize = FALSE)
  s4 = susie(x,y,intercept = FALSE,standardize = TRUE)

  expect_equal(predict(s2),predict(s2,x))
  expect_equal(predict(s4),predict(s4,x))
  expect_equal(predict(s1),predict(s1,x))
  expect_equal(predict(s3),predict(s3,x))

  expect_equal(s2$intercept, 0)
  expect_equal(s4$intercept, 0)

})
