
n=70
p=100

X= matrix(rnorm(n*p),ncol=p)
y= X[,1]+ rnorm(n)
test_that("check that susie small runs correctly",   {
  res1 = susie(X, y, small = TRUE)
  expect_equal(res1$sets$cs[[1]],1)
})
