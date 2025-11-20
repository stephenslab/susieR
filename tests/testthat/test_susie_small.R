context("test_susie_small.R")

test_that(paste("check that susie-small runs correctly in a data set",
                "with 1 causal SNP"),{
  skip()
  set.seed(1)
  n <- 70
  p <- 100
  X <- matrix(rnorm(n*p),nrow = n,ncol = p)
  X <- cbind(X[,1],X)
  y <-  X[,1] + rnorm(n)
  capture.output(res1 <- susie(X,y,L = 1,small = FALSE,max_iter = 10,
                               min_abs_corr = 0,verbose = TRUE))
  capture.output(res2 <- susie(X,y,L = 1,small = TRUE,max_iter = 10,
                               min_abs_corr = 0,verbose = TRUE))
  capture.output(res3 <- susie(X,y,small = FALSE,max_iter = 10,
                               verbose = TRUE))
  suppressWarnings(capture.output(
    res4 <- susie(X,y,small = TRUE,max_iter = 10,verbose = TRUE)))
  expect_equal(res1$sets,res2$sets,scale = 1,tolerance = 1e-4)
  expect_equal(res3$sets,res4$sets,scale = 1,tolerance = 1e-4)
})
