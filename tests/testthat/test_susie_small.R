context("test_susie_small.R")

test_that(paste("check that susie-small runs correctly in a data set",
                "with 1 causal SNP"),{
  set.seed(1)
  n <- 70
  p <- 100
  X <- matrix(rnorm(n*p),nrow = n,ncol = p)
  X <- cbind(X[,1],X)
  y <-  X[,1] + rnorm(n)
  res1 <- susie(X,y,L = 1,small = FALSE,max_iter = 10,min_abs_corr = 0,
                verbose = TRUE)
  res2 <- susie(X,y,L = 1,small = TRUE,max_iter = 10,min_abs_corr = 0,
                verbose = TRUE)
  expect_equal(res1$sets,res2$sets,scale = 1,tolerance = 1e-4)
})
