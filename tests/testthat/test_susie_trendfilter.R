context("test_susie_trendfilter.R")

#check if the susie_trendfilter produces the same results as the original version
test_that("susie for trend filtering",{
  #order 0 trend filtering
  simulate_tf(0)
  original.s0 = susie(X,y)
  s0 = susie_trendfilter(y,0)
  original.s0.ns = susie(X,y,standardize=FALSE,intercept=FALSE)
  s0.ns = susie_trendfilter(y,0,normalize=FALSE)
  #order 1 trend filtering
  simulate_tf(1)
  original.s1 = susie(X,y)
  s1 = susie_trendfilter(y,1)
  original.s1.ns = susie(X,y,standardize=FALSE,intercept=FALSE)
  s1.ns = susie_trendfilter(y,1,normalize=FALSE)
  #order 2 trend filtering
  simulate_tf(2)
  original.s2 = susie(X,y)
  s2 = susie_trendfilter(y,2)
  original.s2.ns = susie(X,y,standardize=FALSE,intercept=FALSE)
  s2.ns = susie_trendfilter(y,2,normalize=FALSE)
  
  expect_equal_susie(s0, original.s0)
  expect_equal_susie(s1, original.s1)
  expect_equal_susie(s2, original.s2)
  expect_equal_susie(s0.ns, original.s0.ns)
  expect_equal_susie(s1.ns, original.s1.ns)
  expect_equal_susie(s2.ns, original.s2.ns)
})

#check if susie_trendfilter can fit different arguments
test_that("susie_trendfilter interface",{
  simulate_tf(0)
  expect_error(susie_trendfilter(y,0,estimate_prior_variance=TRUE), NA)
  expect_error(susie_trendfilter(y,0,null_weight=1/(n+1)), NA)
  expect_error(susie_trendfilter(y,0,estimate_prior_variance=TRUE, null_weight=1/(n+1)), NA)
  expect_error(susie_trendfilter(y,0,intercept=TRUE, standardize=FALSE), NA)
  expect_error(susie_trendfilter(y,0,intercept=FALSE, standardize=FALSE), NA)
  expect_error(susie_trendfilter(y,0,intercept=FALSE, standardize=TRUE), NA)
  expect_error(susie_trendfilter(y,0,compute_univariate_zscore =TRUE), NA)
  
})