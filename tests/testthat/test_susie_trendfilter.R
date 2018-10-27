context("test_susie_trendfilter.R")

#check if the susie_trendfilter produces the same results as the original version
test_that("susie for trend filtering",{
  #order 0 trend filtering
  simulate_tf(0)
  original.s0 = susie(X,y)
  s0 = susie_trendfilter(y,0)
  #order 1 trend filtering
  simulate_tf(1)
  original.s1 = susie(X,y)
  s1 = susie_trendfilter(y,1)
  #order 2 trend filtering
  simulate_tf(2)
  original.s2 = susie(X,y)
  s2 = susie_trendfilter(y,2)
  
  expect_equal_susie(s0, original.s0)
  expect_equal_susie(s1, original.s1)
  expect_equal_susie(s2, original.s2)
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