context("test_susie_trendfilter.R")

#check if the susie_trendfilter produces the same results as the original version
test_that("susie for trend filtering",{
  original.s0  = readRDS('tf_original_s0.rds')
  original.s1 = readRDS('tf_original_s1.rds')
  original.s2 = readRDS('tf_original_s2.rds')
  #order 0 trend filtering
  n = 1000
  set.seed(1)
  beta = c(rep(0,100),rep(1,100),rep(3,100),rep(-2,100),rep(0,600))
  y0 = beta + rnorm(n)
  s0 = susie_trendfilter(y0,0)
  #order 1 trend filtering
  set.seed(1)
  beta = numeric(n)
  for (i in 1:n){
    if (i <= 100){
      beta[i] = 0.001*i + 2
    } else if (i <= 300){
      beta[i] = 5*0.001*i + 1.6
    } else{
      beta[i] = 6.1 - 10*0.001*i
    }
  }
  y1 = beta + rnorm(n)
  s1 = susie_trendfilter(y1,1)
  #order 2 trend filtering
  set.seed(1)
  beta = numeric(n)
  for (i in 1:n){
    if (i <= 100){
      beta[i] = (0.001*i)^2
    } else if (i <= 700){
      beta[i] = -5*(0.001*i)^2 + 0.06
    } else{
      beta[i] = 3*(0.001*i)^2 - 3.86
    }
  }
  y2 = beta + rnorm(n)
  s2 = susie_trendfilter(y2,2)
  
  expect_equal_susie(s0, original.s0)
  expect_equal_susie(s1, original.s1)
  expect_equal_susie(s2, original.s2)
})

#check if susie_trendfilter can fit different arguments
test_that("susie_trendfilter interface",{
  n = 1000
  set.seed(1)
  beta = c(rep(0,100),rep(1,100),rep(3,100),rep(-2,100),rep(0,600))
  y0 = beta + rnorm(n)
  
  expect_error(susie_trendfilter(y0,0,estimate_prior_variance=TRUE), NA)
  expect_error(susie_trendfilter(y0,0,null_weight=1/(n+1)), NA)
  expect_error(susie_trendfilter(y0,0,estimate_prior_variance=TRUE, null_weight=1/(n+1)), NA)
  expect_error(susie_trendfilter(y0,0,intercept=TRUE, standardize=FALSE), NA)
  expect_error(susie_trendfilter(y0,0,intercept=FALSE, standardize=FALSE), NA)
  expect_error(susie_trendfilter(y0,0,intercept=FALSE, standardize=TRUE), NA)
  expect_error(susie_trendfilter(y0,0,compute_univariate_zscore =TRUE), NA)
  
})