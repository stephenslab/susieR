context("test_susie_trendfilter.R")

#check if the susie_trendfilter produces the same results as the original version
test_that("susie for trend filtering",{
  #order 0 trend filtering
  with(simulate_tf(0), {
    original.s0 = susie(X,y,estimate_prior_variance = FALSE)
    s0 = susie_trendfilter(y,0,estimate_prior_variance = FALSE, standardize = TRUE, use_mad=FALSE)
    original.s0.ns.ni = susie(X,y,standardize=FALSE,intercept=FALSE,estimate_prior_variance = FALSE)
    s0.ns.ni = susie_trendfilter(y,0,standardize=FALSE,intercept=FALSE,estimate_prior_variance = FALSE, use_mad=FALSE)
    original.s0.ns = susie(X,y,standardize=FALSE,intercept=TRUE,estimate_prior_variance = FALSE)
    s0.ns = susie_trendfilter(y,0,standardize=FALSE,intercept=TRUE,estimate_prior_variance = FALSE, use_mad=FALSE)
    original.s0.ni = susie(X,y,standardize=TRUE,intercept=FALSE,estimate_prior_variance = FALSE)
    s0.ni = susie_trendfilter(y,0,standardize=TRUE,intercept=FALSE,estimate_prior_variance = FALSE, use_mad=FALSE)

    expect_equal_susie(s0, original.s0)
    expect_equal_susie(s0.ns.ni, original.s0.ns.ni)
    expect_equal_susie(s0.ns, original.s0.ns)
    expect_equal_susie(s0.ni, original.s0.ni)
  })
  #order 1 trend filtering
  with(simulate_tf(1), {
    suppressWarnings({
      original.s1 <- susie(X,y,estimate_prior_variance = FALSE)
      s1 <- suppressWarnings(susie_trendfilter(y,1,estimate_prior_variance = FALSE, standardize = TRUE, use_mad=FALSE))
      original.s1.ns.ni <- susie(X,y,standardize=FALSE,intercept=FALSE,estimate_prior_variance = FALSE)
      s1.ns.ni <- suppressWarnings(susie_trendfilter(y,1,standardize=FALSE,intercept=FALSE,estimate_prior_variance = FALSE, use_mad=FALSE))
      original.s1.ns <- susie(X,y,standardize=FALSE,intercept=TRUE,estimate_prior_variance = FALSE)
      s1.ns <- suppressWarnings(susie_trendfilter(y,1,standardize=FALSE,intercept=TRUE,estimate_prior_variance = FALSE, use_mad=FALSE))
      original.s1.ni <- susie(X,y,standardize=TRUE,intercept=FALSE,estimate_prior_variance = FALSE)
      s1.ni <- suppressWarnings(susie_trendfilter(y,1,standardize=TRUE,intercept=FALSE,estimate_prior_variance = FALSE, use_mad=FALSE))
    })

    expect_equal_susie(s1, original.s1)
    expect_equal_susie(s1.ns.ni, original.s1.ns.ni)
    expect_equal_susie(s1.ns, original.s1.ns)
    expect_equal_susie(s1.ni, original.s1.ni)
  })
  #order 2 trend filtering
  with(simulate_tf(2), {
    original.s2 = susie(X,y,estimate_prior_variance = FALSE)
    s2 = suppressWarnings(susie_trendfilter(y,2,estimate_prior_variance = FALSE, standardize = TRUE, use_mad=FALSE))
    original.s2.ns.ni = susie(X,y,standardize=FALSE,intercept=FALSE,estimate_prior_variance = FALSE)
    s2.ns.ni = suppressWarnings(susie_trendfilter(y,2,standardize=FALSE,intercept=FALSE,estimate_prior_variance = FALSE, use_mad=FALSE))
    original.s2.ns = susie(X,y,standardize=FALSE,intercept=TRUE,estimate_prior_variance = FALSE)
    s2.ns = suppressWarnings(susie_trendfilter(y,2,standardize=FALSE,intercept=TRUE,estimate_prior_variance = FALSE, use_mad=FALSE))
    original.s2.ni = susie(X,y,standardize=TRUE,intercept=FALSE,estimate_prior_variance = FALSE)
    s2.ni = suppressWarnings(susie_trendfilter(y,2,standardize=TRUE,intercept=FALSE,estimate_prior_variance = FALSE, use_mad=FALSE))

    expect_equal_susie(s2, original.s2)
    expect_equal_susie(s2.ns.ni, original.s2.ns.ni)
    expect_equal_susie(s2.ns, original.s2.ns)
    expect_equal_susie(s2.ni, original.s2.ni)
  })
})

#check if susie_trendfilter can fit different arguments
test_that("susie_trendfilter interface",{
  with(simulate_tf(0), {
    n = length(y)
    expect_error(susie_trendfilter(y,0,standardize = TRUE, estimate_prior_variance=TRUE, use_mad=FALSE), NA)
    expect_error(susie_trendfilter(y,0,standardize = TRUE, null_weight=1/(n+1),estimate_prior_variance = FALSE, use_mad=FALSE), NA)
    expect_error(susie_trendfilter(y,0,standardize = TRUE, estimate_prior_variance=TRUE, null_weight=1/(n+1), use_mad=FALSE), NA)
    expect_error(susie_trendfilter(y,0,standardize=FALSE,intercept=FALSE,estimate_prior_variance = FALSE, use_mad=FALSE), NA)
    expect_error(susie_trendfilter(y,0,standardize=FALSE,intercept=TRUE,estimate_prior_variance = FALSE, use_mad=FALSE), NA)
    expect_error(susie_trendfilter(y,0,standardize=TRUE,intercept=FALSE,estimate_prior_variance = FALSE, use_mad=FALSE), NA)
    expect_error(suppressWarnings(susie_trendfilter(y,0,standardize = TRUE, compute_univariate_zscore =TRUE,estimate_prior_variance = FALSE, use_mad=FALSE)), NA)
  })
})
