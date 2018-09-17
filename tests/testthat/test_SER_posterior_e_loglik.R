test_that("sparse version SER_posterior_e_loglik",{
  original.res = load_data('SER_original_res.rds')
  simulate(sparse=T)
  Eb = rep(1, p)
  Eb2 = rep(1, p)
  s2 = residual_variance
  
  scaledX = susieR:::safe_colScale(X)
  scaledX.sparse = susieR:::safe_colScale(X.sparse)
  
  dense.res = susieR:::SER_posterior_e_loglik(scaledX,y,s2,Eb,Eb2)
  sparse.res = susieR:::SER_posterior_e_loglik(scaledX.sparse,y,s2,Eb,Eb2)
  
  expect_equal(dense.res, original.res)
  expect_equal(sparse.res, original.res)
  
})
