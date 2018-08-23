create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)     
}

test_that("sparse version safe_colScale",{
  set.seed(1)
  n = 1000
  p = 10000
  beta = rep(0,p)
  beta[1]    = 10 
  beta[300]  = 10
  beta[400]  = 10
  beta[1000] = 10
  X.dense = create_sparsity_mat(0.99,n,p)
  X.sparse = as(X.dense,'dgCMatrix')
  
  dense.res = susieR:::safe_colScale(X.dense)
  sparse.res = susieR:::safe_colScale(X.sparse)
  
  dense.faceX = dense.res
  attributes(dense.faceX) = NULL
  dense.scaledX = attr(dense.res, 'scaled.X')
  attributes(dense.scaledX) = NULL
  
  sparse.faceX = sparse.res
  attributes(sparse.faceX) = NULL
  attributes(X.sparse) = NULL
  sparse.scaledX = attr(sparse.res, 'scaled.X')
  attributes(sparse.scaledX) = NULL
  
  expect_equal(sparse.faceX, X.sparse)
  expect_equal(sparse.scaledX, dense.faceX)
  expect_equal(attr(dense.res, 'scaled:center'), attr(sparse.res, 'scaled:center'))
  expect_equal(attr(dense.res, 'scaled:scale'), attr(sparse.res, 'scaled:scale'))
})
