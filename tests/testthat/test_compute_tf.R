context("test_compute_tf.R")

expect_equal_compute_tf = function(order){
  simulate_tf(order)
  scaled_X = safe_colScale(X)
  set.seed(1)
  b = rnorm(length(y))
  Xb = compute_tf_Xb(order, b=b) 
  expect_equal(as.matrix(Xb), X%*%b)
  
  Xty = compute_tf_Xty(order,y=y)
  expect_equal(as.matrix(Xty), t(X)%*%y)
  
  d = compute_tf_d(order, n=length(y))
  expect_equal(d, colSums(X*X))
  
  cm = compute_tf_cm(order, n=length(y))
  expect_equal(cm, attr(scaled_X, 'scaled:center'))
  
  csd = compute_tf_csd(order, n=length(y))
  expect_equal(csd, attr(scaled_X, 'scaled:scale'))
  
  std_d = compute_tf_std_d(order, n=length(y))
  expect_equal(std_d, colSums(scaled_X*scaled_X))
  
  set.seed(2)
  n = length(y)
  Eb2 = rnorm(n)
  M <- Matrix::sparseMatrix(i=NULL,j=NULL,dims=c(n,n))
  attr(M, "matrix.type") = "tfmatrix"
  attr(M, "order") = order
  attr(M, "d") <- compute_tf_std_d(order, n)
  attr(M, "scaled:center") <- compute_tf_cm(order, n)
  attr(M, "scaled:scale") <- compute_tf_csd(order, n)
  X2tEb2 = compute_tf_X2tEb2(M, Eb2=Eb2)
  expect_equal(X2tEb2, sum(t(scaled_X * scaled_X) * Eb2))
  
  attr(M, "d") <- compute_tf_d(order, n)
  attr(M, "scaled:center") <- rep(0,n)
  attr(M, "scaled:scale") <- rep(1,n)
  X2tEb2 = compute_tf_X2tEb2(M, Eb2=Eb2)
  expect_equal(X2tEb2, sum(t(X * X) * Eb2))
}

test_that("computation about trend filtering",{
  expect_equal_compute_tf(order=0)
  expect_equal_compute_tf(order=1)
  expect_equal_compute_tf(order=2)
})
