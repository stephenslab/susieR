context("test_compute_tf.R")

test_that("computation about trend filtering",{
  #order 0 trend filtering
  simulate_tf(order=0)
  scaled_X = safe_colScale(X)
  set.seed(1)
  b = rnorm(length(y))
  Xb = compute_tf_Xb(order=0, b=b) 
  expect_equal(as.matrix(Xb), X%*%b)
  
  Xty = compute_tf_Xty(order=0,y=y)
  expect_equal(as.matrix(Xty), t(X)%*%y)
  
  d = compute_tf_d(order=0, n=length(y))
  expect_equal(d, colSums(X*X))
  
  cm = compute_tf_cm(order=0, n=length(y))
  expect_equal(cm, attr(scaled_X, 'scaled:center'))
  
  csd = compute_tf_csd(order=0, n=length(y))
  expect_equal(csd, attr(scaled_X, 'scaled:scale'))
  
  std_d = compute_tf_std_d(order=0, n=length(y))
  expect_equal(std_d, colSums(scaled_X*scaled_X))
  
  set.seed(2)
  Eb2 = rnorm(length(y))
  X2tEb2 = compute_tf_X2tEb2(order=0, n=length(y), Eb2=Eb2)
  expect_equal(X2tEb2, sum(t(scaled_X * scaled_X) * Eb2))
  
  
  #order 1 trend filtering
  simulate_tf(order=1)
  scaled_X = safe_colScale(X)
  set.seed(1)
  b = rnorm(length(y))
  Xb = compute_tf_Xb(order=1, b=b) 
  expect_equal(as.matrix(Xb), X%*%b)
  
  Xty = compute_tf_Xty(order=1,y=y)
  expect_equal(as.matrix(Xty), t(X)%*%y)
  
  d = compute_tf_d(order=1, n=length(y))
  expect_equal(d, colSums(X*X))
  
  cm = compute_tf_cm(order=1, n=length(y))
  expect_equal(cm, attr(scaled_X, 'scaled:center'))
  
  csd = compute_tf_csd(order=1, n=length(y))
  expect_equal(csd, attr(scaled_X, 'scaled:scale'))
  
  std_d = compute_tf_std_d(order=1, n=length(y))
  expect_equal(std_d, colSums(scaled_X*scaled_X))
  
  set.seed(2)
  Eb2 = rnorm(length(y))
  X2tEb2 = compute_tf_X2tEb2(order=1, n=length(y), Eb2=Eb2)
  expect_equal(X2tEb2, sum(t(scaled_X * scaled_X) * Eb2))
  
  
  #order 2 trend filtering
  simulate_tf(order=2)
  scaled_X = safe_colScale(X)
  set.seed(1)
  b = rnorm(length(y))
  Xb = compute_tf_Xb(order=2, b=b) 
  expect_equal(as.matrix(Xb), X%*%b)
  
  Xty = compute_tf_Xty(order=2,y=y)
  expect_equal(as.matrix(Xty), t(X)%*%y)
  
  d = compute_tf_d(order=2, n=length(y))
  expect_equal(d, colSums(X*X))
  
  cm = compute_tf_cm(order=2, n=length(y))
  expect_equal(cm, attr(scaled_X, 'scaled:center'))
  
  csd = compute_tf_csd(order=2, n=length(y))
  expect_equal(csd, attr(scaled_X, 'scaled:scale'))
  
  std_d = compute_tf_std_d(order=2, n=length(y))
  expect_equal(std_d, colSums(scaled_X*scaled_X))
  
  set.seed(2)
  Eb2 = rnorm(length(y))
  X2tEb2 = compute_tf_X2tEb2(order=2, n=length(y), Eb2=Eb2)
  expect_equal(X2tEb2, sum(t(scaled_X * scaled_X) * Eb2))

})
