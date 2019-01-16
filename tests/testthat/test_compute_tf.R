context("test_compute_tf.R")

expect_equal_compute_tf = function(order){
  with(simulate_tf(order), {
  scaled_X = set_X_attributes(X)
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
  })
}

test_that("computation about trend filtering",{
  expect_equal_compute_tf(order=0)
  expect_equal_compute_tf(order=1)
  expect_equal_compute_tf(order=2)
})
