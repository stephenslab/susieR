context("test_tfg.R")

#check the tfg  matrix multiplication match the original tf when data are used  as breaks
test_that("tfg matches tf",{
  set.seed(1)
  n <- 100
  t = rnorm(n)
  X_tf = make_tf_matrix(n,0)
  X_tfg =  make_tfg_matrix(sort(t),sort(t),0)
  b <- rnorm(n)
  res_tfg = compute_tfg_Xb(X_tfg,b)
  res_tg = compute_tf_Xb(0,b)
  expect_equal(res_tfg,res_tg)

  X_tfg2 =  make_tfg_matrix(t,sort(t),0)
  res_tfg2 = compute_tfg_Xb(X_tfg2,b)
  expect_equal(res_tfg2[order(t)],  res_tfg)

  y = rnorm(n)
  res_tfg = compute_tfg_Xty(X_tfg,y)
  res_tg = compute_tf_Xty(0,y)
  expect_equal(res_tfg[1:n],res_tg) # note that res_tfg produces one extra point at end

})
