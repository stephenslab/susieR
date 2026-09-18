test_that("the repository entry point is the slider and uses its own engine", {
  expect_identical(environmentName(environment(susieSlide::susie)),"susieSlide")
  expect_identical(environmentName(environment(susieSlide::susie_workhorse)),"susieSlide")
  expect_false("susieR" %in% names(getNamespaceImports("susieSlide")))
  expect_identical(getFromNamespace(".engine","susieSlide")("susie_workhorse"),
                   susieSlide::susie_workhorse)
  set.seed(441)
  X <- genotypes(600,12)
  y <- 1.5*(X[,2]==2)-.9*(X[,9]>=1)+rnorm(600,sd=.4)
  fit <- quiet(susieSlide::susie,X,y,L=2,residual_variance=.16,
    estimate_residual_variance=FALSE,estimate_prior_variance=FALSE,tol=1e-8)
  expect_s3_class(fit,"susie_slide")
  expect_true(fit$converged)
  expect_lt(fit$delta[which.max(fit$alpha[,2]),2],-.8)
  expect_gt(fit$delta[which.max(fit$alpha[,9]),9],.8)
  expect_equal(susieSlide::predict.susie(fit,newx=X),predict(fit,newx=X),tolerance=1e-10)
  expect_equal(susieSlide::coef.susie(fit),coef(fit),tolerance=1e-10)
  expect_equal(summary(fit)$delta,fit$delta_cs)
})

test_that("the copied additive engine agrees with the separately installed upstream", {
  skip_if_not_installed("susieR","0.16.6")
  set.seed(545)
  X <- genotypes(500,20)
  y <- .8*X[,2]-.7*X[,14]+rnorm(500,sd=.5)
  args <- list(X=X,y=y,L=2,max_iter=500,tol=1e-8)
  original <- quiet(do.call,susieR::susie,args)
  copied <- quiet(do.call,susieSlide::susie_additive,args)
  fixed <- quiet(do.call,susieSlide::susie,c(args,list(delta=0)))
  for(field in c("alpha","mu","mu2","lbf_variable","V","sigma2","fitted","pip","elbo")) {
    expect_equal(copied[[field]],original[[field]],tolerance=2e-6,info=field)
    expect_equal(fixed[[field]],original[[field]],tolerance=2e-6,info=field)
  }
})
