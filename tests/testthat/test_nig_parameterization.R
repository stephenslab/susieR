test_that("NIG agrees with numerical integration of the joint density", {
  set.seed(71)
  X = matrix(rnorm(27),9,3)
  y = X[,1] + rnorm(9)
  V = 0.37
  kappa = 1.4
  lambda = 0.8
  prior = c(0.1,0.3,0.6)
  # Independently integrate b after integrating sigma2 in the joint density.
  # This deliberately uses residual vectors, not the implementation's formulas.
  a = (kappa + length(y) + 1)/2
  constant = -(length(y)+1)/2*log(2*pi) - log(V)/2 +
    kappa/2*log(lambda/2) - lgamma(kappa/2) + lgamma(a)
  reference = sapply(seq_len(ncol(X)),function(j) {
    rate = function(b) (lambda + sum((y-X[,j]*b)^2) + b^2/V)/2
    density = Vectorize(function(b) exp(constant-a*log(rate(b))))
    integral = function(f) integrate(f,-Inf,Inf,rel.tol=1e-10,
                                    abs.tol=1e-20)$value
    z = integral(density)
    c(z=z,mu=integral(function(b) b*density(b))/z,
      mu2=integral(function(b) b^2*density(b))/z,
      em=integral(Vectorize(function(b) b^2*a/rate(b)*density(b)))/z)
  })
  weights = prior*reference["z",]
  alpha = weights/sum(weights)
  fit = single_effect_regression_nig(drop(crossprod(X,y)),colSums(X^2),
           sum(y^2),9,V,prior,alpha0=kappa,beta0=lambda)
  expect_equal(fit$alpha,alpha,tolerance=1e-8)
  expect_equal(fit$mu,unname(reference["mu",]),tolerance=1e-8)
  expect_equal(fit$mu2,unname(reference["mu2",]),tolerance=1e-8)
  expect_equal(fit$loglik,log(sum(weights)),tolerance=1e-8)
  em = single_effect_regression_nig(drop(crossprod(X,y)),colSums(X^2),
          sum(y^2),9,V,prior,"EM",alpha0=kappa,beta0=lambda)
  expect_equal(em$V,sum(alpha*reference["em",]),tolerance=1e-8)
  expect_gte(em$loglik,fit$loglik-1e-10)
  refit = single_effect_regression_nig(drop(crossprod(X,y)),colSums(X^2),
            sum(y^2),9,em$V,prior,alpha0=kappa,beta0=lambda)
  expect_equal(em$mu2,refit$mu2)
})

test_that("NIG scale and standardization invariances hold for dense and sparse data", {
  set.seed(4)
  X = matrix(rnorm(600),30,20)
  X[,2] = 0.8*X[,1] + 0.6*X[,2]
  y = X[,1]+rnorm(30)
  for (L in c(1,3)) {
    f = function(X,y,...) suppressWarnings(susie(X,y,L=L,small=TRUE,
                    max_iter=300,min_abs_corr=0,estimate_prior_variance=FALSE,...))
    a = f(X,y)
    for (scale in c(0.01,10,-2)) {
      b = f(X,scale*y)
      expect_equal(a$alpha,b$alpha,tolerance=1e-9)
      expect_equal(a$V,b$V)
      expect_equal(a$mu,b$mu/scale,tolerance=1e-9)
      expect_equal(a$mu2,b$mu2/scale^2,tolerance=1e-9)
      expect_equal(a$sigma2,b$sigma2/scale^2,tolerance=1e-9)
    }
    b = f(sweep(X,2,seq(0.1,5,length.out=20),"*"),y)
    expect_equal(a$alpha,b$alpha,tolerance=1e-9)
    expect_equal(a$mu,b$mu,tolerance=1e-9)
    b = f(as(X,"CsparseMatrix"),y)
    expect_equal(a$alpha,b$alpha,tolerance=1e-9)
    expect_equal(a$mu2,b$mu2,tolerance=1e-9)
  }
  a = susie(X,y,small=TRUE,L=1,alpha0=2,beta0=3)
  b = susie(X,10*y,small=TRUE,L=1,alpha0=2,beta0=300)
  expect_equal(a$alpha,b$alpha,tolerance=1e-8)
  expect_equal(a$V,b$V,tolerance=1e-8)
})

test_that("NIG honors zero weights, null effects, and zero-information columns", {
  fit = single_effect_regression_nig(c(1,0,0),c(2,0,0),5,6,0,
                                     c(0.2,0,0.8),alpha0=1,beta0=1)
  expect_equal(fit$alpha,c(0.2,0,0.8))
  expect_equal(fit$lbf,c(0,0,0))
  expect_equal(fit$mu2,c(0,0,0))
  fit = single_effect_regression_nig(c(1,0,0),c(2,0,0),5,6,0.3,
                                     c(0.2,0,0.8),alpha0=1,beta0=1)
  expect_identical(fit$alpha[2],0)
  expect_equal(fit$lbf[2:3],c(0,0))
  expect_true(all(is.finite(fit$mu2)))
  set.seed(8)
  X = matrix(rnorm(150),30,5)
  f = susie(X,rnorm(30),L=1,small=TRUE,null_weight=0.8,
            scaled_prior_variance=0,estimate_prior_variance=FALSE)
  expect_equal(f$alpha[1,6],0.8)
  expect_true(all(f$pip == 0))
})

test_that("NIG honors intercept and effective residual dimension", {
  set.seed(10)
  X = matrix(rnorm(120,2),30,4)
  y = 3+X[,1]+rnorm(30)
  f = susie(X,y,L=1,small=TRUE,intercept=FALSE,standardize=FALSE,
             scaled_prior_variance=0.4,estimate_prior_variance=FALSE,
             alpha0=2,beta0=3)
  r = single_effect_regression_nig(drop(crossprod(X,y)),colSums(X^2),
                                   sum(y^2),30,0.4,alpha0=2,beta0=3)
  expect_equal(drop(f$alpha),r$alpha)
  expect_equal(drop(f$mu2),r$mu2)
  f = susie(X,y,L=1,small=TRUE,alpha0=2,beta0=3)
  expect_equal(f$nig$shape,(29+2)/2)
  shifted = susie(X+5,y-7,L=1,small=TRUE,alpha0=2,beta0=3)
  expect_equal(f$alpha,shifted$alpha,tolerance=1e-8)
})

test_that("NIG optim and EM fit the NIG likelihood and monitor all parameters", {
  set.seed(7)
  x = matrix(rnorm(30),30,1)
  y = 2*x[,1]+rnorm(30)
  args = list(xty=drop(crossprod(x,y)),d=sum(x^2),yty=sum(y^2),n=30,
              V=0.2,alpha0=2,beta0=3)
  fit = do.call(single_effect_regression_nig,c(args,list(optimize_V="optim")))
  # Closed-form ML estimate for a single covariate (no mixture).
  explained = args$xty^2/args$d
  expected = max(0,((args$n+args$alpha0-1)*explained/
                   (args$beta0+args$yty-explained)-1)/args$d)
  expect_equal(fit$V,expected,tolerance=1e-4)
  e = susie(x,y,L=1,small=TRUE,scaled_prior_variance=1e-4,
            estimate_prior_method="EM",max_iter=1000,tol_small=1e-7)
  expect_true(e$converged)
  expect_gt(e$niter,3)
  expect_true(all(diff(e$elbo) >= -1e-8))
  expect_lt(tail(e$cv_criterion,1),1e-7)
  one = suppressWarnings(susie(x,y,L=1,small=TRUE,max_iter=1))
  expect_length(one$cv_criterion,1)
  expect_false(one$converged)
  fixed = susie(x,y,L=1,small=TRUE,scaled_prior_variance=2,
                estimate_prior_variance=FALSE)
  expect_equal(fixed$V,2)
})

test_that("NIG exposes evidence only where defined and uses Student posterior summaries", {
  set.seed(9)
  X = matrix(rnorm(150),30,5)
  y = X[,1]+rnorm(30)
  f = susie(X,y,L=1,small=TRUE)
  expect_equal(susie_get_objective(f),f$nig$loglik[1])
  expect_true(all(is.na(f$KL)))
  expected = sum(f$alpha*stats::pt(-abs(f$mu)/sqrt(f$nig$scale2),
                                   df=f$nig$shape))
  expect_equal(susie_get_lfsr(f),expected)
  set.seed(11)
  samples = susie_get_posterior_samples(f,10000)$b
  expect_equal(rowMeans(samples),susie_get_posterior_mean(f),tolerance=0.03)
  g = susie(X,y,L=2,small=TRUE,estimate_prior_variance=FALSE)
  expect_true(all(is.na(g$elbo)))
  expect_error(susie_get_objective(g),"No joint evidence")
  g = susie(X,y,L=1,small=TRUE,alpha0=0,beta0=0)
  expect_true(all(is.finite(g$alpha)))
  expect_error(susie_get_objective(g),"No joint evidence")
  expect_error(susie(X,y,small=TRUE,refine=TRUE),"no joint ELBO")
  expect_error(susie(X,y,small=TRUE,alpha0=-1),"nonnegative")
  expect_error(susie(X[1:3,],y[1:3],small=TRUE,alpha0=0,beta0=0),
                "second moments")
})
