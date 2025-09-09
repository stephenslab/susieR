# This R file contains my very simple implementations of the single
# effect regression (SER) models.

# Compute the log-normalizing factor for the IG(a,b) distribution.
inv_gamma_factor <- function (a, b)
  a*log(b) - lgamma(a)

# Safely compute probabilities from log-weights.
normalizelogweights <- function (logw) {
  x <- max(logw)
  w <- exp(logw - x)
  return(w/sum(w))
}

# For each column of X, compute the least-squares estimate of the
# regression coefficient (bhat) and its variance (shat). Here, s is
# the residual variance.
least_squares_coefs <- function (X, y, s) {
  xx   <- colSums(X^2)
  xy   <- drop(crossprod(X,y))
  return(list(bhat = xy/xx,shat = s/xx))
}

# Compute the (log) Bayes factors given the least-squares estimates
# (bhat) and their corresponding variances (shat). Here, s0 is the prior
# variance of the regression coefficient.
compute_bayes_factors <- function (bhat, shat, s0)
  dnorm(bhat,0,sqrt(s0 + shat),log = TRUE) -
  dnorm(bhat,0,sqrt(shat),log = TRUE)

# Compute the (log) Bayes factors and additional statistics under the
# normal-inverse-gamma (NIG) prior given the least-squares estimates
# (bhat) and their corresponding variances (shat). Additional outputs
# include the compute the posterior moments, mean (b1), second moment
# (b2) and variance (s1), and the posterior mode of the residual
# variance (rv). Here, s0 is the prior variance of the regression
# coefficient, n is the sample size, yy = sum(y^2), sxy = cor(x,y),
# and the residual variance has an IG(a0/2,b0/2) prior. I perform all
# these computations at once because it is just easier that way.
compute_stats_NIG <- function (n, xx, xy, yy, sxy, s0, a0, b0) {
  r0   <- s0/(s0 + 1/xx)
  rss  <- yy*(1 - r0*sxy^2)
  a1   <- a0 + n
  b1   <- b0 + rss
  lbf  <- -(log(1 + s0*xx) + a1*log(b1/(b0 + yy)))/2
  # Note that the line above is the same as:
  # lbf = (log(1 - r0) - a1*log(b1/(b0 + yy)))/2
  bhat <- xy/xx
  mu1  <- r0*bhat
  s1   <- b1/(a1 - 2)*r0/xx
  rv   <- (b1/2)/(a1/2 - 1)
  return(list(lbf = lbf,
              b1  = mu1,
              b2  = s1 + mu1^2,
              s1  = s1,
              rv  = rv))
}

# Compute the posterior moments: mean (b1), second moment (b2),
# variance (s1). The inputs are the least-squares estimates (bhat),
# their variances (shat), and the prior variance (s0).
compute_posterior_moments <- function (bhat, shat, s0) {
  s1 <- 1/(1/s0 + 1/shat)
  b1 <- s1*bhat/shat
  b2 <- b1^2 + s1
  return(list(b1 = b1,b2 = b2,s1 = s1))
}

# Compute the log-likelihood under the "null" model.
compute_null_loglik <- function (y, s)
  sum(dnorm(y,0,sqrt(s),log = TRUE))

# Compute the log-likelihood under the "null" model under the
# normal-inverse-gamma (NIG) prior.
compute_null_loglik_NIG <- function (n, yy, a0, b0) {
  return(-n*log(2*pi)/2 +
           inv_gamma_factor(a0/2,b0/2) -
           inv_gamma_factor((a0 + n)/2,(b0 + yy)/2))
}

# Compute the expected residual sum of squares. The inputs are the
# data (X, y); the posterior inclusion probabilities (pip), and the
# posterior second moments, the mean (b1) and the second moment (b2).
compute_erss <- function (X, y, pip, b1, b2) {
  xx  <- colSums(X^2)
  bb  <- pip * b1
  bb2 <- pip * b2
  xb  <- drop(X %*% bb)
  r   <- drop(y - xb)
  return(sum(r^2) - sum(xb^2) + sum(xx * bb2))
}

# Compute the ELBO (which is also the marginal log-likelihood since
# the ELBO is exact). The inputs are: ll0, the log-likelihood under
# the "null" model; and lbf, the log-Bayes factors.
compute_elbo <- function (ll0, lbf) {
  x <- max(lbf)
  bf <- exp(lbf - x)
  return(ll0 + x + log(mean(bf)))
}

# The EM update for the prior variance is the posterior second moment
# of the coefficient. The inputs are the posterior inclusion
# probabilities (pip) and the posterior second moments (b2).
update_prior_variance <- function (pip, b2)
  sum(pip * b2)

# Update the prior variance under the NIG prior by maximizing the
# marginal log-likelihood.
update_prior_variance_NIG <- function (n, xx, xy, yy, sxy, s0, a0, b0) {
  ll0 <- compute_null_loglik_NIG(n,yy,a0,b0)
  f <- function (s0) {
    lbf <- compute_stats_NIG(n,xx,xy,yy,sxy,s0,a0,b0)$lbf
    return(-compute_elbo(ll0,lbf))
  }
  out <- optim(s0,f,method = "Brent",lower = 0,upper = 10)
  return(out$par)
}

# This is an experimental implementation to check that my EM update is
# correct.
update_prior_variance_NIG2 <- function (n, xx, xy, yy, sxy, pip,
                                        s0, a0, b0) {
  p    <- length(pip)
  r0   <- s0/(s0 + 1/xx)
  rss  <- yy*(1 - r0*sxy^2)
  a1   <- a0 + n
  b1   <- b0 + rss
  bhat <- xy/xx
  mu1  <- r0*bhat
  s1   <- r0/xx
  ns  <- 100
  out <- 0
  for (i in 1:ns) {
    # ss = rinvgamma(p,a1/2,b1/2)
    ss <- 1/rgamma(p,a1/2,b1/2)
    bb  <- rnorm(p,mu1/sqrt(ss),sqrt(s1))
    out <- out + sum(pip * bb^2)
  }
  return(out/ns)
}

# My simple implementation of EM for the SER models.
# For the "normal-inverse-gamma" (NIG) prior, in which the residual
# variance has an IG(a0/2,b0/2) prior. For this prior, note that
# prior mean = (b0/2)/(a0/2 - 1),
# prior mode = (b0/2)/((a0/2) + 1)
simple_susie <- function (X, y, s0 = 1, s = 1, numiter = 10,
                          prior = c("normal","NIG"),
                          a0 = 0.01, b0 = 0.01, verbose = TRUE) {
  prior <- match.arg(prior)

  # Center X and y which is needed to account for an intercept.
  X <- scale(X,center = TRUE,scale = FALSE)
  y <- drop(scale(y,center = TRUE,scale = FALSE))

  # Iterate the EM updates.
  n <- length(y)
  elbo <- rep(0,numiter)
  for (i in 1:numiter) {
    if (verbose)
      cat("(",i,") s=",s,", s0=",s0,"\n",sep = "")

    # Compute the Bayes factors and posterior inclusion probabilities
    # (PIPs).
    if (prior == "normal") {
      dat  <- least_squares_coefs(X,y,s)
      bhat <- dat$bhat
      shat <- dat$shat
      lbf  <- compute_bayes_factors(bhat,shat,s0)
    } else if (prior == "NIG") {
      dat <- list(xx  = colSums(X^2),
                  xy  = drop(crossprod(X,y)),
                  yy  = sum(y^2),
                  sxy = drop(cor(X,y)))
      lbf <- with(dat,compute_stats_NIG(n,xx,xy,yy,sxy,s0,a0,b0))$lbf
    }
    pip <- normalizelogweights(lbf)

    # Compute the ELBO.
    if (prior == "normal") {
      ll0 <- compute_null_loglik(y,s)
    } else if (prior == "NIG") {
      ll0 <- compute_null_loglik_NIG(n,dat$yy,a0,b0)
    }
    elbo[i] <- compute_elbo(ll0,lbf)

    # Compute the posterior means (b1), posterior variances (s1) and
    # posterior second moments (b2) of the coefficients.
    if (prior == "normal") {
      out <- compute_posterior_moments(bhat,shat,s0)
    } else if (prior == "NIG") {
      out <- with(dat,compute_stats_NIG(n,xx,xy,yy,sxy,s0,a0,b0))
      rv <- out$rv
    }
    b1 <- out$b1
    b2 <- out$b2
    s1 <- out$s1

    if (prior == "normal") {

      # The EM update for the residual variance is ERSS/n, where "ERSS" is
      # the expected residual sum of squares.
      erss <- compute_erss(X,y,pip,b1,b2)
      s <- erss/n
    } else if (prior == "NIG") {

      # Obtain the posterior mean of the residual variance parameter.
      s <- sum(pip * rv)
    }

    # The EM update for the prior variance is the posterior second
    # moment of the coefficient.
    if (prior == "normal") {
      s0 <- update_prior_variance(pip,b2)
    } else if (prior == "NIG") {

      # This update uses optimize().
      # s0 <- with(dat,update_prior_variance_NIG(n,xx,xy,yy,sxy,s0,a0,b0))

      # This is a stochastic EM update.
      s0 <- with(dat,update_prior_variance_NIG2(n,xx,xy,yy,sxy,pip,s0,a0,b0))
    }
  }

  return(list(s = s,s0 = s0,pip = pip,b1 = b1,b2 = b2,s1 = s1,elbo = elbo))
}



library(susieR)

data(data_small)
y <- data_small$y
X <- data_small$X
dim(X)
tt0= simple_susie(X,y)
plot(tt0$elbo, col="green", main="matching ELBO")
tt= simple_susie(X,y, prior="NIG")
points(tt$elbo, col="blue")

res0=susie(X,y, L=1)
points(res0$elbo, pch=19, col="green")
res1=susie(X,y, L=1, small=TRUE)
points(res1$elbo, pch=19, col="blue")


plot( tt0$b1, res0$mu, ) # susie Gaussian effect estimate
#vs susie simple estimate with Gaussian prior
abline(a=0,b=1)
points( tt$b1, res1$mu, pch=19, col="lightgreen")
# susie NIG effect estimate
#vs susie simple estimate
plot (res1$pip)
abline(a=0,b=1)


#problem when L>1

resL2=susie(X,y, L=2, small=TRUE)
plot(resL2$pip)
resL2$sets$cs
resL6=susie(X,y, L=6, small=TRUE)
plot(resL6$pip)
resL6$sets$cs
plot(resL6$alpha[1,])

