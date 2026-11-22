context("Corrected NIG single-effect regression")

test_that("NIG reproduces the corrected simulation fits", {
  reference <- dget(test_path("fixtures", "nig_simulation.R"))
  set.seed(1122)
  X <- matrix(rnorm(24 * 12), 24, 12)
  X[, 2] <- 0.9 * X[, 1] + 0.1 * X[, 2]
  y <- 2 * X[, 1] - 1.5 * X[, 4] + rnorm(24)
  cases <- list(defaults_l1 = list(L = 1), defaults_l3 = list(L = 3),
                em = list(L = 3, estimate_prior_method = "EM"),
                custom = list(L = 3, alpha0 = 4, beta0 = 3,
                              scaled_prior_variance = 0.7))
  for (name in names(cases)) {
    fit <- do.call(susie, c(list(X = X, y = y, prior = "nig",
                                 coverage = NULL), cases[[name]]))
    for (field in names(reference[[name]]))
      expect_equal(unname(fit[[field]]), unname(reference[[name]][[field]]),
                   tolerance = 1e-7, info = paste(name, field))
  }
})

test_that("NIG defaults, explicit settings and verbose hints are respected", {
  set.seed(1122)
  X <- matrix(rnorm(20 * 30), 20, 30)
  y <- X[, 1] + rnorm(20)
  init <- function(...) susie(X, y, prior = "nig", init_only = TRUE, ...)$params

  expect_no_message(p <- init())
  expect_equal(p$max_iter, 300)
  expect_equal(p$alpha0, 2)
  expect_equal(p$beta0, 2 * var(y))
  expect_equal(p$n_effective, 19)
  expect_equal(p$tol, 1e-4)
  expect_equal(p$estimate_prior_method, "optim")
  expect_equal(init(max_iter = 17)$max_iter, 17)
  expect_equal(init(alpha0 = 4)$beta0, 4 * var(y))
  expect_equal(init(beta0 = 3)$beta0, 3)
  expect_equal(init(intercept = FALSE)$beta0, 2 * mean(y^2))
  expect_equal(init(estimate_prior_variance = FALSE)$estimate_prior_method, "none")
  expect_message(init(verbose = TRUE), "GIBSS.*Default parameterization")
  msgs <- character()
  withCallingHandlers(invisible(init(verbose = TRUE, alpha0 = 4, beta0 = 3)),
                      message = function(m) {
                        msgs <<- c(msgs, conditionMessage(m))
                        invokeRestart("muffleMessage")
                      })
  expect_true(any(grepl("GIBSS", msgs)))
  expect_false(any(grepl("Default parameterization", msgs)))
  expect_no_message(init(verbose = FALSE, alpha0 = 4, beta0 = 3))
  expect_error(init(estimate_residual_method = "MLE"), "different priors")
  expect_error(init(refine = TRUE), "no joint ELBO")
  expect_equal(susie(X, y, init_only = TRUE)$params$max_iter, 100)
  expect_equal(susie(X, y, estimate_residual_method = "NIG",
                     init_only = TRUE)$params, p)

  # Exercise the actual loop limit, including explicitly supplied limits.
  for (limit in c(300, 5)) {
    args <- list(X = X, y = y, L = 1, prior = "nig", tol = 0,
                 estimate_prior_variance = FALSE, coverage = NULL)
    if (limit == 5) args$max_iter <- limit
    expect_warning(fit <- suppressMessages(do.call(susie, args)), "did not converge")
    expect_equal(fit$niter, limit)
  }
})

test_that("NIG Bayes factors agree with numerical marginalization", {
  x <- c(-2, -1, 0, 1, 2)
  y <- c(-1.5, -0.8, 0.3, 0.6, 1.4)
  y <- y - mean(y)
  n <- length(y) - 1
  xx <- sum(x^2); xy <- sum(x * y); yy <- sum(y^2)
  V <- 0.7; alpha0 <- 2; beta0 <- 2 * var(y)
  evidence <- function(v) integrate(function(precision) {
    omega <- v / (1 + v * xx)
    exp(n / 2 * log(precision / (2 * pi)) - log1p(v * xx) / 2 -
          precision * (yy - omega * xy^2) / 2 +
          dgamma(precision, alpha0 / 2, rate = beta0 / 2, log = TRUE))
  }, 0, Inf, rel.tol = 1e-10)$value
  expect_equal(compute_lbf_NIG(n, xx, xy, yy, NULL, V, alpha0, beta0),
               log(evidence(V) / evidence(0)), tolerance = 1e-8)

  # A zero-information column retains its prior variance, and V=0 is exact.
  post <- compute_posterior_moments_NIG(n, c(xx, 0), c(xy, 0), yy,
                                        NULL, V, alpha0, beta0)
  expect_equal(post$post_mean[2], 0)
  expect_equal(post$post_mean2[2], V * (beta0 + yy) / (n + alpha0 - 2))
  expect_equal(compute_lbf_NIG(n, c(xx, 0), c(xy, 0), yy, NULL,
                              0, alpha0, beta0), c(0, 0))
})

test_that("individual, sparse and centered sufficient-statistic NIG fits agree", {
  set.seed(1123)
  X <- matrix(rnorm(20 * 30), 20, 30)
  X <- cbind(X, 0)
  y <- 2 * X[, 1] - X[, 4] + rnorm(20)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  fit <- susie(X, y, L = 3, prior = "nig", coverage = NULL)
  sparse <- susie(as(X, "dgCMatrix"), y, L = 3, prior = "nig", coverage = NULL)
  ss <- susie_ss(crossprod(Xc), drop(crossprod(Xc, yc)), sum(yc^2), n = 20,
                 L = 3, prior = "nig", coverage = NULL)
  for (field in c("alpha", "mu", "mu2", "V", "lbf", "sigma2", "niter")) {
    expect_equal(unname(sparse[[field]]), unname(fit[[field]]), tolerance = 1e-7)
    expect_equal(unname(ss[[field]]), unname(fit[[field]]), tolerance = 1e-7)
  }
  expect_true(all(is.na(fit$elbo)))
  expect_identical(fit$V_parameterization, "variance_ratio")
  expect_true(all(is.finite(fit$mu2)))

  # User-specified zero prior weight must remain exactly zero.
  weights <- c(0, rep(1/30, 30))
  weighted <- susie(X, y, L = 1, prior = "nig", prior_weights = weights,
                    null_weight = 0.2, coverage = NULL)
  expect_identical(weighted$alpha[1, 1], 0)
})

test_that("NIG uses a scale-invariant variance ratio and Student t summaries", {
  set.seed(1124)
  X <- matrix(rnorm(12), 12, 1)
  y <- 2 * X[, 1] + rnorm(12)
  fit <- susie(X, y, L = 1, prior = "nig", scaled_prior_variance = 0.3,
               estimate_prior_variance = FALSE, coverage = NULL)
  scaled <- susie(X, 7 * y, L = 1, prior = "nig", scaled_prior_variance = 0.3,
                  estimate_prior_variance = FALSE, coverage = NULL)
  expect_equal(fit$V, 0.3)
  expect_equal(scaled$V, fit$V)
  expect_equal(scaled$alpha, fit$alpha)
  expect_equal(scaled$mu, 7 * fit$mu)
  expect_equal(scaled$mu2, 49 * fit$mu2)
  expect_equal(scaled$sigma2, 49 * fit$sigma2)
  expect_equal(susie_get_lfsr(fit),
               as.numeric(pt(-abs(fit$mu) / sqrt(fit$nig$scale2),
                              df = 2 * fit$nig$shape)))
  expect_true(all(is.finite(fit$elbo)))

  # Stable inclusion probabilities alone must not stop EM (p = 1).
  em <- susie(X, y, L = 1, prior = "nig", estimate_prior_method = "EM",
              coverage = NULL)
  expect_gt(em$niter, 2)
})
