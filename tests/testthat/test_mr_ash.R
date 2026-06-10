context("mr.ash family")

# Drives mr.ash() option branches, coef.mr.ash, predict.mr.ash,
# get.full.posterior, var.n, and remove_covariate (R/univariate_regression.R).
# Exercising mr.ash also runs the C++ in caisa.cpp / mr_ash.h via caisa_cpp.
# devtools::load_all is active so functions are called directly, and mr.ash()
# is always called with verbose = FALSE. Test data comes from mr_ash_data()
# in helper_mr_ash.R. The benign "mixture proportion associated with the
# largest prior variance" warning is wrapped in suppressWarnings() except in
# the test where that warning is the assertion.

# ---- Input validation ------------------------------------------------------

test_that("mr.ash rejects malformed sa2 and unsupported covariates", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 1)

  # A negative prior-variance entry is rejected.
  bad_neg <- d$sa2; bad_neg[2] <- -1
  expect_error(mr.ash(d$X, d$y, sa2 = bad_neg, verbose = FALSE),
               "non-negative")

  # The first prior-variance entry must be exactly 0 (the null component).
  bad_first <- d$sa2; bad_first[1] <- 0.5
  expect_error(mr.ash(d$X, d$y, sa2 = bad_first, verbose = FALSE),
               "sa2\\[1\\] must be 0")

  # Covariates Z are not implemented.
  Z <- matrix(rnorm(d$n * 2), d$n, 2)
  expect_error(mr.ash(d$X, d$y, Z = Z, verbose = FALSE),
               "covariates are not currently fully implemented")
})

# ---- beta.init x standardize -----------------------------------------------

test_that("mr.ash beta.init runs under standardize TRUE and FALSE without mutating the input", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 4)
  for (st in c(TRUE, FALSE)) {
    binit <- rep(0, d$p); binit[1:3] <- 0.5
    fit <- suppressWarnings(
      mr.ash(d$X, d$y, beta.init = binit, standardize = st,
             max.iter = 50, verbose = FALSE))
    expect_s3_class(fit, "mr.ash")
    expect_length(fit$beta, d$p)
    expect_true(all(is.finite(fit$beta)),
                info = sprintf("standardize=%s", st))
    # The caller's beta.init must not be modified in place.
    expect_equal(binit[1], 0.5)
  }
})

# ---- Default initialisation paths ------------------------------------------

test_that("mr.ash runs end to end with all defaults (sigma2/sa2/pi/update.order NULL)", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 6)
  fit <- suppressWarnings(mr.ash(d$X, d$y, max.iter = 100, verbose = FALSE))
  expect_s3_class(fit, "mr.ash")
  expect_length(fit$beta, d$p)
  expect_length(fit$pi, 25)               # default sa2 grid has 25 entries
  expect_true(is.finite(fit$sigma2) && fit$sigma2 > 0)
  expect_true(all(is.finite(fit$beta)))
})

test_that("mr.ash builds Phi from beta.init when pi is NULL", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 7)
  binit <- rep(0, d$p); binit[1:4] <- 1
  # sa2 supplied so S/Phi is well-defined; pi NULL triggers the Phi-from-beta
  # construction.
  fit <- suppressWarnings(
    mr.ash(d$X, d$y, sa2 = d$sa2, pi = NULL, beta.init = binit,
           standardize = FALSE, max.iter = 50, verbose = FALSE))
  expect_s3_class(fit, "mr.ash")
  expect_length(fit$pi, d$K)
  expect_equal(sum(fit$pi), 1, tolerance = 1e-6)
})

test_that("mr.ash keeps supplied mixture proportions when update.pi = FALSE", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 8)
  fit <- mr.ash(d$X, d$y, sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
                update.pi = FALSE, update.sigma2 = FALSE,
                max.iter = 50, verbose = FALSE)
  expect_s3_class(fit, "mr.ash")
  expect_equal(fit$pi, d$pi0, tolerance = 1e-6)
})

# ---- update.order branches -------------------------------------------------

test_that("mr.ash honors numeric update.order; 1:p matches the default order", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 9)
  # update.order = 1:p is the sequential default; results must match exactly.
  fit_uo  <- suppressWarnings(
    mr.ash(d$X, d$y, update.order = 1:d$p, max.iter = 50, verbose = FALSE))
  fit_def <- suppressWarnings(
    mr.ash(d$X, d$y, max.iter = 50, verbose = FALSE))
  expect_s3_class(fit_uo, "mr.ash")
  expect_equal(fit_uo$beta, fit_def$beta, tolerance = 1e-10)

  # A non-trivial permutation also runs and stays finite.
  set.seed(91)
  fit_perm <- suppressWarnings(
    mr.ash(d$X, d$y, update.order = sample(d$p), max.iter = 50,
           verbose = FALSE))
  expect_length(fit_perm$beta, d$p)
  expect_true(all(is.finite(fit_perm$beta)))
})

test_that("mr.ash accepts update.order = 'random'", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 10)
  set.seed(99)
  fit <- suppressWarnings(
    mr.ash(d$X, d$y, update.order = "random", max.iter = 50, verbose = FALSE))
  expect_s3_class(fit, "mr.ash")
  expect_length(fit$beta, d$p)
  expect_true(all(is.finite(fit$beta)))
})

# ---- standardize rescaling + the update.pi warning -------------------------

test_that("mr.ash with standardize = TRUE returns finite coefficients on the original scale", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 11)
  fit <- suppressWarnings(
    mr.ash(d$X, d$y, standardize = TRUE, max.iter = 80, verbose = FALSE))
  expect_s3_class(fit, "mr.ash")
  expect_length(fit$beta, d$p)
  expect_true(all(is.finite(fit$beta)))
})

test_that("mr.ash warns when mass piles on the largest prior-variance component", {
  # Strong effects relative to noise push EM mass onto the top sa2 component
  # (out$pi[K] > 1/K), which triggers the largest-prior-variance warning.
  set.seed(2024)
  n <- 100; p <- 16
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = FALSE)
  beta_true <- rep(0, p); beta_true[1:3] <- c(40, -35, 30)
  y <- c(X %*% beta_true + rnorm(n, sd = 0.1)); y <- y - mean(y)

  expect_warning(
    mr.ash(X, y, update.pi = TRUE, max.iter = 1000, verbose = FALSE),
    "mixture proportion associated with the")
})

# ---- coef / predict --------------------------------------------------------

test_that("coef.mr.ash returns c(intercept, beta) of length p + 1", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 12)
  fit <- suppressWarnings(mr.ash(d$X, d$y, max.iter = 50, verbose = FALSE))
  cf <- coef(fit)
  expect_length(cf, d$p + 1)
  expect_type(cf, "double")
  expect_equal(cf[1], fit$intercept, tolerance = 1e-6)
  expect_equal(cf[-1], fit$beta, tolerance = 1e-6)
})

test_that("predict.mr.ash type = 'coefficients' equals coef and rejects newx", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 13)
  fit <- suppressWarnings(mr.ash(d$X, d$y, max.iter = 50, verbose = FALSE))
  pc <- predict(fit, type = "coefficients")
  expect_equal(pc, coef(fit), tolerance = 1e-6)
  expect_length(pc, d$p + 1)
  expect_error(predict(fit, newx = d$X, type = "coefficients"),
               "Do not supply newx")
})

test_that("predict.mr.ash without newx and type = 'response' returns object$fitted", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 15)
  fit <- suppressWarnings(mr.ash(d$X, d$y, max.iter = 50, verbose = FALSE))
  # mr.ash does not set $fitted, so this returns NULL.
  expect_null(predict(fit, type = "response"))
})

test_that("predict.mr.ash with newx returns intercept + newx %*% beta", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 16)
  fit <- suppressWarnings(mr.ash(d$X, d$y, max.iter = 50, verbose = FALSE))
  set.seed(500)
  newx <- matrix(rnorm(10 * d$p), 10, d$p)
  pr <- predict(fit, newx = newx)
  expect_length(pr, 10)
  expect_true(all(is.finite(pr)))
  expect_equal(pr, drop(fit$intercept + newx %*% coef(fit)[-1]),
               tolerance = 1e-6)
})

test_that("predict.mr.ash errors for a non-intercept covariate Z", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 17)
  fit <- suppressWarnings(mr.ash(d$X, d$y, max.iter = 50, verbose = FALSE))
  # mr.ash never produces Z != 1; patch the fit to hit the covariate guard.
  fit$data$Z <- matrix(2, d$n, 1)
  set.seed(501)
  newx <- matrix(rnorm(5 * d$p), 5, d$p)
  expect_error(predict(fit, newx = newx),
               "not implemented for covariates")
})

# ---- get.full.posterior ----------------------------------------------------

test_that("get.full.posterior returns list(phi, m, s2) with p x K dims", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 18)
  fit <- suppressWarnings(
    mr.ash(d$X, d$y, sa2 = d$sa2, max.iter = 80, verbose = FALSE))
  fp <- get.full.posterior(fit)
  expect_type(fp, "list")
  expect_named(fp, c("phi", "m", "s2"))
  expect_equal(dim(fp$phi), c(d$p, d$K))
  expect_equal(dim(fp$m), c(d$p, d$K))
  expect_equal(dim(fp$s2), c(d$p, d$K))
  expect_equal(rowSums(fp$phi), rep(1, d$p), tolerance = 1e-6)
  expect_true(all(is.finite(fp$m)))
  expect_true(all(fp$s2 >= 0))
})

# ---- var.n -----------------------------------------------------------------

test_that("var.n equals sum((x - mean)^2) / length", {
  set.seed(20)
  x <- rnorm(20)
  expect_equal(var.n(x), sum((x - mean(x))^2) / length(x), tolerance = 1e-12)
})

# ---- remove_covariate (R/univariate_regression.R) --------------------------

test_that("mr.ash exercises remove_covariate over standardize x intercept combos", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 21)
  for (st in c(TRUE, FALSE)) {
    for (ic in c(TRUE, FALSE)) {
      fit <- suppressWarnings(
        mr.ash(d$X, d$y, standardize = st, intercept = ic,
               max.iter = 40, verbose = FALSE))
      info <- sprintf("standardize=%s intercept=%s", st, ic)
      expect_s3_class(fit, "mr.ash")
      expect_length(fit$beta, d$p)
      expect_true(all(is.finite(fit$beta)), info = info)
    }
  }
})

test_that("remove_covariate with intercept = FALSE, Z = NULL takes the early return", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 22)
  out <- remove_covariate(d$X, d$y, Z = NULL,
                          standardize = FALSE, intercept = FALSE)
  expect_named(out, c("X", "y", "Z", "ZtZiZX", "ZtZiZy"))
  expect_null(out$Z)
  expect_equal(out$ZtZiZX, rep(0, d$p))
  expect_equal(out$ZtZiZy, 0)
  expect_equal(out$X, d$X)   # X and y untouched in this branch
})

test_that("remove_covariate with intercept = TRUE, Z = NULL centers y", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 23)
  out <- remove_covariate(d$X, d$y, Z = NULL,
                          standardize = FALSE, intercept = TRUE)
  expect_equal(dim(out$Z), c(d$n, 1))
  expect_equal(mean(out$y), 0, tolerance = 1e-10)
  expect_equal(out$ZtZiZy, mean(d$y), tolerance = 1e-10)
})

test_that("remove_covariate with standardize = TRUE scales X columns to unit sd", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 24)
  out <- remove_covariate(d$X, d$y, Z = NULL,
                          standardize = TRUE, intercept = TRUE)
  expect_equal(apply(out$X, 2, sd), rep(1, d$p), tolerance = 1e-8)
  expect_false(is.null(attr(out$X, "scaled:scale")))
})

test_that("remove_covariate residualizes against a multi-column Z", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 25)
  set.seed(250)
  Z <- cbind(1, rnorm(d$n))   # ncol(Z) > 1 path
  out <- remove_covariate(d$X, d$y, Z = Z,
                          standardize = FALSE, intercept = FALSE)
  # Residualized y and X are orthogonal to Z.
  expect_equal(c(crossprod(Z, out$y)), rep(0, 2), tolerance = 1e-8)
  expect_equal(max(abs(crossprod(Z, out$X))), 0, tolerance = 1e-8)
})

test_that("remove_covariate with intercept = TRUE and non-NULL Z prepends an intercept column", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 26)
  set.seed(260)
  Z <- matrix(rnorm(d$n), d$n, 1)   # non-NULL single column -> Z <- cbind(1, Z)
  out <- remove_covariate(d$X, d$y, Z = Z,
                          standardize = FALSE, intercept = TRUE)
  expect_equal(ncol(out$Z), 2L)
  expect_equal(out$Z[, 1], rep(1, d$n))
})

# ---- C++ (caisa.cpp) paths reached through mr.ash() / random_order ---------

test_that("random_order returns a length p*numiter 0-based permutation", {
  set.seed(27)
  o <- random_order(5L, 3L)
  expect_length(o, 15L)
  expect_true(all(o >= 0L & o <= 4L))
  expect_setequal(sort(o[1:5]), 0:4)   # first block is a permutation of 0:4
})

test_that("mr.ash method_q = 'sigma_indep_q' updates sigma2 to a positive value", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 28)
  fit <- suppressWarnings(
    mr.ash(d$X, d$y, method_q = "sigma_indep_q",
           update.sigma2 = TRUE, update.pi = FALSE,
           max.iter = 60, verbose = FALSE))
  expect_s3_class(fit, "mr.ash")
  expect_true(is.finite(fit$sigma2) && fit$sigma2 > 0)
})

test_that("mr.ash verbose = TRUE prints the termination line", {
  d <- mr_ash_data(n = 80, p = 16, k = 3, seed = 29)
  out <- suppressWarnings(capture.output(
    fit <- mr.ash(d$X, d$y, max.iter = 20, verbose = TRUE)))
  expect_true(any(grepl("Mr\\.ASH terminated", out)))
  expect_s3_class(fit, "mr.ash")
})
