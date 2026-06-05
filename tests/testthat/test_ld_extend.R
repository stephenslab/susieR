context("LD extension of credible sets (ld_extend_threshold)")

# =============================================================================
# extend_cs_by_ld() helper (internal)
# =============================================================================

test_that("extend_cs_by_ld adds tight-LD proxies via the Xcorr branch", {
  R <- diag(5)
  R[1, 4] <- R[4, 1] <- 0.995          # variant 4 is a near-perfect proxy of 1
  cs <- list(c(1L))

  ext <- extend_cs_by_ld(cs, X = NULL, Xcorr = R, threshold = 0.99)
  expect_true(setequal(ext[[1]], c(1, 4)))

  # Higher threshold than the proxy correlation -> no extension.
  ext_none <- extend_cs_by_ld(cs, X = NULL, Xcorr = R, threshold = 0.999)
  expect_true(setequal(ext_none[[1]], 1))
})

test_that("extend_cs_by_ld X branch equals the Xcorr branch (the oracle)", {
  set.seed(101)
  X <- matrix(rnorm(80 * 6), 80, 6)
  X[, 4] <- X[, 1] + rnorm(80, sd = 1e-6)   # |cor(1,4)| > 0.99
  R <- cor(X)

  e_x <- extend_cs_by_ld(list(c(1L)), X = X,    Xcorr = NULL, threshold = 0.99)
  e_r <- extend_cs_by_ld(list(c(1L)), X = NULL, Xcorr = R,    threshold = 0.99)
  expect_identical(lapply(e_x, sort), lapply(e_r, sort))
  expect_true(setequal(e_x[[1]], c(1, 4)))
})

test_that("extend_cs_by_ld is a no-op for NULL threshold or empty cs", {
  R <- diag(3); R[1, 2] <- R[2, 1] <- 0.999
  cs <- list(c(1L))
  expect_identical(extend_cs_by_ld(cs, NULL, R, threshold = NULL), cs)
  expect_identical(extend_cs_by_ld(list(), NULL, R, threshold = 0.99), list())
  # No correlation info at all -> unchanged.
  expect_identical(extend_cs_by_ld(cs, NULL, NULL, threshold = 0.99), cs)
})

test_that("extend_cs_by_ld never pulls in the null index", {
  R <- diag(5)
  R[1, 4] <- R[4, 1] <- 0.995
  ext <- extend_cs_by_ld(list(c(1L)), X = NULL, Xcorr = R,
                         threshold = 0.99, null_index = 4)
  expect_true(setequal(ext[[1]], 1))     # variant 4 would be added but is the null
})

test_that("extend_cs_by_ld skips (with a hint) a sparse X instead of densifying", {
  set.seed(102)
  X <- matrix(rnorm(60 * 5), 60, 5)
  X[, 4] <- X[, 1] + rnorm(60, sd = 1e-6)
  spX <- Matrix::Matrix(X, sparse = TRUE)

  expect_message(
    res <- extend_cs_by_ld(list(c(1L)), X = spX, Xcorr = NULL, threshold = 0.99),
    "sparse"
  )
  expect_identical(res, list(c(1L)))     # unchanged because extension was skipped
})

# =============================================================================
# susie_get_cs(): extension wiring, default OFF, and the inline guard
# =============================================================================

# Minimal hand-built fit with deterministic credible sets:
#   effect 1 -> {1}, effect 2 -> {3} at coverage 0.95.
make_fit <- function() {
  res <- list(
    alpha = rbind(c(0.97, 0.01, 0.005, 0.01, 0.005),
                  c(0.005, 0.01, 0.97, 0.01, 0.005)),
    V = c(1, 1),
    null_index = 0
  )
  res
}

test_that("susie_get_cs extends with a threshold and not by default", {
  res <- make_fit()
  R <- diag(5); R[1, 4] <- R[4, 1] <- 0.995

  off <- susie_get_cs(res, Xcorr = R)                              # default NULL
  on  <- susie_get_cs(res, Xcorr = R, ld_extend_threshold = 0.99)

  off_sets <- lapply(off$cs, sort)
  on_sets  <- lapply(on$cs,  sort)
  expect_false(any(vapply(off_sets, setequal, logical(1), c(1, 4))))
  expect_true(any(vapply(off_sets, setequal, logical(1), 1)))      # un-extended {1}
  expect_true(any(vapply(on_sets,  setequal, logical(1), c(1, 4)))) # extended {1,4}
})

test_that("susie_get_cs X and Xcorr pathways agree end-to-end with extension on", {
  set.seed(103)
  dat <- simulate_regression(n = 200, p = 30, k = 2, signal_sd = 2)
  X <- dat$X
  X[, 2] <- X[, 1] + rnorm(nrow(X), sd = 0.05)      # tight-LD proxy
  fit <- suppressWarnings(susie(X, dat$y, L = 5, verbose = FALSE))

  a <- susie_get_cs(fit, X = X,         ld_extend_threshold = 0.99, n_purity = ncol(X))
  b <- susie_get_cs(fit, Xcorr = cor(X), ld_extend_threshold = 0.99, n_purity = ncol(X))
  expect_identical(lapply(a$cs, sort), lapply(b$cs, sort))
})

test_that("susie_get_cs validates ld_extend_threshold (inline guard)", {
  res <- make_fit()
  R <- diag(5)
  for (bad in list(1.5, -0.1, NA_real_, Inf, c(0.9, 0.99), "x")) {
    expect_error(
      susie_get_cs(res, Xcorr = R, ld_extend_threshold = bad),
      "ld_extend_threshold must be NULL or a single numeric value"
    )
  }
  # Boundary and NULL values are accepted.
  for (ok in list(NULL, 0, 1, 0.99))
    expect_silent(susie_get_cs(res, Xcorr = R, ld_extend_threshold = ok))
})

# =============================================================================
# validate_and_override_params(): the new validation block (via a constructor,
# which assembles a complete params list before validating)
# =============================================================================

test_that("constructor validation rejects bad ld_extend_threshold and accepts good", {
  dat <- simulate_regression(n = 100, p = 10, k = 2)
  expect_error(
    individual_data_constructor(dat$X, dat$y, ld_extend_threshold = 1.5),
    "ld_extend_threshold must be NULL or a single numeric value in \\[0, 1\\]"
  )
  expect_error(
    individual_data_constructor(dat$X, dat$y, ld_extend_threshold = c(0.5, 0.6)),
    "ld_extend_threshold must be NULL"
  )
  expect_error(
    individual_data_constructor(dat$X, dat$y, ld_extend_threshold = "x"),
    "ld_extend_threshold must be NULL"
  )
  expect_equal(
    individual_data_constructor(dat$X, dat$y, ld_extend_threshold = 0)$params$ld_extend_threshold, 0
  )
  expect_equal(
    individual_data_constructor(dat$X, dat$y, ld_extend_threshold = 1)$params$ld_extend_threshold, 1
  )
})

# =============================================================================
# Parameter threading through every constructor
# =============================================================================

test_that("all constructors thread ld_extend_threshold into params", {
  set.seed(104)
  dat <- simulate_regression(n = 200, p = 20, k = 2, signal_sd = 2)
  X <- dat$X; y <- dat$y
  XtX <- crossprod(X); Xty <- as.vector(crossprod(X, y)); yty <- sum(y^2)
  R <- cor(X); z <- Xty / sqrt(diag(XtX)); n <- nrow(X)

  ind <- individual_data_constructor(X, y, ld_extend_threshold = 0.99)
  expect_equal(ind$params$ld_extend_threshold, 0.99)

  ss <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = XtX,
                                     ld_extend_threshold = 0.7)
  expect_equal(ss$params$ld_extend_threshold, 0.7)

  rl <- suppressWarnings(rss_lambda_constructor(z = z, R = R, n = n, lambda = 0.5,
                                                ld_extend_threshold = 0.8))
  expect_equal(rl$params$ld_extend_threshold, 0.8)

  mix <- suppressWarnings(ss_mixture_constructor(z = z, R = list(R, R), n = n,
                                                 ld_extend_threshold = 0.6))
  expect_equal(mix$params$ld_extend_threshold, 0.6)

  # summary_stats_constructor forwards to sufficient_stats (single panel)...
  sm1 <- suppressWarnings(summary_stats_constructor(z = z, R = R, n = n,
                                                    ld_extend_threshold = 0.95))
  expect_equal(sm1$params$ld_extend_threshold, 0.95)
  # ...and to ss_mixture (multi panel).
  sm2 <- suppressWarnings(summary_stats_constructor(z = z, R = list(R, R), n = n,
                                                    ld_extend_threshold = 0.55))
  expect_equal(sm2$params$ld_extend_threshold, 0.55)
})

# =============================================================================
# Public fitters: accept the argument and reach the get_cs.* call sites
# =============================================================================

test_that("susie / susie_ss / susie_rss / susie_rss_lambda accept ld_extend_threshold", {
  set.seed(105)
  dat <- simulate_regression(n = 200, p = 20, k = 2, signal_sd = 2)
  X <- dat$X; y <- dat$y
  XtX <- crossprod(X); Xty <- as.vector(crossprod(X, y)); yty <- sum(y^2)
  R <- cor(X); z <- Xty / sqrt(diag(XtX)); n <- nrow(X)

  f1 <- suppressWarnings(susie(X, y, L = 5, ld_extend_threshold = 0.99))
  expect_s3_class(f1, "susie")

  f2 <- suppressWarnings(susie_ss(XtX = XtX, Xty = Xty, yty = yty, n = n,
                                  ld_extend_threshold = 0.99))
  expect_s3_class(f2, "susie")

  f3 <- suppressWarnings(susie_rss(z = z, R = R, n = n, ld_extend_threshold = 0.99))
  expect_s3_class(f3, "susie")

  f4 <- suppressWarnings(susie_rss_lambda(z = z, R = R, n = n, lambda = 0.5,
                                          ld_extend_threshold = 0.99))
  expect_s3_class(f4, "susie")

  # Invalid value is rejected at the constructor/validation layer.
  expect_error(susie(X, y, L = 5, ld_extend_threshold = 1.5),
               "ld_extend_threshold must be NULL or a single numeric value")
})

test_that("susie_rss low-rank X path threads ld_extend_threshold (get_cs.ss X branch)", {
  set.seed(106)
  dat <- simulate_regression(n = 150, p = 12, k = 2, signal_sd = 2)
  X <- dat$X; y <- dat$y
  z <- as.vector(crossprod(X, y)) / sqrt(diag(crossprod(X)))
  fit <- suppressWarnings(susie_rss(z = z, X = X, n = nrow(X),
                                    ld_extend_threshold = 0.99))
  expect_s3_class(fit, "susie")
})
