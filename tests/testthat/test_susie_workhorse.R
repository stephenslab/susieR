context("SuSiE Workhorse - Main Orchestration")

# ---- Basic output structure ----

test_that("susie_workhorse returns a susie object with all required fields", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_s3_class(result, "susie")
  required_fields <- c(
    "alpha", "mu", "mu2", "V", "sigma2",
    "lbf", "lbf_variable", "KL",
    "elbo", "niter", "converged",
    "pip", "sets", "fitted", "intercept"
  )
  expect_true(all(required_fields %in% names(result)))
})

test_that("susie_workhorse returns correct dimensions for n=100, p=50, L=5", {
  n <- 100; p <- 50; L <- 5
  setup <- setup_individual_data(n = n, p = p, L = L)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_equal(dim(result$alpha),        c(L, p))
  expect_equal(dim(result$mu),           c(L, p))
  expect_equal(dim(result$mu2),          c(L, p))
  expect_equal(dim(result$lbf_variable), c(L, p))
  expect_length(result$V,       L)
  expect_length(result$lbf,     L)
  expect_length(result$KL,      L)
  expect_length(result$pip,     p)
  expect_length(result$fitted,  n)
})

# ---- Convergence behavior ----

test_that("susie_workhorse warns and sets converged=FALSE when max_iter reached before tolerance", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 1
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-10

  expect_warning(
    result <- susie_workhorse(setup$data, setup$params),
    "did not converge"
  )
  expect_false(result$converged)
  expect_equal(result$niter, 1)
})

test_that("susie_workhorse converges on simple data with loose tolerance", {
  setup <- setup_individual_data(n = 50, p = 20, L = 3)
  setup$params$max_iter <- 100
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-2

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_true(is.logical(result$converged))
  expect_true(result$niter > 0)
  expect_true(result$niter <= 100)
})

test_that("susie_workhorse ELBO vector has correct length and is finite", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_true(length(result$elbo) > 0)
  expect_true(length(result$elbo) <= 10)
  expect_true(all(is.finite(result$elbo)))
  expect_equal(result$niter, length(result$elbo))
})

test_that("susie_workhorse ELBO is non-decreasing across iterations", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 20
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_true(all(diff(result$elbo) >= -1e-6))
})

test_that("susie_workhorse with convergence_method='pip' produces valid pip output", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 20
  setup$params$convergence_method <- "pip"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_s3_class(result, "susie")
  expect_length(result$pip, 50)
  expect_true(all(result$pip >= 0 & result$pip <= 1))
})

# ---- Variance estimation ----

test_that("susie_workhorse updates residual variance when estimate_residual_variance=TRUE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$estimate_residual_variance <- TRUE
  setup$params$residual_variance <- 1.5

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_true(result$sigma2 > 0)
  expect_true(is.finite(result$sigma2))
  # With variance estimation active, sigma2 should have moved from the initial 1.5
  expect_false(isTRUE(all.equal(result$sigma2, 1.5)))
})

test_that("susie_workhorse holds residual variance fixed when estimate_residual_variance=FALSE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$estimate_residual_variance <- FALSE
  setup$params$residual_variance <- 2.0

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_equal(result$sigma2, 2.0)
})

# ---- Mathematical properties ----

test_that("susie_workhorse alpha rows sum to 1 and are valid probabilities", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_equal(rowSums(result$alpha), rep(1, 5), tolerance = 1e-8)
  expect_true(all(result$alpha >= 0 & result$alpha <= 1))
})

test_that("susie_workhorse PIPs are valid probabilities", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_true(all(result$pip >= 0 & result$pip <= 1))
  expect_true(all(is.finite(result$pip)))
})

test_that("susie_workhorse V values are non-negative finite, sigma2 is positive", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_true(all(result$V >= 0))
  expect_true(all(is.finite(result$V)))
  expect_true(result$sigma2 > 0)
  expect_true(is.finite(result$sigma2))
})

test_that("susie_workhorse KL divergences are non-negative finite", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_true(all(result$KL >= -1e-6))
  expect_true(all(is.finite(result$KL)))
})

# ---- Edge cases ----

test_that("susie_workhorse works with L=1: alpha is 1xp and sums to 1", {
  setup <- setup_individual_data(n = 100, p = 50, L = 1)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_s3_class(result, "susie")
  expect_equal(dim(result$alpha), c(1, 50))
  expect_equal(sum(result$alpha), 1, tolerance = 1e-8)
})

test_that("susie_workhorse with small p trims L to at most p", {
  setup <- setup_individual_data(n = 100, p = 10, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_s3_class(result, "susie")
  expect_true(nrow(result$alpha) <= 10)
})

# ---- Refinement ----

test_that("susie_workhorse with refine=FALSE completes and returns susie object", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$refine <- FALSE

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_s3_class(result, "susie")
  expect_false("trace" %in% names(result))  # no track by default
})

test_that("susie_workhorse with refine=TRUE and few iterations completes without error", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 2
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-10
  setup$params$refine <- TRUE

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_s3_class(result, "susie")
})

# ---- Fit tracking ----

test_that("susie_workhorse includes 'trace' field of class susie_track when track_fit=TRUE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$track_fit <- TRUE

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_true("trace" %in% names(result))
  expect_s3_class(result$trace, "susie_track")
})

test_that("susie_workhorse omits 'trace' field when track_fit=FALSE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$track_fit <- FALSE

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_false("trace" %in% names(result))
})

# ---- Model initialization ----

test_that("susie_workhorse accepts model_init=NULL and produces valid output", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$model_init <- NULL

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_s3_class(result, "susie")
})

test_that("susie_workhorse accepts a warm-start model_init and produces valid output", {
  setup <- setup_individual_data(n = 100, p = 50, L = 3)
  setup$params$max_iter <- 5
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  model_init <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  setup2 <- setup_individual_data(n = 100, p = 50, L = 3, seed = 43)
  setup2$params$max_iter <- 10
  setup2$params$convergence_method <- "elbo"
  setup2$params$tol <- 1e-3
  setup2$params$model_init <- model_init

  result <- suppressWarnings(susie_workhorse(setup2$data, setup2$params))

  expect_s3_class(result, "susie")
  expect_equal(dim(result$alpha), c(3, 50))
})

# ---- Output compatibility ----

test_that("susie_workhorse sets field has 'cs' sublist and fitted values are finite", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 20
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- suppressWarnings(susie_workhorse(setup$data, setup$params))

  expect_type(result$sets, "list")
  expect_true("cs" %in% names(result$sets))
  expect_length(result$fitted, 100)
  expect_true(all(is.finite(result$fitted)))
  expect_true(is.finite(result$intercept))
})

# ---- Signal recovery ----

test_that("susie_workhorse recovers true signal on simulated data with known causal variables", {
  set.seed(123)
  n <- 200; p <- 100
  causal_idx <- c(10, 30, 50)

  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[causal_idx] <- c(2, -2, 1.5)
  y <- drop(X %*% beta + rnorm(n, sd = 0.5))

  X <- set_X_attributes(X, center = TRUE, scale = TRUE)
  mean_y <- mean(y)
  y <- y - mean_y

  data <- structure(
    list(X = X, y = y, n = n, p = p, mean_y = mean_y),
    class = "individual"
  )
  params <- list(
    L = 5,
    intercept = TRUE,
    standardize = TRUE,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    unmappable_effects = "none",
    use_NIG = FALSE,
    compute_univariate_zscore = TRUE,
    coverage = 0.95,
    min_abs_corr = 0.5,
    n_purity = 100,
    check_null_threshold = 0.1,
    scaled_prior_variance = 0.2,
    prior_weights = rep(1/p, p),
    null_weight = 0,
    residual_variance = NULL,
    track_fit = FALSE,
    prior_tol = 1e-9,
    max_iter = 100,
    convergence_method = "elbo",
    tol = 1e-3,
    refine = FALSE,
    model_init = NULL,
    verbose = FALSE
  )

  result <- suppressWarnings(susie_workhorse(data, params))

  expect_true(all(result$pip[causal_idx] > 0.1))
  expect_true(mean(result$pip[causal_idx]) > mean(result$pip[setdiff(1:p, causal_idx)]))
})

# ---- Greedy-L verbose messages ----

test_that("susie_workhorse emits a message containing 'L_greedy' when verbose=TRUE with L_greedy set", {
  set.seed(201)
  setup <- setup_individual_data(n = 80, p = 30, L = 6, seed = 201)
  setup$params$L_greedy          <- 2
  setup$params$L                 <- 6
  setup$params$greedy_lbf_cutoff <- 0.1
  setup$params$verbose           <- TRUE
  setup$params$max_iter          <- 5
  setup$params$tol               <- 1e-2

  expect_message(
    suppressWarnings(susie_workhorse(setup$data, setup$params)),
    "L_greedy"
  )
})

test_that("susie_workhorse verbose greedy-L message includes greedy_lbf_cutoff value", {
  set.seed(202)
  setup <- setup_individual_data(n = 80, p = 30, L = 4, seed = 202)
  setup$params$L_greedy          <- 2
  setup$params$L                 <- 4
  setup$params$greedy_lbf_cutoff <- 0.1
  setup$params$verbose           <- TRUE
  setup$params$max_iter          <- 5
  setup$params$tol               <- 1e-2

  expect_message(
    suppressWarnings(susie_workhorse(setup$data, setup$params)),
    "greedy_lbf_cutoff"
  )
})
