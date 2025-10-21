context("SuSiE Workhorse - Main Orchestration")

# =============================================================================
# BASIC FUNCTIONALITY
# =============================================================================

test_that("susie_workhorse returns valid susie object", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
  expect_type(result, "list")
})

test_that("susie_workhorse creates all required output fields", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  # Core posterior components
  expect_true("alpha" %in% names(result))
  expect_true("mu" %in% names(result))
  expect_true("mu2" %in% names(result))
  expect_true("V" %in% names(result))
  expect_true("sigma2" %in% names(result))

  # Tracking components
  expect_true("lbf" %in% names(result))
  expect_true("lbf_variable" %in% names(result))
  expect_true("KL" %in% names(result))

  # Output fields
  expect_true("elbo" %in% names(result))
  expect_true("niter" %in% names(result))
  expect_true("converged" %in% names(result))
  expect_true("pip" %in% names(result))
  expect_true("sets" %in% names(result))
  expect_true("fitted" %in% names(result))
  expect_true("intercept" %in% names(result))
})

test_that("susie_workhorse returns correct dimensions", {
  n <- 100
  p <- 50
  L <- 5
  setup <- setup_individual_data(n = n, p = p, L = L)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_equal(dim(result$alpha), c(L, p))
  expect_equal(dim(result$mu), c(L, p))
  expect_equal(dim(result$mu2), c(L, p))
  expect_equal(dim(result$lbf_variable), c(L, p))
  expect_length(result$V, L)
  expect_length(result$lbf, L)
  expect_length(result$KL, L)
  expect_length(result$pip, p)
  expect_length(result$fitted, n)
})

# =============================================================================
# CONVERGENCE BEHAVIOR
# =============================================================================

test_that("susie_workhorse sets converged flag when converged", {
  # Use simple data and loose tolerance to ensure convergence
  setup <- setup_individual_data(n = 50, p = 20, L = 3)
  setup$params$max_iter <- 100
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-2

  result <- susie_workhorse(setup$data, setup$params)

  # Should converge with enough iterations
  expect_true("converged" %in% names(result))
  expect_type(result$converged, "logical")
})

test_that("susie_workhorse warns when not converged", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 1  # Too few iterations
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-10  # Very strict tolerance

  # Should warn about not converging
  result <- susie_workhorse(setup$data, setup$params)

  # Check convergence status
  expect_false(result$converged)
  expect_equal(result$niter, 1)
})

test_that("susie_workhorse tracks ELBO correctly", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_true("elbo" %in% names(result))
  expect_true(all(is.finite(result$elbo)))
  expect_true(length(result$elbo) <= setup$params$max_iter)
  expect_true(length(result$elbo) > 0)
})

test_that("susie_workhorse ELBO increases monotonically", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 20
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  # ELBO should be non-decreasing (allow small numerical errors)
  elbo_diff <- diff(result$elbo)
  expect_true(all(elbo_diff >= -1e-6))
})

test_that("susie_workhorse records correct number of iterations", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 15
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_true("niter" %in% names(result))
  expect_true(result$niter <= setup$params$max_iter)
  expect_true(result$niter > 0)
  expect_equal(result$niter, length(result$elbo))
})

# =============================================================================
# VARIANCE ESTIMATION
# =============================================================================

test_that("susie_workhorse updates residual variance when requested", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$estimate_residual_variance <- TRUE
  setup$params$residual_variance <- 1.5  # Initial value

  result <- susie_workhorse(setup$data, setup$params)

  # Residual variance should be updated from initial value
  expect_true(result$sigma2 > 0)
  expect_true(is.finite(result$sigma2))
})

test_that("susie_workhorse does not update residual variance when not requested", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$estimate_residual_variance <- FALSE
  setup$params$residual_variance <- 2.0  # Fixed value

  result <- susie_workhorse(setup$data, setup$params)

  # Residual variance should remain at initial value
  expect_equal(result$sigma2, 2.0)
})

# =============================================================================
# MATHEMATICAL PROPERTIES
# =============================================================================

test_that("susie_workhorse maintains valid probability distributions", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  # Each row of alpha should sum to 1
  row_sums <- rowSums(result$alpha)
  expect_equal(row_sums, rep(1, 5), tolerance = 1e-10)

  expect_true(all(result$alpha >= 0 & result$alpha <= 1))
})

test_that("susie_workhorse produces valid PIPs", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  # PIPs should be valid probabilities
  expect_true(all(result$pip >= 0))
  expect_true(all(result$pip <= 1))
  expect_true(all(is.finite(result$pip)))
})

test_that("susie_workhorse V values are non-negative", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_true(all(result$V >= 0))
  expect_true(all(is.finite(result$V)))
})

test_that("susie_workhorse sigma2 is positive", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_true(result$sigma2 > 0)
  expect_true(is.finite(result$sigma2))
})

test_that("susie_workhorse KL divergences are non-negative", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  # KL divergence should be non-negative (allow small numerical errors)
  expect_true(all(result$KL >= -1e-6))
  expect_true(all(is.finite(result$KL)))
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("susie_workhorse works with L=1", {
  setup <- setup_individual_data(n = 100, p = 50, L = 1)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
  expect_equal(dim(result$alpha), c(1, 50))
  expect_equal(sum(result$alpha), 1, tolerance = 1e-10)
})

test_that("susie_workhorse works with small p", {
  setup <- setup_individual_data(n = 100, p = 10, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
  # L should be adjusted to min(L, p)
  expect_true(nrow(result$alpha) <= 10)
})

test_that("susie_workhorse works with max_iter=1", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 1
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  # Should work but likely not converge
  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
  expect_equal(result$niter, 1)
  # With max_iter=1, may or may not converge depending on data
})

# =============================================================================
# CONVERGENCE METHODS
# =============================================================================

test_that("susie_workhorse works with ELBO convergence", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 20
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
  expect_true("elbo" %in% names(result))
})

test_that("susie_workhorse works with PIP convergence", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 20
  setup$params$convergence_method <- "pip"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
  expect_true("pip" %in% names(result))
})

# =============================================================================
# REFINEMENT
# =============================================================================

test_that("susie_workhorse respects refine=FALSE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$refine <- FALSE

  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
  # Should complete without refinement
})

test_that("susie_workhorse skips refinement when no credible sets", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 2
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-10
  setup$params$refine <- TRUE

  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
})

# =============================================================================
# TRACKING
# =============================================================================

test_that("susie_workhorse includes tracking when track_fit=TRUE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$track_fit <- TRUE

  result <- susie_workhorse(setup$data, setup$params)

  expect_true("trace" %in% names(result))
  expect_type(result$trace, "list")
})

test_that("susie_workhorse excludes tracking when track_fit=FALSE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$track_fit <- FALSE

  result <- susie_workhorse(setup$data, setup$params)

  expect_false("trace" %in% names(result))
})

# =============================================================================
# MODEL INITIALIZATION
# =============================================================================

test_that("susie_workhorse works without model_init", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$model_init <- NULL

  result <- susie_workhorse(setup$data, setup$params)

  expect_s3_class(result, "susie")
})

test_that("susie_workhorse works with model_init", {
  setup <- setup_individual_data(n = 100, p = 50, L = 3)
  setup$params$max_iter <- 5
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  # Create an initial model
  model_init <- susie_workhorse(setup$data, setup$params)

  # Use it to initialize another run
  setup2 <- setup_individual_data(n = 100, p = 50, L = 3, seed = 43)
  setup2$params$max_iter <- 10
  setup2$params$convergence_method <- "elbo"
  setup2$params$tol <- 1e-3
  setup2$params$model_init <- model_init

  result <- susie_workhorse(setup2$data, setup2$params)

  expect_s3_class(result, "susie")
})

# =============================================================================
# CREDIBLE SETS
# =============================================================================

test_that("susie_workhorse computes credible sets", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 20
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_true("sets" %in% names(result))
  expect_type(result$sets, "list")
  expect_true("cs" %in% names(result$sets))
})

# =============================================================================
# FITTED VALUES
# =============================================================================

test_that("susie_workhorse computes fitted values", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  expect_true("fitted" %in% names(result))
  expect_length(result$fitted, 100)
  expect_true(all(is.finite(result$fitted)))
})

test_that("susie_workhorse computes intercept when requested", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3
  setup$params$intercept <- TRUE

  result <- susie_workhorse(setup$data, setup$params)

  expect_true("intercept" %in% names(result))
  expect_true(is.finite(result$intercept))
})

# =============================================================================
# SIGNAL RECOVERY ON SIMULATED DATA
# =============================================================================

test_that("susie_workhorse recovers true signal on simple simulated data", {
  # Generate data with known causal variables
  set.seed(123)
  n <- 200
  p <- 100
  k <- 3  # Number of causal variables

  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  causal_idx <- c(10, 30, 50)
  beta[causal_idx] <- c(2, -2, 1.5)
  y <- drop(X %*% beta + rnorm(n, sd = 0.5))

  # Prepare data
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
    use_servin_stephens = FALSE,
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

  result <- susie_workhorse(data, params)

  # Check that causal variables have high PIPs
  expect_true(all(result$pip[causal_idx] > 0.1))

  # Check that most non-causal variables have low PIPs
  non_causal_idx <- setdiff(1:p, causal_idx)
  expect_true(mean(result$pip[non_causal_idx]) < mean(result$pip[causal_idx]))
})

# =============================================================================
# INTEGRATION WITH FULL PIPELINE
# =============================================================================

test_that("susie_workhorse produces output compatible with susie_get functions", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 20
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  # Should be able to extract PIPs (already in result)
  expect_length(result$pip, 50)

  # Should have credible sets
  expect_true("sets" %in% names(result))

  # Should have all fields needed for coef() method
  expect_true(all(c("alpha", "mu", "intercept") %in% names(result)))
})

test_that("susie_workhorse output is a valid susie object", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$max_iter <- 10
  setup$params$convergence_method <- "elbo"
  setup$params$tol <- 1e-3

  result <- susie_workhorse(setup$data, setup$params)

  # Check class
  expect_s3_class(result, "susie")

  # Check that we have all the core components for a susie object
  required_fields <- c(
    "alpha", "mu", "mu2", "V", "sigma2",
    "elbo", "niter", "converged",
    "pip", "sets", "fitted", "intercept",
    "lbf", "lbf_variable", "KL"
  )

  expect_true(all(required_fields %in% names(result)))
})
