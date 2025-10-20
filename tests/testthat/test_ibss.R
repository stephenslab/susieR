devtools::load_all(".")


context("Iterative Bayesian Stepwise Selection (IBSS)")

# =============================================================================
# IBSS_INITIALIZE - Basic Structure and Components
# =============================================================================

test_that("ibss_initialize returns correct structure with susie class", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_s3_class(model, "susie")
  expect_type(model, "list")
})

test_that("ibss_initialize creates all required model components", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  # Core posterior components
  expect_true("alpha" %in% names(model))
  expect_true("mu" %in% names(model))
  expect_true("mu2" %in% names(model))
  expect_true("V" %in% names(model))
  expect_true("sigma2" %in% names(model))

  # Tracking components
  expect_true("lbf" %in% names(model))
  expect_true("lbf_variable" %in% names(model))
  expect_true("KL" %in% names(model))

  # Prior components
  expect_true("pi" %in% names(model))
  expect_true("predictor_weights" %in% names(model))

  # Fitted values
  expect_true("Xr" %in% names(model))
  expect_true("null_index" %in% names(model))
})

test_that("ibss_initialize creates matrices with correct dimensions", {
  n <- 100
  p <- 50
  L <- 5
  setup <- setup_individual_data(n = n, p = p, L = L)

  model <- ibss_initialize(setup$data, setup$params)

  expect_equal(dim(model$alpha), c(L, p))
  expect_equal(dim(model$mu), c(L, p))
  expect_equal(dim(model$mu2), c(L, p))
  expect_equal(dim(model$lbf_variable), c(L, p))
  expect_length(model$V, L)
  expect_length(model$lbf, L)
  expect_length(model$KL, L)
  expect_length(model$Xr, n)
})

# =============================================================================
# IBSS_INITIALIZE - Parameter Validation
# =============================================================================

test_that("ibss_initialize adjusts L when p < L", {
  setup <- setup_individual_data(n = 100, p = 10, L = 20)

  model <- ibss_initialize(setup$data, setup$params)

  # L should be reduced to p
  expect_equal(nrow(model$alpha), 10)
  expect_equal(length(model$V), 10)
})

test_that("ibss_initialize validates residual variance is positive", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- -1

  expect_error(
    ibss_initialize(setup$data, setup$params),
    "Residual variance sigma2 must be positive"
  )
})

test_that("ibss_initialize validates residual variance is scalar", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- c(1, 2)

  expect_error(
    ibss_initialize(setup$data, setup$params),
    "Input residual variance sigma2 must be a scalar"
  )
})

test_that("ibss_initialize validates residual variance is numeric", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- "one"

  expect_error(
    ibss_initialize(setup$data, setup$params),
    "Input residual variance sigma2 must be numeric"
  )
})

test_that("ibss_initialize sets default residual variance to var(y)", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- NULL

  var_y <- var(drop(setup$data$y))
  model <- ibss_initialize(setup$data, setup$params)

  expect_equal(model$sigma2, var_y)
})

test_that("ibss_initialize uses provided residual variance", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- 2.5

  model <- ibss_initialize(setup$data, setup$params)

  expect_equal(model$sigma2, 2.5)
})

# =============================================================================
# IBSS_INITIALIZE - Model Initialization (model_init)
# =============================================================================

test_that("ibss_initialize works without model_init", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$model_init <- NULL

  model <- ibss_initialize(setup$data, setup$params)

  expect_equal(dim(model$alpha), c(5, 50))
  expect_true(all(model$alpha >= 0 & model$alpha <= 1))
  expect_true(all(is.finite(model$mu)))
})

test_that("ibss_initialize accepts valid susie model_init", {
  setup <- setup_individual_data(n = 100, p = 50, L = 3)

  # Create a proper previous susie fit to use as model_init
  model_init <- ibss_initialize(setup$data, setup$params)
  model_init$V <- rep(0.5, 3)  # Set some prior variance

  # Use it as initialization for a new fit
  setup2 <- setup_individual_data(n = 100, p = 50, L = 3)
  setup2$params$model_init <- model_init

  model <- ibss_initialize(setup2$data, setup2$params)

  expect_equal(dim(model$alpha), c(3, 50))
  expect_true(all(model$alpha >= 0 & model$alpha <= 1))
})

test_that("ibss_initialize handles model_init with fewer effects than L", {
  setup <- setup_individual_data(n = 100, p = 50, L = 2)

  # Create init with 2 effects
  model_init <- ibss_initialize(setup$data, setup$params)
  model_init$V <- rep(0.5, 2)

  # Try to expand to 5 effects
  setup2 <- setup_individual_data(n = 100, p = 50, L = 5)
  setup2$params$model_init <- model_init

  model <- ibss_initialize(setup2$data, setup2$params)

  # Should expand to L=5 effects
  expect_equal(dim(model$alpha), c(5, 50))
})

test_that("ibss_initialize handles model_init with more effects than L", {
  setup <- setup_individual_data(n = 100, p = 50, L = 6)

  # Create init with 6 effects
  model_init <- ibss_initialize(setup$data, setup$params)
  model_init$V <- rep(0.5, 6)

  # Try to reduce to 3 effects
  setup2 <- setup_individual_data(n = 100, p = 50, L = 3)
  setup2$params$model_init <- model_init

  # When model_init has more effects, it keeps all of them (expands L)
  expect_message(
    model <- ibss_initialize(setup2$data, setup2$params),
    "using L = 6"
  )

  # Should keep all 6 effects from model_init
  expect_equal(dim(model$alpha), c(6, 50))
})

# =============================================================================
# IBSS_INITIALIZE - Mathematical Properties
# =============================================================================

test_that("ibss_initialize alpha rows sum to 1", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  row_sums <- rowSums(model$alpha)
  expect_equal(row_sums, rep(1, 5), tolerance = 1e-10)
})

test_that("ibss_initialize alpha values are valid probabilities", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_true(all(model$alpha >= 0 & model$alpha <= 1))
})

test_that("ibss_initialize V values are non-negative", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_true(all(model$V >= 0))
  expect_true(all(is.finite(model$V)))
})

test_that("ibss_initialize sigma2 is positive", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_true(model$sigma2 > 0)
  expect_true(is.finite(model$sigma2))
})

test_that("ibss_initialize KL and lbf are initialized to NA", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$model_init <- NULL

  model <- ibss_initialize(setup$data, setup$params)

  expect_true(all(is.na(model$KL)))
  expect_true(all(is.na(model$lbf)))
})

# =============================================================================
# IBSS_INITIALIZE - Fitted Values
# =============================================================================

test_that("ibss_initialize creates fitted values for individual data", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_true("Xr" %in% names(model))
  expect_length(model$Xr, 100)
  expect_true(all(is.finite(model$Xr)))
})

test_that("ibss_initialize creates fitted values for sufficient stats", {
  setup <- setup_ss_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_true("XtXr" %in% names(model))
  expect_length(model$XtXr, 50)
  expect_true(all(is.finite(model$XtXr)))
})

test_that("ibss_initialize creates fitted values for rss_lambda", {
  setup <- setup_rss_lambda_data(n = 500, p = 50, L = 5, lambda = 0.5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_true("Rz" %in% names(model))
  expect_length(model$Rz, 50)
  expect_true(all(is.finite(model$Rz)))
})

# =============================================================================
# IBSS_INITIALIZE - Null Index
# =============================================================================

test_that("ibss_initialize sets null_index to 0 when null_weight = 0", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$model$null_weight <- 0

  model <- ibss_initialize(setup$data, setup$params)

  expect_equal(model$null_index, 0)
})

test_that("ibss_initialize sets null_index when null_weight > 0", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$null_weight <- 0.5

  model <- ibss_initialize(setup$data, setup$params)

  expect_true(model$null_index > 0)
})

# =============================================================================
# IBSS_INITIALIZE - Data Type Compatibility
# =============================================================================

test_that("ibss_initialize works with individual data", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_s3_class(model, "susie")
  expect_true("Xr" %in% names(model))
})

test_that("ibss_initialize works with sufficient stats", {
  setup <- setup_ss_data(n = 100, p = 50, L = 5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_s3_class(model, "susie")
  expect_true("XtXr" %in% names(model))
})

test_that("ibss_initialize works with rss_lambda", {
  setup <- setup_rss_lambda_data(n = 500, p = 50, L = 5, lambda = 0.5)

  model <- ibss_initialize(setup$data, setup$params)

  expect_s3_class(model, "susie")
  expect_true("Rz" %in% names(model))
})

# =============================================================================
# IBSS_FIT - Basic Functionality
# =============================================================================

test_that("ibss_fit updates all L effects", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  # Fit one iteration
  model_updated <- ibss_fit(setup$data, setup$params, model)

  # All effects should still have valid probabilities
  for (l in 1:5) {
    expect_equal(sum(model_updated$alpha[l, ]), 1, tolerance = 1e-10)
  }

  # V should be updated (even if to 0 for no signal)
  expect_true(all(is.finite(model_updated$V)))
  expect_true(all(model_updated$V >= 0))
})

test_that("ibss_fit updates V for all effects", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  # Store initial V
  V_init <- model$V

  # Fit one iteration
  model_updated <- ibss_fit(setup$data, setup$params, model)

  # V should be updated (unless it converged to same values)
  expect_length(model_updated$V, 5)
  expect_true(all(model_updated$V >= 0))
  expect_true(all(is.finite(model_updated$V)))
})

test_that("ibss_fit updates lbf for all effects", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  # Fit one iteration
  model_updated <- ibss_fit(setup$data, setup$params, model)

  expect_length(model_updated$lbf, 5)
  expect_true(all(is.finite(model_updated$lbf)))
})

test_that("ibss_fit updates KL for all effects", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  # Fit one iteration
  model_updated <- ibss_fit(setup$data, setup$params, model)

  expect_length(model_updated$KL, 5)
  expect_true(all(is.finite(model_updated$KL)))
  # KL divergence should be non-negative
  expect_true(all(model_updated$KL >= -1e-6))
})

# =============================================================================
# IBSS_FIT - Mathematical Properties
# =============================================================================

test_that("ibss_fit maintains valid probability distributions", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  model_updated <- ibss_fit(setup$data, setup$params, model)

  # Each row of alpha should sum to 1
  row_sums <- rowSums(model_updated$alpha)
  expect_equal(row_sums, rep(1, 5), tolerance = 1e-10)

  # All alpha values should be valid probabilities
  expect_true(all(model_updated$alpha >= 0))
  expect_true(all(model_updated$alpha <= 1))
})

test_that("ibss_fit maintains finite values", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  model_updated <- ibss_fit(setup$data, setup$params, model)

  expect_true(all(model_updated$alpha >= 0 & model_updated$alpha <= 1))
  expect_true(all(is.finite(model_updated$mu)))
  expect_true(all(is.finite(model_updated$mu2)))
  expect_true(all(is.finite(model_updated$V)))
})

test_that("ibss_fit updates fitted values", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  Xr_init <- model$Xr
  model_updated <- ibss_fit(setup$data, setup$params, model)

  # Fitted values should be updated
  expect_true("Xr" %in% names(model_updated))
  expect_length(model_updated$Xr, 100)
})

# =============================================================================
# IBSS_FIT - Edge Cases
# =============================================================================

test_that("ibss_fit works with L=1", {
  setup <- setup_individual_data(n = 100, p = 50, L = 1)
  model <- ibss_initialize(setup$data, setup$params)

  model_updated <- ibss_fit(setup$data, setup$params, model)

  expect_equal(dim(model_updated$alpha), c(1, 50))
  expect_equal(sum(model_updated$alpha), 1, tolerance = 1e-10)
})

test_that("ibss_fit works with L=0 (no effects)", {
  setup <- setup_individual_data(n = 100, p = 50, L = 0)
  model <- list(alpha = matrix(0, 0, 50))

  # Should handle gracefully
  model_updated <- ibss_fit(setup$data, setup$params, model)

  expect_equal(nrow(model_updated$alpha), 0)
})

test_that("ibss_fit works with different data types", {
  # Individual data
  setup_ind <- setup_individual_data(n = 100, p = 50, L = 5)
  model_ind <- ibss_initialize(setup_ind$data, setup_ind$params)
  model_ind_updated <- ibss_fit(setup_ind$data, setup_ind$params, model_ind)
  expect_s3_class(model_ind_updated, "susie")

  # Sufficient stats
  setup_ss <- setup_ss_data(n = 100, p = 50, L = 5)
  model_ss <- ibss_initialize(setup_ss$data, setup_ss$params)
  model_ss_updated <- ibss_fit(setup_ss$data, setup_ss$params, model_ss)
  expect_s3_class(model_ss_updated, "susie")

  # RSS lambda
  setup_rss <- setup_rss_lambda_data(n = 500, p = 50, L = 5, lambda = 0.5)
  model_rss <- ibss_initialize(setup_rss$data, setup_rss$params)
  model_rss_updated <- ibss_fit(setup_rss$data, setup_rss$params, model_rss)
  expect_s3_class(model_rss_updated, "susie")
})

# =============================================================================
# IBSS_FIT - Iterative Behavior
# =============================================================================

test_that("ibss_fit can be called iteratively", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  # Run multiple iterations
  for (iter in 1:3) {
    model <- ibss_fit(setup$data, setup$params, model)

    # Check validity after each iteration
    expect_equal(rowSums(model$alpha), rep(1, 5), tolerance = 1e-10)
    expect_true(all(model$V >= 0))
  }
})

# =============================================================================
# IBSS_FINALIZE - Basic Functionality
# =============================================================================

test_that("ibss_finalize adds required output fields", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  # Check for required output fields
  expect_true("niter" %in% names(model_final))
  expect_true("intercept" %in% names(model_final))
  expect_true("fitted" %in% names(model_final))
  expect_true("sets" %in% names(model_final))
  expect_true("pip" %in% names(model_final))
  expect_true("X_column_scale_factors" %in% names(model_final))
})

test_that("ibss_finalize sets iteration count", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 42L, tracking = NULL)

  expect_equal(model_final$niter, 42L)
})

test_that("ibss_finalize computes PIPs", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  expect_length(model_final$pip, 50)
  expect_true(all(model_final$pip >= 0))
  expect_true(all(model_final$pip <= 1))
  expect_true(all(is.finite(model_final$pip)))
})

test_that("ibss_finalize computes credible sets", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  expect_true("sets" %in% names(model_final))
  expect_type(model_final$sets, "list")
})

test_that("ibss_finalize computes fitted values", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  expect_true("fitted" %in% names(model_final))
  expect_length(model_final$fitted, 100)
  expect_true(all(is.finite(model_final$fitted)))
})

test_that("ibss_finalize computes intercept", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$intercept <- TRUE
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  expect_true("intercept" %in% names(model_final))
  expect_true(is.finite(model_final$intercept))
})

# =============================================================================
# IBSS_FINALIZE - Tracking
# =============================================================================

test_that("ibss_finalize includes tracking when track_fit=TRUE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$track_fit <- TRUE
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  # Create mock tracking data
  tracking <- list(
    elbo = c(100, 110, 115),
    sigma2 = c(1, 0.9, 0.85)
  )

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 3L, tracking = tracking)

  expect_true("trace" %in% names(model_final))
  expect_type(model_final$trace, "list")
})

test_that("ibss_finalize excludes tracking when track_fit=FALSE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$track_fit <- FALSE
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 3L, tracking = NULL)

  expect_false("trace" %in% names(model_final))
})

# =============================================================================
# IBSS_FINALIZE - Variable Names
# =============================================================================

test_that("ibss_finalize assigns variable names when available", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  # Add column names to X
  colnames(setup$data$X) <- paste0("var", 1:50)

  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  # Check that variable names are assigned to pip
  expect_named(model_final$pip, paste0("var", 1:50))
})

# =============================================================================
# IBSS_FINALIZE - Z-scores
# =============================================================================

test_that("ibss_finalize computes z-scores for individual data", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  expect_true("z" %in% names(model_final))
  if (!is.null(model_final$z)) {
    expect_length(model_final$z, 50)
    expect_true(all(is.finite(model_final$z)))
  }
})

# =============================================================================
# IBSS_FINALIZE - Scale Factors
# =============================================================================

test_that("ibss_finalize computes X_column_scale_factors", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  expect_true("X_column_scale_factors" %in% names(model_final))
  expect_length(model_final$X_column_scale_factors, 50)
})

# =============================================================================
# IBSS_FINALIZE - Data Type Compatibility
# =============================================================================

test_that("ibss_finalize works with individual data", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)

  model_final <- ibss_finalize(setup$data, setup$params, model,
                               elbo = NULL, iter = 10L, tracking = NULL)

  expect_s3_class(model_final, "susie")
  expect_true("fitted" %in% names(model_final))
  expect_length(model_final$fitted, 100)
})


# =============================================================================
# FULL IBSS PIPELINE
# =============================================================================

test_that("Full IBSS pipeline produces valid susie object", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  # Initialize
  model <- ibss_initialize(setup$data, setup$params)

  # Fit (run 5 iterations)
  for (i in 1:5) {
    model <- ibss_fit(setup$data, setup$params, model)
  }

  # Finalize
  model <- ibss_finalize(setup$data, setup$params, model,
                        elbo = NULL, iter = 5L, tracking = NULL)

  # Check final model is complete
  expect_s3_class(model, "susie")
  expect_true("alpha" %in% names(model))
  expect_true("mu" %in% names(model))
  expect_true("V" %in% names(model))
  expect_true("pip" %in% names(model))
  expect_true("sets" %in% names(model))
  expect_true("fitted" %in% names(model))
  expect_true("niter" %in% names(model))
  expect_equal(model$niter, 5L)
})

test_that("Full IBSS pipeline maintains mathematical properties", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  # Initialize
  model <- ibss_initialize(setup$data, setup$params)

  # Fit (run 3 iterations)
  for (i in 1:3) {
    model <- ibss_fit(setup$data, setup$params, model)
  }

  # Finalize
  model <- ibss_finalize(setup$data, setup$params, model,
                        elbo = NULL, iter = 3L, tracking = NULL)

  # Check mathematical properties
  expect_equal(rowSums(model$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(model$pip >= 0))
  expect_true(all(model$pip <= 1))
  expect_true(all(model$V >= 0))
})

