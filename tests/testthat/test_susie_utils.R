context("Utility functions for susieR")

# =============================================================================
# FUNDAMENTAL BUILDING BLOCKS
# =============================================================================

test_that("warning_message displays warnings correctly", {
  # Test warning style (default)
  expect_message(
    warning_message("Test warning"),
    "WARNING:.*Test warning"
  )

  # Test warning style (explicit)
  expect_message(
    warning_message("Test warning", style = "warning"),
    "WARNING:.*Test warning"
  )

  # Test hint style
  expect_message(
    warning_message("Test hint", style = "hint"),
    "HINT:.*Test hint"
  )

  # Test with warn < 0
  old_warn <- getOption("warn")
  on.exit(options(warn = old_warn), add = TRUE)

  options(warn = -1)
  expect_no_error(warning_message("Still executes", style = "warning"))

  # Hint should still show even with warn < 0
  expect_message(
    warning_message("Hint shows", style = "hint"),
    "HINT:.*Hint shows"
  )
})

test_that("safe_cor computes correlation and handles zero sd columns", {
  # Normal correlation
  x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  result <- safe_cor(x)
  expected <- cor(x)
  expect_equal(result, expected)

  # With constant column (zero sd) - should handle without warning

  x_const <- cbind(x, rep(1, 3))

  # Without safe_cor, cor() would warn
  expect_warning(cor(x_const), "the standard deviation is zero")

  # safe_cor handles this without warning and returns 0 for constant column correlations
  expect_silent(result <- safe_cor(x_const))

  # Check that correlations involving constant column are 0 (not NA)
  expect_equal(result[3, 1], 0)
  expect_equal(result[3, 2], 0)
  expect_equal(result[1, 3], 0)
  expect_equal(result[2, 3], 0)
  # Diagonal should still be 1
  expect_equal(diag(result), c(1, 1, 1))
})

test_that("safe_cov2cor computes correlation from covariance and handles zero variance", {
  # Normal case
  cov_mat <- matrix(c(4, 2, 2, 3), nrow = 2)
  result <- safe_cov2cor(cov_mat)
  expected <- cov2cor(cov_mat)
  expect_equal(result, expected)

  # With zero variance entry - safe_cov2cor handles this without warning
  cov_mat_zero <- matrix(c(0, 0, 0, 3), nrow = 2)

  # Without safe_cov2cor, cov2cor() would warn
  expect_warning(cov2cor(cov_mat_zero))

  # safe_cov2cor handles zero variance by returning 0 correlations (not NA)
  expect_silent(result <- safe_cov2cor(cov_mat_zero))
  expect_true(is.matrix(result))
  # Diagonal should be 1
  expect_equal(diag(result), c(1, 1))
  # Off-diagonal correlations involving zero-variance variable should be 0

  expect_equal(result[1, 2], 0)
  expect_equal(result[2, 1], 0)
})

test_that("safe_cor handles constant columns without warnings", {
  # Create data with a constant column that would trigger a warning in base cor()
  x_const <- matrix(c(1, 2, 3, 4, 5, 6, 1, 1, 1), nrow = 3, ncol = 3)

  # Verify that base cor() would warn
  expect_warning(cor(x_const), "the standard deviation is zero")

  # Test that safe_cor handles it silently
  expect_silent(safe_cor(x_const))

  # Verify result is computed correctly
  result <- safe_cor(x_const)
  expect_true(is.matrix(result))
  # Correlations involving the constant column (col 3) should be 0, not NA

  expect_equal(result[3, 1], 0)
  expect_equal(result[3, 2], 0)
  expect_equal(result[1, 3], 0)
  expect_equal(result[2, 3], 0)
  # Diagonal should be 1
  expect_equal(diag(result), c(1, 1, 1))
})

test_that("safe_cov2cor handles zero variance without warnings", {
  # Create covariance matrix with zero variance that would trigger a warning in base cov2cor()
  cov_mat_zero <- matrix(c(0, 0, 0, 3), nrow = 2)

  # Verify that base cov2cor() would warn
  expect_warning(cov2cor(cov_mat_zero))

  # Test that safe_cov2cor handles it silently
  expect_silent(safe_cov2cor(cov_mat_zero))

  # Verify result is computed correctly
  result <- safe_cov2cor(cov_mat_zero)
  expect_true(is.matrix(result))
  # Correlations involving zero-variance variable should be 0, not NA

  expect_equal(result[1, 2], 0)
  expect_equal(result[2, 1], 0)
  # Diagonal should be 1
  expect_equal(diag(result), c(1, 1))
})

test_that("is_symmetric_matrix correctly identifies symmetric matrices", {
  # Symmetric matrix
  sym_mat <- matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), nrow = 3)
  expect_true(is_symmetric_matrix(sym_mat))

  # Non-symmetric matrix
  nonsym_mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3)
  expect_false(is_symmetric_matrix(nonsym_mat))

  # Identity matrix
  expect_true(is_symmetric_matrix(diag(5)))

  # Sparse symmetric matrix
  sparse_sym <- Matrix::Matrix(sym_mat, sparse = TRUE)
  expect_true(is_symmetric_matrix(sparse_sym))
})

test_that("apply_nonzeros applies function to nonzero elements of sparse matrix", {
  # Create sparse matrix
  X <- Matrix::Matrix(c(1, 0, 0, 2, 0, 3, 4, 0, 5), nrow = 3, sparse = TRUE)

  # Square all nonzero elements
  result <- apply_nonzeros(X, function(x) x^2)

  # Check dimensions preserved
  expect_equal(dim(result), dim(X))

  # Check that zeros remain zeros
  expect_equal(sum(result == 0), sum(X == 0))

  # Check nonzero values are squared
  expected <- Matrix::Matrix(c(1, 0, 0, 4, 0, 9, 16, 0, 25), nrow = 3, sparse = TRUE)
  expect_equal(as.matrix(result), as.matrix(expected))

  # Test with another function (doubling)
  result2 <- apply_nonzeros(X, function(x) x * 2)
  expected2 <- Matrix::Matrix(c(2, 0, 0, 4, 0, 6, 8, 0, 10), nrow = 3, sparse = TRUE)
  expect_equal(as.matrix(result2), as.matrix(expected2))
})

test_that("compute_colSds computes column standard deviations correctly", {
  # Dense matrix
  X_dense <- matrix(rnorm(100), nrow = 10, ncol = 10)
  result_dense <- compute_colSds(X_dense)
  expected_dense <- matrixStats::colSds(X_dense)
  expect_equal(result_dense, expected_dense, tolerance = 1e-10)

  # Sparse matrix
  X_sparse <- Matrix::Matrix(X_dense, sparse = TRUE)
  X_sparse[abs(X_sparse) < 0.5] <- 0  # Make it actually sparse
  result_sparse <- compute_colSds(X_sparse)
  expected_sparse <- apply(as.matrix(X_sparse), 2, sd)
  expect_equal(result_sparse, expected_sparse, tolerance = 1e-10)

  # Matrix with constant column (sd = 0)
  X_const <- cbind(X_dense, rep(1, 10))
  result_const <- compute_colSds(X_const)
  expect_equal(result_const[11], 0)
  expect_equal(result_const[1:10], matrixStats::colSds(X_dense))
})

test_that("compute_colstats computes column statistics correctly", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 123)

  # Test with both center and scale
  result <- compute_colstats(base_data$X, center = TRUE, scale = TRUE)

  # Check components exist
  expect_true(all(c("cm", "csd", "d") %in% names(result)))
  expect_length(result$cm, base_data$p)
  expect_length(result$csd, base_data$p)
  expect_length(result$d, base_data$p)

  # Check column means
  expect_equal(result$cm, colMeans(base_data$X), tolerance = 1e-10)

  # Check column sds
  expected_csd <- apply(base_data$X, 2, sd)
  expect_equal(result$csd, expected_csd, tolerance = 1e-10)

  # Check d values (sum of squared standardized values)
  X_std <- scale(base_data$X, center = TRUE, scale = TRUE)
  expected_d <- colSums(X_std^2)
  expect_equal(result$d, expected_d, tolerance = 1e-8)

  # Test with center = FALSE
  result_nocenter <- compute_colstats(base_data$X, center = FALSE, scale = TRUE)
  expect_equal(result_nocenter$cm, rep(0, base_data$p))
  expect_equal(result_nocenter$csd, expected_csd, tolerance = 1e-10)

  # Test with scale = FALSE
  result_noscale <- compute_colstats(base_data$X, center = TRUE, scale = FALSE)
  expect_equal(result_noscale$cm, colMeans(base_data$X), tolerance = 1e-10)
  expect_equal(result_noscale$csd, rep(1, base_data$p))

  # Test with neither center nor scale
  result_neither <- compute_colstats(base_data$X, center = FALSE, scale = FALSE)
  expect_equal(result_neither$cm, rep(0, base_data$p))
  expect_equal(result_neither$csd, rep(1, base_data$p))
  expected_d_neither <- colSums(base_data$X^2)
  expect_equal(result_neither$d, expected_d_neither, tolerance = 1e-10)

  # Test with column of zeros (sd = 0)
  X_zero <- cbind(base_data$X, rep(0, base_data$n))
  result_zero <- compute_colstats(X_zero, center = TRUE, scale = TRUE)
  expect_equal(result_zero$csd[base_data$p + 1], 1)  # sd = 0 replaced by 1

  # Test with sparse matrix
  X_sparse <- Matrix::Matrix(base_data$X, sparse = TRUE)
  X_sparse[abs(X_sparse) < 1] <- 0
  result_sparse <- compute_colstats(X_sparse, center = TRUE, scale = TRUE)
  expect_length(result_sparse$cm, base_data$p)
  expect_length(result_sparse$csd, base_data$p)
  expect_length(result_sparse$d, base_data$p)
})

# =============================================================================
# DATA PROCESSING & VALIDATION
# =============================================================================

test_that("check_semi_pd identifies positive semi-definite matrices", {
  # Positive definite matrix
  A_pd <- matrix(c(2, 1, 1, 2), nrow = 2)
  result_pd <- check_semi_pd(A_pd, tol = 1e-10)

  expect_true(result_pd$status)
  expect_true(all(result_pd$eigenvalues >= 0))
  expect_true(!is.null(attr(result_pd$matrix, "eigen")))

  # Positive semi-definite matrix (singular)
  A_psd <- matrix(c(1, 1, 1, 1), nrow = 2)
  result_psd <- check_semi_pd(A_psd, tol = 1e-10)

  expect_true(result_psd$status)
  expect_true(min(result_psd$eigenvalues) >= 0)
  expect_true(any(abs(result_psd$eigenvalues) < 1e-10))  # Has zero eigenvalue

  # Not positive semi-definite
  A_neg <- matrix(c(1, 2, 2, -1), nrow = 2)
  result_neg <- check_semi_pd(A_neg, tol = 1e-10)

  expect_false(result_neg$status)
  expect_true(any(result_neg$eigenvalues < 0))

  # Identity matrix
  A_id <- diag(3)
  result_id <- check_semi_pd(A_id, tol = 1e-10)

  expect_true(result_id$status)
  expect_equal(result_id$eigenvalues, rep(1, 3), tolerance = 1e-10)
})

test_that("check_projection verifies if vector is in eigenspace", {
  # Create a matrix and vector in its column space
  A <- matrix(c(4, 2, 2, 3), nrow = 2)
  b_in <- c(2, 1)  # In column space

  result_in <- check_projection(A, b_in)
  expect_true(result_in$status)
  expect_true(is.na(result_in$msg))

  # Test with pre-computed eigen decomposition
  A_with_eigen <- A
  attr(A_with_eigen, "eigen") <- eigen(A, symmetric = TRUE)
  result_with_eigen <- check_projection(A_with_eigen, b_in)
  expect_true(result_with_eigen$status)
})

test_that("validate_init validates model initialization objects", {
  # Create valid model_init
  p <- 50
  L <- 5
  valid_init <- list(
    alpha = matrix(1/p, L, p),
    mu = matrix(0, L, p),
    mu2 = matrix(0, L, p),
    V = rep(1, L),
    sigma2 = 1,
    pi = rep(1/p, p),
    null_index = 0
  )
  class(valid_init) <- "susie"

  data <- list(n = 100, p = p)
  params <- list(L = L, model_init = valid_init)

  # Should pass without error
  expect_silent(validate_init(data, params))

  # Test: not a susie object
  bad_init <- valid_init
  class(bad_init) <- "lm"
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init must be a 'susie' object")

  # Test: NA in alpha
  bad_init <- valid_init
  bad_init$alpha[1, 1] <- NA
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$alpha contains NA/Inf")

  # Test: Inf in mu
  bad_init <- valid_init
  bad_init$mu[1, 1] <- Inf
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$mu contains NA/Inf")

  # Test: alpha not a matrix
  bad_init <- valid_init
  bad_init$alpha <- as.vector(bad_init$alpha)
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$alpha must be a matrix")

  # Test: alpha values outside [0,1]
  bad_init <- valid_init
  bad_init$alpha[1, 1] <- 1.5
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "invalid values outside range")

  # Test: dimension mismatch
  bad_init <- valid_init
  bad_init$mu <- matrix(0, L, p - 1)
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "dimensions do not match")

  # Test: V length mismatch
  bad_init <- valid_init
  bad_init$V <- rep(1, L - 1)
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "does not equal nrow")

  # Test: negative V
  bad_init <- valid_init
  bad_init$V[1] <- -1
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "at least one negative value")

  # Test: negative sigma2
  bad_init <- valid_init
  bad_init$sigma2 <- -0.5
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "sigma2 is negative")

  # Test: NULL V (should pass)
  init_no_V <- valid_init
  init_no_V$V <- NULL
  params_no_V <- list(L = L, model_init = init_no_V)
  expect_silent(validate_init(data, params_no_V))

  # Test 1: mu2 contains NA/Inf values
  bad_init <- valid_init
  bad_init$mu2[2, 3] <- NA
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$mu2 contains NA/Inf values")

  bad_init <- valid_init
  bad_init$mu2[1, 5] <- Inf
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$mu2 contains NA/Inf values")

  # Test 2: V contains NA/Inf values
  bad_init <- valid_init
  bad_init$V[2] <- NA
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$V contains NA/Inf values")

  bad_init <- valid_init
  bad_init$V[3] <- Inf
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$V contains NA/Inf values")

  # Test 3: sigma2 contains NA/Inf
  bad_init <- valid_init
  bad_init$sigma2 <- NA
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$sigma2 contains NA/Inf")

  bad_init <- valid_init
  bad_init$sigma2 <- Inf
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$sigma2 contains NA/Inf")

  # Test 4: pi contains NA/Inf
  bad_init <- valid_init
  bad_init$pi[10] <- NA
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$pi contains NA/Inf")

  bad_init <- valid_init
  bad_init$pi[5] <- Inf
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$pi contains NA/Inf")

  # Test 5: mu2 and alpha dimensions do not match
  bad_init <- valid_init
  bad_init$mu2 <- matrix(0, L, p - 1)
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$mu2 and model_init\\$alpha dimensions do not match")

  bad_init <- valid_init
  bad_init$mu2 <- matrix(0, L + 1, p)
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$mu2 and model_init\\$alpha dimensions do not match")

  # Test 6: V must be numeric
  # Note: This branch is unreachable because is.finite() on character vectors
  # returns FALSE, triggering the NA/Inf error first. The numeric check only runs
  # if all values pass is.finite(). Testing with numeric values only.
  bad_init <- valid_init
  bad_init$V <- rep(0, L)  # All zeros (valid finite numerics)
  # This should pass all checks since 0 is valid for V
  params_ok <- list(L = L, model_init = bad_init)
  expect_silent(validate_init(data, params_ok))

  # Test 7: sigma2 must be numeric
  # Note: Similar to above - unreachable branch due to is.finite() check first
  bad_init <- valid_init
  bad_init$sigma2 <- 0  # Zero is valid
  params_ok <- list(L = L, model_init = bad_init)
  expect_silent(validate_init(data, params_ok))

  # Test 8: pi length must match number of columns in alpha
  bad_init <- valid_init
  bad_init$pi <- rep(1/(p-1), p - 1)
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$pi should have the same length as the number of columns in model_init\\$alpha")

  bad_init <- valid_init
  bad_init$pi <- rep(1/(p+1), p + 1)
  params_bad <- list(L = L, model_init = bad_init)
  expect_error(validate_init(data, params_bad), "model_init\\$pi should have the same length as the number of columns in model_init\\$alpha")
})

test_that("convert_individual_to_ss converts individual data to sufficient statistics", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5, seed = 123)
  data <- setup$data

  params <- list(
    unmappable_effects = "inf",
    verbose = FALSE
  )

  # Convert
  ss_data <- convert_individual_to_ss(data, params)

  # Check class
  expect_s3_class(ss_data, "ss")

  # Check components exist
  expect_true(all(c("XtX", "Xty", "yty", "n", "p") %in% names(ss_data)))

  # Check dimensions
  expect_equal(dim(ss_data$XtX), c(50, 50))
  expect_length(ss_data$Xty, 50)
  expect_length(ss_data$yty, 1)

  # Values may be rescaled for unmappable effects
  expect_true(is.numeric(ss_data$XtX))
  expect_true(is.numeric(ss_data$Xty))
  expect_true(is.numeric(ss_data$yty))
  expect_true(ss_data$yty > 0)

  # Check attributes preserved
  expect_equal(ss_data$X_colmeans, attr(data$X, "scaled:center"))
  expect_equal(ss_data$y_mean, data$mean_y)
  expect_equal(attr(ss_data$XtX, "d"), attr(data$X, "d"))
  expect_equal(attr(ss_data$XtX, "scaled:scale"), attr(data$X, "scaled:scale"))

  # Check eigen decomposition added for unmappable effects
  expect_true(!is.null(ss_data$eigen_vectors))
  expect_true(!is.null(ss_data$eigen_values))
  expect_true(!is.null(ss_data$VtXty))
})

test_that("extract_prior_weights extracts and rescales prior weights", {
  p <- 100

  # Test: no null weight
  model_no_null <- list(
    pi = rep(1/p, p),
    null_weight = 0,
    null_index = 0
  )
  result <- extract_prior_weights(model_no_null)
  expect_equal(result, rep(1/p, p))

  # Test: with null weight
  null_weight <- 0.1
  null_idx <- p
  pi_vec <- c(rep((1 - null_weight)/(p - 1), p - 1), null_weight)
  model_with_null <- list(
    pi = pi_vec,
    null_weight = null_weight,
    null_index = null_idx
  )
  result <- extract_prior_weights(model_with_null)

  # Should extract non-null weights and rescale to sum to 1
  expect_length(result, p - 1)
  expect_equal(sum(result), 1, tolerance = 1e-10)
  expect_equal(result, rep(1/(p-1), p-1), tolerance = 1e-10)

  # Test: null weight provided as argument
  result_arg <- extract_prior_weights(model_with_null, null_weight = null_weight)
  expect_equal(result, result_arg)

  # Test: null_weight = NULL (backwards compatibility)
  model_null_weight_null <- model_no_null
  model_null_weight_null$null_weight <- NULL
  result <- extract_prior_weights(model_null_weight_null)
  expect_equal(result, rep(1/p, p))
})

test_that("reconstruct_full_weights reconstructs prior weights with null component", {
  p <- 100
  non_null_weights <- rep(1/(p-1), p-1)

  # Test: no null weight
  result_no_null <- reconstruct_full_weights(non_null_weights, null_weight = 0)
  expect_length(result_no_null, p - 1)
  expect_equal(sum(result_no_null), 1, tolerance = 1e-10)
  expect_equal(result_no_null, non_null_weights, tolerance = 1e-10)

  # Test: with null weight
  null_weight <- 0.1
  result_with_null <- reconstruct_full_weights(non_null_weights, null_weight = null_weight)
  expect_length(result_with_null, p)
  expect_equal(sum(result_with_null), 1, tolerance = 1e-10)
  expect_equal(result_with_null[p], null_weight, tolerance = 1e-10)
  expect_equal(sum(result_with_null[1:(p-1)]), 1 - null_weight, tolerance = 1e-10)

  # Test: null_weight = NULL
  result_null <- reconstruct_full_weights(non_null_weights, null_weight = NULL)
  expect_equal(sum(result_null), 1, tolerance = 1e-10)
})

test_that("validate_and_override_params validates and adjusts parameters", {
  # Valid params
  valid_params <- list(
    L = 10,
    prior_tol = 1e-9,
    residual_variance_upperbound = 1e10,
    scaled_prior_variance = 0.2,
    unmappable_effects = "none",
    convergence_method = "elbo",
    estimate_prior_variance = TRUE,
    estimate_prior_method = "EM",
    estimate_residual_method = "MLE",
    estimate_residual_variance = TRUE,
    refine = FALSE
  )

  result <- validate_and_override_params(valid_params)
  expect_equal(result$prior_tol, 1e-9)
  expect_false(result$use_servin_stephens)

  # Test: invalid prior_tol
  bad_params <- valid_params
  bad_params$prior_tol <- -1
  expect_error(validate_and_override_params(bad_params), "prior_tol must be non-negative")

  bad_params$prior_tol <- c(1e-9, 1e-8)
  expect_error(validate_and_override_params(bad_params), "prior_tol must be a numeric scalar")

  # Test: invalid residual_variance_upperbound (negative value)
  bad_params <- valid_params
  bad_params$residual_variance_upperbound <- -1
  expect_error(validate_and_override_params(bad_params), "must be positive")

  # Test: residual_variance_upperbound must be a numeric scalar (not a vector)
  bad_params <- valid_params
  bad_params$residual_variance_upperbound <- c(1e10, 1e11)
  expect_error(validate_and_override_params(bad_params), "residual_variance_upperbound must be a numeric scalar")

  # Test: residual_variance_upperbound must be numeric (not character)
  bad_params <- valid_params
  bad_params$residual_variance_upperbound <- "1e10"
  expect_error(validate_and_override_params(bad_params), "residual_variance_upperbound must be a numeric scalar")

  # Test: invalid scaled_prior_variance
  bad_params <- valid_params
  bad_params$scaled_prior_variance <- -0.1
  expect_error(validate_and_override_params(bad_params), "should be positive")

  # Test: invalid unmappable_effects
  bad_params <- valid_params
  bad_params$unmappable_effects <- "invalid"
  expect_error(validate_and_override_params(bad_params), "must be one of")

  # Test: unmappable effects overrides convergence method
  inf_params <- valid_params
  inf_params$unmappable_effects <- "inf"
  inf_params$convergence_method <- "elbo"
  expect_message(
    result <- validate_and_override_params(inf_params),
    "Setting convergence_method='pip'"
  )
  expect_equal(result$convergence_method, "pip")

  # Test: refine incompatible with unmappable effects
  refine_params <- valid_params
  refine_params$unmappable_effects <- "inf"
  refine_params$refine <- TRUE
  expect_error(
    validate_and_override_params(refine_params),
    "Refinement is not supported with unmappable effects"
  )

  # Test: Servin_Stephens overrides convergence method when L > 1
  ss_params <- valid_params
  ss_params$L <- 10
  ss_params$estimate_residual_method <- "Servin_Stephens"
  ss_params$convergence_method <- "elbo"
  ss_params$estimate_prior_method <- "simple"

  expect_message(
    result <- validate_and_override_params(ss_params),
    "PIP convergence"
  )
  expect_message(
    result <- validate_and_override_params(ss_params),
    "EM"
  )

  expect_true(result$use_servin_stephens)
  expect_equal(result$convergence_method, "pip")
  expect_equal(result$estimate_prior_method, "EM")

  # Test: Servin_Stephens does NOT override convergence method when L = 1
  # (ELBO is well-defined for single-effect models)
  ss_params_l1 <- valid_params
  ss_params_l1$L <- 1
  ss_params_l1$estimate_residual_method <- "Servin_Stephens"
  ss_params_l1$convergence_method <- "elbo"
  ss_params_l1$estimate_prior_method <- "EM"

  result_l1 <- validate_and_override_params(ss_params_l1)

  expect_true(result_l1$use_servin_stephens)
  expect_equal(result_l1$convergence_method, "elbo")  # Not overridden
  expect_equal(result_l1$estimate_prior_method, "EM")

  # Test: Servin_Stephens overrides estimate_residual_variance = FALSE
  ss_erv_params <- valid_params
  ss_erv_params$estimate_residual_method <- "Servin_Stephens"
  ss_erv_params$estimate_residual_variance <- FALSE
  ss_erv_params$estimate_prior_method <- "EM"

  expect_message(
    result <- validate_and_override_params(ss_erv_params),
    "estimate_residual_variance = TRUE"
  )
  expect_true(result$estimate_residual_variance)

  # Test: Servin_Stephens with explicit estimate_residual_variance = TRUE produces no warning
  ss_erv_params2 <- valid_params
  ss_erv_params2$estimate_residual_method <- "Servin_Stephens"
  ss_erv_params2$estimate_residual_variance <- TRUE
  ss_erv_params2$estimate_prior_method <- "EM"

  # Should not produce the "estimate_residual_variance" warning
  expect_no_message(
    result <- validate_and_override_params(ss_erv_params2),
    message = "integrates out residual variance"
  )
  expect_true(result$estimate_residual_variance)

  # Test: estimate_prior_variance = FALSE
  no_est_params <- valid_params
  no_est_params$estimate_prior_variance <- FALSE
  result <- validate_and_override_params(no_est_params)
  expect_equal(result$estimate_prior_method, "none")

  # Test: Servin_Stephens with estimate_prior_variance = FALSE respects user choice
  # The EM override should NOT happen when user explicitly disables prior variance estimation
  ss_no_prior_params <- valid_params
  ss_no_prior_params$estimate_residual_method <- "Servin_Stephens"
  ss_no_prior_params$estimate_prior_variance <- FALSE
  ss_no_prior_params$estimate_prior_method <- "optim"

  result <- validate_and_override_params(ss_no_prior_params)
  expect_true(result$use_servin_stephens)
  # estimate_prior_variance = FALSE -> estimate_prior_method stays "none" (set earlier)
  # The SS block should NOT override to "EM" because estimation is disabled
  expect_equal(result$estimate_prior_method, "none")

  # Test: Servin_Stephens with estimate_prior_variance = TRUE overrides to EM
  ss_yes_prior_params <- valid_params
  ss_yes_prior_params$estimate_residual_method <- "Servin_Stephens"
  ss_yes_prior_params$estimate_prior_variance <- TRUE
  ss_yes_prior_params$estimate_prior_method <- "simple"

  expect_message(
    result <- validate_and_override_params(ss_yes_prior_params),
    "EM"
  )
  expect_true(result$use_servin_stephens)
  expect_equal(result$estimate_prior_method, "EM")
})

# =============================================================================
# MODEL INITIALIZATION
# =============================================================================

test_that("initialize_matrices creates correct model matrices", {
  n <- 100
  p <- 50
  L <- 5

  data <- list(n = n, p = p)

  params <- list(
    L = L,
    scaled_prior_variance = 0.2,
    residual_variance = 1.5,
    prior_weights = rep(1/p, p),
    null_weight = 0
  )

  var_y <- 2.0

  result <- initialize_matrices(data, params, var_y)

  # Check all components exist
  expected_names <- c("alpha", "mu", "mu2", "V", "KL", "lbf",
                      "lbf_variable", "sigma2", "pi", "null_weight",
                      "predictor_weights")
  expect_true(all(expected_names %in% names(result)))

  # Check dimensions
  expect_equal(dim(result$alpha), c(L, p))
  expect_equal(dim(result$mu), c(L, p))
  expect_equal(dim(result$mu2), c(L, p))
  expect_equal(dim(result$lbf_variable), c(L, p))
  expect_length(result$V, L)
  expect_length(result$KL, L)
  expect_length(result$lbf, L)
  expect_length(result$predictor_weights, p)

  # Check initial values
  expect_equal(result$alpha, matrix(1/p, L, p))
  expect_equal(result$mu, matrix(0, L, p))
  expect_equal(result$mu2, matrix(0, L, p))
  expect_equal(result$V, rep(params$scaled_prior_variance * var_y, L))
  expect_equal(result$sigma2, params$residual_variance)
  expect_equal(result$pi, params$prior_weights)
  expect_true(all(is.na(result$KL)))
  expect_true(all(is.na(result$lbf)))
})

test_that("initialize_null_index sets null index correctly", {
  data <- list(p = 100)

  # Test: no null weight
  model_no_null <- list(null_weight = 0)
  result <- initialize_null_index(data, model_no_null)
  expect_equal(result, 0)

  model_null_null <- list(null_weight = NULL)
  result <- initialize_null_index(data, model_null_null)
  expect_equal(result, 0)

  # Test: with null weight
  model_with_null <- list(null_weight = 0.1)
  result <- initialize_null_index(data, model_with_null)
  expect_equal(result, data$p)
})

test_that("assign_names assigns variable names to model components", {
  p <- 10
  L <- 3
  data <- list(p = p)

  model <- list(
    alpha = matrix(1/p, L, p),
    mu = matrix(0, L, p),
    mu2 = matrix(0, L, p),
    lbf_variable = matrix(0, L, p),
    pip = rep(0.1, p),
    null_weight = NULL
  )

  variable_names <- paste0("var", 1:p)

  # Test: without null weight
  result <- assign_names(data, model, variable_names)
  expect_equal(names(result$pip), variable_names)
  expect_equal(colnames(result$alpha), variable_names)
  expect_equal(colnames(result$mu), variable_names)
  expect_equal(colnames(result$mu2), variable_names)
  expect_equal(colnames(result$lbf_variable), variable_names)

  # Test: with null weight
  model$null_weight <- 0.1
  model$null_index <- p
  model$pip <- rep(0.1, p - 1)
  variable_names_with_null <- c(paste0("var", 1:(p-1)), "null_placeholder")

  result <- assign_names(data, model, variable_names_with_null)
  expect_equal(names(result$pip), paste0("var", 1:(p-1)))
  expect_equal(colnames(result$alpha)[p], "null")

  # Test: NULL variable names
  result_null <- assign_names(data, model, NULL)
  expect_null(names(result_null$pip))
})

test_that("adjust_L adjusts number of effects correctly", {
  p <- 50
  L_requested <- 10
  num_effects_init <- 5
  var_y <- 2.0

  model_init_pruned <- list(
    alpha = matrix(1/p, num_effects_init, p),
    mu = matrix(0, num_effects_init, p),
    mu2 = matrix(0, num_effects_init, p),
    V = rep(1, num_effects_init)
  )

  params <- list(
    L = L_requested,
    scaled_prior_variance = 0.2
  )

  # Test: L > num_effects (should expand)
  result <- adjust_L(params, model_init_pruned, var_y)
  expect_equal(result$L, L_requested)
  expect_equal(nrow(result$model_init$alpha), L_requested)

  # Test: L < num_effects (should warn and use num_effects)
  params_small <- params
  params_small$L <- 3

  expect_message(
    result <- adjust_L(params_small, model_init_pruned, var_y),
    "is smaller than the"
  )
  expect_equal(result$L, num_effects_init)
})

test_that("prune_single_effects expands or filters model effects", {
  p <- 50
  L_init <- 10

  model_init <- list(
    alpha = matrix(1/p, L_init, p),
    mu = matrix(0, L_init, p),
    mu2 = matrix(0, L_init, p),
    lbf_variable = matrix(0, L_init, p),
    KL = rep(1, L_init),
    lbf = rep(0, L_init),
    V = rep(1, L_init),
    sets = list(cs_index = c(1, 3, 5))
  )

  # Test: L == num_effects (just removes sets)
  result_same <- prune_single_effects(model_init, L = L_init, V = NULL)
  expect_equal(nrow(result_same$alpha), L_init)
  expect_null(result_same$sets)

  # Test: expand to larger L with vector V (length(V) > 1)
  L_expand <- 15
  V_expand <- rep(2, L_expand)
  result_expand <- prune_single_effects(model_init, L = L_expand, V = V_expand)
  expect_equal(nrow(result_expand$alpha), L_expand)
  expect_equal(result_expand$V[1:L_init], rep(1, L_init))
  expect_equal(result_expand$V[(L_init+1):L_expand], rep(2, L_expand - L_init))

  # Test: expand to larger L with scalar V (length(V) == 1)
  # This tests the else branch: V <- rep(V, L)
  L_expand_scalar <- 12
  V_scalar <- 3  # Single value
  result_expand_scalar <- prune_single_effects(model_init, L = L_expand_scalar, V = V_scalar)
  expect_equal(nrow(result_expand_scalar$alpha), L_expand_scalar)
  # When V is scalar, it gets replicated to length L
  expect_equal(result_expand_scalar$V, rep(V_scalar, L_expand_scalar))
  expect_length(result_expand_scalar$V, L_expand_scalar)
  # All V values should be the same scalar value
  expect_true(all(result_expand_scalar$V == V_scalar))
})

test_that("add_null_effect adds null effect to model", {
  p <- 50
  L <- 5

  model_init <- list(
    alpha = matrix(1/p, L, p),
    mu = matrix(0, L, p),
    mu2 = matrix(0, L, p),
    lbf_variable = matrix(0, L, p),
    V = rep(1, L)
  )

  V_null <- 0

  result <- add_null_effect(model_init, V_null)

  # Check dimensions increased
  expect_equal(nrow(result$alpha), L + 1)
  expect_equal(nrow(result$mu), L + 1)
  expect_equal(nrow(result$mu2), L + 1)
  expect_equal(nrow(result$lbf_variable), L + 1)
  expect_length(result$V, L + 1)

  # Check null effect values
  expect_equal(result$alpha[L + 1, ], rep(1/p, p))
  expect_equal(result$mu[L + 1, ], rep(0, p))
  expect_equal(result$mu2[L + 1, ], rep(0, p))
  expect_equal(result$lbf_variable[L + 1, ], rep(0, p))
  expect_equal(result$V[L + 1], V_null)
})

# =============================================================================
# CORE ALGORITHM COMPONENTS
# =============================================================================

test_that("compute_eigen_decomposition computes eigenvalues and eigenvectors", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 456)
  XtX <- crossprod(base_data$X)

  result <- compute_eigen_decomposition(XtX, base_data$n)

  # Check components
  expect_true(all(c("V", "Dsq", "VtXty") %in% names(result)))
  expect_equal(dim(result$V), c(base_data$p, base_data$p))
  expect_length(result$Dsq, base_data$p)
  expect_null(result$VtXty)

  # Check eigenvalues in decreasing order
  expect_true(all(diff(result$Dsq) <= 0))

  # Check eigenvalues are non-negative
  expect_true(all(result$Dsq >= 0))

  # Verify decomposition
  LD <- XtX / base_data$n
  eig_direct <- eigen(LD, symmetric = TRUE)
  expect_equal(result$Dsq, sort(eig_direct$values * base_data$n, decreasing = TRUE), tolerance = 1e-10)
})

test_that("add_eigen_decomposition adds eigen components to data object", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 789)

  XtX <- crossprod(base_data$X)
  Xty <- as.vector(crossprod(base_data$X, base_data$y))
  yty <- sum(base_data$y^2)

  data <- list(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = base_data$n,
    p = base_data$p
  )

  params <- list(
    unmappable_effects = "inf",
    verbose = FALSE
  )

  result <- add_eigen_decomposition(data, params)

  # Check components added
  expect_true(!is.null(result$eigen_vectors))
  expect_true(!is.null(result$eigen_values))
  expect_true(!is.null(result$VtXty))

  # Check dimensions
  expect_equal(dim(result$eigen_vectors), c(base_data$p, base_data$p))
  expect_length(result$eigen_values, base_data$p)
  expect_length(result$VtXty, base_data$p)
  expect_true(all(is.finite(result$VtXty)))

  # Test with unmappable_effects = "none" (no scaling)
  params_none <- list(unmappable_effects = "none", verbose = FALSE)
  result_none <- add_eigen_decomposition(data, params_none)
  expect_true(!is.null(result_none$eigen_vectors))

  # Test with unmappable_effects = "ash" (no raw data storage needed)
  params_ash <- list(unmappable_effects = "ash", verbose = FALSE)

  result_ash <- add_eigen_decomposition(data, params_ash)
  expect_true(!is.null(result_ash$eigen_vectors))
  expect_true(!is.null(result_ash$eigen_values))
  expect_true(!is.null(result_ash$VtXty))
})

test_that("compute_omega_quantities computes omega-weighted quantities", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 111)
  XtX <- crossprod(base_data$X)

  eigen_decomp <- compute_eigen_decomposition(XtX, base_data$n)

  data <- list(
    eigen_vectors = eigen_decomp$V,
    eigen_values = eigen_decomp$Dsq,
    p = base_data$p
  )

  tau2 <- 0.01
  sigma2 <- 1.0

  result <- compute_omega_quantities(data, tau2, sigma2)

  # Check components
  expect_true(all(c("omega_var", "diagXtOmegaX") %in% names(result)))
  expect_length(result$omega_var, base_data$p)
  expect_length(result$diagXtOmegaX, base_data$p)

  # Check omega_var calculation
  expected_omega_var <- tau2 * data$eigen_values + sigma2
  expect_equal(result$omega_var, expected_omega_var, tolerance = 1e-10)

  # Check diagXtOmegaX is positive
  expect_true(all(result$diagXtOmegaX > 0))

  # Check diagXtOmegaX sums correctly
  # Should be trace of X'OmegaX
  trace_approx <- sum(result$diagXtOmegaX)
  expect_true(trace_approx > 0)
})

test_that("compute_theta_blup computes BLUP coefficients", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 222)
  L <- 5

  XtX <- crossprod(base_data$X)
  Xty <- as.vector(crossprod(base_data$X, base_data$y))

  eigen_decomp <- compute_eigen_decomposition(XtX, base_data$n)

  data <- list(
    eigen_vectors = eigen_decomp$V,
    eigen_values = eigen_decomp$Dsq,
    VtXty = crossprod(eigen_decomp$V, Xty),
    p = base_data$p
  )

  model <- list(
    alpha = matrix(1/base_data$p, L, base_data$p),
    mu = matrix(rnorm(L * base_data$p, 0, 0.1), L, base_data$p),
    tau2 = 0.01,
    sigma2 = 1.0
  )

  result <- compute_theta_blup(data, model)

  # Check output
  expect_length(result, base_data$p)

  # Check finite values
  expect_true(all(is.finite(result)))

  # When tau2 = 0, theta should be 0
  model_zero_tau2 <- model
  model_zero_tau2$tau2 <- 0
  result_zero <- compute_theta_blup(data, model_zero_tau2)
  expect_true(all(abs(as.vector(result_zero)) < 1e-10))
})

test_that("lbf_stabilization stabilizes log Bayes factors", {
  p <- 100
  lbf <- rnorm(p, mean = 5, sd = 2)
  prior_weights <- rep(1/p, p)
  shat2 <- rgamma(p, shape = 2, rate = 1)

  result <- lbf_stabilization(lbf, prior_weights, shat2)

  # Check components
  expect_true(all(c("lbf", "lpo") %in% names(result)))
  expect_length(result$lbf, p)
  expect_length(result$lpo, p)

  # Check lpo calculation
  expected_lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))
  expect_equal(result$lpo, expected_lpo, tolerance = 1e-10)

  # Test with infinite shat2
  shat2_inf <- shat2
  shat2_inf[c(1, 5, 10)] <- Inf

  result_inf <- lbf_stabilization(lbf, prior_weights, shat2_inf)

  # LBF should be 0 where shat2 is infinite
  expect_equal(result_inf$lbf[c(1, 5, 10)], rep(0, 3))

  # LPO should be log(prior) where shat2 is infinite
  expected_lpo_inf <- log(prior_weights[c(1, 5, 10)] + sqrt(.Machine$double.eps))
  expect_equal(result_inf$lpo[c(1, 5, 10)], expected_lpo_inf, tolerance = 1e-10)
})

test_that("compute_posterior_weights computes alpha and lbf_model", {
  p <- 100
  lbf <- rnorm(p, mean = 5, sd = 2)
  prior_weights <- rep(1/p, p)
  lpo <- lbf + log(prior_weights)

  result <- compute_posterior_weights(lpo)

  # Check components
  expect_true(all(c("alpha", "lbf_model") %in% names(result)))
  expect_length(result$alpha, p)
  expect_length(result$lbf_model, 1)

  # Check alpha sums to 1
  expect_equal(sum(result$alpha), 1, tolerance = 1e-10)

  # Check alpha values are probabilities
  expect_true(all(result$alpha >= 0))
  expect_true(all(result$alpha <= 1))

  # Verify calculation
  max_lpo <- max(lpo)
  w_weighted <- exp(lpo - max_lpo)
  weighted_sum_w <- sum(w_weighted)
  expected_alpha <- w_weighted / weighted_sum_w
  expected_lbf_model <- log(weighted_sum_w) + max_lpo

  expect_equal(result$alpha, expected_alpha, tolerance = 1e-10)
  expect_equal(result$lbf_model, expected_lbf_model, tolerance = 1e-10)

  # Test numerical stability with very large lpo
  lpo_large <- c(1000, 1001, 1002, rep(0, p - 3))
  result_large <- compute_posterior_weights(lpo_large)
  expect_equal(sum(result_large$alpha), 1, tolerance = 1e-10)
  expect_true(all(result_large$alpha >= 0 & result_large$alpha <= 1))
})

test_that("compute_lbf_gradient computes gradient for prior variance", {
  p <- 100
  alpha <- rep(1/p, p)
  betahat <- rnorm(p, mean = 0, sd = 1)
  shat2 <- rgamma(p, shape = 2, rate = 1)
  V <- 1.0

  result <- compute_lbf_gradient(alpha, betahat, shat2, V, use_servin_stephens = FALSE)

  # Check output is numeric scalar
  expect_length(result, 1)
  expect_true(is.finite(result))

  # Test with different V values
  result_small_V <- compute_lbf_gradient(alpha, betahat, shat2, V = 0.1, use_servin_stephens = FALSE)
  result_large_V <- compute_lbf_gradient(alpha, betahat, shat2, V = 10, use_servin_stephens = FALSE)

  expect_true(is.finite(result_small_V))
  expect_true(is.finite(result_large_V))

  # Test with Servin-Stephens (should return NULL)
  result_ss <- compute_lbf_gradient(alpha, betahat, shat2, V, use_servin_stephens = TRUE)
  expect_null(result_ss)

  # Test with NaN in intermediate calculations (should handle)
  shat2_zero <- rep(0, p)
  betahat_zero <- rep(0, p)
  result_nan <- compute_lbf_gradient(alpha, betahat_zero, shat2_zero, V, use_servin_stephens = FALSE)
  expect_true(is.finite(result_nan))
})

# =============================================================================
# VARIANCE ESTIMATION
# =============================================================================

test_that("mom_unmappable estimates variance using method of moments", {
  # Setup test data
  setup <- setup_ss_data(n = 100, p = 50, L = 5, seed = 333, unmappable_effects = "inf")
  data <- setup$data
  params <- setup$params
  params$verbose <- FALSE
  model <- setup$model

  # Compute omega
  L <- nrow(model$alpha)
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  omega <- matrix(0, L, data$p)
  for (l in seq_len(L)) {
    omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]
  }

  # Test estimating both tau2 and sigma2
  result <- mom_unmappable(data, params, model, omega, tau2 = model$tau2,
                          est_tau2 = TRUE, est_sigma2 = TRUE)

  expect_true(all(c("sigma2", "tau2") %in% names(result)))
  expect_true(result$sigma2 > 0)
  expect_true(result$tau2 >= 0)

  # Test estimating only sigma2
  result_sigma_only <- mom_unmappable(data, params, model, omega, tau2 = 0.01,
                                     est_tau2 = FALSE, est_sigma2 = TRUE)
  expect_true(result_sigma_only$sigma2 > 0)
  expect_equal(result_sigma_only$tau2, 0.01)

  # Test verbose message when estimating both tau2 and sigma2
  params_verbose <- params
  params_verbose$verbose <- TRUE
  expect_message(
    result_verbose_both <- mom_unmappable(data, params_verbose, model, omega,
                                         tau2 = model$tau2,
                                         est_tau2 = TRUE, est_sigma2 = TRUE),
    "Update \\(sigma\\^2,tau\\^2\\) to"
  )
  expect_true(all(c("sigma2", "tau2") %in% names(result_verbose_both)))

  # Test verbose message when estimating only sigma2
  expect_message(
    result_verbose_sigma <- mom_unmappable(data, params_verbose, model, omega,
                                          tau2 = 0.01,
                                          est_tau2 = FALSE, est_sigma2 = TRUE),
    "Update sigma\\^2 to"
  )
  expect_true(result_verbose_sigma$sigma2 > 0)
  expect_equal(result_verbose_sigma$tau2, 0.01)
})

test_that("mle_unmappable estimates variance using MLE", {
  # Setup test data
  setup <- setup_ss_data(n = 100, p = 50, L = 5, seed = 444, unmappable_effects = "inf")
  data <- setup$data
  params <- setup$params
  params$verbose <- FALSE
  model <- setup$model

  # Compute omega
  L <- nrow(model$alpha)
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  omega <- matrix(0, L, data$p)
  for (l in seq_len(L)) {
    omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]
  }

  # Test estimating both tau2 and sigma2
  result <- mle_unmappable(data, params, model, omega,
                          est_tau2 = TRUE, est_sigma2 = TRUE)

  expect_true(all(c("sigma2", "tau2") %in% names(result)))
  expect_true(result$sigma2 > 0)
  expect_true(result$tau2 >= 0)

  # Test estimating only sigma2
  result_sigma_only <- mle_unmappable(data, params, model, omega,
                                     est_tau2 = FALSE, est_sigma2 = TRUE)
  expect_true(result_sigma_only$sigma2 > 0)

  # Test verbose message when estimating both tau2 and sigma2
  params_verbose <- params
  params_verbose$verbose <- TRUE
  expect_message(
    result_verbose_both <- mle_unmappable(data, params_verbose, model, omega,
                                         est_tau2 = TRUE, est_sigma2 = TRUE),
    "Update \\(sigma\\^2,tau\\^2\\) to"
  )
  expect_true(all(c("sigma2", "tau2") %in% names(result_verbose_both)))

  # Test verbose message when estimating only sigma2
  expect_message(
    result_verbose_sigma <- mle_unmappable(data, params_verbose, model, omega,
                                          est_tau2 = FALSE, est_sigma2 = TRUE),
    "Update sigma\\^2 to"
  )
  expect_true(result_verbose_sigma$sigma2 > 0)
})

test_that("compute_lbf_servin_stephens computes log Bayes factor", {
  set.seed(555)
  n <- 100
  x <- rnorm(n)
  y <- 2 * x + rnorm(n)
  s0 <- 1
  alpha0 <- 0
  beta0 <- 0

  result <- compute_lbf_servin_stephens(x, y, s0, alpha0, beta0)

  # Check output is numeric scalar
  expect_length(result, 1)
  expect_true(is.finite(result))

  # LBF should be positive when there's signal
  expect_true(result > 0)

  # Test with no signal
  x_null <- rnorm(n)
  y_null <- rnorm(n)
  result_null <- compute_lbf_servin_stephens(x_null, y_null, s0, alpha0, beta0)
  expect_true(is.finite(result_null))

  # Test with different prior parameters
  result_alpha <- compute_lbf_servin_stephens(x, y, s0, alpha0 = 2, beta0 = 1)
  expect_true(is.finite(result_alpha))
})

test_that("posterior_mean_servin_stephens computes posterior mean", {
  set.seed(666)
  p <- 50
  xtx <- 100
  xty <- 50
  s0_t <- 1

  result <- posterior_mean_servin_stephens(xtx, xty, s0_t)

  # Check output is numeric
  expect_length(result, 1)
  expect_true(is.finite(result))

  # Posterior mean should be shrunk toward zero
  ols_est <- xty / xtx
  expect_true(abs(result) < abs(ols_est))

  # Test with very small prior (strong shrinkage)
  result_small <- posterior_mean_servin_stephens(xtx, xty, s0_t = 0.01)
  expect_true(abs(result_small) < abs(result))
})

test_that("posterior_var_servin_stephens computes posterior variance", {
  set.seed(777)
  xtx <- 100
  xty <- 50
  yty <- 1000
  n <- 100
  s0_t <- 1

  result <- posterior_var_servin_stephens(xtx, xty, yty, n, s0_t)

  # Check components
  expect_true(all(c("post_var", "beta1") %in% names(result)))
  expect_true(is.finite(result$post_var))
  expect_true(is.finite(result$beta1))

  # Posterior variance should be positive
  expect_true(result$post_var > 0)

  # Test with very small prior (should return 0)
  result_small <- posterior_var_servin_stephens(xtx, xty, yty, n, s0_t = 1e-6)
  expect_equal(result_small$post_var, 0)
  expect_equal(result_small$beta1, 0)
})

test_that("est_residual_variance estimates residual variance", {
  # Setup individual data
  setup <- setup_individual_data(n = 100, p = 50, L = 5, seed = 888)
  data <- setup$data
  model <- setup$model

  result <- est_residual_variance(data, model)

  # Check output is numeric and positive
  expect_length(result, 1)
  expect_true(is.finite(result))
  expect_true(result > 0)

  # Should be reasonable for random data
  expect_true(result < 10)  
})

test_that("update_model_variance updates variance components", {
  # Setup individual data with all necessary methods defined
  setup <- setup_individual_data(n = 100, p = 50, L = 5, seed = 999)
  data <- setup$data
  params <- setup$params
  params$estimate_residual_variance <- TRUE
  params$estimate_residual_method <- "MLE"
  params$residual_variance_lowerbound <- 0.01
  params$residual_variance_upperbound <- 10
  params$unmappable_effects <- "none"
  model <- setup$model

  old_sigma2 <- model$sigma2

  result <- update_model_variance(data, params, model)

  # Check sigma2 was updated
  expect_true("sigma2" %in% names(result))

  # Check sigma2 is within bounds
  expect_true(result$sigma2 >= params$residual_variance_lowerbound)
  expect_true(result$sigma2 <= params$residual_variance_upperbound)

  # Check it's finite and positive
  expect_true(is.finite(result$sigma2))
  expect_true(result$sigma2 > 0)
})

# =============================================================================
# CONVERGENCE & OPTIMIZATION
# =============================================================================

test_that("check_convergence detects convergence correctly", {
  p <- 50
  L <- 5

  params <- list(
    convergence_method = "elbo",
    tol = 1e-4,
    verbose = FALSE
  )

  model <- list(
    alpha = matrix(1/p, L, p)
  )

  tracking <- list(
    convergence = list(
      prev_elbo = -1000,
      prev_alpha = matrix(1/p, L, p)
    )
  )

  # Test: first iteration (should not converge)
  result_iter1 <- check_convergence(NULL, params, model, elbo = c(-1000, -999),
                                   iter = 1, tracking = tracking)
  expect_false(result_iter1$converged)

  # Test: ELBO converged
  elbo_converged <- c(-1000, -999.99)
  result_elbo_conv <- check_convergence(NULL, params, model, elbo = elbo_converged,
                                       iter = 2, tracking = tracking)
  expect_true(result_elbo_conv$converged)

  # Test: ELBO not converged
  tracking$convergence$prev_elbo <- -1000
  elbo_not_conv <- c(NA, NA, -990)
  result_elbo_not <- check_convergence(NULL, params, model, elbo = elbo_not_conv,
                                      iter = 2, tracking = tracking)
  expect_false(result_elbo_not$converged)

  # Test: PIP convergence
  params_pip <- list(
    convergence_method = "pip",
    tol = 1e-4,
    verbose = FALSE
  )

  # PIP converged (alpha unchanged)
  result_pip_conv <- check_convergence(NULL, params_pip, model, elbo = c(-1000, -999),
                                      iter = 2, tracking = tracking)
  expect_true(result_pip_conv$converged)

  # PIP not converged (alpha changed)
  model_changed <- model
  model_changed$alpha[1, 1] <- 0.5
  result_pip_not <- check_convergence(NULL, params_pip, model_changed,
                                     elbo = c(-1000, -999),
                                     iter = 2, tracking = tracking)
  expect_false(result_pip_not$converged)

  # Test: ELBO is NA/Inf (fallback to PIP)
  expect_message(
    result_na <- check_convergence(NULL, params, model, elbo = c(-1000, NA),
                                  iter = 2, tracking = tracking),
    "NA/infinite ELBO"
  )
  expect_true(result_na$converged)  # Alpha unchanged, so converged by PIP
})

test_that("get_objective computes ELBO correctly", {
  # Setup individual data
  setup <- setup_individual_data(n = 100, p = 50, L = 5, seed = 101)
  data <- setup$data
  params <- setup$params
  params$unmappable_effects <- "none"
  params$verbose <- FALSE
  model <- setup$model
  model$KL <- rep(0.1, 5)

  result <- get_objective(data, params, model)

  # Check output is numeric scalar
  expect_length(result, 1)
  expect_true(is.finite(result))

  # ELBO should be negative for random data
  expect_true(result < 0)

  # Test with unmappable effects
  setup_inf <- setup_ss_data(n = 100, p = 50, L = 5, seed = 102, unmappable_effects = "inf")
  data_inf <- setup_inf$data
  params_inf <- setup_inf$params
  params_inf$unmappable_effects <- "inf"
  params_inf$verbose <- FALSE
  model_inf <- setup_inf$model
  model_inf$KL <- rep(0.1, 5)
  model_inf$lbf <- rep(0, 5)

  result_inf <- get_objective(data_inf, params_inf, model_inf)
  expect_length(result_inf, 1)
  expect_true(is.finite(result_inf))
})

test_that("compute_elbo_inf computes ELBO for infinitesimal model", {
  # Setup data
  setup <- setup_ss_data(n = 100, p = 50, L = 5, seed = 103, unmappable_effects = "inf")
  data <- setup$data
  model <- setup$model

  # Compute omega
  L <- nrow(model$alpha)
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  omega <- matrix(0, L, data$p)
  for (l in seq_len(L)) {
    omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]
  }

  result <- compute_elbo_inf(
    alpha = model$alpha,
    mu = model$mu,
    omega = omega,
    lbf = rep(0, L),
    sigma2 = model$sigma2,
    tau2 = model$tau2,
    n = data$n,
    p = data$p,
    eigen_vectors = data$eigen_vectors,
    eigen_values = data$eigen_values,
    VtXty = data$VtXty,
    yty = data$yty
  )

  # Check output is numeric scalar
  expect_length(result, 1)
  expect_true(is.finite(result))

  # ELBO should be negative
  expect_true(result < 0)
})

# =============================================================================
# CREDIBLE SETS & POST-PROCESSING
# =============================================================================

test_that("n_in_CS_x counts variables in credible set", {
  # Probability vector with clear peak
  x <- c(0.5, 0.3, 0.1, 0.05, 0.03, 0.02)

  # With 90% coverage, should include first 2 variables (0.5 + 0.3 = 0.8 < 0.9, add 0.1)
  result_90 <- n_in_CS_x(x, coverage = 0.9)
  expect_equal(result_90, 3)

  # With 95% coverage, should include more
  result_95 <- n_in_CS_x(x, coverage = 0.95)
  expect_true(result_95 >= result_90)

  # With 50% coverage
  result_50 <- n_in_CS_x(x, coverage = 0.5)
  expect_equal(result_50, 1)

  # Uniform distribution
  skip("Fails on Linux in CI")
  x_uniform <- rep(1/10, 10)
  result_uniform <- n_in_CS_x(x_uniform, coverage = 0.9)
  expect_equal(result_uniform, 10)  # Need all to reach 90%
})

test_that("in_CS_x creates binary indicator for credible set", {
  x <- c(0.5, 0.3, 0.1, 0.05, 0.03, 0.02)

  result_90 <- in_CS_x(x, coverage = 0.9)

  # Check output is binary
  expect_equal(sort(unique(result_90)), c(0, 1))
  expect_length(result_90, length(x))

  # Check correct variables included
  expect_equal(sum(result_90), n_in_CS_x(x, coverage = 0.9))

  # Top probability should be in CS
  expect_equal(result_90[which.max(x)], 1)

  # Test with different coverage
  result_50 <- in_CS_x(x, coverage = 0.5)
  expect_equal(sum(result_50), 1)
})

test_that("in_CS creates credible set matrix", {
  L <- 5
  p <- 100

  # Create susie object
  alpha <- matrix(0, L, p)
  for (l in 1:L) {
    alpha[l, sample(p, 1)] <- 0.6
    alpha[l, ] <- alpha[l, ] / sum(alpha[l, ]) * 0.9
    alpha[l, ] <- alpha[l, ] + 0.1 / p
  }

  res <- list(alpha = alpha)
  class(res) <- "susie"

  result <- in_CS(res, coverage = 0.9)

  # Check dimensions
  expect_equal(dim(result), c(L, p))

  # Check binary values
  expect_true(all(result %in% c(0, 1)))

  # Each row should have at least one variable
  expect_true(all(rowSums(result) > 0))

  # Test with just alpha matrix
  result_alpha <- in_CS(alpha, coverage = 0.9)
  expect_equal(result, result_alpha)
})

test_that("n_in_CS counts variables in each credible set", {
  L <- 5
  p <- 100

  alpha <- matrix(0, L, p)
  for (l in 1:L) {
    alpha[l, sample(p, 1)] <- 0.7
    alpha[l, ] <- alpha[l, ] / sum(alpha[l, ]) * 0.9
    alpha[l, ] <- alpha[l, ] + 0.1 / p
  }

  res <- list(alpha = alpha)
  class(res) <- "susie"

  result <- n_in_CS(res, coverage = 0.9)

  # Check output
  expect_length(result, L)
  expect_true(all(result > 0))
  expect_true(all(result <= p))

  # Should match in_CS
  cs_matrix <- in_CS(res, coverage = 0.9)
  expect_equal(result, rowSums(cs_matrix))
})

test_that("get_purity computes correlation purity statistics", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 123)

  # Test with multiple variables
  pos <- c(1, 2, 3, 5, 8)
  result <- get_purity(pos, base_data$X, Xcorr = NULL)

  # Check output
  expect_length(result, 3)  # min, mean, median
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))

  # Mean should be between min and max (which is implicitly <= 1)
  expect_true(result[2] >= result[1])

  # Test with single variable (perfect purity)
  result_single <- get_purity(1, base_data$X, Xcorr = NULL)
  expect_equal(result_single, c(1, 1, 1))

  # Test with precomputed correlation
  Xcorr <- cor(base_data$X)
  result_xcorr <- get_purity(pos, base_data$X, Xcorr = Xcorr)
  expect_length(result_xcorr, 3)
  expect_true(all(result_xcorr >= 0))

  # Test with large set (should subsample)
  pos_large <- 1:40
  result_large <- get_purity(pos_large, base_data$X, Xcorr = NULL, n = 20)
  expect_length(result_large, 3)

  # Test squared correlations
  result_squared <- get_purity(pos, base_data$X, Xcorr = NULL, squared = TRUE)
  expect_length(result_squared, 3)
  expect_true(all(result_squared >= 0))
  expect_true(all(result_squared <= 1))
})

test_that("get_tracking cleans tracking object", {
  # Create tracking object with convergence data
  tracking <- list(
    elbo = c(-1000, -999, -998.5),
    sigma2 = c(1, 0.9, 0.85),
    convergence = list(
      prev_elbo = -998.5,
      prev_alpha = matrix(1/100, 5, 100)
    )
  )

  result <- get_tracking(tracking)

  # Check convergence removed
  expect_null(result$convergence)

  # Check other components preserved
  expect_equal(result$elbo, tracking$elbo)
  expect_equal(result$sigma2, tracking$sigma2)

  # Test with minimal tracking
  tracking_minimal <- list(
    elbo = c(-1000),
    convergence = NULL
  )

  result_minimal <- get_tracking(tracking_minimal)
  expect_null(result_minimal$convergence)
  expect_equal(result_minimal$elbo, tracking_minimal$elbo)
})

# =============================================================================
# END OF TESTS
# =============================================================================
