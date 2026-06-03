context("Utility functions for susieR")

# ---- Fundamental building blocks ----

test_that("warning_message emits correct prefix styles and respects warn option", {
  expect_message(warning_message("Test warning"),
                 "WARNING:.*Test warning")
  expect_message(warning_message("Test warning", style = "warning"),
                 "WARNING:.*Test warning")
  expect_message(warning_message("Test hint", style = "hint"),
                 "HINT:.*Test hint")

  old_warn <- getOption("warn")
  on.exit(options(warn = old_warn), add = TRUE)
  options(warn = -1)

  expect_no_error(warning_message("Still executes", style = "warning"))
  expect_message(warning_message("Hint shows", style = "hint"),
                 "HINT:.*Hint shows")
})

test_that("safe_cor returns correct correlation and zeroes out constant columns", {
  x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  expect_equal(safe_cor(x), cor(x))

  x_const <- cbind(x, rep(1, 3))
  expect_warning(cor(x_const), "the standard deviation is zero")
  expect_silent(result <- safe_cor(x_const))

  expect_true(is.matrix(result))
  expect_equal(diag(result), c(1, 1, 1))
  expect_equal(result[3, 1], 0)
  expect_equal(result[3, 2], 0)
  expect_equal(result[1, 3], 0)
  expect_equal(result[2, 3], 0)
})

test_that("safe_cov2cor returns correct correlation and zeroes out zero-variance rows/cols", {
  cov_mat <- matrix(c(4, 2, 2, 3), nrow = 2)
  expect_equal(safe_cov2cor(cov_mat), cov2cor(cov_mat))

  cov_mat_zero <- matrix(c(0, 0, 0, 3), nrow = 2)
  expect_warning(cov2cor(cov_mat_zero))
  expect_silent(result <- safe_cov2cor(cov_mat_zero))

  expect_true(is.matrix(result))
  expect_equal(diag(result), c(1, 1))
  expect_equal(result[1, 2], 0)
  expect_equal(result[2, 1], 0)
})

test_that("is_symmetric_matrix correctly identifies symmetric matrices", {
  sym_mat <- matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), nrow = 3)
  expect_true(is_symmetric_matrix(sym_mat))

  nonsym_mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3)
  expect_false(is_symmetric_matrix(nonsym_mat))

  expect_true(is_symmetric_matrix(diag(5)))

  sparse_sym <- Matrix::Matrix(sym_mat, sparse = TRUE)
  expect_true(is_symmetric_matrix(sparse_sym))
})

test_that("apply_nonzeros applies function to nonzero elements of sparse matrix", {
  X <- Matrix::Matrix(c(1, 0, 0, 2, 0, 3, 4, 0, 5), nrow = 3, sparse = TRUE)

  result <- apply_nonzeros(X, function(x) x^2)
  expect_equal(dim(result), dim(X))
  expect_equal(sum(result == 0), sum(X == 0))
  expected <- Matrix::Matrix(c(1, 0, 0, 4, 0, 9, 16, 0, 25), nrow = 3, sparse = TRUE)
  expect_equal(as.matrix(result), as.matrix(expected))

  result2 <- apply_nonzeros(X, function(x) x * 2)
  expected2 <- Matrix::Matrix(c(2, 0, 0, 4, 0, 6, 8, 0, 10), nrow = 3, sparse = TRUE)
  expect_equal(as.matrix(result2), as.matrix(expected2))
})

test_that("compute_colSds handles dense, sparse, and zero-variance columns", {
  set.seed(42)
  X_dense <- matrix(rnorm(100), nrow = 10, ncol = 10)
  expect_equal(compute_colSds(X_dense),
               matrixStats::colSds(X_dense), tolerance = 1e-10)

  X_sparse <- Matrix::Matrix(X_dense, sparse = TRUE)
  X_sparse[abs(X_sparse) < 0.5] <- 0
  expect_equal(compute_colSds(X_sparse),
               apply(as.matrix(X_sparse), 2, sd), tolerance = 1e-10)

  # Zero-variance column
  X_const <- cbind(X_dense, rep(1, 10))
  result_const <- compute_colSds(X_const)
  expect_equal(result_const[11], 0)
  expect_equal(result_const[1:10], matrixStats::colSds(X_dense))
})

test_that("compute_colstats returns correct stats for all center/scale combinations", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 123)

  result <- compute_colstats(base_data$X, center = TRUE, scale = TRUE)
  expect_true(all(c("cm", "csd", "d") %in% names(result)))
  expect_equal(result$cm, colMeans(base_data$X), tolerance = 1e-10)
  expected_csd <- apply(base_data$X, 2, sd)
  expect_equal(result$csd, expected_csd, tolerance = 1e-10)
  X_std <- scale(base_data$X, center = TRUE, scale = TRUE)
  expect_equal(result$d, colSums(X_std^2), tolerance = 1e-8)

  result_nocenter <- compute_colstats(base_data$X, center = FALSE, scale = TRUE)
  expect_equal(result_nocenter$cm, rep(0, base_data$p))
  expect_equal(result_nocenter$csd, expected_csd, tolerance = 1e-10)

  result_noscale <- compute_colstats(base_data$X, center = TRUE, scale = FALSE)
  expect_equal(result_noscale$cm, colMeans(base_data$X), tolerance = 1e-10)
  expect_equal(result_noscale$csd, rep(1, base_data$p))

  result_neither <- compute_colstats(base_data$X, center = FALSE, scale = FALSE)
  expect_equal(result_neither$cm, rep(0, base_data$p))
  expect_equal(result_neither$csd, rep(1, base_data$p))
  expect_equal(result_neither$d, colSums(base_data$X^2), tolerance = 1e-10)

  # sd = 0 column -> csd replaced by 1
  X_zero <- cbind(base_data$X, rep(0, base_data$n))
  result_zero <- compute_colstats(X_zero, center = TRUE, scale = TRUE)
  expect_equal(result_zero$csd[base_data$p + 1], 1)

  # Sparse matrix
  X_sparse <- Matrix::Matrix(base_data$X, sparse = TRUE)
  X_sparse[abs(X_sparse) < 1] <- 0
  result_sparse <- compute_colstats(X_sparse, center = TRUE, scale = TRUE)
  expect_length(result_sparse$cm, base_data$p)
  expect_length(result_sparse$csd, base_data$p)
  expect_length(result_sparse$d, base_data$p)
})

# ---- Data processing & validation ----

test_that("check_semi_pd identifies positive (semi-)definite and indefinite matrices", {
  result_pd <- check_semi_pd(matrix(c(2, 1, 1, 2), nrow = 2), tol = 1e-10)
  expect_true(result_pd$status)
  expect_true(all(result_pd$eigenvalues >= 0))
  expect_false(is.null(attr(result_pd$matrix, "eigen")))

  result_psd <- check_semi_pd(matrix(c(1, 1, 1, 1), nrow = 2), tol = 1e-10)
  expect_true(result_psd$status)
  expect_true(any(abs(result_psd$eigenvalues) < 1e-10))

  result_neg <- check_semi_pd(matrix(c(1, 2, 2, -1), nrow = 2), tol = 1e-10)
  expect_false(result_neg$status)
  expect_true(any(result_neg$eigenvalues < 0))

  result_id <- check_semi_pd(diag(3), tol = 1e-10)
  expect_true(result_id$status)
  expect_equal(result_id$eigenvalues, rep(1, 3), tolerance = 1e-10)
})

test_that("check_projection verifies whether a vector lies in the eigenspace", {
  A <- matrix(c(4, 2, 2, 3), nrow = 2)
  b_in <- c(2, 1)

  result_in <- check_projection(A, b_in)
  expect_true(result_in$status)
  expect_true(is.na(result_in$msg))

  A_with_eigen <- A
  attr(A_with_eigen, "eigen") <- eigen(A, symmetric = TRUE)
  expect_true(check_projection(A_with_eigen, b_in)$status)
})

test_that("validate_init rejects invalid model_init and accepts valid ones", {
  p <- 50; L <- 5
  valid_init <- list(
    alpha  = matrix(1/p, L, p),
    mu     = matrix(0, L, p),
    mu2    = matrix(0, L, p),
    V      = rep(1, L),
    sigma2 = 1,
    pi     = rep(1/p, p),
    null_index = 0
  )
  class(valid_init) <- "susie"
  data <- list(n = 100, p = p)

  # Valid init passes
  expect_silent(validate_init(data, list(L = L, model_init = valid_init)))

  # Wrong class
  bad <- valid_init; class(bad) <- "lm"
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init must be a 'susie' object")

  # NA in alpha
  bad <- valid_init; bad$alpha[1, 1] <- NA
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$alpha contains NA/Inf")

  # Inf in mu
  bad <- valid_init; bad$mu[1, 1] <- Inf
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$mu contains NA/Inf")

  # alpha not a matrix
  bad <- valid_init; bad$alpha <- as.vector(bad$alpha)
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$alpha must be a matrix")

  # alpha outside [0,1]
  bad <- valid_init; bad$alpha[1, 1] <- 1.5
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "invalid values outside range")

  # mu dimension mismatch
  bad <- valid_init; bad$mu <- matrix(0, L, p - 1)
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "dimensions do not match")

  # V length mismatch
  bad <- valid_init; bad$V <- rep(1, L - 1)
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "does not equal nrow")

  # Negative V
  bad <- valid_init; bad$V[1] <- -1
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "at least one negative value")

  # Negative sigma2
  bad <- valid_init; bad$sigma2 <- -0.5
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "sigma2 is negative")

  # NULL V passes
  init_no_V <- valid_init; init_no_V$V <- NULL
  expect_silent(validate_init(data, list(L = L, model_init = init_no_V)))

  # mu2 NA/Inf
  bad <- valid_init; bad$mu2[2, 3] <- NA
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$mu2 contains NA/Inf values")

  bad <- valid_init; bad$mu2[1, 5] <- Inf
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$mu2 contains NA/Inf values")

  # V NA/Inf
  bad <- valid_init; bad$V[2] <- NA
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$V contains NA/Inf values")

  bad <- valid_init; bad$V[3] <- Inf
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$V contains NA/Inf values")

  # sigma2 NA/Inf
  bad <- valid_init; bad$sigma2 <- NA
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$sigma2 contains NA/Inf")

  bad <- valid_init; bad$sigma2 <- Inf
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$sigma2 contains NA/Inf")

  # pi NA/Inf
  bad <- valid_init; bad$pi[10] <- NA
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$pi contains NA/Inf")

  bad <- valid_init; bad$pi[5] <- Inf
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$pi contains NA/Inf")

  # mu2 / alpha dimension mismatch
  bad <- valid_init; bad$mu2 <- matrix(0, L, p - 1)
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$mu2 and model_init\\$alpha dimensions do not match")

  bad <- valid_init; bad$mu2 <- matrix(0, L + 1, p)
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$mu2 and model_init\\$alpha dimensions do not match")

  # V = 0 is valid
  bad <- valid_init; bad$V <- rep(0, L)
  expect_silent(validate_init(data, list(L = L, model_init = bad)))

  # sigma2 = 0 is valid
  bad <- valid_init; bad$sigma2 <- 0
  expect_silent(validate_init(data, list(L = L, model_init = bad)))

  # pi length mismatch
  bad <- valid_init; bad$pi <- rep(1/(p-1), p - 1)
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$pi should have the same length as the number of columns in model_init\\$alpha")

  bad <- valid_init; bad$pi <- rep(1/(p+1), p + 1)
  expect_error(validate_init(data, list(L = L, model_init = bad)),
               "model_init\\$pi should have the same length as the number of columns in model_init\\$alpha")
})

test_that("validate_init rejects non-numeric V and sigma2", {
  # logical type passes is.finite() but must fail is.numeric()
  p <- 6; L <- 2
  valid_init <- list(
    alpha = matrix(1/p, L, p), mu = matrix(0, L, p), mu2 = matrix(0, L, p),
    V = rep(1, L), sigma2 = 1, pi = rep(1/p, p), null_index = 0)
  class(valid_init) <- "susie"
  data <- list(n = 50, p = p)

  bad_V <- valid_init; bad_V$V <- c(TRUE, FALSE)
  expect_error(validate_init(data, list(L = L, model_init = bad_V)),
               "model_init\\$V must be numeric")

  bad_s2 <- valid_init; bad_s2$sigma2 <- TRUE
  expect_error(validate_init(data, list(L = L, model_init = bad_s2)),
               "model_init\\$sigma2 must be numeric")
})

test_that("convert_individual_to_ss returns a valid ss object with correct structure", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5, seed = 123)
  data <- setup$data
  params <- list(unmappable_effects = "inf", verbose = FALSE)

  ss_data <- convert_individual_to_ss(data, params)

  expect_s3_class(ss_data, "ss")
  expect_true(all(c("XtX", "Xty", "yty", "n", "p") %in% names(ss_data)))
  expect_equal(dim(ss_data$XtX), c(50, 50))
  expect_length(ss_data$Xty, 50)
  expect_length(ss_data$yty, 1)
  expect_true(ss_data$yty > 0)

  expect_equal(ss_data$X_colmeans, attr(data$X, "scaled:center"))
  expect_equal(ss_data$y_mean, data$mean_y)
  expect_equal(attr(ss_data$XtX, "d"), attr(data$X, "d"))
  expect_equal(attr(ss_data$XtX, "scaled:scale"), attr(data$X, "scaled:scale"))

  # Eigen decomposition attached for unmappable effects
  expect_false(is.null(ss_data$eigen_vectors))
  expect_false(is.null(ss_data$eigen_values))
  expect_false(is.null(ss_data$VtXty))
})

test_that("extract_prior_weights extracts and rescales prior weights", {
  p <- 100

  model_no_null <- list(pi = rep(1/p, p), null_weight = 0, null_index = 0)
  expect_equal(extract_prior_weights(model_no_null), rep(1/p, p))

  null_weight <- 0.1
  pi_vec <- c(rep((1 - null_weight)/(p - 1), p - 1), null_weight)
  model_with_null <- list(pi = pi_vec, null_weight = null_weight, null_index = p)
  result <- extract_prior_weights(model_with_null)
  expect_length(result, p - 1)
  expect_equal(sum(result), 1, tolerance = 1e-10)
  expect_equal(result, rep(1/(p-1), p-1), tolerance = 1e-10)

  # Explicit null_weight argument matches stored one
  expect_equal(extract_prior_weights(model_with_null, null_weight = null_weight),
               result)

  # null_weight = NULL (backwards compatibility)
  model_null_weight_null <- model_no_null
  model_null_weight_null$null_weight <- NULL
  expect_equal(extract_prior_weights(model_null_weight_null), rep(1/p, p))
})

test_that("reconstruct_full_weights reconstructs prior weights with null component", {
  p <- 100
  non_null_weights <- rep(1/(p-1), p-1)

  result_no_null <- reconstruct_full_weights(non_null_weights, null_weight = 0)
  expect_length(result_no_null, p - 1)
  expect_equal(sum(result_no_null), 1, tolerance = 1e-10)
  expect_equal(result_no_null, non_null_weights, tolerance = 1e-10)

  null_weight <- 0.1
  result_with_null <- reconstruct_full_weights(non_null_weights, null_weight = null_weight)
  expect_length(result_with_null, p)
  expect_equal(sum(result_with_null), 1, tolerance = 1e-10)
  expect_equal(result_with_null[p], null_weight, tolerance = 1e-10)
  expect_equal(sum(result_with_null[1:(p-1)]), 1 - null_weight, tolerance = 1e-10)

  result_null <- reconstruct_full_weights(non_null_weights, null_weight = NULL)
  expect_equal(sum(result_null), 1, tolerance = 1e-10)
})

test_that("validate_and_override_params validates and adjusts parameters", {
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
    refine = FALSE,
    alpha0 = 0.1,
    beta0 = 0.1,
    n = 100,
    L_greedy = NULL,
    greedy_lbf_cutoff = 0.1
  )

  result <- validate_and_override_params(valid_params)
  expect_equal(result$prior_tol, 1e-9)
  expect_equal(result$tol, 1e-4)
  expect_false(result$use_NIG)

  # NIG auto-tightens tol only when tol is omitted
  nig_tol_params <- valid_params
  nig_tol_params$estimate_residual_method <- "NIG"
  nig_tol_params$estimate_prior_method <- "EM"
  nig_tol_params$L <- 1
  expect_message(result <- validate_and_override_params(nig_tol_params), "tol = 1e-6")
  expect_equal(result$tol, 1e-6)

  nig_tol_params$tol <- 1e-4
  expect_no_message(result <- validate_and_override_params(nig_tol_params))
  expect_equal(result$tol, 1e-4)

  expect_equal(result$alpha0, 1 / sqrt(valid_params$n))
  expect_equal(result$beta0, 1 / sqrt(valid_params$n))

  # Invalid prior_tol
  bad <- valid_params; bad$prior_tol <- -1
  expect_error(validate_and_override_params(bad), "prior_tol must be non-negative")
  bad$prior_tol <- c(1e-9, 1e-8)
  expect_error(validate_and_override_params(bad), "prior_tol must be a numeric scalar")

  # Invalid residual_variance_upperbound
  bad <- valid_params; bad$residual_variance_upperbound <- -1
  expect_error(validate_and_override_params(bad), "must be positive")
  bad <- valid_params; bad$residual_variance_upperbound <- c(1e10, 1e11)
  expect_error(validate_and_override_params(bad), "residual_variance_upperbound must be a numeric scalar")
  bad <- valid_params; bad$residual_variance_upperbound <- "1e10"
  expect_error(validate_and_override_params(bad), "residual_variance_upperbound must be a numeric scalar")

  # Invalid scaled_prior_variance
  bad <- valid_params; bad$scaled_prior_variance <- -0.1
  expect_error(validate_and_override_params(bad), "should be positive")

  # Invalid unmappable_effects
  bad <- valid_params; bad$unmappable_effects <- "invalid"
  expect_error(validate_and_override_params(bad), "must be one of")

  # unmappable_effects overrides convergence method
  inf_params <- valid_params
  inf_params$unmappable_effects <- "inf"
  inf_params$estimate_residual_method <- "MoM"
  inf_params$convergence_method <- "elbo"
  expect_message(result <- validate_and_override_params(inf_params),
                 "Switching to PIP-based convergence")
  expect_equal(result$convergence_method, "pip")

  # refine incompatible with unmappable effects
  refine_params <- valid_params
  refine_params$unmappable_effects <- "inf"
  refine_params$refine <- TRUE
  expect_error(validate_and_override_params(refine_params),
               "Refinement is not supported with unmappable effects")

  # NIG overrides convergence method when L > 1
  nig_params <- valid_params
  nig_params$L <- 10
  nig_params$estimate_residual_method <- "NIG"
  nig_params$convergence_method <- "elbo"
  nig_params$estimate_prior_method <- "simple"
  expect_message(result <- validate_and_override_params(nig_params), "PIP-based convergence")
  expect_message(result <- validate_and_override_params(nig_params), "EM")
  expect_true(result$use_NIG)
  expect_equal(result$convergence_method, "pip")
  expect_equal(result$estimate_prior_method, "EM")

  # NIG does NOT override convergence method when L = 1
  nig_params_l1 <- valid_params
  nig_params_l1$L <- 1
  nig_params_l1$estimate_residual_method <- "NIG"
  nig_params_l1$convergence_method <- "elbo"
  nig_params_l1$estimate_prior_method <- "EM"
  result_l1 <- validate_and_override_params(nig_params_l1)
  expect_true(result_l1$use_NIG)
  expect_equal(result_l1$convergence_method, "elbo")
  expect_equal(result_l1$estimate_prior_method, "EM")

  # NIG overrides estimate_residual_variance = FALSE
  nig_erv_params <- valid_params
  nig_erv_params$estimate_residual_method <- "NIG"
  nig_erv_params$estimate_residual_variance <- FALSE
  nig_erv_params$estimate_prior_method <- "EM"
  expect_message(result <- validate_and_override_params(nig_erv_params),
                 "estimate_residual_variance = TRUE")
  expect_true(result$estimate_residual_variance)

  # NIG with explicit estimate_residual_variance = TRUE produces no warning
  nig_erv_params2 <- valid_params
  nig_erv_params2$estimate_residual_method <- "NIG"
  nig_erv_params2$estimate_residual_variance <- TRUE
  nig_erv_params2$estimate_prior_method <- "EM"
  expect_no_message(result <- validate_and_override_params(nig_erv_params2),
                    message = "integrates out residual variance")
  expect_true(result$estimate_residual_variance)

  # estimate_prior_variance = FALSE -> method becomes "none"
  no_est_params <- valid_params; no_est_params$estimate_prior_variance <- FALSE
  result <- validate_and_override_params(no_est_params)
  expect_equal(result$estimate_prior_method, "none")

  # NIG with estimate_prior_variance = FALSE: EM override must NOT happen
  nig_no_prior_params <- valid_params
  nig_no_prior_params$estimate_residual_method <- "NIG"
  nig_no_prior_params$estimate_prior_variance <- FALSE
  nig_no_prior_params$estimate_prior_method <- "optim"
  result <- validate_and_override_params(nig_no_prior_params)
  expect_true(result$use_NIG)
  expect_equal(result$estimate_prior_method, "none")

  # NIG with estimate_prior_variance = TRUE overrides to EM
  nig_yes_prior_params <- valid_params
  nig_yes_prior_params$estimate_residual_method <- "NIG"
  nig_yes_prior_params$estimate_prior_variance <- TRUE
  nig_yes_prior_params$estimate_prior_method <- "simple"
  expect_message(result <- validate_and_override_params(nig_yes_prior_params), "EM")
  expect_true(result$use_NIG)
  expect_equal(result$estimate_prior_method, "EM")

  # NIG rejects invalid alpha0/beta0 combinations
  for (bad_pair in list(
    list(alpha0 = 0,    beta0 = 0.5),
    list(alpha0 = 0,    beta0 = 0),
    list(alpha0 = -0.5, beta0 = 1),
    list(alpha0 = 1,    beta0 = -0.5),
    list(alpha0 = Inf,  beta0 = 0.1),
    list(alpha0 = NA_real_, beta0 = 0.1),
    list(alpha0 = c(0.1, 0.2), beta0 = 0.1)
  )) {
    p_bad <- valid_params
    p_bad$estimate_residual_method <- "NIG"
    p_bad$alpha0 <- bad_pair$alpha0
    p_bad$beta0  <- bad_pair$beta0
    expect_error(validate_and_override_params(p_bad), "alpha0 > 0 and beta0 > 0")
  }

  # NIG fills NULL alpha0/beta0 from n
  nig_null <- valid_params
  nig_null$estimate_residual_method <- "NIG"
  nig_null$alpha0 <- NULL; nig_null$beta0 <- NULL
  result <- validate_and_override_params(nig_null)
  expect_equal(result$alpha0, 1 / sqrt(valid_params$n))
  expect_equal(result$beta0,  1 / sqrt(valid_params$n))

  # Non-NIG path does NOT validate alpha0/beta0
  no_nig_bad <- valid_params
  no_nig_bad$estimate_residual_method <- "MLE"
  no_nig_bad$alpha0 <- 0; no_nig_bad$beta0 <- 0
  result <- validate_and_override_params(no_nig_bad)
  expect_false(result$use_NIG)
  expect_null(result$alpha0)
  expect_null(result$beta0)

  # NIG requires a valid sample size
  nig_needs_n <- valid_params
  nig_needs_n$estimate_residual_method <- "NIG"
  for (bad_n in list(NULL, 0, -5, c(100, 200))) {
    nig_needs_n$n <- bad_n
    expect_error(validate_and_override_params(nig_needs_n), "requires a valid sample size")
  }
  nig_needs_n$n <- 100
  result <- suppressMessages(validate_and_override_params(nig_needs_n))
  expect_true(result$use_NIG)
})

test_that("validate_and_override_params rejects an invalid explicit tol", {
  base <- list(
    L = 10, prior_tol = 1e-9, residual_variance_upperbound = 1e10,
    scaled_prior_variance = 0.2, unmappable_effects = "none",
    convergence_method = "elbo", estimate_prior_variance = TRUE,
    estimate_prior_method = "EM", estimate_residual_method = "MLE",
    estimate_residual_variance = TRUE, refine = FALSE,
    alpha0 = 0.1, beta0 = 0.1, n = 100,
    L_greedy = NULL, greedy_lbf_cutoff = 0.1)

  for (bad in list(-1, NA_real_, Inf, c(1e-4, 1e-5))) {
    p_bad <- base; p_bad$tol <- bad
    expect_error(validate_and_override_params(p_bad),
                 "tol must be a non-negative numeric scalar")
  }
})

test_that("validate_and_override_params validates greedy-L parameters", {
  base <- list(
    L = 10, prior_tol = 1e-9, residual_variance_upperbound = 1e10,
    scaled_prior_variance = 0.2, unmappable_effects = "none",
    convergence_method = "elbo", estimate_prior_variance = TRUE,
    estimate_prior_method = "EM", estimate_residual_method = "MLE",
    estimate_residual_variance = TRUE, refine = FALSE,
    alpha0 = 0.1, beta0 = 0.1, n = 100,
    L_greedy = NULL, greedy_lbf_cutoff = 0.1
  )

  p_noninteger <- base; p_noninteger$L_greedy <- 2.5
  expect_error(validate_and_override_params(p_noninteger),
               "L_greedy must be NULL or a positive integer")

  p_zero <- base; p_zero$L_greedy <- 0
  expect_error(validate_and_override_params(p_zero),
               "L_greedy must be NULL or a positive integer")

  # L_greedy > L -> hint and clamp to L
  p_big <- base; p_big$L_greedy <- 20
  expect_message(res_big <- validate_and_override_params(p_big), "L_greedy is greater than L")
  expect_equal(res_big$L_greedy, base$L)

  p_ok <- base; p_ok$L_greedy <- 5L
  expect_equal(validate_and_override_params(p_ok)$L_greedy, 5L)

  p_bad_cut <- base; p_bad_cut$greedy_lbf_cutoff <- "x"
  expect_error(validate_and_override_params(p_bad_cut),
               "greedy_lbf_cutoff must be a numeric scalar")

  p_na_cut <- base; p_na_cut$greedy_lbf_cutoff <- NA_real_
  expect_error(validate_and_override_params(p_na_cut),
               "greedy_lbf_cutoff must be a numeric scalar")
})

test_that("validate_and_override_params auto-creates slot_prior and hints for ash", {
  base <- list(
    L = 10, prior_tol = 1e-9, residual_variance_upperbound = 1e10,
    scaled_prior_variance = 0.2, unmappable_effects = "none",
    convergence_method = "pip", estimate_prior_variance = TRUE,
    estimate_prior_method = "EM", estimate_residual_method = "MLE",
    estimate_residual_variance = TRUE, refine = FALSE,
    alpha0 = 0.1, beta0 = 0.1, n = 100,
    L_greedy = NULL, greedy_lbf_cutoff = 0.1, slot_prior = NULL
  )

  p_ash <- base; p_ash$unmappable_effects <- "ash"
  expect_message(res_ash <- validate_and_override_params(p_ash),
                 "slot_prior was not specified")
  expect_false(is.null(res_ash$slot_prior))
  expect_s3_class(res_ash$slot_prior, "slot_prior_betabinom")

  # Explicit non-default Beta-Binomial -> no default hint
  p_bb <- base; p_bb$slot_prior <- slot_prior_betabinom(2, 3)
  expect_no_message(validate_and_override_params(p_bb),
                    message = "slot_prior was not specified")

  # Poisson slot_prior with default nu -> nu hint
  p_pois <- base; p_pois$slot_prior <- slot_prior_poisson(C = 3)
  expect_message(validate_and_override_params(p_pois),
                 "Overdispersion parameter nu not specified")
})

test_that("validate_and_override_params rejects inf + MLE and overrides inf residual variance", {
  base <- list(
    L = 10, prior_tol = 1e-9, residual_variance_upperbound = 1e10,
    scaled_prior_variance = 0.2, unmappable_effects = "none",
    convergence_method = "pip", estimate_prior_variance = TRUE,
    estimate_prior_method = "EM", estimate_residual_method = "MLE",
    estimate_residual_variance = TRUE, refine = FALSE,
    alpha0 = 0.1, beta0 = 0.1, n = 100,
    L_greedy = NULL, greedy_lbf_cutoff = 0.1
  )

  # inf + MLE is unsupported
  p_inf_mle <- base
  p_inf_mle$unmappable_effects <- "inf"
  p_inf_mle$estimate_residual_method <- "MLE"
  expect_error(suppressMessages(validate_and_override_params(p_inf_mle)),
               "'MLE' is not supported with .*unmappable_effects = 'inf'")

  # inf with estimate_residual_variance = FALSE -> hint + override to TRUE
  p_inf_erv <- base
  p_inf_erv$unmappable_effects <- "inf"
  p_inf_erv$estimate_residual_method <- "MoM"
  p_inf_erv$estimate_residual_variance <- FALSE
  expect_message(res_inf <- validate_and_override_params(p_inf_erv),
                 "requires estimate_residual_variance = TRUE")
  expect_true(res_inf$estimate_residual_variance)
})

test_that("validate_and_override_params rejects wrong-length scaled_prior_variance", {
  base_params <- list(
    prior_tol = 1e-9,
    residual_variance_upperbound = 1e4,
    scaled_prior_variance = c(0.1, 0.2, 0.3),
    L = 5,
    unmappable_effects = "none",
    slot_prior = NULL,
    L_greedy = NULL,
    greedy_lbf_cutoff = 0.1
  )
  expect_error(validate_and_override_params(base_params),
               "scalar or a vector of length L")
})

# ---- Model initialization ----

test_that("initialize_matrices creates correct model matrices", {
  n <- 100; p <- 50; L <- 5
  data <- list(n = n, p = p)
  params <- list(
    L = L,
    scaled_prior_variance = 0.2,
    residual_variance = 1.5,
    prior_weights = rep(1/p, p),
    null_weight = 0
  )
  var_y <- 2.0

  result <- initialize_matrices.default(data, params, var_y)

  expected_names <- c("alpha", "mu", "mu2", "V", "KL", "lbf",
                      "lbf_variable", "sigma2", "pi", "null_weight",
                      "predictor_weights")
  expect_true(all(expected_names %in% names(result)))
  expect_equal(dim(result$alpha), c(L, p))
  expect_equal(dim(result$mu), c(L, p))
  expect_equal(dim(result$mu2), c(L, p))
  expect_equal(dim(result$lbf_variable), c(L, p))
  expect_length(result$V, L)
  expect_length(result$KL, L)
  expect_length(result$lbf, L)
  expect_length(result$predictor_weights, p)

  expect_equal(result$alpha, matrix(1/p, L, p))
  expect_equal(result$mu, matrix(0, L, p))
  expect_equal(result$mu2, matrix(0, L, p))
  expect_equal(result$V, rep(params$scaled_prior_variance * var_y, L))
  expect_equal(result$sigma2, params$residual_variance)
  expect_equal(result$pi, params$prior_weights)
  expect_true(all(is.na(result$KL)))
  expect_true(all(is.na(result$lbf)))
})

test_that("initialize_matrices handles vector scaled_prior_variance of length L", {
  # rep(vec * var_y, L) bug would produce length L*L
  n <- 100; p <- 50; L <- 5
  spv <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  params <- list(L = L, scaled_prior_variance = spv, residual_variance = 1.5,
                 prior_weights = rep(1/p, p), null_weight = 0)
  result <- initialize_matrices.default(list(n = n, p = p), params, var_y = 2.0)
  expect_length(result$V, L)
  expect_equal(result$V, spv * 2.0)
})

test_that("initialize_matrices keeps NIG prior variance on scaled sigma units", {
  n <- 100; p <- 50; L <- 5
  params <- list(L = L, scaled_prior_variance = 0.2, residual_variance = 1.5,
                 prior_weights = rep(1/p, p), null_weight = 0, use_NIG = TRUE)
  result <- initialize_matrices.default(list(n = n, p = p), params, var_y = 20)
  expect_equal(result$V, rep(0.2, L))
})

test_that("initialize_matrices preserves vector NIG scaled prior variance", {
  n <- 100; p <- 50; L <- 5
  spv <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  params <- list(L = L, scaled_prior_variance = spv, residual_variance = 1.5,
                 prior_weights = rep(1/p, p), null_weight = 0, use_NIG = TRUE)
  result <- initialize_matrices.default(list(n = n, p = p), params, var_y = 20)
  expect_equal(result$V, spv)
})

test_that("NIG initial log BF uses scaled prior variance, not var_y-scaled V", {
  set.seed(946)
  n <- 40; p <- 10
  X <- scale(matrix(rnorm(n * p), n, p), center = TRUE, scale = FALSE)
  y <- drop(scale(10 * rnorm(n), center = TRUE, scale = FALSE))
  data <- list(n = n, p = p, X = X, y = y)
  class(data) <- "individual"
  params <- list(L = 1, scaled_prior_variance = 0.2, residual_variance = var(y),
                 prior_weights = rep(1/p, p), null_weight = 0, use_NIG = TRUE)
  model <- initialize_matrices.default(data, params, var_y = var(y))
  model$predictor_weights <- colSums(X^2)
  model$raw_residuals <- y
  model$residuals <- drop(crossprod(X, y))
  nig_ss <- get_nig_sufficient_stats(data, model)

  lbf <- compute_lbf_NIG(n, model$predictor_weights, model$residuals,
                         nig_ss$yy, nig_ss$sxy, model$V[1],
                         a0 = 0.1, b0 = 0.1)
  expected_lbf <- compute_lbf_NIG(n, model$predictor_weights, model$residuals,
                                  nig_ss$yy, nig_ss$sxy, s0 = 0.2,
                                  a0 = 0.1, b0 = 0.1)
  old_lbf <- compute_lbf_NIG(n, model$predictor_weights, model$residuals,
                             nig_ss$yy, nig_ss$sxy, s0 = 0.2 * var(y),
                             a0 = 0.1, b0 = 0.1)

  expect_equal(lbf, expected_lbf, tolerance = 1e-14)
  expect_gt(max(abs(old_lbf - expected_lbf)), 0.1)
})

test_that("NIG final model reports prior variance in y units", {
  set.seed(947)
  n <- 40; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n, sd = 3)
  scaled_prior_variance <- 0.2

  fit <- susie(X, y, L = 1, estimate_residual_method = "NIG",
               estimate_prior_variance = FALSE,
               scaled_prior_variance = scaled_prior_variance,
               verbose = FALSE)

  expect_true("rv" %in% names(fit))
  expect_equal(fit$V / fit$rv, scaled_prior_variance, tolerance = 1e-12)
  expect_gt(abs(fit$V - scaled_prior_variance), 1)
})

test_that("expand_scaled_prior_variance recycles scalar and preserves vector", {
  expect_equal(expand_scaled_prior_variance(0.2, 2.0, 5), rep(0.4, 5))
  expect_equal(expand_scaled_prior_variance(c(0.1, 0.2, 0.3, 0.4, 0.5), 2.0, 5),
               c(0.2, 0.4, 0.6, 0.8, 1.0))
})

test_that("susie with vector scaled_prior_variance runs end-to-end", {
  # Before fix, rep(vec * var_y, L) produced length L*L and poisoned downstream state
  set.seed(1)
  n <- 200; p <- 100
  beta <- rep(0, p); beta[1:4] <- 1
  X <- scale(matrix(rnorm(n * p), nrow = n), center = TRUE, scale = TRUE)
  y <- drop(X %*% beta + rnorm(n))

  fit <- susie(X, y, L = 10, estimate_prior_variance = FALSE,
               scaled_prior_variance = rep(1, 10))
  expect_length(fit$V, 10)
  expect_true(all(is.finite(fit$V)))
})

test_that("vector scaled_prior_variance composes with model_init L expansion", {
  set.seed(2)
  n <- 200; p <- 80
  beta <- rep(0, p); beta[1:3] <- 1
  X <- scale(matrix(rnorm(n * p), nrow = n), center = TRUE, scale = TRUE)
  y <- drop(X %*% beta + rnorm(n))

  init <- susie(X, y, L = 2, estimate_prior_variance = TRUE)
  spv <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  fit <- susie(X, y, L = 5, estimate_prior_variance = FALSE,
               scaled_prior_variance = spv, model_init = init)
  expect_length(fit$V, 5)
  expect_true(all(is.finite(fit$V)))
})

test_that("initialize_null_index sets null index correctly", {
  data <- list(p = 100)

  expect_equal(initialize_null_index(data, list(null_weight = 0)), 0)
  expect_equal(initialize_null_index(data, list(null_weight = NULL)), 0)
  expect_equal(initialize_null_index(data, list(null_weight = 0.1)), data$p)
})

test_that("assign_names assigns variable names to model components", {
  p <- 10; L <- 3
  data <- list(p = p)
  model <- list(
    alpha       = matrix(1/p, L, p),
    mu          = matrix(0, L, p),
    mu2         = matrix(0, L, p),
    lbf_variable = matrix(0, L, p),
    pip         = rep(0.1, p),
    null_weight = NULL
  )
  variable_names <- paste0("var", 1:p)

  result <- assign_names(data, model, variable_names)
  expect_equal(names(result$pip), variable_names)
  expect_equal(colnames(result$alpha), variable_names)
  expect_equal(colnames(result$mu), variable_names)
  expect_equal(colnames(result$mu2), variable_names)
  expect_equal(colnames(result$lbf_variable), variable_names)

  # With null weight
  model$null_weight <- 0.1
  model$null_index <- p
  model$pip <- rep(0.1, p - 1)
  variable_names_with_null <- c(paste0("var", 1:(p-1)), "null_placeholder")
  result <- assign_names(data, model, variable_names_with_null)
  expect_equal(names(result$pip), paste0("var", 1:(p-1)))
  expect_equal(colnames(result$alpha)[p], "null")

  # NULL names
  result_null <- assign_names(data, model, NULL)
  expect_null(names(result_null$pip))
})

test_that("adjust_L expands or warns when num_effects exceeds requested L", {
  p <- 50; var_y <- 2.0
  model_init_pruned <- list(
    alpha = matrix(1/p, 5, p), mu = matrix(0, 5, p),
    mu2   = matrix(0, 5, p), V = rep(1, 5))
  params <- list(L = 10, scaled_prior_variance = 0.2)

  result <- adjust_L(params, model_init_pruned, var_y)
  expect_equal(result$L, 10)
  expect_equal(nrow(result$model_init$alpha), 10)

  params_small <- list(L = 3, scaled_prior_variance = 0.2)
  expect_message(result <- adjust_L(params_small, model_init_pruned, var_y),
                 "is smaller than the")
  expect_equal(result$L, 5)
})

test_that("prune_single_effects expands or filters model effects", {
  p <- 50; L_init <- 10
  model_init <- list(
    alpha        = matrix(1/p, L_init, p),
    mu           = matrix(0, L_init, p),
    mu2          = matrix(0, L_init, p),
    lbf_variable = matrix(0, L_init, p),
    KL   = rep(1, L_init),
    lbf  = rep(0, L_init),
    V    = rep(1, L_init),
    sets = list(cs_index = c(1, 3, 5))
  )

  # L == num_effects removes sets
  result_same <- prune_single_effects(model_init, L = L_init, V = NULL)
  expect_equal(nrow(result_same$alpha), L_init)
  expect_null(result_same$sets)

  # Expand with vector V: new rows use V_expand values
  L_expand <- 15; V_expand <- rep(2, L_expand)
  result_expand <- prune_single_effects(model_init, L = L_expand, V = V_expand)
  expect_equal(nrow(result_expand$alpha), L_expand)
  expect_equal(result_expand$V[1:L_init], rep(1, L_init))
  expect_equal(result_expand$V[(L_init+1):L_expand], rep(2, L_expand - L_init))

  # Expand with scalar V: replicated to length L
  V_scalar <- 3
  result_scalar <- prune_single_effects(model_init, L = 12, V = V_scalar)
  expect_equal(nrow(result_scalar$alpha), 12)
  expect_equal(result_scalar$V, rep(V_scalar, 12))
})

test_that("add_null_effect adds null effect row to model matrices", {
  p <- 50; L <- 5
  model_init <- list(
    alpha        = matrix(1/p, L, p),
    mu           = matrix(0, L, p),
    mu2          = matrix(0, L, p),
    lbf_variable = matrix(0, L, p),
    V            = rep(1, L))

  result <- add_null_effect(model_init, 0)

  expect_equal(nrow(result$alpha), L + 1)
  expect_equal(nrow(result$mu), L + 1)
  expect_equal(nrow(result$mu2), L + 1)
  expect_equal(nrow(result$lbf_variable), L + 1)
  expect_length(result$V, L + 1)
  expect_equal(result$alpha[L + 1, ], rep(1/p, p))
  expect_equal(result$mu[L + 1, ], rep(0, p))
  expect_equal(result$mu2[L + 1, ], rep(0, p))
  expect_equal(result$lbf_variable[L + 1, ], rep(0, p))
  expect_equal(result$V[L + 1], 0)
})

# ---- Core algorithm components ----

test_that("compute_eigen_decomposition computes eigenvalues and eigenvectors", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 456)
  XtX <- crossprod(base_data$X)

  result <- compute_eigen_decomposition(XtX, base_data$n)
  expect_true(all(c("V", "Dsq", "VtXty") %in% names(result)))
  expect_equal(dim(result$V), c(base_data$p, base_data$p))
  expect_length(result$Dsq, base_data$p)
  expect_null(result$VtXty)
  expect_true(all(diff(result$Dsq) <= 0))
  expect_true(all(result$Dsq >= 0))

  LD <- XtX / base_data$n
  eig_direct <- eigen(LD, symmetric = TRUE)
  expect_equal(result$Dsq, sort(eig_direct$values * base_data$n, decreasing = TRUE),
               tolerance = 1e-10)
})

test_that("compute_eigen_decomposition handles X-SVD, Cholesky+SVD, and legacy paths", {
  set.seed(1005)

  # Thin SVD path: X supplied directly (n x p with n < p)
  n <- 30; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  res_a <- compute_eigen_decomposition(NULL, n, X = X)
  expect_equal(dim(res_a$V), c(p, min(n, p)))
  expect_length(res_a$Dsq, min(n, p))
  expect_null(res_a$VtXty)
  expect_true(all(diff(res_a$Dsq) <= 1e-9))
  expect_true(all(res_a$Dsq >= 0))

  # Pivoted Cholesky + SVD path: rank-deficient XtX with n + 1 < p
  XtX <- crossprod(X)
  res_b <- suppressWarnings(compute_eigen_decomposition(XtX, n))
  expect_equal(nrow(res_b$V), p)
  expect_null(res_b$VtXty)
  ev_direct <- sort(eigen(XtX, symmetric = TRUE)$values, decreasing = TRUE)
  expect_equal(res_b$Dsq[1:5], ev_direct[1:5], tolerance = 1e-6)
  expect_true(all(res_b$Dsq >= 0))

  # Legacy full eigen path: p <= n + 1
  n2 <- 60; p2 <- 20
  X2 <- matrix(rnorm(n2 * p2), n2, p2)
  XtX2 <- crossprod(X2)
  res_c <- compute_eigen_decomposition(XtX2, n2)
  expect_equal(dim(res_c$V), c(p2, p2))
  expect_null(res_c$VtXty)
  ev2 <- sort(eigen(XtX2 / n2, symmetric = TRUE)$values * n2, decreasing = TRUE)
  expect_equal(res_c$Dsq, ev2, tolerance = 1e-8)
})

test_that("add_eigen_decomposition adds eigen components to data object", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 789)
  XtX <- crossprod(base_data$X)
  Xty <- as.vector(crossprod(base_data$X, base_data$y))
  data <- list(XtX = XtX, Xty = Xty, yty = sum(base_data$y^2),
               n = base_data$n, p = base_data$p)

  result <- add_eigen_decomposition(data, list(unmappable_effects = "inf", verbose = FALSE))
  expect_false(is.null(result$eigen_vectors))
  expect_false(is.null(result$eigen_values))
  expect_false(is.null(result$VtXty))
  expect_equal(dim(result$eigen_vectors), c(base_data$p, base_data$p))
  expect_length(result$eigen_values, base_data$p)
  expect_length(result$VtXty, base_data$p)
  expect_true(all(is.finite(result$VtXty)))

  result_none <- add_eigen_decomposition(data, list(unmappable_effects = "none", verbose = FALSE))
  expect_false(is.null(result_none$eigen_vectors))

  result_ash <- add_eigen_decomposition(data, list(unmappable_effects = "ash", verbose = FALSE))
  expect_false(is.null(result_ash$eigen_vectors))
  expect_false(is.null(result_ash$VtXty))
})

test_that("compute_omega_quantities computes omega-weighted quantities", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 111)
  XtX <- crossprod(base_data$X)
  eigen_decomp <- compute_eigen_decomposition(XtX, base_data$n)
  data <- list(eigen_vectors = eigen_decomp$V, eigen_values = eigen_decomp$Dsq,
               p = base_data$p)

  result <- compute_omega_quantities(data, tau2 = 0.01, sigma2 = 1.0)
  expect_true(all(c("omega_var", "diagXtOmegaX") %in% names(result)))
  expect_length(result$omega_var, base_data$p)
  expect_length(result$diagXtOmegaX, base_data$p)
  expect_equal(result$omega_var, 0.01 * data$eigen_values + 1.0, tolerance = 1e-10)
  expect_true(all(result$diagXtOmegaX > 0))
})

test_that("compute_theta_blup computes BLUP coefficients", {
  set.seed(222)
  base_data <- generate_base_data(n = 100, p = 50, seed = 222)
  L <- 5
  XtX <- crossprod(base_data$X)
  Xty <- as.vector(crossprod(base_data$X, base_data$y))
  eigen_decomp <- compute_eigen_decomposition(XtX, base_data$n)
  data <- list(eigen_vectors = eigen_decomp$V, eigen_values = eigen_decomp$Dsq,
               VtXty = crossprod(eigen_decomp$V, Xty), p = base_data$p)
  model <- list(alpha = matrix(1/base_data$p, L, base_data$p),
                mu    = matrix(rnorm(L * base_data$p, 0, 0.1), L, base_data$p),
                tau2 = 0.01, sigma2 = 1.0)

  result <- compute_theta_blup(data, model)
  expect_length(result, base_data$p)
  expect_true(all(is.finite(result)))

  model$tau2 <- 0
  expect_true(all(abs(as.vector(compute_theta_blup(data, model))) < 1e-10))
})

test_that("lbf_stabilization stabilizes log Bayes factors and zeroes infinite shat2", {
  set.seed(42)
  p <- 100
  lbf <- rnorm(p, mean = 5, sd = 2)
  prior_weights <- rep(1/p, p)
  shat2 <- rgamma(p, shape = 2, rate = 1)

  result <- lbf_stabilization(lbf, prior_weights, shat2)
  expect_true(all(c("lbf", "lpo") %in% names(result)))
  expected_lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))
  expect_equal(result$lpo, expected_lpo, tolerance = 1e-10)

  shat2_inf <- shat2; shat2_inf[c(1, 5, 10)] <- Inf
  result_inf <- lbf_stabilization(lbf, prior_weights, shat2_inf)
  expect_equal(result_inf$lbf[c(1, 5, 10)], rep(0, 3))
  expect_equal(result_inf$lpo[c(1, 5, 10)],
               log(prior_weights[c(1, 5, 10)] + sqrt(.Machine$double.eps)),
               tolerance = 1e-10)
})

test_that("compute_posterior_weights returns valid alpha summing to 1", {
  set.seed(42)
  p <- 100
  lbf <- rnorm(p, mean = 5, sd = 2)
  prior_weights <- rep(1/p, p)
  lpo <- lbf + log(prior_weights)

  result <- compute_posterior_weights(lpo)
  expect_true(all(c("alpha", "lbf_model") %in% names(result)))
  expect_length(result$alpha, p)
  expect_equal(sum(result$alpha), 1, tolerance = 1e-10)
  expect_true(all(result$alpha >= 0 & result$alpha <= 1))

  max_lpo <- max(lpo)
  w <- exp(lpo - max_lpo)
  expect_equal(result$alpha, w / sum(w), tolerance = 1e-10)
  expect_equal(result$lbf_model, log(sum(w)) + max_lpo, tolerance = 1e-10)

  # Numerical stability with very large lpo
  lpo_large <- c(1000, 1001, 1002, rep(0, p - 3))
  result_large <- compute_posterior_weights(lpo_large)
  expect_equal(sum(result_large$alpha), 1, tolerance = 1e-10)
  expect_true(all(result_large$alpha >= 0 & result_large$alpha <= 1))
})

test_that("compute_lbf_gradient returns finite scalar or NULL for NIG", {
  set.seed(42)
  p <- 100
  alpha <- rep(1/p, p)
  betahat <- rnorm(p)
  shat2 <- rgamma(p, shape = 2, rate = 1)

  result <- compute_lbf_gradient(alpha, betahat, shat2, V = 1.0, use_NIG = FALSE)
  expect_length(result, 1)
  expect_true(is.finite(result))
  expect_true(is.finite(compute_lbf_gradient(alpha, betahat, shat2, V = 0.1, use_NIG = FALSE)))
  expect_true(is.finite(compute_lbf_gradient(alpha, betahat, shat2, V = 10,  use_NIG = FALSE)))

  expect_null(compute_lbf_gradient(alpha, betahat, shat2, V = 1.0, use_NIG = TRUE))

  # Degenerate input (shat2 = 0) must not crash
  result_nan <- compute_lbf_gradient(alpha, rep(0, p), rep(0, p), V = 1.0, use_NIG = FALSE)
  expect_true(is.finite(result_nan))
})

# ---- Variance estimation ----

test_that("mom_unmappable estimates sigma2 and tau2 via method of moments", {
  setup <- setup_ss_data(n = 100, p = 50, L = 5, seed = 333, unmappable_effects = "inf")
  data <- setup$data; params <- setup$params; model <- setup$model
  params$verbose <- FALSE

  L <- nrow(model$alpha)
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  omega <- matrix(0, L, data$p)
  for (l in seq_len(L)) omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]

  result <- mom_unmappable(data, params, model, omega, tau2 = model$tau2,
                           est_tau2 = TRUE, est_sigma2 = TRUE)
  expect_true(all(c("sigma2", "tau2") %in% names(result)))
  expect_true(result$sigma2 > 0)
  expect_true(result$tau2 >= 0)

  result_sigma_only <- mom_unmappable(data, params, model, omega, tau2 = 0.01,
                                      est_tau2 = FALSE, est_sigma2 = TRUE)
  expect_true(result_sigma_only$sigma2 > 0)
  expect_equal(result_sigma_only$tau2, 0.01)

  params_verbose <- params; params_verbose$verbose <- TRUE
  expect_message(
    mom_unmappable(data, params_verbose, model, omega, tau2 = model$tau2,
                   est_tau2 = TRUE, est_sigma2 = TRUE),
    "Update \\(sigma\\^2,tau\\^2\\) to")
  expect_message(
    mom_unmappable(data, params_verbose, model, omega, tau2 = 0.01,
                   est_tau2 = FALSE, est_sigma2 = TRUE),
    "Update sigma\\^2 to")
})

test_that("mle_unmappable estimates sigma2 and tau2 via MLE", {
  setup <- setup_ss_data(n = 100, p = 50, L = 5, seed = 444, unmappable_effects = "inf")
  data <- setup$data; params <- setup$params; model <- setup$model
  params$verbose <- FALSE

  L <- nrow(model$alpha)
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  omega <- matrix(0, L, data$p)
  for (l in seq_len(L)) omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]

  result <- mle_unmappable(data, params, model, omega, est_tau2 = TRUE, est_sigma2 = TRUE)
  expect_true(all(c("sigma2", "tau2") %in% names(result)))
  expect_true(result$sigma2 > 0)
  expect_true(result$tau2 >= 0)

  result_sigma_only <- mle_unmappable(data, params, model, omega,
                                      est_tau2 = FALSE, est_sigma2 = TRUE)
  expect_true(result_sigma_only$sigma2 > 0)

  params_verbose <- params; params_verbose$verbose <- TRUE
  expect_message(
    mle_unmappable(data, params_verbose, model, omega, est_tau2 = TRUE, est_sigma2 = TRUE),
    "Update \\(sigma\\^2,tau\\^2\\) to")
  expect_message(
    mle_unmappable(data, params_verbose, model, omega, est_tau2 = FALSE, est_sigma2 = TRUE),
    "Update sigma\\^2 to")
})

test_that("mle_unmappable warns and retains previous params when optim fails", {
  # Degenerate data (yty = 0) collapses the sigma^2 search range, triggering
  # the non-convergence warning branch that retains previous parameters.
  set.seed(1007)
  setup <- setup_ss_data(n = 100, p = 50, L = 5, seed = 444, unmappable_effects = "inf")
  data <- setup$data; params <- setup$params; model <- setup$model
  params$verbose <- FALSE

  L <- nrow(model$alpha)
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  omega <- matrix(0, L, data$p)
  for (l in seq_len(L)) omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]

  data_bad <- data; data_bad$yty <- 0
  expect_message(
    result <- mle_unmappable(data_bad, params, model, omega,
                             est_tau2 = TRUE, est_sigma2 = TRUE),
    "failed to converge")
  expect_equal(result$sigma2, model$sigma2)
  expect_equal(result$tau2, model$tau2)
})

test_that("compute_lbf_NIG_univariate returns finite positive LBF for strong signal", {
  set.seed(555)
  n <- 100
  x <- rnorm(n)
  y <- 2 * x + rnorm(n)

  result <- compute_lbf_NIG_univariate(x, y, s0 = 1, alpha0 = 0, beta0 = 0)
  expect_length(result, 1)
  expect_true(is.finite(result))
  expect_true(result > 0)

  # Null data -> finite result
  expect_true(is.finite(compute_lbf_NIG_univariate(rnorm(n), rnorm(n), s0 = 1,
                                                    alpha0 = 0, beta0 = 0)))
  # Different prior parameters
  expect_true(is.finite(compute_lbf_NIG_univariate(x, y, s0 = 1, alpha0 = 2, beta0 = 1)))
})

test_that("get_nig_sufficient_stats matches cor() formula after intercept centering", {
  set.seed(556)
  n <- 40; p <- 5
  X <- scale(matrix(rnorm(n * p, mean = 2), n, p), center = TRUE, scale = FALSE)
  r <- drop(scale(3 + rnorm(n), center = TRUE, scale = FALSE))
  model <- list(raw_residuals = r, residuals = drop(crossprod(X, r)),
                predictor_weights = colSums(X^2))

  out <- get_nig_sufficient_stats(list(X = X), model)
  expected <- drop(crossprod(X, r)) / sqrt(colSums(X^2) * sum(r^2))
  expect_equal(out$sxy, expected, tolerance = 1e-14)
  expect_equal(out$sxy, drop(cor(X, r)), tolerance = 1e-14)
})

test_that("get_nig_sufficient_stats differs from cor() when intercept is not centered", {
  set.seed(557)
  n <- 40; p <- 5
  X <- matrix(rnorm(n * p, mean = 2), n, p)
  r <- 3 + rnorm(n)
  model <- list(raw_residuals = r, residuals = drop(crossprod(X, r)),
                predictor_weights = colSums(X^2))

  out <- get_nig_sufficient_stats(list(X = X), model)
  expected <- drop(crossprod(X, r)) / sqrt(colSums(X^2) * sum(r^2))
  old_cor_formula <- drop(cor(X, r))

  expect_equal(out$sxy, expected, tolerance = 1e-14)
  expect_gt(max(abs(out$sxy - old_cor_formula)), 0.1)

  new_lbf <- compute_lbf_NIG(n, colSums(X^2), drop(crossprod(X, r)),
                             sum(r^2), out$sxy, s0 = 1, a0 = 0.1, b0 = 0.1)
  old_lbf <- compute_lbf_NIG(n, colSums(X^2), drop(crossprod(X, r)),
                             sum(r^2), old_cor_formula, s0 = 1, a0 = 0.1, b0 = 0.1)
  expected_lbf <- compute_lbf_NIG(n, colSums(X^2), drop(crossprod(X, r)),
                                  sum(r^2), expected, s0 = 1, a0 = 0.1, b0 = 0.1)
  expect_equal(new_lbf, expected_lbf, tolerance = 1e-14)
  expect_gt(max(abs(old_lbf - expected_lbf)), 1)
})

test_that("NIG split helpers match bundled formula reference", {
  compute_stats_NIG_reference <- function(n, xx, xy, yy, sxy, s0, a0, b0, tau = 1) {
    r0   <- s0 / (s0 + tau / xx)
    rss  <- yy * (1 - r0 * sxy^2)
    a1   <- a0 + n
    b1   <- b0 + rss
    lbf  <- -(log(1 + s0 * xx / tau) + a1 * log(b1 / (b0 + yy))) / 2
    bhat <- xy / xx
    post_mean <- r0 * bhat
    post_var  <- b1 / (a1 - 2) * r0 * tau / xx
    rv        <- (b1 / 2) / (a1 / 2 - 1)
    list(lbf = lbf, post_mean = post_mean, post_mean2 = post_var + post_mean^2,
         post_var = post_var, rv = rv)
  }

  set.seed(666)
  n <- 40; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  xx <- colSums(X^2); xy <- drop(crossprod(X, y)); yy <- sum(y^2)
  sxy <- xy / sqrt(xx * yy)
  s0 <- 0.8; a0 <- 0.1; b0 <- 0.2
  tau <- seq(1, 1.5, length.out = p)

  ref     <- compute_stats_NIG_reference(n, xx, xy, yy, sxy, s0, a0, b0, tau)
  lbf     <- compute_lbf_NIG(n, xx, xy, yy, sxy, s0, a0, b0, tau)
  moments <- compute_posterior_moments_NIG(n, xx, xy, yy, sxy, s0, a0, b0, tau)

  expect_equal(lbf,              ref$lbf,        tolerance = 1e-14)
  expect_equal(moments$post_mean,  ref$post_mean,  tolerance = 1e-14)
  expect_equal(moments$post_mean2, ref$post_mean2, tolerance = 1e-14)
  expect_equal(moments$post_var,   ref$post_var,   tolerance = 1e-14)
  expect_equal(moments$rv,         ref$rv,         tolerance = 1e-14)
})

test_that("est_residual_variance returns a finite positive estimate", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5, seed = 888)
  result <- est_residual_variance(setup$data, setup$model)
  expect_length(result, 1)
  expect_true(is.finite(result))
  expect_true(result > 0)
})

test_that("est_residual_variance stops on negative estimate", {
  # Degenerate SS (yty = 0 with large cross term) drives get_ER2 negative
  setup <- setup_ss_data(n = 100, p = 20, L = 3, seed = 7, unmappable_effects = "none")
  data <- setup$data; model <- setup$model

  data$yty <- 0; data$Xty[1] <- 1e6
  model$alpha <- matrix(0, 3, 20); model$alpha[1, 1] <- 1
  model$mu    <- matrix(0, 3, 20); model$mu[1, 1]    <- 1
  model$mu2   <- matrix(0, 3, 20); model$mu2[1, 1]   <- 1

  expect_error(est_residual_variance(data, model),
               "est_residual_variance\\(\\) failed: the estimated value is negative")
})

test_that("update_model_variance updates sigma2 within bounds", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5, seed = 999)
  data <- setup$data; params <- setup$params; model <- setup$model
  params$estimate_residual_variance <- TRUE
  params$estimate_residual_method <- "MLE"
  params$residual_variance_lowerbound <- 0.01
  params$residual_variance_upperbound <- 10
  params$unmappable_effects <- "none"

  result <- update_model_variance.default(data, params, model)
  expect_true("sigma2" %in% names(result))
  expect_true(result$sigma2 >= params$residual_variance_lowerbound)
  expect_true(result$sigma2 <= params$residual_variance_upperbound)
  expect_true(is.finite(result$sigma2) && result$sigma2 > 0)
})

# ---- NIG helpers - NULL and active branches ----

test_that("compute_null_loglik_NIG returns NULL for non-NIG and correct value for NIG", {
  expect_null(compute_null_loglik_NIG(n = 100, yy = 5, a0 = 0.1, b0 = 0.1, use_NIG = FALSE))

  val <- compute_null_loglik_NIG(n = 100, yy = 5, a0 = 0.2, b0 = 0.3, use_NIG = TRUE)
  expect_length(val, 1)
  expect_true(is.finite(val))
  expected <- -100 * log(2 * pi) / 2 +
    inv_gamma_factor(0.2 / 2, 0.3 / 2) -
    inv_gamma_factor((0.2 + 100) / 2, (0.3 + 5) / 2)
  expect_equal(val, expected, tolerance = 1e-12)
})

test_that("compute_marginal_loglik returns NULL for non-NIG and lbf + null for NIG", {
  expect_null(compute_marginal_loglik(lbf_model = 2, n = 100, yy = 5,
                                      a0 = 0.1, b0 = 0.1, use_NIG = FALSE))

  ll0 <- compute_null_loglik_NIG(n = 80, yy = 4, a0 = 0.1, b0 = 0.1, use_NIG = TRUE)
  val <- compute_marginal_loglik(lbf_model = 3, n = 80, yy = 4,
                                 a0 = 0.1, b0 = 0.1, use_NIG = TRUE)
  expect_equal(val, 3 + ll0, tolerance = 1e-12)
})

# ---- Convergence & optimization ----

test_that("check_convergence detects ELBO and PIP convergence correctly", {
  p <- 50; L <- 5
  params <- list(convergence_method = "elbo", tol = 1e-4, verbose = FALSE)
  model <- list(
    alpha = matrix(1/p, L, p),
    runtime = list(prev_elbo = -1000, prev_alpha = matrix(1/p, L, p),
                   prev_pip_diff = NULL)
  )

  # First iteration must not converge
  expect_false(check_convergence.default(NULL, params, model,
                                         elbo = c(-1000, -999), iter = 1)$converged)

  # ELBO converged
  expect_true(check_convergence.default(NULL, params, model,
                                        elbo = c(-1000, -999.99), iter = 2)$converged)

  # ELBO not converged
  model$runtime$prev_elbo <- -1000
  expect_false(check_convergence.default(NULL, params, model,
                                         elbo = c(NA, NA, -990), iter = 2)$converged)

  # PIP converged (alpha unchanged)
  params_pip <- list(convergence_method = "pip", tol = 1e-4, verbose = FALSE)
  expect_true(check_convergence.default(NULL, params_pip, model,
                                        elbo = c(-1000, -999), iter = 2)$converged)

  # PIP not converged (alpha changed)
  model_changed <- model; model_changed$alpha[1, 1] <- 0.5
  expect_false(check_convergence.default(NULL, params_pip, model_changed,
                                         elbo = c(-1000, -999), iter = 2)$converged)

  # NA ELBO -> fallback to PIP
  expect_message(
    result_na <- check_convergence.default(NULL, params, model,
                                           elbo = c(-1000, NA), iter = 2),
    "NA/infinite ELBO")
  expect_true(result_na$converged)
})

test_that("PIP convergence detects and averages short alpha cycles", {
  alpha_a <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
  alpha_b <- matrix(c(0.1, 0.9, 0.8, 0.2), nrow = 2, byrow = TRUE)
  model <- list(
    alpha   = alpha_a,
    runtime = list(prev_alpha = alpha_b,
                   alpha_history = list(alpha_a, alpha_b),
                   pip_history   = list(susie_get_pip(alpha_a), susie_get_pip(alpha_b)))
  )
  params <- list(tol = 1e-4, pip_stall_window = 5, prior_tol = 1e-9)

  result <- check_alpha_pip_cycle_convergence(NULL, params, model)
  expect_true(result$converged)
  expect_equal(result$convergence_reason, "alpha_pip_cycle_2")
  expect_equal(result$alpha, (alpha_a + alpha_b) / 2)
})

test_that("get_objective computes finite ELBO for individual and unmappable data", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5, seed = 101)
  data <- setup$data; params <- setup$params; model <- setup$model
  params$unmappable_effects <- "none"; params$verbose <- FALSE
  model$KL <- rep(0.1, 5)

  result <- get_objective(data, params, model)
  expect_length(result, 1)
  expect_true(is.finite(result) && result < 0)

  setup_inf <- setup_ss_data(n = 100, p = 50, L = 5, seed = 102, unmappable_effects = "inf")
  data_inf <- setup_inf$data; params_inf <- setup_inf$params; model_inf <- setup_inf$model
  params_inf$unmappable_effects <- "inf"; params_inf$verbose <- FALSE
  model_inf$KL <- rep(0.1, 5); model_inf$lbf <- rep(0, 5)

  result_inf <- get_objective(data_inf, params_inf, model_inf)
  expect_length(result_inf, 1)
  expect_true(is.finite(result_inf))
})

test_that("compute_elbo_inf computes finite negative ELBO for infinitesimal model", {
  setup <- setup_ss_data(n = 100, p = 50, L = 5, seed = 103, unmappable_effects = "inf")
  data <- setup$data; model <- setup$model

  L <- nrow(model$alpha)
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  omega <- matrix(0, L, data$p)
  for (l in seq_len(L)) omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]

  result <- compute_elbo_inf(
    alpha = model$alpha, mu = model$mu, omega = omega, lbf = rep(0, L),
    sigma2 = model$sigma2, tau2 = model$tau2, n = data$n, p = data$p,
    eigen_vectors = data$eigen_vectors, eigen_values = data$eigen_values,
    VtXty = data$VtXty, yty = data$yty)

  expect_length(result, 1)
  expect_true(is.finite(result) && result < 0)
})

# ---- Credible sets & post-processing ----

test_that("n_in_CS_x counts variables in credible set", {
  x <- c(0.5, 0.3, 0.1, 0.05, 0.03, 0.02)
  expect_equal(n_in_CS_x(x, coverage = 0.9), 3)
  expect_true(n_in_CS_x(x, coverage = 0.95) >= n_in_CS_x(x, coverage = 0.9))
  expect_equal(n_in_CS_x(x, coverage = 0.5), 1)

  skip("Fails on Linux in CI")
  expect_equal(n_in_CS_x(rep(1/10, 10), coverage = 0.9), 10)
})

test_that("in_CS_x creates binary indicator for credible set", {
  x <- c(0.5, 0.3, 0.1, 0.05, 0.03, 0.02)

  result_90 <- in_CS_x(x, coverage = 0.9)
  expect_equal(sort(unique(result_90)), c(0, 1))
  expect_length(result_90, length(x))
  expect_equal(sum(result_90), n_in_CS_x(x, coverage = 0.9))
  expect_equal(result_90[which.max(x)], 1)
  expect_equal(sum(in_CS_x(x, coverage = 0.5)), 1)
})

test_that("in_CS creates binary credible set matrix from susie object or alpha", {
  set.seed(42)
  L <- 5; p <- 100
  alpha <- matrix(0, L, p)
  for (l in 1:L) {
    alpha[l, sample(p, 1)] <- 0.6
    alpha[l, ] <- alpha[l, ] / sum(alpha[l, ]) * 0.9 + 0.1 / p
  }
  res <- structure(list(alpha = alpha), class = "susie")

  result <- in_CS(res, coverage = 0.9)
  expect_equal(dim(result), c(L, p))
  expect_true(all(result %in% c(0, 1)))
  expect_true(all(rowSums(result) > 0))

  # alpha matrix dispatch
  expect_equal(in_CS(alpha, coverage = 0.9), result)
})

test_that("n_in_CS counts variables per credible set", {
  set.seed(42)
  L <- 5; p <- 100
  alpha <- matrix(0, L, p)
  for (l in 1:L) {
    alpha[l, sample(p, 1)] <- 0.7
    alpha[l, ] <- alpha[l, ] / sum(alpha[l, ]) * 0.9 + 0.1 / p
  }
  res <- structure(list(alpha = alpha), class = "susie")

  result <- n_in_CS(res, coverage = 0.9)
  expect_length(result, L)
  expect_true(all(result > 0 & result <= p))
  expect_equal(result, rowSums(in_CS(res, coverage = 0.9)))
})

test_that("get_purity computes valid purity statistics", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 123)
  pos <- c(1, 2, 3, 5, 8)

  result <- get_purity(pos, base_data$X, Xcorr = NULL)
  expect_length(result, 3)
  expect_true(all(result >= 0 & result <= 1))
  expect_true(result[2] >= result[1])

  # Single variable -> perfect purity
  expect_equal(get_purity(1, base_data$X, Xcorr = NULL), c(1, 1, 1))

  # Precomputed correlation
  Xcorr <- cor(base_data$X)
  result_xcorr <- get_purity(pos, base_data$X, Xcorr = Xcorr)
  expect_length(result_xcorr, 3)
  expect_true(all(result_xcorr >= 0))

  # Large set with subsampling
  result_large <- get_purity(1:40, base_data$X, Xcorr = NULL, n = 20)
  expect_length(result_large, 3)

  # Squared correlations
  result_sq <- get_purity(pos, base_data$X, Xcorr = NULL, squared = TRUE)
  expect_length(result_sq, 3)
  expect_true(all(result_sq >= 0 & result_sq <= 1))
})

test_that("get_purity errors when correlations contain NaN/NA", {
  Xc <- matrix(0.5, 5, 5); diag(Xc) <- 1
  Xc[1, 2] <- NA; Xc[2, 1] <- NA
  expect_error(get_purity(c(1, 2, 3), X = NULL, Xcorr = Xc),
               "get_purity returned NaN/NA")
})

# ---- Small matrix helpers ----

test_that("compute_Rv dispatches across X / XtX / R / Rv_matrix and errors", {
  set.seed(1001)
  n <- 50; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  v <- rnorm(p)
  XtX <- crossprod(X); R <- cor(X)

  expect_equal(compute_Rv(list(X = X), v),
               as.vector(crossprod(X, X %*% v)), tolerance = 1e-10)
  expect_equal(compute_Rv(list(X = NULL, XtX = XtX), v),
               as.vector(XtX %*% v), tolerance = 1e-10)
  expect_equal(compute_Rv(list(X = NULL, XtX = NULL, R = R), v),
               as.vector(R %*% v), tolerance = 1e-10)

  M <- matrix(rnorm(n * p), n, p)
  expect_equal(compute_Rv(list(X = X), v, Rv_matrix = M),
               as.vector(crossprod(M, M %*% v)), tolerance = 1e-10)

  expect_error(compute_Rv(list(a = 1), v), "No predictor matrix available")
})

test_that("compute_BR dispatches across X / XtX and errors", {
  set.seed(1002)
  n <- 40; p <- 8; L <- 3
  X <- matrix(rnorm(n * p), n, p)
  B <- matrix(rnorm(L * p), L, p)
  XtX <- crossprod(X)

  expect_equal(compute_BR(list(X = X), B), (B %*% t(X)) %*% X, tolerance = 1e-10)
  expect_equal(compute_BR(list(X = NULL, XtX = XtX), B), B %*% XtX, tolerance = 1e-10)
  expect_error(compute_BR(list(a = 1), B),
               "No predictor matrix available for compute_BR")
})

test_that("compute_XtXv_eigen reconstructs XtX %*% v from eigen factors", {
  set.seed(1003)
  n <- 60; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  XtX <- crossprod(X)
  v <- rnorm(p)
  ev <- eigen(XtX, symmetric = TRUE)
  data <- list(eigen_vectors = ev$vectors, eigen_values = ev$values)

  result <- compute_XtXv_eigen(data, v)
  expect_length(result, p)
  expect_equal(result, as.vector(XtX %*% v), tolerance = 1e-8)
})

test_that("scale_design_matrix centers and scales with NULL defaults", {
  set.seed(1004)
  n <- 30; p <- 6
  X <- matrix(rnorm(n * p), n, p)

  expect_equal(scale_design_matrix(X), X, tolerance = 1e-12)

  cm  <- colMeans(X)
  csd <- apply(X, 2, sd)
  ref_cs <- scale(X, center = cm, scale = csd)
  attributes(ref_cs) <- list(dim = dim(ref_cs))
  expect_equal(scale_design_matrix(X, center = cm, scale = csd), ref_cs, tolerance = 1e-10)

  ref_c <- scale(X, center = cm, scale = FALSE)
  attributes(ref_c) <- list(dim = dim(ref_c))
  expect_equal(scale_design_matrix(X, center = cm), ref_c, tolerance = 1e-10)
})

test_that("format_V_summary formats non-zero, zero, and NA entries", {
  expect_equal(format_V_summary(c(0.123, 0, 0, NA, 0.05)),
               "[1.23e-01, 5.00e-02, 0 x 2, NA x 1]")
  expect_equal(format_V_summary(c(0.1, 0.2)), "[1.00e-01, 2.00e-01]")
  expect_equal(format_V_summary(c(0, 0, 0)), "[0 x 3]")
})

# ---- get_xcorr - correlation-matrix dispatch and caching ----

test_that("get_xcorr dispatches across cache / SS / individual paths and errors", {
  set.seed(1006)
  n <- 60; p <- 10
  X <- matrix(rnorm(n * p), n, p)

  # Individual path: safe_cor(X), result cached
  res_ind <- get_xcorr(list(X = X))
  expect_equal(dim(res_ind$Xcorr), c(p, p))
  expect_equal(res_ind$Xcorr, safe_cor(X), tolerance = 1e-10)
  expect_false(is.null(res_ind$data$Xcorr_cache))

  # Cached path
  res_cached <- get_xcorr(list(X = X, Xcorr_cache = diag(p)))
  expect_equal(res_cached$Xcorr, diag(p))

  # SS path with covariance diagonal
  XtX <- crossprod(scale(X, center = TRUE, scale = FALSE))
  res_ss <- get_xcorr(list(XtX = XtX))
  expect_equal(res_ss$Xcorr, safe_cov2cor(XtX), tolerance = 1e-10)

  # SS path when XtX is already a correlation matrix
  R <- safe_cov2cor(XtX)
  expect_equal(get_xcorr(list(XtX = R))$Xcorr, R)

  # Error: neither XtX nor X
  expect_error(get_xcorr(list(a = 1)), "neither XtX nor X")
})

# ---- SuSiE-ash end-to-end ----

test_that("susie ash mode (default path) runs and produces valid structure", {
  set.seed(1008)
  sim <- simulate_regression(n = 120, p = 40, k = 3, signal_sd = 3)

  fit <- suppressWarnings(susie(
    sim$X, sim$y, L = 8, slot_prior = slot_prior_betabinom(),
    unmappable_effects = "ash", max_iter = 30, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_false(is.null(fit$alpha))
  expect_false(is.null(fit$theta))
  expect_false(is.null(fit$tau2))
  expect_length(fit$theta, sim$p)
  expect_true(all(is.finite(fit$theta)))
  expect_true(is.finite(fit$sigma2) && fit$sigma2 > 0)
  expect_false(is.null(fit$pip))
  expect_length(fit$pip, sim$p)
})

test_that("susie ash mode on unmappable_data drives masking decision tiers", {
  # Real unmappable data with genuine LD runs many iterations, exercising
  # cumulative masking state inside update_ash_variance_components().
  data(unmappable_data)
  set.seed(1009)
  X <- unmappable_data$X[, 1:300]
  y <- as.vector(unmappable_data$y)

  fit <- suppressWarnings(susie(
    X, y, L = 10, slot_prior = slot_prior_betabinom(),
    unmappable_effects = "ash", max_iter = 30, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_length(fit$theta, ncol(X))
  expect_true(all(is.finite(fit$theta)))
  expect_true(fit$niter >= 2)
  expect_true(is.finite(fit$tau2))
  expect_true(all(abs(rowSums(fit$alpha) - 1) < 1e-6))
})

test_that("susie ash_filter_archived mode exercises compute_ash_masking", {
  data(unmappable_data)
  set.seed(1010)
  X <- unmappable_data$X[, 1:300]
  y <- as.vector(unmappable_data$y)

  fit <- suppressWarnings(susie(
    X, y, L = 10, unmappable_effects = "ash_filter_archived",
    max_iter = 30, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_false(is.null(fit$theta))
  expect_length(fit$theta, ncol(X))
  expect_true(all(is.finite(fit$theta)))
  expect_true(is.finite(fit$tau2))
})

test_that("susie_ss ash_filter_archived mode runs through summary-stat masking", {
  # SS path: get_xcorr derives correlation from XtX and masking dispatches to
  # compute_ash_from_summary_stats via update_ash_variance_components_filter_archived.
  set.seed(1011)
  n <- 150; p <- 40
  Z <- matrix(rnorm(n * 8), n, 8)
  X <- matrix(0, n, p)
  for (b in 1:8) {
    idx <- ((b - 1) * 5 + 1):(b * 5)
    X[, idx] <- Z[, b] + matrix(rnorm(n * 5, sd = 0.3), n, 5)
  }
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p); beta[c(3, 12, 23)] <- 4
  y <- drop(X %*% beta + rnorm(n)); y <- y - mean(y)

  fit <- suppressWarnings(susie_ss(
    XtX = crossprod(X), Xty = as.vector(crossprod(X, y)),
    yty = sum(y^2), n = n, L = 8,
    unmappable_effects = "ash_filter_archived", max_iter = 30, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_false(is.null(fit$theta))
  expect_length(fit$theta, p)
  expect_true(all(is.finite(fit$theta)))
})

# ---- SuSiE-ash masking decision tiers - direct internal-function drivers ----
#
# The end-to-end ash fits above exercise the masking machinery, but specific
# sub-branches (wait-then-expose, oscillation reversal, second-chance restore,
# lazy ash_iter init) require precise per-effect states. We drive the internal
# updaters directly with hand-crafted alpha matrices and an LD structure that
# guarantees those states, using options(susie.skip_mrash = TRUE) so
# branch-hitting is deterministic and fast.

# Build a small SS data object with two tight-LD blocks plus a moderate-LD
# partner:
#   - cols 1,2: tight LD (|r| ~ 0.999) -> high-purity CS (CASE 3)
#   - col 3:    moderate LD with block A -> CS purity ~ 0.35 (CASE 2)
#   - cols 7,8: tight LD (|r| ~ 0.999) -> second clean signal / collision
.make_ash_collision_ss <- function(n = 400, p = 12, seed = 1) {
  set.seed(seed)
  Z1 <- rnorm(n); Z2 <- rnorm(n)
  X <- matrix(rnorm(n * p), n, p)
  X[, 1] <- Z1 + rnorm(n, sd = 0.02)
  X[, 2] <- Z1 + rnorm(n, sd = 0.02)
  X[, 3] <- 0.35 * Z1 + rnorm(n, sd = sqrt(1 - 0.35^2))
  X[, 7] <- Z2 + rnorm(n, sd = 0.02)
  X[, 8] <- Z2 + rnorm(n, sd = 0.02)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p); beta[1] <- 5; beta[7] <- 5
  y <- drop(X %*% beta + rnorm(n)); y <- y - mean(y)
  data <- structure(
    list(XtX = crossprod(X), X = NULL, Xty = as.vector(crossprod(X, y)),
         yty = sum(y^2), n = n, p = p, X_colmeans = rep(0, p), y_mean = 0),
    class = "ss")
  list(data = data, Xcorr = stats::cor(X), p = p)
}

# An alpha row putting `mass` on `idx` and spreading the rest uniformly.
.ash_alpha_row <- function(idx, mass, p) {
  a <- rep((1 - sum(mass)) / (p - length(idx)), p)
  a[idx] <- mass
  a
}

# Carry the mutable masking-state fields from a result back onto the model so
# the next call sees accumulated state (mimics the per-iteration loop).
.ash_carry_state <- function(model, res) {
  for (nm in c("ash_iter", "prev_case", "prev_sentinel", "ever_diffuse",
               "diffuse_iter_count", "masked", "ever_unmasked",
               "unmask_candidate_iters", "force_exposed_iter",
               "second_chance_used")) {
    model[[nm]] <- res[[nm]]
  }
  model$theta  <- res$theta
  model$sigma2 <- res$sigma2
  model$tau2   <- res$tau2
  model$ash_pi <- res$ash_pi
  model
}

test_that("update_ash_variance_components lazily initializes ash_iter", {
  # The defensive guard fires only when a caller passes a model without the field
  d <- .make_ash_collision_ss(seed = 21)
  L <- 2; p <- d$p
  model <- list(
    alpha = matrix(1/p, L, p), mu = matrix(0, L, p), mu2 = matrix(0, L, p),
    V = rep(1, L), lbf = rep(0, L), sigma2 = 1, theta = rep(0, p),
    ash_pi = NULL, tau2 = 0)   # NOTE: no ash_iter field
  model$alpha[1, ] <- .ash_alpha_row(c(1, 2), c(0.5, 0.49), p)
  params <- list(verbose = FALSE, estimate_residual_variance = TRUE, slot_prior = NULL)

  old <- getOption("susie.skip_mrash", FALSE)
  on.exit(options(susie.skip_mrash = old), add = TRUE)
  options(susie.skip_mrash = TRUE)

  res <- update_ash_variance_components(d$data, model, params)
  expect_equal(res$ash_iter, 1L)
})

test_that("update_ash_variance_components drives CASE 2 wait-then-expose and second-chance", {
  # A persistent moderate-purity effect (CASE 2) accumulates diffuse_iter_count;
  # after enough stable iters its neighborhood is exposed (force_exposed_iter set),
  # and later the second-chance window restores masking.
  d <- .make_ash_collision_ss(seed = 22)
  L <- 3; p <- d$p
  model <- list(
    alpha = matrix(1/p, L, p), mu = matrix(0.5, L, p), mu2 = matrix(0.5, L, p),
    V = rep(1, L), lbf = rep(0, L), sigma2 = 1, theta = rep(0, p),
    ash_pi = NULL, tau2 = 0, ash_iter = 0L)
  params <- list(verbose = FALSE, estimate_residual_variance = TRUE, slot_prior = NULL)

  old <- getOption("susie.skip_mrash", FALSE)
  on.exit(options(susie.skip_mrash = old), add = TRUE)
  options(susie.skip_mrash = TRUE)

  exposed_seen <- FALSE; restored_seen <- FALSE
  for (it in 1:7) {
    model$alpha[1, ] <- .ash_alpha_row(c(1, 3), c(0.5, 0.45), p)
    model$alpha[2, ] <- .ash_alpha_row(c(7, 8), c(0.5, 0.49), p)
    res <- update_ash_variance_components(d$data, model, params)
    if (any(res$force_exposed_iter > 0)) exposed_seen <- TRUE
    if (any(res$second_chance_used))     restored_seen <- TRUE
    model <- .ash_carry_state(model, res)
  }

  expect_true(exposed_seen)
  expect_true(restored_seen)
  expect_true(all(is.finite(model$theta)))
})

test_that("update_ash_variance_components reverses subtraction on CASE 2->3 oscillation", {
  # When a slot flips CASE 2 -> CASE 3 it is flagged unstable and the
  # c_hat-weighted subtraction added this iter is reversed.
  d <- .make_ash_collision_ss(seed = 23)
  L <- 3; p <- d$p
  model <- list(
    alpha = matrix(1/p, L, p), mu = matrix(2, L, p), mu2 = matrix(4, L, p),
    V = rep(1, L), lbf = rep(0, L), sigma2 = 1, theta = rep(0, p),
    ash_pi = NULL, tau2 = 0, ash_iter = 0L)
  params <- list(verbose = FALSE, estimate_residual_variance = TRUE, slot_prior = NULL)

  old <- getOption("susie.skip_mrash", FALSE)
  on.exit(options(susie.skip_mrash = old), add = TRUE)
  options(susie.skip_mrash = TRUE)

  model$alpha[1, ] <- .ash_alpha_row(c(1, 3), c(0.5, 0.45), p)
  res <- update_ash_variance_components(d$data, model, params)
  expect_equal(res$prev_case[1], 2L)
  model <- .ash_carry_state(model, res)

  model$alpha[1, ] <- .ash_alpha_row(c(1, 2), c(0.5, 0.49), p)
  res <- update_ash_variance_components(d$data, model, params)
  expect_equal(res$prev_case[1], 3L)
  expect_gt(res$ever_diffuse[1], 0)
})

test_that("compute_ash_masking drives all decision tiers (archived classifier)", {
  d <- .make_ash_collision_ss(seed = 24)
  L <- 3; p <- d$p
  params <- list(verbose = FALSE)

  base_model <- function() {
    m <- list(alpha = matrix(1/p, L, p), mu = matrix(2, L, p),
              mu2 = matrix(4, L, p), V = rep(1, L), lbf = rep(0, L), sigma2 = 1)
    init_ash_fields_filter_archived(m, n = d$data$n, p = p, L = L, is_individual = FALSE)
  }

  # wait-then-expose + second-chance: persistent CASE 2 slot 1
  m <- base_model()
  exposed_seen <- FALSE; restored_seen <- FALSE
  for (it in 1:7) {
    m$alpha[1, ] <- .ash_alpha_row(c(1, 3), c(0.5, 0.45), p)
    m$alpha[2, ] <- .ash_alpha_row(c(7, 8), c(0.5, 0.49), p)
    r <- compute_ash_masking(d$Xcorr, m, params)
    m <- r$model
    if (any(m$force_exposed_iter > 0)) exposed_seen <- TRUE
    if (any(m$second_chance_used))     restored_seen <- TRUE
  }
  expect_true(exposed_seen)
  expect_true(restored_seen)

  # oscillation CASE 2 -> CASE 3
  m2 <- base_model()
  m2$alpha[1, ] <- .ash_alpha_row(c(1, 3), c(0.5, 0.45), p)
  r1 <- compute_ash_masking(d$Xcorr, m2, params)
  expect_equal(r1$current_case[1], 2)
  m2 <- r1$model
  m2$alpha[1, ] <- .ash_alpha_row(c(1, 2), c(0.5, 0.49), p)
  r2 <- compute_ash_masking(d$Xcorr, m2, params)
  expect_equal(r2$current_case[1], 3)
  expect_gt(r2$model$ever_diffuse[1], 0)
})

test_that("update_ash_variance_components_filter_archived runs in skip_mrash mode", {
  d <- .make_ash_collision_ss(seed = 25)
  L <- 2; p <- d$p
  model <- list(alpha = matrix(1/p, L, p), mu = matrix(1, L, p),
                mu2 = matrix(1, L, p), V = rep(1, L), lbf = rep(0, L), sigma2 = 1)
  model <- init_ash_fields_filter_archived(model, n = d$data$n, p = p, L = L,
                                           is_individual = FALSE)
  params <- list(verbose = FALSE, estimate_residual_variance = TRUE, slot_prior = NULL)

  old <- getOption("susie.skip_mrash", FALSE)
  on.exit(options(susie.skip_mrash = old), add = TRUE)
  options(susie.skip_mrash = TRUE)

  for (it in 1:6) {
    model$alpha[1, ] <- .ash_alpha_row(c(1, 3), c(0.5, 0.45), p)
    model$alpha[2, ] <- .ash_alpha_row(c(7, 8), c(0.5, 0.49), p)
    res <- update_ash_variance_components_filter_archived(d$data, model, params)
    model <- .ash_carry_state(model, res)
  }
  expect_equal(model$ash_iter, 6)
  expect_true(all(is.finite(model$theta)))
  expect_true(is.finite(model$tau2))
})

test_that("update_ash_variance_components: collision resets diffuse_iter_count", {
  # CASE 2 path with current_collision[l] = TRUE:
  # both sentinels highly correlated forces collision; ever_diffuse > 0 forces CASE 2.
  set.seed(1100)
  n <- 200; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  X[, 2] <- X[, 1] + rnorm(n, 0, 0.1)
  X[, 3] <- X[, 1] * 0.5 + rnorm(n, 0, sqrt(1 - 0.25))
  X[, 4] <- X[, 2] * 0.5 + rnorm(n, 0, sqrt(1 - 0.25))
  X <- scale(X)

  data <- structure(list(X = X, n = n, p = p, XtX = crossprod(X), y_mean = 0),
                    class = c("individual", "list"))

  L <- 2
  alpha <- matrix(0, L, p)
  alpha[1, 1] <- 0.55; alpha[1, 3] <- 0.44; alpha[1, 5:p] <- 0.01 / (p - 4)
  alpha[2, 2] <- 0.55; alpha[2, 4] <- 0.44; alpha[2, 5:p] <- 0.01 / (p - 4)
  alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
  alpha[2, ] <- alpha[2, ] / sum(alpha[2, ])

  options(susie.skip_mrash = TRUE)
  on.exit(options(susie.skip_mrash = NULL), add = TRUE)

  model <- list(
    alpha               = alpha,
    mu                  = matrix(0.05, L, p),
    sigma2              = 1,
    V                   = rep(0.1, L),
    ash_iter            = 1L,
    prev_case           = rep(2L, L),
    prev_sentinel       = rep(0L, L),
    ever_diffuse        = rep(1L, L),
    diffuse_iter_count  = rep(3L, L),
    masked              = rep(FALSE, p),
    ever_unmasked       = rep(FALSE, p),
    unmask_candidate_iters = rep(0L, p),
    force_exposed_iter  = rep(0L, p),
    second_chance_used  = rep(FALSE, p),
    slot_weights        = rep(1, L),
    theta               = rep(0, p),
    tau2                = 0,
    ash_pi              = rep(1/p, p),
    ash_s0              = rep(0.1, p),
    X_theta             = rep(0, n)
  )

  result <- update_ash_variance_components(data, model, list())

  # After a collision both slots remain CASE 2 and diffuse_iter_count resets to 0
  expect_equal(result$diffuse_iter_count, rep(0L, L))
  expect_equal(result$prev_case, rep(2L, L))
})
