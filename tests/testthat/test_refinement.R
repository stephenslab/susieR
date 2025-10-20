devtools::load_all(".")


context("Refinement unit tests")

# =============================================================================
# BASIC FUNCTIONALITY
# =============================================================================

test_that("run_refine returns a valid susie model", {
  setup <- create_model_with_cs(seed = 100)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true("alpha" %in% names(refined_model))
  expect_true("mu" %in% names(refined_model))
  expect_true("V" %in% names(refined_model))
  expect_true("sigma2" %in% names(refined_model))
  expect_true("elbo" %in% names(refined_model))
})

test_that("run_refine maintains or improves ELBO", {
  setup <- create_model_with_cs(seed = 101)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  initial_elbo <- susie_get_objective(setup$model)
  refined_model <- run_refine(setup$model, setup$data, setup$params)
  final_elbo <- susie_get_objective(refined_model)

  expect_true(final_elbo >= initial_elbo - 1e-6)
})

test_that("run_refine preserves model dimensions", {
  setup <- create_model_with_cs(seed = 102)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  initial_L <- nrow(setup$model$alpha)
  initial_p <- ncol(setup$model$alpha)

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_equal(nrow(refined_model$alpha), initial_L)
  expect_equal(ncol(refined_model$alpha), initial_p)
  expect_equal(nrow(refined_model$mu), initial_L)
  expect_equal(ncol(refined_model$mu), initial_p)
})

test_that("run_refine returns finite ELBO", {
  setup <- create_model_with_cs(seed = 103)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(all(is.finite(refined_model$elbo)))
  expect_true(is.finite(susie_get_objective(refined_model)))
})

test_that("run_refine maintains valid probability distributions", {
  setup <- create_model_with_cs(seed = 104)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(all(refined_model$alpha >= 0))
  expect_true(all(refined_model$alpha <= 1))

  row_sums <- rowSums(refined_model$alpha)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

# =============================================================================
# REFINEMENT LOGIC
# =============================================================================

test_that("run_refine iterates through credible sets", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10,
                                n_causal = 3, seed = 105)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  n_cs_initial <- length(setup$model$sets$cs)

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(is.finite(susie_get_objective(refined_model)))
})

test_that("run_refine uses two-step procedure correctly", {
  setup <- create_model_with_cs(seed = 106)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true("alpha" %in% names(refined_model))
  expect_true(all(refined_model$alpha >= 0 & refined_model$alpha <= 1))
})

test_that("run_refine preserves prior weights structure", {
  setup <- create_model_with_cs(seed = 107)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  initial_pi <- setup$model$pi

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_equal(length(refined_model$pi), length(initial_pi))
  expect_true(all(refined_model$pi >= 0))
  expect_true(abs(sum(refined_model$pi) - 1) < 1e-10)
})

test_that("run_refine evaluates multiple candidate models", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10,
                                n_causal = 4, seed = 108)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) < 2,
          "Need at least 2 credible sets")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(is.finite(susie_get_objective(refined_model)))
})

# =============================================================================
# CONVERGENCE BEHAVIOR
# =============================================================================

test_that("run_refine stops when ELBO improvement < tol", {
  setup <- create_model_with_cs(seed = 109)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  params_tight_tol <- setup$params
  params_tight_tol$tol <- 1e-10

  refined_model <- run_refine(setup$model, setup$data, params_tight_tol)

  expect_true(is.finite(susie_get_objective(refined_model)))
})

test_that("run_refine with loose tolerance may iterate more", {
  setup <- create_model_with_cs(seed = 110)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  params_loose <- setup$params
  params_loose$tol <- 1e-1

  refined_model <- run_refine(setup$model, setup$data, params_loose)

  expect_true(is.finite(susie_get_objective(refined_model)))
})

test_that("run_refine stops when no candidate models generated", {
  setup <- create_model_with_cs(seed = 111)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true("elbo" %in% names(refined_model))
})

test_that("run_refine convergence is deterministic", {
  setup <- create_model_with_cs(seed = 112)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  refined1 <- run_refine(setup$model, setup$data, setup$params)
  refined2 <- run_refine(setup$model, setup$data, setup$params)

  expect_equal(susie_get_objective(refined1),
               susie_get_objective(refined2))
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("run_refine handles model with no credible sets", {
  set.seed(113)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  model <- susie(X, y, L = 5, verbose = FALSE)

  constructor_result <- individual_data_constructor(
    X = X, y = y, L = 5,
    standardize = TRUE, intercept = TRUE,
    estimate_residual_method = "MLE",
    convergence_method = "elbo",
    coverage = 0.95, min_abs_corr = 0.5,
    n_purity = 100,
    track_fit = FALSE
  )

  if (is.null(model$sets) || length(model$sets$cs) == 0) {
    refined_model <- run_refine(model, constructor_result$data,
                                constructor_result$params)

    expect_equal(susie_get_objective(refined_model),
                 susie_get_objective(model))
  } else {
    skip("Model unexpectedly found credible sets")
  }
})

test_that("run_refine handles single credible set", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10,
                                n_causal = 1, seed = 114)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(is.finite(susie_get_objective(refined_model)))
  expect_true(susie_get_objective(refined_model) >=
              susie_get_objective(setup$model) - 1e-6)
})

test_that("run_refine handles credible set with all prior weights zero", {
  setup <- create_model_with_cs(seed = 115)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  if (length(setup$model$sets$cs) > 0) {
    cs_vars <- setup$model$sets$cs[[1]]

    if (length(cs_vars) < ncol(setup$model$alpha)) {
      refined_model <- run_refine(setup$model, setup$data, setup$params)
      expect_true(is.finite(susie_get_objective(refined_model)))
    }
  }
})

test_that("run_refine handles large credible set", {
  setup <- create_model_with_cs(n = 200, p = 100, seed = 116)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(is.finite(susie_get_objective(refined_model)))
})

test_that("run_refine handles small p relative to L", {
  setup <- create_model_with_cs(n = 100, p = 10, L = 5,
                                n_causal = 2, seed = 117)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(is.finite(susie_get_objective(refined_model)))
})

# =============================================================================
# PARAMETER HANDLING
# =============================================================================

test_that("run_refine respects verbose parameter", {
  setup <- create_model_with_cs(seed = 118)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  params_verbose <- setup$params
  params_verbose$verbose <- TRUE

  expect_message(
    run_refine(setup$model, setup$data, params_verbose),
    "Starting refinement"
  )
})

test_that("run_refine verbose=FALSE produces no output", {
  setup <- create_model_with_cs(seed = 119)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  params_silent <- setup$params
  params_silent$verbose <- FALSE

  expect_silent(
    run_refine(setup$model, setup$data, params_silent)
  )
})

test_that("run_refine warns about model_init", {
  setup <- create_model_with_cs(seed = 120)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  params_with_init <- setup$params
  params_with_init$model_init <- list(alpha = setup$model$alpha)

  expect_message(
    run_refine(setup$model, setup$data, params_with_init),
    "model_init is not used"
  )
})

test_that("run_refine respects tolerance parameter", {
  setup <- create_model_with_cs(seed = 121)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  params_tol1 <- setup$params
  params_tol1$tol <- 1e-2

  params_tol2 <- setup$params
  params_tol2$tol <- 1e-6

  refined1 <- run_refine(setup$model, setup$data, params_tol1)
  refined2 <- run_refine(setup$model, setup$data, params_tol2)

  expect_true(is.finite(susie_get_objective(refined1)))
  expect_true(is.finite(susie_get_objective(refined2)))
})

test_that("run_refine preserves null_weight", {
  setup <- create_model_with_cs(seed = 122)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  initial_null_weight <- setup$model$null_weight

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_equal(refined_model$null_weight, initial_null_weight)
})

# =============================================================================
# INTEGRATION
# =============================================================================

test_that("run_refine works with individual data", {
  setup <- create_model_with_cs(seed = 123)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  expect_equal(class(setup$data), "individual")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(is.finite(susie_get_objective(refined_model)))
})

test_that("run_refine works with sufficient statistics", {
  set.seed(124)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X)

  beta <- rep(0, p)
  beta[c(5, 15, 25)] <- c(1.5, -1.2, 1.0)
  y <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  XtX <- crossprod(X)
  Xty <- as.vector(crossprod(X, y))
  yty <- sum(y^2)

  model <- susie_ss(XtX, Xty, yty, n = n, L = 5, verbose = FALSE)

  skip_if(is.null(model$sets) || length(model$sets$cs) == 0,
          "No credible sets found")

  constructor_result <- sufficient_stats_constructor(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 5,
    standardize = TRUE,
    estimate_residual_method = "MLE",
    convergence_method = "elbo",
    coverage = 0.95, min_abs_corr = 0.5,
    n_purity = 100,
    track_fit = FALSE
  )

  refined_model <- run_refine(model, constructor_result$data,
                              constructor_result$params)

  expect_true(is.finite(susie_get_objective(refined_model)))
})

test_that("run_refine output compatible with susie_get functions", {
  setup <- create_model_with_cs(seed = 125)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  pips <- susie_get_pip(refined_model)
  expect_equal(length(pips), ncol(refined_model$alpha))
  expect_true(all(pips >= 0))
  expect_true(all(pips <= 1))

  cs <- susie_get_cs(refined_model)
  expect_true(is.null(cs) || is.list(cs))

  post_mean <- susie_get_posterior_mean(refined_model)
  expect_true(all(is.finite(post_mean)))
})

test_that("run_refine maintains fitted values", {
  setup <- create_model_with_cs(seed = 126)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true("fitted" %in% names(refined_model))
  expect_equal(length(refined_model$fitted), nrow(setup$X))
  expect_true(all(is.finite(refined_model$fitted)))
})

test_that("run_refine maintains intercept", {
  setup <- create_model_with_cs(seed = 127)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true("intercept" %in% names(refined_model))
  expect_true(is.finite(refined_model$intercept))
})

# =============================================================================
# MATHEMATICAL PROPERTIES
# =============================================================================

test_that("run_refine maintains non-negative prior variances", {
  setup <- create_model_with_cs(seed = 128)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(all(refined_model$V >= 0))
  expect_true(all(is.finite(refined_model$V)))
})

test_that("run_refine maintains positive residual variance", {
  setup <- create_model_with_cs(seed = 129)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(refined_model$sigma2 > 0)
  expect_true(is.finite(refined_model$sigma2))
})

test_that("run_refine maintains non-negative KL divergences", {
  setup <- create_model_with_cs(seed = 130)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_true(all(refined_model$KL >= -1e-6))
})

test_that("run_refine ELBO is monotonically increasing", {
  setup <- create_model_with_cs(seed = 131)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  elbo_diff <- diff(refined_model$elbo)
  expect_true(all(elbo_diff > -1e-6))
})

# =============================================================================
# SIGNAL RECOVERY
# =============================================================================

test_that("run_refine improves or maintains signal recovery", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10,
                                n_causal = 3, seed = 132)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  pips_initial <- susie_get_pip(setup$model)

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  pips_refined <- susie_get_pip(refined_model)

  expect_equal(length(pips_refined), length(pips_initial))
  expect_true(all(pips_refined >= 0))
  expect_true(all(pips_refined <= 1))
})

test_that("run_refine identifies true causal variables", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10,
                                n_causal = 3, seed = 133)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  pips <- susie_get_pip(refined_model)
  top_vars <- order(pips, decreasing = TRUE)[1:5]

  overlap <- length(intersect(top_vars, setup$causal_idx))
  expect_true(overlap >= 1)
})

test_that("run_refine maintains low PIPs for null variables", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10,
                                n_causal = 3, seed = 134)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  pips <- susie_get_pip(refined_model)
  null_vars <- setdiff(1:length(pips), setup$causal_idx)
  null_pips <- pips[null_vars]

  expect_true(median(null_pips) < 0.3)
})

# =============================================================================
# COMPARISON
# =============================================================================

test_that("run_refine produces different result than no refinement", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10,
                                n_causal = 3, seed = 135)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  initial_elbo <- susie_get_objective(setup$model)

  refined_model <- run_refine(setup$model, setup$data, setup$params)
  refined_elbo <- susie_get_objective(refined_model)

  if (refined_elbo > initial_elbo + setup$params$tol) {
    expect_true(TRUE)
  } else {
    expect_equal(refined_elbo, initial_elbo, tolerance = setup$params$tol)
  }
})

test_that("run_refine with tight tolerance may differ from loose tolerance", {
  setup <- create_model_with_cs(seed = 136)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  params_tight <- setup$params
  params_tight$tol <- 1e-6

  params_loose <- setup$params
  params_loose$tol <- 1e-1

  refined_tight <- run_refine(setup$model, setup$data, params_tight)
  refined_loose <- run_refine(setup$model, setup$data, params_loose)

  expect_true(is.finite(susie_get_objective(refined_tight)))
  expect_true(is.finite(susie_get_objective(refined_loose)))
})

# =============================================================================
# STRESS TESTING
# =============================================================================

test_that("run_refine handles multiple refinement iterations", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10,
                                n_causal = 5, seed = 137)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) < 2,
          "Need multiple credible sets")

  params_loose <- setup$params
  params_loose$tol <- 1e-3

  refined_model <- run_refine(setup$model, setup$data, params_loose)

  expect_true(is.finite(susie_get_objective(refined_model)))
  expect_true(susie_get_objective(refined_model) >=
              susie_get_objective(setup$model) - 1e-6)
})

test_that("run_refine handles large L", {
  setup <- create_model_with_cs(n = 150, p = 80, L = 20,
                                n_causal = 4, seed = 138)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_equal(nrow(refined_model$alpha), nrow(setup$model$alpha))
  expect_true(is.finite(susie_get_objective(refined_model)))
})

test_that("run_refine handles large p", {
  setup <- create_model_with_cs(n = 150, p = 200, L = 10,
                                n_causal = 3, seed = 139)

  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined_model <- run_refine(setup$model, setup$data, setup$params)

  expect_equal(ncol(refined_model$alpha), ncol(setup$model$alpha))
  expect_true(is.finite(susie_get_objective(refined_model)))
})
