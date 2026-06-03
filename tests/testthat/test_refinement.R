context("Refinement unit tests")

# ---- Basic properties ----

test_that("run_refine returns a valid susie model with correct structure and probability constraints", {
  setup <- create_model_with_cs(seed = 100)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  refined <- run_refine(setup$model, setup$data, setup$params)

  # Required fields present
  expect_true("alpha" %in% names(refined))
  expect_true("mu"    %in% names(refined))
  expect_true("V"     %in% names(refined))
  expect_true("sigma2" %in% names(refined))
  expect_true("elbo"  %in% names(refined))

  # Dimensions preserved
  expect_equal(nrow(refined$alpha), nrow(setup$model$alpha))
  expect_equal(ncol(refined$alpha), ncol(setup$model$alpha))
  expect_equal(nrow(refined$mu),    nrow(setup$model$mu))
  expect_equal(ncol(refined$mu),    ncol(setup$model$mu))

  # Valid probability distributions
  expect_true(all(refined$alpha >= 0))
  expect_true(all(refined$alpha <= 1))
  row_sums <- rowSums(refined$alpha)
  expect_true(all(abs(row_sums - 1) < 1e-10))

  # Finite ELBO
  expect_true(all(is.finite(refined$elbo)))
  expect_true(is.finite(susie_get_objective(refined)))
})

test_that("run_refine maintains or improves ELBO", {
  setup <- create_model_with_cs(seed = 101)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  initial_elbo <- susie_get_objective(setup$model)
  refined <- run_refine(setup$model, setup$data, setup$params)

  expect_gte(susie_get_objective(refined), initial_elbo - 1e-6)
})

# ---- Convergence behavior ----

test_that("run_refine convergence is deterministic", {
  setup <- create_model_with_cs(seed = 112)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found in initial model")

  refined1 <- run_refine(setup$model, setup$data, setup$params)
  refined2 <- run_refine(setup$model, setup$data, setup$params)

  expect_equal(susie_get_objective(refined1), susie_get_objective(refined2),
               tolerance = 1e-8)
})

test_that("run_refine ELBO is monotonically non-decreasing within a run", {
  setup <- create_model_with_cs(seed = 131)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined <- run_refine(setup$model, setup$data, setup$params)

  elbo_diff <- diff(refined$elbo)
  expect_true(all(elbo_diff > -1e-6))
})

# ---- Edge cases ----

test_that("run_refine with no credible sets returns model unchanged", {
  set.seed(113)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  model <- suppressWarnings(susie(X, y, L = 5, verbose = FALSE))

  constructor_result <- individual_data_constructor(
    X = X, y = y, L = 5,
    standardize = TRUE, intercept = TRUE,
    estimate_residual_method = "MLE",
    convergence_method = "elbo",
    coverage = 0.95, min_abs_corr = 0.5,
    n_purity = 100,
    track_fit = FALSE
  )

  # If susie happens to find CS with this pure-noise data, re-seed to guarantee no CS
  if (!is.null(model$sets) && length(model$sets$cs) > 0) {
    set.seed(9999)
    X2 <- matrix(rnorm(n * p), n, p)
    y2 <- rnorm(n)
    model <- suppressWarnings(susie(X2, y2, L = 5, verbose = FALSE))
    constructor_result <- individual_data_constructor(
      X = X2, y = y2, L = 5,
      standardize = TRUE, intercept = TRUE,
      estimate_residual_method = "MLE",
      convergence_method = "elbo",
      coverage = 0.95, min_abs_corr = 0.5,
      n_purity = 100,
      track_fit = FALSE
    )
  }

  skip_if(!is.null(model$sets) && length(model$sets$cs) > 0,
          "Could not produce a no-CS model")

  initial_elbo <- susie_get_objective(model)
  refined <- run_refine(model, constructor_result$data, constructor_result$params)

  expect_equal(susie_get_objective(refined), initial_elbo, tolerance = 1e-8)
})

test_that("run_refine with single credible set maintains or improves ELBO", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10, n_causal = 1, seed = 114)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  initial_elbo <- susie_get_objective(setup$model)
  refined <- run_refine(setup$model, setup$data, setup$params)

  expect_true(is.finite(susie_get_objective(refined)))
  expect_gte(susie_get_objective(refined), initial_elbo - 1e-6)
})

# ---- Parameter handling ----

test_that("run_refine with verbose=TRUE prints Block ascent progress messages", {
  setup <- create_model_with_cs(seed = 118)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  params_verbose <- setup$params
  params_verbose$verbose <- TRUE

  expect_message(
    run_refine(setup$model, setup$data, params_verbose),
    "Block ascent iter"
  )
})

test_that("run_refine with verbose=FALSE produces no output", {
  setup <- create_model_with_cs(seed = 119)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  params_silent <- setup$params
  params_silent$verbose <- FALSE

  expect_silent(run_refine(setup$model, setup$data, params_silent))
})

test_that("run_refine emits hint message when model_init is supplied", {
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

test_that("run_refine preserves null_weight unchanged", {
  setup <- create_model_with_cs(seed = 122)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  initial_null_weight <- setup$model$null_weight
  refined <- run_refine(setup$model, setup$data, setup$params)

  expect_equal(refined$null_weight, initial_null_weight)
})

# ---- Mathematical properties ----

test_that("run_refine maintains valid prior variances, residual variance, prior weights, and KL", {
  setup <- create_model_with_cs(seed = 128)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined <- run_refine(setup$model, setup$data, setup$params)

  # Prior variances non-negative and finite
  expect_true(all(refined$V >= 0))
  expect_true(all(is.finite(refined$V)))

  # Residual variance positive and finite
  expect_gt(refined$sigma2, 0)
  expect_true(is.finite(refined$sigma2))

  # Prior weights valid distribution
  expect_equal(length(refined$pi), length(setup$model$pi))
  expect_true(all(refined$pi >= 0))
  expect_lt(abs(sum(refined$pi) - 1), 1e-10)

  # KL non-negative (within numerical tolerance)
  expect_true(all(refined$KL >= -1e-6))
})

# ---- Integration ----

test_that("run_refine works with individual data and output is compatible with susie_get functions", {
  setup <- create_model_with_cs(seed = 123)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  expect_equal(class(setup$data), "individual")
  refined <- run_refine(setup$model, setup$data, setup$params)

  # PIPs valid
  pips <- susie_get_pip(refined)
  expect_equal(length(pips), ncol(refined$alpha))
  expect_true(all(pips >= 0 & pips <= 1))

  # CS is NULL or list
  cs <- susie_get_cs(refined)
  expect_true(is.null(cs) || is.list(cs))

  # Posterior mean finite
  expect_true(all(is.finite(susie_get_posterior_mean(refined))))

  # Fitted values present, correct length, finite
  expect_true("fitted" %in% names(refined))
  expect_equal(length(refined$fitted), nrow(setup$X))
  expect_true(all(is.finite(refined$fitted)))

  # Intercept present and finite
  expect_true("intercept" %in% names(refined))
  expect_true(is.finite(refined$intercept))
})

test_that("run_refine works with sufficient statistics data", {
  set.seed(124)
  n <- 100; p <- 50
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

  refined <- run_refine(model, constructor_result$data, constructor_result$params)

  expect_true(is.finite(susie_get_objective(refined)))
  expect_gte(susie_get_objective(refined), susie_get_objective(model) - 1e-6)
})

# ---- Signal recovery ----

test_that("run_refine maintains valid PIPs and recovers causal signal", {
  setup <- create_model_with_cs(n = 200, p = 100, L = 10, n_causal = 3, seed = 132)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined <- run_refine(setup$model, setup$data, setup$params)
  pips <- susie_get_pip(refined)

  expect_equal(length(pips), ncol(refined$alpha))
  expect_true(all(pips >= 0 & pips <= 1))

  # At least one of the top-5 PIPs overlaps with a true causal variable
  top_vars <- order(pips, decreasing = TRUE)[1:5]
  expect_gte(length(intersect(top_vars, setup$causal_idx)), 1L)

  # Null variables have low median PIP
  null_pips <- pips[setdiff(seq_along(pips), setup$causal_idx)]
  expect_lt(median(null_pips), 0.3)
})

# ---- Stress testing ----

test_that("run_refine handles varied dimensions: large L, large p, multiple signals", {
  configs <- list(
    list(n = 150, p = 80,  L = 20, n_causal = 4, seed = 138),
    list(n = 150, p = 200, L = 10, n_causal = 3, seed = 139),
    list(n = 200, p = 100, L = 10, n_causal = 5, seed = 137)
  )

  for (cfg in configs) {
    setup <- create_model_with_cs(n = cfg$n, p = cfg$p, L = cfg$L,
                                  n_causal = cfg$n_causal, seed = cfg$seed)
    skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
            "No credible sets found")

    params_loose <- setup$params
    params_loose$tol <- 1e-3

    refined <- run_refine(setup$model, setup$data, params_loose)

    expect_equal(nrow(refined$alpha), nrow(setup$model$alpha))
    expect_equal(ncol(refined$alpha), ncol(setup$model$alpha))
    expect_true(is.finite(susie_get_objective(refined)))
    expect_gte(susie_get_objective(refined), susie_get_objective(setup$model) - 1e-6)
  }
})

# ---- block_coordinate_ascent direct tests ----

test_that("block_coordinate_ascent returns converged=TRUE when step_fn signals converged", {
  setup <- create_model_with_cs(seed = 140)

  done_step <- function(model, data, iter)
    list(model = model, converged = TRUE)

  result <- block_coordinate_ascent(setup$model, setup$data, done_step,
                                    max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_true(result$converged)
})

test_that("block_coordinate_ascent reverts to previous model when ELBO decreases", {
  setup <- create_model_with_cs(seed = 141)
  orig_elbo <- susie_get_objective(setup$model)

  worse_step <- function(model, data, iter) {
    bad <- model
    bad$elbo <- c(bad$elbo, susie_get_objective(model) - 5)
    list(model = bad)
  }

  result <- suppressWarnings(
    block_coordinate_ascent(setup$model, setup$data, worse_step,
                            max_iter = 5, tol = 1e-3, verbose = FALSE)
  )

  expect_true(result$converged)
  expect_equal(susie_get_objective(result), orig_elbo, tolerance = 1e-8)
})

test_that("block_coordinate_ascent prints revert message when verbose and ELBO decreases", {
  setup <- create_model_with_cs(seed = 142)

  worse_step <- function(model, data, iter) {
    bad <- model
    bad$elbo <- c(bad$elbo, susie_get_objective(model) - 5)
    list(model = bad)
  }

  expect_message(
    suppressWarnings(
      block_coordinate_ascent(setup$model, setup$data, worse_step,
                              max_iter = 5, tol = 1e-3, verbose = TRUE)
    ),
    "did not improve ELBO"
  )
})

test_that("block_coordinate_ascent converges when ELBO change falls below tolerance", {
  setup <- create_model_with_cs(seed = 143)

  noop_step <- function(model, data, iter)
    list(model = model)

  result <- block_coordinate_ascent(setup$model, setup$data, noop_step,
                                    max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_true(result$converged)
})

test_that("block_coordinate_ascent prints log_msg in verbose mode", {
  setup <- create_model_with_cs(seed = 144)

  step_fn <- function(model, data, iter) {
    out <- model
    out$elbo <- c(out$elbo, susie_get_objective(model) + 1e-6)
    list(model = out, converged = TRUE, log_msg = "custom-log")
  }

  expect_message(
    block_coordinate_ascent(setup$model, setup$data, step_fn,
                            max_iter = 5, tol = 1e-3, verbose = TRUE),
    "custom-log"
  )
})

test_that("block_coordinate_ascent emits did-not-converge message when max_iter exhausted", {
  setup <- create_model_with_cs(seed = 145)

  improving_step <- function(model, data, iter) {
    out <- model
    out$elbo <- c(out$elbo, susie_get_objective(model) + 10)
    list(model = out, converged = FALSE)
  }

  expect_message(
    result <- block_coordinate_ascent(setup$model, setup$data, improving_step,
                                      max_iter = 3, tol = 1e-12, verbose = FALSE),
    "did not converge"
  )
  expect_false(result$converged)
})

test_that("block_coordinate_ascent accepts and forwards updated data from step_fn", {
  setup <- create_model_with_cs(seed = 146)

  step_fn <- function(model, data, iter) {
    out <- model
    out$elbo <- c(out$elbo, susie_get_objective(model))
    list(model = out, data = data)
  }

  result <- block_coordinate_ascent(setup$model, setup$data, step_fn,
                                    max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_true(result$converged)
})

# ---- Candidate generation and acceptance ----

test_that("run_refine accepts best candidate when multiple signals allow ELBO improvement", {
  setup <- create_model_with_cs(n = 300, p = 150, L = 12, n_causal = 5, seed = 147)
  skip_if(is.null(setup$model$sets) || length(setup$model$sets$cs) == 0,
          "No credible sets found")

  refined <- run_refine(setup$model, setup$data, setup$params)

  expect_gte(susie_get_objective(refined), susie_get_objective(setup$model) - 1e-6)
  expect_true(is.finite(susie_get_objective(refined)))
  expect_equal(dim(refined$alpha), dim(setup$model$alpha))
})

test_that("run_refine handles all-zero pw_cs (break) and empty candidates (short-circuit) branches", {
  set.seed(148)
  n <- 100; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(3, p)
  y <- as.vector(X %*% b + rnorm(n, sd = 0.1))

  pw <- rep(0, p)
  pw[1] <- 1

  fit <- suppressWarnings(suppressMessages(
    susie(X, y, L = 3, prior_weights = pw, refine = TRUE,
          max_iter = 10, verbose = FALSE)
  ))

  expect_s3_class(fit, "susie")
  expect_gt(fit$sigma2, 0)
})

test_that("run_refine accepts a strictly higher-ELBO candidate when the greedy fit is trapped", {
  # A decoy predictor correlated with the sum of the two causal variables traps
  # greedy IBSS in a worse local optimum; refinement escapes it, so refine = TRUE
  # ends at a STRICTLY higher ELBO than refine = FALSE (the accept-candidate path).
  set.seed(14)
  n  <- 200
  s1 <- rnorm(n); s2 <- rnorm(n)
  decoy <- 0.7 * s1 + 0.7 * s2 + rnorm(n) * 0.25
  X <- scale(cbind(s1, s2, decoy, matrix(rnorm(n * 2), n, 2)))
  y <- as.vector(1.6 * s1 + 1.6 * s2 + rnorm(n))

  fit_norefine <- suppressWarnings(susie(X, y, L = 3, refine = FALSE, verbose = FALSE))
  fit_refine   <- suppressWarnings(susie(X, y, L = 3, refine = TRUE,  verbose = FALSE))

  expect_s3_class(fit_refine, "susie")
  expect_gt(susie_get_objective(fit_refine),
            susie_get_objective(fit_norefine) + 1e-6)
})
