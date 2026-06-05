context("model methods")

# ---- Helper: workhorse model state ----
# Build a post-ibss_fit model with runtime convergence state.

build_workhorse_model <- function(n = 100, p = 50, L = 5, seed = 42) {
  setup <- setup_individual_data(n = n, p = p, L = L, seed = seed)
  model <- ibss_initialize(setup$data, setup$params)
  model$runtime <- list(prev_elbo = -Inf, prev_alpha = model$alpha)
  model <- ibss_fit(setup$data, setup$params, model)
  list(data = setup$data, params = setup$params, model = model)
}

# ---- Posterior Mean Accessors ----

test_that("get_posterior_mean_l returns alpha * mu for every effect", {
  wh    <- build_workhorse_model(seed = 1)
  model <- wh$model

  for (l in seq_len(nrow(model$alpha))) {
    pm <- get_posterior_mean_l(model, l)
    expect_equal(pm, model$alpha[l, ] * model$mu[l, ], tolerance = 1e-12)
  }
})

test_that("get_posterior_mean_l.default matches the generic", {
  wh    <- build_workhorse_model(seed = 2)
  model <- wh$model

  expect_equal(get_posterior_mean_l.default(model, 1),
               get_posterior_mean_l(model, 1),
               tolerance = 1e-12)
})

test_that("get_posterior_mean_sum equals colSums(alpha * mu)", {
  wh    <- build_workhorse_model(seed = 3)
  model <- wh$model

  expect_equal(get_posterior_mean_sum(model),
               colSums(model$alpha * model$mu),
               tolerance = 1e-12)
})

test_that("get_posterior_mean_sum.default matches the generic", {
  wh    <- build_workhorse_model(seed = 4)
  model <- wh$model

  expect_equal(get_posterior_mean_sum.default(model),
               get_posterior_mean_sum(model),
               tolerance = 1e-12)
})

test_that("get_posterior_mean_sum equals row-sum of per-effect posterior means", {
  wh    <- build_workhorse_model(seed = 5)
  model <- wh$model
  L     <- nrow(model$alpha)

  per_effect <- sapply(seq_len(L), function(l) get_posterior_mean_l(model, l))
  expect_equal(get_posterior_mean_sum(model), rowSums(per_effect), tolerance = 1e-12)
})

# ---- get_objective ----

test_that("get_objective returns a finite scalar ELBO", {
  wh  <- build_workhorse_model(seed = 6)
  obj <- get_objective(wh$data, wh$params, wh$model)

  expect_length(obj, 1)
  expect_true(is.finite(obj))
})

test_that("get_objective.default matches the generic", {
  wh <- build_workhorse_model(seed = 7)

  expect_equal(get_objective.default(wh$data, wh$params, wh$model),
               get_objective(wh$data, wh$params, wh$model),
               tolerance = 1e-10)
})

test_that("get_objective errors when KL contains Inf", {
  wh          <- build_workhorse_model(n = 60, p = 30, L = 3, seed = 8)
  model       <- wh$model
  model$KL[1] <- Inf

  expect_error(get_objective(wh$data, wh$params, model), "infinite ELBO value")
})

test_that("get_objective NIG model returns finite ELBO distinct from default ELBO", {
  set.seed(9)
  n <- 80; p <- 20
  X    <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[3] <- 2
  y    <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  fit_default <- suppressWarnings(susie(X, y, L = 2,
                                        estimate_residual_method = "MLE",
                                        max_iter = 5, verbose = FALSE))
  fit_nig     <- suppressWarnings(susie(X, y, L = 2,
                                        estimate_residual_method = "NIG",
                                        alpha0 = 2, beta0 = 1,
                                        max_iter = 5, verbose = FALSE))

  expect_true(is.finite(tail(fit_default$elbo, 1)))
  expect_true(is.finite(tail(fit_nig$elbo, 1)))
  # NIG and default use different objectives; ELBO values need not be equal.
  expect_true(!isTRUE(all.equal(tail(fit_default$elbo, 1),
                                tail(fit_nig$elbo, 1),
                                tolerance = 1e-6)))
})

# ---- format_extra_diag ----

test_that("format_extra_diag returns empty string when lambda_bias is NULL", {
  model <- structure(list(), class = "susie")
  expect_equal(format_extra_diag(model), "")
})

test_that("format_extra_diag formats a scalar lambda_bias", {
  model <- structure(list(lambda_bias = 1.5), class = "susie")
  expect_true(grepl("lambda_infl=", format_extra_diag(model)))
})

test_that("format_extra_diag appends B_eff when B_corrected is constant", {
  model <- structure(list(lambda_bias = 1.5, B_corrected = rep(2, 5)), class = "susie")
  out   <- format_extra_diag(model)
  expect_true(grepl("lambda_infl=", out))
  expect_true(grepl("B_eff=", out))
})

test_that("format_extra_diag sanitizes non-finite lambda_bias to zero", {
  model <- structure(list(lambda_bias = Inf), class = "susie")
  out   <- format_extra_diag(model)
  expect_true(grepl("lambda_infl=", out))
  expect_true(grepl("0", out))
})

test_that("format_extra_diag errors when lambda_bias is not a scalar", {
  model <- structure(list(lambda_bias = c(1, 2)), class = "susie")
  expect_error(format_extra_diag(model), "lambda_bias must be a scalar")
})

# ---- check_convergence: first iteration ----

test_that("check_convergence does not converge on iter=1", {
  wh     <- build_workhorse_model(seed = 10)
  params <- wh$params
  params$convergence_method <- "elbo"

  model <- check_convergence(wh$data, params, wh$model, c(-Inf, -90, -89, -88.5), iter = 1)

  expect_false(model$converged)
})

test_that("check_convergence prints tabular header when verbose on iter=1", {
  wh     <- build_workhorse_model(seed = 11)
  params <- wh$params
  params$convergence_method <- "elbo"
  params$verbose            <- TRUE

  expect_message(
    check_convergence(wh$data, params, wh$model, c(-Inf, -90, -89, -88.5), iter = 1),
    "iter"
  )
})

# ---- check_convergence: ELBO path ----

test_that("check_convergence (elbo) converges when ELBO change < tol", {
  wh     <- build_workhorse_model(seed = 12)
  params <- wh$params
  params$convergence_method <- "elbo"
  params$tol                <- 1e-3

  model                    <- wh$model
  model$runtime$prev_elbo  <- -88.6

  model <- check_convergence(wh$data, params, model, c(-Inf, -90, -89, -88.5999), iter = 3)
  expect_true(model$converged)
})

test_that("check_convergence (elbo) does not converge when ELBO change > tol", {
  wh     <- build_workhorse_model(seed = 13)
  params <- wh$params
  params$convergence_method <- "elbo"
  params$tol                <- 1e-3

  model                   <- wh$model
  model$runtime$prev_elbo <- -90

  model <- check_convergence(wh$data, params, model, c(-Inf, -92, -91, -88), iter = 3)
  expect_false(model$converged)
})

test_that("check_convergence (elbo) warns when the ELBO decreases", {
  wh     <- build_workhorse_model(seed = 14)
  params <- wh$params
  params$convergence_method <- "elbo"
  params$tol                <- 1e-3

  model                   <- wh$model
  model$runtime$prev_elbo <- -88

  expect_message(
    model <- check_convergence(wh$data, params, model, c(-Inf, -90, -89, -90), iter = 3),
    "ELBO decreased"
  )
  expect_false(model$converged)
})

test_that("check_convergence (elbo) prints converged message when verbose", {
  wh     <- build_workhorse_model(seed = 15)
  params <- wh$params
  params$convergence_method <- "elbo"
  params$tol                <- 1e-3
  params$verbose            <- TRUE

  model                   <- wh$model
  model$runtime$prev_elbo <- -88.6

  expect_message(
    check_convergence(wh$data, params, model, c(-Inf, -90, -89, -88.5999), iter = 3),
    "converged"
  )
})

test_that("check_convergence (elbo) lambda_bias guard blocks convergence", {
  wh     <- build_workhorse_model(seed = 16)
  params <- wh$params
  params$convergence_method <- "elbo"
  params$tol                <- 1e-3

  model                            <- wh$model
  model$runtime$prev_elbo          <- -88.6
  model$runtime$lambda_bias_diff   <- 0.5

  model <- check_convergence(wh$data, params, model, c(-Inf, -90, -89, -88.5999), iter = 3)
  expect_false(model$converged)
})

# ---- check_convergence: ELBO failure fallback ----

test_that("check_convergence warns and falls back to PIP when ELBO is NA", {
  wh     <- build_workhorse_model(seed = 17)
  params <- wh$params
  params$convergence_method <- "elbo"

  model                   <- wh$model
  model$runtime$prev_elbo <- -88

  expect_message(
    model <- check_convergence(wh$data, params, model, c(-Inf, -90, -89, NA), iter = 3),
    "NA/infinite ELBO"
  )
  expect_false(is.null(model$converged))
})

# ---- check_convergence: PIP path ----

test_that("check_convergence (pip) sets converged flag and pip_diff", {
  wh     <- build_workhorse_model(seed = 18)
  params <- wh$params
  params$convergence_method <- "pip"
  params$tol                <- 1e-3

  model                   <- wh$model
  model$runtime$prev_elbo <- -88

  model <- check_convergence(wh$data, params, model, c(-Inf, -90, -89, -88.5), iter = 3)

  expect_false(is.null(model$converged))
  expect_true(is.numeric(model$runtime$pip_diff))
})

test_that("check_convergence (pip) converges when alpha is at a fixed point", {
  wh     <- build_workhorse_model(seed = 19)
  params <- wh$params
  params$convergence_method <- "pip"
  params$tol                <- 1e-3

  model                      <- wh$model
  model$runtime$prev_alpha   <- model$alpha
  model$runtime$prev_elbo    <- -88

  model <- check_convergence(wh$data, params, model, c(-Inf, -90, -89, -88.5), iter = 3)
  expect_true(model$converged)
})

test_that("check_convergence (pip) prints diagnostics when verbose", {
  wh     <- build_workhorse_model(seed = 20)
  params <- wh$params
  params$convergence_method <- "pip"
  params$tol                <- 1e-3
  params$verbose            <- TRUE

  model                    <- wh$model
  model$runtime$prev_alpha <- model$alpha
  model$runtime$prev_elbo  <- -88

  expect_message(
    check_convergence(wh$data, params, model, c(-Inf, -90, -89, -88.5), iter = 3),
    "max\\|d\\(alpha,PIP\\)\\|"
  )
})

test_that("check_convergence (pip) lambda_bias guard blocks convergence and sets reason", {
  wh     <- build_workhorse_model(seed = 21)
  params <- wh$params
  params$convergence_method <- "pip"
  params$tol                <- 1e-3

  model                          <- wh$model
  model$runtime$prev_alpha       <- model$alpha
  model$runtime$prev_elbo        <- -88
  model$runtime$lambda_bias_diff <- 0.5

  model <- check_convergence(wh$data, params, model, c(-Inf, -90, -89, -88.5), iter = 3)

  expect_false(model$converged)
  expect_true(grepl("lambda_infl_changed", model$convergence_reason))
})

# ---- check_convergence: verbose NA-ELBO on iter=1 ----

test_that("check_convergence.default emits message when ELBO is NA on iter=1 with verbose", {
  set.seed(26)
  setup                             <- setup_individual_data(n = 80, p = 40)
  setup$params$verbose              <- TRUE
  setup$params$convergence_method   <- "pip"

  var_y        <- get_var_y.individual(setup$data)
  model        <- initialize_susie_model.individual(setup$data, setup$params, var_y)
  model$sigma2 <- var_y

  expect_message(
    check_convergence.default(setup$data, setup$params, model,
                              c(-Inf, NA_real_), iter = 1),
    regexp = "iter"
  )
})

# ---- check_convergence: ash final-pass on convergence ----

test_that("check_convergence.default triggers run_final_ash_pass on ELBO convergence with ash", {
  set.seed(1010)
  setup                          <- setup_individual_data(n = 80, p = 30, L = 3)
  setup$params$unmappable_effects <- "ash"
  setup$params$convergence_method <- "elbo"
  setup$params$tol                <- 0.1

  options(susie.skip_mrash = TRUE)
  on.exit(options(susie.skip_mrash = NULL), add = TRUE)

  model         <- ibss_initialize(setup$data, setup$params)
  model         <- init_ash_fields(model, setup$data$n, setup$data$p, setup$params$L,
                                   is_individual = TRUE)
  model$runtime <- list(prev_elbo = -100 + 1e-5, prev_alpha = model$alpha,
                        lambda_bias_diff = 0)

  result <- check_convergence.default(setup$data, setup$params, model,
                                      c(-Inf, -100, -100 + 1e-5), iter = 2)

  expect_true(result$converged)
  expect_s3_class(result, "susie")
})

# ---- check_convergence: end-to-end via susie() ----

test_that("susie converges with convergence_method='elbo'", {
  set.seed(22)
  n <- 120; p <- 60
  X    <- matrix(rnorm(n * p), n, p)
  b    <- rep(0, p); b[c(1, 10, 20)] <- c(1.5, -1.2, 1.0)
  y    <- as.vector(X %*% b + rnorm(n))

  fit <- susie(X, y, L = 6, convergence_method = "elbo", verbose = FALSE)

  expect_true(fit$converged)
  expect_true(all(is.finite(fit$elbo)))
})

test_that("susie converges with convergence_method='pip'", {
  set.seed(23)
  n <- 120; p <- 60
  X    <- matrix(rnorm(n * p), n, p)
  b    <- rep(0, p); b[c(1, 10, 20)] <- c(1.5, -1.2, 1.0)
  y    <- as.vector(X %*% b + rnorm(n))

  fit <- susie(X, y, L = 6, convergence_method = "pip", verbose = FALSE)

  expect_true(fit$converged)
})

test_that("susie verbose elbo path emits per-iteration ELBO output", {
  set.seed(24)
  n <- 100; p <- 50
  X    <- matrix(rnorm(n * p), n, p)
  b    <- rep(0, p); b[c(1, 10)] <- c(1.5, -1.2)
  y    <- as.vector(X %*% b + rnorm(n))

  expect_message(
    susie(X, y, L = 5, convergence_method = "elbo", verbose = TRUE),
    "ELBO"
  )
})

test_that("susie with max_iter=1 does not converge and warns", {
  set.seed(25)
  n <- 100; p <- 50
  X    <- matrix(rnorm(n * p), n, p)
  b    <- rep(0, p); b[c(1, 10, 20)] <- c(1.5, -1.2, 1.0)
  y    <- as.vector(X %*% b + rnorm(n))

  expect_warning(
    fit <- susie(X, y, L = 6, max_iter = 1, convergence_method = "elbo"),
    "did not converge"
  )
  expect_false(fit$converged)
})

# ---- PIP vs ELBO convergence comparison ----

test_that("pip and elbo convergence methods both yield finite ELBO traces", {
  set.seed(30)
  n <- 100; p <- 40
  X    <- matrix(rnorm(n * p), n, p)
  b    <- rep(0, p); b[c(2, 15)] <- c(1.5, -1.5)
  y    <- as.vector(X %*% b + rnorm(n))

  fit_elbo <- suppressWarnings(susie(X, y, L = 4,
                                     convergence_method = "elbo", verbose = FALSE))
  fit_pip  <- suppressWarnings(susie(X, y, L = 4,
                                     convergence_method = "pip",  verbose = FALSE))

  expect_true(all(is.finite(fit_elbo$elbo)))
  expect_true(all(is.finite(fit_pip$elbo)))
})
