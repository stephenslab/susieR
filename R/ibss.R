# This file contains the three core functions for SuSiE (and its modifications):
# 1) ibss_initialize(), 2) ibss_fit(), and 3) ibss_finalize(). All of these
# functions make up the main workhorse "susie_engine()" function.


### Initialization ###

# This function initializes the SuSiE model object

ibss_initialize <- function(data,
                            L                     = min(10, data$p),
                            scaled_prior_variance = 0.2,
                            residual_variance     = NULL,
                            prior_weights         = NULL,
                            null_weight           = NULL,
                            model_init            = NULL){
  # Define p
  p <- data$p

  # Check prior variance
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")

  # Check prior weights
  if (is.null(prior_weights))
    prior_weights <- rep(1 / p, p)

  # Check number of single effects
  if (p < L)
    L = p

  # Set variance of y
  var_y <- if (!is.null(data$y)) var(drop(data$y)) else data$yty / (data$n - 1)

  # Check residual variance
  if (is.null(residual_variance)) {
    residual_variance <- var_y
  }

  # Set alpha, mu, mu^2, V
  alpha  <- matrix(1 / p, L, p)
  mu     <- matrix(0,     L, p)
  mu2    <- matrix(0,     L, p)
  V      <- rep(scaled_prior_variance * var_y, L)

  # Check for pre-initialized model
  if (!missing(model_init) && !is.null(model_init)) {
    if (!inherits(model_init, "susie"))
      stop("model_init must be a 'susie' object")

    if (nrow(model_init$alpha) > L)
      stop("model_init has more effects than requested L")

    if (max(model_init$alpha) > 1 || min(model_init$alpha) < 0)
      stop("model_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")

    # TODO: add more checks

    model_init <- prune_zero_variance_effects(model_init)  # TODO: add generic here to prune effects with zero prior var

    idx <- seq_len(nrow(model_init$alpha))
    alpha[idx, ] <- model_init$alpha
    mu[idx, ]    <- model_init$mu
    mu2[idx, ]   <- model_init$mu2
    V[idx]       <- model_init$V[idx]

    if (!is.numeric(V))
      stop("Input prior variance must be numeric")
    if (!all(V >= 0))
      stop("prior variance must be non-negative")
    if (!all(dim(mu) == dim(mu2)))
      stop("dimension of mu and mu2 in input object do not match")
    if (!all(dim(mu) == dim(alpha)))
      stop("dimension of mu and alpha in input object do not match")
    if (nrow(alpha) != length(V))
      stop("Input prior variance V must have length of nrow of alpha in ",
           "input object")
  }

  fitted <- initialize_fitted(data, alpha, mu)


  # Return Initialized Model
  model <- c(
    list(alpha        = alpha,
         mu           = mu,
         mu2          = mu2,
         KL           = rep(as.numeric(NA), L),
         lbf          = rep(as.numeric(NA), L),
         lbf_variable = matrix(as.numeric(NA), L, p),
         sigma2       = residual_variance,
         V            = V,
         pi           = prior_weights),
         fitted)

  class(model) <- "susie"

  # Set null index (for refine step)
  model$null_index <- if (is.null(null_weight)) 0 else length(prior_weights)

  return(model)
}


### Fitting ###

# This function updates each of the L effects

ibss_fit = function(data, model,
                    estimate_prior_variance = TRUE,
                    estimate_prior_method = c("optim"),
                    check_null_threshold = 0){
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Repeat for each effect to update
  L <- nrow(model$alpha)
  if (L > 0)
    for (l in seq_len(L)){
      model <- single_effect_update(data, model, l,
                                  optimize_V = estimate_prior_method,
                                  check_null_threshold = check_null_threshold)
    }

  return(model)
}

### Finalize ###

# This function takes the final SuSiE model object and appends credible sets,
# posterior inclusion probabilities, and runs other post-processing steps.

ibss_finalize <- function(data,
                          model,
                          coverage               = 0.95,
                          min_abs_corr           = 0.50,
                          median_abs_corr        = NULL,
                          prior_tol              = 1e-9,
                          n_purity               = 100,
                          compute_univariate_zscore = FALSE,
                          intercept              = TRUE,
                          standardize            = TRUE,
                          elbo                   = NULL,
                          iter                   = NA_integer_,
                          null_weight            = NULL,
                          track_fit              = FALSE,
                          tracking               = NULL) {

  # Append ELBO & iteration count to model output
  model$elbo  <- elbo[2:(iter + 1)] # consider changing this to (susie_ss.R L262 format)
                                    # this would instead remove infinite values + NA
  model$niter <- iter

  # Intercept & Fitted Values
  model$X_column_scale_factors <- get_scale_factors(data)
  model$intercept              <- get_intercept(data, model, intercept)
  model$fitted                 <- get_fitted(data, model, intercept)

  # Tracking Across Iterations
  if (track_fit)
    model$trace = tracking

  # Credible Sets
  model$sets <- get_cs(data, model, coverage, min_abs_corr, n_purity) # TODO: Check on median_abs_corr parameter

  # Posterior Inclusion Probabilities
  model$pip <- get_pip(data, model, coverage, min_abs_corr, prior_tol)

  # Assign Variable Names
  model <- get_variable_names(data, model, null_weight)

  # Compute z-scores
  model$z <- get_zscore(data, model, compute_univariate_zscore,
                        intercept, standardize, null_weight)

  return(model)
}

