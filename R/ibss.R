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
  # Define p and var_y
  p <- data$p
  var_y <- get_var_y(data)

  # Check prior variance
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")

  # Check prior weights
  if (is.null(prior_weights))
    prior_weights <- rep(1 / p, p)

  # Check number of single effects
  if (p < L)
    L = p

  # Check residual variance
  if (is.null(residual_variance)) {
    residual_variance <- var_y
  }

  # Check for pre-initialized model
  if (!missing(model_init) && !is.null(model_init)) {

    # Validate the contents of model_init
    validate_init(model_init, L, null_weight)

    # Prune effects with zero prior variance
    model_init <- susie_prune_single_effects(model_init)

    # Adjust the number of effects
    adjustment <- adjust_L(model_init, L, V = rep(scaled_prior_variance * var_y, L))

    # Return adjusted initialized model and number of effects
    model_init <- adjustment$model_init
    L          <- adjustment$L
  }

  # Build blank matrices
  # alpha <- matrix(1 / p, L, p)
  # mu    <- matrix(0,     L, p)
  # mu2   <- matrix(0,     L, p)
  # V     <- rep(scaled_prior_variance * var_y, L)

  mat_init <- initialize_matrices(data, L,
                                  scaled_prior_variance, var_y,
                                  residual_variance, prior_weights)

  # Overwrite blank rows if we used an initialized model
  if (!missing(model_init) && !is.null(model_init)) {
    idx                        <- seq_len(nrow(model_init$alpha))
    mat_init$alpha[idx, ]      <- model_init$alpha
    mat_init$mu[idx, ]         <- model_init$mu
    mat_init$mu2[idx, ]        <- model_init$mu2
    mat_init$V[idx]            <- model_init$V[idx]
    mat_init$residual_variance <- model_init$sigma2
    mat_init$prior_weights     <- model_init$pi
  }

  # Initialize fitted values
  fitted <- initialize_fitted(data, mat_init$alpha, mat_init$mu)

  # Set null index (for refine step)
  null_index <- initialize_null_index(null_weight, p)

  # TODO: Make generic function for nonsparse methods (inf/ash)

  # Return assembled SuSiE object (this can be cleaned up later -- works for now)
  model <- c(
    mat_init,
    list(null_index = null_index),
    fitted
  )

  class(model) <- "susie"

  return(model)
}


### Fitting ###

# This function updates each of the L effects

ibss_fit = function(data, model,
                    estimate_prior_variance = TRUE,
                    estimate_prior_method   = c("optim","EM","simple"),
                    check_null_threshold    = 0,
                    check_prior             = FALSE){
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Repeat for each effect to update
  L <- nrow(model$alpha)
  if (L > 0)
    for (l in seq_len(L)){
      model <- single_effect_update(data, model, l,
                                  optimize_V = estimate_prior_method,
                                  check_null_threshold = check_null_threshold)
    }

  # Validate prior variance is reasonable
  validate_prior(data, model, check_prior)

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
  #model$elbo  <- elbo[2:(iter + 1)] # consider changing this to (susie_ss.R L262 format)
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

