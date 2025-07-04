### SuSiE.ash data backend functions ###

# Initialize Fitted Values
initialize_fitted.individual <- function(data, alpha, mu){
  return(list(Xr = compute_Xb(data$X, colSums(alpha * mu))))
}

# Initialize Matrices
initialize_matrices.individual <- function(data, L, scaled_prior_variance, var_y,
                                           residual_variance, prior_weights, ...){
  p <- data$p

  mat_init <- list(
    alpha        = matrix(1 / p, L, p),
    mu           = matrix(0,     L, p),
    mu2          = matrix(0,     L, p),
    V            = rep(scaled_prior_variance * var_y, L),
    KL           = rep(as.numeric(NA), L),
    lbf          = rep(as.numeric(NA), L),
    lbf_variable = matrix(as.numeric(NA), L, p),
    sigma2       = residual_variance,
    pi           = prior_weights
  )

  return(mat_init)
}

# Get variance of y
get_var_y.individual <- function(data, ...) {
  return(var(drop(data$y)))
}
