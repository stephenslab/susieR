#' @title Compute posterior samples from susie_fit.
#'
#' @param susie_fit A susie fit, the output form `susieR::susie()`
#'
#' @param num_samples the number of samples to be drawn from the posterior distribution
#'
#' @return The return value is a list containing the effect sizes samples and causal status samples
#'    components:
#'
#'    \item{b}{ num_variables x num_samples matrix of effect sizes draw}
#'
#'    \item{gamma}{ num_variables  x num_samples matrix of causal status draw}
#'
#' @export
#'

susie_get_posterior_samples <- function(susie_fit, num_samples){
    # removed effects having estimated prior variance equals zero
    if (is.numeric(susie_fit$V)) {
        include_idx = which(susie_fit$V > 1E-9)
    }
    else {
        include_idx = 1:nrow(susie_fit$alpha)
    }
    
    posterior_mean <- sweep(susie_fit$mu, 2, susie_fit$X_column_scale_factors, "/")
    posterior_sd <- sweep(sqrt(susie_fit$mu2 - (susie_fit$mu)^2), 2, susie_fit$X_column_scale_factors, "/")

    pip <- susie_fit$alpha
    L = nrow(pip)
    num_snps <- ncol(pip)
    b_samples <- matrix(NA, num_snps, num_samples)
    gamma_samples <- matrix(NA, num_snps, num_samples)
    for (sample_i in 1 : num_samples){
        b <- 0
        for (l in include_idx){
            gamma_l <- rmultinom(1, 1, pip[l, ])
            effect_size <- rnorm(1, mean=posterior_mean[l, which(gamma_l != 0)], 
                                    sd=posterior_sd[l, which(gamma_l != 0)])
            b_l <- gamma_l * effect_size
            b <- b + b_l
        }
        b_samples[, sample_i] <- b
        gamma_samples[, sample_i] <- as.numeric(b != 0)
    }
    return(list(b=b_samples,
                gamma=gamma_samples))
}