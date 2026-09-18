#' @export
#' @noRd
initialize_susie_model.slide_individual <- function(data,params,var_y,...) {
  model <- .engine("initialize_susie_model.individual")(data,params,var_y,...)
  model$delta <- matrix(0,nrow(model$alpha),data$p)
  model$slide_s <- matrix(data$xx,nrow(model$alpha),data$p,byrow=TRUE)
  model$component_fitted <- matrix(0,data$n,nrow(model$alpha))
  class(model) <- c("susie_slide","susie")
  model
}
#' @export
#' @noRd
ibss_initialize.slide_individual <- function(data,params) {
  model <- .engine("ibss_initialize.default")(data,params)
  if(!is.null(data$warm)) {
    warm <- data$warm
    if(!inherits(warm,"susie_slide") || !identical(dim(warm$alpha),dim(model$alpha)) ||
       !identical(dim(warm$delta),dim(model$delta)) ||
       !isTRUE(all.equal(unname(warm$X_column_scale_factors),unname(data$scale))))
      stop("model_init must be a susieSlide fit with matching L, SNPs, and genotype scales.")
    for(nm in c("alpha","mu","mu2","V","delta")) model[[nm]] <- warm[[nm]]
    if(isTRUE(params$estimate_residual_variance)) model$sigma2 <- warm$sigma2
    model$delta[,data$forced] <- 0
    if(!is.null(data$fixed_delta)) {
      model$delta <- data$fixed_delta
      model$delta[,data$forced] <- 0
    }
    for(l in seq_len(nrow(model$alpha))) {
      d <- model$delta[l,]
      model$slide_s[l,] <- pmax(0,data$xx+2*d*data$xh+d^2*data$hh)
      model$component_fitted[,l] <- .slide_product(data,model$alpha[l,]*model$mu[l,],
                                                  model$alpha[l,]*model$mu[l,]*d)
    }
    model$Xr <- rowSums(model$component_fitted)
  }
  model
}
#' @export
#' @noRd
compute_residuals.slide_individual <- function(data,params,model,l,...) {
  model$fitted_without_l <- model$Xr-model$component_fitted[,l]
  model$raw_residuals <- data$y-model$fitted_without_l
  model$residuals <- .engine("compute_Xty")(data$X,model$raw_residuals)
  model$hresiduals <- .h_score(data,model$raw_residuals)
  model$residual_variance <- model$sigma2
  model
}
#' @export
#' @noRd
compute_ser_statistics.slide_individual <- function(data,params,model,l,...) {
  stats <- .engine("compute_ser_statistics.individual")(data,params,model,l,...)
  stats$fixed_delta <- if(is.null(data$fixed_delta)) rep(NA_real_,data$p) else data$fixed_delta[l,]
  stats$fixed_delta[data$forced] <- 0
  stats
}
.slide_ser <- function(data,model,V,ser_stats) {
  ans <- slide_ser_native(as.double(data$xx),as.double(data$xh),as.double(data$hh),
    as.double(model$residuals),as.double(model$hresiduals),as.double(V),
    as.double(model$sigma2),as.double(ser_stats$fixed_delta))
  colnames(ans) <- c("delta","lbf","mu","mu2","s","t","variance")
  ans
}
#' @export
#' @noRd
loglik.slide_individual <- function(data,params,model,V,ser_stats,l=NULL,...) {
  ser <- .slide_ser(data,model,V,ser_stats)
  # Use the engine's normalization, including its tiny prior-weight offset,
  # so the all-additive path matches this exact susieR version.
  updated <- .engine("apply_ser_lbf")(model,ser[,"lbf"],model$sigma2/ser[,"s"],l)
  if(is.null(l)) return(updated$lbf_model)
  model <- updated$model
  model$delta[l,] <- ser[,"delta"]
  model$slide_s[l,] <- ser[,"s"]
  model$slide_ser <- ser
  model
}
#' @export
#' @noRd
neg_loglik.slide_individual <- function(data,params,model,V_param,ser_stats,...) {
  V <- if(ser_stats$optim_scale=="log") exp(V_param) else V_param
  -loglik.slide_individual(data,params,model,V,ser_stats)
}
#' @export
#' @noRd
calculate_posterior_moments.slide_individual <- function(data,params,model,V,l,...) {
  model$mu[l,] <- model$slide_ser[,"mu"]
  model$mu2[l,] <- model$slide_ser[,"mu2"]
  model
}
#' @export
#' @noRd
SER_posterior_e_loglik.slide_individual <- function(data,params,model,l) {
  d <- model$delta[l,]
  linear <- sum(model$alpha[l,]*model$mu[l,]*(model$residuals+d*model$hresiduals))
  quadratic <- sum(model$alpha[l,]*model$mu2[l,]*model$slide_s[l,])
  -data$n/2*log(2*pi*model$sigma2) -
    (sum(model$raw_residuals^2)-2*linear+quadratic)/(2*model$sigma2)
}
#' @export
#' @noRd
compute_kl.slide_individual <- function(data,params,model,l) {
  # Same conjugate SER KL identity as SuSiE, conditional on fitted delta.
  linear <- sum(model$alpha[l,]*model$mu[l,]*
                  (model$residuals+model$delta[l,]*model$hresiduals))
  quadratic <- sum(model$alpha[l,]*model$mu2[l,]*model$slide_s[l,])
  model$KL[l] <- -model$lbf[l]+(2*linear-quadratic)/(2*model$sigma2)
  model
}
#' @export
#' @noRd
post_loglik_prior_hook.slide_individual <- function(data,params,model,ser_stats,l,V_init) {
  if(params$estimate_prior_method!="EM") return(list(V=V_init,model=model))
  V <- sum(model$alpha[l,]*model$mu2[l,])
  # Refresh the posterior/evidence/KL at the new variance, rather than leaving
  # the returned evidence and moments at the preceding EM variance.
  model <- loglik.slide_individual(data,params,model,V,ser_stats,l)
  model <- calculate_posterior_moments.slide_individual(data,params,model,V,l)
  model <- compute_kl.slide_individual(data,params,model,l)
  list(V=V,model=model)
}
#' @export
#' @noRd
update_fitted_values.slide_individual <- function(data,params,model,l,...) {
  b <- model$alpha[l,]*model$mu[l,]
  model$component_fitted[,l] <- .slide_product(data,b,b*model$delta[l,])
  model$Xr <- model$fitted_without_l+model$component_fitted[,l]
  model
}
#' @export
#' @noRd
get_ER2.slide_individual <- function(data,model) {
  sum((data$y-model$Xr)^2)+sum(model$alpha*model$mu2*model$slide_s)-
    sum(model$component_fitted^2)
}
#' @export
#' @noRd
Eloglik.slide_individual <- function(data,model) {
  -data$n/2*log(2*pi*model$sigma2)-get_ER2.slide_individual(data,model)/(2*model$sigma2)
}
#' @export
#' @noRd
trim_null_effects.slide_individual <- function(data,params,model) {
  model <- .engine("trim_null_effects.default")(data,params,model)
  null <- which(model$V==0)
  if(length(null)) {
    model$delta[null,] <- 0
    model$component_fitted[,null] <- 0
    model$Xr <- rowSums(model$component_fitted)
  }
  model
}
#' @export
#' @noRd
get_intercept.slide_individual <- function(data,params,model,...) {
  if(!params$intercept) return(0)
  data$mean_y-sum(data$xmean*colSums(model$alpha*model$mu)/data$scale)-
    sum(data$hmean*colSums(model$alpha*model$mu*model$delta)/data$scale)
}
#' @export
#' @noRd
get_cs.slide_individual <- function(data,params,model,...) {
  if(is.null(params$coverage)) return(NULL)
  if(all(model$delta==0) && !is.null(params$min_abs_corr))
    return(.engine("get_cs.individual")(data,params,model,...))
  cs <- list(); purity <- list(); coverage <- numeric(0); indices <- integer(0)
  for(l in which(model$V>params$prior_tol)) {
    order <- order(model$alpha[l,],decreasing=TRUE)
    k <- which(cumsum(model$alpha[l,order])>=params$coverage)[1]
    if(is.na(k)) k <- length(order)
    j <- sort(order[seq_len(k)])
    if((model$null_index>0 && model$null_index %in% j) ||
       any(vapply(cs,identical,logical(1),j))) next
    # Correlations must use THIS component's fitted slider for each SNP.
    jpurity <- j
    cap <- .engine("resolve_n_purity")(params$n_purity,data$n,length(j))
    if(length(jpurity)>cap) jpurity <- sample(jpurity,cap)
    z <- as.matrix(data$X[,jpurity,drop=FALSE])+sweep(as.matrix(.H(data,jpurity)),
                                                  2,model$delta[l,jpurity],"*")
    pur <- .engine("get_purity")(seq_along(jpurity),z,NULL,n=-1)
    keep <- (is.null(params$min_abs_corr) && is.null(params$median_abs_corr)) ||
      (!is.null(params$min_abs_corr) && pur[1]>=params$min_abs_corr) ||
      (!is.null(params$median_abs_corr) && pur[3]>=params$median_abs_corr)
    if(!keep) next
    cs[[length(cs)+1L]] <- j
    purity[[length(purity)+1L]] <- pur
    coverage <- c(coverage,sum(model$alpha[l,j])); indices <- c(indices,l)
  }
  if(!length(cs)) return(list(cs=NULL,coverage=NULL,requested_coverage=params$coverage))
  purity <- as.data.frame(do.call(rbind,purity))
  names(purity) <- c("min.abs.corr","mean.abs.corr","median.abs.corr")
  names(cs) <- rownames(purity) <- paste0("L",indices)
  ordering <- order(purity[,if(is.null(params$min_abs_corr)) 3 else 1],decreasing=TRUE)
  list(cs=cs[ordering],purity=purity[ordering,,drop=FALSE],cs_index=indices[ordering],
       coverage=coverage[ordering],requested_coverage=params$coverage)
}
#' @export
#' @noRd
get_variable_names.slide_individual <- function(data,model,...) {
  model <- .engine("get_variable_names.individual")(data,model,...)
  dimnames(model$delta) <- dimnames(model$alpha)
  model
}
#' @export
#' @noRd
cleanup_model.slide_individual <- function(data,params,model,...) {
  model$min_obs <- data$min_obs
  model$genotype_counts <- data$counts
  rownames(model$genotype_counts) <- colnames(model$alpha)
  model$delta_forced <- setNames(data$forced,colnames(model$alpha))
  model$input_p <- data$input_p
  model$mu_delta <- model$mu*model$delta
  model$expected_squared_residuals <- get_ER2.slide_individual(data,model)
  for(nm in c("component_fitted","slide_ser","slide_s","hresiduals")) model[[nm]] <- NULL
  .engine("cleanup_model.individual")(data,params,model,...)
}
