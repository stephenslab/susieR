# Performs sum of single-effect (SuSiE) linear regression with z
# scores (with lambda). The summary data required are the p by p
# correlation matrix R, the p vector z. The summary stats should come
# from the same individuals. This function fits the regression model z
# = sum_l Rb_l + e, where e is N(0,residual_variance * R + lambda I)
# and the sum_l b_l is a p vector of effects to be estimated. The
# assumption is that each b_l has exactly one non-zero element, with
# all elements equally likely to be non-zero. The prior on the
# non-zero element is N(0,var = prior_variance).
#
#' @importFrom stats optimize
susie_rss_lambda = function(z, R, maf = NULL, maf_thresh = 0,
                            L = 10, lambda = 0,
                            prior_variance = 50, residual_variance = NULL,
                            r_tol = 1e-08, prior_weights = NULL,
                            null_weight = 0,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE,
                            estimate_prior_method = c("optim", "EM", "simple"),
                            check_null_threshold = 0, prior_tol = 1e-9,
                            max_iter = 100, s_init = NULL, intercept_value = 0,
                            coverage = 0.95, min_abs_corr = 0.5,
                            tol = 1e-3, verbose = FALSE, track_fit = FALSE,
                            check_R = TRUE, check_z = FALSE) {

  # Check input R.
  if (nrow(R) != length(z))
    stop(paste0("The dimension of correlation matrix (",nrow(R)," by ",
                ncol(R),") does not agree with expected (",length(z)," by ",
                length(z),")"))
  if (!is_symmetric_matrix(R))
    stop("R is not a symmetric matrix")
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix or a sparse matrix")
  
  # MAF filter.
  if (!is.null(maf)) {
    if (length(maf) != length(z))
      stop(paste0("The length of maf does not agree with expected ",length(z)))
    id = which(maf > maf_thresh)
    R = R[id,id]
    z = z[id]
  }

  if (any(is.infinite(z)))
    stop("z contains infinite values")

  # Check for NAs in R.
  if (anyNA(R))
    stop("R matrix contains missing values")

  # Replace NAs in z with zero.
  if (anyNA(z)) {
    warning_message("NA values in z-scores are replaced with 0")
    z[is.na(z)] = 0
  }

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R)*(1-null_weight),ncol(R)),null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight),null_weight)
    R = cbind(rbind(R,0),0)
    z = c(z,0)
  }

  # Eigen decomposition for R, filter on eigenvalues.
  p = ncol(R)
  attr(R,"eigen") = eigen(R,symmetric = TRUE)
  if (check_R && any(attr(R,"eigen")$values < -r_tol))
    stop(paste0("The correlation matrix (",nrow(R)," by ",ncol(R),
                "is not a positive semidefinite matrix. ",
                "The smallest eigenvalue is ",min(attr(R,"eigen")$values),
                ". You can bypass this by \"check_R = FALSE\" which instead ",
                "sets negative eigenvalues to 0 to allow for continued ",
                "computations."))

  # Check whether z in space spanned by the non-zero eigenvectors of R.
  if (check_z) {
    proj = check_projection(R,z)
    if (!proj$status)
      warning_message("Input z does not lie in the space of non-zero eigenvectors ",
              "of R.")
    else
      message("Input z is in space spanned by the non-zero eigenvectors of ",
              "R.")
  }
  R = set_R_attributes(R,r_tol)
  
  if (lambda == "estimate"){
    colspace = which(attr(R,"eigen")$values > 0)
    if (length(colspace) == length(z))
      lambda = 0
    else {
      znull = crossprod(attr(R,"eigen")$vectors[,-colspace], z) # U2^T z
      lambda = sum(znull^2)/length(znull)
    }
  }

  # Initialize susie fit.
  s = init_setup_rss(p,L,prior_variance,residual_variance,prior_weights,
                     null_weight)
  if (!missing(s_init) && !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")
    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init)
    num_effects = nrow(s_init$alpha)
    if(missing(L)){
      L = num_effects
    }else if(L < num_effects){
      warning_message(paste("Specified number of effects L =",L,
                    "is smaller than the number of effects",num_effects,
                    "in input SuSiE model. It will have",num_effects,"effects."))
      L = num_effects
    }
    # expand s_init if L > num_effects.
    s_init = susie_prune_single_effects(s_init, L, s$V)
    s = modifyList(s,s_init)
    s = init_finalize_rss(s,R = R)
  } else
    s = init_finalize_rss(s)

  s$sigma2 = s$sigma2 - lambda
  estimate_prior_method = match.arg(estimate_prior_method)

  # Intialize elbo to NA.
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()

  attr(R,"lambda") = lambda
  Sigma = update_Sigma(R,s$sigma2,z)  # sigma2*R + lambda*I

  for (i in 1:max_iter) {
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    s = update_each_effect_rss(R,z,s,Sigma,estimate_prior_variance,
                               estimate_prior_method,check_null_threshold)
    if (verbose)
      print(paste0("before estimate sigma2 objective:",
                   get_objective_rss(R,z,s)))

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance, before the update.
    elbo[i+1] = get_objective_rss(R,z,s)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }
    if (estimate_residual_variance) {
      if (lambda == 0) {
        est_sigma2 = (1/sum(attr(R,"eigen")$values != 0))*get_ER2_rss(1,R,z,s)
        if (est_sigma2 < 0)
          stop("Estimating residual variance failed: the estimated value ",
               "is negative")
        if (est_sigma2 > 1)
          est_sigma2 = 1
      } else {
        est_sigma2 = optimize(Eloglik_rss, interval = c(1e-4, 1-lambda),
                              R = R, z = z, s = s, maximum = TRUE)$maximum
        if(Eloglik_rss(est_sigma2, R, z, s) < Eloglik_rss(1-lambda, R, z, s)){
          est_sigma2 = 1-lambda
        }
      }
      s$sigma2 = est_sigma2

      if (verbose)
        print(paste0("after estimate sigma2 objective:",
                     get_objective_rss(R,z,s)))
      Sigma = update_Sigma(R,s$sigma2,z)
    }
  }

  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(i+1)]
  s$elbo = elbo
  s$niter = i
  s$lambda = lambda

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  s$intercept = intercept_value
  s$fitted = s$Rz

  s$X_column_scale_factors = attr(R,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(R)
    s$sets = susie_get_cs(s,coverage = coverage,Xcorr = R,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }
  if (!is.null(names(z))) {
    variable_names = names(z)
    if (!is.null(null_weight))
      variable_names = c("null",variable_names)
    colnames(s$alpha) = variable_names
    colnames(s$mu) = variable_names
    colnames(s$mu2) = variable_names
    colnames(s$lbf_variable) = variable_names
    names(s$pip) = variable_names
  }

  return(s)
}

update_Sigma = function (R, sigma2, z) {
  Sigma = sigma2*R + attr(R,"lambda") * diag(length(z))
  eigenS = attr(R,"eigen")
  eigenS$values = sigma2 * eigenS$values + attr(R,"lambda")

  Dinv = 1/(eigenS$values)
  Dinv[is.infinite(Dinv)] = 0
  attr(Sigma,"eigenS") = eigenS

  # Sigma^(-1) R_j = U (sigma2 D + lambda)^(-1) D U^T e_j
  attr(Sigma,"SinvRj") = eigenS$vectors %*% (Dinv * attr(R,"eigen")$values *
                           t(eigenS$vectors))

  if (attr(R,"lambda") == 0)
    attr(Sigma,"RjSinvRj") = attr(R,"d")/sigma2
  else {
    tmp = t(eigenS$vectors)
    attr(Sigma,"RjSinvRj") =
      colSums(tmp * (Dinv*(attr(R,"eigen")$values^2) * tmp))
  }

  return(Sigma)
}
