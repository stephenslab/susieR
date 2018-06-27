#' @title Local false sign rate (lfsr) for susie confidence sets
#' @details This computes the average lfsr across SNPs for each l, weighted by the
#' posterior inclusion probability alpha
#' @param fitted a susie fit, the output of `susieR::susie()`
#' @return an l vector of lfsr for confidence sets
#' @export
susie_get_lfsr = function(res){
  pos_prob = pnorm(0,mean=t(res$mu),sd=sqrt(res$mu2-res$mu^2))
  neg_prob = 1-pos_prob
  1-rowSums(res$alpha*t(pmax(pos_prob,neg_prob)))
}

#find how many variables in the CS
# x is a probability vector
n_in_CS_x = function(x, coverage = 0.9){
  sum(cumsum(sort(x,decreasing = TRUE))<coverage)+1
}

# return binary vector indicating if each point is in CS
# x is a probability vector
in_CS_x = function(x, coverage = 0.9){
  n = n_in_CS_x(x, coverage)
  o = order(x,decreasing=TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}

# returns an l by p binary matrix
# indicating which variables are in susie confidence sets
in_CS = function(res, coverage = 0.9){
  if (class(res) == "susie")
    res = res$alpha
  t(apply(res,1,function(x) in_CS_x(x, coverage)))
}

n_in_CS = function(res, coverage = 0.9){
  if (class(res) == "susie")
    res = res$alpha
  apply(res,1,function(x) n_in_CS_x(x, coverage))
}

get_purity = function(pos, corr) {
  value = abs(corr[pos, pos])
  c(min(value), mean(value), median(value))
}

#' @title Extract confidence sets from fitted SuSiE model
#' @details It reports indices of variables in each confidence set identified,
#' as well as summaries of correlation between variables within each set.
#' @param fitted a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param X N by P matrix of variables.
#' When provided, correlation between variables will be computed and used to remove
#' confidence sets whose minimum correlation between variables is smaller than `min_abs_corr` (see below).
#' @param Xcorr P by P matrix of correlations between variables.
#' when provided, it will be used to remove confidence sets
#' whose minimum correlation between variables is smaller than `min_abs_corr` (see below).
#' @param coverage coverage of confident sets. Default to 0.9 for 90% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @return a list of `cs`, and additionally `purity` and selected `cs_index` if `X` or `Xcorr` was provided.
#' @export
susie_get_CS = function(fitted,
                        X = NULL, Xcorr = NULL,
                        coverage = 0.9,
                        min_abs_corr = 0.5) {
  if (class(fitted) == "susie")
    fitted = fitted$alpha
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  if (!is.null(X) && ncol(X) != length(fitted[1,])) {
    stop("Column of X matrix does not match length of fitted posterior")
  }
  if (!is.null(Xcorr) && ncol(Xcorr) != length(fitted[1,])) {
    stop("Dimension of Xcorr matrix does not match length of fitted posterior")
  }
  if (!is.null(Xcorr) && !isSymmetric(Xcorr)) {
    stop("Xcorr matrix must be symmetric")
  }
  if (!is.null(X)) Xcorr = cor(X)
  # L by P binary matrix
  status = in_CS(fitted, coverage)
  # an L list of CS positions
  cs = lapply(1:nrow(status), function(i) which(status[i,]!=0))
  cs = cs[lapply(cs, length) > 0]
  # compute and filter by "purity"
  if (is.null(Xcorr)) {
    return(list(cs=cs))
  } else {
    purity = data.frame(do.call(rbind, lapply(1:length(cs), function(i) get_purity(cs[[i]], Xcorr))))
    colnames(purity) = c('min.abs.corr', 'mean.abs.corr', 'median.abs.corr')
    is_pure = which(purity$min.abs.corr > min_abs_corr)
    return(list(cs = cs[is_pure], purity = purity[is_pure,], cs_index = is_pure))
  }
}

#' @title Compute posterior inclusion probability (PIP) for all variables
#' @param fitted a susie fit, the output of `susieR::susie()`, or simply the posterior
#' inclusion probability matrix alpha.
#' @param include_index index of single effect models to consider when calculating PIP. Default to NULL.
#' @return a vector of posterior inclusion probability.
#' @export
susie_get_PIP = function(fitted, include_index = NULL) {
  if (class(fitted) == "susie")
    fitted = fitted$alpha
  if (!is.null(include_index)) {
    fitted = fitted[include_index,,drop=FALSE]
  }
  return(as.vector(1 - apply(1 - fitted, 2, prod)))
}

# compute standard error for regression coef
calc_stderr = function(X, residuals) {
    # S = (X'X)^-1 \Sigma
    sqrt(diag(sum(residuals^2) / (nrow(X) - 2) * chol2inv(chol(t(X) %*% X))))
  }

# univariate regression between each column of X and y
# Remove covariates if Z is not NULL
univariate_regression = function(X, y, Z=NULL, return_residue=FALSE) {
  if (!is.null(Z)) {
    y = .lm.fit(Z, y)$residuals
  }
  output = do.call(rbind,
                   lapply(1:ncol(X), function(i) {
                     g = .lm.fit(cbind(1, X[,i]), y)
                     return(c(coef(g)[2], calc_stderr(cbind(1, X[,i]), g$residuals)[2]))
                   })
                   )
  if (return_residue) {
    return(list(betahat = output[,1], sebetahat = output[,2],
                residuals = y))
  } else {
    return(list(betahat = output[,1], sebetahat = output[,2]))
  }
}

# computes z score (t-statistic) for association between each
# column of X and y
calc_z = function(X,y){
  out = univariate_regression(X,y)
  return(out$betahat/out$sebetahat)
}

# plot p values of data and color in the CSs
# for simulated data, specify b = true effects (highlights in red)
pplot = function(X,y,res,pos=NULL,b=NULL,CSmax = 400,...){
  z = calc_z(X,y)
  zneg = -abs(z)
  logp = log10(pnorm(zneg))
  if(is.null(b)){b = rep(0,ncol(X))}
  if(is.null(pos)){pos = 1:ncol(X)}
  plot(pos,-logp,col="grey",xlab="",ylab="-log10(p)",...)
  points(pos[b!=0],-logp[b!=0],col=2,pch=16)
  for(i in 1:nrow(res$alpha)){
    if(n_in_CS(res)[i]<CSmax)
      points(pos[which(in_CS(res)[i,]>0)],-logp[which(in_CS(res)[i,]>0)],col=i+2)
  }

}


# return residuals from Y after removing susie fit
get_R = function(X,Y,s){
  Y- X %*%  coef(s)
}

# get number of iterations
susie_get_niter <- function(res) {
  return(res$niter)
}
