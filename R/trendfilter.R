#susie_trendfilter

#' @title perform trend filtering using SuSiE
#' @param Y an n vector
#' @param order a scalar for the order of trend filtering
#' @return SuSiE fit for trend filtering
susie_trendfilter = function(Y, order, ...){
  # standardize = FALSE while trend filtering
  X = diag(length(Y))
  attr(X, "d") <- compute_colSumsDinv2(order, Y)
  attr(X, "scaled:center") = rep(0, length(Y))
  attr(X, "scaled:scale") = rep(1, length(Y))
  attr(X, "order") = order
  s = susie(X=X, Y=Y, trendfiltering=TRUE, order=order, ...)
  return(s)
}