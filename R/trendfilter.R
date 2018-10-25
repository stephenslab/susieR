#susie_trendfilter

#' @title perform trend filtering using SuSiE
#' @param Y an n vector
#' @param order a scalar for the order of trend filtering
#' @return SuSiE fit for trend filtering
#' @export
susie_trendfilter = function(Y, order, ...){
  #X = diag(length(Y))
  X = create_Dinv(order, length(Y))
  class(X) = "tfmatrix"
  attr(X, "order") = order
  #attr(X, "d") <- compute_colSumsDinv2(order, Y)
  #attr(X, "scaled:center") = rep(0, length(Y))
  #attr(X, "scaled:scale") = rep(1, length(Y))
  s = susie(X=X, Y=Y, ...)
  return(s)
}

