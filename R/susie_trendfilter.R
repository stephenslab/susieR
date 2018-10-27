#susie_trendfilter

#' @title perform trend filtering using SuSiE
#' @param Y an n vector
#' @param order a scalar for the order of trend filtering
#' @return SuSiE fit for trend filtering
#' @export
susie_trendfilter = function(Y, order, ...){
  X = create_Dinv(order, length(Y))
  class(X) = "tfmatrix"
  attr(X, "order") = order
  s = susie(X=X, Y=Y, ...)
  return(s)
}

