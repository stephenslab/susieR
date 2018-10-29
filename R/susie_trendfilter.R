#susie_trendfilter

#' @title perform trend filtering using SuSiE
#' @param Y an n vector
#' @param order a scalar for the order of trend filtering
#' @return SuSiE fit for trend filtering
#' @examples
#' susie_trendfilter(y, 0)
#' susie_trendfilter(y, 1, L=20)
#' susie_trenndfilter(y, 0, estimate_prior_variance = TRUE)
#' @export
susie_trendfilter = function(Y, order, ...){
  X = create_Dinv(order, length(Y))
  class(X) = "tfmatrix"
  attr(X, "order") = order
  s = susie(X=X, Y=Y, ...)
  return(s)
}

