#' splits a vector into a list of vectors, each containing a specified number
#' of vectors of the orignal matrix
#' @param b vector to be split
#' @param nvar vector of integers saying how many columns of b to include in each
#' @examples
#' b = c(rep(1,3),rep(2,4),rep(3,3))
#' split_vector(b,c(3,4,3))
split_vector = function(b,nvar){
  split(b,rep(1:length(nvar), nvar))
}
