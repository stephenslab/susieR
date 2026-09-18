# Look up the shared engine in this package, including during load_all().
#' @importFrom stats cov2cor uniroot
#' @importFrom utils tail
#' @noRd
.engine <- function(name) get(name,envir=environment(.engine),inherits=FALSE)
