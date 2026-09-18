#' @rdname slider_cs_table
#' @method coef susie_slide
#' @export
coef.susie_slide <- function(object,...) {
  scale <- object$X_column_scale_factors
  ans <- cbind(additive=c(object$intercept,colSums(object$alpha*object$mu)/scale),
               heterozygote=c(0,colSums(object$alpha*object$mu*object$delta)/scale))
  snps <- colnames(object$alpha)
  if(is.null(snps)) snps <- paste0("SNP",seq_len(ncol(object$alpha)))
  rownames(ans) <- c("(Intercept)",snps)
  if(object$null_index>0) ans <- ans[-(object$null_index+1L),,drop=FALSE]
  ans
}
#' @rdname slider_cs_table
#' @method predict susie_slide
#' @export
predict.susie_slide <- function(object,newx=NULL,type=c("response","coefficients"),...) {
  type <- match.arg(type)
  if(type=="coefficients") {
    if(!is.null(newx)) stop("Do not supply newx when requesting coefficients.")
    return(coef(object))
  }
  if(is.null(newx)) return(object$fitted)
  if((!is.matrix(newx) && !inherits(newx,"sparseMatrix")) || ncol(newx)!=object$input_p)
    stop("newx must have the same SNP columns, in the same order, as the training genotypes.")
  values <- if(inherits(newx,"sparseMatrix")) newx@x else newx
  if(any(!is.finite(values)) || any(!values %in% 0:2)) stop("newx must contain finite hard-call genotypes.")
  cf <- coef(object)
  if(!is.null(colnames(newx)) && !is.null(colnames(object$alpha)) &&
     !identical(colnames(newx),rownames(cf)[-1])) stop("newx SNP names/order differ from the training data.")
  result <- as.numeric(cf[1,1]+newx %*% cf[-1,1])
  active <- which(cf[-1,2]!=0)
  if(length(active)) {
    cap <- max(1,floor(64e6/(8*max(1,nrow(newx)))))
    for(k in .chunks(length(active),cap)) {
      j <- active[k]
      result <- result+as.numeric(((newx[,j,drop=FALSE]==1)*1) %*% cf[j+1,2])
    }
  }
  result
}
#' Slider coefficients, predictions, and credible-set output
#'
#' Extracts output using both the additive and heterozygote terms.
#' @param fit,object A susieSlide fit.
#' @param newx Complete hard-call genotypes with the same SNPs and allele
#'   orientation, in the same order. NULL returns training fitted values.
#' @param type Return response predictions or coefficients.
#' @param ... Reserved additional arguments.
#' @return coef returns a two-column matrix, additive and heterozygote,
#'   including an intercept row. Predictions equal the intercept plus X
#'   times additive coefficients plus I(X=1) times heterozygote coefficients.
#'   predict returns predictions or that coefficient matrix. slider_cs_table
#'   returns one row per SNP per reported credible set, with the component
#'   index, delta, alpha, SNP PIP, and count-forced status, matching
#'   fit$delta_cs$summary. The fit also stores delta_cs$delta, a matrix with
#'   one row per reported CS and one column per input SNP (including SNPs
#'   outside that CS, excluding an explicit null column). Rows follow
#'   sets$cs and select component rows using sets$cs_index. summary returns
#'   the usual SuSiE summary with an additional delta table.
#' @importFrom stats setNames
#' @export
slider_cs_table <- function(fit) {
  empty <- data.frame(cs=character(),component=integer(),snp_index=integer(),
                      snp=character(),alpha=numeric(),pip=numeric(),delta=numeric(),
                      delta_forced=logical())
  if(is.null(fit$sets$cs) || !length(fit$sets$cs)) return(empty)
  snps <- colnames(fit$alpha)
  if(is.null(snps)) snps <- paste0("SNP",seq_len(ncol(fit$alpha)))
  rows <- lapply(seq_along(fit$sets$cs),function(k) {
    l <- fit$sets$cs_index[k]; j <- fit$sets$cs[[k]]
    data.frame(cs=names(fit$sets$cs)[k],component=l,snp_index=j,snp=snps[j],
      alpha=fit$alpha[l,j],pip=fit$pip[j],delta=fit$delta[l,j],
      delta_forced=unname(fit$delta_forced[j]),row.names=NULL)
  })
  do.call(rbind,rows)
}

.slider_cs_output <- function(fit) {
  components <- if(length(fit$sets$cs)) fit$sets$cs_index else integer(0)
  snps <- seq_len(ncol(fit$delta))
  if(!is.null(fit$null_index) && fit$null_index>0)
    snps <- snps[snps!=fit$null_index]
  deltas <- fit$delta[components,snps,drop=FALSE]
  rownames(deltas) <- names(fit$sets$cs)
  snp_names <- colnames(fit$alpha)
  if(is.null(snp_names)) snp_names <- paste0("SNP",seq_len(ncol(fit$alpha)))
  colnames(deltas) <- snp_names[snps]
  list(summary=slider_cs_table(fit),delta=deltas)
}

#' @rdname slider_cs_table
#' @method summary susie_slide
#' @export
summary.susie_slide <- function(object,...) {
  # Use the familiar summary layout; avoid the upstream null-column row drop.
  temporary <- object
  temporary$null_index <- 0
  result <- summary.susie(temporary,...)
  result$delta <- slider_cs_table(object)
  result
}
