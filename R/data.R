#' @name N2finemapping
#'
#' @title Simulated Fine-mapping Data with Two Effect Variables
#'
#' @docType data
#'
#' @description This data set contains a genotype matrix for 574
#'   individuals and 1,002 variables. The variables are genotypes after
#'   centering and scaling, and therefore retain the correlation
#'   structure of the original genotype data. Two of the variables have
#'   non-zero effects on the multivariate response. The response data
#'   are generated under a multivariate linear regression model. See
#'   Wang \emph{et al} (2020) for details.
#'
#' @format \code{N2finemapping} is a list with the following elements:
#'
#' \describe{
#'
#'   \item{X}{Centered and scaled genotype data.}
#'
#'   \item{chrom}{Chromomsome of the original data, in hg38 coordinates.}
#'
#'   \item{pos}{Chromomosomal position of the original data, in hg38
#'     coordinates. The information can be used to compare impact of using
#'     other genotype references of the same variables in \code{susie_rss}
#'     application.}
#'
#'   \item{true_coef}{Simulated effect sizes.}
#'
#'   \item{residual_variance}{Simulated residual covariance matrix.}
#'
#'   \item{Y}{Simulated multivariate response.}
#'
#'   \item{allele_freq}{Allele frequencies based on the original
#'     genotype data.}
#'
#'   \item{V}{Suggested prior covariance matrix for effect sizes of
#'      the two non-zero effect variables.}
#' }
#'
#' @keywords data
#'
#' @references
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#'   new approach to variable selection in regression, with application
#'   to genetic fine-mapping. \emph{Journal of the Royal Statistical
#'   Society, Series B} \url{https://doi.org/10.1101/501114}.
#'
#' @examples
#' data(N2finemapping)
NULL

#' @name N3finemapping
#'
#' @title Simulated Fine-mapping Data with Three Effect Variables.
#'
#' @docType data
#'
#' @description The data-set contains a matrix of 574
#' individuals and 1,001 variables. These variables are real-world
#' genotypes centered and scaled, and therefore retains the
#' correlation structure of variables in the original genotype data. 3
#' out of the variables have non-zero effects.  The response data is
#' generated under a multivariate linear regression model.  See Wang
#' \emph{et al} (2020) for more details.
#'
#' @format \code{N3finemapping} is a list with the following elements:
#'
#' \describe{
#'
#'   \item{X}{N by P variable matrix of centered and scaled genotype
#' data.}
#'
#'   \item{chrom}{Chromomsome of the original data, in hg38 coordinate.}
#'
#'   \item{pos}{Chromomosomal positoin of the original data, in hg38
#' coordinate. The information can be used to compare impact of using
#' other genotype references of the same variables in susie_rss
#' application.}
#'
#'   \item{true_coef}{The simulated effect sizes.}
#'
#'   \item{residual_variance}{The simulated residual covariance matrix.}
#'
#'   \item{Y}{The simulated response variables.}
#'
#'   \item{allele_freq}{Allele frequency of the original genotype data.}
#'
#'   \item{V}{Prior covariance matrix for effect size of the three
#' non-zero effect variables.}  }
#'
#' @keywords data
#'
#' @references
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#'   new approach to variable selection in regression, with application
#'   to genetic fine-mapping. \emph{Journal of the Royal Statistical
#'   Society, Series B} \url{https://doi.org/10.1101/501114}.
#'
#' @examples
#' data(N3finemapping)
NULL

#' @name FinemappingConvergence
#'
#' @title Simulated Fine-mapping Data with Convergence Problem.
#'
#' @description The data is simulated using real genotypes
#' from 50000 individuals and 1001 SNPs. Two of the variables have
#' non-zero effects on the multivariate response. The response data
#' are generated under a linear regression model.
#' The genotypes and simulated response are column centered.
#'
#' @format \code{FinemappingConvergence} is a list with the following
#' elements:
#'
#' \describe{
#'
#'   \item{XtX}{XtX is computed using centered and scaled genotype
#' data.}
#'
#'   \item{Xty}{XtX is computed using centered and scaled genotype
#' data and centered simulated response.}
#'
#'   \item{yty}{yty is computed using centered simulated response.}
#'
#'   \item{n}{Sample size.}
#'
#'   \item{true_coef}{The simulated effect sizes.}
#'
#'   \item{z}{Marginal z scores from simple linear regression.}}
#'
#' @docType data
#'
#' @keywords data
#'
#' @examples
#' data(FinemappingConvergence)
NULL
