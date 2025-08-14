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
#'   Society, Series B} \doi{10.1101/501114}.
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
#' generated under a multivariate linear regression model. See Wang
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
#'   Society, Series B} \doi{10.1101/501114}.
#'
#' @examples
#' data(N3finemapping)
NULL

#' @name FinemappingConvergence
#'
#' @title Simulated Fine-mapping Data with Convergence Problem.
#'
#' @description Data simulated using real genotypes from 50,000
#'   individuals and 200 SNPs. Two of the SNPs have non-zero effects
#'   on the multivariate response. The response data are generated under
#'   a linear regression model. The simulated response and the columns
#'   of the genotype matrix are centered.
#'
#' @format \code{FinemappingConvergence} is a list with the following
#' elements:
#'
#' \describe{
#'
#'   \item{XtX}{Summary statistics computed using the centered and
#'     scaled genotype matrix.}
#'
#'   \item{Xty}{Summary statistics computed using the centered and
#'     scaled genotype data, and the centered simulated response.}
#'
#'   \item{yty}{yty is computed using the centered simulated response.}
#'
#'   \item{n}{The sample size (50,000).}
#'
#'   \item{true_coef}{The coefficients used to simulate the responses.}
#'
#'   \item{z}{z-scores from a simple (single-SNP) linear regression.}}
#'
#' @docType data
#'
#' @keywords data
#'
#' @seealso A similar data set with more SNPs is used in the
#'   \dQuote{Refine SuSiE model} vignette.
#'
#' @examples
#' data(FinemappingConvergence)
NULL

#' @name SummaryConsistency
#'
#' @title Simulated Fine-mapping Data with LD matrix From Reference Panel.
#'
#' @description Data simulated using real genotypes from 10,000
#'   individuals and 200 SNPs. One SNP have non-zero effect
#'   on the multivariate response. The response data are generated under
#'   a linear regression model. There is also one SNP with flipped allele
#'   between summary statistics and the reference panel.
#'
#' @format \code{SummaryConsistency} is a list with the following
#' elements:
#'
#' \describe{
#'
#'   \item{z}{z-scores computed by fitting univariate simple regression
#'     variable-by-variable.}
#'
#'   \item{ldref}{LD matrix estimated from the reference panel.}
#'
#'   \item{flip_id}{The index of the SNP with the flipped allele.}
#'
#'   \item{signal_id}{The index of the SNP with the non-zero effect.}}
#'
#' @seealso A similar data set with more samples is used in the
#'   \dQuote{Diagnostic for fine-mapping with summary statistics}
#'   vignette.
#'
#' @docType data
#'
#' @keywords data
#'
#' @examples
#' data(SummaryConsistency)
NULL
