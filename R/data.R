#' @name N2finemapping
#'
#' @title Simulated fine-mapping data with 2 effect variables.
#'
#' @docType data
#' 
#' @description The data-set contains a matrix of about 600 individuals and 1,000 variables. 
#' These variables are real-world genotypes centered and scaled, and therefore retains the correlation
#' structure of variables in the original genotype data. 2 out of the variables have non-zero effects.
#' The response data is generated under a multivariate linear regression model. 
#' See Wang et al (2020) for more details.
#' 
#' @format \code{N2finemapping} is a list with the following elements:
#' 
#' \describe{
#'   \item{X}{N by P variable matrix of centered and scaled genotype data.}
#' 
#'   \item{chrom}{Chromomsome of the original data, in hg38 coordinate.}
#' 
#'   \item{pos}{Chromomosomal positoin of the original data, in hg38 coordinate. The information can be 
#'      used to compare impact of using other genotype references of the same variables in susie_rss application.}
#'
#'   \item{true_coef}{The simulated effect sizes.}
#'
#'   \item{residual_variance}{The simulated residual covariance matrix.}
#'
#'   \item{Y}{The simulated response variables.}
#'
#'   \item{allele_freq}{Allele frequency of the original genotype data.}
#'
#'   \item{V}{Prior covariance matrix for effect size of the two non-zero effect variables.}
#' }
#' 
#' @keywords data
#'
#' @references Wang, G., Sarkar, A., Carbonetto, P., & Stephens, M. (2020). 
#' A simple new approach to variable selection in regression, 
#' with application to genetic fine mapping. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology).
NULL

#' @name N3finemapping
#'
#' @title Simulated fine-mapping data with 3 effect variables.
#'
#' @docType data
#'  
#' @description The data-set contains a matrix of about 600 individuals and 1,000 variables. 
#' These variables are real-world genotypes centered and scaled, and therefore retains the correlation
#' structure of variables in the original genotype data. 3 out of the variables have non-zero effects.
#' The response data is generated under a multivariate linear regression model. 
#' See Wang et al (2020) for more details.
#' 
#' @format \code{N3finemapping} is a list with the following elements:
#' 
#' \describe{
#'   \item{X}{N by P variable matrix of centered and scaled genotype data.}
#' 
#'   \item{chrom}{Chromomsome of the original data, in hg38 coordinate.}
#' 
#'   \item{pos}{Chromomosomal positoin of the original data, in hg38 coordinate. The information can be 
#'      used to compare impact of using other genotype references of the same variables in susie_rss application.}
#'
#'   \item{true_coef}{The simulated effect sizes.}
#'
#'   \item{residual_variance}{The simulated residual covariance matrix.}
#'
#'   \item{Y}{The simulated response variables.}
#'
#'   \item{allele_freq}{Allele frequency of the original genotype data.}
#'
#'   \item{V}{Prior covariance matrix for effect size of the three non-zero effect variables.}
#' }
#' 
#' @keywords data
#' @references Wang, G., Sarkar, A., Carbonetto, P., & Stephens, M. (2020). 
#' A simple new approach to variable selection in regression, 
#' with application to genetic fine mapping. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology).
NULL
