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

#' @name rss_mismatch_cases
#'
#' @title Simulated Summary Statistics with LD Reference Mismatch.
#'
#' @description Five simulated regions, each with one causal variant and
#'   121 variants, used to illustrate how \code{susie_rss} handles a
#'   mismatch between the summary statistics and the LD reference. The
#'   first four cases simulate a phenotype from GTEx genotypes (n = 838);
#'   the last is a meta-analysis of an African-ancestry and a
#'   European-ancestry cohort.
#'
#' @format \code{rss_mismatch_cases} is a list of five cases. Each case is
#' a list with the inputs to \code{susie_rss}:
#'
#' \describe{
#'
#'   \item{z}{z-scores.}
#'
#'   \item{R}{LD matrix; for \code{multi_panel}, a named list of LD
#'     matrices, one per reference panel.}
#'
#'   \item{X}{For \code{ld_mismatch_sketch} only, a genotype sketch
#'     (factor matrix) of the reference panel, used instead of \code{R}.
#'     This case also stores \code{R}, a copy of \code{R_insample}.}
#'
#'   \item{n}{Sample size of the summary statistics.}
#'
#'   \item{B}{Reference panel sample size, or the sketch dimension
#'     (number of rows of \code{X}), passed as \code{R_finite}; one value
#'     per panel for \code{multi_panel}. Absent for \code{allele_flip},
#'     whose LD is in-sample.}
#'
#'   \item{causal}{Index of the true causal variant.}}
#'
#' The cases are:
#'
#' \describe{
#'
#'   \item{allele_flip}{In-sample LD, with the sign of the z-score
#'     reversed at 1\% of the variants. Also contains \code{z_clean}, the
#'     z-scores before flipping, and \code{flip_id}, the indices of the
#'     flipped variants.}
#'
#'   \item{ld_mismatch}{LD from 489 European individuals in 1000 Genomes.}
#'
#'   \item{ld_mismatch_sketch}{The reference panel given as a 500-row
#'     genotype sketch \code{X} of the ADSP panel.}
#'
#'   \item{ser_fallback}{LD from 489 European individuals in 1000 Genomes,
#'     a region where \code{SER_fallback} is \code{TRUE}.}
#'
#'   \item{multi_panel}{Meta-analysis of 400 African-ancestry and 400
#'     European-ancestry samples, with AFR and EUR panel LD in \code{R}.
#'     Also contains \code{prop}, the true proportion of samples from each
#'     ancestry.}}
#'
#' \code{ld_mismatch}, \code{ld_mismatch_sketch} and \code{ser_fallback}
#' also contain \code{z_insample} and \code{R_insample}, the z-scores and
#' the in-sample LD on the same variants.
#'
#' @docType data
#'
#' @keywords data
#'
#' @seealso The \dQuote{Modeling and Accounting for LD Reference Mismatch
#'   in Summary Statistics Fine-mapping} vignette.
#'
#' @examples
#' data(rss_mismatch_cases)
#' names(rss_mismatch_cases)
NULL

#' @name data_small
#'
#' @title Simulated Small-sample eQTL Data.
#'
#' @description A simulated eQTL data set with 47 individuals and 7,430
#'   variables. The response is a simulated gene expression phenotype and
#'   the variables are genotypes. This data set illustrates the small
#'   sample-size setting considered in Denault \emph{et al} (2025).
#'
#' @format \code{data_small} is a list with the following elements:
#'
#' \describe{
#'
#'   \item{y}{Simulated gene expression response.}
#'
#'   \item{X}{Genotype matrix.}}
#'
#' @docType data
#'
#' @keywords data
#'
#' @seealso The \dQuote{Small data example} vignette.
#'
#' @references
#' W. R. P. Denault \emph{et al} (2025). Accounting for uncertainty in
#'   residual variances improves fine-mapping in small sample studies.
#'   \emph{bioRxiv} \doi{10.1101/2025.05.16.654543}.
#'
#' @examples
#' data(data_small)
NULL

#' @name unmappable_data
#'
#' @title Simulated Fine-mapping Data with Sparse, Oligogenic and Polygenic Effects.
#'
#' @description A simulated data set with 1,000 individuals and 5,000
#'   variants, combining 3 sparse, 5 oligogenic and 15 polygenic
#'   non-zero effects. The response is generated under a linear
#'   regression model. This data set illustrates fine-mapping with
#'   SuSiE-ash and SuSiE-inf.
#'
#' @format \code{unmappable_data} is a list with the following elements:
#'
#' \describe{
#'
#'   \item{X}{Centered and scaled genotype matrix.}
#'
#'   \item{y}{Simulated response.}
#'
#'   \item{beta}{Simulated effect sizes.}
#'
#'   \item{h2g}{Total proportion of variance in the response explained
#'     by the simulated effects.}}
#'
#' @docType data
#'
#' @keywords data
#'
#' @seealso The \dQuote{Fine-mapping with SuSiE-ash and SuSiE-inf}
#'   vignette.
#'
#' @examples
#' data(unmappable_data)
NULL
