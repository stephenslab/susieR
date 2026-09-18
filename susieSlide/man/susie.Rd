\name{susie}
\alias{susie}
\title{SuSiE with an empirical-Bayes genotype slider}
\description{Fits individual-level genotype effects using
\eqn{x+\delta I(x=1)}, with \eqn{-1\leq\delta\leq1} for each SNP and
single-effect component. Uses the IBSS engine from susieR.}
\usage{
susie(X, y, L = min(10, ncol(X)), scaled_prior_variance = 0.2,
  residual_variance = NULL, prior_weights = NULL, null_weight = 0,
  standardize = TRUE, intercept = TRUE,
  estimate_residual_variance = TRUE, estimate_prior_variance = TRUE,
  estimate_prior_method = "optim",
  min_obs = 5, delta = NULL, chunk_size = 1000L,
  cache_heterozygotes = FALSE, coverage = 0.95, min_abs_corr = 0.5,
  max_iter = 100, tol = 1e-3, verbose = FALSE, ...)
}
\arguments{
\item{X}{Sample-by-SNP numeric matrix or numeric sparse matrix of complete
hard calls 0, 1, 2. Supply original genotypes, not standardized predictors.}
\item{y}{Numeric phenotype vector.}
\item{L}{Maximum number of single-effect components.}
\item{scaled_prior_variance}{Gaussian coefficient prior variance divided by
the sample variance of y; scalar or length L, as in susieR.}
\item{residual_variance}{Residual variance, initially var(y) if NULL.}
\item{prior_weights}{Prior SNP weights, as in susieR.}
\item{null_weight}{Prior probability of the explicit no-effect column.}
\item{standardize}{Use the original additive genotype standard deviation
to scale both genotype and heterozygote terms.}
\item{intercept}{Include an unpenalized intercept by centering.}
\item{estimate_residual_variance}{Update Gaussian residual variance.}
\item{estimate_prior_variance}{Estimate the effect prior variance.}
\item{estimate_prior_method}{Prior-variance update: optim, EM, or simple.
Set estimate_prior_variance=FALSE to keep the variance fixed.
Each likelihood evaluation uses fitted slider values.
EM refreshes evidence and moments after its variance update.}
\item{min_obs}{Nonnegative integer. If any genotype count is smaller,
force delta=0 for that SNP in every component. Zero disables the rule.}
\item{delta}{NULL to estimate sliders, or fixed values: scalar, length-p
vector, or L-by-p matrix. The min_obs safeguard takes precedence.}
\item{chunk_size}{Maximum SNP block size for heterozygote calculations;
also capped internally to bound temporary block memory.}
\item{cache_heterozygotes}{Cache the full heterozygote matrix; otherwise
reconstruct it in blocks to avoid a second full genotype-sized matrix.}
\item{coverage}{Credible-set probability, or NULL to omit credible sets.}
\item{min_abs_corr}{Minimum absolute correlation for credible-set purity.
Correlations use each component's fitted transformed genotypes.}
\item{max_iter,tol,verbose}{IBSS iteration limit, ELBO convergence tolerance,
and progress output.}
\item{...}{Supported additional susieR options: check_null_threshold,
prior_tol, residual_variance_upperbound, residual_variance_lowerbound,
na.rm (missing y only), n_purity, median_abs_corr,
estimate_residual_method (MoM or MLE), track_fit, and model_init.
Warm starts require a susieSlide fit with matching dimensions and scales.
Other options are rejected rather than silently using additive calculations.}
}
\details{
Let h=I(x=1). With an intercept and standardization, the internal predictor is
\eqn{\{x-\bar{x}+\delta(h-\bar{h})\}/s_x}. The scale s_x does not change
with delta. The Gaussian coefficient prior is on this internal scale;
delta always retains its raw-genotype interpretation. Setting standardize
to FALSE gives a prior on the raw genotype coefficient.

Beta is integrated analytically. Delta maximizes the marginal likelihood:
the compiled solver compares all stationary points and endpoints on [-1,1].
The reported log-BFs are plug-in values conditional on fitted delta, not
evidence integrated over a delta prior. PIPs and credible sets likewise
condition on fitted sliders; frequentist calibration is not established by
the example simulations. Distinct components can select the same SNP and
have different sliders.

The first version supports individual-level hard-call genotypes and Gaussian
single-effect priors. RSS input, genotype missingness/dosages, NIG residual
priors, infinitesimal effects, refinement, greedy L, and slot priors are not
implemented. Ordinary additive summary statistics are insufficient.
}
\value{A list of class c("susie_slide","susie") with usual SuSiE fields
(alpha, mu, mu2, pip, sets, lbf, lbf_variable, V, sigma2, elbo, fitted), plus
delta (L-by-p, matching alpha), mu_delta, genotype_counts, delta_forced,
delta_cs, and expected_squared_residuals. An explicit null column, if used,
appears in component matrices but not SNP PIPs. mu and mu2 use internal
coefficient units. See coef and predict methods for original-scale effects.}
\examples{
set.seed(1)
X <- matrix(rbinom(4000, 2, 0.3), 200, 20)
y <- 1.5 * (X[, 1] == 2) + rnorm(200)
fit <- susie(X, y, L = 1, estimate_prior_variance = FALSE)
fit$delta_cs
head(predict(fit, newx = X))
}
