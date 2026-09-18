# susieSlide

An individual-level SuSiE extension with one heterozygote slider per SNP
and single-effect component. This is a separate package in the `susie_slide`
branch; it reuses `susieR`'s IBSS engine through registered extension methods.
The upstream package is not renamed or patched at runtime.

```r
fit <- susieSlide::susie(X, y, L = 10, min_obs = 5)

fit$pip                         # SNP inclusion probabilities
fit$sets                        # familiar credible-set output
fit$delta                       # L-by-p matrix, aligned with fit$alpha
fit$delta_cs                    # one row per SNP in each reported CS
susieSlide::slider_cs_table(fit)
predict(fit, newx = X_test)      # original 0/1/2 genotype matrix
coef(fit)                       # additive and heterozygote coefficient columns
```

The input must contain original complete hard-call genotypes. Create the
heterozygote indicator before any centering or other transformation. Missing
phenotypes can be removed using `na.rm=TRUE`; genotype imputation is not
performed by this package.

## Model and scaling

For component l and candidate SNP j, the internal predictor is

```
z_j(delta_lj) = ((x_j - mean(x_j)) +
                delta_lj * (I(x_j == 1) - mean(I(x_j == 1)))) / sd(x_j)
```

The original additive genotype SD is fixed across delta values. Both terms
use the same divisor. Delta remains in raw-genotype units: -1 is recessive,
0 additive, and +1 dominant. It does not need to be back-transformed.
`standardize=FALSE` sets the divisor to one; `intercept=FALSE` omits centering.
Constant-column scale factors are set to one, as in SuSiE.

The effect coefficient has the same Gaussian prior convention as SuSiE:
initial variance `scaled_prior_variance * var(y)` on the internal coefficient
scale, optionally learned separately for each effect. Thus a raw-scale
coefficient has variance `V_l / sd(x_j)^2` under standardization. Keeping a
fixed raw-scale prior instead is a different prior specification.

Beta is integrated analytically. Delta is estimated by maximizing the
marginal likelihood in [-1,1]. A compiled solver enumerates every root of
the cubic derivative in that interval, plus the boundaries. It brackets
roots using the quadratic derivative's turning points, avoiding a
unimodality assumption and unstable closed-form cubic divisions.

If any of the counts of genotypes 0, 1, or 2 is below `min_obs`, delta is
fixed at zero, even when a different fixed delta was supplied. An absent
class counts as zero; exactly five observations passes the default rule.
This is an additive fallback, not a statistical test of additivity.

## Output and interpretation

- `delta[l,j]` is conditional on SNP j being the selected SNP for component l.
  Its row corresponds to the same row of `alpha`, `mu`, and `mu2`.
- `sets$cs_index` maps reported credible sets to component rows.
- `mu` and `mu2` are first and second moments on the internal coefficient
  scale. `mu_delta` is `mu * delta`, valid for this plug-in delta model.
- `lbf_variable` contains fixed-delta Gaussian log-BFs evaluated at fitted
  delta, and `lbf` is their SNP-prior-weighted component log evidence score.
  **These are not integrated over a delta prior.**
- `coef(fit)` returns two coefficient columns because one additive vector
  cannot represent a slider prediction. On the original scale,
  `prediction = intercept + X %*% b + I(X == 1) %*% b_heterozygote`.
- Fitted values, residual updates, expected squared residuals, residual
  variance updates, and the conditional ELBO all include both terms.
- Credible-set purity uses each component's fitted transformed genotype
  columns. Calling `susieR::susie_get_cs(fit, X=X)` manually would instead
  apply additive-genotype purity; use the returned `fit$sets`.
- An explicit null column is included in component matrices when requested,
  but excluded from SNP PIPs and coefficient rows.
- Multiple components may select the same SNP with different deltas. The
  resulting summed effect need not itself have one bounded slider.

PIPs, effect moments, and credible sets condition on estimated deltas.
Optimizing many deltas can overfit null data. Example performance does not
establish genome-wide PIP or credible-set calibration; inspect the null
diagnostic alongside the effect-recovery examples.

## Supported interface

Supports dense and numeric sparse genotype matrices, fixed or estimated
Gaussian residual/prior variances, scalar/vector/component-specific fixed
deltas, SNP prior weights, an optional null column, count filtering,
predictions, summaries, and same-dimension/scaling warm starts. The EM option
updates the Gaussian effect-prior variance and refreshes the slider posterior
at the new variance before updating fitted values.

This release does not implement ordinary RSS/summary-statistic input, dosage
or missing-genotype handling, covariate input, NIG priors, infinitesimal/ash
effects, slot priors, greedy-L expansion, or refinement. Unsupported options
raise explicit errors. Covariate-adjusted extensions must construct the
heterozygote indicator before adjustment and adjust both basis matrices.

Genotype geometry and counts are cached. Residual products are recomputed for
each single-effect update. The default reconstructs heterozygotes in memory-
bounded blocks; `cache_heterozygotes=TRUE` trades a second genotype-sized
matrix for faster repeated products. No dense p-by-p LD matrix is created.
The package still takes an in-memory genotype matrix: packed-file streaming
and a full million-variant application are not implemented or benchmarked.

## Build and test

Tested against this repository's `susieR` 0.16.6 at commit
`8e56a8e038e989856d106d9ca5175cc664fea9d2`. Internal extension hooks are used;
the package checks required hooks at load time. Compatibility with arbitrary
future susieR versions is not implied.

From the repository root, using Rtools on Windows:

```sh
R CMD INSTALL .
R CMD INSTALL susieSlide
R CMD build susieSlide
R CMD check --no-manual susieSlide_0.1.0.tar.gz
```

From R in the package source directory:

```r
library(susieSlide)
testthat::test_dir("tests/testthat")
```

The DESCRIPTION maintainer email is an explicit placeholder
(`maintainer@example.org`) for local development; replace it before public
distribution.

## Reproduce comparisons

From the package source directory:

```sh
Rscript inst/examples/compare_models.R
Rscript inst/examples/validation_diagnostics.R
```

`compare_models.R` runs six inheritance scenarios with correlated hard-call
genotypes: additive, recessive, partially recessive, dominant, partially
dominant, and mixed. Ten independent replicates use n=1,000 training and
2,000 test observations, 80 SNPs, and two causal variants. It compares additive
SuSiE, the slider, equal-weight stacked coding, and an additional coding-weight
variational-EM comparator. The main comparison matches effect-prior scales and
uses L=2; an L=4 sensitivity fit permits stacking to use multiple components
per biological SNP.

`validation_diagnostics.R` additionally reproduces the core settings inspected
in the existing susie_mix workhorse: raw 0/1/2, 0/1/1, 0/0/1 columns;
`standardize=FALSE`, `estimate_prior_method="EM"`, L=10, learned residual
variance, and uniform predictor weights. Here EM estimates coefficient-prior
variances, not coding-class weights. Its settings differ from the matched-
prior experiment and are reported separately. This script also runs null
examples and an actual one-million-summary compiled inference benchmark.

All raw results and the interpretation are in `validation/` in the source
checkout. These files are excluded from the installed package. Source scripts
are also available after installation via
`system.file("examples", package="susieSlide")`.
