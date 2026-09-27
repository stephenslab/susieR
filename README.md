# susieRSlidePrior

This repository's **root package is susieRSlidePrior** on branch `susie_slide_prior`.
It integrates a finite heterozygote-slider prior for each SNP and single-effect
component. The default grid has 17 values from -1 to 1, spaced by 0.125, with
fixed uniform probabilities. Probability learning belongs to a separate outer
workflow, as for mixed coding; this fitter never updates `delta_prior`.

The root `DESCRIPTION` says `Package: susieRSlidePrior`.

```r
# From C:/Document/Serieux/Travail/Package/git/susieR:
devtools::document(roclets = c("rd", "collate", "namespace"))
devtools::install(upgrade = "never")
```

Use `susieRSlidePrior::susie()` for the slider. The original additive entry point
is explicitly named `susieRSlidePrior::susie_additive()`. Other inherited interfaces
such as `susie_ss()` and `susie_rss()` remain additive and do not estimate
sliders. Both `susieR` and the original `susieSlide` can be installed alongside
`susieRSlidePrior` for comparison.

```r
fit <- susieRSlidePrior::susie(X, y, L = 10, min_obs = 5)

# Supply different fixed probabilities (zeros are allowed).
w <- rep(1, 17)
w[9] <- 8                       # put more prior mass on additive coding
fit <- susieRSlidePrior::susie(X, y, L = 10, delta_prior = w / sum(w))

fit$pip                         # SNP inclusion probabilities
fit$sets                        # familiar credible-set output
fit$delta                       # L-by-p matrix, aligned with fit$alpha
fit$delta_prior                 # supplied probabilities, fixed throughout fitting
fit$alpha_delta                 # L x p x 17 joint SNP-slider posterior probabilities
fit$delta_weights               # L x 17 posterior slider probabilities
fit$delta_prior_counts          # L x 17 counts for an external learning step
fit$delta_cs$summary             # one row per member SNP in each reported CS
fit$delta_cs$delta               # reported CSs x all input SNPs
susieRSlidePrior::slider_cs_table(fit)
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

Beta is integrated analytically and delta is averaged over its fixed discrete
prior. The compiled kernel reuses the genotype and heterozygote residual scores
for all grid values; it does not build an n-by-17p matrix. Supply `delta_grid`
and matching `delta_prior` for a different finite grid (increasing, bounded by
-1 and 1, including zero). `delta_prior=NULL` explicitly selects the original
continuous plug-in maximization, and `delta=...` overrides the prior with a
fixed scalar, SNP vector or component-by-SNP matrix.

If any of the counts of genotypes 0, 1, or 2 is below `min_obs`, delta is
fixed at zero, even when a different fixed delta was supplied. An absent
class counts as zero; exactly five observations passes the default rule.
This is an additive fallback, not a statistical test of additivity.

## Output and interpretation

- `delta[l,j]` is the posterior mean conditional on SNP j being selected in
  component l. It is a summary, not a plug-in coding used to fit the model.
  Its row corresponds to the same row of `alpha`, `mu`, and `mu2`.
- `sets$cs_index` maps reported credible sets to component rows.
- `delta_cs$summary` is the member-SNP table returned by `slider_cs_table()`.
  `delta_cs$delta` selects the corresponding rows of `delta` using
  `sets$cs_index`, with row names matching `sets$cs`. Its columns contain
  all input SNPs in their original order, including SNPs outside each CS,
  and exclude an explicit null column. If no CS is reported, it is a
  zero-row matrix with one column per input SNP. Deltas are estimated for
  each candidate SNP within each component, rather than shared across a CS.
- `mu` and `mu2` are first and second moments on the internal coefficient
  scale. `mu_delta` is E(beta*delta | SNP, component), which generally differs
  from `mu * delta`. `mu2_delta` and `mu2_delta2` retain the corresponding
  second moments for expected residual sums of squares.
- `alpha_delta[l,j,k]` is the joint posterior mass on SNP j and grid value k.
  Sum over k to recover `alpha`. `mu_grid` and `mu2_grid` store the conditional
  Gaussian beta moments for each pair. Together these permit exact posterior
  sampling within each variational component.
- `delta_weights` sums joint mass over SNPs. `delta_prior_counts` excludes
  count-forced SNPs, explicit null columns and zero-variance components. Its
  rows retain the original components so an outer mixed-coding-style workflow
  can apply its eligibility rules; no purity filter or M-step is run here.
- `lbf_variable` integrates over the fixed slider prior; `lbf` also integrates
  over the SNP prior. Zero prior probabilities receive exactly zero posterior
  mass. Count-forced additive SNPs have a separate point prior at zero, even
  when the supplied grid prior has zero additive mass.
- `coef(fit)` returns two coefficient columns because one additive vector
  cannot represent a slider prediction. On the original scale,
  `prediction = intercept + X %*% b + I(X == 1) %*% b_heterozygote`.
- Fitted values, residual updates, expected squared residuals, residual
  variance updates, and the conditional ELBO all include both terms.
- Credible-set membership uses SNP probabilities marginalized over delta.
  Purity uses each component's posterior-mean slider coding as a representative
  coding diagnostic; it does not average correlations over the posterior.
  Calling `susieR::susie_get_cs(fit, X=X)` manually would instead
  apply additive-genotype purity; use the returned `fit$sets`.
- An explicit null column is included in component matrices when requested,
  but excluded from SNP PIPs and coefficient rows.
- Multiple components may select the same SNP with different deltas. The
  resulting summed effect need not itself have one bounded slider.

PIPs and effect moments integrate over delta conditional on the supplied
probabilities. Calibration remains an empirical property to assess. The
optional legacy plug-in path instead conditions on optimized deltas.

## Supported interface

Supports dense and numeric sparse genotype matrices, fixed or estimated
Gaussian residual/prior variances, scalar/vector/component-specific fixed
deltas, SNP prior weights, an optional null column, count filtering,
predictions, summaries, and matching-grid/dimension/scaling warm starts. New
`prior_weights` and `delta_prior` inputs are respected when warm-starting an
outer coding-prior loop. The `estimate_prior_method="EM"` option updates only
the Gaussian effect-prior variance and refreshes the posterior at that variance.
It does not estimate slider probabilities.

The slider entry point does not implement ordinary RSS/summary-statistic input, dosage
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

The included engine derives from this repository's `susieR` 0.16.6 at commit
`8e56a8e038e989856d106d9ca5175cc664fea9d2`. Slider methods dispatch within this
package. Gaussian priors and the inherited additive engine are retained.

From the repository root, using Rtools on Windows:

```sh
R CMD INSTALL .
R CMD build .
R CMD check --no-manual susieRSlidePrior_0.3.0.tar.gz
```

From R in the package source directory:

```r
devtools::test()
```

Documentation and NAMESPACE are generated from the R sources. Native-code
registration is generated with `cpp11::cpp_register()` (also run by devtools
when compiling). If roxygen2 reports a missing `decor` dependency, install it
with `install.packages("decor")`; this is a development-tool dependency.
Do not rename only DESCRIPTION: native registration and documentation must
use the same package name. They are already configured for susieRSlidePrior here.

## Reproduce comparisons

From the package source directory:

```sh
Rscript inst/examples/compare_models.R
Rscript inst/examples/validation_diagnostics.R
Rscript inst/examples/plot_results.R
Rscript inst/examples/benchmark_fit.R
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
`system.file("examples", package="susieRSlidePrior")`.

The external susieR package is suggested only for the comparison scripts,
which deliberately compare against a separately installed additive package.
The current package vignette is `vignettes/slider-model.Rmd`. Inherited
additive tutorials are preserved in `validation/upstream-vignettes/` for
reference; they are not built as tutorials for the slider model.
