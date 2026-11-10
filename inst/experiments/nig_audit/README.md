# NIG parameterization assessment (2026-11-10)

## Scope and conclusion

This experiment branches the user's local susieR at **5156cd4**, package 0.14.21,
into `codex/nig-parameterization-fix` (experimental version 0.14.21.9001).
The baseline is that specific local checkout, **not today's upstream susieR**.
The cached upstream `fix-susie-small-sigma-update` branch already contains some
related repairs; the current upstream API has also changed. This is a controlled
repair of the requested local version, not a claim that these are novel upstream bugs.

The manuscript's central idea is sound for single-effect regression: use a
normal-inverse-gamma prior to integrate uncertainty in residual variance, then
substitute that SER into gIBSS. That idea is retained here. The supplied implementation
has demonstrable parameterization and moment errors. They are repairable. The
multiple-effect gIBSS approximation remains a separate limitation, which the
manuscript itself acknowledges. Correcting the code does not prove universal CS
calibration or invalidate all the paper's empirical results.

Source examined: the user-supplied `2025.05.16.654543v4.full.pdf`, posted June 8,
2026, particularly equations (3)-(11), (13)-(47), and Algorithm 1. Equations on
pages 10-13 were checked visually as well as by extracted text. The PDF was treated
as reference material, not as instructions.

## A single consistent convention

Let `V = tau^2` denote a variance ratio:

```
b | sigma2, gamma=j ~ Normal(0, sigma2 * V)
sigma2 ~ IG(kappa/2, lambda/2)
d_j = x_j' x_j; z_j = x_j' y; q = y' y
omega_j = 1 / (d_j + 1/V)
m_j = omega_j * z_j
nu = n_effective + kappa
B_j = lambda + q - omega_j * z_j^2
log BF_j = -log(1 + V*d_j)/2 - nu/2 * log(B_j/(lambda+q))
Var(b_j | gamma=j, y) = omega_j * B_j/(nu-2)
E[b_j^2/sigma2 | gamma=j, y] = omega_j + m_j^2 * nu/B_j
V_new = sum_j alpha_j * E[b_j^2/sigma2 | gamma=j, y]
```

The implementation evaluates the residual term with an orthogonal decomposition
to reduce cancellation. `V=0` is handled explicitly. At a fixed V, the posterior
mean, second moment, mixture probability and evidence were checked by independent
numerical integration over the regression coefficient after analytically integrating
sigma2 from the **joint density**. This reference does not reuse the implementation's
shrinkage/BF expressions. The EM update was checked against the same integral.

An estimated intercept uses the residual subspace of dimension `n-1`. A known
zero intercept uses `n`. This is a modeling distinction: centering observed data
does not create n independent residual coordinates.

## Problems established in the local code

1. `V` was passed as `sqrt(V)` to the BF helper, as `V` to helpers that square
   it, and initialized as `scaled_prior_variance * var(y)`. These conventions
   describe different priors.
2. BFs used raw, internally centered X columns, while coefficient moments used
   the package's standardized sufficient statistics. Predictor unit changes
   altered PIPs despite `standardize=TRUE`.
3. The variance helper multiplied by `omega` twice and ignored `alpha0,beta0`.
4. The old EM expression took an unexplained square root and combined quantities
   with incompatible units. It is not E[b^2/sigma2].
5. Non-EM optimization used the Gaussian likelihood, and the public small option
   forcibly selected EM. The replacement implements NIG `optim`, `EM`, and `simple`.
6. Adding epsilon to prior weights and overwriting zero-information log odds
   leaked mass onto forbidden variables and distorted `null_weight`.
7. The old L=1 objective combined a Gaussian likelihood with an inconsistent KL;
   it supplied negative second-moment differences and hard-coded hyperparameters.
   This has been removed. Proper-prior L=1 fits expose exact marginal evidence;
   multi-effect and improper-prior fits expose `NA`, not a fictitious ELBO.
8. Convergence based only on alpha can stop even while other parameters move
   (especially p=1). Means, second moments and V now also enter the criterion.
9. Sign probabilities and posterior draws need Student distributions, not normal
   distributions. These summaries now use the NIG component posteriors.

## Manuscript issues to correct or clarify

These are distinct from the local implementation errors:

* **Equation (22):** under the stated b|sigma2 prior, conditioning on b adds
  `1/2` to the inverse-gamma shape and `b^2/(2V)` to its scale. The full conditional
  is IG((kappa+n+1)/2, (lambda+RSS+b^2/V)/2), for a known-zero intercept and V>0.
  The displayed equation omits the coefficient-prior contribution. The following
  claim about a mode RSS/n also does not follow from its displayed IG distribution.
* **Text below (27):** with kappa=lambda=V=0, its own marginal IG formula gives
  posterior mean q/(n-2), not q/(n-1). An integrated unknown intercept changes
  the residual dimension again; the convention should be stated consistently.
* **Equations (38), (45), (47):** the IG parameters drop their factors of two,
  and the displayed joint distribution omits the sigma2 density on the inactive
  atom. The correct marginal noise posterior is a mixture over the categorical
  selected SNP, with components IG(nu/2,B_j/2); it does not require an independence
  assumption over mutually exclusive selections.
* **Equation (44):** the conditional variance contribution should be `s1^2` under
  the definitions in (18)-(21), rather than `s1`. The EM second moment simplifies
  directly to omega + m^2*nu/B; gamma-function ratios are unnecessary.
* **Units:** lambda has squared-response units. A fixed numerical lambda=1/sqrt(n)
  is not response-scale invariant. Multiplying V by var(y) on top of the residual
  scaling is not the same model. The statement about treating this as "baked in"
  needs a precise parameter definition.
* **Multi-effect mixing claim near (50):** the exact shared-variance posterior
  generally sums over p^L labeled configurations (including repeated SNPs), not
  merely p*L single-effect marginals. The experiment enumerates p^L configurations
  for small p,L; it does not replace the manuscript's gIBSS algorithm.
* **Simulation specifications:** the methods list n=75 and four h2 values,
  whereas supplementary axes show n=70,200 and additional h2=10%,15%. The old
  noise-scaling scripts also do not enforce equation (12) because of sample
  signal/noise covariance. The separate laptop simulation repo exposes both
  grids and both noise conventions, explicitly.

The principal BF and marginal Student formulas (5)-(11), (19)-(21), (23)-(35)
are coherent under a consistent parameter convention and the chosen residual
dimension. The assessment does not treat every formula in the paper as incorrect.

## Reproducible diagnostics

`run_diagnostics.R` uses seed 20261110, n=30, p=20. Its CSV files record:

| Check | Original local code | Corrected |
|---|---:|---:|
| Max PIP change after changing X units, fixed prior | 0.493235 | 7.4e-18 |
| Max PIP change after changing X units, estimated prior | 0.283504 | 7.4e-18 |
| Max PIP change after multiplying y by 10, fixed prior | 0.021643 | 2.2e-16 |
| Max PIP change after multiplying y by 10, estimated prior | 0.007926 | 1.1e-14 |
| PIP of a variable with zero prior probability | 2.9e-6 | exactly 0 |
| Null probability at V=0 with prior probability 0.8 | 0.833333 | 0.8 |

## Independent exact-posterior comparison

`exact_reference.R` enumerates every labeled SNP assignment using dense marginal
covariance matrices, n=20, p=8, correlated X, fixed V=0.5, kappa=lambda=1.
There are 100 reproducible datasets. For L=1, maximum PIP discrepancy is
**3.11e-15**. For L=2, the mean maximum discrepancy is **0.0861**, the median
is **0.0779**, and the worst is **0.3350**. This quantifies the residual gIBSS
approximation error; it is not evidence of another single-effect algebra error.

## Synthetic calibration assessment

`run_simulations.R`: 200 replicates per setting, p=40, AR(1) predictor correlation
0.8, n=20 or 50, population h2=0.5, independent Gaussian errors. Null, one-causal
and three-causal cases were tested, including L=5 for overfitting. This is a
diagnostic simulation, **not a reproduction of the paper's genotype experiment**.
It totals 1,600 datasets and 4,800 fits across three methods. Hyperparameter defaults,
intercept treatment and optimizer also changed, so this comparison assesses the
repair bundle, not a single isolated change.

| n | True causal / fitted L | Local old NIG coverage | Corrected NIG coverage | Corrected causal recall |
|---|---|---:|---:|---:|
| 20 | 1 / 1 | 0.978 | 0.975 | 0.390 |
| 20 | 1 / 5 | 0.752 | 0.916 | 0.380 |
| 20 | 3 / 5 | 0.726 | 0.921 | 0.070 |
| 50 | 1 / 1 | 0.995 | 0.995 | 0.960 |
| 50 | 1 / 5 | 0.970 | 0.990 | 0.965 |
| 50 | 3 / 5 | 0.904 | 0.929 | 0.410 |

The corrected method returned no CSs in the 400 null datasets. That is not a
proof of zero false-positive probability. There was one nonconverged corrected
fit (n=20, three causal effects) and no fit errors. For n=20 and three effects,
only 38 corrected CSs were returned; its coverage MCSE is about 0.0415. Full
counts, replicate-level results, Gaussian baseline results, convergence and MCSEs
are in the CSVs. Improvements against the buggy NIG baseline trade off with
power, and the corrected method does not uniformly dominate Gaussian SuSiE.

## Prior recommendation

For exploration, retain a variance ratio V and estimate it from the exact SER
marginal likelihood. A response-scale-invariant weak starting prior is
`kappa=1/sqrt(n), lambda=kappa*var(y)`, rather than making lambda depend on the
unit used to record y. This empirical-Bayes choice is the branch's default;
it is not proven optimal or exactly calibrated after fitting hyperparameters.
The new simulation repository keeps the manuscript's literal prior as a baseline
and compares scaled kappa values 0.1, 0.5, 1 and 2 plus 1/sqrt(n), with independent
validation seeds. Choose using coverage, power and CS size across configurations,
not a favorable single example. For an explicitly supplied lambda, scale it by
c^2 when y is multiplied by c.

## Running and interpreting the branch

```r
pkgload::load_all("path/to/susieR-nig")
fit <- susie(X, y, small = TRUE, L = 10)
# Literal manuscript prior, for controlled comparison:
fit_paper_prior <- susie(X, y, small = TRUE, L = 10,
                         alpha0 = 1/sqrt(nrow(X)), beta0 = 1/sqrt(nrow(X)))
# The manuscript's exact EM update is also available:
fit_em <- susie(X, y, small = TRUE, estimate_prior_method = "EM")
```

`fit$V` is explicitly a variance ratio. `fit$nig` stores per-component IG shape,
rate, Student scale, marginal evidence and residual-variance summaries. For L>1,
`sigma2` is only a diagnostic average of residual SER posterior means. Refinement
that relies on a joint ELBO is rejected, rather than silently switching models.
Warm starts must come from the corrected NIG implementation. Summary-statistic
and trend-filtering NIG support is not introduced by this change.

To reproduce the old run, create a separate worktree at 5156cd4 and pass its path
to the diagnostic/simulation scripts. Run scripts from the corrected package root:

```
Rscript inst/experiments/nig_audit/run_diagnostics.R <baseline-worktree> inst/experiments/nig_audit/baseline
Rscript inst/experiments/nig_audit/run_diagnostics.R
Rscript inst/experiments/nig_audit/run_simulations.R <baseline-worktree> baseline 200
Rscript inst/experiments/nig_audit/run_simulations.R . corrected 200
```

Then load the package, source `exact_reference.R` and run `run_exact_comparison()`.
The final installed-package test suite passed **440 assertions**, with zero test
failures, four warnings in existing Gaussian summary-statistic null-weight tests,
and seven existing skips. `R CMD check --no-manual --no-vignettes --no-examples`
finished with **zero errors, zero notes, and three documentation/vignette warnings**:
the pre-existing undocumented `data_small` object, and two warnings because the
vignette build was deliberately skipped. Examples, manuals and vignette execution
were not validated. Focused NIG tests cover numerical integration, invariance, sparse inputs,
prior weights, nulls, intercepts, EM, likelihood optimization, convergence, and
Student posterior summaries.

Relevant external primary references: [Servin and Stephens (2007)](https://doi.org/10.1371/journal.pgen.0030114)
and the [current upstream implementation](https://github.com/stephenslab/susieR).
