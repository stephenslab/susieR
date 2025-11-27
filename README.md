# susieR

[![CI](https://github.com/stephenslab/susieR/actions/workflows/ci.yml/badge.svg)](https://github.com/stephenslab/susieR/actions/workflows/ci.yml)
[![CRAN status badge](https://www.r-pkg.org/badges/version/susieR)](https://cran.r-project.org/package=susieR)
[![Codecov test coverage](https://codecov.io/gh/StatFunGen/susieR/graph/badge.svg)](https://app.codecov.io/gh/StatFunGen/susieR)

The `susieR` package implements a simple new way to perform variable
selection in multiple regression ($y=Xb+e$). The methods implemented
here are particularly well-suited to settings where some of the X
variables are highly correlated, and the true effects are highly
sparse (e.g. <20 non-zero effects in the vector $b$).  One example of
this is genetic fine-mapping applications, and this application was a
major motivation for developing these methods. However, the methods
should also be useful more generally.

The methods are based on a new model for sparse multiple regression,
which we call the "Sum of Single Effects" (SuSiE) model.  This model,
which is described in [Wang et al. (2020)](https://doi.org/10.1111/rssb.12388), lends itself to a particularly simple and intuitive fitting 
procedure -- effectively a Bayesian modification of simple forward 
selection, which we call "Iterative Bayesian Step-wise Selection".

The output of the fitting procedure is a number of "Credible Sets"
(CSs), which are each designed to have high probability to contain a
variable with non-zero effect, while at the same time being as small
as possible. You can think of the CSs as being a set of "highly
correlated" variables that are each associated with the response: you
can be confident that one of the variables has a non-zero coefficient,
but they are too correlated to be sure which one.

The package was initially developed by Gao Wang, Peter Carbonetto,
Yuxin Zou, Kaiqian Zhang, and Matthew Stephens from the
[Stephens Lab](https://stephenslab.uchicago.edu) at the University of
Chicago. It was later extended with new methods and implementations by
Alexander McCreight from the [StatFunGen Lab](https://wanggroup.org/) at
Columbia University.

Please
[post issues](https://github.com/stephenslab/susieR/issues) to ask
questions, get our support or provide us feedback; please
[send pull requests](https://github.com/stephenslab/susieR/pulls) if
you have helped fixing bugs or making improvements to the source code.

## Quick Start

Install susieR from [CRAN](https://cran.r-project.org/package=susieR):

```R
install.packages("susieR")
```

Alternatively, install the latest development version of `susieR`
from GitHub:

```R
# install.packages("remotes")
remotes::install_github("stephenslab/susieR")
```

See [here](https://stephenslab.github.io/susieR/articles/mwe.html) for
a brief illustration of `susieR`. For more documentation and examples
please visit https://stephenslab.github.io/susieR

## Citing this work

If you find the `susieR` package or any of the source code in this
repository useful for your work, please cite both:

> Wang, G., Sarkar, A., Carbonetto, P. & Stephens, M. (2020). A
> simple new approach to variable selection in regression, with
> application to genetic fine mapping. *Journal of the Royal
> Statistical Society, Series B* **82**, 1273–1300.
> https://doi.org/10.1111/rssb.12388

> McCreight, A., Cho, Y., Li, R., Nachun, D., Gan, H-Y., Carbonetto, P., Stephens,
> M., Denault, W.R.P. & Wang, G. (2025). SuSiE 2.0:
> improved methods and implementations for genetic fine-mapping and
> phenotype prediction. Submitting to *Genome Biology*.

If you use any of the summary data methods such as `susie_ss` or 
`susie_rss`, please also cite:

> Zou, Y., Carbonetto, P., Wang, G. & Stephens, M. (2022). Fine-mapping
> from summary data with the "Sum of Single Effects" model. *PLoS
> Genetics* **18**, e1010299. https://doi.org/10.1371/journal.pgen.1010299

If you use the Servin-Stephens prior on residual variance estimates
(`estimate_residual_method = "Servin_Stephens"`), please also cite:

> Denault, W.R.P., Carbonetto, P., Li, R., Alzheimer's Disease Functional
> Genomics Consortium, Wang, G. & Stephens, M. (2025). Accounting for
> uncertainty in residual variances improves calibration of the "Sum of
> Single Effects" model for small sample sizes. *bioRxiv*, 2025-05.
> Under review for *Nature Methods*.

If you use infinitesimal effects modeling (`unmappable_effects = "inf"`), 
please also cite:

> Cui, R., Elzur, R.A., Kanai, M. et al. (2024). Improving fine-mapping
> by modeling infinitesimal effects. *Nature Genetics* **56**, 162–169.
> https://doi.org/10.1038/s41588-023-01597-3

## Developer notes


+ The `Makefile` contains various R commands to build and maintain the package. 
For example to build the website via `pkgdown`:

   ```bash
   make pkgdown
   ```

+ When any changes are made to `roxygen2` markup, run
`make document` to update package `NAMESPACE` and documentation
files.


+ To format R codes in the `R` folder,

   ```bash
   for i in `ls R/*.R`; do bash inst/misc/format_r_code.sh $i; done
   ```

[susie-preprint]: https://doi.org/10.1101/501114
