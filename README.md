# susieR

[![CRAN status badge](https://www.r-pkg.org/badges/version/susieR)](https://cran.r-project.org/package=susieR)
[![Travis Build Status](https://travis-ci.org/stephenslab/susieR.svg?branch=master)](https://travis-ci.org/stephenslab/susieR)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/tcgeqxd8q8krija6?svg=true)](https://ci.appveyor.com/project/pcarbo/susier)
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/stephenslab/susieR/tree/master.svg?style=svg)](https://app.circleci.com/pipelines/github/stephenslab/susieR?branch=master)
[![codecov](https://codecov.io/gh/stephenslab/susieR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/stephenslab/susieR)

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
which will be described in a manuscript in preparation (Wang et al),
lends itself to a particularly simple and intuitive fitting procedure
-- effectively a Bayesian modification of simple forward selection,
which we call "Iterative Bayesian Step-wise Selection".

The output of the fitting procedure is a number of "Credible Sets"
(CSs), which are each designed to have high probability to contain a
variable with non-zero effect, while at the same time being as small
as possible. You can think of the CSs as being a set of "highly
correlated" variables that are each associated with the response: you
can be confident that one of the variables has a non-zero coefficient,
but they are too correlated to be sure which one.

The package is developed by Gao Wang, Peter Carbonetto, Yuxin Zou,
Kaiqian Zhang, and Matthew Stephens from the
[Stephens Lab](https://stephenslab.uchicago.edu) at the University of
Chicago.

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
repository useful for your work, please cite:

> G. Wang, G., Sarkar, A., Carbonetto, P. & Stephens, M. (2020). A
> simple new approach to variable selection in regression, with
> application to genetic fine mapping. *Journal of the Royal
> Statistical Society, Series B* **82**, 1273–1300.
> https://doi.org/10.1111/rssb.12388

If you use any of the summary data methods such as `susie_suff_stat`
or `susie_rss`, please also cite:

> Zou, Y., Carbonetto, P., Wang, G. & Stephens, M. (2022). Fine-mapping
> from summary data with the "Sum of Single Effects" model. *PLoS
> Genetics* **18**, e1010299. https://doi.org/10.1371/journal.pgen.1010299

## Developer notes

+ When any changes are made to `roxygen2` markup, run
`devtools::document()` to update package `NAMESPACE` and documentation
files.

+ Run `pkgdown::build_site()` to build the website. Getting `pkgdown`
to work properly can be frustrating due to numerous & fragile dependencies.
If `pkgdown` does not work for you out of the box you can use this `docker`
command to run all vignettes and build the site:

```bash
docker run --rm --security-opt label:disable -t -P -w $PWD -v $PWD:$PWD \
  -u $UID:${GROUPS[0]} -e HOME=/home/$USER -e USER=$USER gaow/susie \
  R --slave -e "pkgdown::build_site(lazy=TRUE, examples=FALSE)"
```

[susie-preprint]: https://doi.org/10.1101/501114
