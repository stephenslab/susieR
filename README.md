# susieR

[![Travis Build Status](https://travis-ci.org/stephenslab/susieR.svg?branch=master)](https://travis-ci.org/stephenslab/susieR)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/tcgeqxd8q8krija6?svg=true)](https://ci.appveyor.com/project/pcarbo/susier)
[![codecov](https://codecov.io/gh/stephenslab/susieR/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/susieR)

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
[Stephens Lab](http://stephenslab.uchicago.edu/) at the University of
Chicago.

This is very much work in progress. Please
[post issues](https://github.com/stephenslab/susieR/issues) to ask
questions, get our support or provide us feedback; please
[send pull requests](https://github.com/stephenslab/susieR/pulls) if
you have helped fixing bugs or making improvements to the source code.

## Setup

To automatically retrieve and install `susieR` from this repository,
with the vignettes included, run:

```R
devtools::install_github("stephenslab/susieR@0.8.0",build_vignettes = TRUE)
```

## Quick Start

See [here](https://stephenslab.github.io/susieR/articles/mwe.html) for
a brief illustration of `susieR`.  For more documentation and examples
please visit https://stephenslab.github.io/susieR

## Developer notes

+ When any changes are made to `roxygen2` markup, simply run 
`devtools::document()` to update package `NAMESPACE`
and documentation files.

+ To install and test the development version of susieR, run the
following commands in the shell:

    ```bash
    R CMD build --resave-data --no-build-vignettes susieR
    R CMD INSTALL susieR_0.8.1.0525.tar.gz
    R CMD check --as-cran --ignore-vignettes susieR_0.8.1.0525.tar.gz
    ```

+ Run `pkgdown::build_site()` to build the website. Getting `pkgdown`
to work properly can be frustrating due to numerous & fragile dependencies. 
If `pkgdown` does not work for you out of the box you can use this `docker`
command to run all vignettes and build the site:

```bash
docker run --rm --security-opt label:disable -t -P -w $PWD -v $PWD:$PWD \
  -u $UID:${GROUPS[0]} -e HOME=/home/$USER -e USER=$USER gaow/susie \
  R --slave -e "pkgdown::build_site(lazy=TRUE, examples=FALSE)"
```
