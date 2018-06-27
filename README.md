# `susieR`

An R package for "sum of single effects" (SuSiE) regression.

SuSiE's single effect model is based on simple univariate Bayesian regression 
applied in the context of multiple regression with the assumption that there is exactly
one variable has a non-zero effect. It extends the single effect assumption to multiple regression
with more than one non-zero effect using a simple and intuitive model fitting procedure based on
variational inference. It provides approximate posterior distributions from which quantities such as
posterior inclusion probability, posterior mean effect size and 95% confidence sets of select variables
can be derived. Particularly, the SuSiE model structure and algorithm results in readily interpretable output 
in the context of genetic fine-mapping studies where variables are highly correlated due to Linkage Disequilibrium.

SuSiE model is developed by [Matthew Stephens and his lab members](http://stephenslab.uchicago.edu/) at the University of Chicago.

This is very much work in progress. Please [post issues](https://github.com/stephenslab/susieR/issues)
to ask questions, get our support or provide us feedback; 
please [send pull requests](https://github.com/stephenslab/susieR/pulls) 
if you have helped fixing bugs or making improvements to the source code.

## Setup

To automatically retrieve and install `susieR` from this repository,

   ```R
   devtools::install_github("stephenslab/susieR")
   ```

## Quick Start

[Here](https://stephenslab.github.io/susieR/articles/mwe.html) is a quick document to show `susieR` in action.
For more documentation and examples please visit: https://stephenslab.github.io/susieR

## Developer notes

+ When any changes are made to `roxygen2` markup, simply run 
`devtools::document()` to update package `NAMESPACE`
and documentation files.

+ Run `pkgdown::build_site()` to build the website. Getting `pkgdown`
to work properly can be frustrating due to numerous & fragile dependencies. 
If `pkgdown` does not work for you out of the box you can use this `docker`
command to run all vignettes and build the site:

```bash
docker run --rm --security-opt label:disable -t -P -w $PWD -v /tmp:/tmp -v $PWD:$PWD \
	-u $UID:${GROUPS[0]} -e HOME=/home/$USER -e USER=$USER gaow/susie \
	R --slave -e "pkgdown::build_site(lazy=TRUE, examples=FALSE)"
```
