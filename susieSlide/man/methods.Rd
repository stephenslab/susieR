\name{slider_cs_table}
\alias{slider_cs_table}
\alias{coef.susie_slide}
\alias{predict.susie_slide}
\alias{summary.susie_slide}
\title{Slider coefficients, predictions, and credible-set output}
\usage{
slider_cs_table(fit)
\method{coef}{susie_slide}(object, ...)
\method{predict}{susie_slide}(object, newx = NULL,
  type = c("response", "coefficients"), ...)
\method{summary}{susie_slide}(object, ...)
}
\arguments{
\item{fit,object}{A susieSlide fit.}
\item{newx}{Complete hard-call genotypes with the same SNPs and allele
orientation, in the same order. NULL returns training fitted values.}
\item{type}{Return response predictions or coefficients.}
\item{...}{Reserved additional arguments.}
}
\description{Extracts slider-aware output rather than treating the fit as
an ordinary additive linear predictor.}
\value{coef returns a two-column matrix, additive and heterozygote, with an
intercept row. Predictions equal the intercept plus X times additive
coefficients plus I(X=1) times heterozygote coefficients. predict returns
predictions or that coefficient matrix. slider_cs_table returns one row per
SNP per reported credible set, including the component index, delta,
conditional component probability alpha, SNP PIP, and count-forced status.
summary returns the usual SuSiE summary with an additional delta table.}
