% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustered_mpple_se_bs.R
\name{clustered_mpple_se_bs}
\alias{clustered_mpple_se_bs}
\title{Obtain coefficient estimates standard error using bootstrap method}
\usage{
clustered_mpple_se_bs(
  data,
  formula1 = y ~ x + Z1 + Z2,
  formula2 = Surv(x, d) ~ Z1 + Z2,
  w = TRUE,
  nboot = 1000,
  x,
  t,
  ...
)
}
\arguments{
\item{data}{a data frame in suitable format.}

\item{formula1}{a formula relating the response for cause 1 \code{y} to a set of covariates.}

\item{formula2}{a formula relating the survival object \code{Surv(x, d)} to a set of covariates.}

\item{w}{logical value: if TRUE, the estimation procedure is weighed by cluster size.}

\item{nboot}{a number of bootstrap samples for estimating variances when bootstrapping methods are used.}

\item{x}{a vector of time points for cumulative incidence function.}

\item{t}{a vector of time points for cumulative residual process.}

\item{...}{for future methods.}
}
\value{
A list with components:
\item{se1}{a list of bootstrap samples for cause 1}
\item{se2}{a list of bootstrap samples for cause 1}
}
\description{
Obtain coefficient estimates standard error for semiparametric regression for competing risks data
with missing cause of failure using bootstrap.
}
\author{
Wenxian Zhou, \email{wz11 at iu dot edu}

Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
}
