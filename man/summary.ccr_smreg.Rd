% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccr_smreg.R
\name{summary.ccr_smreg}
\alias{summary.ccr_smreg}
\title{Summary method for \code{ccr_smreg} object}
\usage{
\method{summary}{ccr_smreg}(object, ...)
}
\arguments{
\item{object}{the object of a fitted ccr_smreg model}

\item{...}{for future methods}
}
\value{
An object of class \code{summary.ccr_smreg} with components:
\item{cause1}{a dataframe for the model summary for cause 1}
\item{cause2}{a dataframe for the model summary for cause 2}
}
\description{
Produce a summary of the fitted semiparametric marginal regression model
}
\examples{
\dontrun{
summary(fit)
}
}
\seealso{
\code{\link[ClusteredMPPLE]{print.summary.ccr_smreg}} for printing the summarized results
}
\author{
Wenxian Zhou, \email{wz11 at iu dot edu}

Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
}
