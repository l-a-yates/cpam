% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_changepoint.R
\name{gcv.gam}
\alias{gcv.gam}
\title{Calculate pointwise GCV for GAM models}
\usage{
gcv.gam(fit)
}
\arguments{
\item{fit}{A fitted GAM object}
}
\value{
Vector of pointwise GCV scores
}
\description{
Calculate pointwise GCV for GAM models
}
\details{
Implements the standard GCV calculation for GAM models on a pointwise basis.
}
\keyword{internal}
