% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIC.R
\name{computeCorrelation}
\alias{computeCorrelation}
\title{MICor - Matrix of implied correlation for a RAM model in OpenMx.}
\usage{
computeCorrelation(model, latents = TRUE, stdize = TRUE)
}
\arguments{
\item{model}{An Mx model in the RAM formulation}

\item{latents}{Include latent variables? (default TRUE)}

\item{stdize}{Return a correlation matrix instead of a covariance? (default TRUE)}
}
\value{
The covariance or correlation matrix implied by this model.
}
\description{
Computes a correlation or covariance matrix for a RAM model.  You can do this
better by using `mxGetExpected(model, "covariance")`.  This function is here
because a reader asked that we make it.
}
\details{
See Paper.
}
