% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rDynMix.R
\name{rDynMix}
\alias{rDynMix}
\title{Simulating a dynamic lognormal-Pareto mixture}
\usage{
rDynMix(nreps, x)
}
\arguments{
\item{nreps}{integer: number of observations sampled from the mixture.}

\item{x}{(6 by 1) numerical vector: values of mu_c, tau, mu, sigma, xi, beta.}
}
\value{
ysim (nreps x 1) vector: nreps random numbers from the lognormal-GPD dynamic mixture.
}
\description{
This function fits a dynamic mixture by Approximate Maximum Likelihood and by standard maximum likelihood.
Currently only implemented for the lognormal - generalized Pareto case.
}
\details{
This function simulates a dynamic lognormal-GPD mixture using
the algorithm of Frigessi et al. (2002, p. 221).
}
\examples{
ysim <- rDynMix(100,c(1,2,0,0.5,0.25,3))
}
\references{
{
\insertRef{fri02}{FitDynMix}
}
}
\keyword{dynamic}
\keyword{mixture;}
\keyword{simulation.}
