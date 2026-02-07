\name{DEGPD3}
\alias{DEGPD3}
\alias{dDEGPD3}
\alias{pDEGPD3}
\alias{qDEGPD3}
\alias{rDEGPD3}
\title{gamlss Family for Discrete EGPD Model 3}
\description{
  The \code{DEGPD3()} function defines a \code{gamlss.family} distribution
  for fitting a discrete Extended Generalized Pareto Distribution with
  model 3 (incomplete beta G-transformation).

  The functions \code{dDEGPD3}, \code{pDEGPD3}, \code{qDEGPD3}, and
  \code{rDEGPD3} define the PMF, distribution function, quantile
  function, and random generation.
}
\usage{
DEGPD3(mu.link = "log", sigma.link = "log", nu.link = "log")

dDEGPD3(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pDEGPD3(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qDEGPD3(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rDEGPD3(n, mu = 1, sigma = 0.5, nu = 1)
}
\arguments{
  \item{mu.link}{link function for mu, default \code{"log"}}
  \item{sigma.link}{link function for sigma, default \code{"log"}}
  \item{nu.link}{link function for nu, default \code{"log"}}
  \item{x, q}{vector of quantiles (non-negative integers)}
  \item{p}{vector of probabilities}
  \item{n}{number of observations}
  \item{mu}{GPD scale parameter (sigma), positive}
  \item{sigma}{GPD shape parameter (xi), positive}
  \item{nu}{G-transformation parameter (delta), positive}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{DEGPD3()} returns an object of class \code{"gamlss.family"}.

  \code{dDEGPD3} gives the PMF, \code{pDEGPD3} gives the distribution
  function, \code{qDEGPD3} gives the quantile function, and \code{rDEGPD3}
  generates random deviates.
}
\references{
Ahmad, T. and Arshad, I. A. (2024).
New flexible versions of extended generalized Pareto model for count data.
\emph{arXiv preprint} arXiv:2409.18719.
\url{https://arxiv.org/abs/2409.18719}
}
\seealso{
  \code{\link{DEGPD1}}, \code{\link{DEGPD4}}, \code{\link{degpd_bamlss}}
}
