\name{EGPD3}
\alias{EGPD3}
\alias{dEGPD3}
\alias{pEGPD3}
\alias{qEGPD3}
\alias{rEGPD3}
\title{gamlss Family for Continuous EGPD Model 3}
\description{
  The \code{EGPD3()} function defines a \code{gamlss.family} distribution
  for fitting a continuous Extended Generalized Pareto Distribution with
  model 3 (incomplete beta G-transformation).

  The functions \code{dEGPD3}, \code{pEGPD3}, \code{qEGPD3}, and
  \code{rEGPD3} define the density, distribution function, quantile
  function, and random generation.
}
\usage{
EGPD3(mu.link = "log", sigma.link = "log", nu.link = "log")

dEGPD3(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pEGPD3(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qEGPD3(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rEGPD3(n, mu = 1, sigma = 0.5, nu = 1)
}
\arguments{
  \item{mu.link}{link function for mu, default \code{"log"}}
  \item{sigma.link}{link function for sigma, default \code{"log"}}
  \item{nu.link}{link function for nu, default \code{"log"}}
  \item{x, q}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observations}
  \item{mu}{GPD scale parameter (sigma), positive}
  \item{sigma}{GPD shape parameter (xi), positive}
  \item{nu}{G-transformation parameter (delta), positive}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{EGPD3()} returns an object of class \code{"gamlss.family"}.

  \code{dEGPD3} gives the density, \code{pEGPD3} gives the distribution
  function, \code{qEGPD3} gives the quantile function, and \code{rEGPD3}
  generates random deviates.
}
\references{
Abbas, A., Ahmad, T. and Ahmad, I. (2025).
Modeling zero-inflated precipitation extremes.
\emph{arXiv preprint} arXiv:2504.11058.
\url{https://arxiv.org/abs/2504.11058}
}
\seealso{
  \code{\link{EGPD1}}, \code{\link{EGPD4}}, \code{\link{egpd_bamlss}}
}
