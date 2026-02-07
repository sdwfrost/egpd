\name{DEGPD4}
\alias{DEGPD4}
\alias{dDEGPD4}
\alias{pDEGPD4}
\alias{qDEGPD4}
\alias{rDEGPD4}
\title{gamlss Family for Discrete EGPD Model 4}
\description{
  The \code{DEGPD4()} function defines a \code{gamlss.family} distribution
  for fitting a discrete Extended Generalized Pareto Distribution with
  model 4 (power-beta G-transformation). This is a 4-parameter family.

  The functions \code{dDEGPD4}, \code{pDEGPD4}, \code{qDEGPD4}, and
  \code{rDEGPD4} define the PMF, distribution function, quantile
  function, and random generation.
}
\usage{
DEGPD4(mu.link = "log", sigma.link = "log", nu.link = "log",
       tau.link = "log")

dDEGPD4(x, mu = 1, sigma = 0.5, nu = 1, tau = 1, log = FALSE)
pDEGPD4(q, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE,
        log.p = FALSE)
qDEGPD4(p, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE,
        log.p = FALSE)
rDEGPD4(n, mu = 1, sigma = 0.5, nu = 1, tau = 1)
}
\arguments{
  \item{mu.link}{link function for mu, default \code{"log"}}
  \item{sigma.link}{link function for sigma, default \code{"log"}}
  \item{nu.link}{link function for nu, default \code{"log"}}
  \item{tau.link}{link function for tau, default \code{"log"}}
  \item{x, q}{vector of quantiles (non-negative integers)}
  \item{p}{vector of probabilities}
  \item{n}{number of observations}
  \item{mu}{GPD scale parameter (sigma), positive}
  \item{sigma}{GPD shape parameter (xi), positive}
  \item{nu}{G-transformation parameter (delta), positive}
  \item{tau}{G-transformation parameter (kappa), positive}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{DEGPD4()} returns an object of class \code{"gamlss.family"}.

  \code{dDEGPD4} gives the PMF, \code{pDEGPD4} gives the distribution
  function, \code{qDEGPD4} gives the quantile function, and \code{rDEGPD4}
  generates random deviates.
}
\references{
Ahmad, T. and Arshad, I. A. (2024).
New flexible versions of extended generalized Pareto model for count data.
\emph{arXiv preprint} arXiv:2409.18719.
\url{https://arxiv.org/abs/2409.18719}
}
\seealso{
  \code{\link{DEGPD1}}, \code{\link{DEGPD3}}, \code{\link{degpd_bamlss}}
}
