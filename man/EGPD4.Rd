\name{EGPD4}
\alias{EGPD4}
\alias{dEGPD4}
\alias{pEGPD4}
\alias{qEGPD4}
\alias{rEGPD4}
\title{gamlss Family for Continuous EGPD Model 4}
\description{
  The \code{EGPD4()} function defines a \code{gamlss.family} distribution
  for fitting a continuous Extended Generalized Pareto Distribution with
  model 4 (power-beta G-transformation). This is a 4-parameter family.

  The functions \code{dEGPD4}, \code{pEGPD4}, \code{qEGPD4}, and
  \code{rEGPD4} define the density, distribution function, quantile
  function, and random generation.
}
\usage{
EGPD4(mu.link = "log", sigma.link = "log", nu.link = "log",
      tau.link = "log")

dEGPD4(x, mu = 1, sigma = 0.5, nu = 1, tau = 1, log = FALSE)
pEGPD4(q, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE,
       log.p = FALSE)
qEGPD4(p, mu = 1, sigma = 0.5, nu = 1, tau = 1, lower.tail = TRUE,
       log.p = FALSE)
rEGPD4(n, mu = 1, sigma = 0.5, nu = 1, tau = 1)
}
\arguments{
  \item{mu.link}{link function for mu, default \code{"log"}}
  \item{sigma.link}{link function for sigma, default \code{"log"}}
  \item{nu.link}{link function for nu, default \code{"log"}}
  \item{tau.link}{link function for tau, default \code{"log"}}
  \item{x, q}{vector of quantiles}
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
  \code{EGPD4()} returns an object of class \code{"gamlss.family"}.

  \code{dEGPD4} gives the density, \code{pEGPD4} gives the distribution
  function, \code{qEGPD4} gives the quantile function, and \code{rEGPD4}
  generates random deviates.
}
\references{
Abbas, A., Ahmad, T. and Ahmad, I. (2025).
Modeling zero-inflated precipitation extremes.
\emph{arXiv preprint} arXiv:2504.11058.
\url{https://arxiv.org/abs/2504.11058}
}
\seealso{
  \code{\link{EGPD1}}, \code{\link{EGPD3}}, \code{\link{egpd_bamlss}}
}
