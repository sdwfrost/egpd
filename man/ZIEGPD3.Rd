\name{ZIEGPD3}
\alias{ZIEGPD3}
\alias{dZIEGPD3}
\alias{pZIEGPD3}
\alias{qZIEGPD3}
\alias{rZIEGPD3}
\title{gamlss Family for Zero-Inflated Continuous EGPD Model 3}
\description{
  The \code{ZIEGPD3()} function defines a \code{gamlss.family} distribution
  for fitting a zero-inflated continuous Extended Generalized Pareto
  Distribution with model 3 (incomplete beta G-transformation).

  The functions \code{dZIEGPD3}, \code{pZIEGPD3}, \code{qZIEGPD3}, and
  \code{rZIEGPD3} define the density, distribution function, quantile
  function, and random generation.
}
\usage{
ZIEGPD3(mu.link = "log", sigma.link = "log", nu.link = "log",
        tau.link = "logit")

dZIEGPD3(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIEGPD3(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
qZIEGPD3(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
rZIEGPD3(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
}
\arguments{
  \item{mu.link}{link function for mu, default \code{"log"}}
  \item{sigma.link}{link function for sigma, default \code{"log"}}
  \item{nu.link}{link function for nu, default \code{"log"}}
  \item{tau.link}{link function for tau, default \code{"logit"}}
  \item{x, q}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observations}
  \item{mu}{GPD scale parameter (sigma), positive}
  \item{sigma}{GPD shape parameter (xi), positive}
  \item{nu}{G-transformation parameter (delta), positive}
  \item{tau}{zero-inflation probability (pi), in (0,1)}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{ZIEGPD3()} returns an object of class \code{"gamlss.family"}.

  \code{dZIEGPD3} gives the density, \code{pZIEGPD3} gives the distribution
  function, \code{qZIEGPD3} gives the quantile function, and \code{rZIEGPD3}
  generates random deviates.
}
\references{
Abbas, A., Ahmad, T. and Ahmad, I. (2025).
Modeling zero-inflated precipitation extremes.
\emph{arXiv preprint} arXiv:2504.11058.
\url{https://arxiv.org/abs/2504.11058}
}
\seealso{
  \code{\link{ZIEGPD1}}, \code{\link{ziegpd_bamlss}}, \code{\link{dziegpd}}
}
