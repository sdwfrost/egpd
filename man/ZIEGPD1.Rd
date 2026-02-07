\name{ZIEGPD1}
\alias{ZIEGPD1}
\alias{dZIEGPD1}
\alias{pZIEGPD1}
\alias{qZIEGPD1}
\alias{rZIEGPD1}
\title{gamlss Family for Zero-Inflated Continuous EGPD Model 1}
\description{
  The \code{ZIEGPD1()} function defines a \code{gamlss.family} distribution
  for fitting a zero-inflated continuous Extended Generalized Pareto
  Distribution with model 1 (power G-transformation).

  The functions \code{dZIEGPD1}, \code{pZIEGPD1}, \code{qZIEGPD1}, and
  \code{rZIEGPD1} define the density, distribution function, quantile
  function, and random generation.
}
\usage{
ZIEGPD1(mu.link = "log", sigma.link = "log", nu.link = "log",
        tau.link = "logit")

dZIEGPD1(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIEGPD1(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
qZIEGPD1(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
rZIEGPD1(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
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
  \item{nu}{G-transformation parameter (kappa), positive}
  \item{tau}{zero-inflation probability (pi), in (0,1)}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{ZIEGPD1()} returns an object of class \code{"gamlss.family"}.

  \code{dZIEGPD1} gives the density, \code{pZIEGPD1} gives the distribution
  function, \code{qZIEGPD1} gives the quantile function, and \code{rZIEGPD1}
  generates random deviates.
}
\references{
Abbas, A., Ahmad, T. and Ahmad, I. (2025).
Modeling zero-inflated precipitation extremes.
\emph{arXiv preprint} arXiv:2504.11058.
\url{https://arxiv.org/abs/2504.11058}
}
\seealso{
  \code{\link{ZIEGPD3}}, \code{\link{ziegpd_bamlss}}, \code{\link{dziegpd}}
}
