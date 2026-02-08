\name{ZIEGPD6}
\alias{ZIEGPD6}
\alias{dZIEGPD6}
\alias{pZIEGPD6}
\alias{qZIEGPD6}
\alias{rZIEGPD6}
\title{gamlss Family for Zero-Inflated Continuous EGPD Model 6}
\description{
  The \code{ZIEGPD6()} function defines a \code{gamlss.family} distribution
  for fitting a zero-inflated continuous Extended Generalized Pareto
  Distribution with model 6 (truncated beta G-transformation).

  The functions \code{dZIEGPD6}, \code{pZIEGPD6}, \code{qZIEGPD6}, and
  \code{rZIEGPD6} define the density, distribution function, quantile
  function, and random generation.
}
\usage{
ZIEGPD6(mu.link = "log", sigma.link = "log", nu.link = "log",
        tau.link = "logit")

dZIEGPD6(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIEGPD6(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
qZIEGPD6(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
rZIEGPD6(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
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
  \code{ZIEGPD6()} returns an object of class \code{"gamlss.family"}.

  \code{dZIEGPD6} gives the density, \code{pZIEGPD6} gives the distribution
  function, \code{qZIEGPD6} gives the quantile function, and \code{rZIEGPD6}
  generates random deviates.
}
\references{
Gamet, P. and Jalbert, J. (2022).
A flexible extended generalized Pareto distribution for tail estimation.
\emph{Environmetrics}, 33(6), e2744.

Abbas, A., Ahmad, T. and Ahmad, I. (2025).
Modeling zero-inflated precipitation extremes.
\emph{arXiv preprint} arXiv:2504.11058.
\url{https://arxiv.org/abs/2504.11058}
}
\seealso{
  \code{\link{ZIEGPD1}}, \code{\link{ZIEGPD5}}, \code{\link{ziegpd_bamlss}},
  \code{\link{dziegpd}}
}
