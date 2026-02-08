\name{ZIEGPD5}
\alias{ZIEGPD5}
\alias{dZIEGPD5}
\alias{pZIEGPD5}
\alias{qZIEGPD5}
\alias{rZIEGPD5}
\title{gamlss Family for Zero-Inflated Continuous EGPD Model 5}
\description{
  The \code{ZIEGPD5()} function defines a \code{gamlss.family} distribution
  for fitting a zero-inflated continuous Extended Generalized Pareto
  Distribution with model 5 (truncated normal G-transformation).

  The functions \code{dZIEGPD5}, \code{pZIEGPD5}, \code{qZIEGPD5}, and
  \code{rZIEGPD5} define the density, distribution function, quantile
  function, and random generation.
}
\usage{
ZIEGPD5(mu.link = "log", sigma.link = "log", nu.link = "log",
        tau.link = "logit")

dZIEGPD5(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIEGPD5(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
qZIEGPD5(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
         log.p = FALSE)
rZIEGPD5(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
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
  \code{ZIEGPD5()} returns an object of class \code{"gamlss.family"}.

  \code{dZIEGPD5} gives the density, \code{pZIEGPD5} gives the distribution
  function, \code{qZIEGPD5} gives the quantile function, and \code{rZIEGPD5}
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
  \code{\link{ZIEGPD1}}, \code{\link{ZIEGPD6}}, \code{\link{ziegpd_bamlss}},
  \code{\link{dziegpd}}
}
