\name{DEGPD5}
\alias{DEGPD5}
\alias{dDEGPD5}
\alias{pDEGPD5}
\alias{qDEGPD5}
\alias{rDEGPD5}
\title{gamlss Family for Discrete EGPD Model 5}
\description{
  The \code{DEGPD5()} function defines a \code{gamlss.family} distribution
  for fitting a discrete Extended Generalized Pareto Distribution with
  model 5 (truncated normal G-transformation).

  The functions \code{dDEGPD5}, \code{pDEGPD5}, \code{qDEGPD5}, and
  \code{rDEGPD5} define the PMF, distribution function, quantile
  function, and random generation.
}
\usage{
DEGPD5(mu.link = "log", sigma.link = "log", nu.link = "log")

dDEGPD5(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pDEGPD5(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qDEGPD5(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rDEGPD5(n, mu = 1, sigma = 0.5, nu = 1)
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
  \item{nu}{G-transformation parameter (kappa), positive}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{DEGPD5()} returns an object of class \code{"gamlss.family"}.

  \code{dDEGPD5} gives the PMF, \code{pDEGPD5} gives the distribution
  function, \code{qDEGPD5} gives the quantile function, and \code{rDEGPD5}
  generates random deviates.
}
\references{
Gamet, P. and Jalbert, J. (2022).
A flexible extended generalized Pareto distribution for tail estimation.
\emph{Environmetrics}, 33(6), e2744.

Ahmad, T. and Arshad, I. A. (2024).
New flexible versions of extended generalized Pareto model for count data.
\emph{arXiv preprint} arXiv:2409.18719.
\url{https://arxiv.org/abs/2409.18719}
}
\seealso{
  \code{\link{DEGPD1}}, \code{\link{DEGPD6}}, \code{\link{degpd_bamlss}},
  \code{\link{ddiscegpd}}
}
\examples{
\dontrun{
  library(gamlss)
  set.seed(1)
  y <- rDEGPD5(500, mu = 3, sigma = 0.15, nu = 2)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = DEGPD5())
}
}
