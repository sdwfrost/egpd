\name{DEGPD6}
\alias{DEGPD6}
\alias{dDEGPD6}
\alias{pDEGPD6}
\alias{qDEGPD6}
\alias{rDEGPD6}
\title{gamlss Family for Discrete EGPD Model 6}
\description{
  The \code{DEGPD6()} function defines a \code{gamlss.family} distribution
  for fitting a discrete Extended Generalized Pareto Distribution with
  model 6 (truncated beta G-transformation).

  The functions \code{dDEGPD6}, \code{pDEGPD6}, \code{qDEGPD6}, and
  \code{rDEGPD6} define the PMF, distribution function, quantile
  function, and random generation.
}
\usage{
DEGPD6(mu.link = "log", sigma.link = "log", nu.link = "log")

dDEGPD6(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pDEGPD6(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qDEGPD6(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rDEGPD6(n, mu = 1, sigma = 0.5, nu = 1)
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
  \code{DEGPD6()} returns an object of class \code{"gamlss.family"}.

  \code{dDEGPD6} gives the PMF, \code{pDEGPD6} gives the distribution
  function, \code{qDEGPD6} gives the quantile function, and \code{rDEGPD6}
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
  \code{\link{DEGPD1}}, \code{\link{DEGPD5}}, \code{\link{degpd_bamlss}},
  \code{\link{ddiscegpd}}
}
\examples{
\dontrun{
  library(gamlss)
  set.seed(1)
  y <- rDEGPD6(500, mu = 3, sigma = 0.15, nu = 2)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = DEGPD6())
}
}
