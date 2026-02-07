\name{DEGPD1}
\alias{DEGPD1}
\alias{dDEGPD1}
\alias{pDEGPD1}
\alias{qDEGPD1}
\alias{rDEGPD1}
\title{gamlss Family for Discrete EGPD Model 1}
\description{
  The \code{DEGPD1()} function defines a \code{gamlss.family} distribution
  for fitting a discrete Extended Generalized Pareto Distribution with
  model 1 (power G-transformation).

  The functions \code{dDEGPD1}, \code{pDEGPD1}, \code{qDEGPD1}, and
  \code{rDEGPD1} define the PMF, distribution function, quantile
  function, and random generation.
}
\usage{
DEGPD1(mu.link = "log", sigma.link = "log", nu.link = "log")

dDEGPD1(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pDEGPD1(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qDEGPD1(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rDEGPD1(n, mu = 1, sigma = 0.5, nu = 1)
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
  \code{DEGPD1()} returns an object of class \code{"gamlss.family"}.

  \code{dDEGPD1} gives the PMF, \code{pDEGPD1} gives the distribution
  function, \code{qDEGPD1} gives the quantile function, and \code{rDEGPD1}
  generates random deviates.
}
\references{
Ahmad, T. and Arshad, I. A. (2024).
New flexible versions of extended generalized Pareto model for count data.
\emph{arXiv preprint} arXiv:2409.18719.
\url{https://arxiv.org/abs/2409.18719}
}
\seealso{
  \code{\link{DEGPD3}}, \code{\link{DEGPD4}}, \code{\link{degpd_bamlss}},
  \code{\link{ddiscegpd}}
}
\examples{
\dontrun{
  library(gamlss)
  set.seed(1)
  y <- rDEGPD1(500, mu = 3, sigma = 0.15, nu = 2)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = DEGPD1())
}
}
