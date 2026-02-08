\name{EGPD6}
\alias{EGPD6}
\alias{dEGPD6}
\alias{pEGPD6}
\alias{qEGPD6}
\alias{rEGPD6}
\title{gamlss Family for Continuous EGPD Model 6}
\description{
  The \code{EGPD6()} function defines a \code{gamlss.family} distribution
  for fitting a continuous Extended Generalized Pareto Distribution with
  model 6 (truncated beta G-transformation).

  The functions \code{dEGPD6}, \code{pEGPD6}, \code{qEGPD6}, and
  \code{rEGPD6} define the density, distribution function, quantile
  function, and random generation.
}
\usage{
EGPD6(mu.link = "log", sigma.link = "log", nu.link = "log")

dEGPD6(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pEGPD6(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qEGPD6(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rEGPD6(n, mu = 1, sigma = 0.5, nu = 1)
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
  \item{nu}{G-transformation parameter (kappa), positive}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{EGPD6()} returns an object of class \code{"gamlss.family"}.

  \code{dEGPD6} gives the density, \code{pEGPD6} gives the distribution
  function, \code{qEGPD6} gives the quantile function, and \code{rEGPD6}
  generates random deviates.
}
\references{
Gamet, P. and Jalbert, J. (2022).
A flexible extended generalized Pareto distribution for tail estimation.
\emph{Environmetrics}, 33(6), e2744.
}
\seealso{
  \code{\link{EGPD1}}, \code{\link{EGPD5}}, \code{\link{egpd_bamlss}},
  \code{\link{pegpd}}
}
\examples{
\dontrun{
  library(gamlss)
  set.seed(1)
  y <- rEGPD6(500, mu = 2, sigma = 0.2, nu = 1.5)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = EGPD6())
}
}
