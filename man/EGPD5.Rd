\name{EGPD5}
\alias{EGPD5}
\alias{dEGPD5}
\alias{pEGPD5}
\alias{qEGPD5}
\alias{rEGPD5}
\title{gamlss Family for Continuous EGPD Model 5}
\description{
  The \code{EGPD5()} function defines a \code{gamlss.family} distribution
  for fitting a continuous Extended Generalized Pareto Distribution with
  model 5 (truncated normal G-transformation:
  \eqn{G(u) = (\Phi(\sqrt{\kappa}(u-1)) - \Phi(-\sqrt{\kappa})) / (0.5 - \Phi(-\sqrt{\kappa}))}).

  The functions \code{dEGPD5}, \code{pEGPD5}, \code{qEGPD5}, and
  \code{rEGPD5} define the density, distribution function, quantile
  function, and random generation.
}
\usage{
EGPD5(mu.link = "log", sigma.link = "log", nu.link = "log")

dEGPD5(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pEGPD5(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qEGPD5(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rEGPD5(n, mu = 1, sigma = 0.5, nu = 1)
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
  \code{EGPD5()} returns an object of class \code{"gamlss.family"}.

  \code{dEGPD5} gives the density, \code{pEGPD5} gives the distribution
  function, \code{qEGPD5} gives the quantile function, and \code{rEGPD5}
  generates random deviates.
}
\references{
Gamet, P. and Jalbert, J. (2022).
A flexible extended generalized Pareto distribution for tail estimation.
\emph{Environmetrics}, 33(6), e2744.
}
\seealso{
  \code{\link{EGPD1}}, \code{\link{EGPD6}}, \code{\link{egpd_bamlss}},
  \code{\link{pegpd}}
}
\examples{
\dontrun{
  library(gamlss)
  set.seed(1)
  y <- rEGPD5(500, mu = 2, sigma = 0.2, nu = 1.5)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = EGPD5())
}
}
