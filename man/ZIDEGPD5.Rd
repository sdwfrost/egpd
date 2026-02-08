\name{ZIDEGPD5}
\alias{ZIDEGPD5}
\alias{dZIDEGPD5}
\alias{pZIDEGPD5}
\alias{qZIDEGPD5}
\alias{rZIDEGPD5}
\title{gamlss Family for Zero-Inflated Discrete EGPD Model 5}
\description{
  The \code{ZIDEGPD5()} function defines a \code{gamlss.family} distribution
  for fitting a zero-inflated discrete Extended Generalized Pareto
  Distribution with model 5 (truncated normal G-transformation).

  The functions \code{dZIDEGPD5}, \code{pZIDEGPD5}, \code{qZIDEGPD5}, and
  \code{rZIDEGPD5} define the PMF, distribution function, quantile
  function, and random generation.
}
\usage{
ZIDEGPD5(mu.link = "log", sigma.link = "log", nu.link = "log",
         tau.link = "logit")

dZIDEGPD5(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIDEGPD5(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
          log.p = FALSE)
qZIDEGPD5(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
          log.p = FALSE)
rZIDEGPD5(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
}
\arguments{
  \item{mu.link}{link function for mu, default \code{"log"}}
  \item{sigma.link}{link function for sigma, default \code{"log"}}
  \item{nu.link}{link function for nu, default \code{"log"}}
  \item{tau.link}{link function for tau, default \code{"logit"}}
  \item{x, q}{vector of quantiles (non-negative integers)}
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
  \code{ZIDEGPD5()} returns an object of class \code{"gamlss.family"}.

  \code{dZIDEGPD5} gives the PMF, \code{pZIDEGPD5} gives the distribution
  function, \code{qZIDEGPD5} gives the quantile function, and \code{rZIDEGPD5}
  generates random deviates.
}
\references{
Gamet, P. and Jalbert, J. (2022).
A flexible extended generalized Pareto distribution for tail estimation.
\emph{Environmetrics}, 33(6), e2744.

Ahmad, T. and Hussain, A. (2025).
Flexible model for varying levels of zeros and outliers in count data.
\emph{arXiv preprint} arXiv:2510.27365.
\url{https://arxiv.org/abs/2510.27365}
}
\seealso{
  \code{\link{ZIDEGPD1}}, \code{\link{ZIDEGPD6}}, \code{\link{zidegpd_bamlss}},
  \code{\link{dzidiscegpd}}
}
\examples{
\dontrun{
  library(gamlss)
  set.seed(1)
  y <- rZIDEGPD5(500, mu = 3, sigma = 0.15, nu = 2, tau = 0.25)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                tau.formula = ~ 1, family = ZIDEGPD5())
}
}
