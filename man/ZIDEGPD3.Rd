\name{ZIDEGPD3}
\alias{ZIDEGPD3}
\alias{dZIDEGPD3}
\alias{pZIDEGPD3}
\alias{qZIDEGPD3}
\alias{rZIDEGPD3}
\title{gamlss Family for Zero-Inflated Discrete EGPD Model 3}
\description{
  The \code{ZIDEGPD3()} function defines a \code{gamlss.family} distribution
  for fitting a zero-inflated discrete Extended Generalized Pareto
  Distribution with model 3 (incomplete beta G-transformation).

  The functions \code{dZIDEGPD3}, \code{pZIDEGPD3}, \code{qZIDEGPD3}, and
  \code{rZIDEGPD3} define the PMF, distribution function, quantile
  function, and random generation.
}
\usage{
ZIDEGPD3(mu.link = "log", sigma.link = "log", nu.link = "log",
         tau.link = "logit")

dZIDEGPD3(x, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, log = FALSE)
pZIDEGPD3(q, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
          log.p = FALSE)
qZIDEGPD3(p, mu = 1, sigma = 0.5, nu = 1, tau = 0.1, lower.tail = TRUE,
          log.p = FALSE)
rZIDEGPD3(n, mu = 1, sigma = 0.5, nu = 1, tau = 0.1)
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
  \item{nu}{G-transformation parameter (delta), positive}
  \item{tau}{zero-inflation probability (pi), in (0,1)}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{ZIDEGPD3()} returns an object of class \code{"gamlss.family"}.

  \code{dZIDEGPD3} gives the PMF, \code{pZIDEGPD3} gives the distribution
  function, \code{qZIDEGPD3} gives the quantile function, and \code{rZIDEGPD3}
  generates random deviates.
}
\references{
Ahmad, T. and Arshad, I. A. (2024).
New flexible versions of extended generalized Pareto model for count data.
\emph{arXiv preprint} arXiv:2409.18719.
\url{https://arxiv.org/abs/2409.18719}

Ahmad, T. and Hussain, A. (2025).
Flexible model for varying levels of zeros and outliers in count data.
\emph{arXiv preprint} arXiv:2510.27365.
\url{https://arxiv.org/abs/2510.27365}
}
\seealso{
  \code{\link{ZIDEGPD1}}, \code{\link{zidegpd_bamlss}},
  \code{\link{dzidiscegpd}}
}
