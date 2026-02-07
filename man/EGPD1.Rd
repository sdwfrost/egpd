\name{EGPD1}
\alias{EGPD1}
\alias{dEGPD1}
\alias{pEGPD1}
\alias{qEGPD1}
\alias{rEGPD1}
\title{gamlss Family for Continuous EGPD Model 1}
\description{
  The \code{EGPD1()} function defines a \code{gamlss.family} distribution
for fitting a continuous Extended Generalized Pareto Distribution with
model 1 (power G-transformation: \eqn{G(u) = u^\kappa}).

  The functions \code{dEGPD1}, \code{pEGPD1}, \code{qEGPD1}, and
  \code{rEGPD1} define the density, distribution function, quantile
  function, and random generation for this distribution using gamlss
  parameter names.
}
\usage{
EGPD1(mu.link = "log", sigma.link = "log", nu.link = "log")

dEGPD1(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)
pEGPD1(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
qEGPD1(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)
rEGPD1(n, mu = 1, sigma = 0.5, nu = 1)
}
\arguments{
  \item{mu.link}{link function for mu, default \code{"log"}}
  \item{sigma.link}{link function for sigma, default \code{"log"}}
  \item{nu.link}{link function for nu, default \code{"log"}}
  \item{x, q}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observations}
  \item{mu}{GPD scale parameter (sigma in EGPD notation), positive}
  \item{sigma}{GPD shape parameter (xi in EGPD notation), positive}
  \item{nu}{G-transformation parameter (kappa in EGPD notation), positive}
  \item{log, log.p}{logical; if TRUE, probabilities/densities given as log}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x)}
}
\value{
  \code{EGPD1()} returns an object of class \code{"gamlss.family"}.

  \code{dEGPD1} gives the density, \code{pEGPD1} gives the distribution
  function, \code{qEGPD1} gives the quantile function, and \code{rEGPD1}
  generates random deviates.
}
\details{
  The gamlss parameters map to EGPD parameters as:
  \code{mu} = sigma (GPD scale), \code{sigma} = xi (GPD shape),
  \code{nu} = kappa (G-transformation power).

  First derivatives are computed numerically via central finite differences;
  second derivatives use the outer product of gradients (OPG) approximation.
}
\references{
Abbas, A., Ahmad, T. and Ahmad, I. (2025).
Modeling zero-inflated precipitation extremes.
\emph{arXiv preprint} arXiv:2504.11058.
\url{https://arxiv.org/abs/2504.11058}
}
\seealso{
  \code{\link{egpd_bamlss}}, \code{\link{degpd_density}}, \code{\link{pegpd}}
}
\examples{
\dontrun{
  library(gamlss)
  set.seed(1)
  y <- rEGPD1(500, mu = 2, sigma = 0.2, nu = 1.5)
  fit <- gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                family = EGPD1())
}
}
