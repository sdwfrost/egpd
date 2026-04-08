\name{fit_aggregate_egpd}
\alias{fit_aggregate_egpd}
\alias{coef.aggregate_egpd}
\alias{print.aggregate_egpd}
\alias{summary.aggregate_egpd}
\alias{print.summary.aggregate_egpd}
\alias{predict.aggregate_egpd}
\alias{plot.aggregate_egpd}
\title{Fit a multiscale aggregated EGPD framework}
\usage{
fit_aggregate_egpd(
  data,
  durations,
  family = c("cpegpd", "cpdegpd"),
  p = 2L,
  q = 2L,
  overlap = TRUE,
  cpegpd.h = 0.2,
  start = NULL,
  kappa_bounds = c(0.05, 5),
  xi_bounds = c(0.01, 1.5),
  optim.method = "Nelder-Mead",
  control = list(),
  keep_data = FALSE
)

\method{coef}{aggregate_egpd}(object, ...)

\method{print}{aggregate_egpd}(x, ...)

\method{summary}{aggregate_egpd}(object, ...)

\method{print}{summary.aggregate_egpd}(x, digits = 4, ...)

\method{predict}{aggregate_egpd}(
  object,
  durations = object$durations,
  type = c("parameter", "quantile"),
  prob = NULL,
  ...
)

\method{plot}{aggregate_egpd}(x, ...)
}
\arguments{
\item{data}{either a numeric vector of non-negative observations at the
finest time scale, or a list of non-negative aggregated samples for each
duration.}

\item{durations}{positive integer aggregation durations.}

\item{family}{either \code{"cpegpd"} for continuous aggregates or
\code{"cpdegpd"} for discrete aggregates.}

\item{p}{degree of the polynomial used for \eqn{\log \sigma_d}.}

\item{q}{degree of the polynomial used for \eqn{\log \lambda_d}.}

\item{overlap}{logical; if \code{TRUE}, rolling sums are used when
\code{data} is a finest-scale series.}

\item{cpegpd.h}{discretization step used by the Panjer recursion for the
continuous compound family.}

\item{start}{optional named list with entries \code{sigma_coef},
\code{lambda_coef}, \code{kappa}, and \code{xi}.}

\item{kappa_bounds}{numeric vector of length two giving lower and upper
bounds for \eqn{\kappa}.}

\item{xi_bounds}{numeric vector of length two giving lower and upper
bounds for \eqn{\xi}.}

\item{optim.method}{optimization method passed to \code{\link{optim}}.}

\item{control}{optional optimization control list.}

\item{keep_data}{logical; should aggregated samples be stored in the fitted
object?}

\item{object}{an object of class \code{aggregate_egpd}.}

\item{x}{an object of class \code{aggregate_egpd} or
\code{summary.aggregate_egpd}.}

\item{digits}{number of digits used when printing summaries.}

\item{type}{either \code{"parameter"} or \code{"quantile"}.}

\item{prob}{probabilities used when \code{type = "quantile"}.}

\item{\dots}{ignored additional arguments.}
}
\value{
\code{fit_aggregate_egpd()} returns an object of class
\code{"aggregate_egpd"} with the fitted polynomial coefficients for
\eqn{\log \sigma_d} and \eqn{\log \lambda_d}, the shared shape
parameters \eqn{\kappa} and \eqn{\xi}, the fitted duration-specific
parameters, and the composite log-likelihood.

\code{predict.aggregate_egpd(type = "parameter")} returns a data frame
with the fitted parameters at each duration.

\code{predict.aggregate_egpd(type = "quantile")} returns a matrix of
positive-part quantiles, with rows corresponding to durations and columns
to probabilities.
}
\description{
Fits the multiscale aggregated EGPD framework of Ailliot, Gaetan, and
Naveau (2026). Positive aggregates at duration \eqn{d} are modeled as a
compound Poisson EGPD law with shared shape parameters and duration-wise
scale and Poisson mean:
\deqn{A_d \sim EGPD(\sigma_d, \kappa, \xi, \lambda_d)}
where \eqn{\log \sigma_d} and \eqn{\log \lambda_d} are polynomials in
\eqn{\log d}. The fit is obtained by maximizing a composite likelihood
across aggregation durations.
}
\details{
This implementation follows the manuscript-level framework for positive
aggregated rainfall amounts. It uses the positive part of the
compound-Poisson EGPD or compound-Poisson discrete EGPD distributions
already implemented in the package. Monotone fitted curves for
\eqn{\sigma_d} and \eqn{\lambda_d} are enforced numerically during
optimization to preserve the non-crossing behavior of fitted return
levels across durations.
}
\examples{
\dontrun{
set.seed(1)
durs <- c(1, 2, 4, 8)
agg <- setNames(lapply(durs, function(d) {
  rcpegpd(200,
    sigma = exp(log(0.6) + 0.12 * log(d)),
    xi = 0.2,
    kappa = 0.35,
    lambda = exp(log(0.8) + 0.25 * log(d)),
    type = 1
  )
}), durs)

fit <- fit_aggregate_egpd(
  agg,
  durations = durs,
  family = "cpegpd",
  p = 1,
  q = 1
)

summary(fit)
predict(fit, type = "parameter")
predict(fit, type = "quantile", prob = c(0.5, 0.9))
plot(fit)
}
}
