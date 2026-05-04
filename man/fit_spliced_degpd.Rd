\name{fit_spliced_degpd}
\alias{fit_spliced_degpd}
\alias{print.spliced_degpd}
\alias{summary.spliced_degpd}
\alias{print.summary.spliced_degpd}
\alias{predict.spliced_degpd}
\title{Fit a spliced discrete EGPD forecasting model}
\usage{
fit_spliced_degpd(
  formula,
  data,
  bulk_probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
  threshold_prob = 0.9,
  fit_tail = TRUE,
  tail_formula = NULL,
  degpd.args = list(m = 1),
  qgam_control = list(),
  egpd_control = list(),
  keep_data = FALSE
)

\method{print}{spliced_degpd}(x, ...)

\method{summary}{spliced_degpd}(object, ...)

\method{print}{summary.spliced_degpd}(x, ...)

\method{predict}{spliced_degpd}(
  object,
  newdata,
  type = c("quantile", "cdf", "pmf", "parameter", "bulk"),
  prob = NULL,
  at = NULL,
  ...
)
}
\arguments{
\item{formula}{a formula for the bulk additive quantile regression.}

\item{data}{a data frame.}

\item{bulk_probs}{probabilities used for the bulk quantile regression fit.}

\item{threshold_prob}{probability threshold where the tail model begins.}

\item{fit_tail}{logical; should a DEGPD tail be fitted above the threshold?}

\item{tail_formula}{NULL, a formula, or a named list of formulae for the tail
model. If a single formula is supplied, it is used for \code{lsigma} and the
shape and kappa submodels default to \code{~1}. If a list is supplied, valid
names are \code{lsigma}, \code{lxi}, and \code{lkappa}.}

\item{degpd.args}{list of DEGPD family arguments for the tail fit. Currently
only \code{m = 1} is supported.}

\item{qgam_control}{optional control list passed to \code{qgam::mqgam()}.}

\item{egpd_control}{optional list of tail optimizer options. Supported
entries are \code{start}, \code{fix.arg}, \code{optim.method}, and
\code{control}.}

\item{keep_data}{logical; should the training data be stored in the fitted
object?}

\item{x}{an object of class \code{spliced_degpd} or
\code{summary.spliced_degpd}.}

\item{object}{an object of class \code{spliced_degpd}.}

\item{newdata}{a data frame of covariates. If omitted, the original training
data are used only when \code{keep_data = TRUE}.}

\item{type}{prediction type: \code{"quantile"}, \code{"cdf"},
\code{"pmf"}, \code{"parameter"}, or \code{"bulk"}.}

\item{prob}{probabilities used for \code{type = "quantile"}.}

\item{at}{non-negative integer support values used for
\code{type = "cdf"} or \code{"pmf"}.}

\item{\dots}{unused additional arguments.}
}
\value{
\code{fit_spliced_degpd()} returns an object of class
\code{"spliced_degpd"} containing the fitted bulk quantile model, the fitted
DEGPD tail model, the quantile threshold, and the modeling metadata.

\code{predict.spliced_degpd(type = "quantile")} returns a data frame of
predicted quantiles.

\code{predict.spliced_degpd(type = "cdf")} returns a data frame of predicted
CDF values on the supplied support.

\code{predict.spliced_degpd(type = "pmf")} returns a data frame of predicted
probability masses on the supplied support.

\code{predict.spliced_degpd(type = "parameter")} returns the threshold shift
and fitted tail parameters.
}
\description{
Fits a hybrid count-distribution model inspired by the X-flexForecast
framework of Maia, Castro-Camilo, and Browell (2026): multiple additive
quantile regressions describe the bulk of the count distribution, while
threshold exceedances are modeled with a type-1 DEGPD tail.
}
\details{
This implementation adapts the paper's bulk-and-tail splicing idea to the
\pkg{egpd} package. The bulk is fitted with \pkg{qgam}'s additive quantile
regression, while the tail is fitted with a native type-1 DEGPD likelihood
using the package's existing discrete distribution functions. By default the
tail optimizer fixes \code{kappa = 1}, which recovers the paper's DGP-style
tail; this can be changed through \code{egpd_control$fix.arg}.

The \code{qgam} package is an optional dependency and must be installed to use
\code{fit_spliced_degpd()}.
}
\examples{
\dontrun{
if (requireNamespace("qgam", quietly = TRUE)) {
  set.seed(1)
  n <- 400
  x <- runif(n, -1, 1)
  y <- rpois(n, lambda = exp(1 + x)) +
    rbinom(n, 1, 0.15) * rdiscegpd(n, sigma = 3, xi = 0.2, kappa = 1.4)
  df <- data.frame(y = y, x = x)

  fit <- fit_spliced_degpd(
    y ~ s(x, k = 6),
    data = df,
    bulk_probs = c(0.1, 0.5, 0.8, 0.9),
    threshold_prob = 0.8,
    tail_formula = list(lsigma = ~x, lxi = ~1, lkappa = ~1)
  )

  predict(fit, newdata = df[1:5, ], type = "quantile", prob = c(0.5, 0.9, 0.99))
  predict(fit, newdata = df[1:2, ], type = "pmf", at = 0:20)
}
}
}
