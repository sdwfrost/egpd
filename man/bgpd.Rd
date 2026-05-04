\name{bgpd}
\alias{pbgpd}
\alias{dbgpd}
\alias{qbgpd}
\alias{rbgpd}
\title{Blended generalized Pareto distribution functions}
\usage{
pbgpd(
  q,
  xi,
  pbulk = pbeta,
  dbulk = dbeta,
  qbulk = qbeta,
  bulk.args = list(shape1 = 2, shape2 = 5),
  pa = 0.8,
  pb = 0.99,
  c1 = 5,
  c2 = 5
)

dbgpd(
  x,
  xi,
  pbulk = pbeta,
  dbulk = dbeta,
  qbulk = qbeta,
  bulk.args = list(shape1 = 2, shape2 = 5),
  pa = 0.8,
  pb = 0.99,
  c1 = 5,
  c2 = 5,
  log = FALSE
)

qbgpd(
  p,
  xi,
  pbulk = pbeta,
  dbulk = dbeta,
  qbulk = qbeta,
  bulk.args = list(shape1 = 2, shape2 = 5),
  pa = 0.8,
  pb = 0.99,
  c1 = 5,
  c2 = 5
)

rbgpd(
  n,
  xi,
  pbulk = pbeta,
  dbulk = dbeta,
  qbulk = qbeta,
  bulk.args = list(shape1 = 2, shape2 = 5),
  pa = 0.8,
  pb = 0.99,
  c1 = 5,
  c2 = 5,
  unifsamp = NULL
)
}
\arguments{
\item{q}{quantiles.}

\item{xi}{GP tail shape parameter.}

\item{pbulk}{bulk CDF supported on \eqn{[0,1]}.}

\item{dbulk}{bulk density supported on \eqn{[0,1]}. It is also used to
numerically validate the blended construction.}

\item{qbulk}{bulk quantile function supported on \eqn{[0,1]}.}

\item{bulk.args}{named list of arguments passed to \code{pbulk},
\code{dbulk}, and \code{qbulk}. The default corresponds to a Beta bulk.}

\item{pa}{lower bulk quantile where blending starts.}

\item{pb}{upper bulk quantile where blending ends and the exact GP tail
begins.}

\item{c1, c2}{shape parameters of the Beta weighting function used for
blending.}

\item{x}{values.}

\item{log}{logical; should log-densities be returned?}

\item{p}{probabilities.}

\item{n}{number of samples.}

\item{unifsamp}{optional uniforms to transform instead of drawing new ones.}
}
\value{
\code{pbgpd()} returns CDF values.

\code{dbgpd()} returns density values.

\code{qbgpd()} returns quantiles.

\code{rbgpd()} returns random draws.
}
\description{
Distribution utilities for the blended generalized Pareto construction of
Majumder and Richards (2026). A finite-support bulk model on \eqn{[0,1]} is
blended with a GP tail using the paper's multiplicative CDF blend and
Beta-CDF weighting function.
}
\details{
The GP threshold and scale are derived from \code{pa} and \code{pb} so that
the blended CDF matches the bulk CDF at the start and end of the blending
interval. For \eqn{y \le a}, the distribution equals the bulk model; for
\eqn{y \ge b}, it equals the exact GP model.

Because the paper notes that some bulk/weight combinations can produce a
negative density in the blending interval, these helpers numerically screen
the blended density over that interval and stop if the supplied specification
appears invalid.
}
\examples{
## Default Beta-bulk blended GP
pbgpd(c(0.2, 0.8, 1.5), xi = 0.2)
dbgpd(c(0.2, 0.8, 1.5), xi = 0.2)
qbgpd(c(0.5, 0.9, 0.99), xi = 0.2)

## Custom bulk distribution on [0, 1]
qbgpd(
  c(0.5, 0.9, 0.99),
  xi = -0.15,
  pbulk = pbeta,
  dbulk = dbeta,
  qbulk = qbeta,
  bulk.args = list(shape1 = 3, shape2 = 7),
  pa = 0.7,
  pb = 0.95
)
}
