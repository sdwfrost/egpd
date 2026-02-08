\name{rqresid}
\alias{rqresid}
\alias{rqresid.egpd}
\title{Randomized Quantile Residuals}
\description{
  Computes randomized quantile residuals (Dunn & Smyth, 1996) for
  discrete or continuous EGPD models. For discrete models, a uniform
  random variate is drawn between the lower and upper CDF bounds at
  each observation and transformed to the normal scale via
  \code{\link{qnorm}}. For continuous models, the CDF value is transformed
  directly.
}
\usage{
rqresid(object, \dots)

\method{rqresid}{egpd}(object, seed = NULL, \dots)
}
\arguments{
  \item{object}{a fitted \code{\link{egpd}} object}
  \item{seed}{optional random seed for reproducibility}
  \item{\dots}{unused}
}
\value{
  A numeric vector of randomized quantile residuals on the
  standard normal scale. Values that are infinite are replaced with
  \code{NA}.
}
\details{
  For a well-specified model the residuals should be approximately
  standard normal. The generic \code{rqresid} dispatches to the
  \code{rqresid.egpd} method for objects of class \code{"egpd"}.

  The method supports all EGPD families (\code{"egpd"}, \code{"degpd"},
  \code{"zidegpd"}) and all model types (1--6).
}
\references{
Dunn, P. K. and Smyth, G. K. (1996).
Randomized quantile residuals.
\emph{Journal of Computational and Graphical Statistics}, 5(3), 236--244.
}
\seealso{
  \code{\link{egpd}}, \code{\link{predict.egpd}}
}
\examples{
\dontrun{
  set.seed(42)
  dat <- data.frame(y = regpd(200, sigma = 1, xi = 0.2, kappa = 2, type = 1))
  fit <- egpd(y ~ 1, data = dat, family = "egpd", egpd = list(m = 1))
  r <- rqresid(fit, seed = 1)
  qqnorm(r); qqline(r)
}
}
