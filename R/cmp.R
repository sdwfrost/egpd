## Conway-Maxwell-Poisson (CMP), mean-parameterised. See src/cmp.cpp for why Z is computed
## exactly (Gaunt's asymptotic is invalid at the nu these data force) and why that is still
## affordable (nu is intercept-only, so the exact sums are gridded once per likelihood call).

#' Conway-Maxwell-Poisson distribution
#'
#' \code{P(Y = y) = lambda^y / ((y!)^nu Z(lambda, nu))}, in the mean parameterisation:
#' \code{lambda} is chosen so that \code{E[Y] = mu}. \code{nu < 1} is overdispersed,
#' \code{nu = 1} is Poisson, \code{nu > 1} underdispersed. The tail decays like
#' \code{exp(-nu y log y)}, lighter than geometric for any \code{nu > 0}, so the
#' extreme-value index is \code{0}: CMP models dispersion, not a heavy tail.
#'
#' @param x,q vector of counts
#' @param p vector of probabilities
#' @param n number of draws
#' @param mu mean, \code{mu > 0}
#' @param nu dispersion, \code{nu > 0}
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf
#' @return a numeric vector
#' @name cmp
NULL

#' @rdname cmp
#' @export
dcmp <- function(x, mu = 1, nu = 1, log = FALSE) {
  stopifnot(mu > 0, nu > 0)
  out <- cmp_logpmf_cpp(as.numeric(x), mu, nu)
  out[x < 0 | x != floor(x)] <- -Inf
  if (log) out else exp(out)
}
#' @rdname cmp
#' @export
pcmp <- function(q, mu = 1, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  hi <- max(0, floor(max(q)))
  cdf <- cumsum(dcmp(0:hi, mu, nu))
  out <- ifelse(q < 0, 0, pmin(cdf[pmax(floor(q), 0) + 1L], 1))
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
#' @rdname cmp
#' @export
qcmp <- function(p, mu = 1, nu = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e6) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  n <- max(32, ceiling(4 * mu))
  repeat {
    cdf <- cumsum(dcmp(0:n, mu, nu))
    if (max(cdf) >= max(p, na.rm = TRUE) - 1e-12 || n >= max.value) break
    n <- min(max.value, n * 4L)
  }
  vapply(p, function(pp) {
    if (is.na(pp)) return(NA_real_)
    i <- which(cdf >= pp - 1e-12)[1]
    if (is.na(i)) NA_real_ else i - 1
  }, 0)
}
#' @rdname cmp
#' @export
rcmp <- function(n, mu = 1, nu = 1) qcmp(stats::runif(n), mu, nu)

.iG_cmp <- function(v, mu, nu) qcmp(v, mu = mu, nu = nu)

.cmp1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for cmp.")
  cmp1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
         likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
         likdata$dupid, likdata$duplicate, likdata$offsets)
}
.cmp1.d12 <- function(pars, likdata) {
  cmp1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
          likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
          likdata$dupid, likdata$duplicate, likdata$offsets)
}
.cmp1fns <- list(d0 = .cmp1.d0, d120 = .cmp1.d12, d340 = NULL, m = 1, iG = .iG_cmp)
