## Compound Poisson-EGPD distribution
## Panjer recursion machinery and d/p/q/r functions

.compound_type <- function(type) {
  if (!is.numeric(type) || length(type) != 1L || !is.finite(type) ||
      !(type %in% 1:6))
    stop("'type' must be an integer in 1:6")
  as.integer(type)
}

## Distribution parameters required by `type` must be supplied: a missing one is
## an error (NA defaults hide missingness from the EGPD functions called here),
## while supplied NA or NaN values give NA or NaN results, as for the EGPD
## d/p/q/r functions. Returns that fill value, or NULL when all are non-missing.
.compound_required <- function(type, env = parent.frame()) {
  need <- c("lambda", "sigma", "xi",
            if (type %in% c(1L, 2L, 3L, 5L, 6L)) "kappa",
            if (type %in% c(4L, 5L, 6L)) "delta",
            if (type == 6L) "prob")
  miss <- need[vapply(need, function(nm) eval(call("missing", as.name(nm)), env),
                      logical(1))]
  if (length(miss) == 1L)
    stop(sprintf("Argument `%s' missing.", miss), call. = FALSE)
  if (length(miss) > 1L)
    stop(sprintf("Arguments %s missing.", paste0("`", miss, "'", collapse = ", ")),
         call. = FALSE)
  vals <- mget(need, envir = env)
  many <- names(vals)[lengths(vals) != 1L]
  if (length(many))
    stop("Compound distribution parameters must be scalar: ",
         paste(many, collapse = ", "), ".", call. = FALSE)
  na <- vapply(vals, function(v) length(v) == 1L && is.na(v), logical(1))
  if (!any(na)) return(NULL)
  if (any(vapply(vals[na], is.nan, logical(1)))) NaN else NA_real_
}

.compound_query_fill <- function(query, fill) {
  out <- rep(fill, length(query))
  out[is.nan(query)] <- NaN
  out
}

.validate_compound_params <- function(lambda, sigma, xi, kappa = NA,
                                      delta = NA, prob = NA, type = 1) {
  type <- .compound_type(type)
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda <= 0)
    stop("'lambda' must be positive and finite")
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) ||
      sigma <= 0)
    stop("'sigma' must be positive and finite")
  if (!is.numeric(xi) || length(xi) != 1L || !is.finite(xi))
    stop("'xi' must be finite")
  if (type %in% c(1L, 2L, 3L, 5L, 6L)) {
    if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
        kappa < 0 || (kappa == 0 && type != 2L))
      stop("'kappa' must be ", if (type == 2L) "non-negative" else "positive", " and finite")
  }
  if (type %in% c(4L, 5L, 6L)) {
    if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
        delta <= 0)
      stop("'delta' must be positive and finite")
  }
  if (type == 6L) {
    if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) ||
        prob < 0 || prob > 1)
      stop("'prob' must be finite and in [0, 1]")
  }
  type
}

.compound_grid_size <- function(K, values = numeric(), h = 1) {
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0)
    stop("'h' must be positive and finite")
  if (is.null(K))
    K <- max(ceiling(c(1, values[is.finite(values)]) / h)) + 100
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 1 || K != floor(K))
    stop("'K' must be a positive finite integer")
  K
}

## Relative tolerance for treating x / h as a lattice point. Gauge totals
## built by floating-point sums (e.g. rolling sums of 0.2 mm tips) carry
## absolute errors far above a few ulps, and ceiling() must not move them
## one bin up.
.compound_grid_tol <- sqrt(.Machine$double.eps)

.compound_lookup <- function(x, g, h = 1, cdf = FALSE, discrete = FALSE,
                             log_p = FALSE) {
  out <- rep(NA_real_, length(x))
  out[is.nan(x)] <- NaN
  zero_val <- if (!cdf && log_p) -Inf else 0
  out[!is.na(x) & x < 0] <- zero_val
  out[is.infinite(x) & x > 0] <- if (cdf) 1 else zero_val
  ok <- is.finite(x) & x >= 0
  if (!cdf && discrete) {
    out[ok & x != floor(x)] <- zero_val
    ok <- ok & x == floor(x)
  }
  scaled <- x[ok] / h
  if (!discrete) {
    nearest <- round(scaled)
    ## Only gauge totals have measurement-grid roundoff.
    near_pos_int <- is.finite(scaled) & nearest >= 1 & abs(scaled - nearest) <=
      .compound_grid_tol * pmax(1, abs(scaled))
    scaled[near_pos_int] <- nearest[near_pos_int]
  }
  bins <- if (cdf) floor(scaled) else if (discrete) round(scaled) else
    pmax(ceiling(scaled), as.numeric(x[ok] > 0))
  if (any(bins >= length(g)))
    stop("Query exceeds the compound search grid; increase 'K'.", call. = FALSE)
  out[ok] <- g[bins + 1L]
  out
}

.compound_grid_exceeded <- function()
  stop("Requested quantile exceeds the compound search grid; increase 'K'.",
       call. = FALSE)

## Evaluate or invert P(S > k), k = 1,...,K. By Poisson
## thinning, any jump > K already contributes 'overflow' to every survivor.
## b[n] is the mass at n with NO such jump, on the same (possibly conditional)
## scale as overflow. Through K these are exactly the original positive masses.
.compound_upper_grid <- function(b, log_f, lambda, overflow, target = NULL) {
  K <- length(b)
  within <- c(rev(cumsum(rev(b)))[-1L], 0)  # mass in (k,K]
  j <- seq_len(K)
  w <- j * exp(log(lambda) + log_f)
  mu <- sum(w)
  active <- which(w > 0)
  weights <- w[active]
  offset <- K - active
  invert <- !is.null(target)
  if (invert) target <- target - overflow
  L <- as.double(K)
  repeat {
    ## For the remaining small-jump tail R_L = sum_{n>L} b[n], sum the
    ## Panjer identity n b[n] = sum_j w[j] b[n-j] over n > L. This gives
    ## (L+1-mu) R_L <= sum_j w[j] sum_{i=L-j+1}^L b[i]. All terms in
    ## the bound and in the survivor are non-negative, even deep in the tail.
    bound <- if (L + 1 > mu) sum(w * cumsum(rev(b))) / (L + 1 - mu) else Inf
    if (invert) {
      if (any(target < 0)) .compound_grid_exceeded()
      k <- vapply(target, function(t) which(within <= t)[1L], integer(1))
      if (all(bound <= target - within[k])) return(k)
    } else if (bound <= .Machine$double.eps / 16) {
      ## The omitted mass is below 1/8 ulp of any CDF in [1/2,1].
      return(overflow + within)
    }

    ## These extra masses refine P(S > K), not the search grid.
    block <- c(b, numeric(K))
    for (i in j)
      block[K + i] <- sum(weights * block[i + offset]) / (L + i)
    b <- block[K + j]
    if (invert) target <- target - sum(b) else overflow <- overflow + sum(b)
    L <- L + K
  }
}

.compound_quantiles <- function(p, lambda, sigma, xi, kappa, delta, prob,
                                type, K, h = NULL, positive = FALSE) {
  grid <- .compound_grid(lambda, sigma, xi, kappa, delta, prob, type, K, h,
                         log_out = positive, upper = TRUE)
  log_norm <- if (positive) .compound_log_nonzero(grid$log_rate) else 0
  if (!is.finite(log_norm)) stop("Positive aggregate probability is numerically zero")
  g <- if (positive) c(0, exp(grid$g[-1L] - log_norm)) else grid$g
  cg <- pmin(1, cumsum(g))
  unit <- if (is.null(h)) 1 else h
  ## Check support from the parameters; every positive grid mass may underflow.
  unbounded <- !is.null(h) || xi >= 0 || sigma > -xi
  out <- vapply(p, function(pi) {
    if (is.na(pi)) return(as.numeric(pi))
    if (!is.finite(pi) || pi < 0 || pi > 1) {
      warning("NaNs produced", call. = FALSE)
      return(NaN)
    }
    if (pi == 0) return(if (positive) unit else 0)
    if (pi == 1) return(if (unbounded) Inf else 0)
    if (pi > 0.5) return(0)  # handle on the survivor scale below
    idx <- which(cg >= pi)[1]
    if (is.na(idx)) .compound_grid_exceeded()
    (idx - 1) * unit
  }, numeric(1))

  upper <- which(p > 0.5 & p < 1)
  if (!positive)  # targets reached at the zero atom need no positive bin
    upper <- upper[(1 - p[upper]) < -expm1(grid$log_g0)]
  if (length(upper)) {
    if (cg[K + 1L] < 0.25) .compound_grid_exceeded()  # well away from every target
    overflow <- exp(grid$log_overflow - log_norm)
    out[upper] <- unit * .compound_upper_grid(g[-1L], grid$log_f[-1L], lambda,
                                              overflow, target = 1 - p[upper])
  }
  out
}

#' Panjer recursion for Poisson case (a=0, b=lambda) with dynamic log-rescaling
#' @param lambda Poisson rate parameter
#' @param f discretized severity PMF vector (f[1] = P(X = 0), f[j+1] = P(X in ((j-1)h, jh]))
#' @param log_out logical; if TRUE, return log-probabilities directly
#' @param log_f original log severity masses, including finite values below the
#'   linear probability range
#' @param log_g0 zero atom on the log scale
#' @return vector of compound distribution (log-)PMF values
#' @noRd
.panjer_poisson <- function(lambda, f, log_out = FALSE, log_f = log(f),
                            log_g0 = -lambda * (1 - f[1L])) {
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
    stop("'lambda' must be positive and finite")
  K <- length(f)
  log_g <- rep(-Inf, K)
  logscale <- log_g0
  log_g[1L] <- logscale
  if (K == 1L) {
    return(if (log_out) log_g else exp(log_g))
  }

  pos_f <- which(is.finite(log_f[-1L]))
  if (length(pos_f) == 0L) {
    return(if (log_out) log_g else exp(log_g))
  }
  log_jf <- log(pos_f) + log_f[pos_f + 1L]
  jf <- seq_len(K - 1L) * f[-1L]
  ## From the first subnormal jump on, a linear convolution need not be
  ## accurate, even when its sum is finite. Log masses distinguish this
  ## case from a genuinely impossible jump (log f_j = -Inf).
  log_from <- match(TRUE, log_f[pos_f + 1L] < log(.Machine$double.xmin))
  log_from <- if (is.na(log_from)) Inf else pos_f[log_from] + 1L

  g_tilde <- numeric(K)
  g_tilde[1L] <- 1

  for (k in 2:K) {
    idx <- seq_len(k - 1L)
    val <- if (k < log_from) (lambda / (k - 1L)) * sum(jf[idx] * g_tilde[k - idx]) else 0
    if (is.finite(val) && val >= 1e-250) {
      g_tilde[k] <- val
      log_g[k] <- logscale + log(val)
      if (val > 1e200) {
        shift <- log(val)
        logscale <- logscale + shift
        g_tilde[1:k] <- g_tilde[1:k] / val
      }
    } else if (log_out || k >= log_from || logscale < -50 || !is.finite(val)) {
      jj <- pos_f[pos_f <= k - 1L]
      if (length(jj) > 0L) {
        terms <- log_jf[seq_along(jj)] + log_g[k - jj]
        m <- max(terms)
        if (is.finite(m)) {
          lgk <- log(lambda) - log(k - 1L) + m + log(sum(exp(terms - m)))
          log_g[k] <- lgk
          if (lgk - logscale > 460) {
            logscale <- lgk
            g_tilde[1:k] <- exp(log_g[1:k] - logscale)
          } else {
            g_tilde[k] <- exp(lgk - logscale)
          }
        }
      }
    } else if (val > 0) {
      g_tilde[k] <- val
      log_g[k] <- logscale + log(val)
    }
  }

  if (log_out) log_g else exp(log_g)
}


#' Discretize EGPD severity on grid {0, h, 2h, ..., K*h}
#'
#' Uses upper-bin discretization per Ailliot et al. (2026) Eq. (F.2):
#'   f_0 = 0
#'   f_j = F(j * h) - F((j - 1) * h) for j = 1, ..., K
#' Preserves the exact continuous compound zero atom P(S = 0) = exp(-lambda).
#'
#' @param sigma GPD scale
#' @param xi GPD shape
#' @param kappa G-transformation parameter
#' @param delta G-transformation parameter
#' @param prob G-transformation parameter
#' @param type integer 1-6
#' @param h bin width
#' @param K number of bins
#' @return numeric vector of length K+1 with log bin probabilities
#' @noRd
.gauge_logpmf <- function(sigma, xi, kappa, delta, prob, type, h, K) {
  breaks <- seq_len(K) * h
  lF <- pegpd(breaks, prob = prob, kappa = kappa, delta = delta,
              sigma = sigma, xi = xi, type = type, log.p = TRUE)
  lf <- lF + .log1mexp(c(-Inf, lF[-K]) - lF)
  lf[lF == -Inf] <- -Inf
  up <- which(lF > -log(2))
  if (length(up)) {
    j <- up[1L]:K
    lS <- pegpd(c(j[1L] - 1L, j) * h, prob = prob, kappa = kappa, delta = delta,
                sigma = sigma, xi = xi, type = type, lower.tail = FALSE, log.p = TRUE)
    lo <- lS[-length(lS)]
    hi <- lS[-1L]
    lf[j] <- ifelse(lo == -Inf, -Inf, lo + .log1mexp(hi - lo))
  }
  c(-Inf, lf)
}

## log P(Poisson(r) > 0), including rates below the linear range.
.compound_log_nonzero <- function(log_rate) {
  if (log_rate < log(.Machine$double.eps)) log_rate else log(-expm1(-exp(log_rate)))
}

.compound_grid <- function(lambda, sigma, xi, kappa, delta, prob, type, K,
                           h = NULL, log_out = FALSE, upper = FALSE) {
  if (is.null(h)) {
    lf <- ddiscegpd(0:K, prob = prob, kappa = kappa, delta = delta,
                    sigma = sigma, xi = xi, type = type, log = TRUE)
    ls <- pegpd(1, prob = prob, kappa = kappa, delta = delta,
                sigma = sigma, xi = xi, type = type, lower.tail = FALSE, log.p = TRUE)
    lg0 <- if (ls > log(.Machine$double.xmin)) -lambda * exp(ls) else -exp(log(lambda) + ls)
  } else {
    lf <- .gauge_logpmf(sigma, xi, kappa, delta, prob, type, h, K)
    ls <- 0
    lg0 <- -lambda
  }
  log_overflow <- if (upper) {
    ## J = ceiling(X/h) for the gauge, floor(X) for integer severity.
    log_big <- pegpd(if (is.null(h)) K + 1 else K * h, prob = prob, kappa = kappa,
                     delta = delta, sigma = sigma, xi = xi, type = type,
                     lower.tail = FALSE, log.p = TRUE)
    .compound_log_nonzero(log(lambda) + log_big)
  }
  list(g = .panjer_poisson(lambda, exp(lf), log_out = log_out, log_f = lf, log_g0 = lg0),
       log_f = lf, log_g0 = lg0, log_rate = log(lambda) + ls, log_overflow = log_overflow)
}

.compound_pmf_grid <- function(...) .compound_grid(...)$g

.compound_cdf_grid <- function(lambda, ...) {
  grid <- .compound_grid(lambda, ..., upper = TRUE)
  cg <- cumsum(grid$g)
  high <- which(cg[-1L] > 0.5)
  if (length(high)) {
    ## Forward mass sums inherit the error in the whole law near one.
    ## Sum the other tail, including mass beyond K (do not normalize the grid).
    upper <- .compound_upper_grid(grid$g[-1L], grid$log_f[-1L], lambda,
                                  exp(grid$log_overflow))
    cg[high + 1L] <- 1 - upper[high]
  }
  cummax(cg)
}


#' Density (PMF) of Compound Poisson-EGPD
#'
#' Lattice probabilities of the gauge-discretised compound Poisson-EGPD law
#' of Ailliot, Gaetan and Naveau (2026, Definition 2.5 and Appendix F).
#'
#' The random variable is a compound Poisson sum of EGPD event amounts, each
#' rounded up to the gauge resolution \eqn{h}:
#' \deqn{S_h = \sum_{i=1}^{N} h \lceil X_i / h \rceil,}
#' where \eqn{N \sim \mathrm{Poisson}(\lambda)} and the event amounts
#' \eqn{X_i} are i.i.d. \eqn{\mathrm{EGPD}(\sigma, \xi, \kappa, \ldots)},
#' independent of \eqn{N}. \eqn{S_h} lives on the lattice \eqn{\{0, h, 2h, \ldots\}}. Each event
#' contributes \eqn{j \ge 1} bins with probability
#' \eqn{f_j = F(jh) - F((j-1)h)} and \eqn{f_0 = 0} (Eq. (F.2)), so the zero
#' atom is exactly \eqn{P(S_h = 0) = \exp(-\lambda)}, and the lattice
#' probabilities follow from the Panjer recursion (Eq. (F.3)). They are exact
#' for \eqn{S_h} up to truncation of the severity at \code{K} bins, and
#' \code{rcpegpd()} with the same \code{h} simulates the same law.
#'
#' This is the model for gauge data recorded in multiples of \code{h} (for
#' example 0.2 mm tipping buckets). For a continuous compound sum
#' \eqn{S = \sum X_i}, \eqn{S_h \ge S} and \eqn{E[S_h] - E[S]} is about
#' \eqn{\lambda h / 2}; use an \code{h} much smaller than typical event
#' amounts, and \code{rcpegpd(..., h = 0)} to simulate \eqn{S} itself.
#'
#' \code{dcpegpd()} returns lattice probabilities, not densities: a positive
#' \eqn{x} is mapped to the bin \eqn{k = \lceil x / h \rceil}, i.e.
#' \eqn{P(S_h = kh)}, after snapping values within a relative
#' \code{sqrt(.Machine$double.eps)} of a lattice point onto it (so
#' floating-point sums of gauge data do not move up one bin). Negative
#' inputs have mass zero, and \code{pcpegpd(q)} is \eqn{P(S_h \le q)}.
#' Explicit grids never extrapolate past \code{K}.
#'
#' The parameters required by \code{type} (\code{lambda}, \code{sigma},
#' \code{xi} and the carrier parameters \code{kappa}, \code{delta},
#' \code{prob} that the type uses) must be supplied: a missing one is an
#' error, whereas a supplied \code{NA} (\code{NaN}) gives \code{NA}
#' (\code{NaN}) results, as for \code{\link{pegpd}}. Missing queries also
#' preserve \code{NA} versus \code{NaN}. Other invalid parameters are errors.
#'
#' @references Ailliot, P., Gaetan, C. and Naveau, P. (2026). A parsimonious
#'   tail compliant multiscale statistical model for aggregated rainfall.
#'   \emph{Advances in Water Resources}, 208, 105216.
#'   \doi{10.1016/j.advwatres.2026.105216}. Definition, equation and appendix
#'   numbers on this page refer to the original-author preprint
#'   arXiv:2601.08350v1.
#'
#' @param x numeric vector of quantiles
#' @param lambda Scalar Poisson rate (> 0). The distribution parameters
#'   describe one common severity and rate; they are not recycled over queries.
#' @param prob G-transformation parameter (type 6 only)
#' @param kappa G-transformation parameter; for type 2, zero is the
#'   uniform-carrier limit (ordinary GPD severity before gauge rounding).
#' @param delta G-transformation parameter
#' @param sigma GPD scale parameter (> 0)
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G-transformation type
#' @param h positive gauge resolution: event amounts are rounded up to
#'   multiples of \code{h} (default 0.2).
#' @param K number of bins (default: covers finite queries with a grid margin).
#'   Quantile functions default to 5000 bins.
#' @param log logical; if TRUE, return log-probabilities
#'
#' @return numeric vector of (log-)probabilities
#'   \eqn{P(S_h = h\lceil x/h \rceil)}
#' @seealso \code{\link{pcpegpd}}, \code{\link{qcpegpd}},
#'   \code{\link{rcpegpd}}, \code{\link{fit_aggregate_egpd}}
#' @export
dcpegpd <- function(x, lambda, prob = NA, kappa = NA, delta = NA,
                     sigma, xi, type = 1, h = 0.2, K = NULL,
                     log = FALSE) {
  type <- .compound_type(type)
  fill <- .compound_required(type)
  K <- .compound_grid_size(K, x, h)
  if (!is.null(fill)) return(.compound_query_fill(x, fill))
  type <- .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                                    kappa = kappa, delta = delta, prob = prob,
                                    type = type)

  g <- .compound_pmf_grid(lambda, sigma, xi, kappa, delta, prob, type, K, h, isTRUE(log))

  .compound_lookup(x, g, h, log_p = isTRUE(log))
}


#' Distribution function of Compound Poisson-EGPD
#'
#' \eqn{P(S_h \le q)} for the gauge-discretised compound sum
#' \eqn{S_h = \sum_{i=1}^{N} h \lceil X_i / h \rceil} described in
#' \code{\link{dcpegpd}} (Ailliot et al. 2026, Appendix F).
#'
#' @inheritParams dcpegpd
#' @param q numeric vector of quantiles
#'
#' @return numeric vector of cumulative probabilities
#' @export
pcpegpd <- function(q, lambda, prob = NA, kappa = NA, delta = NA,
                     sigma, xi, type = 1, h = 0.2, K = NULL) {
  type <- .compound_type(type)
  fill <- .compound_required(type)
  K <- .compound_grid_size(K, q, h)
  if (!is.null(fill)) return(.compound_query_fill(q, fill))
  type <- .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                                    kappa = kappa, delta = delta, prob = prob,
                                    type = type)

  g <- .compound_cdf_grid(lambda, sigma, xi, kappa, delta, prob, type, K, h)
  .compound_lookup(q, g, h, cdf = TRUE)
}


#' Quantile function of Compound Poisson-EGPD
#'
#' Quantiles of the gauge-discretised compound sum
#' \eqn{S_h = \sum_{i=1}^{N} h \lceil X_i / h \rceil} described in
#' \code{\link{dcpegpd}} (Ailliot et al. 2026, Appendix F).
#'
#' @inheritParams dcpegpd
#' @param p numeric vector of probabilities
#' @details Returns the smallest lattice point \eqn{kh} with
#'   \eqn{P(S_h \le kh) \ge p}, and \code{Inf} at probability one for a
#'   nondegenerate compound distribution. Invalid probabilities produce
#'   \code{NaN}. If a probability in \eqn{(0,1)} exceeds the grid CDF, an
#'   error requests a larger \code{K}; the endpoint is never returned as a
#'   false quantile.
#'
#' @return numeric vector of quantiles (multiples of \code{h})
#' @export
qcpegpd <- function(p, lambda, prob = NA, kappa = NA, delta = NA,
                     sigma, xi, type = 1, h = 0.2, K = NULL) {
  type <- .compound_type(type)
  fill <- .compound_required(type)
  K <- .compound_grid_size(if (is.null(K)) 5000 else K, h = h)
  if (!is.null(fill)) return(.compound_query_fill(p, fill))
  type <- .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                                    kappa = kappa, delta = delta, prob = prob,
                                    type = type)

  .compound_quantiles(p, lambda, sigma, xi, kappa, delta, prob, type, K, h)
}


#' Random generation from Compound Poisson-EGPD
#'
#' Draws \eqn{N \sim \mathrm{Poisson}(\lambda)} and \eqn{N} EGPD event
#' amounts, rounds each amount up to a multiple of \code{h} and returns the
#' sum: a draw of \eqn{S_h = \sum_{i=1}^{N} h \lceil X_i / h \rceil}, which
#' has exactly the law of \code{\link{dcpegpd}}, \code{\link{pcpegpd}} and
#' \code{\link{qcpegpd}} for the same \code{h} (Ailliot et al. 2026,
#' Appendix F). Such draws can be passed directly to
#' \code{distfit(family = "cpegpd")} or \code{\link{fit_aggregate_egpd}}
#' with matching \code{cpegpd.h}.
#'
#' \code{h = 0} returns the continuous compound sum
#' \eqn{S = \sum_{i=1}^{N} X_i}, whose law d/p/qcpegpd approximate only as
#' \eqn{h \to 0} (\eqn{S_h \ge S}, with \eqn{E[S_h] - E[S]} about
#' \eqn{\lambda h / 2}).
#'
#' Missing required parameters are an error and supplied \code{NA}
#' parameters give \code{NA} draws with a warning (see \code{\link{dcpegpd}}).
#'
#' @inheritParams dcpegpd
#' @param n number of observations
#' @param h gauge resolution (default 0.2, as in \code{dcpegpd()}); each
#'   event amount is rounded up to a multiple of \code{h}. Use \code{h = 0}
#'   for the continuous sum.
#'
#' @return numeric vector of length \code{n}; multiples of \code{h} when
#'   \code{h > 0}
#' @export
rcpegpd <- function(n, lambda, prob = NA, kappa = NA, delta = NA,
                     sigma, xi, type = 1, h = 0.2) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n))
    stop("'n' must be a non-negative integer")
  type <- .compound_type(type)
  fill <- .compound_required(type)
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h < 0)
    stop("'h' must be zero or positive and finite")
  n <- as.integer(n)
  if (n == 0L) return(numeric(0))
  if (!is.null(fill)) {
    warning("NAs produced", call. = FALSE)
    return(rep(fill, n))
  }
  type <- .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                                    kappa = kappa, delta = delta, prob = prob,
                                    type = type)
  N <- rpois(n, lambda)
  vapply(N, function(ni) {
    if (ni == 0) return(0)
    x <- regpd(ni, prob = prob, kappa = kappa, delta = delta,
               sigma = sigma, xi = xi, type = type)
    if (h > 0) {
      ## X > 0 almost surely, so every event occupies at least bin 1 (f_0 = 0)
      ## even if a tiny amount underflows to zero.
      h * sum(pmax(ceiling(x / h), 1))
    } else {
      sum(x)
    }
  }, numeric(1))
}


############################################
## Compound Poisson-Discrete EGPD        ##
## (exact integer-valued compound sum)    ##
############################################

#' Density (PMF) of Compound Poisson-Discrete EGPD
#'
#' Computes the probability mass function of the Compound Poisson-Discrete
#' EGPD distribution using Panjer recursion on the exact discrete EGPD
#' severity PMF (no discretization step).
#'
#' The model is \eqn{S = X_1 + \cdots + X_N} where
#' \eqn{N \sim \mathrm{Poisson}(\lambda)} and
#' \eqn{X_i \sim \mathrm{Discrete\mbox{-}EGPD}(\sigma, \xi, \kappa, \ldots)}.
#' Since the severity is already integer-valued, the Panjer recursion
#' produces an exact compound distribution with no discretization error.
#' Negative and noninteger inputs have mass zero; negative CDF queries
#' return zero. Finite queries outside an explicit grid fail rather than
#' being clipped to its endpoint. Parameters required by \code{type} must be
#' supplied (a missing one is an error), whereas supplied \code{NA}
#' (\code{NaN}) values give \code{NA} (\code{NaN}) results. Missing queries
#' also preserve \code{NA} versus \code{NaN}. Other invalid parameters are errors.
#'
#' @param x numeric vector of non-negative integer quantiles
#' @param lambda Scalar Poisson rate (> 0). The distribution parameters
#'   describe one common severity and rate; they are not recycled over queries.
#' @param prob G-transformation parameter (type 6 only)
#' @param kappa G-transformation parameter; for type 2, zero is the
#'   uniform-carrier limit (floored GPD severity).
#' @param delta G-transformation parameter
#' @param sigma GPD scale parameter (> 0)
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G-transformation type
#' @param K grid size for Panjer recursion (default: covers finite queries
#'   with a grid margin; the quantile function defaults to 5000).
#' @param log logical; if TRUE, return log-probabilities
#'
#' @return numeric vector of (log-)probabilities
#' @export
#' @seealso \code{\link{pcpdegpd}}, \code{\link{qcpdegpd}},
#'   \code{\link{rcpdegpd}}, \code{\link{fitegpd}}
dcpdegpd <- function(x, lambda, prob = NA, kappa = NA, delta = NA,
                      sigma, xi, type = 1, K = NULL, log = FALSE) {
  type <- .compound_type(type)
  fill <- .compound_required(type)
  K <- .compound_grid_size(K, x)
  if (!is.null(fill)) return(.compound_query_fill(x, fill))
  type <- .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                                    kappa = kappa, delta = delta, prob = prob,
                                    type = type)

  g <- .compound_pmf_grid(lambda, sigma, xi, kappa, delta, prob, type, K, log_out = isTRUE(log))

  ## Map x values to g (R is 1-indexed: g[1] = P(S=0))
  .compound_lookup(x, g, discrete = TRUE, log_p = isTRUE(log))
}


#' Distribution function of Compound Poisson-Discrete EGPD
#'
#' @inheritParams dcpdegpd
#' @param q numeric vector of quantiles
#'
#' @return numeric vector of cumulative probabilities
#' @export
#' @seealso \code{\link{dcpdegpd}}, \code{\link{qcpdegpd}},
#'   \code{\link{rcpdegpd}}
pcpdegpd <- function(q, lambda, prob = NA, kappa = NA, delta = NA,
                      sigma, xi, type = 1, K = NULL) {
  type <- .compound_type(type)
  fill <- .compound_required(type)
  K <- .compound_grid_size(K, q)
  if (!is.null(fill)) return(.compound_query_fill(q, fill))
  type <- .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                                    kappa = kappa, delta = delta, prob = prob,
                                    type = type)

  g <- .compound_cdf_grid(lambda, sigma, xi, kappa, delta, prob, type, K)
  .compound_lookup(q, g, cdf = TRUE, discrete = TRUE)
}


#' Quantile function of Compound Poisson-Discrete EGPD
#'
#' @inheritParams dcpdegpd
#' @param p numeric vector of probabilities
#' @details Invalid probabilities produce \code{NaN}; probability one
#'   returns \code{Inf} for a nondegenerate compound distribution.
#'   Unattained interior probabilities cause an error requesting larger
#'   \code{K}, rather than a misleading grid endpoint.
#'
#' @return numeric vector of quantiles (non-negative integers)
#' @export
#' @seealso \code{\link{dcpdegpd}}, \code{\link{pcpdegpd}},
#'   \code{\link{rcpdegpd}}
qcpdegpd <- function(p, lambda, prob = NA, kappa = NA, delta = NA,
                      sigma, xi, type = 1, K = NULL) {
  type <- .compound_type(type)
  fill <- .compound_required(type)
  K <- .compound_grid_size(if (is.null(K)) 5000 else K)
  if (!is.null(fill)) return(.compound_query_fill(p, fill))
  type <- .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                                    kappa = kappa, delta = delta, prob = prob,
                                    type = type)

  .compound_quantiles(p, lambda, sigma, xi, kappa, delta, prob, type, K)
}


#' Random generation from Compound Poisson-Discrete EGPD
#'
#' Direct simulation: draw \eqn{N \sim \mathrm{Poisson}(\lambda)}, then
#' sum \eqn{N} i.i.d. Discrete-EGPD draws. Missing required parameters are an
#' error; supplied \code{NA} parameters give \code{NA} draws with a warning.
#'
#' @inheritParams dcpdegpd
#' @param n number of observations
#'
#' @return integer vector of length \code{n}
#' @export
#' @seealso \code{\link{dcpdegpd}}, \code{\link{pcpdegpd}},
#'   \code{\link{qcpdegpd}}
rcpdegpd <- function(n, lambda, prob = NA, kappa = NA, delta = NA,
                      sigma, xi, type = 1) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n))
    stop("'n' must be a non-negative integer")
  type <- .compound_type(type)
  fill <- .compound_required(type)
  n <- as.integer(n)
  if (n == 0L) return(integer(0))
  if (!is.null(fill)) {
    warning("NAs produced", call. = FALSE)
    return(rep(NA_integer_, n))
  }
  type <- .validate_compound_params(lambda = lambda, sigma = sigma, xi = xi,
                                    kappa = kappa, delta = delta, prob = prob,
                                    type = type)
  N <- rpois(n, lambda)
  vapply(N, function(ni) {
    if (ni == 0L) return(0L)
    s <- sum(as.numeric(rdiscegpd(ni, prob = prob, kappa = kappa,
                                  delta = delta, sigma = sigma,
                                  xi = xi, type = type)))
    if (!is.finite(s) || s > .Machine$integer.max)
      stop("Simulated compound counts exceed R's 32-bit integer range.")
    as.integer(s)
  }, integer(1))
}
