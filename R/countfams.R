## Generalized Poisson, generalized Waring and Poisson-lognormal: d/p/q/r functions and
## the likelihood bundles for the native GAM fitter. The log pmfs live in
## src/countfams.cpp; see its header for the parameterisations, links and tail behaviour.
##
## The R names avoid masking stats::plnorm (the lognormal cdf): the Poisson-lognormal
## functions are d/p/q/rpoislnorm, while the egpd family string remains "plnorm".

## ---------------------------------------------------------------- generalized Poisson

#' Generalized Poisson distribution
#'
#' Consul & Jain's generalized Poisson in the mean parameterisation:
#' \code{theta = mu (1 - lambda)}, so \code{E[Y] = mu} and
#' \code{Var[Y] = mu / (1 - lambda)^2}. The tail decays geometrically, so this family
#' captures overdispersion but not a heavy tail.
#'
#' @param x,q vector of counts
#' @param p vector of probabilities
#' @param n number of draws; if \code{length(n) != 1}, the length is taken to be
#'   the number required (\code{numeric(0)} gives none), as in \code{\link[stats]{rpois}}
#' @param mu mean, \code{mu > 0}
#' @param lambda dispersion, \code{lambda} in \code{[0, 1)}; \code{lambda = 0} is Poisson
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf; nonattaining
#'   caps return \code{NA} with a warning. Invalid probabilities return \code{NaN},
#'   with endpoints 0 and \code{Inf}.
#' @return a numeric vector
#' @name gpois
NULL

#' @rdname gpois
#' @export
dgpois <- function(x, mu = 1, lambda = 0.5, log = FALSE) {
  .count_d(.gpois_fam, x, list(mu = mu, lambda = lambda), log = log)
}
#' @rdname gpois
#' @export
pgpois <- function(q, mu = 1, lambda = 0.5, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.gpois_fam, q, list(mu = mu, lambda = lambda), lower.tail, log.p)
}
#' @rdname gpois
#' @export
qgpois <- function(p, mu = 1, lambda = 0.5, lower.tail = TRUE, log.p = FALSE,
                   max.value = 1e6) {
  .count_q(.gpois_fam, p, list(mu = mu, lambda = lambda), lower.tail, log.p, max.value)
}
#' @rdname gpois
#' @export
rgpois <- function(n, mu = 1, lambda = 0.5) {
  .count_r(.gpois_fam, n, list(mu = mu, lambda = lambda), max.value = 1e6)
}

.gpois_fam <- list(
  pars = c("mu", "lambda"),
  valid = function(p) is.finite(p$mu) & p$mu > 0 &
    is.finite(p$lambda) & p$lambda >= 0 & p$lambda < 1,
  lpmf = function(n, p) gpois_logpmf_cpp(0:n, p$mu, p$lambda),
  lpmf_at = function(x, p) gpois_logpmf_cpp(as.numeric(x), p$mu, p$lambda),
  rlim = function(p) p$lambda * exp(1 - p$lambda),
  positive_tail = function(n, p) {
    theta <- p$mu * (1 - p$lambda)
    if (theta > .01 || p$lambda < .95) return(NULL)
    ## Divide by theta before subtracting: p_j/theta = Pois(j-1;theta+j*lambda)/j.
    ## S(0)/theta tends to one even if the atom itself rounds to one.
    j <- seq_len(n)
    total <- .expm1_over_x(theta)
    rem <- total - c(0, cumsum(stats::dpois(j-1, theta + j*p$lambda)/j))
    rem[rem < 1e-4*total] <- NA_real_  # use the direct tail for ill-conditioned differences
    log(p$mu) + log1p(-p$lambda) + log(rem)
  },
  init = function(p) if (p$mu*(1-p$lambda) < .01) 32 else
    max(32, ceiling(4 * p$mu / (1 - p$lambda))))

## ---------------------------------------------------------------- generalized Waring

#' Generalized Waring distribution
#'
#' The three-parameter generalized Waring distribution in the mean parameterisation, with
#' \code{a = mu (rho - 1) / k}, so \code{E[Y] = mu}. Its tail is regularly varying,
#' \code{P(Y = y) ~ y^-(rho + 1)}: a power law with extreme-value index \code{xi = 1/rho},
#' placing it in the Frechet domain, unlike the negative binomial, PIG, generalized Poisson
#' or Poisson-lognormal. The mean exists for \code{rho > 1} and the variance for
#' \code{rho > 2}. Setting \code{k = 1} recovers the Waring distribution.
#'
#' As \code{rho -> Inf} with \code{mu} and \code{k} fixed it tends to the negative binomial
#' with \code{size = k}, and as \code{rho} and \code{k} grow with \code{rho >> k} to the
#' Poisson; these functions stay accurate along the way. When fitted by
#' \code{egpd(family = "gwaring")} the likelihood is evaluated with the parameters clamped to
#' \code{mu} in \code{[1e-10, 1e12]}, \code{k} in \code{[1e-8, 1e15]} and \code{rho} in
#' \code{[1 + 1e-6, 1e15]}. At \code{rho = 1e15} the distribution differs from its
#' negative-binomial limit only by \code{O(1/rho)}, and the Poisson limit is approached most
#' closely at \code{k = sqrt(mu rho)}, inside the bound on \code{k}. Beyond a bound the
#' likelihood is flat (a fit to light-tailed data can drift there), and \code{predict()}
#' reports the clamped values.
#'
#' Tail probabilities beyond the run of the probability mass function that the functions sum
#' come from a beta-negative-binomial mixture integral. Against 512-bit references it is
#' accurate to about 1e-11 (relative) throughout and beyond the fitting bounds. Its one known
#' limitation is a shape \code{k > 1e15}, beyond the fitting bound, with the tail evaluated at
#' a count above 1000: \code{P(Y > q)} can then be off by up to about 1e-3 (relative), though by
#' less than 1e-9 for counts up to 2e4, and \code{P(Y <= q)} by up to 1e-3 times
#' \code{P(Y > q)}. When a call uses such a value, \code{pgwaring()}, \code{qgwaring()} and
#' \code{rgwaring()} return it unchanged and warn once, with a warning of class
#' \code{"egpd_gwaring_accuracy"}. Calls that use only lower-tail probabilities on the
#' probability scale (\code{lower.tail = TRUE, log.p = FALSE}, and \code{rgwaring()}, which
#' inverts them) warn only where that \code{P(Y > q)} exceeds 1e-8, since below it they are
#' accurate to about 1e-11.
#' \code{suppressWarnings(expr, classes = "egpd_gwaring_accuracy")} silences the warning.
#' At a CDF jump, \code{qgwaring()} allows for numerical rounding: targets
#' within the comparison tolerance need not give the strict inverse of an
#' arbitrary-precision CDF. An ordinary upper probability above one half
#' is inverted using the smaller lower probability \code{1 - p}.
#'
#' @param x,q vector of counts
#' @param p vector of probabilities
#' @param n number of draws; if \code{length(n) != 1}, the length is taken to be
#'   the number required (\code{numeric(0)} gives none), as in \code{\link[stats]{rpois}}
#' @param mu mean, \code{mu > 0}
#' @param k shape, \code{k > 0}
#' @param rho tail parameter, \code{rho > 1}
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf; nonattaining
#'   caps return \code{NA} with a warning. Invalid probabilities return \code{NaN},
#'   with endpoints 0 and \code{Inf}.
#' @return a numeric vector
#' @name gwaring
NULL

#' @rdname gwaring
#' @export
dgwaring <- function(x, mu = 1, k = 1, rho = 2.5, log = FALSE) {
  .count_d(.gwaring_fam, x, list(mu = mu, k = k, rho = rho), log = log)
}
#' @rdname gwaring
#' @export
pgwaring <- function(q, mu = 1, k = 1, rho = 2.5, lower.tail = TRUE, log.p = FALSE) {
  .gwaring_watch(.count_p(.gwaring_fam, q, list(mu = mu, k = k, rho = rho), lower.tail, log.p),
                 sys.call(), isTRUE(as.logical(lower.tail)) && isFALSE(as.logical(log.p)))
}
#' @rdname gwaring
#' @export
qgwaring <- function(p, mu = 1, k = 1, rho = 2.5, lower.tail = TRUE, log.p = FALSE,
                     max.value = 1e7) {
  .gwaring_watch(.count_q(.gwaring_fam, p, list(mu = mu, k = k, rho = rho), lower.tail, log.p,
                          max.value), sys.call(),
                 isTRUE(as.logical(lower.tail)) && isFALSE(as.logical(log.p)))
}
#' @rdname gwaring
#' @export
rgwaring <- function(n, mu = 1, k = 1, rho = 2.5) {
  .gwaring_watch(.count_r(.gwaring_fam, n, list(mu = mu, k = k, rho = rho), max.value = 1e7),
                 sys.call(), TRUE)
}

## The power-law tail is not summable in practice, so the remainder beyond the ladder is
## the exact beta-negative-binomial mixture integral (src/countfams.cpp, gw_lcdf).
.gwaring_fam <- list(
  pars = c("mu", "k", "rho"),
  valid = function(p) is.finite(p$mu) & p$mu > 0 & is.finite(p$k) & p$k > 0 &
    is.finite(p$rho) & p$rho > 1,
  lpmf = function(n, p) gwaring_logpmf_cpp(0:n, p$mu, p$k, p$rho),
  lpmf_at = function(x, p) gwaring_logpmf_cpp(as.numeric(x), p$mu, p$k, p$rho),
  lsurv = function(k, p) {
    ls <- gwaring_lcdf_cpp(as.numeric(k), p$mu, p$k, p$rho, FALSE)
    if (p$k > .GW_ACC_K && any(k > .GW_ACC_Q)) .gwaring_flag(p$k, k, ls)
    ls
  },
  init = function(p) max(32, ceiling(4 * p$mu)))

## Accuracy limit of that integral, mapped against 512-bit references: for k > 1e15 (beyond the
## fitting bound) and an evaluation point q > 1000, its quadrature can miss the sharp drop of
## P(NB(k, w) > q) where k (1 - w)/w ~ q, and P(Y > q) was off by up to ~1e-3 (relative; below
## 1e-9 for q <= 2e4). Elsewhere it is exact to ~1e-11. The engine evaluates it only beyond the
## pmf ladder, at one count `need` with P(Y > need) < 1e-3, and each P(Y > q) it reports there
## is an exact ladder sum plus P(Y > need). Upper-tail and log-scale values thus inherit the
## relative error, while P(Y <= q) >= 0.999 is off by at most ~1e-3 P(Y > need): below 1e-11
## when P(Y > need) <= 1e-8. lsurv records such evaluations (largest count, largest log S), and
## the p/q/r function whose call used one warns once (class egpd_gwaring_accuracy), unless it
## uses only lower-tail probabilities on the probability scale and every log S is <= log(1e-8).
## Values are unchanged.
.GW_ACC_K <- 1e15
.GW_ACC_Q <- 1000
.GW_ACC_LS <- log(1e-8)
.gwaring_acc <- new.env(parent = emptyenv())
.gwaring_flag <- function(k, q, ls) {
  big <- which(q > .GW_ACC_Q)
  qm <- max(q[big])
  lm <- max(ls[big])                                 # NA/NaN propagate: then always warn
  h <- .gwaring_acc$hit
  if (is.null(h)) {
    h <- c(k = k, q = qm, lS = lm)
  } else {
    if (qm > h[["q"]]) h[c("k", "q")] <- c(k, qm)
    h[["lS"]] <- max(h[["lS"]], lm)
  }
  .gwaring_acc$hit <- h
}
.gwaring_watch <- function(expr, call, lower_p = FALSE) {
  .gwaring_acc$hit <- NULL
  on.exit(.gwaring_acc$hit <- NULL)
  val <- expr
  hit <- .gwaring_acc$hit
  if (!is.null(hit) && !(lower_p && isTRUE(hit[["lS"]] <= .GW_ACC_LS)))
    warning(structure(class = c("egpd_gwaring_accuracy", "warning", "condition"), list(
      message = sprintf(paste0(
        "generalized Waring tail with k = %.3g > 1e15 evaluated at a count of %.0f > 1000: ",
        "beyond the fitting bound on k the tail integral can be off by up to about 1e-3 ",
        "(relative) in P(Y > q) (below 1e-9 for counts up to 2e4), and P(Y <= q) by up to ",
        "1e-3 times P(Y > q), so these results are approximate"), hit[["k"]], hit[["q"]]),
      call = call)))
  val
}

## ---------------------------------------------------------------- Poisson-lognormal

#' Poisson-lognormal distribution
#'
#' \code{Y | Z ~ Poisson(exp(Z))} with \code{Z ~ N(muz, sigma^2)} and
#' \code{muz = log(mu) - sigma^2 / 2}, so \code{E[Y] = mu}. The pmf has no closed form and
#' is evaluated by adaptive Gauss-Kronrod quadrature centred on the mode of the integrand.
#' The tail is lognormal-like: heavier than the negative binomial's, but still in the
#' Gumbel domain (\code{xi = 0}), so this family models a heavy body rather than a heavy
#' tail.
#'
#' @param x,q vector of counts
#' @param p vector of probabilities
#' @param n number of draws; if \code{length(n) != 1}, the length is taken to be
#'   the number required (\code{numeric(0)} gives none), as in \code{\link[stats]{rpois}}
#' @param mu mean, \code{mu > 0}
#' @param sigma standard deviation of the log-intensity, \code{sigma > 0}
#' @param log,log.p logical; return the log density / probability
#' @param lower.tail logical; if TRUE probabilities are \code{P[X <= x]}
#' @param max.value largest count searched when inverting the cdf; nonattaining
#'   caps return \code{NA} with a warning. Invalid probabilities return \code{NaN},
#'   with endpoints 0 and \code{Inf}.
#' @return a numeric vector. As with \code{\link[stats]{rpois}}, \code{rpoislnorm}
#'   returns integers unless a draw exceeds \code{.Machine$integer.max}, when it
#'   returns doubles.
#' @name poislnorm
NULL

#' @rdname poislnorm
#' @export
dpoislnorm <- function(x, mu = 1, sigma = 1, log = FALSE) {
  .count_d(.poislnorm_fam, x, list(mu = mu, sigma = sigma), log = log)
}
#' @rdname poislnorm
#' @export
ppoislnorm <- function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .count_p(.poislnorm_fam, q, list(mu = mu, sigma = sigma), lower.tail, log.p)
}
#' @rdname poislnorm
#' @export
qpoislnorm <- function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE,
                       max.value = 1e6) {
  .count_q(.poislnorm_fam, p, list(mu = mu, sigma = sigma), lower.tail, log.p, max.value)
}
#' @rdname poislnorm
#' @export
rpoislnorm <- function(n, mu = 1, sigma = 1) {
  n <- .count_rn(n)
  if (n == 0L) return(integer(0))
  if (length(mu) == 0L || length(sigma) == 0L) {
    warning("NaNs produced", call. = FALSE)
    return(rep(NaN, n))
  }
  mu <- rep_len(mu, n); sigma <- rep_len(sigma, n)
  st <- .count_state(list(mu = mu, sigma = sigma), .poislnorm_fam$valid)
  out <- rep(NaN, n)
  out[st$na] <- NA_real_
  if (any(st$ok)) {
    m <- sum(st$ok)
    z <- stats::rnorm(m, log(mu[st$ok]) - sigma[st$ok]^2 / 2, sigma[st$ok])
    out[st$ok] <- stats::rpois(m, exp(z))
  }
  ## integer, as stats::rpois, unless a draw exceeds the integer range
  if (all(st$ok) && !any(out > .Machine$integer.max, na.rm = TRUE)) as.integer(out) else out
}

## Lognormal tails are subexponential (summing them to negligible remainder can take
## ~1e11 terms), so the remainder beyond the ladder is the exact mixture integral of the
## Poisson survivor (src/countfams.cpp, pln_lcdf).
.poislnorm_fam <- list(
  pars = c("mu", "sigma"),
  valid = function(p) is.finite(p$mu) & p$mu > 0 & is.finite(p$sigma) & p$sigma > 0,
  lpmf = function(n, p) plnorm_logpmf_cpp(0:n, p$mu, p$sigma),
  lpmf_at = function(x, p) plnorm_logpmf_cpp(as.numeric(x), p$mu, p$sigma),
  lsurv = function(k, p) plnorm_lcdf_cpp(as.numeric(k), p$mu, p$sigma, FALSE),
  init = function(p) max(32, ceiling(4 * p$mu)))

## Shared vectorisation, RNG length, probability conversion, and cdf inversion helpers.
## As in base R, a vector n of length other than 1 gives length(n) draws (numeric(0): none).
.count_rn <- function(n) {
  if (length(n) != 1L) return(length(n))
  if (is.na(n) || !is.finite(n) || n < 0)
    stop("'n' must be a non-negative integer")
  as.integer(floor(n))
}

.count_prob <- function(p, lower.tail, log.p) {
  if (log.p) {
    if (lower.tail) exp(p) else -expm1(p)
  } else {
    if (lower.tail) p else 1 - p
  }
}

.count_groups <- function(pars_list) {
  k <- length(pars_list[[1L]])
  if (k == 0L) return(list())
  if (all(vapply(pars_list, function(v) length(unique(v)) <= 1L, TRUE)))
    return(list(seq_len(k)))
  key <- do.call(paste, c(lapply(pars_list, sprintf, fmt = "%a"), sep = "|"))
  unname(split(seq_len(k), factor(key, levels = unique(key))))
}

.count_logadd <- function(x, y) {
  m <- pmax(x, y)
  out <- m + log1p(exp(pmin(x, y) - m))
  out[!is.na(m) & m == -Inf] <- -Inf
  out
}

.count_cap <- function(max.value) {
  if (length(max.value) != 1L || !is.finite(max.value) ||
      max.value < 0 || max.value != floor(max.value) ||
      max.value > .Machine$integer.max - 1)
    stop("'max.value' must be a non-negative integer within integer range")
  max.value
}

## ---------------------------------------------------------------------------------
## Common d/p/q/r engine for the count families. A family is a list
##   pars   parameter names, in order
##   valid  function(par list) -> logical: which rows are valid
##   lpmf   function(n, p1) -> log p_0, ..., log p_n at one parameter row p1
##   lpmf_at (optional) function(x, p1) -> log pmf at the integers x
##   lsurv  (optional) function(k, p1) -> exact log P(X > k); used for the remainder beyond
##          the ladder (heavy tails), or, with surv_all = TRUE, at the requested k alone
##          (closed-form survivor: no ladder beyond the far lower tail)
##   rlim   (optional) function(p1) -> limit of p_{k+1}/p_k, which bounds the remainder of
##          a light tail beyond the ladder
##   maxN   (optional) function(p1, need) -> longest ladder for summing the tail
##   slow_ladder (optional) TRUE when the ladder costs O(n^2): the survivor beyond the
##          forward sum is then formed only when needed to resolve the requested tail
##   init   function(p1) -> first extent searched by the quantile function
##
## Argument contracts, uniform across families and following stats::dpois:
##   * a parameter that is NA (not NaN) gives NA, silently;
##   * an invalid parameter (NaN, +-Inf or out of range) gives NaN with a warning;
##   * x: NA stays NA and NaN stays NaN, silently, even with an invalid parameter
##     (dpois(NA, -1) is NA); x < 0 and x = +-Inf have mass 0; a non-integer x
##     has mass 0 with a "non-integer x" warning;
##   * q: NA/NaN propagate, q < 0 gives F = 0, q = Inf gives F = 1, otherwise floor(q);
##   * p outside [0, 1] (log.p: outside [-Inf, 0]) gives NaN with a warning; the
##     quantile endpoints are 0 and Inf.
##
## Tails. F(k) is the forward sum of the pmf. Where the survivor S(k) = P(X > k) is below
## 1e-3 it is summed from the far tail (or taken from an exact tail formula), and the upper
## tail and log F = log(1 - S) are reported from it; both tails are then accurate to
## rounding on both the probability and the log scale, where 1 - cumsum(pmf) lost all
## relative accuracy below ~1e-12. If the tail sum cannot be certified within the family's
## summation limit, the upper tail is NaN with a warning while the lower tail keeps the
## forward sum (bounded by one).
## The quantile function inverts exactly these values -- on the survivor scale for upper
## tails -- with a relative fuzz of 64 ulps, as base R's discrete quantile functions long
## used; on the log scale it is relative to 1 - p where log p is near 0.

## `x` (optional) is the first argument of a d/p/q function: an NA or NaN there takes
## precedence over an invalid parameter, as in stats::dpois, and is not flagged as bad.
.count_state <- function(pl, valid, zi = FALSE, x = NULL) {
  na <- Reduce(`|`, lapply(pl, function(v) is.na(v) & !is.nan(v)), FALSE)
  ok <- valid(pl)
  if (zi) ok <- ok & is.finite(pl$pi) & pl$pi >= 0 & pl$pi < 1
  ok[is.na(ok)] <- FALSE
  ok <- ok & !na
  bad <- !ok & !na
  if (!is.null(x)) bad <- bad & !is.na(x)
  if (any(bad)) warning("NaNs produced", call. = FALSE)
  list(ok = ok, na = na, bad = bad)
}

.count_lpmf_at <- function(fam, x, p1) {
  if (!is.null(fam$lpmf_at)) return(fam$lpmf_at(x, p1))
  fam$lpmf(max(x), p1)[x + 1]
}

## log pmf (or pmf) for a family, optionally zero-inflated with probability pi.
.count_d <- function(fam, x, pars, log = FALSE, pi = NULL) {
  zi <- !is.null(pi)
  z <- do.call(.gpig_recycle, c(list(x = x), pars, if (zi) list(pi = pi)))
  xv <- z$x
  n <- length(xv)
  if (n == 0L) return(numeric(0))
  pz <- z[fam$pars]
  st <- .count_state(c(pz, if (zi) list(pi = z$pi)), fam$valid, zi, xv)
  out <- rep(NA_real_, n)
  out[st$bad] <- NaN
  isna <- !st$na & is.na(xv)
  out[isna] <- xv[isna]
  fin <- st$ok & !is.na(xv)
  nonint <- fin & is.finite(xv) & xv != floor(xv)
  for (v in unique(xv[nonint]))
    warning(sprintf("non-integer x = %f", v), call. = FALSE)
  out[fin & (!is.finite(xv) | xv < 0 | nonint)] <- -Inf
  good <- fin & is.finite(xv) & xv >= 0 & !nonint
  if (any(good)) {
    gi <- which(good)
    for (g in .count_groups(lapply(pz, `[`, gi))) {
      ii <- gi[g]
      out[ii] <- .count_lpmf_at(fam, xv[ii], lapply(pz, `[`, ii[1L]))
    }
    if (zi) {
      piv <- z$pi[gi]
      lp <- log1p(-piv) + out[gi]
      zero <- xv[gi] == 0
      if (any(zero)) {
        ## log P(0) = log(pi + (1 - pi) p0) is a log-probability, but the log-sum keeps
        ## only absolute accuracy and can round above 0: where P(Y > 0) =
        ## (1 - pi)(1 - p0) is below 1e-3 it is log1p(-P(Y > 0)) instead
        l0 <- .count_logadd(log(piv[zero]), lp[zero])
        m <- (1 - piv[zero]) * -expm1(out[gi][zero])
        near <- !is.na(m) & m < 1e-3 & piv[zero] > 0
        l0[near] <- log1p(-m[near])
        lp[zero] <- pmin(l0, 0)
      }
      out[gi] <- lp
    }
  }
  if (log) out else exp(out)
}

## log(1 - e^l), l <= 0, accurate on both sides of log(1/2).
.count_l1mexp <- function(l) ifelse(l > -log(2), log(-expm1(l)), log1p(-exp(l)))

## log F(k) and log S(k) = log P(X > k), k = 0..need, for one parameter row. The survivor
## below tail_cut is summed directly, and NaN where it cannot be resolved
## (the bounded forward lF is then kept). For other rows use log(1 - F).
.count_tails <- function(fam, p1, need, tail_cut = 1e-3) {
  lp <- fam$lpmf(need, p1)
  lF <- pmin(count_lcse_cpp(lp, FALSE, -Inf), 0)
  Sc <- -expm1(lF)
  far <- Sc < tail_cut
  lS <- rep(NA_real_, length(Sc))
  lS[!far] <- .count_l1mexp(lF[!far])
  if (any(far)) {
    unresolved <- far
    if (!is.null(fam$positive_tail)) {
      pos <- fam$positive_tail(need, p1)
      if (length(pos)) {
        good <- far & !is.na(pos)
        lS[good] <- pos[good]
        unresolved[good] <- FALSE
      }
    }
    if (any(unresolved)) {
      k0 <- which(unresolved)[1L] - 1L
      lS[unresolved] <- .count_tail_sums(fam, p1, k0, need, lp, Sc[unresolved])
    }
    ## Recover the analytically known zero atom for an unresolved heavy
    ## tail. Do not overwrite a resolved survivor with a rounded log P(0).
    if (is.nan(lS[1L])) lS[1L] <- .count_l1mexp(lp[1L])
    res <- far & !is.nan(lS)
    lF[res] <- log1p(-exp(lS[res]))
  }
  list(lF = lF, lS = lS)
}

## Families with an exact survivor at every k (surv_all): S from lsurv at the requested k
## alone, F = 1 - S, except in the far lower tail (F < 1e-3), where the forward sum of the
## pmf keeps F's relative accuracy -- that ladder stops at the 1e-3 quantile.
.count_tails_surv <- function(fam, p1, k) {
  lS <- fam$lsurv(k, p1)
  lF <- .count_l1mexp(lS)
  low <- !is.na(lS) & lS > log1p(-1e-3)
  if (any(low)) {
    kl <- k[low]
    lF[low] <- count_lcse_cpp(fam$lpmf(max(kl), p1), FALSE, -Inf)[kl + 1]
  }
  list(lF = lF, lS = lS)
}

## log S(k) for k = k0..need, summed from the tail of the ladder; NaN where the sum cannot
## be certified within the summation limit.
.count_tail_sums <- function(fam, p1, k0, need, lp, Sc) {
  if (!is.null(fam$lsurv)) {
    ltail <- fam$lsurv(need, p1)
    rs <- if (need > k0) count_lcse_cpp(lp[(k0 + 2L):(need + 1L)], TRUE, ltail) else
      numeric(0)
    return(c(rs, ltail))
  }
  ## The ladder must reach past need, whatever the limit; a family's maxN (GPIG's ladder
  ## costs O(N^2)) bounds the extra terms, otherwise up to 1e5 more are allowed.
  cap <- if (is.null(fam$maxN)) max(1e7, need + 1e5) else max(fam$maxN(p1, need), need + 1)
  rlim <- if (is.null(fam$rlim)) 0 else fam$rlim(p1)
  ## For an untempered power tail no extension can certify a geometric
  ## remainder. Apply the same unresolved-tail rule without building a ladder.
  if (rlim >= 1) return(ifelse(Sc >= 1e-8, log(Sc), NaN))
  N <- min(cap, max(2 * need, need + 64))
  repeat {
    lpN <- fam$lpmf(N, p1)
    last <- lpN[N + 1L]
    if (last == -Inf) {
      lrem <- -Inf
    } else {
      r <- max(exp(last - lpN[N]), rlim)
      lrem <- if (r < 1) last + log(r) - log1p(-r) else Inf
    }
    out <- count_lcse_cpp(lpN[(k0 + 2L):(N + 1L)], TRUE, -Inf)[seq_len(need - k0 + 1L)]
    if (!(lrem >= out[length(out)] + log(1e-17)) || N >= cap) break
    N <- min(cap, 2 * N)
  }
  ## With the cap reached the partial sums are lower bounds: keep those certified to 1e-12
  ## (a remainder that is exactly zero certifies them all), fall back to 1 - F where that is
  ## still accurate, and leave the rest unresolved.
  bad <- lrem > -Inf & lrem >= out + log(1e-12)
  use_c <- bad & Sc >= 1e-8
  out[use_c] <- log(Sc[use_c])
  out[bad & !use_c] <- NaN
  out
}

## Zero inflation of a tail table: F_zi = pi + (1 - pi) F, S_zi = (1 - pi) S. log F_zi is a
## log-probability: the log-sum can round it above 0 when F is 1 to rounding.
.count_tails_zi <- function(tab, pi) {
  lF <- ifelse(tab$lF > -log(2), log1p((1 - pi) * expm1(tab$lF)),
               .count_logadd(log(pi), log1p(-pi) + tab$lF))
  list(lF = pmin(lF, 0), lS = log1p(-pi) + tab$lS)
}

## Tail probabilities reported by the p-functions, from a tail table.
.count_tail_out <- function(tab, lower.tail, log.p) {
  far <- !is.na(tab$lS) & tab$lS < log(1e-3)
  if (lower.tail) {
    if (log.p) ifelse(far, log1p(-exp(tab$lS)), tab$lF)
    else ifelse(far, -expm1(tab$lS), exp(tab$lF))
  } else {
    if (log.p) tab$lS else exp(tab$lS)
  }
}

.count_p <- function(fam, q, pars, lower.tail = TRUE, log.p = FALSE, pi = NULL) {
  zi <- !is.null(pi)
  z <- do.call(.gpig_recycle, c(list(q = q), pars, if (zi) list(pi = pi)))
  qv <- z$q
  n <- length(qv)
  if (n == 0L) return(numeric(0))
  pz <- z[fam$pars]
  st <- .count_state(c(pz, if (zi) list(pi = z$pi)), fam$valid, zi, qv)
  out <- rep(NA_real_, n)
  out[st$bad] <- NaN
  isna <- !st$na & is.na(qv)
  out[isna] <- qv[isna]
  lo <- if (lower.tail) (if (log.p) -Inf else 0) else (if (log.p) 0 else 1)
  hi <- if (lower.tail) (if (log.p) 0 else 1) else (if (log.p) -Inf else 0)
  out[st$ok & !is.na(qv) & qv < 0] <- lo
  out[st$ok & !is.na(qv) & qv == Inf] <- hi
  good <- st$ok & is.finite(qv) & qv >= 0
  if (any(good)) {
    gi <- which(good)
    k <- floor(qv[gi])
    ## Even a linear lower CDF needs the survivor where rounding in a
    ## slow forward sum can swamp its remaining probability.
    tail_cut <- if (lower.tail && !log.p && isTRUE(fam$slow_ladder))
      sqrt(.Machine$double.eps) else 1e-3
    for (g in .count_groups(lapply(c(pz, if (zi) list(pi = z$pi)), `[`, gi))) {
      ii <- gi[g]
      p1 <- lapply(pz, `[`, ii[1L])
      if (isTRUE(fam$surv_all)) {
        tab <- .count_tails_surv(fam, p1, k[g])
      } else {
        tab <- .count_tails(fam, p1, max(k[g]), tail_cut)
        tab <- list(lF = tab$lF[k[g] + 1], lS = tab$lS[k[g] + 1])
      }
      if (zi) tab <- .count_tails_zi(tab, z$pi[ii[1L]])
      out[ii] <- .count_tail_out(tab, lower.tail, log.p)
    }
    if (!lower.tail && any(is.nan(out[gi])))
      warning("upper tail not resolved within the summation limit; returning NaN",
              call. = FALSE)
  }
  out
}

## Generalised inverse at interior probabilities pv (one parameter row): the smallest k
## with F(k) >= p (lower tail) or S(k) <= p (upper tail), on the values the p-function
## reports, with a relative fuzz of 64 ulps. On the log scale it is the shift
## log(1 -/+ fuzz) of log p for log p <= -1, and a relative change of log p itself above
## -1: near 0, log p ~ -(1 - p) carries 1 - p to full relative accuracy, which a fixed
## shift would swamp (base R's discrete quantiles also scale log p).
.count_qsearch <- function(fam, pv, p1, lower.tail, log.p, max.value, pi = NULL) {
  if (!lower.tail && isTRUE(fam$slow_ladder)) {
    ls0 <- .count_l1mexp(fam$lpmf(0L, p1)[1L])
    if (!is.null(pi)) ls0 <- log1p(-pi) + ls0
    at0 <- (if (log.p) pv else log(pv)) >= ls0
    if (any(at0)) {
      out <- numeric(length(pv))
      if (any(!at0))
        out[!at0] <- .count_qsearch(fam, pv[!at0], p1, FALSE, log.p, max.value, pi)
      return(out)
    }
  }
  if (lower.tail && !log.p) {
    far <- 1 - pv < sqrt(.Machine$double.eps)
    if (any(far)) {
      ## A tolerance relative to p near 1 can exceed the requested
      ## survivor. The log-CDF comparison scales it to the small tail and
      ## still works from the forward sum if a survivor is unavailable.
      out <- numeric(length(pv))
      q <- .count_qsearch(fam, log(pv[far]), p1, TRUE, TRUE, max.value, pi)
      if (!isTRUE(fam$slow_ladder)) {
        ## A rounded-up p(k) can make the log-tail inverse k+1. Honor
        ## that published knot only when both adjacent jumps are resolved;
        ## never use a plateau of rounded CDF values to move the quantile.
        k <- ifelse(is.na(q), as.double(max.value), q - 1)
        i <- which(k >= 0)
        if (length(i)) {
          at <- outer(k[i], -1:1, "+")
          F <- matrix(.count_p(fam, as.vector(at), p1, pi = pi), ncol = 3L)
          gap <- .Machine$double.eps
          tie <- which(F[,2] == pv[far][i] & F[,2] - F[,1] > gap & F[,3] - F[,2] > gap)
          q[i[tie]] <- k[i[tie]]
        }
      }
      out[far] <- q
      if (any(!far))
        out[!far] <- .count_qsearch(fam, pv[!far], p1, TRUE, FALSE, max.value, pi)
      return(out)
    }
  }
  fz <- 64 * .Machine$double.eps
  tgt <- if (log.p) {
    if (lower.tail) pmax(pv + log1p(-fz), pv * (1 + fz)) else
      pmin(pv + log1p(fz), pv * (1 - fz))
  } else if (lower.tail) pv - fz * pmin(pv, 1 - pv) else pv * (1 + fz)
  if (isTRUE(fam$surv_all))
    return(.count_qsearch_surv(fam, tgt, p1, lower.tail, log.p, max.value, pi))
  tail_cut <- if (lower.tail && !log.p && isTRUE(fam$slow_ladder)) 0 else 1e-3
  need <- min(max.value, fam$init(p1))
  repeat {
    tab <- .count_tails(fam, p1, need, tail_cut)
    if (!is.null(pi)) tab <- .count_tails_zi(tab, pi)
    v <- .count_tail_out(tab, lower.tail, log.p)
    end <- v[length(v)]
    done <- if (lower.tail) isTRUE(end >= max(tgt)) else isTRUE(end <= min(tgt))
    if (done || need >= max.value) break
    need <- min(max.value, max(4 * need, need + 16))
  }
  unresolved <- is.na(end)
  v[is.na(v)] <- if (lower.tail) -Inf else Inf
  i <- if (lower.tail) findInterval(tgt, cummax(v), left.open = TRUE) else
    findInterval(-tgt, cummax(-v), left.open = TRUE)
  out <- as.numeric(i)
  out[i > need] <- NA_real_
  if (unresolved && anyNA(out)) attr(out,"unresolved") <- TRUE
  out
}

## The same generalised inverse for surv_all families, by bisection on the pointwise
## values (no table of the support): an exponential search brackets the most extreme
## target, then every target is bisected at once.
.count_qsearch_surv <- function(fam, tgt, p1, lower.tail, log.p, max.value, pi = NULL) {
  val <- function(k) {
    tab <- .count_tails_surv(fam, p1, k)
    if (!is.null(pi)) tab <- .count_tails_zi(tab, pi)
    v <- .count_tail_out(tab, lower.tail, log.p)
    v[is.na(v)] <- if (lower.tail) -Inf else Inf
    v
  }
  hit <- function(v, t) if (lower.tail) v >= t else v <= t
  ext <- if (lower.tail) max(tgt) else min(tgt)
  hi <- min(max.value, fam$init(p1))
  while (!hit(val(hi), ext) && hi < max.value) hi <- min(max.value, 4 * hi + 16)
  n <- length(tgt)
  att <- hit(val(rep(hi, n)), tgt)
  lo <- rep(-1, n)                                  # F(-1) = 0, S(-1) = 1: never a hit
  up <- rep(hi, n)
  repeat {
    act <- which(up - lo > 1)
    if (!length(act)) break
    mid <- floor((lo[act] + up[act]) / 2)
    h <- hit(val(mid), tgt[act])
    up[act[h]] <- mid[h]
    lo[act[!h]] <- mid[!h]
  }
  up[!att] <- NA_real_
  up
}

.count_q <- function(fam, p, pars, lower.tail = TRUE, log.p = FALSE, max.value,
                     pi = NULL) {
  max.value <- .count_cap(max.value)
  zi <- !is.null(pi)
  z <- do.call(.gpig_recycle, c(list(p = p), pars, if (zi) list(pi = pi)))
  pv <- z$p
  n <- length(pv)
  if (n == 0L) return(numeric(0))
  pz <- z[fam$pars]
  st <- .count_state(c(pz, if (zi) list(pi = z$pi)), fam$valid, zi, pv)
  out <- rep(NA_real_, n)
  out[st$bad] <- NaN
  isna <- !st$na & is.na(pv)
  out[isna] <- pv[isna]
  inr <- if (log.p) pv <= 0 else pv >= 0 & pv <= 1
  oor <- st$ok & !is.na(pv) & !inr
  if (any(oor)) {
    warning("NaNs produced", call. = FALSE)
    out[oor] <- NaN
  }
  ok <- st$ok & !is.na(pv) & inr
  if (!is.null(fam$point0) && any(ok)) {             # a point mass at zero: q = 0 for all p
    oi <- which(ok)
    d0 <- fam$point0(lapply(pz, `[`, oi))
    out[oi[d0]] <- 0
    ok[oi[d0]] <- FALSE
  }
  p0 <- if (log.p) pv == -Inf else pv == 0
  p1e <- if (log.p) pv == 0 else pv == 1
  out[ok & (if (lower.tail) p0 else p1e)] <- 0
  out[ok & (if (lower.tail) p1e else p0)] <- Inf
  mid <- ok & !p0 & !p1e
  if (any(mid)) {
    gi <- which(mid)
    unresolved <- short <- FALSE
    groups <- lapply(c(pz, if (zi) list(pi = z$pi)), `[`, gi)
    ## For p > 1/2 the complement is exact in binary arithmetic. Comparing
    ## near-one survivors (even with no fuzz) can merge distinct lower probabilities.
    flip <- !lower.tail & !log.p & pv[gi] > .5
    if (any(flip)) {
      pv[gi[flip]] <- 1 - pv[gi[flip]]
      groups$lower <- as.numeric(flip)
    }
    for (g in .count_groups(groups)) {
      ii <- gi[g]
      q <- .count_qsearch(fam, pv[ii], lapply(pz, `[`, ii[1L]), lower.tail || flip[g[1L]], log.p,
                          max.value, if (zi) z$pi[ii[1L]])
      if (isTRUE(attr(q,"unresolved"))) unresolved <- TRUE else if (anyNA(q)) short <- TRUE
      out[ii] <- q
    }
    if (unresolved) warning("Quantile could not be resolved within the tail summation limit.",
                            call. = FALSE)
    if (short) warning("quantile not attained within 'max.value'", call. = FALSE)
  }
  out
}

## Random generation by inversion; zero-inflated families first draw the structural zero.
.count_r <- function(fam, n, pars, max.value, pi = NULL) {
  n <- .count_rn(n)
  if (n == 0L) return(numeric(0))
  zi <- !is.null(pi)
  if (any(lengths(pars) == 0L) || (zi && length(pi) == 0L)) {
    warning("NaNs produced", call. = FALSE)
    return(rep(NaN, n))
  }
  pars <- lapply(pars, rep_len, n)
  if (!zi) return(.count_q(fam, stats::runif(n), pars, max.value = max.value))
  pi <- rep_len(pi, n)
  st <- .count_state(c(pars, list(pi = pi)), fam$valid, TRUE)
  out <- rep(NaN, n)
  out[st$na] <- NA_real_
  u <- stats::runif(n)
  z0 <- st$ok & u < pi
  out[z0] <- 0
  draw <- st$ok & !z0
  if (any(draw))
    out[draw] <- .count_q(fam, stats::runif(sum(draw)), lapply(pars, `[`, draw),
                          max.value = max.value)
  out
}

## Random generation by a family's exact sampler: `draw(p)` returns one draw per row of the
## valid parameter list p. NA, invalid and zero-length parameters, and the structural zero
## of zero-inflated families, are handled as in .count_r.
.count_r_exact <- function(fam, n, pars, draw, pi = NULL) {
  n <- .count_rn(n)
  if (n == 0L) return(numeric(0))
  zi <- !is.null(pi)
  if (any(lengths(pars) == 0L) || (zi && length(pi) == 0L)) {
    warning("NaNs produced", call. = FALSE)
    return(rep(NaN, n))
  }
  pars <- lapply(pars, rep_len, n)
  pl <- if (zi) c(pars, list(pi = rep_len(pi, n))) else pars
  st <- .count_state(pl, fam$valid, zi)
  out <- rep(NaN, n)
  out[st$na] <- NA_real_
  ok <- st$ok
  if (zi) {
    z0 <- ok & stats::runif(n) < pl$pi
    out[z0] <- 0
    ok <- ok & !z0
  }
  if (any(ok)) out[ok] <- draw(lapply(pars, `[`, ok))
  out
}

## ---------------------------------------------------------------- likelihood bundles
## Each wrapper splits the flat coefficient vector by parameter and calls the Rcpp routine,
## mirroring R/pig_lik.R and R/gpig_lik.R.

.iG_gpois   <- function(v, mu, lambda) qgpois(v, mu = mu, lambda = lambda)
.iG_gwaring <- function(v, mu, k, rho)  qgwaring(v, mu = mu, k = k, rho = rho)
.iG_plnorm  <- function(v, mu, sigma)   qpoislnorm(v, mu = mu, sigma = sigma)

.gpois1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for gpois.")
  gpois1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
           likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
           likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gpois1.d12 <- function(pars, likdata) {
  gpois1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
            likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
            likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gpois1fns <- list(d0 = .gpois1.d0, d120 = .gpois1.d12, d340 = NULL, m = 1, iG = .iG_gpois)

.gwaring1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for gwaring.")
  gwaring1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
             likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
             likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gwaring1.d12 <- function(pars, likdata) {
  gwaring1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
              likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[, 1],
              likdata$dupid, likdata$duplicate, likdata$offsets)
}
.gwaring1fns <- list(d0 = .gwaring1.d0, d120 = .gwaring1.d12, d340 = NULL, m = 1,
                     iG = .iG_gwaring)

.plnorm1.d0 <- function(pars, likdata) {
  if (likdata$censored) stop("Censored likelihoods not available for plnorm.")
  plnorm1d0(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
            likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
            likdata$dupid, likdata$duplicate, likdata$offsets)
}
.plnorm1.d12 <- function(pars, likdata) {
  plnorm1d12(split(pars, factor(likdata$idpars, levels = seq_along(likdata$X))),
             likdata$X[[1]], likdata$X[[2]], likdata$y[, 1],
             likdata$dupid, likdata$duplicate, likdata$offsets)
}
.plnorm1fns <- list(d0 = .plnorm1.d0, d120 = .plnorm1.d12, d340 = NULL, m = 1,
                    iG = .iG_plnorm)
