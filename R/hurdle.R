## Hurdle (zero-adjusted) count models ----------------------------------------
##
## A hurdle model gives the zero probability its own parameter instead of taking
## whatever the count distribution happens to imply:
##
##   P(Y = 0) = pi,    P(Y = y) = (1 - pi) f(y) / (1 - f(0)),  y >= 1
##
## Unlike zero *inflation*, which can only add mass at zero, pi is free in both
## directions. That matters here: generalized Waring systematically over-predicts
## zeros on disease counts, so the correction needed is downwards.
##
## The log-likelihood factorises, and the two blocks share no parameters:
##
##   l = sum_i [1{y=0} log pi + 1{y>0} log(1 - pi)]        <- Bernoulli
##     + sum_{y>0} [log f(y) - log(1 - f(0))]              <- zero-truncated base
##
## so the hurdle half is an ordinary binomial GAM and only the second block needs
## anything new. For that block the negative log-likelihood and its derivatives follow
## from the base family's own d0/d12 evaluated twice -- at y, and at 0 -- with no
## new derivative algebra. With w_j = f0_j / (1 - f0_j) and g0, H0 the base
## gradient and Hessian of -log f0:
##
##   nll_j        = base_nll_j + log(1 - f0_j)
##   dnll_j/deta  = base_g_j(y) + w_j * g0_j
##   d2nll_j      = base_H_j(y) + w_j * H0_j - f0_j / (1 - f0_j)^2 * g0_j g0_j'
##
## For DEGPD model 1, whose base Hessians are exact, the second-order block is the
## exact observed information. The count families of src/countfams.cpp use the
## outer product of the per-observation gradient (BHHH; cf_pack()) throughout, and
## the zero-truncated generalized Waring follows them: with smooths on its weakly
## identified k and rho the exact observed information is indefinite, and the REML
## and edf computations then break down (edf -2937 on one test series), while the
## outer product stays positive semi-definite. Its standard errors are therefore
## BHHH-based, as for egpd(family = "gwaring").

## linear predictors, mirroring cf_eta() in src/countfams.cpp
.zt_eta <- function(pars, likdata) {
  b <- split(pars, factor(likdata$idpars, levels = seq_along(likdata$X)))
  lapply(seq_along(likdata$X), function(i) {
    e <- as.numeric(likdata$X[[i]] %*% b[[i]])
    if (isTRUE(likdata$duplicate == 1)) e <- e[likdata$dupid + 1]
    o <- likdata$offsets[[i]]
    if (length(o) > 0) e <- e + o
    e
  })
}

## The optimiser visits parameter values the density functions reject outright, so
## failures return NA and are turned into an infinite negative log-likelihood by the
## caller rather than an error.
.zt_safe <- function(expr) {
  v <- tryCatch(suppressWarnings(expr), error = function(e) NA_real_)
  if (!is.numeric(v)) NA_real_ else v
}

.zt_degpd1_xi <- function(eta2, likdata = NULL) {
  use_bounded <- !is.null(likdata) &&
    !identical(likdata$bounded.xi, FALSE) &&
    !is.null(likdata$xi.max) && is.finite(likdata$xi.max)
  if (use_bounded) {
    a <- likdata$xi.max
    tt <- exp(eta2 - log(a))
    ifelse(tt > 700, a, a * (-expm1(-tt)))
  } else {
    exp(eta2)
  }
}

## The truncation term needs log(1 - f0) = log(-expm1(log f0)) to full relative
## accuracy: f0 -> 1 as k -> Inf for the Waring and as sigma -> 0 for DEGPD, where
## log1p(-f0) of a linear-scale f0 (or a guard such as f0 < 1 - 1e-12) breaks the
## objective. So log f0 itself is computed without cancellation, vectorised over rows.

## lgamma(z + d) - lgamma(z) and digamma(z + d) - digamma(z) without cancellation:
## a Taylor series in d when d << z (as cf_lgamma_diff() in src/countfams.cpp).
.cf_lgamma_diff <- function(z, d) {
  out <- lgamma(z + d) - lgamma(z)
  s <- which(d < 1e-4 * z)
  if (length(s)) {
    zs <- z[s]; ds <- d[s]; d2 <- ds * ds
    out[s] <- ds * digamma(zs) + 0.5 * d2 * trigamma(zs) +
      (d2 * ds / 6) * psigamma(zs, 2L) + (d2 * d2 / 24) * psigamma(zs, 3L)
  }
  out[which(d == 0)] <- 0
  out
}
.cf_digamma_diff <- function(z, d) {
  out <- digamma(z + d) - digamma(z)
  s <- which(d < 1e-4 * z)
  if (length(s)) {
    zs <- z[s]; ds <- d[s]; d2 <- ds * ds
    out[s] <- ds * trigamma(zs) + 0.5 * d2 * psigamma(zs, 2L) +
      (d2 * ds / 6) * psigamma(zs, 3L) + (d2 * d2 / 24) * psigamma(zs, 4L)
  }
  out[which(d == 0)] <- 0
  out
}

## log f(0) of the generalized Waring and, if `deriv`, its gradient in the linear
## predictors (log mu, log k, log rho). The parameters are clamped exactly as the C++
## likelihood clamps them (gwaring_lpg() in src/countfams.cpp, bounds from
## cf_bounds_cpp(): mu in [1e-10, 1e12], k in [1e-8, 1e15], rho in [1 + 1e-6, 1e15];
## a = mu (rho - 1) / k is not clamped, as in C++), with zero derivative in a clamped
## direction, and
##   log p(0) = lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) - lgamma(a + k + rho)
## is evaluated as in dgwaring(0, log = TRUE): the pair of lgamma differences in the
## smaller of a and k, each without cancellation. The C++ score at y = 0 suffers that
## cancellation (e.g. psi(k + rho) - psi(a + k + rho) as k -> Inf), and the truncation
## term multiplies it by f0 / (1 - f0), so the gradient is evaluated here too.
.zt_gwaring_logf0 <- function(eta, deriv = FALSE) {
  bd <- cf_bounds_cpp()[["gwaring"]]
  raw <- lapply(eta[1:3], exp)
  clamp <- function(v, lo, hi) pmin(pmax(v, lo), hi)
  mu  <- clamp(raw[[1]], bd[["mu_lo"]],  bd[["mu_hi"]])
  k   <- clamp(raw[[2]], bd[["k_lo"]],   bd[["k_hi"]])
  rho <- clamp(raw[[3]], bd[["rho_lo"]], bd[["rho_hi"]])
  a <- mu * ((rho - 1) / k)
  by_a <- !is.na(a) & !is.na(k) & a <= k
  val <- ifelse(by_a,
                .cf_lgamma_diff(rho, a) - .cf_lgamma_diff(rho + k, a),
                .cf_lgamma_diff(rho, k) - .cf_lgamma_diff(rho + a, k))
  if (!deriv) return(val)
  d_a <- -.cf_digamma_diff(a + rho, k)                  # d log f0 / d a
  d_k <- -.cf_digamma_diff(k + rho, a)                  # d log f0 / d k
  d_r <- ifelse(by_a,                                   # d log f0 / d rho, a fixed
                .cf_digamma_diff(rho, a) - .cf_digamma_diff(rho + k, a),
                .cf_digamma_diff(rho, k) - .cf_digamma_diff(rho + a, k))
  ## a = mu (rho - 1) / k: da/dlog mu = a, da/dlog k = -a, da/dlog rho = a rho / (rho - 1)
  grad <- cbind(a * d_a, k * d_k - a * d_a, rho * d_r + a * rho / (rho - 1) * d_a)
  outside <- function(v, lo, hi) !is.na(v) & (v < lo | v > hi)
  grad[outside(raw[[1]], bd[["mu_lo"]],  bd[["mu_hi"]]),  1] <- 0
  grad[outside(raw[[2]], bd[["k_lo"]],   bd[["k_hi"]]),   2] <- 0
  grad[outside(raw[[3]], bd[["rho_lo"]], bd[["rho_hi"]]), 3] <- 0
  list(value = val, grad = grad)
}

## log f(0) = kappa log H(1) = kappa log1p(-S(1)) for DEGPD model 1, with the GPD
## survival S(1) = (1 + xi / sigma)^(-1 / xi) on the log scale (S(1) = 0 when the
## support ends at or below 1). Where xi / sigma is subnormal (a log-link xi below
## ~2e-308 sigma) it has lost its relative precision, and log1p(z) / xi -- a ratio of
## differently rounded subnormals -- was wrong by up to O(1) (the objective up to 50
## below its minimum); there log S(1) = -1 / sigma, the exponential limit, to double
## precision.
.zt_degpd1_logf0 <- function(sigma, xi, kappa) {
  z <- pmax(xi / sigma, -1)
  logS <- ifelse(xi == 0 | abs(z) < .Machine$double.xmin, -1 / sigma, -log1p(z) / xi)
  kappa * .log1mexp(logS)
}

## log f(0) per observation (NA where the parameters are invalid)
.zt_logf0 <- function(eta, family, likdata = NULL) {
  v <- switch(family,
    gwaring  = .zt_gwaring_logf0(eta),
    degpd1   = .zt_degpd1_logf0(exp(eta[[1]]), .zt_degpd1_xi(eta[[2]], likdata),
                                exp(eta[[3]])),
    degpd1id = .zt_degpd1_logf0(exp(eta[[1]]), eta[[2]], exp(eta[[3]])),
    stop("no zero-truncated form for family '", family, "'", call. = FALSE))
  v[is.nan(v)] <- NA_real_
  v
}

## DEGPD rows at the kappa -> 0 limit of the truncation term: where kappa (log kappa
## below -708.4) or log f(0) = kappa log H(1) is subnormal it keeps only a few
## significant bits while log kappa is exact (the objective varied by O(1) there and came
## out up to 0.7 below its minimum), and below |log f(0)| = sqrt(DBL_MIN) the Hessian
## weight f0 / (1 - f0)^2 ~ 1 / log(f0)^2 overflows (the Hessian came out zero from log
## kappa ~ -355 down). For them log(1 - f(0)) = log(kappa) + log(-log H(1)) + O(log f(0))
## is taken from the predictor, as is log f(y) in the native cell
## (log_power_sf_from_pair), and the derivatives are those of this limit.
## NULL if there are none (or for the Waring).
.zt_kappa_limit <- function(eta, family, lf0, likdata = NULL) {
  lH1 <- switch(family,
    degpd1   = .zt_degpd1_logf0(exp(eta[[1]]), .zt_degpd1_xi(eta[[2]], likdata), 1),
    degpd1id = .zt_degpd1_logf0(exp(eta[[1]]), eta[[2]], 1),
    NULL)
  if (is.null(lH1)) return(NULL)
  rows <- which(!is.na(lf0) & !is.na(lH1) & lH1 < 0 &
                (exp(eta[[3]]) < .Machine$double.xmin | abs(lf0) < sqrt(.Machine$double.xmin)))
  if (!length(rows)) return(NULL)
  list(rows = rows, lH1 = lH1[rows])
}

.zt_d0 <- function(pars, likdata, base_d0, family) {
  ld <- likdata
  if (identical(ld$bounded.xi, FALSE)) ld$xi.max <- Inf
  nll <- base_d0(pars, ld)
  if (!is.finite(nll)) return(nll)
  eta <- .zt_eta(pars, ld)
  lf0 <- .zt_logf0(eta, family, ld)
  sb <- .zt_kappa_limit(eta, family, lf0, ld)
  ## f(0) = 1 would leave no mass for the positive counts
  bad <- is.na(lf0) | lf0 >= 0
  if (!is.null(sb)) bad[sb$rows] <- FALSE
  if (any(bad)) return(Inf)
  l1m <- log(-expm1(lf0))
  if (!is.null(sb)) l1m[sb$rows] <- eta[[3]][sb$rows] + log(-sb$lH1)
  nll + sum(l1m)
}

## Per-observation pieces shared by the gradient and Hessian: the base derivatives at
## y and at 0, log f(0), the weight w = f0 / (1 - f0) and c2 = f0 / (1 - f0)^2.
## For the Waring, the derivatives at 0 are the gradient of -log f0 computed above.
## For the rows of .zt_kappa_limit, g01 holds the base derivatives at 0 with
## kappa = 1, i.e. of -log H(1).
.zt_parts <- function(pars, ld, base_d12, family) {
  g  <- base_d12(pars, ld)
  eta <- .zt_eta(pars, ld)
  sb <- NULL
  if (family == "gwaring") {
    lf <- .zt_gwaring_logf0(eta, deriv = TRUE)
    lf0 <- lf$value
    g0 <- -lf$grad
  } else {
    ld0 <- ld
    ld0$y <- matrix(0, nrow(ld$y), ncol(ld$y))
    g0 <- base_d12(pars, ld0)
    lf0 <- .zt_logf0(eta, family, ld)
    sb <- .zt_kappa_limit(eta, family, lf0, ld)
    if (!is.null(sb)) {
      off <- ld0$offsets[[3]]
      if (!length(off)) off <- numeric(length(eta[[3]]))
      off[sb$rows] <- off[sb$rows] - eta[[3]][sb$rows]
      ld0$offsets[[3]] <- off
      sb$g01 <- base_d12(pars, ld0)[sb$rows, , drop = FALSE]
    }
  }
  ok <- !is.na(lf0) & lf0 < 0
  w <- ifelse(ok, 1 / expm1(-lf0), 0)
  list(np = length(ld$X), g = g, g0 = g0, ok = ok, w = w,
       c2 = ifelse(ok, w / -expm1(lf0), 0), sb = sb)
}

.zt_d12 <- function(pars, likdata, base_d12, family) {
  ld <- likdata
  if (identical(ld$bounded.xi, FALSE)) ld$xi.max <- Inf
  ## Note: for degpd1 with bounded.xi = TRUE (xi.max < Inf), .degpd1.d12()
  ## already applies the bounded-link chain factor mp_xi = (dxi/deta2)/xi to
  ## both g (at y) and g0 (at 0) inside degpd1d12().
  pt <- .zt_parts(pars, ld, base_d12, family)
  np <- pt$np
  G <- pt$g[, seq_len(np), drop = FALSE] + pt$w * pt$g0[, seq_len(np), drop = FALSE]
  exact <- family != "gwaring"
  c2 <- pt$c2
  out <- matrix(0, nrow(G), np + np * (np + 1) / 2)
  out[, seq_len(np)] <- G
  col <- np
  for (i in seq_len(np)) for (k in i:np) {
    col <- col + 1
    out[, col] <- if (exact) {
      pt$g[, col] + pt$w * pt$g0[, col] - c2 * pt$g0[, i] * pt$g0[, k]
    } else {
      G[, i] * G[, k]
    }
  }
  ## kappa -> 0 limit rows (.zt_kappa_limit): log(1 - f0) = log(kappa) + log(L), L = -log H(1),
  ## whose derivatives are 1 in log kappa and those of log(L) in the other predictors
  if (!is.null(pt$sb)) {
    s <- pt$sb$rows; L <- -pt$sb$lH1; g1 <- pt$sb$g01; ik <- 3L
    out[s, seq_len(np)] <- pt$g[s, seq_len(np), drop = FALSE] +
      vapply(seq_len(np), function(i) if (i == ik) rep(1, length(s)) else g1[, i] / L,
             numeric(length(s)))
    col <- np
    for (i in seq_len(np)) for (k in i:np) {
      col <- col + 1
      out[s, col] <- pt$g[s, col] +
        (if (i == ik || k == ik) 0 else g1[, col] / L - g1[, i] * g1[, k] / L^2)
    }
  }
  ## rows where the objective is infinite carry zero derivatives (the native
  ## convention; the optimiser rejects such steps on the objective)
  out[!is.finite(out)] <- 0
  out
}

.zt_fns <- function(base_d0, base_d12, family, iG, identity.xi = FALSE, inits = NULL) {
  res <- list(d0   = function(pars, likdata) .zt_d0(pars, likdata, base_d0, family),
              d120 = function(pars, likdata) .zt_d12(pars, likdata, base_d12, family),
              d340 = NULL, m = 1, iG = iG)
  if (isTRUE(identity.xi)) res$identity.xi <- TRUE
  if (!is.null(inits)) res$inits <- inits
  res
}

## The default start of egpd(family = "ztgwaring") (likfns$inits). .setup.inner.inits()
## takes Newton steps on the intercept-only model from the generic count start (the mean,
## k = 1, rho = 2.5). On 222 positive generalized Waring draws (mu = 4, k = 3, rho = 1.6)
## they ran past rho = 1 into the clamp, and the fit ended with convergence 3, 4.27
## log-likelihood units below the interior MLE (rho = 1.59) that fit_hurdle() and `inits`
## reach. So those steps are taken here first, and if they end outside the interior while
## the intercept-only MLE of the zero-truncated law lies inside it, the MLE is the start.
## Otherwise nothing changes: the generic start is returned -- or, if the steps ended by
## the gradient criterion, their end point, from which .setup.inner.inits() stops at once
## at the same point -- so those fits, and the start of fit_hurdle() that is egpd()'s
## default, are as before. (Steps that end inside, short of the MLE, are left to the fit:
## on data with rho = 1.8 they ended 0.05 log-likelihood units short, and the smooth fit
## that continued from there had the best REML criterion of fit_hurdle()'s starts, while
## the one started at the MLE failed.)
## - Interior: inside the box in which the likelihood is evaluated (cf_bounds_cpp()) and
##   off the plateaus of the links by the bound of .plateau_flags(), |log k| and
##   |log rho| <= 9, and log(rho - 1) >= -9. Where the likelihood rises towards rho -> 1
##   the MLE is not interior (a fit started at it stopped at rho = 1 + 1.001e-6, just
##   inside the clamp, escaping its convergence-3 report), and the generic start stays.
## - The MLE is found as fit_hurdle()'s first start is (.hurdle_zt_mle: Nelder-Mead from
##   the same point, with rho = 1 + exp(t) > 1 built in, ignoring offsets and `fixed`),
##   here on the distinct counts weighted by their frequencies, and is then polished by
##   BFGS on the intercept-only model's own objective and gradient: from Nelder-Mead's
##   |gradient| of about 1e-3, Newton's steps with the outer-product Hessian, nearly
##   singular along one direction there, crawled, and the fit ran out of iterations.
.ztgwaring_inits <- function(likdata) {
  y <- likdata$y[, 1]
  ld0 <- .start_model(likdata)
  start <- .start_coefficients(c(log(max(mean(y), 1e-3)), 0, log(2.5)), ld0)
  npar <- length(likdata$X)
  pinned <- vapply(likdata$X, function(X) !ncol(X) || .is_pinned_block(X), logical(1))
  bd <- cf_bounds_cpp()[["gwaring"]]
  interior <- function(b) {
    v <- exp(b)
    length(v) == 3L && all(is.finite(v)) &&
      isTRUE(all((v > c(bd[["mu_lo"]], bd[["k_lo"]], bd[["rho_lo"]]) &
                  v < c(bd[["mu_hi"]], bd[["k_hi"]], bd[["rho_hi"]]) &
                  c(TRUE, abs(b[2:3]) <= 9) &
                  c(TRUE, TRUE, b[3] > 0 && log(expm1(b[3])) >= -9))[!pinned]))
  }
  ## the estimate, on the distinct counts u with their frequencies w
  u <- sort(unique(y))
  w <- tabulate(match(y, u), length(u))
  ldu <- likdata
  ldu$y <- matrix(u, ncol = 1L)
  ldu$X <- rep(list(matrix(1, length(u), 1L)), npar)
  ldu$offsets <- rep(list(numeric(0)), npar)
  ldu$idpars <- seq_len(npar)
  ldu$dupid <- 0
  ldu$duplicate <- 0
  fnu <- function(b) {
    eta <- .zt_eta(b, ldu)
    v <- .zt_safe(-sum(w * (dgwaring(u, exp(eta[[1]]), exp(eta[[2]]), exp(eta[[3]]),
                                     log = TRUE) -
                            log(-expm1(.zt_gwaring_logf0(eta))))))
    if (length(v) == 1L && is.finite(v)) v else 1e20
  }
  gru <- function(b) colSums(w * .zt_d12(b, ldu, .gwaring1.d12, "gwaring")[, 1:3, drop = FALSE])
  nm <- tryCatch(suppressWarnings(stats::optim(
          c(log(mean(y)), 0, log(1.5)),
          function(t) fnu(c(t[1:2], if (t[3] > 30) t[3] else log1p(exp(t[3])))),
          control = list(maxit = 2000, reltol = 1e-10))), error = function(e) NULL)
  if (is.null(nm) || nm$value >= 1e20)
    return(start)
  mle <- c(nm$par[1:2], if (nm$par[3] > 30) nm$par[3] else log1p(exp(nm$par[3])))
  if (!interior(mle))
    return(start)
  ## the intercept-only model of .setup.inner.inits(), and its Newton steps from `start`
  interior_fit <- function(b) {
    eta <- .zt_eta(b, ld0)
    interior(vapply(eta, min, 0)) && interior(vapply(eta, max, 0))
  }
  pre <- tryCatch(.newton_step_inner(replace(start, pinned, 0), .nllh.nopen, .search.nopen,
                                     likdata = ld0, likfns = .ztgwaring1fns,
                                     control = likdata$control$inner),
                  error = function(e) NULL)
  if (!is.null(pre) && interior_fit(pre$pars))
    return(if (isTRUE(pre$gradconv)) pre$pars else start)
  ## polish the estimate, on u and w unless offsets or `fixed` make the rows differ
  nll0 <- function(b) .nllh.nopen(b, ld0, .ztgwaring1fns)
  b0 <- .start_coefficients(mle, ld0)
  plain <- all(vapply(ld0$X, function(z) all(z == 1), logical(1))) &&
    all(vapply(likdata$offsets, function(o) all(o == 0), logical(1)))
  o <- tryCatch(suppressWarnings(
         if (plain) stats::optim(b0, fnu, gru, method = "BFGS",
                                 control = list(maxit = 500, reltol = 1e-15))
         else stats::optim(b0, nll0,
                           function(b) as.vector(.gH.nopen(b, ld0, .ztgwaring1fns, deriv = 1)[[1]]),
                           method = "BFGS", control = list(maxit = 500, reltol = 1e-15))),
         error = function(e) NULL)
  if (!is.null(o) && interior_fit(o$par) && nll0(o$par) <= nll0(b0)) o$par else b0
}

.ztgwaring1fns <- .zt_fns(.gwaring1.d0,   .gwaring1.d12,   "gwaring",  .iG_gwaring,
                          inits = .ztgwaring_inits)
.ztdegpd1fns   <- .zt_fns(.degpd1.d0,     .degpd1.d12,     "degpd1",   .iG1_degpd)
.ztdegpd1idfns <- .zt_fns(.degpd1id.d0,   .degpd1id.d12,   "degpd1id", .iG1_degpd,
                          identity.xi = TRUE)

#' Fit a Hurdle (Zero-Adjusted) Count GAM
#'
#' @description Fits a hurdle model in which the probability of a zero has its own
#'   linear predictor, and the positive counts follow a zero-truncated generalized
#'   Waring or discrete EGPD (model 1). Unlike zero inflation, the zero probability
#'   is free in both directions, so it can correct over- as well as under-prediction
#'   of zeros.
#'
#' @details The hurdle log-likelihood factorises into a Bernoulli block for
#'   zero-versus-positive and a zero-truncated block for the positive counts, and the
#'   two share no parameters. They are therefore fitted separately -- the first with
#'   \code{mgcv::gam(family = binomial)}, the second with \code{\link{egpd}} on the
#'   subset \code{y > 0} -- and their log-likelihoods and degrees of freedom add.
#'   Complete cases across the response and all covariates referenced in
#'   \code{formula} and \code{zero.formula} are selected before splitting, with a
#'   warning when incomplete rows are omitted. This includes missing values
#'   introduced by evaluated terms and offsets, not just missing raw columns.
#'   A response expression is evaluated
#'   once before this selection. If \code{maxdata} truncates the data, both parts,
#'   all starts, the likelihood and \code{nobs()} refer to the same sampled rows.
#'
#'   The zero-truncated fit is started from several points and the best converged fit
#'   is kept, compared as for the restarts of \code{\link{egpd}} (see its argument
#'   \code{restarts}): by the log-likelihood, or, with smooth terms, by the REML
#'   criterion when the smoothing parameters are estimated (the penalised
#'   log-likelihood when \code{sp} is supplied through \code{...}), with fits within
#'   0.01 of the best counted as the same optimum and the highest log-likelihood
#'   kept among those; with a sandwich adjustment (\code{sandwich.args} through
#'   \code{...}), always by the log-likelihood. The starts are the intercept-only
#'   maximum likelihood estimate of the zero-truncated base law (computed directly from
#'   \code{\link{dgwaring}} or \code{\link{ddiscegpd}}, slopes and smooths starting
#'   at zero), a fixed start (for DEGPD, \eqn{\xi = 0.05} on the scale of the
#'   chosen shape link), and the \code{\link{egpd}} default. Supplying \code{inits}
#'   through \code{...} replaces these with that single start. A start counts only
#'   if \code{egpd()} reports convergence and a log-likelihood above the
#'   \eqn{-10^{20}} failure sentinel. If no start qualifies, the start that did not
#'   error and has the highest unpenalised log-likelihood at its returned
#'   coefficients is kept, with one warning giving its convergence code and report;
#'   as for a non-converged \code{egpd()} fit, the hurdle log-likelihood and AIC are
#'   then \code{NA} and \code{converged} is \code{FALSE}. \code{fit_hurdle()} stops
#'   only if every start errors or none has a finite log-likelihood at its
#'   coefficients, listing the reason for each start. Warnings and messages are
#'   reported for the selected start only.
#'
#'   For the zero-truncated DEGPD the Hessian of the positive part, and hence its
#'   standard errors, is the exact observed information. The zero-truncated
#'   generalized Waring, like \code{egpd(family = "gwaring")} and the other count
#'   families, uses the outer product of the per-observation scores (BHHH), which
#'   stays positive semi-definite where the observed information of its weakly
#'   identified \eqn{k} and \eqn{\rho} need not; its standard errors are
#'   BHHH-based, not observed-information standard errors.
#'   The outer product is singular where the scores are linearly dependent -- e.g. an
#'   intercept-only fit whose estimate has \eqn{a = k}, where
#'   \eqn{a = \mu(\rho - 1)/k} -- and then the standard errors of the positive part do
#'   not reflect the uncertainty along the unidentified direction. Even standard
#'   errors for the mean and tail parameters can differ substantially from those
#'   based on observed information near this singularity.
#'   \code{positive$hessian.rank} records the rank, \code{print()} and \code{summary()}
#'   note it, and the warning of class \code{"egpd_singular_hessian"} that
#'   \code{\link{egpd}} gives (see its Value) is reported for the start that is kept.
#'
#' @param formula model formula, or list of formulae, for the positive counts,
#'   exactly as passed to \code{\link{egpd}}: the response may be on any element of a
#'   named list, and on the first element of an unnamed one
#' @param data a data frame
#' @param family base family for the positive counts: \code{"gwaring"} or \code{"degpd"}
#' @param zero.formula one-sided formula for the probability of a positive count;
#'   defaults to the terms of the element of \code{formula} that has the
#'   response, including its offsets. Supply \code{zero.formula} explicitly
#'   if that offset should not enter the log odds of a positive count.
#' @param degpd.args list of arguments for the zero-truncated DEGPD model, as in
#'   \code{\link{egpd}}. Only \code{m = 1} is supported; \code{link = "identity"}
#'   gives the unconstrained shape
#' @param ... further arguments passed to \code{\link{egpd}} (unused by the
#'   \code{print}, \code{logLik} and \code{nobs} methods)
#' @param maxdata maximum number of complete rows for the joint model. A smaller
#'   sample is drawn without replacement before splitting off the positive counts.
#' @param x,object an object of class \code{"hurdle.egpd"}
#'
#' @return \code{fit_hurdle} returns an object of class \code{"hurdle.egpd"}, a
#'   list containing:
#'   \describe{
#'     \item{zero}{the fitted binomial \code{mgcv::gam} object for \eqn{P(Y > 0)}}
#'     \item{positive}{the fitted zero-truncated \code{\link{egpd}} object for
#'       \eqn{Y \mid Y > 0}}
#'     \item{family}{character string naming the positive-count family}
#'     \item{logLik}{combined log-likelihood (\code{NA} when no start of the
#'       zero-truncated fit converged)}
#'     \item{df}{combined effective degrees of freedom: those of the binomial part
#'       plus the positive part's \code{df} (see \code{\link{egpd}}, Value)}
#'     \item{AIC}{Akaike information criterion (\code{NA} when \code{logLik} is)}
#'     \item{converged}{logical; whether a start of the zero-truncated fit
#'       converged}
#'     \item{n}{integer number of complete observations used}
#'     \item{n_zero}{integer number of zero counts among the complete
#'       observations}
#'     \item{response.name}{name of the response variable}
#'     \item{data}{the complete-case data (not stored when \code{removeData = TRUE}
#'       is passed)}
#'   }
#'   Methods are provided for \code{print}, \code{summary}
#'   (\code{\link{summary.hurdle.egpd}}), \code{logLik}, \code{nobs} and
#'   \code{predict} (\code{\link{predict.hurdle.egpd}}).
#'
#'   \code{logLik.hurdle.egpd} returns an object of class \code{"logLik"} with
#'   \code{"df"} and \code{"nobs"} attributes suitable for
#'   \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}};
#'   \code{nobs.hurdle.egpd} returns \code{n}.
#'
#' @seealso \code{\link{egpd}}, \code{\link{summary.hurdle.egpd}},
#'   \code{\link{predict.hurdle.egpd}}
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(100), y = rpois(100, 2))
#' m <- fit_hurdle(y ~ x, data = d, family = "degpd")
#' print(m)
#' AIC(m)
#' BIC(m)
#' summary(m)
#' predict(m, newdata = d[1:3, ], type = "quantile", prob = c(0.5, 0.9))
#' @export
fit_hurdle <- function(formula, data, family = c("gwaring", "degpd"),
                       zero.formula = NULL, degpd.args = list(), ..., maxdata = Inf) {
  family <- match.arg(family)
  if (!is.list(degpd.args))
    stop("`degpd.args` must be a list", call. = FALSE)
  if (!is.null(degpd.args$m)) {
    m_val <- degpd.args$m
    if (!is.numeric(m_val) || is.logical(m_val) || length(m_val) != 1L ||
        !is.finite(m_val) || m_val != floor(m_val) || as.integer(m_val) != 1L) {
      stop("fit_hurdle() supports only degpd.args$m == 1", call. = FALSE)
    }
  }
  if (!is.data.frame(data))
    stop("`data` must be a data frame", call. = FALSE)
  if (!is.numeric(maxdata) || length(maxdata) != 1L || is.na(maxdata) ||
      maxdata < 1 || maxdata != floor(maxdata))
    stop("`maxdata` must be a positive integer or Inf", call. = FALSE)

  fl <- if (inherits(formula, "formula")) list(formula) else formula
  if (!is.list(fl) || !length(fl) || !all(vapply(fl, inherits, "formula", FUN.VALUE = logical(1))))
    stop("`formula` must be a formula or non-empty list of formulae", call. = FALSE)
  ## As in egpd(), any element of a named list may carry the response; an unnamed
  ## list is positional, so there it must be the first.
  has_resp <- vapply(fl, function(f) attr(stats::terms(f), "response") == 1L, logical(1))
  named <- !is.null(names(fl)) && any(nzchar(names(fl)))
  if (!any(has_resp) || (!named && !has_resp[1L]))
    stop("`formula` must include a response variable", call. = FALSE)
  material <- .materialise_response(fl, data, named = named)
  fl <- material$formula
  data <- material$data
  formula <- if (inherits(formula, "formula")) fl[[1L]] else fl
  rf <- fl[[which(has_resp)[1L]]]
  resp <- if (is.symbol(rf[[2L]])) as.character(rf[[2L]]) else character()
  if (length(resp) != 1L || !(resp %in% names(data)))
    stop("response '", paste(resp, collapse = ""), "' not found in `data`", call. = FALSE)

  if (is.null(zero.formula)) {
    zero.formula <- rf
    zero.formula[[2L]] <- NULL  # keep the RHS, offsets and formula environment
  }
  if (!inherits(zero.formula, "formula"))
    stop("`zero.formula` must be a formula or NULL", call. = FALSE)

  y_raw <- data[[resp]]
  if (!is.numeric(y_raw) || is.logical(y_raw))
    stop("response must consist of non-negative integers", call. = FALSE)
  y_ok <- y_raw[!is.na(y_raw)]
  if (any(!is.finite(y_ok)) || any(y_ok < 0) || any(y_ok != floor(y_ok)))
    stop("response must consist of non-negative integers", call. = FALSE)

  cc <- .gam_complete_rows(c(fl, list(zero.formula)), data, required = resp)
  if (any(!cc)) {
    warning(sprintf("omitted %d observation(s) with missing values", sum(!cc)),
            call. = FALSE)
    data <- data[cc, , drop = FALSE]
  }
  dots <- list(...)
  if (nrow(data) > maxdata) {
    data <- data[sort(sample.int(nrow(data), maxdata)), , drop = FALSE]
    if (is.null(dots$trace) || dots$trace >= 0)
      message("Hurdle model using one sample of ", maxdata, " rows for both parts (`maxdata`).")
  }
  y <- data[[resp]]
  if (length(y) == 0L || all(y == 0))
    stop("no positive counts available to fit the zero-truncated model", call. = FALSE)

  ## the hurdle itself: an ordinary binomial GAM for P(Y > 0)
  zdat <- data
  lhs <- .aux_response_name(".positive", c(names(zdat), all.vars(zero.formula)))
  zdat[[lhs]] <- as.integer(y > 0)
  zform <- stats::as.formula(call("~", as.name(lhs), zero.formula[[length(zero.formula)]]),
                             env = environment(zero.formula))
  zform <- .gam_offset_formula(zform, zdat)
  zdat <- .gam_offset_data(zform, zdat)
  zfit <- mgcv::gam(zform,
                    family = stats::binomial(), data = zdat, method = "REML")

  ## The positive counts: the same base family, zero-truncated, fitted from several
  ## starts (see Details). Every candidate is the same model fitted to the same data, so
  ## keeping the best converged fit is plain multi-start optimisation. DEGPD model 1 has
  ## a nearly flat ridge along which sigma -> 0 and kappa -> infinity, and from a poor
  ## start the optimiser can land on it (one test panel reached kappa = 7e7, 2300
  ## log-likelihood units below a sane fit) or slide to xi -> 0.
  ztfam  <- if (family == "gwaring") "ztgwaring" else "ztdegpd"
  dpos   <- data[y > 0, , drop = FALSE]
  ypos   <- y[y > 0]
  if ("inits" %in% names(dots)) {
    starts <- list(dots$inits)
    dots$inits <- NULL
  } else {
    starts <- lapply(.hurdle_starts(ypos, family, degpd.args), .eta_start)
  }

  rng <- .rng_freeze(!is.null(dots$maxspline) && nrow(dpos) > dots$maxspline)
  runs <- lapply(starts, function(s) {
    .rng_restore(rng)
    .hurdle_try_start(c(list(formula = formula, data = dpos, family = ztfam,
                             degpd.args = degpd.args, inits = s, maxdata = Inf), dots))
  })
  ## Prefer a genuine fit, the best of them by the criterion the fits optimised, as for
  ## egpd()'s restarts (.pick_start: with smooth terms, the REML criterion rather than
  ## the log-likelihood, unless a sandwich adjustment is in force). Failing that, keep
  ## the start that did not error and has the highest unpenalised log-likelihood at its
  ## returned coefficients (its logLik is NA, as egpd() reports for a non-converged
  ## fit), with one warning.
  ok <- vapply(runs, function(r) is.null(r$reason), logical(1))
  converged <- any(ok)
  if (converged) {
    best <- .pick_start(lapply(runs, `[[`, "fit"), which(ok))
  } else {
    ll_coef <- vapply(runs, function(r) r$ll_coef, numeric(1))
    if (all(is.na(ll_coef))) {
      reasons <- vapply(runs, `[[`, character(1), "reason")
      stop("the zero-truncated ", family, " fit failed from every start:\n",
           paste0("  start ", seq_along(reasons), ": ", reasons, collapse = "\n"),
           call. = FALSE)
    }
    best <- which.max(ll_coef)
  }
  for (cnd in runs[[best]]$conditions) {
    if (inherits(cnd, "warning")) warning(cnd) else message(cnd)
  }
  pfit <- runs[[best]]$fit
  if (nobs(zfit) != nrow(data) || nobs(pfit) != nrow(dpos))
    stop("The hurdle components did not retain the same joint sample; ",
         "check for missing values in the model terms.", call. = FALSE)
  if (!converged) {
    crep <- if (length(pfit$convergence.report))
      paste0(": ", paste(as.character(unlist(pfit$convergence.report)), collapse = "; "))
    warning(sprintf(paste0(
      "no start of the zero-truncated %s fit converged; returning the start with ",
      "the highest log-likelihood at its coefficients (%.3f; convergence code %s%s). ",
      "Its log-likelihood and AIC are NA. Consider a simpler model for the positive ",
      "counts, or supply `inits`."),
      family, ll_coef[best], format(pfit$convergence), if (is.null(crep)) "" else crep),
      call. = FALSE)
  }

  ll <- if (converged) as.numeric(stats::logLik(zfit)) + as.numeric(pfit$logLik) else NA_real_
  df <- sum(zfit$edf) + as.numeric(pfit$df)
  out <- list(zero = zfit, positive = pfit, family = family,
              logLik = ll, df = df, AIC = -2 * ll + 2 * df, converged = converged,
              n = as.integer(length(y)), n_zero = as.integer(sum(y == 0)),
              response.name = resp)
  if (!isTRUE(dots$removeData)) out$data <- data
  class(out) <- "hurdle.egpd"
  out
}

## Fit the zero-truncated block from one start, collecting its warnings (including
## egpd()'s "egpd_nonconvergence" warnings) and messages; fit_hurdle() re-emits only
## those of the start it keeps. `reason` is NULL for a genuine fit: convergence == 0
## and a finite log-likelihood above the -1e20 sentinel that .nllh.nopen() returns when
## the objective is infinite (finite, so a finiteness check alone would accept a fit
## that never moved). egpd() reports a failed fit with logLik = NA, a non-zero
## convergence code and the optimiser's message in `convergence.report`.
.hurdle_try_start <- function(args) {
  conds <- list()
  fit <- tryCatch(
    withCallingHandlers(
      do.call(egpd, args),
      warning = function(w) {
        conds[[length(conds) + 1L]] <<- w
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        conds[[length(conds) + 1L]] <<- m
        invokeRestart("muffleMessage")
      }),
    error = function(e) e)
  ## unpenalised log-likelihood at the returned coefficients, which ranks the starts
  ## when none converged; egpd() keeps likfns and likdata even with removeData = TRUE
  ll_coef <- if (inherits(fit, "error")) NA_real_ else
    .zt_safe(-fit$likfns$d0(unname(fit$coefficients), fit$likdata))
  if (length(ll_coef) != 1L || !is.finite(ll_coef) || ll_coef <= -1e19) ll_coef <- NA_real_
  reason <- if (inherits(fit, "error")) {
    paste("error:", conditionMessage(fit))
  } else if (!isTRUE(as.numeric(fit$convergence) == 0)) {
    paste0("did not converge (convergence code ", format(fit$convergence),
           if (length(fit$convergence.report))
             paste0(": ", paste(as.character(unlist(fit$convergence.report)),
                                collapse = "; ")), ")")
  } else if (!is.finite(as.numeric(fit$logLik)) || as.numeric(fit$logLik) <= -1e19) {
    paste0("no valid likelihood at the reported solution (log-likelihood ",
           format(as.numeric(fit$logLik)), ")")
  }
  list(fit = fit, reason = reason, conditions = conds, ll_coef = ll_coef)
}

## Candidate starts for the zero-truncated block, one value per linear predictor
## (egpd() expands each to the intercept, with any further coefficients at zero).
## The egpd() default (NULL) always comes last.
.hurdle_starts <- function(ypos, family, degpd.args) {
  if (family == "gwaring") {
    mle <- .hurdle_zt_mle(ypos, "gwaring")
    ## without the MLE, fall back to the gwaring seed of egpd(): mean, k = 1, rho = 2.5
    starts <- list(if (is.null(mle)) c(log(mean(ypos)), 0, log(2.5)) else mle)
  } else {
    link <- if (identical(degpd.args$link, "identity")) "identity" else "log"
    xi.max <- if (is.null(degpd.args$xi.max)) Inf else degpd.args$xi.max
    starts <- list(.hurdle_zt_mle(ypos, "degpd", link, xi.max),
                   c(log(mean(ypos) + 1),
                     .hurdle_xi_link(min(0.05, xi.max / 2), link, xi.max), 0))
  }
  c(Filter(Negate(is.null), starts), list(NULL))
}

## xi on the scale of the fitted shape link: identity, log, or the bounded log link
## xi = xi.max (1 - exp(-exp(eta) / xi.max)) of egpd(degpd.args = list(xi.max = )).
.hurdle_xi_link <- function(xi, link, xi.max = Inf) {
  if (identical(link, "identity")) return(xi)
  if (is.finite(xi.max)) return(log(-xi.max * log1p(-xi / xi.max)))
  log(xi)
}

## Intercept-only maximum likelihood estimate of the zero-truncated base law, on the
## link scale of the GAM, computed directly from the public pmf. NULL if it fails.
.hurdle_zt_mle <- function(ypos, family, link = "log", xi.max = Inf) {
  ypos <- ypos[is.finite(ypos) & ypos >= 1]
  if (length(ypos) < 2L) return(NULL)
  n <- length(ypos)
  zt_nll <- function(logp, logp0) {
    v <- -(sum(logp) - n * log1p(-exp(logp0)))
    if (length(v) == 1L && is.finite(v)) v else 1e10
  }
  best <- NULL
  fit_from <- function(st, nll) {
    o <- tryCatch(stats::optim(st, nll, control = list(maxit = 2000, reltol = 1e-10)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && o$value < 1e10 &&
        (is.null(best) || o$value < best$value)) best <<- o
  }
  if (family == "gwaring") {
    ## rho = 1 + exp(t3) keeps the mean finite, as the C++ clamp does
    nll <- function(th) {
      mu <- exp(th[1]); k <- exp(th[2]); rho <- 1 + exp(th[3])
      v <- .zt_safe(zt_nll(dgwaring(ypos, mu = mu, k = k, rho = rho, log = TRUE),
                           dgwaring(0, mu = mu, k = k, rho = rho, log = TRUE)))
      if (is.finite(v)) v else 1e10
    }
    fit_from(c(log(mean(ypos)), 0, log(1.5)), nll)
    if (is.null(best)) return(NULL)
    t3 <- best$par[3]
    return(c(best$par[1:2], if (t3 > 30) t3 else log1p(exp(t3))))   # log rho
  }
  ## zero-truncated DEGPD model 1, shape on the GAM's link scale
  xi_of <- function(t) {
    if (identical(link, "identity")) return(t)
    if (is.finite(xi.max)) return(-xi.max * expm1(-exp(t) / xi.max))
    exp(t)
  }
  nll <- function(th) {
    s <- exp(th[1]); xi <- xi_of(th[2]); k <- exp(th[3])
    v <- .zt_safe(zt_nll(ddiscegpd(ypos, sigma = s, xi = xi, kappa = k, type = 1, log = TRUE),
                         ddiscegpd(0, sigma = s, xi = xi, kappa = k, type = 1, log = TRUE)))
    if (is.finite(v)) v else 1e10
  }
  for (xi0 in unique(pmin(c(0.05, 0.3), xi.max / 2)))
    fit_from(c(log(mean(ypos)), .hurdle_xi_link(xi0, link, xi.max), 0), nll)
  if (is.null(best)) NULL else best$par
}

## Quantiles of a zero-truncated count law at per-row probabilities u:
## Q(u) = Q_base(f0 + u (1 - f0)), f0 = P_base(0), which is at least 1. `pars` holds the
## response-scale parameters by name, as returned by predict.egpd(). Where the upper
## tail of that linear target, (1 - u) (1 - f0), is below 1e-3 it has lost relative
## accuracy (it is exactly 1 once P_base(Y > 0) < 1e-16, giving Inf), so there Q inverts
## the base survivor instead: the smallest q with
## log S_base(q) <= log(1 - u) + log S_base(0). `l1mu` = log(1 - u), which a caller may
## know more accurately than log1p(-u).
.zt_quantile <- function(family, pars, u, l1mu = log1p(-u)) {
  out <- rep(NA_real_, length(u))
  bad <- !is.na(u) & (u < 0 | u > 1)
  out[bad] <- NaN
  ok <- !is.na(u) & !bad
  if (!any(ok)) return(out)
  if (family %in% c("gwaring", "ztgwaring")) {
    a <- list(mu = pars[["mu"]][ok], k = pars[["k"]][ok], rho = pars[["rho"]][ok])
    dfun <- dgwaring; pfun <- pgwaring; qfun <- qgwaring
  } else if (family %in% c("degpd", "ztdegpd")) {
    a <- list(sigma = pars[["scale"]][ok], xi = pars[["shape"]][ok],
              kappa = pars[["kappa"]][ok], type = 1)
    dfun <- ddiscegpd; pfun <- pdiscegpd; qfun <- qdiscegpd
  } else {
    stop("no zero-truncated quantile for family '", family, "'", call. = FALSE)
  }
  uo <- u[ok]
  lo <- rep_len(l1mu, length(u))[ok]
  f0 <- do.call(dfun, c(list(0), a))
  lin <- !(exp(lo) * (1 - f0) < 1e-3)
  lin[is.na(lin)] <- TRUE
  sub <- function(i) lapply(a, function(v) if (length(v) > 1L) v[i] else v)
  q <- rep(NA_real_, length(uo))
  if (any(lin))
    q[lin] <- do.call(qfun, c(list(f0[lin] + uo[lin] * (1 - f0[lin])), sub(lin)))
  if (any(!lin)) {
    s <- !lin
    ls0 <- do.call(pfun, c(list(0), sub(s), list(lower.tail = FALSE, log.p = TRUE)))
    q[s] <- do.call(qfun, c(list(lo[s] + ls0), sub(s),
                            list(lower.tail = FALSE, log.p = TRUE)))
  }
  out[ok] <- pmax(1, q)
  out
}

#' @rdname fit_hurdle
#' @export
print.hurdle.egpd <- function(x, ...) {
  cat("Hurdle model: binomial zero part + zero-truncated", x$family, "\n")
  cat(sprintf("n = %d (%d zeros, %.1f%%)\n", x$n, x$n_zero, 100 * x$n_zero / x$n))
  cat(sprintf("logLik = %.2f  df = %.2f  AIC = %.2f\n", x$logLik, x$df, x$AIC))
  if (isFALSE(x$converged))
    cat("The zero-truncated fit did not converge from any start.\n")
  note <- .hessian_note(x$positive$hessian.rank, "`positive$hessian.rank`")
  if (!is.null(note)) cat(note, "\n", sep = "")
  invisible(x)
}

#' @rdname fit_hurdle
#' @export
logLik.hurdle.egpd <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- object$df
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

#' @rdname fit_hurdle
#' @export
nobs.hurdle.egpd <- function(object, ...) object$n

#' Summarise a hurdle model
#'
#' Combines \code{summary()} of the binomial zero part (an \code{mgcv::gam}) and of
#' the zero-truncated positive part (an \code{\link{egpd}} fit) with the joint
#' log-likelihood, degrees of freedom and AIC.
#'
#' @param object a \code{"hurdle.egpd"} fit from \code{\link{fit_hurdle}}
#' @param x a \code{"summary.hurdle.egpd"} object
#' @param ... unused
#' @return An object of class \code{"summary.hurdle.egpd"} with components
#'   \code{zero} (\code{summary.gam}), \code{positive} (\code{summary.egpd}),
#'   \code{family}, \code{logLik}, \code{df}, \code{AIC}, \code{n} and \code{n_zero}.
#' @export
summary.hurdle.egpd <- function(object, ...) {
  structure(list(zero = summary(object$zero), positive = summary(object$positive),
                 family = object$family, logLik = object$logLik, df = object$df,
                 AIC = object$AIC, n = object$n, n_zero = object$n_zero),
            class = "summary.hurdle.egpd")
}

#' @rdname summary.hurdle.egpd
#' @export
print.summary.hurdle.egpd <- function(x, ...) {
  cat("Hurdle model: binomial zero part + zero-truncated", x$family, "\n")
  cat(sprintf("n = %d (%d zeros, %.1f%%)\n", x$n, x$n_zero, 100 * x$n_zero / x$n))
  cat("\n==== Zero part: binomial GAM for P(Y > 0) ====\n")
  print(x$zero)
  cat("\n==== Positive part: zero-truncated", x$family, "====\n")
  print(x$positive)
  cat(sprintf("\nJoint logLik = %.2f  df = %.2f  AIC = %.2f\n", x$logLik, x$df, x$AIC))
  invisible(x)
}

#' Predictions from a hurdle model
#'
#' The hurdle distribution is \eqn{P(Y = 0) = \pi_0} and
#' \eqn{P(Y = y) = (1 - \pi_0) f(y) / (1 - f(0))} for \eqn{y \ge 1}, where
#' \eqn{\pi_0} comes from the binomial zero part and \eqn{f} is the base law of the
#' positive part (generalized Waring or DEGPD model 1).
#'
#' @param object a \code{"hurdle.egpd"} fit from \code{\link{fit_hurdle}}
#' @param newdata a data frame holding the variables of both parts; defaults to the
#'   complete-case data used for the fit (not stored if it was fitted with
#'   \code{removeData = TRUE}). A parametric factor level absent from a
#'   component's training rows gives \code{NA} for that component, with a
#'   warning; predictions for other levels are unaffected.
#' @param type \code{"response"} (default): the probability of a zero, column
#'   \code{pzero}, followed by the response-scale parameters of the positive part
#'   as returned by \code{\link{predict.egpd}} (for DEGPD this includes the
#'   \code{xi} alias of \code{shape}); \code{"link"}: the zero-part linear
#'   predictor, column \code{logitppositive} (the logit of \eqn{P(Y > 0)}),
#'   followed by the positive part's linear predictors; \code{"quantile"}:
#'   quantiles of the full hurdle distribution at \code{prob}, i.e. 0 when
#'   \eqn{p \le \pi_0} and otherwise the zero-truncated quantile at
#'   \eqn{(p - \pi_0) / (1 - \pi_0)}. The zero-atom comparison uses the
#'   fitted log odds rather than a rounded response-scale probability.
#' @param prob probabilities for \code{type = "quantile"}; supplying \code{prob}
#'   implies \code{type = "quantile"}
#' @param ... unused
#' @return A data frame with one row per row of \code{newdata}; quantile columns are
#'   named \code{"q:<prob>"}. Select columns by name.
#' @export
predict.hurdle.egpd <- function(object, newdata, type = c("response", "link", "quantile"),
                                prob = NULL, ...) {
  type <- match.arg(type)
  if (!is.null(prob)) type <- "quantile"
  if (type == "quantile" && is.null(prob))
    stop("non-NULL `prob' required if `type = \"quantile\"'", call. = FALSE)
  if (missing(newdata) || is.null(newdata)) {
    newdata <- object$data
    if (is.null(newdata))
      stop("supply `newdata': the data were not stored (removeData = TRUE)", call. = FALSE)
  }
  eta0 <- as.numeric(.gam_predict_rows(object$zero, newdata, "link"))
  if (type == "link") {
    pos <- predict(object$positive, newdata = newdata, type = "link")
    return(data.frame(logitppositive = eta0, pos, check.names = FALSE))
  }
  pzero <- stats::plogis(eta0, lower.tail = FALSE)
  pars <- predict(object$positive, newdata = newdata, type = "response")
  if (type == "response")
    return(data.frame(pzero = pzero, pars, check.names = FALSE))
  logpos <- stats::plogis(eta0, log.p = TRUE)
  out <- vapply(prob, function(p) {
    if (is.na(p)) return(rep(NA_real_, length(pzero)))
    if (p < 0 || p > 1) return(rep(NaN, length(pzero)))
    if (p == 0) return(ifelse(is.na(eta0), NA_real_, 0))
    if (p == 1) return(ifelse(eta0 == -Inf, 0,
      .zt_quantile(object$family, pars, rep(1, length(eta0)), -Inf)))
    ## On the positive branch u = p * (1-exp(-(eta + logit(p)))).
    ## Decide the side before rounding either fitted probability.
    cut <- if (p < .25) log1p(-p) - log(p) else -qlogis(p)
    gap <- eta0 - cut
    near <- which(abs(gap) <= 8*.Machine$double.eps*(1+abs(cut)))
    if (p != .5 && length(near)) gap[near] <- hurdle_logodds_gap_cpp(eta0[near], p)
    u <- p * -expm1(-pmax(gap, 0))
    l1mu <- ifelse(u < .5, log1p(-u), pmin(0, log1p(-p) - logpos))
    ifelse(gap <= 0, 0, .zt_quantile(object$family, pars, u, l1mu))
  }, numeric(length(pzero)))
  out <- as.data.frame(matrix(out, nrow = length(pzero)))
  names(out) <- .egpd_qnames(prob)
  out
}
