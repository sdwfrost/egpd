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
## anything new. For that block the negative log-likelihood and its gradient follow
## from the base family's own d0/d12 evaluated twice -- at y, and at 0 -- with no
## new derivative algebra:
##
##   nll_j       = base_nll_j + log(1 - f0_j)
##   dnll_j/deta = base_g_j(y) + [f0_j / (1 - f0_j)] * base_g_j(0)
##
## The second-order block is the outer product of the gradient, matching cf_pack()
## in src/countfams.cpp.

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

## f(0) on the response scale, per observation. The optimiser visits parameter
## values the density functions reject outright -- dgwaring() requires rho > 1, for
## instance -- so failures return NA and are turned into an infinite negative
## log-likelihood by the caller rather than an error.
.zt_safe <- function(expr) {
  v <- tryCatch(suppressWarnings(expr), error = function(e) NA_real_)
  if (!is.numeric(v)) NA_real_ else v
}

.zt_f0 <- function(eta, family) {
  switch(family,
    gwaring  = mapply(function(a, b, c)
                 .zt_safe(dgwaring(0, mu = exp(a), k = exp(b), rho = exp(c))),
                 eta[[1]], eta[[2]], eta[[3]]),
    degpd1   = .zt_safe(dDEGPD1(0, mu = exp(eta[[1]]), sigma = exp(eta[[2]]),
                                nu = exp(eta[[3]]))),
    degpd1id = .zt_safe(dDEGPD1(0, mu = exp(eta[[1]]), sigma = eta[[2]],
                                nu = exp(eta[[3]]))),
    stop("no zero-truncated form for family '", family, "'", call. = FALSE))
}

.zt_d0 <- function(pars, likdata, base_d0, family) {
  nll <- base_d0(pars, likdata)
  if (!is.finite(nll)) return(nll)
  f0 <- .zt_f0(.zt_eta(pars, likdata), family)
  ## f(0) = 1 would leave no mass for the positive counts
  if (any(!is.finite(f0)) || any(f0 >= 1 - 1e-12) || any(f0 < 0)) return(Inf)
  nll + sum(log1p(-f0))
}

.zt_d12 <- function(pars, likdata, base_d12, family) {
  np <- length(likdata$X)
  g  <- base_d12(pars, likdata)
  ld0 <- likdata
  ld0$y <- matrix(0, nrow(likdata$y), ncol(likdata$y))
  g0 <- base_d12(pars, ld0)

  f0 <- .zt_f0(.zt_eta(pars, likdata), family)
  w  <- ifelse(is.finite(f0) & f0 >= 0 & f0 < 1 - 1e-12, f0 / (1 - f0), 0)

  G <- g[, seq_len(np), drop = FALSE] + w * g0[, seq_len(np), drop = FALSE]
  G[!is.finite(G)] <- 0

  out <- matrix(0, nrow(G), np + np * (np + 1) / 2)
  out[, seq_len(np)] <- G
  col <- np
  for (i in seq_len(np)) for (k in i:np) {
    col <- col + 1
    out[, col] <- G[, i] * G[, k]
  }
  out
}

.zt_fns <- function(base_d0, base_d12, family, iG) {
  list(d0   = function(pars, likdata) .zt_d0(pars, likdata, base_d0, family),
       d120 = function(pars, likdata) .zt_d12(pars, likdata, base_d12, family),
       d340 = NULL, m = 1, iG = iG)
}

.ztgwaring1fns <- .zt_fns(.gwaring1.d0,   .gwaring1.d12,   "gwaring",  .iG_gwaring)
.ztdegpd1fns   <- .zt_fns(.degpd1.d0,     .degpd1.d12,     "degpd1",   .iG1_degpd)
.ztdegpd1idfns <- .zt_fns(.degpd1id.d0,   .degpd1id.d12,   "degpd1id", .iG1_degpd)

#' Fit a Hurdle (Zero-Adjusted) Count GAM
#'
#' @description Fits a hurdle model in which the probability of a zero has its own
#'   linear predictor, and the positive counts follow a zero-truncated generalized
#'   Waring or discrete EGPD. Unlike zero inflation, the zero probability is free in
#'   both directions, so it can correct over- as well as under-prediction of zeros.
#'
#' @details The hurdle log-likelihood factorises into a Bernoulli block for
#'   zero-versus-positive and a zero-truncated block for the positive counts, and the
#'   two share no parameters. They are therefore fitted separately -- the first with
#'   \code{mgcv::gam(family = binomial)}, the second with \code{\link{egpd}} on the
#'   subset \code{y > 0} -- and their log-likelihoods and degrees of freedom add.
#'
#' @param formula model formula, or list of formulae, for the positive counts,
#'   exactly as passed to \code{\link{egpd}}
#' @param data a data frame
#' @param family base family for the positive counts: \code{"gwaring"} or \code{"degpd"}
#' @param zero.formula one-sided formula for the probability of a positive count;
#'   defaults to the terms of the first element of \code{formula}
#' @param degpd.args as in \code{\link{egpd}}; \code{link = "identity"} gives the
#'   unconstrained shape
#' @param ... further arguments passed to \code{\link{egpd}}
#'
#' @return An object of class \code{"hurdle.egpd"}: a list holding the binomial fit
#'   (\code{zero}), the zero-truncated fit (\code{positive}), and the combined
#'   \code{logLik}, \code{df} and \code{AIC}.
#'
#' @examples
#' \dontrun{
#' m <- fit_hurdle(y ~ s(t), data = d, family = "gwaring")
#' m$AIC
#' }
#' @export
fit_hurdle <- function(formula, data, family = c("gwaring", "degpd"),
                       zero.formula = NULL, degpd.args = list(), ...) {
  family <- match.arg(family)
  fl <- if (inherits(formula, "formula")) list(formula) else formula
  resp <- all.vars(fl[[1]])[1]
  y <- data[[resp]]
  if (is.null(y))
    stop("response '", resp, "' not found in `data`", call. = FALSE)
  if (any(y < 0, na.rm = TRUE) || any(y != floor(y), na.rm = TRUE))
    stop("response must be non-negative counts", call. = FALSE)

  ## the hurdle itself: an ordinary binomial GAM for P(Y > 0)
  if (is.null(zero.formula))
    zero.formula <- stats::reformulate(attr(stats::terms(fl[[1]]), "term.labels"))
  zdat <- data
  zdat$.positive <- as.integer(y > 0)
  zfit <- mgcv::gam(stats::update(zero.formula, .positive ~ .),
                    family = stats::binomial(), data = zdat, method = "REML")

  ## The positive counts: the same base family, zero-truncated. DEGPD model 1 has a
  ## nearly flat ridge along which sigma -> 0 and kappa -> infinity, and from the
  ## default start the optimiser can land on it -- on one test panel it reached
  ## kappa = 7e7, costing 2300 log-likelihood units against a sane fit. So try more
  ## than one start and keep whichever attains the highest likelihood; every
  ## candidate lies in the same parameter space, so this is plain multi-start MLE.
  ztfam  <- if (family == "gwaring") "ztgwaring" else "ztdegpd"
  dpos   <- data[y > 0, , drop = FALSE]
  starts <- list(NULL)
  if (family == "degpd")
    starts <- c(starts, list(c(log(mean(y[y > 0]) + 1), 0.05, 0)))

  fits <- lapply(starts, function(s)
    tryCatch(egpd(formula, data = dpos, family = ztfam, degpd.args = degpd.args,
                  inits = s, ...),
             error = function(e) NULL))
  fits <- Filter(function(z) !is.null(z) && is.finite(as.numeric(z$logLik)), fits)
  if (!length(fits))
    stop("the zero-truncated fit failed from every start", call. = FALSE)
  pfit <- fits[[which.max(vapply(fits, function(z) as.numeric(z$logLik), 0))]]

  ll <- as.numeric(stats::logLik(zfit)) + as.numeric(pfit$logLik)
  df <- sum(zfit$edf) + as.numeric(pfit$df)
  out <- list(zero = zfit, positive = pfit, family = family,
              logLik = ll, df = df, AIC = -2 * ll + 2 * df,
              n = length(y), n_zero = sum(y == 0))
  class(out) <- "hurdle.egpd"
  out
}

#' @export
print.hurdle.egpd <- function(x, ...) {
  cat("Hurdle model: binomial zero part + zero-truncated", x$family, "\n")
  cat(sprintf("n = %d (%d zeros, %.1f%%)\n", x$n, x$n_zero, 100 * x$n_zero / x$n))
  cat(sprintf("logLik = %.2f  df = %.2f  AIC = %.2f\n", x$logLik, x$df, x$AIC))
  invisible(x)
}

#' @export
logLik.hurdle.egpd <- function(object, ...) {
  structure(object$logLik, df = object$df, class = "logLik")
}
