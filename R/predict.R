## Inverse links shared by predictions, parameter draws, and confidence intervals.
## A custom family's own inverse link (custom.fns$unlink, by predictor) takes
## precedence; otherwise the link is the one named by the predictor's prefix.
.strip_link <- function(name) sub("^(cloglog|probit|logit|log)", "", name)

.egpd_invlink <- function(object, name, eta, deriv = FALSE) {
  if (identical(object$family, "custom")) {
    i <- match(name, names(object)[seq_len(object$ngam)])
    ul <- object$likfns$unlink
    fn <- if (!is.na(i) && i <= length(ul)) ul[[i]] else NULL
    if (!is.null(fn)) {
      if (!deriv) return(.egpd_apply_unlink(fn, eta, name))
      dfn <- attr(fn, "deriv")
      if (!is.null(dfn)) return(.egpd_apply_unlink(dfn, eta, name))
      ## no derivative supplied: a central difference, with a relative step near
      ## the optimum eps^(1/3) for its truncation and rounding errors
      h <- 6e-6 * pmax(abs(eta), 1)
      return((.egpd_apply_unlink(fn, eta + h, name) -
              .egpd_apply_unlink(fn, eta - h, name)) / (2 * h))
    }
  }
  if (identical(name, "logshape") && isTRUE(object$likdata$bounded.xi) &&
      !isTRUE(object$likfns$identity.xi)) {
    a <- object$likdata$xi.max
    t <- exp(eta) / a
    if (!deriv) return(-a * expm1(-t))
    out <- exp(eta - t)
    out[is.infinite(t) & !is.na(t)] <- 0
    return(out)
  }
  if (startsWith(name, "logit")) {
    p <- plogis(eta)
    return(if (deriv) p * (1 - p) else p)
  }
  if (startsWith(name, "probit"))
    return(if (deriv) dnorm(eta) else pnorm(eta))
  if (startsWith(name, "cloglog")) {
    if (!deriv) return(-expm1(-exp(eta)))
    out <- exp(eta - exp(eta))
    out[is.infinite(eta) & eta > 0] <- 0
    return(out)
  }
  if (startsWith(name, "log")) return(exp(eta))
  if (deriv) return(eta * 0 + 1)
  eta
}

## Apply a user-supplied (custom-family) inverse link, or its derivative, keeping the
## shape of `eta` (a vector, or the observation-by-draw matrix of simulate()).
.egpd_apply_unlink <- function(fn, eta, name) {
  out <- fn(eta)
  if (is.numeric(out) && length(out) == 1L) out <- rep(out, length(eta))
  if (!is.numeric(out) || length(out) != length(eta))
    stop("custom.fns$unlink for '", name, "' must return one number per value of ",
         "its linear predictor", call. = FALSE)
  attributes(out) <- attributes(eta)
  out
}

## Effective response-scale bounds used by the count likelihoods.
.egpd_count_bounds <- function(family) {
  switch(as.vector(family),
         gpig = , gpignat = , zigpig = , zigpignat =  # eta_to_abc in src/gpig.cpp
           c(a_lo = 1e-3, a_hi = 1 - 1e-6, c_lo = 1e-6, c_hi = 1 - 1e-6),
         cmp = cmp_bounds_cpp(), prl = prl_bounds_cpp(),
         gpois = , plnorm = , gwaring = cf_bounds_cpp()[[family]],
         ztgwaring = cf_bounds_cpp()[["gwaring"]], NULL)
}

.egpd_param_bounds <- function(bd, name, link = FALSE) {
  key <- paste0(.strip_link(name), c("_lo", "_hi"))
  if (!all(key %in% names(bd))) return(NULL)
  lim <- unname(bd[key])
  if (!link) return(lim)
  if (startsWith(name, "logit")) qlogis(lim) else log(lim)
}

## Response-scale parameters clamped to the safe interior that the likelihood itself
## uses, so that predict() and simulate() describe the model that was fitted. `pars`
## is any named list of vectors or matrices (a data frame from predict(), or the
## draw matrices of simulate()); entries are matched by name.
.egpd_clamp_response <- function(pars, family, bd = .egpd_count_bounds(family)) {
  for (pn in names(pars)) {
    lim <- .egpd_param_bounds(bd, pn)
    if (length(lim)) pars[[pn]] <- pmin(pmax(pars[[pn]], lim[1]), lim[2])
  }
  pars
}

## Arguments of the EGPD distribution functions (sigma, xi, kappa, delta, prob, type)
## from response-scale parameters, matched by name. Model 2's second power is
## kappa2 = kappa1 + dkappa.
.egpd_dist_args <- function(pars, m) {
  out <- list(sigma = pars[["scale"]], xi = pars[["shape"]],
              type = c(1L, 6L, 4L, 5L, 2L, 3L)[m])
  if (m %in% c(1L, 4L, 5L, 6L)) out$kappa <- pars[["kappa"]]
  if (m %in% c(3L, 4L)) out$delta <- pars[["delta"]]
  if (m == 2L) {
    out$kappa <- pars[["kappa1"]]
    out$delta <- pars[["kappa1"]] + pars[["dkappa"]]
    out$prob  <- pars[["p"]]
  }
  out
}

## Keep familiar short labels unless two distinct probabilities need more
## digits to distinguish them (e.g. the largest double below one and one).
.egpd_qnames <- function(prob) {
  if (!length(prob)) return(character(0))
  p <- unique(prob)
  lab <- vapply(p, format, character(1), digits = 15)
  dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
  if (any(dup)) lab[dup] <- vapply(p[dup], format, character(1), digits = 17)
  paste0("q:", lab[match(prob, p)])
}

## Quantiles at probabilities `prob` for response-scale parameters `pars` (a data
## frame or named list of equal-length vectors, by parameter name, as returned by
## predict(type = "response") before the `xi` alias is added). Returns a
## length(rows) x length(prob) matrix. Shared by predict() and simulate().
.egpd_quantile <- function(object, pars, prob) {
  family <- object$family
  pars <- as.data.frame(pars, optional = TRUE)
  nms <- names(pars)
  ndat <- nrow(pars)
  nprob <- length(prob)
  out <- matrix(NA_real_, ndat, nprob)
  if (family %in% c("custom", "comppareto")) q_fn <- object$likfns$q

  if (family %in% c("egpd", "degpd", "zidegpd", "ztdegpd"))
    q_args <- .egpd_dist_args(pars, object$likfns$m)

  for (j in seq_len(nprob)) {
    pj <- prob[j]

    if (family == "custom") {
      if (length(nms) == 1) {
        out[, j] <- q_fn(pj, pars[,1])
      } else if (length(nms) == 2) {
        out[, j] <- q_fn(pj, pars[,1], pars[,2])
      } else if (length(nms) == 3) {
        out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3])
      } else {
        out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3], pars[,4])
      }
    } else if (family == "comppareto") {
      out[, j] <- do.call(q_fn, c(list(pj), unname(as.list(pars[seq_along(nms)]))))
    } else if (family == "egpd") {
      out[, j] <- do.call(qegpd, c(list(p = rep(pj, ndat)), q_args))
    } else if (family == "degpd") {
      out[, j] <- do.call(qdiscegpd, c(list(p = rep(pj, ndat)), q_args))
    } else if (family == "zidegpd") {
      out[, j] <- do.call(qzidiscegpd, c(list(p = rep(pj, ndat), pi = pars[["pi"]]), q_args))
    } else if (family %in% c("ztdegpd", "ztgwaring")) {
      out[, j] <- .zt_quantile(family, pars, rep(pj, ndat))
    } else if (family == "gpd") {
      out[, j] <- .qgpd(pj, 0, pars[,1], pars[,2])
    } else if (family == "pig") {
      out[, j] <- qpig(rep(pj, ndat), mu = pars$mu, sigma = pars$sigma)
    } else if (family == "zipig") {
      out[, j] <- qzipig(rep(pj, ndat), mu = pars$mu, sigma = pars$sigma, pi = pars$pi)
    } else if (family == "gpig") {
      out[, j] <- qGPIG(rep(pj, ndat), mu = pars$mu, sigma = pars$a, nu = pars$c)
    } else if (family == "gpignat") {
      out[, j] <- qgpig(rep(pj, ndat), a = pars$a, b = pars$b, c = pars$c)
    } else if (family == "zigpig") {
      out[, j] <- qZIGPIG(rep(pj, ndat), mu = pars$mu, sigma = pars$a, nu = pars$c, tau = pars$pi)
    } else if (family == "zigpignat") {
      out[, j] <- qzigpig(rep(pj, ndat), a = pars$a, b = pars$b, c = pars$c, pi = pars$pi)
    } else if (family == "bell") {
      out[, j] <- qBELL(rep(pj, ndat), mu = pars$mu)
    } else if (family == "bellnat") {
      out[, j] <- qbell(rep(pj, ndat), theta = pars$theta)
    } else if (family == "zibell") {
      out[, j] <- qZIBELL(rep(pj, ndat), mu = pars$mu, sigma = pars$pi)
    } else if (family == "zibellnat") {
      out[, j] <- qzibell(rep(pj, ndat), theta = pars$theta, pi = pars$pi)
    } else if (family == "cmp") {
      out[, j] <- vapply(seq_len(ndat), function(i) qcmp(pj, mu = pars$mu[i], nu = pars$nu[i]), numeric(1))
    } else if (family == "prl") {
      out[, j] <- vapply(seq_len(ndat), function(i) qprl(pj, tau = .prl_tau_from_mu(pars$mu[i])), numeric(1))
    } else if (family == "gpois") {
      out[, j] <- vapply(seq_len(ndat), function(i) qgpois(pj, mu = pars$mu[i], lambda = pars$lambda[i]), numeric(1))
    } else if (family == "gwaring") {
      out[, j] <- vapply(seq_len(ndat), function(i) qgwaring(pj, mu = pars$mu[i], k = pars$k[i], rho = pars$rho[i]), numeric(1))
    } else if (family == "plnorm") {
      out[, j] <- vapply(seq_len(ndat), function(i) qpoislnorm(pj, mu = pars$mu[i], sigma = pars$sigma[i]), numeric(1))
    } else {
      stop("invalid family for quantile prediction")
    }
  }
  out
}

#' Predictions from a fitted \code{egpd} object
#'
#' @param object a fitted \code{egpd} object
#' @param newdata a data frame. A missing input to one linear predictor
#'   gives a missing value for that predictor on that row; other rows
#'   retain their own design values and offsets. Evaluated predictor terms
#'   must have one row per prediction row. Supply row-specific covariates
#'   and offsets as columns rather than vectors in the formula environment.
#' @param type a character string: "link", "response", "lpmatrix", or "quantile"
#' @param prob a scalar or vector of probabilities for quantile estimation
#' @param se.fit logical: should standard errors be returned? Defaults to FALSE.
#'   Response-scale delta-method SEs are zero on a flat part of a clamped
#'   inverse link. At the kink there is no derivative, so the SE is
#'   \code{NA} unless the linear predictor has zero uncertainty.
#' @param marginal logical: should uncertainty integrate out smoothing parameter uncertainty? Defaults to TRUE
#' @param trace an integer controlling output verbosity
#' @param ... unused
#'
#' @return Depending on \code{type}:
#'   \describe{
#'     \item{\code{"link"}}{a data frame with one column per linear predictor,
#'       named by its link (e.g. \code{logscale}, \code{logshape}), including
#'       offsets.}
#'     \item{\code{"response"}}{a data frame with one column per distributional
#'       parameter on its natural scale, named by parameter (e.g. \code{scale},
#'       \code{shape}, \code{kappa}, \code{pi}), clamped to the same safe interior
#'       as the likelihood. For the EGPD families (\code{"egpd"}, \code{"degpd"},
#'       \code{"zidegpd"}, \code{"ztdegpd"}) a further column \code{xi} repeats
#'       \code{shape}, the GPD tail index, as a convenience alias; the data frame
#'       therefore has one more column than the model has parameters, so select
#'       columns by name (\code{pars$shape}, \code{pars[["kappa"]]}), not by
#'       position or \code{ncol()}.}
#'     \item{\code{"quantile"}}{a data frame with one column per probability,
#'       named \code{"q:<prob>"} with the probability to 15 significant digits
#'       (e.g. \code{"q:0.95"}, \code{"q:0.9995"}).}
#'     \item{\code{"lpmatrix"}}{a list of design matrices, one per linear
#'       predictor, with the offsets in \code{attr(, "offsets")}. They are in the
#'       basis of the fitted coefficients (the matrix times the predictor's
#'       coefficients, plus its offset, is the linear predictor); for \code{t2()}
#'       terms this is \pkg{mgcv}'s fit basis, not the prediction basis of
#'       \code{mgcv::predict.gam()}.}
#'   }
#'   With \code{se.fit = TRUE} (link or response), a list with elements
#'   \code{fitted} and \code{se.fit} (for the response scale, delta-method
#'   standard errors: the link-scale standard error times the absolute value of the
#'   inverse link's derivative, which is negative for a decreasing
#'   \code{custom.fns$unlink}).
#'
#' @export
predict.egpd <- function(object, newdata, type="link", prob=NULL, se.fit=FALSE, marginal=TRUE,
trace = 0, ...) {

## a few checks
family <- object$family

if (family == "egpd") {
  egpd_m <- object$likfns$m
  egpd_iG <- object$likfns$iG
}
if (family == "degpd") {
  degpd_m <- object$likfns$m
  degpd_iG <- object$likfns$iG
}
if (family == "zidegpd") {
  zidegpd_m <- object$likfns$m
  zidegpd_iG <- object$likfns$iG
}
if (family == "custom") {
  q_fn <- object$likfns$q
}
if (family == "comppareto") {
  q_fn <- object$likfns$q
  p_fn <- object$likfns$p
}

if (!is.null(prob))
  type <- "quantile"
if (type == "quantile" & is.null(prob))
  stop("non-NULL `prob' required if `type = quantile'")

## end checks

## standard error setup
if (se.fit) {
  if (marginal) {
    V.type <- "Vc"
  } else {
    V.type <- "Vp"
  }
  conf.pars <- list(object$coefficients, object[[V.type]], object$idpars)
}

got.newdata <- !missing(newdata)

if (got.newdata) {
  ## Parameters pinned via `fixed = list(...)` are carried by synthetic constant
  ## columns the caller has no reason to know about, so re-inject them rather than
  ## demanding the user reproduce them in newdata. No-op when `fixed` was not used.
  if (length(object$fixed.offsets)) {
    for (cn in names(object$fixed.offsets))
      newdata[[cn]] <- rep(object$fixed.offsets[[cn]], nrow(newdata))
  }
  pred.vars <- object$predictor.names
  missing.vars <- pred.vars[!(pred.vars %in% names(newdata))]
  if (length(missing.vars) > 0)
    stop(paste("Variable(s) '", paste(missing.vars, collapse=", "), "' not supplied to `newdata'.", sep=""))
}

if (got.newdata) {
  ndat <- nrow(newdata)
} else {
  if (is.null(object$data))
    stop("Supply `egpd' with `removeData = FALSE' if not supplying `newdata'.")
  ndat <- nrow(object$data)
}

## X creation
X <- .X.egpd(object, newdata)
nX <- length(X)
nms <- names(object)[seq_len(nX)]

## offsets
if (got.newdata) {
  offsets <- attr(X, "offsets")
} else {
  offsets <- object$likdata$offsets
}

if (type == "lpmatrix") {
  attr(X, "offsets") <- offsets
  return(X)
} else {

out <- lapply(seq_len(nX), function(i) {
  eta <- X[[i]] %*% object[[i]]$coefficients
  if (length(offsets[[i]]) > 0) eta <- eta + offsets[[i]]
  eta
})
names(out) <- names(X)
out <- as.data.frame(lapply(out, function(x) x[,1]))

if (se.fit) {
  if (type != "quantile") {
    std.err <- lapply(seq_len(nX), function(i) sqrt(rowSums(X[[i]] * (X[[i]] %*% object[[i]][[V.type]]))))
    std.err <- as.data.frame(std.err)
    names(std.err) <- nms
  }
}

if (type %in% c("response", "quantile")) {
  bd <- .egpd_count_bounds(family)
  for (i in seq_along(nms)) {
    ## delta method: |d par / d eta| times the link-scale SE (a custom.fns$unlink
    ## may be decreasing, e.g. 1/eta). A pin has SE zero without evaluating
    ## a derivative, whose finite-difference probes could leave its domain.
    val <- .egpd_invlink(object, nms[i], out[, i])
    if (se.fit && type == "response" && !.is_pinned_block(object[[i]]$X)) {
      slope <- .egpd_invlink(object, nms[i], out[, i], deriv = TRUE)
      lim <- .egpd_param_bounds(bd, nms[i], link = TRUE)
      if (length(lim)) {
        eta <- out[, i]
        slope[which(eta < lim[1] | eta > lim[2])] <- 0
        slope[which(eta %in% lim & std.err[, i] > 0)] <- NA_real_
      }
      std.err[, i] <- abs(slope) * std.err[, i]
    }
    out[, i] <- val
  }

  nms <- .strip_link(nms)
  names(out) <- nms

  if (se.fit & type == "response")
    names(std.err) <- nms

  out <- .egpd_clamp_response(out, family, bd)

  if (type == "quantile") {

    if (se.fit) {
      stop("Standard errors not yet available for quantile predictions.")
    }

    out <- .egpd_quantile(object, out, prob)
    out <- as.data.frame(out)
    names(out) <- .egpd_qnames(prob)

  } ## end quantile

  ## Convenience alias: the GPD tail index is returned in the column named
  ## "shape" (from the logshape/shape link). Many users reach for `xi`; expose
  ## it as an alias so `pars[["xi"]]` does not silently return NULL (see #5).
  ## Only for the EGPD families, where "shape" is the tail index: for CompPareto
  ## "shape" is the body (gamma/Weibull) shape.
  if (type == "response" && is.data.frame(out) &&
      family %in% c("egpd", "degpd", "zidegpd", "ztdegpd") &&
      "shape" %in% names(out) && !("xi" %in% names(out))) {
    out$xi <- out[["shape"]]
  }

  if (se.fit) {
    out <- list(fitted = out, se.fit = std.err)
  }

} ## end response/quantile

if (type == "link" & se.fit) {
  out <- list(fitted = out, se.fit = std.err)
}

return(out)

} ## end not lpmatrix

}


#' Randomized Quantile Residuals
#'
#' Computes randomized quantile residuals (Dunn & Smyth, 1996) for
#' discrete or continuous EGPD models. For discrete models, a uniform
#' random variate is drawn between the lower and upper CDF bounds at
#' each observation and transformed to the normal scale via
#' \code{\link{qnorm}}. For continuous models, the CDF value is transformed
#' directly. Extreme upper-tail residuals use log survival probabilities so
#' that a CDF rounded to one does not hide a finite residual.
#'
#' For a well-specified model the residuals should be approximately
#' standard normal. The generic \code{rqresid} dispatches to the
#' \code{rqresid.egpd} method for objects of class \code{"egpd"}.
#'
#' The method supports every \code{\link{egpd}} family except
#' \code{"custom"}: the EGPD families (\code{"egpd"}, \code{"degpd"},
#' \code{"zidegpd"}, \code{"ztdegpd"}) with all model types (1--6), the
#' composite-Pareto family (\code{"comppareto"}), and the count families
#' \code{"pig"}, \code{"zipig"}, \code{"gpig"}, \code{"gpignat"},
#' \code{"zigpig"}, \code{"zigpignat"}, \code{"bell"}, \code{"bellnat"},
#' \code{"zibell"}, \code{"zibellnat"}, \code{"cmp"}, \code{"prl"},
#' \code{"gpois"}, \code{"gwaring"}, \code{"ztgwaring"} and
#' \code{"plnorm"}. The fit must keep its data (the default
#' \code{removeData = FALSE}). A non-\code{NULL} \code{seed} is passed to
#' \code{\link{set.seed}}, so it resets the global random number stream.
#'
#' @param object a fitted \code{\link{egpd}} object
#' @param seed optional random seed for reproducibility
#' @param ... unused
#'
#' @return A numeric vector of randomized quantile residuals on the
#'   standard normal scale. Values that are infinite are replaced with
#'   \code{NA}.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996).
#' Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics}, 5(3), 236--244.
#' \doi{10.1080/10618600.1996.10474708}
#'
#' @seealso \code{\link{egpd}}, \code{\link{predict.egpd}}
#'
#' @examples
#' set.seed(42)
#' dat <- data.frame(y = regpd(200, sigma = 1, xi = 0.2, kappa = 2, type = 1))
#' fit <- egpd(y ~ 1, data = dat, family = "egpd", egpd.args = list(m = 1))
#' r <- rqresid(fit, seed = 1)
#' qqnorm(r); qqline(r)
#' @export
rqresid <- function(object, ...) UseMethod("rqresid")

#' @rdname rqresid
#' @export
rqresid.egpd <- function(object, seed = NULL, ...) {

  if (!is.null(seed)) set.seed(seed)
  if (is.null(object$data))
    stop("rqresid requires the original data; refit with removeData = FALSE", call. = FALSE)

  family <- object$family
  y <- object$data[[object$response.name]]
  pars <- predict(object, type = "response")
  pfun <- switch(family, egpd = pegpd, degpd = , ztdegpd = pdiscegpd, zidegpd = pzidiscegpd,
                 comppareto = pcomppareto, pig = ppig, zipig = pzipig, gpig = pGPIG,
                 gpignat = pgpig, zigpig = pZIGPIG, zigpignat = pzigpig,
                 bell = pBELL, bellnat = pbell, zibell = pZIBELL, zibellnat = pzibell,
                 cmp = pcmp, prl = pprl, gpois = pgpois, gwaring = , ztgwaring = pgwaring,
                 plnorm = ppoislnorm, stop("rqresid not implemented for family '", family, "'"))
  args <- switch(family,
    egpd = , degpd = , zidegpd = , ztdegpd = .egpd_dist_args(pars, object$likfns$m),
    comppareto = c(list(spec = object$likfns$spec),
                   .comppareto_named_pars(object$likfns$spec, as.list(pars[seq_len(object$ngam)]))),
    gpig = list(mu = pars$mu, sigma = pars$a, nu = pars$c),
    zigpig = list(mu = pars$mu, sigma = pars$a, nu = pars$c, tau = pars$pi),
    zibell = list(mu = pars$mu, sigma = pars$pi),
    prl = list(tau = .prl_tau_from_mu(pars$mu)),
    as.list(pars))
  if (family == "zidegpd") args$pi <- pars$pi
  base <- function(q, i, lower) {
    a <- lapply(args, function(v) if (length(v) == 1L) v else v[i])
    do.call(pfun, c(list(q = q, lower.tail = lower, log.p = TRUE), a))
  }
  rows <- seq_along(y)
  cdf <- base
  if (family %in% c("ztdegpd", "ztgwaring")) {
    ## S(k | Y > 0) = S(k)/S(0); do not subtract a zero mass rounded to one.
    ls0 <- base(0, rows, FALSE)
    cdf <- function(q, i, lower) {
      ls <- pmin(0, base(q, i, FALSE) - ls0[i])
      ls[q < 1] <- 0
      if (!lower) return(ls)
      lf <- .log1mexp(ls)
      ## When the conditional CDF is tiny, the two log survivors above
      ## can coincide. Difference the lower CDFs instead if they resolve it.
      j <- which(q >= 1 & lf < log(sqrt(.Machine$double.eps)))
      if (length(j)) {
        hi <- base(q[j], i[j], TRUE)
        lo <- base(0, i[j], TRUE)
        k <- which(lo < hi)
        lf[j[k]] <- hi[k] + .log1mexp(lo[k] - hi[k]) - ls0[i[j[k]]]
      }
      lf
    }
  }

  lp <- cdf(y, rows, TRUE)
  upper <- !is.na(lp) & lp > -sqrt(.Machine$double.eps)
  if (any(upper)) lp[upper] <- cdf(y[upper], rows[upper], FALSE)
  if (!family %in% c("egpd", "comppareto")) {
    prev <- rep(NA_real_, length(y))
    if (any(!upper)) prev[!upper] <- cdf(y[!upper] - 1, rows[!upper], TRUE)
    if (any(upper)) prev[upper] <- cdf(y[upper] - 1, rows[upper], FALSE)
    ## U = (1-w) F(y-1) + w F(y); 1-U is the same mixture of survivors.
    w <- runif(length(y))
    lp <- .logspace_add(log1p(-w) + prev, log(w) + lp)
  }
  r <- qnorm(pmin(lp, 0), log.p = TRUE)
  r[upper] <- -r[upper]
  r[is.infinite(r)] <- NA_real_
  return(r)
}
