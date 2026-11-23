## confint and vcov methods for egpd objects

# Backward-compatible inverse-link detector returning a unary function
.detect_invlink <- function(name, object = NULL) {
  function(eta) .egpd_invlink(object, name, eta)
}

# Resolve parm argument to coefficient indices
.resolve_parm <- function(object, parm) {
  coef_names <- names(object$coefficients)
  ngam <- object$ngam
  link_names <- names(object)[seq_len(ngam)]
  stripped <- vapply(link_names, .strip_link, character(1))

  if (is.null(parm)) {
    # Default: one intercept per distributional parameter
    idx <- integer(0)
    for (i in seq_len(ngam)) {
      sub_coefs <- names(object[[i]]$coefficients)
      intercept_pos <- match("(Intercept)", sub_coefs)
      if (is.na(intercept_pos)) intercept_pos <- match("(fixed)", sub_coefs)
      if (!is.na(intercept_pos)) {
        # Find the global index
        offset <- sum(vapply(seq_len(i - 1), function(j) length(object[[j]]$coefficients), integer(1)))
        idx <- c(idx, offset + intercept_pos)
      }
    }
    return(idx)
  }

  if (is.numeric(parm)) {
    if (any(!is.finite(parm) | parm != trunc(parm) | parm < 1 | parm > length(coef_names)))
      stop("Numeric 'parm' out of range: supply integer coefficient indices")
    return(as.integer(parm))
  }

  if (is.character(parm)) {
    idx <- integer(0)
    for (p in parm) {
      # Try matching against stripped link names
      m <- match(p, stripped)
      if (is.na(m)) {
        aliases <- c(xi = "shape", sigma = "scale", psi = "scale")
        if (p %in% names(aliases)) m <- match(aliases[[p]], stripped)
      }
      if (!is.na(m)) {
        # Return the intercept (or pinned offset-only term) of this parameter
        offset <- sum(vapply(seq_len(m - 1), function(j) length(object[[j]]$coefficients), integer(1)))
        sub_coefs <- names(object[[m]]$coefficients)
        intercept_pos <- match("(Intercept)", sub_coefs)
        if (is.na(intercept_pos)) intercept_pos <- match("(fixed)", sub_coefs)
        if (is.na(intercept_pos))
          stop(sprintf(paste0("Parameter '%s' has no intercept; request its coefficients ",
                              "by name (%s) or supply `newdata`"),
                       p, paste(coef_names[offset + seq_along(sub_coefs)], collapse = ", ")),
               call. = FALSE)
        idx <- c(idx, offset + intercept_pos)
      } else {
        # Try matching against full coefficient names
        m2 <- match(p, coef_names)
        if (!is.na(m2)) {
          idx <- c(idx, m2)
        } else {
          stop(sprintf("Parameter '%s' not found", p))
        }
      }
    }
    return(idx)
  }

  stop("'parm' must be NULL, numeric, or character")
}

# Find a root of root_fn in one direction via bracketing + uniroot
.find_profile_root <- function(root_fn, mle, se, direction, cutoff) {
  # direction: -1 for lower, +1 for upper
  step <- direction * se
  lo <- mle
  hi <- mle + 3 * step

  # Evaluate at MLE (should be -cutoff since deviance = 0 there)
  f_lo <- suppressWarnings(root_fn(mle))
  if (!is.finite(f_lo)) return(NA_real_)

  # Evaluate at initial bracket endpoint, expand if needed (up to 5 doublings,
  # i.e. out to 96 steps: standard errors, or .profile_one_param()'s fallback step)
  f_hi <- suppressWarnings(root_fn(hi))
  f_prev <- NA_real_
  for (j in seq_len(5)) {
    if (is.finite(f_hi) && f_hi > 0) break
    hi <- mle + (hi - mle) * 2
    f_prev <- f_hi
    f_hi <- suppressWarnings(root_fn(hi))
  }

  if (!is.finite(f_hi) || f_hi <= 0) {
    # Still below the cutoff 96 steps out. If the profile deviance has levelled off
    # (no change between 48 and 96 steps), it stays below the cutoff as the coefficient
    # tends to infinity: the interval reaches the boundary of the parameter space,
    # i.e. the link-scale limit, which the inverse link maps to e.g. 0 for log(xi).
    if (is.finite(f_hi) && is.finite(f_prev) && abs(f_hi - f_prev) < 1e-3)
      return(direction * Inf)
    warning("Could not bracket profile likelihood root", call. = FALSE)
    return(NA_real_)
  }

  # Order bounds for uniroot
  bounds <- sort(c(lo, hi))

  res <- tryCatch(
    suppressWarnings(uniroot(root_fn, interval = bounds,
                             tol = .Machine$double.eps^0.25)),
    error = function(e) NULL
  )

  if (!is.null(res)) return(res$root)

  # Fallback: manual bisection (30 iterations gives ~1e-9 precision)
  a <- bounds[1]
  b <- bounds[2]
  fa <- suppressWarnings(root_fn(a))
  fb <- suppressWarnings(root_fn(b))
  # Ensure fa < 0 and fb > 0
  if (fa > 0 && fb < 0) {
    tmp <- a; a <- b; b <- tmp
    tmp <- fa; fa <- fb; fb <- tmp
  }
  for (j in seq_len(30)) {
    mid <- (a + b) / 2
    fmid <- suppressWarnings(root_fn(mid))
    if (!is.finite(fmid)) {
      # Shrink toward the safe side (MLE)
      if (abs(mid - mle) > abs(a - mle)) b <- mid else a <- mid
      next
    }
    if (fmid < 0) a <- mid else b <- mid
  }
  (a + b) / 2
}

## Box for intercept-only nuisance coefficients. Beyond these limits *every* row
## is on the flat side of its clamp; clipping there loses no feasible likelihood,
## even with varying offsets. A fixed parameter is not a free intercept, nor are
## bounds on transformed predictors coefficient bounds when CH mixes coordinates.
.profile_nuisance_box <- function(object, k) {
  bd <- .egpd_count_bounds(object$family)
  if (!length(bd)) return(NULL)
  ld <- object$likdata
  ix <- seq_len(object$ngam)[-k]
  lo <- rep(-Inf, length(ix)); hi <- -lo
  origin <- numeric(length(ix))
  for (j in seq_along(ix)) {
    i <- ix[j]
    nm <- names(object)[i]
    eta <- .egpd_param_bounds(bd, nm, link = TRUE)
    if (!length(eta) || .is_pinned_block(ld$X[[i]])) next
    off <- ld$offsets[[i]]
    if (!length(off)) off <- 0
    a <- 1
    if (isTRUE(ld$adjust == 2L)) {
      ## A bound on i is only a box constraint if moving i moves no other
      ## predictor, and no other coefficient moves predictor i.
      a <- ld$CH[i, i]
      if (a == 0 || any(ld$CH[i, -i] != 0) || any(ld$CH[-i, i] != 0)) next
      off <- off + (1 - a) * ld$compmode[i]
    }
    bounds <- sort((eta - c(max(off), min(off))) / a)
    lo[j] <- bounds[1L]; hi[j] <- bounds[2L]
    origin[j] <- (max(eta[1L] + 1, min(0, eta[2L] - 1)) - mean(off)) / a
  }
  if (any(is.finite(c(lo, hi)))) list(lo = lo, hi = hi, origin = origin)
}

# Profile one coefficient. `monitor`, an environment with `dev` (0 to start with), records
# the lowest deviance of any sound inner fit: one below zero is a point with a higher
# likelihood than the fit's (see .warn_profile_better_fit()).
.profile_one_param <- function(object, k, level, monitor = NULL) {
  beta_hat <- object$coefficients
  likdata <- object$likdata
  likfns <- object$likfns
  cutoff <- qchisq(level, 1)
  ## The fitting objective is the negative log-likelihood times likdata$k, which is
  ## 1/gamma unless a sandwich adjustment is in force (gamma must then be 1, and k is
  ## the adjustment, whose adjusted likelihood ratio is the one to profile). As in
  ## mgcv, `gamma` only steers smoothing-parameter selection, so the profile uses the
  ## unscaled log-likelihood against the chi-squared cutoff: the scaled one would
  ## solve a likelihood ratio of gamma times the cutoff.
  if (!isTRUE(likdata$adjust > 0)) likdata$k <- 1

  nllh_mle <- suppressWarnings(.nllh.nopen(beta_hat, likdata, likfns))
  ## The bracketing step is the coefficient's standard error, unless that gives no usable
  ## scale: where the Hessian behind Vp is singular (`hessian.rank`, see egpd()) its
  ## variances do not describe the likelihood -- at the generalized Waring's a = k, Vp gave
  ## log k a standard error of 1.6e-8, 96 of which never reached the cutoff, so the interval
  ## for k came out as [0, Inf] instead of [1.06, 8.52] -- and a variance that is zero or
  ## negligible (a standard error below sqrt(eps) (|estimate| + 1)) gives none either.
  vk <- diag(object$Vp)[k]
  hr <- object$hessian.rank
  singular <- length(hr) == 2L && !anyNA(hr) && hr[[1L]] < hr[[2L]]
  se_k <- if (!singular && vk > .Machine$double.eps * (abs(beta_hat[k]) + 1)^2) sqrt(vk) else
    abs(beta_hat[k]) * 0.1 + 0.1
  nuisance_start <- beta_hat[-k]
  box <- .profile_nuisance_box(object, k)

  ## A fit with a coefficient beyond +-700 on its link scale, where exp() leaves the normal
  ## floating-point range, gives no value of the likelihood (see `trusted` below).
  sound <- function(par) all(abs(par) <= 700)

  ## L-BFGS-B differences the objective *inside* the feasible box at an edge:
  ## the analytic score of the clamp itself is zero, not the inward derivative.
  bounded_fit <- function(start, v) {
    opt <- function(p) suppressWarnings(tryCatch(
      optim(p, .nllh.nopen, likdata = likdata, likfns = likfns, id = k, vals = v,
            method = "L-BFGS-B", lower = box$lo, upper = box$hi,
            ## Scale the search, not the LR. A step along the *summed* score
            ## can go straight to an invalid corner and stop at its start.
            control = list(fnscale = likdata$nobs, ndeps = rep(1e-5, length(p)),
                           factr = 1e4, maxit = 500)),
      error = function(e) NULL))
    p <- pmax(box$lo, pmin(box$hi, start))
    best <- opt(p)
    edge <- p == box$lo | p == box$hi
    if (any(edge)) {
      ## A limiting law (e.g. Waring at rho = 1e15) may already be flat well
      ## *inside* the box. Compare the limiting and interior branches.
      other <- opt(replace(p, edge, box$origin[edge]))
      if (!is.null(other) && is.finite(other$value) &&
          (is.null(best) || !is.finite(best$value) || other$value < best$value)) best <- other
    }
    best
  }

  inner_fit <- function(start, v) {
    opt <- if (!is.null(box)) bounded_fit(start, v) else suppressWarnings(tryCatch(
      optim(par = start, fn = .nllh.nopen, gr = .grad.nopen,
            likdata = likdata, likfns = likfns,
            id = k, vals = v, method = "BFGS",
            control = list(maxit = 200)),
      error = function(e) NULL
    ))

    if (is.null(box) && (is.null(opt) || !is.finite(opt$value))) {
      # Fallback to Nelder-Mead
      opt <- suppressWarnings(tryCatch(
        optim(par = start, fn = .nllh.nopen,
              likdata = likdata, likfns = likfns,
              id = k, vals = v, method = "Nelder-Mead",
              control = list(maxit = 500)),
        error = function(e) NULL
      ))
    }

    if (is.null(opt) || !is.finite(opt$value)) return(NULL)
    if (!is.null(monitor) && sound(opt$par) && 2 * (opt$value - nllh_mle) < monitor$dev) {
      monitor$dev <- 2 * (opt$value - nllh_mle)
      monitor$coef <- add_in(opt$par, k, v)
      monitor$k <- k
    }
    opt
  }

  ## Each fit starts from the estimates of the other coefficients. That start can lie in the
  ## wrong basin once the profile path of a weakly identified coefficient has moved far from
  ## it: for a continuous EGPD model-3 fit with log delta = -6.76 (standard error 861), the
  ## optimal log delta is about 1.1 at scale 1.93, and the fit from -6.76 ran instead onto
  ## the plateau delta -> Inf, where the likelihood is flat (log delta 170, deviance 69.6
  ## instead of 0.52). A failed fit can only overstate the deviance, so the only harm is a
  ## spurious crossing of the cutoff, which gave the upper scale limit 1.92 instead of 2.385.
  ## So every evaluated point is kept with its solution, and a limit is accepted once the
  ## point just beyond it, refitted from other starts (confirmed()), still exceeds the
  ## cutoff. If it does not, the search continues from that point with every fit also
  ## started from the solution at the nearest point below the cutoff (the profile path
  ## followed by continuation). Values are replaced only where the second fit is lower by
  ## more than `tol_dev`, so where the first fit attains the minimum nothing changes.
  ##
  ## Only trusted fits serve as starts or replacements: a deviance below -tol_dev (a fit
  ## better than the estimate) or a coefficient beyond +-700 on its link scale (where exp()
  ## leaves the normal floating-point range) is no value of the profile. Fits from the
  ## estimate reached both on a continuous EGPD model-5 fit, where log kappa ran to -742
  ## and the log-likelihood, evaluated with a denormal kappa, came out 3 above the maximum;
  ## followed, such points gave a lower shape limit of -0.319 instead of -0.130.
  tol_dev <- 0.01
  trusted <- function(dev, par) is.finite(dev) && dev >= -tol_dev && sound(par)
  pts <- new.env()
  pts$v <- unname(beta_hat[k])
  pts$dev <- 0
  pts$ok <- TRUE
  pts$par <- list(unname(nuisance_start))
  keep <- function(v, dev, par) {
    pts$v <- c(pts$v, v)
    pts$dev <- c(pts$dev, dev)
    pts$ok <- c(pts$ok, trusted(dev, par))
    pts$par[[length(pts$par) + 1L]] <- par
  }
  ## the trusted kept point below the cutoff nearest to v, on v's side of the estimate
  anchor <- function(v) {
    ok <- which(pts$ok & pts$dev < cutoff & (pts$v - beta_hat[k]) * (v - beta_hat[k]) >= 0)
    ok[which.min(abs(pts$v[ok] - v))]
  }

  profile_deviance <- function(v, follow = FALSE) {
    opt <- inner_fit(nuisance_start, v)
    dev <- if (is.null(opt)) Inf else 2 * (opt$value - nllh_mle)
    par <- if (is.null(opt)) NULL else opt$par
    if (follow) {
      if (!is.null(par) && !trusted(dev, par)) {
        dev <- Inf
        par <- NULL
      }
      a <- anchor(v)
      ow <- if (length(a)) inner_fit(pts$par[[a]], v)
      if (!is.null(ow) && trusted(2 * (ow$value - nllh_mle), ow$par) &&
          2 * (ow$value - nllh_mle) < dev - tol_dev) {
        dev <- 2 * (ow$value - nllh_mle)
        par <- ow$par
      }
    }
    if (!is.null(par)) keep(v, dev, par)
    dev
  }

  ## Coordinates with a standard error above 10 on their link scale sit on a plateau at the
  ## estimate, where the derivative-based fits cannot move them: for a DEGPD model-4 fit
  ## with log xi = -13.5 (standard error 711), every fit of the kappa profile kept xi near
  ## 1e-6, past a flat stretch from the interior optimum at xi = 0.12. A restart sets them
  ## to 0.
  weak <- which(diag(object$Vp)[-k] > 100)
  reset_start <- if (length(weak)) replace(nuisance_start, weak, 0)

  ## Model 4 has two limiting carriers even when delta is well identified at the
  ## MLE. Both can compete with the interior optimum elsewhere on the profile.
  ## In terms of the GPD cumulative hazard A, their base survivors are
  ## exp(-A)*(1+A) (delta -> 0) and exp(-A) (delta -> Inf). Seed both sides;
  ## optimise and compare the actual finite-parameter likelihood, not a limit.
  delta_starts <- numeric()
  # ZIDEGPD puts kappa before delta; the other model-4 families reverse them.
  d <- match("logdelta", names(object)[seq_len(object$ngam)])
  if (object$family %in% c("egpd", "degpd", "zidegpd") && isTRUE(likfns$m == 4L) &&
      !is.na(d) && identical(names(object[[d]]$coefficients), "(Intercept)")) {
    delta_idx <- which(likdata$idpars[-k] == d)
    if (length(delta_idx) == 1L) {
      off <- likdata$offsets[[d]]
      delta_starts <- c(1, -1) * log(sqrt(.Machine$double.eps)) -
        (if (length(off)) stats::median(off) else 0)
    }
  }

  ## Is the limit r real? The kept point b beyond it that exceeds the cutoff must still
  ## exceed it when refitted from the solution at the trusted kept point below the cutoff
  ## nearest to r, from `reset_start`, from the origin of the link scales, and at the end
  ## of a short continuation from the estimate to b (K fits, each from the previous
  ## solution, K = 4 for a limit within 4 steps of the estimate). The continuation is what
  ## finds the interior minimum where the fits from the estimate, at b and at its
  ## neighbours alike, stall on a flat direction or end in a degenerate basin with a
  ## nuisance coefficient run off to -26 or -197 on its link scale (DEGPD model 5 and
  ## continuous EGPD model 5 fits: deviances 3.548 and 3.786 where they gave 3.841). The
  ## origin finds a branch the others miss where the profile path switches basins (a DEGPD
  ## model-4 scale profile below 2.45, where log delta drops from 2.4 to -18: deviance 1.31
  ## where they gave 3.84). If b does not exceed the cutoff after all, its corrected value
  ## is kept and b is returned, to continue from.
  confirmed <- function(r, direction) {
    hi <- which((pts$v - r) * direction >= 0 & pts$dev >= cutoff)
    lo <- which(pts$ok & (pts$v - r) * direction <= 0 & pts$dev < cutoff &
                  (pts$v - beta_hat[k]) * direction >= 0)
    if (!length(hi) || !length(lo)) return(TRUE)
    b <- hi[which.min(abs(pts$v[hi] - r))]
    a <- lo[which.min(abs(pts$v[lo] - r))]
    vb <- pts$v[b]
    best <- NULL
    take <- function(opt) {
      if (!is.null(opt) && trusted(2 * (opt$value - nllh_mle), opt$par) &&
          (is.null(best) || opt$value < best$value)) best <<- opt
    }
    take(inner_fit(pts$par[[a]], vb))
    if (!is.null(reset_start)) take(inner_fit(reset_start, vb))
    take(inner_fit(numeric(length(nuisance_start)), vb))
    for (ld in delta_starts)
      take(inner_fit(replace(pts$par[[b]], delta_idx, ld), vb))
    nk <- min(16L, max(4L, ceiling(abs(vb - beta_hat[k]) / se_k)))
    st <- nuisance_start
    for (j in seq_len(nk)) {
      v <- beta_hat[k] + (vb - beta_hat[k]) * j / nk
      opt <- inner_fit(st, v)
      if (is.null(opt) || !trusted(2 * (opt$value - nllh_mle), opt$par)) break
      st <- opt$par
      if (j < nk) keep(v, 2 * (opt$value - nllh_mle), opt$par) else take(opt)
    }
    if (is.null(best)) return(TRUE)
    dev <- 2 * (best$value - nllh_mle)
    if (!(dev < min(cutoff, pts$dev[b] - tol_dev))) return(TRUE)
    pts$dev[b] <- dev
    pts$ok[b] <- TRUE
    pts$par[[b]] <- best$par
    vb
  }

  limit <- function(direction) {
    r <- .find_profile_root(function(v) profile_deviance(v) - cutoff, beta_hat[k], se_k,
                            direction, cutoff)
    for (i in seq_len(10L)) {
      if (!is.finite(r)) return(r)
      from <- confirmed(r, direction)
      if (isTRUE(from)) return(r)
      r <- .find_profile_root(function(v) profile_deviance(v, follow = TRUE) - cutoff, from,
                              se_k, direction, cutoff)
    }
    r
  }

  lower <- limit(-1)
  upper <- limit(+1)

  c(lower, upper)
}

#' Variance-covariance matrix for egpd fits
#'
#' @param object A fitted \code{egpd} object.
#' @param marginal Logical; if \code{TRUE} (the default), return \code{object$Vc},
#'   which also accounts for smoothing-parameter uncertainty (it equals
#'   \code{object$Vp} for models without smooth terms); if \code{FALSE}, return
#'   the conditional Bayesian covariance \code{object$Vp}. The default matches
#'   \code{\link{predict.egpd}}, \code{\link{simulate.egpd}},
#'   \code{\link{confint.egpd}} and \code{\link{predictive_coverage}}.
#' @param ... Not used.
#' @return The variance-covariance matrix with row and column
#'   names matching the coefficient names.
#' @export
vcov.egpd <- function(object, marginal = TRUE, ...) {
  V <- if (isTRUE(marginal) && !is.null(object$Vc)) object$Vc else object$Vp
  nms <- names(object$coefficients)
  rownames(V) <- colnames(V) <- nms
  V
}

## A parameter's response-scale value(s) from its linear predictor, by the fitted
## inverse link. Wald intervals are for the parameter actually predicted.
## A profile plateau can still extend to the statistical boundary, not the
## numerical guard. Respect the Waring rho > 1 domain in both cases.
.confint_response <- function(object, name, eta, effective = TRUE) {
  val <- .egpd_invlink(object, name, eta)
  fam <- object$family
  pn <- .strip_link(name)
  lim <- .egpd_param_bounds(.egpd_count_bounds(fam), pn)
  if (length(lim)) {
    if (effective) val <- pmin(pmax(val, lim[1]), lim[2])
    else if (fam %in% c("gwaring", "ztgwaring") && pn == "rho") val <- pmax(val, lim[1])
  }
  val
}

.check_confint_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level < 0 || level > 1)
    stop("'level' must be one finite probability in [0, 1]", call. = FALSE)
}

## Response-scale limits from `lo` and `hi`, the images of link-scale limits
## lower <= upper under a parameter's (monotone) inverse link. A decreasing one, e.g.
## a custom.fns$unlink exp(-eta) or 1/eta, maps the lower limit above the upper, so
## such pairs are swapped. Where one image is NA (a profile root that could not be
## bracketed), the other is compared with `mid(i)`, the images of the estimates at the
## positions `i` (logical), which lie between the two limits on the link scale.
.confint_orient <- function(lo, hi, mid, check_full = FALSE) {
  dec <- !is.na(lo) & !is.na(hi) & lo > hi
  one <- xor(is.na(lo), is.na(hi))
  if (any(one)) {
    m <- mid(one)
    dec[one] <- (!is.na(lo[one]) & lo[one] > m) | (!is.na(hi[one]) & hi[one] < m)
    dec[is.na(dec)] <- FALSE
  }
  if (any(dec)) {
    tmp <- lo[dec]
    lo[dec] <- hi[dec]
    hi[dec] <- tmp
  }
  if (any(check_full)) {
    i <- rep_len(check_full, length(lo))
    center <- mid(i)
    if (!all(is.finite(center) & !is.na(lo[i]) & !is.na(hi[i]) &
             lo[i] <= center & center <= hi[i]))
      stop("Cannot determine a 100% response-scale interval: the custom inverse link ",
           "at -Inf and Inf does not bound its fitted value.", call. = FALSE)
  }
  list(lower = lo, upper = hi)
}

## Per selected coefficient: its predictor, whether it is a parameter-level quantity
## (the only coefficient of a predictor whose offset is constant, i.e. an
## intercept-only or pinned parameter), that constant offset, and the row label.
.confint_coef_info <- function(object, idx) {
  coef_names <- names(object$coefficients)
  link_names <- names(object)[seq_len(object$ngam)]
  cumlen <- cumsum(vapply(seq_len(object$ngam),
                          function(i) length(object[[i]]$coefficients), integer(1)))
  lapply(idx, function(k) {
    sub_idx <- which(k <= cumlen)[1]
    local_idx <- k - if (sub_idx > 1L) cumlen[sub_idx - 1L] else 0L
    sub_coefs <- names(object[[sub_idx]]$coefficients)
    off <- object$likdata$offsets[[sub_idx]]
    shift <- if (!length(off)) 0 else if (length(unique(off)) == 1L) off[1] else NA_real_
    natural <- length(sub_coefs) == 1L && sub_coefs %in% c("(Intercept)", "(fixed)") &&
      !is.na(shift)
    nm <- link_names[sub_idx]
    list(natural = natural,
         pinned = natural && identical(sub_coefs[local_idx], "(fixed)"),
         shift = shift,
         label = if (natural) .strip_link(nm) else coef_names[k],
         invlink = function(eta, effective) .confint_response(object, nm, eta, effective))
  })
}

#' Confidence intervals for egpd model parameters
#'
#' @param object A fitted \code{egpd} object.
#' @param parm Which parameters to compute CIs for.  \code{NULL} (default)
#'   selects the intercept or fixed offset-only term of each distributional
#'   parameter. Can also be a
#'   character vector of stripped parameter names (e.g. \code{"scale"},
#'   \code{"shape"}; aliases \code{"sigma"}/\code{"psi"} and \code{"xi"}
#'   are accepted) or of full coefficient names (e.g. \code{"logscale_x"}), or an
#'   integer vector of coefficient indices. Requesting a parameter by name whose
#'   predictor has no intercept is an error. With \code{newdata}, a character
#'   vector of parameter names (default: all parameters).
#' @param level Confidence level in \eqn{[0,1]} (default 0.95). At level 1
#'   there is no finite likelihood-ratio cutoff: profile limits are the
#'   link-scale boundaries, transformed for response-scale rows. Parameters
#'   held fixed retain their point intervals. If a custom inverse link has
#'   a restricted domain and its two limits do not bound the fitted value,
#'   a 100\% response-scale interval is not available.
#' @param method Either \code{"wald"} (default) or \code{"profile"}.  Profile
#'   likelihood CIs are only available for intercept-only models (no smooth
#'   terms or covariates beyond an intercept), and not with \code{newdata}. A
#'   profile interval can extend to the boundary of the parameter space (see
#'   Value); an endpoint that cannot be bracketed otherwise is \code{NA}, with a
#'   warning. The profile is that of the log-likelihood itself, compared with the
#'   chi-squared cutoff: a fit made with \code{gamma != 1} (see \code{\link{egpd}})
#'   is profiled without dividing the log-likelihood by \code{gamma}, since, as in
#'   \pkg{mgcv}, \code{gamma} only affects how smoothing parameters are selected. A
#'   sandwich-adjusted fit is profiled on its adjusted log-likelihood. GPIG
#'   and ZIGPIG profiles use the fit's saved likelihood evaluator, not the
#'   current global setting; older fits without these settings must be refitted.
#' @param marginal Logical; for \code{method = "wald"}, use \code{object$Vc}
#'   (accounting for smoothing-parameter uncertainty) when \code{TRUE}, the
#'   default (as in \code{\link{predict.egpd}}); \code{FALSE} uses \code{object$Vp}.
#' @param newdata Optional data frame. If supplied, Wald intervals for the
#'   distributional parameters themselves are returned at each of its rows: the
#'   interval for the linear predictor (including offsets) is transformed by the
#'   fitted inverse link, including count-family clamps. The same clamps
#'   apply to parameter-level Wald intervals without \code{newdata}, but
#'   not to coefficient-scale intervals.
#' @param ... Not used.
#' @return A matrix with columns \code{"lower"} and \code{"upper"}.
#'   \itemize{
#'     \item Without \code{newdata}, a coefficient that alone determines its
#'       parameter -- the intercept of an intercept-only predictor (plus any
#'       constant offset), or a pinned offset-only term -- gives a
#'       \emph{parameter-level} interval on the response (natural) scale, using the
#'       fitted inverse link (including bounded shape links, and
#'       \code{custom.fns$unlink} for \code{family = "custom"}); its row is named by
#'       the parameter (e.g. \code{"scale"}). A decreasing inverse link (a
#'       \code{custom.fns$unlink} such as \code{1/eta} or \code{exp(-eta)}) maps the
#'       upper link-scale limit to the lower response-scale one, so
#'       \code{"lower"} never exceeds \code{"upper"} (here and with
#'       \code{newdata}). Pinned (constant offset-only)
#'       parameters have degenerate intervals at their fixed values. The
#'       generalized Waring \eqn{\rho}, log-linked although its parameter space is
#'       \eqn{\rho > 1}, is evaluated by the likelihood (and reported by
#'       \code{\link{predict.egpd}}) as at least \eqn{1 + 10^{-6}}, so a lower
#'       limit that would fall below that is \eqn{1 + 10^{-6}}.
#'     \item Every other coefficient (slopes, smooth coefficients, and intercepts of
#'       predictors with covariates or varying offsets) gives an interval for the
#'       coefficient itself on the link scale; its row is named by the full
#'       coefficient name (e.g. \code{"logscale_x"}). Transforming such an
#'       interval by the inverse link would not be an interval for any parameter.
#'     \item With \code{newdata}, rows are named \code{"<parameter>[<row>]"} and
#'       all intervals are on the response scale (\eqn{\rho} as above). A row whose
#'       linear predictor involves a coefficient without positive variance (an
#'       unidentified one) has \code{NA} limits, as without \code{newdata}; a
#'       pinned parameter has a degenerate interval at its fixed value.
#'     \item For \code{method = "profile"}, the search for each limit starts from the
#'       estimate in steps of the coefficient's standard error, or, where that is not a
#'       usable scale -- the Hessian behind the standard errors is singular
#'       (\code{hessian.rank}, see \code{\link{egpd}}), or the variance is zero or
#'       negligible -- of \eqn{0.1 (|\hat\beta| + 1)}{0.1 (|estimate| + 1)}. Each
#'       point of the profile is fitted from the estimates of the other coefficients.
#'       Where hard count-family bounds can be imposed on the nuisance intercepts,
#'       they are optimised within those bounds, using inward differences at an edge
#'       and comparing an interior start for a nuisance estimate on a clamp.
#'       On weakly identified fits that start can lie on a plateau or in another basin
#'       far from the estimate, where the fit overstates the deviance and would cut the
#'       interval short, so a finite limit is accepted only once the point just beyond
#'       it still exceeds the cutoff when refitted from the solution at its neighbour
#'       inside the interval, from the estimates with every coefficient whose standard
#'       error exceeds 10 set to 0, from all coefficients 0, and at the end of a short
#'       continuation from the estimate. For a model-4 carrier with a free nuisance
#'       \eqn{\delta}, starts near both limiting carriers are also checked: the better
#'       branch can change between the estimate and the limit. If a checked point is
#'       below the cutoff, the search continues outward, each
#'       point also fitted from the solution at the nearest point below the cutoff,
#'       whose value replaces the first where it is lower by more than 0.01; so where
#'       the fits from the estimates attain the minimum the interval is unchanged. Fits
#'       with a deviance below that of the estimate, or with a coefficient beyond 700 in
#'       absolute value on its link scale (where floating-point evaluation of the
#'       likelihood fails), are not used.
#'       If the profile deviance is still below the chi-squared cutoff 96 steps from the
#'       estimate and has levelled off (it changes by less than 0.001 between 48 and 96
#'       steps), the interval extends to the boundary of the parameter space in that
#'       direction.
#'       The endpoint is then the limit of the inverse link as the link-scale
#'       coefficient tends to \code{-Inf} or \code{Inf}: for example 0 (lower) or
#'       \code{Inf} (upper) for a log-link parameter such as the shape \eqn{\xi},
#'       0 or 1 for logit, probit and cloglog links, and 0 or \code{xi.max} for the
#'       bounded shape link (\eqn{1 + 10^{-6}} for a lower generalized Waring
#'       \eqn{\rho}, as above). Link-scale rows report \code{-Inf} or \code{Inf}.
#'     \item For \code{method = "profile"}, each profile point refits the other
#'       coefficients. If one of these fits reaches a log-likelihood more than 0.05
#'       above the fitted one (a profile deviance below \eqn{-0.1}), the fit is not the
#'       maximum likelihood estimate (as \pkg{MASS}'s \code{profile.glm()} reports, the
#'       profile "has found a better solution"): a warning of class
#'       \code{"egpd_profile_better_fit"}, given once per call, reports the largest
#'       gain found and the coefficient whose profile found it. The intervals are still
#'       those of the profile relative to the fit, which is not optimal; refit the
#'       model, e.g. from other starting values (\code{inits}) or, for families
#'       \code{"egpd"}, \code{"degpd"}, \code{"zidegpd"} and \code{"comppareto"}, with
#'       \code{restarts = TRUE}, and profile the new fit.
#'   }
#' @export
confint.egpd <- function(object, parm = NULL, level = 0.95, method = c("wald", "profile"),
                         marginal = TRUE, newdata = NULL, ...) {
  method <- match.arg(method)
  .check_confint_level(level)
  if (!is.null(newdata))
    return(.confint_newdata(object, parm, level, method, marginal, newdata))
  idx <- .resolve_parm(object, parm)
  info <- .confint_coef_info(object, idx)
  labels <- vapply(info, `[[`, character(1), "label")

  ## response-scale bounds for parameter-level rows, link-scale bounds otherwise
  finish <- function(j, bounds) {
    if (!info[[j]]$natural) return(bounds)
    sh <- info[[j]]$shift
    effective <- method == "wald" || info[[j]]$pinned
    b <- info[[j]]$invlink(bounds + sh, effective)
    if (length(b) != 2L) return(b)  # a pinned parameter's single value
    b <- .confint_orient(b[1L], b[2L],
                         function(i) info[[j]]$invlink(object$coefficients[[idx[j]]] + sh, effective),
                         check_full = level == 1 && object$family == "custom")
    c(b$lower, b$upper)
  }

  ci <- matrix(NA_real_, nrow = length(idx), ncol = 2)
  colnames(ci) <- c("lower", "upper")
  rownames(ci) <- labels

  if (method == "wald") {
    z <- qnorm((1 - level) / 2, lower.tail = FALSE)
    coefs <- object$coefficients[idx]
    V_mat <- if (isTRUE(marginal) && !is.null(object$Vc)) object$Vc else object$Vp
    vars <- diag(V_mat)[idx]

    for (j in seq_along(idx)) {
      if (info[[j]]$pinned) {
        ci[j, ] <- finish(j, 0)
        next
      }
      if (vars[j] <= 0) next  # non-positive variance: CI undefined
      ci[j, ] <- finish(j, coefs[j] + c(-1, 1) * z * sqrt(vars[j]))
    }

    return(ci)
  }

  # Profile method
  if (length(object$gotsmooth) > 0)
    stop("Profile likelihood CIs are only available for intercept-only models (no smooth terms)")

  # Check that each sub-model is intercept-only
  for (i in seq_len(object$ngam)) {
    cn <- names(object[[i]]$coefficients)
    if (length(cn) != 1L || !cn %in% c("(Intercept)", "(fixed)"))
      stop("Profile likelihood CIs are only available for intercept-only models (no covariates beyond intercept)")
  }

  if (object$family %in% .gpig_families && level > 0 && level < 1 &&
      any(!vapply(info, `[[`, logical(1), "pinned"))) {
    if (is.null(object$gpig.control))
      stop("GPIG likelihood settings were not recorded; refit the model before profiling.",
           call. = FALSE)
    old <- do.call(gpig_control, object$gpig.control)
    on.exit(do.call(gpig_control, old), add = TRUE)
  }
  monitor <- new.env()
  monitor$dev <- 0
  for (j in seq_along(idx)) {
    if (info[[j]]$pinned) {
      ci[j, ] <- finish(j, 0)
      next
    }
    k <- idx[j]
    bounds <- if (level == 0) {
      rep(object$coefficients[k], 2L)
    } else if (level == 1) {
      c(-Inf, Inf)
    } else .profile_one_param(object, k, level, monitor)
    ci[j, ] <- finish(j, bounds)
  }
  if (monitor$dev < -0.1) {
    nm <- names(object$coefficients)
    .warn_profile_better_fit(monitor$dev,
                             if (is.null(nm)) paste("coefficient", monitor$k) else nm[monitor$k])
  }

  ci
}

## Profiling refits the other coefficients at each point. A fit that reaches a deviance
## below -0.1 (a log-likelihood more than 0.05 above the fitted one, well beyond optimiser
## tolerance) shows that the fit is not the maximum, as MASS's profile.glm() reports
## ("profiling has found a better solution"). The intervals are still those of the profile
## relative to the fit.
.warn_profile_better_fit <- function(dev, coef_name) {
  msg <- sprintf(paste0(
    "confint(): profiling found a log-likelihood %s higher than the fitted one (in the ",
    "profile of %s), so the fit is not the maximum likelihood estimate, and the profile ",
    "intervals are relative to a fit that is not optimal. Refit the model, e.g. from other ",
    "starting values (`inits`) or, for families \"egpd\", \"degpd\", \"zidegpd\" and \"comppareto\", with ",
    "`restarts = TRUE`."), format(signif(-dev / 2, 3)), coef_name)
  warning(structure(class = c("egpd_profile_better_fit", "warning", "condition"),
                    list(message = msg, call = NULL)))
  invisible(TRUE)
}

## Wald intervals for the distributional parameters at the rows of `newdata`.
.confint_newdata <- function(object, parm, level, method, marginal, newdata) {
  if (method != "wald")
    stop("intervals at `newdata` are available for method = \"wald\" only", call. = FALSE)
  link_names <- names(object)[seq_len(object$ngam)]
  stripped <- vapply(link_names, .strip_link, character(1))
  sel <- seq_len(object$ngam)
  if (!is.null(parm)) {
    if (!is.character(parm))
      stop("with `newdata`, `parm` must name distributional parameters", call. = FALSE)
    aliases <- c(xi = "shape", sigma = "scale", psi = "scale")
    sel <- vapply(parm, function(p) {
      m <- match(p, stripped)
      if (is.na(m) && p %in% names(aliases)) m <- match(aliases[[p]], stripped)
      if (is.na(m)) stop(sprintf("Parameter '%s' not found", p), call. = FALSE)
      m
    }, integer(1))
  }
  X <- predict.egpd(object, newdata, type = "lpmatrix")
  offs <- attr(X, "offsets")
  V.type <- if (isTRUE(marginal) && !is.null(object$Vc)) "Vc" else "Vp"
  z <- qnorm((1 - level) / 2, lower.tail = FALSE)
  out <- lapply(sel, function(i) {
    Xi <- X[[i]]
    eta <- as.numeric(Xi %*% object[[i]]$coefficients)
    if (length(offs[[i]])) eta <- eta + offs[[i]]
    Vi <- object[[i]][[V.type]]
    v <- rowSums(Xi * (Xi %*% Vi))
    se <- sqrt(pmax(v, 0))
    ## As confint() without newdata, a coefficient without positive variance (an
    ## unidentified one) gives no interval: rows that load on one are NA, as are rows
    ## whose variance is not positive although they load on some coefficient. Rows
    ## that load on none, e.g. a pinned parameter (whose column is zero), are exact.
    loads <- Xi != 0
    undefined <- rowSums(loads[, !(diag(Vi) > 0), drop = FALSE]) > 0 |
      (v <= 0 & rowSums(loads) > 0)
    se[which(undefined)] <- NA_real_
    width <- ifelse(se == 0, 0, z * se)
    ci <- .confint_orient(.confint_response(object, link_names[i], eta - width),
                          .confint_response(object, link_names[i], eta + width),
                          function(k) .confint_response(object, link_names[i], eta[k]),
                          check_full = level == 1 & object$family == "custom" & is.infinite(width))
    ci <- cbind(lower = ci$lower, upper = ci$upper)
    rownames(ci) <- paste0(stripped[[i]], "[", seq_len(nrow(ci)), "]")
    ci
  })
  do.call(rbind, out)
}
