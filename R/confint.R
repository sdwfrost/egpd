## confint and vcov methods for egpd objects

# Detect inverse link from coefficient name prefix
.detect_invlink <- function(name) {
  if (grepl("^logit", name))
    return(plogis)
  if (grepl("^probit", name))
    return(qnorm)
  if (grepl("^log", name))
    return(exp)
  identity
}

# Strip link prefix to get parameter name
.strip_link <- function(name) {
  name <- sub("^logit", "", name)
  name <- sub("^probit", "", name)
  name <- sub("^log", "", name)
  name
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
      if (!is.na(intercept_pos)) {
        # Find the global index
        offset <- sum(vapply(seq_len(i - 1), function(j) length(object[[j]]$coefficients), integer(1)))
        idx <- c(idx, offset + intercept_pos)
      }
    }
    return(idx)
  }

  if (is.numeric(parm)) {
    if (any(parm < 1 | parm > length(coef_names)))
      stop("Numeric 'parm' out of range")
    return(as.integer(parm))
  }

  if (is.character(parm)) {
    idx <- integer(0)
    for (p in parm) {
      # Try matching against stripped link names
      m <- match(p, stripped)
      if (!is.na(m)) {
        # Return the intercept index for this parameter
        offset <- sum(vapply(seq_len(m - 1), function(j) length(object[[j]]$coefficients), integer(1)))
        sub_coefs <- names(object[[m]]$coefficients)
        intercept_pos <- match("(Intercept)", sub_coefs)
        if (is.na(intercept_pos)) intercept_pos <- 1L
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

  # Evaluate at initial bracket endpoint, expand if needed (up to 5 doublings)
  f_hi <- suppressWarnings(root_fn(hi))
  for (j in seq_len(5)) {
    if (is.finite(f_hi) && f_hi > 0) break
    hi <- mle + (hi - mle) * 2
    f_hi <- suppressWarnings(root_fn(hi))
  }

  if (!is.finite(f_hi) || f_hi <= 0) {
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

# Profile one coefficient
.profile_one_param <- function(object, k, level) {
  beta_hat <- object$coefficients
  likdata <- object$likdata
  likfns <- object$likfns
  cutoff <- qchisq(level, 1)

  nllh_mle <- suppressWarnings(.nllh.nopen(beta_hat, likdata, likfns))
  vk <- diag(object$Vp)[k]
  se_k <- if (vk > 0) sqrt(vk) else abs(beta_hat[k]) * 0.1 + 0.1
  nuisance_start <- beta_hat[-k]

  profile_deviance <- function(v) {
    # Always start from MLE nuisance values for robustness
    opt <- suppressWarnings(tryCatch(
      optim(par = nuisance_start, fn = .nllh.nopen, gr = .grad.nopen,
            likdata = likdata, likfns = likfns,
            id = k, vals = v, method = "BFGS",
            control = list(maxit = 200)),
      error = function(e) NULL
    ))

    if (is.null(opt) || !is.finite(opt$value)) {
      # Fallback to Nelder-Mead
      opt <- suppressWarnings(tryCatch(
        optim(par = nuisance_start, fn = .nllh.nopen,
              likdata = likdata, likfns = likfns,
              id = k, vals = v, method = "Nelder-Mead",
              control = list(maxit = 500)),
        error = function(e) NULL
      ))
    }

    if (is.null(opt) || !is.finite(opt$value)) return(Inf)

    2 * (opt$value - nllh_mle)
  }

  root_fn <- function(v) profile_deviance(v) - cutoff

  lower <- .find_profile_root(root_fn, beta_hat[k], se_k, -1, cutoff)
  upper <- .find_profile_root(root_fn, beta_hat[k], se_k, +1, cutoff)

  c(lower, upper)
}

#' Variance-covariance matrix for egpd fits
#'
#' @param object A fitted \code{egpd} object.
#' @param ... Not used.
#' @return The variance-covariance matrix \code{object$Vp} with row and column
#'   names matching the coefficient names.
#' @export
vcov.egpd <- function(object, ...) {
  V <- object$Vp
  nms <- names(object$coefficients)
  rownames(V) <- colnames(V) <- nms
  V
}

#' Confidence intervals for egpd model parameters
#'
#' @param object A fitted \code{egpd} object.
#' @param parm Which parameters to compute CIs for.  \code{NULL} (default)
#'   selects the intercept of each distributional parameter.  Can also be a
#'   character vector of stripped parameter names (e.g. \code{"scale"},
#'   \code{"shape"}) or an integer vector of coefficient indices.
#' @param level Confidence level (default 0.95).
#' @param method Either \code{"wald"} (default) or \code{"profile"}.  Profile
#'   likelihood CIs are only available for intercept-only models (no smooth
#'   terms or covariates beyond an intercept).
#' @param ... Not used.
#' @return A matrix with columns \code{"lower"} and \code{"upper"}, on the
#'   response (natural) scale.  Rows are named by parameter.
#' @export
confint.egpd <- function(object, parm = NULL, level = 0.95, method = c("wald", "profile"), ...) {
  method <- match.arg(method)
  idx <- .resolve_parm(object, parm)

  coef_names <- names(object$coefficients)
  link_names <- names(object)[seq_len(object$ngam)]

  # Determine which link name each coefficient belongs to
  param_labels <- character(length(idx))
  invlinks <- vector("list", length(idx))
  for (j in seq_along(idx)) {
    k <- idx[j]
    # Find which sub-model this coefficient belongs to
    cumlen <- cumsum(vapply(seq_len(object$ngam), function(i) length(object[[i]]$coefficients), integer(1)))
    sub_idx <- which(k <= cumlen)[1]
    param_labels[j] <- .strip_link(link_names[sub_idx])
    invlinks[[j]] <- .detect_invlink(link_names[sub_idx])
  }

  if (method == "wald") {
    z <- qnorm(1 - (1 - level) / 2)
    coefs <- object$coefficients[idx]
    vars <- diag(object$Vp)[idx]

    ci <- matrix(NA_real_, nrow = length(idx), ncol = 2)
    colnames(ci) <- c("lower", "upper")
    rownames(ci) <- param_labels

    for (j in seq_along(idx)) {
      if (vars[j] <= 0) next  # non-positive variance: CI undefined
      se_j <- sqrt(vars[j])
      ci_link <- coefs[j] + c(-1, 1) * z * se_j
      ci[j, ] <- invlinks[[j]](ci_link)
    }

    return(ci)
  }

  # Profile method
  if (length(object$gotsmooth) > 0)
    stop("Profile likelihood CIs are only available for intercept-only models (no smooth terms)")

  # Check that each sub-model is intercept-only
  for (i in seq_len(object$ngam)) {
    if (length(object[[i]]$coefficients) > 1)
      stop("Profile likelihood CIs are only available for intercept-only models (no covariates beyond intercept)")
  }

  ci <- matrix(NA_real_, nrow = length(idx), ncol = 2)
  colnames(ci) <- c("lower", "upper")
  rownames(ci) <- param_labels

  for (j in seq_along(idx)) {
    k <- idx[j]
    bounds <- .profile_one_param(object, k, level)
    ci[j, ] <- invlinks[[j]](bounds)
  }

  ci
}
