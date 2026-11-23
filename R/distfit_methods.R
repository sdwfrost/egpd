## S3 methods for fitegpd objects

## Helper: check if family is bivariate (neural Bayes)
.is_bivariate <- function(family) family %in% c("begpd", "bdegpd", "bzidegpd",
                                                  "mdgpd", "zimdgpd")

## Family labels for display
.family_label <- function(family, data_dim = NULL) {
  dim_str <- if (!is.null(data_dim) && data_dim > 2L) paste0(data_dim, "-variate ") else "bivariate "
  switch(family,
    "begpd"    = "bivariate BEGPD",
    "bdegpd"   = "bivariate discrete EGPD (BDEGPD) [Experimental]",
    "bzidegpd" = "zero-inflated bivariate discrete EGPD (BZIDEGPD) [Experimental]",
    "mdgpd"    = paste0(dim_str, "MDGPD (Aka-Kratz-Naveau) [Experimental]"),
    "zimdgpd"  = paste0("zero-inflated ", dim_str, "MDGPD (ZIMDGPD) [Experimental]"),
    family
  )
}

#' @export
print.distfit <- function(x, ...) {
  if (.is_bivariate(x$family)) {
    est_label <- if (x$estimator_type == "npe") "npe" else "nbe"
    dd <- if (is.matrix(x$data)) ncol(x$data) else 2L
    cat("Fitting of ", .family_label(x$family, dd), " by neuralbayes (",
        est_label, ")\n", sep = "")
    if (x$estimator_type == "npe") {
      cat("Posterior median estimates:\n")
    } else {
      cat("Point estimates:\n")
    }
    print(round(x$estimate, 4))
    return(invisible(x))
  }

  cat("Fitting of the distribution '", x$family, "' (type ", x$type,
      ") by ", x$method, "\n", sep = "")
  if (!is.null(x$reduction)) cat(x$reduction$description, "\n")
  cat("Parameters:\n")
  print(round(x$estimate, 4))
  invisible(x)
}

#' @export
summary.distfit <- function(object, ...) {
  if (.is_bivariate(object$family)) {
    est <- object$estimate

    if (object$estimator_type == "npe" && !is.null(object$posterior_samples)) {
      post <- object$posterior_samples  # d x nsamples
      se <- apply(post, 1, sd)
      q025 <- apply(post, 1, quantile, probs = 0.025)
      q975 <- apply(post, 1, quantile, probs = 0.975)
      tab <- cbind(Median = est, Post.SD = se, `2.5%` = q025, `97.5%` = q975)
    } else {
      se <- object$sd
      tab <- cbind(Estimate = est, SE = se)
    }
    rownames(tab) <- names(est)

    dd <- if (is.matrix(object$data)) ncol(object$data) else 2L
    return(structure(list(
      table       = tab,
      fix.tab     = NULL,
      loglik      = object$loglik,
      aic         = object$aic,
      bic         = object$bic,
      n           = object$n,
      npar        = object$npar,
      family      = object$family,
      type        = object$type,
      method      = object$method,
      convergence = object$convergence,
      bernstein.m = NULL,
      cpegpd.h    = NULL,
      estimator_type = object$estimator_type,
      nsamples    = object$nsamples,
      data_dim    = dd
    ), class = "summary.distfit"))
  }

  est <- object$estimate
  se  <- object$sd
  z   <- est / se
  p   <- 2 * pnorm(-abs(z))

  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = p)
  rownames(tab) <- names(est)

  ## Include fixed args in display
  fix_tab <- NULL
  if (length(object$fix.arg) > 0) {
    fix_tab <- cbind(Value = unlist(object$fix.arg))
  }

  structure(list(
    table    = tab,
    fix.tab  = fix_tab,
    loglik   = object$loglik,
    aic      = object$aic,
    bic      = object$bic,
    n        = object$n,
    npar     = object$npar,
    reduction = object$reduction,
    family   = object$family,
    type     = object$type,
    method   = object$method,
    convergence = object$convergence,
    bernstein.m = object$bernstein.m,
    cpegpd.h = object$cpegpd.h
  ), class = "summary.distfit")
}

#' @export
print.summary.distfit <- function(x, digits = 4, ...) {
  if (.is_bivariate(x$family)) {
    est_label <- if (x$estimator_type == "npe") "npe" else "nbe"
    cat("Fitting of ", .family_label(x$family, x$data_dim), "\n", sep = "")
    cat("Method: neuralbayes (", est_label, ")", sep = "")
    if (!is.null(x$nsamples))
      cat("  [", x$nsamples, " posterior samples]", sep = "")
    cat("\n\n")

    if (x$estimator_type == "npe") {
      cat("Posterior summary:\n")
    } else {
      cat("Estimated parameters:\n")
    }
    print(round(x$table, digits))

    cat("\nNote: log-likelihood, AIC, and BIC are not available for neural estimation\n")
    cat("Number of observations: ", x$n, "\n")
    return(invisible(x))
  }

  cat("Fitting of the distribution '", x$family, "' (type ", x$type, ")\n", sep = "")
  cat("Method: ", x$method, sep = "")
  if (!is.null(x$bernstein.m)) cat(" (Bernstein degree = ", x$bernstein.m, ")", sep = "")
  if (!is.null(x$cpegpd.h)) cat(" (cpegpd discretization h = ", x$cpegpd.h, ")", sep = "")
  cat("\n\n")

  cat("Estimated parameters:\n")
  if (!is.null(x$reduction)) cat(x$reduction$description, "\n")
  if (nrow(x$table)) printCoefmat(x$table, digits = digits, na.print = "NA",
                                  P.values = TRUE, has.Pvalue = TRUE, ...)

  if (!is.null(x$fix.tab)) {
    cat("\nFixed parameters:\n")
    print(round(x$fix.tab, digits))
  }

  cat("\nConvergence: ", ifelse(x$convergence == 0, "successful", "FAILED"), "\n")
  cat("Loglikelihood: ", round(x$loglik, 2),
      "  AIC: ", round(x$aic, 2),
      "  BIC: ", round(x$bic, 2), "\n")
  cat("Number of observations: ", x$n, "\n")
  invisible(x)
}

#' @export
coef.distfit <- function(object, ...) {
  object$estimate
}

#' @export
vcov.distfit <- function(object, ...) {
  if (.is_bivariate(object$family)) {
    if (object$estimator_type == "npe" && !is.null(object$posterior_samples)) {
      V <- cov(t(object$posterior_samples))
      rownames(V) <- colnames(V) <- names(object$estimate)
      return(V)
    }
    stop("Variance-covariance matrix not available for NBE estimator (no posterior samples)",
         call. = FALSE)
  }
  if (isTRUE(object$support.ridge))
    stop("Variance-covariance matrix not available: free parameters are not separately identifiable.",
         call. = FALSE)
  if (isTRUE(object$support.change))
    stop("Variance-covariance matrix not available at a change in count support.", call. = FALSE)
  if (isTRUE(object$optim$pi.boundary))
    stop("Variance-covariance matrix not available: zero-inflation MLE on its natural boundary.",
         call. = FALSE)
  if (length(object$optim$carrier.boundary))
    stop("Variance-covariance matrix not available: carrier MLE on its natural boundary.",
         call. = FALSE)
  if (isTRUE(object$optim$bound.active))
    stop("Variance-covariance matrix not available: the MLE is on an optimizer bound.",
         call. = FALSE)
  if (isTRUE(object$optim$stalled))
    stop("Variance-covariance matrix not available: the likelihood optimizer stalled.",
         call. = FALSE)
  if (is.null(object$vcov))
    stop("Variance-covariance matrix not available (hessian was not computed)")
  object$vcov
}

#' @export
logLik.distfit <- function(object, ...) {
  if (.is_bivariate(object$family)) {
    warning("Log-likelihood is not available for neural Bayes estimation",
            call. = FALSE)
    ll <- NA_real_
    attr(ll, "df") <- object$npar
    attr(ll, "nobs") <- object$n
    class(ll) <- "logLik"
    return(ll)
  }
  ll <- object$loglik
  attr(ll, "df") <- object$npar
  attr(ll, "nobs") <- object$n
  class(ll) <- "logLik"
  ll
}

#' @export
nobs.distfit <- function(object, ...) {
  object$n
}

#' Confidence intervals for distfit parameters
#'
#' For maximum likelihood and Bernstein fits, Wald intervals computed on the scale
#' on which each parameter is estimated -- log for positive parameters (e.g.
#' \code{sigma}, \code{kappa}, \code{delta}, \code{mu}, \code{theta}), logit for
#' probabilities (e.g. \code{pi}, \code{prob}, GPIG \code{a} and \code{c}), and
#' identity for the shape \code{xi} -- except the generalised Waring \code{rho},
#' whose interval is formed for \eqn{\log(\rho - 1)} since \eqn{\rho > 1};
#' back-transformed, so they stay inside the parameter space. For neural posterior
#' estimation fits, equal-tailed posterior credible intervals.
#'
#' @param object a \code{distfit} object
#' @param parm the parameters, as names or as indices into \code{coef(object)};
#'   all estimated parameters by default. Unknown names and out-of-range indices
#'   are an error.
#' @param level the confidence level, a single probability in \eqn{[0,1]}
#' @param ... unused
#' @return A matrix with one row per parameter and columns for the lower and
#'   upper limits.
#' @export
confint.distfit <- function(object, parm, level = 0.95, ...) {
  .check_confint_level(level)
  a <- (1 - level) / 2
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), "%")
  cf <- object$estimate
  parm <- .distfit_parm(if (missing(parm)) NULL else parm, names(cf))
  if (.is_bivariate(object$family)) {
    if (object$estimator_type == "npe" && !is.null(object$posterior_samples)) {
      ## Credible intervals from posterior quantiles
      post <- object$posterior_samples  # d x nsamples
      idx <- match(parm, rownames(post))
      ci <- matrix(unlist(lapply(idx, function(i)
        stats::quantile(post[i, ], probs = c(a, 1 - a), names = FALSE))),
        ncol = 2L, byrow = TRUE, dimnames = list(parm, pct))
      return(ci)
    }
    stop("Credible intervals not available for NBE estimator (no posterior samples)",
         call. = FALSE)
  }

  se <- object$sd
  z <- qnorm(a, lower.tail = FALSE)
  ## Wald intervals on the scale the parameter was estimated on (log for positive
  ## parameters, logit for probabilities; log(rho - 1) for the Waring rho > 1),
  ## back-transformed, so they respect the parameter space; the delta-method SE on
  ## that scale is sd / (dpar / dtheta).
  tr <- .distfit_transforms(object)[parm]
  est <- cf[parm]
  lg <- tr == "log"
  lt <- tr == "logit"
  l1 <- tr == "log1"
  th <- est
  jac <- rep(1, length(est))
  th[lg] <- log(est[lg])
  jac[lg] <- est[lg]
  th[lt] <- qlogis(est[lt])
  jac[lt] <- est[lt] * (1 - est[lt])
  th[l1] <- log(est[l1] - 1)
  jac[l1] <- est[l1] - 1
  se_th <- se[parm] / jac
  back <- function(v) {
    v[lg] <- exp(v[lg])
    v[lt] <- plogis(v[lt])
    v[l1] <- 1 + exp(v[l1])
    v
  }
  ci <- cbind(back(th - z * se_th), back(th + z * se_th))
  colnames(ci) <- pct
  rownames(ci) <- parm
  ci
}

## Resolve confint()'s `parm` to parameter names: NULL gives all of them; numbers
## index `nms`, names must be among them (as confint.egpd() requires).
.distfit_parm <- function(parm, nms) {
  if (is.null(parm)) return(nms)
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(parm < 1 | parm > length(nms) | parm != round(parm)))
      stop("Numeric 'parm' out of range: the fit has ", length(nms),
           " parameters (", paste(nms, collapse = ", "), ")", call. = FALSE)
    return(nms[parm])
  }
  if (is.character(parm)) {
    bad <- setdiff(parm, nms)
    if (length(bad))
      stop(sprintf("Parameter%s %s not found; the fit has %s",
                   if (length(bad) > 1L) "s" else "",
                   paste0("'", bad, "'", collapse = ", "), paste(nms, collapse = ", ")),
           call. = FALSE)
    return(parm)
  }
  stop("'parm' must be numeric or character", call. = FALSE)
}

## Scale on which confint.distfit() forms each Wald interval of a univariate fit:
## the estimation scale used by optim() ("log", "logit" or "identity"), except the
## generalised Waring rho > 1, which uses "log1", i.e. log(rho - 1), so that the
## interval stays above 1.
.distfit_transforms <- function(object) {
  spec <- if (identical(object$method, "bernstein"))
    list(sigma = list(transform = "log"), xi = list(transform = "identity"),
         kappa = list(transform = "log"))
  else .distfit_parspec(object$type, object$family)
  if (isTRUE(object$reduction$product) || isTRUE(object$reduction$thinned))
    spec$mu <- list(transform = "log")
  if (isTRUE(object$reduction$bernoulli)) spec$p1 <- list(transform = "logit")
  if (identical(object$family, "gwaring")) spec$rho$transform <- "log1"
  vapply(names(object$estimate), function(nm) spec[[nm]]$transform, character(1))
}


#' Plot diagnostics for a distfit object
#'
#' Produces a 4-panel diagnostic plot. For univariate fits: histogram with
#' fitted density, empirical vs fitted CDF, Q-Q plot, and P-P plot. For
#' bivariate BEGPD fits: observed scatter, simulated scatter, radial Q-Q
#' plot, and posterior marginals or parameter bar chart.
#'
#' @param x a \code{distfit} object
#' @param ... additional graphical parameters
#'
#' @export
plot.distfit <- function(x, ...) {
  ## Dispatch to bivariate-specific plots
  if (x$family == "begpd") {
    return(.plot_begpd(x, ...))
  }
  if (x$family %in% c("bdegpd", "bzidegpd")) {
    return(.plot_bdegpd(x, ...))
  }
  if (x$family %in% c("mdgpd", "zimdgpd")) {
    return(.plot_mdgpd(x, ...))
  }

  obj <- x
  dat <- obj$data
  n <- obj$n
  est <- c(as.list(obj$estimate), as.list(obj$fix.arg), obj$reduction$constants)
  if (isTRUE(obj$reduction$product)) est$b <- est$mu
  is_discrete <- obj$family %in% c("degpd", "zidegpd", "cpdegpd", "pig", "zipig",
                                   "gpig", "zigpig", "bell", "zibell", "gwaring")
  is_zi <- obj$family %in% c("ziegpd", "zidegpd")
  is_pig <- obj$family %in% c("pig", "zipig")
  is_gpig <- obj$family %in% c("gpig", "zigpig")
  is_bell <- obj$family %in% c("bell", "zibell")
  is_gwaring <- obj$family == "gwaring"
  is_cpegpd <- obj$family == "cpegpd"
  is_cpdegpd <- obj$family == "cpdegpd"
  is_bernstein <- obj$method == "bernstein"

  ## Build parameter list for distribution functions
  sigma  <- est[["sigma"]]
  mu     <- if ("mu" %in% names(est)) est[["mu"]] else NA
  xi     <- est[["xi"]]
  kappa  <- if ("kappa" %in% names(est)) est[["kappa"]] else NA
  delta  <- if ("delta" %in% names(est)) est[["delta"]] else NA
  prob   <- if ("prob" %in% names(est)) est[["prob"]] else NA
  pi_val <- if ("pi" %in% names(est)) est[["pi"]] else NA
  lambda <- if ("lambda" %in% names(est)) est[["lambda"]] else NA
  a_gp   <- if ("a" %in% names(est)) est[["a"]] else NA
  b_gp   <- if ("b" %in% names(est)) est[["b"]] else NA
  c_gp   <- if ("c" %in% names(est)) est[["c"]] else NA
  theta_bell <- if ("theta" %in% names(est)) est[["theta"]] else NA
  k_gw   <- if ("k" %in% names(est)) est[["k"]] else NA
  rho_gw <- if ("rho" %in% names(est)) est[["rho"]] else NA

  ## Density, CDF, and quantile wrappers
  if (isTRUE(obj$reduction$bernoulli)) {
    dfun <- function(xv) stats::dbinom(xv, 1, est$p1)
    pfun <- function(qv) stats::pbinom(qv, 1, est$p1)
    qfun <- function(pv) stats::qbinom(pv, 1, est$p1)
  } else if (isTRUE(obj$reduction$thinned)) {
    dfun <- function(xv) stats::dpois(xv, est$mu)
    pfun <- function(qv) stats::ppois(qv, est$mu)
    qfun <- function(pv) stats::qpois(pv, est$mu)
  } else if (is_gwaring) {
    dfun <- function(xv) dgwaring(xv, mu = mu, k = k_gw, rho = rho_gw)
    pfun <- function(qv) pgwaring(qv, mu = mu, k = k_gw, rho = rho_gw)
    qfun <- function(pv) qgwaring(pv, mu = mu, k = k_gw, rho = rho_gw)
  } else if (is_pig) {
    ## The package's native PIG/ZIPIG functions, as in the fit (gamlss.dist's
    ## break down near the Poisson limit sigma -> 0).
    dfun <- function(xv) if (obj$family == "pig")
      dpig(xv, mu = mu, sigma = sigma)
      else dzipig(xv, mu = mu, sigma = sigma, pi = pi_val)
    pfun <- function(qv) if (obj$family == "pig")
      ppig(qv, mu = mu, sigma = sigma)
      else pzipig(qv, mu = mu, sigma = sigma, pi = pi_val)
    qfun <- function(pv) if (obj$family == "pig")
      qpig(pv, mu = mu, sigma = sigma)
      else qzipig(pv, mu = mu, sigma = sigma, pi = pi_val)
  } else if (is_gpig) {
    dfun <- function(xv) if (obj$family == "gpig")
      dgpig(xv, a = a_gp, b = b_gp, c = c_gp)
      else dzigpig(xv, a = a_gp, b = b_gp, c = c_gp, pi = pi_val)
    pfun <- function(qv) if (obj$family == "gpig")
      pgpig(qv, a = a_gp, b = b_gp, c = c_gp)
      else pzigpig(qv, a = a_gp, b = b_gp, c = c_gp, pi = pi_val)
    qfun <- function(pv) if (obj$family == "gpig")
      qgpig(pv, a = a_gp, b = b_gp, c = c_gp)
      else qzigpig(pv, a = a_gp, b = b_gp, c = c_gp, pi = pi_val)
  } else if (is_bell) {
    dfun <- function(xv) if (obj$family == "bell")
      dbell(xv, theta = theta_bell)
      else dzibell(xv, theta = theta_bell, pi = pi_val)
    pfun <- function(qv) if (obj$family == "bell")
      pbell(qv, theta = theta_bell)
      else pzibell(qv, theta = theta_bell, pi = pi_val)
    qfun <- function(pv) if (obj$family == "bell")
      qbell(pv, theta = theta_bell)
      else qzibell(pv, theta = theta_bell, pi = pi_val)
  } else if (is_cpdegpd) {
    dfun <- function(xv) dcpdegpd(xv, lambda = lambda, prob = prob, kappa = kappa,
                                    delta = delta, sigma = sigma, xi = xi,
                                    type = obj$type)
    pfun <- function(qv) pcpdegpd(qv, lambda = lambda, prob = prob, kappa = kappa,
                                    delta = delta, sigma = sigma, xi = xi,
                                    type = obj$type)
    qfun <- function(pv) qcpdegpd(pv, lambda = lambda, prob = prob, kappa = kappa,
                                    delta = delta, sigma = sigma, xi = xi,
                                    type = obj$type)
  } else if (is_cpegpd) {
    cp_h <- obj$cpegpd.h
    dfun <- function(xv) dcpegpd(xv, lambda = lambda, prob = prob, kappa = kappa,
                                   delta = delta, sigma = sigma, xi = xi,
                                   type = obj$type, h = cp_h)
    pfun <- function(qv) pcpegpd(qv, lambda = lambda, prob = prob, kappa = kappa,
                                   delta = delta, sigma = sigma, xi = xi,
                                   type = obj$type, h = cp_h)
    qfun <- function(pv) qcpegpd(pv, lambda = lambda, prob = prob, kappa = kappa,
                                   delta = delta, sigma = sigma, xi = xi,
                                   type = obj$type, h = cp_h)
  } else if (is_bernstein) {
    ## Bernstein model: custom density/CDF/quantile
    wt <- obj$bernstein.weights
    m  <- obj$bernstein.m
    dfun <- function(xv) .bernstein_full_density(xv, sigma, xi, kappa, wt, m)
    pfun <- function(qv) .bernstein_full_cdf(qv, sigma, xi, kappa, wt, m)
    qfun <- function(pv) .bernstein_full_quantile(pv, sigma, xi, kappa, wt, m)
  } else if (is_zi && is_discrete) {
    dfun <- function(xv) dzidiscegpd(xv, pi = pi_val, prob = prob, kappa = kappa,
                                      delta = delta, sigma = sigma, xi = xi, type = obj$type)
    pfun <- function(qv) pzidiscegpd(qv, pi = pi_val, prob = prob, kappa = kappa,
                                      delta = delta, sigma = sigma, xi = xi, type = obj$type)
    qfun <- function(pv) qzidiscegpd(pv, pi = pi_val, prob = prob, kappa = kappa,
                                      delta = delta, sigma = sigma, xi = xi, type = obj$type)
  } else if (is_zi) {
    dfun <- function(xv) dziegpd(xv, pi = pi_val, prob = prob, kappa = kappa,
                                  delta = delta, sigma = sigma, xi = xi, type = obj$type)
    pfun <- function(qv) pziegpd(qv, pi = pi_val, prob = prob, kappa = kappa,
                                  delta = delta, sigma = sigma, xi = xi, type = obj$type)
    qfun <- function(pv) qziegpd(pv, pi = pi_val, prob = prob, kappa = kappa,
                                  delta = delta, sigma = sigma, xi = xi, type = obj$type)
  } else if (is_discrete) {
    dfun <- function(xv) ddiscegpd(xv, prob = prob, kappa = kappa, delta = delta,
                                    sigma = sigma, xi = xi, type = obj$type)
    pfun <- function(qv) pdiscegpd(qv, prob = prob, kappa = kappa, delta = delta,
                                    sigma = sigma, xi = xi, type = obj$type)
    qfun <- function(pv) qdiscegpd(pv, prob = prob, kappa = kappa, delta = delta,
                                    sigma = sigma, xi = xi, type = obj$type)
  } else {
    dfun <- function(xv) degpd_density(xv, prob = prob, kappa = kappa, delta = delta,
                                        sigma = sigma, xi = xi, type = obj$type)
    pfun <- function(qv) pegpd(qv, prob = prob, kappa = kappa, delta = delta,
                                sigma = sigma, xi = xi, type = obj$type)
    qfun <- function(pv) qegpd(pv, prob = prob, kappa = kappa, delta = delta,
                                sigma = sigma, xi = xi, type = obj$type)
  }

  ## Randomized PIT for discrete Q-Q and P-P plots (Dunn & Smyth, 1996)
  pit <- NULL
  if (is_discrete) {
    dat_sort_pit <- sort(dat)
    ## Save and restore RNG state so plotting doesn't affect user's stream
    rng_existed <- exists(".Random.seed", envir = .GlobalEnv)
    if (rng_existed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(1)
    a_pit <- ifelse(dat_sort_pit == 0, 0, pfun(dat_sort_pit - 1))
    b_pit <- pfun(dat_sort_pit)
    pit <- runif(n, a_pit, b_pit)
    pit <- pmax(pmin(pit, 1 - 1e-10), 1e-10)
    if (rng_existed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }

  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(op))

  ## Panel 1: Histogram + density
  if (is_cpegpd) {
    hist(dat, breaks = "FD", freq = FALSE, main = "Histogram and Fitted Density",
         xlab = "x", col = "lightblue", border = "grey")
    xseq <- seq(0, max(dat), by = cp_h)
    dens <- dfun(xseq) / cp_h
    lines(xseq, dens, col = "red", lwd = 2)
  } else if (is_discrete) {
    tab <- table(dat)
    vals <- as.numeric(names(tab))
    ## the fitted pmf goes at the bar centres that barplot() returns, not at 1..K
    bp <- barplot(tab / n, names.arg = vals, main = "Empirical vs Fitted PMF",
                  xlab = "x", ylab = "Probability", col = "lightblue", border = "grey")
    points(bp, dfun(vals), pch = 16, col = "red", type = "b")
  } else {
    hist(dat, breaks = "FD", freq = FALSE, main = "Histogram and Fitted Density",
         xlab = "x", col = "lightblue", border = "grey")
    xseq <- seq(max(min(dat), 1e-6), max(dat), length.out = 200)
    lines(xseq, dfun(xseq), col = "red", lwd = 2)
  }

  ## Panel 2: Empirical vs fitted CDF
  dat_sort <- sort(dat)
  emp_cdf <- (1:n) / (n + 1)

  if (is_discrete) {
    plot(stepfun(dat_sort, c(0, emp_cdf)), main = "Empirical vs Fitted CDF",
         xlab = "x", ylab = "CDF", col = "black", do.points = FALSE)
    xseq <- seq(min(dat), max(dat))
    lines(xseq, pfun(xseq), col = "red", lwd = 2, type = "s")
  } else {
    plot(dat_sort, emp_cdf, pch = ".", main = "Empirical vs Fitted CDF",
         xlab = "x", ylab = "CDF")
    xseq <- seq(min(dat), max(dat), length.out = 200)
    lines(xseq, pfun(xseq), col = "red", lwd = 2)
  }

  ## Panel 3: Q-Q plot
  if (is_discrete) {
    qres <- sort(qnorm(pit))
    theo_norm <- qnorm(emp_cdf)
    plot(theo_norm, qres, pch = 1, cex = 0.5,
         main = "Q-Q Plot (Quantile Residuals)",
         xlab = "Normal Theoretical Quantiles",
         ylab = "Randomized Quantile Residuals")
    abline(0, 1, col = "red", lwd = 2)
  } else {
    theo_q <- qfun(emp_cdf)
    plot(theo_q, dat_sort, pch = 1, cex = 0.5,
         main = "Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
    abline(0, 1, col = "red", lwd = 2)
  }

  ## Panel 4: P-P plot
  if (is_discrete) {
    plot(sort(pit), emp_cdf, pch = 1, cex = 0.5,
         main = "P-P Plot (Randomized PIT)",
         xlab = "Randomized PIT", ylab = "Empirical CDF")
    abline(0, 1, col = "red", lwd = 2)
  } else {
    theo_p <- pfun(dat_sort)
    plot(theo_p, emp_cdf, pch = 1, cex = 0.5,
         main = "P-P Plot", xlab = "Theoretical CDF", ylab = "Empirical CDF")
    abline(0, 1, col = "red", lwd = 2)
  }

  invisible(obj)
}

## ---- deprecated S3 shims -----------------------------------------------------
## The fit object's class was renamed from "fitegpd" to "distfit". New fits carry
## class "distfit"; these shims keep S3 dispatch working for objects saved by
## older versions (class "fitegpd") by relabelling and re-dispatching. The
## fitegpd() function itself (see R/distfit.R) emits the deprecation message.

#' @export
print.fitegpd <- function(x, ...) { class(x) <- "distfit"; print(x, ...) }
#' @export
summary.fitegpd <- function(object, ...) { class(object) <- "distfit"; summary(object, ...) }
#' @export
print.summary.fitegpd <- function(x, ...) { class(x) <- "summary.distfit"; print(x, ...) }
#' @export
coef.fitegpd <- function(object, ...) { class(object) <- "distfit"; coef(object, ...) }
#' @export
vcov.fitegpd <- function(object, ...) { class(object) <- "distfit"; vcov(object, ...) }
#' @export
logLik.fitegpd <- function(object, ...) { class(object) <- "distfit"; logLik(object, ...) }
#' @export
nobs.fitegpd <- function(object, ...) { class(object) <- "distfit"; nobs(object, ...) }
#' @export
confint.fitegpd <- function(object, ...) { class(object) <- "distfit"; confint(object, ...) }
#' @export
plot.fitegpd <- function(x, ...) { class(x) <- "distfit"; plot(x, ...) }
