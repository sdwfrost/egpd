## S3 methods for fitegpd objects

#' @export
print.fitegpd <- function(x, ...) {
  cat("Fitting of the distribution '", x$family, "' (type ", x$type,
      ") by ", x$method, "\n", sep = "")
  cat("Parameters:\n")
  print(round(x$estimate, 4))
  invisible(x)
}

#' @export
summary.fitegpd <- function(object, ...) {
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
    family   = object$family,
    type     = object$type,
    method   = object$method,
    convergence = object$convergence,
    bernstein.m = object$bernstein.m,
    cpegpd.h = object$cpegpd.h
  ), class = "summary.fitegpd")
}

#' @export
print.summary.fitegpd <- function(x, digits = 4, ...) {
  cat("Fitting of the distribution '", x$family, "' (type ", x$type, ")\n", sep = "")
  cat("Method: ", x$method, sep = "")
  if (!is.null(x$bernstein.m)) cat(" (Bernstein degree = ", x$bernstein.m, ")", sep = "")
  if (!is.null(x$cpegpd.h)) cat(" (cpegpd discretization h = ", x$cpegpd.h, ")", sep = "")
  cat("\n\n")

  cat("Estimated parameters:\n")
  printCoefmat(x$table, digits = digits, na.print = "NA",
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
coef.fitegpd <- function(object, ...) {
  object$estimate
}

#' @export
vcov.fitegpd <- function(object, ...) {
  if (is.null(object$vcov))
    stop("Variance-covariance matrix not available (hessian was not computed)")
  object$vcov
}

#' @export
logLik.fitegpd <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- object$npar
  attr(ll, "nobs") <- object$n
  class(ll) <- "logLik"
  ll
}

#' @export
nobs.fitegpd <- function(object, ...) {
  object$n
}

#' @export
confint.fitegpd <- function(object, parm, level = 0.95, ...) {
  cf <- object$estimate
  se <- object$sd
  if (missing(parm)) parm <- names(cf)
  parm <- parm[parm %in% names(cf)]
  a <- (1 - level) / 2
  z <- qnorm(1 - a)
  ci <- cbind(cf[parm] - z * se[parm], cf[parm] + z * se[parm])
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), "%")
  colnames(ci) <- pct
  rownames(ci) <- parm
  ci
}


#' Plot diagnostics for a fitegpd object
#'
#' Produces a 4-panel diagnostic plot: histogram with fitted density,
#' empirical vs fitted CDF, Q-Q plot, and P-P plot.
#'
#' @param x a \code{fitegpd} object
#' @param ... additional graphical parameters
#'
#' @export
plot.fitegpd <- function(x, ...) {
  obj <- x
  dat <- obj$data
  n <- obj$n
  est <- as.list(c(obj$estimate, unlist(obj$fix.arg)))
  is_discrete <- obj$family %in% c("degpd", "zidegpd", "cpdegpd")
  is_zi <- obj$family %in% c("ziegpd", "zidegpd")
  is_cpegpd <- obj$family == "cpegpd"
  is_cpdegpd <- obj$family == "cpdegpd"
  is_bernstein <- obj$method == "bernstein"

  ## Build parameter list for distribution functions
  sigma  <- est[["sigma"]]
  xi     <- est[["xi"]]
  kappa  <- if ("kappa" %in% names(est)) est[["kappa"]] else NA
  delta  <- if ("delta" %in% names(est)) est[["delta"]] else NA
  prob   <- if ("prob" %in% names(est)) est[["prob"]] else NA
  pi_val <- if ("pi" %in% names(est)) est[["pi"]] else NA
  lambda <- if ("lambda" %in% names(est)) est[["lambda"]] else NA

  ## Density, CDF, and quantile wrappers
  if (is_cpdegpd) {
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
    barplot(tab / n, names.arg = vals, main = "Empirical vs Fitted PMF",
            xlab = "x", ylab = "Probability", col = "lightblue", border = "grey")
    points(seq_along(vals), dfun(vals), pch = 16, col = "red", type = "b")
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
  theo_q <- qfun(emp_cdf)
  plot(theo_q, dat_sort, pch = 1, cex = 0.5,
       main = "Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  abline(0, 1, col = "red", lwd = 2)

  ## Panel 4: P-P plot
  theo_p <- pfun(dat_sort)
  plot(theo_p, emp_cdf, pch = 1, cex = 0.5,
       main = "P-P Plot", xlab = "Theoretical CDF", ylab = "Empirical CDF")
  abline(0, 1, col = "red", lwd = 2)

  invisible(obj)
}
