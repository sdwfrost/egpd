## Univariate distribution fitting for EGPD families
## MLE and Bernstein polynomial fitting

#' Fit EGPD distribution to data
#'
#' Maximum likelihood or Bernstein polynomial fitting of EGPD, discrete EGPD,
#' zero-inflated EGPD, or zero-inflated discrete EGPD distributions.
#'
#' @param x numeric vector of observations
#' @param type integer 1-6 specifying the G-transformation type
#' @param family character: "egpd", "degpd", "ziegpd", or "zidegpd"
#' @param method character: "mle" or "bernstein" (Bernstein only for family="egpd")
#' @param start named list of starting values, or NULL for automatic
#' @param fix.arg named list of fixed parameters (not optimized)
#' @param optim.method optimization method passed to \code{\link{optim}}
#' @param hessian logical: compute standard errors via Hessian?
#' @param bernstein.m integer: Bernstein polynomial degree (method="bernstein" only)
#' @param ... additional arguments passed to \code{\link{optim}}
#'
#' @return An object of class \code{"fitegpd"} with components:
#' \describe{
#'   \item{estimate}{named vector of parameter estimates}
#'   \item{sd}{named vector of standard errors (NA if Hessian not computed)}
#'   \item{vcov}{variance-covariance matrix on natural scale}
#'   \item{loglik}{maximized log-likelihood}
#'   \item{aic}{Akaike information criterion}
#'   \item{bic}{Bayesian information criterion}
#'   \item{n}{number of observations}
#'   \item{npar}{number of estimated parameters}
#'   \item{data}{the input data vector}
#'   \item{type}{the G-transformation type}
#'   \item{family}{the distribution family}
#'   \item{method}{the fitting method}
#'   \item{fix.arg}{list of fixed arguments}
#'   \item{convergence}{convergence code from optim (0 = success)}
#'   \item{optim}{full optim output}
#'   \item{call}{the matched call}
#'   \item{bernstein.m}{Bernstein degree (NULL for method="mle")}
#'   \item{bernstein.weights}{Bernstein weights (NULL for method="mle")}
#' }
#'
#' @examples
#' \dontrun{
#' x <- regpd(500, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
#' fit <- fitegpd(x, type = 1)
#' summary(fit)
#' plot(fit)
#' }
#'
#' @export
fitegpd <- function(x, type = 1, family = c("egpd", "degpd", "ziegpd", "zidegpd", "cpegpd"),
                    method = c("mle", "bernstein"), start = NULL, fix.arg = NULL,
                    optim.method = "Nelder-Mead", hessian = TRUE,
                    bernstein.m = 8, cpegpd.h = 0.2, ...) {
  cl <- match.call()
  family <- match.arg(family)
  method <- match.arg(method)

  ## Validate inputs
  if (!is.numeric(x) || length(x) < 2)
    stop("'x' must be a numeric vector with at least 2 observations")
  if (!type %in% 1:6)
    stop("'type' must be an integer from 1 to 6")
  if (method == "bernstein" && family != "egpd")
    stop("Bernstein method is only available for family='egpd'")
  if (family %in% c("degpd", "zidegpd") && any(x != floor(x) | x < 0))
    warning("Discrete family expects non-negative integer data")
  if (family %in% c("ziegpd", "zidegpd") && !any(x == 0))
    warning("Zero-inflated family but no zeros observed in data")
  if (family == "cpegpd" && any(x < 0))
    warning("cpegpd family expects non-negative data")

  ## Dispatch to Bernstein fitting

  if (method == "bernstein") {
    return(.fitegpd_bernstein(x, type = type, start = start, fix.arg = fix.arg,
                              m = bernstein.m, optim.method = optim.method,
                              hessian = hessian, call = cl, ...))
  }

  ## Parameter specification
  spec <- .fitegpd_parspec(type, family)

  ## Remove fixed args from spec
  free_spec <- spec
  if (!is.null(fix.arg)) {
    bad <- setdiff(names(fix.arg), names(spec))
    if (length(bad) > 0)
      stop("Unknown fixed parameters: ", paste(bad, collapse = ", "))
    free_spec <- spec[setdiff(names(spec), names(fix.arg))]
  }
  if (length(free_spec) == 0)
    stop("No free parameters to estimate")

  ## Starting values
  if (is.null(start)) {
    start_vals <- .fitegpd_start(x, type, family, free_spec)
  } else {
    bad <- setdiff(names(start), names(free_spec))
    if (length(bad) > 0)
      stop("Start values for fixed or unknown parameters: ", paste(bad, collapse = ", "))
    start_vals <- .fitegpd_start(x, type, family, free_spec)
    for (nm in names(start)) start_vals[[nm]] <- start[[nm]]
  }

  ## Transform to unconstrained scale
  start_par <- unlist(start_vals)
  theta0 <- .par_to_theta(start_par, free_spec)

  ## Optimize
  opt <- optim(theta0, fn = .fitegpd_nll, x = x, type = type, family = family,
               fix.arg = fix.arg, spec = free_spec,
               h = if (family == "cpegpd") cpegpd.h else NULL,
               method = optim.method, hessian = hessian, ...)

  ## Back-transform estimates
  estimate <- .theta_to_par(opt$par, free_spec)

  ## SE via delta method
  se <- rep(NA_real_, length(estimate))
  names(se) <- names(estimate)
  V_par <- NULL

  if (hessian && !is.null(opt$hessian)) {
    V_theta <- tryCatch(solve(opt$hessian), error = function(e) NULL)
    if (!is.null(V_theta)) {
      jac <- .fitegpd_jacobian(opt$par, free_spec)
      J <- diag(jac, nrow = length(jac))
      V_par <- J %*% V_theta %*% t(J)
      rownames(V_par) <- colnames(V_par) <- names(estimate)
      se <- sqrt(pmax(diag(V_par), 0))
    }
  }

  ## Assemble result
  npar <- length(estimate)
  ll <- -opt$value
  n <- length(x)

  structure(list(
    estimate  = estimate,
    sd        = se,
    vcov      = V_par,
    loglik    = ll,
    aic       = -2 * ll + 2 * npar,
    bic       = -2 * ll + log(n) * npar,
    n         = n,
    npar      = npar,
    data      = x,
    type      = type,
    family    = family,
    method    = method,
    fix.arg   = fix.arg,
    convergence = opt$convergence,
    optim     = opt,
    call      = cl,
    bernstein.m       = NULL,
    bernstein.weights = NULL,
    cpegpd.h  = if (family == "cpegpd") cpegpd.h else NULL
  ), class = "fitegpd")
}


## ---- Internal helpers ----

#' Parameter specification for each type and family
#' @noRd
.fitegpd_parspec <- function(type, family) {
  ## Base GPD parameters
  spec <- list(
    sigma = list(transform = "log"),
    xi    = list(transform = "identity")
  )

  ## G-transformation parameters by type
  g_params <- switch(as.character(type),
    "1" = list(kappa = list(transform = "log")),
    "2" = list(kappa = list(transform = "log")),
    "3" = list(kappa = list(transform = "log")),
    "4" = list(delta = list(transform = "log")),
    "5" = list(delta = list(transform = "log"),
               kappa = list(transform = "log")),
    "6" = list(kappa = list(transform = "log"),
               delta = list(transform = "log"),
               prob  = list(transform = "logit"))
  )

  spec <- c(spec, g_params)

  ## Zero-inflation parameter
  if (family %in% c("ziegpd", "zidegpd")) {
    spec <- c(spec, list(pi = list(transform = "logit")))
  }

  ## Compound Poisson rate parameter
  if (family == "cpegpd") {
    spec <- c(spec, list(lambda = list(transform = "log")))
  }

  spec
}


#' Transform natural-scale parameters to unconstrained scale
#' @noRd
.par_to_theta <- function(par, spec) {
  theta <- numeric(length(par))
  names(theta) <- names(par)
  for (nm in names(par)) {
    theta[nm] <- switch(spec[[nm]]$transform,
      "log"      = log(par[nm]),
      "logit"    = qlogis(par[nm]),
      "identity" = par[nm]
    )
  }
  theta
}


#' Transform unconstrained parameters to natural scale
#' @noRd
.theta_to_par <- function(theta, spec) {
  par <- numeric(length(theta))
  names(par) <- names(theta)
  for (nm in names(theta)) {
    par[nm] <- switch(spec[[nm]]$transform,
      "log"      = exp(theta[nm]),
      "logit"    = plogis(theta[nm]),
      "identity" = theta[nm]
    )
  }
  par
}


#' Jacobian diagonal of back-transformation (for delta method)
#' @noRd
.fitegpd_jacobian <- function(theta, spec) {
  jac <- numeric(length(theta))
  names(jac) <- names(theta)
  for (nm in names(theta)) {
    jac[nm] <- switch(spec[[nm]]$transform,
      "log" = exp(theta[nm]),
      "logit" = {
        p <- plogis(theta[nm])
        p * (1 - p)
      },
      "identity" = 1
    )
  }
  jac
}


#' Automatic starting values
#' @noRd
.fitegpd_start <- function(x, type, family, spec) {
  ## Remove non-finite values for moment estimation
  xclean <- x[is.finite(x)]
  if (length(xclean) < 2) xclean <- c(0.1, 0.2)

  ## Use positive data for GPD moment estimation
  if (family %in% c("degpd", "zidegpd")) {
    xpos <- xclean[xclean > 0] + 0.5
  } else {
    xpos <- xclean[xclean > 0]
  }
  if (length(xpos) < 2) xpos <- abs(xclean) + 0.01

  m <- mean(xpos)
  v <- var(xpos)
  if (!is.finite(m) || !is.finite(v) || v <= 0) { m <- 1; v <- 1 }

  ## For cpegpd, the data is S = X_1 + ... + X_N where N ~ Poisson(lambda).
  ## E[S|S>0] = lambda*E[X]/(1-exp(-lambda)). Adjust moments to target
  ## the individual severity rather than the aggregate.
  if (family == "cpegpd") {
    p0 <- max(min(mean(xclean == 0), 0.99), 0.01)
    lambda_init <- -log(p0)
    ## Rough per-event mean: E[X] ~ E[S] / lambda = mean(xclean) / lambda
    per_event_mean <- mean(xclean) / lambda_init
    m <- max(per_event_mean, 0.01)
    v <- max(v / lambda_init, 0.01)  # Var(S) ~ lambda * E[X^2], rough scale
  }

  ## GPD moment estimates
  xi_start <- 0.5 * (m^2 / v - 1)
  xi_start <- max(min(xi_start, 0.5), -0.5)
  sigma_start <- m * (1 - xi_start)
  sigma_start <- max(sigma_start, 0.01)

  ## Safety: ensure sigma + xi * max(x) > 0
  xmax <- max(xclean)
  if (is.finite(xmax) && xmax > 0 &&
      is.finite(sigma_start + xi_start * xmax) &&
      sigma_start + xi_start * xmax <= 0) {
    xi_start <- max(xi_start, -sigma_start / (2 * xmax))
  }

  start <- list()
  for (nm in names(spec)) {
    start[[nm]] <- switch(nm,
      "sigma" = sigma_start,
      "xi"    = xi_start,
      "kappa" = 1,
      "delta" = 1,
      "prob"  = 0.5,
      "pi"    = max(min(mean(x == 0), 0.99), 0.01),
      "lambda" = -log(max(min(mean(x == 0), 0.99), 0.01)),
      1 # fallback
    )
  }

  start
}


#' Negative log-likelihood for fitegpd
#' @noRd
.fitegpd_nll <- function(theta, x, type, family, fix.arg, spec, h = NULL) {
  ## Back-transform to natural scale
  par <- .theta_to_par(theta, spec)

  ## Merge with fixed args
  all_par <- as.list(c(par, unlist(fix.arg)))

  sigma  <- all_par[["sigma"]]
  xi     <- all_par[["xi"]]
  kappa  <- if ("kappa" %in% names(all_par)) all_par[["kappa"]] else NA
  delta  <- if ("delta" %in% names(all_par)) all_par[["delta"]] else NA
  prob   <- if ("prob" %in% names(all_par)) all_par[["prob"]] else NA
  pi_val <- if ("pi" %in% names(all_par)) all_par[["pi"]] else NA
  lambda <- if ("lambda" %in% names(all_par)) all_par[["lambda"]] else NA

  ## Compute NLL based on family
  nll <- tryCatch({
    switch(family,
      "egpd" = {
        ll <- degpd_density(x, prob = prob, kappa = kappa, delta = delta,
                            sigma = sigma, xi = xi, type = type, log = TRUE)
        -sum(ll)
      },
      "degpd" = {
        pmf <- ddiscegpd(x, prob = prob, kappa = kappa, delta = delta,
                         sigma = sigma, xi = xi, type = type)
        pmf <- pmax(pmf, .Machine$double.xmin)
        -sum(log(pmf))
      },
      "ziegpd" = {
        ll <- dziegpd(x, pi = pi_val, prob = prob, kappa = kappa, delta = delta,
                      sigma = sigma, xi = xi, type = type, log = TRUE)
        -sum(ll)
      },
      "zidegpd" = {
        pmf <- dzidiscegpd(x, pi = pi_val, prob = prob, kappa = kappa, delta = delta,
                           sigma = sigma, xi = xi, type = type)
        pmf <- pmax(pmf, .Machine$double.xmin)
        -sum(log(pmf))
      },
      "cpegpd" = {
        ll <- dcpegpd(x, lambda = lambda, prob = prob, kappa = kappa,
                      delta = delta, sigma = sigma, xi = xi, type = type,
                      h = h, log = TRUE)
        -sum(ll)
      }
    )
  }, error = function(e) Inf)

  if (!is.finite(nll)) return(1e20)
  nll
}
