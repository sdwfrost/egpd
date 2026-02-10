## Bernstein polynomial fitting for continuous EGPD
##
## Replaces the parametric G-transformation with a flexible
## Bernstein polynomial density on [0,1], applied after a power
## transform u^kappa. The full density is:
##   f(x) = b(H(x)^kappa) * kappa * H(x)^(kappa-1) * h(x)
## where H is the GPD CDF, h is the GPD density, and b is the
## Bernstein polynomial density with softmax-parameterized weights.


#' Bernstein polynomial density on [0,1]
#'
#' @param u values in (0,1)
#' @param weights probability weights (length m, sum to 1)
#' @param m Bernstein degree
#' @return density values
#' @noRd
.bernstein_density <- function(u, weights, m) {
  ## b(u) = sum_{j=0}^{m-1} w_j * dbeta(u, j+1, m-j)
  d <- numeric(length(u))
  for (j in seq_len(m)) {
    d <- d + weights[j] * dbeta(u, j, m - j + 1)
  }
  d
}


#' Bernstein polynomial CDF on [0,1]
#'
#' @param u values in (0,1)
#' @param weights probability weights (length m, sum to 1)
#' @param m Bernstein degree
#' @return CDF values
#' @noRd
.bernstein_cdf <- function(u, weights, m) {
  p <- numeric(length(u))
  for (j in seq_len(m)) {
    p <- p + weights[j] * pbeta(u, j, m - j + 1)
  }
  p
}


#' Full density of Bernstein-EGPD model
#'
#' @param x data values (>= 0)
#' @param sigma GPD scale
#' @param xi GPD shape
#' @param kappa power transform parameter
#' @param weights Bernstein weights
#' @param m Bernstein degree
#' @return density values
#' @noRd
.bernstein_full_density <- function(x, sigma, xi, kappa, weights, m) {
  u <- .pgpd_std(x, sigma, xi)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  uk <- u^kappa
  uk <- pmin(pmax(uk, 1e-10), 1 - 1e-10)
  .bernstein_density(uk, weights, m) * kappa * u^(kappa - 1) *
    .dgpd_std(x, sigma, xi)
}


#' Full CDF of Bernstein-EGPD model
#'
#' @param q quantile values
#' @param sigma GPD scale
#' @param xi GPD shape
#' @param kappa power transform parameter
#' @param weights Bernstein weights
#' @param m Bernstein degree
#' @return CDF values
#' @noRd
.bernstein_full_cdf <- function(q, sigma, xi, kappa, weights, m) {
  u <- .pgpd_std(q, sigma, xi)
  u <- pmin(pmax(u, 0), 1)
  uk <- u^kappa
  .bernstein_cdf(uk, weights, m)
}


#' Full quantile of Bernstein-EGPD model (numerical inversion)
#'
#' @param p probability values
#' @param sigma GPD scale
#' @param xi GPD shape
#' @param kappa power transform parameter
#' @param weights Bernstein weights
#' @param m Bernstein degree
#' @return quantile values
#' @noRd
.bernstein_full_quantile <- function(p, sigma, xi, kappa, weights, m) {
  ## Numerically invert the Bernstein CDF on [0,1] first
  ## B(u^kappa) = p => find u^kappa such that B(v) = p => v = B^{-1}(p)
  ## Then u = v^{1/kappa}, x = Q_GPD(u)
  sapply(p, function(pp) {
    if (pp <= 0) return(0)
    if (pp >= 1) return(Inf)
    ## Find v such that bernstein_cdf(v) = pp
    v <- tryCatch(
      uniroot(function(v) .bernstein_cdf(v, weights, m) - pp,
              interval = c(0, 1), tol = 1e-8)$root,
      error = function(e) pp  # fallback
    )
    u <- v^(1 / kappa)
    .qgpd_std(u, sigma, xi)
  })
}


#' Negative log-likelihood for Bernstein-EGPD model
#'
#' @param theta unconstrained parameter vector:
#'   c(log(kappa), log(sigma), xi, alpha_1, ..., alpha_m)
#' @param x data
#' @param m Bernstein degree
#' @param fix.arg named list of fixed parameters
#' @return scalar negative log-likelihood
#' @noRd
.bernstein_nll <- function(theta, x, m, fix.arg) {
  ## Parse theta
  kappa <- exp(theta[1])
  sigma <- exp(theta[2])
  xi    <- theta[3]
  alpha <- theta[4:(3 + m)]

  ## Override with fixed args
  if (!is.null(fix.arg)) {
    if ("kappa" %in% names(fix.arg)) kappa <- fix.arg[["kappa"]]
    if ("sigma" %in% names(fix.arg)) sigma <- fix.arg[["sigma"]]
    if ("xi" %in% names(fix.arg)) xi <- fix.arg[["xi"]]
  }

  ## Softmax weights
  alpha_shifted <- alpha - max(alpha)
  weights <- exp(alpha_shifted) / sum(exp(alpha_shifted))

  ## GPD quantities
  u <- .pgpd_std(x, sigma, xi)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)

  ## Power transform
  uk <- u^kappa
  uk <- pmin(pmax(uk, 1e-10), 1 - 1e-10)

  ## Bernstein density at u^kappa
  b_dens <- .bernstein_density(uk, weights, m)
  b_dens <- pmax(b_dens, .Machine$double.xmin)

  ## GPD log-density
  gpd_logd <- .dgpd_std(x, sigma, xi, log = TRUE)

  ## Full log-likelihood: log(b(u^kappa)) + log(kappa) + (kappa-1)*log(u) + log(h(x))
  ll <- sum(log(b_dens) + log(kappa) + (kappa - 1) * log(u) + gpd_logd)

  if (!is.finite(ll)) return(1e20)
  -ll
}


#' Bernstein polynomial EGPD fitting
#'
#' Two-stage procedure: (1) initial parametric MLE for GPD parameters,
#' (2) joint optimization of GPD + kappa + Bernstein weights.
#'
#' @param x data
#' @param type type for initial parametric fit
#' @param start optional starting values for sigma, xi, kappa
#' @param fix.arg fixed parameters
#' @param m Bernstein degree
#' @param optim.method optimization method
#' @param hessian compute SEs?
#' @param call matched call from fitegpd
#' @param ... extra args to optim
#' @return fitegpd object
#' @noRd
.fitegpd_bernstein <- function(x, type, start, fix.arg, m,
                                optim.method, hessian, call, ...) {
  n <- length(x)

  ## Stage 1: initial parametric MLE for GPD parameters
  init_fit <- tryCatch(
    fitegpd(x, type = type, family = "egpd", method = "mle",
            start = start, fix.arg = fix.arg, hessian = FALSE),
    error = function(e) NULL
  )

  if (!is.null(init_fit) && init_fit$convergence == 0) {
    sigma0 <- init_fit$estimate[["sigma"]]
    xi0    <- init_fit$estimate[["xi"]]
    kappa0 <- if ("kappa" %in% names(init_fit$estimate)) init_fit$estimate[["kappa"]] else 1
  } else {
    ## Fallback starting values
    xpos <- x[x > 0]
    if (length(xpos) < 2) xpos <- abs(x) + 0.01
    m_x <- mean(xpos); v_x <- var(xpos)
    xi0 <- max(min(0.5 * (m_x^2 / v_x - 1), 0.5), -0.5)
    sigma0 <- max(m_x * (1 - xi0), 0.01)
    kappa0 <- 1
  }

  ## Override with user start values
  if (!is.null(start)) {
    if ("sigma" %in% names(start)) sigma0 <- start[["sigma"]]
    if ("xi" %in% names(start)) xi0 <- start[["xi"]]
    if ("kappa" %in% names(start)) kappa0 <- start[["kappa"]]
  }

  ## Initial Bernstein weights: uniform (alpha = 0)
  alpha0 <- rep(0, m)

  ## Build theta0: c(log(kappa), log(sigma), xi, alpha_1..alpha_m)
  theta0 <- c(log(kappa0), log(sigma0), xi0, alpha0)

  ## Stage 2: joint optimization (higher-dimensional, needs more iterations)
  dots <- list(...)
  if (is.null(dots$control)) dots$control <- list()
  if (is.null(dots$control$maxit)) dots$control$maxit <- 10000
  opt <- do.call(optim, c(list(par = theta0, fn = .bernstein_nll, x = x, m = m,
                                fix.arg = fix.arg, method = optim.method,
                                hessian = hessian), dots))

  ## Extract estimates
  kappa_hat <- exp(opt$par[1])
  sigma_hat <- exp(opt$par[2])
  xi_hat    <- opt$par[3]
  alpha_hat <- opt$par[4:(3 + m)]

  ## Softmax weights
  alpha_shifted <- alpha_hat - max(alpha_hat)
  weights <- exp(alpha_shifted) / sum(exp(alpha_shifted))

  estimate <- c(sigma = sigma_hat, xi = xi_hat, kappa = kappa_hat)

  ## SE via delta method for the 3 main params
  se <- rep(NA_real_, 3)
  names(se) <- names(estimate)
  V_par <- NULL

  if (hessian && !is.null(opt$hessian)) {
    V_theta <- tryCatch(solve(opt$hessian), error = function(e) NULL)
    if (!is.null(V_theta)) {
      ## Jacobian for the 3 main params (kappa=log, sigma=log, xi=identity)
      jac3 <- c(kappa_hat, sigma_hat, 1)
      ## Extract 3x3 sub-block
      V3 <- V_theta[1:3, 1:3]
      J3 <- diag(jac3)
      ## Reorder: optim has (kappa, sigma, xi), we want (sigma, xi, kappa)
      reord <- c(2, 3, 1)
      J3 <- J3[reord, reord]
      V3 <- V3[reord, reord]
      V_par <- J3 %*% V3 %*% t(J3)
      rownames(V_par) <- colnames(V_par) <- names(estimate)
      se <- sqrt(pmax(diag(V_par), 0))
    }
  }

  ## Total number of free parameters: 3 (sigma, xi, kappa) + (m-1) Bernstein weights
  ## (softmax has m-1 degrees of freedom since weights sum to 1)
  npar <- 3 + (m - 1)
  ll <- -opt$value

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
    family    = "egpd",
    method    = "bernstein",
    fix.arg   = fix.arg,
    convergence = opt$convergence,
    optim     = opt,
    call      = call,
    bernstein.m       = m,
    bernstein.weights = weights
  ), class = "fitegpd")
}
