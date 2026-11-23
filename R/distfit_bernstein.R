## Bernstein polynomial fitting for continuous EGPD
##
## Replaces the parametric G-transformation with a flexible
## Bernstein polynomial density on [0,1], applied after a power
## transform u^kappa. The full density is:
##   f(x) = b(H(x)^kappa) * kappa * H(x)^(kappa-1) * h(x)
## where H is the GPD CDF, h is the GPD density, and b is the
## Bernstein polynomial density with softmax-parameterized weights.


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


#' Log-density of the Bernstein-EGPD model, with its gradient
#'
#' Evaluated in log space throughout, so that no value is clipped: log H from the
#' GPD log-survival log S as log(1 - exp(log S)), log v = kappa log H, log(1 - v)
#' as log(1 - exp(log v)), and the Bernstein density b(v) as the log-sum-exp of
#' log w_j + log dbeta(v; j, m - j + 1). For a large kappa, kappa S (and so
#' 1 - v) can be of order one however small S is. Where S and kappa S are both
#' below e^-40, log(1 - v) is log(kappa) + log S to relative accuracy S + kappa S
#' (kappa log H may underflow there); where S is subnormal or zero,
#' kappa log H = -kappa S is formed as -exp(log(kappa) + log S).
#'
#' @param x data values, positive and inside the GPD support
#' @param sigma,xi GPD scale and shape
#' @param kappa power transform parameter
#' @param lw log Bernstein weights (length m)
#' @param m Bernstein degree
#' @param deriv logical: also return the derivatives?
#' @return list: \code{ll}, the log-densities, and with \code{deriv = TRUE}
#'   \code{d}, the length(x)-by-(3 + m) matrix of their derivatives with respect to
#'   (log kappa, log sigma, xi, alpha_1, ..., alpha_m), alpha the softmax
#'   log-weights.
#' @noRd
.bernstein_logdens <- function(x, sigma, xi, kappa, lw, m, deriv = FALSE, ratios = FALSE) {
  n <- length(x)
  j <- seq_len(m)
  z <- x / sigma
  u <- xi * z
  lS <- .gpd_logsurv_inside(x, rep(sigma,n), rep(xi,n))
  L <- if (xi == 0) numeric(n) else -xi * lS        # log(1 + xi z)
  lH <- .log1mexp(lS)                              # GPD log-CDF
  tiny <- which(x > 0 & -lS < .Machine$double.xmin)
  if (length(tiny))
    lH[tiny] <- .gpd_logu_small(x[tiny], rep(sigma,length(tiny)), rep(xi,length(tiny)))
  lkS <- log(kappa) + lS                           # log(kappa S)
  lv <- kappa * lH                                 # log v
  lHk <- (kappa - 1) * lH
  ## a subnormal or zero S has lost the digits of log H = -S, while
  ## kappa log H = -kappa S need not be small
  sub <- which(lS < log(.Machine$double.xmin))
  if (length(sub)) {
    lv[sub] <- -exp(lkS[sub])
    lHk[sub] <- lv[sub] - lH[sub]
  }
  l1v <- ifelse(lS < -40 & lkS < -40, lkS, .log1mexp(lv))   # log(1 - v)
  ## log(w_j dbeta(v; j, m - j + 1)); the exponents j - 1 and m - j vanish at the
  ## ends, so those products are left out rather than formed as 0 * log(0)
  M <- matrix((if (ratios) 0 else lw) + log(m) + lchoose(m - 1, j - 1), n, m, byrow = TRUE)
  if (m > 1L) {
    M[, -1L] <- M[, -1L] + outer(lv, j[-1L] - 1)
    M[, -m] <- M[, -m] + outer(l1v, m - j[-m])
  }
  if (ratios) {
    basis <- M
    M <- sweep(M, 2L, lw, "+")
  }
  mx <- M[cbind(seq_len(n), max.col(M, ties.method = "first"))]
  lb <- mx + log(rowSums(exp(M - mx)))
  if (ratios) return(basis - lb)
  ll <- lb + log(kappa) + lHk - log(sigma) + lS - L
  if (!deriv) return(list(ll = ll))
  R <- exp(M - lb)                                 # posterior weights of the j
  ## A = v b'(v) / b(v) = sum_j r_j [(j - 1) - (m - j) v / (1 - v)]
  A <- drop(R %*% (j - 1))
  if (m > 1L)
    A <- A - drop(exp(M[, -m, drop = FALSE] - lb + lv - l1v) %*% (m - j[-m]))
  gap <- 1 + u
  edge <- which(u < -.5)
  if (length(edge))
    gap[edge] <- gpd_endpoint_gap_cpp(x[edge], rep(sigma,length(edge)), rep(xi,length(edge)))
  zt <- z / gap
  ## g = d log S / d xi = log(1 + u) / xi^2 - z / (xi (1 + u)); its series in u
  ## for small |u| avoids the cancellation of the two terms
  g <- if (xi == 0) z^2 / 2 else
    ifelse(abs(u) < 1e-3,
           z^2 * (1/2 + u * (-2/3 + u * (3/4 + u * (-4/5 + u * (5/6 + u * (-6/7 + u * 7/8)))))),
           (L / xi - zt) / xi)
  B <- kappa * A + kappa - 1                       # H d ll / d H
  BSH <- B * exp(lS - lH)                          # B S / H
  if (length(sub))                                 # (A + 1 - 1/kappa) kappa S / H
    BSH[sub] <- (A[sub] + 1 - 1 / kappa) * exp(lkS[sub] - lH[sub])
  dscale <- -BSH * zt - 1 + (1 + xi) * zt
  dshape <- (1 - BSH) * g - zt
  if (length(tiny)) {
    ## S/H alone can overflow although (S/H)*zt and (S/H)*g are finite.
    lz <- log(x[tiny]) - log(sigma)
    t <- sign(xi)*exp(log(abs(xi)) + lz)
    a <- abs(t) < 1e-3
    r <- numeric(length(t))
    r[a] <- .5 + t[a]*(-2/3 + t[a]*(3/4 - 4/5*t[a]))
    r[!a] <- (L[tiny[!a]] - t[!a]/gap[tiny[!a]]) / t[!a]^2
    lratio <- lS[tiny] - lH[tiny]
    dscale[tiny] <- -B[tiny]*exp(lratio + lz - L[tiny]) - 1 + (1 + xi)*zt[tiny]
    dshape[tiny] <- g[tiny] - zt[tiny] -
      sign(B[tiny])*exp(log(abs(B[tiny])) + lratio + 2*lz + log(r))
  }
  d <- cbind(lv * (A + 1) + 1, dscale, dshape,
             R - matrix(exp(lw), n, m, byrow = TRUE), deparse.level=0)
  list(ll = ll, d = d)
}


#' Full density of Bernstein-EGPD model
#'
#' @param x data values (>= 0)
#' @param sigma GPD scale
#' @param xi GPD shape
#' @param kappa power transform parameter
#' @param weights Bernstein weights
#' @param m Bernstein degree
#' @return density values; at x = 0 the limit (infinite for kappa < 1,
#'   m w_1 / sigma for kappa = 1, 0 for kappa > 1)
#' @noRd
.bernstein_full_density <- function(x, sigma, xi, kappa, weights, m) {
  out <- numeric(length(x))
  out[is.na(x)] <- x[is.na(x)]
  ok <- is.finite(x) & x > 0 & is.finite(.gpd_cdf_surv_std(x, sigma, xi)$log_surv)
  if (any(ok))
    out[ok] <- exp(.bernstein_logdens(x[ok], sigma, xi, kappa, log(weights), m)$ll)
  zero <- which(x == 0)
  if (length(zero))
    out[zero] <- if (kappa < 1) Inf else if (kappa == 1) m * weights[1] / sigma else 0
  out
}


#' Full CDF of Bernstein-EGPD model
#'
#' B(v) with v = H^kappa. The power u^kappa of the GPD CDF u = H has relative
#' error about kappa * 1e-16 (near u = 1, u is 1 - S rounded), so for kappa > 1000
#' v is formed as exp(kappa log H), with log H from the GPD log-survival; v is
#' then right however small S is, e.g. exp(-kappa S) where u rounds to 1.
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
  uk <- pegpd(q, sigma=sigma, xi=xi, kappa=kappa, type=1)
  .bernstein_cdf(uk, weights, m)
}


#' Full quantile of Bernstein-EGPD model (numerical inversion)
#'
#' Solves B(v) = p for v = H(x)^kappa, with B the Bernstein CDF, then
#' x = Q_GPD(v^(1/kappa)). The root is found on the log scale of v for p <= 1/2
#' and of t = 1 - v for p > 1/2, using the exact complement
#' 1 - B(1 - t) = sum_j w_j pbeta(t, m - j + 1, j); otherwise the upper tail
#' rounds to v = 1 and the quantile to Inf (or the endpoint) long before p = 1.
#' The GPD quantile is taken from log(1 - u), with u = v^(1/kappa), by the
#' log1mexp switch: log1p(-exp(log u)) for u < 1/2, log(-expm1(log u)) otherwise,
#' so that both a tiny u (lower tail, x ~ sigma u) and u near 1 (upper tail) keep
#' their relative precision.
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
  j <- seq_len(m)
  lw <- log(weights[j])
  lse <- function(a) { mx <- max(a); if (!is.finite(mx)) mx else mx + log(sum(exp(a - mx))) }
  logB <- function(lv, a) {
    ## Beta(a,m-a+1) is a binomial tail. Below the normal range its
    ## first term gives the log CDF to rounding without materializing v.
    lp <- if (lv < log(.Machine$double.xmin)) lchoose(m,a) + a*lv else
      pbeta(exp(lv), a, m-a+1, log.p=TRUE)
    lse(lw + lp)
  }
  upper <- if (xi < 0) -sigma / xi else Inf
  vapply(p, function(pp) {
    if (is.na(pp)) return(pp)
    if (pp < 0 || pp > 1) return(NaN)
    if (pp == 0) return(0)
    if (pp == 1) return(upper)
    lower <- pp <= .5
    target <- if (lower) log(pp) else log1p(-pp)
    a <- if (lower) j else m-j+1L
    ## Every component CDF is at most m*v, so this brackets either tail.
    lv <- uniroot(function(v) logB(v,a)-target, c(target-log(m)-1,0), tol=1e-13)$root
    qegpd(lv, sigma=sigma, xi=xi, kappa=kappa, type=1, lower.tail=lower, log.p=TRUE)
  }, numeric(1))
}


#' Negative log-likelihood for Bernstein-EGPD model, or its gradient
#'
#' @param theta unconstrained parameter vector:
#'   c(log(kappa), log(sigma), xi, alpha_1, ..., alpha_m)
#' @param x data (positive)
#' @param m Bernstein degree
#' @param fix.arg named list of fixed parameters
#' @param deriv logical: return the gradient with respect to theta instead?
#' @return scalar negative log-likelihood, 1e20 where the likelihood is zero or
#'   not finite (the gradient is zero there)
#' @noRd
.bernstein_nll <- function(theta, x, m, fix.arg, deriv = FALSE) {
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
  infeasible <- if (deriv) numeric(length(theta)) else 1e20
  if (!all(is.finite(alpha)) || !is.finite(kappa) || kappa <= 0 ||
      !is.finite(sigma) || sigma <= 0 || !is.finite(xi))
    return(infeasible)
  if (xi < 0) {
    edge <- x[x >= -sigma/xi]
    if (length(edge) && any(gpd_endpoint_gap_cpp(edge, rep(sigma,length(edge)), rep(xi,length(edge))) <= 0))
      return(infeasible)
  }

  ## Softmax log-weights
  alpha_shifted <- alpha - max(alpha)
  lw <- alpha_shifted - log(sum(exp(alpha_shifted)))

  ## log(b(H^kappa)) + log(kappa) + (kappa - 1) log(H) + log(h), in log space
  r <- .bernstein_logdens(x, sigma, xi, kappa, lw, m, deriv = deriv)
  if (deriv) {
    g <- -colSums(r$d)
    return(if (all(is.finite(g))) g else infeasible)
  }
  ll <- sum(r$ll)
  if (!is.finite(ll)) return(1e20)
  -ll
}


#' Bernstein polynomial EGPD fitting
#'
#' Two-stage procedure: (1) initial parametric MLE for GPD parameters,
#' (2) joint optimization of GPD + kappa + Bernstein weights, which for
#' Nelder-Mead is polished by alternating BFGS (exact gradient) and Nelder-Mead
#' restarts until a cycle no longer improves it.
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
.distfit_bernstein <- function(x, type, start, fix.arg, m,
                                optim.method, hessian, call, ...) {
  n <- length(x)
  if (length(m) != 1L || !is.finite(m) || m < 1 || m != as.integer(m))
    stop("Bernstein degree must be a positive integer.", call.=FALSE)
  if (length(fix.arg)) {
    if (is.null(names(fix.arg)) || anyDuplicated(names(fix.arg)) ||
        any(!names(fix.arg) %in% c("sigma", "xi", "kappa")))
      stop("Fixed Bernstein parameters must be uniquely named sigma, xi or kappa.",
           call.=FALSE)
    for (nm in names(fix.arg)) {
      v <- fix.arg[[nm]]
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
          (nm != "xi" && v <= 0))
        stop("Invalid fixed Bernstein parameter: ", nm, call.=FALSE)
    }
  }
  if (length(start)) {
    if (is.null(names(start)) || anyDuplicated(names(start)) ||
        any(!names(start) %in% c("sigma", "xi", "kappa")))
      stop("Bernstein start values must be uniquely named sigma, xi or kappa.",
           call.=FALSE)
    for (nm in names(start)) {
      v <- start[[nm]]
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
          (nm != "xi" && v <= 0))
        stop("Invalid Bernstein start value: ", nm, call.=FALSE)
    }
  }

  main_names <- c("kappa", "sigma", "xi")
  free_names <- c(setdiff(main_names, names(fix.arg)), rep("weight", m - 1L))
  fixed <- as.list(fix.arg)
  if (!is.null(fixed$xi))
    .distfit_check_endpoint(max(x), fixed$xi, fixed$sigma, free_names, list(...), monotone = m == 1L)

  ## Stage 1: initial parametric MLE for GPD parameters
  ## only a source of starting values (its convergence is checked below), so its
  ## non-convergence warning would be noise
  init_fit <- tryCatch(
    suppressWarnings(distfit(x, type = type, family = "egpd", method = "mle",
                             start = start, fix.arg = fix.arg, hessian = FALSE)),
    error = function(e) {
      warning("Parametric initialization failed: ", conditionMessage(e),
              "; using moment starting values.", call.=FALSE)
      NULL
    }
  )

  if (!is.null(init_fit) && init_fit$convergence == 0) {
    initial <- c(init_fit$estimate, unlist(fix.arg))
    sigma0 <- initial[["sigma"]]
    xi0    <- initial[["xi"]]
    kappa0 <- if ("kappa" %in% names(initial)) initial[["kappa"]] else 1
  } else {
    ## Fallback starting values
    xpos <- x[x > 0]
    if (length(xpos) < 2) xpos <- abs(x) + 0.01
    m_x <- mean(xpos); v_x <- var(xpos)
    xi0 <- max(min(0.5 * (1 - m_x^2 / v_x), 0.5), -0.5)
    sigma0 <- max(m_x * (1 - xi0), 0.01)
    kappa0 <- 1
  }

  ## Override with user start values
  if (!is.null(start)) {
    if ("sigma" %in% names(start)) sigma0 <- start[["sigma"]]
    if ("xi" %in% names(start)) xi0 <- start[["xi"]]
    if ("kappa" %in% names(start)) kappa0 <- start[["kappa"]]
  }
  if ("sigma" %in% names(fix.arg)) sigma0 <- fix.arg[["sigma"]]
  if ("xi" %in% names(fix.arg)) xi0 <- fix.arg[["xi"]]
  if ("kappa" %in% names(fix.arg)) kappa0 <- fix.arg[["kappa"]]

  ## Initial Bernstein weights: uniform (alpha = 0)
  alpha0 <- rep(0, m)

  ## Build theta0: c(log(kappa), log(sigma), xi, alpha_1..alpha_m)
  theta0 <- c(log(kappa0), log(sigma0), xi0, alpha0)
  names(theta0) <- c("logkappa", "logsigma", "xi", paste0("alpha", seq_len(m)))
  fixed_idx <- which(main_names %in% names(fix.arg))
  ## One reference log-weight removes the exact softmax shift invariance.
  free_idx <- setdiff(seq_along(theta0), c(fixed_idx, 3L + m))
  if (!length(free_idx))
    stop("No free Bernstein parameters to estimate.", call.=FALSE)
  nll <- function(theta) {
    full <- theta0
    full[free_idx] <- theta
    .bernstein_nll(full, x, m, fix.arg)
  }
  gr <- function(theta) {
    full <- theta0
    full[free_idx] <- theta
    .bernstein_nll(full, x, m, fix.arg, deriv = TRUE)[free_idx]
  }

  ## Stage 2: joint optimization (higher-dimensional, needs more iterations).
  ## As in distfit(), a single free parameter is fitted by BFGS, not Nelder-Mead.
  dots <- list(...)
  optim.method <- .distfit_method(optim.method, length(free_idx), dots)
  if (is.null(dots$control)) dots$control <- list()
  if (is.null(dots$control$maxit)) dots$control$maxit <- 10000
  opt <- do.call(optim, c(list(par = theta0[free_idx], fn = nll,
                                gr = if (optim.method %in% c("BFGS", "CG", "L-BFGS-B")) gr,
                                method = optim.method), dots))
  ## Nelder-Mead's simplex can collapse on the flat ridges of this objective and
  ## report convergence well short of the maximum (by more than one log-likelihood
  ## unit). Its solution is therefore refined by alternating BFGS, with the exact
  ## gradient, and fresh Nelder-Mead simplices until a whole cycle gains nothing.
  if (optim.method == "Nelder-Mead" && is.finite(opt$value) && opt$value < 1e20) {
    better <- function(o, cand)
      if (!is.null(cand) && is.finite(cand$value) && cand$value < o$value) cand else o
    ## A tiny softmax gradient can hide a profitable component with almost
    ## zero weight. At fixed main parameters the mixture likelihood is concave
    ## in the weights: check the best vertex direction in those coordinates.
    weight_step <- function(o) {
      if (m == 1L) return(o)
      full <- theta0
      full[free_idx] <- o$par
      wi <- 3L + seq_len(m)
      lw <- full[wi] - max(full[wi])
      lw <- lw - log(sum(exp(lw)))
      lr <- .bernstein_logdens(x, exp(full[2]), full[3], exp(full[1]), lw, m, ratios = TRUE)
      score <- colSums(expm1(lr))
      j <- which.max(score)
      if (score[j] <= 0) return(o)
      gain <- function(t) sum(.logspace_add(log1p(-t), log(t) + lr[, j]))
      step <- stats::optimize(gain, c(0, 1), maximum = TRUE, tol = 1e-12)
      if (step$objective <= 1e-10 * n) return(o)
      t <- step$maximum
      lw <- lw + log1p(-t)
      lw[j] <- .logspace_add(lw[j], log(t))
      full[wi] <- lw - lw[m]
      par <- full[free_idx]
      value <- nll(par)
      if (value < o$value) { o$par <- par; o$value <- value }
      o
    }
    settled <- FALSE
    for (cycle in seq_len(20L)) {
      v0 <- opt$value
      opt <- better(opt, tryCatch(optim(opt$par, nll, gr, method = "BFGS",
                                        control = list(maxit = 1000L, reltol = 1e-12)),
                                  error = function(e) NULL))
      opt <- better(opt, optim(opt$par, nll, method = "Nelder-Mead",
                               control = dots$control))
      if (v0 - opt$value <= 1e-9 * (abs(opt$value) + 1)) {
        candidate <- weight_step(opt)
        if (candidate$value == opt$value) { settled <- TRUE; break }
        opt <- candidate
      }
    }
    ## converged: a whole cycle left the solution unchanged
    opt$convergence <- if (settled) 0L else 1L
  }
  opt <- .distfit_check_box(opt, optim.method, dots, hessian, nll)
  .distfit_check_optim(opt, optim.method)

  ## Extract estimates
  full <- theta0
  full[free_idx] <- opt$par
  kappa_hat <- exp(full[1])
  sigma_hat <- exp(full[2])
  xi_hat    <- full[3]
  alpha_hat <- full[4:(3 + m)]
  if (is.null(fixed$xi))
    .distfit_check_endpoint(max(x), xi_hat, sigma_hat, free_names, dots, monotone = m == 1L)
  if (hessian && !isTRUE(opt$bound.active))
    opt$hessian <- optimHess(opt$par, nll, gr, control = dots$control)

  ## Softmax weights
  alpha_shifted <- alpha_hat - max(alpha_hat)
  weights <- exp(alpha_shifted) / sum(exp(alpha_shifted))

  estimates <- c(sigma = unname(sigma_hat), xi = unname(xi_hat),
                 kappa = unname(kappa_hat))
  free_main <- setdiff(names(estimates), names(fix.arg))
  estimate <- estimates[free_main]

  ## SE via delta method for the 3 main params
  se <- rep(NA_real_, length(estimate))
  names(se) <- names(estimate)
  V_par <- NULL

  if (hessian && !is.null(opt$hessian)) {
    R <- tryCatch(chol(opt$hessian), error = function(e) NULL)
    if (!is.null(R)) {
      V_theta <- chol2inv(R)
      pos <- match(match(free_main, main_names), free_idx)
      jac <- ifelse(free_main == "xi", 1, estimate)
      J <- diag(jac, length(jac))
      V_par <- J %*% V_theta[pos, pos, drop=FALSE] %*% J
      rownames(V_par) <- colnames(V_par) <- names(estimate)
      se <- sqrt(diag(V_par))
    } else {
      warning("Bernstein Hessian is not positive definite; standard errors unavailable.",
              call.=FALSE)
    }
  }

  ## Total number of free parameters: 3 (sigma, xi, kappa) + (m-1) Bernstein weights
  ## (softmax has m-1 degrees of freedom since weights sum to 1)
  npar <- length(free_idx)
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
  ), class = "distfit")
}
