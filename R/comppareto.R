## Composite Pareto distribution utilities and GAM family hooks

.comppareto_specs <- c("lnorm", "gamma", "weibull", "exp")

.comppareto_exp_link <- function(x) exp(x)
attr(.comppareto_exp_link, "deriv") <- function(x) exp(x)

.comppareto_identity_link <- function(x) x
attr(.comppareto_identity_link, "deriv") <- function(x) rep(1, length(x))

.comppareto_body_info <- function(spec) {
  spec <- match.arg(spec, .comppareto_specs)

  switch(
    spec,
    lnorm = list(
      spec = spec,
      gam_names = c("meanlog", "logsdlog", "logalpha", "logtheta"),
      response_names = c("meanlog", "sdlog", "alpha", "theta"),
      body_names = c("meanlog", "sdlog"),
      d = stats::dlnorm,
      p = stats::plnorm,
      q = stats::qlnorm,
      unlink = list(
        .comppareto_identity_link,
        .comppareto_exp_link,
        .comppareto_exp_link,
        .comppareto_exp_link
      )
    ),
    gamma = list(
      spec = spec,
      gam_names = c("logshape", "logscale", "logalpha", "logtheta"),
      response_names = c("shape", "scale", "alpha", "theta"),
      body_names = c("shape", "scale"),
      d = stats::dgamma,
      p = stats::pgamma,
      q = stats::qgamma,
      unlink = rep(list(.comppareto_exp_link), 4L)
    ),
    weibull = list(
      spec = spec,
      gam_names = c("logshape", "logscale", "logalpha", "logtheta"),
      response_names = c("shape", "scale", "alpha", "theta"),
      body_names = c("shape", "scale"),
      d = stats::dweibull,
      p = stats::pweibull,
      q = stats::qweibull,
      unlink = rep(list(.comppareto_exp_link), 4L)
    ),
    exp = list(
      spec = spec,
      gam_names = c("lograte", "logalpha", "logtheta"),
      response_names = c("rate", "alpha", "theta"),
      body_names = c("rate"),
      d = stats::dexp,
      p = stats::pexp,
      q = stats::qexp,
      unlink = rep(list(.comppareto_exp_link), 3L)
    )
  )
}

.comppareto_named_pars <- function(spec, pars) {
  info <- .comppareto_body_info(spec)
  if (length(pars) != length(info$response_names)) {
    stop("Incorrect number of parameters for CompPareto specification '", spec, "'.")
  }
  names(pars) <- info$response_names
  pars
}

.comppareto_body_args <- function(spec, pars) {
  info <- .comppareto_body_info(spec)
  body <- pars[info$body_names]
  names(body) <- info$body_names
  body
}

.comppareto_recycle_args <- function(n, args) {
  lapply(args, rep_len, length.out = n)
}

.comppareto_subset_args <- function(args, idx) {
  lapply(args, function(x) x[idx])
}

.comppareto_pareto_density <- function(x, alpha, theta) {
  out <- alpha * theta^alpha / (x + theta)^(alpha + 1)
  out[!is.finite(x) | x <= 0] <- 0
  out
}

.comppareto_pareto_cdf <- function(q, alpha, theta) {
  out <- 1 - (theta / (q + theta))^alpha
  out[!is.finite(q)] <- ifelse(q[!is.finite(q)] > 0, 1, NA_real_)
  out[q <= 0] <- 0
  out
}

.comppareto_pareto_quantile <- function(p, alpha, theta) {
  out <- rep(NA_real_, length(p))
  out[p == 1] <- Inf
  idx <- is.finite(p) & p >= 0 & p < 1
  out[idx] <- theta[idx] * ((1 - p[idx])^(-1 / alpha[idx]) - 1)
  out
}

.comppareto_phi <- function(spec, alpha, theta, body_args) {
  info <- .comppareto_body_info(spec)
  G_theta <- do.call(info$p, c(list(theta), body_args))
  g_theta <- do.call(info$d, c(list(theta), body_args))
  G2_theta <- .comppareto_pareto_cdf(theta, alpha = alpha, theta = theta)
  g2_theta <- .comppareto_pareto_density(theta, alpha = alpha, theta = theta)
  tail_surv_theta <- 1 - G2_theta
  phi <- g_theta * tail_surv_theta / (g2_theta * G_theta)
  phi[
    !is.finite(phi) | phi <= 0 |
      !is.finite(G_theta) | G_theta <= 0 |
      !is.finite(g_theta) | g_theta <= 0 |
      !is.finite(tail_surv_theta) | tail_surv_theta <= 0
  ] <- NA_real_
  list(phi = phi, G_theta = G_theta, G2_theta = G2_theta, tail_surv_theta = tail_surv_theta)
}

.comppareto_fd_step <- function(eta) {
  ## Note: only used for body specs without analytic derivatives.
  ## h=1e-4 is unstable because observations near the splice point theta
  ## flip body/tail under +/- h perturbation; h=1e-3 is empirically robust.
  1e-3 * pmax(1, abs(eta))
}

## --------------------------------------------------------------------------
## Analytic gradient/Hessian helpers for the CompPareto likelihood.
##
## Notation:
##   y        observed response
##   theta    splice point (= Pareto scale)
##   alpha    Pareto shape
##   body     parameters of the body distribution g(.; body), G(.; body)
##   phi      = 2 * theta * g(theta; body) / (alpha * G(theta; body))
##   T        = log(phi) - log(1 + phi)
##   log f(y) = s_part(y; ...) + T(...)
## with
##   s_body(y) = log g(y; body) - log g(theta; body) + log alpha - log theta - log 2
##   s_tail(y) = log alpha + alpha * log(2 * theta) - (alpha + 1) * log(y + theta)
##
## Each spec helper returns, for a vector of observations,
##   logf_y         log g(y; body)                                       (n)
##   logg_theta     log g(theta; body)                                   (n)
##   logG_theta     log G(theta; body)                                   (n)
## together with their gradients and Hessians with respect to the linear
## predictors (eta_b, eta_theta).  All terms are stored using the
## convention used elsewhere in egpd: linear-predictor (eta) derivatives,
## upper-triangle Hessian columns laid out as (1,1), (1,2), ..., (p,p).
##
## .comppareto_assemble_d120() then converts those pieces into per-row
## gradient/Hessian contributions of the negative log-likelihood with
## respect to (eta_b, eta_alpha, eta_theta).
## --------------------------------------------------------------------------

.comppareto_m_z <- function(z) {
  ## m(z) = z / (1 - exp(-z)); m'(z) = ((1 - e^{-z}) - z e^{-z}) / (1 - e^{-z})^2
  ## Use a Taylor expansion for small |z| to avoid cancellation.
  m  <- numeric(length(z))
  mp <- numeric(length(z))
  small <- abs(z) < 1e-3
  if (any(small)) {
    zs <- z[small]
    m[small]  <- 1 + zs / 2 + zs^2 / 12 - zs^4 / 720
    mp[small] <- 0.5 + zs / 6 - zs^3 / 180
  }
  if (any(!small)) {
    zb <- z[!small]
    one_emz <- -expm1(-zb)        # 1 - exp(-z)
    emz <- exp(-zb)
    m[!small]  <- zb / one_emz
    mp[!small] <- (one_emz - zb * emz) / one_emz^2
  }
  list(m = m, mp = mp)
}

.comppareto_log1mexp_neg <- function(z) {
  ## log(1 - exp(-z)) for z > 0, numerically stable.
  out <- numeric(length(z))
  big <- z > log(2)
  out[big] <- log1p(-exp(-z[big]))
  out[!big] <- log(-expm1(-z[!big]))
  out
}

## Spec-specific blocks: each returns a list with the body quantities and
## their (eta_b, eta_theta) derivatives.  Body params are passed in their
## linear-predictor (eta) parameterization to make chain rules trivial.

.comppareto_body_block_exp <- function(y, eta_b, eta_theta) {
  eta_lambda <- eta_b[, 1L]
  lambda <- exp(eta_lambda)
  theta  <- exp(eta_theta)
  z      <- lambda * theta            # = exp(eta_lambda + eta_theta)

  logf_y <- eta_lambda - lambda * y
  dlogf_y_b <- matrix(1 - lambda * y, ncol = 1L)
  d2logf_y_b <- matrix(-lambda * y, ncol = 1L)   # only entry of upper triangle

  ## A = log g(theta; lambda) = eta_lambda - z
  logg_theta <- eta_lambda - z
  ## gradient wrt (eta_lambda, eta_theta): (1 - z, -z)
  dlogg_theta <- cbind(1 - z, -z)
  ## Hessian columns (b,b), (b,theta), (theta,theta)
  d2logg_theta <- cbind(-z, -z, -z)

  ## B = log G(theta; lambda) = log(1 - exp(-z))
  q_z <- .comppareto_log1mexp_neg(z)
  logG_theta <- q_z
  ## q'(z) = 1 / (e^z - 1)
  qprime <- 1 / expm1(z)
  ## q''(z) = -e^z / (e^z - 1)^2 = -(qprime + qprime^2) ??? simpler: derivative of 1/(e^z-1) wrt z = -e^z/(e^z-1)^2
  qpp    <- -exp(z) / expm1(z)^2
  ## chain rule: dz/d eta_lambda = dz/d eta_theta = z; d^2 z/(d eta_*)^2 = z; cross d^2 = z
  dlogG_theta <- cbind(qprime * z, qprime * z)
  d_bb_zz <- qpp * z^2 + qprime * z
  d2logG_theta <- cbind(d_bb_zz, d_bb_zz, d_bb_zz)

  ## u = log phi = log 2 + eta_theta + eta_lambda - z - eta_alpha - log(1 - e^{-z})
  ## but phi assembly is generic; expose via direct formula too.
  list(
    p_b = 1L,
    body_link_factor = matrix(lambda, ncol = 1L), # not used yet; kept for symmetry
    logf_y       = logf_y,
    dlogf_y_b    = dlogf_y_b,
    d2logf_y_b   = d2logf_y_b,
    logg_theta   = logg_theta,
    dlogg_theta  = dlogg_theta,
    d2logg_theta = d2logg_theta,
    logG_theta   = logG_theta,
    dlogG_theta  = dlogG_theta,
    d2logG_theta = d2logG_theta
  )
}

.comppareto_body_block <- function(spec, y, eta_b, eta_theta) {
  switch(
    spec,
    exp = .comppareto_body_block_exp(y, eta_b, eta_theta),
    stop("Analytic derivatives for spec '", spec, "' not yet implemented.")
  )
}

.comppareto_has_analytic <- function(spec) spec %in% c("exp")

.comppareto_assemble_d120 <- function(y, eta, spec) {
  n <- nrow(eta)
  blk <- .comppareto_body_block(spec, y, eta[, seq_len(ncol(eta) - 2L), drop = FALSE], eta[, ncol(eta) - 1L])
  ## Wait: eta column layout: [body_etas..., eta_alpha, eta_theta]
  ## Re-extract correctly.
  p_b <- blk$p_b
  eta_b      <- eta[, seq_len(p_b), drop = FALSE]
  eta_alpha  <- eta[, p_b + 1L]
  eta_theta  <- eta[, p_b + 2L]
  ## Recompute block with correct eta_theta in case the call above used wrong columns.
  blk <- .comppareto_body_block(spec, y, eta_b, eta_theta)

  alpha <- exp(eta_alpha)
  theta <- exp(eta_theta)

  P <- p_b + 2L
  ## Indices of (eta_b..., eta_alpha, eta_theta)
  idx_b     <- seq_len(p_b)
  idx_alpha <- p_b + 1L
  idx_theta <- p_b + 2L

  is_body <- y <= theta

  ## ---- u = log phi and its derivatives wrt (eta_b, eta_alpha, eta_theta) ----
  ## u = log 2 + eta_theta + log g(theta; body) - eta_alpha - log G(theta; body)
  u_val <- log(2) + eta_theta + blk$logg_theta - eta_alpha - blk$logG_theta
  phi   <- exp(u_val)
  inv1p <- 1 / (1 + phi)

  ## Gradient of u wrt (eta_b, eta_alpha, eta_theta)
  u_grad <- matrix(0, n, P)
  ## body-eta partials: come only from logg_theta - logG_theta
  if (p_b > 0L) {
    u_grad[, idx_b] <- blk$dlogg_theta[, idx_b, drop = FALSE] - blk$dlogG_theta[, idx_b, drop = FALSE]
  }
  u_grad[, idx_alpha] <- -1
  ## eta_theta column: explicit "+ eta_theta" gives 1, plus theta-partials from logg/logG
  u_grad[, idx_theta] <- 1 + blk$dlogg_theta[, p_b + 1L] - blk$dlogG_theta[, p_b + 1L]

  ## Hessian columns of u in upper-triangle order (1,1),(1,2),...,(P,P)
  n_h <- P * (P + 1L) / 2L
  u_hess <- matrix(0, n, n_h)

  ## Map (i,j) i<=j -> column index
  hcol <- function(i, j) {
    if (i > j) { tmp <- i; i <- j; j <- tmp }
    ## column index of (i,j) in row-major upper triangle laid out as in d120
    ## ordering used in d120 above: pos 1 = (1,1); pos 2 = (1,2); ... (1,P); (2,2); (2,3); ...
    ## Compute: for row k from 1..i-1, the row contributes (P - k + 1) entries.
    base <- 0L
    if (i > 1L) {
      base <- sum((P:1)[seq_len(i - 1L)])
    }
    base + (j - i + 1L)
  }

  ## d2 logg(theta) and d2 logG(theta) cover the (eta_b, eta_theta) block only.
  ## Indexing within blk$d2logg_theta: columns are upper-triangle of (eta_b..., eta_theta) of size (p_b+1).
  pb1 <- p_b + 1L
  bcol <- function(i, j) {
    if (i > j) { tmp <- i; i <- j; j <- tmp }
    base <- 0L
    if (i > 1L) base <- sum((pb1:1)[seq_len(i - 1L)])
    base + (j - i + 1L)
  }

  ## Fill (b,b) entries
  if (p_b > 0L) {
    for (i in seq_len(p_b)) for (j in i:p_b) {
      val <- blk$d2logg_theta[, bcol(i, j)] - blk$d2logG_theta[, bcol(i, j)]
      u_hess[, hcol(i, j)] <- val
    }
    ## (b, theta)
    for (i in seq_len(p_b)) {
      val <- blk$d2logg_theta[, bcol(i, pb1)] - blk$d2logG_theta[, bcol(i, pb1)]
      u_hess[, hcol(i, idx_theta)] <- val
    }
  }
  ## (alpha, *) = 0 (already zeros)
  ## (theta, theta)
  u_hess[, hcol(idx_theta, idx_theta)] <-
    blk$d2logg_theta[, bcol(pb1, pb1)] - blk$d2logG_theta[, bcol(pb1, pb1)]

  ## ---- T = u - log(1 + e^u) derivatives ----
  T_grad <- u_grad * inv1p
  pf <- phi * inv1p^2
  T_hess <- matrix(0, n, n_h)
  for (i in seq_len(P)) for (j in i:P) {
    T_hess[, hcol(i, j)] <- u_hess[, hcol(i, j)] * inv1p - pf * u_grad[, i] * u_grad[, j]
  }

  ## ---- s_part derivatives ----
  s_grad <- matrix(0, n, P)
  s_hess <- matrix(0, n, n_h)

  ## Body case: s_body = log g(y;b) - log g(theta;b) + log alpha - log theta - log 2
  if (any(is_body)) {
    i_idx <- which(is_body)
    if (p_b > 0L) {
      s_grad[i_idx, idx_b] <- blk$dlogf_y_b[i_idx, , drop = FALSE] -
        blk$dlogg_theta[i_idx, idx_b, drop = FALSE]
      ## Body Hessian (b,b) of s_body
      for (i in seq_len(p_b)) for (j in i:p_b) {
        val <- blk$d2logf_y_b[i_idx, bcol(i, j)] - blk$d2logg_theta[i_idx, bcol(i, j)]
        s_hess[i_idx, hcol(i, j)] <- val
      }
      ## (b, theta) of s_body: only -log g(theta;b) contributes via cross derivative
      for (i in seq_len(p_b)) {
        val <- -blk$d2logg_theta[i_idx, bcol(i, pb1)]
        s_hess[i_idx, hcol(i, idx_theta)] <- val
      }
    }
    s_grad[i_idx, idx_alpha] <- 1
    ## eta_theta partial of s_body: -d log g(theta;b)/d eta_theta - 1
    s_grad[i_idx, idx_theta] <- -blk$dlogg_theta[i_idx, pb1] - 1
    s_hess[i_idx, hcol(idx_theta, idx_theta)] <- -blk$d2logg_theta[i_idx, bcol(pb1, pb1)]
  }

  ## Tail case: s_tail = eta_alpha + alpha * (log 2 + eta_theta) - (alpha+1) * log(y + theta)
  if (any(!is_body)) {
    i_idx <- which(!is_body)
    yt   <- y[i_idx]
    al   <- alpha[i_idx]
    th   <- theta[i_idx]
    log_ratio <- log(2 * th) - log(yt + th)
    s_grad[i_idx, idx_alpha] <- 1 + al * log_ratio
    s_grad[i_idx, idx_theta] <- (al * yt - th) / (yt + th)
    s_hess[i_idx, hcol(idx_alpha, idx_alpha)] <- al * log_ratio
    s_hess[i_idx, hcol(idx_alpha, idx_theta)] <- al * yt / (yt + th)
    s_hess[i_idx, hcol(idx_theta, idx_theta)] <- -(al + 1) * th * yt / (yt + th)^2
  }

  ## ---- assemble negative-log-likelihood derivatives ----
  ## ell = -log f = -(s_part + T)
  grad <- -(s_grad + T_grad)
  hess <- -(s_hess + T_hess)

  cbind(grad, hess)
}

.comppareto_eta_matrix <- function(pars, likdata) {
  split_pars <- split(pars, likdata$idpars)
  eta <- lapply(seq_along(split_pars), function(i) {
    val <- drop(likdata$X[[i]] %*% split_pars[[i]])
    if (likdata$duplicate == 1) {
      val <- val[likdata$dupid + 1L]
    }
    off <- likdata$offsets[[i]]
    if (length(off) > 0L) {
      val <- val + off
    }
    val
  })
  do.call(cbind, eta)
}

.comppareto_eta_to_response <- function(spec, eta) {
  info <- .comppareto_body_info(spec)
  pars <- lapply(seq_len(ncol(eta)), function(i) info$unlink[[i]](eta[, i]))
  names(pars) <- info$response_names
  pars
}

.comppareto_nll_obs <- function(y, eta, spec) {
  pars <- .comppareto_eta_to_response(spec, eta)
  dens <- do.call(
    dcomppareto,
    c(list(x = y, spec = spec, log = FALSE), pars)
  )
  dens[!is.finite(dens)] <- 0
  -log(pmax(dens, .Machine$double.xmin))
}

.comppareto_d120 <- function(y, eta, spec) {
  if (.comppareto_has_analytic(spec)) {
    return(.comppareto_assemble_d120(y, eta, spec))
  }
  n <- nrow(eta)
  p <- ncol(eta)
  f0 <- .comppareto_nll_obs(y, eta, spec)

  grad <- matrix(NA_real_, nrow = n, ncol = p)
  hess_terms <- vector("list", length = p * (p + 1L) / 2L)
  pos <- 1L

  for (k in seq_len(p)) {
    hk <- .comppareto_fd_step(eta[, k])

    eta_pk <- eta
    eta_mk <- eta
    eta_pk[, k] <- eta_pk[, k] + hk
    eta_mk[, k] <- eta_mk[, k] - hk

    f_pk <- .comppareto_nll_obs(y, eta_pk, spec)
    f_mk <- .comppareto_nll_obs(y, eta_mk, spec)

    grad[, k] <- (f_pk - f_mk) / (2 * hk)
    hess_terms[[pos]] <- (f_pk - 2 * f0 + f_mk) / (hk^2)
    pos <- pos + 1L

    if (k < p) {
      for (j in (k + 1L):p) {
        hj <- .comppareto_fd_step(eta[, j])

        eta_pp <- eta
        eta_pm <- eta
        eta_mp <- eta
        eta_mm <- eta

        eta_pp[, k] <- eta_pp[, k] + hk
        eta_pp[, j] <- eta_pp[, j] + hj
        eta_pm[, k] <- eta_pm[, k] + hk
        eta_pm[, j] <- eta_pm[, j] - hj
        eta_mp[, k] <- eta_mp[, k] - hk
        eta_mp[, j] <- eta_mp[, j] + hj
        eta_mm[, k] <- eta_mm[, k] - hk
        eta_mm[, j] <- eta_mm[, j] - hj

        f_pp <- .comppareto_nll_obs(y, eta_pp, spec)
        f_pm <- .comppareto_nll_obs(y, eta_pm, spec)
        f_mp <- .comppareto_nll_obs(y, eta_mp, spec)
        f_mm <- .comppareto_nll_obs(y, eta_mm, spec)

        hess_terms[[pos]] <- (f_pp - f_pm - f_mp + f_mm) / (4 * hk * hj)
        pos <- pos + 1L
      }
    }
  }

  cbind(grad, do.call(cbind, hess_terms))
}

.comppareto_initial_values <- function(spec, likdata) {
  y <- likdata$y[, 1]
  y <- y[is.finite(y) & y > 0]
  if (length(y) == 0L) {
    stop("CompPareto models require strictly positive response values.")
  }

  theta0 <- stats::quantile(y, probs = 0.8, names = FALSE)
  theta0 <- max(theta0, .Machine$double.eps)

  switch(
    spec,
    lnorm = {
      logy <- log(y)
      c(
        mean(logy),
        log(max(stats::sd(logy), 0.1)),
        0,
        log(theta0)
      )
    },
    gamma = {
      y_mean <- mean(y)
      y_var <- stats::var(y)
      shape0 <- if (is.finite(y_var) && y_var > 0) y_mean^2 / y_var else 1
      scale0 <- if (is.finite(y_mean) && y_mean > 0) y_var / y_mean else y_mean
      scale0 <- if (!is.finite(scale0) || scale0 <= 0) y_mean else scale0
      c(log(max(shape0, 0.1)), log(max(scale0, 0.1)), 0, log(theta0))
    },
    weibull = c(0, log(max(mean(y), 0.1)), 0, log(theta0)),
    exp = c(log(max(1 / mean(y), 0.1)), 0, log(theta0))
  )
}

.comppareto_q_closure <- function(spec) {
  function(p, ...) {
    pars <- .comppareto_named_pars(spec, list(...))
    do.call(qcomppareto, c(list(p = p, spec = spec), pars))
  }
}

.comppareto_p_closure <- function(spec) {
  function(q, ...) {
    pars <- .comppareto_named_pars(spec, list(...))
    do.call(pcomppareto, c(list(q = q, spec = spec), pars))
  }
}

.comppareto_d_closure <- function(spec) {
  function(x, ...) {
    pars <- .comppareto_named_pars(spec, list(...))
    do.call(dcomppareto, c(list(x = x, spec = spec), pars))
  }
}

.comppareto_family_setup <- function(comppareto) {
  if (is.null(comppareto$spec)) {
    spec <- "lnorm"
  } else {
    spec <- match.arg(comppareto$spec, .comppareto_specs)
  }

  info <- .comppareto_body_info(spec)

  lik.fns <- list(
    d0 = function(pars, likdata) {
      if (likdata$censored) {
        stop("Censored likelihoods are not currently available for CompPareto models.")
      }
      eta <- .comppareto_eta_matrix(pars, likdata)
      sum(.comppareto_nll_obs(likdata$y[, 1], eta, spec))
    },
    d120 = function(pars, likdata) {
      eta <- .comppareto_eta_matrix(pars, likdata)
      .comppareto_d120(likdata$y[, 1], eta, spec)
    },
    d340 = NULL,
    q = .comppareto_q_closure(spec),
    p = .comppareto_p_closure(spec),
    d = .comppareto_d_closure(spec),
    unlink = info$unlink,
    inits = function(likdata) .comppareto_initial_values(spec, likdata),
    spec = spec
  )

  list(
    lik.fns = lik.fns,
    npar = length(info$gam_names),
    nms = info$gam_names,
    nms2 = info$gam_names
  )
}

#' Density of the continuous composite Pareto distribution
#'
#' @param x vector of quantiles.
#' @param spec body distribution specification. Currently one of
#'   \code{"lnorm"}, \code{"gamma"}, \code{"weibull"}, or \code{"exp"}.
#' @param alpha positive Pareto shape parameter.
#' @param theta positive splice point and Pareto scale parameter.
#' @param log logical; return log-density?
#' @param ... body-distribution parameters for the chosen \code{spec}.
#'
#' @return Numeric vector of densities.
#' @export
dcomppareto <- function(x, spec, alpha = 1, theta = 1, log = FALSE, ...) {
  spec <- match.arg(spec, .comppareto_specs)
  info <- .comppareto_body_info(spec)
  n <- length(x)
  alpha <- rep_len(alpha, n)
  theta <- rep_len(theta, n)
  body_args <- .comppareto_recycle_args(n, list(...))
  comp <- .comppareto_phi(spec, alpha = alpha, theta = theta, body_args = body_args)

  out <- rep(NA_real_, length(x))
  invalid <- !is.finite(alpha) | alpha <= 0 | !is.finite(theta) | theta <= 0 |
    !is.finite(comp$phi) | !is.finite(comp$G_theta)
  if (length(invalid) == 1L) {
    invalid <- rep(invalid, length(x))
  }

  idx_lower <- !invalid & is.finite(x) & x > 0 & x <= theta
  idx_upper <- !invalid & is.finite(x) & x > theta

  if (any(idx_lower)) {
    out[idx_lower] <- do.call(info$d, c(list(x[idx_lower]), .comppareto_subset_args(body_args, idx_lower))) /
      ((1 + comp$phi[idx_lower]) * comp$G_theta[idx_lower])
  }
  if (any(idx_upper)) {
    out[idx_upper] <- (comp$phi[idx_upper] / (1 + comp$phi[idx_upper])) *
      .comppareto_pareto_density(x[idx_upper], alpha = alpha[idx_upper], theta = theta[idx_upper]) /
      comp$tail_surv_theta[idx_upper]
  }

  out[!is.finite(x) | x <= 0] <- 0

  if (log) {
    out <- log(out)
  }

  out
}

#' Distribution function of the continuous composite Pareto distribution
#'
#' @inheritParams dcomppareto
#' @param q vector of quantiles.
#' @param lower.tail logical; if \code{FALSE}, return the upper tail.
#' @param log.p logical; return log-probabilities?
#'
#' @return Numeric vector of probabilities.
#' @export
pcomppareto <- function(q, spec, alpha = 1, theta = 1, lower.tail = TRUE, log.p = FALSE, ...) {
  spec <- match.arg(spec, .comppareto_specs)
  info <- .comppareto_body_info(spec)
  n <- length(q)
  alpha <- rep_len(alpha, n)
  theta <- rep_len(theta, n)
  body_args <- .comppareto_recycle_args(n, list(...))
  comp <- .comppareto_phi(spec, alpha = alpha, theta = theta, body_args = body_args)

  out <- rep(NA_real_, length(q))
  invalid <- !is.finite(alpha) | alpha <= 0 | !is.finite(theta) | theta <= 0 |
    !is.finite(comp$phi) | !is.finite(comp$G_theta)
  if (length(invalid) == 1L) {
    invalid <- rep(invalid, length(q))
  }

  idx_lower <- !invalid & is.finite(q) & q > 0 & q <= theta
  idx_upper <- !invalid & is.finite(q) & q > theta

  if (any(idx_lower)) {
    out[idx_lower] <- do.call(info$p, c(list(q[idx_lower]), .comppareto_subset_args(body_args, idx_lower))) /
      ((1 + comp$phi[idx_lower]) * comp$G_theta[idx_lower])
  }
  if (any(idx_upper)) {
    out[idx_upper] <- 1 / (1 + comp$phi[idx_upper]) +
      (comp$phi[idx_upper] / (1 + comp$phi[idx_upper])) *
      (.comppareto_pareto_cdf(q[idx_upper], alpha = alpha[idx_upper], theta = theta[idx_upper]) -
        comp$G2_theta[idx_upper]) / comp$tail_surv_theta[idx_upper]
  }

  out[!is.finite(q)] <- NA_real_
  out[q <= 0] <- 0

  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }

  out
}

#' Quantile function of the continuous composite Pareto distribution
#'
#' @inheritParams dcomppareto
#' @param p vector of probabilities.
#' @param log.p logical; if \code{TRUE}, \code{p} is supplied on the log scale.
#'
#' @return Numeric vector of quantiles.
#' @export
qcomppareto <- function(p, spec, alpha = 1, theta = 1, log.p = FALSE, ...) {
  spec <- match.arg(spec, .comppareto_specs)
  info <- .comppareto_body_info(spec)
  if (log.p) {
    p <- exp(p)
  }
  n <- length(p)
  alpha <- rep_len(alpha, n)
  theta <- rep_len(theta, n)
  body_args <- .comppareto_recycle_args(n, list(...))

  comp <- .comppareto_phi(spec, alpha = alpha, theta = theta, body_args = body_args)
  boundary <- 1 / (1 + comp$phi)

  out <- rep(NA_real_, length(p))
  invalid <- !is.finite(alpha) | alpha <= 0 | !is.finite(theta) | theta <= 0 |
    !is.finite(comp$phi) | !is.finite(comp$G_theta) | !is.finite(p) | p < 0 | p > 1
  if (length(invalid) == 1L) {
    invalid <- rep(invalid, length(p))
  }

  idx_lower <- !invalid & p <= boundary
  idx_upper <- !invalid & p > boundary & p < 1
  idx_one <- !invalid & p == 1

  if (any(idx_lower)) {
    p_body <- p[idx_lower] * (1 + comp$phi[idx_lower]) * comp$G_theta[idx_lower]
    out[idx_lower] <- do.call(info$q, c(list(p_body), .comppareto_subset_args(body_args, idx_lower)))
  }
  if (any(idx_upper)) {
    p_tail <- ((p[idx_upper] * (1 + comp$phi[idx_upper]) - 1) / comp$phi[idx_upper]) *
      comp$tail_surv_theta[idx_upper] + comp$G2_theta[idx_upper]
    out[idx_upper] <- .comppareto_pareto_quantile(
      p_tail,
      alpha = alpha[idx_upper],
      theta = theta[idx_upper]
    )
  }
  out[idx_one] <- Inf

  out
}

#' Random generation from the continuous composite Pareto distribution
#'
#' @inheritParams dcomppareto
#' @param n number of observations.
#'
#' @return Numeric vector of random draws.
#' @export
rcomppareto <- function(n, spec, alpha = 1, theta = 1, ...) {
  n <- as.integer(n[1])
  if (!is.finite(n) || n < 0) {
    stop("`n` must be a non-negative integer.")
  }
  qcomppareto(stats::runif(n), spec = spec, alpha = alpha, theta = theta, ...)
}
