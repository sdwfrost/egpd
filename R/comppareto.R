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
  1e-4 * pmax(1, abs(eta))
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
