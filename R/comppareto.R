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

.comppareto_recycle <- function(primary, alpha, theta, body_args) {
  all_args <- c(list(primary = primary, alpha = alpha, theta = theta), body_args)
  if (any(lengths(all_args) == 0L)) {
    return(list(
      n = 0L,
      primary = numeric(0),
      alpha = numeric(0),
      theta = numeric(0),
      body_args = lapply(body_args, function(x) numeric(0))
    ))
  }
  n <- max(lengths(all_args))
  list(
    n = n,
    primary = rep_len(primary, n),
    alpha = rep_len(alpha, n),
    theta = rep_len(theta, n),
    body_args = lapply(body_args, rep_len, length.out = n)
  )
}

.comppareto_log1mexp <- function(a) {
  out <- rep(NaN, length(a))
  isna <- is.na(a)
  out[isna] <- a[isna]
  out[!isna & a == 0] <- -Inf
  out[!isna & a == -Inf] <- 0
  ok <- is.finite(a) & a < 0
  if (any(ok)) {
    ao <- a[ok]
    out[ok] <- ifelse(ao < -log(2), log1p(-exp(ao)), log(-expm1(ao)))
  }
  out
}

.comppareto_subset_args <- function(args, idx) {
  lapply(args, function(x) x[idx])
}

## log(exp(a) + exp(b)), with -Inf terms allowed
.comppareto_logadd <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[!is.na(m) & m == -Inf] <- -Inf
  out
}

## log(2 theta / (x + theta)), x > theta. Near the splice even a one-ulp
## gap matters when multiplied by a large tail shape.
.comppareto_tail_logratio <- function(x, theta) {
  ifelse(x / theta <= 3, -log1p(0.5 * ((x - theta) / theta)),
         log(2) + log(theta) - log(x) - log1p(theta / x))
}

## Body-region (q <= theta) survivor from the body's own survivor S_G:
##   S(q) = [S_G(q) - S_G(theta) + phi G(theta)] / [(1 + phi) G(theta)],
## a sum of two non-negative terms, with S_G(q) - S_G(theta) = S_G(q)(1 - S_G(theta)/S_G(q))
## and log G(theta) from the body's log CDF. Forming it as 1 - F keeps only absolute
## accuracy, and log G(theta) taken from a linear G(theta) near 1 adds an error of
## S_G(theta): both swamp a survivor of order phi when the splice is far in the body's
## tail. Returns log S(q) for rows with parameters `body_args` (already subset).
.comppareto_body_lsurv <- function(info, q, theta, body_args, lphi, l1p) {
  lGt <- do.call(info$p, c(list(theta), body_args, list(log.p = TRUE)))
  lSt <- do.call(info$p, c(list(theta), body_args, list(lower.tail = FALSE, log.p = TRUE)))
  lSq <- do.call(info$p, c(list(q), body_args, list(lower.tail = FALSE, log.p = TRUE)))
  ldiff <- rep(-Inf, length(q))
  pos <- lSq > -Inf
  ldiff[pos] <- lSq[pos] + .comppareto_log1mexp(pmin(lSt[pos] - lSq[pos], 0))
  .comppareto_logadd(ldiff, lphi + lGt) - l1p - lGt
}

## Body-region quantile at an upper-tail target log S = ls (in the body, so
## ls >= log(phi / (1 + phi))), from the body's own survivor: solving the survivor above,
## S_G(q) = S_G(theta) + G(theta) phi (e^t - 1) with t = log((1 + phi) S / phi) >= 0.
.comppareto_body_qsurv <- function(info, ls, theta, body_args, lphi, lwt) {
  lGt <- do.call(info$p, c(list(theta), body_args, list(log.p = TRUE)))
  lSt <- do.call(info$p, c(list(theta), body_args, list(lower.tail = FALSE, log.p = TRUE)))
  t <- pmax(ls - lwt, 0)
  lem <- ifelse(t > 1, t + log1p(-exp(-t)), log(expm1(t)))   # log(e^t - 1)
  lSq <- pmin(.comppareto_logadd(lSt, lGt + lphi + lem), 0)
  q <- do.call(info$q, c(list(lSq), body_args, list(lower.tail = FALSE, log.p = TRUE)))
  pmin(q, theta)
}

## Rows whose parameters are NA (not NaN): these give NA silently. Any other invalid row
## gives NaN with a warning, as stats::dexp etc. do.
.comppareto_na_par <- function(alpha, theta, body_args) {
  Reduce(`|`, lapply(c(list(alpha, theta), body_args), function(v) is.na(v) & !is.nan(v)),
         FALSE)
}

.comppareto_invalid_out <- function(n, invalid, na_par) {
  out <- rep(NA_real_, n)
  bad <- invalid & !na_par
  if (any(bad)) {
    warning("NaNs produced", call. = FALSE)
    out[bad] <- NaN
  }
  out
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
  G_theta <- suppressWarnings(do.call(info$p, c(list(theta), body_args)))
  g_theta <- suppressWarnings(do.call(info$d, c(list(theta), body_args)))
  tail_surv_theta <- exp(-alpha * log(2))
  G2_theta <- -expm1(-alpha * log(2))
  phi <- 2 * theta * g_theta / (alpha * G_theta)
  lphi <- l1p <- lG <- rep(NA_real_, length(phi))
  pos <- !is.na(phi) & phi > 0
  lphi[pos] <- log(phi[pos])
  l1p[pos] <- log1p(phi[pos])
  gpos <- !is.na(G_theta) & G_theta >= 0
  lG[gpos] <- log(G_theta[gpos])
  ## Where g(theta), G(theta) or phi leave the normal range, the linear ratio has lost its
  ## relative accuracy, or is 0, Inf or NaN, although phi is a valid number: there log phi
  ## is formed from log g(theta) and log G(theta) (the `redo` rows; elsewhere the linear
  ## values and their logs are used unchanged).
  xmin <- .Machine$double.xmin
  redo <- is.finite(alpha) & alpha > 0 & is.finite(theta) & theta > 0 &
    !is.na(g_theta) & !is.na(G_theta) &
    !(g_theta >= xmin & G_theta >= xmin & !is.na(phi) & phi >= xmin & phi < Inf)
  if (any(redo)) {
    ba <- .comppareto_subset_args(body_args, redo)
    lg_r <- suppressWarnings(do.call(info$d, c(list(theta[redo]), ba, list(log = TRUE))))
    lG_r <- suppressWarnings(do.call(info$p, c(list(theta[redo]), ba, list(log.p = TRUE))))
    lphi[redo] <- log(2) + log(theta[redo]) + lg_r - log(alpha[redo]) - lG_r
    l1p[redo] <- pmax(lphi[redo], 0) + log1p(exp(-abs(lphi[redo])))
    lG[redo] <- lG_r
    phi[redo] <- exp(lphi[redo])
    G_theta[redo] <- exp(lG_r)
  }
  bad <- !(is.finite(lphi) & is.finite(lG) & is.finite(l1p))
  phi[bad] <- lphi[bad] <- l1p[bad] <- lG[bad] <- NA_real_
  ## log(phi / (1 + phi)) without subtracting nearly equal log weights.
  list(phi = phi, lphi = lphi, l1p = l1p, lwt = plogis(lphi, log.p = TRUE),
       G_theta = G_theta, lG = lG,
       redo = redo & !bad, G2_theta = G2_theta, tail_surv_theta = tail_surv_theta)
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
##   theta    splice point (= Lomax scale)
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
    log_ratio <- .comppareto_tail_logratio(yt, th)
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
  split_pars <- split(pars, factor(likdata$idpars, levels = seq_along(likdata$X)))
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

## Per-observation negative log-likelihood. The density is zero at y <= 0 (the support
## is x > 0), so such an observation contributes +Inf, as does any parameter value at
## which the density is undefined; the fitter treats a non-finite objective as a failed
## point. (Clamping the density at DBL_MIN instead turned every y <= 0 into a constant
## 708.4 that the analytic score did not see.)
.comppareto_nll_obs <- function(y, eta, spec) {
  pars <- .comppareto_eta_to_response(spec, eta)
  lp <- suppressWarnings(do.call(
    dcomppareto,
    c(list(x = y, spec = spec, log = TRUE), pars)
  ))
  out <- -lp
  out[is.na(out)] <- Inf
  out
}

## Rows outside the support (y <= 0) have an infinite, parameter-free contribution, so
## their gradient and Hessian are zero, consistently in the analytic and FD paths.
.comppareto_zero_offsupport <- function(out, y) {
  off <- !(is.finite(y) & y > 0)
  if (any(off)) out[off, ] <- 0
  out
}

.comppareto_d120 <- function(y, eta, spec) {
  if (.comppareto_has_analytic(spec)) {
    return(.comppareto_zero_offsupport(.comppareto_assemble_d120(y, eta, spec), y))
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

  .comppareto_zero_offsupport(cbind(grad, do.call(cbind, hess_terms)), y)
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

## Further starting values (link scale, one intercept per linear predictor) for the
## restarts of egpd(family = "comppareto"): the default start's body parameters, the
## splice point theta at the 50%, 75%, 90% and 95% quantiles of the response, and alpha
## at the maximum likelihood estimate of the conditional Lomax tail above theta,
## n_t / sum(log((y + theta) / (2 theta))) over the n_t observations above it (1 with
## fewer than two of them), kept in [0.05, 50]. The default start's alpha = 1 is a much
## heavier tail than most data have.
.comppareto_restart_starts <- function(spec, y) {
  y <- as.numeric(y)
  y <- y[is.finite(y) & y > 0]
  if (length(y) < 2L)
    return(list())
  d0 <- .comppareto_initial_values(spec, list(y = matrix(y)))
  nb <- length(d0) - 2L
  lapply(c(0.5, 0.75, 0.9, 0.95), function(q) {
    th <- max(stats::quantile(y, probs = q, names = FALSE), .Machine$double.eps)
    yt <- y[y > th]
    al <- if (length(yt) >= 2L) length(yt) / sum(log((yt + th) / (2 * th))) else 1
    c(d0[seq_len(nb)], log(min(max(al, 0.05), 50)), log(th))
  })
}

## Why egpd(family = "comppareto") should try the restarts' other starts after its
## default fit, or NULL if it should not: every observation lies on one side of the
## fitted splice point. With no observation above theta (or none at or below it) the
## likelihood does not identify alpha and theta, and it keeps increasing towards a
## boundary of the parameter space -- with none above, towards the body distribution
## truncated at the largest observation as alpha -> Inf -- so such a fit is never an
## interior maximum, however small its gradient.
.comppareto_restart_reason <- function(fit) {
  y <- fit$likdata$y[, 1]
  lt <- fit[[length(fit$likdata$X)]]$fitted
  if (!length(y) || length(lt) != length(y) || anyNA(lt))
    return(NULL)
  above <- y > exp(lt)
  if (!any(above))
    return("has no observation above its fitted splice point")
  if (all(above))
    return("has no observation at or below its fitted splice point")
  NULL
}

## Diagnose the split that survived the default search, not every exploratory
## start. With no tail observations, increasing an unpenalised log-alpha intercept
## strictly increases the likelihood at any finite alpha. This conclusion need
## not apply to a pinned alpha or a curvature-transformed objective.
.comppareto_check_split <- function(fit) {
  why <- .comppareto_restart_reason(fit)
  if (is.null(why) || isTRUE(fit$convergence == 2L)) return(fit)
  ld <- fit$likdata
  p <- length(ld$X)
  ia <- which(ld$idpars == p-1L)[match("(Intercept)", colnames(ld$X[[p-1L]]))]
  if (any(ld$y[, 1L] > exp(fit[[p]]$fitted)) || ld$adjust == 2L ||
      !length(ia) || is.na(ia) || (!is.null(ld$S) && any(ld$S[ia, ] != 0))) {
    fit$splice.boundary <- paste("the selected fit", why, "(one side is unobserved)")
    return(fit)
  }
  msg <- paste("the selected fit", why, "and the likelihood increases as alpha tends to infinity; no finite tail-shape maximum at this split")
  fit$splice.boundary <- fit$convergence.report <- msg
  fit$convergence <- 3L
  fit$logLik <- NA_real_
  fit$logLik.adjusted <- NULL
  fit$AIC <- fit$BIC <- 1e20
  fit
}

## For a common splice t = log(theta), C = 2 theta g(theta)/G(theta), and
## A = sum_tail log((y+theta)/(2 theta)), the alpha-dependent part of the full
## log likelihood is -n log(1+C/alpha) - A alpha. At a split with an observed
## tail it has a unique positive maximizer, alpha(alpha+C) = n C/A. Profiling
## this exactly leaves at most two smooth coordinates and no design/Hessian work.
.comppareto_split_point <- function(t, body, ly, info, free, la, polish = FALSE) {
  ib <- seq_along(body)
  lower <- ly <= t
  v <- ly[!lower] - t
  A <- sum(v + log1p(expm1(-v)/2))
  n <- length(ly)
  nt <- length(v)
  yb <- exp(ly[lower])
  th <- exp(t)
  log_alpha <- function(lc) {
    if (!free[length(free)]) return(la)
    q <- log(n) - log(A)
    log(2) + q - .comppareto_logadd(0, .comppareto_logadd(0, log(4) + q - lc)/2)
  }
  if (info$spec == "exp" && free[1L]) {
    ## Jointly concave in z = rate*theta and log(alpha): see the normalizer
    ## in .comppareto_rate_boundary. Root the profile score in z itself; on the
    ## log-rate scale a nearly uniform body can stall far from a finite optimum.
    D <- nt + sum(yb/th)
    score <- function(z) {
      lI <- if (z == 0) 0 else log(-expm1(-z)/z)
      lc <- log(2) - z - lI
      w <- stats::plogis(lc - log_alpha(lc))
      m <- if (z < 1e-4) .5 - z/12 else 1/z - 1/expm1(z)
      n * (w + (1-w)*m) - D
    }
    z <- if (score(0) <= 0) sqrt(.Machine$double.eps) else
      stats::uniroot(score, c(0, 1), extendInt = "downX", tol = 1e-9)$root
    body <- log(z) - t
    polish <- FALSE
  }
  evaluate <- function(b) {
    args <- Map(function(f, x) f(x), info$unlink[ib], b)
    names(args) <- info$body_names  # gamma's second body parameter is scale, not rate
    lG <- do.call(info$p, c(list(th), args, log.p = TRUE))
    lg <- do.call(info$d, c(list(th), args, log = TRUE))
    lc <- log(2) + t + lg - lG
    a <- log_alpha(lc)
    ll <- sum(do.call(info$d, c(list(yb), args, log = TRUE)) - lG) +
      nt * (lc - log(2) - t) - n * .comppareto_logadd(0, lc - a) - (exp(a) + 1)*A
    c(ll = if (.egpd_invalid_objective(-ll)) -Inf else ll, b, a, t)
  }
  best <- evaluate(body)
  i <- which(free[ib])
  if (polish && length(i)) {
    objective <- function(x) -evaluate(replace(body, i, x))[1L]
    o <- try(stats::optim(body[i], objective, method = "BFGS",
                          control = list(maxit = 60L, reltol = 1e-10)), silent = TRUE)
    if (!inherits(o, "try-error") && !.egpd_invalid_objective(o$value) && -o$value > best[1L])
      best <- evaluate(replace(body, i, o$par))
  }
  unname(best)
}

## A bounded scan of actual knots and midpoints in *adjacent open pieces* (not
## the midpoint between two thinned knots, which can cross many observations).
.comppareto_split_grid <- function(breaks, size) {
  u <- sort(unique(breaks))
  if (length(u) < 3L) return(numeric())
  k <- unique(round(seq(1, length(u) - 1L, length.out = min(length(u) - 1L, size))))
  sort(unique(c(u[k], u[k] + (u[k+1L] - u[k])/2)))
}

## Score cheaply; only the best separated peaks with evidence of a material gain
## pay for a full call to egpd(). With covariates, move the theta intercept of the
## *fitted predictors* (so the knots are log(y_i) - eta_theta_i); this is a line
## search, not an enumeration of multidimensional splice regions. The refit uses
## ordinary intercept starts and re-estimates the other coefficients and any sp.
.comppareto_split_starts <- function(fit) {
  ld <- fit$likdata
  X <- ld$X
  p <- length(X)
  intercept <- vapply(X, function(x) ncol(x) > 0L && identical(colnames(x)[1L], "(Intercept)"),
                      logical(1))
  if (!intercept[p]) return(list())
  y <- ld$y[, 1L]
  eta <- .comppareto_eta_matrix(fit$coefficients, ld)
  ll0 <- -fit$likfns$d0(fit$coefficients, ld)
  if (.egpd_invalid_objective(-ll0)) return(list())
  spec <- fit$likfns$spec
  common <- all(vapply(X, ncol, 0L) == 1L) &&
    all(intercept | vapply(X, .is_pinned_block, logical(1))) &&
    all(vapply(ld$offsets, function(o) length(o) < 2L || all(o == o[1L]), logical(1)))
  if (common) {
    ly <- sort(log(y))
    grid <- .comppareto_split_grid(ly, 128L)
    if (!length(grid)) return(list())
    info <- .comppareto_body_info(spec)
    ib <- seq_len(p-2L)
    b0 <- .comppareto_initial_values(spec, ld)[ib]
    b0[!intercept[ib]] <- eta[1L, ib][!intercept[ib]]
    point <- function(t, b, polish = FALSE)
      .comppareto_split_point(t, b, ly, info, intercept[-p], eta[1L, p-1L], polish)
    ## Conditional fits at a few ranks give the scan body shapes appropriate
    ## to different splits. Optimising the nuisance parameters at every rank,
    ## let alone running a full GAM at each, is unnecessary.
    if (spec == "exp" || !any(intercept[ib])) {
      scan <- lapply(grid, point, b = b0)
    } else {
      j <- unique(round(seq(1, length(grid), length.out = min(11L, length(grid)))))
      anchors <- lapply(grid[j], point, b = b0, polish = TRUE)
      body <- vapply(ib, function(i) stats::approx(grid[j], vapply(anchors, `[`, 0, i+1L),
                                                   xout = grid)$y, numeric(length(grid)))
      scan <- lapply(seq_along(grid), function(k) {
        r <- point(grid[k], body[k, ])
        s <- point(grid[k], eta[1L, ib])
        if (s[1L] > r[1L]) s else r
      })
    }
    offset <- vapply(ld$offsets, function(o) if (length(o)) o[1L] else 0, 0)
    start <- function(k) {
      v <- point(grid[k], scan[[k]][1L + ib], polish = TRUE)[-1L] - offset
      replace(v, !intercept, 0)
    }
  } else {
    grid <- .comppareto_split_grid(log(y) - eta[, p], 32L)
    if (!length(grid)) return(list())
    scan <- lapply(grid, function(shift) {
      e <- eta
      e[, p] <- eta[, p] + shift
      -sum(.comppareto_nll_obs(y, e, spec))
    })
    start <- function(k) {
      v <- .intercept_starts(fit)
      v[p] <- v[p] + grid[k]
      replace(v, !intercept, 0)
    }
  }
  ll <- vapply(scan, `[`, 0, 1L)
  ll[!is.finite(ll)] <- -Inf
  peaks <- which(ll > ll0 + 0.01 &
                 ll >= c(-Inf, ll[-length(ll)]) & ll >= c(ll[-1L], -Inf))
  peaks <- peaks[order(ll[peaks], decreasing = TRUE)]
  lapply(utils::head(peaks, 2L), start)
}

## At fixed alpha and theta the exponential-body log likelihood is strictly concave
## in a common positive multiplier of the rates. For z = rate * theta its normalizer
## (apart from parameter-independent factors) is integral_0^1 exp(-z*t) dt +
## 2*exp(-z)/alpha. The sufficient statistic for t is min(y/theta, 1). A non-positive
## score at multiplier zero thus rules out a finite optimum along this direction,
## even where the *log* rate score is numerically zero. Only use this for a free,
## unpenalised log-rate intercept and an unmixed coefficient scale.
.comppareto_rate_boundary <- function(b, ld) {
  ir <- which(ld$idpars == 1L)[match("(Intercept)", colnames(ld$X[[1L]]))]
  if (!length(ir) || is.na(ir) || (!is.null(ld$S) && any(ld$S[ir, ] != 0)))
    return(FALSE)
  eta <- .comppareto_eta_matrix(b, ld)
  lz <- eta[, 1L] + eta[, 3L]
  mean0 <- (1 + stats::plogis(log(2) - eta[, 2L])) / 2
  t <- exp(pmin(0, log(ld$y[, 1L]) - eta[, 3L]))
  isTRUE(sum(exp(lz - max(lz)) * (mean0 - t)) <= 0)
}

## Strict one-sided derivatives *outward* from the knot of the objective being
## minimised. Use two points inside each adjacent open piece; extrapolating their
## secants removes the leading curvature bias (one uphill trial is not a slope test).
## Unresolved signs, or even one invalid trial on either side, certify nothing.
.comppareto_kink_sides <- function(x, breaks, objective) {
  f0 <- objective(x)
  if (.egpd_invalid_objective(f0)) return(FALSE)
  all(vapply(c(-1, 1), function(dir) {
    gap <- dir * (breaks - x)
    h <- min(1e-4, gap[gap > 0] / 2)
    v <- vapply(c(h/2, h), function(s) objective(x + dir*s), 0)
    if (any(vapply(v, .egpd_invalid_objective, logical(1)))) return(FALSE)
    noise <- 32 * .Machine$double.eps * max(1, abs(f0), abs(v))
    isTRUE(v[1L] >= f0 && 4*(v[1L] - f0) - (v[2L] - f0) > noise)
  }, logical(1)))
}

## The ordinary Newton decrement uses a single smooth branch. At a scalar splice
## knot this can predict an impossible improvement across the kink and give code 3.
## Snap to the actual knot, fit/check *all* free tangent coordinates (not only the
## columns that Newton kept), then check both normals. For a varying theta design
## the active surfaces need not be parallel: coordinatewise tests are insufficient.
.comppareto_finish_exp <- function(fit, fn, sfn, likdata, likfns) {
  ld <- likdata
  if (!(fit$convergence %in% c(0L, 1L, 3L)) || ld$adjust == 2L ||
      .egpd_invalid_objective(fit$objective)) return(fit)
  b <- as.vector(fit$pars)
  fail <- function(reason) {
    ## An exhausted search stays an iteration-limit failure unless the knot
    ## passes the entire stationarity check, just as for an early stop there.
    if (fit$convergence == 1L) return(fit)
    fit$convergence <- 3L
    fit$gradconv <- FALSE
    fit$report <- reason
    fit
  }
  boundary <- "exponential body at a boundary: the likelihood increases as rate tends to zero at this splice and tail shape; no finite rate maximum there"
  it <- which(ld$idpars == 3L)
  if (length(it) != 1L || !all(ld$X[[3L]] == 1)) return(fit)
  breaks <- log(ld$y[, 1L])
  if (length(ld$offsets[[3L]])) breaks <- breaks - ld$offsets[[3L]]
  knot <- breaks[which.min(abs(breaks - b[it]))]
  if (abs(b[it] - knot) > 1e-7 * max(1, abs(knot))) {
    if (.comppareto_rate_boundary(b, ld)) return(fail(boundary))
    return(fit)
  }
  if (!any(breaks > knot) || !any(breaks < knot)) return(fit)
  b[it] <- knot
  pinned <- vapply(ld$X, .is_pinned_block, logical(1))
  tangent <- which(ld$idpars != 3L & !pinned[ld$idpars])
  obj <- function(p) as.numeric(fn(p, likdata = ld, likfns = likfns))
  gh <- if (is.null(ld$S)) .gH.nopen else .gH.pen
  deriv <- function(p, order) gh(p, ld, likfns, deriv = order)
  scale <- .coef_scale(ld)[tangent]
  tol <- min(1, ld$k) * min(1e-4, ld$control$inner$gradtol)
  stationary <- function(p) {
    if (!length(tangent)) return(TRUE)
    d <- deriv(p, 2)
    g <- as.vector(d[[1L]])[tangent] / scale
    H <- d[[2L]][tangent, tangent, drop = FALSE]
    all(is.finite(g)) && max(abs(g)) <= tol &&
      !inherits(try(chol(H), silent = TRUE), "try-error")
  }
  if (!stationary(b)) {
    full <- function(v) { p <- b; p[tangent] <- v; p }
    o <- try(stats::optim(b[tangent], function(v) obj(full(v)),
                          function(v) as.vector(deriv(full(v), 1)[[1L]])[tangent],
                          method = "BFGS", control = list(reltol = 1e-15,
                            maxit = min(100L, ld$control$inner$itlim))), silent = TRUE)
    if (!inherits(o, "try-error") && !.egpd_invalid_objective(o$value) && o$value <= obj(b))
      b <- full(o$par)
  }
  if (.comppareto_rate_boundary(b, ld)) return(fail(boundary))
  ## Exact jump in the slope of the *negative* log likelihood: only s_body
  ## versus s_tail changes, not their common smooth normalizer (or a penalty).
  ## exp(log(y)) need not equal y in floating point: correct from the side
  ## actually used by d120 for each tied row.
  eta <- .comppareto_eta_matrix(b, ld)
  tied <- breaks == knot
  used.body <- ld$y[, 1L] <= exp(eta[, 3L])
  jump <- ld$k * ((1 + exp(eta[, 2L]))/2 - exp(eta[, 1L] + eta[, 3L]))
  g <- as.vector(deriv(b, 1)[[1L]])[it]
  left <- g - sum(jump[tied & used.body])
  right <- g + sum(jump[tied & !used.body])
  at <- function(x) { b[it] <- x; obj(b) }
  if (!stationary(b) || !isTRUE(left < 0 && right > 0) ||
      !.comppareto_kink_sides(knot, breaks, at))
    return(fail("splice knot is not a verified minimum of the negative log-likelihood in the smooth and both one-sided directions"))
  if (.egpd_invalid_objective(obj(b)) ||
      obj(b) > fit$objective + 16 * .Machine$double.eps * max(1, abs(fit$objective)))
    return(fail("stationary splice knot not reached at the returned objective"))
  if (fit$convergence == 0L && identical(b[-it], fit$pars[-it]) &&
      abs(b[it] - fit$pars[it]) <= 16 * .Machine$double.eps * max(1, abs(knot)))
    return(fit)
  out <- .refresh_inner(fit, b, fn, sfn, ld, likfns)
  if (!identical(out$pars, b) || .egpd_invalid_objective(out$objective))
    return(fail("could not evaluate derivatives at the stationary splice knot"))
  out$convergence <- 0L
  out$gradconv <- FALSE  # directional stationarity, not a zero ordinary gradient
  out$report <- "stationary splice knot (smooth coefficients and both one-sided objective slopes checked)"
  out
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
    finish_inner = if (spec == "exp") .comppareto_finish_exp,
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
#' The implemented tail is a Lomax (Pareto-II) density proportional to
#' \eqn{(x+\theta)^{-\alpha-1}}, conditioned on \eqn{x>\theta}; its conditional
#' survivor is \eqn{[2\theta/(x+\theta)]^\alpha}. It is not a Pareto-I tail
#' with survivor \eqn{(\theta/x)^\alpha}. The splice weight is
#' \eqn{\phi=2\theta g(\theta)/[\alpha G(\theta)]}; it is formed on the log scale
#' where \eqn{g(\theta)}, \eqn{G(\theta)} or \eqn{\phi} underflow or overflow, so a
#' splice point far in a tail of the body still gives a finite density.
#'
#' The support is \eqn{x > 0}: the density is 0 (log-density \code{-Inf}) at
#' \eqn{x \le 0} and at \eqn{x = \infty}, \code{pcomppareto} is 0 at \eqn{q \le 0}, and
#' \code{qcomppareto(0)} is 0, so a model fitted with \code{family = "comppareto"}
#' requires strictly positive responses. A parameter that is \code{NA} gives \code{NA};
#' an invalid parameter (non-positive \code{alpha} or \code{theta}, invalid body
#' parameters, \code{NaN}) gives \code{NaN} with a warning, as do probabilities outside
#' \eqn{[0, 1]} in \code{qcomppareto}.
#'
#' @param x vector of quantiles.
#' @param spec body distribution specification. Currently one of
#'   \code{"lnorm"}, \code{"gamma"}, \code{"weibull"}, or \code{"exp"}.
#' @param alpha positive Pareto shape parameter.
#' @param theta positive splice point and Lomax scale parameter.
#' @param log logical; return log-density?
#' @param ... body-distribution parameters for the chosen \code{spec}.
#'
#' @return Numeric vector of densities.
#' @export
dcomppareto <- function(x, spec, alpha = 1, theta = 1, log = FALSE, ...) {
  spec <- match.arg(spec, .comppareto_specs)
  info <- .comppareto_body_info(spec)
  rec <- .comppareto_recycle(x, alpha, theta, list(...))
  if (rec$n == 0L) return(numeric(0))
  x <- rec$primary
  alpha <- rec$alpha
  theta <- rec$theta
  body_args <- rec$body_args
  comp <- .comppareto_phi(spec, alpha = alpha, theta = theta, body_args = body_args)

  invalid <- !is.finite(alpha) | alpha <= 0 | !is.finite(theta) | theta <= 0 |
    !is.finite(comp$lphi) | !is.finite(comp$lG)

  out <- .comppareto_invalid_out(rec$n, invalid,
                                 .comppareto_na_par(alpha, theta, body_args))
  isna_x <- !invalid & is.na(x)
  out[isna_x] <- x[isna_x]

  zero_idx <- !invalid & !is.na(x) & (x <= 0 | !is.finite(x))
  out[zero_idx] <- if (log) -Inf else 0

  idx_lower <- !invalid & is.finite(x) & x > 0 & x <= theta
  idx_upper <- !invalid & is.finite(x) & x > theta

  if (any(idx_lower)) {
    ii <- idx_lower & (log | comp$redo)
    if (any(ii)) {
      log_body <- do.call(
        info$d,
        c(list(x[ii]), .comppareto_subset_args(body_args, ii), list(log = TRUE))
      )
      out[ii] <- log_body - comp$l1p[ii] - comp$lG[ii]
      if (!log) out[ii] <- exp(out[ii])
    }
    ii <- idx_lower & !(log | comp$redo)
    if (any(ii)) {
      out[ii] <- do.call(
        info$d,
        c(list(x[ii]), .comppareto_subset_args(body_args, ii))
      ) / ((1 + comp$phi[ii]) * comp$G_theta[ii])
    }
  }
  if (any(idx_upper)) {
    ii <- idx_upper
    log_ratio <- .comppareto_tail_logratio(x[ii], theta[ii])
    x_p_th <- as.double(x[ii]) + theta[ii]
    log_x_p_th <- ifelse(
      is.finite(x_p_th) & x[ii] / theta[ii] <= 3,
      log(x_p_th),
      log(x[ii]) + log1p(theta[ii] / x[ii])
    )
    ## The conditional density can overflow even when its weighted value is finite.
    ld <- comp$lwt[ii] + log(alpha[ii]) - log_x_p_th + alpha[ii] * log_ratio
    out[ii] <- if (log) ld else exp(ld)
  }

  out
}

#' Distribution function of the continuous composite Pareto distribution
#'
#' Below the splice point the upper tail
#' \eqn{S(q) = [S_G(q) - S_G(\theta) + \phi G(\theta)] / [(1 + \phi) G(\theta)]} is computed
#' from the body's own survivor \eqn{S_G} where \eqn{S(q)} is below about \code{1e-3}, as is
#' \code{log.p = TRUE} near 0, so that both keep full relative accuracy when the splice
#' point lies far in the body's tail and the Pareto tail weight \eqn{\phi / (1 + \phi)} is
#' tiny.
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
  rec <- .comppareto_recycle(q, alpha, theta, list(...))
  if (rec$n == 0L) return(numeric(0))
  q <- rec$primary
  alpha <- rec$alpha
  theta <- rec$theta
  body_args <- rec$body_args
  comp <- .comppareto_phi(spec, alpha = alpha, theta = theta, body_args = body_args)

  invalid <- !is.finite(alpha) | alpha <= 0 | !is.finite(theta) | theta <= 0 |
    !is.finite(comp$lphi) | !is.finite(comp$lG)

  out <- .comppareto_invalid_out(rec$n, invalid,
                                 .comppareto_na_par(alpha, theta, body_args))
  isna_q <- !invalid & is.na(q)
  out[isna_q] <- q[isna_q]

  idx_le0 <- !invalid & !is.na(q) & q <= 0
  if (any(idx_le0)) {
    out[idx_le0] <- if (lower.tail) {
      if (log.p) -Inf else 0
    } else {
      if (log.p) 0 else 1
    }
  }

  idx_inf <- !invalid & !is.na(q) & q == Inf
  if (any(idx_inf)) {
    out[idx_inf] <- if (lower.tail) {
      if (log.p) 0 else 1
    } else {
      if (log.p) -Inf else 0
    }
  }

  idx_lower <- !invalid & is.finite(q) & q > 0 & q <= theta
  idx_upper <- !invalid & is.finite(q) & q > theta

  if (any(idx_lower)) {
    lin <- lower.tail && !log.p
    ii <- idx_lower & lin & !comp$redo
    if (any(ii)) {
      out[ii] <- do.call(
        info$p,
        c(list(q[ii]), .comppareto_subset_args(body_args, ii))
      ) / ((1 + comp$phi[ii]) * comp$G_theta[ii])
    }
    ii <- idx_lower & !(lin & !comp$redo)
    if (any(ii)) {
      log_Gq <- do.call(
        info$p,
        c(list(q[ii]), .comppareto_subset_args(body_args, ii), list(log.p = TRUE))
      )
      log_F <- pmin(log_Gq - comp$l1p[ii] - comp$lG[ii], 0)
      out[ii] <- if (lower.tail) {
        if (log.p) log_F else exp(log_F)
      } else if (log.p) {
        .comppareto_log1mexp(log_F)
      } else {
        -expm1(log_F)
      }
      ## a survivor below about 1e-3: the upper tail, and the lower tail on the log
      ## scale, come from the body's own survivor instead of 1 - F
      if (!lower.tail || log.p) {
        jj <- which(ii)[log_F > log1p(-1e-3)]
        if (length(jj)) {
          lS <- .comppareto_body_lsurv(info, q[jj], theta[jj],
                                       .comppareto_subset_args(body_args, jj),
                                       comp$lphi[jj], comp$l1p[jj])
          out[jj] <- if (lower.tail) .comppareto_log1mexp(lS) else if (log.p) lS else exp(lS)
        }
      }
    }
  }
  if (any(idx_upper)) {
    ii <- idx_upper
    log_ratio <- .comppareto_tail_logratio(q[ii], theta[ii])
    log_surv <- comp$lwt[ii] + alpha[ii] * log_ratio
    out[ii] <- if (!lower.tail) {
      if (log.p) log_surv else exp(log_surv)
    } else {
      if (log.p) .comppareto_log1mexp(log_surv) else -expm1(log_surv)
    }
    if (lower.tail && log.p) {
      ## A subnormal log survivor no longer retains the relative accuracy
      ## of the CDF. Add the body and tail increment on their log scales.
      sub <- which(-log_surv < .Machine$double.xmin)
      if (length(sub)) {
        j <- which(ii)[sub]
        inc <- log(alpha[j]) + log(-log_ratio[sub])
        out[j] <- .comppareto_logadd(-comp$l1p[j], comp$lwt[j] + inc)
      }
    }
  }

  out
}

#' Quantile function of the continuous composite Pareto distribution
#'
#' Upper-tail and \code{log.p} targets below the splice point whose body probability is
#' within about \code{1e-3} of 1 are inverted on the body's own survivor, matching
#' \code{\link{pcomppareto}}.
#'
#' @inheritParams dcomppareto
#' @param p vector of probabilities.
#' @param lower.tail logical; if \code{FALSE}, probabilities are \eqn{P[X > x]}.
#' @param log.p logical; if \code{TRUE}, \code{p} is supplied on the log scale.
#'
#' @return Numeric vector of quantiles.
#' @export
qcomppareto <- function(p, spec, alpha = 1, theta = 1, lower.tail = TRUE, log.p = FALSE, ...) {
  spec <- match.arg(spec, .comppareto_specs)
  info <- .comppareto_body_info(spec)
  rec <- .comppareto_recycle(p, alpha, theta, list(...))
  if (rec$n == 0L) return(numeric(0))
  p <- rec$primary
  alpha <- rec$alpha
  theta <- rec$theta
  body_args <- rec$body_args

  comp <- .comppareto_phi(spec, alpha = alpha, theta = theta, body_args = body_args)

  invalid_par <- !is.finite(alpha) | alpha <= 0 | !is.finite(theta) | theta <= 0 |
    !is.finite(comp$lphi) | !is.finite(comp$lG)

  out <- .comppareto_invalid_out(rec$n, invalid_par,
                                 .comppareto_na_par(alpha, theta, body_args))
  isna_p <- !invalid_par & is.na(p)
  out[isna_p] <- p[isna_p]

  valid_p <- if (log.p) {
    !invalid_par & !is.na(p) & p <= 0 & p != Inf
  } else {
    !invalid_par & is.finite(p) & p >= 0 & p <= 1
  }
  oor <- !invalid_par & !is.na(p) & !valid_p
  if (any(oor)) {
    warning("NaNs produced", call. = FALSE)
    out[oor] <- NaN
  }

  if (!any(valid_p)) return(out)

  vp <- which(valid_p)
  pv <- p[vp]
  if (log.p) {
    if (lower.tail) {
      log_cdf <- pv
      log_surv <- .comppareto_log1mexp(pv)
    } else {
      log_surv <- pv
      log_cdf <- .comppareto_log1mexp(pv)
    }
  } else {
    if (lower.tail) {
      log_cdf <- log(pv)
      log_surv <- log1p(-pv)
    } else {
      log_surv <- log(pv)
      log_cdf <- log1p(-pv)
    }
  }

  phi_v <- comp$phi[vp]
  Gth_v <- comp$G_theta[vp]
  lphi_v <- comp$lphi[vp]
  l1p_v <- comp$l1p[vp]
  lwt_v <- comp$lwt[vp]
  lG_v  <- comp$lG[vp]
  redo_v <- comp$redo[vp]
  al_v  <- alpha[vp]
  th_v  <- theta[vp]

  is_zero <- log_cdf == -Inf
  is_inf  <- log_surv == -Inf
  in_body <- !is_zero & !is_inf & (
    if (lower.tail) log_cdf <= -l1p_v else log_surv >= lwt_v
  )
  in_tail <- !is_zero & !is_inf & !in_body

  res <- rep(NA_real_, length(vp))
  res[is_zero] <- 0
  res[is_inf]  <- Inf

  if (any(in_body)) {
    lin <- !log.p && lower.tail
    ib <- which(in_body & lin & !redo_v)
    if (length(ib)) {
      p_body <- pmin(pv[ib] * (1 + phi_v[ib]) * Gth_v[ib], 1)
      res[ib] <- do.call(
        info$q,
        c(list(p_body), .comppareto_subset_args(body_args, vp[ib]))
      )
    }
    ib <- which(in_body & !(lin & !redo_v))
    if (length(ib)) {
      lp_body <- pmin(log_cdf[ib] + l1p_v[ib] + lG_v[ib], 0)
      ## a body target within 1e-3 of G = 1, given through its upper tail or on the log
      ## scale: invert the body's own survivor (lower-tail body quantiles lose their
      ## accuracy there, and 1 - p its relative accuracy)
      near <- !lin & !is.na(lp_body) & lp_body > log1p(-1e-3)
      ibo <- ib[!near]
      if (length(ibo)) {
        res[ibo] <- do.call(
          info$q,
          c(list(lp_body[!near]), .comppareto_subset_args(body_args, vp[ibo]), list(log.p = TRUE))
        )
      }
      ibn <- ib[near]
      if (length(ibn)) {
        res[ibn] <- .comppareto_body_qsurv(info, log_surv[ibn], th_v[ibn],
                                           .comppareto_subset_args(body_args, vp[ibn]),
                                           lphi_v[ibn], lwt_v[ibn])
      }
    }
  }
  if (any(in_tail)) {
    it <- which(in_tail)
    log_tail_surv <- log_surv[it] - lwt_v[it]
    hazard <- pmax(-log_tail_surv / al_v[it], 0)
    if (lower.tail && log.p) {
      sub <- which(-log_surv[it] < .Machine$double.xmin)
      if (length(sub)) {
        j <- it[sub]
        ## Increment = (F-F(theta))/tail weight. At subnormal F this
        ## also equals the conditional cumulative hazard to working precision.
        linc <- log_cdf[j] + .comppareto_log1mexp(pmin(-l1p_v[j] - log_cdf[j], 0))
        hazard[sub] <- exp(linc - lwt_v[j] - log(al_v[j]))
      }
    }
    ## Scale the whole dimensionless inverse, not theta or its increment
    ## separately (both can lose range). For a huge multiplier the -1
    ## in 2*exp(hazard)-1 is negligible.
    factor <- 1 + 2 * expm1(hazard)
    res[it] <- th_v[it] * factor
    big <- which(is.infinite(factor))
    res[it[big]] <- exp((log(th_v[it[big]]) + hazard[big]) + log(2))
  }

  out[vp] <- res
  out
}

#' Random generation from the continuous composite Pareto distribution
#'
#' @inheritParams dcomppareto
#' @param n number of observations; if \code{length(n) != 1}, the length is taken to
#'   be the number required (\code{numeric(0)} gives none), as in \code{\link[stats]{runif}}.
#'
#' @return Numeric vector of random draws.
#' @export
rcomppareto <- function(n, spec, alpha = 1, theta = 1, ...) {
  spec <- match.arg(spec, .comppareto_specs)
  if (length(n) != 1L) {
    n <- length(n)
  } else {
    if (is.na(n) || !is.finite(n) || n < 0) {
      stop("`n` must be a non-negative integer.")
    }
    n <- as.integer(floor(n))
  }
  if (n == 0L) return(numeric(0))
  body_args <- list(...)
  if (length(alpha) == 0L || length(theta) == 0L || any(lengths(body_args) == 0L)) {
    warning("NaNs produced", call. = FALSE)
    return(rep(NaN, n))
  }
  do.call(
    qcomppareto,
    c(
      list(
        p = stats::runif(n),
        spec = spec,
        alpha = rep_len(alpha, n),
        theta = rep_len(theta, n)
      ),
      .comppareto_recycle_args(n, body_args)
    )
  )
}
