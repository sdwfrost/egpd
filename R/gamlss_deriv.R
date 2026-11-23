## Analytic scores for the egpd gamlss families
##
## The score of an observation is the derivative of its log-density (or
## log-PMF) with respect to the natural gamlss parameters
##   mu = GPD scale sigma, sigma = GPD shape xi,
##   nu = kappa (carriers 1, 2, 3) or delta (carriers 4, 5),
##   tau = kappa (carrier 5) or pi (zero-inflated families).
## Scores are assembled from numerically stable closed-form pieces of the
## GPD and of the carrier G, all written in terms of A = -log(1 - F_GPD(y))
## so that u = F_GPD(y) and s = 1 - u are both available without
## cancellation.  For a discrete observation
##   d log P(y) / d theta = [dH(y + 1) / d theta - dH(y) / d theta] / P(y),
## with H = G(F_GPD(.)) and P(y) the stable PMF from ddiscegpd(); each term
## is evaluated as sign * exp(log|term| - log P(y)), so neither underflow
## nor the u^(kappa - 1) * du/dtheta = Inf * 0 limit at y = 0 (where
## F_GPD(0) = 0 and dH(0) / d theta = 0 exactly) can arise.  Small |xi| uses
## series for log1p(w) / w and its xi-derivative, so no switch to finite
## differences is needed near xi = 0.  For truncated-beta shape scores use
## the conditional cell and full normalizer derivatives directly.
##
## Observations with zero (or infinite) likelihood at the current parameters,
## e.g. counts beyond a finite upper endpoint, have no score: NaN is returned
## with a warning rather than a silent zero.  Invalid parameters give NaN
## with a warning and NA inputs give NA, following base R conventions.

.egpd_score_info <- function(key) {
  zi <- startsWith(key, "zi")
  base <- sub("^zi", "", key)
  disc <- startsWith(base, "degpd")
  model <- sub("^d?egpd", "", base)
  type <- switch(model, "1" = 1L, "3" = 4L, "4" = 5L, "5" = 2L, "6" = 3L,
                 stop("Unknown family: ", key, call. = FALSE))
  list(type = type, disc = disc, zi = zi,
       npar = if (type == 5L || zi) 4L else 3L)
}

## log1p(x) / x
.sc_lg1p_over_x <- function(x) {
  out <- log1p(x) / x
  sm <- which(abs(x) < 1e-4)
  if (length(sm)) {
    xs <- x[sm]
    out[sm] <- 1 - xs / 2 + xs^2 / 3 - xs^3 / 4
  }
  out
}

## (log1p(w) - w / (1 + w)) / w^2, so that dA/dxi = -z^2 * phi(w)
.sc_phi <- function(w) {
  out <- (log1p(w) - w / (1 + w)) / w^2
  sm <- which(abs(w) < 1e-3)
  if (length(sm)) {
    ws <- w[sm]
    out[sm] <- 0.5 - 2 * ws / 3 + 3 * ws^2 / 4 - 4 * ws^3 / 5 + 5 * ws^4 / 6
  }
  out
}

## x / expm1(x) - 1
.sc_xm1 <- function(x) {
  out <- x / expm1(x) - 1
  sm <- which(abs(x) < 1e-4)
  if (length(sm)) {
    xs <- x[sm]
    out[sm] <- -xs / 2 + xs^2 / 12 - xs^4 / 720
  }
  out[which(x == Inf)] <- -1
  out
}

## log(1 - (1 + x) exp(-x)) for x >= 0
.sc_log_q <- function(x) {
  out <- rep(0, length(x))
  fin <- which(is.finite(x))
  if (length(fin)) {
    xf <- x[fin]
    val <- log(-expm1(-xf) - xf * exp(-xf))
    sm <- which(xf < 1e-3)
    if (length(sm)) {
      xs <- xf[sm]
      val[sm] <- 2 * log(xs) - log(2) +
        log1p(-2 * xs / 3 + xs^2 / 4 - xs^3 / 15 + xs^4 / 72)
    }
    out[fin] <- val
  }
  out
}

## GPD pieces at t >= 0: A = -log S(t), dA/dscale, dA/dxi and the log-density
## derivatives.  Points at or beyond a finite upper endpoint have A = Inf.
.sc_gpd <- function(t, sc, xi) {
  n <- length(t)
  z <- t / sc
  w <- xi * z
  at0 <- t <= 0
  beyond <- !at0 & xi < 0 & w <= -1
  ok <- which(!at0 & !beyond)
  A <- dAs <- dAx <- dlfs <- dlfx <- numeric(n)
  A[beyond] <- Inf
  dlfs[at0] <- -1 / sc[at0]
  if (length(ok)) {
    zo <- z[ok]; wo <- w[ok]; so <- sc[ok]; xo <- xi[ok]
    A[ok] <- zo * .sc_lg1p_over_x(wo)
    dAs[ok] <- -(zo / so) / (1 + wo)
    ph <- .sc_phi(wo)
    dAx[ok] <- -zo^2 * ph
    dlfs[ok] <- -1 / so + (1 + xo) * (zo / so) / (1 + wo)
    dlfx[ok] <- -A[ok] + (1 + xo) * zo^2 * ph
  }
  list(A = A, dAs = dAs, dAx = dAx, dlfs = dlfs, dlfx = dlfx,
       at0 = at0, beyond = beyond)
}

## Carrier pieces at A (vectors of equal length).  Returns
##   lg   log g(u)
##   dlgs (d log g / du) * s
##   dlgk, dlgd  d log g / d kappa, d log g / d delta  (continuous scores)
##   Gk, Gd      list(sign, lm): dG(u)/d kappa = sign * exp(lm)  (discrete)
.sc_carrier <- function(type, A, k, d, need = c("lg", "dlgs", "dlgk", "dlgd", "Gk", "Gd")) {
  s <- exp(-A)
  u <- -expm1(-A)
  logu <- .egpd_log_u(A)
  out <- list()
  if (type == 1L) {
    if ("lg" %in% need) out$lg <- log(k) + ifelse(k == 1, 0, (k - 1) * logu)
    if ("dlgs" %in% need) out$dlgs <- ifelse(k == 1, 0, (k - 1) / expm1(A))
    if ("dlgk" %in% need) out$dlgk <- 1 / k + logu
    if ("Gk" %in% need) out$Gk <- list(sign = -1, lm = k * logu + .egpd_log_mlogu(A))
  } else if (type == 2L) {
    small <- k < 1e-2
    C <- .egpd_C2(k)
    if ("lg" %in% need) {
      lg <- 0.5 * log(k) - 0.5 * log(2 * pi) - 0.5 * k * s^2 - log(C)
      sm <- which(small)
      if (length(sm)) lg[sm] <- -0.5 * k[sm] * s[sm]^2 - log(.erf_series(k[sm]))
      out$lg <- lg
    }
    if ("dlgs" %in% need) out$dlgs <- k * s^2
    if ("dlgk" %in% need) {
      dser <- -1 / 6 + k / 20 - 3 * k^2 / 336 + 4 * k^3 / 3456 - 5 * k^4 / 42240
      out$dlgk <- ifelse(small, -s^2 / 2 - dser / .erf_series(k),
                         1 / (2 * k) - s^2 / 2 - dnorm(sqrt(k)) / (2 * sqrt(k) * C))
    }
    if ("Gk" %in% need) {
      sk <- sqrt(k)
      W <- (dnorm(sk * s) - dnorm(sk) * .egpd_S2_over_s(s, k)) / (2 * sk * C)
      sW <- sign(W)
      lW <- log(abs(W))
      ## W = phi(sk s) (1 - R) / (2 sk C) with R = exp(-D) S2(s) / s and
      ## D = kappa (1 - s^2) / 2.  dnorm() underflows for large kappa and the
      ## bracket cancels as u -> 0 (S2(s) / s rounds to 1), so evaluate it on
      ## the log scale there, with log(S2(s) / s) = log(1 - G(u)) + A for u < 1/2
      lg <- which(!is.finite(lW) | u < 1e-2)
      if (length(lg)) {
        Ag <- A[lg]; kg <- rep_len(k, length(A))[lg]; ug <- u[lg]
        lS2s <- numeric(length(lg))
        lo <- which(ug < 0.5)
        if (length(lo)) {
          lGg <- .egpd_carrier_log(Ag[lo], 2L, rep(NA_real_, length(lo)), kg[lo],
                                   rep(NA_real_, length(lo)))$lG
          lS2s[lo] <- .log1mexp(lGg) + Ag[lo]
        }
        hi <- which(ug >= 0.5)
        if (length(hi)) lS2s[hi] <- log(.egpd_S2_over_s(exp(-Ag[hi]), kg[hi]))
        lR <- lS2s - kg * (-expm1(-2 * Ag)) / 2
        l1mR <- ifelse(lR < 0, .log1mexp(lR),
                       ifelse(lR > 0, lR + .log1mexp(-lR), -Inf))
        lW[lg] <- -kg * exp(-2 * Ag) / 2 - 0.5 * log(2 * pi) + l1mR -
          log(2 * sqrt(kg) * rep_len(C, length(A))[lg])
        sW[lg] <- -sign(lR)
      }
      out$Gk <- list(sign = -sW, lm = -A + lW)
    }
  } else if (type == 3L) {
    z <- ((15/16) * s)^2
    if ("lg" %in% need) out$lg <- .egpd_tb_logdens(s, k)
    if ("dlgs" %in% need) out$dlgs <- (k - 1) * (2 * z / (1 - z))
    if ("dlgk" %in% need) out$dlgk <- log1p(-z) - .egpd_tb_norm(k)[, 2L]
  } else if (type == 4L) {
    lg4 <- log1p(1 / d) + log(-expm1(-d * A))
    if ("lg" %in% need) out$lg <- lg4
    if ("dlgs" %in% need) out$dlgs <- d / expm1(d * A)
    if ("dlgd" %in% need) out$dlgd <- .sc_xm1(d * A) / d + 1 / (d + 1)
    if ("Gd" %in% need) out$Gd <- list(sign = 1, lm = -A + .sc_log_q(d * A) - 2 * log(d))
  } else {
    na <- rep(NA_real_, length(A))
    cl4 <- .egpd_carrier_log(A, 4L, na, na, d)
    lG4 <- cl4$lG
    lmlG4 <- ifelse(lG4 <= -log(2), suppressWarnings(log(-lG4)), .egpd_log_mlog1m(cl4$lS))
    lg4 <- log1p(1 / d) + log(-expm1(-d * A))
    if ("lg" %in% need) out$lg <- log(k / 2) + ifelse(k == 2, 0, (k / 2 - 1) * lG4) + lg4
    if ("dlgs" %in% need) {
      out$dlgs <- ifelse(k == 2, 0, (k / 2 - 1) * exp(lg4 - A - lG4)) + d / expm1(d * A)
    }
    if ("dlgk" %in% need) out$dlgk <- 1 / k + lG4 / 2
    lq <- -A + .sc_log_q(d * A) - 2 * log(d)
    if ("dlgd" %in% need) {
      ## (dG4 / d delta) / G4 = exp(lq - lG4) tends to 1 / (delta + 1) as A -> 0
      ## (G4 ~ (1 + delta) A^2 / 2), where both logs are -Inf
      rq <- exp(lq - lG4)
      a0 <- which(A == 0)
      if (length(a0)) rq[a0] <- 1 / (rep_len(d, length(A))[a0] + 1)
      out$dlgd <- ifelse(k == 2, 0, (k / 2 - 1) * rq) +
        .sc_xm1(d * A) / d + 1 / (d + 1)
    }
    if ("Gk" %in% need) out$Gk <- list(sign = -1, lm = log(0.5) + (k / 2) * lG4 + lmlG4)
    if ("Gd" %in% need) out$Gd <- list(sign = 1, lm = log(k / 2) + (k / 2 - 1) * lG4 + lq)
  }
  out
}

## Map gamlss (nu, tau) to carrier (kappa, delta)
.sc_shapes <- function(type, nu, tau) {
  if (type %in% 1:3) list(k = nu, d = rep(NA_real_, length(nu)))
  else if (type == 4L) list(k = rep(NA_real_, length(nu)), d = nu)
  else list(k = tau, d = nu)
}

## Stable log-density / log-PMF of the base (non-zero-inflated) model
.sc_logdens <- function(y, sc, xi, k, d, type, disc) {
  if (disc) {
    ddiscegpd(y, sigma = sc, xi = xi, kappa = k, delta = d, type = type, log = TRUE)
  } else {
    degpd_density(y, sigma = sc, xi = xi, kappa = k, delta = d, type = type, log = TRUE)
  }
}

## Scores of the base continuous model for rows with a finite log-density
.sc_score_cont <- function(y, sc, xi, k, d, type, idx) {
  g <- .sc_gpd(y, sc, xi)
  t0 <- g$at0
  if (idx <= 2L) {
    cp <- .sc_carrier(type, g$A, k, d, need = "dlgs")
    dA <- if (idx == 1L) g$dAs else g$dAx
    carrier_term <- ifelse(t0, 0, cp$dlgs * dA)
    return(carrier_term + (if (idx == 1L) g$dlfs else g$dlfx))
  }
  if (type == 5L) {
    cp <- .sc_carrier(type, g$A, k, d, need = if (idx == 4L) "dlgk" else "dlgd")
    return(if (idx == 4L) cp$dlgk else cp$dlgd)
  }
  cp <- .sc_carrier(type, g$A, k, d, need = if (type == 4L) "dlgd" else "dlgk")
  if (type == 4L) cp$dlgd else cp$dlgk
}

## On a fixed interior cell, d log P is the mean continuous score weighted
## by its density. This avoids subtracting two nearly equal endpoint terms.
.sc_cell_score <- function(y, sc, xi, k, d, type, idx) {
  m <- length(y)
  rule <- .egpd_gl20()
  tt <- y + rep((1 + rule$x) / 2, each = m)
  rp <- function(v) rep(v, times = length(rule$x))
  At <- -.gpd_cdf_surv_std(tt, rp(sc), rp(xi))$log_surv
  ld <- matrix(.egpd_carrier_logdens(At, type, NA_real_, rp(k), rp(d)) +
                 .dgpd_std(tt, rp(sc), rp(xi), log = TRUE), nrow = m)
  sk <- matrix(.sc_score_cont(tt, rp(sc), rp(xi), rp(k), rp(d), type, idx), nrow = m)
  fine <- rowSums(!is.finite(ld) | !is.finite(sk)) == 0
  fine[fine] <- (.egpd_row_max(ld[fine, , drop = FALSE]) -
                   .egpd_row_min(ld[fine, , drop = FALSE])) <= 3
  out <- rep(NA_real_, m)
  if (any(fine)) {
    lw <- ld[fine, , drop = FALSE] + rep(log(rule$w), each = sum(fine))
    w <- exp(lw - .egpd_row_max(lw))
    out[fine] <- rowSums(w * sk[fine, , drop = FALSE]) / rowSums(w)
  }
  out
}

## Scores of the base discrete model for rows with log P(y) = logP finite
.sc_score_disc <- function(y, sc, xi, k, d, type, idx, logP) {
  n <- length(y)
  g0 <- .sc_gpd(y, sc, xi)
  g1 <- .sc_gpd(y + 1, sc, xi)
  live0 <- !g0$at0 & !g0$beyond
  live1 <- !g1$at0 & !g1$beyond
  if (type == 3L && idx == 3L)
    return(egpd_tb_moments_cpp(k, exp(-g1$A), exp(-g0$A))[, 2L] - .egpd_tb_norm(k)[, 2L])
  difference <- function(i, lm0, lm1, sign0, sign1) {
    l0 <- ifelse(live0[i] & sign0 != 0, lm0, -Inf)
    l1 <- ifelse(live1[i] & sign1 != 0, lm1, -Inf)
    out <- sign1 * exp(l1 - logP[i]) - sign0 * exp(l0 - logP[i])
    narrow <- which(live0[i] & live1[i] & (!is.finite(out) |
                    (sign0 == sign1 & abs(l0 - l1) <= 1e-5 * pmax(1, abs(l0), abs(l1)))))
    if (length(narrow)) {
      j <- i[narrow]
      avg <- .sc_cell_score(y[j], sc[j], xi[j], k[j], d[j], type, idx)
      ok <- is.finite(avg)
      out[narrow[ok]] <- avg[ok]
    }
    out
  }
  carrier_param <- idx >= 3L
  if (!carrier_param) {
    A0 <- ifelse(live0, g0$A, 1)
    A1 <- ifelse(live1, g1$A, 1)
    lg0 <- .sc_carrier(type, A0, k, d, need = "lg")$lg
    lg1 <- .sc_carrier(type, A1, k, d, need = "lg")$lg
    dA0 <- if (idx == 1L) g0$dAs else g0$dAx
    dA1 <- if (idx == 1L) g1$dAs else g1$dAx
    return(difference(seq_len(n), lg0 - A0 + log(-dA0), lg1 - A1 + log(-dA1), -1, -1))
  }
  which_shape <- if (type == 5L) (if (idx == 4L) "k" else "d") else if (type == 4L) "d" else "k"
  num <- if (type == 2L) k < 1e-3 else rep(FALSE, n)
  out <- numeric(n)
  an <- which(!num)
  if (length(an)) {
    nm <- if (which_shape == "k") "Gk" else "Gd"
    A0 <- ifelse(live0, g0$A, 1)[an]
    A1 <- ifelse(live1, g1$A, 1)[an]
    c0 <- .sc_carrier(type, A0, k[an], d[an], need = nm)[[nm]]
    c1 <- .sc_carrier(type, A1, k[an], d[an], need = nm)[[nm]]
    s0 <- rep_len(c0$sign, length(an)); s1 <- rep_len(c1$sign, length(an))
    out[an] <- difference(an, c0$lm, c1$lm, s0, s1)
  }
  nn <- which(num)
  if (length(nn)) {
    ## d log P / d kappa = int_y^{y+1} f d(log f)/d kappa / int f: where the
    ## density varies slowly over [y, y + 1] this weighted mean of the
    ## continuous carrier-shape score is evaluated by Gauss-Legendre, which
    ## avoids differentiating the far-tail log P numerically
    out[nn] <- NA_real_
    i <- nn[y[nn] >= 1 & live1[nn]]
    if (length(i))
      out[i] <- .sc_cell_score(y[i], sc[i], xi[i], k[i], d[i], type, idx)
    rest <- nn[!is.finite(out[nn])]
    if (length(rest)) {
      yn <- y[rest]; sn <- sc[rest]; xn <- xi[rest]; dn <- d[rest]; kn <- k[rest]
      ldf <- function(p, j) .sc_logdens(yn[j], sn[j], xn[j], p, dn[j], type, TRUE)
      out[rest] <- .egpd_nd_core(ldf, kn, pmax(1e-5 * kn, 1e-8))
    }
  }
  out
}

## Score engine shared by all EGPD/DEGPD/ZIEGPD/ZIDEGPD gamlss families
.egpd_score <- function(y, mu, sigma, nu, tau, key, idx) {
  info <- .egpd_score_info(key)
  if (!(idx %in% seq_len(info$npar))) stop("Unknown family/idx: ", key, "_", idx, call. = FALSE)
  if (is.null(tau)) tau <- numeric(1)
  lens <- c(length(y), length(mu), length(sigma), length(nu), if (info$npar == 4L) length(tau))
  if (min(lens) == 0L) return(numeric(0))
  n <- max(lens)
  y <- rep_len(as.numeric(y), n); mu <- rep_len(as.numeric(mu), n)
  sigma <- rep_len(as.numeric(sigma), n); nu <- rep_len(as.numeric(nu), n)
  tau <- rep_len(as.numeric(tau), n)
  type <- info$type
  sh <- .sc_shapes(type, nu, tau)
  k <- sh$k; d <- sh$d
  pi0 <- if (info$zi) tau else rep(0, n)

  active <- list(y, mu, sigma, nu)
  if (info$npar == 4L) active <- c(active, list(tau))
  has_na <- Reduce(`|`, lapply(active, is.na))
  shape_ok <- if (type %in% 1:3) is.finite(k) & (k > 0 | (type == 2L & k == 0)) else TRUE
  if (type %in% 4:5) shape_ok <- shape_ok & is.finite(d) & d > 0
  if (type == 5L) shape_ok <- shape_ok & is.finite(k) & k > 0
  par_ok <- !has_na & is.finite(mu) & mu > 0 & is.finite(sigma) & shape_ok &
    is.finite(pi0) & pi0 >= 0 & pi0 <= 1
  out <- rep(NaN, n)
  bad_par <- sum(!has_na & !par_ok)

  ok <- which(par_ok)
  undefined <- 0L
  if (length(ok)) {
    yo <- y[ok]; so <- mu[ok]; xo <- sigma[ok]; ko <- k[ok]; dd <- d[ok]; po <- pi0[ok]
    res <- rep(NaN, length(ok))
    zero_obs <- info$zi & yo == 0
    if (info$zi && idx == 4L) {
      if (info$disc) {
        lp0 <- .sc_logdens(rep(0, length(ok)), so, xo, ko, dd, type, TRUE)
        res[zero_obs] <- (-expm1(lp0) / (po + (1 - po) * exp(lp0)))[zero_obs]
      } else {
        res[zero_obs] <- 1 / po[zero_obs]
      }
      pos <- which(!zero_obs)
      if (length(pos)) {
        lp <- .sc_logdens(yo[pos], so[pos], xo[pos], ko[pos], dd[pos], type, info$disc)
        res[pos] <- ifelse(is.finite(lp) & po[pos] < 1, -1 / (1 - po[pos]), NaN)
      }
    } else {
      ## base-model score at y (or at 0 for zero-inflated zeros)
      lp <- .sc_logdens(yo, so, xo, ko, dd, type, info$disc)
      fin <- if (info$zi && !info$disc) !zero_obs & is.finite(lp) else is.finite(lp)
      if (info$zi) fin <- fin & (zero_obs | po < 1)
      f <- which(fin)
      if (length(f)) {
        res[f] <- if (info$disc) {
          .sc_score_disc(yo[f], so[f], xo[f], ko[f], dd[f], type, idx, lp[f])
        } else {
          .sc_score_cont(yo[f], so[f], xo[f], ko[f], dd[f], type, idx)
        }
      }
      if (info$zi) {
        z <- which(zero_obs)
        if (length(z)) {
          if (info$disc) {
            ## d log(pi + (1 - pi) P0) = w0 * d log P0, w0 = (1 - pi) P0 / f(0)
            w0 <- plogis(log1p(-po[z]) + lp[z] - log(po[z]))
            res[z] <- ifelse(w0 == 0, 0, w0 * res[z])
          } else {
            res[z] <- 0
          }
        }
      }
    }
    res[!is.finite(res)] <- NaN
    undefined <- sum(is.nan(res))
    out[ok] <- res
  }
  if (bad_par) warning("NaNs produced", call. = FALSE)
  if (undefined) {
    warning(sprintf(paste0(
      "score undefined for %d observation(s) with zero or infinite likelihood ",
      "at the current parameters (e.g. outside the support); NaN returned"),
      undefined), call. = FALSE)
  }
  .egpd_apply_na_nan(out, active)
}

## ---------------------------------------------------------------
## Dispatchers used by the gamlss family objects
## ---------------------------------------------------------------

#' @keywords internal
.egpd_ad3 <- function(y, mu, sigma, nu, family, idx) {
  .egpd_score(y, mu, sigma, nu, NULL, family, idx)
}

#' @keywords internal
.egpd_ad4 <- function(y, mu, sigma, nu, tau, family, idx) {
  .egpd_score(y, mu, sigma, nu, tau, family, idx)
}

## ---------------------------------------------------------------
## Override helpers: plug analytical derivatives into a gamlss family
## ---------------------------------------------------------------

## 3-parameter family override
.override_derivs3 <- function(fam, family_key) {
  fam$dldm <- eval(bquote(function(y, mu, sigma, nu, ...) {
    egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 1L)
  }))
  fam$dldd <- eval(bquote(function(y, mu, sigma, nu, ...) {
    egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 2L)
  }))
  fam$dldv <- eval(bquote(function(y, mu, sigma, nu, ...) {
    egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 3L)
  }))

  fam$d2ldm2 <- eval(bquote(function(y, mu, sigma, nu, ...) {
    dl <- egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 1L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldd2 <- eval(bquote(function(y, mu, sigma, nu, ...) {
    dl <- egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 2L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldv2 <- eval(bquote(function(y, mu, sigma, nu, ...) {
    dl <- egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 3L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))

  fam$d2ldmdd <- eval(bquote(function(y, mu, sigma, nu, ...) {
    -egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 1L) *
      egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 2L)
  }))
  fam$d2ldmdv <- eval(bquote(function(y, mu, sigma, nu, ...) {
    -egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 1L) *
      egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 3L)
  }))
  fam$d2ldddv <- eval(bquote(function(y, mu, sigma, nu, ...) {
    -egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 2L) *
      egpd::.egpd_ad3(y, mu, sigma, nu, .(family_key), idx = 3L)
  }))

  fam
}

## 4-parameter family override
.override_derivs4 <- function(fam, family_key) {
  fam$dldm <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L)
  }))
  fam$dldd <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L)
  }))
  fam$dldv <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L)
  }))
  fam$dldt <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
  }))

  fam$d2ldm2 <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    dl <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldd2 <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    dl <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldv2 <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    dl <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))
  fam$d2ldt2 <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    dl <- egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
    d2 <- -dl * dl; ifelse(d2 < -1e-15, d2, -1e-15)
  }))

  fam$d2ldmdd <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L)
  }))
  fam$d2ldmdv <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L)
  }))
  fam$d2ldmdt <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 1L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
  }))
  fam$d2ldddv <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L)
  }))
  fam$d2ldddt <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 2L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
  }))
  fam$d2ldvdt <- eval(bquote(function(y, mu, sigma, nu, tau, ...) {
    -egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 3L) *
      egpd::.egpd_ad4(y, mu, sigma, nu, tau, .(family_key), idx = 4L)
  }))

  fam
}
