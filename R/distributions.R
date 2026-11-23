## Distribution functions for EGPD, Discrete EGPD, and Zero-Inflated Discrete EGPD
## Adapted from degpd_functions.R in the degpd-and-zidegpd prototype

## Self-contained GPD CDF/PDF/quantile (replaces extraDistr dependency)

.egpd_restore_dim <- function(out, template) {
  if (!is.null(dim(template)) && length(out) == length(template)) {
    dim(out) <- dim(template)
    dimnames(out) <- dimnames(template)
  } else if (!is.null(names(template)) && length(out) == length(template)) {
    names(out) <- names(template)
  }
  out
}

.egpd_apply_na_nan <- function(out, active) {
  if (!length(out) || !length(active)) return(out)
  n <- length(out)
  has_nan <- rep(FALSE, n)
  has_na  <- rep(FALSE, n)
  for (v in active) {
    vn <- rep_len(as.numeric(v), n)
    inan <- is.nan(vn)
    has_nan <- has_nan | inan
    has_na  <- has_na  | (is.na(vn) & !inan)
  }
  has_na <- has_na & !has_nan
  out[has_nan] <- NaN
  out[has_na]  <- NA_real_
  out
}

.egpd_active_carrier_args <- function(type, prob, kappa, delta) {
  if (type == 0) list()
  else if (type %in% c(1, 2, 3)) list(kappa)
  else if (type == 4) list(delta)
  else if (type == 5) list(kappa, delta)
  else list(prob, kappa, delta)
}

.egpd_carrier_valid <- function(type, prob, kappa, delta, n) {
  ok <- rep(TRUE, n)
  if (type %in% c(1, 3, 5, 6)) {
    ok <- ok & !is.na(kappa) & is.finite(kappa) & kappa > 0
  } else if (type == 2) {
    ok <- ok & !is.na(kappa) & is.finite(kappa) & kappa >= 0
  }
  if (type %in% c(4, 5, 6)) {
    ok <- ok & !is.na(delta) & is.finite(delta) & delta > 0
  }
  if (type == 6) {
    ok <- ok & !is.na(prob) & is.finite(prob) & prob >= 0 & prob <= 1
  }
  ok
}

## Parameters that `type` (and the GPD) need but that were not supplied are an
## error; supplied NA values propagate as NA results instead.
.egpd_missing_args <- function(nms, env = parent.frame()) {
  nms[vapply(nms, function(nm) eval(call("missing", as.name(nm)), env), logical(1))]
}

.egpd_check_required <- function(type, missing_args, gpd = TRUE, pi = FALSE) {
  need <- c(if (gpd) c("sigma", "xi"),
            if (type %in% c(1, 2, 3, 5, 6)) "kappa",
            if (type %in% c(4, 5, 6)) "delta",
            if (type == 6) "prob",
            if (pi) "pi")
  miss <- need[need %in% missing_args]
  if (length(miss) == 1L) {
    stop(sprintf("Argument `%s' missing.", miss), call. = FALSE)
  } else if (length(miss) > 1L) {
    stop(sprintf("Arguments %s missing.", paste0("`", miss, "'", collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

## Number of draws for r*() functions, following base R: length(n) when
## length(n) != 1, otherwise n truncated to an integer (an error if invalid).
.egpd_rn <- function(n) {
  if (length(n) != 1L) return(length(n))
  nn <- suppressWarnings(as.numeric(n))
  if (is.na(nn) || nn < 0 || !is.finite(nn)) stop("invalid arguments", call. = FALSE)
  as.integer(floor(nn))
}

.egpd_parse_prob <- function(p, lower.tail = TRUE, log.p = FALSE) {
  ## Returns both tail probabilities on the linear and the log scale, so that
  ## log-scale inputs far in either tail are never rounded to 0 or 1.
  p <- as.numeric(p)
  if (log.p) {
    valid <- !is.na(p) & p <= 0
    lp <- ifelse(valid, p, NaN)
    if (lower.tail) {
      lp_low <- lp
      lp_upp <- .log1mexp(lp)
    } else {
      lp_upp <- lp
      lp_low <- .log1mexp(lp)
    }
    p_low <- exp(lp_low)
    p_upp <- exp(lp_upp)
  } else {
    valid <- !is.na(p) & p >= 0 & p <= 1
    if (lower.tail) {
      p_low <- p
      p_upp <- 1 - p
    } else {
      p_upp <- p
      p_low <- 1 - p
    }
    lp_low <- suppressWarnings(log(p_low))
    lp_upp <- suppressWarnings(log(p_upp))
  }
  list(p_low = p_low, p_upp = p_upp, lp_low = lp_low, lp_upp = lp_upp, valid = valid)
}

## Zero-inflated quantiles: le_pi flags lower-tail probabilities at or below pi
## (quantile zero).  A log.p probability below the normal double range, whose
## linear value is subnormal or 0, is compared with pi on the log scale and
## flagged `tiny`, so that its zero mass is also removed on the log scale.
.egpd_zi_split <- function(pr, piv, valid, log.p) {
  le_pi <- pr$p_low <= piv
  tiny <- rep(FALSE, length(piv))
  if (log.p) {
    tiny <- valid & pr$p_low < .Machine$double.xmin
    it <- which(tiny)
    if (length(it)) le_pi[it] <- pr$lp_low[it] <= log(piv[it])
  }
  list(le_pi = le_pi, tiny = tiny)
}

## log((p - pi) / (1 - pi)) from lp = log(p) > log(pi)
.egpd_zi_lp_adj <- function(lp, piv) {
  lp + .log1mexp(log(piv) - lp) - log1p(-piv)
}

## Zero-inflated log CDF, lcdf = log(pi + (1 - pi) F), from the base log CDF lb. Where
## the zero-inflated survivor (1 - pi) S is below 1e-3, the log-sum keeps only absolute
## accuracy (~1e-17); log1p(-(1 - pi) S), with S = -expm1(lb), keeps the relative
## accuracy of lb there. Rows with pi = 0 keep lb itself.
.egpd_zi_lcdf_near1 <- function(lcdf, lb, piv) {
  sz <- (1 - piv) * -expm1(lb)
  nr <- which(sz < 1e-3 & piv > 0)
  lcdf[nr] <- log1p(-sz[nr])
  lcdf
}

.erf_series <- function(z) {
  ## sum_{j=0}^5 (-z)^j / (2^j j! (2j+1)) for small |z|
  1 - z / 6 + z^2 / 40 - z^3 / 336 + z^4 / 3456 - z^5 / 42240
}

.expm1_over_x <- function(x) {
  ## (-expm1(-x)) / x for x >= 0
  out <- numeric(length(x))
  small <- !is.na(x) & x < 1e-2
  if (any(small)) {
    xs <- x[small]
    out[small] <- 1 - xs / 2 + xs^2 / 6 - xs^3 / 24 + xs^4 / 120 - xs^5 / 720
  }
  large <- !is.na(x) & !small
  if (any(large)) {
    xl <- x[large]
    out[large] <- ifelse(is.infinite(xl), 0, -expm1(-xl) / xl)
  }
  out
}

.m1_expm1_over_x <- function(x) {
  ## 1 - (-expm1(-x)) / x = (x + expm1(-x)) / x for x >= 0
  out <- numeric(length(x))
  small <- !is.na(x) & x < 1e-2
  if (any(small)) {
    xs <- x[small]
    out[small] <- xs / 2 - xs^2 / 6 + xs^3 / 24 - xs^4 / 120 + xs^5 / 720 - xs^6 / 5040
  }
  large <- !is.na(x) & !small
  if (any(large)) {
    xl <- x[large]
    out[large] <- ifelse(is.infinite(xl), 1, (xl + expm1(-xl)) / xl)
  }
  out
}

.pG4_inner <- function(u, delta) {
  ## Exact closed form G4(u; delta) = 1 - (1 - u) * (1 + (-expm1(delta * log(1 - u))) / delta)
  ## for u in (0, 1) and delta > 0, cancellation-free for both small u and large delta.
  out <- numeric(length(u))
  lo <- u < 0.5
  if (any(lo)) {
    ul <- u[lo]; dl <- delta[lo]
    L <- -log1p(-ul)
    dL <- dl * L
    b <- L * .expm1_over_x(dL)
    out[lo] <- L * (.m1_expm1_over_x(dL) - .m1_expm1_over_x(L)) + ul * b
  }
  hi <- !lo
  if (any(hi)) {
    uh <- u[hi]; dh <- delta[hi]
    L <- -log1p(-uh)
    b <- L * .expm1_over_x(dh * L)
    out[hi] <- 1 - (1 - uh) * (1 + b)
  }
  pmin(pmax(out, 0), 1)
}

.sG4_inner <- function(s, delta) {
  ## Exact survival 1 - G4(1 - s; delta) for s in (0, 1) and delta > 0.
  L <- -log(s)
  pmin(pmax(s * (1 + L * .expm1_over_x(delta * L)), 0), 1)
}

.qG4_pair <- function(p_low, p_upp, delta) {
  ## Returns list(u = Q4(p_low; delta), s = 1 - Q4(p_low; delta)) for valid
  ## p_low, p_upp in [0, 1] (p_low + p_upp = 1) and delta > 0.
  n <- length(p_low)
  u_out <- numeric(n)
  s_out <- numeric(n)
  at0 <- p_low == 0
  at1 <- p_upp == 0
  u_out[at0] <- 0; s_out[at0] <- 1
  u_out[at1] <- 1; s_out[at1] <- 0

  mid <- which(!at0 & !at1)
  if (!length(mid)) return(list(u = u_out, s = s_out))

  ## A = -log(s) has the law E1 + E2/(1+delta), with independent unit
  ## exponentials. Its quantile is between Q_Gamma(2)/(1+delta) and
  ## Q_Gamma(2): below half an ulp the limit is accurate to working precision.
  tiny <- mid[delta[mid] <= .Machine$double.eps / 2]
  if (length(tiny)) {
    lower <- p_low[tiny] <= 0.5
    A <- numeric(length(tiny))
    A[lower] <- qgamma(p_low[tiny[lower]], shape = 2)
    A[!lower] <- qgamma(p_upp[tiny[!lower]], shape = 2, lower.tail = FALSE)
    u_out[tiny] <- -expm1(-A)
    s_out[tiny] <- exp(-A)
    mid <- setdiff(mid, tiny)
  }

  ## B = s^delta ~ Beta(1/delta,2). Obtain the smaller of B and 1-B
  ## directly: a rounded near-one B raised to 1/delta can give zero,
  ## infinity or a badly wrong median. P(B <= .5) = 2^-a (1+a/2).
  fast <- mid[delta[mid] <= 50 & p_low[mid] >= 1e-6 & p_upp[mid] >= 1e-6]
  if (length(fast)) {
    a <- 1 / delta[fast]
    left <- log(p_upp[fast]) <= log1p(a / 2) - a * log(2)
    ls <- numeric(length(fast))
    ls[left] <- a[left] * log(qbeta(p_upp[fast[left]], a[left], 2))
    ls[!left] <- a[!left] * log1p(-qbeta(p_low[fast[!left]], 2, a[!left]))
    s_out[fast] <- exp(ls)
    u_out[fast] <- -expm1(ls)
  }

  ## Bracketed Newton / bisection for large delta or extreme tails
  rem <- setdiff(mid, fast)
  for (i in rem) {
    pl <- p_low[i]; pu <- p_upp[i]; d <- delta[i]
    if (pl < 0.5) {
      lo <- pl
      hi <- min(1, pl + pu / (1 + d))
      v <- sqrt(2 * pl) / sqrt(1 + d)
      if (!is.finite(v) || v <= lo || v >= hi) v <- 0.5 * (lo + hi)
      for (it in 1:60) {
        fv <- .pG4_inner(v, d) - pl
        if (abs(fv) <= 1e-15 * pl || (hi - lo) <= 1e-15 * v) break
        if (fv < 0) lo <- v else hi <- v
        dv <- (1 + 1 / d) * (-expm1(d * log1p(-v)))
        vn <- v - fv / dv
        v <- if (is.finite(vn) && vn > lo && vn < hi) vn else 0.5 * (lo + hi)
      }
      u_out[i] <- v
      s_out[i] <- 1 - v
    } else {
      lo <- (d / (1 + d)) * pu
      hi <- pu
      s <- 0.5 * (lo + hi)
      for (it in 1:60) {
        fs <- .sG4_inner(s, d) - pu
        if (abs(fs) <= 1e-15 * pu || (hi - lo) <= 1e-15 * s) break
        if (fs < 0) lo <- s else hi <- s
        ds <- (1 + 1 / d) * (-expm1(d * log(s)))
        sn <- s - fs / ds
        s <- if (is.finite(sn) && sn > lo && sn < hi) sn else 0.5 * (lo + hi)
      }
      s_out[i] <- s
      u_out[i] <- 1 - s
    }
  }
  list(u = u_out, s = s_out)
}

.pG2_inner <- function(u, kappa) {
  ## Truncated normal carrier CDF on u in (0, 1), kappa >= 0.
  out <- numeric(length(u))
  small <- kappa < 1e-2
  if (any(small)) {
    us <- u[small]; ks <- kappa[small]
    w <- 1 - us; w2 <- w * w
    dS <- -ks * (1 - w2) / 6 + ks^2 * (1 - w2^2) / 40 -
          ks^3 * (1 - w2^3) / 336 + ks^4 * (1 - w2^4) / 3456 -
          ks^5 * (1 - w2^5) / 42240
    out[small] <- us + w * (dS / .erf_series(ks))
  }
  large <- !small
  if (any(large)) {
    ul <- u[large]; kl <- kappa[large]
    sk <- sqrt(kl)
    fmin <- pnorm(-sk)
    out[large] <- (pnorm(sk * (ul - 1)) - fmin) / (0.5 - fmin)
  }
  pmin(pmax(out, 0), 1)
}

.log_dG2_inner <- function(u, kappa) {
  ## Log density log g2(u; kappa) on u in [0, 1], kappa >= 0.
  out <- numeric(length(u))
  small <- kappa < 1e-2
  if (any(small)) {
    us <- u[small]; ks <- kappa[small]
    out[small] <- -0.5 * ks * (1 - us)^2 -
      log1p(-ks / 6 + ks^2 / 40 - ks^3 / 336 + ks^4 / 3456 - ks^5 / 42240)
  }
  large <- !small
  if (any(large)) {
    ul <- u[large]; kl <- kappa[large]
    sk <- sqrt(kl)
    out[large] <- 0.5 * log(kl) - 0.5 * log(2 * pi) - 0.5 * kl * (ul - 1)^2 -
      log(0.5 - pnorm(-sk))
  }
  out
}

.qG2_pair <- function(p_low, p_upp, kappa) {
  ## Returns list(u = Q2(p_low; kappa), s = 1 - Q2(p_low; kappa)) for valid inputs.
  n <- length(p_low)
  u_out <- numeric(n)
  s_out <- numeric(n)
  at0 <- p_low == 0
  at1 <- p_upp == 0
  u_out[at0] <- 0; s_out[at0] <- 1
  u_out[at1] <- 1; s_out[at1] <- 0

  mid <- which(!at0 & !at1)
  if (!length(mid)) return(list(u = u_out, s = s_out))

  small <- mid[kappa[mid] < 1e-3]
  if (length(small)) {
    pl <- p_low[small]; y <- p_upp[small]; k <- kappa[small]
    y2 <- y * y
    poly <- k / 6 + (7 * y2 - 3) * k^2 / 120 +
            (127 * y2^2 - 118 * y2 + 15) * k^3 / 5040 +
            (4369 * y2^3 - 6299 * y2^2 + 2227 * y2 - 105) * k^4 / 362880
    corr <- y * (1 - y2) * poly
    u_out[small] <- pl + corr
    s_out[small] <- y * (1 - (1 - y2) * poly)
  }

  large <- setdiff(mid, small)
  if (length(large)) {
    pl <- p_low[large]; pu <- p_upp[large]; k <- kappa[large]
    sk <- sqrt(k)
    fmin <- pnorm(-sk)
    c_norm <- 0.5 - fmin
    eps <- pu * c_norm
    s_val <- ifelse(
      eps < 1e-5,
      (sqrt(2 * pi) * eps * (1 + (pi / 3) * eps^2 + (7 * pi^2 / 30) * eps^4)) / sk,
      -qnorm(ifelse(pl < 0.5, fmin + pl * c_norm, 0.5 - eps)) / sk
    )
    s_out[large] <- s_val
    u_out[large] <- 1 - s_val
  }
  list(u = pmin(pmax(u_out, 0), 1), s = pmin(pmax(s_out, 0), 1))
}

.qG3_pair <- function(p_low, p_upp, kappa) {
  lower <- 1 / 32
  upper <- 0.5
  hml <- upper - lower
  n <- length(p_low)
  u_out <- numeric(n); s_out <- numeric(n)
  at0 <- p_low == 0; at1 <- p_upp == 0
  u_out[at0] <- 0; s_out[at0] <- 1
  u_out[at1] <- 1; s_out[at1] <- 0
  mid <- which(!at0 & !at1)
  if (!length(mid)) return(list(u = u_out, s = s_out))

  k <- kappa[mid]; pl <- p_low[mid]; pu <- p_upp[mid]
  plo <- pbeta(lower, k, k)
  phi <- pbeta(upper, k, k)
  cn <- phi - plo
  ## the difference cancels as kappa -> 0 (cn is about 1.7 kappa): use the
  ## Beta(1/2, kappa) form of the normaliser there
  ksm <- which(k < 1e-2)
  if (length(ksm)) cn[ksm] <- exp(.egpd_logZ3(k[ksm]))
  delta_p <- pu * cn
  d05 <- dbeta(0.5, k, k)
  ## The expansion parameter is (k - 1) * t^2, not t^2 alone.
  tiny <- delta_p / d05 < 1e-6 * hml / sqrt(pmax(1, k))
  if (any(tiny)) {
    idx <- mid[tiny]
    kt <- k[tiny]
    t1 <- delta_p[tiny] / d05[tiny]
    ## Invert t * (1 - 4*(k-1)/3 * t^2) = t1 => t = t1 * (1 + 4*(k-1)/3 * t1^2)
    t_val <- t1 * (1 + 4 * (kt - 1) * t1^2 / 3)
    st <- t_val / hml
    s_out[idx] <- st
    u_out[idx] <- 1 - st
  }
  if (any(!tiny)) {
    idx <- mid[!tiny]
    kn <- k[!tiny]
    lo_half <- pl[!tiny] <= 0.5
    ## Near the truncation point invert in z. Else invert in (1-2*z)^2:
    ## for concentrated betas even an interior z rounds close to 1/2.
    near_lo <- lo_half & plo[!tiny] >= pl[!tiny] * cn[!tiny]
    il <- which(near_lo)
    if (length(il)) {
      qb <- qbeta(plo[!tiny][il] + pl[!tiny][il] * cn[!tiny][il], kn[il], kn[il])
      u_out[idx[il]] <- (qb - lower) / hml
      s_out[idx[il]] <- (upper - qb) / hml
    }
    iu <- which(!near_lo)
    if (length(iu)) {
      ## P(1/2 - t < X < 1/2) = I_{4 t^2}(1/2, k) / 2 for X ~ Beta(k, k), so the
      ## upper tail is inverted without subtracting from pbeta(1/2)
      w <- numeric(length(iu))
      left <- lo_half[iu]
      if (any(left)) {
        j <- iu[left]
        w[left] <- qbeta(2 * (plo[!tiny][j] + pl[!tiny][j] * cn[!tiny][j]),
                         0.5, kn[j], lower.tail = FALSE)
      }
      if (any(!left)) {
        j <- iu[!left]
        w[!left] <- qbeta(2 * delta_p[!tiny][j], 0.5, kn[j])
      }
      t <- sqrt(w) / 2
      s_out[idx[iu]] <- t / hml
      u_out[idx[iu]] <- 1 - t / hml
    }
  }
  list(u = pmin(pmax(u_out, 0), 1), s = pmin(pmax(s_out, 0), 1))
}

## Log survival at finite points inside the GPD support. When xi*z is
## tiny, dividing a rounded (possibly zero) product by xi loses z.
.gpd_logsurv_inside <- function(q, sigma, xi) {
  z <- q / sigma
  t <- xi * z
  ls <- -z
  reg <- which(is.finite(t) & t > -.5 & abs(t) > .Machine$double.eps)
  ls[reg] <- -log1p(t[reg]) / xi[reg]
  edge <- which(t <= -.5)
  if (length(edge)) {
    gap <- gpd_endpoint_gap_cpp(q[edge], sigma[edge], xi[edge])
    ls[edge] <- -log(pmax(0, gap)) / xi[edge]
  }
  ov <- which(xi > 0 & is.infinite(t))
  if (length(ov)) {
    lt <- log(xi[ov]) + log(q[ov]) - log(sigma[ov])
    ls[ov] <- -(pmax(lt, 0) + log1p(exp(-abs(lt)))) / xi[ov]
  }
  ls
}

## Largest count in the support of floor(X) for strictly negative shapes.
.gpd_last_count <- function(s, xi) {
  end <- -s / xi
  k <- pmax(0, ceiling(end) - 1)
  j <- which(is.finite(end) & end == floor(end) & end < 2^53)
  if (length(j)) {
    i <- j[gpd_endpoint_gap_cpp(end[j], s[j], xi[j]) > 0]
    k[i] <- end[i]
  }
  k
}

## A subnormal base CDF can become an ordinary probability under a power
## carrier. Retain log(u) (u equals the cumulative hazard to this precision).
.gpd_logu_small <- function(q, s, xi) {
  lu <- log(q) - log(s)
  t <- sign(xi)*exp(log(abs(xi)) + lu)
  j <- which(abs(t) > .Machine$double.eps)
  v <- log1p(pmax(-1, t[j]))
  edge <- which(t[j] < -.5)
  if (length(edge)) {
    i <- j[edge]
    v[edge] <- log(gpd_endpoint_gap_cpp(q[i], s[i], xi[i]))
  }
  lu[j] <- log(abs(v)) - log(abs(xi[j]))
  lu
}

.gpd_from_small_logu <- function(lu, s, xi) {
  t <- sign(xi)*exp(log(abs(xi)) + lu)
  adj <- numeric(length(t))
  j <- which(abs(t) > .Machine$double.eps)
  adj[j] <- log(expm1(t[j])/t[j])
  exp(log(s) + lu + adj)
}

.gpd_cdf_surv_std <- function(q, sigma, xi) {
  ## Returns list(cdf, surv, log_surv) for standard GPD(0, sigma, xi).
  if (min(length(q), length(sigma), length(xi)) == 0L) {
    return(list(cdf = numeric(0), surv = numeric(0), log_surv = numeric(0)))
  }
  n <- max(length(q), length(sigma), length(xi))
  qv <- rep_len(as.numeric(q), n)
  sv <- rep_len(as.numeric(sigma), n)
  xv <- rep_len(as.numeric(xi), n)

  cdf <- rep(NaN, n)
  surv <- rep(NaN, n)
  lsurv <- rep(NaN, n)

  valid <- !is.na(qv) & !is.na(sv) & !is.na(xv) & is.finite(sv) & sv > 0 & is.finite(xv)
  below <- which(valid & qv <= 0)
  cdf[below] <- 0; surv[below] <- 1; lsurv[below] <- 0

  above <- which(valid & qv == Inf)
  cdf[above] <- 1; surv[above] <- 0; lsurv[above] <- -Inf

  inside <- which(valid & qv > 0 & is.finite(qv))
  if (length(inside)) {
    ls <- .gpd_logsurv_inside(qv[inside], sv[inside], xv[inside])
    lsurv[inside] <- ls
    cdf[inside]   <- -expm1(ls)
    surv[inside]  <- exp(ls)
  }

  active <- list(q, sigma, xi)
  list(
    cdf      = .egpd_apply_na_nan(cdf, active),
    surv     = .egpd_apply_na_nan(surv, active),
    log_surv = .egpd_apply_na_nan(lsurv, active)
  )
}

.pgpd_std <- function(q, sigma, xi) {
  ## GPD CDF with location=0, scale=sigma, shape=xi
  if (min(length(q), length(sigma), length(xi)) == 0L) return(numeric(0))
  template <- q
  res <- .gpd_cdf_surv_std(q, sigma, xi)$cdf
  .egpd_restore_dim(res, template)
}

.dgpd_std <- function(x, sigma, xi, log = FALSE) {
  ## GPD density with location=0, scale=sigma, shape=xi
  if (min(length(x), length(sigma), length(xi)) == 0L) return(numeric(0))
  template <- x
  n <- max(length(x), length(sigma), length(xi))
  xv <- rep_len(as.numeric(x), n)
  sv <- rep_len(as.numeric(sigma), n)
  xiv <- rep_len(as.numeric(xi), n)
  out <- rep(NaN, n)
  valid <- !is.na(xv) & !is.na(sv) & !is.na(xiv) & is.finite(sv) & sv > 0 & is.finite(xiv)
  out[which(valid)] <- -Inf
  i <- which(valid & xv >= 0 & is.finite(xv))
  if (length(i)) {
    ls <- .gpd_logsurv_inside(xv[i], sv[i], xiv[i])
    out[i] <- -log(sv[i]) + (1 + xiv[i]) * ls
    out[i[xiv[i] < 0 & ls == -Inf]] <- -Inf
  }
  res <- if (log) out else exp(out)
  res <- .egpd_apply_na_nan(res, list(x, sigma, xi))
  .egpd_restore_dim(res, template)
}

.qgpd_from_logsurv <- function(logsurv, sigma, xi, template = logsurv, active = list(logsurv, sigma, xi)) {
  n <- max(length(logsurv), length(sigma), length(xi))
  ls <- rep_len(as.numeric(logsurv), n)
  sv <- rep_len(as.numeric(sigma), n)
  xv <- rep_len(as.numeric(xi), n)
  out <- rep(NaN, n)
  valid <- !is.na(ls) & !is.na(sv) & !is.na(xv) & ls <= 0 & is.finite(sv) & sv > 0 & is.finite(xv)
  out[valid] <- -sv[valid] * ls[valid]
  nz <- which(valid & xv != 0 & ls < 0)
  w <- -xv[nz] * ls[nz]
  i <- nz[abs(w) > .Machine$double.eps]
  w <- w[abs(w) > .Machine$double.eps]
  if (length(i)) {
    e <- expm1(w)
    prod <- sv[i] * e
    val <- prod / xv[i]
    ## Reorder the scaling if the intermediate product lost its range.
    j <- which(!is.finite(prod) | abs(prod) < .Machine$double.xmin)
    if (length(j)) {
      ratio <- e[j] / xv[i[j]]
      val[j] <- sv[i[j]] * ratio
      k <- j[!is.finite(ratio) | ratio < .Machine$double.xmin]
      if (length(k)) {
        le <- log(abs(e[k]))
        le[is.infinite(e[k])] <- w[k][is.infinite(e[k])]
        val[k] <- exp(le + (log(sv[i[k]]) - log(abs(xv[i[k]]))))
      }
    }
    out[i] <- val
  }
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

.qgpd_std <- function(p, sigma, xi) {
  ## GPD quantile with location=0, scale=sigma, shape=xi
  if (min(length(p), length(sigma), length(xi)) == 0L) return(numeric(0))
  template <- p
  n <- max(length(p), length(sigma), length(xi))
  pv <- rep_len(as.numeric(p), n)
  ls <- rep(NaN, n)
  ok <- which(!is.na(pv) & pv >= 0 & pv <= 1)
  ls[ok] <- log1p(-pv[ok])
  .qgpd_from_logsurv(ls, sigma, xi, template = template, active = list(p, sigma, xi))
}

.qG6_logsurv <- function(p, prob, kappa, delta, p_upp = NULL) {
  if (min(length(p), length(prob), length(kappa), length(delta)) == 0L) {
    return(numeric(0))
  }
  find.root <- function(p, pu, prob, kappa, delta) {
    if (is.nan(p) || is.nan(pu) || is.nan(prob) || is.nan(kappa) || is.nan(delta)) {
      return(NaN)
    }
    if (is.na(p) || is.na(pu) || is.na(prob) || is.na(kappa) || is.na(delta)) {
      return(NA_real_)
    }
    if (!is.finite(prob) || prob < 0 || prob > 1 ||
        !is.finite(kappa) || kappa <= 0 ||
        !is.finite(delta) || delta <= 0) {
      return(NaN)
    }
    if (p < 0 || p > 1 || pu < 0 || pu > 1) return(NaN)
    if (p == 0) return(0)
    if (pu == 0) return(-Inf)
    if (prob == 1) return(if (pu < 0.5) log(-expm1(log1p(-pu) / kappa)) else log1p(-p^(1 / kappa)))
    if (prob == 0) return(if (pu < 0.5) log(-expm1(log1p(-pu) / delta)) else log1p(-p^(1 / delta)))
    ## Solve in t = -log(1 - u): absolute tolerances in u lose extreme tails.
    target <- pu
    f <- function(t) {
      lu <- if (t > log(2)) log1p(-exp(-t)) else log(-expm1(-t))
      -prob * expm1(kappa * lu) - (1 - prob) * expm1(delta * lu) - target
    }
    hi <- max(1, -log(target) + log(max(kappa, delta)) + 2)
    iter <- 0L
    while (iter < 1075L) {
      fhi <- f(hi)
      if (!is.finite(fhi) || fhi <= 0) break
      hi <- 2 * hi
      iter <- iter + 1L
    }
    if (!is.finite(f(hi)) || f(hi) > 0) return(-Inf)
    -uniroot(f, c(0, hi), tol = 1e-12)$root
  }
  n <- max(length(p), length(prob), length(kappa), length(delta))
  pv <- rep_len(as.numeric(p), n)
  puv <- if (is.null(p_upp)) 1 - pv else rep_len(as.numeric(p_upp), n)
  mapply(find.root, pv, puv, rep_len(as.numeric(prob), n),
         rep_len(as.numeric(kappa), n), rep_len(as.numeric(delta), n),
         USE.NAMES = FALSE)
}

.qG_logsurv <- function(p_low, p_upp, type, prob, kappa, delta,
                        lp_low = log(p_low), lp_upp = log(p_upp)) {
  ## Computes log(1 - Q_G(p_low)) accurately across both tails for valid inputs.
  ## Targets whose tail probability is tiny (or underflows on the linear scale)
  ## are inverted on the log scale.
  n <- length(p_low)
  out <- numeric(n)
  at0 <- lp_low == -Inf
  at1 <- lp_upp == -Inf
  out[at0] <- 0
  out[at1] <- -Inf
  live <- !at0 & !at1
  if (!any(live)) return(out)

  if (type == 1) {
    k <- kappa
    lo <- which(live & lp_low <= -log(2))
    if (length(lo)) out[lo] <- .log1mexp(lp_low[lo] / k[lo])
    up <- which(live & lp_low > -log(2))
    if (length(up)) {
      lmlu <- .egpd_log_mlog1m(lp_upp[up]) - log(k[up])
      out[up] <- lmlu + log(.expm1_over_x(exp(lmlu)))
      ## p > 1/2 but u = p^(1/kappa) <= 1/2 (kappa < 1): the form above
      ## cancels for large v = -log(u), so take log(1 - exp(-v)) directly
      vb <- which(lmlu >= log(log(2)))
      if (length(vb)) out[up[vb]] <- .log1mexp(-exp(lmlu[vb]))
    }
    return(out)
  }
  ## log-scale root for the mixture carrier, for tiny lower-tail targets of the
  ## truncated normal/beta carriers, for the truncated beta carrier with
  ## kappa -> 0 (whose closed-form inverse cancels), for type-5 targets whose
  ## G4 value G^(2 / kappa) is below the normal double range, and for targets
  ## beyond linear precision
  ## The central qbeta branch uses squared distances down to O(1e-12/k)
  ## before the small-distance series takes over. Bisect in A if those
  ## distances can be subnormal: their square roots need not be tiny.
  beta_extreme <- type == 3 & kappa > 1e-12 / .Machine$double.xmin
  root <- live & (type == 6 | beta_extreme | lp_low < -700 | lp_upp < -700 |
                    (type %in% c(2, 3) & p_low < 0.01) |
                    (type == 3 & kappa < 1e-2) |
                    (type == 5 & (2 / kappa) * lp_low < -700))
  ## lower-tail targets of the truncated normal/beta carriers down to 1e-10:
  ## the closed-form inverse (accurate to ~1e-9 relative there) polished by
  ## Newton steps in log(A); rows that do not converge are bisected
  nw <- which(root & type %in% c(2, 3) & p_low < 0.01 & lp_low >= log(1e-10) &
                lp_upp >= -700 & !beta_extreme)
  if (length(nw)) {
    pair <- if (type == 2) .qG2_pair(p_low[nw], p_upp[nw], kappa[nw])
            else .qG3_pair(p_low[nw], p_upp[nw], kappa[nw])
    A_nw <- .egpd_carrier_newton(log(-log1p(-pair$u)), lp_low[nw], type, prob[nw],
                                 kappa[nw], delta[nw])
    ok_nw <- !is.na(A_nw)
    out[nw[ok_nw]] <- -A_nw[ok_nw]
    root[nw[ok_nw]] <- FALSE
    live[nw[ok_nw]] <- FALSE
  }
  ir <- which(root)
  if (length(ir)) {
    lower <- lp_low[ir] <= -log(2)
    A <- numeric(length(ir))
    il <- which(lower)
    if (length(il)) {
      A[il] <- .egpd_carrier_root(lp_low[ir[il]], TRUE, type, prob[ir[il]],
                                  kappa[ir[il]], delta[ir[il]])
    }
    iu <- which(!lower)
    if (length(iu)) {
      A[iu] <- .egpd_carrier_root(lp_upp[ir[iu]], FALSE, type, prob[ir[iu]],
                                  kappa[ir[iu]], delta[ir[iu]])
    }
    out[ir] <- -A
  }
  mid <- which(live & !root)
  if (!length(mid)) return(out)

  pl <- p_low[mid]; pu <- p_upp[mid]
  if (type == 1) {
    kv <- kappa[mid]
    lu <- ifelse(pu < 0.5, log1p(-pu) / kv, log(pl) / kv)
    out[mid] <- ifelse(lu < -log(2), log1p(-exp(lu)), log(-expm1(lu)))
  } else if (type == 2) {
    pair <- .qG2_pair(pl, pu, kappa[mid])
    out[mid] <- ifelse(pair$s < 0.5, log(pair$s), log1p(-pair$u))
  } else if (type == 3) {
    pair <- .qG3_pair(pl, pu, kappa[mid])
    out[mid] <- ifelse(pair$s < 0.5, log(pair$s), log1p(-pair$u))
  } else if (type == 4) {
    pair <- .qG4_pair(pl, pu, delta[mid])
    out[mid] <- ifelse(pair$s < 0.5, log(pair$s), log1p(-pair$u))
  } else if (type == 5) {
    kv <- kappa[mid]
    lu4 <- ifelse(pu < 0.5, (2 / kv) * log1p(-pu), (2 / kv) * log(pl))
    pl4 <- exp(lu4)
    pu4 <- -expm1(lu4)
    pair <- .qG4_pair(pl4, pu4, delta[mid])
    out[mid] <- ifelse(pair$s < 0.5, log(pair$s), log1p(-pair$u))
  } else {
    out[mid] <- .qG6_logsurv(pl, prob[mid], kappa[mid], delta[mid], p_upp = pu)
  }
  out
}

## ---------------------------------------------------------------------------
## Log-scale carrier functions
##
## Everything below is written in terms of A = -log(1 - F_GPD(y)), so that the
## GPD probability u = 1 - exp(-A) and its complement s = exp(-A) are both
## available without cancellation or underflow.  For each carrier G these give
## log G(u) and log(1 - G(u)) (computed directly in the tail where each is
## small and derived through log1mexp() otherwise), the log carrier density,
## and a vectorised log-scale inverse.
## ---------------------------------------------------------------------------

## log(1 - exp(x)) for x <= 0
.log1mexp <- function(x) {
  out <- rep(-Inf, length(x))
  ok <- which(!is.na(x) & x < 0)
  if (length(ok)) {
    xo <- x[ok]
    out[ok] <- ifelse(xo > -log(2), log(-expm1(xo)), log1p(-exp(xo)))
  }
  out[is.na(x)] <- x[is.na(x)]
  out
}

## log(exp(a) + exp(b))
.logspace_add <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[which(m == -Inf)] <- -Inf
  out[which(m == Inf)] <- Inf
  out
}

## log(u) with u = 1 - exp(-A), A >= 0
.egpd_log_u <- function(A) {
  ifelse(A > log(2), log1p(-exp(-A)), log(-expm1(-A)))
}

## log(-log(u)) with u = 1 - exp(-A), A > 0 (about -A for large A)
.egpd_log_mlogu <- function(A) {
  out <- suppressWarnings(log(-.egpd_log_u(A)))
  big <- which(A > 30)
  if (length(big)) out[big] <- -A[big] + log1p(exp(-A[big]) / 2)
  out
}

## log(-log1p(-x)) given lx = log(x), 0 < x < 1
.egpd_log_mlog1m <- function(lx) {
  out <- suppressWarnings(log(-log1p(-exp(lx))))
  small <- which(lx < -30)
  if (length(small)) out[small] <- lx[small] + exp(lx[small]) / 2
  out
}

## Phi(sqrt(k)) - 1/2, accurate for small k
.egpd_C2 <- function(k) 0.5 * pchisq(k, df = 1)

## S2(s) / s for carrier 2, where S2 = 1 - G2(1 - s); finite as s -> 0
.egpd_S2_over_s <- function(s, k) {
  out <- numeric(length(s))
  small <- k < 1e-2
  if (any(small)) {
    out[small] <- .erf_series(k[small] * s[small]^2) / .erf_series(k[small])
  }
  big <- which(!small)
  if (length(big)) {
    sk <- sqrt(k[big]); sb <- s[big]; z <- sk * sb
    num <- ifelse(z < 1e-2, (sk / sqrt(2 * pi)) * .erf_series(z^2),
                  0.5 * pchisq(z^2, df = 1) / sb)
    out[big] <- num / .egpd_C2(k[big])
  }
  out
}

## Row-wise max / min of a matrix without apply()
.egpd_row_max <- function(M) {
  r <- M[, 1L]
  for (j in seq_len(ncol(M))[-1L]) r <- pmax(r, M[, j])
  r
}
.egpd_row_min <- function(M) {
  r <- M[, 1L]
  for (j in seq_len(ncol(M))[-1L]) r <- pmin(r, M[, j])
  r
}

## 20-point Gauss-Legendre rule on [-1, 1] (Golub-Welsch), computed once
.egpd_gl_cache <- new.env(parent = emptyenv())
.egpd_gl20 <- function() {
  if (is.null(.egpd_gl_cache$rule)) {
    n <- 20L
    i <- seq_len(n - 1L)
    b <- i / sqrt(4 * i^2 - 1)
    J <- matrix(0, n, n)
    J[cbind(i, i + 1L)] <- b
    J[cbind(i + 1L, i)] <- b
    e <- eigen(J, symmetric = TRUE)
    .egpd_gl_cache$rule <- list(x = e$values, w = 2 * e$vectors[1, ]^2)
  }
  .egpd_gl_cache$rule
}

## log of P(lo < X < lo + width) for X ~ Beta(k, k), accurate for narrow intervals
.egpd_log_beta_interval <- function(lo, width, k) {
  out <- numeric(length(lo))
  gl <- width <= lo & abs(k - 1) * width / lo <= 2
  ig <- which(gl)
  if (length(ig)) {
    rule <- .egpd_gl20()
    x <- outer(width[ig] / 2, 1 + rule$x) + lo[ig]
    ld <- dbeta(x, k[ig], k[ig], log = TRUE) + rep(log(rule$w), each = length(ig))
    m <- .egpd_row_max(ld)
    out[ig] <- log(width[ig] / 2) + m + log(rowSums(exp(ld - m)))
  }
  id <- which(!gl)
  if (length(id)) {
    ## F(hi) - F(lo) directly, or via the symmetry identity
    ## P(x < X < 1/2) = I_{(1 - 2x)^2}(1/2, k) / 2 (x <= 1/2); use whichever
    ## difference cancels less
    kk <- k[id]
    hi <- pbeta(lo[id] + width[id], kk, kk, log.p = TRUE)
    lw <- pbeta(lo[id], kk, kk, log.p = TRUE)
    b0 <- pbeta((1 - 2 * lo[id])^2, 0.5, kk, log.p = TRUE)
    b1 <- pbeta(pmax(1 - 2 * (lo[id] + width[id]), 0)^2, 0.5, kk, log.p = TRUE)
    out[id] <- ifelse(lw - hi <= b1 - b0,
                      hi + .log1mexp(lw - hi),
                      log(0.5) + b0 + .log1mexp(b1 - b0))
  }
  out
}

## log normalising constant of carrier 3, P(1/32 < X < 1/2) for X ~ Beta(k, k);
## P(1/2 - t < X < 1/2) = pbeta(4 t^2, 1/2, k) / 2 avoids cancellation
.egpd_logZ3 <- function(k) {
  log(0.5) + pbeta((1 - 2 / 32)^2, 0.5, k, log.p = TRUE)
}

## Beta(1/2,k) at x = (15/16)^2 exp(-2A). Do not round x to a
## subnormal before accounting for k: k*x can still be of order one.
.egpd_log_beta3 <- function(A, k, lower.tail) {
  ls <- log(15/16) - A
  lx <- 2 * ls
  reg <- lx >= log(.Machine$double.xmin)
  res <- numeric(length(A))
  res[reg] <- pbeta(exp(lx[reg]), .5, k[reg], lower.tail = lower.tail, log.p = TRUE)
  lkx <- log(k) + lx
  flat <- !reg & lkx < log(.Machine$double.eps)
  if (any(flat)) {
    lp <- log(2) + ls[flat] - lbeta(.5, k[flat])
    res[flat] <- if (lower.tail) lp else .log1mexp(lp)
  }
  ## Here k exceeds eps/xmin and k*x is bounded by xmax*xmin:
  ## the scaled beta and Gamma(1/2,1) differ below double precision.
  big <- !reg & !flat
  res[big] <- pgamma(exp(lkx[big]), .5, lower.tail = lower.tail, log.p = TRUE)
  res
}

## G4 = (1+1/d) integral_0^A exp(-t) (1-exp(-d*t)) dt. For tiny u,
## replacing exp(-t) by one and A by u costs relative O(u), uniformly in d.
.egpd_logG4_small <- function(lu, d) {
  w <- exp(log(d) + lu)
  lv <- log1p(d) + 2*lu - log(2)
  i <- w >= .Machine$double.eps
  lv[i] <- log1p(1/d[i]) + lu[i] + log(.m1_expm1_over_x(w[i]))
  lv
}

.egpd_logG4 <- function(u, d) {
  small <- u < 1e-100
  out <- numeric(length(u))
  out[!small] <- log(.pG4_inner(u[!small], d[!small]))
  out[small] <- .egpd_logG4_small(log(u[small]), d[small])
  out
}

## Beta constants cancel in the truncated law. Columns are log integral
## of exp((k-1)*L) and its k-derivative, with L = log(1-(15*s/16)^2).
.egpd_tb_norm <- function(k) {
  keys <- unique(k)
  egpd_tb_moments_cpp(keys, 0, 1)[match(k, keys), , drop = FALSE]
}

.egpd_tb_logdens <- function(s, k) {
  (k - 1) * log1p(-((15/16) * s)^2) - .egpd_tb_norm(k)[, 1L]
}

## log G(u) at u = 1 - exp(-A), 0 < A < Inf; relatively accurate wherever G
## is small (in particular as u -> 0)
.egpd_lG_direct <- function(A, type, prob, kappa, delta) {
  logu <- .egpd_log_u(A)
  u <- -expm1(-A)
  if (type == 1L) return(kappa * logu)
  if (type == 2L) {
    out <- numeric(length(A))
    k <- kappa
    unif <- k == 0
    ser <- !unif & k * u <= 1
    is <- which(ser)
    if (length(is)) {
      a <- sqrt(k[is]); x <- a * u[is]
      cj_prev <- 0; cj <- 1; S <- 1; xp <- 1
      for (j in 0:29) {
        cn <- (a * cj - cj_prev) / (j + 1)
        xp <- xp * x
        S <- S + cn * xp / (j + 2)
        cj_prev <- cj; cj <- cn
      }
      out[is] <- dnorm(a, log = TRUE) - log(.egpd_C2(k[is])) + log(a) + logu[is] + log(S)
    }
    id <- which(!ser & !unif)
    if (length(id)) {
      a <- sqrt(k[id])
      hi <- pnorm(-a * exp(-A[id]), log.p = TRUE)
      lw <- pnorm(-a, log.p = TRUE)
      out[id] <- hi + .log1mexp(lw - hi) - log(.egpd_C2(k[id]))
    }
    out[unif] <- logu[unif]
    return(out)
  }
  if (type == 3L) {
    om <- 1 / 32
    hml <- 0.5 - om
    center <- A > 1 & kappa >= 1
    out <- numeric(length(A))
    if (any(!center))
      out[!center] <- .egpd_log_beta_interval(rep(om, sum(!center)), hml*u[!center],
                                              kappa[!center]) - .egpd_logZ3(kappa[!center])
    if (any(center)) {
      ## The left tail can be tiny even where u rounds to one.
      k <- kappa[center]
      hi <- .egpd_log_beta3(A[center], k, FALSE)
      lo <- pbeta((15/16)^2, .5, k, lower.tail = FALSE, log.p = TRUE)
      out[center] <- hi + .log1mexp(lo - hi) - pbeta((15/16)^2, .5, k, log.p = TRUE)
    }
    return(out)
  }
  if (type %in% c(4L, 5L)) {
    lG4 <- .egpd_logG4(u, delta)
    return(if (type == 4L) lG4 else (kappa / 2) * lG4)
  }
  .logspace_add(log(prob) + kappa * logu, log1p(-prob) + delta * logu)
}

## Direct log survivor. For power carriers this also works below A = log(2).
.egpd_lS_direct <- function(A, type, prob, kappa, delta) {
  lS1 <- function(k) {
    lv <- log(k) + .egpd_log_mlogu(A)
    v <- exp(lv)
    ifelse(v < 0.5, lv + log(.expm1_over_x(v)), .log1mexp(-v))
  }
  if (type == 1L) return(lS1(kappa))
  if (type == 6L)
    return(.logspace_add(log(prob) + lS1(kappa), log1p(-prob) + lS1(delta)))
  if (type == 2L) return(-A + log(.egpd_S2_over_s(exp(-A), kappa)))
  if (type == 3L) {
    return(.egpd_log_beta3(A, kappa, TRUE) - pbeta((15/16)^2, .5, kappa, log.p = TRUE))
  }
  w <- delta * A
  lS4 <- -A + log1p(ifelse(w < 1e-2, A * .expm1_over_x(w), -expm1(-w) / delta))
  if (type == 4L) return(lS4)
  if (type == 5L) {
    lv <- log(kappa / 2) + .egpd_log_mlog1m(lS4)
    return(lv + log(.expm1_over_x(exp(lv))))
  }
}

## list(lG = log G(u), lS = log(1 - G(u))) at u = 1 - exp(-A), A in [0, Inf]
.egpd_carrier_log <- function(A, type, prob, kappa, delta) {
  n <- length(A)
  lG <- rep(NaN, n)
  lS <- rep(NaN, n)
  at0 <- which(!is.na(A) & A <= 0)
  at1 <- which(!is.na(A) & A == Inf)
  lG[at0] <- -Inf; lS[at0] <- 0
  lG[at1] <- 0; lS[at1] <- -Inf
  ## compute the smaller of G and 1 - G directly and the other by log1mexp()
  mid <- which(!is.na(A) & A > 0 & A < Inf)
  if (length(mid)) {
    g <- .egpd_lG_direct(A[mid], type, prob[mid], kappa[mid], delta[mid])
    use_g <- !is.na(g) & (g <= -log(2) | (type != 6L & A[mid] < log(2)))
    lG[mid[use_g]] <- g[use_g]
    lS[mid[use_g]] <- .log1mexp(g[use_g])
    hi <- mid[!use_g]
    if (length(hi)) {
      sv <- .egpd_lS_direct(A[hi], type, prob[hi], kappa[hi], delta[hi])
      lS[hi] <- sv
      lG[hi] <- .log1mexp(sv)
    }
  }
  list(lG = lG, lS = lS)
}

## Type-4 log density; factor out delta before its reciprocal can overflow.
.egpd_logg4 <- function(A, delta) {
  w <- delta * A
  ifelse(w < 1e-2, log1p(delta) + log(A) + log(.expm1_over_x(w)),
         ifelse(delta < 1, log1p(delta) - log(delta), log1p(1 / delta)) + log(-expm1(-w)))
}

## log g(u) at u = 1 - exp(-A), 0 < A < Inf
.egpd_carrier_logdens <- function(A, type, prob, kappa, delta) {
  logu <- .egpd_log_u(A)
  if (type == 1L) return(log(kappa) + ifelse(kappa == 1, 0, (kappa - 1) * logu))
  if (type == 2L) {
    s <- exp(-A)
    k <- rep_len(kappa, length(A))
    out <- 0.5 * log(k) - 0.5 * log(2 * pi) - 0.5 * k * s^2 - log(.egpd_C2(k))
    sm <- which(k < 1e-2)
    if (length(sm)) out[sm] <- -0.5 * k[sm] * s[sm]^2 - log(.erf_series(k[sm]))
    return(out)
  }
  if (type == 3L) return(.egpd_tb_logdens(exp(-A), kappa))
  lg4 <- .egpd_logg4(A, delta)
  if (type == 4L) return(lg4)
  if (type == 5L) {
    lG4 <- .egpd_carrier_log(A, 4L, prob, kappa, delta)$lG
    return(log(kappa / 2) + ifelse(kappa == 2, 0, (kappa / 2 - 1) * lG4) + lg4)
  }
  .logspace_add(log(prob * kappa) + ifelse(kappa == 1, 0, (kappa - 1) * logu),
                log((1 - prob) * delta) + ifelse(delta == 1, 0, (delta - 1) * logu))
}

## Newton iterations in x = log(A) for the lower-tail equation log G = lp, with
## d log G / d x = A (1 - u) g(u) / G(u).  Starts x with non-finite entries,
## and rows whose correction does not fall below 1e-9 (after which the
## remaining error is ~1e-18) within maxit steps, are returned as NA.
.egpd_carrier_newton <- function(x, lp, type, prob, kappa, delta, maxit = 6L) {
  n <- length(x)
  done <- rep(FALSE, n)
  act <- which(is.finite(x))
  for (it in seq_len(maxit)) {
    if (!length(act)) break
    A <- exp(x[act])
    lG <- .egpd_carrier_log(A, type, prob[act], kappa[act], delta[act])$lG
    lg <- .egpd_carrier_logdens(A, type, prob[act], kappa[act], delta[act])
    step <- (lG - lp[act]) / exp(lg - A + x[act] - lG)
    step <- pmin(pmax(step, -2), 2)
    x[act] <- x[act] - step
    fin <- !is.na(step) & is.finite(x[act])
    conv <- fin & abs(step) <= 1e-9
    done[act[conv]] <- TRUE
    act <- act[fin & !conv]
  }
  out <- rep(NA_real_, n)
  out[done] <- exp(x[done])
  out
}

## Vectorised log-scale inverse: A with log G = lp (lower = TRUE) or
## log(1 - G) = lp (lower = FALSE), by bisection in log(A) (relative accuracy).
.egpd_carrier_root <- function(lp, lower, type, prob, kappa, delta) {
  n <- length(lp)
  lo <- rep(-760, n)
  hi <- pmax(7, log(2 * abs(lp) + 50))
  f <- function(x, i) {
    cl <- .egpd_carrier_log(exp(x), type, prob[i], kappa[i], delta[i])
    if (lower) cl$lG - lp[i] else lp[i] - cl$lS
  }
  idx <- seq_len(n)
  ## make sure the upper end brackets the root (f increasing in x)
  for (it in 1:20) {
    bad <- which(f(hi, idx) < 0)
    if (!length(bad)) break
    hi[bad] <- hi[bad] + 2
  }
  for (it in 1:70) {
    mid <- (lo + hi) / 2
    fm <- f(mid, idx)
    up <- !is.na(fm) & fm >= 0
    hi[up] <- mid[up]
    lo[!up] <- mid[!up]
  }
  exp((lo + hi) / 2)
}

## Transformation functions G, g, G^{-1} for EGPD types 1-6

#' Transformation CDF for EGPD
#'
#' Near \eqn{u = 0} (carriers 2 and 3), for carrier 3 with \code{kappa < 0.01}
#' and where the carrier-4 value is below the double range (carrier 5) the CDF
#' is computed on the log scale, so it keeps its relative accuracy there.
#'
#' @details Carrier parameters required by \code{type} must be supplied: a missing
#'   one is an error, whereas supplied \code{NA} values give \code{NA} results.
#' @param u values in (0,1)
#' @param type integer 1-6 specifying the transformation type
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter
#' @param delta shape parameter (types 4-6)
#'
#' @return Transformed probabilities
#' @export
p.G <- function(u, type = 1, prob, kappa, delta) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("prob", "kappa", "delta")), gpd = FALSE)
  if (length(u) == 0L ||
      (!missing(prob) && length(prob) == 0L) ||
      (!missing(kappa) && length(kappa) == 0L) ||
      (!missing(delta) && length(delta) == 0L)) {
    return(numeric(0))
  }

  template <- u
  prob_arg  <- if (missing(prob))  NA_real_ else prob
  kappa_arg <- if (missing(kappa)) NA_real_ else kappa
  delta_arg <- if (missing(delta)) NA_real_ else delta

  n <- max(length(u), length(prob_arg), length(kappa_arg), length(delta_arg))
  uv <- rep_len(as.numeric(u), n)
  pv <- rep_len(as.numeric(prob_arg), n)
  kv <- rep_len(as.numeric(kappa_arg), n)
  dv <- rep_len(as.numeric(delta_arg), n)

  c_ok <- .egpd_carrier_valid(type, pv, kv, dv, n)
  out <- rep(NaN, n)
  out[which(c_ok & !is.na(uv) & uv <= 0)] <- 0
  out[which(c_ok & !is.na(uv) & uv >= 1)] <- 1

  inside <- which(c_ok & !is.na(uv) & uv > 0 & uv < 1)
  if (length(inside)) {
    ui <- uv[inside]; pi_v <- pv[inside]; ki <- kv[inside]; di <- dv[inside]
    g4 <- if (type == 5) .pG4_inner(ui, di)
    res <- if (type == 1) {
      ui^ki
    } else if (type == 2) {
      .pG2_inner(ui, ki)
    } else if (type == 3) {
      lower <- 1 / 32; upper <- 0.5
      nf <- pbeta(upper, ki, ki) - pbeta(lower, ki, ki)
      (pbeta((upper - lower) * ui + lower, ki, ki) - pbeta(lower, ki, ki)) / nf
    } else if (type == 4) {
      .pG4_inner(ui, di)
    } else if (type == 5) {
      g4^(ki / 2)
    } else {
      pi_v * ui^ki + (1 - pi_v) * ui^di
    }
    ## the log-scale carrier where the linear forms lose precision: the CDF
    ## differences of carriers 2 and 3 near u = 0 (kappa = 0 is the exact
    ## uniform carrier), the carrier-3 normaliser as kappa -> 0, and G4 below
    ## the normal double range for carrier 5
    lg <- if (type == 2) {
      which(ui < 0.01 & ki > 0)
    } else if (type == 3) {
      which(ui < 0.01 | ki < 1e-2)
    } else if (type == 5) {
      which(g4 < exp(-700))
    } else {
      integer(0)
    }
    if (length(lg)) {
      res[lg] <- exp(.egpd_carrier_log(-log1p(-ui[lg]), type, pi_v[lg], ki[lg], di[lg])$lG)
    }
    out[inside] <- res
  }

  active <- c(list(u), .egpd_active_carrier_args(type, prob_arg, kappa_arg, delta_arg))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' Transformation density for EGPD
#'
#' @details Carrier parameters required by \code{type} must be supplied: a missing
#'   one is an error, whereas supplied \code{NA} values give \code{NA} results.
#' @param u values in (0,1)
#' @param type integer 1-6 specifying the transformation type
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter
#' @param delta shape parameter (types 4-6)
#' @param log logical: return log-density?
#'
#' @return Transformation density values
#' @export
d.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA, log = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("prob", "kappa", "delta")), gpd = FALSE)
  if (min(length(u), length(prob), length(kappa), length(delta)) == 0L) {
    return(numeric(0))
  }

  template <- u
  n <- max(length(u), length(prob), length(kappa), length(delta))
  uv <- rep_len(as.numeric(u), n)
  pv <- rep_len(as.numeric(prob), n)
  kv <- rep_len(as.numeric(kappa), n)
  dv <- rep_len(as.numeric(delta), n)

  c_ok <- .egpd_carrier_valid(type, pv, kv, dv, n)
  out <- rep(NaN, n)
  outside <- which(c_ok & !is.na(uv) & (uv < 0 | uv > 1))
  out[outside] <- if (log) -Inf else 0

  inside <- which(c_ok & !is.na(uv) & uv >= 0 & uv <= 1)
  if (length(inside)) {
    ui <- uv[inside]; pi_v <- pv[inside]; ki <- kv[inside]; di <- dv[inside]
    ld <- if (type == 1) {
      res <- log(ki) + (ki - 1) * log(ui)
      res[which(ki == 1)] <- 0
      res
    } else if (type == 2) {
      .log_dG2_inner(ui, ki)
    } else if (type == 3) {
      .egpd_tb_logdens(1 - ui, ki)
    } else if (type == 4) {
      .egpd_logg4(-log1p(-ui), di)
    } else if (type == 5) {
      res <- numeric(length(ui))
      at0 <- which(ui == 0)
      at1 <- which(ui == 1)
      mid <- which(ui > 0 & ui < 1)
      if (length(at0)) {
        res[at0] <- ifelse(ki[at0] < 1, Inf,
                           ifelse(ki[at0] > 1, -Inf, 0.5 * log((di[at0] + 1) / 2)))
      }
      if (length(at1)) {
        res[at1] <- log(ki[at1] / 2) + .egpd_logg4(Inf, di[at1])
      }
      if (length(mid)) {
        um <- ui[mid]; km <- ki[mid]; dm <- di[mid]
        lG4 <- .egpd_logG4(um, dm)
        pwr <- (km / 2 - 1) * lG4
        pwr[which(km == 2)] <- 0
        res[mid] <- log(km / 2) + pwr + .egpd_logg4(-log1p(-um), dm)
      }
      res
    } else {
      ## log-scale mixture (u^(kappa - 1) can under- or overflow); a
      ## zero-weight component contributes nothing, also at u = 0
      lu <- log(ui)
      l1 <- log(pi_v) + log(ki) + ifelse(ki == 1, 0, (ki - 1) * lu)
      l2 <- log1p(-pi_v) + log(di) + ifelse(di == 1, 0, (di - 1) * lu)
      l1[pi_v == 0] <- -Inf
      l2[pi_v == 1] <- -Inf
      .logspace_add(l1, l2)
    }
    out[inside] <- if (log) ld else exp(ld)
  }

  active <- c(list(u), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' Inverse transformation (quantile) for EGPD
#'
#' Numeric arguments are recycled to their common maximum length. The mixture
#' carrier (type 6), the lower tails of carriers 2 and 3, carrier 3 with
#' very small or extreme positive \code{kappa}, and carrier-5 targets whose intermediate value
#' \eqn{u^{2/\kappa}} is below the double range are inverted by a log-scale
#' root search with relative accuracy, preserving accuracy near both ends of
#' the unit interval.
#'
#' @details Carrier parameters required by \code{type} must be supplied: a missing
#'   one is an error, whereas supplied \code{NA} values give \code{NA} results.
#' @param u values in (0,1)
#' @param type integer 1-6
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter
#' @param delta shape parameter (types 4-6)
#'
#' @return Inverse-transformed values
#' @export
q.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("prob", "kappa", "delta")), gpd = FALSE)
  if (min(length(u), length(prob), length(kappa), length(delta)) == 0L) {
    return(numeric(0))
  }

  template <- u
  n <- max(length(u), length(prob), length(kappa), length(delta))
  uv <- rep_len(as.numeric(u), n)
  pv <- rep_len(as.numeric(prob), n)
  kv <- rep_len(as.numeric(kappa), n)
  dv <- rep_len(as.numeric(delta), n)

  c_ok <- .egpd_carrier_valid(type, pv, kv, dv, n)
  valid <- which(c_ok & !is.na(uv) & uv >= 0 & uv <= 1)
  out <- rep(NaN, n)

  if (length(valid)) {
    ui <- uv[valid]; pi_v <- pv[valid]; ki <- kv[valid]; di <- dv[valid]
    out[valid] <- if (type == 1) {
      ui^(1 / ki)
    } else if (type %in% c(2, 3)) {
      -expm1(.qG_logsurv(ui, 1 - ui, type = type, prob = pi_v, kappa = ki, delta = di))
    } else if (type == 4) {
      .qG4_pair(ui, 1 - ui, di)$u
    } else if (type == 5) {
      u4 <- ui^(2 / ki)
      s4 <- -expm1((2 / ki) * log(ui))
      s4[ui == 0] <- 1
      s4[ui == 1] <- 0
      res <- numeric(length(ui))
      ## G4 targets u^(2 / kappa) below the normal double range are inverted
      ## on the log scale
      deep <- ui > 0 & (2 / ki) * log(ui) < -700
      sh <- which(!deep)
      if (length(sh)) res[sh] <- .qG4_pair(u4[sh], s4[sh], di[sh])$u
      dp <- which(deep)
      if (length(dp)) {
        res[dp] <- -expm1(.qG_logsurv(ui[dp], 1 - ui[dp], type = 5, prob = pi_v[dp],
                                      kappa = ki[dp], delta = di[dp]))
      }
      res
    } else {
      -expm1(.qG_logsurv(ui, 1 - ui, type = type, prob = pi_v, kappa = ki, delta = di))
    }
  }

  active <- c(list(u), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' Random generation from transformation G
#'
#' @details Carrier parameters required by \code{type} must be supplied: a missing
#'   one is an error, whereas supplied \code{NA} values give \code{NA} results.
#'   For the random generators, \code{n} follows base R: if \code{length(n) > 1}
#'   its length is used, so \code{numeric(0)} gives \code{numeric(0)}.
#' @param n number of samples
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter
#' @param delta shape parameter (types 4-6)
#' @param type integer 1-6
#' @param unifsamp optional uniform samples, recycled to length \code{n}. When
#'   supplied the draws are \code{q.G(unifsamp, ...)} (inverse-CDF sampling),
#'   except for \code{type = 6} with \code{direct = FALSE}, where the mixture
#'   component of each draw is still chosen with \code{runif}.
#' @param direct logical for type 6
#'
#' @return Random samples
#' @export
r.G <- function(n, prob = NA, kappa = NA, delta = NA, type = 1, unifsamp = NULL, direct = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("prob", "kappa", "delta")), gpd = FALSE)
  unifsamp <- .egpd_unif(n, unifsamp)
  n <- length(unifsamp)
  if (n == 0L) return(numeric(0))

  prob_r  <- if (missing(prob))  NA else rep_len(prob, n)
  kappa_r <- if (missing(kappa)) NA else rep_len(kappa, n)
  delta_r <- if (missing(delta)) NA else rep_len(delta, n)
  .egpd_rG_core(unifsamp, prob_r, kappa_r, delta_r, type, direct)$u
}

## Uniforms for r.G() and regpd(): n follows base R (.egpd_rn) and supplied
## samples are recycled to it; length 0 when there is nothing to draw.
.egpd_unif <- function(n, unifsamp) {
  n <- .egpd_rn(n)
  if (!is.null(unifsamp)) {
    if (length(unifsamp) == 0L || n == 0L) return(numeric(0))
    return(rep_len(as.numeric(unifsamp), n))
  }
  if (n == 0L) return(numeric(0))
  runif(n)
}

## Carrier draws u = G^{-1}(U) for r.G() and regpd().  For type 6 with
## direct = FALSE the mixture component of each draw is chosen with runif()
## and inverted as a type-1 carrier.  `k1` is the type-1 shape that each draw
## was inverted with (types 1 and 6; NA otherwise).
.egpd_rG_core <- function(U, prob_r, kappa_r, delta_r, type, direct) {
  n <- length(U)
  if (type != 6 || direct) {
    u <- q.G(U, prob = prob_r, kappa = kappa_r, delta = delta_r, type = type)
    k1 <- if (type == 1) suppressWarnings(rep_len(as.numeric(kappa_r), n)) else rep(NA_real_, n)
    return(list(u = u, k1 = k1))
  }

  pv <- as.numeric(prob_r)
  kv <- as.numeric(kappa_r)
  dv <- as.numeric(delta_r)
  u_comp <- runif(n)
  res <- rep(NaN, n)
  k1 <- rep(NA_real_, n)
  c_ok <- .egpd_carrier_valid(6L, pv, kv, dv, n)
  idx1 <- which(c_ok & (u_comp < pv))
  idx2 <- which(c_ok & !(u_comp < pv))
  if (length(idx1) > 0L) {
    res[idx1] <- q.G(U[idx1], kappa = kv[idx1], type = 1)
    k1[idx1] <- kv[idx1]
  }
  if (length(idx2) > 0L) {
    res[idx2] <- q.G(U[idx2], kappa = dv[idx2], type = 1)
    k1[idx2] <- dv[idx2]
  }
  list(u = .egpd_apply_na_nan(res, list(U, prob_r, kappa_r, delta_r)), k1 = k1)
}

############################################
##    Continuous EGPD                     ##
############################################

#' CDF of the Extended GPD
#'
#' The carrier is evaluated on the log scale from the GPD log-survival, so
#' \code{log.p = TRUE} results stay accurate far into either tail, where the
#' linear-scale probability underflows.
#'
#' @details Parameters required by \code{type} (the GPD \code{sigma} and \code{xi} and the
#'   carrier parameters) must be supplied: a missing one is an error, whereas
#'   supplied \code{NA} values give \code{NA} results.
#' @param q quantiles
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le q)},
#'   otherwise \eqn{P(X > q)}
#' @param log.p logical; if TRUE, probabilities are returned on the log scale
#'
#' @return CDF values
#' @export
pegpd <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA,
                  type = 1, lower.tail = TRUE, log.p = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta")))
  if (min(length(q), length(prob), length(kappa), length(delta),
          length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- q
  n <- max(length(q), length(prob), length(kappa), length(delta),
           length(sigma), length(xi))
  qv <- rep_len(as.numeric(q), n)
  pv <- rep_len(as.numeric(prob), n)
  kv <- rep_len(as.numeric(kappa), n)
  dv <- rep_len(as.numeric(delta), n)
  sv <- rep_len(as.numeric(sigma), n)
  xv <- rep_len(as.numeric(xi), n)

  ## A = -log(1 - F_GPD(q)); the carrier is evaluated on the log scale so that
  ## neither tail underflows before log.p is applied
  A <- -.gpd_cdf_surv_std(qv, sv, xv)$log_surv
  c_ok <- .egpd_carrier_valid(type, pv, kv, dv, n) &
          !is.na(sv) & is.finite(sv) & sv > 0 & !is.na(xv) & is.finite(xv)

  out <- rep(NaN, n)
  ok <- which(c_ok & !is.na(A))
  if (length(ok)) {
    cl <- .egpd_carrier_log(A[ok], type, pv[ok], kv[ok], dv[ok])
    lv <- if (lower.tail) cl$lG else cl$lS
    out[ok] <- if (log.p) lv else exp(lv)
  }
  if (type %in% c(1L, 4L, 5L, 6L)) {
    j <- ok[qv[ok] > 0 & A[ok] < .Machine$double.xmin]
    if (length(j)) {
      lu <- .gpd_logu_small(qv[j], sv[j], xv[j])
      lg <- if (type == 6L) .logspace_add(log(pv[j])+kv[j]*lu, log1p(-pv[j])+dv[j]*lu) else
        if (type == 1L) kv[j]*lu else .egpd_logG4_small(lu,dv[j])
      if (type == 5L) lg <- kv[j]*(lg/2)
      val <- if (lower.tail) lg else .log1mexp(pmin(lg, 0))
      out[j] <- if (log.p) val else exp(val)
    }
  }
  if (!lower.tail && !log.p && type == 1L) {
    ## A unit-power carrier has an algebraic GPD survivor. A log/exp
    ## round trip can move even an exact dyadic jump to the wrong side
    ## of a quantile target. Keep the log1p route near shape zero where
    ## adding one would lose the shape increment.
    i <- ok[kv[ok] == 1 & xv[ok] != 0 & A[ok] > 0 & is.finite(A[ok])]
    t <- xv[i] * (qv[i] / sv[i])
    base <- 1 + t
    edge <- which(t < -.5)
    if (length(edge)) {
      k <- i[edge]
      base[edge] <- gpd_endpoint_gap_cpp(qv[k], sv[k], xv[k])
    }
    j <- which(is.finite(base) & base > 0 & t != 0 &
                 (abs(t) >= .5 | base - 1 == t))
    out[i[j]] <- base[j]^(-1 / xv[i[j]])
  }
  ## With xi=-1 the base is uniform. Identity carriers must retain its
  ## exact cell probabilities, e.g. when profiling a zero-inflation atom.
  i <- ok[xv[ok] == -1 & qv[ok] > 0 & qv[ok] < sv[ok]]
  if (length(i)) {
    power <- switch(as.integer(type), kv[i], ifelse(kv[i] == 0, 1, NA_real_),
                    ifelse(kv[i] == 1, 1, NA_real_), ifelse(dv[i] == 1, 2, NA_real_),
                    ifelse(dv[i] == 1, kv[i], NA_real_),
                    ifelse(pv[i] == 1 | kv[i] == dv[i], kv[i],
                           ifelse(pv[i] == 0, dv[i], NA_real_)))
    if (lower.tail && !log.p) {
      ## A dyadic power or mixture can also have an exact zero mass. Use
      ## the direct power away from the endpoints where u loses accuracy.
      u <- qv[i] / sv[i]
      j <- which(u >= .Machine$double.xmin & u <= .5)
      if (type == 6L) {
        out[i[j]] <- pv[i[j]] * u[j]^kv[i[j]] + (1 - pv[i[j]]) * u[j]^dv[i[j]]
      } else {
        j <- j[!is.na(power[j])]
        out[i[j]] <- u[j]^power[j]
      }
    }
    i <- i[which(power == 1)]
    num <- if (lower.tail) qv[i] else sv[i] - qv[i]
    opp <- if (lower.tail) sv[i] - qv[i] else qv[i]
    p <- num / sv[i]
    out[i] <- if (!log.p) p else ifelse(p > .5, log1p(-opp / sv[i]),
      ifelse(p < .Machine$double.xmin, log(num) - log(sv[i]), log(p)))
  }

  active <- c(list(q, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' Density of the Extended GPD
#'
#' The shape-zero case uses the exact exponential GPD limit. Inside the support
#' the carrier density is evaluated from the GPD log-survival, which keeps the
#' far tail accurate when \eqn{1 - F_{GPD}(x)} is below machine precision. Values below zero,
#' at or beyond a finite upper endpoint, or at infinity have density zero
#' (log-density \code{-Inf}).
#'
#' @details Parameters required by \code{type} (the GPD \code{sigma} and \code{xi} and the
#'   carrier parameters) must be supplied: a missing one is an error, whereas
#'   supplied \code{NA} values give \code{NA} results.
#' @param x values
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param log logical: return log-density?
#'
#' @return Density values
#' @export
degpd_density <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, log = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta")))
  if (min(length(x), length(prob), length(kappa), length(delta),
          length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }
  template <- x
  n <- max(length(x), length(prob), length(kappa), length(delta), length(sigma), length(xi))
  xv <- rep_len(as.numeric(x), n)
  sv <- rep_len(as.numeric(sigma), n)
  xiv <- rep_len(as.numeric(xi), n)
  pv <- rep_len(prob, n); kv <- rep_len(kappa, n); dv <- rep_len(delta, n)

  param_ok <- .egpd_carrier_valid(type, as.numeric(pv), as.numeric(kv), as.numeric(dv), n) &
              !is.na(sv) & is.finite(sv) & sv > 0 & !is.na(xiv) & is.finite(xiv)
  outside <- param_ok & !is.na(xv) & (xv < 0 | is.infinite(xv))
  bound <- which(param_ok & is.finite(xv) & xv > 0 & xiv < 0 & xv >= -sv/xiv)
  if (length(bound))
    outside[bound] <- gpd_endpoint_gap_cpp(xv[bound], sv[bound], xiv[bound]) <= 0

  out <- rep(NaN, n)
  ## Strictly inside the support, evaluate the carrier from the GPD
  ## log-survival so that 1 - u is not rounded away in the far tail.
  inner <- which(param_ok & !is.na(xv) & !outside & xv > 0)
  if (length(inner)) {
    A <- -.gpd_cdf_surv_std(xv[inner], sv[inner], xiv[inner])$log_surv
    out[inner] <- .egpd_carrier_logdens(A, type, as.numeric(pv[inner]),
                                        as.numeric(kv[inner]), as.numeric(dv[inner])) +
      .dgpd_std(xv[inner], sigma = sv[inner], xi = xiv[inner], log = TRUE)
    if (type %in% c(1L, 4L, 5L, 6L)) {
      j <- inner[A < .Machine$double.xmin]
      if (length(j)) {
        lu <- .gpd_logu_small(xv[j], sv[j], xiv[j])
        if (type %in% c(4L,5L)) {
          w <- exp(log(dv[j]) + lu)
          lg <- log1p(dv[j]) + lu + log(.expm1_over_x(w))
          if (type == 5L) {
            lG <- .egpd_logG4_small(lu,dv[j])
            lg <- log(kv[j])-log(2) + (lg-lG) + kv[j]*(lG/2)
          }
        } else {
          lg <- log(kv[j]) + (kv[j]-1)*lu
          if (type == 6L) lg <- .logspace_add(log(pv[j])+lg,
            log1p(-pv[j])+log(dv[j])+(dv[j]-1)*lu)
        }
        out[j] <- lg + .dgpd_std(xv[j], sv[j], xiv[j], log=TRUE)
      }
    }
  }
  ## Remaining rows (x = 0, off support, invalid parameters): carrier density
  ## at u = F(x), with off-support rows moved to an interior point to avoid
  ## 0 * Inf before they are set to -Inf.
  rest <- setdiff(seq_len(n), inner)
  if (length(rest)) {
    eval_x <- xv[rest]
    off <- outside[rest]
    sr <- sv[rest]; xr <- xiv[rest]
    eval_x[off] <- 0.5 * sr[off]
    finite_tail <- which(off & xr < 0)
    eval_x[finite_tail] <- pmin(eval_x[finite_tail], -0.5 * sr[finite_tail] / xr[finite_tail])
    out[rest] <- d.G(.pgpd_std(eval_x, sigma = sr, xi = xr), prob = pv[rest], kappa = kv[rest],
                     delta = dv[rest], type = type, log = TRUE) +
      .dgpd_std(eval_x, sigma = sr, xi = xr, log = TRUE)
  }
  out[outside] <- -Inf
  if (!log) out <- exp(out)
  active <- c(list(x, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' Quantile function of the Extended GPD
#'
#' Shape zero uses the exponential limit. Log-scale probabilities
#' (\code{log.p = TRUE}) are inverted on the log scale, so tail probabilities
#' below the double-precision range still give finite quantiles. The mixture
#' carrier (type 6), the lower tails of the truncated-normal and
#' truncated-beta carriers (types 2 and 3), the truncated-beta carrier with
#' very small or extreme positive \code{kappa}, and carrier-5 probabilities whose intermediate value
#' \eqn{p^{2/\kappa}} is below the double range are inverted by a log-scale
#' root search with relative accuracy.
#'
#' @details Parameters required by \code{type} (the GPD \code{sigma} and \code{xi} and the
#'   carrier parameters) must be supplied: a missing one is an error, whereas
#'   supplied \code{NA} values give \code{NA} results.
#' @param p probabilities
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le q)},
#'   otherwise \eqn{P(X > q)}
#' @param log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}
#'
#' @return Quantile values
#' @export
qegpd <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA,
                  type = 1, lower.tail = TRUE, log.p = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta")))
  if (min(length(p), length(prob), length(kappa), length(delta),
          length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- p
  n <- max(length(p), length(prob), length(kappa), length(delta), length(sigma), length(xi))
  pr <- .egpd_parse_prob(rep_len(p, n), lower.tail = lower.tail, log.p = log.p)
  pv <- rep_len(as.numeric(prob), n)
  kv <- rep_len(as.numeric(kappa), n)
  dv <- rep_len(as.numeric(delta), n)
  sv <- rep_len(as.numeric(sigma), n)
  xv <- rep_len(as.numeric(xi), n)

  c_ok <- .egpd_carrier_valid(type, pv, kv, dv, n) & pr$valid
  lsurv <- rep(NaN, n)
  ok <- which(c_ok)
  if (length(ok)) {
    lsurv[ok] <- .qG_logsurv(pr$p_low[ok], pr$p_upp[ok], type = type,
                             prob = pv[ok], kappa = kv[ok], delta = dv[ok],
                             lp_low = pr$lp_low[ok], lp_upp = pr$lp_upp[ok])
  }

  active <- c(list(p, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .qgpd_from_logsurv(lsurv, sv, xv, template = template, active = active)
  if (type %in% c(1L, 4L, 5L, 6L)) {
    j <- which(c_ok & is.finite(sv) & sv > 0 & is.finite(xv) &
                 is.finite(pr$lp_low) & -lsurv < .Machine$double.xmin)
    if (length(j)) {
      if (type == 1L) {
        lu <- pr$lp_low[j]/kv[j]
      } else if (type %in% c(4L,5L)) {
        target <- pr$lp_low[j]
        if (type == 5L) {
          target <- 2*(target/kv[j])
          up <- which(pr$lp_upp[j] < -log(2))
          target[up] <- -exp(.egpd_log_mlog1m(pr$lp_upp[j[up]]) - log(kv[j[up]]) + log(2))
        }
        lu <- (target - log1p(dv[j]) + log(2))/2
        ## The quadratic lower-tail limit is exact to rounding unless d*u
        ## is appreciable. Only those rows need inversion in log(u).
        i <- which(log(dv[j])+lu > log(.Machine$double.eps))
        if (length(i)) {
          lo <- lu[i]; hi <- rep(log(.Machine$double.xmin)+1,length(i))
          for (it in seq_len(64L)) {
            mid <- (lo + hi)/2
            hit <- .egpd_logG4_small(mid, dv[j[i]]) >= target[i]
            hi[hit] <- mid[hit]; lo[!hit] <- mid[!hit]
          }
          lu[i] <- (lo + hi)/2
        }
      } else {
        ## Invert the mixture without materializing either component's u.
        ## Below this bracket even scaling by DBL_MAX gives rounded zero.
        lo <- rep(-1600, length(j))
        hi <- rep(0, length(j))
        lower <- pr$lp_low[j] <= -log(2)
        for (it in seq_len(64L)) {
          mid <- (lo + hi)/2
          a <- kv[j]*mid; b <- dv[j]*mid
          hit <- ifelse(lower,
            .logspace_add(log(pv[j])+a, log1p(-pv[j])+b) >= pr$lp_low[j],
            .logspace_add(log(pv[j])+.log1mexp(a), log1p(-pv[j])+.log1mexp(b)) <= pr$lp_upp[j])
          hi[hit] <- mid[hit]; lo[!hit] <- mid[!hit]
        }
        lu <- (lo + hi)/2
      }
      out[j] <- .gpd_from_small_logu(lu, sv[j], xv[j])
    }
  }
  out
}

#' Random generation from the Extended GPD
#'
#' @details Parameters required by \code{type} (the GPD \code{sigma} and \code{xi} and the
#'   carrier parameters) must be supplied: a missing one is an error, whereas
#'   supplied \code{NA} values give \code{NA} results.
#'   For the random generators, \code{n} follows base R: if \code{length(n) > 1}
#'   its length is used, so \code{numeric(0)} gives \code{numeric(0)}.
#'
#'   Where the carrier draw lies within rounding of 1 -- a power
#'   \eqn{U^{1/\kappa}} with \eqn{\kappa > 10^4} (carriers 1 and 6, and
#'   \eqn{U^{2/\kappa}} with \eqn{\kappa > 2 \times 10^4} for carrier 5) or a
#'   draw that rounded to 1 although \eqn{U < 1} -- the GPD quantile is taken
#'   from the log survival probability computed from \eqn{U}, as in
#'   \code{\link{qegpd}}, so the draws stay finite and inside the support.
#' @param n number of samples
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param unifsamp optional uniform samples, recycled to length \code{n}; when
#'   supplied the draws are the GPD quantiles of
#'   \code{r.G(n, ..., unifsamp = unifsamp)}, i.e. inverse-CDF sampling (for
#'   \code{type = 6} the mixture component is still drawn at random; see
#'   \code{\link{r.G}}).
#'
#' @return Random samples
#' @export
regpd <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  .egpd_check_required(type[1], .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta")))
  type <- type[1]
  U <- .egpd_unif(n, unifsamp)
  m <- length(U)
  if (!m) return(numeric(0))
  pr <- rep_len(prob, m); kr <- rep_len(kappa, m); dr <- rep_len(delta, m)
  sv <- rep_len(sigma, m); xv <- rep_len(xi, m)
  ## The incomplete-beta inverse has both tails. Carry its log survivor
  ## into the GPD inverse instead of reconstructing it from a rounded u.
  if (type %in% c(4L, 5L))
    return(qegpd(U, kappa = kr, delta = dr, sigma = sv, xi = xv, type = type))
  draw <- .egpd_rG_core(U, pr, kr, dr, type, direct = FALSE)
  ug <- draw$u
  out <- .qgpd_std(ug, sigma = sv, xi = xv)

  ## 1 - u is lost where the carrier draw u lies within rounding of 1: take
  ## log(1 - u) from U on the log scale there (a type-1 power U^(1/kappa) with
  ## a tiny exponent, or a draw rounded to 1)
  sw <- which(!is.na(ug) & ((is.finite(draw$k1) & draw$k1 > 1e4) |
                              (ug == 1 & U < 1)))
  if (length(sw)) {
    pp <- .egpd_parse_prob(U[sw])
    ls <- numeric(length(sw))
    t1 <- !is.na(draw$k1[sw])
    i1 <- which(t1)
    if (length(i1)) {
      ls[i1] <- .qG_logsurv(pp$p_low[i1], pp$p_upp[i1], type = 1, prob = NA,
                            kappa = draw$k1[sw[i1]], delta = NA,
                            lp_low = pp$lp_low[i1], lp_upp = pp$lp_upp[i1])
    }
    io <- which(!t1)
    if (length(io)) {
      s <- sw[io]
      num <- function(v) suppressWarnings(as.numeric(v[s]))
      ls[io] <- .qG_logsurv(pp$p_low[io], pp$p_upp[io], type = type, prob = num(pr),
                            kappa = num(kr), delta = num(dr),
                            lp_low = pp$lp_low[io], lp_upp = pp$lp_upp[io])
    }
    out[sw] <- .qgpd_from_logsurv(ls, sv[sw], xv[sw])
  }
  tiny <- which(ug < .Machine$double.xmin & U > 0 & U < 1 &
                  is.finite(draw$k1) & is.finite(sv) & sv > 0 & is.finite(xv))
  if (length(tiny))
    out[tiny] <- .gpd_from_small_logu(log(U[tiny])/draw$k1[tiny], sv[tiny], xv[tiny])
  out
}

############################################
##   Discrete Extended GPD               ##
############################################

#' Density of the Discrete Extended GPD
#'
#' Mass is zero off the non-negative integer lattice. The PMF
#' \eqn{H(x+1) - H(x)} is differenced on the log scale in whichever tail of
#' \eqn{H} both values are small (and integrated by Gauss-Legendre quadrature
#' far in the upper tail), so \code{log = TRUE} stays accurate where the mass
#' underflows.
#'
#' @details Parameters required by \code{type} (the GPD \code{sigma} and \code{xi} and the
#'   carrier parameters) must be supplied: a missing one is an error, whereas
#'   supplied \code{NA} values give \code{NA} results.
#' @param x non-negative integer values
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param log logical: return log-PMF?
#'
#' @return PMF values
#' @export
ddiscegpd <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA,
                      type = 1, log = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta")))
  if (min(length(x), length(prob), length(kappa), length(delta),
          length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- x
  n <- max(length(x), length(prob), length(kappa), length(delta), length(sigma), length(xi))
  xv <- rep_len(as.numeric(x), n)
  pv <- rep_len(as.numeric(prob), n)
  kv <- rep_len(as.numeric(kappa), n)
  dv <- rep_len(as.numeric(delta), n)
  sv <- rep_len(as.numeric(sigma), n)
  xiv <- rep_len(as.numeric(xi), n)

  param_ok <- .egpd_carrier_valid(type, pv, kv, dv, n) &
              !is.na(sv) & is.finite(sv) & sv > 0 & !is.na(xiv) & is.finite(xiv)
  out <- rep(NaN, n)
  off <- which(param_ok & !is.na(xv) & (xv < 0 | xv != floor(xv) | is.infinite(xv)))
  out[off] <- 0

  if (log) out[off] <- -Inf
  on_lat <- which(param_ok & !is.na(xv) & is.finite(xv) & xv >= 0 & xv == floor(xv))
  if (length(on_lat)) {
    xl <- xv[on_lat]; pl <- pv[on_lat]; kl <- kv[on_lat]
    dl <- dv[on_lat]; sl <- sv[on_lat]; xil <- xiv[on_lat]
    ## log P(x) = log(H(x + 1) - H(x)), differenced on the log scale in
    ## whichever tail of H the two values are small
    c0 <- .egpd_carrier_log(-.gpd_cdf_surv_std(xl, sl, xil)$log_surv, type, pl, kl, dl)
    c1 <- .egpd_carrier_log(-.gpd_cdf_surv_std(xl + 1, sl, xil)$log_surv, type, pl, kl, dl)
    use_g <- c1$lG <= -log(2)
    lpmf <- ifelse(use_g,
                   c1$lG + .log1mexp(c0$lG - c1$lG),
                   c0$lS + .log1mexp(c1$lS - c0$lS))
    lpmf[is.nan(lpmf) & !is.na(c0$lG) & !is.na(c1$lG)] <- -Inf
    ## Far in the upper tail P(x) can be a tiny fraction of S(x): the
    ## log-space difference then has relative error about
    ## eps * |log S(x)| / |log S(x + 1) - log S(x)|, so where that ratio is
    ## large integrate the (slowly varying) density over [x, x + 1] with
    ## Gauss-Legendre instead.
    dS <- c1$lS - c0$lS
    gl <- which(!use_g & is.finite(c1$lS) & dS > -1e-6 * pmax(1, -c0$lS) & xl >= 1)
    if (length(gl)) {
      rule <- .egpd_gl20()
      m <- length(gl)
      tt <- xl[gl] + rep((1 + rule$x) / 2, each = m)
      rp <- function(v) rep(v[gl], times = 20L)
      At <- -.gpd_cdf_surv_std(tt, rp(sl), rp(xil))$log_surv
      ld <- matrix(.egpd_carrier_logdens(At, type, rp(pl), rp(kl), rp(dl)) +
                     .dgpd_std(tt, sigma = rp(sl), xi = rp(xil), log = TRUE), nrow = m)
      ok_gl <- which(rowSums(!is.finite(ld)) == 0)
      ok_gl <- ok_gl[(.egpd_row_max(ld[ok_gl, , drop = FALSE]) -
                        .egpd_row_min(ld[ok_gl, , drop = FALSE])) <= 3]
      if (length(ok_gl)) {
        lw <- ld[ok_gl, , drop = FALSE] + rep(log(rule$w / 2), each = length(ok_gl))
        mx <- .egpd_row_max(lw)
        lpmf[gl[ok_gl]] <- mx + log(rowSums(exp(lw - mx)))
      }
    }
    out[on_lat] <- if (log) lpmf else exp(lpmf)
  }

  active <- c(list(x, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' CDF of the Discrete Extended GPD
#'
#' Non-integer arguments are rounded down; negative arguments return zero
#' (or one when \code{lower.tail = FALSE}).
#'
#' @details Parameters required by \code{type} (the GPD \code{sigma} and \code{xi} and the
#'   carrier parameters) must be supplied: a missing one is an error, whereas
#'   supplied \code{NA} values give \code{NA} results.
#' @param q non-negative integer quantiles
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le q)},
#'   otherwise \eqn{P(X > q)}
#' @param log.p logical; if TRUE, probabilities are returned on the log scale
#'
#' @return CDF values
#' @export
pdiscegpd <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA,
                      type = 1, lower.tail = TRUE, log.p = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  .egpd_check_required(type[1], .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta")))
  if (min(length(q), length(prob), length(kappa), length(delta),
          length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }
  template <- q
  res <- pegpd(floor(as.numeric(q)) + 1, prob = prob, kappa = kappa, delta = delta,
               sigma = sigma, xi = xi, type = type,
               lower.tail = lower.tail, log.p = log.p)
  active <- c(list(q, sigma, xi), .egpd_active_carrier_args(type[1], prob, kappa, delta))
  res <- .egpd_apply_na_nan(res, active)
  .egpd_restore_dim(res, template)
}

## Correct a continuous-inverse seed at integer jumps. The two exact probability
## endpoints already have their support values; never move them using a rounded
## CDF. Near one, compare an interior linear lower-tail target by its survivor.
.egpd_refine_count_q <- function(q, idx, p, pr, lower.tail, log.p, cdf) {
  idx <- idx[is.finite(q[idx]) & q[idx] < 2^53 &
               is.finite(pr$lp_low[idx]) & is.finite(pr$lp_upp[idx])]
  if (!length(idx)) return(q)
  target <- rep_len(as.numeric(p), length(q))
  far <- lower.tail & !log.p & pr$p_upp < sqrt(.Machine$double.eps)
  reached <- function(at, i) {
    hi <- far[i]
    ans <- logical(length(i))
    if (any(hi))
      ans[hi] <- cdf(at[hi], i[hi], FALSE, FALSE) <= pr$p_upp[i[hi]]
    if (any(!hi)) {
      v <- cdf(at[!hi], i[!hi], lower.tail, log.p)
      t <- target[i[!hi]]
      ans[!hi] <- if (lower.tail) v >= t else v <= t
    }
    ans
  }
  for (step in 1:3) {
    pos <- idx[q[idx] > 0]
    dec <- pos[which(reached(q[pos] - 1, pos))]
    if (!length(dec)) break
    q[dec] <- q[dec] - 1
  }
  for (step in 1:3) {
    inc <- idx[which(!reached(q[idx], idx))]
    if (!length(inc)) break
    q[inc] <- q[inc] + 1
  }
  q
}

#' Quantile function of the Discrete Extended GPD
#'
#' @details Parameters required by \code{type} (the GPD \code{sigma} and \code{xi} and the
#'   carrier parameters) must be supplied: a missing one is an error, whereas
#'   supplied \code{NA} values give \code{NA} results.
#' @param p probabilities
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le q)},
#'   otherwise \eqn{P(X > q)}
#' @param log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}
#'
#' @return Quantile values (non-negative integers)
#' @export
qdiscegpd <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA,
                      type = 1, lower.tail = TRUE, log.p = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  .egpd_check_required(type[1], .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta")))
  if (min(length(p), length(prob), length(kappa), length(delta),
          length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }
  template <- p
  qc <- qegpd(p, prob = prob, kappa = kappa, delta = delta,
              sigma = sigma, xi = xi, type = type,
              lower.tail = lower.tail, log.p = log.p)
  n <- length(qc)
  q <- ceiling(as.numeric(qc)) - 1
  q[which(!is.na(q) & q < 0)] <- 0

  pr <- .egpd_parse_prob(rep_len(p, n), lower.tail = lower.tail, log.p = log.p)
  pv <- rep_len(prob, n); kv <- rep_len(kappa, n); dv <- rep_len(delta, n)
  sv <- rep_len(sigma, n); xv <- rep_len(xi, n)

  ## Zero probability has quantile 0 (tested on the log scale: a log.p input
  ## below about -745 underflows to 0 on the linear scale)
  q[which(pr$valid & is.finite(q) & pr$lp_low == -Inf)] <- 0
  last <- which(pr$valid & !is.na(qc) & xv < 0 & pr$lp_upp == -Inf)
  q[last] <- .gpd_last_count(sv[last], xv[last])

  eval_p <- function(at, idx, lower, logp) {
    pdiscegpd(at, prob = pv[idx], kappa = kv[idx], delta = dv[idx],
              sigma = sv[idx], xi = xv[idx], type = type, lower.tail = lower, log.p = logp)
  }
  q <- .egpd_refine_count_q(q, which(pr$valid), p, pr, lower.tail, log.p, eval_p)

  active <- c(list(p, sigma, xi), .egpd_active_carrier_args(type[1], prob, kappa, delta))
  q <- .egpd_apply_na_nan(q, active)
  .egpd_restore_dim(q, template)
}

#' Random generation from the Discrete Extended GPD
#'
#' @details Parameters required by \code{type} (the GPD \code{sigma} and \code{xi} and the
#'   carrier parameters) must be supplied: a missing one is an error, whereas
#'   supplied \code{NA} values give \code{NA} results.
#'   For the random generators, \code{n} follows base R: if \code{length(n) > 1}
#'   its length is used, so \code{numeric(0)} gives \code{numeric(0)}.
#' @param n number of samples
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param unifsamp optional uniform samples, passed to \code{\link{regpd}}; the
#'   draws are the \code{floor} of the continuous draws. At a finite upper
#'   endpoint, rounding of the continuous draw cannot create an extra cell.
#'
#' @return Random non-negative integer samples
#' @export
rdiscegpd <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  .egpd_check_required(type[1], .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta")))
  u <- .egpd_unif(n, unifsamp)
  n <- length(u)
  if (!n) return(numeric())
  ## Type 6 samples one of two powers. Keep that choice when checking a
  ## rounded draw against its integer cell edge, without redrawing uniforms.
  if (type[1L] == 6L) {
    draw <- .egpd_rG_core(u, rep_len(prob,n), rep_len(kappa,n), rep_len(delta,n), 6L, FALSE)
    kappa <- draw$k1
    kappa[is.na(draw$u)] <- draw$u[is.na(draw$u)]
    type <- 1L
  }
  x <- regpd(n, prob=prob, kappa=kappa, delta=delta, sigma=sigma, xi=xi, type=type, unifsamp=u)
  k <- floor(x)
  shape <- rep_len(xi,n); scale <- rep_len(sigma,n)
  j <- which(is.finite(k) & k > 0 & k < 2^53 & k == x & shape < 0)
  if (length(j)) {
    ls <- pegpd(k[j], prob=rep_len(prob,n)[j], kappa=rep_len(kappa,n)[j],
                delta=rep_len(delta,n)[j], sigma=scale[j], xi=shape[j], type=type,
                lower.tail=FALSE, log.p=TRUE)
    down <- j[log1p(-u[j]) > ls]
    k[down] <- k[down]-1
    at1 <- j[u[j] == 1]
    k[at1] <- .gpd_last_count(scale[at1], shape[at1])
  }
  k
}

############################################
##   Zero-Inflated Discrete Extended GPD  ##
############################################

#' Density of the Zero-Inflated Discrete Extended GPD
#'
#' Numeric arguments, including carrier parameters, are recycled to their
#' common maximum length.
#'
#' @details Parameters required by \code{type} (\code{pi}, the GPD \code{sigma} and \code{xi}
#'   and the carrier parameters) must be supplied: a missing one is an error,
#'   whereas supplied \code{NA} values give \code{NA} results.
#' @param x non-negative integer values
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param log logical: return log-PMF?
#'
#' @return PMF values
#' @export
dzidiscegpd <- function(x, pi = NA, prob = NA, kappa = NA, delta = NA,
                        sigma = NA, xi = NA, type = 1, log = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta", "pi")),
                       pi = TRUE)
  if (min(length(x), length(pi), length(prob), length(kappa),
          length(delta), length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- x
  ly <- max(length(x), length(pi), length(sigma), length(xi),
            length(prob), length(kappa), length(delta))
  xv <- rep_len(as.numeric(x), ly)
  piv <- rep_len(as.numeric(pi), ly)
  pv <- rep_len(prob, ly); kv <- rep_len(kappa, ly); dv <- rep_len(delta, ly)
  sv <- rep_len(sigma, ly); xiv <- rep_len(xi, ly)

  pi_ok <- !is.na(piv) & is.finite(piv) & piv >= 0 & piv <= 1
  base_ok <- .egpd_carrier_valid(type, as.numeric(pv), as.numeric(kv), as.numeric(dv), ly) &
             !is.na(as.numeric(sv)) & is.finite(as.numeric(sv)) & as.numeric(sv) > 0 &
             !is.na(as.numeric(xiv)) & is.finite(as.numeric(xiv))
  all_ok <- pi_ok & base_ok

  out <- rep(NaN, ly)
  off <- which(all_ok & !is.na(xv) & (xv < 0 | xv != floor(xv) | is.infinite(xv)))
  out[off] <- if (log) -Inf else 0

  pi1 <- which(all_ok & !is.na(xv) & piv == 1 & xv >= 0 & xv == floor(xv) & is.finite(xv))
  if (length(pi1)) {
    out[pi1] <- ifelse(xv[pi1] == 0, if (log) 0 else 1, if (log) -Inf else 0)
  }

  reg <- which(all_ok & !is.na(xv) & piv < 1 & xv >= 0 & xv == floor(xv) & is.finite(xv))
  if (length(reg)) {
    z0 <- reg[xv[reg] == 0]
    if (length(z0)) {
      ld0 <- ddiscegpd(0, prob = pv[z0], kappa = kv[z0], delta = dv[z0],
                       sigma = sv[z0], xi = xiv[z0], type = type, log = TRUE)
      ## a log-probability: rounding in the log-sum can leave it just above 0,
      ## and P(0) = 1 under the base law gives exactly 1
      lp0 <- pmin(.logspace_add(log(piv[z0]), log1p(-piv[z0]) + ld0), 0)
      lp0[ld0 == 0] <- 0
      out[z0] <- if (log) lp0 else exp(lp0)
    }
    zp <- reg[xv[reg] > 0]
    if (length(zp)) {
      if (log) {
        ld <- ddiscegpd(xv[zp], prob = pv[zp], kappa = kv[zp], delta = dv[zp],
                        sigma = sv[zp], xi = xiv[zp], type = type, log = TRUE)
        out[zp] <- log1p(-piv[zp]) + ld
      } else {
        d <- ddiscegpd(xv[zp], prob = pv[zp], kappa = kv[zp], delta = dv[zp],
                       sigma = sv[zp], xi = xiv[zp], type = type, log = FALSE)
        out[zp] <- (1 - piv[zp]) * d
      }
    }
  }

  active <- c(list(x, pi, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' CDF of the Zero-Inflated Discrete Extended GPD
#'
#' @details Parameters required by \code{type} (\code{pi}, the GPD \code{sigma} and \code{xi}
#'   and the carrier parameters) must be supplied: a missing one is an error,
#'   whereas supplied \code{NA} values give \code{NA} results.
#' @param q non-negative integer quantiles
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le q)},
#'   otherwise \eqn{P(X > q)}
#' @param log.p logical; if TRUE, probabilities are returned on the log scale
#'
#' @return CDF values
#' @export
pzidiscegpd <- function(q, pi = NA, prob = NA, kappa = NA, delta = NA,
                        sigma = NA, xi = NA, type = 1,
                        lower.tail = TRUE, log.p = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta", "pi")),
                       pi = TRUE)
  if (min(length(q), length(pi), length(prob), length(kappa),
          length(delta), length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- q
  ly <- max(length(q), length(pi), length(sigma), length(xi),
            length(prob), length(kappa), length(delta))
  qv <- rep_len(as.numeric(q), ly)
  piv <- rep_len(as.numeric(pi), ly)
  pi_ok <- !is.na(piv) & is.finite(piv) & piv >= 0 & piv <= 1

  if (lower.tail) {
    lcdf_base <- pdiscegpd(qv, prob = prob, kappa = kappa, delta = delta,
                           sigma = sigma, xi = xi, type = type,
                           lower.tail = TRUE, log.p = TRUE)
    out <- rep(NaN, ly)
    ok <- which(pi_ok & !is.na(lcdf_base))
    if (length(ok)) {
      ## a log-probability: rounding in the log-sum can leave it just above 0,
      ## and a base CDF of exactly 1 gives exactly 1
      lcdf <- pmin(.logspace_add(log(piv[ok]), log1p(-piv[ok]) + lcdf_base[ok]), 0)
      if (log.p) lcdf <- .egpd_zi_lcdf_near1(lcdf, lcdf_base[ok], piv[ok])
      lcdf[lcdf_base[ok] == 0] <- 0
      lcdf[qv[ok] < 0] <- -Inf
      out[ok] <- if (log.p) lcdf else exp(lcdf)
    }
  } else {
    surv_base <- pdiscegpd(qv, prob = prob, kappa = kappa, delta = delta,
                           sigma = sigma, xi = xi, type = type,
                           lower.tail = FALSE, log.p = log.p)
    out <- rep(NaN, ly)
    ok <- which(pi_ok & !is.na(surv_base))
    if (length(ok)) {
      if (log.p) {
        out[ok] <- ifelse(qv[ok] < 0, 0,
                          ifelse(piv[ok] == 1, -Inf, log1p(-piv[ok]) + surv_base[ok]))
      } else {
        out[ok] <- ifelse(qv[ok] < 0, 1, (1 - piv[ok]) * surv_base[ok])
      }
    }
  }

  active <- c(list(q, pi, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' Quantile function of the Zero-Inflated Discrete Extended GPD
#'
#' Returns the smallest support value with CDF at least \code{p}, without
#' shifting probabilities by a fixed tolerance. With \code{pi = 1}, every valid
#' probability has quantile zero; otherwise \code{p = 1} returns the upper endpoint.
#'
#' @details Parameters required by \code{type} (\code{pi}, the GPD \code{sigma} and \code{xi}
#'   and the carrier parameters) must be supplied: a missing one is an error,
#'   whereas supplied \code{NA} values give \code{NA} results.
#' @param p probabilities
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le q)},
#'   otherwise \eqn{P(X > q)}
#' @param log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}
#'
#' @return Quantile values (non-negative integers)
#' @export
qzidiscegpd <- function(p, pi = NA, prob = NA, kappa = NA, delta = NA,
                        sigma = NA, xi = NA, type = 1,
                        lower.tail = TRUE, log.p = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta", "pi")),
                       pi = TRUE)
  if (min(length(p), length(pi), length(prob), length(kappa),
          length(delta), length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- p
  ly <- max(length(p), length(pi), length(sigma), length(xi),
            length(prob), length(kappa), length(delta))
  pr <- .egpd_parse_prob(rep_len(p, ly), lower.tail = lower.tail, log.p = log.p)
  piv <- rep_len(as.numeric(pi), ly)
  pv <- rep_len(prob, ly); kv <- rep_len(kappa, ly); dv <- rep_len(delta, ly)
  sv <- rep_len(sigma, ly); xv <- rep_len(xi, ly)

  pi_ok <- !is.na(piv) & is.finite(piv) & piv >= 0 & piv <= 1
  base_ok <- .egpd_carrier_valid(type, as.numeric(pv), as.numeric(kv), as.numeric(dv), ly) &
             !is.na(as.numeric(sv)) & is.finite(as.numeric(sv)) & as.numeric(sv) > 0 &
             !is.na(as.numeric(xv)) & is.finite(as.numeric(xv))
  valid <- pr$valid & pi_ok & base_ok

  q <- rep(NaN, ly)
  zs <- .egpd_zi_split(pr, piv, valid, log.p)
  at_zero <- which(valid & (piv == 1 | zs$le_pi))
  q[at_zero] <- 0

  above <- which(valid & piv < 1 & !zs$le_pi)
  if (length(above)) {
    ## Use upper-tail probability p_upp / (1 - pi) when in upper half to avoid 1 - p cancellation
    use_upp <- pr$p_low[above] > 0.5
    if (any(!use_upp)) {
      idx <- above[!use_upp]
      lin <- idx[!zs$tiny[idx]]
      if (length(lin)) {
        p_adj <- (pr$p_low[lin] - piv[lin]) / (1 - piv[lin])
        q[lin] <- qdiscegpd(p_adj, prob = pv[lin], kappa = kv[lin], delta = dv[lin],
                            sigma = sv[lin], xi = xv[lin], type = type,
                            lower.tail = TRUE, log.p = FALSE)
      }
      lgs <- idx[zs$tiny[idx]]
      if (length(lgs)) {
        q[lgs] <- qdiscegpd(.egpd_zi_lp_adj(pr$lp_low[lgs], piv[lgs]), prob = pv[lgs],
                            kappa = kv[lgs], delta = dv[lgs], sigma = sv[lgs],
                            xi = xv[lgs], type = type, lower.tail = TRUE, log.p = TRUE)
      }
    }
    if (any(use_upp)) {
      idx <- above[use_upp]
      lpu_adj <- pmin(pr$lp_upp[idx] - log1p(-piv[idx]), 0)
      q[idx] <- qdiscegpd(lpu_adj, prob = pv[idx], kappa = kv[idx], delta = dv[idx],
                          sigma = sv[idx], xi = xv[idx], type = type,
                          lower.tail = FALSE, log.p = TRUE)
    }

    eval_pz <- function(at, idx, lower, logp) {
      pzidiscegpd(at, pi = piv[idx], prob = pv[idx], kappa = kv[idx],
                  delta = dv[idx], sigma = sv[idx], xi = xv[idx],
                  type = type, lower.tail = lower, log.p = logp)
    }
    q <- .egpd_refine_count_q(q, above, p, pr, lower.tail, log.p, eval_pz)
  }

  active <- c(list(p, pi, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  q <- .egpd_apply_na_nan(q, active)
  .egpd_restore_dim(q, template)
}

#' Random generation from the Zero-Inflated Discrete Extended GPD
#'
#' @details Parameters required by \code{type} (\code{pi}, the GPD \code{sigma} and \code{xi}
#'   and the carrier parameters) must be supplied: a missing one is an error,
#'   whereas supplied \code{NA} values give \code{NA} results.
#'   For the random generators, \code{n} follows base R: if \code{length(n) > 1}
#'   its length is used, so \code{numeric(0)} gives \code{numeric(0)}.
#' @param n number of samples
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param unifsamp optional uniform samples; when supplied the draws are
#'   \code{qzidiscegpd(unifsamp, ...)}, i.e. inverse-CDF sampling of the whole
#'   zero-inflated law, so the zero indicator is also determined by
#'   \code{unifsamp} and no other random numbers are used. \code{unifsamp}
#'   and the parameters are recycled to (or truncated at) length \code{n}.
#'
#' @return Random non-negative integer samples
#' @export
rzidiscegpd <- function(n, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  .egpd_check_required(type[1], .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta", "pi")),
                       pi = TRUE)
  n <- .egpd_rn(n)
  if (!is.null(unifsamp) && length(unifsamp) == 0L) return(numeric(0))
  if (n == 0L) return(numeric(0))
  if (!is.null(unifsamp)) {
    ## inverse-CDF sampling of the whole zero-inflated law: the zero indicator
    ## is determined by unifsamp too, so the draws are reproducible.  As in
    ## base R generators, parameters are recycled to (or truncated at) length n.
    rl <- function(v) rep_len(v, n)
    return(qzidiscegpd(rep_len(as.numeric(unifsamp), n), pi = rl(pi), prob = rl(prob),
                       kappa = rl(kappa), delta = rl(delta), sigma = rl(sigma), xi = rl(xi),
                       type = type))
  }
  z <- rbinom(n, size = 1, prob = pi)
  y <- (1 - z) * rdiscegpd(n, prob = prob, kappa = kappa, delta = delta,
                           sigma = sigma, xi = xi, type = type)
  return(y)
}

############################################
##  Continuous Zero-Inflated EGPD         ##
############################################

#' Density of the Zero-Inflated Extended GPD
#'
#' Numeric arguments, including carrier parameters, are recycled to their
#' common maximum length.
#'
#' Continuous zero-inflated EGPD density:
#' mass \eqn{\pi} at zero and density \eqn{(1-\pi)f_{EGPD}(x)} for \eqn{x > 0}.
#'
#' @details Parameters required by \code{type} (\code{pi}, the GPD \code{sigma} and \code{xi}
#'   and the carrier parameters) must be supplied: a missing one is an error,
#'   whereas supplied \code{NA} values give \code{NA} results.
#' @param x non-negative values
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param log logical: return log-density?
#'
#' @return Density values
#' @export
dziegpd <- function(x, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, log = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta", "pi")),
                       pi = TRUE)
  if (min(length(x), length(pi), length(prob), length(kappa),
          length(delta), length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- x
  ly <- max(length(x), length(pi), length(sigma), length(xi), length(prob), length(kappa), length(delta))
  xv <- rep_len(as.numeric(x), ly)
  piv <- rep_len(as.numeric(pi), ly)
  pv <- rep_len(prob, ly); kv <- rep_len(kappa, ly); dv <- rep_len(delta, ly)
  sv <- rep_len(sigma, ly); xiv <- rep_len(xi, ly)

  pi_ok <- !is.na(piv) & is.finite(piv) & piv >= 0 & piv <= 1
  base_ok <- .egpd_carrier_valid(type, as.numeric(pv), as.numeric(kv), as.numeric(dv), ly) &
             !is.na(as.numeric(sv)) & is.finite(as.numeric(sv)) & as.numeric(sv) > 0 &
             !is.na(as.numeric(xiv)) & is.finite(as.numeric(xiv))
  all_ok <- pi_ok & base_ok

  out <- rep(NaN, ly)
  neg <- which(all_ok & !is.na(xv) & xv < 0)
  out[neg] <- if (log) -Inf else 0

  z0 <- which(all_ok & !is.na(xv) & xv == 0)
  if (length(z0)) {
    out[z0] <- if (log) ifelse(piv[z0] > 0, log(piv[z0]), -Inf) else piv[z0]
  }

  zp <- which(all_ok & !is.na(xv) & xv > 0)
  if (length(zp)) {
    pi1 <- zp[piv[zp] == 1]
    if (length(pi1)) out[pi1] <- if (log) -Inf else 0
    reg <- zp[piv[zp] < 1]
    if (length(reg)) {
      if (log) {
        out[reg] <- log1p(-piv[reg]) +
          degpd_density(xv[reg], prob = pv[reg], kappa = kv[reg],
                        delta = dv[reg], sigma = sv[reg], xi = xiv[reg],
                        type = type, log = TRUE)
      } else {
        out[reg] <- (1 - piv[reg]) *
          degpd_density(xv[reg], prob = pv[reg], kappa = kv[reg],
                        delta = dv[reg], sigma = sv[reg], xi = xiv[reg],
                        type = type, log = FALSE)
      }
    }
  }

  active <- c(list(x, pi, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' CDF of the Zero-Inflated Extended GPD
#'
#' @details Parameters required by \code{type} (\code{pi}, the GPD \code{sigma} and \code{xi}
#'   and the carrier parameters) must be supplied: a missing one is an error,
#'   whereas supplied \code{NA} values give \code{NA} results.
#' @param q quantiles
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le q)},
#'   otherwise \eqn{P(X > q)}
#' @param log.p logical; if TRUE, probabilities are returned on the log scale
#'
#' @return CDF values
#' @export
pziegpd <- function(q, pi = NA, prob = NA, kappa = NA, delta = NA,
                    sigma = NA, xi = NA, type = 1,
                    lower.tail = TRUE, log.p = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta", "pi")),
                       pi = TRUE)
  if (min(length(q), length(pi), length(prob), length(kappa),
          length(delta), length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- q
  ly <- max(length(q), length(pi), length(sigma), length(xi),
            length(prob), length(kappa), length(delta))
  qv <- rep_len(as.numeric(q), ly)
  piv <- rep_len(as.numeric(pi), ly)
  pi_ok <- !is.na(piv) & is.finite(piv) & piv >= 0 & piv <= 1

  if (lower.tail) {
    lcdf_base <- pegpd(qv, prob = prob, kappa = kappa, delta = delta,
                       sigma = sigma, xi = xi, type = type,
                       lower.tail = TRUE, log.p = TRUE)
    out <- rep(NaN, ly)
    ok <- which(pi_ok & !is.na(lcdf_base))
    if (length(ok)) {
      ## a log-probability: rounding in the log-sum can leave it just above 0,
      ## and a base CDF of exactly 1 gives exactly 1
      lcdf <- pmin(.logspace_add(log(piv[ok]), log1p(-piv[ok]) + lcdf_base[ok]), 0)
      if (log.p) lcdf <- .egpd_zi_lcdf_near1(lcdf, lcdf_base[ok], piv[ok])
      lcdf[lcdf_base[ok] == 0] <- 0
      lcdf[qv[ok] < 0] <- -Inf
      out[ok] <- if (log.p) lcdf else exp(lcdf)
    }
  } else {
    surv_base <- pegpd(qv, prob = prob, kappa = kappa, delta = delta,
                       sigma = sigma, xi = xi, type = type,
                       lower.tail = FALSE, log.p = log.p)
    out <- rep(NaN, ly)
    ok <- which(pi_ok & !is.na(surv_base))
    if (length(ok)) {
      if (log.p) {
        out[ok] <- ifelse(qv[ok] < 0, 0,
                          ifelse(piv[ok] == 1, -Inf, log1p(-piv[ok]) + surv_base[ok]))
      } else {
        out[ok] <- ifelse(qv[ok] < 0, 1, (1 - piv[ok]) * surv_base[ok])
      }
    }
  }

  active <- c(list(q, pi, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  out <- .egpd_apply_na_nan(out, active)
  .egpd_restore_dim(out, template)
}

#' Quantile function of the Zero-Inflated Extended GPD
#'
#' Probabilities at or below \code{pi} return zero, with no epsilon subtraction.
#' With \code{pi = 1}, all valid probabilities return zero. Otherwise the upper
#' endpoint is returned at \code{p = 1}, including infinity for unbounded support.
#'
#' @details Parameters required by \code{type} (\code{pi}, the GPD \code{sigma} and \code{xi}
#'   and the carrier parameters) must be supplied: a missing one is an error,
#'   whereas supplied \code{NA} values give \code{NA} results.
#' @param p probabilities
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le q)},
#'   otherwise \eqn{P(X > q)}
#' @param log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}
#'
#' @return Quantile values
#' @export
qziegpd <- function(p, pi = NA, prob = NA, kappa = NA, delta = NA,
                    sigma = NA, xi = NA, type = 1,
                    lower.tail = TRUE, log.p = FALSE) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  .egpd_check_required(type, .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta", "pi")),
                       pi = TRUE)
  if (min(length(p), length(pi), length(prob), length(kappa),
          length(delta), length(sigma), length(xi)) == 0L) {
    return(numeric(0))
  }

  template <- p
  ly <- max(length(p), length(pi), length(sigma), length(xi),
            length(prob), length(kappa), length(delta))
  pr <- .egpd_parse_prob(rep_len(p, ly), lower.tail = lower.tail, log.p = log.p)
  piv <- rep_len(as.numeric(pi), ly)
  pv <- rep_len(prob, ly); kv <- rep_len(kappa, ly); dv <- rep_len(delta, ly)
  sv <- rep_len(sigma, ly); xv <- rep_len(xi, ly)

  pi_ok <- !is.na(piv) & is.finite(piv) & piv >= 0 & piv <= 1
  base_ok <- .egpd_carrier_valid(type, as.numeric(pv), as.numeric(kv), as.numeric(dv), ly) &
             !is.na(as.numeric(sv)) & is.finite(as.numeric(sv)) & as.numeric(sv) > 0 &
             !is.na(as.numeric(xv)) & is.finite(as.numeric(xv))
  valid <- pr$valid & pi_ok & base_ok

  q <- rep(NaN, ly)
  zs <- .egpd_zi_split(pr, piv, valid, log.p)
  at_zero <- which(valid & (piv == 1 | zs$le_pi))
  q[at_zero] <- 0

  above <- which(valid & piv < 1 & !zs$le_pi)
  if (length(above)) {
    use_upp <- pr$p_low[above] > 0.5
    if (any(!use_upp)) {
      idx <- above[!use_upp]
      lin <- idx[!zs$tiny[idx]]
      if (length(lin)) {
        p_adj <- (pr$p_low[lin] - piv[lin]) / (1 - piv[lin])
        q[lin] <- qegpd(p_adj, prob = pv[lin], kappa = kv[lin], delta = dv[lin],
                        sigma = sv[lin], xi = xv[lin], type = type,
                        lower.tail = TRUE, log.p = FALSE)
      }
      lgs <- idx[zs$tiny[idx]]
      if (length(lgs)) {
        q[lgs] <- qegpd(.egpd_zi_lp_adj(pr$lp_low[lgs], piv[lgs]), prob = pv[lgs],
                        kappa = kv[lgs], delta = dv[lgs], sigma = sv[lgs],
                        xi = xv[lgs], type = type, lower.tail = TRUE, log.p = TRUE)
      }
    }
    if (any(use_upp)) {
      idx <- above[use_upp]
      lpu_adj <- pmin(pr$lp_upp[idx] - log1p(-piv[idx]), 0)
      q[idx] <- qegpd(lpu_adj, prob = pv[idx], kappa = kv[idx], delta = dv[idx],
                      sigma = sv[idx], xi = xv[idx], type = type,
                      lower.tail = FALSE, log.p = TRUE)
    }
  }

  active <- c(list(p, pi, sigma, xi), .egpd_active_carrier_args(type, prob, kappa, delta))
  q <- .egpd_apply_na_nan(q, active)
  .egpd_restore_dim(q, template)
}

#' Random generation from the Zero-Inflated Extended GPD
#'
#' @details Parameters required by \code{type} (\code{pi}, the GPD \code{sigma} and \code{xi}
#'   and the carrier parameters) must be supplied: a missing one is an error,
#'   whereas supplied \code{NA} values give \code{NA} results.
#'   For the random generators, \code{n} follows base R: if \code{length(n) > 1}
#'   its length is used, so \code{numeric(0)} gives \code{numeric(0)}.
#' @param n number of samples
#' @param pi zero-inflation probability
#' @param prob mixing probability (type 6)
#' @param kappa shape parameter for G transformation
#' @param delta shape parameter for G transformation (types 4-6)
#' @param sigma GPD scale parameter
#' @param xi GPD shape parameter
#' @param type integer 1-6 specifying G type
#'
#' @return Random non-negative samples
#' @export
rziegpd <- function(n, pi = NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  if (length(type) == 0L || is.na(type[1]) || !type[1] %in% 1:6) {
    stop("Invalid `type' argument")
  }
  .egpd_check_required(type[1], .egpd_missing_args(c("sigma", "xi", "prob", "kappa", "delta", "pi")),
                       pi = TRUE)
  n <- .egpd_rn(n)
  if (n == 0L) return(numeric(0))
  z <- rbinom(n, size = 1, prob = pi)
  y <- (1 - z) * regpd(n, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type)
  return(y)
}
