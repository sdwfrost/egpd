r2_likdata <- function(y, npar, xi.max = Inf) {
  list(
    X = rep(list(matrix(1, length(y), 1)), npar),
    y = cbind(y),
    idpars = seq_len(npar),
    dupid = 0L,
    duplicate = 0L,
    offsets = rep(list(numeric()), npar),
    censored = FALSE,
    xi.max = xi.max
  )
}

r2_richardson_fd12 <- function(f_vec, b, step = 1e-3) {
  one <- function(h) {
    n <- length(b)
    f0 <- f_vec(b)
    shift <- function(i, d) { z <- b; z[i] <- z[i] + d; z }
    g <- vapply(seq_len(n), function(i)
      (f_vec(shift(i, h)) - f_vec(shift(i, -h))) / (2 * h), f0)
    if (length(f0) == 1L) g <- matrix(g, nrow = 1L)
    H <- list()
    for (i in seq_len(n)) for (j in i:n) {
      if (i == j) {
        z <- (f_vec(shift(i, h)) - 2 * f0 + f_vec(shift(i, -h))) / h^2
      } else {
        pp <- pm <- mp <- mm <- b
        pp[c(i, j)] <- pp[c(i, j)] + h
        mm[c(i, j)] <- mm[c(i, j)] - h
        pm[i] <- pm[i] + h; pm[j] <- pm[j] - h
        mp[i] <- mp[i] - h; mp[j] <- mp[j] + h
        z <- (f_vec(pp) - f_vec(pm) - f_vec(mp) + f_vec(mm)) / (4 * h^2)
      }
      H[[length(H) + 1L]] <- z
    }
    cbind(g, do.call(cbind, H))
  }
  (4 * one(step / 2) - one(step)) / 3
}

test_that("R2-1: identity and log-link DEGPD retain exact d0 and d12 when log F(y+1) < -745", {
  skip_on_cran()
  # Exact reproducer: y = 0, sigma = 100, xi = 0.1, kappa = 200 -> NLL = 922.132663388948
  y0 <- c(0, 1)
  cases <- list(
    `1` = c(log(100), 0.1, log(200)),
    `2` = c(log(100), 0.1, log(200), log(50), 0.3),
    `3` = c(log(1e200), 0.1, log(1.5)),
    `4` = c(log(100), 0.1, log(1.5), log(400)),
    `5` = c(log(100), 0.1, log(2000)),
    `6` = c(log(100), 0.1, log(400))
  )

  # Exact m=1 analytical check at y = 0, sigma = 100, xi = 0.1, kappa = 200
  ld1 <- r2_likdata(0, 3L)
  b1 <- cases[["1"]]
  true_nll1 <- -200 * log(-expm1(-log1p(0.1 / 100) / 0.1))
  expect_equal(egpd:::.degpd1idfns$d0(b1, ld1), true_nll1, tolerance = 1e-13)
  expect_equal(egpd:::.degpd1fns$d0(c(b1[1], log(b1[2]), b1[3]), ld1),
               true_nll1, tolerance = 1e-13)

  for (m in 1:6) {
    b_id <- cases[[as.character(m)]]
    npar <- length(b_id)
    fns_id <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    fns_lg <- get(paste0(".degpd", m, "fns"), asNamespace("egpd"))

    for (yy in y0) {
      ld <- r2_likdata(yy, npar)
      d0_id <- fns_id$d0(b_id, ld)
      d12_id <- fns_id$d120(b_id, ld)
      expect_true(is.finite(d0_id) && d0_id > 700,
                  info = sprintf("id m=%d y=%d d0=%g", m, yy, d0_id))
      expect_true(all(is.finite(d12_id)),
                  info = sprintf("id m=%d y=%d d12 non-finite", m, yy))

      want_id <- r2_richardson_fd12(function(p) fns_id$d0(p, ld), b_id, step = 1e-3)
      expect_equal(unname((d12_id - want_id) / (1 + abs(want_id))),
                   matrix(0, 1, ncol(d12_id)), tolerance = 2e-6,
                   info = sprintf("id FD m=%d y=%d", m, yy))

      b_lg <- b_id
      b_lg[2] <- log(b_id[2])
      d0_lg <- fns_lg$d0(b_lg, ld)
      d12_lg <- fns_lg$d120(b_lg, ld)
      expect_equal(d0_lg, d0_id, tolerance = 1e-13,
                   info = sprintf("log vs id d0 m=%d y=%d", m, yy))
      want_lg <- r2_richardson_fd12(function(p) fns_lg$d0(p, ld), b_lg, step = 1e-3)
      expect_equal(unname((d12_lg - want_lg) / (1 + abs(want_lg))),
                   matrix(0, 1, ncol(d12_lg)), tolerance = 2e-6,
                   info = sprintf("log FD m=%d y=%d", m, yy))
    }

    # Also check the exact parameter point (sigma=100, xi=0.1, kappa=200) across all m=1..6
    b_200 <- switch(m,
      c(log(100), 0.1, log(200)),
      c(log(100), 0.1, log(200), log(50), 0.3),
      c(log(100), 0.1, log(200)),
      c(log(100), 0.1, log(1.5), log(200)),
      c(log(100), 0.1, log(200)),
      c(log(100), 0.1, log(200)))
    ld0 <- r2_likdata(0, npar)
    got_200 <- fns_id$d120(b_200, ld0)
    want_200 <- r2_richardson_fd12(function(p) fns_id$d0(p, ld0), b_200, step = 1e-3)
    expect_equal(unname((got_200 - want_200) / (1 + abs(want_200))),
                 matrix(0, 1, ncol(got_200)), tolerance = 2e-6,
                 info = sprintf("m=%d y=0 sigma=100 kappa=200", m))
  }

  # Extreme Newton trial steps in lkappa for quadrature carriers m = 5, 6 never throw
  for (m in c(5L, 6L)) {
    fns_lg <- get(paste0(".degpd", m, "fns"), asNamespace("egpd"))
    ld_ext <- r2_likdata(c(0, 5, 20), 3L)
    for (lk in c(-50, 30, 800)) {
      expect_true(is.finite(fns_lg$d0(c(log(8), -2, lk), ld_ext)))
      expect_true(all(is.finite(fns_lg$d120(c(log(8), -2, lk), ld_ext))))
    }
  }
})

test_that("R2-1b: carrier 3 and 4 avoid small-delta and large-delta cancellation", {
  skip_on_cran()
  y <- c(0, 1, 6, 20)
  for (m in c(3L, 4L)) {
    npar <- if (m == 3L) 3L else 4L
    fns_id <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    ld <- r2_likdata(y, npar)
    for (ldelta in c(-30, -12, -6, 6, 25, 60)) {
      for (xi in c(-0.1, 0, 0.2)) {
        b <- c(log(4), xi, ldelta, if (m == 4L) log(1.8))
        d0 <- fns_id$d0(b, ld)
        d12 <- fns_id$d120(b, ld)
        expect_true(is.finite(d0), info = sprintf("m=%d ldelta=%g xi=%g", m, ldelta, xi))
        expect_true(all(is.finite(d12)), info = sprintf("m=%d ldelta=%g xi=%g", m, ldelta, xi))
        got <- colSums(d12)
        want <- drop(r2_richardson_fd12(function(p) fns_id$d0(p, ld), b, step = 1e-3))
        expect_equal(unname((got - want) / (1 + abs(want))),
                     rep(0, length(got)), tolerance = 5e-6,
                     info = sprintf("m=%d ldelta=%g xi=%g", m, ldelta, xi))
      }
    }
  }
})

test_that("R2-2: continuous EGPD d12 and d34 guard out-of-support and y = 0 cleanly", {
  theta_for <- function(m, xi = -0.25) {
    switch(m,
      c(log(2), xi, log(1.5)),
      c(log(2), xi, log(1.2), log(0.8), qlogis(0.35)),
      c(log(2), xi, log(1.4)),
      c(log(2), xi, log(1.4), log(1.5)),
      c(log(2), xi, log(2.0)),
      c(log(2), xi, log(2.0)))
  }

  # 1. Out-of-support: y < 0 and 1 + xi*y/sigma <= 0 (sigma = 2, xi = -0.25 => upper = 8)
  y_mix <- c(-1.0, 1.5, 8.0, 10.5, 2.5)
  in_sup <- c(FALSE, TRUE, FALSE, FALSE, TRUE)
  for (m in 1:6) {
    b <- theta_for(m, xi = -0.25)
    npar <- length(b)
    ld <- r2_likdata(y_mix, npar)
    ld_in <- r2_likdata(y_mix[in_sup], npar)
    d0_fn <- get(paste0(".egpd", m, ".d0"), asNamespace("egpd"))
    d12_fn <- get(paste0(".egpd", m, ".d12"), asNamespace("egpd"))

    expect_equal(d0_fn(b, ld), d0_fn(b, ld_in) + 3 * 1e20, tolerance = 1e-12)
    M12 <- d12_fn(b, ld)
    expect_true(all(is.finite(M12)), info = paste("d12 finite OOS m=", m))
    expect_equal(M12[!in_sup, , drop = FALSE],
                 matrix(0, sum(!in_sup), ncol(M12)),
                 info = paste("d12 zero OOS m=", m))
    want_in <- r2_richardson_fd12(
      function(p) vapply(y_mix[in_sup], function(yy) d0_fn(p, r2_likdata(yy, npar)), 0),
      b, step = 1e-4)
    expect_equal(unname((M12[in_sup, , drop = FALSE] - want_in) / (1 + abs(want_in))),
                 matrix(0, sum(in_sup), ncol(M12)), tolerance = 2e-5,
                 info = paste("d12 in-support m=", m))

    if (m <= 4) {
      b34 <- if (m == 2) c(log(2), -0.25, log(1.2), log(2.0), qlogis(0.35)) else b
      d34_fn <- get(paste0(".egpd", m, ".d34"), asNamespace("egpd"))
      M34 <- d34_fn(b34, ld)
      expect_true(all(is.finite(M34)), info = paste("d34 finite OOS m=", m))
      expect_equal(M34[!in_sup, , drop = FALSE],
                   matrix(0, sum(!in_sup), ncol(M34)),
                   info = paste("d34 zero OOS m=", m))
    }
  }

  # 2. Boundary y = 0: finite positive density cases vs zero/infinite density penalties
  for (xi in c(-0.2, 0.0, 0.3)) {
    # m = 1, kappa = 1 -> NLL = lpsi
    ld1 <- r2_likdata(0, 3L)
    b1 <- c(log(2.5), xi, 0.0)
    expect_equal(egpd:::.egpd1.d0(b1, ld1), log(2.5))
    expect_equal(as.numeric(egpd:::.egpd1.d12(b1, ld1)), c(1, rep(0, 8)))
    expect_equal(as.numeric(egpd:::.egpd1.d34(b1, ld1)), rep(0, 25))
    # m = 1, kappa != 1 -> clean penalty 1e20, zero derivatives, never NaN
    for (lk in c(log(0.5), log(2.0))) {
      expect_equal(egpd:::.egpd1.d0(c(log(2.5), xi, lk), ld1), 1e20)
      expect_equal(as.numeric(egpd:::.egpd1.d12(c(log(2.5), xi, lk), ld1)), rep(0, 9))
      expect_equal(as.numeric(egpd:::.egpd1.d34(c(log(2.5), xi, lk), ld1)), rep(0, 25))
    }

    # m = 2, kappa1 = 1 -> NLL = lpsi - log(p)
    ld2 <- r2_likdata(0, 5L)
    lp_val <- 0.4; p_prob <- plogis(lp_val); q_prob <- plogis(-lp_val)
    b2 <- c(log(2.5), xi, 0.0, log(0.8), lp_val)
    expect_equal(egpd:::.egpd2.d0(b2, ld2), log(2.5) - log(p_prob), tolerance = 1e-14)
    got2 <- as.numeric(egpd:::.egpd2.d12(b2, ld2))
    want2 <- rep(0, 20); want2[1] <- 1; want2[5] <- -q_prob; want2[20] <- p_prob * q_prob
    expect_equal(got2, want2, tolerance = 1e-14)
    b2_34 <- c(log(2.5), xi, 0.0, log(1.8), lp_val)
    got2_34 <- as.numeric(egpd:::.egpd2.d34(b2_34, ld2))
    expect_true(all(is.finite(got2_34)))
    expect_equal(got2_34[35], p_prob * q_prob * (1 - 2 * p_prob), tolerance = 1e-14)
    expect_equal(got2_34[105], p_prob * q_prob * (1 - 6 * p_prob * q_prob), tolerance = 1e-14)
    # m = 2, kappa1 != 1 -> penalty
    expect_equal(egpd:::.egpd2.d0(c(log(2.5), xi, 0.3, log(0.8), lp_val), ld2), 1e20)
    expect_equal(as.numeric(egpd:::.egpd2.d12(c(log(2.5), xi, 0.3, log(0.8), lp_val), ld2)), rep(0, 20))

    # m = 3 -> g3(0) = 0 -> penalty 1e20 and zero derivatives
    ld3 <- r2_likdata(0, 3L)
    expect_equal(egpd:::.egpd3.d0(c(log(2.5), xi, 0.3), ld3), 1e20)
    expect_equal(as.numeric(egpd:::.egpd3.d12(c(log(2.5), xi, 0.3), ld3)), rep(0, 9))
    expect_equal(as.numeric(egpd:::.egpd3.d34(c(log(2.5), xi, 0.3), ld3)), rep(0, 25))

    # m = 4, kappa = 1 (lkappa = 0) -> NLL = lpsi - 0.5*(log1p(delta) - log(2))
    ld4 <- r2_likdata(0, 4L)
    b4 <- c(log(2.5), xi, 0.6, 0.0)
    q_d <- plogis(0.6)
    expect_equal(egpd:::.egpd4.d0(b4, ld4), log(2.5) - 0.5 * (log1p(exp(0.6)) - log(2)), tolerance = 1e-14)
    got4 <- as.numeric(egpd:::.egpd4.d12(b4, ld4))
    want4 <- rep(0, 14); want4[1] <- 1; want4[3] <- -0.5 * q_d; want4[12] <- -0.5 * q_d * (1 - q_d)
    expect_equal(got4, want4, tolerance = 1e-14)
    got4_34 <- as.numeric(egpd:::.egpd4.d34(b4, ld4))
    expect_equal(got4_34[17], -0.5 * q_d * (1 - q_d) * (1 - 2 * q_d), tolerance = 1e-14)
    expect_equal(got4_34[51], -0.5 * q_d * (1 - q_d) * (1 - 6 * q_d * (1 - q_d)), tolerance = 1e-14)
    # m = 4, kappa != 1 -> penalty
    expect_equal(egpd:::.egpd4.d0(c(log(2.5), xi, 0.6, log(2)), ld4), 1e20)
    expect_equal(as.numeric(egpd:::.egpd4.d12(c(log(2.5), xi, 0.6, log(2)), ld4)), rep(0, 14))
    expect_equal(as.numeric(egpd:::.egpd4.d34(c(log(2.5), xi, 0.6, log(2)), ld4)), rep(0, 55))

    # m = 5 and m = 6 at y = 0: smooth positive density for any kappa > 0
    for (m in c(5L, 6L)) {
      ldm <- r2_likdata(0, 3L)
      bm <- c(log(2.5), xi, log(1.8))
      d0m <- get(paste0(".egpd", m, ".d0"), asNamespace("egpd"))
      d12m <- get(paste0(".egpd", m, ".d12"), asNamespace("egpd"))
      gotm <- d12m(bm, ldm)
      wantm <- r2_richardson_fd12(function(p) d0m(p, ldm), bm, step = 1e-4)
      expect_true(is.finite(d0m(bm, ldm)))
      expect_equal(unname((gotm - wantm) / (1 + abs(wantm))),
                   matrix(0, 1, 9), tolerance = 2e-6,
                   info = sprintf("m=%d y=0 xi=%g", m, xi))
    }
  }

  # 3. Tiny positive y (H(y) near 0) across m = 1..6
  for (m in 1:6) {
    b <- theta_for(m, xi = 0.25)
    ld <- r2_likdata(1e-15, length(b))
    d0_fn <- get(paste0(".egpd", m, ".d0"), asNamespace("egpd"))
    d12_fn <- get(paste0(".egpd", m, ".d12"), asNamespace("egpd"))
    got <- d12_fn(b, ld)
    want <- r2_richardson_fd12(function(p) d0_fn(p, ld), b, step = 1e-4)
    expect_true(all(is.finite(got)), info = paste("tiny y d12 m=", m))
    expect_equal(unname((got - want) / (1 + abs(want))),
                 matrix(0, 1, ncol(got)), tolerance = 5e-5,
                 info = paste("tiny y FD m=", m))
  }
})

test_that("R2-2b: continuous EGPD model 4 avoids pbeta underflow at large delta", {
  y <- c(0.8, 2.0, 5.0)
  ld4 <- r2_likdata(y, 4L)
  ld1 <- r2_likdata(y, 3L)
  for (ldelta in c(4, 6, 12, 40)) {
    b4 <- c(log(1.5), 0.15, ldelta, log(3.4))
    d <- exp(ldelta); k <- exp(b4[4])
    z <- y / exp(b4[1]); A <- log1p(b4[2] * z) / b4[2]
    G3_exact <- 1 - exp(-A) * (1 + -expm1(-d * A) / d)
    logg4_exact <- log1p(exp(-ldelta)) + log(-expm1(-d * A)) +
      b4[4] - log(2) + (k / 2 - 1) * log(G3_exact)
    want_d0 <- sum(b4[1] + (1 + b4[2]) * A - logg4_exact)
    expect_equal(egpd:::.egpd4.d0(b4, ld4), want_d0, tolerance = 1e-12,
                 info = paste("ldelta=", ldelta))
    got12 <- colSums(egpd:::.egpd4.d12(b4, ld4))
    want12 <- drop(r2_richardson_fd12(function(p) egpd:::.egpd4.d0(p, ld4), b4, step = 1e-4))
    expect_equal(unname((got12 - want12) / (1 + abs(want12))),
                 rep(0, 14), tolerance = 2e-5,
                 info = paste("d12 ldelta=", ldelta))
  }
  # As delta -> Inf, EGPD4(sigma, xi, delta, kappa) -> EGPD1(sigma, xi, kappa/2)
  b4_inf <- c(log(1.5), 0.15, 45, log(3.4))
  b1_lim <- c(log(1.5), 0.15, log(3.4 / 2))
  expect_equal(egpd:::.egpd4.d0(b4_inf, ld4), egpd:::.egpd1.d0(b1_lim, ld1), tolerance = 1e-12)
})

test_that("R2-3a: log-link DEGPD (m = 1..6) is exact down to xi = 0 (lxi = -800)", {
  y <- c(0, 1, 5, 18)
  for (m in 1:6) {
    npar <- c(3L, 5L, 3L, 4L, 3L, 3L)[m]
    fns_lg <- get(paste0(".degpd", m, "fns"), asNamespace("egpd"))
    fns_id <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    ld <- r2_likdata(y, npar)

    extra <- switch(m,
      log(1.6),
      c(log(1.3), log(0.9), 0.25),
      log(1.4),
      c(log(1.4), log(1.9)),
      log(1.7),
      log(1.7))

    for (lxi in c(log(0.25), log(1e-4), log(1e-8), log(1e-12), log(1e-16), -100, -800)) {
      b_lg <- c(log(3.5), lxi, extra)
      b_id <- c(log(3.5), exp(lxi), extra)
      d0_lg <- fns_lg$d0(b_lg, ld)
      d0_id <- fns_id$d0(b_id, ld)
      expect_equal(d0_lg, d0_id, tolerance = 1e-12,
                   info = sprintf("m=%d lxi=%g", m, lxi))

      M_lg <- fns_lg$d120(b_lg, ld)
      expect_true(all(is.finite(M_lg)), info = sprintf("m=%d lxi=%g finite", m, lxi))

      if (lxi >= log(1e-12)) {
        want <- r2_richardson_fd12(
          function(p) vapply(y, function(yy) fns_lg$d0(p, r2_likdata(yy, npar)), 0),
          b_lg, step = 1e-3)
        expect_equal(unname((M_lg - want) / (1 + abs(want))),
                     matrix(0, nrow(M_lg), ncol(M_lg)), tolerance = 2e-5,
                     info = sprintf("m=%d lxi=%g FD", m, lxi))
        # Relative sign/accuracy check specifically on the lxi score column!
        M_id <- fns_id$d120(b_id, ld)
        expect_equal(M_lg[, 2], exp(lxi) * M_id[, 2], tolerance = 1e-12,
                     info = sprintf("m=%d lxi=%g exact chain", m, lxi))
      } else if (lxi <= -750) {
        expect_equal(M_lg[, 2], rep(0, length(y)))
      }
    }

    # Bounded xi.max = 0.5 link across extreme predictors eta_xi in [-800, 800]
    ld_xm <- r2_likdata(y, npar, xi.max = 0.5)
    for (eta_xi in c(-800, -25, -1, 1, 800)) {
      b_xm <- c(log(3.5), eta_xi, extra)
      expect_true(is.finite(fns_lg$d0(b_xm, ld_xm)))
      M_xm <- fns_lg$d120(b_xm, ld_xm)
      expect_true(all(is.finite(M_xm)),
                  info = sprintf("bounded m=%d eta=%g", m, eta_xi))
      if (abs(eta_xi) <= 5) {
        got_xm <- colSums(M_xm)
        want_xm <- drop(r2_richardson_fd12(function(p) fns_lg$d0(p, ld_xm), b_xm, step = 1e-3))
        expect_equal(unname((got_xm - want_xm) / (1 + abs(want_xm))),
                     rep(0, length(got_xm)), tolerance = 2e-5,
                     info = sprintf("bounded FD m=%d eta=%g", m, eta_xi))
      }
    }
  }
})

test_that("R2-3b: ZIDEGPD (m = 1..6) is exact for xi -> 0+, large ldelta, and deep-tail y = 0", {
  y <- c(0, 0, 1, 4, 15)
  for (m in 1:6) {
    # Note zidegpd4 parameter order: (lsigma, lxi, lkappa, ldelta, logitpi)!
    extra_zi <- switch(m,
      log(1.6),
      c(log(1.3), log(0.9), 0.25),
      log(1.4),
      c(log(1.9), log(1.4)),
      log(1.7),
      log(1.7))
    degpd_extra <- switch(m,
      log(1.6),
      c(log(1.3), log(0.9), 0.25),
      log(1.4),
      c(log(1.4), log(1.9)),
      log(1.7),
      log(1.7))

    npar <- length(extra_zi) + 3L
    fns_zi <- get(paste0(".zidegpd", m, "fns"), asNamespace("egpd"))
    fns_id <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    ld_zi <- r2_likdata(y, npar)

    for (lxi in c(log(0.25), log(1e-5), log(1e-12), log(1e-16), -800)) {
      logitpi <- -0.7
      b_zi <- c(log(3.5), lxi, extra_zi, logitpi)
      b_id <- c(log(3.5), exp(lxi), degpd_extra)

      # Independent log-space reference from verified DEGPD log-PMFs
      l0 <- vapply(y, function(yy) -fns_id$d0(b_id, r2_likdata(yy, length(b_id))), 0)
      lpi <- -log1p(exp(-logitpi))
      lompi <- -logitpi + lpi
      want_nll <- sum(ifelse(y > 0, -(lompi + l0),
                             -{
                               a <- lpi; b <- lompi + l0; mx <- pmax(a, b)
                               mx + log(exp(a - mx) + exp(b - mx))
                             }))
      expect_equal(fns_zi$d0(b_zi, ld_zi), want_nll, tolerance = 1e-12,
                   info = sprintf("zidegpd%d d0 lxi=%g", m, lxi))

      M_zi <- fns_zi$d120(b_zi, ld_zi)
      expect_true(all(is.finite(M_zi)),
                  info = sprintf("zidegpd%d d12 finite lxi=%g", m, lxi))

      if (lxi >= log(1e-12)) {
        got_zi <- colSums(M_zi)
        want_zi <- drop(r2_richardson_fd12(function(p) fns_zi$d0(p, ld_zi), b_zi, step = 1e-3))
        expect_equal(unname((got_zi - want_zi) / (1 + abs(want_zi))),
                     rep(0, length(got_zi)), tolerance = 2e-5,
                     info = sprintf("zidegpd%d FD lxi=%g", m, lxi))
      }
    }
  }

  # Reproducer: zidegpd3 and zidegpd4 with large ldelta = 5 (delta ~ 148) on counts up to 30
  y_hi <- c(0, 1, 5, 15, 30)
  for (ldelta in c(4.68, 5.0, 12.0)) {
    b4 <- c(1.77, -2.23, 1.10, ldelta, -1.39)
    ld4 <- r2_likdata(y_hi, 5L)
    d0_4 <- egpd:::.zidegpd4fns$d0(b4, ld4)
    M4 <- egpd:::.zidegpd4fns$d120(b4, ld4)
    expect_true(is.finite(d0_4))
    expect_true(all(is.finite(M4)))
    want4 <- drop(r2_richardson_fd12(function(p) egpd:::.zidegpd4fns$d0(p, ld4), b4, step = 1e-3))
    expect_equal(unname((colSums(M4) - want4) / (1 + abs(want4))),
                 rep(0, 20), tolerance = 2e-5,
                 info = paste("zidegpd4 ldelta=", ldelta))

    b3 <- c(1.77, -2.23, ldelta, -1.39)
    ld3 <- r2_likdata(y_hi, 4L)
    M3 <- egpd:::.zidegpd3fns$d120(b3, ld3)
    expect_true(all(is.finite(M3)))
    want3 <- drop(r2_richardson_fd12(function(p) egpd:::.zidegpd3fns$d0(p, ld3), b3, step = 1e-3))
    expect_equal(unname((colSums(M3) - want3) / (1 + abs(want3))),
                 rep(0, 14), tolerance = 2e-5,
                 info = paste("zidegpd3 ldelta=", ldelta))
  }

  # Extreme y = 0 log-sum-exp corner: pi -> 0 (logitpi = -1000) AND l0 = -922.1327 < -800
  ld0_zi1 <- r2_likdata(0, 4L)
  b_deep <- c(log(100), log(0.1), log(200), -1000)
  true_nll1 <- -200 * log(-expm1(-log1p(0.1 / 100) / 0.1))
  expect_equal(egpd:::.zidegpd1fns$d0(b_deep, ld0_zi1), true_nll1, tolerance = 1e-13)
  M_deep <- egpd:::.zidegpd1fns$d120(b_deep, ld0_zi1)
  expect_true(all(is.finite(M_deep)))
  M_degpd1 <- egpd:::.degpd1fns$d120(b_deep[1:3], r2_likdata(0, 3L))
  expect_equal(M_deep[1, 1:3], M_degpd1[1, 1:3], tolerance = 1e-12)
  expect_equal(M_deep[1, 4], 0, tolerance = 1e-15)
})

test_that("R2-4: ginv.egpd and pinv compute Moore-Penrose pseudoinverses via std SVD", {
  check_mp <- function(A, G, tol = 1e-12) {
    expect_equal(dim(G), rev(dim(A)))
    expect_equal(A %*% G %*% A, A, tolerance = tol)
    expect_equal(G %*% A %*% G, G, tolerance = tol)
    expect_equal(t(A %*% G), A %*% G, tolerance = tol)
    expect_equal(t(G %*% A), G %*% A, tolerance = tol)
  }

  A2 <- matrix(c(2, 1, 1, 2), 2, 2)
  expect_equal(ginv.egpd(A2), matrix(c(2/3, -1/3, -1/3, 2/3), 2, 2), tolerance = 1e-14)
  expect_equal(egpd:::pinv(A2), matrix(c(2/3, -1/3, -1/3, 2/3), 2, 2), tolerance = 1e-14)

  mats <- list(
    matrix(c(1, 2, 3, 2, 4, 6), 3, 2),                    # tall rank-1
    matrix(c(1, 0, 2, -1, 3, 1, 5, 0), 2, 4),             # wide rank-2
    tcrossprod(cbind(c(1, -2, 3), c(0, 1, -1))),          # 3x3 symmetric rank-2
    matrix(0, 3, 2)                                       # all-zero rectangular
  )
  for (M in mats) {
    check_mp(M, ginv.egpd(M))
    check_mp(M, egpd:::pinv(M))
  }
})
