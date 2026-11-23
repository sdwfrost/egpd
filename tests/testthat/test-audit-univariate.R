audit_uni_theta <- function(model, xi = 0) {
  switch(model,
    c(log(2), xi, log(1.5)),
    c(log(2), xi, log(1.2), log(1.4), qlogis(.35)),
    c(log(2), xi, log(1.4)),
    c(log(2), xi, log(1.4), log(1.5)),
    c(log(2), xi, log(2)),
    c(log(2), xi, log(2)))
}

audit_uni_native <- function(model, order, theta, y, offsets = NULL, duplicate = FALSE) {
  n <- if (duplicate) 2L else length(y)
  np <- length(theta)
  ld <- list(X = rep(list(matrix(1, n, 1)), np), idpars = seq_len(np),
             y = matrix(y, ncol = 1),
             dupid = if (duplicate) c(0L, 1L, 0L) else 0L,
             duplicate = as.integer(duplicate), censored = FALSE,
             offsets = if (is.null(offsets)) rep(list(numeric()), np) else offsets)
  get(paste0(".egpd", model, ".d", order), asNamespace("egpd"))(theta, ld)
}

audit_uni_ref <- function(model, t, y) {
  z <- y / exp(t[1])
  A <- if (t[2] == 0) z else log1p(t[2] * z) / t[2]
  H <- -expm1(-A)
  k <- exp(t[3])
  carrier <- switch(model,
    k * H^(k - 1),
    plogis(t[5]) * k * H^(k - 1) +
      plogis(-t[5]) * (k + exp(t[4])) * H^(k + exp(t[4]) - 1),
    (1 + k) / k * (1 - exp(-k * A)),
    {
      G <- pbeta(exp(-k * A), 1 / k, 2, lower.tail = FALSE)
      exp(t[4]) / 2 * G^(exp(t[4]) / 2 - 1) *
        (1 + k) / k * (1 - exp(-k * A))
    },
    dnorm(H, 1, 1 / sqrt(k)) / (pnorm(sqrt(k)) - .5),
    (1/2 - 1/32) * dbeta(1/32 + (1/2 - 1/32) * H, k, k) /
      (pbeta(1/2, k, k) - pbeta(1/32, k, k)))
  sum(t[1] + (1 + t[2]) * A - log(carrier))
}

test_that("UNI-001: zero-shape helpers have exact exponential limits", {
  x <- c(0, 1e-10, 2, 10)
  expect_equal(pegpd(x, sigma = 2, xi = 0, kappa = 1), pexp(x, .5),
               tolerance = 1e-14)
  expect_equal(degpd_density(x, sigma = 2, xi = 0, kappa = 1), dexp(x, .5))
  p <- c(0, 1e-10, .5, 1 - 1e-12, 1)
  expect_equal(qegpd(p, sigma = 2, xi = 0, kappa = 1), qexp(p, .5))
  expect_equal(egpd:::.pgpd(2, 0, 2, 0, 0), pexp(2, .5))
  expect_equal(egpd:::.qgpd(.5, 0, 2, 0), qexp(.5, .5))
  expect_equal(egpd:::.qgpd(1e-20, 0, 2, 0), -2 * log1p(-1e-20),
               tolerance = 1e-14)
  L <- -log1p(-c(.2, .5))
  expect_equal(unname(egpd:::.dqgpd(c(.2, .5), log(2), 0)), cbind(2 * L, L^2))
  expect_equal(unname(egpd:::.dqgpd(1, log(2), 0)), matrix(Inf, 1, 2))
  for (xi in c(-1e-12, 1e-12)) {
    expect_equal(pegpd(x, sigma = 2, xi = xi, kappa = 1), pexp(x, .5),
                 tolerance = 1e-10)
    expect_equal(qegpd(.5, sigma = 2, xi = xi, kappa = 1), qexp(.5, .5),
                 tolerance = 1e-10)
  }
  p <- matrix(c(.2, .5, .7, .9), 2, dimnames = list(c("a", "b"), c("c", "d")))
  expect_equal(pegpd(qegpd(p, sigma = 2, xi = 0, kappa = 1),
                     sigma = 2, xi = 0, kappa = 1), p)
  expect_equal(egpd:::.qgpd(p, 0, 2, 0), qexp(p, .5))
})

test_that("UNI-001: native objectives and derivatives extend through zero shape", {
  y <- c(.5, 1, 3)
  for (m in 1:6) {
    for (xi in c(-.2, -1e-8, 0, 1e-8, .2)) {
      t <- audit_uni_theta(m, xi)
      expect_equal(audit_uni_native(m, "0", t, y), audit_uni_ref(m, t, y),
                   tolerance = 1e-10, info = paste(m, xi))
    }
    t <- audit_uni_theta(m)
    np <- length(t)
    h <- 1e-4
    f <- function(t) audit_uni_ref(m, t, y)
    gradient <- vapply(seq_len(np), function(i) {
      e <- replace(numeric(np), i, h)
      (f(t + e) - f(t - e)) / (2 * h)
    }, 0)
    hessian <- unlist(lapply(seq_len(np), function(i) vapply(i:np, function(j) {
      ei <- replace(numeric(np), i, h)
      ej <- replace(numeric(np), j, h)
      (f(t + ei + ej) - f(t + ei - ej) - f(t - ei + ej) + f(t - ei - ej)) / (4 * h^2)
    }, 0)))
    got <- colSums(audit_uni_native(m, "12", t, y))
    expect_equal(got, c(gradient, hessian), tolerance = 2e-5, info = paste("model", m))
    expect_gte(audit_uni_native(m, "0", t, -1), 1e20)
    expect_equal(audit_uni_native(m, "0", t, y, duplicate = TRUE),
                 audit_uni_native(m, "0", t, y))
    offsets <- rep(list(rep(0, length(y))), np)
    offsets[[2]] <- c(-1e-8, 0, 1e-8)
    expect_equal(audit_uni_native(m, "0", t, y, offsets),
                 sum(vapply(seq_along(y), function(i) {
                   ti <- t; ti[2] <- offsets[[2]][i]
                   audit_uni_ref(m, ti, y[i])
                 }, 0)), tolerance = 1e-10)
  }
  for (m in 1:4) {
    expect_true(all(is.finite(audit_uni_native(m, "34", audit_uni_theta(m), 1))))
  }
  higher <- drop(audit_uni_native(1, "34", c(log(2), 0, 0), 1))
  z <- .5
  expect_equal(higher[c(7, 21)], c(2 * z^3 - 1.5 * z^4, 24 * z^5/5 - 6 * z^4),
               tolerance = 5e-4)
  expect_equal(higher[c(10, 25)], rep(-log(-expm1(-z)), 2), tolerance = 5e-4)
  expect_equal(audit_uni_native(1, "0", c(log(2), 0, 0), 0), log(2))
})

test_that("UNI-002: densities mask support before combining carriers", {
  for (type in 1:6) {
    args <- list(prob = .35, kappa = .5, delta = .5, sigma = 2, xi = -.2, type = type)
    expect_equal(do.call(degpd_density, c(list(x = c(-1, 10, 11, Inf)), args)), rep(0, 4))
    expect_equal(do.call(degpd_density, c(list(x = c(-1, 10, 11, Inf), log = TRUE), args)),
                 rep(-Inf, 4))
    expect_identical(do.call(qziegpd, c(list(p = c(0, .1, .2), pi = .2), args)), rep(0, 3))
  }
  expect_equal(d.G(1, type = 4, delta = .5), 3)
  expect_equal(d.G(1, type = 5, kappa = 2, delta = .5), 3)
  expect_equal(d.G(0, type = 5, kappa = 1, delta = .5), sqrt(.75))
})

test_that("UNI-003: discrete distributions live on the nonnegative lattice", {
  x <- c(-Inf, -1, -.5, .5, 1.5, Inf)
  expect_equal(ddiscegpd(x, sigma = 2, xi = .2, kappa = 1), rep(0, length(x)))
  expect_equal(pdiscegpd(c(-Inf, -.5, 0, .5, 1, 1.5, Inf), sigma = 2, xi = .2, kappa = 1),
               c(0, 0, rep(1 - 1.1^-5, 2), rep(1 - 1.2^-5, 2), 1))
  expect_equal(sum(ddiscegpd(0:10, sigma = 2, xi = -.2, kappa = 1)), 1)
})

test_that("UNI-004: mixture inversion preserves extreme finite quantiles and recycles", {
  p <- c(.2, .8, 1 - 1e-8, 1 - 1e-12)
  s <- 2 * (1 - p) / (1.5 + sqrt(2.25 - 2 * (1 - p)))
  want <- 10 * expm1(-.2 * log(s))
  expect_equal(qegpd(p, type = 6, prob = .5, kappa = 1, delta = 2, sigma = 2, xi = .2),
               want, tolerance = 2e-5)
  got <- q.G(.7, type = 6, prob = c(.2, .6), kappa = c(1, 2), delta = c(2, 3))
  expect_equal(p.G(got, type = 6, prob = c(.2, .6), kappa = c(1, 2), delta = c(2, 3)),
               rep(.7, 2), tolerance = 1e-11)
  expect_equal(egpd:::.iG2_egpd(p, 1, 2, .5), q.G(p, type = 6, kappa = 1, delta = 2, prob = .5))
  expect_equal(egpd:::.iG2_zidegpd(p, 1, 2, .5), q.G(p, type = 6, kappa = 1, delta = 2, prob = .5))
})

test_that("UNI-005: zero-inflated quantiles respect atoms and endpoints", {
  args <- list(pi = .2, sigma = 2, xi = .2, kappa = 1)
  p0 <- do.call(pzidiscegpd, c(list(q = 0), args))
  atoms <- do.call(pzidiscegpd, c(list(q = 0:5), args))
  expect_equal(do.call(qzidiscegpd, c(list(p = atoms), args)), 0:5)
  p <- p0 + 1e-8
  q <- do.call(qzidiscegpd, c(list(p = p), args))
  expect_equal(q, 1)
  expect_gte(do.call(pzidiscegpd, c(list(q = q), args)), p)
  expect_lt(do.call(pzidiscegpd, c(list(q = q - 1), args)), p)
  for (fun in list(qziegpd, qzidiscegpd)) {
    expect_equal(do.call(fun, c(list(p = 1), args)), Inf)
    expect_equal(fun(c(0, .5, 1), pi = 1, sigma = 2, xi = .2, kappa = 1), rep(0, 3))
    expect_true(all(is.nan(fun(c(-.1, 1.1), pi = .2, sigma = 2, xi = .2, kappa = 1))))
  }
  p <- c(.2, .2 + 1e-8, .7, 1 - 1e-10)
  expect_equal(do.call(pziegpd, c(list(q = do.call(qziegpd, c(list(p = p), args))), args)), p,
               tolerance = 1e-13)
})

test_that("UNI-006: all zero-inflated functions recycle every parameter", {
  for (type in 1:6) {
    args <- list(pi = c(.1, .2, .3, .4), sigma = 2, xi = .2,
                 kappa = seq(1, 2, length.out = 4), delta = seq(1, 3, length.out = 4),
                 prob = seq(.2, .8, length.out = 4), type = type)
    for (fun in list(dzidiscegpd, pzidiscegpd, qzidiscegpd, dziegpd, pziegpd, qziegpd)) {
      at <- if (identical(fun, qziegpd) || identical(fun, qzidiscegpd)) .7 else 0
      vector <- do.call(fun, c(list(at), args))
      scalar <- vapply(1:4, function(i) {
        a <- lapply(args, function(x) rep_len(x, 4)[i])
        do.call(fun, c(list(at), a))
      }, 0)
      expect_equal(vector, scalar)
      args2 <- args; args2$pi <- .2
      expect_length(do.call(fun, c(list(at), args2)), 4)
    }
  }
  args <- list(pi = c(.1, .2), sigma = 2:4, xi = c(.1, .2, .3, .4),
               kappa = 1:5, delta = c(1, 2), prob = c(.2, .4, .6), type = 6)
  for (fun in list(dzidiscegpd, pzidiscegpd, qzidiscegpd, dziegpd, pziegpd, qziegpd)) {
    want <- vapply(1:5, function(i) {
      do.call(fun, c(list(.7), lapply(args, function(x) rep_len(x, 5)[i])))
    }, 0)
    expect_equal(do.call(fun, c(list(.7), args)), want)
  }
})

test_that("UNI-007: scalar GAMLSS parameters retain every zero-count score", {
  y <- c(0, 1, 0)
  h <- 1e-5
  want <- (dDEGPD1(y, mu = 2, sigma = .2, nu = 1.5 + h, log = TRUE) -
             dDEGPD1(y, mu = 2, sigma = .2, nu = 1.5 - h, log = TRUE)) / (2 * h)
  expect_equal(DEGPD1()$dldv(y, mu = 2, sigma = .2, nu = 1.5), want, tolerance = 1e-8)
})

test_that("UNI-008: no-tail splice closes its upper atom and handles tied knots", {
  ob <- list(tail_fit = NULL, bulk_probs = c(.2, .5, .8),
             threshold_index = 3, threshold_prob = .8)
  cdf <- function(values, at) drop(egpd:::.spliced_degpd_cdf_matrix(
    ob, data.frame(x = 0), at, bulk_mat = matrix(values, nrow = 1)))
  expect_equal(unname(cdf(c(1, 2, 3), 0:4)), c(0, .2, .5, 1, 1))
  expect_equal(unname(cdf(c(1, 3, 3), 0:4)), c(0, .2, .35, 1, 1))
  expect_equal(unname(cdf(c(2, 2, 3), 0:4)), c(0, .1, .5, 1, 1))
  expect_equal(unname(cdf(c(0, 0, 0), 0:2)), rep(1, 3))
  q <- egpd:::.spliced_degpd_bulk_quantile(c(.2, .5, .8), c(1, 3, 3),
                                          c(0, .1, .2, .35, .7, .9, 1))
  expect_equal(q, c(0, .5, 1, 2, 3, 3, 3))
  expect_equal(sum(diff(c(0, cdf(c(1, 3, 3), 0:5)))), 1)
})
