type_args <- function(type) {
  switch(as.character(type),
    "1" = list(kappa = 1.7),
    "2" = list(kappa = 3.2),
    "3" = list(kappa = 2.3),
    "4" = list(delta = 1.4),
    "5" = list(kappa = 1.6, delta = 1.3),
    "6" = list(prob = 0.35, kappa = 1.4, delta = 2.2)
  )
}

numeric_derivative <- function(cdf_fun, x) {
  h <- 1e-6 * pmax(1, x)
  (cdf_fun(x + h) - cdf_fun(pmax(x - h, 0))) / (2 * h)
}

previous_cdf <- function(q, cdf_fun, step = 1) {
  out <- numeric(length(q))
  idx <- q > 0
  if (any(idx)) {
    out[idx] <- cdf_fun(q[idx] - step)
  }
  out
}

test_that("G transformation p/q/d functions are internally consistent", {
  probs <- c(0.15, 0.5, 0.85)
  u <- c(0.2, 0.5, 0.8)

  for (tp in 1:6) {
    args <- c(list(type = tp), type_args(tp))
    tol <- if (tp == 6) 2e-5 else 1e-10

    qu <- do.call(q.G, c(list(u = probs), args))
    expect_equal(do.call(p.G, c(list(u = qu), args)), probs,
                 tolerance = tol, info = paste("type", tp))

    d_num <- numeric_derivative(function(x) do.call(p.G, c(list(u = x), args)), u)
    d_val <- do.call(d.G, c(list(u = u), args))
    expect_true(max(abs(d_num - d_val) / pmax(abs(d_val), 1e-8)) < 1e-5,
                info = paste("type", tp))
  }
})

test_that("EGPD p/q/d functions roundtrip across all supported types", {
  probs <- c(0.15, 0.5, 0.85)

  for (tp in 1:6) {
    args <- c(list(type = tp, sigma = 2, xi = 0.2), type_args(tp))
    tol <- if (tp == 6) 2e-5 else 1e-10

    q <- do.call(qegpd, c(list(p = probs), args))
    expect_equal(do.call(pegpd, c(list(q = q), args)), probs,
                 tolerance = tol, info = paste("type", tp))

    d_num <- numeric_derivative(function(x) do.call(pegpd, c(list(q = x), args)), q)
    d_val <- do.call(degpd_density, c(list(x = q), args))
    expect_true(max(abs(d_num - d_val) / pmax(abs(d_val), 1e-8)) < 1e-5,
                info = paste("type", tp))
  }
})

test_that("Discrete EGPD p/q/d functions roundtrip across all supported types", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  x <- 0:12

  for (tp in 1:6) {
    args <- c(list(type = tp, sigma = 3, xi = 0.2), type_args(tp))
    pmf <- do.call(ddiscegpd, c(list(x = x), args))
    cdf <- do.call(pdiscegpd, c(list(q = x), args))
    cdf_prev <- c(0, cdf[-length(cdf)])

    expect_equal(pmf, cdf - cdf_prev, tolerance = 1e-12,
                 info = paste("type", tp))

    q <- do.call(qdiscegpd, c(list(p = probs), args))
    cdf_at_q <- do.call(pdiscegpd, c(list(q = q), args))
    cdf_before <- previous_cdf(q, function(v) do.call(pdiscegpd, c(list(q = v), args)))

    expect_true(all(cdf_at_q >= probs - 1e-12), info = paste("type", tp))
    expect_true(all(cdf_before < probs + 1e-12), info = paste("type", tp))
  }
})

test_that("Zero-inflated discrete EGPD p/q/d functions roundtrip across all supported types", {
  probs <- c(0.05, 0.2, 0.35, 0.7, 0.9)
  x <- 0:12

  for (tp in 1:6) {
    args <- c(list(type = tp, pi = 0.25, sigma = 3, xi = 0.2), type_args(tp))
    pmf <- do.call(dzidiscegpd, c(list(x = x), args))
    cdf <- do.call(pzidiscegpd, c(list(q = x), args))
    cdf_prev <- c(0, cdf[-length(cdf)])

    expect_equal(pmf, cdf - cdf_prev, tolerance = 1e-12,
                 info = paste("type", tp))

    q <- do.call(qzidiscegpd, c(list(p = probs), args))
    cdf_at_q <- do.call(pzidiscegpd, c(list(q = q), args))
    cdf_before <- previous_cdf(q, function(v) do.call(pzidiscegpd, c(list(q = v), args)))

    expect_true(all(cdf_at_q >= probs - 1e-12), info = paste("type", tp))
    expect_true(all(cdf_before < probs + 1e-12), info = paste("type", tp))
  }
})

test_that("Zero-inflated continuous EGPD p/q/d functions roundtrip across all supported types", {
  probs <- c(0.05, 0.2, 0.35, 0.7, 0.9)
  pi0 <- 0.25

  for (tp in 1:6) {
    args <- c(list(type = tp, pi = pi0, sigma = 2, xi = 0.2), type_args(tp))
    tol <- if (tp == 6) 2e-5 else 1e-6

    q <- do.call(qziegpd, c(list(p = probs), args))
    cdf_at_q <- do.call(pziegpd, c(list(q = q), args))

    expect_equal(q[probs <= pi0], rep(0, sum(probs <= pi0)),
                 info = paste("type", tp))
    expect_equal(cdf_at_q[probs <= pi0], rep(pi0, sum(probs <= pi0)),
                 tolerance = 1e-12, info = paste("type", tp))
    expect_equal(cdf_at_q[probs > pi0], probs[probs > pi0],
                 tolerance = tol, info = paste("type", tp))

    x <- do.call(qziegpd, c(list(p = c(pi0 + 0.1, 0.75, 0.9)), args))
    d_num <- numeric_derivative(function(v) do.call(pziegpd, c(list(q = v), args)), x)
    d_val <- do.call(dziegpd, c(list(x = x), args))
    expect_true(max(abs(d_num - d_val) / pmax(abs(d_val), 1e-8)) < 1e-5,
                info = paste("type", tp))
  }
})

test_that("Compound Poisson EGPD p/q/d functions roundtrip across all supported types", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  h <- 0.25
  K <- 250
  x <- seq(0, 10, by = h)

  for (tp in 1:6) {
    args <- c(list(type = tp, lambda = 1.2, sigma = 1.5, xi = 0.15, h = h, K = K), type_args(tp))
    pmf <- do.call(dcpegpd, c(list(x = x), args))
    cdf <- do.call(pcpegpd, c(list(q = x), args))
    cdf_prev <- c(0, cdf[-length(cdf)])

    expect_equal(pmf, cdf - cdf_prev, tolerance = 1e-12,
                 info = paste("type", tp))

    q <- do.call(qcpegpd, c(list(p = probs), args))
    cdf_at_q <- do.call(pcpegpd, c(list(q = q), args))
    cdf_before <- previous_cdf(q, function(v) do.call(pcpegpd, c(list(q = v), args)), step = h)

    expect_true(all(cdf_at_q >= probs - 1e-12), info = paste("type", tp))
    expect_true(all(cdf_before < probs + 1e-12), info = paste("type", tp))
  }
})

test_that("Compound Poisson discrete EGPD p/q/d functions roundtrip across all supported types", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  x <- 0:25

  for (tp in 1:6) {
    args <- c(list(type = tp, lambda = 1.2, sigma = 2, xi = 0.15, K = 250), type_args(tp))
    pmf <- do.call(dcpdegpd, c(list(x = x), args))
    cdf <- do.call(pcpdegpd, c(list(q = x), args))
    cdf_prev <- c(0, cdf[-length(cdf)])

    expect_equal(pmf, cdf - cdf_prev, tolerance = 1e-12,
                 info = paste("type", tp))

    q <- do.call(qcpdegpd, c(list(p = probs), args))
    cdf_at_q <- do.call(pcpdegpd, c(list(q = q), args))
    cdf_before <- previous_cdf(q, function(v) do.call(pcpdegpd, c(list(q = v), args)))

    expect_true(all(cdf_at_q >= probs - 1e-12), info = paste("type", tp))
    expect_true(all(cdf_before < probs + 1e-12), info = paste("type", tp))
  }
})
