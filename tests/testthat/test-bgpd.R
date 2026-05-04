numeric_derivative_bgpd <- function(cdf_fun, x) {
  h <- 1e-6 * pmax(1, x)
  (cdf_fun(x + h) - cdf_fun(pmax(x - h, 0))) / (2 * h)
}

test_that("blended GP matches the bulk below the blending interval and the GP in the tail", {
  args <- list(
    xi = 0.2,
    pbulk = pbeta,
    dbulk = dbeta,
    qbulk = qbeta,
    bulk.args = list(shape1 = 2.5, shape2 = 6),
    pa = 0.7,
    pb = 0.95,
    c1 = 5,
    c2 = 6
  )

  a <- do.call(qbeta, c(list(p = args$pa), args$bulk.args))
  b <- do.call(qbeta, c(list(p = args$pb), args$bulk.args))
  A <- (1 - args$pa)^(-args$xi)
  B <- (1 - args$pb)^(-args$xi)
  sigma <- args$xi * (a - b) / (A - B)
  u <- a - sigma * (A - 1) / args$xi

  x_bulk <- a / 2
  x_tail <- b + 0.5

  expect_equal(do.call(pbgpd, c(list(q = x_bulk), args)),
               do.call(pbeta, c(list(q = x_bulk), args$bulk.args)))
  expect_equal(do.call(dbgpd, c(list(x = x_bulk), args)),
               do.call(dbeta, c(list(x = x_bulk), args$bulk.args)))

  expect_equal(do.call(pbgpd, c(list(q = x_tail), args)),
               1 - (1 + args$xi * (x_tail - u) / sigma)^(-1 / args$xi))
  expect_equal(do.call(dbgpd, c(list(x = x_tail), args)),
               (1 / sigma) * (1 + args$xi * (x_tail - u) / sigma)^(-1 / args$xi - 1))
})

test_that("blended GP p/q/d functions roundtrip for positive and negative shape", {
  settings <- list(
    list(xi = 0.2, pa = 0.75, pb = 0.96, c1 = 5, c2 = 5),
    list(xi = -0.15, pa = 0.7, pb = 0.92, c1 = 6, c2 = 5)
  )
  probs <- c(0.05, 0.2, 0.5, 0.8, 0.9, 0.97)

  for (cfg in settings) {
    args <- c(
      cfg,
      list(
        pbulk = pbeta,
        dbulk = dbeta,
        qbulk = qbeta,
        bulk.args = list(shape1 = 3, shape2 = 7)
      )
    )
    q <- do.call(qbgpd, c(list(p = probs), args))
    expect_equal(do.call(pbgpd, c(list(q = q), args)), probs, tolerance = 1e-7)
  }
})

test_that("blended GP density matches a numerical derivative of the CDF", {
  args <- list(
    xi = 0.2,
    pbulk = pbeta,
    dbulk = dbeta,
    qbulk = qbeta,
    bulk.args = list(shape1 = 2.5, shape2 = 6),
    pa = 0.7,
    pb = 0.95,
    c1 = 5,
    c2 = 6
  )

  a <- do.call(qbeta, c(list(p = args$pa), args$bulk.args))
  b <- do.call(qbeta, c(list(p = args$pb), args$bulk.args))
  x <- c(a / 2, (a + b) / 2, b + 0.35)

  d_num <- numeric_derivative_bgpd(function(v) do.call(pbgpd, c(list(q = v), args)), x)
  d_val <- do.call(dbgpd, c(list(x = x), args))

  expect_true(max(abs(d_num - d_val) / pmax(abs(d_val), 1e-8)) < 1e-5)
})

test_that("blended GP random generation respects the finite upper endpoint when xi is negative", {
  args <- list(
    xi = -0.2,
    pbulk = pbeta,
    dbulk = dbeta,
    qbulk = qbeta,
    bulk.args = list(shape1 = 3, shape2 = 6),
    pa = 0.75,
    pb = 0.95,
    c1 = 5,
    c2 = 5
  )

  endpoint <- do.call(qbgpd, c(list(p = 1), args))
  samp <- do.call(rbgpd, c(list(n = 200, unifsamp = seq(0.01, 0.99, length.out = 200)), args))

  expect_true(all(samp >= 0))
  expect_true(all(samp <= endpoint + 1e-10))
  expect_equal(do.call(pbgpd, c(list(q = endpoint), args)), 1)
})
