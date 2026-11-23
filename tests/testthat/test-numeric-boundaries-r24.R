test_that("zero inflation distinguishes an exact tie from a one-ulp interior fit", {
  y <- rep(0:3, each = 3)
  fit <- function(s, ...) distfit(y, family = "zidegpd",
    fix.arg = list(sigma = s, xi = -1, kappa = 1), ...)
  expect_warning(f <- fit(4), "natural boundary")
  expect_identical(unname(coef(f)), 0)
  expect_true(all(is.na(c(f$sd, confint(f)))))
  expect_error(vcov(f), "natural boundary")
  expect_warning(box <- fit(4, lower = qlogis(.1)), "optimizer bound")
  expect_equal(coef(box), c(pi = .1))

  s <- 4 + 2^-50
  interior <- fit(s)
  # For a floor-uniform base, P(0) = 1/s exactly.
  want <- (s - 4) / (4 * (s - 1))
  expect_equal(unname(coef(interior)) / want, 1, tolerance = 1e-13)
  expect_true(is.finite(interior$sd))
  # The squared uniform CDF has an exact mass of 1/16 at count zero.
  y2 <- c(0, rep(1:3, 5))
  expect_warning(square <- distfit(y2, family = "zidegpd",
    fix.arg = list(sigma = 4, xi = -1, kappa = 2)), "natural boundary")
  expect_identical(unname(coef(square)), 0)
})

test_that("the blended quantile inverts on the scale of its matching knots", {
  p <- c(.85, .9, .95)
  args <- list(xi = 0, bulk.args = list(shape1 = 1, shape2 = 1e14))
  q <- do.call(qbgpd, c(list(p), args))
  expect_equal(do.call(pbgpd, c(list(q), args)), p, tolerance = 1e-12)
})

test_that("a steep but finite matched GP tail is not an infinite quantile", {
  p <- c(.75, .75025, .751)
  for (xi in c(1023, 1025, 1100)) {
    args <- list(xi = xi, pa = .5, pb = .75, pbulk = punif, dbulk = dunif,
                 qbulk = qunif, bulk.args = list())
    # At these shapes exp(-xi*log(2)) is negligible. Use the ratio of
    # increments, not two separately enormous powers.
    ref <- .5 + .25 * exp(xi * log1p((p - .75) / (1 - p)))
    got <- do.call(qbgpd, c(list(p), args))
    expect_equal(got / ref, rep(1, length(p)), tolerance = 5e-13)
    expect_equal(do.call(pbgpd, c(list(ref), args)), p, tolerance = 1e-13)
    density <- (1 - p) / (xi * (ref - .5))
    expect_equal(do.call(dbgpd, c(list(ref), args)) / density, rep(1, length(p)),
                 tolerance = 1e-12)
    expect_equal(do.call(rbgpd, c(list(n = length(p), unifsamp = p), args)), got)
  }
})

test_that("inverse composite-Pareto tails preserve finite scaled answers", {
  th <- c(1e308, 1e-305)
  alpha <- c(2, .001)
  p <- c(.75, .51)
  # Conditional Lomax inverse using 256-bit arithmetic at these inputs.
  ref <- c(1.426122638850535e308, 53949.81981983951)
  q <- function(p, ...) qcomppareto(p, "exp", theta = th, alpha = alpha, rate = 1/th, ...)
  for (got in list(q(p), q(1 - p, lower.tail = FALSE),
                   q(log1p(-p), lower.tail = FALSE, log.p = TRUE)))
    expect_equal(got / ref, c(1, 1), tolerance = 3e-12)
  expect_equal(pcomppareto(ref, "exp", theta = th, alpha = alpha, rate = 1/th), p)
})

test_that("multivariate GPD draws keep their finite scaled magnitude and support", {
  tiny <- .Machine$double.xmin * .Machine$double.eps
  set.seed(24010)
  expect_true(all(rmdgpd(64, sigma = tiny, xi = 0, rate = 0) == 1))
  set.seed(241)
  expect_equal(as.vector(rmdgpd(1, sigma = 1e-307, xi = 3900, rate = 0)), c(5, 5))
  set.seed(241)
  radial <- exp(799 * -log1p(-runif(1)) + log(1e-307) - log(799))
  set.seed(241)
  z <- rbegpd(1, kappa = 1, sigma = 1e-307, xi = 799, thL = 3, thU = 3, thw = .2)
  expect_equal(sum(z), radial, tolerance = 5e-13)
})
