test_that("CMP inversion attains high-dispersion-parameter requested means", {
  for (p in list(c(5, 100), c(1000, 10), c(10, 1000), c(5.3, 100))) {
    x <- 0:ceiling(2 * p[1] + 30)
    mass <- dcmp(x, mu = p[1], nu = p[2])
    expect_equal(sum(mass), 1, tolerance = 1e-8)
    expect_equal(sum(x * mass), p[1], tolerance = 1e-8)
  }
  expect_error(dcmp(0, mu = 1e12, nu = 1), "numerical limits")
  expect_equal(dcmp(1, mu = 1e-100, nu = 1, log = TRUE),
               dpois(1, 1e-100, log = TRUE), tolerance = 1e-10)
})

test_that("public GPIG survives the seed underflow using exact probabilities", {
  x <- c(0, 1, 999, 1000, 1001)
  expect_equal(dgpig(x, 1, 1000, 1, log = TRUE),
               dpois(x, 1000, log = TRUE), tolerance = 1e-10)
  mu <- 2000; sigma <- .001
  # Independent Poisson/inverse-Gaussian mixing integral in log mixing scale.
  ref <- integrate(function(t) {
    z <- exp(t)
    exp(dpois(2000, mu * z, log = TRUE) -
          .5 * log(2 * pi * sigma) - .5 * t - (z - 1)^2 / (2 * sigma * z))
  }, -.5, .5, rel.tol = 1e-10)$value
  b <- sqrt(1 + 2 * sigma * mu) / sigma
  cc <- 2 * sigma * mu / (1 + 2 * sigma * mu)
  expect_equal(dgpig(2000, .5, b, cc), ref, tolerance = 1e-9)
  expect_equal(dgpig(0:1, .5, b, cc, log = TRUE),
               c(b * (sqrt(1 - cc) - 1),
                 b * (sqrt(1 - cc) - 1) + log(.5 * b * cc)), tolerance = 1e-10)
  # Far beyond the mode a common rescaling cannot hold the full ladder.
  mu <- 20 * .6 / (2 * sqrt(.4)); sigma <- 1 / (20 * sqrt(.4))
  X <- matrix(1, 1, 1)
  tail <- -egpd:::pig1d0(list(log(mu), log(sigma)), X, X, 2000, integer(), 0L,
                        rep(list(numeric()), 2))
  expect_equal(dgpig(2000, .5, 20, .6, log = TRUE), tail, tolerance = 1e-10)
  expect_true(is.finite(dgpig(2, .5, 1, 1e-300, log = TRUE)))
})

test_that("PIG and ZIPIG retain their Poisson-limit likelihood and derivatives", {
  X <- matrix(1, 1, 1)
  for (sigma in c(1e-8, 1e-12, 1e-20)) for (y in c(0, 2, 5)) {
    pars <- list(log(2), log(sigma))
    lp <- -egpd:::pig1d0(pars, X, X, y, integer(), 0L, rep(list(numeric()), 2))
    expect_equal(lp, dpois(y, 2, log = TRUE), tolerance = 2e-7)
    expect_equal(dpig(y, 2, sigma, log = TRUE), lp, tolerance = 1e-12)
    expect_equal(ppig(y, 2, sigma), ppois(y, 2), tolerance = 2e-7)
    deriv <- egpd:::pig1d12(pars, X, X, y, integer(), 0L, rep(list(numeric()), 2))
    expect_true(all(is.finite(deriv)))
    expect_equal(unname(deriv[1, c(1, 2, 3)]), c(2 - y, 0, 2), tolerance = 2e-7)
    zi <- c(pars, list(qlogis(.2)))
    zlp <- -egpd:::zipig1d0(zi, X, X, X, y, integer(), 0L, rep(list(numeric()), 3))
    ref <- log(.8 * dpois(y, 2) + .2 * (y == 0))
    expect_equal(zlp, ref, tolerance = 2e-7)
    expect_equal(dzipig(y, 2, sigma, .2, log = TRUE), ref, tolerance = 2e-7)
  }
  p <- c(.1, .5, .9)
  expect_equal(qpig(p, 2, 1e-20), qpois(p, 2))
  expect_equal(qzipig(p, 2, 1e-20, .2),
               qpois(pmax(0, (p - .2) / .8), 2))
  ## draws in the Poisson limit follow Poisson(2) (rpig samples the Poisson-inverse
  ## Gaussian mixture exactly, so it is no longer qpois() of the same uniforms)
  gof <- function(x, pr) {
    obs <- tabulate(pmin(x, length(pr)) + 1, length(pr) + 1)
    e <- length(x) * c(pr, 1 - sum(pr))
    pchisq(sum((obs - e)^2 / e), length(pr), lower.tail = FALSE)
  }
  set.seed(841)
  expect_gt(gof(rpig(2e4, 2, 1e-20), dpois(0:7, 2)), 1e-3)
  set.seed(842)
  expect_gt(gof(rzipig(2e4, 2, 1e-20, .2), .2 * (0:7 == 0) + .8 * dpois(0:7, 2)), 1e-3)
})

test_that("normalization sums the selected GPIG lattice including zero", {
  old <- gpig_control()
  on.exit(do.call(gpig_control, old))
  for (method in c("saddlepoint", "hybrid")) {
    gpig_control(method, normalize = FALSE, yswitch = 5, ymax = 100, nquad = 32)
    raw <- vapply(0:200, function(y) egpd:::gpig_logp_cpp(y, .5, 1, .5)[["logp"]], 0)
    expected <- log(sum(exp(raw)))
    expect_equal(egpd:::gpig_logS_cpp(.5, 1, .5), expected, tolerance = 1e-9)
    gpig_control(method, normalize = TRUE, yswitch = 5, ymax = 100, nquad = 32)
    normalized <- vapply(0:200, function(y) egpd:::gpig_logp_cpp(y, .5, 1, .5)[["logp"]], 0)
    expect_equal(normalized, raw - expected, tolerance = 1e-8)
    expect_equal(sum(exp(normalized)), 1, tolerance = 1e-9)
    h <- 1e-5
    lp <- function(b) egpd:::gpig_logp_cpp(0L, .5, b, .5)
    expect_equal(lp(1)[["db"]], (lp(1 + h)[["logp"]] - lp(1 - h)[["logp"]]) / (2 * h),
                 tolerance = 1e-6)
  }
  gpig_control("saddlepoint", normalize = TRUE)
  expect_error(egpd:::gpig_logS_cpp(.5, 1e6, .6), "did not converge")
})

test_that("scores respect active clamps rather than differentiating through them", {
  h <- 1e-5
  prl <- function(e) egpd:::prl_grad_cpp(0, e)
  expect_equal(prl(log(2))[["g"]], 0)
  expect_equal((prl(log(2) + h)[["logp"]] - prl(log(2) - h)[["logp"]]) / (2 * h), 0)
  for (case in list(list("gwaring", c(log(10), log(2), log(.5)), 3),
                    list("gwaring", c(log(1e-20), log(2), log(2)), 1),
                    list("gpois", c(log(1e20), 0), 1),
                    list("gpois", c(0, -100), 2),
                    list("plnorm", c(0, log(1e-10)), 2))) {
    fn <- function(e) egpd:::cf_grad_cpp(case[[1]], 0, e)
    e <- case[[2]]; j <- case[[3]]
    plus <- minus <- e; plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
    expect_equal(unname(fn(e)[j + 1]), 0)
    expect_equal(unname((fn(plus)[1] - fn(minus)[1]) / (2 * h)), 0)
  }
})

test_that("quantile caps never fabricate values and endpoints are mathematical", {
  skip_on_cran()
  for (fn in list(qbell, qgpig, qzibell, qzigpig)) {
    expect_equal(fn(c(0, 1)), c(0, Inf))
    expect_warning(expect_true(all(is.na(fn(c(NA, NaN, -.1, 1.1))))), "NaNs produced")
    expect_error(fn(.5, max.value = NA), "max.value")
    expect_error(fn(.5, max.value = 2.5), "max.value")
  }
  expect_warning(expect_true(is.na(qbell(.99, theta = 2, max.value = 2))), "not attained")
  expect_warning(expect_true(is.na(qgpig(.9, a = .2, b = 1, c = 1, max.value = 100))), "not attained")
  expect_equal(qgpig(c(0, .5, 1), c = 0), c(0, 0, 0))
  for (fn in list(qzibell, qzigpig))
    expect_warning(expect_true(is.na(fn(-.1, pi = .9))), "NaNs produced")
  expect_warning(expect_true(is.na(qbell(.5, theta = -1))), "NaNs produced")
  expect_warning(expect_true(is.na(qgpig(.5, a = 2))), "NaNs produced")
})

test_that("CompPareto large alpha agrees with conditional Lomax oracles", {
  for (alpha in c(54, 60, 1000, 10000)) {
    phi <- 4 * dexp(2) / (alpha * pexp(2))
    expect_equal(dcomppareto(1, "exp", alpha, 2, rate = 1),
                 dexp(1) / ((1 + phi) * pexp(2)), tolerance = 1e-12)
    x <- 2 + 1 / alpha
    expect_equal(dcomppareto(x, "exp", alpha, 2, rate = 1),
                 phi / (1 + phi) * alpha / (x + 2) * (4 / (x + 2))^alpha,
                 tolerance = 1e-10)
    p <- c(.25, .75, .99999)
    q <- qcomppareto(p, "exp", alpha, 2, rate = 1)
    expect_equal(pcomppareto(q, "exp", alpha, 2, rate = 1), p, tolerance = 1e-10)
    expect_equal(pcomppareto(Inf, "exp", alpha, 2, rate = 1), 1)
  }
})

test_that("Bell log probabilities do not exponentiate before logging", {
  y <- 500; theta <- .1
  k <- 1:1000
  terms <- y * log(k) - lgamma(k + 1)
  logbell <- max(terms) + log(sum(exp(terms - max(terms)))) - 1
  ref <- y * log(theta) - expm1(theta) + logbell - lgamma(y + 1)
  expect_equal(dbell(y, theta, log = TRUE), ref, tolerance = 1e-10)
  expect_equal(dzibell(y, theta, .2, log = TRUE), ref + log(.8), tolerance = 1e-10)
  expect_equal(dbell(y, theta), 0)
  expect_equal(dzibell(0, 10, pi = 0, log = TRUE), dbell(0, 10, log = TRUE))
  expect_equal(dzigpig(0, 1, 1000, 1, pi = 0, log = TRUE), -1000)
  expect_equal(dzipig(0, 1000, 1e-20, pi = 0, log = TRUE), -1000)
})

test_that("count CDF infinity and missing endpoints propagate in both tails", {
  for (fn in list(pbell, pgpig, pcmp, pgpois, pgwaring, ppoislnorm, pprl, ppig, pzipig,
                  pzibell, pzigpig)) {
    expect_equal(fn(c(-Inf, Inf, NA_real_)), c(0, 1, NA_real_))
    expect_equal(fn(c(-Inf, Inf), lower.tail = FALSE, log.p = TRUE), c(0, -Inf))
    expect_equal(fn(numeric()), numeric())
  }
})
