## Round-5 audit regressions (countr): count-engine log.p quantiles and zero-inflated
## log-probabilities, ZI-DEGPD log CDF near 0, CompPareto body upper tail, and the gamlss
## start cache and scale clamp. Expected values come from base R, closed forms, or 400-bit
## Rmpfr values hard-coded below.

test_that("VB5-05: lower-tail log.p count quantiles near 0 follow base R", {
  lp <- c(-1e-10, -1e-15, -1e-20, -1e-40)
  ## generalized Poisson with lambda = 0 is Poisson(mu)
  expect_identical(qgpois(lp, 3, 0, log.p = TRUE), qpois(lp, 3, log.p = TRUE))
  expect_identical(qpois(lp, 3, log.p = TRUE), c(19, 25, 30, 48))
  ## GPIG with a = 1 is Poisson(b c); its zero-inflated version moves 1 - F by 1 - pi
  expect_identical(qgpig(lp, a = 1, b = 5, c = 0.6, log.p = TRUE), c(19, 25, 30, 48))
  zref <- vapply(lp, function(l) {
    k <- 0
    while (log1p(-0.7 * ppois(k, 3, lower.tail = FALSE)) < l) k <- k + 1
    k
  }, numeric(1))
  expect_identical(qzigpig(lp, a = 1, b = 5, c = 0.6, pi = 0.3, log.p = TRUE), zref)
  ## PIG(3, 0.5): smallest q with log1p(-S(q)) >= lp, S summed from gamlss.dist::dPIG; the
  ## upper-tail form of the same levels has the same brute-force answer
  expect_identical(qpig(c(-1e-15, -1e-20), 3, 0.5, log.p = TRUE), c(107, 146))
  expect_identical(qpig(log(c(1e-15, 1e-20)), 3, 0.5, lower.tail = FALSE, log.p = TRUE),
                   c(107, 146))
})

test_that("VB5-05: upper-tail log.p counts near 0 (tiny lower tail) follow base R", {
  ## log P(X > k) for Poisson(50) at small k is -F(k), far below the old fixed fuzz
  ## (values compared as ratios: all.equal() tolerances are absolute below the tolerance)
  expect_equal(pgpois(c(0:3, 10), 50, 0, lower.tail = FALSE, log.p = TRUE) /
                 ppois(c(0:3, 10), 50, lower.tail = FALSE, log.p = TRUE), rep(1, 5),
               tolerance = 1e-12)
  ulp <- c(-1e-15, -1e-18, -1e-20, -1e-21)
  expect_identical(qgpois(ulp, 50, 0, lower.tail = FALSE, log.p = TRUE),
                   qpois(ulp, 50, lower.tail = FALSE, log.p = TRUE))
  expect_identical(qpois(ulp, 50, lower.tail = FALSE, log.p = TRUE), c(6, 3, 2, 1))
})

test_that("VB5-05: log.p round trips hold where log F or log S is near 0", {
  k <- 0:50
  expect_identical(qgpois(pgpois(k, 3, 0, log.p = TRUE), 3, 0, log.p = TRUE), as.numeric(k))
  k <- 0:120
  expect_identical(qpig(ppig(k, 3, 0.5, log.p = TRUE), 3, 0.5, log.p = TRUE), as.numeric(k))
  k <- 0:12
  expect_identical(qgpois(pgpois(k, 50, 0, lower.tail = FALSE, log.p = TRUE), 50, 0,
                          lower.tail = FALSE, log.p = TRUE), as.numeric(k))
  ## linear-scale quantiles keep base R's results
  expect_identical(qgpois(c(0.1, 0.5, 0.9, 1 - 1e-10), 3, 0), qpois(c(0.1, 0.5, 0.9, 1 - 1e-10), 3))
})

test_that("VB5-07: zero-inflated count log-probabilities are never positive", {
  ## c = 0 is a point mass at 0, so P(Y = 0) = 1 exactly
  pis <- seq(0.001, 0.999, by = 0.001)
  expect_identical(dzigpig(0, a = 0.5, b = 1, c = 0, pi = pis, log = TRUE), rep(0, length(pis)))
  ## PIG: log p0 = -2 mu / (1 + sqrt(1 + 2 sigma mu)) = -1e-20, so
  ## log P(0) = log1p(-(1 - pi)(1 - p0)) = log1p(-0.95e-20)
  ## (compared as ratios: all.equal() tolerances are absolute below the tolerance)
  expect_equal(dzipig(0, 1e-20, 1, 0.05, log = TRUE) / log1p(-0.95e-20), 1, tolerance = 1e-12)
  ## Bell: p0 = exp(1 - e^theta)
  expect_equal(dzibell(0, 1e-20, 0.05, log = TRUE) / log1p(-0.95 * -expm1(-expm1(1e-20))), 1,
               tolerance = 1e-12)
  expect_equal(dzibell(0, 1e-5, 0.3, log = TRUE) / log1p(-0.7 * -expm1(-expm1(1e-5))), 1,
               tolerance = 1e-14)
  ## a CDF whose base survivor is not resolved: log F_zi is in [-1e-18, 0], and the
  ## quantile round trip is defined
  lp <- suppressWarnings(pzigpig(100, a = 0.99, b = 1e-18, c = 1, pi = 0.05, log.p = TRUE))
  expect_lte(lp, 0)
  expect_gte(lp, -1e-18)
  expect_false(is.nan(expect_silent(qzigpig(lp, a = 0.99, b = 1e-18, c = 1, pi = 0.05,
                                            log.p = TRUE))))
})

test_that("VB5-08: ZI-DEGPD lower-tail log CDF near 0 keeps relative accuracy", {
  ## F(y) = pi + (1 - pi)[1 - (1 + xi (y + 1)/sigma)^(-1/xi)]^kappa; 400-bit Rmpfr values
  ref <- c(-1.0011661561915850e-15, -1.0006040822352972e-15,
           -1.0000423869074633e-15, -9.9948106991055304e-16)
  ## (compared as ratios: all.equal() tolerances are absolute below the tolerance)
  expect_equal(pzidiscegpd(8892:8895, pi = 0.2, kappa = 0.7, sigma = 2, xi = 0.2,
                           log.p = TRUE) / ref, rep(1, 4), tolerance = 1e-10)
  ## smallest y with log F(y) >= lp (400-bit brute force)
  expect_identical(qzidiscegpd(c(-1e-14, -1e-15, -1e-20), pi = 0.2, kappa = 0.7, sigma = 2,
                               xi = 0.2, log.p = TRUE), c(5608, 8895, 89040))
  expect_equal(pziegpd(c(1e3, 1e6), pi = 0.2, kappa = 0.7, sigma = 2, xi = 0.2, log.p = TRUE) /
                 c(-5.3282078508157913e-11, -5.5997200083998199e-26), c(1, 1), tolerance = 1e-10)
  ## pi = 0 keeps the base log CDF, kappa log1p(-(1 + xi (y + 1)/sigma)^(-1/xi))
  y <- 8892:8895
  expect_equal(pzidiscegpd(y, pi = 0, kappa = 0.7, sigma = 2, xi = 0.2, log.p = TRUE) /
                 (0.7 * log1p(-(1 + 0.2 * (y + 1) / 2)^(-5))), rep(1, 4), tolerance = 1e-12)
})

## CompPareto law (documented): for q <= theta,
##   S(q) = [S_G(q) - S_G(theta) + phi G(theta)] / [(1 + phi) G(theta)],
## phi = 2 theta g(theta) / (alpha G(theta)); references in 400-bit Rmpfr from the closed-form
## body survivors (gamma(2, 1): e^-x (1 + x); weibull(2, 1): e^-x^2; exp(1): e^-x;
## lnorm(0, 1): pnorm(log x, lower.tail = FALSE)), quantiles by bisection on them.
test_that("VB5-06: CompPareto body upper tail comes from the body survivor", {
  ## splice far in the gamma tail (g(theta) underflows): log S was -Inf
  expect_equal(pcomppareto(c(799, 800, 801), "gamma", alpha = 3, theta = 800, shape = 2,
                           rate = 1, lower.tail = FALSE, log.p = TRUE),
               c(-787.03302738937767, -787.03624165277231, -787.03811606707884), tolerance = 1e-13)
  expect_equal(qcomppareto(c(-700, -787.03, -787.0362), "gamma", alpha = 3, theta = 800, shape = 2,
                           rate = 1, lower.tail = FALSE, log.p = TRUE),
               c(706.561825010874, 798.531414874622, 799.978027590362), tolerance = 1e-12)
  ## an ordinary row with tail weight 3.8e-19: S was 118% off and its quantile non-monotone
  wS <- exp(c(-41.003603955774922, -42.418418651624963))
  expect_equal(pcomppareto(c(6.4, 6.5), "weibull", alpha = 200, theta = 6.5, shape = 2, scale = 1,
                           lower.tail = FALSE) / wS, c(1, 1), tolerance = 1e-12)
  ## ... and the lower tail on the log scale, log F = log(1 - S) = -S
  expect_equal(pcomppareto(c(6.4, 6.5), "weibull", alpha = 200, theta = 6.5, shape = 2, scale = 1,
                           log.p = TRUE) / -wS, c(1, 1), tolerance = 1e-12)
  expect_equal(qcomppareto(c(8e-19, 5e-19, 3.788e-19), "weibull", alpha = 200, theta = 6.5,
                           shape = 2, scale = 1, lower.tail = FALSE),
               c(6.44875786896218, 6.48148866043365, 6.49992588912777), tolerance = 1e-12)
  ## lower-tail body quantiles at log p near 0 (qgamma loses accuracy there)
  ref_d <- 506.228962504222
  expect_equal(qcomppareto(-500, "gamma", alpha = 3, theta = 650, shape = 2, rate = 1,
                           lower.tail = FALSE, log.p = TRUE), ref_d, tolerance = 1e-12)
  expect_equal(qcomppareto(-exp(-500), "gamma", alpha = 3, theta = 650, shape = 2, rate = 1,
                           log.p = TRUE), ref_d, tolerance = 1e-12)
  ## exponential and lognormal bodies
  q <- c(30, 35, 39.9, 40)
  ref_e <- c(-29.998230968398764, -34.766684413330222, -36.308494723425435, -36.311120545886064)
  expect_equal(pcomppareto(q, "exp", alpha = 2, theta = 40, rate = 1, lower.tail = FALSE,
                           log.p = TRUE), ref_e, tolerance = 1e-13)
  expect_equal(qcomppareto(ref_e, "exp", alpha = 2, theta = 40, rate = 1, lower.tail = FALSE,
                           log.p = TRUE), q, tolerance = 1e-12)
  ref_l <- c(-27.380851895104831, -32.776670242731888, -32.918938533204677)
  expect_equal(pcomppareto(exp(c(7, 7.9, 8)), "lnorm", alpha = 2, theta = exp(8), meanlog = 0,
                           sdlog = 1, lower.tail = FALSE, log.p = TRUE), ref_l, tolerance = 1e-13)
  expect_equal(log(qcomppareto(ref_l, "lnorm", alpha = 2, theta = exp(8), meanlog = 0, sdlog = 1,
                               lower.tail = FALSE, log.p = TRUE)), c(7, 7.9, 8), tolerance = 1e-12)
})

test_that("VB5-02: the EGPD gamlss start cache is keyed on the full response", {
  skip_if_not_installed("gamlss")
  ctl <- gamlss::gamlss.control(trace = FALSE)
  ## A and B share n, sum, sum of squares, sum(i y), and first and last values
  common <- rep(0:3, c(10, 30, 30, 10))
  A <- c(common, rep(c(0, 3, 0, 0, 3, 0), 10))
  B <- c(common, rep(c(0, 1, 1, 4, 0, 0), 10))
  assign("last", NULL, envir = .egpd_start_cache)
  gamlss::gamlss(y ~ 1, family = DEGPD1(sigma.link = "identity"), data = data.frame(y = A),
                 control = ctl)
  fB <- gamlss::gamlss(y ~ 1, family = DEGPD1(sigma.link = "identity"), data = data.frame(y = B),
                       control = ctl)
  ## independent MLE of B: closed-form pmf H(y + 1) - H(y),
  ## H(x) = [1 - (1 + xi x / s)^(-1/xi)]^k, maximised by optim
  expect_true(fB$converged)
  expect_equal(as.numeric(logLik(fB)), -204.275405945949, tolerance = 1e-9)
})

test_that("VB5-02: the GPIG start cache is keyed on the full response", {
  skip_if_not_installed("gamlss")
  ctl <- gamlss::gamlss.control(trace = FALSE)
  ## GPIG: independent MLE from the pgf recursion log G(s) = b[(1 - c)^a - (1 - cs)^a]
  Ag <- rep(c(0, 3, 0, 0, 3, 0), 10)
  Bg <- rep(c(0, 1, 1, 4, 0, 0), 10)
  assign("last", NULL, envir = .gpig_start_cache)
  suppressWarnings(gamlss::gamlss(y ~ 1, family = GPIG(), data = data.frame(y = Ag), control = ctl))
  gB <- gamlss::gamlss(y ~ 1, family = GPIG(), data = data.frame(y = Bg), control = ctl)
  expect_equal(as.numeric(logLik(gB)), -83.0917681961, tolerance = 1e-9)
})

test_that("VB5-10: gamlss starts for data on a large scale keep a bounded tail's endpoint", {
  skip_if_not_installed("gamlss")
  set.seed(5)
  y <- 2e8 * rbeta(300, 2, 1.5)
  assign("last", NULL, envir = .egpd_start_cache)
  st <- .egpd_gamlss_start(y, "EGPD1", sigma_link = "identity")
  expect_lt(max(y), -st[["mu"]] / st[["sigma"]])
  f <- gamlss::gamlss(y ~ 1, family = EGPD1(sigma.link = "identity"), data = data.frame(y = y),
                      control = gamlss::gamlss.control(trace = FALSE))
  ## independent MLE of the closed-form EGPD1 likelihood (kappa H^(kappa - 1) h, H the GPD
  ## cdf), maximised by optim: scale 1.3624e8, xi -0.67993, kappa 2.1012
  expect_true(f$converged)
  expect_equal(as.numeric(logLik(f)), -5700.943874341949, tolerance = 1e-10)
})
