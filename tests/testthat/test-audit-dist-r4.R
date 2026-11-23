## Round-4 distribution fixes (VB-01, 06-10, 14-17, 20, 23). Expected values
## are closed forms or high-precision (Rmpfr) values computed independently
## and hard-coded here.

## Relative comparison: expect_equal() compares absolutely when the expected
## value is below its tolerance, which would accept 0 for 1e-46.
expect_rel <- function(object, expected, tol) {
  expect_lte(max(abs(object / expected - 1)), tol)
}

test_that("VB-01: discrete and ZI quantiles invert lower log.p below -745", {
  ## carrier 1, xi = 0: H(y) = (1 - exp(-y / sigma))^kappa and X = floor(Y),
  ## so F_X(j) = H(j + 1); the quantile is the smallest j with
  ## kappa * log1p(-exp(-(j + 1) / sigma)) >= lp
  s <- 10; k <- 1000
  ref_disc <- function(lp) {
    j <- 0
    while (k * log1p(-exp(-(j + 1) / s)) < lp) j <- j + 1
    j
  }
  for (lp in c(-745.2, -800, -1000, -2000)) {
    expect_equal(qdiscegpd(lp, sigma = s, xi = 0, kappa = k, log.p = TRUE), ref_disc(lp))
    expect_equal(qzidiscegpd(lp, pi = 0, sigma = s, xi = 0, kappa = k, log.p = TRUE), ref_disc(lp))
    expect_equal(qziegpd(lp, pi = 0, sigma = s, xi = 0, kappa = k, log.p = TRUE),
                 -s * log1p(-exp(lp / k)), tolerance = 1e-13)
  }
  expect_equal(qdiscegpd(-1000, sigma = s, xi = 0, kappa = k, log.p = TRUE), 4)
  expect_equal(qziegpd(-1000, pi = 0, sigma = s, xi = 0, kappa = k, log.p = TRUE),
               -10 * log1p(-exp(-1)), tolerance = 1e-14)
  ## subnormal linear probability (exp(-744) = 9.9e-324)
  expect_equal(qziegpd(-744, pi = 0, sigma = s, xi = 0, kappa = k, log.p = TRUE),
               -10 * log1p(-exp(-0.744)), tolerance = 1e-12)
  ## round trip through the log CDF
  lF4 <- pdiscegpd(4, sigma = s, xi = 0, kappa = k, log.p = TRUE)
  expect_equal(lF4, k * log1p(-exp(-5 / s)), tolerance = 1e-13)
  expect_equal(qdiscegpd(lF4, sigma = s, xi = 0, kappa = k, log.p = TRUE), 4)
  ## pi = 1e-322 (subnormal) with lp = -740: log((p - pi) / (1 - pi)) =
  ## -740.2690589684954516 (Rmpfr); continuous quantile 6.4814620728934033
  expect_equal(qziegpd(-740, pi = 1e-322, sigma = s, xi = 0, kappa = k, log.p = TRUE),
               6.4814620728934033, tolerance = 1e-12)
  expect_equal(qzidiscegpd(-740, pi = 1e-322, sigma = s, xi = 0, kappa = k, log.p = TRUE),
               ref_disc(-740.26905896849545))
  ## below pi on the log scale the quantile is zero
  expect_equal(qziegpd(-745, pi = 1e-322, sigma = s, xi = 0, kappa = k, log.p = TRUE), 0)
})

test_that("VB-06: type-1 qegpd is accurate for small kappa with p > 1/2", {
  ## x = -sigma * log1p(-p^(1 / kappa)) for xi = 0
  expect_rel(qegpd(0.9, kappa = 0.001, sigma = 1, xi = 0), -log1p(-0.9^1000), 1e-13)
  expect_rel(qegpd(0.6, kappa = 0.02, sigma = 1, xi = 0), -log1p(-0.6^50), 1e-13)
  expect_rel(qegpd(0.6, kappa = 0.005, sigma = 1, xi = 0), -log1p(-0.6^200), 1e-13)
  ## xi != 0: sigma / xi * ((1 - u)^(-xi) - 1)
  u <- 0.9^200
  expect_rel(qegpd(0.9, kappa = 0.005, sigma = 1, xi = 0.2), expm1(-0.2 * log1p(-u)) / 0.2, 1e-13)
  ## the same quantile from the upper tail on the log scale
  expect_rel(qegpd(log(0.1), kappa = 0.001, sigma = 1, xi = 0, lower.tail = FALSE, log.p = TRUE),
             -log1p(-0.9^1000), 1e-12)
  expect_equal(pegpd(qegpd(0.9, kappa = 0.001, sigma = 1, xi = 0), kappa = 0.001, sigma = 1, xi = 0),
               0.9, tolerance = 1e-12)
})

test_that("VB-07: carrier-5 quantiles and CDF below the double range of G4", {
  ## delta = 1: G4(u) = u^2, so G5(u) = u^kappa and Q5(p) = p^(1 / kappa)
  ## (the quantiles amplify the rounding of log(p) by 1 / kappa)
  expect_rel(q.G(1e-50, type = 5, kappa = 0.2, delta = 1), 1e-250, 1e-12)
  expect_rel(qegpd(1e-20, type = 5, kappa = 0.1, delta = 1, sigma = 1, xi = 0), 1e-200, 1e-12)
  expect_rel(qegpd(log(1e-20), type = 5, kappa = 0.1, delta = 1, sigma = 1, xi = 0, log.p = TRUE),
             1e-200, 1e-12)
  expect_rel(qegpd(0.5, type = 5, kappa = 0.0015, delta = 1, sigma = 1, xi = 0),
             -log1p(-0.5^(1 / 0.0015)), 1e-12)
  expect_rel(qegpd(0.6, type = 5, kappa = 0.001, delta = 1, sigma = 1, xi = 0),
             -log1p(-exp(log(0.6) / 0.001)), 1e-12)
  expect_rel(p.G(1e-200, type = 5, kappa = 0.1, delta = 1), 1e-20, 1e-13)
  expect_rel(p.G(1e-250, type = 5, kappa = 0.2, delta = 1), 1e-50, 1e-13)
})

test_that("VB-08: p.G carriers 2 and 3 are relatively accurate near u = 0", {
  ## carrier 3 with kappa = 1 is the uniform: G3(u) = u
  expect_rel(p.G(1e-20, type = 3, kappa = 1), 1e-20, 1e-13)
  ## first-order closed forms (the O(u^2) terms are 1e-20 relative):
  ## G2(u) ~ u sqrt(k) phi(sqrt(k)) / (Phi(sqrt(k)) - 1/2) and
  ## G3(u) ~ (15 / 32) u f(1/32) / (1/2 - F(1/32)) with F, f the Beta(2, 2) CDF and density
  expect_rel(p.G(1e-20, type = 2, kappa = 1), 1e-20 * dnorm(1) / (pnorm(1) - 0.5), 1e-13)
  F2 <- function(x) 3 * x^2 - 2 * x^3
  f2 <- function(x) 6 * x * (1 - x)
  expect_rel(p.G(1e-20, type = 3, kappa = 2), 1e-20 * (15 / 32) * f2(1 / 32) / (0.5 - F2(1 / 32)),
             1e-13)
  ## kappa < 0.01 (series branch of carrier 2): G2(u) ~ u * (1 - kappa / 3)
  k <- 0.005
  expect_rel(p.G(1e-20, type = 2, kappa = k),
             1e-20 * sqrt(k) * dnorm(sqrt(k)) / (pnorm(sqrt(k)) - 0.5), 1e-13)
  ## u = 1e-12: the linear difference keeps only about 4 digits
  expect_rel(p.G(1e-12, type = 2, kappa = 1), 1e-12 * dnorm(1) / (pnorm(1) - 0.5), 1e-11)
})

test_that("VB-09: carrier 3 stays accurate as kappa -> 0", {
  ## kappa -> 0 limit: G3(u) = log(x / (1 - x) * 31) / log(31) and
  ## g3(u) = (15 / 32) / (x (1 - x) log(31)), x = 1/32 + 15 u / 32 (error O(kappa))
  expect_equal(p.G(0.5, type = 3, kappa = 1e-15), log(527 / 47) / log(31), tolerance = 1e-13)
  expect_equal(d.G(0.5, type = 3, kappa = 1e-15, log = TRUE),
               log(15 / 32) - log((17 / 64) * (47 / 64)) - log(log(31)), tolerance = 1e-13)
  u <- 1 - 1.1^-5
  x <- 1 / 32 + 15 / 32 * u
  G <- (log(x / (1 - x)) + log(31)) / log(31)
  expect_equal(pegpd(0.5, type = 3, kappa = 1e-15, sigma = 1, xi = 0.2), G, tolerance = 1e-13)
  expect_equal(qegpd(G, type = 3, kappa = 1e-15, sigma = 1, xi = 0.2), 0.5, tolerance = 1e-12)
  expect_equal(qegpd(G, type = 3, kappa = 1e-10, sigma = 1, xi = 0.2), 0.5, tolerance = 1e-9)
  expect_equal(q.G(log(527 / 47) / log(31), type = 3, kappa = 1e-15), 0.5, tolerance = 1e-12)
})

test_that("VB-10: carrier 6 with a zero-weight component, and its log density far below 1", {
  ## prob = 0: G(u) = u^delta; prob = 1: G(u) = u^kappa
  expect_equal(d.G(0, type = 6, prob = 0, kappa = 0.5, delta = 1), 1)
  expect_equal(degpd_density(0, type = 6, prob = 0, kappa = 0.5, delta = 1, sigma = 2, xi = 0.1), 0.5)
  expect_equal(degpd_density(0, type = 6, prob = 1, kappa = 2, delta = 0.5, sigma = 1, xi = 0.1), 0)
  expect_equal(degpd_density(0, type = 6, prob = 1, kappa = 1, delta = 0.5, sigma = 4, xi = 0.1,
                             log = TRUE), -log(4))
  ## g(u) = prob kappa u^(kappa - 1) + (1 - prob) delta u^(delta - 1)
  expect_equal(d.G(1e-200, type = 6, prob = 0.5, kappa = 3, delta = 4, log = TRUE),
               log(1.5) - 400 * log(10), tolerance = 1e-14)
})

test_that("VB-14: zero-inflated log-probabilities never exceed 0", {
  ## xi = -0.3, sigma = 2: the base law ends at 6.67, so the CDF is 1 from 6 on
  expect_identical(pzidiscegpd(c(6, 12, 1e6, Inf), pi = 0.05, kappa = 1, sigma = 2, xi = -0.3,
                               log.p = TRUE), rep(0, 4))
  expect_identical(pziegpd(c(10, Inf), pi = 0.05, kappa = 1, sigma = 2, xi = -0.3, log.p = TRUE),
                   c(0, 0))
  ## sigma = 0.1: all mass at 0
  expect_identical(dzidiscegpd(0, pi = 0.05, kappa = 1, sigma = 0.1, xi = -0.3, log = TRUE), 0)
  pis <- seq(0.001, 0.999, by = 0.001)
  expect_true(all(pzidiscegpd(1e6, pi = pis, kappa = 1, sigma = 2, xi = 0.3, log.p = TRUE) <= 0))
  expect_true(all(pziegpd(Inf, pi = pis, kappa = 1, sigma = 2, xi = 0.3, log.p = TRUE) == 0))
  ## the round trip returns the largest support value 6
  lp <- pzidiscegpd(12, pi = 0.05, kappa = 1, sigma = 2, xi = -0.3, log.p = TRUE)
  expect_equal(qzidiscegpd(lp, pi = 0.05, kappa = 1, sigma = 2, xi = -0.3, log.p = TRUE), 6)
  lpc <- pziegpd(10, pi = 0.05, kappa = 1, sigma = 2, xi = -0.3, log.p = TRUE)
  expect_equal(qziegpd(lpc, pi = 0.05, kappa = 1, sigma = 2, xi = -0.3, log.p = TRUE), 2 / 0.3)
})

test_that("VB-16: regpd draws stay accurate and inside the support for huge kappa", {
  u <- c(0.1, 0.5, 0.9)
  ## inverse CDF of carrier 1 with xi = 0: -log(1 - u^(1 / kappa))
  ref <- -log(-expm1(log(u) / 1e16))
  expect_equal(regpd(3, kappa = 1e16, sigma = 1, xi = 0, unifsamp = u), ref, tolerance = 1e-13)
  expect_equal(regpd(3, kappa = 1e12, sigma = 1, xi = 0, unifsamp = u), -log(-expm1(log(u) / 1e12)),
               tolerance = 1e-13)
  ## carrier 6 with both components at kappa = 1e16, carrier 5 with delta = 1
  ## (G5(u) = u^kappa)
  expect_equal(regpd(3, type = 6, prob = 0.5, kappa = 1e16, delta = 1e16, sigma = 1, xi = 0,
                     unifsamp = u), ref, tolerance = 1e-13)
  expect_equal(regpd(3, type = 5, kappa = 1e16, delta = 1, sigma = 1, xi = 0, unifsamp = u), ref,
               tolerance = 1e-13)
  ## a draw that rounds to 1: u = 1 - 2^-53, kappa = 2 gives -log(1 - u^(1/2)) = 54 log(2)
  expect_equal(regpd(1, kappa = 2, sigma = 1, xi = 0, unifsamp = 1 - 2^-53), 54 * log(2),
               tolerance = 1e-12)
  set.seed(418)
  r <- regpd(2000, kappa = 1e16, sigma = 1, xi = 0)
  expect_true(all(is.finite(r)))
  ## xi = -0.5, sigma = 1: Y < 2 and P(X = 0) = 0.75^kappa = 0, so every draw is 1
  set.seed(418)
  expect_true(all(rdiscegpd(1000, kappa = 1e16, sigma = 1, xi = -0.5) == 1))
})

test_that("VB-17: log tails stay finite when q / sigma overflows", {
  ## kappa = 1 is the GPD: log S = -log1p(xi q / sigma) / xi
  lz <- log(0.2) + 400 * log(10)
  expect_equal(pegpd(1e200, sigma = 1e-200, xi = 0.2, kappa = 1, lower.tail = FALSE, log.p = TRUE),
               -lz / 0.2, tolerance = 1e-14)
  expect_equal(degpd_density(1e200, sigma = 1e-200, xi = 0.2, kappa = 1, log = TRUE),
               -log(1e-200) - 1.2 * lz / 0.2, tolerance = 1e-14)
  expect_equal(pegpd(1e200, sigma = 1e-200, xi = 0.2, kappa = 1), 1)
})

test_that("VB-15: a logit nu link starts inside (0, 1)", {
  skip_if_not_installed("gamlss")
  set.seed(2)
  y <- rdiscegpd(300, kappa = 2.5, sigma = 2, xi = 0.2)
  st <- .egpd_gamlss_start(y, "DEGPD1", nu_link = "logit")
  expect_true(st[["nu"]] > 0 && st[["nu"]] < 1)
  fit <- suppressWarnings(gamlss::gamlss(y ~ 1, family = DEGPD1(nu.link = "logit"),
                                         control = gamlss::gamlss.control(trace = FALSE)))
  ## the constrained maximum lies at the boundary kappa = 1, where the
  ## closed-form DEGPD1 likelihood maximised over (sigma, xi) is -750.157015
  ## (independent optim; H(y) = (1 - (1 + xi y / sigma)^(-1 / xi))^kappa)
  H <- function(t, s, x) ifelse(t <= 0, 0, 1 - (1 + x * t / s)^(-1 / x))
  nll <- function(th) -sum(log(H(y + 1, exp(th[1]), exp(th[2])) - H(y, exp(th[1]), exp(th[2]))))
  ref <- -optim(c(log(2), log(0.2)), nll, control = list(reltol = 1e-12))$value
  expect_gt(fitted(fit, "nu")[1], 0.99)
  expect_equal(as.numeric(logLik(fit)), ref, tolerance = 1e-5)
})

test_that("VB-20: DEGPD5 kappa-score for large kappa and for scale >> count", {
  ## Rmpfr (300-bit) central differences of the closed-form log-PMF
  ## log[G2(F(y + 1)) - G2(F(y))], F the GPD CDF, G2 the truncated-normal carrier
  fam <- DEGPD5()
  expect_equal(fam$dldv(1, mu = 10, sigma = 0.2, nu = 2500), -0.3379818482719184, tolerance = 1e-12)
  expect_equal(fam$dldv(0, mu = 1e10, sigma = 0.2, nu = 2), -0.3731479490481578, tolerance = 1e-12)
  expect_equal(fam$dldv(0, mu = 1e17, sigma = 0.2, nu = 2), -0.3731479490981577, tolerance = 1e-12)
  expect_equal(fam$dldv(0, mu = 50, sigma = 0.3, nu = 3000), -0.4806186920521228, tolerance = 1e-11)
  zfam <- ZIDEGPD5()
  expect_equal(zfam$dldv(1, mu = 10, sigma = 0.2, nu = 2500, tau = 0.3), -0.3379818482719184,
               tolerance = 1e-12)
})

test_that("VB-23: EGPD4 delta-score at y = 0 with tau = 1", {
  ## log f(0) = log(1 / sigma) + log(g5(0)), g5(0) = sqrt((1 + delta) / 2) for
  ## kappa = 1, so d log f(0) / d delta = 1 / (2 (1 + delta))
  fam <- EGPD4()
  expect_no_warning(v <- fam$dldv(0, mu = 2, sigma = 0.2, nu = 1.5, tau = 1))
  expect_equal(v, 1 / (2 * (1 + 1.5)), tolerance = 1e-14)
  expect_equal(fam$dldv(c(0, 0), mu = 3, sigma = 0.1, nu = c(0.5, 4), tau = 1), 1 / (2 * c(1.5, 5)),
               tolerance = 1e-14)
})
