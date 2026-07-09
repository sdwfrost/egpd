## Generalized Poisson, generalized Waring and Poisson-lognormal (src/countfams.cpp).
## Every pmf is checked against an INDEPENDENT computation, never against itself:
##   gpois    vs gamlss.dist::dGPO
##   gwaring  vs its beta-negative-binomial representation, integrated numerically
##   plnorm   vs direct quadrature of the Poisson-normal mixture

test_that("all three pmfs sum to one and have mean mu", {
  y <- 0:200000
  for (a in list(list(f = dgpois,      mu = 20,   p = list(lambda = 0.6)),
                 list(f = dgpois,      mu = 500,  p = list(lambda = 0.9)),
                 list(f = dpoislnorm,  mu = 20,   p = list(sigma = 1.2)),
                 list(f = dpoislnorm,  mu = 5000, p = list(sigma = 0.8)))) {
    p <- do.call(a$f, c(list(y, mu = a$mu), a$p))
    expect_equal(sum(p), 1, tolerance = 1e-5)
    expect_equal(sum(y * p), a$mu, tolerance = 1e-4)
  }
})

test_that("generalized Waring sums to one; its mean converges slowly because the tail is a power law", {
  # E[Y] needs a large grid precisely BECAUSE P(Y=y) ~ y^-(rho+1). That slow convergence is
  # the property the family was added for, not a defect.
  y <- 0:10000000
  p <- dgwaring(y, mu = 20, k = 1, rho = 2.5)
  expect_equal(sum(p), 1, tolerance = 1e-8)
  expect_equal(sum(y * p), 20, tolerance = 1e-4)
})

test_that("generalized Waring has a regularly varying tail with index -(rho+1)", {
  for (rho in c(1.5, 2.5, 4)) {
    y  <- c(1e5, 2e5)
    lp <- dgwaring(y, mu = 20, k = 1, rho = rho, log = TRUE)
    expect_equal(as.numeric(diff(lp) / diff(log(y))), -(rho + 1), tolerance = 1e-3)
  }
})

test_that("generalized Waring equals the beta-negative-binomial BNB(r=k, alpha=rho, beta=a)", {
  mu <- 20; k <- 1; rho <- 2.5
  a  <- mu * (rho - 1) / k
  for (y in c(0, 1, 5, 50, 500)) {
    num <- stats::integrate(function(pp) stats::dnbinom(y, size = k, prob = pp) *
                              stats::dbeta(pp, rho, a), 0, 1, rel.tol = 1e-12)$value
    expect_equal(dgwaring(y, mu = mu, k = k, rho = rho), num, tolerance = 1e-9)
  }
})

test_that("generalized Poisson matches gamlss.dist::dGPO", {
  skip_if_not_installed("gamlss.dist")
  # gamlss GPO has Var = mu (1 + sigma mu)^2, so lambda = sigma mu / (1 + sigma mu).
  for (pr in list(c(20, 0.05), c(500, 0.01))) {
    mu <- pr[1]; sg <- pr[2]; lam <- sg * mu / (1 + sg * mu)
    for (y in c(0, 3, 20, 200))
      expect_equal(dgpois(y, mu = mu, lambda = lam),
                   gamlss.dist::dGPO(y, mu = mu, sigma = sg), tolerance = 1e-10)
  }
})

test_that("generalized Poisson variance is mu/(1-lambda)^2", {
  y <- 0:300000
  p <- dgpois(y, mu = 20, lambda = 0.6); m <- sum(y * p)
  expect_equal(sum((y - m)^2 * p), 20 / (1 - 0.6)^2, tolerance = 1e-6)
})

test_that("Poisson-lognormal matches direct quadrature, including where integrate() fails", {
  mu <- 5000; sg <- 0.8; muz <- log(mu) - sg^2 / 2
  # The Poisson kernel is a spike of width ~1/sqrt(y). Integrating over muz +/- 12 sd
  # misses it entirely and returns ~1e-53 at y = 15000; the adaptive rule must not.
  for (y in c(2500, 5000, 15000)) {
    zpk   <- log(y)
    tight <- stats::integrate(function(z) stats::dpois(y, exp(z)) * stats::dnorm(z, muz, sg),
                              zpk - 0.3, zpk + 0.3, rel.tol = 1e-12)$value
    expect_equal(dpoislnorm(y, mu = mu, sigma = sg), tight, tolerance = 1e-6)
  }
  # and P(Y=0) at parameters where zeros actually occur
  for (pr in list(c(15, 1.5), c(5, 1.0), c(2, 0.8))) {
    m <- pr[1]; s <- pr[2]; mz <- log(m) - s^2 / 2
    int <- stats::integrate(function(z) exp(-exp(z)) * stats::dnorm(z, mz, s),
                            mz - 15 * s, mz + 8 * s, rel.tol = 1e-13)$value
    expect_equal(dpoislnorm(0, mu = m, sigma = s), int, tolerance = 1e-5)
  }
})

test_that("Poisson-lognormal variance is mu + mu^2 (exp(sigma^2) - 1)", {
  y <- 0:300000
  p <- dpoislnorm(y, mu = 20, sigma = 1.2); m <- sum(y * p)
  expect_equal(sum((y - m)^2 * p), 20 + 400 * (exp(1.44) - 1), tolerance = 1e-5)
})

test_that("analytic scores match Richardson-extrapolated finite differences", {
  # Plain central differences are limited by roundoff here (the error GROWS as h -> 0, and a
  # near-zero gradient component inflates any relative measure), so the reference is a
  # Richardson-extrapolated derivative, accurate to ~1e-7.
  rich <- function(fam, y, eta, i) {
    f <- function(h) {
      e1 <- e2 <- eta; e1[i] <- e1[i] + h; e2[i] <- e2[i] - h
      (egpd:::cf_grad_cpp(fam, y, e1)[1] - egpd:::cf_grad_cpp(fam, y, e2)[1]) / (2 * h)
    }
    (4 * f(1e-4) - f(2e-4)) / 3
  }
  cases <- list(list(fam = "gpois",   eta = c(log(500), stats::qlogis(0.7))),
                list(fam = "plnorm",  eta = c(log(500), log(0.9))),
                list(fam = "gwaring", eta = c(log(500), log(1.3), log(2.2))))
  for (cs in cases) for (y in c(0, 1, 7, 200, 5000)) {
    g <- egpd:::cf_grad_cpp(cs$fam, y, cs$eta)[-1]
    r <- vapply(seq_along(cs$eta), function(i) rich(cs$fam, y, cs$eta, i), 0)
    expect_equal(unname(g), unname(r), tolerance = 1e-5)
  }
})

test_that("q/p functions invert one another", {
  for (cs in list(list(q = qgpois,     p = pgpois,     a = list(mu = 30, lambda = 0.5)),
                  list(q = qgwaring,   p = pgwaring,   a = list(mu = 30, k = 1, rho = 3)),
                  list(q = qpoislnorm, p = ppoislnorm, a = list(mu = 30, sigma = 0.9)))) {
    for (pr in c(0.1, 0.5, 0.9)) {
      x <- do.call(cs$q, c(list(pr), cs$a))
      expect_gte(do.call(cs$p, c(list(x), cs$a)), pr - 1e-9)
      if (x > 0) expect_lt(do.call(cs$p, c(list(x - 1), cs$a)), pr)
    }
  }
})

test_that("each family fits a GAM through egpd() and predicts on the response scale", {
  set.seed(7)
  n <- 400
  d <- data.frame(t = seq(0, 1, length.out = n))
  d$y_gp <- rgpois(n, mu = 30, lambda = 0.5)
  d$y_gw <- rgwaring(n, mu = 30, k = 1, rho = 3)
  d$y_pl <- rpoislnorm(n, mu = 30, sigma = 0.8)

  f1 <- egpd(list(lmu = y_gp ~ 1, logitlambda = ~1), data = d, family = "gpois")
  f2 <- egpd(list(lmu = y_gw ~ 1, lk = ~1, lrho = ~1), data = d, family = "gwaring")
  f3 <- egpd(list(lmu = y_pl ~ 1, lsigma = ~1), data = d, family = "plnorm")
  for (f in list(f1, f2, f3)) {
    expect_true(is.finite(f$AIC) && f$AIC > 0)
    expect_lte(as.numeric(f$logLik), 0)              # a pmf log-likelihood must be <= 0
  }
  pr <- as.data.frame(predict(f1, newdata = d[1, , drop = FALSE], type = "response"))
  expect_named(pr, c("mu", "lambda"))
  pr <- as.data.frame(predict(f2, newdata = d[1, , drop = FALSE], type = "response"))
  expect_named(pr, c("mu", "k", "rho"))
  expect_gt(pr$rho, 1)                                # rho > 1 is enforced by the clamp
  pr <- as.data.frame(predict(f3, newdata = d[1, , drop = FALSE], type = "response"))
  expect_named(pr, c("mu", "sigma"))
})

test_that("predict() reports the same parameters the likelihood used (clamps applied)", {
  # rho enters the gwaring likelihood only through its clamp, so beyond the bound the
  # objective is flat and the optimiser leaves eta3 unbounded. Unclamped, predict() reported
  # rho = Inf (and 4.8e220) on real Tycho fits while the model had used rho = 1e6 -- a
  # boundary fit that looked converged. The bounds come from C++ so they cannot drift.
  bd <- egpd:::cf_bounds_cpp()
  expect_gt(bd$gwaring[["rho_lo"]], 1)               # rho > 1 keeps the mean finite
  expect_true(all(vapply(bd, function(b) all(is.finite(b)), TRUE)))

  set.seed(11)
  # Poisson data has no heavy tail, so gwaring must escape to its rho -> Inf (NB) limit.
  d <- data.frame(y = stats::rpois(400, 30))
  f <- egpd(list(lmu = y ~ 1, lk = ~1, lrho = ~1), data = d, family = "gwaring")
  pr <- as.data.frame(predict(f, newdata = d[1, , drop = FALSE], type = "response"))
  expect_true(is.finite(pr$rho))
  expect_lte(pr$rho, bd$gwaring[["rho_hi"]])
  expect_gte(pr$rho, bd$gwaring[["rho_lo"]])
  expect_true(is.finite(1 / pr$rho))                 # xi is then reportable, and ~0
})

test_that("gwaring collapses onto the negative binomial as rho -> Inf", {
  # BNB(r=k, alpha=rho, beta=a) with a = mu(rho-1)/k tends to NB(size=k, mu) as rho -> Inf.
  # This is the limit the optimiser escapes to on light-tailed data, so it must be right.
  mu <- 30; k <- 2
  for (y in c(0, 5, 30, 200))
    expect_equal(dgwaring(y, mu = mu, k = k, rho = 1e7),
                 stats::dnbinom(y, size = k, mu = mu), tolerance = 1e-4)
})
