## Poisson Ramos-Louzada (src/prl.cpp), Alkhairy (2023) MBE 20(8):14061-14080.
##
## The pmf is checked against INDEPENDENT statements of the same distribution:
##   * the paper's own closed-form survival function (their Eq. 4), which the pmf code never
##     uses, so agreement is a genuine cross-check;
##   * the defining Poisson-RL mixture (their Eq. 1), integrated numerically;
##   * the moment identities derived from the pmf.

test_that("PRL sums to one and matches its derived moments", {
  # E[X] = tau^2/(tau-1);  Var[X] = tau^2 (tau^2 + tau - 3)/(tau-1)^2
  for (tau in c(2.0001, 2.5, 5, 20, 200)) {
    x <- 0:2e6
    p <- dprl(x, tau); m <- sum(x * p)
    expect_equal(sum(p), 1, tolerance = 1e-8)
    expect_equal(m, tau^2 / (tau - 1), tolerance = 1e-5)
    expect_equal(sum((x - m)^2 * p), tau^2 * (tau^2 + tau - 3) / (tau - 1)^2, tolerance = 1e-4)
  }
})

test_that("PRL is the Poisson mixture of a Ramos-Louzada intensity", {
  # X | theta ~ Poisson(theta), theta ~ RL(tau) with density (their Eq. 1)
  #   f(theta; tau) = (tau^2 - 2 tau + theta) exp(-theta/tau) / (tau^2 (tau - 1))
  for (tau in c(2.5, 5, 20)) {
    for (x in c(0, 1, 5, 40)) {
      num <- stats::integrate(function(th)
        stats::dpois(x, th) * (tau^2 - 2 * tau + th) * exp(-th / tau) / (tau^2 * (tau - 1)),
        0, Inf, rel.tol = 1e-12)$value
      expect_equal(dprl(x, tau), num, tolerance = 1e-8)
    }
  }
})

test_that("the pmf agrees with the paper's closed-form survival function", {
  # S(x) = P(X > x) = (1+1/tau)^-x tau (x + tau^2) / ((tau-1)(1+tau)^2)   [their Eq. 3/4]
  # It is the STRICT survivor, so P(X = x) = S(x-1) - S(x) and F(x) = 1 - S(x).
  S <- function(x, tau) (1 + 1 / tau)^(-x) * tau * (x + tau^2) / ((tau - 1) * (1 + tau)^2)
  expect_equal(S(-1, 2.5), 1)                       # P(X > -1) = 1, a check on the formula
  for (tau in c(2.5, 5, 50)) for (x in c(0, 1, 5, 40)) {
    # an identity the pmf code never uses
    expect_equal(dprl(x, tau), S(x - 1, tau) - S(x, tau), tolerance = 1e-10)
    expect_equal(pprl(x, tau), 1 - S(x, tau), tolerance = 1e-10)
  }
})

test_that("PRL's mean cannot go below 4, and its dispersion is pinned to the mean", {
  # mu(tau) = tau^2/(tau-1) has d mu/d tau = tau(tau-2)/(tau-1)^2, zero at tau = 2.
  mu <- function(tau) tau^2 / (tau - 1)
  expect_equal(mu(2), 4)
  expect_true(all(mu(seq(2, 50, length.out = 200)) >= 4 - 1e-12))
  expect_true(all(diff(mu(seq(2.001, 50, length.out = 200))) > 0))    # increasing after tau=2
  # the dispersion index is a function of the mean alone -- no free dispersion parameter
  DI <- function(tau) (tau^2 + tau - 3) / (tau - 1)
  expect_equal(DI(2), 3)
  for (tau in c(3, 10, 100)) expect_lt(abs(DI(tau) / mu(tau) - 1), 0.35)
})

test_that("the analytic score matches Richardson-extrapolated finite differences", {
  rich <- function(y, eta) {
    f <- function(h) (egpd:::prl_grad_cpp(y, eta + h)[["logp"]] -
                      egpd:::prl_grad_cpp(y, eta - h)[["logp"]]) / (2 * h)
    (4 * f(1e-4) - f(2e-4)) / 3
  }
  for (eta in log(c(4.5, 10, 500, 5000))) for (y in c(0, 1, 7, 200, 5000)) {
    g <- egpd:::prl_grad_cpp(y, eta)[["g"]]
    expect_equal(unname(g), rich(y, eta), tolerance = 1e-5)
  }
})

test_that("q/p invert one another", {
  for (tau in c(2.5, 10)) for (pr in c(0.1, 0.5, 0.9)) {
    x <- qprl(pr, tau)
    expect_gte(pprl(x, tau), pr - 1e-9)
    if (x > 0) expect_lt(pprl(x - 1, tau), pr)
  }
})

test_that("PRL fits a GAM through egpd(), and predict() respects the mu >= 4 clamp", {
  set.seed(3)
  d <- data.frame(y = rprl(400, tau = 8))
  f <- egpd(list(lmu = y ~ 1), data = d, family = "prl")
  expect_true(is.finite(f$AIC) && f$AIC > 0)
  expect_lte(as.numeric(f$logLik), 0)                 # a pmf log-likelihood must be <= 0
  pr <- as.data.frame(predict(f, newdata = d[1, , drop = FALSE], type = "response"))
  expect_named(pr, "mu")
  expect_gte(pr$mu, egpd:::prl_bounds_cpp()[["mu_lo"]] - 1e-12)

  # Data with a mean well below PRL's floor of 4. The fit must stay feasible and not diverge.
  # The optimum need NOT sit on the floor: PRL's shape is not driven by its mean alone, and on
  # Poisson(0.7) data the MLE is interior at mu = 4.256 (confirmed against a 1-D optimize()).
  d2 <- data.frame(y = stats::rpois(300, 0.7))
  f2 <- egpd(list(lmu = y ~ 1), data = d2, family = "prl")
  pr2 <- as.data.frame(predict(f2, newdata = d2[1, , drop = FALSE], type = "response"))
  mu_lo <- egpd:::prl_bounds_cpp()[["mu_lo"]]
  expect_true(is.finite(pr2$mu))
  expect_gte(pr2$mu, mu_lo - 1e-12)
  nll <- function(mu) { tau <- (mu + sqrt(mu^2 - 4 * mu)) / 2; -sum(dprl(d2$y, tau, log = TRUE)) }
  expect_equal(pr2$mu, stats::optimize(nll, c(mu_lo, 20))$minimum, tolerance = 1e-3)
})

test_that("the likelihood is clamped, not extrapolated, below mu = 4", {
  # Below the floor the pmf is evaluated at mu_lo, so log p must be exactly flat there --
  # otherwise the optimiser would be reading a distribution that does not exist.
  mu_lo <- egpd:::prl_bounds_cpp()[["mu_lo"]]
  for (y in c(0, 1, 9)) {
    ref <- egpd:::prl_grad_cpp(y, log(mu_lo))[["logp"]]
    for (mu in c(3.9, 2.0, 0.5)) expect_equal(egpd:::prl_grad_cpp(y, log(mu))[["logp"]], ref)
  }
})
