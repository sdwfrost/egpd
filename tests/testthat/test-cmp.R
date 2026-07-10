## Conway-Maxwell-Poisson (src/cmp.cpp), mean-parameterised.
##
## Every quantity is checked against an INDEPENDENT computation:
##   * Z, E[Y], Var[Y] against a brute-force sum in R;
##   * nu = 1 against dpois, and nu -> 0 against the geometric limit;
##   * the analytic score against Richardson-extrapolated finite differences;
##   * and, most importantly, the GRID path against the exact per-observation path -- the grid
##     is the entire reason CMP is affordable here, so it must not change the answer.

test_that("the exact sums match a brute-force computation", {
  brute <- function(ll, nu, J = 3e5) {
    j <- 0:J
    f <- j * ll - nu * lgamma(j + 1)
    mx <- max(f); w <- exp(f - mx); S <- sum(w)
    EY <- sum(j * w) / S
    list(logZ = mx + log(S), EY = EY, VarY = sum(j^2 * w) / S - EY^2)
  }
  for (pr in list(c(0.5, 1), c(2, 0.5), c(-1, 0.2), c(1e-4, 8.3e-5))) {
    ll <- pr[1]; nu <- pr[2]
    a <- egpd:::cmp_sums_cpp(ll, nu); b <- brute(ll, nu)
    expect_equal(a[["logZ"]], b$logZ, tolerance = 1e-8)
    expect_equal(a[["EY"]],   b$EY,   tolerance = 1e-6)
    expect_equal(a[["VarY"]], b$VarY, tolerance = 1e-5)
  }
})

test_that("nu = 1 is exactly the Poisson", {
  for (mu in c(0.5, 3, 50, 500)) {
    y <- 0:min(2000, ceiling(mu + 12 * sqrt(mu)))
    expect_equal(dcmp(y, mu = mu, nu = 1), stats::dpois(y, mu), tolerance = 1e-7)
  }
})

test_that("nu -> 0 approaches the geometric with the same mean", {
  # As nu -> 0 (with lambda < 1) CMP -> geometric, whose dispersion index is 1 + mu.
  for (mu in c(2, 10)) {
    y <- 0:200
    expect_equal(dcmp(y, mu = mu, nu = 1e-6), stats::dgeom(y, prob = 1 / (1 + mu)),
                 tolerance = 1e-3)
  }
})

test_that("the mean parameterisation delivers E[Y] = mu", {
  for (pr in list(c(5, 1), c(5, 0.3), c(50, 0.1), c(200, 0.02))) {
    mu <- pr[1]; nu <- pr[2]
    y <- 0:2e5
    p <- dcmp(y, mu = mu, nu = nu)
    expect_equal(sum(p), 1, tolerance = 1e-5)
    expect_equal(sum(y * p), mu, tolerance = 1e-3 * mu)
  }
})

test_that("the tail is lighter than geometric for every nu > 0", {
  # log p ~ -nu y log y, so the log-pmf must fall away faster than any straight line.
  mu <- 20; nu <- 0.3
  y <- c(200, 400, 800)
  lp <- dcmp(y, mu = mu, nu = nu, log = TRUE)
  s1 <- (lp[2] - lp[1]) / (y[2] - y[1])
  s2 <- (lp[3] - lp[2]) / (y[3] - y[2])
  expect_lt(s2, s1)                                   # slope steepens: sub-geometric decay
})

test_that("the analytic score matches Richardson-extrapolated finite differences", {
  rich <- function(y, mu, nu, i) {
    f <- function(h) {
      e <- c(log(mu), log(nu)); e1 <- e2 <- e; e1[i] <- e1[i] + h; e2[i] <- e2[i] - h
      (egpd:::cmp_grad_cpp(y, exp(e1[1]), exp(e1[2]))[["logp"]] -
       egpd:::cmp_grad_cpp(y, exp(e2[1]), exp(e2[2]))[["logp"]]) / (2 * h)
    }
    (4 * f(1e-4) - f(2e-4)) / 3
  }
  for (pr in list(c(10, 0.6), c(60, 0.2))) for (y in c(0, 1, 9, 80)) {
    mu <- pr[1]; nu <- pr[2]
    g <- egpd:::cmp_grad_cpp(y, mu, nu)
    expect_equal(unname(g[["g_logmu"]]), rich(y, mu, nu, 1), tolerance = 1e-4)
    expect_equal(unname(g[["g_lognu"]]), rich(y, mu, nu, 2), tolerance = 1e-4)
  }
})

test_that("the grid path reproduces the exact per-observation path", {
  # The load-bearing check: nu being intercept-only lets us grid log lambda once per
  # likelihood evaluation instead of solving per observation. If the grid changed the
  # likelihood by more than rounding, every CMP AIC in the comparison would be wrong.
  set.seed(5)
  n  <- 300
  y  <- as.numeric(stats::rpois(n, 20))
  X1 <- cbind(1, as.numeric(scale(seq_len(n))))        # mu varies across observations
  X2 <- matrix(1, n, 1)                                # nu intercept-only -> grid path
  pars <- list(c(log(20), 0.15), log(0.6))
  off  <- list(numeric(0), numeric(0))
  nll_grid <- egpd:::cmp1d0(pars, X1, X2, y, as.integer(0), 0L, off)

  e1 <- as.vector(X1 %*% pars[[1]]); nu <- exp(pars[[2]])
  nll_slow <- -sum(vapply(seq_len(n), function(i)
    egpd:::cmp_logpmf_cpp(y[i], exp(e1[i]), nu), 0))
  expect_equal(nll_grid, nll_slow, tolerance = 1e-5)
})

test_that("CMP fits a GAM through egpd() and recovers the Poisson case", {
  set.seed(9)
  d <- data.frame(y = stats::rpois(500, 12))
  f <- egpd(list(lmu = y ~ 1, lnu = ~1), data = d, family = "cmp")
  expect_true(is.finite(f$AIC) && f$AIC > 0)
  expect_lte(as.numeric(f$logLik), 0)
  pr <- as.data.frame(predict(f, newdata = d[1, , drop = FALSE], type = "response"))
  expect_named(pr, c("mu", "nu"))
  expect_equal(pr$mu, mean(d$y), tolerance = 0.05 * mean(d$y))
  expect_equal(pr$nu, 1, tolerance = 0.3)             # Poisson data => nu near 1
})

test_that("CMP recovers an overdispersed fit", {
  set.seed(4)
  d <- data.frame(y = stats::rnbinom(600, size = 2, mu = 15))
  f <- egpd(list(lmu = y ~ 1, lnu = ~1), data = d, family = "cmp")
  pr <- as.data.frame(predict(f, newdata = d[1, , drop = FALSE], type = "response"))
  expect_lt(pr$nu, 1)                                 # overdispersion => nu < 1
  expect_equal(pr$mu, mean(d$y), tolerance = 0.2 * mean(d$y))
})
