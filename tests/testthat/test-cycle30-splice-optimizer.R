test_that("a Lomax splice keeps its density under changes of units", {
  a <- 3
  tau <- 1.25
  r <- .5
  t <- c(2.5, 2.875, 3.5)
  odds <- 2*tau*r/(a*expm1(r*tau))
  unit <- odds/(1+odds) * a/(t+tau) * (2*tau/(t+tau))^a
  for (b in c(1, 2^1022)) {
    d <- function(log) dcomppareto(t*b, "exp", theta = tau*b, rate = r/b,
                                   alpha = a, log = log)
    expect_lt(max(abs(d(TRUE) - (log(unit) - log(b)))), 3e-12)
    expect_equal(d(FALSE)/(unit/b), rep(1, length(t)), tolerance = 3e-12)
  }
  th <- 1200000000L
  xx <- 1500000000L
  d <- function(x, theta) dcomppareto(x, "exp", theta = theta, rate = .5/th, log = TRUE)
  expect_identical(d(xx, th), d(as.double(xx), as.double(th)))

  # In these small units only the *unweighted* conditional density overflows.
  b <- 2^-1019
  t <- 1 + 2^-16
  a <- 256
  w <- 2/(2 + a*expm1(1))
  unit <- w*a/(t+1) * (2/(t+1))^a
  expect_equal(b*dcomppareto(t*b, "exp", theta = b, rate = 1/b, alpha = a),
               unit, tolerance = 3e-12)
})

test_that("POT L-BFGS-B does not report a support-blocked line search as a maximum", {
  n <- c(16, 48, 36)
  x <- 11 + rep(0:2, n)
  fit <- function(...) fit_degpd_pot(x, 11, fix.kappa = FALSE, fix.arg = list(xi = -1),
                                      start = list(sigma = 2.625, kappa = 1.75),
                                      optim.method = "L-BFGS-B", ...)
  # F(0)=1/sigma^kappa and F(1)/F(0)=2^kappa recover both parameters.
  f <- fit()
  expect_equal(f$convergence, 0)
  expect_equal(coef(f), c(sigma = 2.5, kappa = 2), tolerance = 5e-5)
  expect_lt(abs(f$loglik - sum(n*log(n/sum(n)))), 1e-6)

  # The improving-step check cannot cross an explicit lower scale bound.
  expect_warning(g <- fit(lower = log(c(2.56, .75)), upper = log(c(2.8, 3))),
                 "optimizer bound")
  expect_equal(unname(coef(g)["sigma"]), 2.56)
  expect_true(all(is.na(g$sd)))
})

test_that("a feasible retry recovers an active box with the default start", {
  set.seed(42)
  x <- 2*((1-runif(200)^(1/1.8))^(-.15)-1)/.15
  expect_warning(f <- distfit(x, lower = c(.5, -.2, .8), upper = c(10, 1, 1.5)),
                 "optimizer bound")
  # A feasible reference on the lower power bound, from the continuous
  # GPD density times the derivative of its power carrier.
  s <- 1.805911; xi <- .155329; k <- exp(.8)
  logS <- -log1p(xi*x/s)/xi
  ref <- sum(log(k/s) + (1+xi)*logS + (k-1)*log(-expm1(logS)))
  expect_equal(f$convergence, 0)
  expect_gte(f$loglik, ref-1e-5)
  expect_equal(unname(coef(f)["kappa"]), k)
  expect_true(all(is.na(f$sd)))
})

test_that("a decreasing continuous beta-carrier likelihood has no interior maximum", {
  x <- seq(.01, .12, length.out = 50)
  fit <- function(x, ...) distfit(x, type = 3, fix.arg = list(sigma = 1, xi = -1), ...)
  # At k=0, logit(z) is uniform on [-log(31),0]. The natural score is
  # sum(T-E[T]) for T=log(z*(1-z)) and decreases strictly with k.
  z <- 1/32 + 15*x/32
  m <- integrate(function(v) log(plogis(v)*plogis(-v)), -log(31), 0)$value/log(31)
  expect_lt(mean(log(z*(1-z))), m)
  expect_error(fit(x), "No interior type-3 kappa")
  expect_warning(b <- fit(x, lower = log(.2), upper = log(4)), "optimizer bound")
  expect_equal(unname(coef(b)), .2)
  g <- fit(seq(.15, .55, length.out = 55))
  expect_equal(unname(coef(g)), .77602152, tolerance = 5e-4)
  expect_true(is.finite(g$sd))
})
