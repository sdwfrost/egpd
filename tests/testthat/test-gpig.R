## Tests for the generalised Poisson-inverse Gaussian (GPIG; Zhu & Joe 2009)
## and zero-inflated GPIG (ZIGPIG), across three categories:
##   (1) comparison against reference/special-case implementations,
##   (2) mutual d/p/q consistency,
##   (3) analytic (C++) derivatives vs finite differences.

expit <- function(x) 1 / (1 + exp(-x))

# ---- (1) reference / special-case comparisons --------------------------------

test_that("GPIG reduces to Poisson(bc) at a = 1", {
  for (par in list(c(3, 0.6), c(1.5, 0.9), c(5, 0.2))) {
    b <- par[1]; c <- par[2]
    expect_equal(egpd::dgpig(0:15, a = 1, b = b, c = c),
                 stats::dpois(0:15, b * c), tolerance = 1e-10)
  }
})

test_that("GPIG reduces to PIG at a = 1/2 (matches gamlss.dist::dPIG)", {
  skip_if_not_installed("gamlss.dist")
  for (par in list(c(2.4, 0.7), c(1.8, 0.5), c(3.0, 0.85))) {
    b <- par[1]; c <- par[2]
    mu_pig  <- b * c / (2 * sqrt(1 - c))
    sig_pig <- 1 / (b * sqrt(1 - c))
    expect_equal(egpd::dgpig(0:40, a = 0.5, b = b, c = c),
                 gamlss.dist::dPIG(0:40, mu = mu_pig, sigma = sig_pig),
                 tolerance = 1e-10)
  }
})

test_that("closed-form mean and variance match the recursion", {
  a <- 0.35; b <- 2.4; c <- 0.9
  p <- egpd::dgpig(0:5000, a = a, b = b, c = c)
  mu_cf  <- a * b * c / (1 - c)^(1 - a)
  var_cf <- a * b * c * (1 - a * c) / (1 - c)^(2 - a)
  m1 <- sum((0:5000) * p); m2 <- sum((0:5000)^2 * p)
  expect_equal(m1, mu_cf, tolerance = 1e-6)
  expect_equal(m2 - m1^2, var_cf, tolerance = 1e-5)
})

test_that("mean-parameterisation dGPIG has the requested mean", {
  p <- egpd::dGPIG(0:20000, mu = 3, sigma = 0.45, nu = 0.8)
  expect_equal(sum((0:20000) * p), 3, tolerance = 1e-3)
})

# ---- (2) d/p/q consistency ---------------------------------------------------

test_that("GPIG pmf sums to 1 and CDF = cumsum(pmf)", {
  for (par in list(c(0.4, 2.4, 0.9), c(0.6, 1.5, 0.5), c(0.3, 3, 0.95))) {
    a <- par[1]; b <- par[2]; c <- par[3]
    d <- egpd::dgpig(0:2000, a, b, c)
    expect_equal(sum(d), 1, tolerance = 1e-8)
    expect_equal(egpd::pgpig(0:2000, a, b, c), cumsum(d), tolerance = 1e-10)
  }
})

test_that("qgpig inverts pgpig", {
  a <- 0.4; b <- 2.4; c <- 0.9
  k <- 0:25
  expect_equal(egpd::qgpig(egpd::pgpig(k, a, b, c), a, b, c), k)
})

test_that("ZIGPIG pmf sums to 1, mixes correctly, and q inverts p", {
  a <- 0.4; b <- 2.4; c <- 0.9; pv <- 0.25
  dz <- egpd::dzigpig(0:2000, a, b, c, pi = pv)
  expect_equal(sum(dz), 1, tolerance = 1e-8)
  expect_equal(dz[1], pv + (1 - pv) * egpd::dgpig(0, a, b, c), tolerance = 1e-12)
  expect_equal(egpd::pzigpig(0:2000, a, b, c, pi = pv), cumsum(dz), tolerance = 1e-10)
  k <- 0:25
  expect_equal(egpd::qzigpig(egpd::pzigpig(k, a, b, c, pi = pv), a, b, c, pi = pv), k)
})

test_that("log-density options agree with log of the density", {
  a <- 0.5; b <- 2; c <- 0.7
  expect_equal(egpd::dgpig(0:10, a, b, c, log = TRUE),
               log(egpd::dgpig(0:10, a, b, c)))
  expect_equal(egpd::dzigpig(0:10, a, b, c, pi = 0.2, log = TRUE),
               log(egpd::dzigpig(0:10, a, b, c, pi = 0.2)))
})

# ---- (3) analytic C++ derivatives vs finite differences ----------------------

## Call the native-fitter Rcpp stubs with an intercept-only design.
.d0 <- function(fn, e, y, np)
  do.call(fn, c(list(as.list(e)),
                rep(list(matrix(1, 1, 1)), np),
                list(y, as.integer(0), 0L, as.list(rep(list(numeric(0)), np)))))
.grad_fd <- function(d0fn, e, y) {
  h <- 1e-6
  sapply(seq_along(e), function(i) {
    ep <- e; em <- e; ep[i] <- ep[i] + h; em[i] <- em[i] - h
    (d0fn(ep, y) - d0fn(em, y)) / (2 * h)
  })
}

test_that("GPIG d12 gradient matches finite differences (both parameterisations)", {
  cfgs <- list(
    list(d0 = egpd:::gpig1d0,    d12 = egpd:::gpig1d12,    e = c(log(3), qlogis(0.4), qlogis(0.75))),
    list(d0 = egpd:::gpignat1d0, d12 = egpd:::gpignat1d12, e = c(qlogis(0.4), log(2.1), qlogis(0.75)))
  )
  for (cf in cfgs) for (y in c(0, 1, 5, 20)) {
    row <- .d0(cf$d12, cf$e, y, 3)
    g <- row[1:3]
    f <- .grad_fd(function(e, yy) .d0(cf$d0, e, yy, 3), cf$e, y)
    expect_equal(as.numeric(g), as.numeric(f), tolerance = 1e-5)
    ## BHHH Hessian is the outer product of the gradient
    expect_equal(row[4:9],
                 c(g[1]*g[1], g[1]*g[2], g[1]*g[3], g[2]*g[2], g[2]*g[3], g[3]*g[3]),
                 tolerance = 1e-12)
  }
})

test_that("ZIGPIG d12 gradient matches finite differences (both parameterisations)", {
  cfgs <- list(
    list(d0 = egpd:::zigpig1d0,    d12 = egpd:::zigpig1d12,    e = c(log(3), qlogis(0.4), qlogis(0.75), qlogis(0.2))),
    list(d0 = egpd:::zigpignat1d0, d12 = egpd:::zigpignat1d12, e = c(qlogis(0.4), log(2.1), qlogis(0.75), qlogis(0.2)))
  )
  for (cf in cfgs) for (y in c(0, 1, 5, 20)) {
    row <- .d0(cf$d12, cf$e, y, 4)
    g <- row[1:4]
    f <- .grad_fd(function(e, yy) .d0(cf$d0, e, yy, 4), cf$e, y)
    expect_equal(as.numeric(g), as.numeric(f), tolerance = 1e-5)
  }
})

# ---- fitting surfaces --------------------------------------------------------

test_that("native GAM fit matches the direct optim MLE (both parameterisations)", {
  skip_on_cran()
  set.seed(7)
  y <- egpd:::rgpig(2000, a = 0.45, b = 2.2, c = 0.75)
  dat <- data.frame(y = y)
  fm <- egpd(y ~ 1, data = dat, family = "gpig", trace = 0)
  fn <- egpd(y ~ 1, data = dat, family = "gpignat", trace = 0)
  ## the two parameterisations describe the same fitted distribution
  expect_equal(as.numeric(logLik(fm)), as.numeric(logLik(fn)), tolerance = 1e-3)
  ## and they reach the direct MLE (up to the REML/Laplace constant)
  nll <- function(par) -sum(egpd:::dgpig(y, expit(par[1]), exp(par[2]), expit(par[3]), log = TRUE))
  o <- optim(c(0, log(mean(y)), 0), nll, method = "BFGS")
  k <- 0.5 * 3 * log(2 * pi)   # Laplace constant for 3 parameters
  expect_equal(as.numeric(logLik(fm)) + k, -o$value, tolerance = 1e-2)
})

test_that("fitegpd recovers GPIG / ZIGPIG parameters and detects the PIG case", {
  skip_on_cran()
  set.seed(3)
  y <- egpd:::rgpig(3000, a = 0.45, b = 2.2, c = 0.8)
  f <- fitegpd(y, family = "gpig")
  expect_equal(unname(f$estimate["a"]), 0.45, tolerance = 0.1)
  expect_equal(unname(f$estimate["c"]), 0.8,  tolerance = 0.1)

  skip_if_not_installed("gamlss.dist")
  set.seed(5)
  yp <- gamlss.dist::rPIG(3000, mu = 3, sigma = 0.8)
  fp <- fitegpd(yp, family = "gpig")
  expect_equal(unname(fp$estimate["a"]), 0.5, tolerance = 0.08)  # a = 1/2 <=> PIG
})

test_that("gamlss and bamlss GPIG family constructors are well-formed", {
  fam <- egpd::GPIG()
  expect_s3_class(fam, "gamlss.family")
  expect_equal(fam$nopar, 3L)
  zfam <- egpd::ZIGPIG()
  expect_equal(zfam$nopar, 4L)

  gb <- egpd::gpig_bamlss()
  expect_s3_class(gb, "family.bamlss")
  expect_equal(gb$names, c("mu", "sigma", "nu"))
  ## d() delegates to the mean-parameterisation density
  b <- 3 * (1 - 0.8)^(1 - 0.45) / (0.45 * 0.8)
  expect_equal(gb$d(3, list(mu = 3, sigma = 0.45, nu = 0.8)),
               egpd::dgpig(3, a = 0.45, b = b, c = 0.8))
})
