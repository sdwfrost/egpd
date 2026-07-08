## Tests for the Bell distribution (Castellares, Ferrari & Lemonte 2018) and
## zero-inflated Bell (ZIBell), across three categories:
##   (1) comparison against the paper's closed-form quantities,
##   (2) mutual d/p/q consistency,
##   (3) analytic (C++) score AND observed-information Hessian vs finite diffs.

expit <- function(x) 1 / (1 + exp(-x))

# ---- (1) reference comparisons -----------------------------------------------

test_that("Bell pmf matches the paper's explicit small-y formulas", {
  ## Pr(Y=y) = theta^y e^{1-e^theta} B_y/y!, with B_0..B_5 = 1,1,2,5,15,52.
  for (theta in c(0.4, 0.7, 1.3)) {
    e <- exp(1 - exp(theta))
    ref <- c(e, theta * e, theta^2 * e, 5 * theta^3 / 6 * e,
             15 * theta^4 / 24 * e, 52 * theta^5 / 120 * e)
    expect_equal(egpd::dbell(0:5, theta = theta), ref, tolerance = 1e-12)
  }
})

test_that("Bell pmf sums to 1", {
  for (theta in c(0.3, 0.8, 1.5, 2.2))
    expect_equal(sum(egpd::dbell(0:5000, theta = theta)), 1, tolerance = 1e-9)
})

test_that("closed-form mean and variance match the recursion", {
  theta <- 1.3
  p <- egpd::dbell(0:6000, theta = theta)
  m1 <- sum((0:6000) * p); m2 <- sum((0:6000)^2 * p)
  expect_equal(m1, theta * exp(theta), tolerance = 1e-8)          # E[Y] = theta e^theta
  expect_equal(m2 - m1^2, theta * (1 + theta) * exp(theta),       # Var = theta(1+theta)e^theta
               tolerance = 1e-6)
})

test_that("mean-parameterisation dBELL has the requested mean", {
  p <- egpd::dBELL(0:30000, mu = 4)
  expect_equal(sum((0:30000) * p), 4, tolerance = 1e-4)
})

test_that("native theta MLE equals W0(ybar)", {
  set.seed(11)
  y <- egpd::rbell(4000, theta = 0.9)
  nll <- function(lt) -sum(egpd::dbell(y, theta = exp(lt), log = TRUE))
  o <- optimize(nll, c(-6, 3))
  expect_equal(exp(o$minimum), as.numeric(egpd:::bell_W0_cpp(mean(y))), tolerance = 1e-4)
})

# ---- (2) d/p/q consistency ---------------------------------------------------

test_that("Bell CDF = cumsum(pmf) and qbell inverts pbell", {
  for (theta in c(0.4, 0.9, 1.5)) {
    d <- egpd::dbell(0:2000, theta)
    expect_equal(egpd::pbell(0:2000, theta), cumsum(d), tolerance = 1e-10)
    ## round-trip over the well-resolved support (1 - F(k) > 1e-8); beyond that
    ## the light Bell tail underflows the CDF resolution.
    k <- 0:(sum(egpd::pbell(0:80, theta) < 1 - 1e-8) - 1)
    expect_equal(egpd::qbell(egpd::pbell(k, theta), theta), k)
  }
})

test_that("ZIBell pmf sums to 1, mixes correctly, and q inverts p", {
  theta <- 0.9; pv <- 0.3
  dz <- egpd::dzibell(0:2000, theta, pi = pv)
  expect_equal(sum(dz), 1, tolerance = 1e-8)
  expect_equal(dz[1], pv + (1 - pv) * egpd::dbell(0, theta), tolerance = 1e-12)
  expect_equal(egpd::pzibell(0:2000, theta, pi = pv), cumsum(dz), tolerance = 1e-10)
  k <- 0:(sum(egpd::pzibell(0:80, theta, pi = pv) < 1 - 1e-8) - 1)
  expect_equal(egpd::qzibell(egpd::pzibell(k, theta, pi = pv), theta, pi = pv), k)
})

test_that("log-density options agree with log of the density", {
  expect_equal(egpd::dbell(0:10, 0.8, log = TRUE), log(egpd::dbell(0:10, 0.8)))
  expect_equal(egpd::dzibell(0:10, 0.8, pi = 0.2, log = TRUE),
               log(egpd::dzibell(0:10, 0.8, pi = 0.2)))
})

# ---- (3) analytic C++ derivatives vs finite differences ----------------------

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
.hess_fd <- function(d0fn, e, y) {
  h <- 1e-4; n <- length(e); H <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) {
    ep <- e; ep[i] <- ep[i] + h; ep[j] <- ep[j] + h
    em <- e; em[i] <- em[i] + h; em[j] <- em[j] - h
    ea <- e; ea[i] <- ea[i] - h; ea[j] <- ea[j] + h
    eb <- e; eb[i] <- eb[i] - h; eb[j] <- eb[j] - h
    H[i, j] <- (d0fn(ep, y) - d0fn(em, y) - d0fn(ea, y) + d0fn(eb, y)) / (4 * h * h)
  }
  H
}

test_that("Bell d12 gradient and exact Hessian match FD (both parameterisations)", {
  cfgs <- list(
    list(d0 = egpd:::bell1d0,    d12 = egpd:::bell1d12,    e = log(3)),
    list(d0 = egpd:::bellnat1d0, d12 = egpd:::bellnat1d12, e = log(0.8))
  )
  for (cf in cfgs) for (y in c(0, 1, 5, 20)) {
    row <- .d0(cf$d12, cf$e, y, 1)
    g <- row[1]; H <- row[2]
    gf <- .grad_fd(function(e, yy) .d0(cf$d0, e, yy, 1), cf$e, y)
    Hf <- .hess_fd(function(e, yy) .d0(cf$d0, e, yy, 1), cf$e, y)[1, 1]
    expect_equal(as.numeric(g), as.numeric(gf), tolerance = 1e-5)
    expect_equal(as.numeric(H), as.numeric(Hf), tolerance = 1e-4)  # observed info, not BHHH
  }
})

test_that("ZIBell d12 gradient and exact Hessian match FD (both parameterisations)", {
  cfgs <- list(
    list(d0 = egpd:::zibell1d0,    d12 = egpd:::zibell1d12,    e = c(log(3),   qlogis(0.25))),
    list(d0 = egpd:::zibellnat1d0, d12 = egpd:::zibellnat1d12, e = c(log(0.8), qlogis(0.25)))
  )
  for (cf in cfgs) for (y in c(0, 1, 5)) {
    row <- .d0(cf$d12, cf$e, y, 2)
    g <- row[1:2]; H <- row[3:5]
    gf <- .grad_fd(function(e, yy) .d0(cf$d0, e, yy, 2), cf$e, y)
    Hm <- .hess_fd(function(e, yy) .d0(cf$d0, e, yy, 2), cf$e, y)
    expect_equal(as.numeric(g), as.numeric(gf), tolerance = 1e-5)
    expect_equal(as.numeric(H), c(Hm[1, 1], Hm[1, 2], Hm[2, 2]), tolerance = 1e-4)
  }
})

# ---- fitting surfaces --------------------------------------------------------

test_that("native GAM fit matches the direct optim MLE (both parameterisations)", {
  skip_on_cran()
  set.seed(7)
  y <- egpd:::rbell(2000, theta = 0.8)
  dat <- data.frame(y = y)
  fm <- egpd(y ~ 1, data = dat, family = "bell", trace = 0)
  fn <- egpd(y ~ 1, data = dat, family = "bellnat", trace = 0)
  ## the two parameterisations describe the same fitted distribution
  expect_equal(as.numeric(logLik(fm)), as.numeric(logLik(fn)), tolerance = 1e-3)
  ## and they reach the direct MLE (up to the REML/Laplace constant)
  nll <- function(lt) -sum(egpd:::dbell(y, theta = exp(lt), log = TRUE))
  o <- optimize(nll, c(-5, 3))
  k <- 0.5 * 1 * log(2 * pi)   # Laplace constant for 1 parameter
  expect_equal(as.numeric(logLik(fn)) + k, -o$objective, tolerance = 1e-2)
})

test_that("fitegpd recovers Bell / ZIBell parameters", {
  skip_on_cran()
  set.seed(3)
  y <- egpd:::rbell(3000, theta = 0.9)
  f <- fitegpd(y, family = "bell")
  expect_equal(unname(f$estimate["theta"]), 0.9, tolerance = 0.08)

  set.seed(5)
  yz <- egpd:::rzibell(3000, theta = 0.9, pi = 0.3)
  fz <- fitegpd(yz, family = "zibell")
  expect_equal(unname(fz$estimate["theta"]), 0.9, tolerance = 0.12)
  expect_equal(unname(fz$estimate["pi"]),    0.3, tolerance = 0.08)
})

test_that("gamlss and bamlss Bell family constructors are well-formed", {
  fam <- egpd::BELL()
  expect_s3_class(fam, "gamlss.family")
  expect_equal(fam$nopar, 1L)
  zfam <- egpd::ZIBELL()
  expect_equal(zfam$nopar, 2L)

  gb <- egpd::bell_bamlss()
  expect_s3_class(gb, "family.bamlss")
  expect_equal(gb$names, "mu")
  ## d() delegates to the mean-parameterisation density (theta = W0(mu))
  expect_equal(gb$d(3, list(mu = 4)),
               egpd::dbell(3, theta = as.numeric(egpd:::bell_W0_cpp(4))))

  zb <- egpd::zibell_bamlss()
  expect_equal(zb$names, c("mu", "pi"))
})
