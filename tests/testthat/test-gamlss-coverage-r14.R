test_that("held-out exponential and Gamma mean curvature follows response units", {
  skip_if_not_installed("gamlss")
  n <- 20L
  u <- -log1p(-(seq_len(n) - .5)/n)
  nd <- data.frame(x = 0)
  nsim <- 2000L
  for (m in c(2.2, .0022, .0011, 2.2e5)) {
    d <- data.frame(y = m*u/mean(u), x = seq_len(n))
    f <- gamlss::gamlss(y ~ 1, data = d,
                        family = egpd::EGPD1(mu.link = "identity", sigma.link = "identity"),
                        mu.start = m, sigma.start = 0, sigma.fix = TRUE,
                        nu.start = 1, nu.fix = TRUE, trace = FALSE)
    ## xi = 0, kappa = 1: the exponential sample mean has information n/m^2.
    expect_equal(sqrt(.gamlss_vcov(f)[1])/m, 1/sqrt(n), tolerance = 2e-7)
    set.seed(84)
    r <- predictive_coverage(f, m, nd, data_fit = d, method = "parametric", nsim = nsim)
    set.seed(84)
    mu <- rnorm(nsim, 1, 1/sqrt(n))
    expect_true(all(mu > 0))
    ref <- quantile(-mu * log1p(-runif(nsim)), c(.025, .975), names = FALSE)
    expect_equal(c(r$L, r$U)/m, ref, tolerance = 2e-7)
  }
  d$y <- .0011*u/mean(u)
  g <- gamlss::gamlss(y ~ 1, data = d, family = gamlss.dist::GA(mu.link = "identity"),
                      mu.start = mean(d$y), sigma.start = .5, sigma.fix = TRUE, trace = FALSE)
  expect_equal(sqrt(.gamlss_vcov(g)[1])/mean(d$y), .5/sqrt(n), tolerance = 2e-7)
})

test_that("rescaled free location and scale use joint weighted observed information", {
  skip_if_not_installed("gamlss")
  set.seed(547)
  n <- 96L
  x <- seq(-.9, .9, length.out = n)
  z <- runif(n, -1, 1)
  o <- sin(seq_len(n))/4
  w <- rep(c(.5, 2, 1.4), length.out = n)
  w[7] <- 0
  y <- 1.4 - .7*x + o + rnorm(n, sd = .75 + .2*z)
  X <- cbind(1, x)
  Z <- cbind(1, z)
  for (s in c(.001, 1e3)) {
    d <- data.frame(y = s*y, o = s*o, x = x/s, z = z*s, w = w)
    f <- gamlss::gamlss(y ~ x + offset(o), sigma.formula = ~ z, weights = w, data = d,
                        family = gamlss.dist::NO(sigma.link = "identity"),
                        mu.start = s*1.4, sigma.start = s*.75, trace = FALSE)
    ## Differentiate sum w*(log(sig) + (y-mu)^2/(2*sig^2)) in original units.
    r <- y - f$mu.fv/s
    sig <- f$sigma.fv/s
    H11 <- crossprod(X, X*(w/sig^2))
    H12 <- crossprod(X, Z*(2*w*r/sig^3))
    H22 <- crossprod(Z, Z*(w*(3*r^2/sig^4 - 1/sig^2)))
    V <- solve(rbind(cbind(H11, H12), cbind(t(H12), H22)))
    units <- c(s, s^2, s, 1)
    expect_equal(unname(.gamlss_vcov(f))/outer(units, units), unname(V), tolerance = 1e-6)
  }
})

test_that("both held-out methods reject non-estimable Poisson rows in either term order", {
  skip_if_not_installed("gamlss")
  x <- rep(0:1, each = 60)
  d <- data.frame(x = x, z = x, y = 10*x + rep(0:4, 24))
  nd <- data.frame(x = c(0, .5, 1), z = c(0, .5, 1))
  bad <- rbind(nd, c(0, 1))
  nsim <- 500L
  for (form in list(y ~ x + z, y ~ z + x)) {
    f <- gamlss::gamlss(form, data = d, family = gamlss.dist::PO(), trace = FALSE)
    for (m in c("parametric", "plug-in"))
      expect_error(predictive_coverage(f, rep(8, 4), bad, data_fit = d, method = m, nsim = 5),
                   "Non-estimable 'mu'")
    set.seed(91)
    a <- predictive_coverage(f, rep(8, 3), nd, data_fit = d, method = "plug-in", nsim = nsim)
    set.seed(91)
    q <- apply(matrix(rpois(3*nsim, rep(2*6^nd$x, each = nsim)), nsim), 2,
               quantile, c(.025, .975), type = 1, names = FALSE)
    expect_identical(a$L, q[1, ])
    expect_identical(a$U, q[2, ])
  }
})

test_that("plug-in checks aliased dispersion without needing curvature of a smooth mean", {
  skip_if_not_installed("gamlss")
  ## gamlss resolves its smoother by name.
  withr::local_package("gamlss")
  set.seed(43)
  n <- 84L
  x <- rep(0:1, each = n/2)
  d <- data.frame(x = x, z = x, t = seq(0, 1, length.out = n))
  d$y <- sin(4*d$t) + rnorm(n, sd = .5 + .3*x)
  f <- gamlss::gamlss(y ~ pb(t), sigma.formula = ~ x + z, data = d,
                      family = gamlss.dist::NO(), trace = FALSE)
  nd <- data.frame(t = c(.1, .4), x = c(0, .5), z = c(0, .5))
  plug <- function(nd) predictive_coverage(f, c(0, 1), nd, data_fit = d, method = "plug-in",
                                          nsim = 50)
  expect_true(all(is.finite(plug(nd)$U)))
  nd$z[2] <- 1
  expect_error(plug(nd), "Non-estimable 'sigma'")
})

test_that("a non-positive information matrix is not repaired or required for plug-in", {
  skip_if_not_installed("gamlss")
  ## Symmetric Cauchy responses: at location zero the score vanishes but
  ## the observed location information 2*(1-y^2)/(1+y^2)^2 is negative.
  d <- expand.grid(y = c(-2, 2), x = rep(0:1, each = 10))
  d$z <- 1e-9*d$x
  f <- gamlss::gamlss(y ~ x + z, data = d, family = gamlss.dist::TF(), mu.start = 0,
                      sigma.start = 1, sigma.fix = TRUE, nu.start = 1, nu.fix = TRUE,
                      trace = FALSE)
  nd <- data.frame(x = .5, z = .5e-9)
  cover <- function(m) predictive_coverage(f, 0, nd, data_fit = d, nsim = 50, method = m)
  expect_error(cover("parametric"), "observed information")
  expect_true(is.finite(cover("plug-in")$U))
  nd$z <- 0
  expect_error(cover("plug-in"), "Non-estimable 'mu'")
})
