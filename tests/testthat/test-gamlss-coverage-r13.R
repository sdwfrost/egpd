test_that("fixed-shape GPD information and predictions do not depend on the fixed link", {
  skip_if_not_installed("gamlss")
  u <- (seq_len(16) - .5)/16
  nsim <- 4000L
  nd <- data.frame(x = 42)
  for (xi in c(.3, 1.8)) {
    d <- data.frame(y = 1.4*expm1(-xi*log1p(-u))/xi, x = seq_along(u))
    for (link in c("identity", "log")) {
      ## A local family object, not a bare constructor on the search path.
      fam <- egpd::EGPD1(sigma.link = link)
      f <- gamlss::gamlss(y ~ 1, data = d, family = fam, mu.start = 1.4,
                          sigma.fix = TRUE, sigma.start = xi, nu.fix = TRUE, nu.start = 1,
                          trace = FALSE)
      eta <- log(f$mu.fv[1])
      t <- d$y/exp(eta)
      expect_true(f$converged)
      expect_equal(as.numeric(logLik(f)), -sum(eta + (1+1/xi)*log1p(xi*t)))
      se <- 1/sqrt(sum((1+xi)*t/(1+xi*t)^2))
      expect_equal(sqrt(.gamlss_vcov(f)[1]), se, tolerance = 1e-7)
      for (m in c("plug-in", "parametric")) {
        set.seed(871)
        r <- predictive_coverage(f, 1, nd, data_fit = d, method = m, nsim = nsim)
        set.seed(871)
        scale <- if (m == "parametric") exp(rnorm(nsim, eta, se)) else exp(eta)
        ref <- quantile(scale*expm1(-xi*log1p(-runif(nsim)))/xi, c(.025, .975), names = FALSE)
        expect_equal(c(r$L, r$U), ref, tolerance = 2e-7)
      }
    }
  }
})

test_that("held-out grouped probit uses trial sizes, case weights, offsets and estimable columns", {
  skip_if_not_installed("gamlss")
  set.seed(107)
  n <- 42L
  d <- data.frame(x = 1e4*seq(-.8, .8, length.out = n), t = runif(n, .8, 1.4),
                  bd = sample(2:14, n, TRUE), w = rep(c(.5, 2, 1.4), length.out = n))
  d$x2 <- 2*d$x
  d$w[4] <- 0
  d$y <- rbinom(n, d$bd, pnorm(.1 - .7*d$x/1e4 + log(d$t)))
  f <- gamlss::gamlss(cbind(y, bd-y) ~ x + x2 + offset(log(t)), weights = w, data = d,
                      family = gamlss.dist::BI(mu.link = "probit"), trace = FALSE)
  expect_true(is.na(coef(f)["x2"]))
  ## -d^2/deta^2 [y log Phi(eta) + (bd-y) log Phi(-eta)], not Fisher weights.
  eta <- f$mu.lp
  p <- pnorm(eta)
  a <- dnorm(eta)
  info <- a^2*(d$y/p^2 + (d$bd-d$y)/(1-p)^2) + eta*a*(d$y-d$bd*p)/(p*(1-p))
  X <- cbind(1, d$x/1e4)
  V <- solve(crossprod(X, X*(d$w*info)))
  Vfit <- .gamlss_vcov(f)
  expect_identical(colnames(Vfit), c("mu.(Intercept)", "mu.x"))
  expect_equal(unname(Vfit)*outer(c(1, 1e4), c(1, 1e4)), unname(V), tolerance = 1e-6)

  i <- c(5, 22, 38)
  nd <- d[i, ]
  sizes <- c(29, 0, 7)
  nsim <- 2500L
  se <- sqrt(rowSums((X[i, ] %*% V) * X[i, ]))
  set.seed(819)
  r <- predictive_coverage(f, c(18, 0, 4), nd, bd = sizes, data_fit = d, nsim = nsim)
  set.seed(819)
  prob <- pnorm(rnorm(3*nsim, rep(eta[i], each = nsim), rep(se, each = nsim)))
  yy <- qbinom(runif(3*nsim), rep(sizes, each = nsim), prob)
  ref <- apply(matrix(yy, nsim), 2, quantile, c(.025, .975), type = 1, names = FALSE)
  expect_equal(r$L, ref[1, ])
  expect_equal(r$U, ref[2, ])
  nd$x2 <- 0
  expect_error(predictive_coverage(f, c(18, 0, 4), nd, bd = sizes, data_fit = d, nsim = 5),
               "Non-estimable")
})

test_that("information expands at the accepted fitted predictor after damping", {
  skip_if_not_installed("gamlss")
  d <- data.frame(x = seq(-1, 1, length.out = 49), t = seq(1, 1.3, length.out = 49),
                  y = c(rep(c(0:2, 2, 1), 7), 1:14))
  f <- gamlss::gamlss(y ~ x + offset(log(t)), data = d, family = gamlss.dist::PO(),
                      mu.start = 2, control = gamlss::gamlss.control(mu.step = .15, c.crit = .2,
                                                                     trace = FALSE),
                      i.control = gamlss::glim.control(cyc = 2))
  ## GAMLSS 5.5-0 saves the last working coefficients, not the damped ones.
  ## Use the observed Poisson curvature at the actual returned rates.
  X <- cbind(1, d$x)
  V <- solve(crossprod(X, X*f$mu.fv))
  expect_equal(unname(.gamlss_vcov(f)), unname(V), tolerance = 1e-6)
  nd <- d[c(5, 25), ]
  z <- suppressWarnings(predict(f, what = "mu", newdata = nd, data = d, type = "link"))
  expect_equal(.gamlss_newdata_se(f, "mu", nd, .gamlss_vcov(f), z),
               sqrt(rowSums((X[c(5, 25), ] %*% V) * X[c(5, 25), ])), tolerance = 1e-6)
  ## Do not report a covariance for an incompatible reconstructed likelihood.
  f$G.deviance <- f$G.deviance + 1
  expect_error(.gamlss_vcov(f), "does not reproduce the fitted weighted log likelihood")
})

test_that("boundary-fixed Bell inflation adds no coordinate under a non-default free link", {
  skip_if_not_installed("gamlss")
  d <- data.frame(y = rep(c(0:4, 9), 9))
  f <- gamlss::gamlss(y ~ 1, family = egpd::ZIBELL(mu.link = "sqrt"), data = d,
                      sigma.fix = TRUE, sigma.start = 0, trace = FALSE)
  ## At the sample mean the Bell score vanishes and the observed information
  ## on sqrt(mean) is 4 n / (1 + W(mean)). The fixed logit(0) is not finite.
  mu <- mean(d$y)
  theta <- uniroot(function(t) t*exp(t) - mu, c(0, 5), tol = 1e-12)$root
  expect_equal(unname(f$mu.fv[1]), mu, tolerance = 1e-7)
  V <- .gamlss_vcov(f)
  expect_identical(colnames(V), "mu.(Intercept)")
  expect_equal(V[1], (1+theta)/(4*nrow(d)), tolerance = 1e-6)
})

test_that("a parameter with only an offset has no artificial uncertainty", {
  skip_if_not_installed("gamlss")
  d <- data.frame(x = seq(0, 1, length.out = 12), y = rep(1:4, 3))
  f <- gamlss::gamlss(y ~ offset(x) - 1, family = gamlss.dist::PO(), data = d, trace = FALSE)
  nd <- data.frame(x = .4)
  set.seed(87); r <- predictive_coverage(f, 2, nd, data_fit = d, nsim = 500)
  set.seed(87); q <- quantile(rpois(500, exp(.4)), c(.025, .975), type = 1, names = FALSE)
  expect_equal(c(r$L, r$U), q)
})
