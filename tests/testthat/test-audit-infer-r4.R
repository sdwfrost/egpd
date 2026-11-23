## Round-4 regressions for inference outputs: predictive coverage, confidence
## intervals, simulate/predict naming and the custom-family inverse links.

## ---- VC-01: gamlss parametric coverage at newdata ----------------------------

test_that("VC-01: gamlss parametric coverage at newdata keeps parameter uncertainty", {
  skip_if_not_installed("gamlss")
  ## Intercept-only Poisson: SE(log mu) = 1/sqrt(sum(y)) = 1/sqrt(2), so the
  ## parametric predictive law at any row is int ppois(k, mu_hat e^{z/sqrt(2)}) dPhi(z):
  ## F(2) = 0.95178 < 0.975 <= F(3) = 0.98340, so its 97.5% quantile is 3, while the
  ## plug-in law ppois(k, 0.5) gives 2.
  y <- c(0, 0, 1, 1)
  d <- data.frame(y = y)
  Fpar <- function(k) integrate(function(z)
    ppois(k, mean(y) * exp(z / sqrt(sum(y)))) * dnorm(z), -12, 12, rel.tol = 1e-12)$value
  q_ref <- min(which(vapply(0:10, Fpar, numeric(1)) >= 0.975)) - 1
  expect_identical(q_ref, 3)
  expect_identical(qpois(0.975, mean(y)), 2)

  g <- gamlss::gamlss(y ~ 1, family = gamlss.dist::PO(), data = d, trace = FALSE)
  set.seed(11)
  r <- predictive_coverage(g, 0, newdata = data.frame(y = 0), data_fit = d,
                           nsim = 4e4, method = "parametric")
  expect_identical(r$method, "parametric")
  expect_identical(r$U, q_ref)
})

test_that("VC-01: gamlss new-data standard errors match the Poisson information", {
  skip_if_not_installed("gamlss")
  ## Poisson log link: the information is X' diag(mu) X, so the SE of the linear
  ## predictor at x0 is sqrt(x0' (X' W X)^{-1} x0); with a factor and an offset.
  set.seed(3)
  n <- 300
  dd <- data.frame(x = runif(n), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                   t = runif(n, 1, 2))
  dd$y <- rpois(n, dd$t * exp(0.5 + 1.2 * dd$x + 0.3 * (dd$f == "b")))
  g <- gamlss::gamlss(y ~ x + f + offset(log(t)), family = gamlss.dist::PO(),
                      data = dd, trace = FALSE)
  X <- cbind(1, dd$x, dd$f == "b", dd$f == "c")
  mu <- as.numeric(fitted(g, "mu"))
  Vref <- solve(crossprod(X * mu, X))
  nd <- data.frame(x = c(0.1, 0.9, 0.5), f = factor(c("c", "a", "b"), levels = c("a", "b", "c")),
                   t = c(1.2, 1.9, 1))
  X0 <- cbind(1, nd$x, nd$f == "b", nd$f == "c")
  se_ref <- sqrt(rowSums((X0 %*% Vref) * X0))
  eta <- drop(X0 %*% coef(g, "mu")) + log(nd$t)
  se <- .gamlss_newdata_se(g, "mu", nd, .gamlss_vcov(g), eta)
  expect_equal(se, se_ref, tolerance = 1e-4)
  ## a design that does not reproduce gamlss's prediction is an error, not a guess
  expect_error(.gamlss_newdata_se(g, "mu", nd, .gamlss_vcov(g), eta + 1), "does not reproduce")
})

test_that("VC-01: gamlss smoothers at newdata are an error, not a plug-in interval", {
  skip_if_not_installed("gamlss")
  skip_on_cran()
  ## gamlss resolves smoothers such as pb() on the search path
  if (!"package:gamlss" %in% search()) {
    suppressPackageStartupMessages(library(gamlss))
    on.exit(detach("package:gamlss", unload = FALSE, character.only = TRUE), add = TRUE)
  }
  set.seed(4)
  dd <- data.frame(x = runif(150))
  dd$y <- rpois(150, exp(0.5 + dd$x))
  gs <- gamlss::gamlss(y ~ pb(x), family = gamlss.dist::PO(), data = dd, trace = FALSE)
  nd <- data.frame(x = c(0.2, 0.8))
  expect_error(predictive_coverage(gs, c(1, 2), newdata = nd, data_fit = dd, nsim = 50,
                                   method = "parametric"),
               "additive smoothing terms")
  r <- predictive_coverage(gs, c(1, 2), newdata = nd, data_fit = dd, nsim = 50,
                           method = "plug-in")
  expect_identical(r$method, "plug-in")
})

## ---- VC-02: custom-family inverse links --------------------------------------

## Poisson with a log link, written as an evgam-style custom family; the predictor is
## named `mu` (no link prefix), so only `unlink` says that it is log-linked.
.custom_pois <- function(unlink = list(exp)) {
  out <- list(
    d0 = function(pars, likdata) {
      eta <- as.vector(likdata$X[[1]] %*% pars)
      y <- likdata$y[, 1]
      sum(exp(eta) - y * eta + lgamma(y + 1))
    },
    d120 = function(pars, likdata) {
      eta <- as.vector(likdata$X[[1]] %*% pars)
      y <- likdata$y[, 1]
      cbind(exp(eta) - y, exp(eta))
    },
    d340 = NULL,
    q = function(p, mu) qpois(p, mu))
  out$unlink <- unlink
  out
}

test_that("VC-02: simulate(), confint() and predict() apply custom.fns$unlink", {
  set.seed(4)
  y <- rpois(150, 0.6)
  d <- data.frame(y = y)
  f <- egpd(list(mu = y ~ 1), data = d, family = "custom", custom.fns = .custom_pois())
  ## closed forms: mu_hat = ybar, SE(log mu_hat) = 1 / sqrt(sum(y))
  mu_hat <- mean(y)
  se_eta <- 1 / sqrt(sum(y))
  wald <- exp(log(mu_hat) + c(-1, 1) * qnorm(0.975) * se_eta)

  lk <- simulate(f, nsim = 5, seed = 9, type = "link")[["mu"]]
  rs <- simulate(f, nsim = 5, seed = 9, type = "response")[["mu"]]
  expect_equal(rs, exp(lk), tolerance = 1e-14)
  qs <- simulate(f, nsim = 5, seed = 9, probs = 0.9)[["q:0.9"]]
  expect_identical(c(qs), qpois(0.9, exp(c(lk))))

  expect_equal(unname(confint(f)["mu", ]), wald, tolerance = 1e-6)
  expect_equal(unname(confint(f, newdata = d[1:2, , drop = FALSE])),
               unname(rbind(wald, wald)), tolerance = 1e-6)

  ## without a "deriv" attribute, predict(se.fit = TRUE) still gives the
  ## delta-method SE on the response scale, mu_hat * SE(log mu_hat)
  pr <- predict(f, type = "response", se.fit = TRUE)
  expect_equal(pr$fitted$mu[1], mu_hat, tolerance = 1e-7)
  expect_equal(pr$se.fit$mu[1], mu_hat * se_eta, tolerance = 1e-7)

  ## no unlink: the link named by the predictor prefix, the same in every method
  f2 <- egpd(list(logmu = y ~ 1), data = d, family = "custom",
             custom.fns = .custom_pois(unlink = NULL))
  expect_equal(predict(f2, type = "response")$mu[1], mu_hat, tolerance = 1e-7)
  expect_equal(unname(confint(f2)["mu", ]), wald, tolerance = 1e-6)
})

## ---- VC-03: profile intervals at gamma != 1 ----------------------------------

test_that("VC-03: profile intervals use the unscaled likelihood whatever gamma is", {
  ## Bell(theta): l(theta) = S log(theta) - n e^theta + const, maximised where
  ## theta e^theta = ybar; the 95% profile interval solves 2 (l(theta_hat) - l(theta))
  ## = qchisq(0.95, 1) -- the raw likelihood, not gamma times it. The data come from
  ## qbell(runif()) rather than rbell(), whose exact sampler draws differently.
  set.seed(2)
  y <- qbell(runif(80), theta = 0.9)
  d <- data.frame(y = y)
  n <- length(y)
  S <- sum(y)
  th <- uniroot(function(t) t * exp(t) - mean(y), c(1e-6, 10), tol = 1e-14)$root
  ll <- function(t) S * log(t) - n * exp(t)
  lr <- function(t) 2 * (ll(th) - ll(t)) - qchisq(0.95, 1)
  ref <- c(uniroot(lr, c(1e-8, th), tol = 1e-14)$root,
           uniroot(lr, c(th, 50), tol = 1e-14)$root)
  expect_equal(ref, c(0.81825244, 1.00979889), tolerance = 1e-8)
  for (g in c(1, 4)) {
    f <- egpd(y ~ 1, d, family = "bellnat", gamma = g)
    expect_equal(unname(confint(f, method = "profile")[1, ]), ref, tolerance = 1e-4,
                 label = paste("profile interval at gamma =", g))
  }
})

## ---- VC-13: newdata intervals for unidentified parameters --------------------

test_that("VC-13: confint(newdata) gives NA, not a point, for zero-variance parameters", {
  ## intercept-only, so the parameter at any row is the global one: a coefficient
  ## without positive variance has no interval in either form
  set.seed(12)
  d <- data.frame(y = rdiscegpd(300, sigma = 3, xi = 0.2, kappa = 1.2, type = 1))
  f <- egpd(list(y ~ 1, ~ 1, ~ 1), data = d, family = "degpd", restarts = FALSE)
  f$Vc[2, ] <- f$Vc[, 2] <- 0
  f[[2]]$Vc[] <- 0
  ci0 <- confint(f)
  cin <- confint(f, newdata = data.frame(row = 1))
  expect_true(all(is.na(ci0["shape", ])))
  expect_true(all(is.na(cin["shape[1]", ])))
  expect_equal(unname(cin[c("scale[1]", "kappa[1]"), ]), unname(ci0[c("scale", "kappa"), ]))

  ## the reviewers' case: a CompPareto fit whose alpha and theta are unidentified. It is
  ## the default start's fit, which has every observation below its splice point; the
  ## default restarts refit it (to the maximum, with both identified), so it is
  ## requested with restarts = FALSE.
  set.seed(2)
  yc <- rcomppareto(120, spec = "lnorm", meanlog = -0.5, sdlog = 0.8, alpha = 1.5, theta = 2)
  dfc <- data.frame(y = yc)
  fit <- suppressWarnings(egpd(list(y ~ 1, ~ 1, ~ 1, ~ 1), data = dfc, family = "comppareto",
                               comppareto.args = list(spec = "lnorm"), restarts = FALSE))
  if (all(diag(fit$Vc)[3:4] <= 0)) {
    ci0 <- confint(fit)
    cin <- confint(fit, newdata = dfc[1, , drop = FALSE])
    expect_true(all(is.na(ci0[c("alpha", "theta"), ])))
    expect_true(all(is.na(cin[c("alpha[1]", "theta[1]"), ])))
    expect_equal(unname(cin[c("meanlog[1]", "sdlog[1]"), ]),
                 unname(ci0[c("meanlog", "sdlog"), ]))
  }

  ## rows that load on no coefficient keep their exact interval: a pinned parameter,
  ## and a row at which a no-intercept predictor is exactly 0
  set.seed(8)
  dd <- data.frame(x = runif(300))
  dd$y <- rdiscegpd(300, sigma = exp(2 * dd$x), xi = 0.15, kappa = 1.5, type = 1)
  f <- egpd(list(y ~ x - 1, ~ 1, ~ 1), data = dd, family = "degpd",
            fixed = list(kappa = 1.5), restarts = FALSE)
  ci <- confint(f, newdata = data.frame(x = c(0, 0.5)))
  expect_equal(unname(ci["kappa[1]", ]), c(1.5, 1.5))
  expect_equal(unname(ci["scale[1]", ]), c(1, 1))
  expect_true(all(is.finite(ci["scale[2]", ])) && ci["scale[2]", 1] < ci["scale[2]", 2])
})

## ---- VC-14: plot() without smooth terms --------------------------------------

test_that("VC-14: plot() of a fit without smooths is a clear error", {
  set.seed(1)
  f <- egpd(y ~ 1, data.frame(y = rbell(50, theta = 1)), family = "bellnat")
  expect_error(plot(f), "No smooth terms")
})

## ---- VC-18: unique quantile names --------------------------------------------

test_that("VC-18: quantile columns are named by the full probability", {
  set.seed(12)
  d <- data.frame(y = rdiscegpd(300, sigma = 3, xi = 0.2, kappa = 1.2, type = 1))
  f <- egpd(list(y ~ 1, ~ 1, ~ 1), data = d, family = "degpd", restarts = FALSE)
  p <- c(0.5, 0.95, 0.999, 0.9995, 0.9999)
  want <- c("q:0.5", "q:0.95", "q:0.999", "q:0.9995", "q:0.9999")
  nd <- data.frame(row = 1)
  pr <- predict(f, newdata = nd, prob = p)
  expect_named(pr, want)
  pe <- predict(f, newdata = nd, type = "response")
  expect_identical(pr[["q:0.9999"]],
                   qdiscegpd(0.9999, sigma = pe$scale, xi = pe$shape, kappa = pe$kappa, type = 1))
  expect_named(simulate(f, newdata = nd, nsim = 2, seed = 1, probs = p), want)
})

## ---- VC-08 (follow-on): generalized Waring rho stays above 1 -----------------

test_that("VC-08: confint() for a generalized Waring rho stays inside rho > 1", {
  set.seed(2)
  x <- rgwaring(150, mu = 4, k = 3, rho = 1.6)
  f <- suppressWarnings(egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = x), family = "gwaring"))
  skip_if_not(f$convergence == 0)
  z <- qnorm(0.975)
  eta <- f$coefficients[[3]]
  se <- sqrt(f$Vc[3, 3])
  ## the log-link Wald interval reaches below rho = 1, outside rho > 1; the fitted
  ## likelihood evaluates rho as at least 1 + 1e-6 (see ?dgwaring), where it ends
  expect_lt(exp(eta - z * se), 1)
  want <- c(1 + 1e-6, exp(eta + z * se))
  expect_equal(unname(confint(f)["rho", ]), want)
  expect_equal(unname(confint(f, newdata = data.frame(row = 1))["rho[1]", ]), want)
  expect_gt(confint(f)["rho", "lower"], 1)
})
