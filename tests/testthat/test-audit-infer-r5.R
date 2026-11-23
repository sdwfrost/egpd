## Round-5 regressions for inference outputs: decreasing custom inverse links in
## predict()/confint() (VC5-01) and gamlss fits with a fixed parameter in
## predictive_coverage() (VC5-07, VB5-11).

## ---- VC5-01: decreasing custom inverse links --------------------------------

## Poisson whose mean is exp(-eta), a decreasing inverse link
.custom_pois_dec <- function(unlink) list(
  d0 = function(pars, likdata) {
    e <- drop(likdata$X[[1]] %*% pars)
    y <- likdata$y[, 1]
    sum(y * e + exp(-e) + lgamma(y + 1))
  },
  d120 = function(pars, likdata) {
    e <- drop(likdata$X[[1]] %*% pars)
    y <- likdata$y[, 1]
    cbind(y - exp(-e), exp(-e))
  },
  d340 = NULL, unlink = list(unlink), q = function(p, mu) qpois(p, mu))

test_that("VC5-01: a decreasing custom inverse link gives non-negative SEs and ordered intervals", {
  set.seed(11)
  y <- rpois(60, 2.3)
  d <- data.frame(y = y)
  n <- length(y)
  S <- sum(y)
  ybar <- mean(y)
  z <- qnorm(0.975)
  ## closed forms: mu_hat = ybar and SE(eta_hat) = 1/sqrt(S), so the delta-method SE
  ## of mu is |d mu/d eta| SE = ybar/sqrt(S), and the Wald interval is the image of
  ## eta_hat -/+ z SE under exp(-eta): ybar exp(-/+ z/sqrt(S)), lower limit first
  wald <- ybar * exp(c(-1, 1) * z / sqrt(S))
  ## profile interval: the Poisson likelihood-ratio roots
  ll <- function(mu) S * log(mu) - n * mu
  lr <- function(mu) 2 * (ll(ybar) - ll(mu)) - qchisq(0.95, 1)
  prof <- c(uniroot(lr, c(1e-3, ybar), tol = 1e-12)$root,
            uniroot(lr, c(ybar, 10 * ybar), tol = 1e-12)$root)
  ul <- function(e) exp(-e)
  ul_d <- ul
  attr(ul_d, "deriv") <- function(e) -exp(-e)
  for (u in list(ul, ul_d)) {
    f <- egpd(list(mu = y ~ 1), d, family = "custom", custom.fns = .custom_pois_dec(u),
              inits = -log(ybar))
    pr <- predict(f, type = "response", se.fit = TRUE)
    expect_equal(pr$fitted$mu[1], ybar, tolerance = 1e-8)
    expect_equal(pr$se.fit$mu, rep(ybar / sqrt(S), n), tolerance = 1e-7)
    expect_equal(unname(confint(f)["mu", ]), wald, tolerance = 1e-7)
    expect_equal(unname(confint(f, newdata = d[1:2, , drop = FALSE])),
                 rbind(wald, wald, deparse.level = 0), tolerance = 1e-7)
    expect_equal(unname(confint(f, method = "profile")["mu", ]), prof, tolerance = 1e-4)
    expect_equal(unname(confint(f, method = "profile", level = 1)["mu", ]), c(0, Inf))
  }
})

test_that("VC5-01: the inverse link mean = 1/eta agrees with glm(Gamma('inverse'))", {
  set.seed(12)
  x <- rexp(40, rate = 1 / 3)
  d <- data.frame(y = x)
  fns <- list(
    d0 = function(pars, likdata) {
      e <- drop(likdata$X[[1]] %*% pars)
      if (any(e <= 0)) return(1e20)
      sum(likdata$y[, 1] * e - log(e))
    },
    d120 = function(pars, likdata) {
      e <- drop(likdata$X[[1]] %*% pars)
      cbind(likdata$y[, 1] - 1 / e, 1 / e^2)
    },
    d340 = NULL, unlink = list(function(e) 1 / e), q = function(p, mu) qexp(p, 1 / mu))
  f <- egpd(list(mu = y ~ 1), d, family = "custom", custom.fns = fns, inits = 1 / mean(x))
  ## the exponential law is the gamma GLM at dispersion 1, whose canonical link is 1/mu
  ## (stats::Gamma: bamlss attaches distributions3, which masks it)
  g <- stats::glm(y ~ 1, family = stats::Gamma("inverse"), data = d)
  eta <- unname(coef(g))
  se <- sqrt(vcov(g, dispersion = 1)[1, 1])
  z <- qnorm(0.975)
  ci <- 1 / (eta + c(1, -1) * z * se)
  pr <- predict(f, type = "response", se.fit = TRUE)
  expect_equal(pr$fitted$mu[1], 1 / eta, tolerance = 1e-6)
  expect_equal(pr$se.fit$mu[1], se / eta^2, tolerance = 1e-5)
  expect_equal(unname(confint(f)["mu", ]), ci, tolerance = 1e-5)
  ## 1/eta is only monotone on its valid half-line: +/-Inf both map to 0.
  ## Without that custom domain, do not report [0,0].
  expect_error(confint(f, method = "profile", level = 1), "custom inverse link")
  expect_error(confint(f, newdata = d[1, , drop = FALSE], level = 1), "custom inverse link")
  expect_equal(unname(confint(f, newdata = d[1, , drop = FALSE])[1, ]), ci, tolerance = 1e-5)
  ## and the closed forms: SE xbar/sqrt(n), interval xbar/(1 +/- z/sqrt(n))
  expect_equal(pr$se.fit$mu[1], mean(x) / sqrt(40), tolerance = 1e-7)
  expect_equal(unname(confint(f)["mu", ]), mean(x) / (1 + c(1, -1) * z / sqrt(40)),
               tolerance = 1e-7)
})

test_that("VC5-01: an unbracketed profile limit stays NA on its own side", {
  ## Negative log-likelihood eta^2/2 for eta <= 0 and c log(1 + eta^2/(2c)) above, c =
  ## 0.15: the estimate is 0 with SE 1; the profile deviance reaches qchisq(0.95, 1) at
  ## eta = -1.959964, but above 0 it is still below the cutoff and rising 96 SEs out
  ## (2.68 at 48 SEs, 3.10 at 96), so that link-scale limit is NA. Under mean = exp(-eta)
  ## it is the lower limit of the mean; under exp(eta), the upper.
  cc <- 0.15
  f0 <- function(e) ifelse(e <= 0, e^2 / 2, cc * log1p(e^2 / (2 * cc)))
  f1 <- function(e) ifelse(e <= 0, e, e / (1 + e^2 / (2 * cc)))
  f2 <- function(e) ifelse(e <= 0, 1, (1 - e^2 / (2 * cc)) / (1 + e^2 / (2 * cc))^2)
  fns <- function(unlink) list(
    d0 = function(pars, likdata) f0(drop(likdata$X[[1]] %*% pars))[1],
    d120 = function(pars, likdata) {
      e <- drop(likdata$X[[1]] %*% pars)
      cbind(f1(e), f2(e)) / length(e)
    },
    d340 = NULL, unlink = list(unlink), q = function(p, mu) qexp(p, 1 / mu))
  d <- data.frame(y = rep(1, 10))
  root <- sqrt(qchisq(0.95, 1))
  f <- egpd(list(mu = y ~ 1), d, family = "custom", custom.fns = fns(function(e) exp(-e)),
            inits = 0.3)
  expect_warning(ci <- confint(f, method = "profile"), "Could not bracket")
  expect_true(is.na(ci["mu", "lower"]))
  expect_equal(ci["mu", "upper"], exp(root), tolerance = 1e-4)
  f <- egpd(list(mu = y ~ 1), d, family = "custom", custom.fns = fns(function(e) exp(e)),
            inits = 0.3)
  expect_warning(ci <- confint(f, method = "profile"), "Could not bracket")
  expect_equal(ci["mu", "lower"], exp(-root), tolerance = 1e-4)
  expect_true(is.na(ci["mu", "upper"]))
})

## ---- VC5-07: gamlss fits with a fixed parameter ------------------------------

## Equal-tailed discrete intervals of an nsim x n matrix of draws (type-1 quantiles)
.eq_tail <- function(ys, level = 0.95) {
  a <- (1 - level) / 2
  list(L = apply(ys, 2, quantile, probs = a, type = 1, names = FALSE),
       U = apply(ys, 2, quantile, probs = 1 - a, type = 1, names = FALSE))
}

test_that("VC5-07: predictive_coverage() holds a gamlss parameter fixed in the fit", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("numDeriv")
  set.seed(21)
  d <- data.frame(x = runif(120))
  d$y <- rnbinom(120, mu = exp(0.5 + d$x), size = 2)
  g <- gamlss::gamlss(y ~ x, sigma.fix = TRUE, sigma.start = 0.5,
                      family = gamlss.dist::NBI(), data = d, trace = FALSE)
  nd <- d[1:10, , drop = FALSE]
  nsim <- 400
  ## The references simulate the NBI law with sigma = 0.5 directly. Plug-in: mu from
  ## gamlss's own predictions, in sample and at newdata.
  mu_in <- as.numeric(fitted(g, "mu"))
  mu_nd <- as.numeric(predict(g, what = "mu", newdata = nd, type = "response", data = d))
  set.seed(5)
  r <- predictive_coverage(g, d$y, nsim = nsim, method = "plug-in")
  set.seed(5)
  ref <- .eq_tail(sapply(mu_in, function(m) gamlss.dist::rNBI(nsim, mu = m, sigma = 0.5)))
  expect_identical(r[c("L", "U")], ref)
  set.seed(6)
  r <- predictive_coverage(g, nd$y, newdata = nd, data_fit = d, nsim = nsim, method = "plug-in")
  set.seed(6)
  ref <- .eq_tail(sapply(mu_nd, function(m) gamlss.dist::rNBI(nsim, mu = m, sigma = 0.5)))
  expect_identical(r[c("L", "U")], ref)
  ## Parametric: log mu ~ N(eta, se) and sigma = 0.5 in every draw. In sample, se is
  ## gamlss's predict(se.fit = TRUE); at newdata, from the observed information of the
  ## NBI log-likelihood in mu's coefficients with sigma held at 0.5.
  pr <- predict(g, what = "mu", type = "link", se.fit = TRUE)
  X <- cbind(1, d$x)
  ll <- function(b) sum(gamlss.dist::dNBI(d$y, mu = exp(drop(X %*% b)), sigma = 0.5, log = TRUE))
  V <- solve(-numDeriv::hessian(ll, coef(g, "mu")))
  Xn <- X[1:10, ]
  draw <- function(m, s) {
    eta <- matrix(rnorm(nsim * length(m), rep(m, each = nsim), rep(s, each = nsim)), nsim)
    .eq_tail(sapply(seq_along(m), function(j) gamlss.dist::rNBI(nsim, exp(eta[, j]), 0.5)))
  }
  set.seed(7)
  r <- predictive_coverage(g, d$y, nsim = nsim, method = "parametric")
  set.seed(7)
  expect_identical(r[c("L", "U")], draw(as.numeric(pr$fit), as.numeric(pr$se.fit)))
  set.seed(8)
  r <- predictive_coverage(g, nd$y, newdata = nd, data_fit = d, nsim = nsim,
                           method = "parametric")
  set.seed(8)
  expect_identical(r[c("L", "U")], draw(drop(Xn %*% coef(g, "mu")),
                                        sqrt(rowSums((Xn %*% V) * Xn))))

  ## fixed at observation-specific values: used in sample, unknown at newdata
  g2 <- gamlss::gamlss(y ~ x, sigma.fix = TRUE, sigma.start = rep(c(0.3, 0.6), 60),
                       family = gamlss.dist::NBI(), data = d, trace = FALSE)
  set.seed(9)
  r <- predictive_coverage(g2, d$y, nsim = nsim, method = "plug-in")
  set.seed(9)
  ref <- .eq_tail(sapply(seq_len(120), function(j)
    gamlss.dist::rNBI(nsim, mu = fitted(g2, "mu")[j], sigma = c(0.3, 0.6)[2 - j %% 2])))
  expect_identical(r[c("L", "U")], ref)
  expect_error(predictive_coverage(g2, nd$y, newdata = nd, data_fit = d, nsim = 50,
                                   method = "plug-in"),
               "vary between observations")
})

test_that("VC5-07: gamlss new-data standard errors are conditional on a fixed parameter", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("numDeriv")
  ## DEGPD1 with the tail index (gamlss's sigma) fixed at 0.2. vcov() includes it as if
  ## it were estimated, which inflates the scale's standard errors (strongly correlated
  ## with the tail index); the reference is the inverse observed information of an
  ## independent DEGPD likelihood in the free coefficients only, P(Y = y) = H(y + 1) -
  ## H(y) with H(x) = (1 - (1 + xi x/s)^(-1/xi))^kappa.
  set.seed(7)
  d <- data.frame(x = runif(250))
  d$y <- rdiscegpd(250, sigma = exp(1 + d$x), xi = 0.2, kappa = 1.5, type = 1)
  g <- suppressWarnings(gamlss::gamlss(
    y ~ x, sigma.fix = TRUE, sigma.start = 0.2, nu.formula = ~ 1, data = d,
    family = DEGPD1(), trace = FALSE,
    control = gamlss::gamlss.control(n.cyc = 200, trace = FALSE)))
  H <- function(x, s, k) (1 - (1 + 0.2 * x / s)^(-1 / 0.2))^k
  X <- cbind(1, d$x)
  ll <- function(th) {
    s <- exp(drop(X %*% th[1:2]))
    sum(log(H(d$y + 1, s, exp(th[3])) - H(d$y, s, exp(th[3]))))
  }
  V <- solve(-numDeriv::hessian(ll, c(coef(g, "mu"), coef(g, "nu"))))
  nd <- data.frame(x = c(0.1, 0.5, 0.9))
  Xn <- cbind(1, nd$x)
  se_ref <- sqrt(rowSums((Xn %*% V[1:2, 1:2]) * Xn))
  Vg <- .gamlss_vcov(g)
  se <- .gamlss_newdata_se(g, "mu", nd, Vg, drop(Xn %*% coef(g, "mu")))
  expect_equal(se, se_ref, tolerance = 1e-5)
  ## All free coordinates, including their cross-parameter covariances.
  expect_equal(unname(Vg), unname(V), tolerance = 1e-5)
  expect_identical(colnames(Vg), c("mu.(Intercept)", "mu.x", "nu.(Intercept)"))
  expect_equal(.gamlss_newdata_se(g, "nu", nd, Vg, rep(coef(g, "nu"), 3)),
               rep(sqrt(V[3, 3]), 3), tolerance = 1e-5)
  set.seed(3)
  r <- predictive_coverage(g, c(2, 5, 9), newdata = nd, data_fit = d, nsim = 200,
                           method = "parametric")
  expect_length(r$U, 3)
})
