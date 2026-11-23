## Round-3 audit regression tests for the core fitting machinery (core-fit area):
## optimiser failure detection and convergence reporting (E-04, E-08, S1), restarts
## (E-08, E-11, S2), parameter names (E-01), formula environments (E-13), maxspline
## (E-07), unscaled log-likelihoods (E-12), the REML Hessian (E-03), finite-difference
## tolerances (P3), negREML and censored responses (P3), and CompPareto support (C-6).
##
## References are closed-form probability mass functions written out here, optim() on
## those closed forms, numDeriv, or constants from independent reviews.

## ---- independent closed forms ----------------------------------------------------

## discrete EGPD, carrier G(v) = v^kappa: P(Y = y) = G(H(y + 1)) - G(H(y)), H the GPD cdf
ldegpd1_ref <- function(y, sigma, xi, kappa) {
  H <- function(t) 1 - (1 + xi * t / sigma)^(-1 / xi)
  log(H(y + 1)^kappa - H(y)^kappa)
}
lzidegpd1_ref <- function(y, pi, sigma, xi, kappa) {
  p <- exp(ldegpd1_ref(y, sigma, xi, kappa))
  log(ifelse(y == 0, pi + (1 - pi) * p, (1 - pi) * p))
}
lztdegpd1_ref <- function(y, sigma, xi, kappa)
  ldegpd1_ref(y, sigma, xi, kappa) - log1p(-exp(ldegpd1_ref(0, sigma, xi, kappa)))
## generalized Waring, mean parameterisation a = mu (rho - 1) / k
lgw_ref <- function(y, mu, k, rho) {
  a <- mu * (rho - 1) / k
  lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) - lgamma(a) - lgamma(k) +
    lgamma(a + y) + lgamma(k + y) - lgamma(a + k + rho + y) - lgamma(y + 1)
}
lztgw_ref <- function(y, mu, k, rho)
  lgw_ref(y, mu, k, rho) - log1p(-exp(lgw_ref(0, mu, k, rho)))

## ---- E-04: zero-truncated Waring and invalid starts ---------------------------------

test_that("E-04: ztgwaring leaves its start and reaches the zero-truncated MLE", {
  set.seed(1)
  y <- rgwaring(1500, 4, 2, 4)
  y <- y[y > 0][1:300]
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "ztgwaring")
  b <- unname(fit$coefficients)
  expect_equal(fit$convergence, 0L)
  expect_equal(fit$logLik, sum(lztgw_ref(y, exp(b[1]), exp(b[2]), exp(b[3]))),
               tolerance = 1e-10)
  nll <- function(th) -sum(lztgw_ref(y, exp(th[1]), exp(th[2]), exp(th[3])))
  ref <- optim(c(log(4), 0, log(3)), nll, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 1000))
  expect_lt(abs(fit$logLik + ref$value), 1e-4)
  expect_equal(fit$logLik, -721.983, tolerance = 1e-6)
})

test_that("E-04: starting values with a non-finite likelihood are an error", {
  set.seed(1)
  y <- rdiscegpd(300, sigma = 3, xi = 0.2, kappa = 1.5)
  ## xi = -2 (identity link) bounds the support below most of the counts
  expect_error(egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "degpd",
                    degpd.args = list(m = 1, link = "identity"), inits = c(log(3), -2, 0)),
               "non-finite log-likelihood")
})

test_that("E-04: a start outside a count family's clamp range is an error", {
  ## beyond a clamp the likelihood is evaluated at the bound with a zero score, so the
  ## parameter never moved: gwaring from rho = exp(0) = 1 reported convergence at
  ## logLik -764.03 against an MLE of -719.88
  set.seed(1)
  y <- rgwaring(1500, 4, 2, 4)
  yp <- y[y > 0][1:300]
  rho_lo <- egpd:::cf_bounds_cpp()$gwaring[["rho_lo"]]
  for (fam in c("ztgwaring", "gwaring")) {
    d <- data.frame(y = if (fam == "gwaring") y[1:300] else yp)
    expect_error(egpd(list(y ~ 1, ~1, ~1), d, family = fam, inits = c(0, 0, 0)),
                 "put rho at 1, at or below its lower bound .*outside the parameter space \\(rho > 1\\)")
  }
  ## the bounds are the C++ clamps, on either side
  set.seed(2)
  yg <- rgpois(300, mu = 5, lambda = 0.3)
  expect_error(egpd(list(y ~ 1, ~1), data.frame(y = yg), family = "gpois", inits = c(log(5), -30)),
               "put lambda at .*outside \\[1e-10, ")
  set.seed(4)
  expect_error(egpd(y ~ 1, data.frame(y = rprl(300, 8)), family = "prl", inits = log(3)),
               "outside the parameter space \\(mu > 4\\)")
  ## a start inside the box is unaffected and reaches the zero-truncated MLE
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = yp), family = "ztgwaring",
              inits = c(log(4), 0, log(2.5)))
  b <- unname(fit$coefficients)
  expect_gt(exp(b[3]), rho_lo)
  expect_equal(fit$logLik, sum(lztgw_ref(yp, exp(b[1]), exp(b[2]), exp(b[3]))), tolerance = 1e-10)
  expect_equal(fit$logLik, -721.983, tolerance = 1e-6)
})

test_that("E-04: a fit left on the generalized Waring's rho bound is not reported as converged", {
  ## let the out-of-range start through to show what the end-of-fit check sees
  local_mocked_bindings(.check_count_box = function(...) invisible(NULL), .package = "egpd")
  set.seed(1)
  y <- rgwaring(1500, 4, 2, 4)[1:300]
  expect_warning(fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "gwaring",
                             inits = c(0, 0, 0)),
                 class = "egpd_nonconvergence")
  expect_equal(fit$convergence, 3L)
  expect_true(is.na(fit$logLik))
  expect_match(fit$convergence.report, "rho reached its lower bound")
})

test_that("E-04: the optimisers never report an invalid objective as convergence", {
  ctl <- list(steptol = 1e-12, itlim = 50, fntol = 1e-8, gradtol = 1e-6, stepmax = 3)
  ## the package's objectives return 1e20 with a zero gradient at invalid points
  f <- function(p) 1e20
  s <- function(p, kept) structure(0 * p, gradient = 0 * p, Hessian = diag(length(p)))
  fit <- egpd:::.newton_step_inner(c(1, 2), f, s, control = ctl)
  expect_equal(fit$convergence, 2L)
  expect_false(fit$gradconv)
  fit2 <- egpd:::.BFGS(c(1, 2), f, function(p) 0 * p, control = ctl)
  expect_equal(fit2$convergence, 2L)
})

test_that("E-04: a positive log-likelihood for a count family is reported as failure", {
  set.seed(1)
  y <- rgwaring(600, 4, 2, 4)
  y <- y[y > 0][1:100]
  ## a likelihood that has broken down numerically: P(y) > 1 everywhere
  flat <- function(pars, likdata) {
    out <- matrix(0, nrow(likdata$y), 9)
    out[, c(4, 7, 9)] <- 1
    out
  }
  broken <- list(d0 = function(pars, likdata) -5, d120 = flat, d340 = NULL, m = 1,
                 iG = function(...) NA)
  real <- egpd:::.setup.family.egpd
  local_mocked_bindings(.setup.family.egpd = function(...) {
    out <- real(...)
    out$lik.fns <- broken
    out
  }, .package = "egpd")
  expect_warning(fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "ztgwaring"),
                 class = "egpd_nonconvergence")
  expect_equal(fit$convergence, 2L)
  expect_true(is.na(fit$logLik))
})

## ---- E-08: convergence warnings, ZIDEGPD starts and restarts -------------------------

test_that("E-08: a fit stopped by the iteration limit warns and reports no logLik", {
  set.seed(615)
  y <- rdiscegpd(200, sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  expect_warning(
    fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "degpd",
                restarts = FALSE, control = list(inner = list(itlim = 1))),
    class = "egpd_nonconvergence")
  expect_equal(fit$convergence, 1L)
  expect_true(is.na(fit$logLik))
  expect_match(fit$convergence.report, "iteration limit")
})

test_that("E-08: the default ZIDEGPD fit with covariates reaches the MLE", {
  set.seed(3)
  n <- 500
  x <- runif(n)
  z <- rnorm(n)
  y <- rzidiscegpd(n, pi = plogis(-1 + 0.8 * z), sigma = exp(1 + 0.3 * z), xi = 0.2,
                   kappa = 1.5, type = 1)
  d <- data.frame(y = y, z = z)
  ll <- function(b) sum(lzidegpd1_ref(y, plogis(b[5] + b[6] * z), exp(b[1] + b[2] * z),
                                      exp(b[3]), exp(b[4])))
  fit <- egpd(list(y ~ z, ~1, ~1, ~z), data = d, family = "zidegpd")
  b <- unname(fit$coefficients)
  expect_equal(fit$convergence, 0L)
  expect_equal(fit$logLik, ll(b), tolerance = 1e-10)
  ref <- optim(b, function(bb) -ll(bb), method = "BFGS",
               control = list(reltol = 1e-14, maxit = 500))
  expect_lt(-ref$value - fit$logLik, 1e-4)
  ## reviewer's value from good starting values; the default start gave logLik NA
  expect_equal(fit$logLik, -985.8561, tolerance = 1e-7)
  ## the single (warm) start now gets there too
  fit1 <- egpd(list(y ~ z, ~1, ~1, ~z), data = d, family = "zidegpd", restarts = FALSE)
  expect_equal(fit1$logLik, fit$logLik, tolerance = 1e-7)
})

test_that("E-08: intercept-only ZIDEGPD models 5 and 6 converge at the kappa -> 0 boundary", {
  skip_on_cran()
  set.seed(7)
  n <- 500
  x <- runif(n)
  sig <- exp(1 + 0.8 * sin(2 * pi * x))
  yc <- sig / 0.2 * ((1 - runif(n)^(1 / 1.5))^(-0.2) - 1)
  yz <- ifelse(runif(n) < 0.25, 0, floor(yc))
  d <- data.frame(y = yz)
  ## optima from 54-start Nelder-Mead + BFGS searches of each likelihood (both on the
  ## kappa -> 0 boundary); before, m5 stopped at the iteration limit and m6 at -1229.87
  best <- c(`5` = -1173.4734, `6` = -1172.1850)
  for (m in 5:6) {
    fit <- egpd(list(y ~ 1, ~1, ~1, ~1), data = d, family = "zidegpd",
                zidegpd.args = list(m = m))
    expect_equal(fit$convergence, 0L)
    expect_gt(fit$logLik, best[[as.character(m)]] - 1e-3)
    p <- predict(fit, type = "response")[1, ]
    direct <- sum(log(dzidiscegpd(yz, pi = p$pi, sigma = p$scale, xi = p$shape,
                                  kappa = p$kappa, type = m - 3L)))
    expect_equal(fit$logLik, direct, tolerance = 1e-8)
  }
})

test_that("S2: the ZIDEGPD warm start does not count offsets twice", {
  set.seed(21)
  n <- 400
  lo <- log(runif(n, 20, 80))
  y <- rzidiscegpd(n, pi = 0.2, sigma = exp(-2 + lo), xi = 0.2, kappa = 1.5, type = 1)
  d <- data.frame(y = y, lo = lo)
  dfit <- egpd(list(y ~ offset(lo), ~1, ~1), d, family = "degpd", restarts = FALSE)
  seen <- NULL
  real <- egpd:::.setup.inner.inits
  local_mocked_bindings(.setup.inner.inits = function(inits, likdata, likfns, npar, family, ...) {
    if (identical(as.vector(family), "zidegpd")) seen <<- inits
    real(inits, likdata, likfns, npar, family, ...)
  }, .package = "egpd")
  suppressWarnings(egpd(list(y ~ offset(lo), ~1, ~1, ~1), d, family = "zidegpd",
                        restarts = FALSE))
  ## the intercept start is the DEGPD intercept, not intercept + mean(offset)
  expect_equal(seen[1:3], unname(dfit$coefficients[1:3]), tolerance = 1e-6)
})

test_that("E-08: a Newton step clipped into an ascent direction is rescaled instead", {
  ## quadratic with an exact Hessian whose Newton step s has g's > 0, while clipping
  ## each component to stepmax = 100 gives g'clip(s) < 0 (as on a ZT-DEGPD fit)
  s <- c(-3.44, 245.9, 300)
  g <- c(39.6, 1.12, 0.016)
  H <- tcrossprod(g) / sum(g * s) + 0.5 * (diag(3) - tcrossprod(s) / sum(s^2))
  p0 <- c(1, 2, 3)
  pstar <- p0 - s
  f <- function(p) 0.5 * sum((p - pstar) * (H %*% (p - pstar)))
  sfn <- function(p, kept) {
    gr <- as.vector(H %*% (p - pstar))
    structure(solve(H, gr), gradient = gr, Hessian = H)
  }
  ctl <- list(steptol = 1e-12, itlim = 100, fntol = 1e-14, gradtol = 1e-9, stepmax = 100)
  fit <- egpd:::.newton_step_inner(p0, f, sfn, control = ctl)
  expect_equal(fit$pars, pstar, tolerance = 1e-8)
  expect_equal(fit$convergence, 0L)
})

test_that("E-08: the default ZT-DEGPD start no longer stalls far from the MLE", {
  set.seed(1)
  y <- rdiscegpd(900, sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  y <- y[y > 0][1:300]
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "ztdegpd")
  b <- unname(fit$coefficients)
  nll <- function(th) -sum(lztdegpd1_ref(y, exp(th[1]), exp(th[2]), exp(th[3])))
  ref <- optim(c(0.9, -1.3, 0.6), nll, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 1000))
  expect_equal(fit$convergence, 0L)
  expect_equal(fit$logLik, -nll(b), tolerance = 1e-10)
  ## before: "step tolerance reached" at logLik -760.28 with a gradient of 40
  expect_lt(abs(fit$logLik + ref$value), 1e-4)
})

test_that("E-08: a final fit that stops with a real gradient is not reported as converged", {
  ## the objective ignores p[2] but the gradient claims a slope there: no step can
  ## reduce it, so the search ends on its step tolerance with a large Newton decrement,
  ## and a gradient search cannot reduce it either
  fn <- function(pars, likdata, likfns, newton = TRUE) (pars[1] - 1)^2
  sfn <- function(pars, likfns, likdata, kept, newton = TRUE) {
    g <- c(2 * (pars[1] - 1), 0.01)
    if (!newton) return(structure(g, gradient = g))
    H <- diag(c(2, 1e-10))
    structure(solve(H, g), gradient = g, Hessian = H)
  }
  ctl <- list(steptol = 1e-12, itlim = 50, fntol = 1e-8, gradtol = 1e-4, stepmax = 100)
  fit <- egpd:::.fit_inner(c(0, 0), fn, sfn, likdata = list(control = list(inner = ctl)),
                           likfns = NULL)
  expect_equal(fit$pars[1], 1, tolerance = 1e-6)
  expect_equal(fit$convergence, 3L)
  expect_match(fit$report, "predicted by the Newton step")
  expect_match(fit$report, "the objective and its gradient disagree")
})

## Synthetic objectives with exact gradients for the check .fit_inner() makes before it
## reports a stalled fit (code 3) or one out of iterations (code 1).
fit_inner_ctl <- function(...)
  list(control = list(inner = modifyList(list(steptol = 1e-12, itlim = 100, fntol = 1e-8,
                                              gradtol = 1e-4, stepmax = 100), list(...))))
fit_inner_fns <- function(f, grad, hess) list(
  fn = function(pars, likdata, likfns, newton = TRUE) f(pars),
  sfn = function(pars, likfns, likdata, kept, newton = TRUE) {
    g <- grad(pars)
    if (!newton) return(structure(g, gradient = g))
    H <- hess(pars)
    structure(solve(H, g), gradient = g, Hessian = H)
  })

test_that("E-08: a decrease predicted only by an outer-product-like Hessian is not non-convergence", {
  ## converged as far as the objective can tell (1e9 + a remaining 2e-8, below its
  ## rounding), but the Hessian underestimates the curvature along p2 by 2e10, so the
  ## Newton step predicts a decrease of 0.005; this used to be reported as code 3
  q <- fit_inner_fns(function(p) 1e9 + sum((p - 1)^2), function(p) 2 * (p - 1),
                     function(p) diag(c(2, 1e-10)))
  p0 <- c(1 + 1.5e-4, 1 + 5e-7)
  fit <- egpd:::.fit_inner(p0, q$fn, q$sfn, fit_inner_ctl(), NULL)
  expect_equal(fit$convergence, 0L)
  expect_match(fit$report, "not realisable")
  expect_equal(fit$pars, p0, tolerance = 1e-12)
  ## a real remaining decrease (Newton stopped on a coarse function tolerance) is
  ## found by the gradient search, which then meets the gradient criterion
  q <- fit_inner_fns(function(p) sum((p - 1)^2), function(p) 2 * (p - 1),
                     function(p) diag(c(2, 1e-4)))
  fit <- egpd:::.fit_inner(c(3, 2), q$fn, q$sfn, fit_inner_ctl(fntol = 1e-3), NULL)
  expect_equal(fit$convergence, 0L)
  expect_match(fit$report, "gradient tolerance reached by a gradient search")
  expect_equal(fit$pars, c(1, 1), tolerance = 1e-6)
})

test_that("E-08: after a gradient search gains, Newton resumes from the better point", {
  ## quartic in p2 with its exact Hessian and a coarse function tolerance: Newton stops
  ## with a real decrease left; a short (two-iteration) search gains part of it and
  ## Newton resumes from there
  q <- fit_inner_fns(function(p) (p[1] - 1)^2 + (p[2] - 1)^4,
                     function(p) c(2 * (p[1] - 1), 4 * (p[2] - 1)^3),
                     function(p) diag(c(2, max(12 * (p[2] - 1)^2, 1e-12))))
  real <- egpd:::.polish_inner
  calls <- 0
  local_mocked_bindings(.polish_inner = function(fit, fn, sfn, likdata, likfns, gradtol, maxit = 100) {
    calls <<- calls + 1
    real(fit, fn, sfn, likdata, likfns, gradtol, maxit = 2)
  }, .package = "egpd")
  fit <- egpd:::.fit_inner(c(1, 2), q$fn, q$sfn, fit_inner_ctl(fntol = 0.2, gradtol = 1e-8), NULL)
  expect_equal(calls, 1)
  expect_equal(fit$convergence, 0L)
  expect_match(fit$report, "Newton resumed after a gradient search")
  ## before: code 3 at objective 0.0066
  expect_lt(as.numeric(fit$objective), 1e-3)
})

test_that("E-08: a ridge along which the objective keeps decreasing stays code 3", {
  ## constant slope 0.01 in p2, whose (zero) curvature the Hessian puts at 1e-6; Newton
  ## stops on a coarse function tolerance, and each gradient search gains 0.01 more
  r <- fit_inner_fns(function(p) (p[1] - 1)^2 - 0.01 * p[2], function(p) c(2 * (p[1] - 1), -0.01),
                     function(p) diag(c(2, 1e-6)))
  real <- egpd:::.polish_inner
  calls <- 0
  local_mocked_bindings(.polish_inner = function(...) {
    calls <<- calls + 1
    real(...)
  }, .package = "egpd")
  fit <- egpd:::.fit_inner(c(1, 0), r$fn, r$sfn, fit_inner_ctl(fntol = 2), NULL)
  expect_equal(calls, 2)
  expect_equal(fit$convergence, 3L)
  expect_match(fit$report, "keep decreasing the objective")
  ## the better point is kept
  expect_lt(as.numeric(fit$objective), -2)
})

test_that("E-08: fits that run out of Newton iterations next to the MLE converge", {
  ## generalized Waring, outer-product Hessian: Newton converges linearly and used to stop
  ## at the iteration limit (logLik NA) within 3e-4 of the MLE
  for (s in 19:20) {
    set.seed(s)
    n <- 120
    x1 <- runif(n, -1, 1)
    y <- rgwaring(n, mu = exp(1 + 0.3 * x1), k = 2, rho = 4)
    fit <- egpd(list(y ~ x1, ~1, ~1), data.frame(y = y, x1 = x1), family = "gwaring")
    b <- unname(fit$coefficients)
    nll <- function(th) -sum(lgw_ref(y, exp(th[1] + th[2] * x1), exp(th[3]), exp(th[4])))
    ref <- optim(b, nll, method = "BFGS", control = list(reltol = 1e-15, maxit = 2000))
    expect_equal(fit$convergence, 0L)
    expect_equal(fit$logLik, -nll(b), tolerance = 1e-10)
    expect_lt(abs(fit$logLik + ref$value), 1e-6)
  }
})

test_that("E-08: ZT-DEGPD fits are reported as converged only where no real decrease remains", {
  ## R2-HURDLE generator. With the outer-product Hessian the default start can stall at
  ## kappa -> 0, where the zero-truncated likelihood loses precision: seed 30 within
  ## 1e-3 of the optimum there (reported as code 3 before), seed 5 three units short of
  ## an interior optimum, which must not be reported as converged
  lzt <- function(y, sigma, xi, kappa) {
    H <- function(t) -expm1(-log1p(xi * t / sigma) / xi)
    a <- log(H(y + 1))
    kappa * a + log(-expm1(kappa * (log(H(y)) - a))) - log(-expm1(kappa * log(H(1))))
  }
  for (s in c(5, 30)) {
    set.seed(s)
    n <- 120
    x1 <- runif(n, -1, 1)
    y <- ifelse(runif(n) < 0.35, 0L,
                rdiscegpd(n, sigma = exp(1 + 0.3 * x1), xi = 0.15, kappa = 1.2) + 1L)
    d <- data.frame(y = y, x1 = x1)[y > 0, ]
    fit <- suppressWarnings(egpd(list(y ~ x1, ~1, ~1), d, family = "ztdegpd"))
    b <- unname(fit$coefficients)
    nll <- function(bb) {
      v <- -sum(lzt(d$y, exp(bb[1] + bb[2] * d$x1), exp(bb[3]), exp(bb[4])))
      if (is.finite(v)) v else 1e10
    }
    o <- optim(b, nll, method = "BFGS", control = list(reltol = 1e-15, maxit = 5000))
    o <- optim(o$par, nll, method = "Nelder-Mead", control = list(reltol = 1e-15, maxit = 20000))
    o <- optim(o$par, nll, method = "BFGS", control = list(reltol = 1e-15, maxit = 5000))
    if (s == 30)
      expect_equal(fit$convergence, 0L)
    if (fit$convergence == 0L) {
      expect_equal(fit$logLik, -nll(b), tolerance = 1e-5)
      expect_lt(nll(b) - o$value, 5e-3)
    }
  }
})

test_that("E-08: the stall check does not run for fits that converge", {
  calls <- 0
  real <- egpd:::.polish_inner
  local_mocked_bindings(.polish_inner = function(...) {
    calls <<- calls + 1
    real(...)
  }, .package = "egpd")
  set.seed(615)
  y <- rdiscegpd(200, sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  f1 <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "degpd")
  set.seed(1)
  x1 <- runif(400, -1, 1)
  yg <- rgwaring(400, mu = exp(1 + 0.3 * x1), k = 2, rho = 4)
  f2 <- egpd(list(y ~ x1, ~1, ~1), data.frame(y = yg, x1 = x1), family = "gwaring")
  f3 <- egpd(list(y ~ x1, ~1), data.frame(y = rgpois(400, mu = exp(1 + 0.3 * x1), lambda = 0.3),
                                          x1 = x1), family = "gpois")
  expect_equal(c(f1$convergence, f2$convergence, f3$convergence), c(0L, 0L, 0L))
  expect_equal(calls, 0)
})

test_that("E-08: ZIDEGPD restarts run only when the default fit fails, except for m = 3, 4", {
  count_fits <- function(expr) {
    calls <- 0
    real <- egpd:::.setup.family.egpd
    local_mocked_bindings(.setup.family.egpd = function(...) {
      calls <<- calls + 1
      real(...)
    }, .package = "egpd")
    force(expr)
    calls
  }
  ## well-behaved data: the default fit converges in the interior, so the four
  ## extra starts are not fitted (outer call, DEGPD warm start, primary fit)
  set.seed(1001)
  n <- 400
  z <- rnorm(n)
  y <- rzidiscegpd(n, pi = 0.2, sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  d <- data.frame(y = y)
  nfit <- count_fits(fit <- egpd(list(y ~ 1, ~1, ~1, ~1), d, family = "zidegpd"))
  expect_lte(nfit, 3)
  nll <- function(b) -sum(lzidegpd1_ref(y, plogis(b[4]), exp(b[1]), exp(b[2]), exp(b[3])))
  ref <- optim(unname(fit$coefficients), nll, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 500))
  expect_equal(fit$convergence, 0L)
  expect_lt(abs(fit$logLik + ref$value), 1e-6)

  ## the default fit used to stop at the iteration limit here (logLik NA), which
  ## triggered the restarts; the gradient search after Newton now reaches the MLE
  set.seed(1078)
  z <- rnorm(n)
  x <- runif(n)
  y <- rzidiscegpd(n, pi = plogis(-1 + 0.8 * z), sigma = exp(1 + 0.3 * z), xi = 0.2,
                   kappa = 1.5, type = 1)
  d <- data.frame(y = y, z = z)
  fl <- list(y ~ z, ~1, ~1, ~z)
  expect_no_warning(fit <- egpd(fl, d, family = "zidegpd", restarts = FALSE))
  ll <- function(b) sum(lzidegpd1_ref(y, plogis(b[5] + b[6] * z), exp(b[1] + b[2] * z),
                                      exp(b[3]), exp(b[4])))
  ref <- optim(unname(fit$coefficients), function(b) -ll(b), method = "BFGS",
               control = list(reltol = 1e-14, maxit = 500))
  expect_equal(fit$convergence, 0L)
  expect_equal(fit$logLik, ll(unname(fit$coefficients)), tolerance = 1e-10)
  expect_lt(-ref$value - fit$logLik, 1e-4)

  ## a default fit on a boundary (no zero inflation: logit(pi) -> -Inf) triggers them
  set.seed(3)
  y <- rdiscegpd(n, sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  expect_message(fit <- egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = y), family = "zidegpd",
                             trace = 1),
                 "ZIDEGPD restarts: the primary fit has logitpi on a boundary")
  ## the supremum is the DEGPD likelihood (pi -> 0)
  nll <- function(b) -sum(ldegpd1_ref(y, exp(b[1]), exp(b[2]), exp(b[3])))
  ref <- optim(unname(fit$coefficients[1:3]), nll, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 500))
  expect_equal(fit$convergence, 0L)
  expect_lt(abs(fit$logLik + ref$value), 1e-5)

  ## carrier 3: the warm start converges to an inferior interior optimum, so all five
  ## starts are always fitted
  set.seed(3008)
  z <- rnorm(n)
  x <- runif(n)
  y3 <- rzidiscegpd(n, pi = 0.2, sigma = 3, xi = 0.2, delta = 2, type = 4)
  d3 <- data.frame(y = y3)
  one <- egpd(list(y ~ 1, ~1, ~1, ~1), d3, family = "zidegpd",
              zidegpd.args = list(m = 3), restarts = FALSE)
  nfit <- count_fits(all5 <- egpd(list(y ~ 1, ~1, ~1, ~1), d3, family = "zidegpd",
                                   zidegpd.args = list(m = 3)))
  expect_equal(nfit, 7)
  expect_gt(all5$logLik, one$logLik + 3)
  p <- predict(all5, type = "response")[1, ]
  expect_equal(all5$logLik, sum(log(dzidiscegpd(y3, pi = p$pi, sigma = p$scale, xi = p$shape,
                                                delta = p$delta, type = 4))),
               tolerance = 1e-8)
})

## ---- E-01: parameter names ------------------------------------------------------------

test_that("E-01: every family accepts natural, internal and long parameter names in any order", {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  fams <- list(
    list("egpd", list(m = 1)), list("egpd", list(m = 2)), list("egpd", list(m = 3)),
    list("egpd", list(m = 4)), list("degpd", NULL, list(m = 1)),
    list("degpd", NULL, list(m = 1, link = "identity")), list("degpd", NULL, list(m = 2)),
    list("degpd", NULL, list(m = 4)), list("zidegpd", NULL, NULL, list(m = 1)),
    list("zidegpd", NULL, NULL, list(m = 4)), list("pig"), list("zipig"), list("gpig"),
    list("gpignat"), list("zigpig"), list("zigpignat"), list("bell"), list("bellnat"),
    list("zibell"), list("zibellnat"), list("cmp"), list("prl"), list("gpois"),
    list("gwaring"), list("ztgwaring"), list("ztdegpd"),
    list("ztdegpd", NULL, list(m = 1, link = "identity")), list("plnorm"),
    list("comppareto", NULL, NULL, NULL, list(spec = "lnorm")),
    list("comppareto", NULL, NULL, NULL, list(spec = "gamma")),
    list("comppareto", NULL, NULL, NULL, list(spec = "weibull")),
    list("comppareto", NULL, NULL, NULL, list(spec = "exp")))
  natural <- list(
    egpd = list(c("psi", "xi", "kappa"), c("psi", "xi", "kappa1", "dkappa", "p"),
                c("psi", "xi", "delta"), c("psi", "xi", "delta", "kappa")),
    comppareto = list(c("meanlog", "sdlog", "alpha", "theta"),
                      c("shape", "scale", "alpha", "theta"),
                      c("shape", "scale", "alpha", "theta"), c("rate", "alpha", "theta")))
  for (a in fams) {
    a <- c(a, rep(list(NULL), 5 - length(a)))
    fi <- egpd:::.setup.family.egpd(a[[1]], a[[2]] %||% list(), a[[3]] %||% list(),
                                    a[[4]] %||% list(), a[[5]] %||% list(), NULL, list())
    np <- fi$npar
    fl <- lapply(seq_len(np), function(i)
      if (i == 1) y ~ x else stats::as.formula(paste0("~ z", i)))
    expect <- c("y~x", paste0("~z", seq_len(np))[-1])
    for (spelling in list(fi$nat, fi$nms, fi$nms2)) {
      names(fl) <- spelling
      r <- egpd:::.reorder_named_formulae(rev(fl), fi$ptab)
      got <- gsub(" ", "", vapply(r, function(f) paste(deparse(f), collapse = ""), ""))
      expect_equal(unname(got), expect, info = paste(a[[1]], paste(spelling, collapse = ",")))
    }
  }
  ## the natural names are the ones used in the papers and in predict() output
  fi <- egpd:::.setup.family.egpd("comppareto", list(), list(), list(),
                                  list(spec = "lnorm"), NULL, list())
  expect_equal(fi$nat, natural$comppareto[[1]])
  fi <- egpd:::.setup.family.egpd("gpois", list(), list(), list(), list(), NULL, list())
  expect_equal(fi$nat, c("mu", "lambda"))
})

test_that("E-01: gpois lambda can be named, in either order", {
  set.seed(2)
  y <- rgpois(300, mu = 4, lambda = 0.3)
  d <- data.frame(y = y)
  base <- egpd(list(y ~ 1, ~1), d, family = "gpois")
  for (fl in list(list(mu = y ~ 1, lambda = ~1), list(lambda = ~1, mu = y ~ 1))) {
    fit <- egpd(fl, d, family = "gpois")
    expect_equal(fit$coefficients, base$coefficients)
  }
})

test_that("E-01: CompPareto parameters can be fixed by their natural names", {
  set.seed(5)
  yc <- rcomppareto(300, spec = "lnorm", alpha = 2, theta = 3, meanlog = 0.5, sdlog = 0.6)
  fit <- egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = yc), family = "comppareto",
              fixed = list(sdlog = 0.6, theta = 3))
  p <- predict(fit, type = "response")[1, ]
  expect_equal(p$sdlog, 0.6)
  expect_equal(p$theta, 3)
  ## names are no longer mangled by gsub(): "logsdlog" -> "sdlog", not "sd"
  expect_equal(names(fit$call), c("meanlog", "sdlog", "alpha", "theta"))
})

test_that("E-01: legacy lpi and lxi are accepted with a deprecation warning", {
  set.seed(616)
  y <- rzidiscegpd(300, pi = 0.2, sigma = 3, xi = 0.25, kappa = 1.5, type = 1)
  d <- data.frame(y = y)
  ref <- egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, logitpi = ~1), d,
              family = "zidegpd", restarts = FALSE)
  expect_warning(
    fit <- egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, lpi = ~1), d,
                family = "zidegpd", restarts = FALSE),
    class = "deprecatedWarning")
  expect_equal(fit$coefficients, ref$coefficients)

  yd <- rdiscegpd(300, sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  dd <- data.frame(y = yd)
  ref <- egpd(list(y ~ 1, ~1, ~1), dd, family = "degpd", restarts = FALSE,
              degpd.args = list(link = "identity"))
  expect_warning(
    fit <- egpd(list(lkappa = ~1, lxi = ~1, lsigma = y ~ 1), dd, family = "degpd",
                restarts = FALSE, degpd.args = list(link = "identity")),
    class = "deprecatedWarning")
  expect_equal(fit$coefficients, ref$coefficients)
})

test_that("E-01: unknown or repeated names are clear errors", {
  d <- data.frame(y = rep(0:5, 10))
  expect_error(egpd(list(sigma = y ~ 1, tail = ~1, kappa = ~1), d, family = "degpd",
                    restarts = FALSE),
               "valid names: sigma \\(or lsigma\\), xi \\(or lxi\\), kappa \\(or lkappa\\)")
  expect_error(egpd(list(y ~ 1, ~1, ~1), d, family = "degpd", restarts = FALSE,
                    fixed = list(kappa = 1, lkappa = 2)),
               "more than once")
})

## ---- E-11: restarts with a named list ----------------------------------------------

test_that("E-11: DEGPD restarts run when the response formula is not listed first", {
  set.seed(12)
  n <- 300
  x <- runif(n)
  y <- rdiscegpd(n, sigma = exp(1 + x), xi = 0.3, kappa = 1.5, type = 1)
  d <- data.frame(y = y, x = x)
  canonical <- egpd(list(lsigma = y ~ x, lxi = ~1, lkappa = ~1), d, family = "degpd")
  calls <- 0
  real <- egpd:::.setup.family.egpd
  local_mocked_bindings(.setup.family.egpd = function(...) {
    calls <<- calls + 1
    real(...)
  }, .package = "egpd")
  expect_no_warning(
    fit <- egpd(list(lxi = ~1, lsigma = y ~ x, lkappa = ~1), d, family = "degpd"))
  ## one setup per candidate start: the five restarts did run
  expect_gte(calls, 5)
  expect_equal(fit$logLik, canonical$logLik)
})

## ---- E-13: formula environments --------------------------------------------------

test_that("E-13: smooths may use the caller's local variables in any slot", {
  set.seed(1)
  n <- 300
  d <- data.frame(x = runif(n))
  d$y <- rdiscegpd(n, sigma = exp(1 + sin(2 * pi * d$x)), xi = 0.2, kappa = 1.5, type = 1)
  run <- function(fl) {
    kk <- 6
    egpd(eval(fl), data = d, family = "degpd", restarts = FALSE,
         fixed = list(kappa = 1.5))
  }
  a <- run(quote(list(y ~ 1, ~ s(x, k = kk), ~1)))
  expect_equal(a[[2]]$smooth[[1]]$bs.dim, 6)
  b <- run(quote(list(lkappa = y ~ 1, lsigma = ~ s(x, k = kk), lxi = ~1)))
  expect_equal(b[[1]]$smooth[[1]]$bs.dim, 6)
})

## ---- E-07: maxspline -----------------------------------------------------------------

test_that("E-07: maxspline subsampling fits and predicts", {
  set.seed(9)
  n <- 800
  x <- round(runif(n), 2)
  f <- factor(sample(c("a", "b", "c"), n, TRUE))
  y <- rdiscegpd(n, sigma = exp(1 + sin(2 * pi * x) + 0.3 * (f == "b")), xi = 0.2,
                 kappa = 1.5, type = 1)
  d <- data.frame(y = y, x = x, f = f)
  fit <- egpd(list(y ~ s(x, k = 6) + f, ~1, ~1), d, family = "degpd", restarts = FALSE,
              maxspline = 300)
  expect_equal(fit$convergence, 0L)
  expect_equal(nobs(fit), n)
  p <- predict(fit, type = "response")
  expect_equal(fit$logLik, sum(ldegpd1_ref(y, p$scale, p$shape, p$kappa)),
               tolerance = 1e-10)
  pn <- predict(fit, newdata = d[1:5, ], type = "response")
  expect_equal(unname(as.matrix(pn[, c("scale", "shape", "kappa")])),
               unname(as.matrix(p[1:5, c("scale", "shape", "kappa")])))
  ## parametric terms have no data-dependent basis, so subsampling cannot change them
  g1 <- egpd(list(y ~ x + f, ~1, ~1), d, family = "degpd", restarts = FALSE,
             maxspline = 300)
  g0 <- egpd(list(y ~ x + f, ~1, ~1), d, family = "degpd", restarts = FALSE)
  expect_equal(g1$coefficients, g0$coefficients, tolerance = 1e-8)
})

## ---- E-12: unscaled log-likelihood -------------------------------------------------

test_that("E-12: logLik is the data log-likelihood under sandwich adjustment and gamma", {
  set.seed(31)
  n <- 600
  x <- runif(n)
  grp <- rep(1:60, each = 10)
  u <- rnorm(60, 0, 0.3)[grp]
  y <- rdiscegpd(n, sigma = exp(1 + x + u), xi = 0.2, kappa = 1.5, type = 1)
  d <- data.frame(y = y, x = x, grp = grp)
  fm <- egpd(list(y ~ x, ~1, ~1), d, family = "degpd", restarts = FALSE,
             sandwich.args = list(id = "grp", method = "magnitude"))
  b <- unname(fm$coefficients)
  manual <- sum(ldegpd1_ref(y, exp(b[1] + b[2] * x), exp(b[3]), exp(b[4])))
  ## the reviewer saw -1765.21 reported against a true -1842.23
  expect_equal(fm$logLik, manual, tolerance = 1e-10)
  expect_equal(fm$logLik.adjusted, fm$likdata$k * manual, tolerance = 1e-10)
  expect_equal(AIC(fm), -2 * manual + 2 * fm$df, tolerance = 1e-10)

  fg <- egpd(list(y ~ s(x, k = 5), ~1, ~1), d, family = "degpd", restarts = FALSE,
             gamma = 2)
  p <- predict(fg, type = "response")
  expect_equal(fg$logLik, sum(ldegpd1_ref(y, p$scale, p$shape, p$kappa)),
               tolerance = 1e-10)
})

## ---- E-03: REML Hessian --------------------------------------------------------------

test_that("E-03: smoothing-parameter covariance is reliable without third derivatives", {
  set.seed(7)
  n <- 600
  x <- runif(n)
  x2 <- runif(n)
  reml_hessian <- function(fit) {
    ld <- fit$likdata
    ld$control$inner$gradtol <- 1e-10
    f <- function(r) as.numeric(egpd:::.reml0(structure(r, beta = fit$coefficients),
                                              fit$likfns, ld, fit$Sdata))
    numDeriv::hessian(f, log(fit$sp), method.args = list(d = 0.1, r = 4))
  }
  set.seed(11)
  y <- rzidiscegpd(n, pi = 0.2, sigma = exp(1 + sin(2 * pi * x)), xi = 0.2, kappa = 1.5,
                   type = 1)
  fz <- egpd(list(y ~ s(x, k = 8), ~1, ~1, ~1), data.frame(y = y, x = x),
             family = "zidegpd", restarts = FALSE)
  ## REML curvature 2.20 from an independent LAML implementation (core review);
  ## the package used to return -0.196 here and hence Vlsp = 1.229
  expect_equal(as.numeric(fz$Vlsp), 1 / 2.20, tolerance = 0.02)
  skip_if_not_installed("numDeriv")
  expect_equal(as.vector(fz$Vlsp), as.vector(solve(reml_hessian(fz))), tolerance = 1e-3)

  set.seed(11)
  ## PIG data by inversion, the draws rpig() gave before it sampled the mixture exactly:
  ## this curvature comparison was calibrated on them (with other draws the s(x2) smooth
  ## is penalised far harder and numDeriv's step is too coarse for 1e-4)
  yp <- qpig(runif(n), mu = exp(1 + sin(2 * pi * x)), sigma = exp(-1 + x2))
  fp <- egpd(list(y ~ s(x, k = 8), ~ s(x2, k = 6)), data.frame(y = yp, x = x, x2 = x2),
             family = "pig")
  expect_equal(as.vector(fp$Vlsp), as.vector(solve(reml_hessian(fp))), tolerance = 1e-4)
})

test_that("E-03: a REML Hessian that is not positive definite is not inverted blindly", {
  V <- egpd:::.reml_vcov(matrix(c(2, 0, 0, 1e-9), 2))
  expect_equal(unname(V[1:2, 1:2]), diag(c(0.5, 0)))
  expect_warning(V <- egpd:::.reml_vcov(matrix(c(2, 0, 0, -0.2), 2)), "not positive definite")
  expect_equal(unname(V[1:2, 1:2]), diag(c(0.5, 0)))
  ## curvature below twice its own estimated error is not inverted either
  H <- structure(matrix(c(2, 0, 0, -1), 2), error = diag(c(0.001, 0.6)))
  expect_warning(V <- egpd:::.reml_vcov(H), "not smooth enough")
  expect_equal(unname(V[1:2, 1:2]), diag(c(0.5, 0)))
  ## ... but a small curvature measured accurately is, whatever the error of the
  ## strongly curved direction (a single error bound for the matrix discarded it)
  H <- structure(diag(c(2, 0.001)), error = diag(c(0.0017, 1e-6)))
  expect_no_warning(V <- egpd:::.reml_vcov(H))
  expect_equal(unname(V[1:2, 1:2]), diag(c(0.5, 1000)))
})

test_that("E-03: small but well-measured REML curvatures keep their variance", {
  skip_if_not_installed("numDeriv")
  set.seed(1)
  n <- 400
  x <- runif(n)
  z <- runif(n)
  w <- runif(n)
  y <- rpig(n, mu = exp(1 + sin(2 * pi * x) + 0.5 * z), sigma = 0.5)
  expect_no_warning(
    fit <- egpd(list(y ~ s(x, k = 8) + s(z, k = 6) + s(w, k = 6), ~ s(x, k = 5)),
                data.frame(y, x, z, w), family = "pig"))
  ld <- fit$likdata
  ld$control$inner$gradtol <- 1e-10
  f <- function(r) as.numeric(egpd:::.reml0(structure(r, beta = fit$coefficients),
                                            fit$likfns, ld, fit$Sdata))
  Hn <- numDeriv::hessian(f, log(fit$sp), method.args = list(d = 0.1, r = 6))
  ## eigenvalues 2.04, 0.186, 0.00169, 0.00100: the two small ones used to be dropped
  ## (variances ~1e-8 instead of ~600 and ~1000)
  expect_equal(attr(fit$Vlsp, "eigen"), eigen(Hn)$values, tolerance = 0.01)
  expect_equal(diag(as.matrix(fit$Vlsp)), diag(solve(Hn)), tolerance = 0.01)
})

## ---- S1: failed refits inside the outer optimiser ----------------------------------

test_that("S1: BFGS rejects a trial point whose gradient fails instead of aborting", {
  ctl <- list(steptol = 1e-12, itlim = 100, fntol = 1e-14, gradtol = 1e-8, stepmax = 3)
  f <- function(p) sum((p - c(1, -2))^2)
  calls <- 0
  g <- function(p) {
    calls <<- calls + 1
    if (calls == 2) stop("refit failed")
    2 * (p - c(1, -2))
  }
  fit <- egpd:::.BFGS(c(0, 0), f, g, control = ctl)
  expect_equal(fit$par, c(1, -2), tolerance = 1e-6)
  expect_equal(fit$convergence, 0L)
})

curvature_data <- function() {
  set.seed(904)
  d <- data.frame(x = seq(0, 1, length.out = 240))
  scale <- exp(1 + 1.5 * sin(2 * pi * d$x))
  d$y <- floor(scale / 0.3 * ((1 - runif(nrow(d)))^(-0.3) - 1))
  d
}

test_that("S1: a failed profile refit at a trial smoothing parameter is not fatal", {
  args <- list(formula = list(y ~ s(x, k = 7), ~1, ~1), data = curvature_data(),
               family = "degpd", fixed = list(xi = 0.3, kappa = 1), restarts = FALSE,
               outer = "fd")
  ref <- do.call(egpd, args)
  real <- egpd:::.reml1.fd
  calls <- 0
  local_mocked_bindings(.reml1.fd = function(...) {
    calls <<- calls + 1
    if (calls == 2)
      stop("Inner optimization failed while differentiating profile REML.", call. = FALSE)
    real(...)
  }, .package = "egpd")
  fit <- do.call(egpd, args)
  expect_gt(calls, 2)
  expect_equal(fit$convergence, 0L)
  expect_equal(fit$negREML, ref$negREML, tolerance = 1e-6)
})

## ---- P3: finite-difference tolerances, negREML, censoring ----------------------------

test_that("P3: an explicit inner gradtol is respected by finite-difference REML", {
  args <- list(formula = list(y ~ s(x, k = 7), ~1, ~1), data = curvature_data(),
               family = "degpd", fixed = list(xi = 0.3, kappa = 1), restarts = FALSE,
               outer = "fd")
  tols <- numeric(0)
  real <- egpd:::.reml0
  local_mocked_bindings(.reml0 = function(pars, likfns, likdata, Sdata, beta = NULL) {
    tols <<- c(tols, likdata$control$inner$gradtol)
    real(pars, likfns, likdata, Sdata, beta)
  }, .package = "egpd")
  fit <- do.call(egpd, c(args, list(control = list(inner = list(gradtol = 1e-5)))))
  expect_equal(unique(tols), 1e-5)
  tols <- numeric(0)
  fit <- do.call(egpd, args)
  ## defaults: the line search keeps 1e-4, only the differences are tightened
  expect_true(1e-4 %in% tols)
  expect_true(any(tols < 1e-4))
})

test_that("P3: negREML for fixed smoothing parameters is the REML criterion", {
  d <- curvature_data()
  fit <- egpd(list(y ~ s(x, k = 7), ~1, ~1), d, family = "degpd",
              fixed = list(xi = 0.3, kappa = 1), restarts = FALSE)
  fsp <- egpd(list(y ~ s(x, k = 7), ~1, ~1), d, family = "degpd",
              fixed = list(xi = 0.3, kappa = 1), restarts = FALSE, sp = fit$sp)
  expect_equal(fsp$negREML, fit$negREML, tolerance = 1e-8)
})

test_that("P3: censored responses fail early with a clear error", {
  d <- data.frame(l = 1:5, r = 2:6)
  expect_error(egpd(list(cens(l, r) ~ 1, ~1, ~1), d, family = "degpd", restarts = FALSE),
               "censored responses.*not supported")
})

## ---- sandwich adjustment with force = TRUE ----------------------------------------

test_that("sandwich force = TRUE reaches its pseudo-inverse branches instead of failing", {
  set.seed(31)
  n <- 600
  x <- runif(n)
  grp <- rep(1:60, each = 10)
  u <- rnorm(60, 0, 0.3)[grp]
  y <- rdiscegpd(n, sigma = exp(1 + x + u), xi = 0.2, kappa = 1.5, type = 1)
  d <- data.frame(y = y, x = x, x2 = 2 * x, grp = grp, g2 = rep(1:2, length.out = n))
  nll <- function(b) -sum(ldegpd1_ref(y, exp(b[1] + b[2] * x), exp(b[3]), exp(b[4])))
  ref <- optim(c(1, 1, log(0.2), log(1.5)), nll, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 1000))
  ## collinear design: the Hessian of the unpenalised MLE is singular (and so, from
  ## round 6, egpd() warns that the standard errors cannot show it)
  expect_warning(expect_message(
    fa <- egpd(list(y ~ x + x2, ~1, ~1), d, family = "degpd", restarts = FALSE,
               sandwich.args = list(id = "grp", method = "magnitude", force = TRUE)),
    "pseudo-inverse for the sandwich adjustment"),
    class = "egpd_singular_hessian")
  expect_equal(fa$convergence, 0L)
  expect_true(is.finite(fa$likdata$k) && fa$likdata$k > 0)
  expect_equal(fa$logLik, -ref$value, tolerance = 1e-8)
  ## two clusters for four coefficients: the sandwich variance is singular, and so is the
  ## adjusted Hessian (round 6: egpd() warns)
  expect_warning(expect_message(
    fb <- egpd(list(y ~ x, ~1, ~1), d, family = "degpd", restarts = FALSE,
               sandwich.args = list(id = "g2", method = "curvature", force = TRUE)),
    "pseudo-inverse for the curvature adjustment"),
    class = "egpd_singular_hessian")
  expect_equal(fb$convergence, 0L)
  expect_equal(unname(fb$coefficients), ref$par, tolerance = 1e-5)
  expect_equal(fb$logLik, -ref$value, tolerance = 1e-8)
})

## ---- C-6: CompPareto support ------------------------------------------------------------

test_that("C-6: comppareto rejects responses that are not strictly positive", {
  set.seed(5)
  yc <- rcomppareto(100, spec = "lnorm", alpha = 2, theta = 3, meanlog = 0.5, sdlog = 0.6)
  for (bad in list(c(0, yc), c(-1, yc), c(Inf, yc)))
    expect_error(egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = bad), family = "comppareto"),
                 "strictly positive")
})

## ---- bounded DEGPD shape (degpd.args$xi.max) ------------------------------------------

## xi under the bounded link xi = a (1 - exp(-exp(eta) / a))
xi_bounded <- function(eta, a) a * -expm1(-exp(eta) / a)

test_that("xi.max: bounded-shape DEGPD fits reach an interior MLE, not a link asymptote", {
  fl <- list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1)
  for (s in c(4, 7)) {
    set.seed(s)
    y <- rdiscegpd(500, sigma = 3, xi = 0.1, kappa = 2, type = 1)
    d <- data.frame(y = y)
    f2 <- egpd(fl, d, family = "degpd", degpd.args = list(m = 1, xi.max = 0.2))
    f3 <- egpd(fl, d, family = "degpd", degpd.args = list(m = 1, xi.max = 0.3))
    b <- unname(f2$coefficients)
    xi <- xi_bounded(b[2], 0.2)
    nll <- function(th) -sum(ldegpd1_ref(y, exp(th[1]), exp(th[2]), exp(th[3])))
    ref <- optim(c(1.2, log(0.05), 0.6), nll, method = "BFGS",
                 control = list(reltol = 1e-14, maxit = 1000))
    expect_equal(f2$convergence, 0L)
    expect_equal(f2$logLik, -nll(c(b[1], log(xi), b[3])), tolerance = 1e-10)
    ## before: xi -> 0 (eta -72 and -101) at logLik -1292.455 and -1267.944
    expect_lt(abs(f2$logLik - f3$logLik), 1e-3)
    expect_lt(abs(f2$logLik + ref$value), 1e-6)
    expect_gt(xi, 0.01)
    expect_lt(xi, 0.19)
    ## the single default start now starts inside the cap and gets there as well
    f1 <- egpd(fl, d, family = "degpd", degpd.args = list(m = 1, xi.max = 0.2), restarts = FALSE)
    expect_equal(f1$logLik, f2$logLik, tolerance = 1e-8)
  }
})

test_that("xi.max: a best fit on a plateau of the bounded link is restarted from the log-link fit", {
  set.seed(4)
  y <- rdiscegpd(500, sigma = 3, xi = 0.1, kappa = 2, type = 1)
  d <- data.frame(y = y)
  fl <- list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1)
  fu <- egpd(fl, d, family = "degpd", degpd.args = list(m = 1))
  ## put every ordinary start deep on the plateau at the cap (eta = 30, where the
  ## likelihood does not depend on eta at all), so that each one stops there
  real <- egpd:::.shape_eta
  local_mocked_bindings(.shape_eta = function(xi, identity = FALSE, xi.max = Inf, cap = 0.5)
    if (cap == 0.5 && is.finite(xi.max)) 30 else real(xi, identity, xi.max, cap),
    .package = "egpd")
  f0 <- suppressMessages(egpd(fl, d, family = "degpd", degpd.args = list(m = 1, xi.max = 0.2),
                              restarts = FALSE))
  expect_lt(f0$logLik, fu$logLik - 1)
  expect_message(f <- egpd(fl, d, family = "degpd", degpd.args = list(m = 1, xi.max = 0.2),
                           trace = 1),
                 "on the xi -> .* asymptote of the bounded link; trying the log-link fit")
  expect_equal(f$convergence, 0L)
  expect_equal(f$logLik, fu$logLik, tolerance = 1e-8)
  expect_equal(xi_bounded(unname(f$coefficients[2]), 0.2), exp(unname(fu$coefficients[2])),
               tolerance = 1e-4)
})

test_that("xi.max: a genuine xi -> 0 fit stays on the boundary, without warnings", {
  ## a discretised exponential (geometric) sample, whose MLE is at xi -> 0
  set.seed(1)
  y <- floor(rexp(500, 1 / 3))
  fl <- list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1)
  expect_no_warning(f <- egpd(fl, data.frame(y = y), family = "degpd",
                              degpd.args = list(m = 1, xi.max = 0.2)))
  expect_equal(f$convergence, 0L)
  expect_lt(xi_bounded(unname(f$coefficients[2]), 0.2), 1e-4)
  ## the xi -> 0 limit of DEGPD model 1, H(y) = 1 - exp(-y / sigma)
  nll0 <- function(th) -sum(log((-expm1(-(y + 1) / exp(th[1])))^exp(th[2]) -
                                  (-expm1(-y / exp(th[1])))^exp(th[2])))
  ref <- optim(c(1, 0), nll0, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  expect_lt(abs(f$logLik + ref$value), 1e-5)
})
