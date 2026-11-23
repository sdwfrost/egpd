## Tests for distfit() Bernstein polynomial fitting

## Independent Bernstein-EGPD log-likelihood (base R only):
## f(x) = sum_j w_j dbeta(H^kappa; j, m - j + 1) kappa H^(kappa - 1) h(x), with H and
## h the GPD CDF and density; as a function of the optimisation parameters
## (log kappa, log sigma, xi, log(w_j / w_m), j < m).
bern_ll <- function(x, sigma, xi, kappa, w) {
  m <- length(w)
  t <- 1 + xi * x / sigma
  logS <- if (xi == 0) -x / sigma else -log(t) / xi
  H <- -expm1(logS)
  b <- rowSums(vapply(seq_len(m), function(j) w[j] * dbeta(H^kappa, j, m - j + 1),
                      numeric(length(x))))
  sum(log(b) + log(kappa) + (kappa - 1) * log(H) - log(sigma) + logS - log(t))
}
bern_ll_theta <- function(th, x, m) {
  a <- c(th[-(1:3)], 0)
  w <- exp(a - max(a))
  bern_ll(x, exp(th[2]), th[3], exp(th[1]), w / sum(w))
}

test_that("distfit Bernstein method converges on type 1 data", {
  set.seed(42)
  x <- regpd(500, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- distfit(x, type = 1, method = "bernstein", bernstein.m = 8, hessian = FALSE)

  expect_s3_class(fit, "distfit")
  expect_equal(fit$convergence, 0)
  expect_equal(fit$method, "bernstein")
  expect_equal(fit$family, "egpd")
  expect_equal(fit$bernstein.m, 8)
  expect_equal(length(fit$bernstein.weights), 8)
  ## the reported log-likelihood is that of the reported estimate, and the estimate
  ## is a stationary point of the independently computed log-likelihood
  est <- fit$estimate
  w <- fit$bernstein.weights
  th <- c(log(est[["kappa"]]), log(est[["sigma"]]), est[["xi"]], log(w[-8] / w[8]))
  expect_equal(fit$loglik, bern_ll_theta(th, x, 8), tolerance = 1e-10)
  g <- numDeriv::grad(function(t) bern_ll_theta(t, x, 8), th)
  expect_lt(max(abs(g)), 1e-3)
})

test_that("distfit Bernstein with m = 1 is the type-1 EGPD", {
  ## m = 1 gives b(v) = dbeta(v; 1, 1) = 1, i.e. the density kappa H^(kappa - 1) h
  ## of the type-1 EGPD, so both methods must reach the same maximum
  set.seed(123)
  x <- regpd(300, sigma = 1.5, xi = 0.2, kappa = 2, type = 1)
  fit_b <- distfit(x, type = 1, method = "bernstein", bernstein.m = 1)
  fit_m <- distfit(x, type = 1)

  expect_equal(fit_b$npar, 3)
  expect_equal(fit_b$loglik, fit_m$loglik, tolerance = 1e-8)
  expect_equal(fit_b$estimate, fit_m$estimate, tolerance = 1e-3)
  expect_equal(fit_b$sd, fit_m$sd, tolerance = 1e-2)
})

test_that("a stalled logit does not strand a beneficial Bernstein component", {
  set.seed(31337)
  x <- regpd(1000, kappa = 1.8, sigma = 2, xi = .15, type = 1)
  fit <- function(z) distfit(z, method = "bernstein", bernstein.m = 8, hessian = FALSE)
  a <- fit(x)
  scaled <- 1e-3*x
  b <- fit(scaled)
  expect_equal(b$convergence, 0)
  expect_lt(abs(b$loglik - length(x)*log(1000) - a$loglik), .001)
  z <- coef(b); w <- b$bernstein.weights
  # At fixed GPD parameters, moving a little mass toward any one component
  # cannot improve the fitted log likelihood materially.
  gain <- vapply(seq_along(w), function(j) {
    v <- .999*w; v[j] <- v[j]+.001
    bern_ll(scaled, z["sigma"], z["xi"], z["kappa"], v) - b$loglik
  }, numeric(1))
  expect_lt(max(gain), 1e-5)
})

test_that("distfit Bernstein S3 methods work", {
  set.seed(42)
  x <- regpd(300, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- distfit(x, type = 1, method = "bernstein", bernstein.m = 8, hessian = FALSE)

  ## coef
  expect_named(coef(fit), c("sigma", "xi", "kappa"))

  ## logLik
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")

  ## npar should account for Bernstein degrees of freedom
  ## 3 params (sigma, xi, kappa) + (m-1) Bernstein weights = 3 + 7 = 10
  expect_equal(fit$npar, 10)

  ## AIC
  expect_true(is.numeric(AIC(fit)))

  ## summary
  s <- summary(fit)
  expect_s3_class(s, "summary.distfit")
  expect_output(print(s))
})

test_that("distfit Bernstein plot does not error", {
  set.seed(42)
  x <- regpd(200, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  fit <- distfit(x, type = 1, method = "bernstein", bernstein.m = 6, hessian = FALSE)

  expect_no_error(plot(fit))
})

test_that("distfit Bernstein with different degrees converges", {
  set.seed(42)
  x <- regpd(300, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)

  fit4 <- distfit(x, type = 1, method = "bernstein", bernstein.m = 4, hessian = FALSE)
  fit12 <- distfit(x, type = 1, method = "bernstein", bernstein.m = 12, hessian = FALSE)

  expect_equal(fit4$convergence, 0)
  expect_equal(fit12$convergence, 0)
  expect_equal(fit4$bernstein.m, 4)
  expect_equal(fit12$bernstein.m, 12)
  ## the fitted density is a density
  for (f in list(fit4, fit12)) {
    dens <- function(z) egpd:::.bernstein_full_density(z, f$estimate[["sigma"]],
      f$estimate[["xi"]], f$estimate[["kappa"]], f$bernstein.weights, f$bernstein.m)
    expect_equal(integrate(dens, 0, Inf, rel.tol = 1e-10)$value, 1, tolerance = 1e-7)
  }
})

test_that("distfit Bernstein rejects non-egpd family", {
  x <- c(1, 2, 3, 4, 5)
  expect_error(
    distfit(x, type = 1, method = "bernstein", family = "degpd"),
    "Bernstein method is only available"
  )
})

test_that("AIC comparison between MLE and Bernstein works", {
  set.seed(42)
  x <- regpd(300, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)

  fit_mle <- distfit(x, type = 1, method = "mle")
  fit_bern <- distfit(x, type = 1, method = "bernstein", bernstein.m = 8, hessian = FALSE)

  ## uniform weights give b = 1, the type-1 EGPD from which the Bernstein fit starts,
  ## so its maximum cannot be lower; it has 3 + (m - 1) parameters
  expect_gte(fit_bern$loglik, fit_mle$loglik - 1e-8)
  expect_equal(AIC(fit_mle), -2 * fit_mle$loglik + 2 * 3)
  expect_equal(AIC(fit_bern), -2 * fit_bern$loglik + 2 * (3 + 8 - 1))
})
