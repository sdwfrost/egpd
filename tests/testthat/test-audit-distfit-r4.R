## Round-4 regressions for distfit(), its methods, and the compound/aggregate EGPD
## (verify-C clusters VC-04 to VC-31). Expected values are closed forms, finite
## sums or independent implementations, never the package compared with itself.

## ---- VC-04: PIG at the Poisson limit ----------------------------------------

test_that("VC-04: distfit(family = 'pig') attains the Poisson-limit maximum", {
  set.seed(5)
  x <- rpois(25, 2)
  ## The sample is underdispersed, so the PIG likelihood is maximised at the
  ## Poisson limit sigma -> 0, whose maximum is the closed form at mu = mean(x).
  ll_pois <- sum(stats::dpois(x, mean(x), log = TRUE))
  expect_equal(ll_pois, -41.39568, tolerance = 1e-6)
  f <- distfit(x, family = "pig")
  expect_lte(f$loglik, 0)
  expect_lt(abs(f$loglik - ll_pois), 1e-3)
  expect_equal(f$estimate[["mu"]], mean(x), tolerance = 1e-3)
  expect_lt(f$estimate[["sigma"]], 1e-3)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(f))
})

test_that("VC-04: a positive discrete log-likelihood is an error, not a fit", {
  ## A pmf is at most one; a likelihood that says otherwise is numerically broken.
  local_mocked_bindings(dbell = function(x, theta, log = FALSE)
    rep(if (log) log(2) else 2, length(x)))
  expect_error(distfit(c(0, 1, 2, 3, 1), family = "bell"), "positive")
})

## ---- VC-07: exact zeros in the continuous EGPD ------------------------------

test_that("VC-07: distfit(family = 'egpd') rejects exact zeros (MLE and Bernstein)", {
  ## f(0) = kappa H(0)^(kappa - 1) h(0) with H(0) = 0 is infinite for kappa < 1 and
  ## 0 for kappa > 1, so one zero used to pin kappa-hat = 1 with SE 1e-13.
  set.seed(77)
  x <- c(0, rexp(100, 0.7))
  expect_error(distfit(x, family = "egpd"), "positive data.*ziegpd")
  expect_error(distfit(x, family = "egpd", method = "bernstein"), "positive data")
  ## zeros remain valid where the law has mass at 0
  expect_s3_class(suppressWarnings(distfit(x, family = "ziegpd")), "distfit")
})

## ---- VC-05 / VC-06: Bernstein likelihood and optimiser ----------------------

## Independent Bernstein-EGPD log-likelihood, from base R only:
## f(x) = sum_j w_j dbeta(H^kappa; j, m - j + 1) kappa H^(kappa - 1) h(x).
bern_ll_ref <- function(x, sigma, xi, kappa, w) {
  m <- length(w)
  t <- 1 + xi * x / sigma
  logS <- if (xi == 0) -x / sigma else -log(t) / xi
  H <- -expm1(logS)
  logh <- -log(sigma) + logS - log(t)
  b <- rowSums(vapply(seq_len(m), function(j) w[j] * stats::dbeta(H^kappa, j, m - j + 1),
                      numeric(length(x))))
  sum(log(b) + log(kappa) + (kappa - 1) * log(H) + logh)
}

test_that("VC-05: the Bernstein fit is not stopped short by Nelder-Mead", {
  skip_on_cran()
  set.seed(101)
  x <- regpd(500, sigma = 2, xi = 0.15, kappa = 1.4, type = 1)
  fb <- distfit(x, method = "bernstein", bernstein.m = 8)
  expect_equal(fb$convergence, 0L)
  ## the reported value is the likelihood at the reported estimate ...
  est <- fb$estimate
  expect_equal(fb$loglik, bern_ll_ref(x, est[["sigma"]], est[["xi"]], est[["kappa"]],
                                      fb$bernstein.weights), tolerance = 1e-10)
  ## ... and it reaches the maximum found by repeated independent optimisation of
  ## bern_ll_ref (-1036.954), not the -1038.293 where Nelder-Mead used to stop
  expect_gte(fb$loglik, -1037.06)
})

test_that("VC-06: Bernstein likelihood and density are not clipped at tiny x", {
  ## m = 1 (b = 1), sigma = 1, xi = 0 fixed: f(x) = kappa H^(kappa - 1) h with
  ## H = 1 - exp(-x), so kappa-hat = -n / sum(log(H)) in closed form.
  set.seed(1)
  x <- c(rexp(50), 1e-20)
  f <- distfit(x, method = "bernstein", bernstein.m = 1,
               fix.arg = list(sigma = 1, xi = 0))
  expect_equal(f$estimate[["kappa"]], -length(x) / sum(log(-expm1(-x))),
               tolerance = 1e-5)
  ## density at x = 1e-20 for kappa = 1/2: (1/2) H^(-1/2) h = 0.5e10
  expect_equal(egpd:::.bernstein_full_density(1e-20, sigma = 1, xi = 0, kappa = 0.5,
                                              weights = 1, m = 1),
               0.5 * (-expm1(-1e-20))^(-0.5) * exp(-1e-20), tolerance = 1e-12)
})

## ---- VC-08 / VC-09: confint.distfit() ---------------------------------------

test_that("VC-08: the generalised Waring rho interval stays above 1", {
  set.seed(2)
  x <- rgwaring(150, mu = 4, k = 3, rho = 1.6)
  f <- distfit(x, family = "gwaring")
  ci <- confint(f)
  r <- f$estimate[["rho"]]
  s <- f$sd[["rho"]]
  ## Wald interval for log(rho - 1), whose delta-method SE is sd(rho) / (rho - 1)
  expect_equal(unname(ci["rho", ]),
               1 + exp(log(r - 1) + c(-1, 1) * stats::qnorm(0.975) * s / (r - 1)),
               tolerance = 1e-12)
  expect_equal(ci["rho", 1], 1.067139, tolerance = 1e-5)   # was 0.2355
  expect_true(all(is.finite(dgwaring(0:2, mu = f$estimate[["mu"]], k = f$estimate[["k"]],
                                     rho = ci["rho", 1]))))
})

test_that("VC-09: confint.distfit() takes indices and rejects unknown parameters", {
  set.seed(5)
  f <- distfit(2 * rexp(150), type = 1)
  ci <- confint(f)
  expect_identical(confint(f, parm = 2), ci["xi", , drop = FALSE])
  expect_identical(confint(f, parm = c(3, 1)), ci[c("kappa", "sigma"), ])
  expect_error(confint(f, parm = "foo"), "'foo' not found")
  expect_error(confint(f, parm = c("xi", "foo")), "'foo' not found")
  expect_error(confint(f, parm = 4), "out of range")
  ## neural posterior fits: equal-tailed posterior quantiles, any number of rows
  post <- rbind(a = stats::qnorm(ppoints(400)), b = stats::qexp(ppoints(400)))
  npe <- structure(list(family = "begpd", estimator_type = "npe",
                        estimate = c(a = 0, b = log(2)), posterior_samples = post),
                   class = c("distfit_npe", "distfit"))
  expect_equal(unname(confint(npe, parm = "b")),
               matrix(stats::quantile(post["b", ], c(0.025, 0.975), names = FALSE), 1))
  expect_equal(confint(npe, parm = 1), confint(npe)["a", , drop = FALSE])
  expect_error(confint(npe, parm = "c"), "'c' not found")
})

## ---- VC-10: fix.arg / start validation --------------------------------------

test_that("VC-10: unnamed or invalid fix.arg/start are errors, not ignored", {
  set.seed(5)
  x <- 2 * rexp(100)
  expect_error(distfit(x, fix.arg = list(1.5)), "named list")
  expect_error(distfit(x, start = list(2, 0.1, 1.5)), "named list")
  expect_error(distfit(x, start = list(sigma = -1)), "positive")
  expect_error(distfit(x, fix.arg = list(kappa = 1.5, kappa = 2)), "duplicated")
  expect_error(distfit(rpois(50, 3), family = "gwaring", start = list(rho = 0.9)),
               "above 1")
  expect_error(distfit(x, method = "bernstein", start = list(2, 0.1)), "named")
  ## a named fixed value is honoured
  expect_named(distfit(x, fix.arg = list(kappa = 1))$estimate, c("sigma", "xi"))
})

## ---- VC-11: plot.distfit() PMF panel ----------------------------------------

test_that("VC-11: the fitted pmf is drawn at the bar centres", {
  set.seed(1)
  x <- c(rpois(150, 3), 15, 25)
  f <- distfit(x, family = "bell")
  vals <- sort(unique(x))
  pts <- NULL
  local_mocked_bindings(points = function(x, y = NULL, ...) {
    if (is.null(pts)) pts <<- list(x = as.numeric(x), y = y)
    invisible(NULL)
  })
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(f)
  ## barplot() with its default width 1 and space 0.2 centres bar k at 0.7 + 1.2 (k - 1)
  expect_equal(pts$x, 0.7 + 1.2 * (seq_along(vals) - 1))
  expect_equal(pts$y, dbell(vals, theta = f$estimate[["theta"]]))
})

## ---- VC-20: Julia set-up ----------------------------------------------------

test_that("VC-20: a failed Julia package load gives the full install command", {
  ## `using BSON` fails unless BSON is a direct dependency of the Julia environment
  local_mocked_bindings(
    .check_julia_deps = function() invisible(NULL),
    .julia_eval = function(...)
      stop("ArgumentError: Package BSON not found in current path.")
  )
  env <- egpd:::.begpd_env
  old <- env$julia_initialized
  env$julia_initialized <- FALSE
  on.exit(env$julia_initialized <- old, add = TRUE)
  expect_error(egpd:::.init_julia_begpd(),
               'BSON not found.*Pkg.add\\(\\["NeuralEstimators", "Flux", "BSON", "Optimisers"\\]\\)')
})

## ---- VC-15: npe.control ranges -----------------------------------------------

test_that("VC-15: npe.control values are range-checked", {
  lr_fn <- function(theta) -0.5 * colSums(((theta - c(0.5, 0.05, -0.5)) / 0.2)^2)
  is_ctl <- function(ctl) egpd:::.npe_importance_sample(lr_fn, "mdgpd", 2L, 200L,
                                                         control = ctl)
  expect_error(is_ctl(list(defensive = 1)), "defensive.*\\[0, 1\\)")
  expect_error(is_ctl(list(defensive = 1.5)), "defensive")
  expect_error(is_ctl(list(df = 0)), "df")
  expect_error(is_ctl(list(inflate = -1)), "inflate")
  expect_error(is_ctl(list(max_stages = -1)), "max_stages")
  expect_error(is_ctl(list(ess_floor = 0)), "ess_floor")
  ## a defensive weight that rounds to the whole stage draws from the prior: the
  ## proposal density stays defined (it was NaN, so every later stage had ESS 0)
  set.seed(1)
  r <- is_ctl(list(defensive = 0.9995, K = 500L, max_stages = 2L))
  expect_true(all(is.finite(r$diagnostics$path$ess) & r$diagnostics$path$ess > 0))
})

## ---- VC-16 / VC-17: compound Poisson-EGPD tails and tiny rates --------------

test_that("VC-16: cpegpd bin masses keep their precision where F is close to 1", {
  ## Exp(1) amounts (type 1, kappa = 1, xi = 0) on h = 1: each event fills a
  ## positive geometric number of bins, f_j = (1 - r) r^(j - 1) with r = exp(-1), so
  ## P(S = k) = sum_n dpois(n, lambda) choose(k - 1, n - 1) (1 - r)^n r^(k - n).
  r <- exp(-1)
  exact_log <- function(k, lambda) {
    n <- seq_len(k)
    a <- stats::dpois(n, lambda, log = TRUE) + lchoose(k - 1, n - 1) +
      n * log1p(-r) + (k - n) * log(r)
    max(a) + log(sum(exp(a - max(a))))
  }
  ks <- c(10, 36, 40, 50, 70)
  for (lambda in c(0.01, 0.5)) {
    got <- dcpegpd(ks, lambda = lambda, sigma = 1, xi = 0, kappa = 1, h = 1,
                   K = 120, log = TRUE)
    expect_equal(got, vapply(ks, exact_log, numeric(1), lambda = lambda),
                 tolerance = 1e-10)
  }
  expect_equal(dcpegpd(50, 0.01, sigma = 1, xi = 0, kappa = 1, h = 1, K = 120,
                       log = TRUE), -53.6798, tolerance = 1e-6)    # was -55.2712
})

test_that("VC-17: positive-part and p = 1 quantiles survive a rate below 1e-13", {
  ## As lambda -> 0 the positive part is one gauge-rounded event:
  ## h ceiling(Q(p) / h) with the GPD quantile Q(p) = ((1 - p)^-xi - 1) / xi
  p <- c(0.5, 0.9, 0.99)
  single <- 0.2 * ceiling(((1 - p)^(-0.2) - 1) / 0.2 / 0.2)
  expect_equal(single, c(0.8, 3.0, 7.6))
  for (lambda in c(1e-14, 1e-17)) {
    q <- egpd:::.aggregate_positive_quantile(p, "cpegpd", sigma = 1, kappa = 1,
                                            xi = 0.2, lambda = lambda,
                                            cpegpd.h = 0.2, K = 2000)
    expect_equal(q, single)
  }
  ## discrete: Y | Y >= 1 for Y = floor(X), P(Y >= y) = (1 + 0.2 y)^-5
  sv <- function(y) (1 + 0.2 * y)^(-5)
  single_d <- vapply(p, function(pp) {
    y <- 1
    while (1 - sv(y + 1) / sv(1) < pp) y <- y + 1
    y
  }, numeric(1))
  expect_equal(egpd:::.aggregate_positive_quantile(p, "cpdegpd", sigma = 1, kappa = 1,
                                                  xi = 0.2, lambda = 1e-17, K = 2000),
               single_d)
  ## the support is unbounded for every lambda > 0, even when exp(-lambda) == 1
  expect_identical(qcpegpd(1, lambda = 1e-17, sigma = 1, xi = 0.2, kappa = 1, K = 200), Inf)
  expect_identical(qcpdegpd(1, lambda = 1e-17, sigma = 1, xi = 0.2, kappa = 1, K = 200), Inf)
})
