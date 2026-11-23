## Round-4 count-area regressions (VB-02, VB-03, VB-04, VB-05, VB-11, VB-12, VB-13,
## VB-18, VB-19, VB-21, VB-22). Every expected value comes from an independent
## reference: a closed form, a brute-force sum, base R, or an independent recursion.

## ---- VB-13: zero-length n ------------------------------------------------------------

test_that("count and CompPareto generators give no draws for a zero-length n (VB-13)", {
  expect_length(stats::rpois(numeric(0), 2), 0L)        # base R
  rfns <- list(
    function(n) rgpois(n, mu = 2, lambda = 0.2),
    function(n) rgwaring(n, mu = 2, k = 2, rho = 3),
    function(n) rpoislnorm(n, mu = 2, sigma = 0.5),
    function(n) rprl(n, tau = 3),
    function(n) rcmp(n, mu = 2, nu = 1),
    function(n) rpig(n, mu = 2, sigma = 0.5),
    function(n) rzipig(n, mu = 2, sigma = 0.5, pi = 0.2),
    function(n) rbell(n, theta = 0.8),
    function(n) rzibell(n, theta = 0.8, pi = 0.2),
    function(n) rgpig(n, a = 0.5, b = 2, c = 0.6),
    function(n) rzigpig(n, a = 0.5, b = 2, c = 0.6, pi = 0.2),
    function(n) rBELL(n, mu = 2),
    function(n) rZIBELL(n, mu = 2, sigma = 0.2),
    function(n) rGPIG(n, mu = 2, sigma = 0.5, nu = 0.6),
    function(n) rZIGPIG(n, mu = 2, sigma = 0.5, nu = 0.6, tau = 0.2),
    function(n) rcomppareto(n, "exp", alpha = 1.5, theta = 2, rate = 1))
  for (rfn in rfns) expect_length(rfn(numeric(0)), 0L)
})

## ---- VB-11: rpoislnorm beyond the integer range ----------------------------------------

test_that("rpoislnorm returns doubles, not NA, above the integer range (VB-11)", {
  ## the documented mixture, built from base R: Z ~ N(log mu - sigma^2/2, sigma), Y ~ Pois(e^Z)
  set.seed(73)
  ref <- stats::rpois(3, exp(stats::rnorm(3, log(3e9) - 0.001^2 / 2, 0.001)))
  set.seed(73)
  expect_warning(got <- rpoislnorm(3, mu = 3e9, sigma = 0.001), NA)
  expect_identical(got, ref)
  expect_true(all(got > .Machine$integer.max))
  ## ordinary draws keep their values and integer type, as stats::rpois
  set.seed(1)
  ref <- stats::rpois(5, exp(stats::rnorm(5, log(3) - 0.5^2 / 2, 0.5)))
  set.seed(1)
  expect_identical(rpoislnorm(5, mu = 3, sigma = 0.5), ref)
})

## ---- VB-18: NA observation with an invalid parameter --------------------------------

test_that("an NA or NaN count with an invalid parameter is returned silently (VB-18)", {
  ## base R: the first argument's NA/NaN takes precedence, without a warning
  expect_warning(expect_identical(stats::dpois(NA_real_, -1), NA_real_), NA)
  expect_warning(expect_true(is.nan(stats::dpois(NaN, -1))), NA)
  isNA <- function(v) is.na(v) && !is.nan(v)
  expect_warning(expect_true(isNA(dpig(NA_real_, mu = -1))), NA)
  expect_warning(expect_true(isNA(dbell(NA_real_, theta = -1))), NA)
  expect_warning(expect_true(isNA(pgpois(NA_real_, mu = -1))), NA)
  expect_warning(expect_true(isNA(qcmp(NA_real_, mu = -1))), NA)
  expect_warning(expect_true(isNA(dzipig(NA_real_, mu = 1, pi = 2))), NA)
  expect_warning(expect_true(isNA(pgwaring(NA_real_, rho = 0.5))), NA)
  expect_warning(expect_true(is.nan(dgpig(NaN, a = 2))), NA)
  ## a non-NA count with the invalid parameter still gives NaN and the warning
  expect_warning(v <- dpig(c(NA, 1), mu = c(-1, -1)), "NaNs produced")
  expect_true(isNA(v[1]) && is.nan(v[2]))
})

## ---- VB-04: bamlss count families forward lower.tail, log.p and max.value -------------

test_that("bamlss count families honour lower.tail, log.p and max.value (VB-04)", {
  ## Bell(theta = 1), i.e. mu = e: P(Y = y) = e^(1 - e) B_y / y! (Castellares et al. 2018)
  B <- c(1, 1, 2, 5, 15, 52, 203)
  pb <- exp(1 - exp(1)) * B / factorial(0:6)
  par <- list(mu = exp(1))
  fb <- bell_bamlss()
  expect_equal(fb$p(0, par, lower.tail = FALSE), 1 - pb[1], tolerance = 1e-13)
  expect_equal(fb$p(0, par, log.p = TRUE), 1 - exp(1), tolerance = 1e-13)
  expect_equal(fb$p(4, par, lower.tail = FALSE, log.p = TRUE), log(1 - sum(pb[1:5])),
               tolerance = 1e-12)
  ## S(5) = 0.1226 > 0.1 >= S(6) = 0.0720
  expect_identical(fb$q(0.1, par, lower.tail = FALSE), 6)
  expect_identical(fb$q(log(0.1), par, lower.tail = FALSE, log.p = TRUE), 6)
  ## the 0.9 quantile is 6 (F(5) = 0.877 < 0.9 <= F(6) = 0.928): not attained within 5
  expect_warning(v <- fb$q(0.9, par, max.value = 5), "max.value")
  expect_true(is.na(v))
  fz <- zibell_bamlss()
  pz <- list(mu = exp(1), pi = 0.2)
  expect_equal(fz$p(0, pz, lower.tail = FALSE), 0.8 * (1 - pb[1]), tolerance = 1e-13)
  expect_equal(fz$p(2, pz, log.p = TRUE), log(0.2 + 0.8 * sum(pb[1:3])), tolerance = 1e-13)
  expect_identical(fz$q(0.1, pz, lower.tail = FALSE), 5)   # 0.8 S(5) = 0.098 <= 0.1
  ## GPIG from the Zhu-Joe pgf exp{b[(1-c)^a - (1-cs)^a]}: (n+1) p_{n+1} = sum g_k p_{n-k}
  gpig_pmf <- function(N, mu, a, c) {
    b <- mu * (1 - c)^(1 - a) / (a * c)
    k <- 0:N
    g <- exp(log(a) + log(b) + (k + 1) * log(c) + lgamma(k + 1 - a) - lgamma(1 - a) -
               lgamma(k + 1))
    p <- numeric(N + 1)
    p[1] <- exp(b * ((1 - c)^a - 1))
    for (n in 0:(N - 1)) p[n + 2] <- sum(g[1:(n + 1)] * p[(n + 1):1]) / (n + 1)
    p
  }
  pg <- gpig_pmf(10, 3, 0.4, 0.8)
  fg <- gpig_bamlss()
  pr <- list(mu = 3, sigma = 0.4, nu = 0.8)
  expect_equal(fg$p(2, pr, lower.tail = FALSE), 1 - sum(pg[1:3]), tolerance = 1e-12)
  expect_equal(fg$p(2, pr, log.p = TRUE), log(sum(pg[1:3])), tolerance = 1e-12)
  q <- fg$q(0.3, pr, lower.tail = FALSE)
  expect_identical(q, which(1 - cumsum(pg) <= 0.3)[1] - 1)
  fzg <- zigpig_bamlss()
  prz <- c(pr, list(pi = 0.25))
  expect_equal(fzg$p(2, prz, lower.tail = FALSE), 0.75 * (1 - sum(pg[1:3])), tolerance = 1e-12)
  expect_identical(fzg$q(0.3, prz, lower.tail = FALSE),
                   which(0.75 * (1 - cumsum(pg)) <= 0.3)[1] - 1)
})

## ---- VB-12: exact PIG and Bell samplers; max.value for GPIG ---------------------------

## Pearson chi-square of draws x against probabilities pr of the values 0..K-1 (and >= K).
count_gof_p <- function(x, pr) {
  K <- length(pr)
  obs <- tabulate(pmin(x, K) + 1, K + 1)
  e <- length(x) * c(pr, 1 - sum(pr))
  stats::pchisq(sum((obs - e)^2 / e), K, lower.tail = FALSE)
}

test_that("rpig and rzipig draw exactly from the PIG mixture, without a cap (VB-12)", {
  ## dpig's parameterisation: Poisson(mu Z) with Z inverse Gaussian, mean 1, variance sigma
  dig <- function(z, s) sqrt(1 / (2 * pi * s * z^3)) * exp(-(z - 1)^2 / (2 * s * z))
  ref <- vapply(0:6, function(y) stats::integrate(function(z) stats::dpois(y, 3 * z) * dig(z, 0.8),
                                                  0, Inf, rel.tol = 1e-12)$value, 0)
  expect_equal(dpig(0:6, mu = 3, sigma = 0.8), ref, tolerance = 1e-9)
  set.seed(20261117)
  x <- rpig(1e5, mu = 3, sigma = 0.8)
  expect_gt(count_gof_p(x, dpig(0:24, 3, 0.8)), 1e-3)
  ## mean mu and variance mu + sigma mu^2 (5 and 3 standard errors respectively)
  expect_lt(abs(mean(x) - 3), 5 * sqrt((3 + 0.8 * 9) / 1e5))
  expect_equal(stats::var(x), 3 + 0.8 * 9, tolerance = 0.03)
  set.seed(20261118)
  z <- rzipig(1e5, mu = 3, sigma = 0.8, pi = 0.3)
  expect_gt(count_gof_p(z, 0.3 * (0:24 == 0) + 0.7 * dpig(0:24, 3, 0.8)), 1e-3)
  ## means far beyond the old 1e5 inversion cap
  set.seed(1)
  expect_warning(big <- rpig(5, mu = 2e5, sigma = 1), NA)
  expect_true(all(is.finite(big)))
  set.seed(2)
  big <- rpig(2e4, mu = 2e5, sigma = 1)
  expect_lt(abs(mean(big) / 2e5 - 1), 5 * sqrt((2e5 + 4e10) / 2e4) / 2e5)
  ## NA and invalid parameters as before
  expect_warning(v <- rpig(3, mu = c(1, NA, -1)), "NaNs produced")
  expect_true(is.finite(v[1]) && is.na(v[2]) && !is.nan(v[2]) && is.nan(v[3]))
})

test_that("rbell and rzibell draw exactly (Neyman type A), without a cap (VB-12)", {
  ## Bell(theta): P(Y = y) = theta^y e^(1 - e^theta) B_y / y!, B_y from the Bell triangle
  bell_numbers <- function(n) {
    B <- numeric(n + 1); B[1] <- 1; row <- 1
    for (i in seq_len(n)) {
      new <- numeric(i + 1); new[1] <- row[i]
      for (j in seq_len(i)) new[j + 1] <- new[j] + row[j]
      row <- new; B[i + 1] <- row[1]
    }
    B
  }
  th <- 1.2
  pr <- exp(1 - exp(th) + (0:24) * log(th) + log(bell_numbers(24)) - lgamma(1:25))
  set.seed(20261119)
  x <- rbell(1e5, theta = th)
  expect_gt(count_gof_p(x, pr), 1e-3)
  ## mean theta e^theta and variance theta (1 + theta) e^theta
  set.seed(5)
  x <- rbell(2e5, theta = 2.5)
  expect_lt(abs(mean(x) - 2.5 * exp(2.5)), 5 * sqrt(2.5 * 3.5 * exp(2.5) / 2e5))
  expect_equal(stats::var(x), 2.5 * 3.5 * exp(2.5), tolerance = 0.02)
  set.seed(20261120)
  z <- rzibell(1e5, theta = th, pi = 0.25)
  expect_gt(count_gof_p(z, 0.25 * (0:24 == 0) + 0.75 * pr), 1e-3)
  ## theta = 10 has mean 10 e^10 = 220265, far beyond the old 1e5 inversion cap
  set.seed(1)
  expect_warning(big <- rbell(5, theta = 10), NA)
  expect_true(all(is.finite(big)))
  set.seed(4)
  big <- rbell(2e4, theta = 10)
  expect_lt(abs(mean(big) - 10 * exp(10)), 5 * sqrt(110 * exp(10) / 2e4))
  expect_warning(v <- rzibell(3, theta = c(1, NA, -1)), "NaNs produced")
  expect_true(is.finite(v[1]) && is.na(v[2]) && !is.nan(v[2]) && is.nan(v[3]))
})

test_that("rgpig and rzigpig pass max.value to the inversion (VB-12)", {
  ## the draws are the generalised inverse of uniforms, as before
  set.seed(9)
  u <- stats::runif(20)
  set.seed(9)
  expect_identical(rgpig(20, a = 0.4, b = 2, c = 0.97), qgpig(u, a = 0.4, b = 2, c = 0.97))
  ## a draw above max.value is NA with a warning
  set.seed(9)
  expect_warning(v <- rgpig(20, a = 0.4, b = 2, c = 0.97, max.value = 10), "max.value")
  ref <- qgpig(u, a = 0.4, b = 2, c = 0.97)
  expect_identical(is.na(v), ref > 10)
  expect_identical(v[!is.na(v)], ref[ref <= 10])
  set.seed(9)
  expect_warning(v <- rZIGPIG(10, mu = 50, sigma = 0.4, nu = 0.97, tau = 0.2, max.value = 20),
                 "max.value")
  expect_true(anyNA(v) && all(v[!is.na(v)] <= 20))
})

## ---- VB-22: fit_hurdle with the response formula not first ----------------------------

test_that("fit_hurdle finds the response on any element of a named list (VB-22)", {
  set.seed(1)
  d <- data.frame(x = stats::runif(100), y = stats::rpois(100, 2))
  f <- fit_hurdle(list(lkappa = ~1, lsigma = y ~ x, lxi = ~1), data = d, family = "degpd")
  expect_equal(deparse(f$zero$formula), ".positive ~ x")
  ## independent MLE: binomial GLM for P(Y > 0) plus the zero-truncated DEGPD (type 1),
  ## P(Y = y) = H(y + 1) - H(y), H(x) = (1 - (1 + xi x / sigma)^(-1/xi))^kappa,
  ## log sigma = b0 + b1 x, divided by 1 - P(Y = 0) = 1 - H(1)
  z <- stats::glm(I(y > 0) ~ x, family = stats::binomial, data = d)
  pos <- d[d$y > 0, ]
  lH <- function(x, s, xi, k) k * log1p(-(1 + xi * x / s)^(-1 / xi))
  nll <- function(th) {
    s <- exp(th[1] + th[2] * pos$x); xi <- exp(th[3]); k <- exp(th[4])
    l1 <- lH(pos$y + 1, s, xi, k); l0 <- lH(pos$y, s, xi, k); lz <- lH(1, s, xi, k)
    v <- -sum(l1 + log1p(-exp(l0 - l1)) - log1p(-exp(lz)))
    if (is.finite(v)) v else 1e10
  }
  best <- Inf
  for (st in list(c(log(2), 0, log(0.1), 0), c(0, 0, -3, 2))) {
    o <- stats::optim(st, nll, method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
    best <- min(best, stats::optim(o$par, nll, control = list(maxit = 2e4, reltol = 1e-15))$value)
  }
  expect_equal(as.numeric(logLik(f)), as.numeric(stats::logLik(z)) - best, tolerance = 1e-6)
  ## an unnamed list is positional, as in egpd(): the first element must carry the response
  expect_error(fit_hurdle(list(~1, y ~ x, ~1), data = d, family = "degpd"),
               "must include a response")
})

## ---- VB-03: hurdle quantiles on the log-survivor scale ---------------------------------

test_that("zero-truncated quantiles stay finite when P(Y > 0) is below 1e-16 (VB-03)", {
  zq <- egpd:::.zt_quantile
  ## xi = 0, kappa = 1, sigma = 0.01: geometric with ratio r = exp(-100), so
  ## P(Y > 1 | Y > 0) = r = 3.7e-44 and every quantile with u < 1 - r is 1
  expect_identical(zq("degpd", list(scale = rep(0.01, 3), shape = rep(0, 3), kappa = rep(1, 3)),
                      c(0.2, 0.8, 0.999)), c(1, 1, 1))
  ## kappa = 1: S(q) = (1 + xi (q + 1) / sigma)^(-1/xi); the smallest q with
  ## S(q) / S(0) <= 1 - u, for the double u = 1 - 1e-15 (1 - u = 9.992e-16), is 10991
  u <- 1 - 1e-15
  lS <- function(q) -log1p(0.2 * (q + 1) / 2) / 0.2
  expect_true(lS(10991) - lS(0) <= log1p(-u) && lS(10990) - lS(0) > log1p(-u))
  expect_identical(zq("degpd", list(scale = 2, shape = 0.2, kappa = 1), u), 10991)
  expect_identical(zq("degpd", list(scale = 2, shape = 0.2, kappa = 1), 1), Inf)
  ## generalized Waring (a, k, rho) = (mu (rho - 1) / k, 1.5, 3) with mu = 1e-10 (the
  ## fitting clamp): P(Y > 0) = 6.2e-11; quantiles from the closed-form pmf
  ## Gamma(a+rho) Gamma(k+rho) / (Gamma(rho) Gamma(a+k+rho)) (a)_y (k)_y / ((a+k+rho)_y y!)
  a <- 1e-10 * 2 / 1.5
  y <- 1:2e5
  p <- exp(lgamma(a + 3) + lgamma(4.5) - lgamma(3) - lgamma(a + 4.5) + lgamma(a + y) - lgamma(a) +
             lgamma(1.5 + y) - lgamma(1.5) - (lgamma(a + 4.5 + y) - lgamma(a + 4.5)) - lgamma(y + 1))
  S <- c(rev(cumsum(rev(p)))[-1], 0) / sum(p)            # S(q) / S(0), q = 1, 2, ...
  ref <- vapply(c(0.3, 0.9, 1 - 1e-6), function(u) which(S <= 1 - u)[1], 0)
  expect_identical(ref, c(1, 3, 209))
  expect_identical(zq("gwaring", list(mu = rep(1e-10, 3), k = rep(1.5, 3), rho = rep(3, 3)),
                      c(0.3, 0.9, 1 - 1e-6)), c(1, 3, 209))
})

test_that("hurdle predictions extrapolated to a degenerate base law are 1, not Inf (VB-03)", {
  set.seed(11)
  n <- 400
  x <- stats::runif(n)
  pos <- stats::rbinom(n, 1, 0.6)
  ypos <- ifelse(x > 0.6, 1, 1 + stats::rpois(n, 4 * (1 - x)))
  d <- data.frame(x = x, y = pos * ypos)
  f <- suppressWarnings(fit_hurdle(list(y ~ x, ~1, ~1), data = d, family = "degpd"))
  nd <- data.frame(x = c(2, 2.5, 3))
  pr <- predict(f, newdata = nd)
  ## closed form, kappa-power of the GPD: log S(q) = log(1 - H(q + 1)),
  ## H(x) = (1 - t)^kappa, log t = -log1p(xi x / sigma) / xi
  lS <- function(q) {
    lt <- -log1p(pr$shape * (q + 1) / pr$scale) / pr$shape
    log(-expm1(pr$kappa * log1p(-exp(lt))))
  }
  expect_true(all(lS(1) - lS(0) < log1p(-0.9)))           # P(Y > 1 | Y > 0) < 0.1
  q <- predict(f, newdata = nd, prob = c(0.5, 0.9))
  expect_true(all(pr$pzero < 0.5))
  expect_identical(unname(as.matrix(q)), matrix(1, 3, 2))
})

## ---- VB-21: spliced tail fit ------------------------------------------------------------

test_that("the spliced tail fit stops at zero likelihood and reports non-convergence (VB-21)", {
  skip_if_not_installed("qgam")
  set.seed(42)
  x <- stats::runif(200, -1, 1)
  y <- stats::rpois(200, exp(1.2 + 0.5 * x)) +
    stats::rbinom(200, 1, 0.3) * rdiscegpd(200, sigma = 8, xi = 0.3, kappa = 1, type = 1)
  df <- data.frame(y = y, x = x)
  fit <- function(...) fit_spliced_degpd(y ~ x, data = df, bulk_probs = c(0.5, 0.8),
                                          threshold_prob = 0.8, ...)
  ## xi = -2 with sigma = exp(-5) ends the support at sigma / 2 < 1: every exceedance >= 1
  ## has probability zero, so no fit exists from that start
  expect_error(fit(egpd_control = list(start = list(beta = -5, xi = -2))),
               "positive likelihood")
  expect_warning(fb <- fit(egpd_control = list(optim.method = "Nelder-Mead",
                                               control = list(maxit = 2))),
                 class = "egpd_nonconvergence")
  expect_identical(fb$tail_fit$convergence, 1L)
  expect_true(is.na(fb$tail_fit$loglik))
  expect_warning(fd <- fit(), NA)
  expect_identical(fd$tail_fit$convergence, 0L)
  expect_true(is.finite(fd$tail_fit$loglik))
})

## ---- VB-19: CompPareto splice weight on the log scale -----------------------------------

test_that("CompPareto stays finite when g(theta) or G(theta) underflow (VB-19)", {
  ## documented law: f = g(x) / ((1 + phi) G(theta)) for x <= theta, else
  ## phi / (1 + phi) alpha (2 theta)^alpha / (x + theta)^(alpha + 1), with
  ## log phi = log(2 theta / alpha) + log g(theta) - log G(theta) (here from dlnorm/plnorm logs)
  lphi <- log(2) + stats::dlnorm(1, 10, 0.2, log = TRUE) - log(2) -
    stats::plnorm(1, 10, 0.2, log.p = TRUE)
  lf <- lphi - log1p(exp(lphi)) + log(2) + 2 * log(2) - 3 * log(20001)
  expect_warning(v <- dcomppareto(20000, "lnorm", alpha = 2, theta = 1, meanlog = 10, sdlog = 0.2),
                 NA)
  expect_equal(v, exp(lf), tolerance = 1e-12)
  expect_equal(dcomppareto(20000, "lnorm", alpha = 2, theta = 1, meanlog = 10, sdlog = 0.2,
                           log = TRUE), lf, tolerance = 1e-12)
  lS <- lphi - log1p(exp(lphi)) + 2 * (log(2) - log(20001))
  expect_equal(pcomppareto(20000, "lnorm", alpha = 2, theta = 1, meanlog = 10, sdlog = 0.2,
                           lower.tail = FALSE, log.p = TRUE), lS, tolerance = 1e-12)
  expect_equal(qcomppareto(lS, "lnorm", alpha = 2, theta = 1, meanlog = 10, sdlog = 0.2,
                           lower.tail = FALSE, log.p = TRUE), 20000, tolerance = 1e-10)
  ## Weibull body with the splice point far in its tail: log phi = -892.5, so phi
  ## underflows and f(0.5) = g(0.5) / G(30)
  ref <- stats::dweibull(0.5, 2, 1) / stats::pweibull(30, 2, 1)
  expect_warning(v <- dcomppareto(0.5, "weibull", alpha = 2, theta = 30, shape = 2, scale = 1), NA)
  expect_equal(v, ref, tolerance = 1e-13)
  expect_equal(pcomppareto(0.5, "weibull", alpha = 2, theta = 30, shape = 2, scale = 1),
               stats::pweibull(0.5, 2, 1), tolerance = 1e-13)
  expect_equal(qcomppareto(stats::pweibull(0.5, 2, 1), "weibull", alpha = 2, theta = 30,
                           shape = 2, scale = 1), 0.5, tolerance = 1e-12)
  ## exponential body, rate 10, theta = 80 (g(theta) = 10 e^-800 underflows): the GAM
  ## objective is -log f(1) = -(log 10 - 10), as at theta = 40
  eta <- matrix(c(log(10), log(2), log(80)), nrow = 1)
  expect_equal(egpd:::.comppareto_nll_obs(1, eta, "exp"), 10 - log(10), tolerance = 1e-13)
})

## ---- VB-05: GPIG()/ZIGPIG() gamlss starts ----------------------------------------------

test_that("gamlss with GPIG() and ZIGPIG() reaches the MLE on heavy-tailed data (VB-05)", {
  skip_if_not_installed("gamlss")
  skip_on_cran()
  ## independent GPIG pmf from the Zhu-Joe pgf exp{b[(1-c)^a - (1-cs)^a]}
  gpig_pmf <- function(N, mu, a, c) {
    b <- mu * (1 - c)^(1 - a) / (a * c)
    k <- 0:N
    g <- exp(log(a) + log(b) + (k + 1) * log(c) + lgamma(k + 1 - a) - lgamma(1 - a) -
               lgamma(k + 1))
    p <- numeric(N + 1)
    p[1] <- exp(b * ((1 - c)^a - 1))
    for (n in 0:(N - 1)) p[n + 2] <- sum(g[1:(n + 1)] * p[(n + 1):1]) / (n + 1)
    p
  }
  ctrl <- gamlss::gamlss.control(n.cyc = 20, trace = FALSE)
  ## maximum log-likelihoods of these samples (optim on the recursion above; round-4 verifier)
  set.seed(101)
  y <- rGPIG(300, 6, 0.4, 0.95)
  m <- gamlss::gamlss(y ~ 1, family = GPIG(), control = ctrl)
  expect_true(m$converged)
  p <- gpig_pmf(max(y), fitted(m, "mu")[1], fitted(m, "sigma")[1], fitted(m, "nu")[1])
  expect_lt(abs(sum(log(p[y + 1])) - -863.899310634678), 1e-3)
  set.seed(31)
  y <- rZIGPIG(300, 6, 0.4, 0.95, 0.2)
  m <- gamlss::gamlss(y ~ 1, family = ZIGPIG(), control = ctrl)
  expect_true(m$converged)
  p <- gpig_pmf(max(y), fitted(m, "mu")[1], fitted(m, "sigma")[1], fitted(m, "nu")[1])
  pz <- fitted(m, "tau")[1]
  expect_lt(abs(sum(log(pz * (y == 0) + (1 - pz) * p[y + 1])) - -769.602838887151), 1e-3)
  set.seed(5)
  x <- stats::runif(300)
  y <- rGPIG(300, mu = exp(1.2 + 0.8 * x), sigma = 0.4, nu = 0.95)
  m <- gamlss::gamlss(y ~ x, family = GPIG(), control = gamlss::gamlss.control(n.cyc = 30, trace = FALSE))
  expect_true(m$converged)
  expect_lt(abs(as.numeric(logLik(m)) - -786.336419213823), 5e-3)
})

## ---- VB-02: CMP grid validated against the exact law -----------------------------------

test_that("the CMP likelihood is exact for strongly underdispersed laws (VB-02)", {
  X <- matrix(1, 1, 1)
  off <- list(numeric(0), numeric(0))
  lp <- function(y, mu, nu)
    vapply(y, function(v) -egpd:::cmp1d0(list(log(mu), log(nu)), X, X, v, 0L, 0L, off), 0)
  ## independent enumeration (Shmueli et al. 2005): lambda^j / (j!)^nu / Z with E[Y] = mu
  cmp_ref <- function(mu, nu, J) {
    j <- 0:J
    lf <- function(ll) j * ll - nu * lgamma(j + 1)
    mean_at <- function(ll) { f <- lf(ll); w <- exp(f - max(f)); sum(j * w) / sum(w) }
    ll <- stats::uniroot(function(l) mean_at(l) - mu, c(-50, nu * log(mu + 2) + 5),
                         tol = 1e-13)$root
    f <- lf(ll)
    f - (max(f) + log(sum(exp(f - max(f)))))
  }
  ## log p(100 | mu = 100.2, nu = 1000), from a 200-bit Rmpfr enumeration
  expect_equal(lp(100, 100.2, 1000), -0.22352145341103595, tolerance = 1e-9)
  for (cs in list(c(100.2, 1000, 250), c(3, 32, 40), c(10, 100, 40), c(3, 100, 40))) {
    ref <- cmp_ref(cs[1], cs[2], cs[3])
    got <- lp(0:cs[3], cs[1], cs[2])
    expect_true(all(is.finite(got)))
    expect_equal(sum(exp(got)), 1, tolerance = 1e-8)
    bulk <- ref > log(1e-10)
    expect_lt(max(abs(got[bulk] - ref[bulk])), 1e-7)
  }
  ## at fixed nu, CMP is an exponential family in log lambda, so the MLE has mean
  ## ybar = 100.2; the log-likelihood there is -50.055762133188662 (200-bit Rmpfr)
  d <- data.frame(y = c(rep(100, 80), rep(101, 20)))
  f <- egpd(list(mu = y ~ 1, nu = ~1), data = d, family = "cmp", fixed = list(nu = 1000),
            inits = c(log(100), log(1000)))
  p <- predict(f, newdata = d[1, , drop = FALSE], type = "response")
  expect_equal(p$mu, 100.2, tolerance = 1e-7)
  expect_equal(as.numeric(logLik(f)), -50.055762133188662, tolerance = 1e-9)
  ## free nu: the independent MLE of these counts is mu = 100.055 (= ybar), nu = 279.48,
  ## log-likelihood -179.6556851 (optim on the enumeration above)
  d <- data.frame(y = c(rep(99, 28), rep(100, 135), rep(101, 35), rep(102, 2)))
  f <- egpd(list(mu = y ~ 1, nu = ~1), data = d, family = "cmp")
  p <- predict(f, newdata = d[1, , drop = FALSE], type = "response")
  expect_equal(p$mu, 100.055, tolerance = 1e-6)
  expect_equal(p$nu, 279.4823, tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f)), -179.6556851015, tolerance = 1e-8)
})
