## Round-4 fixes to the outer (smoothing-parameter) optimisers and the REML criterion:
## VA-01 (unavailable third derivatives at trial points), VA-03 (REML gradient under
## the curvature sandwich adjustment), VA-10 (fits ending on the xi < -1 likelihood
## singularity), VA-14 (multi-penalty log|S|+) and VA-16 (negREML constants).

## ---- VA-01: egpd_d34_unavailable at a trial point reaches .outer() ------------------

d34_condition <- function()
  structure(class = c("egpd_d34_unavailable", "error", "condition"),
            list(message = "analytic DEGPD d34 is unstable near zero shape", call = NULL))

test_that("VA-01: .BFGS passes egpd_d34_unavailable at a trial point on to .outer()", {
  ctl <- list(steptol = 1e-12, itlim = 100, fntol = 1e-14, gradtol = 1e-8, stepmax = 3)
  f <- function(p) sum((p - c(1, -2))^2)
  calls <- 0
  g <- function(p) {
    calls <<- calls + 1
    if (calls == 2) stop(d34_condition())
    2 * (p - c(1, -2))
  }
  expect_error(egpd:::.BFGS(c(0, 0), f, g, control = ctl), class = "egpd_d34_unavailable")
  ## any other failure at a trial point still only rejects that point
  calls <- 0
  g2 <- function(p) {
    calls <<- calls + 1
    if (calls == 2) stop("refit failed")
    2 * (p - c(1, -2))
  }
  fit <- egpd:::.BFGS(c(0, 0), f, g2, control = ctl)
  expect_equal(fit$par, c(1, -2), tolerance = 1e-6)
})

test_that("VA-01: the Newton outer search passes egpd_d34_unavailable on as well", {
  ctl <- list(steptol = 1e-12, itlim = 100, fntol = 1e-14, gradtol = 1e-8, stepmax = 3)
  f <- function(p) sum((p - c(1, -2))^2)
  calls <- 0
  search <- function(fail) function(p, kept) {
    calls <<- calls + 1
    if (calls == 2) fail()
    g <- 2 * (p - c(1, -2))
    structure(g / 2, gradient = g, Hessian = diag(2, 2))
  }
  expect_error(egpd:::.newton_step_inner(c(0, 0), f, search(function() stop(d34_condition())),
                                         control = ctl),
               class = "egpd_d34_unavailable")
  calls <- 0
  fit <- egpd:::.newton_step_inner(c(0, 0), f, search(function() stop("refit failed")),
                                   control = ctl)
  expect_equal(fit$pars, c(1, -2), tolerance = 1e-6)
})

## DEGPD with the identity shape link and true xi = 0.06, just outside the |xi| < 0.05
## band where the analytic d34 is unavailable (verifier probe v01_bfgs_d34.R).
va01_data <- function() {
  set.seed(5)
  n <- 600
  x <- runif(n)
  runif(n)
  u <- runif(n)^(1 / 1.5)
  data.frame(yd = floor(exp(1 + sin(2 * pi * x)) / 0.06 * ((1 - u)^(-0.06) - 1)), x = x)
}

test_that("VA-01: an analytic outer search that reaches the d34 band retries with finite differences", {
  d <- va01_data()
  args <- list(formula = list(yd ~ s(x, k = 8), ~1, ~1), data = d, family = "degpd",
               degpd.args = list(link = "identity"), restarts = FALSE)
  fd <- do.call(egpd, c(args, list(outer = "fd")))
  for (outer in c("bfgs", "newton")) {
    expect_warning(fit <- do.call(egpd, c(args, list(outer = outer))),
                   "finite-difference outer optimization")
    expect_equal(fit$convergence, 0L)
    ## The minimiser of an independent Laplace REML criterion (own identity-link DEGPD1
    ## pmf, own Newton inner fit, eigenvalue log-determinants) is sp = 0.08211, with
    ## xi = 0.03562 there; before the fix both searches stopped on the band edge,
    ## at sp 0.5949 and xi 0.0500, 7.96 REML units above that minimum.
    expect_equal(unname(fit$sp), 0.08211, tolerance = 5e-3)
    expect_equal(unname(fit$coefficients["shape_(Intercept)"]), 0.03562, tolerance = 2e-3)
    expect_lt(abs(fit$negREML - fd$negREML), 1e-3)
  }
})

test_that("VA-01: identity-link fits stopping on the d34 band edge now match finite differences", {
  ## verifier probe v01d_freq.R, true xi = 0.06: these seeds stopped on the band edge
  ## (xi = 0.050) 3.81 and 1.04 REML units above the finite-difference optimum
  for (seed in c(2, 12)) {
    set.seed(seed)
    n <- 400
    x <- runif(n)
    u <- runif(n)^(1 / 1.5)
    d <- data.frame(yd = floor(exp(1 + sin(2 * pi * x)) / 0.06 * ((1 - u)^(-0.06) - 1)), x = x)
    args <- list(formula = list(yd ~ s(x, k = 8), ~1, ~1), data = d, family = "degpd",
                 degpd.args = list(link = "identity"), restarts = FALSE)
    expect_warning(fit <- do.call(egpd, args), "finite-difference outer optimization")
    fd <- do.call(egpd, c(args, list(outer = "fd")))
    expect_equal(unname(fit$sp), unname(fd$sp), tolerance = 1e-3)
    expect_lt(abs(fit$negREML - fd$negREML), 1e-4)
    expect_lt(abs(fit$coefficients["shape_(Intercept)"] - 0.05), 0.04)
    expect_gt(abs(fit$coefficients["shape_(Intercept)"] - 0.05), 1e-3)
  }
})

test_that("VA-01: the fully default DEGPD call no longer stops on the guard band", {
  skip_on_cran()
  ## verifier probe v01b_default.R: log shape link, smooths on scale and shape, default
  ## outer and restarts. Before the fix the fit stopped with min xi = 0.001 (the band
  ## edge), sp (0.458, 1.006), 3.93 REML units above the finite-difference fit and
  ## 3.98 units above it by an independent LAML.
  set.seed(6)
  n <- 800
  x <- runif(n)
  x2 <- runif(n)
  xi <- 0.35 * exp(-7 * x2)
  u <- runif(n)^(1 / 1.5)
  d <- data.frame(yd = floor(exp(1 + sin(2 * pi * x)) / xi * ((1 - u)^(-xi) - 1)), x = x, x2 = x2)
  form <- list(yd ~ s(x, k = 8), ~ s(x2, k = 6), ~1)
  expect_warning(fit <- egpd(form, data = d, family = "degpd"),
                 "finite-difference outer optimization")
  fd <- suppressWarnings(egpd(form, data = d, family = "degpd", outer = "fd"))
  expect_equal(unname(fit$sp), unname(fd$sp), tolerance = 1e-3)
  expect_lt(abs(fit$negREML - fd$negREML), 1e-3)
})

## ---- VA-03: REML gradient under the curvature sandwich adjustment -------------------

## EGPD m1 with a smooth on the scale, 40 clusters of 6 with a cluster effect on the log
## scale (verifier probe v03_curv.R); the curvature adjustment's CH is far from diagonal.
va03_data <- function() {
  set.seed(31)
  ncl <- 40
  n <- ncl * 6
  id <- rep(seq_len(ncl), each = 6)
  x <- runif(n)
  ue <- rnorm(ncl, sd = 0.4)[id]
  u <- runif(n)^(1 / 1.5)
  data.frame(y = exp(0.5 + sin(2 * pi * x) + ue) / 0.2 * ((1 - u)^(-0.2) - 1), x = x, id = id)
}

test_that("VA-03: the analytic REML derivatives are right under the curvature sandwich adjustment", {
  skip_if_not_installed("numDeriv")
  fit <- egpd(list(y ~ s(x, k = 6), ~1, ~1), data = va03_data(), family = "egpd",
              sandwich.args = list(id = "id", method = "curvature"))
  expect_gt(max(abs(fit$likdata$CH - diag(diag(fit$likdata$CH)))), 0.5)
  ## minimiser of an independent adjusted LAML (own EGPD1 density at
  ## compmode + CH (beta - compmode), own Newton inner fit); before: sp 0.2050
  expect_equal(unname(fit$sp), 0.197551, tolerance = 5e-3)
  ld <- egpd:::.tighten_inner(fit$likdata, 1e-12)
  rho <- log(0.197303) + 0.7
  b <- attr(egpd:::.reml0(rho, fit$likfns, ld, fit$Sdata, beta = fit$coefficients), "beta")
  crit <- function(r) as.numeric(egpd:::.reml0(r, fit$likfns, ld, fit$Sdata, beta = b))
  g <- egpd:::.reml1(rho, fit$likfns, ld, fit$Sdata, beta = b)
  ## Richardson derivative of the independent LAML: 0.7569229 (before: 0.6910398)
  expect_equal(g, 0.7569229, tolerance = 1e-5)
  expect_equal(g, numDeriv::grad(crit, rho), tolerance = 1e-6)
  ## the REML Hessian (differences of the analytic gradient) against second differences
  ## of the criterion: 0.998829 (before: 0.963839)
  h <- egpd:::.reml12(structure(rho, beta = b), fit$likfns, ld, fit$Sdata)[[2]]
  expect_equal(c(h), numDeriv::hessian(crit, rho)[1, 1], tolerance = 1e-4)
})

## ---- VA-14: multi-penalty log|S|+ with widely spread smoothing parameters -----------

test_that("VA-14: log|S|+ keeps the directions penalised only by small-lambda terms", {
  ## exact: S1 = diag(1, 0), S2 = diag(0, 1), so log|S| = rho1 + rho2 (before: d0 20, d1 (0, 1))
  r <- egpd:::.logdetS_multi(c(-10, 20), list(S = list(diag(c(1, 0)), diag(c(0, 1)))), deriv = 2)
  expect_equal(r$d0, 10, tolerance = 1e-12)
  expect_equal(r$d1, c(1, 1), tolerance = 1e-12)
  expect_equal(r$d2, matrix(0, 2, 2), tolerance = 1e-12)
  ## dense blocks with exact zeros: S1 = bdiag(A, D, 0), S2 = bdiag(B, 0, 0), so that
  ## log|S|+ = log|l1 A + l2 B| + 2 rho1 + log|D| over the 5-dimensional range
  set.seed(7)
  A <- crossprod(matrix(rnorm(9), 3))
  B <- crossprod(matrix(rnorm(9), 3))
  D <- crossprod(matrix(rnorm(4), 2))
  S1 <- S2 <- matrix(0, 6, 6)
  S1[1:3, 1:3] <- A
  S1[4:5, 4:5] <- D
  S2[1:3, 1:3] <- B
  for (rho in list(c(0, 0), c(-10, 20), c(-15, 15), c(20, -10))) {
    l <- exp(rho)
    iM <- solve(l[1] * A + l[2] * B)
    r <- egpd:::.logdetS_multi(rho, list(S = list(S1, S2)), deriv = 1)
    expect_equal(r$d0, as.numeric(determinant(l[1] * A + l[2] * B)$modulus) + 2 * rho[1] +
                   as.numeric(determinant(D)$modulus), tolerance = 1e-10)
    expect_equal(r$d1, c(l[1] * sum(iM * A) + 2, l[2] * sum(iM * B)), tolerance = 1e-10)
  }
})

test_that("VA-14: tensor-product log|S|+ matches mgcv's stable evaluation at any spread", {
  skip_if_not(exists("gam.reparam", asNamespace("mgcv")))
  set.seed(3)
  d <- data.frame(x = runif(300), z = runif(300), w = runif(300), y = rnorm(300))
  ## mgcv's own log|S|+ (Wood 2011, Appendix B, in C) on the penalty square roots
  ## projected onto the range space of the total penalty, as in mgcv::gam
  mgcv_ldetS <- function(Sl, rho) {
    St <- Reduce("+", lapply(Sl, function(S) S / sqrt(sum(S * S))))
    es <- eigen(St, symmetric = TRUE)
    Y <- es$vectors[, es$values > max(es$values) * .Machine$double.eps^.66, drop = FALSE]
    rp <- utils::getFromNamespace("gam.reparam", "mgcv")(
      lapply(Sl, function(S) t(Y) %*% mgcv::mroot(S)), rho, 2)
    list(d0 = rp$det, d1 = rp$det1, d2 = rp$det2)
  }
  te2 <- mgcv::gam(y ~ te(x, z, k = c(5, 5)), data = d, fit = FALSE)$smooth[[1]]$S
  te3 <- mgcv::gam(y ~ te(x, z, w, k = c(4, 4, 4)), data = d, fit = FALSE)$smooth[[1]]$S
  ## before: (-15, 15) gave d0 218.09 for 124.60 and d1 (0, 15) for (6, 15); the
  ## three-penalty case failed outright ("number of rows of matrices must match")
  for (case in list(list(te2, c(-2, 5)), list(te2, c(-10, 10)), list(te2, c(-15, 15)),
                    list(te2, c(-10, 20)), list(te3, c(-3, 2, 6)), list(te3, c(-12.5, 0, 12.5)),
                    list(te3, c(15, 0, -15)))) {
    got <- egpd:::.logdetS_multi(case[[2]], list(S = case[[1]]), deriv = 2)
    ref <- mgcv_ldetS(case[[1]], case[[2]])
    expect_equal(got$d0, ref$d0, tolerance = 1e-9)
    expect_equal(got$d1, ref$d1, tolerance = 1e-9)
    expect_equal(got$d2, ref$d2, tolerance = 1e-8)
  }
})

## ---- VA-16: negREML is the REML criterion with mgcv's constants ---------------------

## a Poisson log-likelihood (with log y!) through family = "custom", as in mgcv's poisson()
va16_pois <- list(
  d0 = function(pars, likdata) {
    eta <- as.vector(likdata$X[[1]] %*% pars)
    sum(exp(eta) - likdata$y[, 1] * eta + lgamma(likdata$y[, 1] + 1))
  },
  d120 = function(pars, likdata) {
    mu <- exp(as.vector(likdata$X[[1]] %*% pars))
    cbind(mu - likdata$y[, 1], mu)
  },
  d340 = function(pars, likdata) {
    mu <- exp(as.vector(likdata$X[[1]] %*% pars))
    cbind(mu, mu)
  })

va16_data <- function() {
  set.seed(16)
  n <- 400
  x <- runif(n)
  z <- runif(n)
  data.frame(y = rpois(n, exp(1 + sin(2 * pi * x) + 0.5 * z)), x = x, z = z)
}

test_that("VA-16: negREML equals mgcv's REML score at the same smoothing parameters", {
  d <- va16_data()
  ## supplied smoothing parameters, so no optimiser tolerance is involved
  ## (before: 19.40 above mgcv's score)
  form <- y ~ s(x, k = 8) + s(z, k = 6)
  fe <- egpd(list(form), data = d, family = "custom", custom.fns = va16_pois, sp = c(0.3, 50))
  g <- mgcv::gam(form, data = d, family = poisson, sp = c(0.3, 50), method = "REML")
  expect_equal(fe$negREML, as.numeric(g$gcv.ubre), tolerance = 1e-10)
  ## estimated smoothing parameters, with a parametric term and with a tensor product
  ## (before: 18.43 and 7.35 above)
  for (form in list(y ~ z + s(x, k = 10), y ~ te(x, z, k = c(5, 5)))) {
    fe <- egpd(list(form), data = d, family = "custom", custom.fns = va16_pois)
    g <- mgcv::gam(form, data = d, family = poisson, sp = fe$sp, method = "REML")
    expect_equal(fe$negREML, as.numeric(g$gcv.ubre), tolerance = 1e-9)
  }
})

test_that("log|S|+ is taken from the right smooth when an unpenalised smooth comes first", {
  ## (found with VA-16) .logdetS indexed the smooths by the position of their penalty
  ## group, so s(z, fx = TRUE) before s(x) gave s(x)'s log|S|+ the rank of s(z):
  ## sp 0.0535 against mgcv's REML optimum 0.1145
  d <- va16_data()
  form <- y ~ s(z, k = 5, fx = TRUE) + s(x, k = 8)
  fe <- egpd(list(form), data = d, family = "custom", custom.fns = va16_pois)
  ref <- mgcv::gam(form, data = d, family = poisson, method = "REML")
  expect_equal(unname(fe$sp), unname(ref$sp), tolerance = 2e-3)
  g <- mgcv::gam(form, data = d, family = poisson, sp = fe$sp, method = "REML")
  expect_equal(fe$negREML, as.numeric(g$gcv.ubre), tolerance = 1e-9)
})

test_that("VA-16: negREML of a three-parameter family is the Laplace REML criterion", {
  ## The verifier's independent LAML for the VA-01 fit at its smoothing parameter (own
  ## identity-link DEGPD1 pmf, own Newton inner fit, numDeriv Hessian, eigenvalue
  ## log|S|+) is 1424.417768 without the 2 pi term; the improper prior on the
  ## Mp = 10 - 6 = 4 unpenalised coefficients adds -2 log(2 pi). (Before: 1437.045770.)
  fit <- egpd(list(yd ~ s(x, k = 8), ~1, ~1), data = va01_data(), family = "degpd",
              degpd.args = list(link = "identity"), restarts = FALSE, outer = "fd")
  expect_equal(unname(fit$sp), 0.082138, tolerance = 1e-4)
  expect_equal(fit$negREML, 1424.417768 - 2 * log(2 * pi), tolerance = 1e-8)
})

## ---- VA-10: fits ending on the xi < -1 likelihood singularity -------------------------

## min(1 + xi y / sigma) at a continuous EGPD fit, from its fitted linear predictors
## (log scale; identity-link shape)
va10_zmin <- function(fit, y) min(1 + fit[[2]]$fitted * y / exp(fit[[1]]$fitted))

test_that("VA-10: a coefficient search that stops on an unbounded objective is not converged", {
  ## f -> -Inf as p1 -> 1 (invalid beyond), concave in p1: the perturbed Newton step,
  ## and the decrease it predicts, are negligible at the edge, which used to pass as
  ## code 0 ("step tolerance reached") with a gradient of 1.8e13
  fn <- function(pars, likdata, likfns, newton = TRUE)
    if (pars[1] < 1) 0.01 * log1p(-pars[1]) + (pars[2] - 1)^2 else 1e20
  sfn <- function(pars, likfns, likdata, kept, newton = TRUE) {
    g <- c(if (pars[1] < 1) -0.01 / (1 - pars[1]) else 0, 2 * (pars[2] - 1))
    if (!newton) return(structure(g, gradient = g))
    egpd:::.search.dir(g, diag(c(if (pars[1] < 1) -0.01 / (1 - pars[1])^2 else 0, 2)), kept)
  }
  ctl <- list(steptol = 1e-12, itlim = 100, fntol = 1e-8, gradtol = 1e-4, stepmax = 100)
  ## Exercise the curvature check, not the separate large-decrement stopping rule.
  dectol <- 0.1
  fit <- egpd:::.fit_inner(c(0, 0), fn, sfn, likdata = list(control = list(inner = ctl)),
                           likfns = NULL, dectol = dectol)
  expect_gt(fit$pars[1], 0.99)
  expect_lt(0.5 * fit$decrement, dectol)
  expect_equal(fit$convergence, 3L)
  expect_match(fit$report, "not convex along its gradient")
})

test_that("VA-10: EGPD fits that end on the xi < -1 singularity report non-convergence", {
  ## verifier probe v10_singular.R: for xi < -1 the density grows without bound as
  ## 1 + xi y / sigma -> 0+, so there is no MLE; these fits stopped with one observation
  ## at the upper end of the support (and a gradient of order 1e11-1e13) and reported
  ## convergence 0 with a finite logLik and no warning
  set.seed(2)
  xp <- runif(40)
  d <- data.frame(y = regpd(40, sigma = exp(0.5 + 0.6 * xp), xi = 0.15, kappa = 1.5, type = 1),
                  x1 = xp, x2 = xp^2, x3 = xp^3, x4 = sin(2 * pi * xp))
  expect_warning(fit <- suppressMessages(egpd(list(y ~ x1 + x2 + x3 + x4, ~1, ~ x1 + x2), data = d,
                                              family = "egpd")),
                 class = "egpd_nonconvergence")
  expect_lt(min(fit[[2]]$fitted), -1)
  expect_lt(va10_zmin(fit, d$y), 1e-8)
  expect_equal(fit$convergence, 3L)
  expect_true(is.na(fit$logLik))
  expect_match(fit$convergence.report, "unbounded")

  set.seed(26)
  xs <- runif(50)
  ds <- data.frame(y = regpd(50, sigma = exp(0.5 + 0.6 * sin(2 * pi * xs)), xi = 0.15, kappa = 1.5,
                             type = 1), x = xs)
  expect_warning(fit <- suppressMessages(egpd(list(y ~ s(x, k = 6), ~1, ~ s(x, k = 5)), data = ds,
                                              family = "egpd")),
                 class = "egpd_nonconvergence")
  expect_lt(min(fit[[2]]$fitted), -1)
  expect_lt(va10_zmin(fit, ds$y), 1e-8)
  expect_equal(fit$convergence, 3L)
  expect_true(is.na(fit$logLik))
  expect_match(fit$convergence.report, "at the returned smoothing parameters")

  ## a regular fit of the same design is unaffected
  set.seed(3)
  x <- runif(40)
  d3 <- data.frame(y = regpd(40, sigma = exp(0.5 + 0.6 * x), xi = 0.15, kappa = 1.5, type = 1),
                   x1 = x, x2 = x^2, x3 = x^3, x4 = sin(2 * pi * x))
  expect_silent(fit <- egpd(list(y ~ x1 + x2 + x3 + x4, ~1, ~ x1 + x2), data = d3, family = "egpd"))
  expect_equal(fit$convergence, 0L)
  expect_gt(va10_zmin(fit, d3$y), 0.5)
})
