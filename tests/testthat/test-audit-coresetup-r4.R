## Round-4 regression tests, coresetup (egpd() set-up and fitting core).
## Expected values come from independent implementations (own densities, numDeriv,
## brute-force optimisation) or closed forms; see files/egpd-fix-r4/coresetup/NOTES.md.

## ---- independent reference densities ------------------------------------------------

## continuous EGPD log-densities: H is the GPD(sigma, xi) cdf, h its density
ref_lgpd <- function(y, sigma, xi) {
  lz <- log1p(xi * y / sigma)
  list(lS = -lz / xi, lh = -log(sigma) - (1 / xi + 1) * lz)
}
## model 1 (Naveau et al. 2016, G(v) = v^kappa)
ref_legpd1 <- function(y, sigma, xi, kappa) {
  g <- ref_lgpd(y, sigma, xi)
  log(kappa) + (kappa - 1) * log(-expm1(g$lS)) + g$lh
}
## model 3 (G(v) = 1 - D_delta{(1 - v)^delta}, D_delta the Beta(1/delta, 2) cdf):
## g(v) = (1 + delta) / delta * (1 - (1 - v)^delta)
ref_legpd3 <- function(y, sigma, xi, delta) {
  g <- ref_lgpd(y, sigma, xi)
  log1p(delta) - log(delta) + log(-expm1(delta * g$lS)) + g$lh
}

## ---- VA-02: compact = TRUE with continuous EGPD smooths ----------------------------

va02_rgen <- function(n, x, xi = 0.2, kap = 1.5) {
  u <- runif(n)^(1 / kap)
  exp(0.5 + sin(2 * pi * x)) / xi * ((1 - u)^(-xi) - 1)
}

test_that("VA-02: compact continuous EGPD smooth fits equal the non-compact fits", {
  ## 25 distinct covariate values, n = 300: the compact design has 25 rows, and the
  ## analytic d340 (300 rows) used to be combined with it ("dims do not match")
  set.seed(11)
  n <- 300
  x <- sample(seq(0, 1, length.out = 25), n, replace = TRUE)
  d <- data.frame(y = va02_rgen(n, x), x = x)
  fl <- list(y ~ s(x, k = 6), ~1, ~1)
  for (m in c(1, 3)) {
    fc <- egpd(fl, d, family = "egpd", egpd.args = list(m = m), compact = TRUE,
               restarts = FALSE)
    fn <- egpd(fl, d, family = "egpd", egpd.args = list(m = m), restarts = FALSE)
    expect_equal(fc$convergence, 0L)
    expect_equal(fc$sp, fn$sp, tolerance = 1e-3)
    expect_equal(fc$coefficients, fn$coefficients, tolerance = 1e-4)
    expect_equal(fc$logLik, fn$logLik, tolerance = 1e-6)
    ## the reported log-likelihood is that of the data at the fitted coefficients
    p <- predict(fc, type = "response")
    ref <- if (m == 1) ref_legpd1(d$y, p$scale, p$shape, p$kappa) else
      ref_legpd3(d$y, p$scale, p$shape, p$delta)
    expect_equal(fc$logLik, sum(ref), tolerance = 1e-10)
  }
})

test_that("VA-02: a compact fit with commensurate row counts uses the right REML gradient", {
  ## two smooths and every covariate row exactly twice: the row counts were then
  ## commensurate (n = n_unique * nsp), so the compact fit ran silently with a wrong
  ## REML gradient and stopped at sp 82.8 against 80.0
  set.seed(12)
  xu <- seq(0, 1, length.out = 150)
  x <- rep(xu, 2)
  d2 <- data.frame(y = va02_rgen(length(x), x), x = x)
  f2 <- list(y ~ s(x, k = 6), ~1, ~ s(x, k = 5))
  ## Compare gradients at the same point, not the stopping points on a flat REML curve.
  sp <- c(0.15, 2)
  ctl <- list(inner = list(gradtol = 1e-8))
  fc <- egpd(f2, d2, family = "egpd", compact = TRUE, restarts = FALSE, sp = sp, control = ctl)
  fn <- egpd(f2, d2, family = "egpd", restarts = FALSE, sp = sp, control = ctl)
  expect_null(fc$likfns$d340)
  gc <- egpd:::.reml1.fd(log(sp), fc$likfns, fc$likdata, fc$Sdata, beta = fc$coefficients)
  gn <- egpd:::.reml1(log(sp), fn$likfns, fn$likdata, fn$Sdata, beta = fn$coefficients)
  expect_lt(max(abs(gc - gn)), 1e-3)
  expect_equal(fc$logLik, fn$logLik, tolerance = 1e-6)
})

## ---- VA-04: `fixed` with the sandwich adjustment ------------------------------------

test_that("VA-04: the sandwich adjustment works with `fixed` and counts free coefficients", {
  skip_if_not_installed("numDeriv")
  set.seed(41)
  ncl <- 60
  m <- 8
  n <- ncl * m
  id <- rep(seq_len(ncl), each = m)
  x <- runif(n)
  ue <- rnorm(ncl, sd = 0.3)[id]
  sig <- exp(1 + 0.5 * x + ue)
  u <- runif(n)
  y <- floor(sig / 0.2 * ((1 - u)^(-0.2) - 1))      # kappa = 1: a discretised GPD
  d <- data.frame(y = y, x = x, id = id)
  form <- list(y ~ x, ~1, ~1)
  ## independent: own discrete GPD log-pmf, own MLE, numDeriv information and
  ## cluster-summed scores; k = p / tr(H^-1 J) over the p = 3 free coefficients
  lp <- function(b) {
    s <- exp(b[1] + b[2] * x)
    xi <- exp(b[3])
    lS <- function(t) -log1p(xi * t / s) / xi
    a <- lS(y)
    a + log(-expm1(lS(y + 1) - a))
  }
  nll <- function(b) -sum(lp(b))
  bh <- optim(c(1, 0, log(0.2)), nll, method = "BFGS",
              control = list(reltol = 1e-14, maxit = 1000))$par
  for (i in 1:3) bh <- bh - solve(numDeriv::hessian(nll, bh), numDeriv::grad(nll, bh))
  Hn <- numDeriv::hessian(nll, bh)
  Jc <- crossprod(rowsum(numDeriv::jacobian(lp, bh), id))
  k_ref <- 3 / sum(diag(solve(Hn, Jc)))
  ## headtree: "Hessian of unpenalised MLE not positive definite" (and with force = TRUE,
  ## k was 4/3 of this, so the standard errors shrank instead of growing)
  fm <- egpd(form, d, family = "degpd", fixed = list(kappa = 1), restarts = FALSE,
             sandwich.args = list(id = "id", method = "magnitude"))
  expect_equal(fm$likdata$k, k_ref, tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(fm$Vp))[1:3]), sqrt(diag(solve(k_ref * Hn))),
               tolerance = 1e-4)
  expect_equal(fm$logLik, -nll(bh), tolerance = 1e-8)
  ## the curvature adjustment makes Vp the sandwich covariance H^-1 J H^-1
  fc <- egpd(form, d, family = "degpd", fixed = list(kappa = 1), restarts = FALSE,
             sandwich.args = list(id = "id", method = "curvature"))
  iH <- solve(Hn)
  expect_equal(unname(fc$Vp[1:3, 1:3]), iH %*% Jc %*% iH, tolerance = 1e-4)
  expect_equal(unname(fc$Vp[4, ]), rep(0, 4))
})

## ---- VA-05 / VC-12: other side effects of `fixed` -----------------------------------

test_that("VA-05: holding every parameter fixed is a clear error", {
  set.seed(1)
  yp <- rpois(60, 8)
  expect_error(egpd(list(y ~ 1), data.frame(y = yp), family = "bell", fixed = list(mu = 8)),
               "nothing to estimate")
  expect_error(egpd(list(y ~ 1, ~1, ~1), data.frame(y = yp), family = "degpd",
                    fixed = list(sigma = 5, xi = 0.1, kappa = 2)),
               "nothing to estimate")
})

test_that("VA-05/VC-12: a fixed parameter's inert coefficient is 0, not a start value", {
  set.seed(2)
  yd <- rdiscegpd(200, sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  dd <- data.frame(y = yd, x = seq(0, 1, length.out = 200))
  ## the inert coefficient multiplies an all-zero column: it is not estimated, so it
  ## must not show up as one ("(fixed) -0.3246, z = -Inf, p = 0" on headtree)
  f <- egpd(list(y ~ x, ~1, ~1), data = dd, family = "degpd", fixed = list(kappa = 2),
            restarts = FALSE)
  expect_equal(unname(f$coefficients[f$likdata$idpars == 3]), 0)
  s <- summary(f)[[1]][[3]]
  expect_true(all(is.na(unlist(s[, c("z value", "Pr(>|z|)")]))))
  expect_equal(unique(predict(f, type = "response")$kappa), 2)
  ## the fit is that of the kappa = 2 model: own DEGPD1 log-pmf at the coefficients
  b <- unname(f$coefficients)
  lH <- function(t, s, xi) log(-expm1(-log1p(xi * t / s) / xi))
  s1 <- exp(b[1] + b[2] * dd$x)
  a <- 2 * lH(yd + 1, s1, exp(b[3]))
  bb <- ifelse(yd > 0, 2 * lH(yd, s1, exp(b[3])), -Inf)
  expect_equal(f$logLik, sum(a + log(-expm1(bb - a))), tolerance = 1e-10)
})

test_that("VA-05: count-family clamp boxes are checked at `fixed` values", {
  set.seed(3)
  yo <- rnbinom(200, mu = 10, size = 0.8)
  d <- data.frame(y = yo)
  ## rho = 0.8 is outside the parameter space rho > 1: headtree reported a
  ## non-converged fit from the default start, and blamed `inits` when they were given
  expect_error(egpd(list(y ~ 1, ~1, ~1), d, family = "gwaring", fixed = list(rho = 0.8)),
               "`fixed` put rho at 0.8")
  expect_error(egpd(list(y ~ 1, ~1, ~1), d, family = "gwaring", fixed = list(rho = 0.8),
                    inits = c(log(10), 0, 0)),
               "`fixed` put rho at 0.8")
  ## a value beyond a numerical bound of the box (gpois lambda >= 1e-10) is not
  ## silently evaluated at the bound
  expect_error(egpd(list(y ~ 1, ~1), d, family = "gpois", fixed = list(lambda = 5e-11)),
               "`fixed` put lambda at 5e-11")
})

## ---- VA-17: custom family with a bare formula ----------------------------------------

test_that("VA-17: a custom family given a bare formula has one parameter", {
  set.seed(17)
  n <- 500
  x <- runif(n)
  y <- rpois(n, exp(1 + x))
  d <- data.frame(y = y, x = x)
  pois <- list(
    d0 = function(pars, likdata) {
      eta <- as.vector(likdata$X[[1]] %*% pars)
      sum(exp(eta) - likdata$y[, 1] * eta + lgamma(likdata$y[, 1] + 1))
    },
    d120 = function(pars, likdata) {
      mu <- exp(as.vector(likdata$X[[1]] %*% pars))
      cbind(mu - likdata$y[, 1], mu)
    },
    d340 = NULL)
  ## headtree: three parameters were set up and Armadillo failed ("Mat::col(): index
  ## out of bounds")
  f <- egpd(y ~ x, data = d, family = "custom", custom.fns = pois)
  ref <- glm(y ~ x, family = poisson, data = d)
  expect_equal(unname(f$coefficients), unname(coef(ref)), tolerance = 1e-6)
  expect_equal(f$logLik, as.numeric(logLik(ref)), tolerance = 1e-10)
})

## ---- VA-18: .shape_plateau() and a fixed shape --------------------------------------

test_that("VA-18: a fixed shape is never on the bounded link's plateau", {
  set.seed(18)
  n <- 300
  x <- runif(n)
  yd <- rdiscegpd(n, sigma = exp(1 + 0.5 * x), xi = 0.1, kappa = 1.5, type = 1)
  d <- data.frame(y = yd, x = x)
  f1 <- egpd(list(y ~ x, ~1, ~1), data = d, family = "degpd",
             degpd.args = list(xi.max = 0.2), fixed = list(xi = 5e-5), restarts = FALSE)
  expect_null(egpd:::.shape_plateau(f1, 0.2))
  ## so the default restarts do not add the log-link refit and its mapped start
  ## (headtree: 7 model set-ups)
  calls <- 0
  real <- egpd:::.setup.data
  local_mocked_bindings(.setup.data = function(...) {
    calls <<- calls + 1
    real(...)
  }, .package = "egpd")
  f2 <- egpd(list(y ~ x, ~1, ~1), data = d, family = "degpd",
             degpd.args = list(xi.max = 0.2), fixed = list(xi = 5e-5))
  expect_lte(calls, 5)
  expect_equal(f2$logLik, f1$logLik, tolerance = 1e-8)
})

## ---- VA-12 / VA-13: subsampled set-ups ------------------------------------------------

test_that("VA-12: maxspline keeps every factor level of the model", {
  set.seed(12)
  n <- 2000
  x <- runif(n)
  f <- factor(ifelse(seq_len(n) <= 6, "rare", sample(c("a", "b"), n, TRUE)))
  u <- runif(n)^(1 / 1.5)
  y <- exp(0.5 + sin(2 * pi * x) + 0.3 * (f == "b")) / 0.2 * ((1 - u)^(-0.2) - 1)
  d <- data.frame(y = y, x = x, f = f)
  full <- egpd(list(y ~ s(x, k = 6) + f, ~1, ~1), d, family = "egpd", restarts = FALSE)
  mm <- model.matrix(~ f, d)
  ## seeds whose 200-row subsample has no "rare" row: headtree failed ("number of items
  ## to replace is not a multiple of replacement length")
  for (s in c(1, 2, 4)) {
    set.seed(s)
    fit <- egpd(list(y ~ s(x, k = 6) + f, ~1, ~1), d, family = "egpd", restarts = FALSE,
                maxspline = 200)
    expect_equal(fit$convergence, 0L)
    expect_equal(unname(fit[[1]]$X[, colnames(mm)]), unname(mm[, colnames(mm)]))
    ## only the spline basis comes from the subsample
    expect_equal(fit$logLik, full$logLik, tolerance = 1e-5)
  }
})

test_that("VA-13: one-column data frames can be subsampled", {
  set.seed(13)
  d <- data.frame(y = rdiscegpd(100, sigma = 3, xi = 0.2, kappa = 1.5, type = 1))
  ## headtree: "argument is of length zero" (maxdata) and "'data' must be a data.frame,
  ## environment, or list" (maxspline), as data[id, ] dropped to a vector
  set.seed(1)
  f1 <- suppressMessages(egpd(list(y ~ 1, ~1, ~1), d, family = "degpd", restarts = FALSE,
                              maxdata = 50))
  expect_equal(nobs(f1), 50)
  ## the fit is the MLE on the 50 retained rows (own DEGPD1 log-pmf)
  ys <- f1$data$y
  nll <- function(b) {
    lH <- function(t) log(-expm1(-log1p(exp(b[2]) * t / exp(b[1])) / exp(b[2])))
    a <- exp(b[3]) * lH(ys + 1)
    bb <- ifelse(ys > 0, exp(b[3]) * lH(ys), -Inf)
    -sum(a + log(-expm1(bb - a)))
  }
  ref <- optim(unname(f1$coefficients), nll, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 1000))
  expect_lt(abs(f1$logLik + ref$value), 1e-6)
  f2 <- egpd(list(y ~ 1, ~1, ~1), d, family = "degpd", restarts = FALSE, maxspline = 50)
  f0 <- egpd(list(y ~ 1, ~1, ~1), d, family = "degpd", restarts = FALSE)
  expect_equal(nobs(f2), 100)
  expect_equal(f2$logLik, f0$logLik, tolerance = 1e-10)
})

## ---- VA-15: expressions in formulae --------------------------------------------------

## own DEGPD1 log-pmf
ref_ldegpd1 <- function(y, sigma, xi, kappa) {
  lH <- function(t) log(-expm1(-log1p(xi * t / sigma) / xi))
  a <- kappa * lH(y + 1)
  b <- ifelse(y > 0, kappa * lH(y), -Inf)
  a + log(-expm1(b - a))
}

va15_data <- function() {
  set.seed(15)
  n <- 300
  x <- runif(n, 0.5, 2.5)
  g <- rep(1:3, length.out = n)
  u <- runif(n)^(1 / 1.5)
  y <- floor(exp(0.8 + 0.6 * log(x) + 0.3 * (g == 2)) / 0.2 * ((1 - u)^(-0.2) - 1))
  data.frame(y = y, x = x, g = g)
}

test_that("VA-15: parametric terms may be expressions of data columns", {
  d <- va15_data()
  ## headtree: "Variable(s) 'log(x), factor(g)' not supplied to `data'."
  fit <- egpd(list(y ~ log(x) + factor(g), ~1, ~1), d, family = "degpd", restarts = FALSE)
  expect_equal(fit$convergence, 0L)
  expect_setequal(fit$predictor.names, c("x", "g"))
  ## independent: base R design, own log-pmf, own optimiser started at the package's
  ## estimate cannot improve on it
  X <- model.matrix(~ log(x) + factor(g), d)
  nll <- function(b) -sum(ref_ldegpd1(d$y, exp(drop(X %*% b[1:4])), exp(b[5]), exp(b[6])))
  b <- unname(fit$coefficients)
  expect_equal(fit$logLik, -nll(b), tolerance = 1e-10)
  ref <- optim(b, nll, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  expect_lt(nll(b) - ref$value, 1e-6)
  ## predict() on newdata holding a single level of g rebuilds factor(g) correctly
  nd <- d[d$g == 2, ][1:4, ]
  pn <- predict(fit, newdata = nd, type = "response")
  expect_equal(unname(pn$scale), unname(exp(drop(X[rownames(nd), ] %*% b[1:4]))),
               tolerance = 1e-12)
  ## and so do I() and poly()
  fp <- egpd(list(y ~ poly(x, 2) + I(x > 1.5), ~1, ~1), d, family = "degpd", restarts = FALSE)
  Xp <- cbind(1, d$x, d$x^2, d$x > 1.5)
  nllp <- function(b) -sum(ref_ldegpd1(d$y, exp(drop(Xp %*% b[1:4])), exp(b[5]), exp(b[6])))
  bp <- optim(c(log(mean(d$y) + 1), 0, 0, 0, log(0.2), 0), nllp, method = "BFGS",
              control = list(reltol = 1e-14, maxit = 5000))
  bp <- optim(bp$par, nllp, method = "BFGS", control = list(reltol = 1e-14, maxit = 5000))
  expect_lt(abs(fp$logLik + bp$value), 1e-5)
  expect_equal(unname(as.matrix(predict(fp, newdata = d[5:8, ]))),
               unname(as.matrix(predict(fp)[5:8, ])), tolerance = 1e-12)
  ## a constant held in a local variable is found, as by mgcv::gam()
  deg <- 2
  fl <- egpd(list(y ~ poly(x, deg) + I(x > 1.5), ~1, ~1), d, family = "degpd",
             restarts = FALSE)
  expect_lt(abs(fl$logLik + bp$value), 1e-5)
  expect_error(egpd(list(y ~ log(zz), ~1, ~1), d, family = "degpd", restarts = FALSE),
               "'zz' not supplied")
})

test_that("VA-15: the response may be an expression", {
  d <- va15_data()
  ## headtree: "undefined columns selected"
  fit <- egpd(list(I(y + 1) ~ x, ~1, ~1), d, family = "degpd", restarts = FALSE)
  b <- unname(fit$coefficients)
  nll <- function(b) -sum(ref_ldegpd1(d$y + 1, exp(b[1] + b[2] * d$x), exp(b[3]), exp(b[4])))
  expect_equal(fit$logLik, -nll(b), tolerance = 1e-10)
  ref <- optim(b, nll, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  expect_lt(nll(b) - ref$value, 1e-6)
  expect_equal(nrow(predict(fit, newdata = d[1:3, "x", drop = FALSE])), 3)
})

## ---- VC-03: gamma only affects smoothing-parameter selection (mgcv semantics) -------

test_that("VC-03: with gamma, Vp is the covariance of the unscaled penalised likelihood", {
  skip_if_not_installed("numDeriv")
  set.seed(5)
  n <- 300
  x <- runif(n)
  y <- rdiscegpd(n, sigma = exp(1 + sin(2 * pi * x)), xi = 0.2, kappa = 1.5, type = 1)
  d <- data.frame(y = y, x = x)
  fit <- egpd(list(y ~ s(x, k = 8), ~1, ~1), d, family = "degpd", gamma = 4,
              restarts = FALSE)
  ## own penalised log-likelihood l(b) - b'S b / 2, S the penalty at the reported
  ## smoothing parameters (those of the penalty on the unscaled likelihood)
  X1 <- fit$likdata$X[[1]]
  S <- Reduce("+", Map("*", attr(fit$Sdata, "Sl"), fit$sp))
  npl <- function(b) -sum(ref_ldegpd1(y, exp(drop(X1 %*% b[1:8])), exp(b[9]), exp(b[10]))) +
    0.5 * sum(b * (S %*% b))
  b <- unname(fit$coefficients)
  ## the fit maximises it ...
  expect_lt(max(abs(numDeriv::grad(npl, b))), 1e-3)
  ## ... and Vp is the inverse of its negative Hessian (headtree: 4 times too large)
  Hn <- numDeriv::hessian(npl, b)
  expect_equal(unname(fit$Vp), solve(Hn), tolerance = 1e-5)
  ## Vc adds smoothing-parameter uncertainty to that
  expect_true(all(diag(fit$Vc) >= diag(fit$Vp) * (1 - 1e-8)))
  ## supplying the reported smoothing parameters gives the same fit, whatever gamma is
  f1 <- egpd(list(y ~ s(x, k = 8), ~1, ~1), d, family = "degpd", sp = fit$sp,
             restarts = FALSE)
  expect_equal(f1$coefficients, fit$coefficients, tolerance = 1e-6)
  expect_equal(f1$Vp, fit$Vp, tolerance = 1e-5)
})

test_that("VC-03: gamma follows mgcv for smoothing parameters and Vp", {
  pois <- list(
    d0 = function(pars, likdata) {
      eta <- as.vector(likdata$X[[1]] %*% pars)
      sum(exp(eta) - likdata$y[, 1] * eta + lgamma(likdata$y[, 1] + 1))
    },
    d120 = function(pars, likdata) {
      mu <- exp(as.vector(likdata$X[[1]] %*% pars))
      cbind(mu - likdata$y[, 1], mu)
    },
    d340 = NULL)
  set.seed(3)
  n <- 300
  x <- runif(n)
  d <- data.frame(y = rpois(n, exp(0.5 + sin(2 * pi * x))), x = x)
  f <- egpd(list(y ~ s(x, k = 10)), d, family = "custom", custom.fns = pois, gamma = 4)
  m <- mgcv::gam(y ~ s(x, k = 10), family = poisson, data = d, method = "REML", gamma = 4)
  ## headtree: sp 0.109 against mgcv's 0.438, and Vp 4 times mgcv's
  expect_equal(unname(f$sp), unname(m$sp), tolerance = 5e-3)
  expect_lt(max(abs(f$Vp - m$Vp)) / max(abs(m$Vp)), 5e-3)
})

test_that("VC-03: without smooths gamma leaves the covariance alone", {
  ## Bell(theta), log link: theta-hat solves theta exp(theta) = ybar, and the
  ## information for log(theta) is n exp(theta) theta (1 + theta)
  set.seed(2)
  y <- rbell(80, theta = 0.9)
  th <- uniroot(function(t) t * exp(t) - mean(y), c(1e-6, 10), tol = 1e-14)$root
  info <- length(y) * exp(th) * th * (1 + th)
  f <- egpd(y ~ 1, data.frame(y = y), family = "bellnat", gamma = 4)
  expect_equal(exp(unname(f$coefficients)), th, tolerance = 1e-6)
  expect_equal(as.numeric(f$Vp), 1 / info, tolerance = 1e-6)
  expect_equal(summary(f)[[1]][[1]][1, "Std. Error"], 1 / sqrt(info), tolerance = 1e-6)
  expect_error(egpd(y ~ 1, data.frame(y = y), family = "bellnat", gamma = 0),
               "`gamma` must be a single positive number")
})

## ---- VA-11: restarts fit one data set ------------------------------------------------

test_that("VA-11: the data are evaluated once and every restart sees the same rows", {
  own_mle <- function(y) {
    nll <- function(b) -sum(ref_ldegpd1(y, exp(b[1]), exp(b[2]), exp(b[3])))
    o <- optim(c(log(mean(y) + 1), log(0.2), 0), nll, method = "BFGS",
               control = list(reltol = 1e-15, maxit = 5000))
    o <- optim(o$par, nll, method = "Nelder-Mead", control = list(reltol = 1e-15, maxit = 5000))
    optim(o$par, nll, method = "BFGS", control = list(reltol = 1e-15, maxit = 5000))$value
  }
  u <- seq(0.005, 0.995, length.out = 120)
  y0 <- floor(3 / 0.2 * ((1 - u^(1 / 1.5))^(-0.2) - 1))
  calls <- 0L
  make_data <- function() {
    calls <<- calls + 1L
    data.frame(y = y0 * calls)
  }
  ## headtree evaluated make_data() for every start (6 times) and returned the fit to
  ## the sixth data set, y0 * 6
  f <- egpd(list(y ~ 1, ~1, ~1), data = make_data(), family = "degpd")
  expect_equal(calls, 1L)
  expect_lt(abs(f$logLik + own_mle(y0)), 1e-5)

  set.seed(7)
  n <- 400
  yy <- rdiscegpd(n, sigma = 4, xi = 0.25, kappa = 1.5, type = 1)
  d <- data.frame(y = yy, id = seq_len(n))
  rows <- list()
  real <- egpd:::.setup.data
  local_mocked_bindings(.setup.data = function(...) {
    out <- real(...)
    rows[[length(rows) + 1L]] <<- out$data$id
    out
  }, .package = "egpd")
  set.seed(99)
  ## headtree: five different subsamples, and the fit to the easiest one was returned
  fm <- suppressMessages(egpd(list(y ~ 1, ~1, ~1), d, family = "degpd", maxdata = 150))
  expect_length(rows, 5)
  expect_length(unique(rows), 1)
  expect_lt(abs(fm$logLik + own_mle(yy[rows[[1]]])), 1e-5)
})

## ---- VA-19: restarts for the continuous EGPD -------------------------------------------

## own continuous EGPD model 2 log-likelihood: G(u) = p u^k1 + (1 - p) u^k2, k2 = k1 + dk
ref_negll_egpd2 <- function(t, y) {
  s <- exp(t[1]); xi <- t[2]; k1 <- exp(t[3]); k2 <- k1 + exp(t[4]); p <- plogis(t[5])
  lz <- log1p(xi * y / s)
  lH <- log(-expm1(-lz / xi))
  -sum(log(p * k1 * exp((k1 - 1) * lH) + (1 - p) * k2 * exp((k2 - 1) * lH)) -
         log(s) - (1 / xi + 1) * lz)
}

test_that("VA-19: continuous EGPD fits are restarted by default", {
  set.seed(42)
  n <- 400
  u <- runif(n)^(1 / 1.5)
  d <- data.frame(y = 3 / 0.2 * ((1 - u)^(-0.2) - 1))
  fl <- list(y ~ 1, ~1, ~1, ~1, ~1)
  f <- suppressMessages(egpd(fl, d, family = "egpd", egpd.args = list(m = 2)))
  f0 <- suppressMessages(egpd(fl, d, family = "egpd", egpd.args = list(m = 2),
                              restarts = FALSE))
  expect_equal(f$logLik, -ref_negll_egpd2(unname(f$coefficients), d$y), tolerance = 1e-10)
  ## -1023.248077: the best of 40 random starts of BFGS on the likelihood above
  ## (verify-A v19_m2.R); the default start alone stops 1.65 units lower (headtree's fit)
  expect_equal(f$logLik, -1023.248077, tolerance = 1e-8)
  expect_equal(f0$logLik, -1024.896745, tolerance = 1e-8)

  ## where the default start reaches the optimum, the fit is the single-start fit itself
  set.seed(11)
  x <- runif(300)
  u <- runif(300)^(1 / 1.5)
  d1 <- data.frame(y = exp(0.5 + sin(2 * pi * x)) / 0.2 * ((1 - u)^(-0.2) - 1), x = x)
  calls <- 0
  real <- egpd:::.setup.data
  local_mocked_bindings(.setup.data = function(...) {
    calls <<- calls + 1
    real(...)
  }, .package = "egpd")
  g <- egpd(list(y ~ s(x, k = 6), ~1, ~1), d1, family = "egpd")
  expect_equal(calls, 5)
  calls <- 0
  g0 <- egpd(list(y ~ s(x, k = 6), ~1, ~1), d1, family = "egpd", restarts = FALSE)
  expect_equal(calls, 1)
  expect_identical(g$coefficients, g0$coefficients)
  expect_identical(g$sp, g0$sp)
  expect_identical(g$Vp, g0$Vp)
  ## starts that differ only in a parameter held by `fixed` are fitted once
  calls <- 0
  gf <- egpd(list(y ~ s(x, k = 6), ~1, ~1), d1, family = "egpd", fixed = list(kappa = 1.5))
  expect_equal(calls, 3)
})
