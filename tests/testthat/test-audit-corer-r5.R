## Round-5 fixes to the model set-up (corer): VA5-09 (t2() smooths and mgcv's fit basis),
## VA5-06 (maxspline level coverage of factor expressions in smooths), VA5-07
## (removeData = TRUE with a finite xi.max), VA5-03/04/05 (expressions in formulae) and
## VB5-09 (a named formula list with the response on a later element).

## a unit-variance Gaussian log-likelihood through family = "custom", the model of
## mgcv::gam(family = gaussian, scale = 1)
r5_gauss <- list(
  d0 = function(pars, likdata)
    -sum(dnorm(likdata$y[, 1], as.vector(likdata$X[[1]] %*% pars), 1, log = TRUE)),
  d120 = function(pars, likdata)
    cbind(as.vector(likdata$X[[1]] %*% pars) - likdata$y[, 1], rep(1, nrow(likdata$y))),
  d340 = function(pars, likdata) matrix(0, nrow(likdata$y), 2))

va509_data <- function(n = 180) {
  set.seed(817)
  d <- data.frame(x = runif(n), z = runif(n))
  d$y <- sin(6 * d$x) + 0.5 * d$z + rnorm(n)
  d
}

## ---- VA5-09: t2() terms are set up in mgcv's fit basis --------------------------------

test_that("VA5-09: t2() terms predict on newdata as mgcv does", {
  d <- va509_data()
  nd <- data.frame(x = c(0.15, 0.4, 0.7), z = c(0.5, 0.2, 0.9))
  z <- qnorm(0.975)
  for (form in list(y ~ t2(x, z, k = c(4, 4)), y ~ s(x, k = 5) + t2(x, z, k = c(4, 4)))) {
    m <- mgcv::gam(form, data = d, scale = 1, method = "REML")
    fit <- egpd(list(form), d, family = "custom", custom.fns = r5_gauss, sp = m$sp)
    ## the training fit already matched mgcv
    expect_equal(predict(fit)[[1]], as.vector(fitted(m)), tolerance = 1e-10)
    ## headtree: off by up to 1.87 (the design on newdata was built in mgcv's
    ## prediction basis, the coefficients are in its fit basis)
    pm <- predict(m, nd, se.fit = TRUE)
    pe <- predict(fit, newdata = nd, se.fit = TRUE, marginal = FALSE)
    expect_equal(pe$fitted[[1]], as.vector(pm$fit), tolerance = 1e-10)
    expect_equal(pe$se.fit[[1]], as.vector(pm$se.fit), tolerance = 1e-8)
    expect_equal(predict(fit, newdata = d)[[1]], as.vector(fitted(m)), tolerance = 1e-10)
    ## Wald intervals at newdata use the same design: mgcv's fit -/+ z se
    ci <- confint(fit, newdata = nd, marginal = FALSE)
    expect_equal(unname(ci[, "lower"]), as.vector(pm$fit - z * pm$se.fit), tolerance = 1e-8)
    expect_equal(unname(ci[, "upper"]), as.vector(pm$fit + z * pm$se.fit), tolerance = 1e-8)
  }
})

test_that("VA5-09: compact fits with a t2() term are mgcv's fits", {
  d <- va509_data()
  form <- y ~ t2(x, z, k = c(4, 4))
  m <- mgcv::gam(form, data = d, scale = 1, method = "REML")
  fit <- egpd(list(form), d, family = "custom", custom.fns = r5_gauss, sp = m$sp,
              compact = TRUE)
  ## headtree: negREML 282.842150568 and logLik -263.370360418, fitted values 0.69 off
  expect_equal(fit$negREML, as.numeric(m$gcv.ubre), tolerance = 1e-10)
  expect_equal(fit$logLik, sum(dnorm(d$y, fitted(m), 1, log = TRUE)), tolerance = 1e-10)
  expect_equal(predict(fit)[[1]], as.vector(fitted(m)), tolerance = 1e-10)
})

test_that("VA5-09: a maxspline subsample sets t2() terms up as mgcv does", {
  d <- va509_data(400)
  form <- y ~ t2(x, z, k = c(4, 4))
  set.seed(2)
  fit <- egpd(list(form), d, family = "custom", custom.fns = r5_gauss, maxspline = 380,
              sp = c(1e3, 0.3, 0.01))
  expect_equal(fit$convergence, 0L)
  ## the basis is set up on the subsample, and on its rows the design is mgcv's
  ## (headtree: up to 5.1 off, which paired a prediction-basis design with fit-basis
  ## penalties)
  set.seed(2)
  id <- sample(nrow(d), 380)
  G <- mgcv::gam(form, data = d[id, ], fit = FALSE, method = "REML")
  expect_equal(unname(fit[[1]]$X[id, ]), unname(G$X), tolerance = 1e-10)
})

test_that("VA5-09: plot() gives mgcv's plotting code prediction-basis coefficients", {
  d <- va509_data()
  form <- y ~ t2(x, z, k = c(4, 4))
  m <- mgcv::gam(form, data = d, scale = 1, method = "REML")
  fit <- egpd(list(form), d, family = "custom", custom.fns = r5_gauss, sp = m$sp)
  got <- NULL
  local_mocked_bindings(plot.gam = function(x, ...) got <<- x, .package = "mgcv")
  plot(fit, onepage = FALSE, ask = FALSE)
  ## as mgcv::gam() reports them, beta = P beta_fit and Vp = P Vp_fit P' (headtree:
  ## the fit-basis coefficients, and a surface up to 2.15 off)
  expect_equal(unname(got$coefficients), unname(coef(m)), tolerance = 1e-8)
  expect_equal(unname(got$Vp), unname(m$Vp), tolerance = 1e-8)
})

## ---- VA5-06: maxspline covers factor expressions in smooths --------------------------

## a Poisson log-likelihood (with log y!) through family = "custom"
r5_pois <- list(
  d0 = function(pars, likdata) {
    eta <- as.vector(likdata$X[[1]] %*% pars)
    sum(exp(eta) - likdata$y[, 1] * eta + lgamma(likdata$y[, 1] + 1))
  },
  d120 = function(pars, likdata) {
    mu <- exp(as.vector(likdata$X[[1]] %*% pars))
    cbind(mu - likdata$y[, 1], mu)
  },
  d340 = NULL)

test_that("VA5-06: a level that the maxspline subsample misses keeps its by-smooth", {
  ## level g = 3 has 6 of 2000 rows, and the 200-row subsample drawn after set.seed(1)
  ## misses all of them
  set.seed(12)
  n <- 2000
  x <- runif(n)
  g <- ifelse(seq_len(n) <= 6, 3L, sample(1:2, n, replace = TRUE))
  eta <- 1.2 + 0.8 * sin(2 * pi * x) * (g == 1) - 1.1 * cos(2 * pi * x) * (g == 3)
  d <- data.frame(y = rpois(n, exp(eta)), x = x, g = g)
  set.seed(1)
  expect_false(any(sample(n, 200) <= 6))
  set.seed(1)
  fit <- egpd(y ~ s(x, by = factor(g), k = 6), data = d, family = "custom",
              custom.fns = r5_pois, maxspline = 200)
  ## mgcv sets the model up with one smooth per level, 16 coefficients (headtree: 2
  ## smooths and 11 columns, no warning, and flat predictions at g = 3)
  mg <- mgcv::gam(y ~ s(x, by = factor(g), k = 6), family = poisson, data = d,
                  method = "REML")
  expect_equal(vapply(fit[[1]]$smooth, `[[`, "", "label"),
               vapply(mg$smooth, `[[`, "", "label"))
  expect_equal(ncol(fit[[1]]$X), ncol(model.matrix(mg)))
  s3 <- fit[[1]]$smooth[[3]]
  expect_equal(unname(which(rowSums(abs(fit[[1]]$X[, s3$first.para:s3$last.para])) > 0)),
               which(d$g == 3))
  ## predictions at g = 3 follow x, close to mgcv's full-data fit (headtree: 1.207 at
  ## every x, against mgcv's 0.652, 1.786, 0.757)
  nd <- data.frame(x = c(0.1, 0.5, 0.9), g = 3L)
  expect_equal(predict(fit, newdata = nd)[[1]], as.vector(predict(mg, nd)), tolerance = 0.1)
  ## and they are those of the same model with the factor as a column
  d$gf <- factor(d$g)
  set.seed(1)
  ref <- egpd(y ~ s(x, by = gf, k = 6), data = d, family = "custom", custom.fns = r5_pois,
              maxspline = 200)
  expect_equal(fit$logLik, ref$logLik, tolerance = 1e-10)
  nd$gf <- factor(3L, levels = 1:3)
  expect_equal(predict(fit, newdata = nd)[[1]], predict(ref, newdata = nd)[[1]],
               tolerance = 1e-10)
})

## ---- VA5-07: removeData = TRUE with a finite xi.max ----------------------------------

test_that("VA5-07: the shape-plateau restart works without the data", {
  set.seed(1)
  dp <- data.frame(y = rpois(150, 4), x = runif(150))
  ## headtree: "Supply `egpd' with `removeData = FALSE' if not supplying `newdata'."
  expect_message(
    fit <- egpd(list(y ~ x, ~1, ~1), data = dp, family = "degpd",
                degpd.args = list(xi.max = 0.3), removeData = TRUE, trace = 1),
    "trying the log-link fit as a start")
  expect_equal(fit$convergence, 0L)
  expect_null(fit$data)
  ## own bounded-link DEGPD1 log-likelihood at the coefficients (the shape is near 0,
  ## so in a form that keeps its precision there)
  b <- unname(fit$coefficients)
  sig <- exp(b[1] + b[2] * dp$x)
  xi <- -0.3 * expm1(-exp(b[3]) / 0.3)
  kap <- exp(b[4])
  lH <- function(t) log(-expm1(-log1p(xi * t / sig) / xi))
  a <- kap * lH(dp$y + 1)
  bb <- ifelse(dp$y > 0, kap * lH(dp$y), -Inf)
  expect_equal(fit$logLik, sum(a + log(-expm1(bb - a))), tolerance = 1e-10)
  ## removing the data does not change the fit
  keep <- suppressMessages(egpd(list(y ~ x, ~1, ~1), data = dp, family = "degpd",
                                degpd.args = list(xi.max = 0.3)))
  expect_equal(fit$coefficients, keep$coefficients, tolerance = 1e-12)
})

## ---- VA5-03: smooths of expressions ----------------------------------------------------

va503_data <- function() {
  set.seed(5)
  n <- 300
  d <- data.frame(x = runif(n, 1, 10), z = runif(n), g = sample(c("a", "b", "c"), n, TRUE))
  d$y <- regpd(n, sigma = exp(0.3 + sin(log(d$x))), xi = 0.1, kappa = 1.5)
  d
}

## own continuous EGPD (carrier 1) log-density
ref_legpd1 <- function(y, sigma, xi, kappa) {
  t <- 1 + xi * y / sigma
  log(kappa) - log(sigma) - (1 / xi + 1) * log(t) + (kappa - 1) * log1p(-t^(-1 / xi))
}

test_that("VA5-03: smooths of expressions are fitted, as mgcv sets them up", {
  d <- va503_data()
  d$lx <- log(d$x)
  d$gf <- factor(d$g)
  nd <- data.frame(x = c(2, 7), z = c(0.3, 0.8), g = c("a", "c"))
  nd$lx <- log(nd$x)
  nd$gf <- factor(nd$g, levels = c("a", "b", "c"))
  forms <- list(list(y ~ s(log(x), k = 6), y ~ s(lx, k = 6)),
                list(y ~ te(log(x), z, k = 4), y ~ te(lx, z, k = 4)),
                list(y ~ s(x, factor(g), bs = "fs", k = 4), y ~ s(x, gf, bs = "fs", k = 4)))
  for (f in forms) {
    ## headtree: "undefined columns selected", after fitting
    fit <- egpd(list(f[[1]], ~1, ~1), d, family = "egpd", restarts = FALSE)
    expect_equal(fit$convergence, 0L)
    ## the design is mgcv's, and logLik is the log-likelihood at the coefficients
    G <- mgcv::gam(f[[1]], data = d, fit = FALSE, method = "REML")
    expect_equal(unname(fit[[1]]$X), unname(G$X), tolerance = 1e-10)
    b <- fit$coefficients
    p <- ncol(G$X)
    expect_equal(fit$logLik, sum(ref_legpd1(d$y, exp(drop(G$X %*% b[1:p])), b[p + 1],
                                            exp(b[p + 2]))), tolerance = 1e-10)
    ## it is the fit on the precomputed columns
    ref <- egpd(list(f[[2]], ~1, ~1), d, family = "egpd", restarts = FALSE)
    expect_equal(fit$logLik, ref$logLik, tolerance = 1e-10)
    expect_equal(fit$sp, ref$sp, tolerance = 1e-8)
    expect_equal(as.matrix(predict(fit, newdata = nd)), as.matrix(predict(ref, newdata = nd)),
                 tolerance = 1e-8)
  }
  ## plotdata holds the expressions' values
  expect_equal(fit$plotdata[[1]][["factor(g)"]], factor(d$g))
  ## headtree: six complete fits under the default restarts, then the same error
  fr <- egpd(list(y ~ s(log(x), k = 6), ~1, ~1), d, family = "egpd")
  expect_equal(fr$convergence, 0L)
  expect_equal(fr$plotdata[[1]][["log(x)"]], log(d$x))
})

test_that("VA5-03: when every restart fails, the warning says so", {
  d <- va503_data()
  ## headtree: "restarts requested but no restart converged; ..."
  expect_warning(
    expect_error(egpd(list(y ~ zz, ~1, ~1), d, family = "egpd"), "'zz' not supplied"),
    "every start failed with an error")
})

## ---- VA5-04: a response expression is evaluated once, in its own environment --------

test_that("VA5-04: a response expression with a local constant", {
  d <- data.frame(y = rep(c(0, 1, 2, 3, 5, 8, 10), 10))
  rf <- local({
    shift_local <- 1
    I(y + shift_local) ~ 1
  })
  ## independent MLE of the discrete GPD (DEGPD carrier 1 at kappa = 1, xi = 0.2) of y + 1
  z <- d$y + 1
  nll <- function(ls) {
    s <- exp(ls)
    -sum(log((1 + 0.2 * z / s)^-5 - (1 + 0.2 * (z + 1) / s)^-5))
  }
  ref <- -optimize(nll, c(-5, 10), tol = 1e-12)$objective
  expect_equal(ref, -195.353231649, tolerance = 1e-11)
  ## headtree: "object 'shift_local' not found" (and, for the second list, "formula has
  ## no response")
  for (fl in list(list(rf, ~1, ~1), list(lsigma = ~1, lxi = rf, lkappa = ~1),
                  list(lxi = ~1, lsigma = rf, lkappa = ~1))) {
    fit <- egpd(fl, d, family = "degpd", fixed = list(xi = 0.2, kappa = 1))
    expect_equal(fit$logLik, ref, tolerance = 1e-9)
    expect_equal(fit$response.name, "I(y + shift_local)")
  }
  ## the expression is evaluated once, by whatever number of fits (restarts) and formulae
  calls <- 0
  shift <- function(v) {
    calls <<- calls + 1
    v + 1
  }
  fit <- egpd(list(shift(y) ~ 1, ~1, ~1), d, family = "degpd", fixed = list(xi = 0.2, kappa = 1))
  expect_equal(fit$logLik, ref, tolerance = 1e-9)
  expect_equal(calls, 1)
})

## ---- VA5-05: only a name refers to a column -----------------------------------------------

test_that("VA5-05: a response expression is evaluated even if a column is named like it", {
  d <- data.frame(y = 1:20)
  d[["I(y + 1)"]] <- d$y + 101
  ## headtree: the column (mean 111.5) was the response
  fit <- egpd(I(y + 1) ~ 1, data = d, family = "bell")
  ## the Bell MLE has mean ybar = 11.5: theta e^theta = 11.5, and the log-likelihood
  ## sum(y log(theta) + 1 - e^theta + log(B_y) - log(y!)), with B_y the Bell numbers
  expect_equal(predict(fit, type = "response")$mu, rep(11.5, 20), tolerance = 1e-6)
  bell <- 1
  for (m in 0:20) bell[m + 2] <- sum(choose(m, 0:m) * bell[1:(m + 1)])
  theta <- uniroot(function(t) t * exp(t) - 11.5, c(0, 5), tol = 1e-14)$root
  y <- d$y + 1
  expect_equal(fit$logLik, sum(y * log(theta) + 1 - exp(theta) + log(bell[y + 1]) -
                                 lgamma(y + 1)), tolerance = 1e-9)
})

test_that("VA5-05: a parametric expression reads its variables, not a same-named column", {
  set.seed(3)
  d <- data.frame(y = rpois(40, 3), x = rep(2:21, 2))
  d[["log(x)"]] <- 0
  fit <- egpd(y ~ log(x), data = d, family = "bell")
  ## headtree: predictor.names "log(x)", and predict() required that column in newdata
  expect_equal(fit$predictor.names, "x")
  nd <- data.frame(x = c(3, 10))
  expect_equal(predict(fit, newdata = nd)[[1]],
               drop(cbind(1, log(nd$x)) %*% fit$coefficients), tolerance = 1e-12)
  ## headtree: compact keyed the rows by the constant column, a fit with slope 0
  fc <- egpd(y ~ log(x), data = d, family = "bell", compact = TRUE)
  expect_equal(fc$coefficients, fit$coefficients, tolerance = 1e-10)
  expect_equal(fc$logLik, fit$logLik, tolerance = 1e-10)
  ## headtree: a missing x was not seen ("incompatible matrix dimensions")
  d$x[3] <- NA
  expect_warning(fn <- egpd(y ~ log(x), data = d, family = "bell"), "Omitting 1 rows")
  expect_equal(nobs(fn), 39)
  expect_equal(fn$logLik, egpd(y ~ log(x), data = d[-3, ], family = "bell")$logLik,
               tolerance = 1e-12)
})

## ---- VB5-09: a named list with the response on a later element ----------------------

test_that("VB5-09: the response may be on a later element of a canonically ordered list", {
  set.seed(11)
  n <- 400
  x <- runif(n)
  pos <- rbinom(n, 1, 0.7)
  y <- ifelse(pos == 1, 1 + rgeom(n, 1 / (1 + 3 * exp(0.5 * x))), 0)
  d <- data.frame(y = y, x = x)
  dp <- d[d$y > 0, ]
  ## headtree: "formula has no response"
  fz <- egpd(list(lsigma = ~x, lxi = y ~ 1, lkappa = ~1), dp, family = "ztdegpd")
  expect_equal(fz$convergence, 0L)
  ## own zero-truncated DEGPD1 log-likelihood: the fit is its maximum
  nll <- function(b) {
    sig <- exp(b[1] + b[2] * dp$x)
    xi <- exp(b[3])
    ka <- exp(b[4])
    H <- function(t) 1 - (1 + xi * t / sig)^(-1 / xi)
    -sum(log(H(dp$y + 1)^ka - H(dp$y)^ka) - log1p(-H(1)^ka))
  }
  b <- unname(fz$coefficients)
  expect_equal(fz$logLik, -nll(b), tolerance = 1e-10)
  opt <- optim(b, nll, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  expect_lt(nll(b) - opt$value, 1e-6)
  ## headtree: "the zero-truncated degpd fit failed from every start"; the zero part of
  ## the hurdle is intercept-only (the terms of the response's formula), a closed form
  h <- fit_hurdle(list(lsigma = ~x, lxi = y ~ 1, lkappa = ~1), d, family = "degpd")
  n1 <- sum(d$y > 0)
  expect_equal(h$logLik, n1 * log(n1 / n) + (n - n1) * log(1 - n1 / n) - opt$value,
               tolerance = 1e-8)
  hg <- fit_hurdle(list(lmu = ~x, lk = y ~ 1, lrho = ~1), d, family = "gwaring")
  hr <- fit_hurdle(list(lk = y ~ 1, lmu = ~x, lrho = ~1), d, family = "gwaring")
  expect_equal(hg$logLik, hr$logLik, tolerance = 1e-10)
})
