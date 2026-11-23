## Round-5 fixes (corediag): effective degrees of freedom where the log-likelihood
## Hessian is indefinite (VA5-12), the plateau flag (VA5-11), and the default start of
## the zero-truncated generalized Waring (VA5-10).

## ---- VA5-12: edf by mgcv's rule where H0 is indefinite ------------------------------

## mgcv:::gam.fit5.post.proc's clamp loop, applied to H0 itself: eigenvalues of H0 below
## tol * (its largest) are raised to that level, tol = 0, 1e-7, 1.1e-6, ..., until the
## clamped H0 plus S has full rank; the edf are diag((Hc + S)^-1 Hc).
mgcv_rule_edf <- function(H0, S) {
  eh <- eigen(H0, symmetric = TRUE)
  mev <- max(eh$values)
  tol <- 0
  dtol <- 1e-7
  repeat {
    v <- eh$values
    v[v < tol * mev] <- tol * mev
    Hc <- crossprod(sqrt(v) * t(eh$vectors))
    if (attr(suppressWarnings(chol(Hc + S, pivot = TRUE)), "rank") == ncol(S)) break
    tol <- tol + dtol
    dtol <- dtol * 10
  }
  diag(solve(Hc + S, Hc))
}

test_that("VA5-12: edf, df and AIC follow mgcv's rule where the Hessian is indefinite", {
  skip_on_cran()
  ## the verifier's zero-truncated DEGPD hurdle, seed 11: smooths on all three parameters
  set.seed(11)
  n <- 250
  x <- runif(n)
  y <- rdiscegpd(n, sigma = exp(0.3 + sin(2 * pi * x)), xi = 0.2, kappa = 1.2, type = 1)
  h <- suppressWarnings(fit_hurdle(y ~ s(x, k = 8), data = data.frame(y = y, x = x),
                                   family = "degpd"))
  f <- h$positive
  expect_identical(f$convergence, 0L)
  H0 <- egpd:::.gH.nopen(unname(f$coefficients), f$likdata, f$likfns)[[2]]
  S <- f$likdata$S
  ## H0 is indefinite at this valid penalised optimum (H0 + S is positive definite), so
  ## tr((H0 + S)^-1 H0) is unbounded below: it was -51.5 (and -81.0 with REML restarts)
  ev <- eigen(H0, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(min(ev), -1)
  expect_gt(min(eigen(H0 + S, symmetric = TRUE, only.values = TRUE)$values), 0)
  edf <- mgcv_rule_edf(H0, S)
  expect_equal(f$df, sum(edf), tolerance = 1e-8)
  expect_equal(unname(unlist(lapply(1:3, function(i) f[[i]]$edf))), edf, tolerance = 1e-6)
  expect_gt(f$df, 0)
  expect_equal(attr(logLik(f), "df"), f$df)
  expect_equal(f$AIC, -2 * f$logLik + 2 * f$df)
  expect_gte(f$AIC, -2 * f$logLik)
  expect_equal(h$df, sum(h$zero$edf) + f$df)
  expect_gte(h$AIC, -2 * h$logLik)
  ## Vp is still the inverse of the exact penalised Hessian
  expect_equal(f$Vp %*% (H0 + S), diag(ncol(S)), tolerance = 1e-6)

  ## Check the rule at this fit's smoothing parameters, not a previous optimum.
  dpos <- data.frame(y = y, x = x)[y > 0, ]
  f1 <- suppressWarnings(egpd(y ~ s(x, k = 8), data = dpos, family = "ztdegpd"))
  H1 <- egpd:::.gH.nopen(unname(f1$coefficients), f1$likdata, f1$likfns)[[2]]
  expect_identical(f1$convergence, 0L)
  expect_lt(min(eigen(H1, symmetric = TRUE, only.values = TRUE)$values), -1)
  expect_equal(f1$df, sum(mgcv_rule_edf(H1, f1$likdata$S)), tolerance = 1e-8)
  f2 <- suppressWarnings(egpd(y ~ s(x, k = 8), data = dpos, family = "ztdegpd", gamma = 1.4))
  H2 <- 1.4 * egpd:::.gH.nopen(unname(f2$coefficients), f2$likdata, f2$likfns)[[2]]
  expect_lt(min(eigen(H2, symmetric = TRUE, only.values = TRUE)$values), -0.1)
  expect_equal(f2$df, sum(mgcv_rule_edf(H2, 1.4 * f2$likdata$S)), tolerance = 1e-8)

  ## where H0 is positive definite the edf are the plain trace, as before
  set.seed(1)
  x <- runif(300)
  d <- data.frame(x = x, y = rdiscegpd(300, sigma = exp(0.5 + sin(2 * pi * x)), xi = 0.2,
                                       kappa = 1.2, type = 1))
  g <- egpd(list(y ~ s(x, k = 6), ~ 1, ~ 1), data = d, family = "degpd")
  G0 <- egpd:::.gH.nopen(unname(g$coefficients), g$likdata, g$likfns)[[2]]
  expect_gt(min(eigen(G0, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_equal(g$df, sum(diag(solve(G0 + g$likdata$S, G0))), tolerance = 1e-10)
})

## ---- VA5-11: fit$plateau ---------------------------------------------------------------

test_that("VA5-11: a predictor on a plateau of its link is flagged in fit$plateau, without a warning", {
  ## the verifier's continuous carrier-5 fit: kappa -> 0 (carrier 5 -> GPD), log kappa -63
  set.seed(1)
  n <- 200
  d <- data.frame(x = runif(n), z = runif(n))
  d$yp <- rnbinom(n, mu = exp(1 + 0.8 * d$x), size = 3)
  d$yd <- rdiscegpd(n, sigma = exp(0.5 + 0.8 * d$x), xi = 0.2, kappa = 1.2)
  d$yc <- regpd(n, sigma = exp(0.5 + 0.8 * d$x), xi = 0.15, kappa = 1.5)
  ws <- character(0)
  f <- withCallingHandlers(
    egpd(list(yc ~ s(x, k = 6), ~ 1, ~ 1), data = d, family = "egpd",
         egpd.args = list(m = 5), restarts = FALSE),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_lt(f$coefficients[["logkappa_(Intercept)"]], -9)
  expect_identical(f$plateau, c(logscale = FALSE, shape = FALSE, logkappa = TRUE))
  expect_false(any(grepl("plateau", ws)))
  expect_output(print(f), "logkappa on a plateau")
  ## an ordinary fit of the same data (carrier 1)
  g <- egpd(list(yc ~ s(x, k = 6), ~ 1, ~ 1), data = d, family = "egpd")
  expect_identical(g$plateau, c(logscale = FALSE, shape = FALSE, logkappa = FALSE))
  expect_false(any(grepl("plateau", capture.output(print(g)))))
  ## a parameter held by `fixed` is never flagged, wherever it is held
  h <- egpd(list(yd ~ 1, ~ 1, ~ 1), data = d, family = "degpd",
            fixed = list(kappa = 1e-5), restarts = FALSE)
  expect_false(h$plateau[["logkappa"]])
  ## families the criterion was not designed for are not assessed
  p <- egpd(list(yp ~ x, ~ 1), data = d, family = "pig")
  expect_identical(p$plateau, c(logmu = NA, logsigma = NA))
})

## ---- VA5-10: the default start of the zero-truncated generalized Waring -----------------

test_that("VA5-10: egpd(family = 'ztgwaring') reaches the interior MLE from its default start", {
  ## C-opus's data: generalized Waring draws (mu = 4, k = 3, rho = 1.6), positive part
  set.seed(3)
  xz <- rgwaring(400, mu = 4, k = 3, rho = 1.6)
  xz <- xz[xz > 0]
  ## independent reference: the zero-truncated log-likelihood from Irwin's pmf in lgamma,
  ## maximised by BFGS with rho = 1 + exp(t3)
  lgw <- function(y, mu, k, rho) {
    a <- mu * (rho - 1) / k
    lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) - lgamma(a + k + rho + y) +
      lgamma(a + y) - lgamma(a) + lgamma(k + y) - lgamma(k) - lgamma(y + 1)
  }
  ztll <- function(th) {
    mu <- exp(th[1]); k <- exp(th[2]); rho <- 1 + exp(th[3])
    sum(lgw(xz, mu, k, rho)) - length(xz) * log(-expm1(lgw(0, mu, k, rho)))
  }
  o <- optim(c(log(4), log(3), log(0.6)), function(t) -ztll(t), method = "BFGS",
             control = list(maxit = 10000, reltol = 1e-15))
  expect_equal(-o$value, -563.632783664, tolerance = 1e-10)
  ## the estimate has a = k, where the outer-product Hessian is singular (round 6)
  expect_warning(f <- egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = xz), family = "ztgwaring"),
                 class = "egpd_singular_hessian")
  expect_identical(f$convergence, 0L)
  expect_equal(f$logLik, -o$value, tolerance = 1e-10)
  pr <- predict(f, newdata = data.frame(row = 1), type = "response")
  expect_equal(c(pr$mu, pr$k, pr$rho), c(exp(o$par[1:2]), 1 + exp(o$par[3])),
               tolerance = 1e-4)
  ## data whose likelihood rises towards rho -> 1 (infinite mean) keep the generic start,
  ## and the fit still reports that it stopped on the clamp
  set.seed(5)
  yh <- pmin(ceiling(runif(300)^(-1 / 0.3)), 1e9)
  expect_warning(fh <- egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = yh),
                            family = "ztgwaring"),
                 class = "egpd_nonconvergence")
  expect_identical(fh$convergence, 3L)
  expect_match(fh$convergence.report, "rho reached its lower bound")
})
