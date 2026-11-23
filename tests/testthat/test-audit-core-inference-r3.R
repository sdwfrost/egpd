# Round-3 core-inference fixes: hurdle models, prediction, intervals, coverage,
# summaries and distfit. References are independent closed forms, finite
# differences of them, or hand computations, never the package compared with itself.

## ---- independent references ------------------------------------------------

## Generalized Waring log-pmf, mean parameterisation a = mu (rho - 1) / k
gw_lpmf <- function(y, mu, k, rho) {
  ## Beyond rho ~ 1e8 the lgamma differences below cancel catastrophically, and the
  ## generalized Waring is within O(1/rho) of its negative-binomial limit NB(k, mu)
  ## (fits on light-tailed samples run to the C++ bound rho = 1e15).
  if (rho > 1e8) return(dnbinom(y, size = k, mu = mu, log = TRUE))
  a <- mu * (rho - 1) / k
  lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) - lgamma(a) - lgamma(k) +
    lgamma(a + y) + lgamma(k + y) - lgamma(a + k + rho + y) - lgamma(y + 1)
}

## GPD survival S(t), vectorised over all arguments
gpd_surv <- function(t, s, xi) {
  n <- max(length(t), length(s), length(xi))
  t <- rep_len(t, n); s <- rep_len(s, n); xi <- rep_len(xi, n)
  ifelse(abs(xi) < 1e-12, exp(-t / s), pmax(1 + xi * t / s, 0)^(-1 / xi))
}

## DEGPD model 1: F(t) = (1 - S(t))^kappa, P(Y = y) = F(y + 1) - F(y)
dgpd_cdf1 <- function(t, s, xi, k) exp(k * log1p(-gpd_surv(t, s, xi)))

## log P(Y = y) without cancellation:
## F(y + 1) - F(y) = F(y) expm1(kappa [log1p(-S(y + 1)) - log1p(-S(y))])
degpd1_lpmf <- function(y, s, xi, k) {
  a <- log1p(-gpd_surv(y, s, xi))
  b <- log1p(-gpd_surv(y + 1, s, xi))
  ifelse(y == 0, k * b, k * a + log(expm1(k * (b - a))))
}

## zero-truncated negative log-likelihoods on the natural scale
zt_gw_nll  <- function(y, mu, k, rho) -sum(gw_lpmf(y, mu, k, rho) - log1p(-exp(gw_lpmf(0, mu, k, rho))))
zt_deg_nll <- function(y, s, xi, k) -sum(degpd1_lpmf(y, s, xi, k) - log1p(-dgpd_cdf1(1, s, xi, k)))

## central-difference Hessian with one Richardson step (error O(h^4))
fd_hessian <- function(f, x, h = 1e-3) {
  k <- length(x)
  cd <- function(h) {
    H <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in i:k) {
      ei <- replace(numeric(k), i, h)
      ej <- replace(numeric(k), j, h)
      H[i, j] <- H[j, i] <- (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) +
                               f(x - ei - ej)) / (4 * h^2)
    }
    H
  }
  (4 * cd(h / 2) - cd(h)) / 3
}

gw_hurdle_data <- function() {
  set.seed(1)
  y <- rgwaring(1500, mu = 4, k = 2, rho = 4)
  y[y > 0][1:300]
}

## ---- E-04: failed hurdle fits are rejected; ztgwaring predictions work -------

test_that("E-04: fit_hurdle(gwaring) reaches the zero-truncated MLE, not the -1e20 sentinel", {
  ypos <- gw_hurdle_data()
  dz <- data.frame(y = c(rep(0L, 120), ypos))
  h <- fit_hurdle(y ~ 1, data = dz, family = "gwaring")

  o <- optim(c(log(4), 0, log(3)),
             function(th) zt_gw_nll(ypos, exp(th[1]), exp(th[2]), 1 + exp(th[3])),
             control = list(maxit = 5000, reltol = 1e-12))
  expect_equal(h$positive$convergence, 0)
  expect_equal(as.numeric(h$positive$logLik), -o$value, tolerance = 1e-6)
  pr <- predict(h$positive, type = "response")[1, ]
  expect_equal(c(pr$mu, pr$k, pr$rho), c(exp(o$par[1:2]), 1 + exp(o$par[3])),
               tolerance = 1e-3)
  ## the joint log-likelihood adds the binomial part: n0 log p0 + n1 log(1 - p0)
  p0 <- 120 / 420
  expect_equal(h$logLik, 120 * log(p0) + 300 * log1p(-p0) - o$value, tolerance = 1e-6)
})

test_that("E-04: fit_hurdle keeps only converged fits above the -1e20 sentinel", {
  set.seed(2)
  dz <- data.frame(y = c(rep(0L, 30), rgwaring(100, mu = 3, k = 2, rho = 4) + 1L))
  sentinel <- function(...) structure(list(logLik = -1e20, convergence = 0, df = 3),
                                      class = "egpd")
  local_mocked_bindings(egpd = sentinel, .package = "egpd")
  expect_error(fit_hurdle(y ~ 1, data = dz, family = "gwaring"),
               "failed from every start")

  ## a non-converged candidate is never preferred, even with a higher likelihood
  calls <- 0L
  mixed <- function(data, ...) {
    calls <<- calls + 1L
    if (calls == 1L)
      structure(list(logLik = -500, convergence = 1, df = 3, nobs = nrow(data)), class = "egpd")
    else
      structure(list(logLik = -700, convergence = 0, df = 3, nobs = nrow(data)), class = "egpd")
  }
  local_mocked_bindings(egpd = mixed, .package = "egpd")
  h <- fit_hurdle(y ~ 1, data = dz, family = "gwaring")
  expect_equal(as.numeric(h$positive$logLik), -700)
  expect_gte(calls, 2L)
})

test_that("E-04: ztgwaring predict() and rqresid() clamp rho exactly as the likelihood does", {
  ypos <- gw_hurdle_data()
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = ypos), family = "ztgwaring",
              inits = c(log(4), 0, log(2.5)))
  pr <- predict(fit, type = "response")[1, ]

  ## zero-truncated quantile from the closed form: smallest y >= 1 with
  ## (F(y) - f0) / (1 - f0) >= p
  zt_q <- function(p, mu, k, rho) {
    cdf <- cumsum(exp(gw_lpmf(0:5000, mu, k, rho)))
    f0 <- cdf[1]
    max(1, which(cdf >= f0 + p * (1 - f0) - 1e-12)[1] - 1)
  }
  q <- predict(fit, type = "quantile", prob = c(0.1, 0.5, 0.9))[1, ]
  expect_equal(unname(unlist(q)),
               vapply(c(0.1, 0.5, 0.9), zt_q, 0, mu = pr$mu, k = pr$k, rho = pr$rho))
  expect_true(all(is.finite(rqresid(fit, seed = 1))))

  ## a fit stuck at the start rho = exp(0) = 1 is clamped to the C++ bound 1 + 1e-6
  stuck <- fit
  stuck$logrho$coefficients[1] <- 0
  expect_equal(predict(stuck, type = "response")$rho[1], 1 + 1e-6)
  qs <- predict(stuck, type = "quantile", prob = 0.5)[1, 1]
  expect_equal(qs, zt_q(0.5, pr$mu, pr$k, 1 + 1e-6))
  expect_true(all(is.finite(rqresid(stuck, seed = 1))))
  rho_hi <- cf_bounds_cpp()[["gwaring"]][["rho_hi"]]
  stuck$logrho$coefficients[1] <- log(rho_hi) + 5
  expect_equal(predict(stuck, type = "response")$rho[1], rho_hi)
})

test_that("E-04: non-converged egpd() fits (logLik NA) are dropped with their warnings", {
  set.seed(2)
  dz <- data.frame(y = c(rep(0L, 30), rgwaring(100, mu = 3, k = 2, rho = 4) + 1L))
  nonconv <- function(msg)
    structure(class = c("egpd_nonconvergence", "warning", "condition"),
              list(message = msg, call = NULL))
  calls <- 0L
  mock <- function(data, ...) {
    calls <<- calls + 1L
    if (calls == 1L) {
      warning(nonconv("egpd: the optimiser did not converge"))
      structure(list(logLik = NA_real_, convergence = 2L, df = 3,
                     convergence.report = "step tolerance reached"), class = "egpd")
    } else {
      warning("a warning of the kept start")
      structure(list(logLik = -700, convergence = 0L, df = 3, nobs = nrow(data)), class = "egpd")
    }
  }
  local_mocked_bindings(egpd = mock, .package = "egpd")
  ws <- list()
  h <- withCallingHandlers(fit_hurdle(y ~ 1, data = dz, family = "gwaring"),
                           warning = function(w) {
                             ws[[length(ws) + 1L]] <<- w
                             invokeRestart("muffleWarning")
                           })
  expect_equal(as.numeric(h$positive$logLik), -700)
  expect_false(any(vapply(ws, inherits, logical(1), "egpd_nonconvergence")))
  expect_true(any(vapply(ws, function(w) grepl("kept start", conditionMessage(w)), logical(1))))

  ## every start fails: an error carrying egpd()'s convergence report, and no warnings
  local_mocked_bindings(egpd = function(...) {
    warning(nonconv("egpd: the optimiser did not converge"))
    structure(list(logLik = NA_real_, convergence = 3L, df = 3,
                   convergence.report = "iteration limit reached"), class = "egpd")
  }, .package = "egpd")
  ws <- list()
  expect_error(withCallingHandlers(fit_hurdle(y ~ 1, data = dz, family = "gwaring"),
                                   warning = function(w) {
                                     ws[[length(ws) + 1L]] <<- w
                                     invokeRestart("muffleWarning")
                                   }),
               "iteration limit reached")
  expect_length(ws, 0)
})

test_that("E-04: with no converged start, the best start is returned with one warning", {
  set.seed(2)
  dz <- data.frame(y = c(rep(0L, 30), rgwaring(100, mu = 3, k = 2, rho = 4) + 1L))
  run <- function(ll, code, warn) {
    fit <- structure(list(logLik = NA_real_, convergence = code, df = 3, tag = ll,
                          nobs = sum(dz$y > 0),
                          convergence.report = "inner Newton stopped, predicted decrease 0.01"),
                     class = "egpd")
    list(fit = fit, reason = sprintf("did not converge (convergence code %d)", code),
         conditions = list(simpleWarning(warn)), ll_coef = ll)
  }
  collect <- function(expr) {
    ws <- character()
    val <- withCallingHandlers(expr, warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = val, warnings = ws)
  }

  ## all three starts end with code 3: the highest log-likelihood at the coefficients wins
  res <- list(run(-210, 3L, "w1"), run(-207.8, 3L, "w2"), run(-208, 3L, "w3"))
  i <- 0L
  local_mocked_bindings(.hurdle_try_start = function(args) { i <<- i + 1L; res[[i]] },
                        .package = "egpd")
  out <- collect(fit_hurdle(y ~ 1, data = dz, family = "degpd"))
  h <- out$value
  expect_equal(i, 3L)
  expect_equal(h$positive$tag, -207.8)
  expect_false(h$converged)
  expect_true(is.na(h$logLik) && is.na(h$AIC) && is.na(as.numeric(logLik(h))))
  expect_equal(sum(grepl("no start of the zero-truncated degpd fit converged", out$warnings)), 1)
  expect_true(any(grepl("convergence code 3: inner Newton", out$warnings)))
  expect_true("w2" %in% out$warnings)
  expect_false(any(c("w1", "w3") %in% out$warnings))
  expect_output(print(h), "did not converge")

  ## a converged start is preferred to a higher non-converged one, without the warning
  conv <- run(-209, 0L, "wc")
  conv$fit$logLik <- -209
  conv$reason <- NULL
  res <- list(run(-207, 3L, "w1"), conv, run(-208, 3L, "w3"))
  i <- 0L
  out <- collect(fit_hurdle(y ~ 1, data = dz, family = "degpd"))
  expect_equal(out$value$positive$tag, -209)
  expect_true(out$value$converged)
  expect_false(any(grepl("no start", out$warnings)))

  ## an error only when every start errored or has no finite log-likelihood
  res <- list(list(fit = simpleError("boom"), reason = "error: boom", conditions = list(),
                   ll_coef = NA_real_),
              run(NA_real_, 3L, "w2"), run(NA_real_, 2L, "w3"))
  i <- 0L
  expect_error(collect(fit_hurdle(y ~ 1, data = dz, family = "degpd")),
               "failed from every start")
})

test_that("E-04: seed-25 DEGPD hurdle returns a fit whether or not a start converges", {
  skip_on_cran()
  set.seed(25)
  n <- 120
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)
  y <- ifelse(runif(n) < 0.35, 0L,
              rdiscegpd(n, sigma = exp(1 + 0.3 * x1), xi = 0.15, kappa = 1.2) + 1L)
  ws <- character()
  h <- withCallingHandlers(fit_hurdle(y ~ x1, data = data.frame(y, x1, x2), family = "degpd"),
                           warning = function(w) {
                             ws <<- c(ws, conditionMessage(w))
                             invokeRestart("muffleWarning")
                           })
  expect_s3_class(h, "hurdle.egpd")
  expect_true(is.finite(-h$positive$likfns$d0(unname(h$positive$coefficients),
                                              h$positive$likdata)))
  if (isTRUE(h$converged)) {
    expect_true(is.finite(h$logLik))
  } else {
    expect_true(is.na(h$logLik))
    expect_equal(sum(grepl("no start of the zero-truncated", ws)), 1)
  }
})

test_that("E-04: the zero-truncated objective and gradient stay accurate as f(0) -> 1", {
  ## log f(0) of the generalized Waring at (mu, k, rho), 256-bit references (Rmpfr)
  pts <- list(c(3, 5e7, 2), c(1e-3, 1e7, 3), c(0.5, 1e3, 1e5), c(2, 9e7, 50))
  ref <- c(-1.0382849543367559e-6, -3.0390623131640657e-9, -0.49739032262544589,
           -1.5694516079989598e-5)
  expect_equal(vapply(pts, function(p) egpd:::.zt_gwaring_logf0(as.list(log(p))), 0),
               ref, tolerance = 1e-9)

  ## zero-truncated Waring negative log-likelihood and gradient at k -> Inf, where
  ## 1 - f(0) = 3.0e-9, 8.8e-9 and 3.0e-13 (Rmpfr references)
  y <- c(1, 1, 2, 3, 5, 10, 40)
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "ztgwaring",
              inits = c(log(4), 0, log(2.5)))
  etas <- list(c(log(1e-3), log(1e7), log(3)), c(log(0.05), log(5e7), log(1.5)),
               c(log(1e-7), log(1e7), log(3)))
  nll_ref <- c(28.439589051629185, 29.504073319022401, 28.439589043344656)
  grad_ref <- rbind(c(8.28535705625e-09, 0.460649686087, -0.545782194659),
                    c(2.50339339343e-08, 0.395678553363, -0.554822480789),
                    c(8.28535705064e-13, 0.460649693671, -0.545782206276))
  for (i in seq_along(etas)) {
    expect_equal(fit$likfns$d0(etas[[i]], fit$likdata), nll_ref[i], tolerance = 1e-7)
    expect_equal(colSums(fit$likfns$d120(etas[[i]], fit$likdata)[, 1:3]), grad_ref[i, ],
                 tolerance = 1e-5)
  }

  ## zero-truncated DEGPD with 1 - f(0) = 9e-14 (sigma = 0.005, xi = 0.1, kappa = 1.5)
  fd <- suppressMessages(suppressWarnings(
    egpd(list(y ~ 1, ~1, ~1), data.frame(y = c(1, 1, 2)), family = "ztdegpd",
         inits = c(0, log(0.1), log(1.5)))))
  expect_equal(fd$likfns$d0(c(log(0.005), log(0.1), log(1.5)), fd$likdata),
               6.7119790105821925, tolerance = 1e-9)
})

## ---- U3-06: the DEGPD hurdle start respects the shape link -------------------

test_that("U3-06: fit_hurdle(degpd) finds the zero-truncated MLE under every shape link", {
  set.seed(11)
  y <- rdiscegpd(500, kappa = 1.5, sigma = 3, xi = 0.2)
  ypos <- y[y > 0]
  ref <- optim(c(log(3), log(0.2), log(1.5)),
               function(th) zt_deg_nll(ypos, exp(th[1]), exp(th[2]), exp(th[3])),
               control = list(maxit = 5000, reltol = 1e-12))
  xi_ref <- exp(ref$par[2])
  expect_equal(-ref$value, -1048.777, tolerance = 1e-6)

  h <- fit_hurdle(y ~ 1, data = data.frame(y = y), family = "degpd")
  expect_equal(as.numeric(h$positive$logLik), -ref$value, tolerance = 1e-6)
  expect_equal(predict(h$positive, type = "response")$shape[1], xi_ref, tolerance = 1e-3)

  hid <- fit_hurdle(y ~ 1, data = data.frame(y = y), family = "degpd",
                    degpd.args = list(link = "identity"))
  expect_equal(as.numeric(hid$positive$logLik), -ref$value, tolerance = 1e-6)
  expect_equal(unname(coef(hid$positive)[2]), xi_ref, tolerance = 1e-3)

  hb <- fit_hurdle(y ~ 1, data = data.frame(y = y), family = "degpd",
                   degpd.args = list(xi.max = 0.5))
  expect_equal(as.numeric(hb$positive$logLik), -ref$value, tolerance = 1e-6)
  expect_equal(predict(hb$positive, type = "response")$shape[1], xi_ref, tolerance = 1e-3)
})

## ---- U3-07: exact zero-truncated Hessian; hurdle methods ---------------------

test_that("U3-07: zero-truncated DEGPD Hessian is the observed information; Waring is BHHH", {
  ## per-observation Hessian columns (upper triangle, row-major) vs differences of the
  ## independent zero-truncated negative log-likelihood in the linear predictors
  ut <- function(H) H[upper.tri(H, diag = TRUE)][c(1, 2, 4, 3, 5, 6)]
  check_rows <- function(fit, nll_eta, rows) {
    beta <- fit$coefficients
    d12 <- fit$likfns$d120(beta, fit$likdata)
    for (j in rows) {
      H <- fd_hessian(function(eta) nll_eta(fit$likdata$y[j, 1], eta), beta)
      expect_equal(d12[j, 4:9], ut(H), tolerance = 1e-6)
    }
  }
  set.seed(11)
  yd <- rdiscegpd(500, kappa = 1.5, sigma = 3, xi = 0.2)
  yd <- yd[yd > 0]
  rows <- c(1, 2, which.max(yd))
  fl <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = yd), family = "ztdegpd",
             inits = c(1, log(0.2), 0.25))
  check_rows(fl, function(y, e) zt_deg_nll(y, exp(e[1]), exp(e[2]), exp(e[3])), rows)
  fi <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = yd), family = "ztdegpd",
             inits = c(1, 0.2, 0.25), degpd.args = list(link = "identity"))
  check_rows(fi, function(y, e) zt_deg_nll(y, exp(e[1]), e[2], exp(e[3])), rows)
  fb <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = yd), family = "ztdegpd",
             inits = c(1, -1.35, 0.25), degpd.args = list(xi.max = 0.5))
  check_rows(fb, function(y, e) zt_deg_nll(y, exp(e[1]), -0.5 * expm1(-exp(e[2]) / 0.5),
                                           exp(e[3])), rows)
  ## and the positive-part covariance is the inverse observed information
  H <- fd_hessian(function(b) zt_deg_nll(yd, exp(b[1]), exp(b[2]), exp(b[3])),
                  fl$coefficients)
  expect_equal(unname(fl$Vp), unname(solve(H)), tolerance = 1e-5)

  ## zero-truncated Waring: documented BHHH, the outer product of the exact score
  ypos <- gw_hurdle_data()
  fg <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = ypos), family = "ztgwaring",
             inits = c(log(4), 0, log(2.5)))
  d12 <- fg$likfns$d120(fg$coefficients, fg$likdata)
  for (j in c(1, 2, which.max(ypos))) {
    f <- function(e) zt_gw_nll(ypos[j], exp(e[1]), exp(e[2]), exp(e[3]))
    g <- vapply(1:3, function(i) {
      h <- replace(numeric(3), i, 1e-5)
      (f(fg$coefficients + h) - f(fg$coefficients - h)) / 2e-5
    }, 0)
    expect_equal(d12[j, 1:3], g, tolerance = 1e-6)
    expect_equal(d12[j, 4:9], ut(outer(g, g)), tolerance = 1e-6)
  }
})

test_that("U3-07: hurdle predict/summary/print/logLik/nobs methods", {
  set.seed(5)
  n <- 400
  x <- runif(n)
  ypos <- rgwaring(n, mu = exp(1 + 0.5 * x), k = 2, rho = 4) + 1L
  y <- ifelse(runif(n) < plogis(-0.5 + x), 0L, ypos)
  d <- data.frame(y = y, x = x)
  h <- fit_hurdle(list(y ~ x, ~1, ~1), data = d, family = "gwaring")

  ll <- logLik(h)
  expect_equal(as.numeric(ll), as.numeric(logLik(h$zero)) + as.numeric(h$positive$logLik))
  expect_identical(nobs(h), as.integer(n))
  expect_equal(attr(ll, "df"), sum(h$zero$edf) + h$positive$df)

  nd <- data.frame(x = c(0.1, 0.5, 0.9))
  bz <- coef(h$zero)
  p0 <- plogis(bz[1] + bz[2] * nd$x, lower.tail = FALSE)
  bp <- coef(h$positive)
  mu <- exp(bp[1] + bp[2] * nd$x); k <- exp(bp[3]); rho <- exp(bp[4])

  r <- predict(h, newdata = nd)
  expect_equal(r$pzero, unname(p0), tolerance = 1e-10)
  expect_equal(r$mu, unname(mu), tolerance = 1e-10)
  lk <- predict(h, newdata = nd, type = "link")
  expect_equal(lk$logitppositive, unname(bz[1] + bz[2] * nd$x), tolerance = 1e-10)

  ## hurdle quantile from the closed-form cdf F(0) = p0, F(y) = p0 + (1 - p0) F_ZT(y)
  hq <- function(p, i) {
    pm <- exp(gw_lpmf(0:5000, mu[i], k, rho))
    cdf <- p0[i] + (1 - p0[i]) * cumsum(c(0, pm[-1])) / (1 - pm[1])
    which(cdf >= p - 1e-12)[1] - 1
  }
  probs <- c(0.05, 0.3, 0.5, 0.9, 0.99)
  q <- predict(h, newdata = nd, prob = probs)
  expect_equal(names(q), paste("q", probs, sep = ":"))
  for (i in seq_len(nrow(nd)))
    expect_equal(unname(unlist(q[i, ])), vapply(probs, hq, 0, i = i))
  expect_equal(nrow(predict(h)), n)

  s <- summary(h)
  expect_s3_class(s, "summary.hurdle.egpd")
  expect_s3_class(s$zero, "summary.gam")
  expect_s3_class(s$positive, "summary.egpd")
  expect_output(print(s), "Zero part")
  expect_output(print(h), "Hurdle model")
})

## ---- E-05: two-sided parametric p-values ------------------------------------

test_that("E-05: summary.egpd reports two-sided Wald p-values (independent Hessian)", {
  set.seed(5)
  y <- rdiscegpd(400, sigma = 3, xi = 0.1, kappa = 1.5)
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "degpd")
  b <- fit$coefficients
  nll <- function(b) -sum(degpd1_lpmf(y, exp(b[1]), exp(b[2]), exp(b[3])))
  z <- b / sqrt(diag(solve(fd_hessian(nll, b))))
  tab <- do.call(rbind, summary(fit)[[1]])
  expect_named(tab, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(tab[["z value"]], unname(z), tolerance = 1e-4)
  ## |z| = 1.63 and 1.84 here: two-sided p = 0.103 and 0.066 (one-sided 0.052 / 0.033)
  expect_equal(tab[["Pr(>|z|)"]], unname(2 * pnorm(-abs(z))), tolerance = 1e-4)
  expect_true(any(tab[["Pr(>|z|)"]] > 0.05 & tab[["Pr(>|z|)"]] < 0.11))
})

## ---- E-06: confint link handling -------------------------------------------

test_that("E-06: confint reports slopes on the link scale and parameters only where defined", {
  set.seed(50)
  d <- data.frame(y = rdiscegpd(200, sigma = 3, xi = 0.2, kappa = 2, type = 1),
                  x = rnorm(200))
  fit <- egpd(list(lsigma = y ~ x, lxi = ~1, lkappa = ~1), data = d, family = "degpd")
  b <- fit$coefficients
  V <- fit$Vc
  z <- qnorm(0.975)
  wald <- function(k) unname(b[k] + c(-1, 1) * z * sqrt(V[k, k]))

  ## a slope: link scale, labelled by its coefficient name (not exp() and "scale")
  ci <- confint(fit, parm = "logscale_x")
  expect_identical(rownames(ci), "logscale_x")
  expect_equal(unname(ci[1, ]), wald(2))
  expect_equal(unname(confint(fit, parm = 2)), unname(ci))

  ## default: the intercept of a predictor with covariates is a coefficient, the
  ## intercept-only parameters are natural-scale
  ci0 <- confint(fit)
  expect_identical(rownames(ci0), c("logscale_(Intercept)", "shape", "kappa"))
  expect_equal(unname(ci0[1, ]), wald(1))
  expect_equal(unname(ci0[2, ]), exp(wald(3)))
  expect_equal(unname(ci0[3, ]), exp(wald(4)))

  ## parameter-level intervals at newdata: exp(eta +- z se(eta))
  nd <- data.frame(x = c(-1, 0, 2))
  cin <- confint(fit, newdata = nd, parm = c("scale", "xi"))
  expect_identical(rownames(cin), c(paste0("scale[", 1:3, "]"), paste0("shape[", 1:3, "]")))
  Xs <- cbind(1, nd$x)
  eta <- drop(Xs %*% b[1:2])
  se <- sqrt(rowSums((Xs %*% V[1:2, 1:2]) * Xs))
  expect_equal(unname(cin[1:3, "lower"]), exp(eta - z * se))
  expect_equal(unname(cin[1:3, "upper"]), exp(eta + z * se))
  expect_equal(unname(cin[4:6, "upper"]), rep(exp(wald(3))[2], 3))
  expect_error(confint(fit, newdata = nd, method = "profile"), "wald")

  ## a predictor without an intercept has no parameter-level interval
  fni <- egpd(list(lsigma = y ~ x - 1, lxi = ~1, lkappa = ~1), data = d, family = "degpd")
  expect_error(confint(fni, parm = "scale"), "no intercept")
  expect_identical(rownames(confint(fni, parm = "logscale_x")), "logscale_x")
})

test_that("E-06: logit-link slope intervals are not passed through plogis", {
  set.seed(3)
  n <- 500
  z <- rnorm(n)
  y <- rzidiscegpd(n, pi = plogis(-1 + 0.8 * z), sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  fit <- egpd(list(y ~ 1, ~1, ~1, ~z), data = data.frame(y = y, z = z), family = "zidegpd",
              inits = c(1, log(0.2), log(1.5), -1))
  k <- which(names(fit$coefficients) == "logitpi_z")
  ci <- confint(fit, parm = k)
  expect_identical(rownames(ci), "logitpi_z")
  expect_equal(unname(ci[1, ]),
               unname(fit$coefficients[k] + c(-1, 1) * qnorm(0.975) * sqrt(fit$Vc[k, k])))
  expect_lt(ci[1, "lower"], 0.8)
  expect_gt(ci[1, "upper"], 0.8)
})

## ---- CI-1: profile intervals that reach the parameter-space boundary ---------

test_that("CI-1: a profile that levels off below the cutoff gives the link-scale limit", {
  cutoff <- qchisq(0.95, 1)
  plateau <- function(v) 2 * (1 - exp(-abs(v))) - cutoff  # deviance -> 2 < cutoff
  expect_identical(.find_profile_root(plateau, 0, 1, -1, cutoff), -Inf)
  expect_identical(.find_profile_root(plateau, 0, 1, +1, cutoff), Inf)
  ## still rising 96 SEs out: not a boundary, so the old warning and NA remain
  rising <- function(v) 3.5 * abs(v) / 96 - cutoff
  expect_warning(r <- .find_profile_root(rising, 0, 1, +1, cutoff), "Could not bracket")
  expect_identical(r, NA_real_)
  expect_equal(.find_profile_root(function(v) v^2 - cutoff, 0, 1, +1, cutoff),
               sqrt(cutoff), tolerance = 1e-4)
})

test_that("CI-1: profile confint reaches xi = 0, pi = 0 and xi.max instead of returning NA", {
  skip_on_cran()
  cutoff <- qchisq(0.95, 1)
  fit_m1 <- function(y, ...) suppressMessages(egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1),
                                                   data = data.frame(y = y), family = "degpd", ...))
  ## parameter-coverage vignette settings (xi = 0.1, n = 500): xi_hat = 0.052
  set.seed(5)
  y <- rdiscegpd(500, sigma = 3, xi = 0.1, kappa = 2, type = 1)
  fit <- fit_m1(y, degpd.args = list(m = 1))
  expect_no_warning(ci <- confint(fit, parm = "shape", method = "profile"))
  expect_identical(ci[1, "lower"], 0)
  expect_true(is.finite(ci[1, "upper"]) && ci[1, "upper"] > exp(fit$coefficients[[2]]))
  ## the xi -> 0 limit is inside the confidence set: profile deviance at xi = xi_hat * e^-40
  b <- fit$coefficients
  opt <- optim(b[-2], .nllh.nopen, likdata = fit$likdata, likfns = fit$likfns,
               id = 2L, vals = b[[2]] - 40, method = "BFGS")
  expect_lt(2 * (opt$value - .nllh.nopen(b, fit$likdata, fit$likfns)), cutoff)

  ## logit link: no zero inflation in the data, so the interval for pi reaches 0
  set.seed(1)
  y <- rdiscegpd(500, sigma = 3, xi = 0.1, kappa = 2, type = 1)
  fit <- suppressMessages(egpd(list(lsigma = y ~ 1, lxi = ~1, lkappa = ~1, logitpi = ~1),
                               data = data.frame(y = y), family = "zidegpd",
                               zidegpd.args = list(m = 1)))
  expect_no_warning(ci <- confint(fit, parm = "pi", method = "profile"))
  expect_identical(ci[1, "lower"], 0)
  expect_true(ci[1, "upper"] > 0 && ci[1, "upper"] < 0.1)

  ## bounded shape link: an interior xi_hat whose interval reaches both 0 and xi.max
  set.seed(2)
  y <- rdiscegpd(500, sigma = 3, xi = 0.1, kappa = 2, type = 1)
  fit <- fit_m1(y, degpd.args = list(m = 1, xi.max = 0.2))
  expect_true(fit$coefficients[[2]] > -3 && fit$coefficients[[2]] < -1)
  expect_no_warning(ci <- confint(fit, parm = "shape", method = "profile"))
  expect_identical(unname(ci[1, ]), c(0, 0.2))
})

## ---- E-09 / E-10: gamlss and bamlss predictive coverage ---------------------

test_that("E-09: gamlss coverage intervals are integer for every discrete family type", {
  skip_if_not_installed("gamlss")
  set.seed(21)
  n <- 300
  x <- runif(n)
  y <- rpois(n, exp(1 + x))
  g <- gamlss::gamlss(y ~ x, family = gamlss.dist::PO(), data = data.frame(y = y, x = x),
                      trace = FALSE)
  for (m in c("plug-in", "parametric"))
    for (wrong in list(rep(y, 2), y[-1L]))
      expect_error(predictive_coverage(g, y = wrong, nsim = 5, method = m),
                   "length\\(y\\) must match")
  set.seed(1)
  r <- predictive_coverage(g, y = y, nsim = 2000, method = "plug-in")
  expect_true(all(r$L == floor(r$L)) && all(r$U == floor(r$U)))
  ## empirical type-1 quantiles of 2000 Poisson draws: mostly the exact quantiles
  mu <- as.numeric(fitted(g, "mu"))
  expect_gt(mean(r$L == qpois(0.025, mu)), 0.8)
  expect_gt(mean(r$U == qpois(0.975, mu)), 0.8)
  rp <- predictive_coverage(g, y = y, nsim = 500, method = "parametric")
  expect_true(all(rp$L == floor(rp$L)) && all(rp$U == floor(rp$U)))

  ## the package's own discrete gamlss family (not named D*/ZID*)
  x2 <- x[1:100]
  y2 <- rpig(100, mu = exp(1 + x2), sigma = 0.5)
  g2 <- gamlss::gamlss(y ~ x, family = PIG(), data = data.frame(y = y2, x = x2), trace = FALSE)
  r2 <- predictive_coverage(g2, y = y2, nsim = 400, method = "plug-in")
  expect_true(all(r2$L == floor(r2$L)) && all(r2$U == floor(r2$U)))
})

test_that("E-10: bamlss parametric coverage draws from the posterior samples", {
  skip_if_not_installed("bamlss")
  skip_on_cran()
  ## bamlss() looks its default optimizer and sampler up on the search path
  if (!"package:bamlss" %in% search()) {
    suppressPackageStartupMessages(library(bamlss))
    on.exit(detach("package:bamlss", unload = FALSE, character.only = TRUE), add = TRUE)
  }
  set.seed(8)
  n <- 150
  x <- runif(n)
  y <- rdiscegpd(n, sigma = exp(1 + x), xi = 0.2, kappa = 1.5, type = 1)
  d <- data.frame(y = y, x = x)
  fam <- degpd_bamlss(m = 1)
  b0 <- suppressWarnings(bamlss::bamlss(list(y ~ x, ~1, ~1), data = d, family = fam,
                                        sampler = FALSE, verbose = FALSE))
  expect_error(predictive_coverage(b0, y = y, nsim = 100, method = "parametric"),
               "posterior samples")

  b <- suppressWarnings(bamlss::bamlss(list(y ~ x, ~1, ~1), data = d, family = fam,
                                       n.iter = 300, burnin = 100, thin = 2,
                                       verbose = FALSE))
  S <- as.matrix(b$samples)
  for (m in c("plug-in", "parametric")) {
    expect_error(predictive_coverage(b, y = rep(y, 2), nsim = 5, method = m),
                 "length\\(y\\) must match")
    expect_length(predictive_coverage(b, y = y, nsim = 5, method = m)$covered, n)
  }
  nd <- data.frame(x = c(0, 0.5, 1))
  ## one draw per stored sample: log sigma reproduces b0 + b1 x for every sample
  eta <- S[, "sigma.p.(Intercept)"] + outer(S[, "sigma.p.x"], nd$x)
  ## drawn row by row through the recorded sample indices, so a different pairing of
  ## the samples across parameters would fail
  set.seed(1)
  dr <- egpd:::.bamlss_param_draws(b, nd, fam$names, nrow(S))
  di <- attr(dr, "draws")
  expect_setequal(di, seq_len(nrow(S)))
  for (j in seq_len(nrow(nd)))
    expect_equal(log(dr$sigma[, j]), unname(eta[di, j]), tolerance = 1e-8)
  expect_equal(log(dr$xi[, 1]), unname(S[di, "xi.p.(Intercept)"]), tolerance = 1e-8)
  expect_equal(log(dr$kappa[, 1]), unname(S[di, "kappa.p.(Intercept)"]), tolerance = 1e-8)

  set.seed(2)
  r <- predictive_coverage(b, y = c(5, 10, 20), newdata = nd, nsim = 400,
                           method = "parametric")
  expect_true(all(is.finite(c(r$L, r$U))) && all(r$L <= r$U))
  expect_true(all(r$L == floor(r$L)))
})

test_that("P3: predictive_coverage on a removeData fit asks for newdata", {
  set.seed(31)
  d <- data.frame(y = rdiscegpd(200, sigma = 3, xi = 0.2, kappa = 1.5), x = runif(200))
  fit <- egpd(list(y ~ x, ~1, ~1), data = d, family = "degpd", removeData = TRUE)
  expect_error(predictive_coverage(fit, y = d$y, nsim = 50), "removeData")
  set.seed(1)
  r <- predictive_coverage(fit, y = d$y[1:5], newdata = d[1:5, ], nsim = 200)
  expect_length(r$covered, 5)
})

## ---- distfit: support, infeasible fits, convergence, intervals ----------------

test_that("distfit: out-of-support observations have zero likelihood, not exp(-708)", {
  spec <- egpd:::.distfit_parspec(1, "degpd")
  th <- egpd:::.par_to_theta(c(sigma = 2, xi = -0.5, kappa = 1), spec)
  ## xi = -0.5, sigma = 2: the discrete support is {0, 1, 2, 3}
  expect_gte(egpd:::.distfit_nll(th, x = c(1, 10), type = 1, family = "degpd",
                                 fix.arg = NULL, spec = spec), 1e20)
  spec_zi <- egpd:::.distfit_parspec(1, "zidegpd")
  th_zi <- egpd:::.par_to_theta(c(sigma = 2, xi = -0.5, kappa = 1, pi = 0.2), spec_zi)
  expect_gte(egpd:::.distfit_nll(th_zi, x = c(0, 10), type = 1, family = "zidegpd",
                                 fix.arg = NULL, spec = spec_zi), 1e20)
  ## inside the support the value is the closed-form negative log-likelihood
  expect_equal(egpd:::.distfit_nll(th, x = c(1, 3), type = 1, family = "degpd",
                                   fix.arg = NULL, spec = spec),
               -sum(degpd1_lpmf(c(1, 3), 2, -0.5, 1)))
  ## no feasible parameter value at all: an error, not a fit at log-likelihood -1e20
  expect_error(distfit(c(0, 1, 2, 10), family = "degpd", fix.arg = list(sigma = 2, xi = -0.5)),
               "positive likelihood")
  ## a start whose support excludes data is first moved into the support
  set.seed(12)
  x <- rdiscegpd(300, sigma = 2, xi = 0.1, kappa = 1.5)
  f0 <- distfit(x, family = "degpd")
  f1 <- distfit(x, family = "degpd", start = list(sigma = 2, xi = -0.2, kappa = 1.5))
  expect_gt(max(x), 2 / 0.2)
  expect_equal(f1$loglik, f0$loglik, tolerance = 1e-6)
})

test_that("distfit: continuous families reject negative data like the discrete ones", {
  expect_error(distfit(c(-1, 2, 3), family = "egpd"), "non-negative")
  expect_error(distfit(c(-1, 0, 2, 3), family = "ziegpd"), "non-negative")
  expect_error(distfit(c(-1, 0, 2, 3), family = "cpegpd"), "non-negative")
  expect_error(distfit(c(0, 1.5, 2), family = "degpd"), "non-negative integer")
})

test_that("distfit: an optimiser stopped by its iteration limit is reported", {
  set.seed(3)
  x <- rgwaring(200, mu = 4, k = 2, rho = 4)
  expect_warning(distfit(x, family = "gwaring", control = list(maxit = 5)),
                 "did not converge")
})

test_that("distfit: the GPD start is the method-of-moments estimator", {
  ## GPD: m = sigma / (1 - xi), v = sigma^2 / ((1 - xi)^2 (1 - 2 xi)), so
  ## xi = (1 - m^2 / v) / 2 and sigma = m (1 - xi)
  mom <- function(m, v) {
    xi <- min(max((1 - m^2 / v) / 2, -0.4), 0.45)
    c(xi = xi, sigma = m * (1 - xi))
  }
  start_of <- function(x, fam) {
    st <- egpd:::.distfit_start(x, 1, fam, egpd:::.distfit_parspec(1, fam))
    c(xi = st$xi, sigma = st$sigma)
  }
  x <- c(1, 1, 1, 1, 16)                       # m = 4, v = 45: heavy-tailed, xi > 0
  expect_equal(start_of(x, "egpd"), mom(4, 45))
  expect_equal(start_of(x, "egpd")[["xi"]], (1 - 16 / 45) / 2)
  ## counts: X = Y + U, U uniform: mean + 1/2, variance + 1/12; DEGPD keeps its zeros
  y <- c(0, 0, 1, 2, 5, 12)
  expect_equal(start_of(y, "degpd"), mom(mean(y) + 0.5, var(y) + 1 / 12))
  ## ZIDEGPD: positive counts are X | X >= 1, whose excess over 1 is GPD(sigma + xi, xi)
  z <- c(0, 0, 0, 1, 2, 5, 12)
  e <- mom(mean(z[z > 0]) - 0.5, var(z[z > 0]) + 1 / 12)
  expect_equal(start_of(z, "zidegpd"), c(xi = e[["xi"]], sigma = e[["sigma"]] - e[["xi"]]))
  ## ZIEGPD: the positive values are draws of X
  w <- c(0, 0, 1, 1, 1, 1, 16)
  expect_equal(start_of(w, "ziegpd"), mom(4, 45))
  ## xi is clipped to [-0.4, 0.45], keeping sigma = m (1 - xi) positive
  expect_equal(start_of(c(rep(1, 99), 1e5), "egpd")[["xi"]], 0.45)
  lt <- start_of(c(9.9, 10, 10.1, 10), "egpd")
  expect_equal(lt[["xi"]], -0.4)
  expect_equal(lt[["sigma"]], 10 * 1.4)
})

test_that("distfit: confint is a back-transformed Wald interval inside the parameter space", {
  set.seed(1)
  x <- regpd(30, sigma = 1, xi = 0.3, kappa = 0.5)
  f <- distfit(x, family = "egpd")
  ci <- confint(f)
  z <- qnorm(0.975)
  est <- f$estimate
  sd <- f$sd
  ## log-scale parameters: est * exp(+- z sd / est); xi is on its own (identity) scale
  for (nm in c("sigma", "kappa"))
    expect_equal(unname(ci[nm, ]), unname(est[nm] * exp(c(-1, 1) * z * sd[nm] / est[nm])))
  expect_equal(unname(ci["xi", ]), unname(est["xi"] + c(-1, 1) * z * sd["xi"]))
  ## the natural-scale Wald lower limit of sigma would be negative here
  expect_lt(est[["sigma"]] - z * sd[["sigma"]], 0)
  expect_true(all(ci[c("sigma", "kappa"), ] > 0))

  set.seed(2)
  xz <- rzidiscegpd(300, pi = 0.3, sigma = 2, xi = 0.2, kappa = 1.2)
  fz <- distfit(xz, family = "zidegpd")
  p <- fz$estimate[["pi"]]
  expect_equal(unname(confint(fz, "pi")[1, ]),
               plogis(qlogis(p) + c(-1, 1) * z * fz$sd[["pi"]] / (p * (1 - p))))
})

test_that("distfit Bernstein: quantiles keep their relative accuracy in both tails", {
  ## uniform weights make the Bernstein density 1, so the model is EGPD type 1 with
  ## Q(p) = Q_GPD(p^(1/kappa)); 256-bit references (Rmpfr) for sigma = 2, kappa = 1.3
  ps <- c(1e-30, 1e-20, 1e-15, 1e-12, 0.5, 1 - 1e-9, 1 - 1e-12)
  ref <- list(
    `0.2` = c(1.6753552801365867e-23, 8.2492527658027153e-16, 5.7885322494435595e-12,
              1.1756032148695963e-9, 1.9330997869630961, 654.94951395561584,
              2637.2233930359502),
    `0` = c(1.6753552801365867e-23, 8.2492527658027150e-16, 5.7885322494418841e-12,
            1.1756032148004941e-9, 1.7673094062810598, 41.971260259160899,
            55.786815004721469),
    `-0.25` = c(1.6753552801365867e-23, 8.2492527658027146e-16, 5.7885322494397900e-12,
                1.1756032147141164e-9, 1.5857128732768179, 7.9578687692688387,
                7.9925079313692950))
  w <- rep(1 / 8, 8)
  for (xi in c(0.2, 0, -0.25)) {
    q <- egpd:::.bernstein_full_quantile(ps, 2, xi, 1.3, w, 8)
    expect_lt(max(abs(q / ref[[as.character(xi)]] - 1)), 1e-12)
  }
  expect_equal(egpd:::.bernstein_full_quantile(c(0, 1, NA), 2, -0.25, 1.3, w, 8),
               c(0, 8, NA))
  ## non-uniform weights: the quantile inverts the cdf, also deep in the lower tail
  set.seed(3)
  wr <- exp(rnorm(8)); wr <- wr / sum(wr)
  pm <- c(1e-20, 0.01, 0.2, 0.5, 0.8, 0.99)
  qm <- egpd:::.bernstein_full_quantile(pm, 2, 0.2, 1.3, wr, 8)
  expect_lt(max(abs(egpd:::.bernstein_full_cdf(qm, 2, 0.2, 1.3, wr, 8) / pm - 1)), 1e-10)
})

## ---- P3: simulate(), fitted(), plot(), print.summary() -----------------------

set_coefs <- function(fit, values) {
  for (i in seq_len(fit$ngam)) {
    fit[[i]]$coefficients[] <- values[[i]]
    fit$coefficients[fit$idpars == i] <- values[[i]]
  }
  fit
}

test_that("P3: simulate(type = 'response') applies the predict() clamps", {
  ypos <- gw_hurdle_data()
  ## the estimate has a = k, where the outer-product Hessian is singular (round 6)
  expect_warning(fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = ypos), family = "gwaring",
                             inits = c(log(3), 0, log(4))),
                 class = "egpd_singular_hessian")
  rho_hi <- cf_bounds_cpp()[["gwaring"]][["rho_hi"]]
  fit <- set_coefs(fit, list(log(3), log(2), log(rho_hi) + 5))
  fit$Vp <- fit$Vc <- diag(1e-8, 3)
  nd <- data.frame(row = 1:2)
  draws <- simulate(fit, nsim = 20, seed = 1, newdata = nd, type = "response")
  ## the likelihood caps rho at its C++ bound (cf_bounds_cpp), and so does predict()
  expect_equal(predict(fit, newdata = nd, type = "response")$rho, c(rho_hi, rho_hi))
  expect_true(all(draws$rho == rho_hi))
  expect_equal(unname(draws$mu), matrix(3, 2, 20), tolerance = 1e-3)
})

test_that("P3: simulate(probs =) returns quantiles of each coefficient draw", {
  set.seed(7)
  y <- rdiscegpd(300, sigma = 3, xi = 0.2, kappa = 1.5)
  fit <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "degpd")
  nd <- data.frame(row = 1:3)
  pars <- simulate(fit, nsim = 25, seed = 3, newdata = nd, type = "response")
  q <- simulate(fit, nsim = 25, seed = 3, newdata = nd, probs = c(0.5, 0.99),
                threshold = 10)
  expect_named(q, c("q:0.5", "q:0.99"))
  ## DEGPD quantile: the smallest y with F(y + 1) >= p, F(t) = H(t)^kappa, so
  ## y = ceiling(H^{-1}(p^(1/kappa))) - 1
  qref <- function(p, s, xi, k) ceiling(s / xi * ((1 - p^(1 / k))^(-xi) - 1)) - 1
  for (p in c(0.5, 0.99)) {
    ref <- qref(p, pars$scale, pars$shape, pars$kappa)
    expect_equal(unname(q[[paste0("q:", p)]]), unname(ref) + 10)
  }
  expect_equal(dim(q[[1]]), c(3L, 25L))
})

test_that("P3: fitted() returns link (default) or response-scale parameters", {
  set.seed(7)
  d <- data.frame(y = rdiscegpd(200, sigma = 3, xi = 0.2, kappa = 1.5), x = runif(200))
  fit <- egpd(list(y ~ x, ~1, ~1), data = d, family = "degpd")
  fl <- fitted(fit)
  expect_named(fl, c("logscale", "logshape", "logkappa"))
  b <- fit$coefficients
  expect_equal(fl$logscale, unname(b[1] + b[2] * d$x))
  fr <- fitted(fit, type = "response")
  expect_equal(fr$scale, exp(fl$logscale))
  expect_equal(fr$shape, exp(fl$logshape))
  expect_equal(fr$kappa, exp(fl$logkappa))
  expect_equal(fr$xi, fr$shape)
})

test_that("P3: plot.egpd and print.summary.egpd run (smoke)", {
  set.seed(9)
  d <- data.frame(x = runif(200))
  d$y <- rdiscegpd(200, sigma = exp(1 + sin(2 * pi * d$x)), xi = 0.2, kappa = 1.5)
  fit <- egpd(list(y ~ s(x, k = 5), ~1, ~1), data = d, family = "degpd")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, onepage = FALSE, which = 1, main = "s(x)"))
  s <- summary(fit)
  expect_output(print(s), "Smooth terms")
  expect_named(s[[2]][[1]], c("edf", "max.df", "Chi.sq", "p-value"))
  expect_true(s[[2]][[1]][["p-value"]] < 0.01)
})

## ---- rqresid for further families (deterministic, same uniforms) -------------

## rqresid(seed = s) draws u = F(y - 1) + U (F(y) - F(y - 1)) with U = runif(n)
## after set.seed(s); rebuild it from independent cdfs
rq_ref <- function(Fhi, Flo, seed, n) {
  set.seed(seed)
  qnorm(Flo + (Fhi - Flo) * runif(n))
}

test_that("rqresid: DEGPD, ZIDEGPD, Waring, zero-truncated Waring and Poisson-lognormal", {
  set.seed(12)
  n <- 120
  ## DEGPD and ZIDEGPD
  y <- rdiscegpd(n, sigma = 3, xi = 0.2, kappa = 1.5)
  f <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "degpd")
  p <- predict(f, type = "response")
  Fd <- function(t) dgpd_cdf1(t + 1, p$scale, p$shape, p$kappa)
  expect_equal(rqresid(f, seed = 4), rq_ref(Fd(y), ifelse(y == 0, 0, Fd(y - 1)), 4, n),
               tolerance = 1e-8)

  yz <- rzidiscegpd(n, pi = 0.3, sigma = 3, xi = 0.2, kappa = 1.5)
  fz <- egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = yz), family = "zidegpd",
             inits = c(1, log(0.2), log(1.5), qlogis(0.3)))
  pz <- predict(fz, type = "response")
  Fz <- function(t) pz$pi + (1 - pz$pi) * dgpd_cdf1(t + 1, pz$scale, pz$shape, pz$kappa)
  expect_equal(rqresid(fz, seed = 4), rq_ref(Fz(yz), ifelse(yz == 0, 0, Fz(yz - 1)), 4, n),
               tolerance = 1e-8)

  ## generalized Waring and its zero-truncated version
  yw <- rgwaring(n, mu = 4, k = 2, rho = 4)
  ## the estimate has a = k, where the outer-product Hessian is singular (round 6)
  expect_warning(fw <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = yw), family = "gwaring"),
                 class = "egpd_singular_hessian")
  pw <- predict(fw, type = "response")[1, ]
  cw <- cumsum(exp(gw_lpmf(0:20000, pw$mu, pw$k, pw$rho)))
  expect_equal(rqresid(fw, seed = 5),
               rq_ref(cw[yw + 1], ifelse(yw == 0, 0, cw[pmax(yw, 1)]), 5, n), tolerance = 1e-7)

  yt <- yw[yw > 0]
  ft <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = yt), family = "ztgwaring",
             inits = c(log(4), 0, log(2.5)))
  pt <- predict(ft, type = "response")[1, ]
  ct <- cumsum(exp(gw_lpmf(0:20000, pt$mu, pt$k, pt$rho)))
  Ft <- function(t) (ct[t + 1] - ct[1]) / (1 - ct[1])
  expect_equal(rqresid(ft, seed = 6),
               rq_ref(Ft(yt), ifelse(yt == 1, 0, Ft(pmax(yt - 1, 0))), 6, length(yt)),
               tolerance = 1e-7)

  ## Poisson-lognormal: P(Y = y) = E[dpois(y, exp(Z))], Z ~ N(log mu - s^2 / 2, s^2)
  yp <- rpoislnorm(n, mu = 3, sigma = 0.6)
  fp <- egpd(list(y ~ 1, ~1), data.frame(y = yp), family = "plnorm")
  pp <- predict(fp, type = "response")[1, ]
  pmf <- vapply(0:max(yp), function(k) integrate(function(z)
    dpois(k, exp(z)) * dnorm(z, log(pp$mu) - pp$sigma^2 / 2, pp$sigma),
    -Inf, Inf, rel.tol = 1e-12)$value, 0)
  cp <- cumsum(pmf)
  expect_equal(rqresid(fp, seed = 7),
               rq_ref(cp[yp + 1], ifelse(yp == 0, 0, cp[pmax(yp, 1)]), 7, n), tolerance = 1e-5)
})

## ---- xi alias: EGPD families only; consumers select by name -----------------

test_that("xi alias: only EGPD families get it; parameters are consumed by name", {
  set.seed(5)
  yc <- rcomppareto(300, spec = "weibull", shape = 1.5, scale = 2, alpha = 2, theta = 3)
  fc <- egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = yc), family = "comppareto",
             comppareto.args = list(spec = "weibull"))
  pr <- predict(fc, type = "response")
  ## the Weibull body shape is not the GPD tail index, so no `xi` alias
  expect_named(pr, c("shape", "scale", "alpha", "theta"))
  Fw <- pcomppareto(yc, spec = "weibull", shape = pr$shape, scale = pr$scale,
                    alpha = pr$alpha, theta = pr$theta)
  expect_equal(rqresid(fc), qnorm(Fw))
  q <- predict(fc, type = "quantile", prob = c(0.5, 0.9))[1, ]
  expect_equal(unname(unlist(q)),
               qcomppareto(c(0.5, 0.9), spec = "weibull", shape = pr$shape[1],
                           scale = pr$scale[1], alpha = pr$alpha[1], theta = pr$theta[1]))

  set.seed(5)
  y <- rdiscegpd(200, sigma = 3, xi = 0.2, kappa = 1.5)
  f <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "degpd")
  p <- predict(f, type = "response")
  expect_named(p, c("scale", "shape", "kappa", "xi"))
  expect_identical(p$xi, p$shape)
})
