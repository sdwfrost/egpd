## Round 5: the restarts of egpd() and the starts of fit_hurdle() keep the fit that is
## best by the criterion the fits optimise -- the REML criterion when smoothing
## parameters are estimated, the penalised log-likelihood when they are supplied, the
## log-likelihood without smooth terms -- not the highest log-likelihood.
##
## The reference is independent of the selection code: each documented start is
## refitted on its own (restarts = FALSE) and the documented rule (see ?egpd,
## `restarts`) is applied to the refits' criteria.

## egpd()'s documented starts for DEGPD and the continuous EGPD: the package default,
## then shapes 0.05 and 0.3 (log xi for DEGPD, xi for the continuous EGPD) with the other
## parameters filled with 0 and log 2, and the scale at log(mean(y) + 1) (DEGPD) or
## log(mean(y)) (continuous)
r5_starts <- function(y, npar, continuous = FALSE) {
  scale0 <- if (continuous) log(mean(y)) else log(mean(y) + 1)
  c(list(NULL), unlist(lapply(c(0.05, 0.3), function(v) lapply(c(0, log(2)), function(cf)
    c(scale0, if (continuous) v else log(v), rep(cf, npar - 2L)))), recursive = FALSE))
}

## the documented choice among the converged fits, given their criterion values `crit`
## (larger is better) and log-likelihoods `ll`: fits within 0.01 of the best count as
## the same optimum; of those, the default start's (first) fit when `prefer_first`
## (continuous EGPD), and otherwise the highest log-likelihood
r5_choice <- function(crit, ll, prefer_first = FALSE) {
  near <- which(max(crit) - crit <= 0.01)
  if (prefer_first && 1L %in% near) return(1L)
  near[which.max(ll[near])]
}

## the penalised log-likelihood: the log-likelihood less 0.5 sum_j sp_j b' S_j b, with
## S_j the penalty matrices of mgcv's own smooth set-up, on each smooth's coefficients
r5_penll <- function(fit, sp) {
  b <- fit$coefficients
  pen <- 0
  j <- 0
  for (i in seq_len(fit$ngam)) {
    for (sm in fit[[i]]$smooth) {
      j <- j + 1
      idx <- names(b)[fit$idpars == i][sm$first.para:sm$last.para]
      pen <- pen + sp[j] * drop(crossprod(b[idx], sm$S[[1]] %*% b[idx]))
    }
  }
  fit$logLik - pen / 2
}

r5_cont <- function(seed, n = 150) {
  set.seed(seed)
  x <- runif(n)
  u <- runif(n)^(1 / 1.5)
  data.frame(y = exp(0.5 + sin(2 * pi * x)) / 0.2 * ((1 - u)^(-0.2) - 1), x = x)
}

test_that("DEGPD restarts keep the lowest REML criterion, not the highest log-likelihood", {
  set.seed(6)
  n <- 150
  x <- runif(n)
  d <- data.frame(y = rdiscegpd(n, sigma = exp(0.5 + 1.2 * sin(2 * pi * x)), xi = 0.15,
                                kappa = 2, type = 1), x = x)
  fl <- list(y ~ s(x, k = 8), ~1, ~1, ~1)
  fits <- lapply(r5_starts(d$y, 4), function(s)
    egpd(fl, d, family = "degpd", degpd.args = list(m = 4), inits = s, restarts = FALSE))
  expect_true(all(vapply(fits, function(f) f$convergence == 0, NA)))
  ll <- vapply(fits, function(f) f$logLik, 0)
  nr <- vapply(fits, function(f) f$negREML, 0)
  ## the start with the highest log-likelihood (headtree's choice) is a wigglier fit
  ## whose REML criterion is 0.18 worse than the best
  expect_gt(nr[which.max(ll)] - min(nr), 0.1)
  b <- r5_choice(-nr, ll)
  f <- egpd(fl, d, family = "degpd", degpd.args = list(m = 4))
  expect_identical(f$coefficients, fits[[b]]$coefficients)
  expect_identical(f$sp, fits[[b]]$sp)
  expect_lt(f$negREML - min(nr), 0.01)
  expect_gt(f$sp, fits[[which.max(ll)]]$sp)
})

test_that("continuous EGPD restarts: the default fit stands unless another start's REML is better by > 0.01", {
  fl <- list(y ~ s(x, k = 8), ~1, ~1, ~1)
  ## seed 7: another start raises the log-likelihood by 0.27 (> 0.01, so headtree took
  ## it), but its REML criterion is 1.46 worse; the default start's REML is within 0.01
  ## of the best, so the default fit is kept
  d <- r5_cont(7)
  fits <- lapply(r5_starts(d$y, 4, continuous = TRUE), function(s)
    egpd(fl, d, family = "egpd", egpd.args = list(m = 4), inits = s, restarts = FALSE))
  ll <- vapply(fits, function(f) f$logLik, 0)
  nr <- vapply(fits, function(f) f$negREML, 0)
  expect_gt(max(ll) - ll[1], 0.01)
  expect_gt(nr[which.max(ll)] - nr[1], 1)
  expect_identical(r5_choice(-nr, ll, prefer_first = TRUE), 1L)
  f <- egpd(fl, d, family = "egpd", egpd.args = list(m = 4))
  expect_identical(f$coefficients, fits[[1]]$coefficients)

  ## seed 6: two other starts lower the REML criterion by 0.37 (smoother fits, larger
  ## sp, lower log-likelihood); headtree kept the default fit, as no start raised the
  ## log-likelihood by more than 0.01
  d <- r5_cont(6)
  fits <- lapply(r5_starts(d$y, 4, continuous = TRUE), function(s)
    egpd(fl, d, family = "egpd", egpd.args = list(m = 4), inits = s, restarts = FALSE))
  ll <- vapply(fits, function(f) f$logLik, 0)
  nr <- vapply(fits, function(f) f$negREML, 0)
  expect_lt(max(ll) - ll[1], 0.01)
  expect_gt(nr[1] - min(nr), 0.3)
  b <- r5_choice(-nr, ll, prefer_first = TRUE)
  expect_false(b == 1L)
  f <- egpd(fl, d, family = "egpd", egpd.args = list(m = 4))
  expect_identical(f$coefficients, fits[[b]]$coefficients)
  expect_lt(f$negREML - min(nr), 0.01)
})

test_that("with sp supplied, restarts keep the highest penalised log-likelihood", {
  set.seed(43)
  n <- 150
  x <- runif(n)
  d <- data.frame(y = rdiscegpd(n, sigma = exp(0.5 + 1.2 * sin(2 * pi * x)), xi = 0.15,
                                kappa = 2, type = 1), x = x)
  fl <- list(y ~ s(x, k = 8), ~1, ~1, ~1)
  fits <- lapply(r5_starts(d$y, 4), function(s)
    egpd(fl, d, family = "degpd", degpd.args = list(m = 4), sp = 0.02, inits = s,
         restarts = FALSE))
  expect_true(all(vapply(fits, function(f) f$convergence == 0, NA)))
  ll <- vapply(fits, function(f) f$logLik, 0)
  pl <- vapply(fits, r5_penll, 0, sp = 0.02)
  ## the default start reaches an optimum whose log-likelihood is 0.0099 lower than the
  ## other starts' (headtree's choice) but whose penalty is 0.032 smaller
  expect_gt(max(pl) - pl[which.max(ll)], 0.015)
  b <- r5_choice(pl, ll)
  f <- egpd(fl, d, family = "degpd", degpd.args = list(m = 4), sp = 0.02)
  expect_identical(f$coefficients, fits[[b]]$coefficients)
  expect_equal(r5_penll(f, 0.02), max(pl))
})

test_that("fits without smooths are still chosen by the log-likelihood", {
  set.seed(1)
  n <- 200
  x <- runif(n)
  d <- data.frame(y = rdiscegpd(n, sigma = exp(1 + sin(2 * pi * x)), xi = 0.2, kappa = 1.5,
                                type = 1), x = x)
  fl <- list(y ~ x, ~1, ~1)
  fits <- lapply(r5_starts(d$y, 3), function(s)
    egpd(fl, d, family = "degpd", inits = s, restarts = FALSE))
  ll <- vapply(fits, function(f) f$logLik, 0)
  f <- egpd(fl, d, family = "degpd")
  expect_identical(f$coefficients, fits[[which.max(ll)]]$coefficients)

  d <- r5_cont(7)
  fl <- list(y ~ x, ~1, ~1, ~1)
  fits <- lapply(r5_starts(d$y, 4, continuous = TRUE), function(s)
    egpd(fl, d, family = "egpd", egpd.args = list(m = 4), inits = s, restarts = FALSE))
  ll <- vapply(fits, function(f) f$logLik, 0)
  f <- egpd(fl, d, family = "egpd", egpd.args = list(m = 4))
  expect_identical(f$coefficients, fits[[r5_choice(ll, ll, prefer_first = TRUE)]]$coefficients)
})

test_that("REML is not compared across fits that ran to a boundary of a link", {
  ## Every start ends with xi -> 0 (log xi < -9), where the likelihood is flat in log xi
  ## and the REML criterion falls by half a unit per unit of log xi the fit drifted.
  ## Its ranking can coincide with the likelihood ranking, so check the criterion itself.
  set.seed(13)
  n <- 150
  x <- runif(n)
  d <- data.frame(y = rdiscegpd(n, sigma = exp(0.5 + 1.2 * sin(2 * pi * x)), xi = 0.15,
                                kappa = 2, type = 1), x = x)
  fl <- list(y ~ s(x, k = 8), ~1, ~1)
  ## these degenerate fits warn that the REML curvature cannot be resolved
  fits <- lapply(r5_starts(d$y, 3), function(s) suppressWarnings(
    egpd(fl, d, family = "degpd", inits = s, restarts = FALSE)))
  ll <- vapply(fits, function(f) f$logLik, 0)
  expect_true(all(vapply(fits, function(f) all(f$logshape$fitted < -9), NA)))
  expect_identical(egpd:::.start_criterion(fits, ll), ll)
  f <- suppressWarnings(egpd(fl, d, family = "degpd"))
  expect_identical(f$coefficients, fits[[which.max(ll)]]$coefficients)
})

test_that("sandwich-adjusted fits are compared by the log-likelihood", {
  ## The sandwich adjustment is estimated from each start's own unpenalised fit, so its
  ## magnitude k (logLik.adjusted / logLik) differs between starts, and the adjusted
  ## penalised log-likelihood and REML criterion are on different scales: they preferred
  ## the start with the smallest k, whatever its fit. With sp = 0.05 (seed 5) the default
  ## start has k = 0.47 and a log-likelihood 0.61 lower than starts with k = 1.39, and
  ## standard errors 3 times larger; with sp estimated (seed 56) the default start has
  ## k = 0.98, a log-likelihood 0.39 lower and a standard error of 2.1 instead of 0.25,
  ## and a REML criterion 12 units lower, 10 of them from k.
  fl <- list(y ~ s(x, k = 8), ~1, ~1, ~1)
  ## Some starts also warn about the rank of the cluster scores; here compare
  ## the scalar-adjusted fits, letting any unrelated warning through.
  fit <- function(...) withCallingHandlers(egpd(...),
    egpd_cluster_rank = function(w) invokeRestart("muffleWarning"))
  for (case in list(list(seed = 5, sp = 0.05), list(seed = 56, sp = NULL))) {
    set.seed(case$seed)
    n <- 200
    x <- runif(n)
    d <- data.frame(y = rdiscegpd(n, sigma = exp(0.5 + 1.2 * sin(2 * pi * x)), xi = 0.15,
                                  kappa = 2, type = 1), x = x, g = rep(1:20, each = 10))
    ## a start whose unpenalised Hessian is not positive definite cannot be adjusted
    fits <- lapply(r5_starts(d$y, 4), function(s) tryCatch(
      fit(fl, d, family = "degpd", degpd.args = list(m = 4), sp = case$sp,
           sandwich.args = list(id = "g"), inits = s, restarts = FALSE),
      error = function(e) NULL))
    ok <- !vapply(fits, is.null, NA)
    ll <- vapply(fits[ok], function(f) f$logLik, 0)
    k <- vapply(fits[ok], function(f) f$logLik.adjusted / f$logLik, 0)
    expect_gt(diff(range(k)), 0.02)
    f <- fit(fl, d, family = "degpd", degpd.args = list(m = 4), sp = case$sp,
              sandwich.args = list(id = "g"))
    expect_identical(f$coefficients, fits[ok][[which.max(ll)]]$coefficients)
    expect_equal(f$logLik, max(ll))
  }
})

test_that("fit_hurdle() keeps the start with the lowest REML criterion", {
  skip_on_cran()
  set.seed(15)
  n <- 250
  x <- runif(n)
  d <- data.frame(y = rgwaring(n, mu = exp(1 + sin(2 * pi * x)), k = 2, rho = 4), x = x)
  dpos <- d[d$y > 0, ]
  starts <- egpd:::.hurdle_starts(dpos$y, "gwaring", list())
  ## the REML curvature of these fits is not resolved in two directions (a warning)
  fits <- lapply(starts, function(s) suppressWarnings(
    egpd(y ~ s(x, k = 8), dpos, family = "ztgwaring", inits = s)))
  ll <- vapply(fits, function(f) f$logLik, 0)
  nr <- vapply(fits, function(f) f$negREML, 0)
  ## the highest log-likelihood (headtree's choice) has a REML criterion 0.07 worse
  expect_gt(nr[which.max(ll)] - min(nr), 0.05)
  b <- r5_choice(-nr, ll)
  h <- suppressWarnings(fit_hurdle(y ~ s(x, k = 8), d, family = "gwaring"))
  expect_identical(h$positive$coefficients, fits[[b]]$coefficients)
  expect_equal(h$logLik, as.numeric(logLik(h$zero)) + ll[b])
})
