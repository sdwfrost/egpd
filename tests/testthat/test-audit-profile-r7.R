## Round 7 (profile): profile limits are no longer cut short by inner fits that stall on a
## flat direction or end in another basin. Each fit of the other coefficients starts from
## their estimates; on weakly identified fits that start lies on a plateau or in the wrong
## basin far along the profile path, and the overstated deviance there crossed the cutoff
## early. Every finite limit is now checked by refitting the point beyond it from other
## starts, and the search continues along the profile path where that point is in fact
## below the cutoff (the generalized Waring a = k case of test-audit-hesswarn-r6.R keeps
## its interval [1.063, 8.519]).
##
## References are independent of the package's profile code and optimisers: each
## log-likelihood is written out from the model definition (and agrees with the package's
## at its estimates to 1e-12), each profile point is minimised by brute force -- an outer
## grid over one nuisance coefficient, at each grid value an inner grid refined by
## optimize() over the other, then optimize() around the best outer value; for the
## four-parameter DEGPD model 4, multi-start BFGS from 21 starts plus Nelder-Mead along a
## 40-point path with a backward pass -- and the limit is the first crossing of the cutoff
## along a grid from the estimate, refined by uniroot() (tolerance 1e-10).

## Naveau et al.'s model 3 (G(v) = 1 - D_delta{(1 - v)^delta}, D_delta the Beta(1/delta, 2)
## cdf): f(x) = (1 + 1/delta) (1 - S^delta) h(x), with S and h the GPD survival function
## and density; coefficients (log sigma, xi, log delta)
pr7_nll_m3 <- function(b, y) {
  sig <- exp(b[1]); xi <- b[2]; del <- exp(b[3])
  A <- log1p(xi * y / sig) / xi
  -sum(log1p(1 / del) + log(-expm1(-del * A)) - log(sig) - A - log1p(xi * y / sig))
}

test_that("EGPD model 3 with log delta on its plateau: profile limits are the likelihood's", {
  set.seed(11)
  y <- regpd(400, sigma = 2, xi = 0.15, kappa = 1.5)
  f <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "egpd", egpd.args = list(m = 3))
  expect_equal(unname(f$hessian.rank), c(2, 3))
  ## the reference model is the package's (log delta = -6.76, standard error 861)
  expect_equal(-pr7_nll_m3(unname(f$coefficients), y), as.numeric(logLik(f)), tolerance = 1e-10)
  ## a weakly identified fit that is the maximum: profiling finds nothing better
  expect_no_warning(ci <- confint(f, method = "profile"))
  ## from the estimate, the fits ran onto the plateau delta -> Inf (deviance 69.6 for 0.52)
  ## or stayed on delta -> 0: before, the upper scale limit was 1.957 and the lower shape
  ## limit 0.1499
  expect_equal(ci["scale", "upper"], 2.3848248953, tolerance = 1e-3)
  expect_equal(ci["shape", "lower"], 0.0599200219, tolerance = 1e-3)
  ## unchanged
  expect_equal(ci["scale", "lower"], 0.8682787543, tolerance = 1e-3)
  expect_equal(ci["shape", "upper"], 0.3354277961, tolerance = 1e-3)
  expect_equal(ci["delta", ], c(lower = 0, upper = 6.6701591816), tolerance = 1e-3)
})

test_that("DEGPD model 3 with log delta on its plateau: the upper scale limit is not cut short", {
  ## discrete model 3: P(Y = y) = S_F(y) - S_F(y + 1), S_F = S {1 + (1 - S^delta) / delta}
  set.seed(11)
  y <- rdiscegpd(400, sigma = 3, xi = 0.15, kappa = 1.5, type = 1)
  f <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "degpd", degpd.args = list(m = 3))
  expect_equal(unname(f$hessian.rank), c(2, 3))
  ci <- confint(f, parm = "scale", method = "profile")
  expect_equal(ci["scale", "upper"], 3.6491909063, tolerance = 1e-3)  # before: 2.189
})

test_that("a full-rank fit whose profile dips back below the cutoff keeps its lower limit", {
  ## the profile deviance of log sigma rises to 3.38, falls to 3.09 and crosses the cutoff
  ## 3.84 at sigma = 0.925; the fits from the estimate gave 3.84 at sigma = 1.49, where it
  ## is 3.33
  set.seed(3)
  y <- regpd(400, sigma = 2, xi = 0.15, kappa = 1.5)
  f <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "egpd", egpd.args = list(m = 3))
  expect_equal(f$hessian.rank[[1]], 3)
  ci <- confint(f, parm = "scale", method = "profile")
  expect_equal(ci["scale", "lower"], 0.9253801268, tolerance = 1e-3)  # before: 1.491
})

test_that("fits whose kappa underflows are not followed (EGPD model 5)", {
  ## g(v) = sqrt(kappa) phi(sqrt(kappa) (v - 1)) / (Phi(sqrt(kappa)) - 1/2). Fits from the
  ## estimate run log kappa to -742, where the log-likelihood evaluated with a denormal
  ## kappa exceeds the maximum by 3; continuing the profile path from such a point gave a
  ## lower shape limit of -0.319. Before round 7 the limit was -0.1270, cut short by a fit
  ## stalled on the plateau kappa -> 0.
  set.seed(5)
  y <- regpd(150, sigma = 2, xi = 0.15, kappa = 1.5)
  f <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y), family = "egpd", egpd.args = list(m = 5))
  ci <- confint(f, parm = "shape", method = "profile")
  expect_equal(ci["shape", "lower"], -0.1295366275, tolerance = 1e-3)
})

test_that("a DEGPD fit with xi on its boundary plateau: the profile leaves the plateau", {
  ## zero-truncated DEGPD model 1 on negative binomial counts: log xi = -15.6 (standard
  ## error 2544), where the derivative-based fits cannot move it; the profile optima of
  ## scale and kappa away from the estimate have xi > 0.
  ## P(Y = y | Y > 0) = {H(y + 1)^kappa - H(y)^kappa} / {1 - H(1)^kappa}
  set.seed(8)
  y <- rnbinom(300, mu = 5, size = 1.5)
  f <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y[y > 0]), family = "ztdegpd")
  expect_equal(unname(f$hessian.rank), c(2, 3))
  ci <- confint(f, parm = c("scale", "kappa"), method = "profile")
  expect_equal(ci["scale", "lower"], 2.5815057103, tolerance = 1e-3)  # before: 3.262
  expect_equal(ci["kappa", "upper"], 3.4398687402, tolerance = 1e-3)  # before: 2.864
})

test_that("a profile path that switches basins is followed (DEGPD model 4)", {
  ## F(x) = {1 - S3(x)}^(kappa / 2), S3 = S {1 + (1 - S^delta) / delta}. Below scale 2.47
  ## the profile optimum leaves log delta = 2.5 for a branch that runs to delta -> 0, along
  ## which the deviance stays near 1.15 down to scale 2.06. The fits from the estimate, from
  ## the estimate with log xi (standard error 3195) reset, and a short continuation all
  ## stayed on the first branch (deviance 3.84 at scale 1.85, where it is 1.31), giving
  ## 1.850; before round 7 the lower limit was 2.677.
  set.seed(3)
  y <- rdiscegpd(150, sigma = 3, xi = 0.15, kappa = 1.5, type = 1)
  f <- egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = y), family = "degpd", degpd.args = list(m = 4))
  ci <- confint(f, parm = "scale", method = "profile")
  expect_equal(ci["scale", "lower"], 1.2785063264, tolerance = 1e-3)
})

test_that("a model-4 scale profile includes the small-delta branch even with a full-rank fit", {
  set.seed(6)
  y <- regpd(150, sigma = 2, xi = 0.15, kappa = 1.5)
  fit <- egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y), family = "egpd", egpd.args = list(m = 4))
  lower <- unname(confint(fit, parm = "scale", method = "profile")[1, "lower"])

  ## Here the lower crossing is on delta -> 0: for GPD cumulative hazard A,
  ## G = 1 - (1 + A) exp(-A), and the conditional MLE of kappa/2 is -n/sum(log G).
  ## This leaves a one-dimensional likelihood, independent of egpd's optimizer.
  dev <- function(s) {
    ll <- function(xi) {
      A <- log1p(xi * y / s) / xi
      logG <- pgamma(A, shape = 2, log.p = TRUE)
      rate <- -length(y) / sum(logG)
      sum(log(rate) + (rate - 1) * logG + log(A) - log(s) - (1 + xi) * A)
    }
    2 * (as.numeric(logLik(fit)) - optimize(ll, c(0, 1.5), maximum = TRUE)$objective)
  }
  ref <- uniroot(function(s) dev(s) - qchisq(.95, 1), c(.45, .5), tol = 1e-9)$root
  expect_equal(lower, ref, tolerance = 1e-4)  # 0.48243, not the interior-only 0.76324
  expect_equal(dev(lower), qchisq(.95, 1), tolerance = 1e-4)
})

## Composite lognormal-Pareto: phi = 2 theta g(theta) / (alpha G(theta)), with g and G the
## lognormal density and cdf; f(x) = g(x) / {(1 + phi) G(theta)} for 0 < x <= theta, and
## phi / (1 + phi) alpha / (x + theta) {2 theta / (x + theta)}^alpha above theta;
## coefficients (meanlog, log sdlog, log alpha, log theta)
pr7_nll_cp <- function(b, y) {
  mu <- b[1]; s <- exp(b[2]); al <- exp(b[3]); te <- exp(b[4])
  lGt <- plnorm(te, mu, s, log.p = TRUE)
  lphi <- log(2) + b[4] + dlnorm(te, mu, s, log = TRUE) - b[3] - lGt
  l1p <- max(lphi, 0) + log1p(exp(-abs(lphi)))
  lo <- y <= te
  -sum(dlnorm(y[lo], mu, s, log = TRUE) - l1p - lGt) -
    sum(lphi - l1p + b[3] - log(y[!lo] + te) + al * (log(2) + b[4] - log(y[!lo] + te)))
}

test_that("a profile that finds a higher likelihood than the fit's warns that it is no maximum", {
  set.seed(11)
  y <- regpd(400, sigma = 2, xi = 0.15, kappa = 1.5)
  ## restarts = FALSE keeps the stranded single-start fit (the default restarts now refit
  ## it to the maximum), which is the non-optimal fit this test needs
  expect_warning(
    f <- egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y = y), family = "comppareto",
              restarts = FALSE),
    class = "egpd_singular_hessian")
  b <- unname(f$coefficients)
  expect_equal(-pr7_nll_cp(b, y), as.numeric(logLik(f)), tolerance = 1e-10)
  ## independently of the package: Nelder-Mead and BFGS on pr7_nll_cp, from the best point
  ## that profiling theta reached, end at a log-likelihood 8.73 above the fit's
  better <- c(1.457168619, 0.3182270118, 1.207972505, 1.220507835)
  gain <- pr7_nll_cp(b, y) - pr7_nll_cp(better, y)
  expect_gt(gain, 8.7)
  w <- NULL
  withCallingHandlers(confint(f, parm = "meanlog", method = "profile"),
    egpd_profile_better_fit = function(cnd) {
      w <<- cnd
      invokeRestart("muffleWarning")
    })
  expect_s3_class(w, "egpd_profile_better_fit")
  found <- if (is.null(w)) NA else
    as.numeric(sub(".*log-likelihood ([0-9.]+) higher.*", "\\1", conditionMessage(w)))
  expect_true(isTRUE(found > 0.05 && found <= gain))
})

test_that("profiles of ordinary fits do not warn that the fit is not the maximum", {
  set.seed(11)
  yc <- regpd(400, sigma = 2, xi = 0.15, kappa = 1.5)
  yd <- rdiscegpd(400, sigma = 3, xi = 0.2, kappa = 1.5, type = 1)
  ynb <- rnbinom(400, mu = 4, size = 2)
  d <- data.frame(yc = yc, yd = yd, ynb = ynb)
  expect_no_warning(confint(egpd(list(yc ~ 1, ~1, ~1), d, family = "egpd"), method = "profile"))
  expect_no_warning(confint(egpd(list(yd ~ 1, ~1, ~1), d, family = "degpd"), method = "profile"))
  expect_no_warning(confint(egpd(ynb ~ 1, d, family = "bell"), method = "profile"))
  expect_no_warning(confint(egpd(list(ynb ~ 1, ~1), d, family = "pig"), method = "profile"))
})
