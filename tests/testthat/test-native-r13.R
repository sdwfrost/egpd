test_that("a free GPIG nuisance at its floor can leave it on the profile path", {
  set.seed(7)
  y <- rnbinom(300, mu = 4, size = 2)
  y[sample(300, 30)] <- 0L
  d <- data.frame(y, o = 1.5)
  for (zi in c(FALSE, TRUE)) {
    form <- list(y ~ 1 + offset(o), ~1, ~1)
    if (zi) form[[4L]] <- ~1
    fam <- if (zi) "zigpignat" else "gpignat"
    expect_warning(f <- egpd(form, d, family = fam), class = "egpd_singular_hessian")
    expect_equal(predict(f, type = "response")$a[1], .001)
    ci <- confint(f, "b", method = "profile")
    ## Direct LR inversion of the coefficients of the probability generating
    ## function, maximised over a >= .001, c and (for ZI) pi at each b.
    expect_lt(abs(ci[1] - if (zi) 4.460916 else 6.323878), .002)
    expect_lt(abs(ci[2] - if (zi) 2573.161 else 1522.136), .2)
  }
  ## A likelihood-ratio profile works here; a regular sandwich inverse does not.
  d$row <- seq_along(y)
  expect_error(egpd(list(y ~ 1, ~1, ~1), d, family = "gpignat",
                    sandwich.args = list(id = "row", force = TRUE)),
               "observed information.*singular")
  f <- egpd(list(y ~ 1, ~1, ~1), d, family = "gpignat", fixed = list(a = .25, c = .75))
  ci <- confint(f, "b", method = "profile")
  dev <- vapply(ci, function(b) 2 * (f$logLik - sum(dgpig(y, .25, b, .75, log = TRUE))), 0)
  expect_lt(max(abs(dev - qchisq(.95, 1))), .002)

  ## Well-identified interior MLE: even the wider upper branch is retained.
  set.seed(8)
  y <- rgpig(1000, a = .5, b = 30, c = .9)
  interior <- egpd(list(y ~ 1, ~1, ~1), data.frame(y), family = "gpignat")
  ci <- confint(interior, "b", method = "profile")
  expect_equal(as.numeric(ci), c(27.475, 190.578), tolerance = 2e-4)
})

test_that("Waring k profiles compare the interior with the NB limit", {
  for (seed in 102:103) {
    set.seed(seed)
    y <- rnbinom(300, mu = 4, size = if (seed == 102) .7 else 5)
    f <- egpd(list(y ~ 1, ~1, ~1), data.frame(y), family = "gwaring")
    ci <- confint(f, "k", method = "profile")
    if (seed == 102) {
      ## Reference from a bounded fit of dgwaring at each k.
      expect_equal(unname(ci[2]), .93032, tolerance = 2e-4)
    } else {
      expect_identical(unname(ci[2]), Inf)
    }
  }
})

test_that("exposure and dispersion offsets do not change the PL-normal fit", {
  set.seed(43)
  n <- 400
  x <- runif(n); E <- exp(runif(n, log(1e3), log(1e5)))
  d <- data.frame(x, y = rnbinom(n, mu = exp(1 + x), size = 2), o = log(1e5), s = -8)
  form <- list(y ~ x + offset(o), ~1 + offset(s))
  f <- egpd(form, d, family = "plnorm")
  g <- egpd(form, d, family = "plnorm", inits = c(1 - d$o[1], -.4 - d$s[1]))
  expect_equal(f$convergence, 0L)
  expect_equal(coef(f), coef(g), tolerance = 1e-5)
  ## Direct standard-normal mixture integral; the trapped fit was -1154.08.
  expect_lt(abs(f$logLik + 993.058501), 1e-5)
  expect_equal(predict(f, type = "response")$sigma[1], .64476, tolerance = 1e-5)
  expect_equal(as.numeric(predict(f, data.frame(x = .5, o = log(1e5), s = -8),
                                  prob = c(.05, .95))), c(0, 11))

  ## Same exposure, but heterogeneous across rows and on the *second* predictor.
  vary <- transform(d, y = rnbinom(n, mu = exp(-7 + x) * E, size = 2), o = log(E))
  g <- egpd(list(y ~ 1, ~x + offset(o), ~1), vary, family = "gpignat")
  expect_equal(g$convergence, 0L)
  expect_lt(abs(g$logLik + 1549.318), .002)  # direct bounded optimisation of dgpig

  set.seed(1)
  d$y[sample(n, 40)] <- 0L
  d$o <- log(1e3)
  z <- withCallingHandlers(egpd(list(y ~ x + offset(o), ~1, ~1, ~1), d, family = "zigpig"),
                           egpd_singular_hessian = function(w) invokeRestart("muffleWarning"))
  expect_equal(z$convergence, 0L)
  expect_lt(abs(z$logLik + 972.710), .002)  # direct bounded optimisation of dZIGPIG
})

test_that("zero-truncated and continuous EGPD starts use total predictor levels", {
  set.seed(35)
  x <- runif(240)
  d <- data.frame(x, y = pmax(1, rdiscegpd(240, sigma = exp(1 + x), xi = .2, kappa = 1.5)),
                  o = log(1e5))
  f <- egpd(list(y ~ x + offset(o), ~1, ~1), d, family = "ztdegpd")
  g <- egpd(list(y ~ x, ~1, ~1), d, family = "ztdegpd")
  expect_equal(f$convergence, 0L)
  expect_lt(abs(f$logLik - g$logLik), 1e-5)
  ## Here a *restart* of continuous model 2 beats the primary fit. Check
  ## the offset parametrisation against an independently feasible mixture.
  set.seed(11)
  d <- data.frame(y = regpd(400, sigma = 2, xi = .15, kappa = 1.5), o = log(1e5))
  mix <- egpd(list(y ~ offset(o), ~1, ~1, ~1, ~1), d, family = "egpd", egpd.args = list(m = 2))
  b <- c(.6460399, .00407764, .4830601, 8.269654, 3.846758)
  A <- log1p(b[2] * d$y / exp(b[1])) / b[2]; F <- -expm1(-A)
  k1 <- exp(b[3]); k2 <- k1 + exp(b[4]); p <- plogis(b[5])
  ll <- sum(log(p*k1*F^(k1-1) + (1-p)*k2*F^(k2-1)) - b[1] - (1 + b[2])*A)
  expect_equal(mix$convergence, 0L)
  expect_gte(mix$logLik, ll - 1e-5)
})

test_that("starts on a no-intercept predictor use its design and leave pins alone", {
  x <- rep(c(-1, 0, 1), each = 80)
  set.seed(24)
  d <- data.frame(x, o = 1 + x / 5, y = rpois(length(x), exp(1 + .6*x)))
  f <- egpd(list(y ~ 0 + x + offset(o), ~1), d, family = "plnorm",
            fixed = list(sigma = .3), compact = TRUE)
  nll <- function(b) -sum(dpoislnorm(d$y, exp(b*x + d$o), .3, log = TRUE))
  expect_equal(unname(coef(f)[1]), optimize(nll, c(-1, 2))$minimum, tolerance = 1e-4)
  expect_identical(unname(coef(f)[2]), 0)
  expect_equal(predict(f, type = "response")$sigma, rep(.3, length(x)))
})

## Independent log-PMF scores, with the supplied design matrix for log(mu).
r13_gpois_scores <- function(b, X, y) {
  p <- ncol(X)
  mu <- exp(drop(X %*% b[seq_len(p)]))
  l <- plogis(b[p + 1L]); th <- mu * (1 - l)
  t <- (y - 1) / (th + l*y) - 1
  cbind(X * (1 + th*t), -l + l*(1 - l)*(y - mu)*t)
}

r13_info <- function(score, b) {
  h <- 2e-5
  H <- -vapply(seq_along(b), function(i) {
    e <- replace(numeric(length(b)), i, h)
    colSums(score(b + e) - score(b - e)) / (2*h)
  }, numeric(length(b)))
  (H + t(H)) / 2
}

test_that("both sandwich adjustments use observed information and cluster score sums", {
  set.seed(213)
  n <- 800; x <- runif(n, -1, 1); s <- .95*(x > .1)
  d <- data.frame(y = rpois(n, exp(1 + .7*x + rnorm(n, -s^2/2, s))), x,
                  row = seq_len(n), group = rep(1:80, each = 10))
  X <- model.matrix(~x, d)
  for (fam in c("gpois", "bell")) {
    gpois <- fam == "gpois"
    form <- if (gpois) list(y ~ x, ~1) else list(y ~ x)
    score <- if (gpois) function(b) r13_gpois_scores(b, X, d$y) else function(b) {
      mu <- exp(drop(X %*% b))
      theta <- vapply(mu, function(m) uniroot(function(t) t*exp(t) - m, c(0, 20),
                                              tol = 1e-12)$root, 0)
      X * ((d$y - mu)/(1 + theta))
    }
    for (id in c("row", "group")) {
      fit <- function(method) egpd(form, d, family = fam,
                                   sandwich.args = list(id = id, method = method))
      fc <- fit("curvature"); fm <- fit("magnitude")
      b <- coef(fc)
      H <- r13_info(score, b); iH <- solve(H)
      J <- crossprod(rowsum(score(b), d[[id]]))
      k <- length(b) / sum(diag(iH %*% J))
      V <- iH %*% J %*% iH
      expect_equal(fc$convergence, 0L)
      expect_equal(unname(fc$Vp), unname(V), tolerance = 5e-5)
      expect_equal(unname(vcov(fc)), unname(V), tolerance = 5e-5)
      expect_equal(fm$likdata$k, k, tolerance = 5e-5)
      expect_equal(unname(fm$Vp), unname(iH / k), tolerance = 5e-5)
      expect_equal(fm$Vc, fm$Vp)
    }
  }
  ## Likelihood curvature can be invertible while the sandwich itself is not.
  d$two <- rep(1:2, length.out = n)
  expect_error(egpd(list(y ~ x, ~1), d, family = "gpois",
                    sandwich.args = list(id = "two", method = "curvature", force = TRUE)),
               "cluster score covariance.*singular")
  expect_warning(egpd(list(y ~ x, ~1), d, family = "gpois", sandwich.args = list(id = "two")),
                 "cluster score covariance is rank deficient")
  fc <- egpd(list(y ~ x, ~1), d, family = "gpois", fixed = list(lambda = .35),
             sandwich.args = list(id = "group", method = "curvature"))
  S <- function(b) r13_gpois_scores(c(b, qlogis(.35)), X, d$y)[, 1:2]
  b <- coef(fc)[1:2]
  iH <- solve(r13_info(S, b))
  V <- iH %*% crossprod(rowsum(S(b), d$group)) %*% iH
  expect_equal(unname(fc$Vp[1:2, 1:2]), unname(V), tolerance = 1e-5)
  expect_equal(unname(fc$Vp[3, ]), c(0, 0, 0))
})

test_that("penalised sandwich covariance is curvature at the adjusted fit, plus S", {
  set.seed(732)
  x <- rep(seq(0, 1, length.out = 40), 10)
  d <- data.frame(y = rnbinom(length(x), mu = exp(1 + sin(2*pi*x)), size = 2),
                  x, grp = rep(seq_len(80), each = 5))
  for (method in c("magnitude", "curvature")) {
    f <- egpd(list(y ~ s(x, k = 5), ~1), d, family = "gpois", compact = TRUE,
              sandwich.args = list(id = "grp", method = method))
    expect_equal(f$convergence, 0L)
    ld <- f$likdata; b <- coef(f); X <- f[[1L]]$X
    score <- function(z) r13_gpois_scores(z, X, d$y)
    at <- drop(ld$compmode + ld$CH %*% (b - ld$compmode))
    H0 <- ld$k * crossprod(ld$CH, r13_info(score, at) %*% ld$CH)
    V <- solve(H0 + ld$S)
    expect_equal(unname(f$Vp), unname(V), tolerance = 1e-6)
    db <- V %*% ld$S %*% b  # minus derivative w.r.t. the single log(sp)
    expect_gt(drop(f$Vlsp), 0)
    expect_equal(as.numeric(vcov(f)), as.numeric(V + drop(f$Vlsp) * tcrossprod(db)),
                 tolerance = 1e-6)
    expect_equal(f$df, sum(diag(V %*% H0)), tolerance = 1e-6)
  }
})
