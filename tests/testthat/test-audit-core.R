test_that("pivoted inverses retain parameter order and rank-deficient draws retain shape", {
  H <- matrix(c(1, .9, .1, .9, 1, .1, .1, .1, 1), 3)
  expect_false(identical(attr(chol(H, pivot=TRUE), "pivot"), 1:3))
  expect_equal(egpd:::.solve_egpd(H), solve(H), tolerance=1e-12,
               ignore_attr=TRUE)
  expect_equal(egpd:::.Hdata(H)$iH, solve(H), tolerance=1e-12,
               ignore_attr=TRUE)

  set.seed(611)
  rank1 <- egpd:::.pivchol_rmvn(8, c(2, 3, 4), diag(c(0, 1, 0)))
  expect_equal(dim(rank1), c(3L, 8L))
  expect_equal(unname(rank1[c(1, 3), ]), rbind(rep(2, 8), rep(4, 8)))
  fixed <- egpd:::.pivchol_rmvn(8, c(2, 3), matrix(0, 2, 2))
  expect_equal(unname(fixed), rbind(rep(2, 8), rep(3, 8)))
})

test_that("optimizers move with zero step coordinates and report exhausted budgets", {
  ctl <- list(steptol=1e-12, itlim=100, fntol=1e-10,
              gradtol=1e-8, stepmax=3)
  f <- function(p) sum((p - c(1, 0))^2)
  g <- function(p) 2 * (p - c(1, 0))
  fit <- egpd:::.BFGS(c(0, 0), f, g, control=ctl)
  expect_equal(fit$par, c(1, 0), tolerance=1e-8)
  expect_equal(fit$convergence, 0L)

  ctl$itlim <- 1
  ctl$stepmax <- .1
  f <- function(p) sum(p^2)
  g <- function(p) 2 * p
  search <- function(p, kept) {
    step <- p
    attr(step, "gradient") <- 2 * p
    attr(step, "Hessian") <- diag(2, length(p))
    step
  }
  for (fit in list(egpd:::.BFGS(c(10, 10), f, g, control=ctl),
                   egpd:::.newton_step_inner(c(10, 10), f, search, control=ctl))) {
    par <- if (!is.null(fit$pars)) fit$pars else fit$par
    expect_equal(par, c(9.9, 9.9))
    expect_equal(as.vector(fit$gradient), 2 * par)
    expect_equal(fit$convergence, 1L)
    expect_match(fit$report, "iteration limit")
    expect_lte(fit$iterations, ctl$itlim)
  }
})

test_that("fixed response parameters use their actual link", {
  form <- list(y ~ 1, ~ 1, ~ 1, ~ 1)
  out <- egpd:::.setup.fixed.egpd(list(xi=.3, pi=.2), form,
                                  data.frame(y=1:4),
                                  c("lsigma", "lxi", "lkappa", "logitpi"),
                                  xi.max=.5)
  eta <- out$offsets[[".fixedpar_lxi"]]
  expect_equal(.5 * -expm1(-exp(eta) / .5), .3)
  expect_equal(plogis(out$offsets[[".fixedpar_logitpi"]]), .2)
  expect_error(egpd:::.setup.fixed.egpd(list(xi=.5), form, data.frame(y=1:4),
                                        c("lsigma", "lxi", "lkappa", "logitpi"),
                                        xi.max=.5), "below xi.max")
  scale <- egpd:::.setup.fixed.egpd(list(sigma=2), form, data.frame(y=1:4),
                                    c("lsigma", "lxi", "lkappa", "logitpi"))
  expect_equal(all.vars(scale$formula[[1]])[1], "y")
})

test_that("compact fits preserve linear covariates and offsets", {
  set.seed(612)
  d <- data.frame(t=rep(1:20, each=4), x=rep(c(0, 1), 40),
                  z=rep(c(0, .2, .4, .6), 20))
  d$y <- rpois(nrow(d), exp(.8 + .02 * d$t + .6 * d$x + d$z))
  args <- list(formula=list(y ~ s(t, k=4) + x + offset(z), ~1, ~1),
               data=d, family="degpd", fixed=list(xi=.3, kappa=1),
               restarts=FALSE, sp=.2)
  compact <- do.call(egpd, c(args, list(compact=TRUE)))
  full <- do.call(egpd, args)
  expect_equal(unname(as.matrix(predict(compact))),
               unname(as.matrix(predict(full))), tolerance=1e-6)
  expect_equal(as.numeric(logLik(compact)), as.numeric(logLik(full)), tolerance=1e-6)
  expect_equal(compact$convergence, 0L)

  intercept <- egpd(list(y~1, ~1, ~1), d, family="degpd",
                    fixed=list(xi=.3, kappa=1), compact=TRUE, restarts=FALSE)
  expect_true(is.finite(as.numeric(logLik(intercept))))
})

test_that("all parameter models share complete rows and reject fractional counts", {
  d <- data.frame(y=rep(0:5, 10), x=seq_len(60))
  d$x[5] <- NA
  expect_warning(f <- egpd(list(y~x, ~1, ~1), d, family="degpd",
                            fixed=list(xi=.3, kappa=1), restarts=FALSE),
                  "Omitting 1 rows")
  expect_equal(nobs(f), 59)
  expect_equal(nrow(predict(f)), 59)
  expect_error(egpd(list(y~1, ~1, ~1), transform(d, y=y+.25),
                    family="degpd", restarts=FALSE),
                "non-negative integer")
  expect_error(distfit(c(0, .5, 2), family="degpd"), "non-negative integer")
})

test_that("REML curvature follows refitted coefficients", {
  set.seed(904)
  d <- data.frame(x=seq(0, 1, length.out=240))
  scale <- exp(1 + 1.5 * sin(2 * pi * d$x))
  d$y <- floor(scale / .3 * ((1 - runif(nrow(d)))^(-.3) - 1))
  fit <- egpd(list(y~s(x, k=7), ~1, ~1), d, family="degpd",
               fixed=list(xi=.3, kappa=1), restarts=FALSE,
               control=list(inner=list(gradtol=1e-8),
                            outer=list(gradtol=1e-6)))
  rho <- log(fit$sp)
  objective <- function(r) as.numeric(egpd:::.reml0(
    r, fit$likfns, fit$likdata, fit$Sdata, beta=fit$coefficients))
  h <- .01
  reference <- (objective(rho+h) - 2*objective(rho) + objective(rho-h)) / h^2
  analytic <- egpd:::.reml12(rho, fit$likfns, fit$likdata, fit$Sdata,
                             beta=fit$coefficients)[[2]][1, 1]
  numeric <- egpd:::.reml2.fdfd(rho, fit$likfns, fit$likdata, fit$Sdata,
                                beta=fit$coefficients)[1, 1]
  expect_equal(analytic, reference, tolerance=.002)
  expect_equal(numeric, reference, tolerance=.02)
})

test_that("failed profile refits cannot become flat finite-difference gradients", {
  local_mocked_bindings(
    .reml0=function(...) structure(1e20, inner.convergence=1L),
    .package="egpd")
  expect_error(egpd:::.reml1.fd(0, NULL, NULL, NULL),
                "Inner optimization failed")
})

test_that("Bernstein fitting excludes fixed parameters from optimization and IC", {
  set.seed(613)
  x <- -2 * log(runif(100))
  for (fixed in list(list(kappa=2), list(sigma=2, xi=0))) {
    fit <- distfit(x, method="bernstein", bernstein.m=3,
                    fix.arg=fixed, hessian=FALSE)
    expect_false(any(names(fixed) %in% names(coef(fit))))
    expect_equal(fit$npar, 3 - length(fixed) + 2)
    pars <- as.list(c(coef(fit), unlist(fit$fix.arg)))
    got <- sum(log(egpd:::.bernstein_full_density(
      x, pars$sigma, pars$xi, pars$kappa, fit$bernstein.weights, 3)))
    expect_equal(as.numeric(logLik(fit)), got, tolerance=1e-8)
  }
})

test_that("Bernstein uncertainty requires identifiable positive curvature", {
  set.seed(614)
  x <- -2 * log(runif(100))
  fit <- distfit(x, method="bernstein", bernstein.m=1,
                  fix.arg=list(sigma=2, xi=0), optim.method="BFGS")
  exact <- -length(x) / sum(log(-expm1(-x/2)))
  expect_equal(unname(coef(fit)), exact, tolerance=1e-4)
  expect_equal(unname(fit$sd), exact / sqrt(length(x)), tolerance=1e-4)
  expect_equal(egpd:::.bernstein_full_quantile(1, sigma = 2, xi = -0.25, kappa = 1,
                                               weights = 1, m = 1), 8)
})

test_that("restarts preserve winning-fit warnings and validate degpd$m", {
  set.seed(615)
  d <- data.frame(y = rdiscegpd(60, sigma = 3, xi = 0.2, kappa = 1, type = 1),
                  x = seq_len(60))
  d$x[5] <- NA
  expect_warning(
    f <- egpd(list(y ~ x, ~1, ~1), d, family = "degpd", fixed = list(kappa = 1)),
    "Omitting 1 rows"
  )
  expect_equal(nobs(f), 59)
  expect_error(egpd(y ~ 1, d[1:10, ], family = "degpd", degpd.args = list(m = 99)),
               "degpd\\$m must be 1, 2, 3, 4, 5, or 6")
})

test_that("named formula lists reorder to canonical family slots including zidegpd m=4", {
  set.seed(616)
  n <- 300
  x1 <- seq(-1, 1, length.out = n)
  y4 <- rzidiscegpd(n, pi = 0.2, sigma = 3, xi = 0.25,
                    delta = exp(-0.5 + 0.5 * x1), kappa = 3, type = 5)
  df4 <- data.frame(y = y4, x1 = x1)
  fit <- suppressMessages(egpd(
    list(lsigma = y ~ 1, lxi = ~1, ldelta = ~x1, lkappa = ~1, logitpi = ~1),
    data = df4, family = "zidegpd", zidegpd.args = list(m = 4),
    fixed = list(kappa = 3), restarts = FALSE
  ))
  expect_true("logdelta_x1" %in% names(coef(fit)))
  expect_false("logkappa_x1" %in% names(coef(fit)))
  pars <- predict(fit, type = "response")[1, ]
  expect_equal(pars$kappa, 3)

  ## zidegpd quantiles match qzidiscegpd at prob = 1 and near F(0)
  F0 <- pzidiscegpd(0, pi = pars$pi, sigma = pars$scale, xi = pars$shape,
                    kappa = pars$kappa, delta = pars$delta, type = 5)
  probs <- c(0, F0 + 1e-9, 0.9, 1)
  q_pred <- as.numeric(predict(fit, type = "quantile", prob = probs)[1, ])
  q_ref <- qzidiscegpd(probs, pi = pars$pi, sigma = pars$scale, xi = pars$shape,
                       kappa = pars$kappa, delta = pars$delta, type = 5)
  expect_equal(q_pred, q_ref)
})

test_that("zero-truncated families reject zeros and support quantiles/rqresid", {
  d_bad <- data.frame(y = c(0, 1, 2, 3, 4))
  expect_error(egpd(list(y ~ 1, ~1, ~1), d_bad, family = "ztdegpd", restarts = FALSE),
               "positive integer responses")
  expect_error(egpd(list(y ~ 1, ~1, ~1), d_bad, family = "ztgwaring", restarts = FALSE),
               "positive integer responses")

  set.seed(617)
  y_pos <- pmax(1L, rdiscegpd(80, sigma = 3, xi = 0.2, kappa = 1, type = 1))
  f_zt <- egpd(list(y ~ 1, ~1, ~1), data.frame(y = y_pos), family = "ztdegpd",
               fixed = list(kappa = 1), restarts = FALSE)
  expect_lte(as.numeric(logLik(f_zt)), 0)
  expect_equal(predict(f_zt, type = "quantile", prob = 0)[1, 1], 1)
  expect_true(all(is.finite(rqresid(f_zt, seed = 1))))

  f_prl <- egpd(y ~ 1, data.frame(y = rprl(60, 3)), family = "prl")
  expect_true(all(is.finite(predict(f_prl, type = "quantile", prob = 0.5)[, 1])))
  expect_true(all(is.finite(rqresid(f_prl, seed = 1))))
})

test_that("distfit warns on indefinite MLE Hessian, normalizes cpegpd by h, and plots gwaring", {
  set.seed(618)
  x <- regpd(30, sigma = 2, xi = 0.2, kappa = 1.5, delta = 3, prob = 0.4, type = 6)
  expect_warning(
    fit_bad <- distfit(x, type = 6, family = "egpd",
                       start = list(sigma = 5, xi = -0.05, kappa = 0.2, delta = 0.21, prob = 0.5),
                       control = list(maxit = 2)),
    "not positive definite"
  )
  expect_null(fit_bad$vcov)
  expect_true(all(is.na(fit_bad$sd)))

  x_cp <- rcpegpd(120, lambda = 1.5, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
  f_h1 <- distfit(x_cp, family = "cpegpd", cpegpd.h = 0.05, hessian = FALSE,
                  fix.arg = list(sigma = 2, xi = 0.1, kappa = 1.5))
  f_h2 <- distfit(x_cp, family = "cpegpd", cpegpd.h = 0.025, hessian = FALSE,
                  fix.arg = list(sigma = 2, xi = 0.1, kappa = 1.5))
  expect_lt(abs(AIC(f_h1) - AIC(f_h2)) / length(x_cp), 0.3)

  xw <- rgwaring(80, mu = 2, k = 1.5, rho = 3)
  fw <- distfit(xw, family = "gwaring", hessian = FALSE)
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(plot(fw))
})

test_that("exported utility helpers and marginal vcov/confint work cleanly", {
  A <- matrix(c(2, 1, 1, 2), 2, 2)
  expect_equal(pinv(A), solve(A), tolerance = 1e-12)
  expect_equal(ginv.egpd(A), solve(A), tolerance = 1e-12)
  expect_equal(seq_between(c(2, 5), length = 4), c(2, 3, 4, 5))
  expect_equal(seq_between(c(2, 4)), 2:4)
  d1 <- data.frame(a = 1:2, d = as.Date("2025-01-01") + 0:1,
                   t = as.POSIXct("2025-01-01", tz = "UTC") + 0:1)
  out_df <- dfbind(list(d1, d1))
  expect_equal(nrow(out_df), 4L)
  expect_s3_class(out_df$d, "Date")

  set.seed(619)
  ds <- data.frame(x = seq(0, 1, length.out = 80))
  ds$y <- rdiscegpd(80, sigma = exp(0.8 + 0.5 * sin(2 * pi * ds$x)), xi = 0.2, kappa = 1, type = 1)
  fs <- egpd(list(y ~ s(x, k = 5), ~1, ~1), ds, family = "degpd",
             fixed = list(kappa = 1), restarts = FALSE)
  expect_equal(unname(vcov(fs, marginal = FALSE)), unname(fs$Vp))
  expect_equal(unname(vcov(fs, marginal = TRUE)), unname(fs$Vc))
  expect_equal(dim(confint(fs, parm = "shape", marginal = TRUE)), c(1L, 2L))
})
