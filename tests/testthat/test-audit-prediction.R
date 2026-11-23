audit_prediction_fit <- function(family = "degpd", m = 1L, offset = FALSE,
                                 xi.max = Inf, identity = FALSE, fixed = list()) {
  set.seed(901)
  d <- data.frame(y = floor(2 / .3 * ((1 - runif(160))^(-.3) - 1)),
                  z = seq(-.5, .5, length.out = 160))
  if (family == "egpd") d$y <- d$y + .25
  np <- c(3L, 5L, 3L, 4L, 3L, 3L)[m] + (family == "zidegpd")
  fl <- c(list(if (offset) y ~ offset(z) else y ~ 1), rep(list(~1), np - 1L))
  args <- list(formula = fl, data = d, family = family, fixed = fixed,
               restarts = FALSE)
  args[[paste0(family, ".args")]] <- if (family == "degpd")
    list(m = m, xi.max = xi.max, link = if (identity) "identity" else "log") else
    list(m = m)
  if (family == "zidegpd")
    args$inits <- c(log(2), log(.3), rep(0, np - 3L), qlogis(.2))
  ## A fixed negative shape bounds the support at sigma / |xi|; start sigma beyond
  ## the largest count, where the likelihood is finite (egpd() rejects starts at
  ## which it is not).
  if (!is.null(fixed$xi) && fixed$xi < 0)
    args$inits <- c(log(max(d$y) + 1), rep(0, np - 1L))
  suppressMessages(suppressWarnings(do.call(egpd, args)))
}

audit_prediction_set <- function(fit, values) {
  for (i in seq_len(fit$ngam)) {
    nm <- names(fit)[i]
    fit[[i]]$coefficients[] <- values[[nm]]
    fit$coefficients[fit$idpars == i] <- values[[nm]]
  }
  fit
}

test_that("parameter draws include varying offsets and the bounded inverse link", {
  fit <- audit_prediction_fit(offset = TRUE, xi.max = .5)
  fit <- audit_prediction_set(fit, c(logscale = log(2), logshape = log(.3),
                                    logkappa = log(1.2)))
  fit$Vp <- fit$Vc <- diag(.0004, 3)
  nd <- data.frame(z = c(-1, .7))
  eta <- simulate(fit, nsim = 10000, seed = 12, newdata = nd, type = "link")
  expect_equal(unname(rowMeans(eta$scale)), log(2) + nd$z, tolerance = .001)
  draws <- simulate(fit, nsim = 10000, seed = 12, newdata = nd, type = "response")
  expect_equal(draws$scale, exp(eta$scale))
  expect_equal(draws$shape, .5 * (1 - exp(-exp(eta$shape) / .5)))
  expect_true(all(draws$shape > 0 & draws$shape <= .5))
  inside <- simulate(fit, nsim = 10000, seed = 12, type = "link")
  expect_equal(unname(rowMeans(inside$scale)), log(2) + fit$data$z, tolerance = .001)
})

test_that("bounded shape confidence intervals and delta-method SE use the fitted link", {
  fit <- audit_prediction_fit(xi.max = .5)
  fit <- audit_prediction_set(fit, c(logscale = log(2), logshape = log(.3),
                                    logkappa = log(1.2)))
  fit$Vp <- fit$Vc <- diag(c(.01, .04, .01))
  fit[[2]]$Vp <- fit[[2]]$Vc <- matrix(.04, 1, 1)
  point <- .5 * (1 - exp(-.3 / .5))
  expect_equal(unname(confint(fit, parm = "shape", level = 0)),
               matrix(point, 1, 2))
  expected <- .5 * (1 - exp(-exp(log(.3) + c(-1, 1) * qnorm(.975) * .2) / .5))
  expect_equal(unname(confint(fit, parm = "xi")), matrix(expected, 1, 2))
  pred <- predict(fit, type = "response", se.fit = TRUE)
  expect_equal(pred$fitted$shape, rep(point, nrow(fit$data)))
  expect_equal(pred$se.fit$shape, rep(.3 * exp(-.3 / .5) * .2, nrow(fit$data)))
  local_mocked_bindings(.profile_one_param = function(...) c(log(.2), log(.8)),
                        .package = "egpd")
  expect_equal(unname(confint(fit, parm = "shape", method = "profile")),
               matrix(.5 * (1 - exp(-c(.2, .8) / .5)), 1, 2))
})

test_that("constant offset parameters have natural-scale point intervals and draws", {
  # Test zero-valued fixed predictions independently of native xi=0 derivatives.
  fit <- audit_prediction_fit(identity = TRUE, fixed = list(xi = -.5, kappa = 2))
  fit$likdata$offsets[[2]][] <- 0
  fit[[2]]$offset[] <- 0
  fit$fixed.offsets$.fixedpar_xi <- 0
  expect_true(all(c("shape", "kappa") %in% rownames(confint(fit))))
  for (method in c("wald", "profile")) {
    ci <- confint(fit, parm = c("shape", "kappa"), method = method)
    expect_equal(unname(ci["shape", ]), c(0, 0))
    expect_equal(unname(ci["kappa", ]), c(2, 2))
  }
  expect_equal(unname(confint(fit, parm = "xi")), matrix(0, 1, 2))
  expect_equal(predict(fit, type = "response", se.fit = TRUE)$se.fit$shape,
               rep(0, nrow(fit$data)))
  draws <- simulate(fit, nsim = 4, newdata = data.frame(z = 1:2), type = "response")
  expect_equal(unname(draws$shape), matrix(0, 2, 4))
  expect_equal(unname(draws$kappa), matrix(2, 2, 4))
})

test_that("logit parameter draws remain stable at extreme finite predictors", {
  fit <- audit_prediction_fit(family = "zidegpd")
  fit <- audit_prediction_set(fit, c(logscale = log(2), logshape = log(.3),
                                    logkappa = 0, logitpi = 1000))
  fit$Vp <- fit$Vc <- diag(.0004, 4)
  draws <- simulate(fit, nsim = 5, seed = 12, type = "response")
  expect_true(all(is.finite(draws$pi) & draws$pi == 1))
})

test_that("bounded fixed shape uses the inverse bounded link", {
  fit <- audit_prediction_fit(xi.max = .5, fixed = list(xi = .3, kappa = 1))
  expect_equal(predict(fit, type = "response")$shape, rep(.3, nrow(fit$data)))
  expect_equal(unname(confint(fit, parm = "shape")), matrix(.3, 1, 2))
  draws <- simulate(fit, nsim = 4, seed = 12, newdata = data.frame(z = c(0, .1)),
                     type = "response")
  expect_equal(unname(draws$shape), matrix(.3, 2, 4))
})

test_that("shape cap validation rejects invalid scalars and identity-link caps", {
  d <- data.frame(y = rep(0:4, 10))
  for (cap in list(0, -1, NA_real_, NaN, c(.3, .5), "0.5")) {
    expect_error(egpd(list(y ~ 1, ~1, ~1), d, family = "degpd",
                      degpd.args = list(xi.max = cap)), "positive scalar")
  }
  expect_error(egpd(list(y ~ 1, ~1, ~1), d, family = "degpd",
                    degpd.args = list(xi.max = .5, link = "identity")),
                "finite.*xi.max.*identity")
})

test_that("coverage preserves nonpositive shape and respects a bounded endpoint", {
  fit <- audit_prediction_fit(identity = TRUE, fixed = list(xi = -.5, kappa = 1))
  fit[[1]]$coefficients[] <- log(4)
  fit$coefficients[fit$idpars == 1] <- log(4)
  set.seed(4)
  ans <- predictive_coverage(fit, y = 0, newdata = data.frame(z = 0),
                             nsim = 5000, method = "plug-in")
  expect_lte(ans$U, 7)
  expect_gte(ans$U, 5)
  for (shape in c(-.5, 0)) {
    fit$likdata$offsets[[2]][] <- shape
    fit[[2]]$offset[] <- shape
    fit$fixed.offsets$.fixedpar_xi <- shape
    local_mocked_bindings(rdiscegpd = function(n, xi, ...) {
      expect_true(all(xi == shape))
      rep(0, n)
    }, .package = "egpd")
    expect_equal(predictive_coverage(fit, y = 0, newdata = data.frame(z = 0),
                                     nsim = 20, method = "plug-in")$U, 0)
  }
})

test_that("non-finite predictive parameter and response draws are rejected", {
  for (nm in c("sigma", "xi", "kappa")) {
    pars <- list(sigma = 2, xi = -.5, kappa = 1)
    pars[[nm]] <- Inf
    expect_error(egpd:::.sim_responses(pars, rdiscegpd, 1, 10),
                  paste0("Non-finite parameter draws for '", nm, "'"))
  }
  expect_error(egpd:::.pi_from_sims(matrix(c(0, Inf), 2), 0, .95, FALSE, TRUE),
                "Non-finite simulated responses")
})

test_that("coverage models five and six use the correct distribution type and carrier", {
  for (family in c("degpd", "zidegpd")) for (m in 5:6) {
    fit <- audit_prediction_fit(family = family, m = m)
    fit <- audit_prediction_set(fit, c(logscale = log(2), logshape = log(.3),
                                      logkappa = log(1.7), logitpi = qlogis(.2)))
    set.seed(77)
    sim <- if (family == "degpd")
      rdiscegpd(2000, sigma = 2, xi = .3, kappa = 1.7, type = m - 3L) else
      rzidiscegpd(2000, sigma = 2, xi = .3, kappa = 1.7, pi = .2, type = m - 3L)
    expected <- quantile(sim, c(.025, .975), type = 1, names = FALSE)
    set.seed(77)
    ans <- predictive_coverage(fit, y = 0, newdata = data.frame(z = 0),
                               nsim = 2000, method = "plug-in")
    expect_equal(c(ans$L, ans$U), expected)
  }
})

test_that("integer coverage is exactly membership in the reported interval", {
  ans <- egpd:::.pi_from_sims(matrix(0, 100, 40), rep(0, 40), .95, FALSE, TRUE)
  expect_equal(ans$L, rep(0, 40))
  expect_equal(ans$U, rep(0, 40))
  expect_true(all(ans$covered))
  ys <- cbind(rep(0:9, each = 10), rep(10:19, each = 10))
  for (shortest in c(FALSE, TRUE)) {
    ans <- egpd:::.pi_from_sims(ys, c(0, 20), .75, shortest, TRUE)
    expect_equal(ans$covered, c(0, 20) >= ans$L & c(0, 20) <= ans$U)
  }
})

test_that("shortest intervals minimize contiguous integer width, not an HDR envelope", {
  tie <- egpd:::.pi_from_sims(matrix(c(5, 0), 2), 1, .5, TRUE, TRUE)
  expect_equal(c(tie$L, tie$U), c(0, 0))
  ys <- matrix(c(rep(0, 25), rep(10, 30), rep(11, 25), rep(12, 20)), ncol = 1)
  ans <- egpd:::.pi_from_sims(ys, 11, .75, TRUE, TRUE)
  expect_equal(c(ans$L, ans$U), c(10, 12))
  expect_true(ans$covered)
})

test_that("shortest intervals count observations and preserve exact numeric endpoints", {
  for (x in list(as.double(1:2000), 1e15 + (0:255)/8, 1:100)) {
    level <- switch(as.character(length(x)), `2000` = .95, `256` = .75, `100` = .14)
    h <- round(level * length(x))
    r <- egpd:::.pi_from_sims(as.matrix(rev(x)), x[h+1], level, TRUE, FALSE)
    expect_equal(c(r$L, r$U), as.numeric(x[c(1, h)]), tolerance = 0)
    expect_false(r$covered)
    expect_equal(sum(x >= r$L & x <= r$U), h)
  }
})

test_that("model four residuals use named delta and kappa across all EGPD families", {
  for (family in c("egpd", "degpd", "zidegpd")) {
    fit <- audit_prediction_fit(family = family, m = 4L)
    fit <- audit_prediction_set(fit, c(logscale = log(2), shape = .3,
      logshape = log(.3), logdelta = log(.5), logkappa = log(2), logitpi = qlogis(.2)))
    y <- fit$data$y
    cdf <- switch(family, egpd = pegpd, degpd = pdiscegpd, zidegpd = pzidiscegpd)
    args <- list(sigma = 2, xi = .3, delta = .5, kappa = 2, type = 5)
    if (family == "zidegpd") args$pi <- .2
    upper <- do.call(cdf, c(list(q = y), args))
    if (family == "egpd") {
      expected <- qnorm(upper)
    } else {
      lower <- do.call(cdf, c(list(q = y - 1), args))
      lower[y == 0] <- 0
      set.seed(9)
      expected <- qnorm(runif(length(y), lower, upper))
    }
    expected[is.infinite(expected)] <- NA
    expect_equal(rqresid(fit, seed = 9), expected)
  }
})
