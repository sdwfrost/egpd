d34_data <- function(y) {
  list(X = rep(list(matrix(1, length(y), 1)), 3),
       y = cbind(y), idpars = 1:3, dupid = 0L, duplicate = 0L,
       offsets = rep(list(numeric()), 3), censored = FALSE, xi.max = Inf)
}

test_that("DEGPD d34 agrees with derivatives of d12 in the packed index order", {
  ind <- egpd:::.indices(3)
  for (identity in c(FALSE, TRUE)) {
    fns <- if (identity) egpd:::.degpd1idfns else egpd:::.degpd1fns
    for (xi in if (identity) c(0.3, -0.3, 0.06, -0.06) else c(0.3, 0.01, 0.002)) {
      p <- c(log(2), if (identity) xi else log(xi), log(1.5))
      ld <- d34_data(if (xi < 0) c(0, 1, 5, 6, 7) else c(0, 1, 5))
      got <- fns$d340(p, ld)
      expect_true(all(is.finite(got)))
      for (k in 1:3) {
        h <- 1e-4
        pp <- pm <- p
        pp[k] <- p[k] + h
        pm[k] <- p[k] - h
        d3 <- (fns$d120(pp, ld) - fns$d120(pm, ld)) / (2 * h)
        d4 <- (fns$d340(pp, ld) - fns$d340(pm, ld)) / (2 * h)
        for (i in 1:3) for (j in i:3) {
          expect_equal(got[, ind$i3[i, j, k]], d3[, ind$i2[i, j]],
                       tolerance = 2e-5)
          for (l in j:3)
            expect_equal(got[, ind$i4[i, j, k, l]], d4[, ind$i3[i, j, l]],
                         tolerance = 2e-4)
        }
      }
    }
  }
})

test_that("DEGPD d34 expands duplicate designs and applies all offsets", {
  for (fns in list(egpd:::.degpd1fns, egpd:::.degpd1idfns)) {
    p <- c(log(2), 0.3, log(1.5))
    ld <- d34_data(c(0, 1, 5))
    ld$X <- rep(list(matrix(1, 2, 1)), 3)
    ld$duplicate <- 1L
    ld$dupid <- c(0L, 1L, 0L)
    ld$offsets <- list(c(0.1, 0.2, 0.3), c(0.02, 0.03, 0.04),
                       c(-0.1, 0, 0.1))
    got <- fns$d340(p, ld)
    want <- do.call(rbind, lapply(1:3, function(i) {
      pi <- p + vapply(ld$offsets, `[`, numeric(1), i)
      fns$d340(pi, d34_data(ld$y[i, 1]))
    }))
    expect_equal(got, want)
  }
})

test_that("unsafe analytic d34 values request the finite-difference solver", {
  ld <- d34_data(c(0, 1, 5))
  for (xi in c(0, 1e-6, -1e-6))
    expect_error(egpd:::.degpd1idfns$d340(c(log(2), xi, log(1.5)), ld),
                 class = "egpd_d34_unavailable")
  expect_error(egpd:::.degpd1fns$d340(c(log(2), log(1e-15), log(1.5)), ld),
               class = "egpd_d34_unavailable")
  ld$xi.max <- 0.5
  expect_error(egpd:::.degpd1fns$d340(c(log(2), log(0.3), log(1.5)), ld),
               class = "egpd_d34_unavailable")
})

test_that("smooth fits retain finite-difference support for capped and small shape", {
  set.seed(12)
  d <- data.frame(y = rnbinom(80, size = 2, mu = 5), x = seq_len(80))
  form <- list(y ~ s(x, k = 5), ~ 1, ~ 1)
  capped <- egpd(form, d, family = "degpd", restarts = FALSE,
                 degpd.args = list(m = 1, xi.max = 0.5))
  expect_true(is.finite(as.numeric(logLik(capped))))

  args <- list(formula = form, data = d, family = "degpd", restarts = FALSE,
               degpd.args = list(m = 1, link = "identity"),
               fixed = list(xi = 0.01, kappa = 1.5))
  expect_warning(fallback <- do.call(egpd, args), "finite-difference outer")
  fd <- do.call(egpd, c(args, list(outer = "fd")))
  expect_equal(as.numeric(logLik(fallback)), as.numeric(logLik(fd)),
               tolerance = 1e-7)
})

test_that("analytic Newton smoothing agrees with finite differences", {
  set.seed(31)
  d <- data.frame(x = seq(0, 1, length.out = 240))
  d$y <- rDEGPD1(240, mu = exp(1 + 1.5 * sin(2 * pi * d$x)),
                  sigma = 0.3, nu = 1.5)
  args <- list(formula = list(y ~ s(x, k = 7), ~ 1, ~ 1),
               data = d, family = "degpd", restarts = FALSE,
               fixed = list(xi = 0.3, kappa = 1.5),
               control = list(inner = list(gradtol = 1e-6),
                              outer = list(gradtol = 1e-6)))
  analytic <- do.call(egpd, c(args, list(outer = "newton")))
  fd <- do.call(egpd, c(args, list(outer = "fd")))
  expect_equal(as.numeric(logLik(analytic)), as.numeric(logLik(fd)),
               tolerance = 1e-5)
  expect_equal(unname(analytic$sp), unname(fd$sp), tolerance = 1e-3)
})

test_that("optimizer control overrides reach the inner and outer iterations", {
  set.seed(11)
  d <- data.frame(y = rnbinom(80, size = 2, mu = 5), x = seq_len(80))
  fit <- egpd(list(y ~ s(x, k = 5), ~ 1, ~ 1), d, family = "degpd",
              degpd.args = list(xi.max = 0.5), restarts = FALSE,
              control = list(inner = list(gradtol = 1e-5),
                             outer = list(gradtol = 1e-3)))
  expect_equal(fit$likdata$control$inner$gradtol, 1e-5)
  expect_equal(fit$likdata$control$outer$gradtol, 1e-3)
  expect_equal(fit$likdata$control$inner$itlim, 100)
  expect_equal(fit$likdata$control$outer$stepmax, 3)

  for (control in list(list(list(gradtol = 1e-5)),
                       list(outer = list(unknown = 1)),
                       list(inner = list(gradtol = NULL)),
                       list(inner = list(gradtol = NA_real_)),
                       list(inner = list(1e-5)))) {
    expect_error(egpd(list(y ~ 1, ~ 1, ~ 1), d, family = "degpd",
                      restarts = FALSE, control = control), "control")
  }
})
