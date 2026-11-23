test_that("fit_spliced_degpd fits and predicts a spliced count model", {
  skip_on_cran()
  skip_if_not_installed("qgam")

  set.seed(101)
  n <- 300
  x <- runif(n, -1, 1)
  y <- rpois(n, lambda = exp(1 + 0.8 * x)) +
    rbinom(n, 1, 0.12) * rdiscegpd(n, sigma = exp(0.8 + 0.5 * x), xi = 0.15, kappa = 1.3)
  df <- data.frame(y = y, x = x)

  fit <- suppressWarnings(
    fit_spliced_degpd(
      y ~ s(x, k = 6),
      data = df,
      bulk_probs = c(0.1, 0.5, 0.75, 0.9),
      threshold_prob = 0.75,
      tail_formula = list(lsigma = ~x, lxi = ~1, lkappa = ~1),
      keep_data = TRUE
    )
  )

  expect_s3_class(fit, "spliced_degpd")
  expect_s3_class(summary(fit), "summary.spliced_degpd")
  expect_equal(fit$threshold_prob, 0.75)
  expect_true(fit$n_exceed > 0)

  qhat <- predict(fit, type = "quantile", prob = c(0.25, 0.75, 0.95))
  expect_equal(dim(qhat), c(nrow(df), 3L))
  expect_true(all(qhat[[2]] >= qhat[[1]]))
  expect_true(all(qhat[[3]] >= qhat[[2]]))

  pars <- predict(fit, type = "parameter")
  expect_true(all(c("threshold", "threshold_quantile", "scale", "shape", "kappa") %in% names(pars)))
  expect_true(all(pars$threshold >= 0))

  cdf <- predict(fit, newdata = df[1:3, ], type = "cdf", at = 0:40)
  expect_equal(dim(cdf), c(3L, 41L))
  expect_true(all(as.matrix(cdf) >= 0))
  expect_true(all(as.matrix(cdf) <= 1))
  expect_true(all(apply(as.matrix(cdf), 1, function(z) all(diff(z) >= -1e-10))))

  pmf <- predict(fit, newdata = df[1:3, ], type = "pmf", at = 0:60)
  expect_equal(dim(pmf), c(3L, 61L))
  expect_true(all(as.matrix(pmf) >= -1e-12))
  expect_true(all(rowSums(as.matrix(pmf)) > 0.85))
  expect_true(all(rowSums(as.matrix(pmf)) <= 1 + 1e-8))
})

test_that("fit_spliced_degpd works without a tail model", {
  skip_if_not_installed("qgam")

  set.seed(102)
  x <- runif(150, -1, 1)
  y <- rpois(150, lambda = exp(0.6 + 0.4 * x))
  df <- data.frame(y = y, x = x)

  fit <- fit_spliced_degpd(
    y ~ x,
    data = df,
    bulk_probs = c(0.2, 0.5, 0.8),
    threshold_prob = 0.8,
    fit_tail = FALSE,
    keep_data = TRUE
  )

  expect_null(fit$tail_fit)

  bulk <- predict(fit, type = "bulk")
  expect_equal(dim(bulk), c(nrow(df), 3L))

  cdf <- predict(fit, newdata = df[1:2, ], type = "cdf", at = 0:20)
  expect_equal(dim(cdf), c(2L, 21L))
})

splice_geom_data <- function() {
  set.seed(81)
  n <- 320
  d <- data.frame(x = runif(n, -1, 1), z = rnorm(n), off = rep(c(0, log(4)), n / 2))
  d$y <- rpois(n, exp(1 + d$x)) +
    rbinom(n, 1, .25) * rgeom(n, plogis(-1 - .65 * d$z - .4 * d$x^2))
  d
}

splice_geom_fit <- function(d, form, probs = c(.2, .6)) {
  set.seed(140)
  fit_spliced_degpd(y ~ x, d, bulk_probs = probs, threshold_prob = .6, tail_formula = form,
                    egpd_control = list(fix.arg = list(xi = 0, kappa = 1)), keep_data = TRUE)
}

splice_time_limit <- function(expr) {
  setTimeLimit(elapsed = 10, transient = TRUE)
  on.exit(setTimeLimit(transient = TRUE))
  force(expr)
}

test_that("tail rank uses factor levels observed in the training data", {
  skip_if_not_installed("qgam")
  d <- splice_geom_data()
  d$g <- ifelse(d$x > 0, "east", "west")
  d$g[1:12] <- "dry"
  d$y[1:12] <- 0L
  for (g in list(d$g, factor(d$g))) {
    d$g <- g
    expect_error(splice_geom_fit(d, ~g), "rank deficient.*complete exceedances")
  }
  ## A declared but unused level does not make the fitted regression deficient.
  d$g <- factor(ifelse(d$x > 0, "east", "west"), levels = c("east","west","dry"))
  f <- splice_geom_fit(d, ~g)
  expect_true(all(is.finite(predict(f, type = "parameter")$scale)))
})

test_that("the spliced likelihood uses the same complete exceedances as the scale design", {
  skip_if_not_installed("qgam")
  d <- splice_geom_data()
  f <- splice_geom_fit(d, ~z)
  y <- d$y - predict(f, type = "parameter")$threshold
  idx <- which(y >= 0)
  ll <- function(fit, exc, z) {
    b <- unname(fit$tail_fit$coefficients)
    sum(dgeom(exc, prob = -expm1(-exp(-b[1] - b[2] * z)), log = TRUE))
  }
  expect_equal(f$tail_fit$loglik, ll(f, y[idx], d$z[idx]), tolerance = 1e-10)

  miss <- idx[round(length(idx) / 2)]
  d$z[miss] <- NA
  expect_warning(fm <- splice_geom_fit(d, ~z), "Omitted 1 incomplete tail observation")
  keep <- setdiff(idx, miss)
  expect_equal(fm$n_exceed, length(keep))
  expect_identical(fm$tail_fit$convergence, 0L)
  expect_equal(fm$tail_fit$loglik, ll(fm, y[keep], d$z[keep]), tolerance = 1e-10)
  recycled <- ll(fm, y[idx], rep(d$z[keep], length.out = length(idx)))
  expect_gt(abs(recycled - fm$tail_fit$loglik), .01)
  expect_identical(which(is.na(predict(fm, type = "parameter")$scale)), miss)

  tail_data <- d[idx, ]
  tail_data$exceedance <- y[idx]
  tail_data$z[-seq_len(9)] <- NA
  expect_warning(expect_error(egpd:::.fit_spliced_degpd_tail(fm$tail_formula, tail_data),
                              "Not enough complete exceedances"), "9 exceedances remain")

  nd <- d[10, , drop = FALSE]
  par <- predict(f, nd, type = "parameter")
  k <- par$threshold + c(ceiling(44 * par$scale), 0, 1, 0)
  mass <- as.numeric(predict(f, nd, type = "pmf", at = k))
  want <- (1 - par$alpha_tilde) * dgeom(k - par$threshold, -expm1(-1 / par$scale))
  expect_equal(mass / want, rep(1, length(k)), tolerance = 1e-10)
  expect_equal(as.numeric(predict(f, nd, type = "cdf", at = k[1] - 0:1)), c(1, 1))

  ## An ordinary published CDF boundary and the next doubles on either side.
  kk <- par$threshold + 2
  pk <- predict(f, nd, type = "cdf", at = kk)[[1]]
  expect_true(.5 < pk && pk < 1)
  expect_equal(as.numeric(predict(f, nd, type = "quantile", prob = pk + (-1:1) * 2^-53)),
               c(kk, kk, kk + 1))
  p <- 1 - c(2^-52, 2^-53)
  want <- par$threshold + ceiling(-par$scale * log((1 - p) / (1 - par$alpha_tilde))) - 1
  expect_equal(as.numeric(predict(f, nd, type = "quantile", prob = c(p, 1))), c(want, Inf))
})

test_that("a one-knot splice preserves the scale formula and trained prediction basis", {
  skip_if_not_installed("qgam")
  d <- splice_geom_data()
  d$grp <- factor(ifelse(d$z > 0, "A", "B"))
  form <- local({
    base_off <- log(1.25)
    ~0 + grp + poly(x, 2) + offset(off + base_off)
  })
  fit <- splice_geom_fit(d, form, probs = .6)
  expect_identical(environment(fit$tail_formula$lsigma), environment(form))
  expect_equal(dim(predict(fit, type = "bulk")), c(nrow(d), 1L))
  y <- d$y - predict(fit, type = "parameter")$threshold
  train <- d[y >= 0, ]
  X <- model.matrix(form, train)
  b <- fit$tail_fit$coefficients
  expect_identical(names(b), colnames(X))
  eta <- drop(X %*% b) + train$off + log(1.25)
  expect_equal(fit$tail_fit$loglik,
               sum(dgeom(y[y >= 0], prob = -expm1(-exp(-eta)), log = TRUE)), tolerance = 1e-10)

  nd <- data.frame(x = .1, grp = "B", off = 0, row.names = "alone")
  basis <- predict(poly(train$x, 2), nd$x)
  want <- exp(b[["grpB"]] + drop(basis %*% b[3:4]) + log(1.25))
  one <- predict(fit, nd, type = "parameter")
  expect_identical(rownames(one), "alone")
  expect_equal(one$scale, want, ignore_attr = TRUE)
  nd2 <- rbind(nd, transform(nd, off = log(4)), transform(nd, x = -.7, grp = "A"))
  many <- predict(fit, nd2, type = "parameter")
  expect_equal(many$scale[1:2] / one$scale, c(1, 4), tolerance = 1e-12)
  expect_error(predict(fit, transform(nd, grp = "new"), type = "parameter"), "new level")
})

test_that("bounded and heavy fitted splice percentiles invert the survivor, including endpoints", {
  skip_if_not_installed("qgam")
  set.seed(173)
  d <- data.frame(x = runif(250))
  d$y <- rpois(nrow(d), exp(1 + .2 * d$x))
  nd <- d[1, , drop = FALSE]
  for (xi in c(-.05, 2)) {
    set.seed(140)
    f <- fit_spliced_degpd(y ~ x, d, bulk_probs = c(.2, .6), threshold_prob = .6,
                           egpd_control = list(fix.arg = list(xi = xi, kappa = 1)))
    expect_identical(f$tail_fit$convergence, 0L)
    par <- predict(f, nd, type = "parameter")
    p <- c(1 - 1e-10, 1 - 2^-53)
    surv <- (1 - p) / (1 - par$alpha_tilde)
    want <- par$threshold + ceiling(par$scale / xi * expm1(-xi * log(surv))) - 1
    last <- if (xi < 0) par$threshold + ceiling(-par$scale / xi) - 1 else Inf
    got <- splice_time_limit(as.numeric(predict(f, nd, type = "quantile", prob = c(0, p, 1))))
    expect_equal(got[c(1, 4)], c(0, last))
    if (xi < 0) {
      expect_equal(got[2:3], want)
    } else {
      expect_true(all(is.finite(got[2:3]) & got[2:3] > 2^53))
      expect_equal(got[2:3] / want, c(1, 1), tolerance = 1e-13)
    }
  }
})

test_that("bulk inversion makes progress on widely spaced representable counts", {
  big <- .Machine$double.xmax
  q <- splice_time_limit(egpd:::.spliced_degpd_row_quantile(c(.25, 1), big, .5, 1))
  expect_identical(q, c(big / 2, big))
})

test_that("a one-ulp move into the tail survives rescaling at the splice", {
  alpha <- .1
  p <- alpha + 2^-56
  expect_identical((1 - p) / (1 - alpha), 1)
  scale <- 1e18
  q <- egpd:::.spliced_degpd_row_quantile(c(alpha, p), 1, alpha, alpha, has_tail = TRUE,
                                          shift = 2, alpha_tilde = alpha, scale = scale)
  want <- 2 + ceiling(-scale * log1p(-(p - alpha) / (1 - alpha))) - 1
  expect_equal(q, c(1, want))
})
