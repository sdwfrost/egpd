test_that("grouped BI coverage uses fitted and held-out trial sizes in row order", {
  skip_if_not_installed("gamlss")
  d <- data.frame(b = c(20, 5, 10, 20, 5, 10), y = c(11, 2, 5, 13, 4, 7), x = 1:6)
  f <- gamlss::gamlss(cbind(y, b - y) ~ 1, data = d, family = gamlss.dist::BI(),
                      trace = FALSE)
  p <- sum(d$y) / sum(d$b)
  se <- 1 / sqrt(sum(d$b) * p * (1 - p))
  nsim <- 4000L
  for (m in c("plug-in", "parametric")) for (insample in c(TRUE, FALSE)) {
    b <- if (insample) d$b else c(42, 0, 9, 2)
    y <- if (insample) d$y else c(25, 0, 8, 1)
    nd <- if (!insample) data.frame(x = seq_along(b))
    set.seed(12)
    r <- predictive_coverage(f, y, nd, bd = if (!insample) b, data_fit = d,
                             method = m, nsim = nsim)
    ## Independent binomial inverse, using the same uniforms and the
    ## intercept-only binomial information for the normal logit perturbation.
    set.seed(12)
    prob <- if (m == "parametric") plogis(rnorm(nsim * length(b), qlogis(p), se)) else p
    sim <- matrix(qbinom(runif(nsim * length(b)), rep(b, each = nsim), prob), nsim)
    ref <- apply(sim, 2, quantile, probs = c(.025, .975), type = 1, names = FALSE)
    expect_equal(r$L, ref[1, ])
    expect_equal(r$U, ref[2, ])
    expect_equal(r$covered, y >= ref[1, ] & y <= ref[2, ])
  }
})

test_that("held-out binomial denominators are explicit and counts are validated", {
  skip_if_not_installed("gamlss")
  d <- data.frame(y = rep(0:1, 4), b = 1, x = 1:8)
  nd <- data.frame(x = 4:5)
  fit <- function(grouped) {
    gamlss::gamlss(if (grouped) cbind(y, b - y) ~ 1 else y ~ 1, data = d,
                   family = gamlss.dist::BI(), trace = FALSE)
  }
  f <- fit(TRUE)  # two-column training response, even though all sizes are one
  for (m in c("parametric", "plug-in"))
    expect_error(predictive_coverage(f, c(0, 1), nd, method = m), "Supply `bd`")
  for (bad in list(1:3, .5, -1, NA_real_, Inf))
    expect_error(predictive_coverage(f, c(0, 1), nd, bd = bad), "non-negative integer")
  for (bad in list(c(0, 3), c(1, -.5), c(NA_real_, 0)))
    expect_error(predictive_coverage(f, bad, nd, bd = 2), "success counts")
  expect_error(predictive_coverage(f, d$y, bd = 2), "only with `newdata`")
  f$bd <- f$bd[-1L]
  expect_error(predictive_coverage(f, d$y), "per fitted row")

  f <- fit(FALSE)
  set.seed(11)
  r <- predictive_coverage(f, 0:1, nd, data_fit = d, method = "plug-in", nsim = 200)
  expect_identical(r$L, c(0, 0))
  expect_identical(r$U, c(1, 1))
  g <- gamlss::gamlss(y ~ 1, data = d, family = gamlss.dist::PO(), trace = FALSE)
  expect_error(predictive_coverage(g, 0:1, nd, bd = 4), "does not accept")
})

test_that("beta-binomial coverage uses its own law, not BI or the RNG's default size", {
  skip_if_not_installed("gamlss")
  d <- data.frame(b = c(17, 4, 9, 27, 8, 2), y = c(9, 1, 4, 16, 5, 0))
  f <- gamlss::gamlss(cbind(y, b - y) ~ 1, data = d, family = gamlss.dist::BB(),
                      sigma.fix = TRUE, sigma.start = .25, trace = FALSE)
  nsim <- 300L
  set.seed(25)
  r <- predictive_coverage(f, d$y, nsim = nsim, method = "plug-in")
  set.seed(25)
  ref <- vapply(seq_along(d$b), function(j) {
    sim <- gamlss.dist::qBB(runif(nsim), mu = f$mu.fv[j], sigma = .25, bd = d$b[j])
    quantile(sim, c(.025, .975), type = 1, names = FALSE)
  }, numeric(2))
  expect_equal(r$L, ref[1, ])
  expect_equal(r$U, ref[2, ])
})

test_that("zero-adjusted binomial families require a possible positive count", {
  skip_if_not_installed("gamlss")
  d <- data.frame(y = c(0, 3, 5, 0, 2, 4), b = c(7, 6, 9, 2, 5, 8))
  for (fam in list(gamlss.dist::ZABI(), gamlss.dist::ZABB())) {
    f <- gamlss::gamlss(cbind(y, b - y) ~ 1, data = d, family = fam, trace = FALSE,
                        sigma.fix = TRUE, sigma.start = .2, nu.fix = TRUE, nu.start = .3)
    expect_error(predictive_coverage(f, 0, newdata = d[1, ], bd = 0), "positive trial")
  }
})
