## Round-5 coordinator regression tests (integration follow-ups).

test_that("summary() reports t2() fits' parametric terms in mgcv's prediction basis", {
  skip_if_not_installed("mgcv")
  ## egpd keeps t2() coefficients in mgcv's fit basis, which mixes the intercept into
  ## the smooth; summary.gam reports the prediction basis (beta = P beta_fit). With a
  ## unit-variance Gaussian likelihood at mgcv's own smoothing parameters, the
  ## intercept row must equal summary.gam's.
  set.seed(817)
  n <- 180
  d <- data.frame(x = runif(n), z = runif(n))
  d$y <- sin(6 * d$x) + 0.5 * d$z + rnorm(n)
  gauss <- list(
    d0 = function(pars, likdata)
      -sum(dnorm(likdata$y[, 1], as.vector(likdata$X[[1]] %*% pars), 1, log = TRUE)),
    d120 = function(pars, likdata)
      cbind(as.vector(likdata$X[[1]] %*% pars) - likdata$y[, 1], rep(1, nrow(likdata$y))),
    d340 = function(pars, likdata) matrix(0, nrow(likdata$y), 2))
  f <- y ~ t2(x, z, k = c(4, 4))
  m <- mgcv::gam(f, data = d, scale = 1, method = "REML")
  e <- egpd(list(f), d, family = "custom", custom.fns = gauss, sp = m$sp)
  ref <- summary(m)$p.table
  got <- summary(e)[[1]][[1]]
  expect_equal(unname(got[1, 1]), unname(ref[1, 1]), tolerance = 1e-6)
  expect_equal(unname(got[1, 2]), unname(ref[1, 2]), tolerance = 1e-6)
})
