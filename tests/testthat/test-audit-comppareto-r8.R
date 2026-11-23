## An independent log likelihood, from the two conditional densities and their
## mixing weights. Do not use the fitter's density, derivatives or profile here.
cp8_ll <- function(b, y, spec) {
  p <- length(b)
  th <- exp(b[p])
  a <- exp(b[p-1L])
  body <- switch(spec,
    exp = list(d = dexp, p = pexp, args = list(rate = exp(b[1L]))),
    lnorm = list(d = dlnorm, p = plnorm, args = list(meanlog = b[1L], sdlog = exp(b[2L]))),
    gamma = list(d = dgamma, p = pgamma, args = list(shape = exp(b[1L]), scale = exp(b[2L]))),
    weibull = list(d = dweibull, p = pweibull, args = list(shape = exp(b[1L]), scale = exp(b[2L]))))
  lg <- function(v) do.call(body$d, c(list(v), body$args, log = TRUE))
  G <- do.call(body$p, c(list(th), body$args))
  phi <- 2 * th * exp(lg(th)) / (a * G)
  lo <- y <= th
  sum(lg(pmin(y, th)) - log(G) - log1p(phi)) -
    (a+1) * sum(log1p((y[!lo] - th)/(2*th)))
}

test_that("a resolvable tail score is not rounded away at an exponential splice", {
  y <- 1 + .Machine$double.eps
  b <- c(0, log(1e16), 0)
  a <- exp(b[2L])
  ## At theta = rate = 1 the tail NLL depends on log(alpha) only through
  ## log(1 + C/alpha) + (alpha+1) * log(1 + (y-1)/2).
  w <- 2 / (2 + a * expm1(1))
  tail <- a * log1p((y-1)/2)
  d <- egpd:::.comppareto_d120(y, matrix(b, nrow = 1), "exp")
  expect_equal(d[1L, 2L], tail - w, tolerance = 1e-13)
  expect_equal(d[1L, 7L], tail + w*(1-w), tolerance = 1e-13)
})

test_that("a finite exponential splice kink passes both sides, unlike nearby false knots", {
  set.seed(54002)
  y <- rcomppareto(400, "exp", rate = 1, alpha = 3, theta = 2)
  ## Isolate the convergence check: even the single start gets to this kink.
  f <- suppressWarnings(egpd(list(y ~ 1, ~1, ~1), data.frame(y), family = "comppareto",
                             comppareto.args = list(spec = "exp"), restarts = FALSE))
  ll <- function(p) cp8_ll(p, y, "exp")
  b <- unname(f$coefficients)
  expect_equal(f$convergence, 0)
  expect_equal(ll(b), as.numeric(logLik(f)), tolerance = 1e-12)
  expect_lt(abs(ll(b) + 444.093698), 1e-5)
  h <- 1e-6
  delta <- diag(h, 3)
  g <- vapply(1:2, function(i) (ll(b+delta[i, ]) - ll(b-delta[i, ]))/(2*h), 0)
  expect_lt(max(abs(g)), 1e-4)
  outward <- function(p) (c(ll(p-delta[3, ]), ll(p+delta[3, ])) - ll(p))/h
  expect_lt(max(outward(b)), -0.1)

  ## At one nearby observation moving right increases the likelihood; at
  ## another it is moving left. First optimize both nuisance parameters there
  ## so their scores, proximity to a knot, or checking only one side cannot pass.
  knots <- sort(log(y))
  k <- which.min(abs(knots - b[3L]))
  for (d in c(-3L, 3L)) {
    x <- knots[k + d]
    o <- stats::optim(b[1:2], function(v) -ll(c(v, x)), method = "BFGS",
                      control = list(reltol = 1e-14))
    p <- c(o$par, x)
    slopes <- outward(p)
    expect_gt(slopes[if (d < 0L) 2L else 1L], 0.3)
    bad <- list(pars = p, objective = o$value, convergence = 3L, kept = rep(TRUE, 3))
    checked <- egpd:::.comppareto_finish_exp(bad, egpd:::.nllh.nopen, egpd:::.search.nopen,
                                             f$likdata, f$likfns)
    expect_equal(checked$convergence, 3L)
  }
  ## A sentinel on the other branch is not evidence that moving there is uphill.
  expect_false(egpd:::.comppareto_kink_sides(0, -1:1, function(x) if (x > 0) 1e20 else abs(x)))
  ## Here both sides of the raw likelihood still pass, but the actual penalised
  ## objective decreases to the right. It must not be certified using raw scores.
  ld <- f$likdata
  ld$S <- diag(c(0, 0, 50))
  obj <- function(p) egpd:::.nllh.pen(p, ld, f$likfns)
  expect_lt(obj(b+delta[3, ]), obj(b))
  bad <- list(pars = b, objective = obj(b), convergence = 3L, kept = rep(TRUE, 3))
  checked <- egpd:::.comppareto_finish_exp(bad, egpd:::.nllh.pen, egpd:::.search.pen, ld, f$likfns)
  expect_equal(checked$convergence, 3L)
})

test_that("a uniform-body boundary is not reported as a finite exponential kink optimum", {
  set.seed(11501)
  y <- regpd(150, sigma = 2, xi = 0.15, kappa = 1.5)
  f <- suppressWarnings(egpd(list(y ~ 1, ~1, ~1), data.frame(y), family = "comppareto",
                             comppareto.args = list(spec = "exp")))
  expect_gt(f$convergence, 0L)
  expect_match(f$convergence.report, "boundary.*rate tends to zero")
  expect_true(is.na(f$logLik))
  b <- unname(f$coefficients)
  expect_gt(cp8_ll(c(-20, b[-1]), y, "exp"), cp8_ll(c(-5, b[-1]), y, "exp"))
})

test_that("an iteration limit at a splice needs the same stationarity certificate", {
  set.seed(31503)
  y <- regpd(150, sigma = 3, xi = .05, kappa = 2)
  ## Directional stationarity at a kink does not imply a positive-definite
  ## ordinary Hessian. Do not mask an optimizer warning.
  f <- withCallingHandlers(egpd(list(y ~ 1, ~1, ~1), data.frame(y), family = "comppareto",
    comppareto.args = list(spec = "exp"), restarts = FALSE),
    egpd_singular_hessian = function(w) invokeRestart("muffleWarning"))
  b <- unname(f$coefficients)
  ll <- function(p) cp8_ll(p, y, "exp")
  expect_equal(f$convergence, 0)
  expect_lt(abs(ll(b) + 372.3665333), 1e-5)
  expect_gt(ll(b), max(ll(b+c(0,0,1e-6)), ll(b-c(0,0,1e-6))))

  ## Budget exhaustion away from the verified knot is not convergence.
  early <- suppressWarnings(egpd(list(y ~ 1, ~1, ~1), data.frame(y), family = "comppareto",
    comppareto.args = list(spec = "exp"), restarts = FALSE, control = list(inner = list(itlim = 1))))
  expect_identical(early$convergence, 1L)
})

test_that("the split screen finds better interior optima, including between observations", {
  ## Round-7 independent multi-start reference values (not package fits). The
  ## original defaults have observations on *both* sides, so the four stranded
  ## starts alone do nothing for these losses of 4.435, 1.260 and 3.017.
  for (spec in c("lnorm", "gamma", "weibull")) {
    set.seed(switch(spec, lnorm = 2150, gamma = 1150, weibull = 1400))
    y <- switch(spec,
      lnorm = rcomppareto(150, spec, meanlog = .2, sdlog = .35, alpha = 1.5, theta = 2),
      gamma = rcomppareto(150, spec, shape = 2, scale = 1, alpha = 2, theta = 4),
      weibull = rcomppareto(400, spec, shape = 1.5, scale = 2, alpha = 2, theta = 3))
    ref <- switch(spec, lnorm = -258.6029, gamma = -325.7198, weibull = -871.9676)
    f <- egpd(list(y ~ 1, ~1, ~1, ~1), data.frame(y), family = "comppareto",
              comppareto.args = list(spec = spec))
    b <- unname(f$coefficients)
    expect_equal(f$convergence, 0)
    expect_lt(abs(cp8_ll(b, y, spec) - ref), 0.01)
    expect_equal(cp8_ll(b, y, spec), f$logLik, tolerance = 1e-12)
    if (spec == "weibull") expect_gt(min(abs(log(y) - b[4L])), 0.01)
  }
})

test_that("splice searches respect offsets, pins, omitted rows and supplied starts", {
  set.seed(2150)
  y <- rcomppareto(150, "lnorm", meanlog = .2, sdlog = .35, alpha = 1.5, theta = 2)
  d <- data.frame(y = c(y, max(y)*10), off = c(rep(.7, 150), NA))
  fit <- function(...) suppressWarnings(egpd(list(y ~ 1, ~1, ~1, ~1+offset(off)), d,
    family = "comppareto", fixed = list(sdlog = .35), compact = TRUE, removeData = TRUE, ...))
  f <- fit()
  expect_equal(f$nobs, 150)
  b <- vapply(1:4, function(i) f[[i]]$fitted[1L], 0)
  expect_equal(b[2L], log(.35))
  expect_equal(b[4L], unname(f$coefficients[4L]) + .7)
  expect_equal(cp8_ll(b, y, "lnorm"), f$logLik, tolerance = 1e-12)
  ## Specified intercepts, including an explicit offset on theta, skip the scan.
  init <- c(.5, 0, log(3), log(1.5) - .7)
  expect_identical(fit(inits = init)$coefficients,
                   fit(inits = init, restarts = FALSE)$coefficients)
})

test_that("a selected all-body split distinguishes a free tail boundary from a pin", {
  set.seed(498)
  y <- runif(50, .1, 2)
  fit <- function(fixed) suppressWarnings(egpd(list(y ~ 1, ~1, ~1), data.frame(y),
    family = "comppareto", comppareto.args = list(spec = "exp"), fixed = fixed))
  f <- fit(list(theta = 3, rate = .8))
  expect_gt(f$convergence, 0L)
  expect_match(f$convergence.report, "no finite tail-shape maximum")
  expect_true(is.na(f$logLik))
  ## With alpha pinned it is no longer an ascent direction.
  g <- fit(list(theta = 3, alpha = 1.5))
  expect_equal(g$convergence, 0)
  expect_true(is.finite(g$logLik))
})

test_that("a boundary supersedes optimizer diagnostics but not an omitted-row warning", {
  set.seed(77)
  y <- rgamma(60, 2, 1)
  d <- data.frame(y = c(y, NA), x = rep(0:1, length.out = 61))
  d$z <- d$x
  ## Pin the splice above every observation: there can be no finite free alpha.
  ## The unclassified fit either hits its iteration limit or has a singular Hessian.
  for (limited in c(TRUE, FALSE)) {
    body <- if (limited) y ~ 1 else y ~ x + z
    ws <- list()
    f <- withCallingHandlers(egpd(list(body, ~1, ~1, ~1), d, family = "comppareto",
      comppareto.args = list(spec = "gamma"), fixed = list(theta = max(y) + 1),
      control = if (limited) list(inner = list(itlim = 1))),
      warning = function(w) { ws[[length(ws) + 1L]] <<- w; invokeRestart("muffleWarning") })
    expect_identical(vapply(ws, function(w) class(w)[1L], ""),
                     c("simpleWarning", "egpd_nonconvergence"))
    expect_match(conditionMessage(ws[[1L]]), "Omitting 1")
    expect_match(conditionMessage(ws[[length(ws)]]), f$convergence.report, fixed = TRUE)
  }
})
