## G1: mgcv's penalty matrices need not have independent multipliers.
g1_data <- function() {
  set.seed(317)
  d <- data.frame(x = runif(260), z = runif(260))
  d$y <- 1.4 * sin(4*pi*d$x) + .7 * cos(2*pi*d$z) + rnorm(nrow(d))
  d$w <- runif(nrow(d))
  d
}

g1_fit <- function(f, d, ...) {
  eta <- function(b, ld) drop(ld$X[[1L]] %*% b)
  gauss <- list(d0 = function(b, ld) -sum(dnorm(ld$y[, 1L], eta(b, ld), 1, log = TRUE)),
                d120 = function(b, ld) cbind(eta(b, ld) - ld$y[, 1L], 1),
                d340 = function(b, ld) matrix(0, nrow(ld$y), 2))
  egpd(list(mu = f), d, family = "custom", custom.fns = gauss,
       control = list(inner = list(gradtol = 1e-10),
                      outer = list(gradtol = 1e-7, fntol = 1e-12)), ...)
}

## Only the basis comes from mgcv. At *specified full multipliers*, solve
## the Gaussian normal equations and evaluate the restricted likelihood
## with a dense pseudo-determinant, independently of egpd's fit and map.
g1_dense <- function(G, sp, gamma = 1) {
  X <- G$X
  XX <- crossprod(X)
  S <- matrix(0, ncol(X), ncol(X))
  for (j in seq_along(G$S)) {
    i <- G$off[j] - 1L + seq_len(nrow(G$S[[j]]))
    S[i, i] <- S[i, i] + sp[j] * G$S[[j]]
  }
  H <- XX + S
  V <- solve(H)
  b <- drop(solve(H, crossprod(X, G$y)))
  d <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  d <- d[d > max(d) * 1e-10]
  Mp <- ncol(X) - length(d)
  nll <- -sum(dnorm(G$y, X %*% b, 1, log = TRUE))
  reml <- (nll + .5 * sum(b * (S %*% b))) / gamma +
    .5 * (as.numeric(determinant(H/gamma)$modulus) - sum(log(d/gamma)) - Mp * log(2*pi))
  list(b = b, S = S, Vp = V, edf = diag(V %*% XX), reml = reml, Mp = Mp)
}

test_that("a term-fixed smooth is penalised even with no outer parameters", {
  d <- g1_data()
  f <- y ~ s(x, k = 8, sp = 5)
  G <- mgcv::gam(f, data = d, fit = FALSE)
  ref <- g1_dense(G, 5, gamma = 3)
  fit <- g1_fit(f, d, gamma = 3)
  top <- g1_fit(y ~ s(x, k = 8), d, sp = 5)
  expect_length(fit$sp, 0)
  expect_equal(unname(fit$mu$full.sp), 5)
  expect_equal(unname(fit$full.sp), unname(top$full.sp))
  expect_lt(max(abs(coef(fit) - ref$b)), 1e-8)
  expect_gt(max(abs(coef(fit) - qr.coef(qr(G$X), d$y))), 1)
  expect_lt(max(abs(fit$mu$edf - ref$edf)), 1e-8)
  expect_lt(max(abs(fit$Vp - ref$Vp)), 1e-8)
  expect_identical(fit$Vc, fit$Vp)
  expect_lt(abs(fit$negREML - ref$reml), 1e-8)
  expect_lt(max(abs(coef(top) - ref$b)), 1e-8)
  expect_error(g1_fit(f, d, sp = 4), "conflicts")
})

test_that("REML and covariance use the independent parameters, including with a fixed block", {
  d <- g1_data()
  ## Shared id; fixed + shared at gamma != 1; overlapping tensor penalties.
  cases <- list(
    list(y ~ s(x, k = 8, id = 1) + s(z, k = 8, id = 1), 1, function(s) c(s, s)),
    list(y ~ s(x, k = 8, id = 1) + s(z, k = 8, id = 1) + s(w, k = 5, sp = .7),
         2.4, function(s) c(s, s, .7)),
    list(y ~ te(x, z, k = 4, sp = c(.5, -1)), 1, function(s) c(.5, s)))
  for (case in cases) {
    f <- case[[1L]]; g <- case[[2L]]; weights <- case[[3L]]
    G <- mgcv::gam(f, data = d, fit = FALSE)
    ref <- function(r) g1_dense(G, weights(exp(r)), g)
    score <- function(r) ref(r)$reml
    rstar <- optimize(score, c(-10, 5), tol = 1e-9)$minimum
    m <- mgcv::gam(f, data = d, method = "REML", scale = 1, gamma = g)
    fit <- g1_fit(f, d, gamma = g)
    expect_length(fit$sp, 1L)
    r <- log(fit$sp)
    expect_lt(abs(r - rstar), 2e-4)
    expect_lt(abs(log(m$sp) - rstar), 2e-4)
    expect_equal(unname(fit$mu$full.sp), weights(fit$sp), tolerance = 1e-12)
    a <- ref(r)
    expect_lt(max(abs(coef(fit) - a$b)), 2e-6)
    expect_lt(abs(fit$negREML - a$reml), 1e-8)
    expect_lt(abs(fit$negREML + a$Mp * log(g) - m$gcv.ubre), 1e-6)
    expect_lt(max(abs(fit$mu$edf - a$edf)), 1e-8)
    expect_lt(max(abs(fit$Vp - a$Vp)), 1e-8)
    ## Differentiate the independent profile and coefficients, not egpd's
    ## gradient/Hessian. Fixed contributions affect H but have no variance.
    h <- .001
    plus <- ref(r + h); minus <- ref(r - h)
    v <- h^2 / (plus$reml - 2*a$reml + minus$reml)
    J <- (plus$b - minus$b) / (2*h)
    expect_lt(abs(as.numeric(fit$Vlsp) - v), 5e-5)
    expect_lt(max(abs(fit$Vc - (a$Vp + v * tcrossprod(J)))), 1e-5)
    expect_gt(max(abs(fit$Vc - fit$Vp)), 1e-4)
    ## Check the sign and magnitude of the analytic score away from its zero.
    r1 <- r + .4
    grad <- egpd:::.reml1(r1 - log(g), fit$likfns, fit$likdata, fit$Sdata, beta = ref(r1)$b)
    expect_lt(abs(grad - (score(r1 + h) - score(r1 - h))/(2*h)), 1e-5)
    replay <- g1_fit(f, d, sp = fit$full.sp)
    expect_lt(max(abs(coef(replay) - a$b)), 1e-7)
    expect_identical(replay$Vc, replay$Vp)
  }
})

test_that("fixed overrides respect ids and term sp; exact zero removes rank", {
  d <- g1_data()
  f <- y ~ s(x, k = 8, id = 1) + s(z, k = 8, id = 1)
  expect_error(g1_fit(f, d, sp = c(.1, .2)), "conflicts")
  expect_error(g1_fit(f, d, sp = rep(.1, 3)), "Expected 1 independent")
  expect_error(g1_fit(f, d, rho0 = c(0, 0)), "per independent")
  expect_error(g1_fit(y ~ s(x, id = 1) + s(z, id = 1, sp = 2), d), "first smooth")
  for (sp in list(c(0, 0), c(0, .4))) {
    f <- y ~ te(x, z, k = 4)
    a <- g1_dense(mgcv::gam(f, data = d, fit = FALSE), sp)
    fit <- g1_fit(f, d, sp = if (all(sp == 0)) 0 else sp)
    expect_equal(unname(fit$full.sp), sp)
    expect_lt(max(abs(coef(fit) - a$b)), 1e-8)
    expect_lt(abs(fit$df - sum(a$edf)), 1e-8)
    expect_lt(abs(fit$negREML - a$reml), 1e-8)
    expect_identical(fit$Vc, fit$Vp)
  }
  ## mgcv's *term-level* zero is a tiny positive multiplier; its determinant
  ## and rank must use the same value as the objective and reported full.sp.
  f <- y ~ s(x, k = 8, sp = 0)
  m <- mgcv::gam(f, data = d, method = "REML", scale = 1)
  fit <- g1_fit(f, d)
  expect_equal(log(unname(fit$full.sp)), log(unname(m$full.sp)), tolerance = 1e-12)
  expect_lt(abs(fit$negREML - m$gcv.ubre), 1e-8)
})

test_that("a built-in factor-by shares within scale, not with another predictor's id", {
  set.seed(7)
  n <- 600
  fac <- factor(sample(c("a", "b", "c"), n, TRUE))
  x <- runif(n)
  sig <- exp(1 + .8 * sin(2*pi*x) + .3 * (fac == "b"))
  d <- data.frame(y = rdiscegpd(n, sigma = sig, xi = .15, kappa = 1.2, type = 1), x, fac)
  shape <- local({ kk <- 5; ~ s(x, k = kk, id = 1) })
  model <- function(scale, ...) egpd(list(scale, shape, ~ 1), d, family = "degpd",
                                     restarts = FALSE, ...)
  fit <- model(y ~ fac + s(x, by = fac, k = 6, id = 1))
  expect_equal(fit$convergence, 0L)
  expect_length(fit$sp, 2L)
  expect_equal(unname(fit$logscale$sp), fit$sp[1L])
  expect_equal(unname(fit$logshape$sp), fit$sp[2L])
  expect_equal(unname(fit$logscale$full.sp), rep(fit$sp[1L], 3))
  expect_equal(unname(fit$full.sp), fit$sp[c(1, 1, 1, 2)])
  untied <- y ~ fac + s(x, by = fac, k = 6)
  imposed <- model(untied, sp = fit$full.sp)
  expect_lt(max(abs(coef(fit) - coef(imposed))), 1e-4)
  expect_lt(max(abs(fit$Vp - imposed$Vp)), 1e-4)
  expect_lt(abs(fit$df - imposed$df), 1e-4)
  expect_lt(abs(fit$negREML - imposed$negREML), 1e-5)
  ## No free scale sp precedes the shape sp here, but the three scale
  ## penalties and their labels must survive in the predictor's report.
  lambda <- fit$sp[1L]
  fixed <- model(y ~ fac + s(x, by = fac, k = 6, id = 1, sp = lambda), sp = fit$sp[2L])
  expect_length(fixed$logscale$sp, 0)
  expect_equal(fixed$logscale$full.sp, fit$logscale$full.sp)
  expect_equal(fixed$logshape$sp, fit$logshape$sp)
  expect_lt(max(abs(coef(fixed) - coef(imposed))), 1e-4)
  free <- model(untied)
  expect_gt(max(abs(as.matrix(fitted(fit)) - as.matrix(fitted(free)))), .01)
})
