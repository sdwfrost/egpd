## Round 6 (hesswarn): a singular Hessian behind the standard errors is recorded on the fit
## (`hessian.rank`), noted by print() and summary(), and warned about (condition class
## "egpd_singular_hessian") where the standard errors cannot show it.
##
## References are independent of the package: the generalized Waring log pmf and its
## scores are written out from lgamma()/digamma() (Irwin's law; a = mu (rho - 1) / k), the
## outer-product (BHHH) matrix is formed from those scores, and design ranks come from qr().

hw_gw_logpmf <- function(y, mu, k, rho) {
  a <- mu * (rho - 1) / k
  lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) + lgamma(a + y) - lgamma(a) +
    lgamma(k + y) - lgamma(k) - lgamma(a + k + rho + y) - lgamma(y + 1)
}

## per-observation scores of the log pmf with respect to (log mu, log k, log rho)
hw_gw_scores <- function(y, mu, k, rho) {
  a <- mu * (rho - 1) / k
  da <- digamma(a + rho) + digamma(a + y) - digamma(a) - digamma(a + k + rho + y)
  dk <- digamma(k + rho) + digamma(k + y) - digamma(k) - digamma(a + k + rho + y)
  dr <- digamma(a + rho) + digamma(k + rho) - digamma(rho) - digamma(a + k + rho + y)
  cbind(a * da, -a * da + k * dk, a * rho / (rho - 1) * da + rho * dr)
}

## the same for the zero-truncated law: log f(y) - log(1 - f(0))
hw_zt_scores <- function(y, mu, k, rho) {
  f0 <- exp(hw_gw_logpmf(0, mu, k, rho))
  hw_gw_scores(y, mu, k, rho) +
    f0 / (1 - f0) * matrix(hw_gw_scores(0, mu, k, rho), length(y), 3, byrow = TRUE)
}

hw_conditions <- function(expr) {
  ws <- list()
  ms <- list()
  val <- withCallingHandlers(expr,
    warning = function(w) { ws[[length(ws) + 1L]] <<- w; invokeRestart("muffleWarning") },
    message = function(m) { ms[[length(ms) + 1L]] <<- m; invokeRestart("muffleMessage") })
  list(value = val, warnings = ws, messages = ms,
       classes = vapply(ws, function(w) class(w)[1L], ""))
}

test_that("the generalized Waring at a = k warns, and its outer-product Hessian is singular", {
  expect_equal(hw_gw_logpmf(0:40, 3, 1.5, 2.2),
               dgwaring(0:40, mu = 3, k = 1.5, rho = 2.2, log = TRUE), tolerance = 1e-12)
  set.seed(1)
  y <- rgwaring(500, mu = 4, k = 2, rho = 2.5)
  r <- hw_conditions(egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = y), family = "gwaring"))
  f <- r$value
  expect_identical(r$classes, "egpd_singular_hessian")
  expect_identical(f$convergence, 0L)
  expect_identical(f$hessian.rank, c(rank = 2L, npar = 3L))
  msg <- conditionMessage(r$warnings[[1L]])
  expect_match(msg, "singular at the estimates: numerical rank 2 of 3", fixed = TRUE)
  expect_match(msg, "outer product of the per-observation scores (BHHH)", fixed = TRUE)
  expect_match(msg, "mainly in logk_(Intercept)", fixed = TRUE)
  expect_match(msg, "confint(method = \"profile\")", fixed = TRUE)

  b <- unname(f$coefficients)
  mu <- exp(b[1L]); k <- exp(b[2L]); rho <- exp(b[3L])
  ## the estimate sits on the symmetry line a = k of the pmf
  expect_equal(mu * (rho - 1) / k, k, tolerance = 1e-5)
  G <- hw_gw_scores(y, mu, k, rho)
  expect_lt(max(abs(colSums(G))), 1e-3)
  ## the log-k score vanishes for every observation (the a and k scores coincide), so the
  ## outer product of the scores has rank 2
  expect_lt(sqrt(sum(G[, 2L]^2)) / sqrt(sum(G[, 1L]^2)), 1e-5)
  B <- crossprod(G)
  ev <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(ev[3L] / ev[1L], 1e-12)
  expect_gt(ev[2L] / ev[1L], 1e-2)
  ## the fit's Hessian is that outer product, and Vp, its pseudo-inverse, gives the null
  ## direction either no variance or an enormous one: neither is its uncertainty
  expect_equal(egpd:::.gH.nopen(b, f$likdata, f$likfns)[[2L]], B, tolerance = 1e-8)
  v <- eigen(B, symmetric = TRUE)$vectors[, 3L]
  vv <- drop(crossprod(v, f$Vp %*% v))
  s <- max(diag(f$Vp)[-2L])
  expect_true(abs(vv) < 1e-10 * s || vv > 1e6 * s)
})

test_that("near a = k the warning fires although the pseudo-inverse keeps the direction", {
  set.seed(3)
  y <- rgwaring(500, mu = 4, k = 2, rho = 2.5)
  expect_warning(f <- egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = y), family = "gwaring"),
                 class = "egpd_singular_hessian")
  expect_identical(f$hessian.rank, c(rank = 2L, npar = 3L))
  b <- unname(f$coefficients)
  B <- crossprod(hw_gw_scores(y, exp(b[1L]), exp(b[2L]), exp(b[3L])))
  e <- eigen(B, symmetric = TRUE)
  expect_lt(e$values[3L] / e$values[1L], sqrt(.Machine$double.eps))
  ## a variance along the null direction of 1e5 and more times any other
  vv <- drop(crossprod(e$vectors[, 3L], f$Vp %*% e$vectors[, 3L]))
  expect_gt(vv, 1e5 * max(diag(f$Vp)[-2L]))
})

test_that("fit_hurdle() reports the positive part's warning once, for the start it keeps", {
  set.seed(3)
  y <- rgwaring(400, mu = 4, k = 3, rho = 1.6)
  n_w <- 0
  h <- withCallingHandlers(fit_hurdle(y ~ 1, data = data.frame(y = y), family = "gwaring"),
                           egpd_singular_hessian = function(w) {
                             n_w <<- n_w + 1
                             invokeRestart("muffleWarning")
                           })
  expect_equal(n_w, 1)
  p <- h$positive
  expect_identical(p$convergence, 0L)
  expect_identical(p$hessian.rank, c(rank = 2L, npar = 3L))
  b <- unname(p$coefficients)
  G <- hw_zt_scores(y[y > 0], exp(b[1L]), exp(b[2L]), exp(b[3L]))
  expect_lt(max(abs(colSums(G))), 1e-3)
  expect_lt(sqrt(sum(G[, 2L]^2)) / sqrt(sum(G[, 1L]^2)), 1e-4)
  ev <- eigen(crossprod(G), symmetric = TRUE, only.values = TRUE)$values
  expect_lt(ev[3L] / ev[1L], 1e-12)
  out <- capture.output(print(summary(h)))
  expect_true(any(grepl("Hessian behind the standard errors is singular (numerical rank 2 of 3",
                        out, fixed = TRUE)))
  ## print() of the hurdle fit has the same note, pointing at the positive part's rank
  out <- capture.output(print(h))
  expect_true(any(grepl(paste("Hessian behind the standard errors is singular (numerical rank",
                              "2 of 3 estimated coefficients; `positive$hessian.rank`"),
                        out, fixed = TRUE)))
  ## and an ordinary hurdle fit (its estimate has a = 7.5, k = 0.83) has none
  set.seed(3)
  y <- rgwaring(400, mu = 4, k = 2, rho = 2.5)
  r <- hw_conditions(fit_hurdle(y ~ 1, data = data.frame(y = y), family = "gwaring"))
  expect_length(r$warnings, 0L)
  expect_identical(r$value$positive$hessian.rank, c(rank = 3L, npar = 3L))
  expect_false(any(grepl("Hessian", capture.output(print(r$value)))))
})

test_that("collinear terms warn once, also through the restarts of egpd()", {
  set.seed(9)
  n <- 200
  x <- runif(n)
  d <- data.frame(y = rnbinom(n, mu = exp(1 + x), size = 3), x = x, x2 = 2 * x)
  X <- model.matrix(~ x + x2, d)
  expect_identical(qr(X)$rank, 2L)
  r <- hw_conditions(egpd(list(y ~ x + x2, ~ 1), data = d, family = "pig"))
  expect_identical(r$classes, "egpd_singular_hessian")
  ## the observed information is exactly singular: Vp gives the aliased direction no variance
  expect_match(conditionMessage(r$warnings[[1L]]), "(Vp gives them zero variance)", fixed = TRUE)
  ## the warning replaces the message that a Hessian failing chol() used to give
  expect_length(r$messages, 0L)
  expect_identical(r$value$hessian.rank, c(rank = 3L, npar = 4L))
  v <- c(0, 2, -1, 0) / sqrt(5)
  expect_lt(abs(drop(crossprod(v, r$value$Vp %*% v))), 1e-10)
  ## DEGPD fits five starts, each with the same aliased design: one warning, the kept fit's
  rd <- hw_conditions(egpd(list(y ~ x + x2, ~ 1, ~ 1), data = d, family = "degpd"))
  expect_identical(sum(rd$classes == "egpd_singular_hessian"), 1L)
  expect_identical(rd$value$hessian.rank[["npar"]], 5L)
  expect_lt(rd$value$hessian.rank[["rank"]], 5L)
})

test_that("ordinary fits do not warn: other estimates, covariate units, large penalties", {
  set.seed(2)
  y <- rgwaring(500, mu = 4, k = 2, rho = 2.5)
  r <- hw_conditions(egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = y), family = "gwaring"))
  expect_length(r$warnings, 0L)
  expect_identical(r$value$hessian.rank, c(rank = 3L, npar = 3L))
  b <- unname(r$value$coefficients)
  ev <- eigen(crossprod(hw_gw_scores(y, exp(b[1L]), exp(b[2L]), exp(b[3L]))),
              symmetric = TRUE, only.values = TRUE)$values
  expect_gt(ev[3L] / ev[1L], 1e-3)
  expect_false(any(grepl("Hessian", capture.output(print(r$value)))))
  expect_false(any(grepl("Hessian", capture.output(print(summary(r$value))))))

  ## a covariate in large units (calendar years): the Hessian of the coefficients has an
  ## eigenvalue ratio far below sqrt(eps), but in units of the covariate's spread it is
  ## well conditioned
  set.seed(5)
  n <- 300
  year <- rep(2001:2020, length.out = n)
  dy <- data.frame(y = rnbinom(n, mu = exp(1 + 0.05 * (year - 2010)), size = 3), year = year)
  r <- hw_conditions(egpd(list(y ~ year, ~ 1), data = dy, family = "pig"))
  expect_length(r$warnings, 0L)
  expect_identical(r$value$hessian.rank, c(rank = 3L, npar = 3L))
  ev <- eigen(solve(r$value$Vp), symmetric = TRUE, only.values = TRUE)$values
  expect_lt(ev[3L] / ev[1L], sqrt(.Machine$double.eps))

  ## a smooth with a very large smoothing parameter: H = H0 + S has an eigenvalue ratio
  ## below sqrt(eps), yet every direction is determined by the data or the penalty
  set.seed(101)
  n <- 400
  x <- runif(n)
  ds <- data.frame(y = rgwaring(n, mu = exp(1 + 0.8 * x), k = 2, rho = 3), x = x)
  r <- hw_conditions(egpd(list(y ~ s(x, k = 10), ~ 1), data = ds, family = "pig", sp = 1e8))
  expect_length(r$warnings, 0L)
  expect_identical(r$value$hessian.rank, c(rank = 11L, npar = 11L))
  ev <- eigen(solve(r$value$Vp), symmetric = TRUE, only.values = TRUE)$values
  expect_lt(ev[11L] / ev[1L], sqrt(.Machine$double.eps))
})

test_that("print() and summary() note a singular Hessian; a plateau is not warned about", {
  set.seed(6)
  y <- rgwaring(500, mu = 4, k = 2, rho = 2.5)
  expect_warning(f <- egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = y), family = "gwaring"),
                 class = "egpd_singular_hessian")
  note <- "Note: the Hessian behind the standard errors is singular (numerical rank 2 of 3"
  expect_true(any(grepl(note, capture.output(print(f)), fixed = TRUE)))
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl(note, out, fixed = TRUE)))
  ## under the coefficient tables
  expect_gt(max(grep(note, out, fixed = TRUE)), max(grep("(Intercept)", out, fixed = TRUE)))
  expect_identical(summary(f)$hessian.rank, f$hessian.rank)

  ## negative binomial data drive rho past the generalized Waring's range: log rho is on a
  ## plateau of its link (`plateau`), already reported there and by print(), so the
  ## singular Hessian is noted but not warned about
  set.seed(21)
  n <- 300
  x <- runif(n)
  dn <- data.frame(y = rnbinom(n, mu = exp(1 + 0.8 * x), size = 2), x = x)
  r <- hw_conditions(egpd(list(y ~ x, ~ 1, ~ 1), data = dn, family = "gwaring"))
  expect_length(r$warnings, 0L)
  g <- r$value
  expect_identical(unname(g$plateau), c(FALSE, FALSE, TRUE))
  ## log rho is beyond the range in which the likelihood is evaluated, so its score, its
  ## row of the Hessian and its variance in Vp are zero
  b <- unname(g$coefficients)
  expect_gt(exp(b[4L]), egpd:::cf_bounds_cpp()$gwaring[["rho_hi"]])
  expect_identical(g$Vp[4L, 4L], 0)
  expect_identical(g$hessian.rank, c(rank = 3L, npar = 4L))
  expect_true(any(grepl("Hessian behind the standard errors is singular (numerical rank 3 of 4",
                        capture.output(print(g)), fixed = TRUE)))
})

test_that("a fit that did not converge warns only that it did not converge", {
  ## the verifier's EGPD fit on the xi < -1 singularity (test-audit-coreopt-r4.R, VA-10):
  ## its Hessian is not positive definite, which used to add a message as well
  set.seed(2)
  xp <- runif(40)
  d <- data.frame(y = regpd(40, sigma = exp(0.5 + 0.6 * xp), xi = 0.15, kappa = 1.5, type = 1),
                  x1 = xp, x2 = xp^2, x3 = xp^3, x4 = sin(2 * pi * xp))
  r <- hw_conditions(egpd(list(y ~ x1 + x2 + x3 + x4, ~1, ~ x1 + x2), data = d, family = "egpd"))
  expect_identical(r$value$convergence, 3L)
  expect_true("egpd_nonconvergence" %in% r$classes)
  expect_false("egpd_singular_hessian" %in% r$classes)
  expect_length(r$messages, 0L)
  H <- egpd:::.gH.nopen(unname(r$value$coefficients), r$value$likdata, r$value$likfns)[[2L]]
  expect_lt(min(eigen(H, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_lt(r$value$hessian.rank[["rank"]], r$value$hessian.rank[["npar"]])
})

test_that("profile intervals at a = k are those of the likelihood, not of the singular Hessian", {
  set.seed(1)
  y <- rgwaring(500, mu = 4, k = 2, rho = 2.5)
  expect_warning(f <- egpd(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = y), family = "gwaring"),
                 class = "egpd_singular_hessian")
  ## independent profile likelihood of log k: the generalized Waring log-likelihood written
  ## stably for any a (lbeta for log f(0) and the finite product
  ## prod_{j < y} (a + j) (k + j) / ((a + k + rho + j) (j + 1)) for the rest), maximised
  ## over the nuisance parameters (log mu, log(rho - 1)) at each fixed log k
  tab <- table(y)
  yv <- as.integer(names(tab))
  nv <- as.vector(tab)
  loglik <- function(mu, k, rho, yv, nv) {
    a <- mu * (rho - 1) / k
    j <- seq_len(max(yv)) - 1
    cum <- c(0, cumsum(log1p(-(k + rho) / (a + k + rho + j)) + log(k + j) - log1p(j)))
    sum(nv * (lbeta(a + rho, k) - lbeta(rho, k) + cum[yv + 1]))
  }
  expect_equal(vapply(0:40, function(v) loglik(3, 1.5, 2.2, v, 1), 0),
               dgwaring(0:40, mu = 3, k = 1.5, rho = 2.2, log = TRUE), tolerance = 1e-12)
  nll <- function(th) {
    v <- -loglik(exp(th[1]), exp(th[2]), 1 + exp(th[3]), yv, nv)
    if (is.finite(v)) v else 1e100
  }
  o <- optim(c(log(mean(y)), log(2), log(1.5)), nll, control = list(reltol = 1e-12, maxit = 5000))
  o <- optim(o$par, nll, method = "BFGS", control = list(reltol = 1e-15, maxit = 5000))
  expect_equal(-o$value, f$logLik, tolerance = 1e-8)
  dev <- function(lk) {
    p <- optim(o$par[-2], function(t) nll(c(t[1], lk, t[2])), method = "BFGS",
               control = list(reltol = 1e-15, maxit = 5000))
    2 * (p$value - o$value) - qchisq(0.95, 1)
  }
  ref <- exp(c(uniroot(dev, o$par[2] + c(-1, 0), tol = 1e-10)$root,
               uniroot(dev, o$par[2] + c(0, 1.5), tol = 1e-10)$root))
  expect_equal(round(ref, 3), c(1.063, 8.519))
  ## the package's interval was [0, Inf]: its bracketing step was the standard error of
  ## log k from the singular Hessian, 1.6e-8
  ci <- confint(f, method = "profile")
  expect_equal(unname(ci["k", ]), ref, tolerance = 1e-3)
  ## a fit without `hessian.rank` (as from an older version) takes the same step, as the
  ## variance of log k is negligible
  f0 <- f
  f0$hessian.rank <- NULL
  expect_identical(confint(f0, parm = "k", method = "profile"), ci["k", , drop = FALSE])
})
