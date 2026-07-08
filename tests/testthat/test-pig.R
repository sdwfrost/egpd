## Tests for the Poisson-inverse Gaussian (PIG) and zero-inflated PIG (ZIPIG)
## families across the native GAM fitter (egpd()), the distfit() MLE frontend,
## the re-exported gamlss families, and the bamlss constructors.
##
## The native GAM d12 supplies the exact score AND the exact observed Hessian
## (analytic, verified against finite differences below).

skip_if_not_installed("gamlss.dist")

## coefficient-space gradient and Hessian assembled from d120, vs numerical
## differences of d0 (mirrors tests/testthat/test-gam-d12-fd.R).
fd_oracle_pig <- function(fns, npar, y, p) {
  n <- length(y)
  X <- c(list(cbind(1, seq(-1, 1, length.out = n))),
         replicate(npar - 1, matrix(1, n, 1), simplify = FALSE))
  nbk <- vapply(X, ncol, 1L)
  ld <- list(X = X, y = cbind(y), dupid = 0L, duplicate = 0L,
             offsets = replicate(npar, numeric(0), simplify = FALSE),
             idpars = rep(seq_len(npar), nbk), censored = FALSE)
  cg <- function(pp) {
    M <- fns$d120(pp, ld)
    unlist(lapply(seq_len(npar), function(i) colSums(M[, i] * X[[i]])))
  }
  col_ij <- function(a, b) {            # packed upper-tri column for 0-based (a,b)
    off <- 0L; if (a > 0) for (r in 0:(a - 1)) off <- off + (npar - r)
    npar + off + (b - a) + 1L
  }
  ch <- function(pp) {
    M <- fns$d120(pp, ld); st <- c(0, cumsum(nbk)); nb <- sum(nbk)
    H <- matrix(0, nb, nb)
    for (a in 0:(npar - 1)) for (b in a:(npar - 1)) {
      cc <- col_ij(a, b); Hab <- t(X[[a + 1]]) %*% (X[[b + 1]] * M[, cc])
      H[(st[a + 1] + 1):st[a + 2], (st[b + 1] + 1):st[b + 2]] <- Hab
      if (a != b) H[(st[b + 1] + 1):st[b + 2], (st[a + 1] + 1):st[a + 2]] <- t(Hab)
    }
    H
  }
  ng <- vapply(seq_along(p), function(i) {
    h <- 1e-6; pp <- p; pm <- p; pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
    (fns$d0(pp, ld) - fns$d0(pm, ld)) / (2 * h)
  }, numeric(1))
  nh <- matrix(0, length(p), length(p))
  for (i in seq_along(p)) {
    h <- 1e-5; pp <- p; pm <- p; pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
    nh[, i] <- (cg(pp) - cg(pm)) / (2 * h)
  }
  nh <- (nh + t(nh)) / 2
  list(grad = max(abs(cg(p) - ng)), hess = max(abs(ch(p) - nh)))
}

ydisc <- c(0, 0, 1, 2, 3, 5, 8, 12, 4, 1)

test_that("native PIG d0 matches gamlss.dist::dPIG log-likelihood", {
  fns <- get(".pig1fns", envir = asNamespace("egpd"))
  n <- length(ydisc)
  X <- list(matrix(1, n, 1), matrix(1, n, 1))
  ld <- list(X = X, y = cbind(ydisc), dupid = 0L, duplicate = 0L,
             offsets = replicate(2, numeric(0), simplify = FALSE),
             idpars = 1:2, censored = FALSE)
  lmu <- 0.6; lsig <- -0.3
  nll <- fns$d0(c(lmu, lsig), ld)
  ref <- -sum(gamlss.dist::dPIG(ydisc, mu = exp(lmu), sigma = exp(lsig), log = TRUE))
  expect_equal(nll, ref, tolerance = 1e-8)
})

test_that("native PIG/ZIPIG analytic gradient AND observed Hessian match FD", {
  pigfns <- get(".pig1fns", envir = asNamespace("egpd"))
  zipfns <- get(".zipig1fns", envir = asNamespace("egpd"))
  rp <- fd_oracle_pig(pigfns, 2, ydisc, c(0.7, 0.3, -0.3))
  expect_lt(rp$grad, 1e-4); expect_lt(rp$hess, 1e-2)
  rz <- fd_oracle_pig(zipfns, 3, ydisc, c(0.7, 0.3, -0.3, 0.2))
  expect_lt(rz$grad, 1e-4); expect_lt(rz$hess, 1e-2)
})

test_that("native egpd() PIG recovers parameters and matches direct MLE", {
  skip_on_cran()
  set.seed(42)
  y <- gamlss.dist::rPIG(4000, mu = 3.5, sigma = 0.8)
  fit <- egpd(list(y ~ 1, ~ 1), data = data.frame(y = y), family = "pig", trace = 0)
  co <- coef(fit)
  expect_equal(exp(co[[1]]), 3.5, tolerance = 0.1)
  expect_equal(exp(co[[2]]), 0.8, tolerance = 0.1)
  ## coefficients coincide with the direct gamlss.dist MLE
  nll <- function(t) -sum(gamlss.dist::dPIG(y, mu = exp(t[1]), sigma = exp(t[2]), log = TRUE))
  o <- optim(c(log(mean(y)), 0), nll, method = "BFGS")
  expect_equal(unname(co), unname(o$par), tolerance = 1e-2)
})

test_that("native egpd() ZIPIG recovers parameters", {
  skip_on_cran()
  set.seed(11)
  yz <- gamlss.dist::rZIPIG(4000, mu = 4, sigma = 0.7, nu = 0.3)
  fz <- egpd(list(yz ~ 1, ~ 1, ~ 1), data = data.frame(yz = yz), family = "zipig", trace = 0)
  cz <- coef(fz)
  expect_equal(exp(cz[[1]]), 4,   tolerance = 0.15)
  expect_equal(exp(cz[[2]]), 0.7, tolerance = 0.15)
  expect_equal(plogis(cz[[3]]), 0.3, tolerance = 0.06)
})

test_that("distfit() PIG MLE matches gamlss.dist and yields working methods", {
  skip_on_cran()
  set.seed(7)
  y <- gamlss.dist::rPIG(3000, mu = 3, sigma = 0.9)
  f <- distfit(y, family = "pig")
  expect_s3_class(f, "distfit")
  expect_named(f$estimate, c("mu", "sigma"))
  o <- optim(c(log(mean(y)), 0),
             function(t) -sum(gamlss.dist::dPIG(y, exp(t[1]), exp(t[2]), log = TRUE)))
  expect_equal(f$loglik, -o$value, tolerance = 0.05)
  expect_equal(as.numeric(logLik(f)), f$loglik, tolerance = 1e-8)
  expect_equal(AIC(f), -2 * f$loglik + 2 * 2, tolerance = 1e-6)
})

test_that("distfit() ZIPIG MLE recovers the zero-inflation probability", {
  skip_on_cran()
  set.seed(3)
  yz <- gamlss.dist::rZIPIG(3000, mu = 3.5, sigma = 0.8, nu = 0.25)
  fz <- distfit(yz, family = "zipig")
  expect_named(fz$estimate, c("mu", "sigma", "pi"))
  expect_equal(fz$estimate[["pi"]], 0.25, tolerance = 0.06)
})

test_that("gamlss PIG/ZIPIG families are re-exported", {
  expect_s3_class(PIG(), "gamlss.family")
  expect_s3_class(ZIPIG(), "gamlss.family")
})

test_that("pig_bamlss()/zipig_bamlss() build valid family.bamlss objects", {
  fp <- pig_bamlss(); fz <- zipig_bamlss()
  expect_s3_class(fp, "family.bamlss")
  expect_s3_class(fz, "family.bamlss")
  expect_equal(fp$names, c("mu", "sigma"))
  expect_equal(fz$names, c("mu", "sigma", "pi"))
  expect_equal(fp$d(2, list(mu = 3, sigma = 0.9)),
               gamlss.dist::dPIG(2, mu = 3, sigma = 0.9))
  expect_equal(fz$d(0, list(mu = 3, sigma = 0.9, pi = 0.2)),
               gamlss.dist::dZIPIG(0, mu = 3, sigma = 0.9, nu = 0.2))
})

## ----------------------------------------------------------------------------
## Distribution-function checks:
##   (a) the native C++ likelihood implements a proper PMF (sums to 1) and
##       matches the reference gamlss.dist implementation across a grid;
##   (b) the d / p / q functions the frontends rely on are mutually consistent;
##   (c) analytic derivatives match finite differences across parameter regimes.
## ----------------------------------------------------------------------------

## native single-observation PMF from d0 (= exp(-nll) for each count y)
pmf_native <- function(fns, npar, ys, p) {
  X <- replicate(npar, matrix(1, 1, 1), simplify = FALSE)
  base <- list(X = X, dupid = 0L, duplicate = 0L,
               offsets = replicate(npar, numeric(0), simplify = FALSE),
               idpars = seq_len(npar), censored = FALSE)
  vapply(ys, function(y) exp(-fns$d0(p, utils::modifyList(base, list(y = cbind(y))))),
         numeric(1))
}

pgrid <- expand.grid(mu = c(0.5, 2, 8), sigma = c(0.2, 1, 5))

test_that("native PIG likelihood is a proper PMF and matches gamlss.dist::dPIG", {
  fns <- get(".pig1fns", envir = asNamespace("egpd"))
  for (i in seq_len(nrow(pgrid))) {
    mu <- pgrid$mu[i]; sg <- pgrid$sigma[i]; p <- c(log(mu), log(sg))
    ## match the reference implementation over the low counts
    expect_equal(pmf_native(fns, 2, 0:49, p),
                 gamlss.dist::dPIG(0:49, mu, sg), tolerance = 1e-9)
  }
  ## normalisation: the native PMF sums to 1 (independent of gamlss.dist)
  for (i in c(1, 5, 9)) {
    mu <- pgrid$mu[i]; sg <- pgrid$sigma[i]
    expect_equal(sum(pmf_native(fns, 2, 0:6000, c(log(mu), log(sg)))), 1,
                 tolerance = 1e-6)
  }
})

test_that("native ZIPIG likelihood is a proper PMF and matches gamlss.dist::dZIPIG", {
  fns <- get(".zipig1fns", envir = asNamespace("egpd"))
  for (pi0 in c(0.1, 0.3, 0.6)) {
    mu <- 3; sg <- 0.8; p <- c(log(mu), log(sg), qlogis(pi0))
    expect_equal(pmf_native(fns, 3, 0:49, p),
                 gamlss.dist::dZIPIG(0:49, mu, sg, nu = pi0), tolerance = 1e-9)
    expect_equal(sum(pmf_native(fns, 3, 0:6000, p)), 1, tolerance = 1e-6)
  }
})

test_that("PIG/ZIPIG d, p, q are mutually consistent (as used by the frontends)", {
  mu <- 3.2; sg <- 0.7; pi0 <- 0.25; qs <- 0:80
  ## CDF equals the cumulative sum of the PMF
  expect_equal(gamlss.dist::pPIG(qs, mu, sg),
               cumsum(gamlss.dist::dPIG(qs, mu, sg)), tolerance = 1e-8)
  expect_equal(gamlss.dist::pZIPIG(qs, mu, sg, nu = pi0),
               cumsum(gamlss.dist::dZIPIG(qs, mu, sg, nu = pi0)), tolerance = 1e-8)
  ## quantile function inverts the CDF: q(p(k)) == k
  expect_equal(gamlss.dist::qPIG(gamlss.dist::pPIG(qs, mu, sg), mu, sg), qs)
  expect_equal(gamlss.dist::qZIPIG(gamlss.dist::pZIPIG(qs, mu, sg, nu = pi0), mu, sg, nu = pi0),
               qs)
  ## the native prediction/quantile helper (.iG_pig) agrees with qPIG
  iG <- get(".iG_pig", envir = asNamespace("egpd"))
  u <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  expect_equal(iG(u, mu = rep(mu, length(u)), sigma = rep(sg, length(u))),
               gamlss.dist::qPIG(u, mu, sg))
})

test_that("dpig/ppig/qpig/rpig wrappers delegate to gamlss.dist and are consistent", {
  mu <- 2.6; sg <- 0.8
  expect_equal(dpig(0:40, mu, sg), gamlss.dist::dPIG(0:40, mu, sg))
  expect_equal(ppig(0:40, mu, sg), gamlss.dist::pPIG(0:40, mu, sg))
  expect_equal(qpig(c(.1, .5, .9), mu, sg), gamlss.dist::qPIG(c(.1, .5, .9), mu, sg))
  expect_equal(dpig(3, mu, sg, log = TRUE), log(dpig(3, mu, sg)))
  expect_equal(ppig(5, mu, sg, lower.tail = FALSE), 1 - ppig(5, mu, sg))
  ## d / p / q mutually consistent through the egpd-native names
  expect_equal(ppig(0:40, mu, sg), cumsum(dpig(0:40, mu, sg)), tolerance = 1e-8)
  expect_equal(qpig(ppig(0:40, mu, sg), mu, sg), 0:40)
  set.seed(1); expect_length(rpig(50, mu, sg), 50)
})

test_that("dzipig/pzipig/qzipig/rzipig wrappers delegate and are consistent", {
  mu <- 3; sg <- 0.7; pv <- 0.3
  expect_equal(dzipig(0:40, mu, sg, pi = pv), gamlss.dist::dZIPIG(0:40, mu, sg, nu = pv))
  expect_equal(pzipig(0:40, mu, sg, pi = pv), gamlss.dist::pZIPIG(0:40, mu, sg, nu = pv))
  expect_equal(dzipig(0, mu, sg, pi = pv), pv + (1 - pv) * dpig(0, mu, sg))  # zero mixture
  expect_equal(pzipig(0:40, mu, sg, pi = pv), cumsum(dzipig(0:40, mu, sg, pi = pv)),
               tolerance = 1e-8)
  expect_equal(qzipig(pzipig(0:40, mu, sg, pi = pv), mu, sg, pi = pv), 0:40)
  set.seed(1); expect_length(rzipig(50, mu, sg, pi = pv), 50)
})

test_that("native PIG/ZIPIG derivatives match FD across parameter regimes", {
  pigfns <- get(".pig1fns", envir = asNamespace("egpd"))
  zipfns <- get(".zipig1fns", envir = asNamespace("egpd"))
  ## (intercept, slope, log-sigma[, logit-pi]) covering small/large mu & sigma
  pig_pts <- list(c(0.7, 0.3, -0.3), c(2.0, 0.2, -1.5),
                  c(-0.5, 0.4, 1.5), c(1.5, -0.3, 0.5))
  for (p in pig_pts) {
    r <- fd_oracle_pig(pigfns, 2, ydisc, p)
    expect_lt(r$grad, 1e-4); expect_lt(r$hess, 1e-2)
  }
  zip_pts <- list(c(0.7, 0.3, -0.3, 0.2), c(1.8, 0.2, -1.2, -0.8),
                  c(-0.4, 0.4, 1.2, 0.9))
  for (p in zip_pts) {
    r <- fd_oracle_pig(zipfns, 3, ydisc, p)
    expect_lt(r$grad, 1e-4); expect_lt(r$hess, 1e-2)
  }
})
