## Run separately from testthat.R, which attaches egpd and would mask lookup failures.
if (requireNamespace("gamlss", quietly = TRUE)) local({
  stopifnot(!"package:egpd" %in% search())
  d <- data.frame(y = rep(c(0, 1, 1, 2, 3, 7), 12), x = rep(c(0, 0, 0, 1, 1, 1), 12))
  bell <- gamlss::gamlss(y ~ x, data = d, family = egpd::BELL(), trace = FALSE)
  disc <- gamlss::gamlss(y ~ 1, data = d, family = egpd::DEGPD1(), trace = FALSE,
                         sigma.fix = TRUE, sigma.start = .2, nu.fix = TRUE, nu.start = 1)
  equal <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-9))
  for (f in list(bell, disc)) {
    name <- f$family[1L]
    density <- getExportedValue("egpd", paste0("d", name))
    cdf <- getExportedValue("egpd", paste0("p", name))
    par <- lapply(stats::setNames(f$parameters, f$parameters), function(p) f[[paste0(p, ".fv")]])
    stopifnot(f$converged,
              equal(deviance(f), -2 * sum(do.call(density, c(list(d$y, log = TRUE), par)))))
    nd <- data.frame(x = c(1, .3, 0))
    pred <- predict(f, newdata = nd, data = d, type = "response")
    X <- model.matrix(delete.response(f$mu.terms), nd)
    stopifnot(equal(pred, exp(X %*% coef(f))))

    ## Re-evaluate the *stored* rqres expression after fitting.
    set.seed(87)
    z <- gamlss::get.rqres(f, howmany = 1)
    set.seed(87)
    lo <- do.call(cdf, c(list(d$y - 1), par))
    hi <- do.call(cdf, c(list(d$y), par))
    stopifnot(equal(pnorm(z), runif(nrow(d), lo, hi)))
  }

  ## Exact exponential and floored exponential tails. At these last rows a
  ## linear CDF is already one; clipping cannot give the correct residual.
  ed <- data.frame(y = c(-1.6*log1p(-((1:140 - .5)/140)), 100, 145))
  ef <- gamlss::gamlss(y ~ 1, data = ed, family = egpd::EGPD1(sigma.link = "identity"),
                       sigma.fix = TRUE, sigma.start = 0, nu.fix = TRUE, nu.start = 1,
                       mu.start = mean(ed$y), trace = FALSE)
  ez <- qnorm(-ed$y/ef$mu.fv, lower.tail = FALSE, log.p = TRUE)
  stopifnot(ef$converged, equal(ef$residuals, ez), equal(gamlss::get.rqres(ef, 1), ez))

  dd <- data.frame(y = c(rep(0:4, c(60, 40, 15, 12, 8)), 120, 155))
  gf <- gamlss::gamlss(y ~ 1, data = dd, family = egpd::DEGPD1(sigma.link = "identity"),
                       sigma.fix = TRUE, sigma.start = 0, nu.fix = TRUE, nu.start = 1,
                       mu.start = 1/log1p(1/mean(dd$y)), trace = FALSE)
  ls <- -dd$y/gf$mu.fv
  lower <- qnorm(ls, lower.tail = FALSE, log.p = TRUE)
  upper <- qnorm(ls - 1/gf$mu.fv, lower.tail = FALSE, log.p = TRUE)
  stopifnot(gf$converged, all(gf$residuals >= lower & gf$residuals <= upper))
  for (seed in c(28, 73)) {
    set.seed(seed); z <- gamlss::get.rqres(gf, 1)
    set.seed(seed); w <- runif(nrow(dd))
    stopifnot(equal(z, qnorm(ls + log1p(w * expm1(-1/gf$mu.fv)), lower.tail = FALSE,
                             log.p = TRUE)))
  }

  md <- data.frame(y = c(0, ed$y, 0))
  mf <- gamlss::gamlss(y ~ 1, data = md, family = egpd::ZIEGPD1(sigma.link = "identity"),
                       sigma.fix = TRUE, sigma.start = 0, nu.fix = TRUE, nu.start = 1,
                       tau.fix = TRUE, tau.start = .25, mu.start = mean(ed$y), trace = FALSE)
  mz <- qnorm(log(.75) - md$y/mf$mu.fv, lower.tail = FALSE, log.p = TRUE)
  atom <- md$y == 0
  stopifnot(mf$converged, equal(mf$residuals[!atom], mz[!atom]),
            all(pnorm(mf$residuals[atom]) > 0 & pnorm(mf$residuals[atom]) < .25))
  set.seed(75); mz[atom] <- qnorm(.25 * runif(nrow(md))[atom])
  set.seed(75); stopifnot(equal(gamlss::get.rqres(mf, 1), mz))

  ## True continuous endpoints remain infinite; a tiny nonzero lower tail
  ## need not fit in a double on probability scale to give a finite quantile.
  rz <- eval(egpd::EGPD1()$rqres, list(y = c(0, 1e-250, 4), mu = 1, sigma = -.25, nu = 2),
             asNamespace("gamlss"))
  stopifnot(equal(rz, c(-Inf, qnorm(2*log(1e-250), log.p = TRUE), Inf)))

  ## Coverage no longer depends on upstream vcov's bare density/constructor
  ## lookup. Even the boundary xi=0 needs no finite default-link coordinate.
  for (m in c("parametric", "plug-in")) for (fit in list(bell, ef)) {
    dat <- if (identical(fit$family[1], "BELL")) d else ed
    r <- egpd::predictive_coverage(fit, y = 2, newdata = dat[1, , drop = FALSE],
                                   data_fit = dat, method = m, nsim = 100)
    stopifnot(is.finite(r$L), is.finite(r$U))
  }

  ## Exercise the other generated callback shapes in the evaluator's namespace
  ## (not in the family's closure). In particular a mixed-family zero must be
  ## randomized over its atom; it is not an ordinary continuous observation.
  for (fam in list(egpd::ZIBELL(), egpd::GPIG(), egpd::ZIGPIG(), egpd::EGPD1(),
                   egpd::EGPD4(), egpd::ZIEGPD1(), egpd::ZIDEGPD5())) {
    y <- if (fam$type == "Continuous") c(.5, 1, 3) else c(0, 1, 3)
    par <- list(mu = 2, sigma = .3, nu = .6, tau = .2)[names(fam$parameters)]
    par <- lapply(par, rep, length(y))
    vals <- c(list(y = y), par)
    foreign_eval <- function(expr) eval(expr, vals, asNamespace("gamlss"))
    density <- getExportedValue("egpd", paste0("d", fam$family[1]))
    cdf <- getExportedValue("egpd", paste0("p", fam$family[1]))
    stopifnot(equal(foreign_eval(body(fam$G.dev.incr)),
                    -2 * do.call(density, c(list(y, log = TRUE), par))),
              all(is.finite(foreign_eval(body(fam$dldm)))))
    set.seed(11)
    z <- foreign_eval(fam$rqres)
    set.seed(11)
    hi <- do.call(cdf, c(list(y), par))
    u <- switch(fam$type,
                Continuous = hi,
                Mixed = ifelse(y == 0, runif(length(y), 0, par$tau), hi),
                Discrete = runif(length(y), do.call(cdf, c(list(y - 1), par)), hi))
    stopifnot(equal(pnorm(z), u))
  }
  stopifnot(!"package:egpd" %in% search())
})
