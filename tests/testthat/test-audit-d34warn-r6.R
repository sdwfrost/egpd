## Round-6 (d34warn): the analytic DEGPD third/fourth derivatives behind the REML
## gradient lose accuracy in cells [y, y + 1] too narrow for the generated forms, i.e.
## at very large counts or scales, where (sigma + xi*y)/|xi| exceeds 1e-5/eps
## (about 4.5e10). Such fits now fall back to finite-difference smoothing-parameter
## selection with a specific warning of class "egpd_d34_fallback"; ordinary fits keep
## the analytic derivatives and do not warn.

d6_ld <- function(eta, y) {
  list(X = rep(list(matrix(1, length(y), 1)), length(eta)), y = cbind(y),
       idpars = seq_along(eta), dupid = 0L, duplicate = 0L,
       offsets = lapply(eta, rep, length(y)), censored = FALSE, xi.max = Inf)
}
d6_fns <- function(link) get(if (link == "log") ".degpd1fns" else ".degpd1idfns",
                             envir = asNamespace("egpd"))
d6_sim <- function(n, lsig0, xi, kappa, seed = 3) {
  set.seed(seed)
  x <- runif(n)
  sig <- exp(lsig0 + sin(2 * pi * x))
  u <- runif(n)^(1 / kappa)
  data.frame(y = floor(sig / xi * ((1 - u)^(-xi) - 1)), x = x)
}

test_that("narrow DEGPD cells (very large counts or scales) request finite differences", {
  limit <- 1e-5 / .Machine$double.eps                   # (sigma + xi*y)/|xi| = 4.5e10
  d34 <- function(link, sigma, xi, y)
    d6_fns(link)$d340(c(0, 0, 0), d6_ld(c(log(sigma), if (link == "log") log(xi) else xi,
                                         log(1.5)), y))
  for (link in c("log", "identity")) {
    ## ordinary counts and scales
    expect_true(all(is.finite(d34(link, 1e8, 0.2, c(0, 3, 1e8)))))
    ## the scale: (sigma + 0.2)/0.2 just below and above the limit (y = 1)
    expect_true(all(is.finite(d34(link, 0.98 * 0.2 * limit, 0.2, 1))))
    expect_error(d34(link, 1.02 * 0.2 * limit, 0.2, 1), class = "egpd_d34_unavailable",
                 regexp = "very large counts or scales")
    ## the count: (1e3 + 0.5 y)/0.5 below and above the limit
    expect_true(all(is.finite(d34(link, 1e3, 0.5, floor(0.98 * limit)))))
    expect_error(d34(link, 1e3, 0.5, floor(1.02 * limit)), class = "egpd_d34_unavailable",
                 regexp = sprintf("%s shape link", link))
  }
  ## negative shape: the distance to the upper endpoint sigma/|xi| counts
  expect_error(d34("identity", 1e11, -0.5, 1), class = "egpd_d34_unavailable")
  expect_true(all(is.finite(d34("identity", 1e11, -0.5, 2e11 - 1e3))))
})

test_that("DEGPD fits at very large counts fall back to finite-difference smoothing", {
  ## counts of 1e10 to 3e13 (sigma about 5e11): the analytic REML gradient was off by
  ## up to 5e-2 at such scales (sp 5% off), more than the finite-difference one
  d <- d6_sim(120, 27, 0.5, 3)
  ## REML Hessian at the fitted smoothing parameter: differences of the REML gradient
  ## with third derivatives from Richardson differences of the exact d120
  ref_H <- c(identity = 1.209434661, log = 1.24135387)
  for (link in c("identity", "log")) {
    args <- list(list(y ~ s(x, k = 6), ~ 1, ~ 1), data = d, family = "degpd",
                 restarts = FALSE, degpd.args = list(m = 1, link = link))
    expect_warning(fit <- do.call(egpd, args), class = "egpd_d34_fallback",
                   regexp = "very large counts or scales.*finite-difference outer")
    fd <- do.call(egpd, c(args, list(outer = "fd")))
    expect_identical(fit$sp, fd$sp)
    expect_identical(fit$coefficients, fd$coefficients)
    expect_identical(fit$Vc, fd$Vc)
    ## the REML Hessian behind Vc (correctV) no longer comes from differences of the
    ## analytic gradient (5.24 and 1.92 here) but from second differences of the criterion
    rho <- log(fd$sp / fd$likdata$gamma)
    expect_error(egpd:::.reml12(rho, fd$likfns, fd$likdata, fd$Sdata, beta = fd$coefficients),
                 class = "egpd_d34_unavailable")
    expect_equal(c(egpd:::.reml2.fdfd(rho, fd$likfns, fd$likdata, fd$Sdata, beta = fd$coefficients)),
                 ref_H[[link]], tolerance = 1e-4)
  }
  ## an ordinary smooth DEGPD fit keeps the analytic derivatives and does not warn
  set.seed(31)
  d <- data.frame(x = seq(0, 1, length.out = 240))
  d$y <- rDEGPD1(240, mu = exp(1 + 1.5 * sin(2 * pi * d$x)), sigma = 0.3, nu = 1.5)
  for (link in c("log", "identity"))
    expect_no_warning(egpd(list(y ~ s(x, k = 7), ~ 1, ~ 1), data = d, family = "degpd",
                           restarts = FALSE, degpd.args = list(m = 1, link = link)))
})
