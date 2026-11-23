test_that("GPIG scores follow the exact PGF down to tiny means and near-unit c", {
  a <- .37
  c <- c(.44, .99989, 1 - 2^-52)
  mu <- c(1.1e-8, 1.2, 1.2)
  K <- ((1-c)^(1-a) - (1-c))/(a*c)
  dK <- (1 - (1-a*c)*(1-c)^(-a))/(a*c^2)
  ## [z^1]G(z) = exp(-K*mu)*mu*(1-c)^(1-a).
  for (fam in list(GPIG(), ZIGPIG())) {
    p <- list(y = 1, mu = mu, sigma = a, nu = c)
    if (fam$nopar == 4L) p$tau <- .2
    sm <- do.call(fam$dldm, p)
    sc <- do.call(fam$dldv, p)
    expect_equal(mu*sm, 1 - K*mu, tolerance = 1e-12)
    expect_equal(c*(1-c)*sc, -mu*c*(1-c)*dK - (1-a)*c, tolerance = 1e-12)
    expect_equal(do.call(fam$d2ldmdv, p), -sm*sc)
  }
  ## The zero's structural-component score, including the right derivative at pi=0.
  m <- c(1e-16, 3, 5)
  pi <- c(.7, 0, .3)
  lp0 <- -K[1]*m
  expect_equal(ZIGPIG()$dldt(0, m, a, c[1], pi) / m,
               -expm1(lp0) / (pi + (1-pi)*exp(lp0)) / m, tolerance = 1e-13)
})

test_that("positive-count GPIG scores do not switch to the fitting approximation", {
  old <- gpig_control(method = "saddlepoint", order = 1)
  on.exit(do.call(gpig_control, old))
  f <- GPIG()
  y <- c(0, 1, 5, 107)
  mu <- 2.4; a <- .5; c <- .83
  eta <- c(log(mu), qlogis(a), qlogis(c))
  jac <- c(mu, a*(1-a), c*(1-c))
  ll <- function(e) dGPIG(y, exp(e[1]), plogis(e[2]), plogis(e[3]), log = TRUE)
  for (i in 1:3) {
    h <- (1:3 == i)*1e-5
    score <- f[[c("dldm", "dldd", "dldv")[i]]](y, mu, a, c)
    expect_equal(score*jac[i], (ll(eta+h) - ll(eta-h))/(2e-5), tolerance = 1e-8)
  }
  expect_identical(f$dldm(NA_real_, 1, a, c), NA_real_)
  expect_warning(bad <- f$dldm(-1, 1, a, c), "undefined")
  expect_true(is.nan(bad))
})

test_that("BAMLSS count residuals use the support, not the held-out minimum", {
  skip_if_not_installed("bamlss")
  skip_on_cran()
  if (!"package:bamlss" %in% search()) {
    suppressPackageStartupMessages(library(bamlss))
    on.exit(detach("package:bamlss", unload = FALSE, character.only = TRUE))
  }
  d <- data.frame(y = rep(c(0, 0, 0, 0, 1, 1, 2, 7), 12))
  fit <- suppressWarnings(bamlss::bamlss(y ~ 1, family = egpd::bell_bamlss(), data = d,
                                          sampler = FALSE, verbose = FALSE))
  mu <- as.numeric(predict(fit, type = "parameter"))[1]
  th <- uniroot(function(t) t*exp(t) - mu, c(0, 10), tol = 1e-14)$root
  p0 <- exp(-expm1(th))
  nd <- data.frame(y = rep(1, 8))
  set.seed(11); r <- as.numeric(residuals(fit, newdata = nd, type = "quantile"))
  set.seed(11); expect_equal(pnorm(r), runif(8, p0, p0*(1+th)), tolerance = 1e-12)
  expect_equal(as.numeric(residuals(fit, type = "response")), d$y - mu)
})

test_that("mixed and count BAMLSS PITs keep row parameters and log tails", {
  p <- list(sigma = c(1, 2, 1, 1), xi = 0, kappa = 1, pi = c(.4, .1, .2, .3))
  y <- c(0, 2, 0, 1000)
  set.seed(37); z <- ziegpd_bamlss()$rqres(y, p)
  set.seed(37); expect_equal(pnorm(z[c(1, 3)]), runif(2)*p$pi[c(1, 3)])
  expect_equal(z[c(2, 4)], qnorm(log1p(-p$pi[c(2, 4)]) - c(1, 1000),
                                  log.p = TRUE, lower.tail = FALSE))
  set.seed(42); z <- zidegpd_bamlss()$rqres(y, p)
  set.seed(42); w <- runif(4)
  expect_equal(z[4], qnorm(log1p(-p$pi[4]) - 1000 + log1p(w[4]*expm1(-1)),
                           log.p = TRUE, lower.tail = FALSE))
  expect_equal(pnorm(z[c(1, 3)]), w[c(1, 3)] *
                 (1 - (1-p$pi[c(1, 3)])*exp(-1/p$sigma[c(1, 3)])))
  for (ctor in list(pig_bamlss, zipig_bamlss, gpig_bamlss, zigpig_bamlss,
                    bell_bamlss, zibell_bamlss)) {
    f <- ctor()
    par <- list(mu = c(1, 1.5, 2), sigma = .37, nu = .6, pi = c(.7, .3, .2))[f$names]
    target <- if (is.null(par$pi)) par$mu else (1-par$pi)*par$mu
    expect_equal(f$mu(par), target)
    set.seed(20); z <- f$rqres(1:3, par)
    set.seed(20); expect_equal(pnorm(z), runif(3, f$p(0:2, par), f$p(1:3, par)))
  }
  ## residuals.bamlss(FUN = identity) supplies one column per posterior draw.
  m <- matrix(c(.4, 1, 2, .8), 2)
  p <- list(sigma = m, xi = 0, kappa = 1)
  set.seed(2); z <- degpd_bamlss()$rqres(c(1, 2), p)
  set.seed(2); w <- matrix(runif(4), 2)
  expect_equal(pnorm(z), 1 - exp(-c(1, 2)/m)*(1 + w*expm1(-1/m)))
  expect_equal(degpd_bamlss()$mu(p), 1/expm1(1/m))
})

test_that("BAMLSS response means integrate EGPD and sum floored EGPD laws", {
  ## A finite endpoint permits independent finite moment integrals and PMF sums.
  p <- list(sigma = 1.5, xi = -.25, kappa = 1.4, delta = 2.7, prob = .35, pi = .4)
  for (m in 1:6) {
    f <- egpd_bamlss(m); fd <- degpd_bamlss(m)
    par <- p[f$names]
    mc <- integrate(function(x) x*f$d(x, par), 0, -p$sigma/p$xi, rel.tol = 1e-9)$value
    y <- 0:ceiling(-p$sigma/p$xi)
    md <- sum(y*fd$d(y, par))
    expect_equal(f$mu(par), mc, tolerance = 1e-8)
    expect_equal(fd$mu(par), md, tolerance = 1e-8)
    expect_equal(ziegpd_bamlss(m)$mu(c(par, p["pi"])), (1-p$pi)*mc, tolerance = 1e-8)
    expect_equal(zidegpd_bamlss(m)$mu(c(par, p["pi"])), (1-p$pi)*md, tolerance = 1e-8)
  }
  ## Geometric, and S_X(n) = (1+n)^-2 with an infinite sum known in closed form.
  f <- degpd_bamlss()
  expect_equal(f$mu(list(sigma = c(3, .5), kappa = 1, xi = c(0, .5))),
               c(1/expm1(1/3), pi^2/6-1), tolerance = 1e-7)
  expect_warning(m <- ziegpd_bamlss()$mu(list(sigma = 1, kappa = 1, xi = 1, pi = c(.5, 1))),
                 "infinite")
  expect_equal(m, c(Inf, 0))
})

test_that("TF2 predictive coverage supports the default shifted df link", {
  skip_if_not_installed("gamlss")
  set.seed(1432)
  d <- data.frame(y = 2.3 + 1.5*rt(120, 4), x = 1:120)
  f <- gamlss::gamlss(y ~ 1, family = gamlss.dist::TF2(), data = d, trace = FALSE)
  n <- 500L
  nd <- data.frame(x = 50)
  set.seed(95)
  ans <- predictive_coverage(f, 2, newdata = nd, data_fit = d, nsim = n, method = "parametric")
  ## Exact Hessian of the variance-standardised t log likelihood in
  ## (location, log SD, log(df-2)); upstream vcov() is itself numerical.
  m <- coef(f, "mu")[1]; s <- exp(coef(f, "sigma")[1]); k <- exp(coef(f, "nu")[1])
  v <- 2+k; r <- d$y-m; a <- k*s*s; den <- a+r*r
  e <- r*r/den; b <- a/den; t <- r/den
  H <- matrix(0, 3, 3)
  H[1,1] <- sum((v+1)*(r*r-a)/den^2)
  H[1,2] <- H[2,1] <- sum(-2*(v+1)*t*b)
  H[1,3] <- H[3,1] <- sum(k*t-(v+1)*t*b)
  H[2,2] <- sum(-2*(v+1)*e*b)
  H[2,3] <- H[3,2] <- sum(k*e-(v+1)*e*b)
  H[3,3] <- sum(k/2*(digamma((v+1)/2)-digamma(v/2)) +
                k*k/4*(trigamma((v+1)/2)-trigamma(v/2)) -
                k/2*log1p(r*r/a)+k*e-(v+1)*e*b/2)
  se <- sqrt(diag(solve(-H)))
  expect_equal(unname(sqrt(diag(.gamlss_vcov(f)))), se, tolerance = 2e-8)
  set.seed(95)
  loc <- rnorm(n, coef(f, "mu"), se[1])
  sig <- exp(rnorm(n, coef(f, "sigma"), se[2]))
  df <- 2 + exp(rnorm(n, coef(f, "nu"), se[3]))
  q <- quantile(loc + sig*sqrt((df-2)/df)*qt(runif(n), df), c(.025, .975), names = FALSE)
  expect_lt(max(abs(c(ans$L, ans$U) - q)), 1e-8)
  ## Plug-in requires only the fitted response-scale parameters.
  f$nu.link <- NULL
  set.seed(19)
  ans <- predictive_coverage(f, d$y, nsim = n, method = "plug-in")
  set.seed(19)
  q <- quantile(f$mu.fv[1] + f$sigma.fv[1]*sqrt(1-2/f$nu.fv[1])*qt(runif(n), f$nu.fv[1]),
                c(.025, .975), names = FALSE)
  expect_equal(c(ans$L[1], ans$U[1]), q)
})
