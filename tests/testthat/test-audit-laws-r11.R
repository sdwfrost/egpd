test_that("small-delta quantiles and inverse draws obey the finite carrier CDF", {
  p <- c(1e-7, 0.1, 0.5, 0.9, 1-1e-7)
  for (delta in c(1, 1e-12, 1e-16, 1e-20, 1e-310)) for (type in 4:5) {
    k <- if (type == 4) 2 else 3
    a <- list(sigma = 10, xi = 0.1, delta = delta, kappa = k, type = type)
    ## A has the law Exp(1) + Exp(1+delta), independently. Its exact CDF
    ## below avoids a beta quantile; at subnormal delta the gamma limit is
    ## exact to working precision.
    F <- function(y) {
      A <- log1p(.01*y)/.1
      g <- if (delta < .Machine$double.xmin) pgamma(A, shape = 2) else
        -expm1(-A + log1p(-expm1(-delta*A)/delta))
      g^(k/2)
    }
    ref <- vapply(p, function(t) uniroot(function(x) F(x)-t, c(0, 1000), tol=1e-12)$root, 0)
    refd <- ceiling(ref) - 1
    expect_equal(do.call(qegpd, c(list(p), a)), ref, tolerance=1e-10)
    expect_equal(do.call(qdiscegpd, c(list(p), a)), refd)
    expect_equal(do.call(regpd, c(list(n=length(p), unifsamp=p), a)), ref, tolerance=1e-10)
    expect_equal(do.call(rdiscegpd, c(list(n=length(p), unifsamp=p), a)), refd)
    expect_equal(do.call(qzidiscegpd, c(list(.25 + .75*p, pi=.25), a)), refd)
  }
})

test_that("GPIG lower probabilities stay bounded and near-one quantiles use the survivor", {
  ## For a=1/2 this is also PIG(mu=50*sqrt(2), sigma=1/(100*sqrt(2))).
  expect_equal(pgpig(500, .5, 200, .5), 1, tolerance=0)
  expect_equal(pgpig(500, .5, 200, .5, lower.tail=FALSE) / 1.23704787587e-100,
               1, tolerance=1e-9)
  ## 120-decimal PGF coefficient recursion: S(793)=5.27954e-15,
  ## S(794)=4.85948e-15. Both bracket the complement of the binary64 target.
  p <- 1 - 5e-15
  expect_identical(qgpig(c(.9, p), .3, 200, .9, max.value=1000), c(329, 794))
  expect_identical(qzigpig(.1+.9*p, .3, 200, .9, pi=.1, max.value=1000), 794)
  expect_warning(q <- qgpig(p, .3, 200, .9, max.value=793), "not attained")
  expect_true(is.na(q))
  expect_identical(qgpig(1, .3, 200, .9), Inf)
  ## The zero atom alone resolves these; a power-law tail need not be summed.
  expect_equal(pgpig(0, .3, 1e-18, 1, lower.tail=FALSE) / -expm1(-1e-18), 1)
  expect_identical(qgpig(c(.9, p), .3, 1e-18, 1, max.value=8), c(0, 0))
  ## With c=1 the upper tail need not be summable, but F(0)=exp(-b)
  ## and F(1)=F(0)*(1+a*b) bracket this target in both mixtures.
  t <- 1 - 7.5e-10
  for (pi in c(0, .1)) {
    F0 <- pi + (1-pi) * exp(-1e-9)
    F1 <- pi + (1-pi) * exp(-1e-9) * (1 + .5e-9)
    expect_true(F0 < t && t < F1)
    q <- if (pi == 0) qgpig(t, .5, 1e-9, 1, max.value=2) else
      qzigpig(t, .5, 1e-9, 1, pi=pi, max.value=2)
    expect_identical(q, 1)
  }
})

test_that("duration polynomial identifiability counts distinct fitted scales", {
  dat <- list(`1`=1:10, `4`=2:11)
  expect_error(fit_aggregate_egpd(dat, c(1,4)), "2 distinct durations.*at most 1")
  expect_error(fit_aggregate_egpd(dat, c(1,4,4), p=1, q=2), "2 distinct durations")
  start <- list(sigma_coef=c(log(2),.2), lambda_coef=log(.7), kappa=1, xi=.2)
  fit <- fit_aggregate_egpd(dat, c(1,4), p=1, q=0, start=start, control=list(maxit=0))
  expect_equal(predict(fit, durations=2)$sigma, 2*2^.2)
})

test_that("continuous BEGPD preserves both small beta shares without creating zeros", {
  n <- 4000; th <- .11
  set.seed(294)
  y <- rbegpd(n, kappa=2, sigma=1, xi=.1, thL=th, thU=th, thw=.35)
  set.seed(294)
  u <- runif(n)
  radius <- expm1(-.1*log1p(-sqrt(u)))/.1
  expect_true(all(y > 0) && all(is.finite(y)))
  expect_equal(unname(rowSums(y)), radius, tolerance=1e-12)
  ## In either unmixed radial region each share has the same Beta(th,th) law.
  pure <- u < .35 | u > .65
  prob <- pbeta(.01, th, th)
  obs <- colMeans(y[pure, ] / radius[pure] < .01)
  expect_true(all(abs(obs-prob) < 6*sqrt(prob*(1-prob)/sum(pure))))
})
