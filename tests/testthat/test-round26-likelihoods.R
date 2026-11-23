test_that("native DEGPD fits retain a positive last cell near an integer endpoint", {
  s <- 76; xi <- -s/3
  # For this binary xi, 76+3*xi is exactly 2^-48. Multiplying xi*3
  # first would round that residual away.
  S <- c(1,(1+xi/s)^(-1/xi),(1+2*xi/s)^(-1/xi),(2^-48/s)^(-1/xi),0)
  for (n3 in c(0L,136L)) {
    counts <- c(10L,16L,439L,n3)
    ll <- function(k) sum(counts*log(diff((1-S)^k)))
    ref <- optimize(ll, c(.15,2), maximum=TRUE)
    d <- data.frame(y=rep(0:3,counts))
    f <- egpd(list(y~1,~1,~1),d,family="degpd",degpd.args=list(m=1,link="identity"),
              fixed=list(sigma=s,xi=xi),inits=c(0,0,0),restarts=FALSE,trace=-1)
    k <- predict(f,d[1,,drop=FALSE],type="response")$kappa
    expect_equal(f$convergence,0L)
    expect_equal(as.numeric(logLik(f)),ll(k),tolerance=1e-10)
    expect_equal(k,ref$maximum,tolerance=1e-4)
  }
})

test_that("compound fits withhold Wald uncertainty at a jump-support kink", {
  fit <- function(m, h) distfit(rep(0:1,c(30-m,m)),family="cpdegpd",
    fix.arg=list(lambda=.5,xi=-2,kappa=1),start=list(sigma=4),
    optim.method="Brent",lower=log(2.05),upper=log(6),control=list(ndeps=h))
  for (h in c(.01,.0001)) {
    expect_warning(f <- fit(18,h),"change in count support")
    expect_equal(unname(coef(f)),4,tolerance=1e-6)
    expect_true(all(is.na(c(f$sd,confint(f)))))
    expect_error(vcov(f),"change in count support")
  }
  g <- fit(9,.0001)
  expect_equal(unname(coef(g)),3.125,tolerance=1e-6)
  expect_equal(unname(g$sd),1.171875,tolerance=1e-5)
})

test_that("degree-one Bernstein likelihood keeps the exponential shape limit", {
  x <- rep(c(.51,.55,.62,.7,.8,.91,1.2,1.35,1.44),4)
  k <- -length(x)/sum(log(-expm1(-x)))
  for (xi in c(0,2^-1074,-2^-1074)) {
    f <- distfit(x,method="bernstein",bernstein.m=1,fix.arg=list(sigma=1,xi=xi),
                 start=list(kappa=1))
    expect_equal(unname(coef(f)),k,tolerance=1e-6)
    expect_equal(f$loglik,sum(log(k)+(k-1)*log(-expm1(-x))-x),tolerance=1e-10)
  }
})

test_that("scaled type-5 lower tails preserve the type-4 transformation", {
  s <- 2^900; k <- 1/512
  for (d in c(1,2^1023)) {
    x <- if (d==1) 2^-636 else 2^-124
    # At d=1 the carrier is u^k. In the second case d*u=1/2,
    # log G4 = log(u) + log(1 - (1-exp(-1/2))/(1/2)) to rounding.
    p <- if (d==1) 1/8 else exp((k/2)*(-1024*log(2) + log1p(2*expm1(-.5))))
    a <- list(sigma=s,xi=-1,kappa=k,delta=d,type=5)
    expect_equal(do.call(pegpd,c(list(x),a)),p,tolerance=1e-12)
    expect_equal(do.call(qegpd,c(list(log1p(-p), lower.tail=FALSE,log.p=TRUE),a))/x,
                 1,tolerance=1e-12)
    expect_equal(do.call(regpd,c(list(n=1,unifsamp=p),a))/x,1,tolerance=1e-12)
    if (d==1) expect_equal(do.call(degpd_density,c(list(x),a))/(k*p/x),1,tolerance=1e-12)
  }
})
