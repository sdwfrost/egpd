test_that("scalar Brent evaluates the named likelihood and returns usable curvature", {
  y <- rep(0:3, c(17,23,11,9))
  k <- uniroot(function(t) t*exp(t)-mean(y), c(.1,3), tol=1e-13)$root
  ll <- sum(y)*log(k) - length(y)*expm1(k) + sum(log(c(1,1,1,5/6)[1+y]))
  for (h in c(TRUE,FALSE)) {
    f <- distfit(y, family="bell", optim.method="Brent", lower=log(.1), upper=log(3), hessian=h)
    expect_equal(coef(f), c(theta=k), tolerance=1e-6)
    expect_equal(as.numeric(logLik(f)), ll, tolerance=1e-12)
    if (h) expect_equal(unname(f$sd), 1/sqrt(sum(y)/k^2+length(y)*exp(k)), tolerance=1e-5)
  }
  g <- distfit(y, family="degpd", fix.arg=list(xi=0,kappa=1),
               optim.method="Brent", lower=log(.1), upper=log(10))
  expect_equal(coef(g), c(sigma=1/log1p(1/mean(y))), tolerance=1e-6)
})

test_that("compound two-cell jumps identify only the effective Poisson rate", {
  y <- rep(0:5, c(22,38,25,10,4,1)); n <- length(y); m <- mean(y)
  fixed <- list(sigma=1.5,xi=-1); ll <- sum(dpois(y,m,log=TRUE))
  for (st in list(NULL,list(kappa=4,lambda=m/(1-(2/3)^4)),list(mu=.5))) {
    f <- distfit(y, family="cpdegpd", fix.arg=fixed, start=st)
    expect_equal(coef(f), c(mu=m))
    expect_equal(attr(logLik(f),"df"), 1)
    expect_equal(as.numeric(logLik(f)), ll)
    expect_equal(AIC(f), 2-2*ll)
    expect_equal(BIC(f), log(n)-2*ll)
    expect_equal(f$sd, c(mu=sqrt(m/n)))
    expect_equal(unname(confint(f)), matrix(m*exp(c(-1,1)*qnorm(.975)/sqrt(n*m)),1L))
  }
  g <- distfit(y, family="cpdegpd", fix.arg=c(fixed,list(kappa=1)))
  expect_equal(coef(g), c(lambda=3*m), tolerance=1e-5)
  expect_equal(AIC(g), 2-2*ll, tolerance=1e-9)
  expect_error(distfit(c(0,0),family="cpdegpd",fix.arg=fixed), "all aggregate counts are zero")
})

test_that("hurdle prediction inverts the fitted positive probability, not one minus its complement", {
  pos <- rep(c(9,23,31,151),each=19)
  off <- c(0, qlogis(1.625*2^-53), qlogis(2^-59))
  d <- data.frame(y=c(pos,0*pos,0,0), off=c(rep(0,2*length(pos)),off[-1]))
  f <- suppressWarnings(fit_hurdle(list(y~1,~1,~1),d, family="degpd",
    zero.formula=~offset(off),degpd.args=list(link="identity"),fixed=list(xi=0,kappa=1),trace=-1))
  expect_true(f$converged)
  nd <- data.frame(off=off)
  prob <- c(0,.75,1-2^-53,1)
  lp <- plogis(predict(f,nd,type="link")$logitppositive,log.p=TRUE)
  s <- predict(f,nd,type="response")$scale
  expected <- pmax(ceiling(s*outer(lp,log1p(-prob),"-")),0)
  got <- predict(f,nd,prob=prob)
  expect_equal(unname(as.matrix(got)), expected)
  expect_equal(as.numeric(sub("^q:", "", names(got))), prob, tolerance=0)
})
