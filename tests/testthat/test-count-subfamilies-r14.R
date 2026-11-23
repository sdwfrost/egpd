test_that("a two-cell zero-inflated power carrier fits just the Bernoulli probability", {
  y <- rep(0:1,c(68,32))
  fix <- list(sigma=.8,xi=-.5)
  for (st in list(list(pi=.2,kappa=3),list(pi=.5,kappa=7),list(p1=.75))) {
    f <- distfit(y, family="zidegpd", fix.arg=fix, start=st)
    expect_equal(coef(f), c(p1=.32))
    expect_equal(attr(logLik(f),"df"), 1)
    expect_equal(AIC(f), 2-2*sum(dbinom(y,1,.32,log=TRUE)))
    expect_equal(BIC(f), log(100)-2*sum(dbinom(y,1,.32,log=TRUE)))
    expect_equal(f$sd, c(p1=sqrt(.32*.68/100)))
    expect_equal(unname(confint(f)), matrix(plogis(qlogis(.32)+c(-1,1)*qnorm(.975)/
                   sqrt(100*.32*.68)),1L))
  }
  expect_error(distfit(c(y,2),family="zidegpd",fix.arg=fix), "outside.*two-cell support")
  expect_error(distfit(rep(0,20),family="zidegpd",fix.arg=fix), "no interior|No interior")
  one <- distfit(rep(0,20),family="zidegpd",fix.arg=list(sigma=.4,xi=-.5))
  expect_length(coef(one),0)
  expect_equal(AIC(one),0)
  expect_error(distfit(y,type=6,family="zidegpd",fix.arg=c(fix,list(prob=0))),
               "Only P\\(Y = 1\\) is identifiable")
  f <- distfit(rep(0:2,c(55,30,15)),family="zidegpd",fix.arg=list(sigma=2.5,xi=-1))
  expect_named(coef(f),c("kappa","pi"),ignore.order=TRUE)
  expect_equal(attr(logLik(f),"df"),2)
})

test_that("large exponential scales retain the hurdle's positive normalization", {
  d <- data.frame(y=rep(c(0,2^28,2^29,2^30),c(18,11,17,13)))
  y <- d$y[d$y>0]
  for(s in 2^c(38,56)) {
    f <- fit_hurdle(list(y~1,~1,~1),d, family="degpd", degpd.args=list(link="identity"),
                    fixed=list(shape=0,scale=s),trace=-1)
    k <- predict(f,newdata=d[1,,drop=FALSE],type="response")$kappa
    h <- -expm1(-y/s); h1 <- -expm1(-1/s)
    delta <- log1p(exp(-y/s)*h1/h)
    ll <- function(k) sum(k*log(h)+log(expm1(k*delta))-log(-expm1(k*log(h1)))) +
      sum(dbinom(as.integer(d$y>0),1,mean(d$y>0),log=TRUE))
    ref <- optimize(ll, c(1e-5,1), maximum=TRUE, tol=1e-12)
    expect_true(f$converged)
    expect_equal(as.numeric(logLik(f)),ll(k),tolerance=1e-10)
    expect_equal(k,ref$maximum,tolerance=1e-5)
  }
})
