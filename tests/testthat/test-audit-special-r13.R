test_that("a splice requires independent scale columns on the actual tail rows", {
  skip_if_not_installed("qgam")
  mid <- function(n) (seq_len(n)-.5)/n
  d <- data.frame(z=rep(0:1,c(160,80)), y=c(5+qgeom(mid(160),.23),rep(0:3,each=20)))
  fit <- function(dat) {
    set.seed(777)
    fit_spliced_degpd(y~1, dat, bulk_probs=c(.25,.58), threshold_prob=.58,
                     tail_formula=~z, egpd_control=list(fix.arg=list(xi=0,kappa=1)))
  }
  expect_error(fit(d), "rank deficient.*complete exceedances.*z")
  d$y[d$z==1] <- 6+qgeom(mid(80),.16)
  f <- fit(d)
  nd <- data.frame(z=0:1)
  p <- predict(f,nd,type="parameter")
  excess <- d$y - p$threshold[1]
  mean_ex <- vapply(0:1,function(z) mean(excess[d$z==z & excess>=0]),0)
  ref <- p$threshold + qgeom((.95-p$alpha_tilde)/(1-p$alpha_tilde),1/(1+mean_ex))
  expect_equal(predict(f,nd,prob=.95)[[1]],ref)
})

test_that("automatic BFGS starts separate the two type-6 powers", {
  mid <- function(n) (seq_len(n)-.5)/n
  check <- function(data, family, fixed, loglik, gen) {
    f <- distfit(data, type=6, family=family, fix.arg=fixed, optim.method="BFGS")
    b <- as.list(coef(f))
    expect_equal(f$convergence,0L)
    expect_equal(as.numeric(logLik(f)),do.call(loglik,b),tolerance=1e-8)
    expect_gte(as.numeric(logLik(f)),do.call(loglik,gen)-1e-3)
  }
  t <- c(qexp(mid(64),rate=.6),qexp(mid(116),rate=8))
  s <- 2.75
  x <- -s*log(-expm1(-t))
  ll <- function(kappa,delta,prob)
    sum(log(prob*dexp(t,kappa)+(1-prob)*dexp(t,delta))+t-x/s-log(s))
  check(x,"egpd",list(sigma=s,xi=0),ll,list(.6,8,64/180))
  k <- 1/mean(t)
  f <- distfit(x,type=6,fix.arg=list(sigma=s,xi=0),optim.method="BFGS",
               start=list(kappa=k,delta=k,prob=.5),hessian=FALSE)
  expect_gte(as.numeric(logLik(f)),ll(.6,8,64/180)-1e-3)

  t <- c(qexp(mid(54),.45),qexp(mid(106),6))
  s <- 12; xi <- .15
  y <- floor(s*expm1(-xi*log1p(-exp(-t)))/xi)
  H <- function(k) -expm1(-log1p(xi*k/s)/xi)
  cell <- function(k) H(y+1)^k-H(y)^k
  ll <- function(kappa,delta,prob) sum(log(prob*cell(kappa)+(1-prob)*cell(delta)))
  check(y,"degpd",list(sigma=s,xi=xi),ll,list(.45,6,54/160))
})

test_that("distfit reports identifiable subfamily coefficients and information criteria", {
  y <- rep(0:5,c(11,24,17,13,16,10)); m <- mean(y); n <- length(y)
  exact <- sum(dpois(y,m,log=TRUE))
  for (st in list(list(b=m/.28,c=.28),list(mu=m))) {
    f <- distfit(y,family="gpig",fix.arg=list(a=1),start=st)
    expect_equal(coef(f),c(mu=m))
    expect_equal(unname(f$sd),sqrt(m/n),tolerance=1e-4)
    expect_equal(attr(logLik(f),"df"),1)
    expect_equal(AIC(f),2-2*exact)
    expect_equal(BIC(f),log(n)-2*exact)
    expect_equal(rownames(confint(f)),"mu")
  }
  expect_error(distfit(y,family="gpig",fix.arg=list(a=1),start=list(mu=2,c=.4)),
               "mean start.*or its factors")
  f1 <- distfit(y,family="degpd",fix.arg=list(xi=0,kappa=1))
  for (w in c(.2,.8)) {
    f <- distfit(y,family="degpd",type=6,fix.arg=list(xi=0,kappa=1,delta=1),
                 start=list(prob=w))
    expect_named(coef(f),"sigma")
    expect_equal(AIC(f),AIC(f1),tolerance=1e-7)
    expect_equal(unname(f$sd),unname(f1$sd),tolerance=1e-4)
  }
  one <- distfit(y,family="degpd",type=6,fix.arg=list(xi=0,prob=0,delta=1))
  expect_named(coef(one),"sigma")
  expect_equal(AIC(one),AIC(f1),tolerance=1e-7)
  p0 <- distfit(rep(0,12),family="zidegpd",fix.arg=list(pi=1))
  expect_length(coef(p0),0)
  expect_equal(dim(vcov(p0)),c(0L,0L))
  expect_equal(AIC(p0),0)
})
