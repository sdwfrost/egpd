test_that("a bounded support must identify all free masses or positive-jump rates", {
  y <- rep(0:2, c(29,53,43))
  fix <- list(sigma=2.5, xi=-1)
  st <- list(kappa=1, delta=3, prob=.5)
  expect_error(distfit(y, type=6, family="degpd", fix.arg=fix, start=st),
               "3 cells: at most 2")
  expect_error(distfit(rep(0:5,12), type=6, family="cpdegpd", fix.arg=fix),
               "jump support has 3 cells")
  fit <- distfit(y, type=6, family="degpd", fix.arg=c(fix,list(prob=.5)),
                 start=st[c("kappa","delta")], optim.method="BFGS")
  expect_equal(coef(fit), c(kappa=1,delta=3), tolerance=1e-5)
  expect_equal(as.numeric(logLik(fit)), sum(c(29,53,43)*log(c(29,53,43)/125)))
  expect_true(all(is.finite(vcov(fit))))
  # With free scale the same local support is not a reduction of the whole family.
  expect_warning(local_fit <- distfit(y, type=6, family="degpd", fix.arg=list(xi=-1,prob=.5),
                 start=list(sigma=2.5,kappa=1,delta=3), optim.method="BFGS"),
                 "Fitted support has 3 cell")
  expect_true(isTRUE(local_fit$support.ridge))
  expect_true(all(is.na(confint(local_fit))))
})

test_that("equivalent continuous carriers agree on the excluded scale maximum", {
  x <- c(.4,.8,1.2,1.4,1.7,2,2.1,2.5,2.8,3)
  fixed <- list(list(kappa=1),list(kappa=0),list(kappa=1),list(delta=1),
                list(delta=1,kappa=1))
  for (t in seq_along(fixed)) {
    fix <- c(fixed[[t]],list(xi=-1))
    expect_error(distfit(x,type=t,fix.arg=fix),"Nonregular continuous endpoint")
  }
  f <- distfit(x,type=3,fix.arg=list(xi=-1,kappa=1),optim.method="Brent",
               lower=log(3.2),upper=log(8),hessian=FALSE)
  expect_equal(unname(coef(f)),3.2)
  g <- distfit(x,type=3,fix.arg=list(xi=0,kappa=1))
  expect_equal(unname(g$sd),mean(x)/sqrt(length(x)),tolerance=1e-4)
  # Not every type-5 density is monotone in scale at xi=-1.
  z <- rep(1,10)
  h <- distfit(z,type=5,fix.arg=list(xi=-1,delta=3,kappa=.2),optim.method="Brent",
               lower=0,upper=log(10))
  s <- unname(coef(h))
  ll <- function(s) sum(degpd_density(z,sigma=s,xi=-1,delta=3,kappa=.2,type=5,log=TRUE))
  expect_true(s>1 && s<10 && all(ll(s)>c(ll(.999*s),ll(1.001*s))))
})

test_that("free zero inflation attains the likelihood maximum including its natural face", {
  y <- rep(0:3,c(1,12,15,12)); n <- length(y)
  for (fam in c("zibell","zidegpd","zipig","zigpig")) {
    fixed <- switch(fam,zibell=list(theta=1),zidegpd=list(sigma=1.6,xi=0,kappa=1),
                    zipig=list(mu=2.2,sigma=.4),zigpig=list(a=1,b=2.2,c=1))
    b <- switch(fam,zibell=exp(1-exp(1)),zidegpd=-expm1(-1/1.6),
                zipig=exp((1-sqrt(1+2*.4*2.2))/.4),zigpig=exp(-2.2))
    fit <- function(x,...) distfit(x,family=fam,fix.arg=fixed,...)
    expect_warning(f <- fit(y),"natural boundary")
    expect_equal(coef(f),c(pi=0))
    expect_true(all(is.na(c(f$sd,confint(f)))))
    expect_error(vcov(f),"natural boundary")
    expect_warning(g <- fit(y,lower=qlogis(.1)), "optimizer bound")
    expect_equal(coef(g),c(pi=.1))
    h <- fit(rep(0:3,c(22,5,7,6)))
    expect_equal(coef(h), c(pi=(.55-b)/(1-b)))
    expect_equal(unname(h$sd),sqrt(.55*.45/n)/(1-b),tolerance=2e-5)
  }
  # Profile a nuisance mean, not just pi with a fixed base.
  expect_warning(f <- distfit(y,family="zigpig",fix.arg=list(a=1,c=.5)), "natural boundary")
  expect_equal(unname(coef(f)),c(2*mean(y),0),tolerance=1e-5)
  expect_equal(as.numeric(logLik(f)),sum(dpois(y,mean(y),log=TRUE)),tolerance=1e-8)
  expect_error(distfit(rep(0,20),family="zibell",fix.arg=list(theta=1)),"pi = 1")
})
