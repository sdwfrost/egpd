test_that("resolved count survivors are not overwritten by rounded zero masses", {
  m <- 3e-18
  expect_equal(pcmp(0:2, m, 1, lower.tail=FALSE) / ppois(0:2, m, lower.tail=FALSE),
               rep(1,3), tolerance=1e-9)
  expect_identical(qcmp(2e-20,m,1,lower.tail=FALSE), 1)
  ## Integrating the lognormal survivor over the gamma(k+1,1) arrival time.
  ref <- c(4.23934596674e-23, 2.14316771238e-23, 1.60812828174e-23, 1.33960099531e-23)
  expect_equal(ppoislnorm(0:3,10,20,lower.tail=FALSE)/ref, rep(1,4), tolerance=1e-9)
  expect_equal(ppoislnorm(0:3,10,20,log.p=TRUE)/(-ref), rep(1,4), tolerance=1e-9)
  expect_equal(qpoislnorm(1.5e-23,10,20,lower.tail=FALSE), 3)
})

test_that("Poisson-lognormal zero mass and scores survive wide log-intensity mixtures", {
  ## S(0) <= E[sqrt(Lambda)] = sqrt(mu)*exp(-sigma^2/8), so both
  ## zero-count log likelihoods round to zero even with log-scale output.
  expect_equal(dpoislnorm(0,10,c(100,400),log=TRUE), c(0,0), tolerance=0)
  expect_equal(ppoislnorm(0,10,400), 1)
  expect_equal(qpoislnorm(.5,10,400), 0)
  expect_true(is.finite(dpoislnorm(1,10,400,log=TRUE)))
  expect_equal(dpoislnorm(0,10,20,log=TRUE)/(-4.23934596674e-23), 1, tolerance=1e-9)
  X <- matrix(1,1,1); offs <- list(numeric(),numeric())
  nll <- function(eta) egpd:::plnorm1d0(as.list(eta),X,X,0,integer(),0L,offs)
  eta <- log(c(10,20))
  M <- egpd:::plnorm1d12(as.list(eta),X,X,0,integer(),0L,offs)
  fd <- apply(diag(2)*1e-6,1,function(h) (nll(eta+h)-nll(eta-h))/(2e-6))
  expect_equal(as.numeric(M[1,1:2])/fd, c(1,1), tolerance=1e-7)
  expect_equal(nll(log(c(10,400))), 0, tolerance=0)
})

test_that("compound likelihood retains subnormal and underflowed jump probabilities", {
  ## For exponential severity condition on the number of positive
  ## geometric jumps; their sum minus that number is negative binomial.
  exact <- function(k, gauge) {
    n <- seq_len(k); r <- 1e-6 * if (gauge) 1 else exp(-1)
    a <- dpois(n,r,log=TRUE) + dnbinom(k-n,n,-expm1(-1),log=TRUE)
    max(a)+log(sum(exp(a-max(a))))
  }
  for (k in c(700,800)) {
    expect_equal(dcpdegpd(k,1e-6,kappa=1,sigma=1,xi=0,log=TRUE),exact(k,FALSE))
    expect_equal(dcpegpd(k,1e-6,kappa=1,sigma=1,xi=0,h=1,log=TRUE),exact(k,TRUE))
  }
  ## Only the 1-jump term contributes at S=1. The individual jump itself
  ## is below the linear range, not outside the distribution's support.
  lh <- 500*log(-expm1(-log1p(.2*(1:2)/10)/.2))
  g1 <- -2*(-expm1(lh[1]))+log(2)+lh[2]+log(-expm1(lh[1]-lh[2]))
  expect_equal(dcpdegpd(1,2,kappa=500,sigma=10,xi=.2,log=TRUE),g1)
  ## Underflow in P(S>0) must not be interpreted as impossible positive
  ## aggregates: conditional on positivity the first bin here has mass 1
  ## to working precision.
  expect_equal(egpd:::.aggregate_positive_logdens(1:2,"cpdegpd",.001,1,0,1),
               c(0,-1000+log(1.5)), tolerance=1e-12)
})

test_that("only gauge CDFs snap near lattice points", {
  q <- c(.7,1-3e-9,1,2-4e-9,2)
  ## Uniform(0,2) floored jumps thin the Poisson rate by one half.
  expect_equal(pcpdegpd(q,1.6,sigma=2,xi=-1,kappa=1),ppois(floor(q),.8))
  expect_equal(pcpegpd(q[2:3],1.6,sigma=1,xi=-1,kappa=1,h=1),rep(ppois(1,1.6),2))
})

test_that("blended log densities preserve both the beta bulk and pure GP tails", {
  a <- qbeta(.8,2,5); b <- qbeta(.99,2,5)
  s <- (b-a)/log(.2/.01)
  expect_equal(dbgpd(c(50,80),xi=0,log=TRUE),log(.2)-log(s)-(c(50,80)-a)/s)
  x <- c(1e-18,1e-100)
  expect_equal(dbgpd(x,xi=.2,bulk.args=list(shape1=10,shape2=5),log=TRUE),dbeta(x,10,5,log=TRUE))
  mid <- seq(a,b,length.out=7)
  expect_equal(dbgpd(mid,xi=.2,log=TRUE),log(dbgpd(mid,xi=.2)),tolerance=1e-12)
})
