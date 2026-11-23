test_that("the GPIG mean conversion uses the final observation and draw rows", {
  m <- c(.8,5.6); nu <- c(.25,.55,.72)
  x <- c(0,2,4,0,1,3)
  prob <- c(.19,.38,.53,.69,.78,.85)
  mu <- rep_len(m,6); c <- rep_len(nu,6)
  ## a=1/2 is PIG in mean/dispersion coordinates.
  dispersion <- c/(2*mu*(1-c))
  pi <- c(.15,.25)
  expect_equal(dGPIG(x,m,.5,nu,log=TRUE),dpig(x,mu,dispersion,log=TRUE))
  expect_equal(pGPIG(x,m,.5,nu),ppig(x,mu,dispersion))
  expect_equal(qGPIG(prob,m,.5,nu),qpig(prob,mu,dispersion))
  expect_equal(dZIGPIG(x,m,.5,nu,pi),dzipig(x,mu,dispersion,pi))
  expect_equal(pZIGPIG(x,m,.5,nu,pi),pzipig(x,mu,dispersion,pi))
  expect_equal(qZIGPIG(prob,m,.5,nu,pi),qzipig(prob,mu,dispersion,pi))
  ## The zero inflation alone can determine the output length.
  pis <- seq(.1,.6,length.out=6)
  expect_equal(pZIGPIG(0,m,.5,nu,pis),pzipig(0,mu,dispersion,pis))
  for (zi in c(FALSE,TRUE)) {
    rng <- if(zi) rZIGPIG else rGPIG
    r <- function(m,nu) do.call(rng,c(list(n=60,mu=m,sigma=.5,nu=nu),if(zi) list(tau=pi)))
    set.seed(213); y <- r(m,nu)
    set.seed(213); ref <- r(rep_len(m,60),rep_len(nu,60))
    expect_identical(y,ref)
  }
  expect_length(rGPIG(numeric(),m,.5,nu),0)
})

test_that("a blended GP matches the GP at both knots even at small nonzero shape", {
  a <- .65; b <- .94; L <- log1p(-a)-log1p(-b)
  args <- list(bulk.args=list(shape1=1,shape2=1),pa=a,pb=b)
  for (xi in c(-6e-9,6e-9,1.01e-8,0)) {
    s <- if (xi==0) (b-a)/L else (b-a)*xi/expm1(xi*L)
    x <- 2400
    ref <- log1p(-a)-log(s)-if (xi==0) (x-a)/s else (1+1/xi)*log1p(xi*(x-a)/s)
    expect_equal(do.call(dbgpd,c(list(x=x,xi=xi,log=TRUE),args)),ref,tolerance=1e-12)
    if (xi<0) {
      ep <- a-s/xi; at <- ep-.15
      expect_equal(do.call(qbgpd,c(list(p=1,xi=xi),args)),ep,tolerance=1e-14)
      lp <- log1p(-a)-log(s)-(1+1/xi)*log((ep-at)/(ep-a))
      expect_equal(do.call(dbgpd,c(list(x=at,xi=xi,log=TRUE),args)),lp,tolerance=1e-8)
    }
    p <- c(b,.9999999)
    xq <- a+if(xi==0) -s*(log1p(-p)-log1p(-a)) else s*expm1(-xi*(log1p(-p)-log1p(-a)))/xi
    expect_equal(do.call(qbgpd,c(list(p=p,xi=xi),args)),xq,tolerance=1e-12)
    expect_equal(do.call(pbgpd,c(list(q=xq,xi=xi),args)),p,tolerance=1e-12)
  }
})

test_that("compound type 2 accepts the uniform carrier in both severity conventions", {
  a <- list(lambda=1.4,sigma=2,xi=-1,kappa=0,type=2)
  k <- 0:4; p <- c(.2,.6,.8)
  ## Floor(Uniform(0,2)) thins the event rate; the h=2 gauge jump is always 2.
  expect_equal(do.call(dcpdegpd,c(list(k),a)),dpois(k,.7))
  expect_equal(do.call(pcpdegpd,c(list(k),a)),ppois(k,.7))
  expect_equal(do.call(qcpdegpd,c(list(p,K=20),a)),qpois(p,.7))
  h <- c(a,list(h=2))
  expect_equal(do.call(dcpegpd,c(list(2*k),h)),dpois(k,1.4))
  expect_equal(do.call(pcpegpd,c(list(2*k),h)),ppois(k,1.4))
  expect_equal(do.call(qcpegpd,c(list(p,K=20),h)),2*qpois(p,1.4))
  for (gauge in c(FALSE,TRUE)) {
    w <- if(gauge) h else a; rng <- if(gauge) rcpegpd else rcpdegpd
    set.seed(63); r <- do.call(rng,c(list(100),w))
    w$type <- 1; w$kappa <- 1
    set.seed(63); expect_equal(r,do.call(rng,c(list(100),w)))
  }
  expect_error(dcpdegpd(1,1,sigma=1,xi=0,kappa=-.2,type=2),"non-negative")
  expect_error(dcpdegpd(1,1,sigma=1,xi=0,kappa=0,type=1),"positive")
})
