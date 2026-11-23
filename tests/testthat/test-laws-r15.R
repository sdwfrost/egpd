test_that("large adjacent geometric cells retain scale and shape scores", {
  f <- DEGPD1(sigma.link="identity")
  for (s in c(8,2^45,2^49)) {
    y <- s*c(0,1.125,2.5,4); z <- y/s; h <- 1/s
    expect_equal(s*f$dldm(y,s,0,1), z-h/expm1(h), tolerance=1e-12)
    expect_equal(f$dldd(y,s,0,1), (z*z-(2*z+h)*(h/expm1(h)))/2, tolerance=1e-12)
    H0 <- -expm1(-z[-1]); v <- exp(-z[-1])*(-expm1(-h))/H0
    expect_equal(f$dldv(y[-1],s,0,1), log(H0)+(1+v)*(log1p(v)/v), tolerance=1e-12)
  }
})

test_that("power-beta lower-tail inverses satisfy the integer crossing", {
  p <- .1; d <- 1e16; k <- .12; scale <- 1e18
  F <- function(y) {
    u <- -expm1(-log1p(.2*y/scale)/.2)
    mass <- u*(1+1/d)*integrate(function(t) -expm1(d*log1p(-u*t)),0,1)$value
    mass^(k/2)
  }
  q <- qEGPD4(p, mu=scale, sigma=.2, nu=d, tau=k)
  j <- qDEGPD4(p, mu=scale, sigma=.2, nu=d, tau=k)
  expect_equal(F(q), p, tolerance=1e-12)
  expect_lt(F(j), p)
  expect_gte(F(j+1), p)
  expect_equal(regpd(1,sigma=scale,xi=.2,kappa=k,delta=d,type=5,unifsamp=p),q)
})

test_that("truncated-beta densities and scores cancel the common normalizer", {
  f <- EGPD6(); fd <- DEGPD6()
  a <- 15/16
  L <- function(s) log1p(-(a*s)^2)
  # Shape derivative is E[L | observed cell] - E[L] on the carrier scale.
  for (k in c(1.5,1e-14)) {
    moment <- function(lo,hi) integrate(function(s) L(s)*exp((k-1)*L(s)),lo,hi,
                                          rel.tol=1e-12)$value /
      integrate(function(s) exp((k-1)*L(s)),lo,hi,rel.tol=1e-12)$value
    surv <- function(y) (1+.2*y/3)^(-5)
    y <- c(0,2)
    expect_equal(f$dldv(y,3,.2,k), L(surv(y))-moment(0,1), tolerance=1e-9)
    expect_equal(fd$dldv(y,3,.2,k), vapply(y,function(j)
      moment(surv(j+1),surv(j))-moment(0,1),0), tolerance=1e-9)
  }
  k <- 1e16; x <- 2e8
  # Rescale z=a*sqrt(k)*s; the omitted integral beyond 14 is negligible.
  w <- function(z) exp((k-1)*log1p(-z*z/k))
  Z <- integrate(w,0,14,rel.tol=1e-12)$value
  density <- function(t) a*sqrt(k)*w(a*sqrt(k)/(1+t))/(Z*(1+t)^2)
  expect_equal(dEGPD6(x,1,1,k)/density(x), 1, tolerance=1e-9)
  expect_equal(dDEGPD6(x,1,1,k)/integrate(function(t)density(x+t),0,1)$value,
               1, tolerance=1e-9)
  E <- integrate(function(z) k*log1p(-z*z/k)*w(z),0,14,rel.tol=1e-12)$value/Z
  expect_equal(k*f$dldv(x,1,1,k), k*L(1/(1+x))-E, tolerance=1e-9)
})

test_that("pinv distinguishes an explicit absolute zero cutoff from the default", {
  A <- rbind(c(0,2^56),c(1,0),c(0,0))
  inv <- rbind(c(0,1,0),c(2^-56,0,0))
  expect_equal(pinv(A,tol=0), inv)
  expect_equal(pinv(A,tol=.5), inv)
  expect_equal(drop(pinv(A)%*%c(0,1,0)), c(0,0))
  expect_equal(pinv(diag(c(2,0)),tol=0), diag(c(.5,0)))
})
