test_that("a Bernstein fit retains the lower tail before scaling", {
  x <- 2^(-230:-210)
  for (s in 2^c(400,900)) {
    lh <- log(x)-log(s)
    k <- -length(x)/sum(lh)
    f <- distfit(x, method="bernstein", bernstein.m=1,
                 fix.arg=list(sigma=s,xi=-1), start=list(kappa=k))
    expect_equal(unname(coef(f))/k, 1, tolerance=1e-7)
    expect_equal(as.numeric(logLik(f)), sum(log(k)+(k-1)*lh-log(s)), tolerance=1e-12)
    d <- k*exp(k*lh)/x
    expect_equal(egpd:::.bernstein_full_density(x,s,-1,k,1,1)/d, rep(1,length(x)),
                 tolerance=1e-12)
    p <- egpd:::.bernstein_full_cdf(x,s,-1,k,1,1)
    expect_equal(p,exp(k*lh),tolerance=1e-12)
    expect_equal(egpd:::.bernstein_full_quantile(p,s,-1,k,1,1)/x,rep(1,length(x)),
                 tolerance=1e-10)
  }
})

test_that("BAMLSS coverage uses the offset in the fitted likelihood", {
  skip_if_not_installed("bamlss")
  skip_on_cran()
  if (!"package:bamlss" %in% search()) {
    suppressPackageStartupMessages(library(bamlss))
    on.exit(detach("package:bamlss", unload=FALSE, character.only=TRUE), add=TRUE)
  }
  set.seed(2751); n <- 120
  d <- data.frame(x=runif(n,.25,2))
  d$y <- floor(exp(1+.4*d$x+log(d$x))/.3 * expm1(-.3*log1p(-runif(n)^.8)))
  fit <- function(ff,...) suppressWarnings(bamlss::bamlss(list(ff,~1,~1), data=d,
    family=degpd_bamlss(),sampler=FALSE,verbose=FALSE,...))
  a <- fit(y~x+offset(log(x)))
  b <- fit(y~x, offset=list(sigma=log(d$x)))
  coverage <- function(f,y,...) {
    set.seed(2761)
    predictive_coverage(f,y,method="plug-in",nsim=256,...)
  }
  expect_warning(ai <- coverage(a,d$y),"saved fitting offsets")
  expect_warning(an <- coverage(a,d$y,newdata=d),"saved fitting offsets")
  expect_identical(an$U,ai$U)
  bi <- coverage(b,d$y)
  bn <- coverage(b,d$y,newdata=d,offset=log(d$x))
  expect_identical(bn$U,bi$U)
  expect_error(coverage(b,d$y,newdata=d),"varying BAMLSS offset")
  cf <- unlist(coef(b))
  xi <- exp(cf[["xi.p.(Intercept)"]]); kap <- exp(cf[["kappa.p.(Intercept)"]])
  s <- exp(cf[["sigma.p.(Intercept)"]] + cf[["sigma.p.x"]]*d$x[1])*d$x[1]
  set.seed(2761)
  sample <- floor(s/xi*expm1(-xi*log1p(-runif(256)^(1/kap))))
  expect_equal(bi$U[1],as.numeric(quantile(sample,.975,type=1)))
})
