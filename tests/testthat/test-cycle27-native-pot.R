test_that("a mixed tensor and random effect has usable smooth tests", {
  set.seed(2759); n <- 256L
  d <- data.frame(x=runif(n),z=runif(n),g=factor(rep(letters[1:8],32)))
  d$y <- .5*sin(d$x*6)+.4*d$z+rnorm(n)
  ff <- y~t2(x,z,k=4)+s(g,bs="re")
  sp <- c(.8,.9,1.1,2)
  cf <- list(d0=function(b,ld) -sum(dnorm(ld$y,drop(ld$X[[1]]%*%b),1,log=TRUE)),
             d120=function(b,ld) cbind(drop(ld$X[[1]]%*%b)-ld$y[,1],1))
  ref <- mgcv::gam(ff,data=d,sp=sp,scale=1,method="REML")
  for (ns in c(1e20,160)) {
    f <- egpd(list(mu=ff),d,family="custom",custom.fns=cf,sp=sp,trace=-1,
              maxspline=ns,removeData=TRUE)
    tab <- summary(f)
    expect_equal(f$convergence,0L)
    expect_true(all(is.finite(as.matrix(tab[[2]][[1]]))))
    if (ns > n) expect_equal(unname(unlist(tab[[1]][[1]][1,1:2])),
                             unname(summary(ref)$p.table[1,1:2]),tolerance=1e-6)
  }
})

test_that("vector formula expressions may not expand a prediction batch", {
  d <- data.frame(x=rep(0:1,each=25),y=rep(c(0:4,3:7),each=5))
  exposure <- seq(-.1,.1,length.out=nrow(d))
  f <- egpd(list(y~I(x+exposure),~1),d,family="cmp",fixed=list(nu=1),trace=-1)
  nd <- data.frame(x=c(0,1))
  expect_error(predict(f,nd,type="response"),"newdata row")
  expect_error(simulate(f,nsim=2,newdata=nd,type="response"),"newdata row")
  expect_error(confint(f,newdata=nd),"newdata row")
  nd$exposure <- c(.03,.07)
  expect_equal(nrow(predict(f,nd,type="response")),2L)
})

test_that("POT validates counts before filtering and honors named fixed parameters", {
  set.seed(2757)
  y <- floor(rexp(260,1/3))
  for (bad in c(-1,-999,.5,Inf)) {
    expect_error(fit_degpd_pot(c(bad,y),2),"non-negative integer")
    expect_error(degpd_pot_stability(c(bad,y)),"non-negative integer")
  }
  f <- fit_degpd_pot(c(NA_real_,y),2,fix.arg=list(xi=0))
  z <- y[y>=2]-2
  expect_equal(unname(coef(f)),1/log1p(1/mean(z)),tolerance=1e-4)
  expect_equal(as.numeric(logLik(f)),sum(dgeom(z,1/(1+mean(z)),log=TRUE)),tolerance=1e-8)
  expect_equal(f$prop.exceed,mean(y>=2))
  st <- degpd_pot_stability(y,c(1,2),fix.arg=list(xi=0))
  expect_equal(st$xi,c(0,0))
  expect_true(all(is.na(st$xi_lo) & is.na(st$xi_hi)))
  expect_error(fit_degpd_pot(y,2,fix.arg=list(kappa=2)),"conflicts")
  b <- fit_degpd_pot(rep(0:3,each=15),0,fix.kappa=FALSE,fix.arg=list(sigma=4,xi=-1))
  expect_equal(b$upper.endpoint,3)
})
