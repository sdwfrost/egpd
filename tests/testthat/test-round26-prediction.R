test_that("prediction treats another model's frame as covariate data", {
  d <- data.frame(x=rep(0:1,each=50),y=rep(c(0:4,3:7),each=10))
  f <- egpd(list(y~x,~1),d,family="cmp",fixed=list(nu=1),trace=-1)
  nd <- data.frame(x=c(0,.5,1))
  mf <- model.frame(~x,nd)
  expect_equal(predict(f,mf,type="response"),predict(f,nd,type="response"))
  expect_equal(predict(f,mf,type="response")$mu,c(2,sqrt(10),5),tolerance=1e-6)
})

test_that("BAMLSS held-out intervals use training polynomial and factor coordinates", {
  skip_if_not_installed("bamlss")
  skip_on_cran()
  if (!"package:bamlss" %in% search()) {
    suppressPackageStartupMessages(library(bamlss))
    on.exit(detach("package:bamlss", unload=FALSE, character.only=TRUE), add=TRUE)
  }
  set.seed(26981); n <- 72L
  d <- data.frame(x=runif(n,-1,1),g=factor(rep(c("a","b","c"),n/3)))
  d$y <- floor(exp(1+d$x+.4*d$x^2+c(0,.35,.85)[d$g])/.2*expm1(-.2*log1p(-runif(n))))
  f <- suppressWarnings(bamlss::bamlss(list(y~poly(x,2)+g,xi~1,kappa~s(x,k=4)),data=d,
    family=degpd_bamlss(),n.iter=160,burnin=60,thin=1,verbose=FALSE))
  S <- as.matrix(bamlss::samples(f))
  nd <- data.frame(x=.6,g=factor("c",levels=c("c","b","a")))
  b <- as.vector(predict(poly(d$x,2),nd$x))
  eta <- S[,"sigma.p.(Intercept)"] + S[,"sigma.p.poly(x, 2)1"]*b[1] +
    S[,"sigma.p.poly(x, 2)2"]*b[2] + S[,"sigma.p.gc"]
  set.seed(26)
  draws <- egpd:::.bamlss_param_draws(f,nd,f$family$names,nrow(S))
  expect_equal(unname(log(draws$sigma[,1])),unname(eta[attr(draws,"draws")]),tolerance=1e-10)
  many <- nd[c(1,1,1),,drop=FALSE]; many$x <- c(.6,-.9,-.85)
  many$g <- factor(c("c","a","b"),levels=levels(d$g))
  smooth <- as.matrix(predict(f,many,model="kappa",type="parameter",FUN=identity))
  draw3 <- egpd:::.bamlss_param_draws(f,many,f$family$names,nrow(S))
  expect_equal(unname(draw3$kappa),unname(t(smooth[,attr(draw3,"draws"),drop=FALSE])))
  for (method in c("parametric","plug-in")) {
    calc <- function(z) {
      set.seed(267)
      res <- predictive_coverage(f,rep(3,nrow(z)),z,method=method,nsim=150)
      c(res$L[1],res$U[1])
    }
    expect_identical(calc(nd),calc(many))
  }
})
