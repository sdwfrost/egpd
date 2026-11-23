test_that("held-out intervals keep the training spline knot as the batch changes", {
  skip_if_not_installed("gamlss")
  d <- data.frame(x = rep((0:4)/4, each = 20))
  rate <- 4^(4*abs(d$x-.5))
  d$y <- rate + rep(c(-1, 1), 50)
  f <- gamlss::gamlss(y ~ splines::bs(x, df=2, degree=1), data=d,
                       family=gamlss.dist::PO(), trace=FALSE)
  expect_equal(as.numeric(f$mu.fv), rate, tolerance=1e-8)
  X <- cbind(1, d$x, pmax(0, d$x-.5))
  v <- c(1, .5, 0)
  se <- sqrt(drop(v %*% solve(crossprod(X, rate*X), v)))
  nsim <- 800L
  for (x in list(.5, c(.5, rep(.25, 120)))) for (m in c("parametric", "plug-in")) {
    set.seed(418)
    a <- predictive_coverage(f, rep(4,length(x)), data.frame(x=x), nsim=nsim, method=m)
    set.seed(418)
    mu <- if (m == "plug-in") 1 else exp(se*rnorm(nsim*length(x))[seq_len(nsim)])
    ref <- quantile(rpois(nsim, mu), c(.025,.975), type=1, names=FALSE)
    expect_equal(c(a$L[1], a$U[1]), ref)
  }
})

test_that("held-out GAMLSS intervals sum all formula offsets", {
  skip_if_not_installed("gamlss")
  d <- expand.grid(a=log(1:2), b=log(c(1,3)), i=1:20)
  rate <- 2*exp(d$a+d$b)
  d$y <- round(rate) + rep(c(-1,1), each=4, times=10)
  f <- gamlss::gamlss(y ~ offset(a)+offset(b), data=d, family=gamlss.dist::PO(), trace=FALSE)
  for (m in c("plug-in","parametric")) {
    set.seed(108)
    a <- predictive_coverage(f, 12, data.frame(a=log(2), b=log(3)), nsim=1500, method=m)
    set.seed(108)
    rate_new <- if (m == "plug-in") 12 else exp(rnorm(1500, log(12), 1/sqrt(sum(rate))))
    ref <- quantile(rpois(1500, rate_new), c(.025,.975), type=1, names=FALSE)
    expect_equal(c(a$L,a$U), ref)
  }
})

test_that("plug-in with a backfitted smooth keeps the other fitted terms", {
  skip_if_not_installed("gamlss")
  withr::local_package("gamlss")
  set.seed(406)
  d <- data.frame(x=rep((0:4)/4, each=20), t=runif(100), a=runif(100), b=runif(100))
  d$ab <- d$a+d$b
  d$y <- sin(6*d$t) + 8*abs(d$x-.5) + d$ab + rnorm(100, sd=.2)
  fit <- function(form) gamlss::gamlss(form, data=d, family=gamlss.dist::NO(),
                                        sigma.start=.2, sigma.fix=TRUE, trace=FALSE)
  f <- fit(y ~ pb(t) + splines::bs(x,df=2,degree=1) + offset(a) + offset(b))
  g <- fit(y ~ pb(t) + splines::bs(x,knots=.5,degree=1,Boundary.knots=c(0,1)) + offset(ab))
  expect_equal(as.numeric(f$mu.fv), as.numeric(g$mu.fv), tolerance=1e-8)
  nd <- d[c(44,rep(23,110)), , drop=FALSE]
  cover <- function(f) {
    set.seed(158)
    r <- predictive_coverage(f,nd$y,nd,data_fit=d,method="plug-in",nsim=100)
    c(r$L,r$U)
  }
  expect_equal(cover(f), cover(g), tolerance=1e-8)
})
