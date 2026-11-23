test_that("automatic POT thresholds retain large integer-valued counts", {
  x <- floor(8*qexp(((1:600)-.5)/600))
  eval_grid <- function(x, ...) degpd_pot_stability(x, hessian=FALSE, ...)
  lo <- eval_grid(x)
  hi <- eval_grid(x+3e9)
  expect_equal(hi$threshold-3e9, lo$threshold)
  expect_equal(hi[c("sigma","xi","n_exceed")], lo[c("sigma","xi","n_exceed")])
  expect_equal(eval_grid(x+3e9,thresholds=lo$threshold+3e9), hi)
})

test_that("hurdle quantiles distinguish either side of the zero atom", {
  set.seed(2241)
  d <- data.frame(x=runif(200,-1,1),h=0)
  d$y <- ifelse(runif(200)<plogis(d$x),1+rgeom(200,.4),0)
  f <- fit_hurdle(list(y~1,~1,~1),d,zero.formula=~0+x+offset(h),
      family="degpd",degpd.args=list(m=1,link="identity"),fixed=list(xi=0,kappa=1),trace=-1)
  # At x=0 the reported link is exactly the offset, with no fitting error.
  h <- c(-0x1.193ea7aad030bp+0,-0x1.193ea7aad030ap+0,-0x1.193ea7aad0308p+0)
  nd <- data.frame(x=0,h=h)
  a <- predict(f,nd,prob=c(.75-1e-8,.75,.75+1e-8))
  expect_equal(unname(as.matrix(a)), cbind(c(0,0,0),c(0,1,1),c(1,1,1)))
  # A rounded qlogis equality can fall on either side of the true atom.
  nd <- data.frame(x=0,h=c(0x1.3553e218d610fp+0,0x1.35891deb9a3bp-2,0))
  expect_equal(predict(f,nd[1,,drop=FALSE],prob=.23)[[1]],1)
  expect_equal(predict(f,nd[2,,drop=FALSE],prob=.425)[[1]],0)
  expect_equal(predict(f,nd[3,,drop=FALSE],prob=.5)[[1]],0)
})
