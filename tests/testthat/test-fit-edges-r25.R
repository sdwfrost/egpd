test_that("zero-inflation inference respects the exact conditional atom", {
  y <- rep(0:3,each=3)
  fit <- function(s, ...) distfit(y, family="zidegpd",
    fix.arg=list(sigma=s,xi=1,kappa=1), ...)
  expect_warning(f <- fit(3), "natural boundary")
  expect_equal(unname(coef(f)), 0)
  expect_true(all(is.na(c(f$sd, confint(f)))))
  s <- 3 + 2^-51
  g <- fit(s)
  expect_equal(unname(coef(g))/((s-3)/(4*s)), 1, tolerance=1e-13)
  expect_true(is.finite(g$sd))
  expect_warning(box <- fit(3, lower=qlogis(.1)), "optimizer bound")
  expect_equal(unname(coef(box)), .1)
  # Rational equality not representable as a dyadic zero probability.
  y2 <- c(0,rep_len(1:6,342))
  expect_warning(f2 <- distfit(y2, family="zidegpd",
    fix.arg=list(sigma=7,xi=-1,kappa=3)), "natural boundary")
  expect_equal(unname(coef(f2)), 0)
})

test_that("POT support changes withhold Wald inference without losing the fit", {
  x <- 16 + c(0,0,0,2,2,6,8,10,16,16,17)
  for (h in c(.001, .00001)) {
    expect_warning(f <- fit_degpd_pot(x,16,control=list(ndeps=c(h,h))), "change in count support")
    expect_equal(f$convergence,0L)
    expect_equal(-unname(coef(f)["sigma"]/coef(f)["xi"]),18,tolerance=1e-6)
    expect_true(all(is.na(c(f$sd, f$xi.ci, confint(f)))))
    expect_error(vcov(f),"change in count support")
  }
  expect_warning(st <- degpd_pot_stability(x, thresholds=16), "change in count support")
  expect_true(is.finite(st$xi) && is.na(st$xi_lo) && is.na(st$xi_hi))
  p <- (seq_len(300)-.5)/300
  regular <- fit_degpd_pot(floor(20*expm1(-.2*log1p(-p))),0)
  expect_true(all(is.finite(regular$xi.ci)))
})

test_that("inline NA replacement uses the same splice sample as a precomputed offset", {
  skip_if_not_installed("qgam")
  set.seed(25112)
  d <- data.frame(x=runif(250,-1,1))
  d$y <- rnbinom(nrow(d),2,mu=6*exp(d$x/3))
  d$off <- ifelse(d$y > 12, NA_real_, 0)
  d$filled <- rep(0, nrow(d))
  fit <- function(bulk,tail) {
    set.seed(25113)
    fit_spliced_degpd(bulk,d,bulk_probs=c(.3,.55,.75),threshold_prob=.75,
                      tail_formula=tail, egpd_control=list(fix.arg=list(xi=0,kappa=1)))
  }
  ref <- fit(y~x+offset(filled),~offset(filled))
  a <- fit(y~x+offset(ifelse(is.na(off),0,off)),~offset(ifelse(is.na(off),0,off)))
  expect_equal(c(a$n,a$n_exceed), c(ref$n,ref$n_exceed))
  expect_equal(a$tail_fit$loglik,ref$tail_fit$loglik)
  nd <- d[which(is.na(d$off))[1:4],]
  expect_equal(predict(a,nd,type="parameter"),predict(ref,nd,type="parameter"))
})
