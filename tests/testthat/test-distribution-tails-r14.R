test_that("a concentrated beta carrier inverts the truncated law, including count cells", {
  p <- c(.02, .05, .5, .99)
  k <- 1e12
  # For t = 1/2 - X, X ~ Beta(k,k), 4*t^2 ~ Beta(1/2,k).
  z <- pbeta((15/16)^2, .5, k)
  s <- (16/15) * sqrt(qbeta((1-p)*z, .5, k))
  expect_equal(qEGPD6(p, mu=1, sigma=1, nu=k), 1/s-1, tolerance=1e-10)
  expect_equal(qDEGPD6(p, mu=1, sigma=.5, nu=k), ceiling(2*(1/sqrt(s)-1))-1)
  set.seed(55082)
  xx <- rEGPD6(5000, mu=1, sigma=1, nu=k)
  expect_lt(abs(mean(xx <= 1/s[2]-1) - .05), .02)
})

test_that("concentrated beta tails retain their transition at large hazards", {
  for (lk in c(640, 709)) {
    k <- exp(lk); A <- lk/2 + c(-2, 0, 2, 8)
    # k times a Beta(1/2,k) converges to Gamma(1/2,1) at these shapes.
    x <- (15/16)^2 * exp(lk - 2*A)
    for (lower in c(TRUE,FALSE)) {
      lp <- pegpd(A, sigma=1, xi=0, type=3, kappa=k, lower.tail=lower, log.p=TRUE)
      expect_equal(lp, pgamma(x, .5, lower.tail=!lower, log.p=TRUE), tolerance=1e-11)
    }
    p <- c(.01,.5,.75,.99)
    target <- log(15/16) + (lk - log(qgamma(p,.5,lower.tail=FALSE)))/2
    expect_equal(qegpd(p, sigma=1, xi=0, type=3, kappa=k), target, tolerance=1e-11)
  }
  expect_equal(pegpd(800, sigma=1, xi=0, type=3, kappa=1, lower.tail=FALSE, log.p=TRUE),
               -800, tolerance=1e-13)
})

test_that("tiny type-4 and type-5 CDFs retain large delta times the hazard", {
  q <- c(1e-110, 1e-110, 1e-110, 1e-308)
  d <- c(1e240, 1e110, 1, 1e308)
  # For d*q=1 the leading coefficient of q is exp(-1); d=1 gives u^2 exactly.
  lg <- c(log(q[1]), log(q[2])-1, 2*log(-expm1(-q[3])), log(q[4])-1)
  for (type in 4:5) {
    pow <- if (type==4L) 1 else .5
    got <- pegpd(q, kappa=1, delta=d, sigma=1, xi=0, type=type, log.p=TRUE)
    expect_equal(got, pow*lg, tolerance=1e-13)
    expect_equal(p.G(-expm1(-q), type=type, kappa=1, delta=d), exp(pow*lg),
                 tolerance=1e-12, scale=max(exp(pow*lg)))
  }
})

test_that("composite Pareto accumulates rather than loses the mass at the splice", {
  x <- c(1,1.01,1.05)
  for (alpha in c(1, 1e-15, 1e-14)) {
    body <- alpha * expm1(1) / (2 + alpha * expm1(1))
    F <- body - (1-body)*expm1(-alpha*log1p((x-1)/2))
    got <- pcomppareto(x, "exp", alpha=alpha, theta=1, rate=1)
    expect_true(all(diff(got)>0))
    expect_equal(got/F, rep(1,3), tolerance=1e-13)
    for (lower in c(TRUE,FALSE)) {
      lp <- if (lower) log(F) else log1p(-F)
      expect_equal(qcomppareto(lp, "exp", alpha=alpha, theta=1, rate=1,
                               lower.tail=lower, log.p=TRUE), x, tolerance=1e-11)
    }
  }
  ## Here ordinary F has only a few representable values, but log F does not.
  alpha <- .Machine$double.xmin * .Machine$double.eps
  lp <- log(alpha) + log(expm1(1)/2 + log1p((x-1)/2))
  expect_equal(pcomppareto(x, "exp", alpha=alpha, theta=1, rate=1, log.p=TRUE),
               lp, tolerance=1e-13)
  expect_equal(qcomppareto(lp, "exp", alpha=alpha, theta=1, rate=1, log.p=TRUE),
               x, tolerance=1e-11)
})

test_that("small unequal powers retain mixture survivors and cell mass on the log scale", {
  w <- .3
  for (k in c(.1, 1e-17, 1e-18)) {
    y <- c(.5,1)
    surv <- function(x) {
      lH <- log(-expm1(-x))
      -w*expm1(k*lH)-(1-w)*expm1(2*k*lH)
    }
    ref <- log(surv(y))
    expect_equal(pegpd(y, w,k,2*k,1,0,6, lower.tail=FALSE, log.p=TRUE), ref,
                 tolerance=1e-12)
    expect_equal(qegpd(ref, w,k,2*k,1,0,6, lower.tail=FALSE, log.p=TRUE), y,
                 tolerance=1e-12)
    expect_equal(ddiscegpd(1, w,k,2*k,5,0,6, log=TRUE),
                 log(surv(.2)-surv(.4)), tolerance=1e-12)
  }
})
