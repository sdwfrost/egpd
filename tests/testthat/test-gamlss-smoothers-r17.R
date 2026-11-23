test_that("plug-in coverage delegates formula smooths, but checks real aliases", {
  skip_if_not_installed("gamlss")
  withr::local_package("gamlss")
  d <- expand.grid(temp = seq(0,1,length.out = 40), grp = 0:1, rep = 1:7)
  d$z <- d$grp
  d$y <- qpois((d$rep-.5)/7, exp(1 + .3*d$grp + .55*sin(5*d$temp)))
  nd <- data.frame(temp = c(.16,.42,.78), grp = c(0,1,1), z = c(0,1,1))
  for (fo in list(y ~ lo(~temp), y ~ lo(~temp) + grp + z)) {
    f <- gamlss::gamlss(fo, data = d, family = gamlss.dist::PO(), trace = FALSE)
    expect_true(f$converged)
    means <- predict(f, newdata = nd, type = "response", data = d)
    set.seed(8171)
    ref <- vapply(means, function(m) quantile(rpois(2000,m), c(.025,.975), type = 1,
                                               names = FALSE), numeric(2))
    cover <- function(nd) predictive_coverage(f, rep(4,nrow(nd)), nd, data_fit = d,
                                               nsim = 2000, method = "plug-in")
    set.seed(8171)
    ans <- cover(nd)
    expect_equal(as.numeric(ans$L), as.numeric(ref[1,]))
    expect_equal(as.numeric(ans$U), as.numeric(ref[2,]))
  }
  expect_error(cover(transform(nd, z = 2)), "Non-estimable 'mu'")
})
