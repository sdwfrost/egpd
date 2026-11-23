test_that("continuous free endpoints below shape -1 cannot give regular Wald inference", {
  x <- (1:60)/10
  fix <- list(xi = -2, kappa = 1)
  ## f(x) = 1/(s*sqrt(1-2*x/s)); likelihood diverges as s decreases to 12.
  for (method in c("mle","bernstein"))
    expect_error(distfit(x, family = "egpd", method = method, fix.arg = fix,
                         start = list(sigma = 14)), "Unbounded continuous likelihood")
  expect_error(distfit(c(0,x), family = "ziegpd", fix.arg = c(fix,list(pi = .2))),
               "Unbounded continuous likelihood")
  ## A fixed endpoint or an explicit bound excluding the singularity is
  ## not the same problem; nor is a floored density.
  interior <- distfit(x, family = "egpd", fix.arg = list(xi = -2, sigma = 14))
  expect_true(is.finite(as.numeric(logLik(interior))))
  bounded <- distfit(x, family = "egpd", fix.arg = fix, optim.method = "Brent",
                     lower = log(15), upper = log(25), hessian = FALSE)
  expect_equal(unname(coef(bounded)), 15, tolerance = 1e-4)
  count <- suppressWarnings(distfit(floor(x), family = "degpd", fix.arg = fix))
  expect_true(is.finite(as.numeric(logLik(count))))
  expfit <- distfit(x, family = "egpd", fix.arg = list(kappa = 1, xi = 0))
  expect_equal(unname(expfit$sd), mean(x)/sqrt(length(x)), tolerance = 1e-5)
})

test_that("an all-zero Bell sample has no interior MLE", {
  for (family in c("bell","zibell"))
    expect_error(distfit(rep(0,100), family = family), "No interior Bell maximum")
  y <- c(rep(0,85),rep(1,14),2)
  f <- distfit(y, family = "bell")
  theta <- unname(coef(f))
  expect_equal(theta*exp(theta), mean(y), tolerance = 1e-6)
  expect_equal(unname(f$sd), (sum(y)/theta^2 + length(y)*exp(theta))^-.5, tolerance = 1e-5)
})
