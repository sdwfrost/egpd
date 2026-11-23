test_that("a near-bound start can leave the logit plateau", {
  skip_on_cran()
  set.seed(505)
  x <- rcpdegpd(400, lambda = 1, sigma = 2, xi = .2, kappa = .8)
  fit <- function(st = NULL, maxit = 500) fit_aggregate_egpd(x,c(1,2,4),
    family = "cpdegpd", p = 1, q = 1, start = st, control = list(maxit = maxit))
  interior <- fit()
  st <- list(sigma_coef = interior$sigma_coef, lambda_coef = interior$lambda_coef,
             kappa = 4.999999, xi = interior$xi)
  f <- fit(st)
  expect_equal(f$convergence,0L)
  expect_gte(f$composite_loglik, interior$composite_loglik - .01)
  expect_lt(f$kappa,4.99)
})

test_that("aggregate bootstrap refits can move off the original shape bound", {
  skip_on_cran()
  set.seed(11)
  x <- rcpdegpd(100, lambda = .8, sigma = 2, xi = .2, kappa = .9)
  f <- fit_aggregate_egpd(x,c(1,2,4),family = "cpdegpd", p = 1, q = 1, keep_data = TRUE)
  b <- bootstrap_aggregate_egpd(f,nboot = 4,block_length = 8,seed = 42)
  expect_identical(b$n_failed,0L)
  expect_identical(b$n_nonconverged,0L)
  expect_gt(unname(b$se["kappa"]),.1)

  ## On the sixth block sample L-BFGS-B probes just below xi = .01.
  set.seed(42)
  for (i in 1:6) starts <- sample.int(93L,13L,replace = TRUE)
  idx <- as.vector(outer(0:7, starts, "+"))[seq_along(x)]
  fresh <- fit_aggregate_egpd(x[idx],c(1,2,4),family = "cpdegpd",p = 1,q = 1)
  expect_equal(fresh$convergence,0L)
  expect_gte(fresh$xi,.01)
  expect_true(is.finite(fresh$composite_loglik))
})
