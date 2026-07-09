## Randomized quantile residuals (rqresid) for the count families added to
## egpd(): PIG/ZIPIG, GPIG/ZIGPIG and Bell/ZIBell. For a correctly specified
## discrete model the residuals should be approximately standard normal. We
## check structural correctness plus approximate N(0, 1) behaviour averaged
## over several seeds (a single heavy-tailed realisation can have an inflated
## sample sd, so we average rather than assert on one fit).

rqr_summary <- function(family, gen, forms, S = 6) {
  m <- s <- numeric(S)
  for (i in seq_len(S)) {
    set.seed(200 + i)
    y <- gen()
    fit <- egpd(forms, data = data.frame(y = y), family = family, trace = 0)
    r <- rqresid(fit, seed = i)
    expect_length(r, length(y))
    ## An extreme heavy-tail draw can saturate the CDF (F(y) rounds to 1),
    ## giving qnorm(1) = Inf which rqresid maps to NA; almost all are finite.
    expect_gt(mean(is.finite(r)), 0.95)
    m[i] <- mean(r, na.rm = TRUE); s[i] <- sd(r, na.rm = TRUE)
  }
  list(mean = mean(m), sd = mean(s))
}

test_that("rqresid is approximately N(0,1) for PIG and ZIPIG", {
  skip_on_cran()
  a <- rqr_summary("pig",  function() egpd:::rpig(1500, mu = 3, sigma = 0.8),
                   list(lmu = y ~ 1, lsigma = ~ 1))
  expect_lt(abs(a$mean), 0.1); expect_equal(a$sd, 1, tolerance = 0.1)
  b <- rqr_summary("zipig", function() egpd:::rzipig(1500, mu = 3, sigma = 0.8, pi = 0.2),
                   list(lmu = y ~ 1, lsigma = ~ 1, logitpi = ~ 1))
  expect_lt(abs(b$mean), 0.1); expect_equal(b$sd, 1, tolerance = 0.1)
})

test_that("rqresid is approximately N(0,1) for GPIG and ZIGPIG", {
  skip_on_cran()
  a <- rqr_summary("gpig",  function() egpd:::rgpig(1500, a = 0.45, b = 2.2, c = 0.8),
                   list(lmu = y ~ 1, logita = ~ 1, logitc = ~ 1))
  expect_lt(abs(a$mean), 0.1); expect_equal(a$sd, 1, tolerance = 0.12)
  b <- rqr_summary("zigpig", function() egpd:::rzigpig(1500, a = 0.45, b = 2.2, c = 0.8, pi = 0.2),
                   list(lmu = y ~ 1, logita = ~ 1, logitc = ~ 1, logitpi = ~ 1))
  expect_lt(abs(b$mean), 0.1); expect_equal(b$sd, 1, tolerance = 0.12)
})

test_that("rqresid is approximately N(0,1) for Bell and ZIBell", {
  skip_on_cran()
  a <- rqr_summary("bell",  function() egpd:::rbell(1500, theta = 0.9),
                   list(lmu = y ~ 1))
  expect_lt(abs(a$mean), 0.1); expect_equal(a$sd, 1, tolerance = 0.1)
  b <- rqr_summary("zibell", function() egpd:::rzibell(1500, theta = 0.9, pi = 0.2),
                   list(lmu = y ~ 1, logitpi = ~ 1))
  expect_lt(abs(b$mean), 0.1); expect_equal(b$sd, 1, tolerance = 0.1)
})
