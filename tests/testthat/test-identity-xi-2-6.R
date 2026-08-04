## Identity link on xi for DEGPD carriers 2-6.
##
## Only model 1 could represent a bounded/short tail (xi < 0); models 2-6 forced
## xi = exp(eta) > 0, so on bounded count data they had nothing to say. Rather than
## re-derive five carriers symbolically (as derive_degpd1_identity.py did for model
## 1), the C++ is handed xi directly -- legitimate because in every carrier the
## scalar lxi appears ONLY inside exp(lxi), so it is used purely as a value and the
## expressions are analytic for xi < 0. What comes back is still differentiated
## w.r.t. log xi, and .identity_xi_chain() inverts that.
##
## The whole construction rests on that chain rule being right, which is what the
## finite-difference tests below check -- at NEGATIVE xi, the case that did not
## exist before.

testthat::skip_on_cran()

EXTRA <- list(`2` = c("lkappa1", "ldkappa", "logitp"), `3` = "ldelta",
              `4` = c("ldelta", "lkappa"), `5` = "lkappa", `6` = "lkappa")
NP <- c(`2` = 5L, `3` = 3L, `4` = 4L, `5` = 3L, `6` = 3L)

bounded_data <- function() {
  set.seed(3)
  data.frame(cases = rbinom(150, 25, 0.35))
}

id_fit <- function(m, d) {
  fl <- c(list(lsigma = cases ~ 1, xi = ~1),
          setNames(rep(list(~1), length(EXTRA[[m]])), EXTRA[[m]]))
  suppressWarnings(egpd(fl, data = d, family = "degpd",
                        degpd.args = list(m = as.integer(m), link = "identity"),
                        restarts = FALSE, trace = 0))
}

test_that("identity-link gradient matches finite differences, including xi < 0", {
  d <- bounded_data()
  for (m in names(EXTRA)) {
    fit <- id_fit(m, d)
    ld <- fit$likdata; np <- NP[[m]]
    f <- function(b) fit$likfns$d0(b, ld)
    for (xiv in c(-0.25, -0.08, 0.15, 0.4)) {
      b <- c(log(8), xiv, rep(0, np - 2L))
      analytic <- colSums(fit$likfns$d120(b, ld)[, seq_len(np), drop = FALSE])
      numeric <- vapply(seq_len(np), function(k) {
        h <- 1e-6 * max(1, abs(b[k])); bp <- bm <- b
        bp[k] <- b[k] + h; bm[k] <- b[k] - h
        (f(bp) - f(bm)) / (2 * h)
      }, 0)
      expect_equal(analytic, numeric, tolerance = 1e-5,
                   info = sprintf("m = %s, xi = %g", m, xiv))
    }
  }
})

test_that("identity-link Hessian matches finite differences, including xi < 0", {
  d <- bounded_data()
  for (m in names(EXTRA)) {
    fit <- id_fit(m, d)
    ld <- fit$likdata; np <- NP[[m]]
    f <- function(b) fit$likfns$d0(b, ld)
    for (xiv in c(-0.25, 0.3)) {
      b <- c(log(8), xiv, rep(0, np - 2L))
      a <- fit$likfns$d120(b, ld)
      analytic <- colSums(a[, (np + 1L):ncol(a), drop = FALSE])
      numeric <- c()
      for (i in seq_len(np)) for (j in i:np) {
        h <- 1e-4; bpp <- bpm <- bmp <- bmm <- b
        bpp[i] <- bpp[i] + h; bpp[j] <- bpp[j] + h
        bpm[i] <- bpm[i] + h; bpm[j] <- bpm[j] - h
        bmp[i] <- bmp[i] - h; bmp[j] <- bmp[j] + h
        bmm[i] <- bmm[i] - h; bmm[j] <- bmm[j] - h
        numeric <- c(numeric, (f(bpp) - f(bpm) - f(bmp) + f(bmm)) / (4 * h * h))
      }
      ## finite-difference truncation dominates here, so compare on a relative scale
      s <- max(1, max(abs(numeric)))
      expect_equal(analytic / s, numeric / s, tolerance = 1e-4,
                   info = sprintf("m = %s, xi = %g", m, xiv))
    }
  }
})

test_that("the two links give the same likelihood at the same positive xi", {
  ## the identity path must not perturb the log-link answer where both are valid:
  ## same model, same xi, same data => same negative log-likelihood
  d <- bounded_data()
  for (m in names(EXTRA)) {
    np <- NP[[m]]
    fl_log <- c(list(lsigma = cases ~ 1, lxi = ~1),
                setNames(rep(list(~1), length(EXTRA[[m]])), EXTRA[[m]]))
    f_log <- suppressWarnings(egpd(fl_log, data = d, family = "degpd",
               degpd.args = list(m = as.integer(m)), restarts = FALSE, trace = 0))
    f_id <- id_fit(m, d)
    xi <- 0.3
    expect_equal(f_log$likfns$d0(c(log(8), log(xi), rep(0, np - 2L)), f_log$likdata),
                 f_id$likfns$d0(c(log(8), xi, rep(0, np - 2L)), f_id$likdata),
                 tolerance = 1e-10, info = sprintf("m = %s", m))
  }
})

test_that("every carrier can reach a bounded tail on bounded data", {
  ## the point of the exercise: xi < 0 was unreachable for m = 2-6
  d <- bounded_data()
  for (m in as.character(1:6)) {
    ex <- if (m == "1") "lkappa" else EXTRA[[m]]
    fl <- c(list(lsigma = cases ~ 1, xi = ~1),
            setNames(rep(list(~1), length(ex)), ex))
    fit <- suppressWarnings(egpd(fl, data = d, family = "degpd",
             degpd.args = list(m = as.integer(m), link = "identity"),
             restarts = FALSE, trace = 0))
    xi <- as.data.frame(predict(fit, newdata = d[1, , drop = FALSE],
                                type = "response"))$shape[1]
    expect_true(is.finite(xi) && xi < 0,
                info = sprintf("m = %s gave xi-hat = %g on binomial counts", m, xi))
  }
})

test_that("an unrecognised link errors instead of being silently ignored", {
  ## it used to fall through to the log link, so a caller who asked for a bounded
  ## tail got a strictly positive xi with no indication anything was wrong
  d <- bounded_data()
  expect_error(
    egpd(list(lsigma = cases ~ 1, lxi = ~1, ldelta = ~1), data = d,
         family = "degpd", degpd.args = list(m = 3L, link = "probit"),
         restarts = FALSE, trace = 0),
    "must be")
})
