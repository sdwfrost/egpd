## Issue #6. egpd() builds its DEGPD multi-start by editing match.call() and
## re-evaluating it. match.call() records argument *expressions*, so the re-evaluation
## has to happen in the caller's frame -- and parent.frame() called from inside the
## lapply() over starts resolved to lapply()'s frame instead. Any argument the caller
## passed as a local variable was then invisible, every restart died with "object not
## found", tryCatch turned that into NULL, and the fit silently fell through to a
## single start.
##
## It was invisible because nothing failed loudly: the answer was just worse. On a
## Tycho measles series the single-start optimum sat 2,416 log-likelihood units below
## the multi-start one, with xi-hat 0.024 against 0.326 -- a different tail, not a
## slightly worse fit.

testthat::skip_on_cran()

## every argument inlined: this path always worked
fit_inlined <- function(y) {
  d <- data.frame(cases = y)
  egpd(list(lsigma = cases ~ 1, lxi = ~1, lkappa = ~1), data = d, family = "degpd",
       degpd.args = list(m = 1), restarts = TRUE)
}

## the pattern that used to break: arguments arrive as the caller's local variables
fit_from_locals <- function(y) {
  d  <- data.frame(cases = y)
  fl <- list(lsigma = cases ~ 1, lxi = ~1, lkappa = ~1)
  mm <- 1L
  egpd(fl, data = d, family = "degpd", degpd.args = list(m = mm), restarts = TRUE)
}

test_that("restarts survive arguments passed as local variables (#6)", {
  set.seed(4)
  y <- rnbinom(400, size = 1.5, mu = 30)

  a <- suppressWarnings(fit_inlined(y))
  b <- suppressWarnings(fit_from_locals(y))

  ## Same data, same model, same restarts: the optimum must not depend on whether the
  ## caller happened to inline its arguments. Before the fix these differed, because
  ## (b) quietly used one start.
  expect_equal(as.numeric(a$logLik), as.numeric(b$logLik), tolerance = 1e-6)
})

test_that("supplying inits still skips restarts without warning (#6)", {
  ## restarts are skipped by design when the caller supplies inits, so the new
  ## warning must not fire on that path
  set.seed(4)
  y <- rnbinom(200, size = 1.5, mu = 20)
  d <- data.frame(cases = y)
  expect_no_warning(
    egpd(list(lsigma = cases ~ 1, lxi = ~1, lkappa = ~1), data = d, family = "degpd",
         degpd.args = list(m = 1), inits = c(log(mean(y) + 1), log(0.1), 0),
         restarts = TRUE))
})
