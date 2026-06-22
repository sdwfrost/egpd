## Finite-difference verification of the GAM-likelihood analytic/bounded derivatives:
## the per-observation d12 (d120) gradient and Hessian, assembled to coefficient space,
## must match numerical differences of d0. Covers DEGPD models 1-6 (including the bounded
## xi.max link), ZIDEGPD models 5-6, and EGPD models 5-6.

## Assemble coef-space gradient/Hessian from d120 and compare to numerical d0 derivatives.
fd_oracle <- function(fns, npar, y, p, xi.max = Inf) {
  n <- length(y)
  X <- c(list(cbind(1, seq(-1, 1, length.out = n))),
         replicate(npar - 1, matrix(1, n, 1), simplify = FALSE))
  nbk <- vapply(X, ncol, 1L)
  ld <- list(X = X, y = cbind(y), dupid = 0L, duplicate = 0L,
             offsets = replicate(npar, numeric(0), simplify = FALSE),
             idpars = rep(seq_len(npar), nbk), censored = FALSE, xi.max = xi.max)
  cg <- function(pp) {
    M <- fns$d120(pp, ld)
    unlist(lapply(seq_len(npar), function(i) colSums(M[, i] * X[[i]])))
  }
  col_ij <- function(a, b) {            # packed upper-tri column (1-based) for 0-based (a,b)
    off <- 0L; if (a > 0) for (r in 0:(a - 1)) off <- off + (npar - r)
    npar + off + (b - a) + 1L
  }
  ch <- function(pp) {
    M <- fns$d120(pp, ld); st <- c(0, cumsum(nbk)); nb <- sum(nbk)
    H <- matrix(0, nb, nb)
    for (a in 0:(npar - 1)) for (b in a:(npar - 1)) {
      cc <- col_ij(a, b); Hab <- t(X[[a + 1]]) %*% (X[[b + 1]] * M[, cc])
      H[(st[a + 1] + 1):st[a + 2], (st[b + 1] + 1):st[b + 2]] <- Hab
      if (a != b) H[(st[b + 1] + 1):st[b + 2], (st[a + 1] + 1):st[a + 2]] <- t(Hab)
    }
    H
  }
  ng <- vapply(seq_along(p), function(i) {
    h <- 1e-6; pp <- p; pm <- p; pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
    (fns$d0(pp, ld) - fns$d0(pm, ld)) / (2 * h)
  }, numeric(1))
  nh <- matrix(0, length(p), length(p))
  for (i in seq_along(p)) {
    h <- 1e-5; pp <- p; pm <- p; pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
    nh[, i] <- (cg(pp) - cg(pm)) / (2 * h)
  }
  nh <- (nh + t(nh)) / 2
  list(grad = max(abs(cg(p) - ng)), hess = max(abs(ch(p) - nh)))
}

ydisc <- c(0, 0, 1, 2, 3, 5, 8, 12, 4, 1)            # discrete counts (DEGPD / ZIDEGPD)
ycont <- c(0.5, 1, 2, 3, 5, 8, 0.2, 4, 1.5, 2.5)     # positive continuous (EGPD)
npar_degpd <- c(`1` = 3, `2` = 5, `3` = 3, `4` = 4, `5` = 3, `6` = 3)
mkp <- function(np) c(0.7, 0.3, -0.3, if (np > 2) rep(-0.4, np - 2) else numeric(0))  # sum(nbk) = np + 1

for (m in 1:6) {
  np <- npar_degpd[[as.character(m)]]
  fns <- get(paste0(".degpd", m, "fns"), envir = asNamespace("egpd"))
  test_that(sprintf("DEGPD model %d: analytic d12 matches finite differences", m), {
    for (xm in c(Inf, 0.5)) {        # Inf = log link, 0.5 = bounded shape link
      r <- fd_oracle(fns, np, ydisc, mkp(np), xi.max = xm)
      expect_lt(r$grad, 1e-4)
      expect_lt(r$hess, 1e-2)
    }
  })
}

for (m in c(5, 6)) {
  fnsz <- get(paste0(".zidegpd", m, "fns"), envir = asNamespace("egpd"))
  test_that(sprintf("ZIDEGPD model %d: analytic d12 matches finite differences", m), {
    r <- fd_oracle(fnsz, 4, ydisc, c(0.6, 0.3, -0.2, -0.4, 0.2))
    expect_lt(r$grad, 1e-4)
    expect_lt(r$hess, 1e-2)
  })
  fnse <- get(paste0(".egpd", m, "fns"), envir = asNamespace("egpd"))
  test_that(sprintf("EGPD model %d: analytic d12 matches finite differences", m), {
    r <- fd_oracle(fnse, 3, ycont, c(0.5, 0.2, 0.25, 0.3))
    expect_lt(r$grad, 1e-4)
    expect_lt(r$hess, 1e-2)
  })
}
