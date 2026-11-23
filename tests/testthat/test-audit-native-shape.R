audit_native_data <- function(y, n = 3L) {
  list(X = rep(list(matrix(1, length(y), 1)), n), y = cbind(y),
       idpars = seq_len(n), dupid = 0L, duplicate = 0L,
       offsets = rep(list(numeric()), n), censored = FALSE, xi.max = Inf)
}

# Independent cell integration of the continuous density, not a difference of
# native CDFs or a differentiation of the implementation's Taylor polynomial.
audit_native_reference <- function(b, y, m = 1L) {
  sigma <- exp(b[1]); xi <- b[2]; k <- exp(b[3])
  vapply(y, function(yy) {
    upper <- if (xi < 0) min(yy + 1, -sigma/xi) else yy + 1
    if (upper <= yy) return(1e20)
    logdensity <- function(t) {
      e <- if (xi == 0) -t/sigma else -log1p(xi*t/sigma)/xi
      H <- -expm1(e)
      carrier <- switch(m,
        k * H^(k-1),
        {
          k2 <- k + exp(b[4]); p <- plogis(b[5])
          p*k*H^(k-1) + (1-p)*k2*H^(k2-1)
        },
        (k+1)/k * -expm1(k*e),
        {
          G <- H + expm1(k*e)*exp(e)/k
          a <- exp(b[4])/2
          a * G^(a-1) * (k+1)/k * -expm1(k*e)
        },
        sqrt(k)*dnorm(-sqrt(k)*exp(e))/(0.5-pnorm(-sqrt(k))),
        (15/32)*dbeta(0.5-(15/32)*exp(e), k, k)/
          (pbeta(0.5, k, k)-pbeta(1/32, k, k)))
      log(carrier) + e - log(sigma) - log1p(xi*t/sigma)
    }
    scale <- logdensity(if (yy == 0) (yy+upper)/2 else yy)
    -scale-log(integrate(function(t) exp(logdensity(t)-scale),
                        yy, upper, rel.tol = 2e-12, subdivisions = 500L)$value)
  }, numeric(1))
}

audit_native_fd <- function(f, b, step = 0.001) {
  one <- function(h) {
    n <- length(b); f0 <- f(b)
    shift <- function(i, d) { z <- b; z[i] <- z[i] + d; z }
    g <- vapply(seq_len(n), function(i)
      (f(shift(i, h))-f(shift(i, -h)))/(2*h), f0)
    H <- list()
    for (i in seq_len(n)) for (j in i:n) {
      if (i == j) {
        z <- (f(shift(i, h))-2*f0+f(shift(i, -h)))/h^2
      } else {
        pp <- pm <- mp <- mm <- b
        pp[c(i,j)] <- pp[c(i,j)] + h
        mm[c(i,j)] <- mm[c(i,j)] - h
        pm[i] <- pm[i]+h; pm[j] <- pm[j]-h
        mp[i] <- mp[i]-h; mp[j] <- mp[j]+h
        z <- (f(pp)-f(pm)-f(mp)+f(mm))/(4*h^2)
      }
      H[[length(H)+1L]] <- z
    }
    cbind(g, do.call(cbind, H))
  }
  (4*one(step/2)-one(step))/3
}

test_that("identity native likelihood and d12 retain the true exponential limit", {
  y <- c(0, 1, 7, 15)
  # z = xi*(y+1)/sigma crosses the log1p(z)/z Taylor boundary.
  shapes <- c(0, -1e-6, 1e-6, -0.02, 0.02, -0.01875-1e-9,
              -0.01875+1e-9, 0.01875-1e-9, 0.01875+1e-9, -0.15, 0.4)
  for (m in 1:6) {
    n <- c(3L,5L,3L,4L,3L,3L)[m]
    fns <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    ld <- audit_native_data(y, n)
    for (xi in shapes) {
      b <- c(log(3), xi, log(1.5), if (m == 2) c(log(0.8), 0.3)
             else if (m == 4) log(1.7))
      ref <- function(z) audit_native_reference(z, y, m)
      expect_equal(fns$d0(b, ld), sum(ref(b)), tolerance = 2e-10,
                   info = paste(m, xi))
      want <- audit_native_fd(ref, b)
      got <- fns$d120(b, ld)
      expect_equal(unname((got-want)/(1+abs(want))), matrix(0, nrow(got), ncol(got)),
                   tolerance = 3e-5, info = paste(m, xi))
    }
  }
  ld <- audit_native_data(y)
  h <- colSums(egpd:::.degpd1idfns$d120(c(log(3), 1e-6, log(1.5)), ld))[7]
  expect_equal(h, 68.89984, tolerance = 2e-6)
})

test_that("identity support, zero and endpoint cells have consistent derivatives", {
  for (m in 1:6) {
    n <- c(3L,5L,3L,4L,3L,3L)[m]
    fns <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    b <- c(log(1.25), -0.4, log(1.5), if (m == 2) c(log(0.8), 0.3)
           else if (m == 4) log(1.7))
    y <- c(0,1,3)
    ld <- audit_native_data(y, n)
    ref <- function(z) audit_native_reference(z, y, m)
    expect_equal(fns$d0(b, ld), sum(ref(b)), tolerance = 2e-9)
    got <- fns$d120(b, ld); want <- audit_native_fd(ref, b, 0.0002)
    expect_equal(unname((got-want)/(1+abs(want))), matrix(0, nrow(got), ncol(got)),
                 tolerance = 5e-5, info = paste("endpoint carrier", m))
    expect_equal(fns$d0(b, audit_native_data(4, n)), 1e20)
    expect_equal(as.numeric(fns$d120(b, audit_native_data(4, n))), rep(0, n+n*(n+1)/2))
    b[1] <- log(0.2)
    expect_equal(fns$d0(b, audit_native_data(0, n)), 0, tolerance = 1e-12)
    expect_equal(as.numeric(fns$d120(b, audit_native_data(0, n))),
                 rep(0, n+n*(n+1)/2), tolerance = 1e-10)
  }
})

test_that("identity native predictors expand duplicate designs before offsets", {
  for (m in 1:6) {
    n <- c(3L,5L,3L,4L,3L,3L)[m]
    fns <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    b <- c(log(3), 0, rep(log(1.5), n-2))
    ld <- audit_native_data(c(0,1,7,15), n)
    ld$X <- rep(list(matrix(1, 2, 1)), n)
    ld$duplicate <- 1L; ld$dupid <- c(0L,1L,0L,1L)
    ld$offsets <- rep(list(c(0.1,-0.1,0.2,-0.2)), n)
    ld$offsets[[2]] <- c(0,1e-6,-1e-6,0)
    got <- fns$d120(b, ld)
    want <- do.call(rbind, lapply(1:4, function(i)
      fns$d120(b+vapply(ld$offsets, `[`, numeric(1), i),
                audit_native_data(ld$y[i,1], n))))
    expect_equal(got, want)
  }
})

test_that("identity high derivatives request a real fallback near zero", {
  ld <- audit_native_data(c(0,1,7,15))
  for (xi in c(0, -1e-6, 1e-6, -0.04999, 0.04999)) {
    b <- c(log(3), xi, log(1.5))
    expect_error(egpd:::.degpd1idfns$d340(b, ld), class = "egpd_d34_unavailable")
    expect_error(egpd:::degpd1id34(as.list(b), ld$X[[1]], ld$X[[2]], ld$X[[3]],
                                  ld$y[,1], 0L, 0L, ld$offsets, Inf),
                 "unstable near zero")
  }
})

test_that("identity model 1 retains finite log mass deep in exponential tails", {
  for (xi in c(0, -1e-6, 1e-6)) {
    b <- c(log(3), xi, 0)
    y <- c(100, 1000, 3000)
    ld <- audit_native_data(y)
    e <- function(t) if (xi == 0) -t/3 else -log1p(xi*t/3)/xi
    want <- -sum(e(y) + log(-expm1(e(y+1)-e(y))))
    expect_equal(egpd:::.degpd1idfns$d0(b, ld), want, tolerance = 1e-11)
    expect_true(all(is.finite(egpd:::.degpd1idfns$d120(b, ld))))
  }
})

test_that("all carriers retain near-zero finite-endpoint log mass and curvature", {
  for (m in 1:6) for (xi in c(-1e-4, -1e-6)) {
    n <- c(3L,5L,3L,4L,3L,3L)[m]
    fns <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    b <- c(log(-3.125*xi), xi, log(1.5), if (m == 2) c(log(0.8), 0.3)
           else if (m == 4) log(1.7))
    # At this endpoint exp(e_lo) underflows, so log mass equals e_lo plus
    # log G'(1) to far better than double precision (no CDF subtraction).
    ref <- function(z) {
      k <- exp(z[3])
      coeff <- switch(m, k,
        plogis(z[5])*k+(1-plogis(z[5]))*(k+exp(z[4])),
        (k+1)/k, (k+1)/k*exp(z[4])/2,
        sqrt(k)*dnorm(0)/(0.5-pnorm(-sqrt(k))),
        (15/32)*dbeta(0.5,k,k)/(pbeta(0.5,k,k)-pbeta(1/32,k,k)))
      log1p(z[2]*3/exp(z[1]))/z[2]-log(coeff)
    }
    ld <- audit_native_data(3, n)
    expect_equal(fns$d0(b, ld), ref(b), tolerance = 2e-12)
    got <- fns$d120(b, ld)
    expect_true(all(is.finite(got)))
    h <- abs(xi)*2e-5; bp <- bm <- b
    bp[2] <- bp[2]+h; bm[2] <- bm[2]-h
    g <- (ref(bp)-ref(bm))/(2*h)
    H <- (ref(bp)-2*ref(b)+ref(bm))/h^2
    # h_xixi is column 2n+1 in the packed layout.
    expect_equal(got[1,2]/g, 1, tolerance = 2e-6)
    expect_equal(got[1,2*n+1]/H, 1, tolerance = 2e-6)
  }
})

test_that("all identity carriers retain deep-tail finite probabilities", {
  for (m in 1:6) for (xi in c(0, -1e-6, 1e-6)) {
    n <- c(3L,5L,3L,4L,3L,3L)[m]
    fns <- get(paste0(".degpd", m, "idfns"), asNamespace("egpd"))
    b <- c(log(3), xi, log(1.5), if (m == 2) c(log(0.8), 0.3)
           else if (m == 4) log(1.7))
    y <- c(100,1000,3000)
    ld <- audit_native_data(y, n)
    expect_equal(fns$d0(b, ld), sum(audit_native_reference(b, y, m)),
                 tolerance = 2e-11)
    expect_true(all(is.finite(fns$d120(b, ld))))
  }
})

test_that("DEGPD ordered-mixture inversion delegates and recycles scalar parameters", {
  v <- c(0, 1e-14, 0.2, 1-1e-12, 1)
  k2 <- c(0.6, 0.8, 1.3, 1.5, 2)
  prob <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  want <- mapply(function(p, delta, prob)
    egpd:::q.G(p, type = 6, kappa = 0.4, delta = delta, prob = prob),
    v, k2, prob, USE.NAMES = FALSE)
  expect_equal(egpd:::.iG2_degpd(v, 0.4, k2, prob), want)
  expect_equal(egpd:::.iG2i_degpd(v[4], 0.4, k2[4], prob[4]), want[4])
  expect_true(all(is.finite(want)))
})
