## Functions for calculating log(|S|)
## Adapted from evgam by Ben Youngman

.logdetS <- function(smth, rho, deriv=0) {
id <- rep(seq_along(smth), sapply(smth, function(x) length(x$S)))
rho <- split(rho, id)
out <- ldSi <- list()
for (i in seq_along(rho)) {
  if (length(rho[[i]]) == 1) {
    ldSi[[i]] <- .logdetS_single(rho[[i]], smth[[i]], deriv)
  } else
    ldSi[[i]] <- .logdetS_multi(rho[[i]], smth[[i]], deriv)
}
out$d0 <- Reduce("+", lapply(ldSi, function(x) x$d0))
if (deriv > 0) {
  out$d1 <- unlist(lapply(ldSi, function(x) x$d1))
  if (deriv > 1)
    out$d2 <- .bdiag(lapply(ldSi, function(x) x$d2))
}
out
}

.logdetS_single <- function(rho, smthi, deriv) {
out <- list(d0 = rho * smthi$rank)
if (deriv > 0) {
  out$d1 <- smthi$rank
  if (deriv > 1)
    out$d2 <- matrix(0, length(rho), length(rho))
}
out
}

.rank <- function(x) {
ev <- eigen(x, symmetric=TRUE, only.values=TRUE)
sum(ev$values > max(ev$values) * .Machine$double.eps^.8)
}

.logdetS_multi <- function(rho, smthi, deriv) {

Sl <- smthi$S

lambda <- exp(rho)
q <- nrow(Sl[[1]])
eps <- .Machine$double.eps^(1/3)

S <- Reduce("+", Map("*", Sl, lambda))

if (.rank(S) < q) {
  tS <- Sl
  for (i in seq_along(tS))
    tS[[i]] <- tS[[i]] / norm(tS[[i]], "F")
  S <- Reduce("+", Map("*", tS, lambda))
  eS <- eigen(S, symmetric=TRUE)
  okay <- eS$values > max(eS$values) * .Machine$double.eps^.8
  Up <- eS$vectors[, okay, drop=FALSE]
  Sl <- lapply(Sl, function(x) crossprod(Up, crossprod(x, Up)))
  q <- nrow(Sl[[1]])
}

S <- Reduce("+", Map("*", Sl, lambda))
M <- length(Sl)

Q <- q
K <- 0
gamma <- seq_len(M)
Sb <- Sl
cond <- TRUE

while (cond) {
  norms <- sapply(Sb[gamma], norm, "F")
  Omega <- try(norms * lambda[gamma])
  if (inherits(Omega, "try-error")) stop("Error computing Omega in logdetS_multi")
  alpha <- gamma[Omega >= eps * max(Omega)]
  gammap <- gamma[Omega < eps * max(Omega)]
  Sa <- Reduce("+", Map("*", Sb[alpha], lambda[alpha]))
  r <- .rank(Sa)
  if (r != Q) {
    eS <- eigen(Sa, symmetric=TRUE)
    U <- eS$vectors
    Ur <- U[,seq_len(r), drop=FALSE]
    Un <- U[,-seq_len(r), drop=FALSE]
    Tg <- .bdiag(list(diag(K), U))
    Sp <- try(crossprod(Tg, crossprod(S, Tg)), silent= TRUE)
    Ta <- cbind(.bdiag(list(diag(K), Ur)), 0 * Un)
    for (i in alpha)
      Sl[[i]] <- crossprod(Ta, crossprod(Sl[[i]], Ta))
    for (i in gammap)
      Sl[[i]] <- crossprod(Tg, crossprod(Sl[[i]], Tg))
    for (i in gammap)
      Sb[[i]] <- crossprod(Un, crossprod(Sb[[i]], Un))
    K <- K + r
    Q <- Q - r
    S <- Sp
    gamma <- gammap
  } else {
    cond <- FALSE
  }
}

cS <- chol(S)
out <- list(d0 = 2 * sum(log(diag(cS))))

if (deriv > 0) {
  iS <- chol2inv(cS)
  if (deriv > 1) {
    Sl <- lapply(Sl, function(x) crossprod(iS, x))
    out$d1 <- lambda * sapply(Sl, function(x) sum(diag(x)))
    out$d2 <- diag(out$d1, M)
    for (i in seq_len(M)) {
      for (j in seq_len(i)) {
        out$d2[i, j] <- out$d2[i, j] - lambda[i] * lambda[j] * sum(t(Sl[[i]]) * Sl[[j]])
        if (i != j)
          out$d2[j, i] <- out$d2[i, j]
      }
    }
  } else {
    out$d1 <- lambda * sapply(Sl, function(x) sum(iS * x))
  }
}

out

}
