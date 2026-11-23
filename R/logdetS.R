## Functions for calculating log(|S|)
## Adapted from evgam by Ben Youngman

.logdetS <- function(smth, rho, deriv=0) {
map <- attr(smth, "sp.index")
nsp <- length(rho)
rho <- .full_lsp(smth, rho)
id <- rep(seq_along(smth), sapply(smth, function(x) length(x$S)))
## A top-level sp of exactly zero removes a penalty from both |S|+ and
## its structural rank, including when other penalties share its block.
keep <- rho != -Inf
if (!all(keep)) {
  masks <- split(keep, factor(id, levels = seq_along(smth)))
  smth <- Map(function(s, k) {
    s$S <- s$S[k]
    s$rank <- s$rank[k]
    s
  }, smth, masks)
  map <- map[keep]
  id <- id[keep]
  rho <- rho[keep]
}
out <- if (deriv == 0) .logdetS_const(smth) else list()
rho <- split(rho, id)
## split() drops the smooths without a penalty (fx = TRUE), so take each group's smooth
## from the group's name: by position, an unpenalised smooth before a penalised one
## gave the penalised one's log|S| the other smooth's rank.
smth <- smth[as.integer(names(rho))]
ldSi <- list()
for (i in seq_along(rho)) {
  if (length(rho[[i]]) == 1) {
    ldSi[[i]] <- .logdetS_single(rho[[i]], smth[[i]], deriv)
  } else
    ldSi[[i]] <- .logdetS_multi(rho[[i]], smth[[i]], deriv)
}
out$d0 <- Reduce("+", lapply(ldSi, function(x) x$d0), init = 0)
if (deriv > 0) {
  out$d1 <- as.numeric(unlist(lapply(ldSi, function(x) x$d1)))
  if (deriv > 1)
    out$d2 <- if (length(ldSi)) .bdiag(lapply(ldSi, function(x) x$d2)) else matrix(0, 0, 0)
  if (!identical(map, seq_len(nsp))) {
    L <- outer(map, seq_len(nsp), "==")
    out$d1 <- drop(crossprod(L, out$d1))
    if (deriv > 1) out$d2 <- crossprod(L, out$d2 %*% L)
  }
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

## For the constants of the REML criterion: the part of log|S|+ that .logdetS leaves
## out, `ldev`, and the rank of the total penalty, `rank`. .logdetS gives a
## single-penalty smooth's log|lambda S|+ as rank * log(lambda), without the sum of the
## logs of S's own nonzero eigenvalues; a multi-penalty smooth's log|S|+ is complete
## there, and its rank is that of the sum of its normalised penalties (its structural
## rank, as in .logdetS_multi).
.logdetS_const <- function(smth) {
ldev <- 0
rank <- 0
for (smthi in smth) {
  if (length(smthi$S) == 1) {
    ev <- eigen(smthi$S[[1]], symmetric=TRUE, only.values=TRUE)$values
    ldev <- ldev + sum(log(ev[seq_len(smthi$rank)]))
    rank <- rank + smthi$rank
  } else if (length(smthi$S) > 1) {
    rank <- rank + .rank(Reduce("+", lapply(smthi$S, function(x) x / norm(x, "F"))))
  }
}
list(ldev=ldev, rank=rank)
}

.logdetS_multi <- function(rho, smthi, deriv) {

Sl <- smthi$S

lambda <- exp(rho)
q <- nrow(Sl[[1]])
eps <- .Machine$double.eps^(1/3)

S <- Reduce("+", Map("*", Sl, lambda))

## Restrict the penalties to their structural range space. That space does not depend
## on lambda: it is the range of the sum of the normalised penalties (as in mgcv's
## totalPenaltySpace). The lambda-weighted sum has the same range in exact arithmetic,
## and its eigenvectors are used when they resolve all of it, but with widely spread
## smoothing parameters (log-sp spreads of about 25 or more) the directions penalised
## only by the small-lambda terms fall below its rank threshold; they were dropped as
## null space, which removed their contribution to log|S|+ and its derivatives.
if (.rank(S) < q) {
  tS <- Sl
  for (i in seq_along(tS))
    tS[[i]] <- tS[[i]] / norm(tS[[i]], "F")
  S <- Reduce("+", Map("*", tS, lambda))
  eS <- eigen(S, symmetric=TRUE)
  okay <- eS$values > max(eS$values) * .Machine$double.eps^.8
  eS0 <- eigen(Reduce("+", tS), symmetric=TRUE)
  okay0 <- eS0$values > max(eS0$values) * .Machine$double.eps^.8
  if (sum(okay) < sum(okay0)) {
    eS <- eS0
    okay <- okay0
  }
  if (!all(okay)) {
    Up <- eS$vectors[, okay, drop=FALSE]
    Sl <- lapply(Sl, function(x) crossprod(Up, crossprod(x, Up)))
    q <- nrow(Sl[[1]])
  }
}

S <- Reduce("+", Map("*", Sl, lambda))
M <- length(Sl)

## Similarity transforms that separate the dominant penalties from the rest (Wood,
## 2011, Appendix B). The rank and range of each dominant group come from the sum of
## its normalised penalties, which does not depend on how far apart the group's
## smoothing parameters are; and S is assembled from the transformed penalties, in
## which the dominant terms are confined to their own range, rather than by
## transforming the lambda-weighted sum, where their rounding error swamps the
## dominated terms.
Q <- q
K <- 0
gamma <- seq_len(M)
Sb <- Sl
cond <- TRUE
transformed <- FALSE

while (cond) {
  norms <- sapply(Sb[gamma], norm, "F")
  Omega <- try(norms * lambda[gamma])
  if (inherits(Omega, "try-error")) stop("Error computing Omega in logdetS_multi")
  dominant <- Omega >= eps * max(Omega)
  alpha <- gamma[dominant]
  gammap <- gamma[!dominant]
  eS <- eigen(Reduce("+", Map("/", Sb[alpha], norms[dominant])), symmetric=TRUE)
  r <- sum(eS$values > max(eS$values) * .Machine$double.eps^.8)
  if (r != Q) {
    U <- eS$vectors
    Ur <- U[,seq_len(r), drop=FALSE]
    Un <- U[,-seq_len(r), drop=FALSE]
    Tg <- .bdiag(list(diag(K), U))
    ## zero columns for the directions outside the dominant terms' range (K + Q rows:
    ## `0 * Un` had only Q, so a second transformation, K > 0, failed)
    Ta <- cbind(.bdiag(list(diag(K), Ur)), matrix(0, K + Q, ncol(Un)))
    for (i in alpha)
      Sl[[i]] <- crossprod(Ta, crossprod(Sl[[i]], Ta))
    for (i in gammap)
      Sl[[i]] <- crossprod(Tg, crossprod(Sl[[i]], Tg))
    for (i in gammap)
      Sb[[i]] <- crossprod(Un, crossprod(Sb[[i]], Un))
    K <- K + r
    Q <- Q - r
    gamma <- gammap
    transformed <- TRUE
  } else {
    cond <- FALSE
  }
}

if (transformed)
  S <- Reduce("+", Map("*", Sl, lambda))
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
