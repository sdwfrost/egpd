## Generic/matrix utility functions
## Adapted from evgam by Ben Youngman

.pivchol_rmvn <- function(n, mu, Sig) {
R <- suppressWarnings(chol(Sig, pivot = TRUE))
piv <- order(attr(R, "pivot"))
r <- attr(R, "rank")
V <- R[seq_len(r), piv, drop=FALSE]
Y <- crossprod(V, matrix(rnorm(n * r), nrow=r, ncol=n))
Y + as.vector(mu)
}

.precondition <- function(H) {
d <- 1/sqrt(abs(pmax(diag(H), 1e-8)))
out <- d * t(t(H) * d)
attr(out, "d") <- d
out
}

.perturb <- function(A) {
d <- attr(A, "d")
eps <- 1e-16
cholA <- suppressWarnings(chol(A, pivot=TRUE))
while(attr(cholA, "rank") < nrow(A)) {
diag(A) <- diag(A) + eps
cholA <- suppressWarnings(chol(A, pivot=TRUE))
eps <- 1e2 * eps
}
attr(A, "chol") <- cholA
return(A)
}

.solve_egpd <- function(H) {
H2 <- .precondition(H)
H2 <- .perturb(H2)
R <- attr(H2, "chol")
attr(R, "d") <- attr(H2, "d")
.precond_solve(R, diag(nrow(H)))
}

.precond_solve <- function(L, x) {
d <- attr(L, "d")
x <- d * as.matrix(x)
piv <- ipiv <- attr(L, "pivot")
ipiv[piv] <- seq_len(nrow(L))
d * backsolve(L, backsolve(L, x[piv, , drop=FALSE], upper.tri=TRUE, transpose=TRUE))[ipiv, , drop=FALSE]
}

.fdHess <- function(pars, gr, ..., eps=1e-4) {
np <- length(pars)
g0 <- gr(pars, ...)
out <- vapply(seq_len(np), function(i) gr(replace(pars, i, pars[i] + eps), ...) - g0, double(np))
out <- out + t(out)
out <- .5 * out / eps
out
}

.bdiag <- function(x) {
dims <- sapply(x, dim)
x <- x[dims[1,] > 0 | dims[2,] > 0]
dims <- sapply(x, dim)
row_ends <- cumsum(dims[1,])
row_starts <- c(1, row_ends + 1)
col_ends <- cumsum(dims[2,])
col_starts <- c(1, col_ends + 1)
out <- matrix(0, tail(row_ends, 1), tail(col_ends, 1))
for (i in seq_along(x))
  out[row_starts[i]:row_ends[i], col_starts[i]:col_ends[i]] <- x[[i]]
out
}

## smoothing matrix manipulation functions

.makeS <- function(lst, sp, spSl = Map("*", attr(lst, "Sl"), sp)) {
S0 <- attr(lst, "S0")
if (is.null(S0)) Reduce("+", spSl) else Reduce("+", spSl, init = S0)
}

## S(rho) = S0 + sum_j exp(rho_j) Sl[[j]], on the tempered likelihood.
## Accumulate tied penalties before fitting, but keep the original small
## smooth blocks and their map for the stable blockwise log determinant.
.joinSmooth <- function(lst, gamma = 1) {
  nb <- sum(vapply(lst, "[[", 0L, "nb"))
  nsp <- sum(vapply(lst, function(g) length(g$sp), 0L))
  mats <- vector("list", 1L + nsp)  # fixed, then one per independent log-sp
  out <- list()
  index <- integer(0)
  lsp0 <- numeric(0)
  b <- p <- 0L
  for (g in lst) {
    id <- g$sp.index
    id[id > 0L] <- id[id > 0L] + p
    off <- g$lsp0
    off[id == 0L] <- off[id == 0L] - log(gamma)
    index <- c(index, id)
    lsp0 <- c(lsp0, off)
    k <- 0L
    for (sm in g$smooth) {
      sm <- sm[c("S", "first.para", "last.para", "rank", "null.space.dim")]
      sm$first.para <- sm$first.para + b
      sm$last.para <- sm$last.para + b
      if (!length(sm$null.space.dim)) sm$null.space.dim <- 0
      out[[length(out) + 1L]] <- sm
      ii <- sm$first.para:sm$last.para
      for (S in sm$S) {
        k <- k + 1L
        j <- 1L + id[k]
        if (is.null(mats[[j]])) mats[[j]] <- matrix(0, nb, nb)
        mats[[j]][ii, ii] <- mats[[j]][ii, ii] + exp(off[k]) * S
      }
    }
    b <- b + g$nb
    p <- p + length(g$sp)
  }
  structure(out, nb = nb, Sl = mats[-1L], S0 = mats[[1L]],
            sp.index = index, lsp0 = unname(lsp0))
}

.full_lsp <- function(Sdata, rho) {
  attr(Sdata, "lsp0") + c(0, rho)[1L + attr(Sdata, "sp.index")]
}

## Read top-level sp on the unscaled likelihood. An expanded vector (e.g.
## from full.sp) is allowed, but cannot untie an id or override a fixed term.
.fixed_rho <- function(sp, Sdata, gamma) {
  if (!is.numeric(sp) || !is.null(dim(sp)) || any(!is.finite(sp) | sp < 0))
    stop("`sp` must be a finite non-negative numeric vector.", call. = FALSE)
  n <- length(attr(Sdata, "Sl"))
  id <- attr(Sdata, "sp.index")
  rho <- log(as.numeric(sp)) - log(gamma)
  if (length(sp) == n) return(rho)
  if (length(sp) == 1L && n > 0L) return(rep(rho, n))
  if (length(sp) != length(id))
    stop(sprintf("Expected %d independent smoothing parameters or %d penalty multipliers in `sp`.",
                 n, length(id)), call. = FALSE)
  first <- match(seq_len(n), id)
  pars <- rho[first] - attr(Sdata, "lsp0")[first]
  full <- .full_lsp(Sdata, pars)
  if (!all(full == rho | abs(full - rho) <= 1e-12))
    stop("`sp` conflicts with a term-level fixed or shared smoothing parameter.", call. = FALSE)
  pars
}

.grad.R <- function(pars, Sdata, R, eps, likfns, likdata, H0=NULL) {
S <- .makeS(Sdata, exp(pars))
if (is.null(H0)) H0 <- .hess.nopen(likdata$beta, likdata, likfns)
H <- H0 + S
H <- H[likdata$LAid, likdata$LAid]
cholH <- try(chol(H), silent=TRUE)
if (inherits(cholH, "try-error")) {
iH <- pinv(H)
} else {
iH <- chol2inv(cholH)
}
choliH <- try(chol(iH), silent=TRUE)
if (inherits(choliH, "try-error")) choliH <- attr(.perturb(iH), "chol")
(choliH - R) / eps
}

## other stuff

.aux_response_name <- function(name, used) {
  while (name %in% used) name <- paste0(name, ".")
  name
}

.padWithVec <- function(df, vec, nm) {
df <- df[1,]
id <- cbind(1, seq_along(vec))
df <- df[id[,1],]
df[,nm] <- vec
df
}

.addNodes <- function(df1, df2, id) {
df1 <- lapply(df1, rep, nrow(df2))
for (i in seq_along(id)) df1[[id[i]]] <- df2[,i]
as.data.frame(df1)
}

.outer.list <- function(lst) {
out <- lst[[1]]
if (length(lst) > 1) {
for (i in 2:length(lst)) out <- out %o% lst[[i]]
}
c(out)
}

.updateControl <- function(lst0, lst) {
for (i in seq_along(lst0)) {
nmi <- names(lst0)[i]
i2 <- which(names(lst) == nmi)
lsti <- lst[[i2]]
for (j in seq_along(lst0[[i]])) {
nmj <- names(lst0[[i]])[j]
j2 <- which(names(lsti) == nmj)
lsti[[j2]] <- lst0[[i]][[j]]
}
lst[[i2]] <- lsti
}
lst
}
