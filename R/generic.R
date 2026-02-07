## Generic/matrix utility functions
## Adapted from evgam by Ben Youngman

.pivchol_rmvn <- function(n, mu, Sig) {
R <- suppressWarnings(chol(Sig, pivot = TRUE))
piv <- order(attr(R, "pivot"))
r <- attr(R, "rank")
V <- R[1:r, piv]
Y <- crossprod(V, matrix(rnorm(n * r), r))
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
D <- attr(H2, "d")
D <- diag(D, length(D))
R <- attr(H2, "chol")
crossprod(backsolve(R, D, transpose=TRUE))
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

.makeS <- function(lst, sp) {
Reduce("+", Map("*", attr(lst, "Sl"), sp))
}

.joinSmooth <- function(lst) {
nms <- c("S", "first.para", "last.para", "rank", "null.space.dim")
nbi <- sapply(lst, function(x) x$nb)
starts <- cumsum(c(0, nbi))
if (length(lst) == 1) {
  lst <- lst[[1]]$smooth
} else {
  lst <- lapply(lst, function(x) x$smooth)
  smooths <- sapply(lst, length) > 0
  for (i in which(smooths)) {
    for (j in seq_along(lst[[i]])) {
      lst[[i]][[j]]$first.para <- lst[[i]][[j]]$first.para + starts[i]
      lst[[i]][[j]]$last.para <- lst[[i]][[j]]$last.para + starts[i]
    }
  }
lst <- unlist(lst, recursive=FALSE)
}
out <- lapply(lst, function(x) subset(x, names(x) %in% nms))
lst2 <- list()
nb <- tail(starts, 1)
for (i in seq_along(lst)) {
  temp <- list()
  ind <- lst[[i]]$first.para:lst[[i]]$last.para
  for (j in seq_along(lst[[i]]$S)) {
    temp2 <- matrix(0, nb, nb)
      temp2[ind, ind] <- lst[[i]]$S[[j]]
    temp[[j]] <- temp2
  }
  lst2[[i]] <- temp
}
lst2 <- unlist(lst2, recursive=FALSE)
for (i in seq_along(out)) {
  if (length(out[[i]]$null.space.dim) == 0)
    out[[i]]$null.space.dim <- 0
}
attr(out, "Sl") <- lst2
attr(out, "nb") <- nb
out
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
