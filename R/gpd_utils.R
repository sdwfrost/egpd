## GPD utility functions needed for quantile predictions
## Adapted from evgam by Ben Youngman

.qgpd <- function(p, loc, scale, shape, zeta=1, theta=1, m=1) {
if (min(length(p), length(loc), length(scale), length(shape),
        length(zeta), length(theta), length(m)) == 0L)
  return(numeric(0))
template <- p
n <- max(length(p), length(loc), length(scale), length(shape),
         length(zeta), length(theta), length(m))
p <- rep_len(p, n); loc <- rep_len(loc, n); scale <- rep_len(scale, n)
shape <- rep_len(shape, n); zeta <- rep_len(zeta, n)
theta <- rep_len(theta, n); m <- rep_len(m, n)
v <- log(p) / (m * theta)
logsurv <- ifelse(v < -log(2), log1p(-exp(v)), log(-expm1(v))) - log(zeta)
L <- -logsurv
q <- scale * expm1(shape * L) / shape
zero <- which(shape == 0)
q[zero] <- scale[zero] * L[zero]
.egpd_restore_dim(loc + q, template)
}

.pgpd <- function(x, loc, scale, shape, tau, NAOK=FALSE, log=FALSE) {
if (min(length(x), length(loc), length(scale), length(shape), length(tau)) == 0L)
  return(numeric(0))
out <- tau + (1 - tau) * .pgpd_std(x - loc, scale, shape)
below <- rep_len(x < loc, length(out))
if (NAOK) {
  outside <- rep_len(1 + shape * (x - loc) / scale < 0, length(out))
  out[which(outside)] <- NaN
}
if (log) out <- log(out)
if (NAOK) out[below] <- NA
out
}

.dqgpd <- function(p, lscale, shape) {
if (min(length(p), length(lscale), length(shape)) == 0L)
  return(matrix(numeric(0), nrow = 0L, ncol = 2L, dimnames = list(NULL, c("d1", "d2"))))
n <- max(length(p), length(lscale), length(shape))
p <- rep_len(p, n); lscale <- rep_len(lscale, n); shape <- rep_len(shape, n)
L <- -log1p(-p)
d1 <- .qgpd_std(p, exp(lscale), shape)
d2 <- exp(lscale) * (shape * L * exp(shape * L) - expm1(shape * L)) / shape^2
small <- which(is.finite(L) & abs(shape * L) < 1e-4)
z <- shape[small] * L[small]
d2[small] <- exp(lscale[small]) * L[small]^2 *
  (1/2 + z/3 + z^2/8 + z^3/30 + z^4/144)
zero <- which(shape == 0)
d2[zero] <- exp(lscale[zero]) * L[zero]^2 / 2
bounded <- which(p == 1 & shape < 0)
d2[bounded] <- exp(lscale[bounded]) / shape[bounded]^2
d2[which(p == 1 & shape >= 0)] <- Inf
cbind(d1, d2)
}
