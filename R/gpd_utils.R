## GPD utility functions needed for quantile predictions
## Adapted from evgam by Ben Youngman

.qgpd <- function(p, loc, scale, shape, zeta=1, theta=1, m=1) {
shape <- sign(shape) * pmax(abs(shape), 1e-6)
out <- 1 - p ^ (1 / (m * theta))
loc + scale * ((out / zeta)^(-shape) - 1) / shape
}

.pgpd <- function(x, loc, scale, shape, tau, NAOK=FALSE, log=FALSE) {
below <- x < loc
x <- pmax(x, loc)
temp <- 1 + shape * (x - loc) / scale
if (!NAOK) temp <- pmax(temp, 0)
out <- temp ^ (-1/shape)
out <- 1 - (1 - tau) * out
if (log) out <- log(out)
if (NAOK) out[below] <- NA
out
}

.dqgpd <- function(p, lscale, shape) {
shape <- sign(shape) * pmax(abs(shape), 1e-6)
.e1 <- 1 - p
.e2 <- .e1^shape
.e4 <- 1/.e2 - 1
.e5 <- exp(lscale)
d1 <- .e4 * .e5/shape
d2 <- -((.e4/shape + log(.e1)/.e2) * .e5/shape)
cbind(d1, d2)
}
