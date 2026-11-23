## cdf(q, i, lower) returns log probabilities for the indicated rows.
.rqres_logtails <- function(y, cdf, support, uniforms) {
  i <- seq_along(y)
  lp <- cdf(y, i, TRUE)
  upper <- which(lp > -sqrt(.Machine$double.eps))
  if (length(upper)) lp[upper] <- cdf(y[upper], upper, FALSE)

  jump <- switch(support, count = which(is.finite(y) & y >= 0 & y == floor(y)),
                 mixed = which(y == 0), continuous = integer(0))
  if (length(jump)) {
    ## Either tail is linear in the same uniform jump weight. Use F(y-1)
    ## even for the smallest observed count; it need not be the support minimum.
    hi <- jump %in% upper
    prev <- rep(-Inf, length(jump))
    prev[hi] <- 0  # F(0-) = 0, S(0-) = 1
    if (support == "count") {
      for (lower in c(TRUE, FALSE)) {
        k <- which(y[jump] > 0 & hi == !lower)
        if (length(k)) prev[k] <- cdf(y[jump[k]] - 1, jump[k], lower)
      }
    }
    w <- if (missing(uniforms)) stats::runif(length(jump)) else uniforms[jump]
    lp[jump] <- .logspace_add(log1p(-w) + prev, log(w) + lp[jump])
  }
  out <- stats::qnorm(pmin(0, lp), log.p = TRUE)
  out[upper] <- -out[upper]
  out
}
