## H is the GPD cdf. The median of a power component solves H(m)^k = 1/2.
## Pretend the bottom p and top 1-p of the data came from different components
## and match their medians. These hard allocations only propose starts: score
## the *mixture*, not its classified likelihood, and leave all fitting to egpd().
.m2_restart_starts <- function(y, cont, zi, idlink, xm, fit) {
  np <- 5L + zi
  offset <- numeric(np)
  pinned <- rep(FALSE, np)
  held <- numeric(np)
  if (!is.null(fit)) {
    ld <- fit$likdata
    y <- ld$y[, 1L]  # after omission of incomplete rows and maxdata subsampling
    for (i in seq_len(np)) {
      if (length(ld$offsets[[i]])) offset[i] <- mean(ld$offsets[[i]])
      pinned[i] <- !ncol(ld$X[[i]]) || .is_pinned_block(ld$X[[i]])
      if (pinned[i]) held[i] <- .egpd_invlink(fit, names(fit)[i], offset[i])
    }
  }
  value <- function(i, proposed) if (pinned[i]) held[i] else proposed
  y <- sort(y, na.last = NA)
  if (length(y) < 6L || all(pinned[3:5]) || any(!is.finite(y)) ||
      y[1L] < 0 || (cont && y[1L] == 0) || y[1L] == y[length(y)])
    return(list())
  ## Equally weighted central ranks of equal-sized strata; no random subsample
  ## and no deliberately isolating the sample maximum as a component.
  n <- min(512L, length(y))
  y <- y[ceiling((seq_len(n) - 0.5) * length(y) / n)]
  pi0 <- if (zi) value(6L, max(0.01, mean(y == 0) / 2)) else 0
  tail <- 2^-(1:6)
  split <- sort(unique(c(tail, 1 - tail)))
  split <- split[pmin(split, 1 - split) * (1 - pi0) * n >= 3]
  if (!length(split)) return(list())
  ## For counts only the initial medians use midpoints, not the likelihood.
  mid <- if (cont) y else y + 0.5
  q <- function(p) stats::quantile(mid, pi0 + (1 - pi0) * p, names = FALSE)
  m1 <- q(split / 2)
  m2 <- q((1 + split) / 2)
  scales <- value(1L, q(0.5) * c(0.25, 0.5, 0.75, 1, 1.5, 2))
  shapes <- value(2L, unique(pmin(c(0.025, 0.2, 0.6), 0.5 * xm)))
  ## One representative per weight range. Ranking all starts together tends to
  ## spend the entire budget in one basin, especially near a one-power fit.
  band <- 1L + (split > 1/2) + (split > 15/16)
  best <- vector("list", 3L)
  score <- rep(-Inf, length(best))
  for (xi in shapes) for (s in scales) {
    if (!is.finite(s) || s <= 0 || !is.finite(xi)) next
    logH <- function(x) .log1mexp(.gpd_cdf_surv_std(x, s, xi)$log_surv)
    a <- log(0.5) / logH(m1)
    b <- log(0.5) / logH(m2)
    if (cont) {
      lo <- logH(y)
      logh <- .dgpd_std(y, s, xi, log = TRUE)
      component <- function(k) log(k) + (k - 1) * lo + logh
    } else {
      hi <- logH(y + 1)
      log_ratio <- logH(y) - hi
      component <- function(k) k * hi + .log1mexp(k * log_ratio)
    }
    for (j in seq_along(split)) {
      k1 <- value(3L, a[j])
      dk <- value(4L, b[j] - k1)
      ## A bound on *proposals*, not on estimation. Don't seed a spike at the
      ## endpoint or a pair already collapsed to the same numerical power.
      k2 <- k1 + dk
      if (!is.finite(k2) || k1 <= 0 || k2 <= k1 || k2 > 1e6) next
      p <- value(5L, split[j])
      ll <- .logspace_add(log(p) + component(k1), log1p(-p) + component(k2))
      if (zi) {
        ll <- log1p(-pi0) + ll
        ll[y == 0] <- .logspace_add(log(pi0), ll[y == 0])
      }
      total <- sum(ll)
      g <- band[j]
      if (is.finite(total) && total > score[g]) {
        eta <- c(log(s), .shape_eta(xi, idlink, xm), log(k1), log(dk),
                 stats::qlogis(p), if (zi) stats::qlogis(pi0))
        best[[g]] <- .eta_start(replace(eta, pinned, 0))
        score[g] <- total
      }
    }
  }
  unique(Filter(Negate(is.null), best))
}
