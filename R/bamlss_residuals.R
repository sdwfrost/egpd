## residuals.bamlss supplies *back-transformed* parameter columns. With
## posterior summaries these can also be observation-by-draw matrices.
.bamlss_residual_rows <- function(par, y = NULL) {
  par <- lapply(par, function(v) if (is.data.frame(v)) as.matrix(v) else v)
  cols <- c(par, if (!is.null(y)) list(y))
  longest <- which(lengths(cols) == max(lengths(cols)))
  mat <- longest[vapply(cols[longest], is.matrix, logical(1))]
  template <- cols[[if (length(mat)) mat[1L] else longest[1L]]]
  n <- if (any(lengths(cols) == 0L)) 0L else length(template)
  list(par = lapply(par, function(v) rep_len(as.numeric(v), n)),
       y = if (!is.null(y)) rep_len(as.numeric(y), n), template = template)
}

.bamlss_rqres <- function(y, par, pfun, support) {
  rows <- .bamlss_residual_rows(par, y)
  y <- rows$y
  cdf <- function(q, i, lower) {
    pfun(q, lapply(rows$par, `[`, i), lower.tail = lower, log.p = TRUE)
  }
  .egpd_restore_dim(.rqres_logtails(y, cdf, support), rows$template)
}

.bamlss_egpd_hazard <- function(x, p) {
  t <- -.gpd_cdf_surv_std(x, p$sigma, p$xi)$log_surv
  if (p$xi < 0) t[x >= -p$sigma/p$xi] <- Inf
  t[x == 0] <- 0
  t
}

.bamlss_egpd_log_surv <- function(t, p, type) {
  n <- length(t)
  v <- .egpd_carrier_log(t, type, rep_len(p$prob, n), rep_len(p$kappa, n),
                         rep_len(p$delta, n))$lS
  ## A log survivor of -Inf is a zero contribution only off support, not
  ## an acceptable quadrature sample at a finite hazard.
  if (anyNA(v) || any(v > 1e-12) || any(v == -Inf & is.finite(t)))
    stop("EGPD mean: could not resolve the carrier survivor.", call. = FALSE)
  v
}

## log(I(from)/sigma), where I is the remaining survivor integral and
## 'from' is in response units. In GPD hazard t, dx = sigma*exp(xi*t) dt.
.bamlss_egpd_log_tail_per_scale <- function(from, p, type) {
  xi <- p$xi
  A <- .bamlss_egpd_hazard(from, p)
  if (A == Inf) return(-Inf)
  k <- p$kappa; d <- p$delta
  rate <- 1 - xi
  if (type == 1L && k == 1) return(-rate*A - log1p(-xi))
  ls <- function(t) .bamlss_egpd_log_surv(t, p, type)
  log_integral <- function(logf, lo, hi) {
    shift <- logf((lo + hi)/2)
    val <- tryCatch(stats::integrate(function(t) exp(logf(t) - shift), lo, hi,
                                     rel.tol = 1e-9, abs.tol = 0)$value,
                    error = function(e) stop("EGPD mean: survivor integration failed: ",
                                             conditionMessage(e), call. = FALSE))
    if (!is.finite(val) || val <= 0)
      stop("EGPD mean: could not resolve the survivor integral.", call. = FALSE)
    shift + log(val)
  }
  if (xi <= -1) {
    ## u = -xi*(t-A) is minus the log fraction of support remaining.
    ## Unlike a linear endpoint coordinate, it resolves the survivor's
    ## narrow transition there without losing very short supports. Since
    ## S decreases, the omitted integral of exp(-u)*S beyond 40 is at
    ## most 1/expm1(40) times the integral retained.
    bulk <- log_integral(function(u) ls(A + u/(-xi)) - u, 0, 40)
    return(-log(-xi) + xi*A + bulk)
  }
  ## Beyond B the carrier survivor is linear in the GPD survivor to
  ## roundoff (even after multiplication by a power-carrier exponent).
  ## For the beta-power carrier it is linear in S4, not necessarily in
  ## exp(-t): delta*t can still be small, so integrate S4's tail exactly.
  start <- 48 + log(max(1, k, if (type == 6L) d, na.rm = TRUE))
  B <- max(A, start)
  if (type %in% c(4L, 5L)) {
    ## S4 = exp(-t) * (1 + (1-exp(-delta*t))/delta).
    logratio <- log1p(B * .expm1_over_x(d*B) + exp(-d*B)/(rate + d))
    if (type == 5L) logratio <- logratio + log(k) - log(2)
  } else {
    ## Evaluate S(B)/exp(-B) at the onset of the linear tail, not at a
    ## possibly enormous A where the addition would cancel its logarithm.
    logratio <- ls(start) + start
  }
  result <- logratio - rate*B - log(rate)
  if (A < B) {
    ## Bound panel lengths so centering the log integrand also works for
    ## large powers (B can exceed 700) and shapes close to one.
    cuts <- seq(A, B, length.out = 1 + ceiling((B-A)/32))
    logf <- function(t) ls(t) + xi*t
    for (j in seq_len(length(cuts) - 1L))
      result <- .logspace_add(result, log_integral(logf, cuts[j], cuts[j+1L]))
  }
  result
}

.bamlss_egpd_log_mean_per_scale <- function(p, type, discrete) {
  if (!discrete) return(.bamlss_egpd_log_tail_per_scale(0, p, type))
  ## E[floor X] = sum_{j>=1} S_X(j), not floor(E[X]). For a decreasing
  ## survivor the remainder after N lies in [I(N)-S(N), I(N)], where I
  ## is the tail integral. Bound the error rather than silently truncating
  ## a heavy tail at an arbitrary largest count.
  ## Accumulate log(E[floor X]/sigma), *not* log(E[floor(X/sigma)]).
  logscale <- log(p$sigma)
  total <- lower <- -Inf
  lo <- 1L
  for (n in 2^(5:19)) {
    s <- .bamlss_egpd_log_surv(.bamlss_egpd_hazard(lo:n, p), p, type)
    if (s[1L] == -Inf) return(total)
    total <- .logspace_add(total, s[1L] - logscale + log(sum(exp(s - s[1L]))))
    width <- s[length(s)] - logscale
    if (width == -Inf) return(total)
    ## Include the known tail in the lower bound for broad distributions:
    ## even S(N) near one can be negligible relative to the response mean.
    if (lo == 1L || width <= log(1e-7) + max(total, lower)) {
      tail <- .bamlss_egpd_log_tail_per_scale(n, p, type)
      upper <- .logspace_add(total, tail)
      ## Bracket width min(I(N), S(N)) as a fraction of its upper end.
      frac <- exp(min(tail, width) - upper)
      lower <- upper + log1p(-frac)
      if (frac <= 1e-7 * (1 - frac))
        return(upper + log1p(-frac/2))
    }
    lo <- n + 1L
  }
  stop("DEGPD mean: could not bound the survivor-sum error to relative tolerance 1e-7.",
       call. = FALSE)
}

.bamlss_egpd_mean <- function(par, type, discrete) {
  n <- length(par$sigma)
  for (nm in c("kappa", "delta", "prob"))
    if (is.null(par[[nm]])) par[[nm]] <- rep(NA_real_, n)
  pi <- if (is.null(par$pi)) rep(0, n) else par$pi
  par$pi <- NULL
  valid <- .egpd_carrier_valid(type, par$prob, par$kappa, par$delta, n) &
    is.finite(par$sigma) & par$sigma > 0 & is.finite(par$xi) &
    is.finite(pi) & pi >= 0 & pi <= 1
  out <- rep(NA_real_, n)
  out[valid & pi == 1] <- 0
  infinite <- which(valid & pi < 1 & par$xi >= 1)
  if (length(infinite)) {
    warning("EGPD mean is infinite for xi >= 1; response residual undefined.", call. = FALSE)
    out[infinite] <- Inf
  }
  ok <- which(valid & pi < 1 & par$xi < 1)
  for (g in .count_groups(lapply(par, `[`, ok))) {
    i <- ok[g]
    logm <- .bamlss_egpd_log_mean_per_scale(lapply(par, `[`, i[1L]), type, discrete)
    ## Combine the weight with the standardized moment before restoring
    ## scale: the unweighted response mean may not fit in a double.
    out[i] <- exp(log(par$sigma[i]) + (logm + log1p(-pi[i])))
  }
  out
}

.bamlss_residual_family <- function(family, support, dtype = NULL) {
  pfun <- family$p
  family$rqres <- function(y, par, ...) .bamlss_rqres(y, par, pfun, support)
  family$mu <- function(par, ...) {
    rows <- .bamlss_residual_rows(par)
    p <- rows$par
    m <- if (!is.null(dtype)) .bamlss_egpd_mean(p, dtype, support == "count") else
      if (is.null(p$pi)) p$mu else (1 - p$pi) * p$mu
    .egpd_restore_dim(m, rows$template)
  }
  class(family) <- "family.bamlss"
  family
}
