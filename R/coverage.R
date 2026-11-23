## Predictive coverage for egpd, gamlss, and bamlss fits
## -------------------------------------------------------
## S3 generic + methods for computing simulation-based prediction
## intervals and their empirical coverage.

## ---------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------

## Simulate responses from an EGPD-family distribution.
## `par_list` is a named list; each element is either an nsim x n
## matrix (one column per observation, one row per simulation draw)
## or a length-n vector (plug-in case).
## `rfun` is one of regpd / rdiscegpd / rzidiscegpd / rziegpd,
## `dtype` is the integer type code (1-6), and `nsim` is the
## number of draws per observation.
## Returns an nsim x n matrix.
.sim_responses <- function(par_list, rfun, dtype, nsim) {
  for (nm in names(par_list)) {
    v <- par_list[[nm]]
    if (any(!is.finite(v)))
      stop("Non-finite parameter draws for '", nm, "'.", call. = FALSE)
    if (nm %in% c("sigma", "kappa", "delta") && any(v <= 0))
      stop("Parameter draws for '", nm, "' must be positive.", call. = FALSE)
    if (nm %in% c("prob", "pi") && any(v < 0 | v > 1))
      stop("Parameter draws for '", nm, "' must be in [0, 1].", call. = FALSE)
  }
  first <- par_list[[1]]
  n <- if (is.matrix(first)) ncol(first) else length(first)
  ysim <- matrix(NA_real_, nsim, n)
  for (j in seq_len(n)) {
    args <- list(n = nsim, type = dtype)
    for (nm in names(par_list)) {
      v <- par_list[[nm]]
      args[[nm]] <- if (is.matrix(v)) v[, j] else v[j]
    }
    ysim[, j] <- do.call(rfun, args)
  }
  ysim
}

## Compute prediction intervals from an nsim x n simulation matrix.
## Returns list(L, U, covered) given observed y and level.
.pi_from_sims <- function(ysim, y, level, use_shortest, is_discrete) {
  n <- length(y)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level < 0 || level > 1)
    stop("'level' must be in [0, 1].")
  if (any(!is.finite(ysim)))
    stop("Non-finite simulated responses; prediction intervals are undefined.")
  alpha <- 1 - level
  if (!use_shortest) {
    qtype <- if (is_discrete) 1 else 7
    L <- apply(ysim, 2, quantile, probs = alpha / 2, type = qtype,
               names = FALSE)
    U <- apply(ysim, 2, quantile, probs = 1 - alpha / 2, type = qtype,
               names = FALSE)
  } else {
    N <- nrow(ysim)
    if (!N) stop("Cannot construct an interval without simulated responses.")
    ## Compare empirical fractions, not accumulated masses or a rounded-up
    ## product (e.g. .14 * 100). At zero mass select the lowest singleton.
    h <- max(1, floor(level * N))
    if (h / N < level) h <- h + 1
    left <- seq_len(N - h + 1)
    right <- left + (h - 1)
    L <- U <- numeric(n)
    for (j in seq_len(n)) {
      x <- sort(ysim[, j])
      width <- x[right] - x[left]
      if (all(is.infinite(width))) width <- x[right]/2 - x[left]/2
      i <- which.min(width)
      L[j] <- x[left[i]]; U[j] <- x[right[i]]
    }
  }
  if (is_discrete) { L <- ceiling(L); U <- floor(U) }
  covered <- (y >= L) & (y <= U)
  list(L = L, U = U, covered = covered)
}

## ---------------------------------------------------------------
## S3 generic
## ---------------------------------------------------------------

#' Predictive coverage for EGPD models
#'
#' Compute simulation-based prediction intervals and their empirical
#' coverage for fitted EGPD-family models.  Methods are provided for
#' objects of class \code{"egpd"}, \code{"gamlss"}, and
#' \code{"bamlss"}.
#'
#' The algorithm has three steps:
#' \enumerate{
#'   \item \strong{Parameter draws.}
#'     \describe{
#'       \item{\code{method = "plug-in"}}{Fitted parameters are held
#'         fixed at their point estimates.}
#'       \item{\code{method = "parametric"}}{Each simulated response uses its
#'         own parameter draw.  For \code{egpd} fits this is a
#'         multivariate-normal draw of the coefficients from their joint
#'         posterior (\code{\link{simulate.egpd}}, using \code{Vc}); for
#'         \code{gamlss} fits, independent normal perturbations of each
#'         linear predictor scaled by its standard error, back-transformed.
#'         In sample these are the standard errors of
#'         \code{predict(fit, se.fit = TRUE)}; \pkg{gamlss} gives none at
#'         \code{newdata}. There the covariance is the inverse observed
#'         information of the free coefficients: the fitted, case-weighted log
#'         likelihood is differentiated in their joint predictor directions
#'         at steps calibrated to changes in log density, keeping fixed
#'         parameters on the response scale and grouped binomial trial sizes
#'         at their training values. The appropriate
#'         covariance block and new design row give each marginal standard error;
#'         predictor draws are still independent, not joint coefficient draws.
#'         This requires saved parametric designs, a density that reproduces
#'         the fitted likelihood, and finite positive-definite information
#'         that is stable under step refinement.
#'         It does not call \code{vcov.gamlss} or require attaching \pkg{egpd}.
#'         A parameter with additive smoothing terms
#'         (\code{pb()}, \code{cs()}, \code{lo()}, \code{random()}, ...) is an
#'         error with \code{newdata} (use \code{method = "plug-in"} there).
#'         A \code{gamlss} parameter fixed in the fit (\code{<p>.fix = TRUE}) is
#'         held at its fitted value \code{fit$<p>.fv} by both methods (at
#'         \code{newdata} it must be a single value), and the standard errors at
#'         \code{newdata} are conditional on it.
#'         For \code{bamlss} fits, the stored posterior (MCMC) samples: the
#'         distributional parameters are predicted for each sample and one
#'         response is simulated per sample (samples are drawn without
#'         replacement when \code{nsim} does not exceed their number, and
#'         with replacement otherwise).  A \code{bamlss} fit without
#'         samples (e.g. \code{sampler = FALSE}) is an error.}
#'     }
#'   \item \strong{Response simulation.}
#'     \code{nsim} responses are drawn at each observation from the
#'     fitted distribution (using the appropriate \code{r*} function
#'     from \pkg{egpd} or, for external GAMLSS families, \pkg{gamlss.dist}).
#'   \item \strong{Intervals and coverage.}
#'     Equal-tailed quantile intervals are computed from the simulated
#'     responses by default. For discrete, continuous, and mixed responses,
#'     \code{use_shortest = TRUE} instead gives the shortest closed interval
#'     attaining the requested simulated mass (ties select the lowest interval).
#'     Its endpoints are simulated response values, hence integers for counts.
#'     Discrete equal-tailed intervals use empirical integer quantiles.
#'     Discreteness is taken from the family: \code{"degpd"}/\code{"zidegpd"}
#'     for \code{egpd} and \code{bamlss} fits, and the family \code{type}
#'     (\code{"Discrete"}) for \code{gamlss} fits.
#'     Coverage is membership in the returned closed interval, without jitter;
#'     discreteness can therefore make coverage exceed the nominal level.
#' }
#' Non-finite parameter draws or simulated responses cause an error rather
#' than being silently replaced or omitted. Identity-linked shape draws may
#' be zero or negative.
#' For a \code{gamlss} fit with aliased coefficients, both methods require
#' the same training-column dependencies to hold at \code{newdata}; otherwise
#' the predictor is not estimable. Parametric terms are evaluated in the
#' fitted basis (including saved knots and factor contrasts), with all
#' formula offsets added. Plug-in prediction does not require likelihood
#' curvature. With backfitted terms, a column named \code{source} need not
#' be renamed or removed from \code{newdata}, including when selected by
#' name in \code{ri(x.vars = ...)}.
#' At \code{newdata}, both \code{bamlss} methods reconstruct the parametric
#' design on the training basis and with the training factor contrasts. An
#' unavailable or unreproducible saved design is an error. Smooth contributions
#' are predicted using the fitted \code{bamlss} smooths. With posterior samples,
#' the plug-in value is the mean on the parameter scale.
#' BAMLSS coverage uses offsets actually supplied to the fitter through its
#' \code{offset} argument. A formula \code{offset()} that the fitter did not
#' use is not introduced at prediction time.
#'
#' @param fit a fitted model object (class \code{"egpd"},
#'   \code{"gamlss"}, or \code{"bamlss"}). For native \code{egpd}
#'   fits the supported families are \code{"egpd"}, \code{"degpd"} and
#'   \code{"zidegpd"}.
#' @param y numeric vector of observed responses, one per prediction row
#'   in the same order: the rows used to fit the model when \code{newdata}
#'   is \code{NULL}, or the rows of \code{newdata} otherwise.
#' @param newdata optional data frame of covariate values at which to
#'   compute predictions.  If \code{NULL} (the default), in-sample
#'   predictions are used, which needs the fitting data (an \code{egpd}
#'   fit made with \code{removeData = TRUE} requires \code{newdata}).
#' @param level nominal coverage level (default 0.95)
#' @param nsim number of parameter/response simulation draws (default
#'   2000)
#' @param method \code{"parametric"} (default) for parameter
#'   uncertainty, or \code{"plug-in"} for fixed-parameter prediction
#' @param use_shortest logical: use the shortest empirical interval instead
#'   of the equal-tailed interval? Applies also to continuous and mixed
#'   responses. Default \code{FALSE}.
#' @param data_fit (gamlss method only) optional data frame passed as
#'   \code{data} to \code{predict()}; needed when the original
#'   fitting data cannot be recovered from the fit object.
#' @param bd (gamlss binomial families only) number of trials at \code{newdata}:
#'   one non-negative integer per row, or one common size (strictly positive
#'   for the zero-adjusted \code{ZABI} and \code{ZABB} families). Supply this
#'   for grouped fits even if the fitted sizes were constant. For an ungrouped
#'   Bernoulli response the default at new rows is one trial. In sample,
#'   omit \code{bd}: the fitted sizes are used in fitting-row order. Supply
#'   \code{y} as success counts, not proportions or a two-column response.
#' @param ... additional arguments (currently unused)
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{family}}{character string identifying the distribution family}
#'   \item{\code{method}}{the method used (\code{"parametric"} or \code{"plug-in"})}
#'   \item{\code{level}}{the nominal coverage level}
#'   \item{\code{coverage}}{the empirical coverage (proportion of \code{y} falling within the prediction interval)}
#'   \item{\code{covered}}{logical vector indicating which observations are covered}
#'   \item{\code{L}}{lower prediction interval bounds}
#'   \item{\code{U}}{upper prediction interval bounds}
#' }
#'
#' @examples
#' \dontrun{
#' y <- rdiscegpd(500, sigma = 3, xi = 0.15, kappa = 2, type = 1)
#' df <- data.frame(y = y)
#'
#' # egpd fit
#' fit_e <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
#'               data = df, family = "degpd", degpd.args = list(m = 1))
#' predictive_coverage(fit_e, y)
#' }
#'
#' @export
predictive_coverage <- function(fit, ...) UseMethod("predictive_coverage")


## ---------------------------------------------------------------
## egpd method
## ---------------------------------------------------------------

#' @rdname predictive_coverage
#' @export
predictive_coverage.egpd <- function(fit, y, newdata = NULL, level = 0.95,
                                     nsim = 2000,
                                     method = c("parametric", "plug-in"),
                                     use_shortest = FALSE, ...) {
  method <- match.arg(method)
  family <- fit$family                          # "egpd", "degpd", "zidegpd"
  m <- fit$likfns$m
  dtype <- c(1L, 6L, 4L, 5L, 2L, 3L)[m]
  is_discrete <- family %in% c("degpd", "zidegpd")

  ## choose rng function
  rfun <- switch(family,
    egpd     = regpd,
    degpd    = rdiscegpd,
    zidegpd  = rzidiscegpd,
    stop("predictive_coverage not implemented for egpd family '", family, "'")
  )

  if (is.null(newdata) && is.null(fit$data))
    stop("the fit was made with removeData = TRUE: supply `newdata` with the ",
         "covariates of the observations in `y`", call. = FALSE)
  n <- if (!is.null(newdata)) nrow(newdata) else nrow(fit$data)
  if (length(y) != n) stop("length(y) must match the number of observations.")

  ## --- parameter draws (nsim x n matrices) -----------------------
  if (method == "parametric") {
    ## Use simulate.egpd infrastructure: multivariate normal on coefficients
    if (!is.null(newdata)) {
      sim <- simulate(fit, nsim = nsim, newdata = newdata, type = "response")
    } else {
      sim <- simulate(fit, nsim = nsim, type = "response")
    }
    ## sim is a named list of n x nsim matrices; transpose to nsim x n
    par_mat <- lapply(sim, function(mat) t(as.matrix(mat)))
  } else {
    ## Plug-in: fixed at MLE
    if (!is.null(newdata)) {
      pars <- predict(fit, newdata = newdata, type = "response")
    } else {
      pars <- predict(fit, type = "response")
    }
    ## the fitted parameters only (predict() adds an `xi` alias of `shape`)
    pars <- pars[names(pars) != "xi" | !("shape" %in% names(pars))]
    par_mat <- lapply(pars, function(v) {
      matrix(rep(as.numeric(v), each = nsim), nrow = nsim, ncol = n)
    })
  }

  ## Map parameter names to distribution argument names.
  ## After link-stripping, simulate.egpd / predict.egpd names are like
  ## "sigma"/"scale"/"psi", "xi"/"shape", "kappa", "delta", "pi", "dkappa", "p".
  ## We need to build a named par_list with names matching rfun args.
  nms <- names(par_mat)

  .find_par <- function(patterns) {
    for (pat in patterns) {
      idx <- grep(pat, nms, ignore.case = TRUE)
      if (length(idx)) return(idx[1])
    }
    NA_integer_
  }

  sigma_idx <- .find_par(c("^sigma$", "^scale$", "^psi$"))
  xi_idx    <- .find_par(c("^xi$", "^shape$"))

  par_list <- list(sigma = par_mat[[sigma_idx]], xi = par_mat[[xi_idx]])

  if (m %in% c(1, 5, 6)) {
    par_list$kappa <- par_mat[[.find_par(c("^kappa$"))]]
  } else if (m == 2) {
    k1_idx <- .find_par(c("^kappa1$", "^kappa$"))
    dk_idx <- .find_par(c("^dkappa$"))
    par_list$kappa <- par_mat[[k1_idx]]
    par_list$delta <- par_mat[[k1_idx]] + par_mat[[dk_idx]]
    par_list$prob  <- par_mat[[.find_par(c("^p$", "^prob$"))]]
  } else if (m == 3) {
    par_list$delta <- par_mat[[.find_par(c("^delta$"))]]
  } else if (m == 4) {
    par_list$delta <- par_mat[[.find_par(c("^delta$"))]]
    par_list$kappa <- par_mat[[.find_par(c("^kappa$"))]]
  }

  if (family == "zidegpd") {
    par_list$pi <- par_mat[[.find_par(c("^pi$"))]]
  }

  ## --- simulate responses ----------------------------------------
  ysim <- .sim_responses(par_list, rfun, dtype, nsim)

  ## --- intervals and coverage ------------------------------------
  res <- .pi_from_sims(ysim, y, level, use_shortest, is_discrete)

  list(
    family   = paste0(family, m),
    method   = method,
    level    = level,
    coverage = mean(res$covered),
    covered  = res$covered,
    L        = res$L,
    U        = res$U
  )
}


## ---------------------------------------------------------------
## gamlss method
## ---------------------------------------------------------------

.gamlss_family_name <- function(fit) {
  if (is.list(fit$family)) fit$family$family[1L] else fit$family[1L]
}

.gamlss_dist_ns <- function(fit) {
  rng <- paste0("r", .gamlss_family_name(fit))
  for (pkg in c("egpd", "gamlss.dist")) {
    ns <- asNamespace(pkg)
    if (exists(rng, ns, mode = "function", inherits = FALSE)) return(ns)
  }
  stop("No RNG found for family '", .gamlss_family_name(fit),
       "' (looked in egpd and gamlss.dist).", call. = FALSE)
}

## Prefer the fitted link function if retained. Standard fits store only
## the name; egpd families use stats links, gamlss.dist uses its own links.
.gamlss_link_function <- function(fit, p, egpd_family, which = "linkinv") {
  field <- function(suffix) {
    key <- paste0(p, suffix)
    if (!is.null(fit[[key]])) fit[[key]] else if (is.list(fit$family)) fit$family[[key]]
  }
  fn <- field(paste0(".", which))
  if (is.function(fn)) return(fn)
  ln <- field(".link")
  if (!is.character(ln) || length(ln) != 1L)
    stop("Could not determine the '", p, "' link for this gamlss fit.")
  link <- tryCatch(if (egpd_family) stats::make.link(ln) else gamlss.dist::make.link.gamlss(ln),
                   error = function(e) NULL)
  if (!is.function(link[[which]]))
    stop("Cannot reconstruct fitted ", p, " link '", ln, "' (", which, ").", call. = FALSE)
  link[[which]]
}

## Does gamlss parameter `p` have additive smoothing terms (pb(), cs(), lo(),
## random(), ...)? They are fitted by backfitting and have no coefficient
## representation, so no design matrix at new data.
.gamlss_has_smoother <- function(fit, p) {
  if (!is.null(fit[[paste0(p, ".s")]])) return(TRUE)
  tt <- fit[[paste0(p, ".terms")]]
  !is.null(tt) && length(unlist(attr(tt, "specials"))) > 0L
}

.gamlss_linear_block <- function(fit, p, egpd_family) {
  if (.gamlss_has_smoother(fit, p))
    stop("The formula for ", p, " has additive smoothing terms; use method = \"plug-in\" ",
         "at `newdata` or assess coverage in sample.", call. = FALSE)
  cf <- fit[[paste0(p, ".coefficients")]]
  X <- fit[[paste0(p, ".x")]]
  n <- length(fit$y)
  if (!is.matrix(X) || nrow(X) != n || is.null(cf) || ncol(X) != length(cf) ||
      !identical(colnames(X), names(cf)) || (length(cf) && is.null(names(cf))) ||
      anyDuplicated(names(cf)) || any(!is.finite(X)))
    stop("Need the saved, named coefficient design of '", p, "' at the fitting rows.",
         call. = FALSE)
  keep <- !is.na(cf)
  X <- X[, keep, drop = FALSE]
  b <- cf[keep]
  inv <- .gamlss_link_function(fit, p, egpd_family)
  fv <- as.numeric(fit[[paste0(p, ".fv")]])
  eta <- as.numeric(fit[[paste0(p, ".lp")]])
  agrees <- function(x) length(x) == n && length(fv) == n &&
    all(is.finite(x) & is.finite(fv) &
          abs(x - fv) <= 1e-10 * pmax(abs(fv), .Machine$double.xmin))
  if (!agrees(inv(eta))) {
    eta <- .gamlss_link_function(fit, p, egpd_family, "linkfun")(fv)
    if (!agrees(inv(eta)))
      stop("The saved link for '", p, "' does not reproduce its fitted values.", call. = FALSE)
  }
  off <- fit[[paste0(p, ".offset")]]
  if (is.null(off)) off <- 0
  if (!length(off) %in% c(1L, n) || any(!is.finite(off)) || any(!is.finite(eta)))
    stop("Non-finite or incomplete fitted predictor/offset for '", p, "'.", call. = FALSE)
  ## Step halving can leave the coefficients of the last working regression
  ## in the fit, not the coefficients of the accepted predictor (including offset).
  represents_eta <- function(b) all(is.finite(b)) &&
    all(abs(drop(X %*% b) + off - eta) <= 1e-8 * (1 + abs(eta)))
  if (!represents_eta(b)) {
    b <- stats::lm.wfit(X, eta - off, fit$weights)$coefficients
    if (!represents_eta(b))
      stop("The parametric design of '", p, "' does not reproduce its fitted predictor.",
           call. = FALSE)
  }
  coord <- if (length(b)) paste(p, names(b), sep = ".") else character()
  list(X = X, eta = eta, inv = inv, b = b, names = coord)
}

.gamlss_free_blocks <- function(fit, ns) {
  free <- Filter(function(p) !isTRUE(fit[[paste0(p, ".fix")]]), fit$parameters)
  lapply(stats::setNames(free, free), function(p)
    .gamlss_linear_block(fit, p, identical(ns, asNamespace("egpd"))))
}

## Choose a coefficient step by the *symmetric* change in individual log
## densities. This signal is quadratic in the step locally and dimensionless:
## changing response units can add a constant to log densities but not to
## their differences. The design column supplies only the initial probe.
.gamlss_fd_step <- function(change, zero, j, h) {
  ## For a fourth-order curvature estimate balance h^4 against eps/h^2:
  ## in local likelihood units the desired second difference is O(eps^(1/3)).
  target <- .Machine$double.eps^(1/3)
  lo <- 0
  hi <- Inf
  for (i in seq_len(100L)) {
    if (!is.finite(h) || h <= 0) break
    d <- zero
    d[j] <- h
    signal <- max(abs(change(d) + change(-d)))
    if (is.finite(signal) && signal >= target/4 && signal <= 4*target) return(h)
    if (!is.finite(signal) || signal > target) hi <- h else lo <- h
    if (lo > 0 && is.finite(hi)) {
      ## A support boundary can intervene before the target is reached.
      ## Use the interior step; the curvature must still pass refinement.
      if (hi/lo < 1.1) return(lo)
      h <- sqrt(lo) * sqrt(hi)
    } else {
      ratio <- if (is.finite(signal)) sqrt(target/signal) else 0
      h <- h * min(100, max(.01, ratio))
    }
  }
  stop("no resolved two-sided likelihood step for coefficient ", j)
}

.gamlss_fd_vcov <- function(change, w, scale) {
  k <- length(scale)
  zero <- numeric(k)
  h <- vapply(seq_len(k), function(j)
    .gamlss_fd_step(change, zero, j, .001/scale[j]), numeric(1))

  ## Information in units of h. Difference before summing over cases;
  ## use all four corners for cross-block as well as within-block terms.
  information <- function(r) {
    H <- matrix(0, k, k)
    for (i in seq_len(k)) {
      a <- zero
      a[i] <- r*h[i]
      H[i, i] <- sum(w * (change(a) + change(-a)))
      for (j in seq_len(i - 1L)) {
        b <- zero
        b[j] <- r*h[j]
        H[i, j] <- H[j, i] <- sum(w * ((change(a+b) - change(a-b)) -
                                          (change(-a+b) - change(-a-b)))) / 4
      }
    }
    H / r^2
  }
  coarse <- information(1)
  for (i in seq_len(10L)) {
    fine <- information(2^(-i))
    E <- fine - coarse
    coarse <- fine
    if (any(!is.finite(E))) next
    ## Central differences have O(h^2) bias. Require agreement in the
    ## metric of the extrapolated information, including its weak directions.
    R <- tryCatch(chol(fine + E/3), error = function(e) NULL)
    if (is.null(R)) next
    D <- backsolve(R, E, transpose = TRUE)
    D <- backsolve(R, t(D), transpose = TRUE)
    if (all(is.finite(D)) && norm(D, "I") < 1e-4) {
      V <- tcrossprod(h * backsolve(R, diag(k)))
      if (any(!is.finite(V)) || any(diag(V) <= 0)) stop("non-finite or unresolved covariance")
      return(V)
    }
  }
  stop("no stable positive-definite Hessian across two-sided likelihood steps")
}

## Differentiate only the free coordinates. gen.likelihood()/vcov.gamlss
## creates coordinates for fixed parameters using the default family links;
## conditioning its inverse cannot correct a Hessian evaluated at the wrong law.
.gamlss_vcov <- function(fit, ns = .gamlss_dist_ns(fit), blocks = .gamlss_free_blocks(fit, ns)) {
  fam <- .gamlss_family_name(fit)
  density <- get0(paste0("d", fam), ns, mode = "function", inherits = FALSE)
  if (is.null(density))
    stop("No density for the fitted ", fam, " likelihood.", call. = FALSE)
  y <- fit$y
  w <- fit$weights
  if (!is.null(dim(y)) || length(w) != length(y) || any(!is.finite(w) | w < 0) || !any(w > 0))
    stop("Need the fitting responses and non-negative case weights for ", fam, ".", call. = FALSE)
  active <- w > 0
  w <- w[active]
  args <- c(list(y[active]), lapply(stats::setNames(fit$parameters, fit$parameters), function(p) {
    fv <- fit[[paste0(p, ".fv")]]
    if (!length(fv) %in% c(1L, length(y)))
      stop("Incomplete fitted '", p, "' values in the ", fam, " likelihood.", call. = FALSE)
    rep_len(fv, length(y))[active]
  }), list(log = TRUE))
  if ("bd" %in% names(formals(density)))
    args$bd <- .gamlss_trials(fit, fam, y, NULL, NULL)[active]
  blocks <- lapply(blocks, function(block) {
    block$X <- block$X[active, , drop = FALSE]
    block$eta <- block$eta[active]
    block
  })
  coord <- unlist(lapply(blocks, `[[`, "names"), use.names = FALSE)
  idx <- lapply(blocks, function(block) match(block$names, coord))
  logdens <- function(delta) {
    for (p in names(blocks)) {
      block <- blocks[[p]]
      args[[p]] <- block$inv(block$eta + drop(block$X %*% delta[idx[[p]]]))
    }
    do.call(density, args)
  }
  zero <- numeric(length(coord))
  ll0 <- tryCatch(logdens(zero), error = function(e)
    stop("Cannot evaluate the fitted ", fam, " likelihood: ", conditionMessage(e), call. = FALSE))
  dev <- -2 * sum(w * ll0)
  if (length(ll0) != length(w) || !is.finite(dev) ||
      !isTRUE(all.equal(dev, fit$G.deviance, tolerance = 1e-8)))
    stop("The ", fam, " density at the expansion point does not reproduce the fitted ",
         "weighted log likelihood.", call. = FALSE)
  if (!length(coord)) return(matrix(numeric(0), 0L, 0L))
  scale <- unlist(lapply(blocks, function(block) apply(abs(block$X), 2L, max)),
                  use.names = FALSE)
  if (any(scale == 0))
    stop("No positive-weight information for ", paste(coord[scale == 0], collapse = ", "),
         ".", call. = FALSE)
  change <- function(delta) {
    ## An exploratory stencil may leave the support or link domain.
    ll <- tryCatch(logdens(delta), error = function(e) NA_real_, warning = function(w) NA_real_)
    if (length(ll) != length(ll0)) return(NA_real_)
    ll0 - ll
  }
  tryCatch({
    V <- .gamlss_fd_vcov(change, w, scale)
    dimnames(V) <- list(coord, coord)
    V
  }, error = function(e) stop("Could not invert the conditional observed information for ",
    fam, " at the fitted values: ", conditionMessage(e), call. = FALSE))
}

.gamlss_newdata_design <- function(fit, p, newdata, omit = character()) {
  cf <- fit[[paste0(p, ".coefficients")]]
  Terms <- stats::delete.response(fit[[paste0(p, ".terms")]])
  train <- fit[[paste0(p, ".x")]]
  if (length(omit)) {
    Terms <- stats::drop.terms(Terms, match(omit, attr(Terms, "term.labels")))
    train <- train[, !names(cf) %in% omit, drop = FALSE]
    cf <- cf[!names(cf) %in% omit]
  }
  mf <- stats::model.frame(Terms, newdata, xlev = fit[[paste0(p, ".xlevels")]],
                           na.action = stats::na.pass)
  X <- stats::model.matrix(Terms, mf,
                           contrasts.arg = attr(fit[[paste0(p, ".x")]], "contrasts"))
  if (nrow(X) != nrow(newdata) || !identical(colnames(X), names(cf)))
    stop("could not rebuild the design matrix of '", p, "' at `newdata`", call. = FALSE)
  keep <- !is.na(cf)
  Xfree <- X[, keep, drop = FALSE]
  if (any(!keep)) {
    ## An aliased training column may only be dropped at new rows where
    ## the same linear dependence holds.
    lost <- X[, !keep, drop = FALSE]
    if (!is.matrix(train) || !identical(colnames(train), names(cf)) ||
        nrow(train) != length(fit$weights))
      stop("Need the saved coefficient design of '", p, "' to check estimability.",
           call. = FALSE)
    active <- fit$weights > 0
    alias <- train[active, !keep, drop = FALSE]
    restored <- size <- 0
    if (any(keep)) {
      basis <- train[active, keep, drop = FALSE]
      map <- stats::lm.wfit(basis, alias, fit$weights[active])$coefficients
      restored <- Xfree %*% map
      size <- abs(Xfree) %*% abs(map)
    }
    ## Relative tolerance in each dropped column's units, plus roundoff
    ## from fitting and applying the relation (also when lost is zero).
    eps <- .Machine$double.eps
    tol <- sweep(1e-7 * abs(lost) + ncol(Xfree)*eps*size, 2L,
                 nrow(alias)*eps*apply(abs(alias), 2L, max), "+")
    if (!isTRUE(all(abs(restored - lost) <= tol)))
      stop("Non-estimable '", p, "' predictor at `newdata` (aliased training columns).",
           call. = FALSE)
  }
  list(X = Xfree, offset = stats::model.offset(mf))
}

.gamlss_design_lp <- function(design, block) {
  drop(design$X %*% block$b) + (design$offset %||% 0)
}

## Backfitted smooths need predict.gamlss(). It rebuilds the formula on the
## combined training and prediction rows and accepts only one offset. Pass
## it the fitted calls for the parametric terms and the sum of the offsets.
.gamlss_smooth_prediction_fit <- function(fit, p) {
  key <- paste0(p, ".terms")
  tt <- stats::delete.response(fit[[key]])
  vars <- as.list(attr(tt, "variables"))[-1L]
  pred <- as.list(attr(tt, "predvars"))[-1L]
  off <- attr(tt, "offset")
  if (!length(pred)) pred <- vars
  if (length(off) < 2L && identical(vars, pred)) return(fit)
  smooth <- colnames(fit[[paste0(p, ".s")]])
  freeze <- function(expr) {
    i <- which(vapply(vars, identical, logical(1), expr))
    if (length(i)) return(pred[[i[1L]]])
    if (!is.call(expr)) return(expr)
    as.call(c(list(expr[[1L]]), lapply(as.list(expr)[-1L], freeze)))
  }
  terms <- lapply(attr(tt, "term.labels"), function(label) {
    expr <- str2lang(label)
    if (label %in% smooth) expr else freeze(expr)
  })
  plus <- function(a, b) call("+", a, b)
  if (length(off))
    terms <- c(terms, list(call("offset", Reduce(plus, lapply(pred[off], `[[`, 2L)))))
  rhs <- Reduce(plus, terms, init = as.integer(attr(tt, "intercept")))
  form <- stats::as.formula(call("~", rhs), env = environment(tt))
  fit[[paste0(p, ".formula")]] <- form
  fit[[key]] <- stats::terms(form)
  fit
}

## Rename the variable, not a call head, argument tag, arbitrary string or
## component name (as in settings$source). ri(x.vars=) selects columns by
## name or position: the fitted ridge coefficients record the resolved names.
.gamlss_alias_source <- function(expr, alias, smoothers) {
  if (identical(expr, quote(source))) return(alias)
  if (!is.call(expr)) return(expr)
  fun <- expr[[1L]]
  if (identical(fun, quote(ri)) || identical(fun, quote(gamlss::ri))) {
    cols <- rownames(smoothers[[deparse1(expr)]][["coef"]])
    if ("source" %in% cols) {
      cl <- match.call(gamlss::ri, expr)
      if (!is.null(cl$x.vars)) {
        cols[cols == "source"] <- as.character(alias)
        cl$x.vars <- as.call(c(list(quote(c)), as.list(cols)))
        expr <- cl
      }
    }
  }
  i <- seq_along(expr)[-1L]
  if (is.symbol(fun)) {
    if (as.character(fun) %in% c("::", ":::", "quote")) return(expr)
    if (as.character(fun) %in% c("$", "@")) i <- 2L
  }
  expr[i] <- lapply(as.list(expr)[i], .gamlss_alias_source, alias, smoothers)
  expr
}

## predict.gamlss() subsets the training frame to names(newdata), then
## adds a column called 'source' to distinguish the two sets of rows.
.gamlss_smooth_prediction_args <- function(fit, p, newdata, data_fit) {
  fit <- .gamlss_smooth_prediction_fit(fit, p)
  key <- paste0(p, ".terms")
  tt <- stats::delete.response(fit[[key]])
  form <- stats::formula(tt)
  used <- c(names(newdata), names(data_fit), all.names(form))
  alias <- ".src"
  while (alias %in% used) alias <- paste0(alias, "_")
  smoothers <- fit[[paste0(p, ".coefSmo")]]
  labs <- colnames(fit[[paste0(p, ".s")]])
  if (length(smoothers) == length(labs) && length(labs))
    names(smoothers) <- vapply(labs, function(s) deparse1(str2lang(s)), "")
  rhs <- .gamlss_alias_source(form[[2L]], as.name(alias), smoothers)
  if (!identical(form[[2L]], rhs)) {
    if (!"source" %in% names(newdata)) stop("Missing `source` in `newdata`.")
    if (is.null(data_fit)) data_fit <- eval(fit$call$data, environment(tt))
    if (!"source" %in% names(data_fit)) stop("Missing `source` in the fitting data.")
    names(newdata)[names(newdata) == "source"] <- alias
    data_fit[[alias]] <- data_fit[["source"]]
    form[[2L]] <- rhs
    tt1 <- stats::terms(form)

    ## Smooth contributions and factor levels are looked up by model-frame
    ## label (e.g. "pb(source)"), not by the name of the raw column.
    before <- rownames(attr(tt, "factors"))
    after <- rownames(attr(tt1, "factors"))
    relabel <- function(x) {
      i <- match(x, before)
      x[!is.na(i)] <- after[i[!is.na(i)]]
      x
    }
    skey <- paste0(p, ".s")
    xkey <- paste0(p, ".xlevels")
    if (!is.null(fit[[skey]])) colnames(fit[[skey]]) <- relabel(colnames(fit[[skey]]))
    if (length(fit[[xkey]])) names(fit[[xkey]]) <- relabel(names(fit[[xkey]]))
    fit[[paste0(p, ".formula")]] <- form
    fit[[key]] <- tt1
  } else {
    newdata[["source"]] <- NULL
  }
  list(object = fit, newdata = newdata, data = data_fit)
}

## Marginal link SE from the named block of the *joint free* covariance.
.gamlss_newdata_se <- function(fit, p, newdata, V, eta = NULL, block = NULL,
                               design = .gamlss_newdata_design(fit, p, newdata)) {
  if (isTRUE(fit[[paste0(p, ".fix")]])) return(rep(0, nrow(newdata)))
  if (is.null(block))
    block <- .gamlss_linear_block(fit, p, identical(.gamlss_dist_ns(fit), asNamespace("egpd")))
  X <- design$X
  if (!is.null(eta)) {
    lp <- .gamlss_design_lp(design, block)
    if (length(eta) != length(lp) || !isTRUE(all(abs(lp - eta) <= 1e-5 * (1 + abs(eta)))))
      stop("the design matrix of '", p, "' rebuilt at `newdata` does not reproduce ",
           "the supplied link prediction", call. = FALSE)
  }
  if (!ncol(X)) return(rep(0, nrow(X)))  # offset-only predictor
  Vb <- V[block$names, block$names, drop = FALSE]
  unname(sqrt(rowSums((X %*% Vb) * X)))
}

## bd is observed data, not a distributional parameter or a fitted weight.
.gamlss_trials <- function(fit, family, y, newdata, bd) {
  fitted_bd <- fit[["bd"]]
  valid_sizes <- function(x) is.numeric(x) && is.null(dim(x)) &&
    all(is.finite(x) & x >= 0 & x == floor(x))
  if (!valid_sizes(fitted_bd) || length(fitted_bd) != NROW(fit$y))
    stop("Expected one non-negative integer trial count in `fit$bd` per fitted row.",
         call. = FALSE)
  if (is.null(newdata)) {
    if (!is.null(bd))
      stop("Supply `bd` only with `newdata`; in sample the fitted sizes are used.",
           call. = FALSE)
    bd <- fitted_bd
  } else {
    if (is.null(bd)) {
      ## A two-column response is grouped even if every training size was one.
      response_class <- unname(attr(fit$mu.terms, "dataClasses")[1L])
      if (any(fitted_bd != 1) || identical(response_class, "nmatrix.2"))
        stop("Supply `bd` (number of trials) for grouped binomial `newdata`.",
             call. = FALSE)
      bd <- 1L
    }
    if (!valid_sizes(bd) || !length(bd) %in% c(1L, length(y)))
      stop("`bd` must be a non-negative integer per new row or one common size.",
           call. = FALSE)
    bd <- rep_len(bd, length(y))
  }
  ## The positive-count component of a hurdle law cannot exist with no trials.
  if (family %in% c("ZABI", "ZABB") && (any(fitted_bd == 0) || any(bd == 0)))
    stop("Zero-adjusted ", family, " requires positive trial counts (`bd`).", call. = FALSE)
  if (!(is.numeric(y) || is.logical(y)) || !is.null(dim(y)) ||
      any(!is.finite(y) | y != floor(y) | y < 0 | y > bd))
    stop("Binomial `y` must contain integer success counts between 0 and `bd`.",
         call. = FALSE)
  bd
}

#' @rdname predictive_coverage
#' @export
predictive_coverage.gamlss <- function(fit, y, newdata = NULL, level = 0.95,
                                       nsim = 2000,
                                       method = c("parametric", "plug-in"),
                                       use_shortest = FALSE,
                                       data_fit = NULL, bd = NULL, ...) {
  method <- match.arg(method)

  fam_name <- .gamlss_family_name(fit)
  ns <- .gamlss_dist_ns(fit)
  egpd_family <- identical(ns, asNamespace("egpd"))
  rng_fun <- get(paste0("r", fam_name), ns, inherits = FALSE)

  ## Discrete or continuous from the gamlss family type ("Discrete", "Continuous",
  ## "Mixed"), which fits carry as fit$type; a name pattern missed PO, NBI and the
  ## package's own PIG/GPIG/BELL families.
  fam_type <- fit$type
  if (is.null(fam_type)) {
    ctor <- get0(fam_name, ns, mode = "function", inherits = FALSE)
    fam_type <- tryCatch(ctor()$type, error = function(e) NULL)
  }
  if (!is.character(fam_type) || length(fam_type) != 1L)
    stop("Could not determine whether gamlss family '", fam_name,
         "' is discrete (no family type).", call. = FALSE)
  is_discrete <- identical(fam_type, "Discrete")

  params <- fit$parameters
  if (!length(params) || anyDuplicated(params) || !all(params %in% c("mu", "sigma", "nu", "tau")))
    stop("Expected named mu/sigma/nu/tau parameters in the gamlss fit.")

  n <- if (!is.null(newdata)) nrow(newdata) else NROW(fit$y)
  if (length(y) != n) stop("length(y) must match the number of observations.")
  if ("bd" %in% names(formals(rng_fun))) {
    bd <- .gamlss_trials(fit, fam_name, y, newdata, bd)
  } else if (!is.null(fit[["bd"]]) || !is.null(bd)) {
    stop("The RNG for this family does not accept binomial trial counts (`bd`).",
         call. = FALSE)
  }

  ## A parameter fixed in the fit (gamlss(<p>.fix = TRUE)) is a known constant without
  ## a link, coefficients or standard error: its fitted values <p>.fv. gamlss cannot
  ## predict it at new data, where it must therefore be a single value.
  fixed <- vapply(params, function(p) isTRUE(fit[[paste0(p, ".fix")]]), logical(1))
  free <- params[!fixed]
  .fixed_value <- function(p) {
    v <- as.numeric(fit[[paste0(p, ".fv")]])
    if (!is.null(newdata)) {
      if (length(unique(v)) != 1L)
        stop("parameter '", p, "' is fixed at values that vary between observations, ",
             "which are unknown at `newdata`", call. = FALSE)
      return(rep(v[1L], n))
    }
    if (length(v) == 1L) return(rep(v, n))
    if (length(v) != n) stop("Wrong number of fixed values for '", p, "'.", call. = FALSE)
    v
  }

  fixed_values <- lapply(params[fixed], .fixed_value)
  names(fixed_values) <- params[fixed]

  ## Use the upstream link method for in-sample SEs or backfitted smooths;
  ## it does not need to re-evaluate the original family constructor.
  .predict_param <- function(p, se.fit = FALSE) {
    tryCatch({
      if (is.null(newdata))
        return(suppressWarnings(predict(fit, what = p, type = "link", se.fit = se.fit,
                                        data = data_fit)))
      args <- .gamlss_smooth_prediction_args(fit, p, newdata, data_fit)
      do.call(predict, c(args, list(what = p, type = "link", se.fit = se.fit)))
    }, error = function(e) {
      hint <- if (!is.null(newdata) && is.null(data_fit)) " Supply the original `data_fit`."
      stop("Cannot predict '", p, "' from the gamlss fit: ", conditionMessage(e), hint,
           call. = FALSE)
    })
  }

  ## Parameter draws on response scale (nsim x n)
  par_draws <- switch(method,
    "plug-in" = {
      out <- list()
      for (p in params) {
        v <- if (fixed[[p]]) fixed_values[[p]] else if (is.null(newdata)) {
          as.numeric(fit[[paste0(p, ".fv")]])
        } else if (.gamlss_has_smoother(fit, p)) {
          cf <- fit[[paste0(p, ".coefficients")]]
          if (anyNA(cf)) {
            ## lo(), re(), etc. contribute a fitted smooth, not a
            ## parametric coefficient. Their all-zero design placeholders
            ## must not be evaluated outside the smoother prediction.
            lost <- names(cf)[is.na(cf)]
            sm <- intersect(lost, colnames(fit[[paste0(p, ".s")]]))
            X <- fit[[paste0(p, ".x")]]
            placeholders <- sm[vapply(sm, function(s) all(X[, s] == 0), logical(1))]
            if (length(setdiff(lost, placeholders)))
              .gamlss_newdata_design(fit, p, newdata, omit = placeholders)
          }
          inv <- .gamlss_link_function(fit, p, egpd_family)
          as.numeric(inv(.predict_param(p)))
        } else {
          block <- .gamlss_linear_block(fit, p, egpd_family)
          block$inv(.gamlss_design_lp(.gamlss_newdata_design(fit, p, newdata), block))
        }
        if (length(v) != n) stop("Wrong number of predictions for '", p, "'.", call. = FALSE)
        out[[p]] <- matrix(rep(v, each = nsim), nrow = nsim, ncol = n)
      }
      out
    },
    "parametric" = {
      if (!is.null(newdata) && length(free)) {
        blocks <- .gamlss_free_blocks(fit, ns)
        Vg <- .gamlss_vcov(fit, ns, blocks)
      }
      draws <- list()
      for (p in params) {
        if (fixed[[p]]) {
          draws[[p]] <- matrix(rep(fixed_values[[p]], each = nsim), nrow = nsim, ncol = n)
          next
        }
        if (is.null(newdata)) {
          pred <- .predict_param(p, se.fit = TRUE)
          if (!is.list(pred) || is.null(pred$se.fit))
            stop("gamlss returned no standard errors for the '", p, "' linear predictor.",
                 call. = FALSE)
          m <- as.numeric(pred$fit)
          s <- as.numeric(pred$se.fit)
        } else {
          design <- .gamlss_newdata_design(fit, p, newdata)
          m <- .gamlss_design_lp(design, blocks[[p]])
          s <- .gamlss_newdata_se(fit, p, newdata, Vg, block = blocks[[p]], design = design)
        }
        if (length(m) != n || length(s) != n)
          stop("Wrong number of predictions or standard errors for '", p, "'.", call. = FALSE)

        if (any(!is.finite(m)) || any(!is.finite(s) | s < 0))
          stop("Non-finite link predictions or invalid standard errors for '", p, "'.")

        eta <- matrix(
          rnorm(nsim * n, mean = rep(m, each = nsim), sd = rep(s, each = nsim)),
          nrow = nsim, ncol = n
        )
        invlink <- .gamlss_link_function(fit, p, egpd_family)
        val <- invlink(eta)
        if (any(!is.finite(val)))
          stop("Non-finite parameter draws for '", p, "'.")
        draws[[p]] <- matrix(val, nrow = nsim, ncol = n)
      }
      draws
    }
  )

  for (p in names(par_draws)) {
    if (any(!is.finite(par_draws[[p]])))
      stop("Non-finite parameter draws for '", p, "'.")
  }

  ## Simulate responses
  ysim <- matrix(NA_real_, nsim, n)
  for (j in seq_len(n)) {
    args <- list(n = nsim)
    if (!is.null(bd)) args$bd <- bd[j]
    if ("mu"    %in% params) args$mu    <- par_draws$mu[, j]
    if ("sigma" %in% params) args$sigma <- par_draws$sigma[, j]
    if ("nu"    %in% params) args$nu    <- par_draws$nu[, j]
    if ("tau"   %in% params) args$tau   <- par_draws$tau[, j]
    ysim[, j] <- do.call(rng_fun, args)
  }

  ## Intervals and coverage
  res <- .pi_from_sims(ysim, y, level, use_shortest, is_discrete)

  list(
    family   = fam_name,
    method   = method,
    level    = level,
    coverage = mean(res$covered),
    covered  = res$covered,
    L        = res$L,
    U        = res$U
  )
}


## ---------------------------------------------------------------
## bamlss method
## ---------------------------------------------------------------

#' @rdname predictive_coverage
#' @param offset (bamlss method only) link-scale offsets at \code{newdata}
#'   for parameters with a saved fitting offset. Supply a named list of
#'   finite numeric vectors (one per parameter, of length one or the number
#'   of rows), or a numeric vector for a single offset. A constant fitted
#'   offset is reused if omitted; a varying fitted offset must be supplied.
#'   Omit this argument for in-sample coverage.
#' @export
predictive_coverage.bamlss <- function(fit, y, newdata = NULL, level = 0.95,
                                       nsim = 2000,
                                       method = c("parametric", "plug-in"),
                                       use_shortest = FALSE, offset = NULL, ...) {
  method <- match.arg(method)

  fam <- fit$family
  fam_name   <- fam$family        # e.g. "egpd1", "degpd1", "zidegpd1"
  param_names <- fam$names        # e.g. c("sigma","xi","kappa")

  ## Determine family type, rng function, and discrete flag
  base_fam <- sub("[0-9]+$", "", fam_name)   # "egpd", "degpd", "zidegpd", "ziegpd"
  m_val <- as.integer(sub("^[a-z]+", "", fam_name))
  dtype <- c(1L, 6L, 4L, 5L, 2L, 3L)[m_val]
  is_discrete <- base_fam %in% c("degpd", "zidegpd")

  rfun <- switch(base_fam,
    egpd     = regpd,
    degpd    = rdiscegpd,
    zidegpd  = rzidiscegpd,
    ziegpd   = rziegpd,
    stop("predictive_coverage not implemented for bamlss family '", fam_name, "'")
  )

  n <- if (!is.null(newdata)) nrow(newdata) else NROW(fit$y)
  if (length(y) != n) stop("length(y) must match the number of observations.")
  offsets <- .bamlss_coverage_offsets(fit, newdata, offset)

  if (method == "plug-in") {
    ## predict(type="parameter") returns named list of vectors on response scale
    pars <- if (is.null(newdata) && !length(offsets)) predict(fit, type = "parameter") else
      lapply(.bamlss_design_draws(fit, newdata, param_names, offsets), rowMeans)
    par_mat <- lapply(param_names, function(nm) {
      v <- as.numeric(pars[[nm]])
      matrix(rep(v, each = nsim), nrow = nsim)
    })
    names(par_mat) <- param_names
  } else {
    if (is.null(fit$samples))
      stop("method = \"parametric\" needs the posterior samples of the bamlss fit, ",
           "and this fit has none (was it fitted with sampler = FALSE?); refit with ",
           "a sampler or use method = \"plug-in\".", call. = FALSE)
    par_mat <- .bamlss_param_draws(fit, newdata, param_names, nsim, offsets)
  }
  if (any(vapply(par_mat, ncol, integer(1)) != n))
    stop("The bamlss predictions do not match the number of observations.", call. = FALSE)

  ## Simulate responses
  ysim <- .sim_responses(par_mat, rfun, dtype, nsim)

  ## Intervals and coverage
  res <- .pi_from_sims(ysim, y, level, use_shortest, is_discrete)

  list(
    family   = fam_name,
    method   = method,
    level    = level,
    coverage = mean(res$covered),
    covered  = res$covered,
    L        = res$L,
    U        = res$U
  )
}

## Posterior parameter draws of a bamlss fit as nsim x n matrices, one stored MCMC
## sample per row: predict(FUN = identity) returns each parameter for every sample
## (an n x S table), and rows are sampled from the S samples -- without replacement
## when nsim <= S. The sample indices are attached as attr(, "draws").
.bamlss_param_draws <- function(fit, newdata, param_names, nsim, offsets = NULL) {
  pp <- if (is.null(newdata) && !length(offsets)) predict(fit, type = "parameter", FUN = identity) else
    .bamlss_design_draws(fit, newdata, param_names, offsets)
  S <- ncol(as.matrix(pp[[param_names[1L]]]))
  if (!is.finite(S) || S < 1L)
    stop("the bamlss fit returned no posterior parameter samples", call. = FALSE)
  idx <- if (nsim <= S) sample.int(S, nsim) else sample.int(S, nsim, replace = TRUE)
  out <- lapply(param_names, function(nm) {
    M <- as.matrix(pp[[nm]])
    if (ncol(M) != S)
      stop("posterior samples of '", nm, "' do not match those of '",
           param_names[1L], "'", call. = FALSE)
    t(M[, idx, drop = FALSE])
  })
  names(out) <- param_names
  attr(out, "draws") <- idx
  out
}

.bamlss_coverage_offsets <- function(fit, newdata, offset) {
  if (is.null(newdata) && !is.null(offset))
    stop("Omit 'offset' in sample: the saved fitting offsets are used.", call. = FALSE)
  if (any(vapply(fit$x, function(b) length(attr(b$terms, "offset")) > 0L, TRUE)))
    warning("BAMLSS coverage uses saved fitting offsets, not offset() terms in the formula.",
            call. = FALSE)
  saved <- fit$model.frame[["(offset)"]]
  if (is.null(saved)) {
    if (!is.null(offset)) stop("The BAMLSS fit has no saved fitting offset.", call. = FALSE)
    return(NULL)
  }
  saved <- as.list(as.data.frame(saved))
  if (!all(names(saved) %in% fit$family$names) || anyDuplicated(names(saved)))
    stop("Cannot match the saved BAMLSS offsets to the fitted parameters.", call. = FALSE)
  if (is.null(newdata)) return(saved)
  if (is.numeric(offset) && length(saved) == 1L) offset <- setNames(list(offset), names(saved))
  if (length(offset) && (!is.list(offset) || is.null(names(offset)) || anyDuplicated(names(offset)) ||
                         !all(names(offset) %in% names(saved))))
    stop("Supply 'offset' by fitted BAMLSS parameter name.", call. = FALSE)
  for (nm in names(saved)) {
    v <- offset[[nm]]
    if (is.null(v)) {
      v <- unique(saved[[nm]])
      if (length(v) != 1L)
        stop("Supply the varying BAMLSS offset for '", nm, "' at newdata.", call. = FALSE)
    }
    if (!is.numeric(v) || !is.null(dim(v)) || !length(v) ||
        !length(v) %in% c(1L, nrow(newdata)) || any(!is.finite(v)))
      stop("Invalid BAMLSS offset for '", nm, "' at newdata.", call. = FALSE)
    saved[[nm]] <- v
  }
  saved
}

## predict.bamlss rebuilds parametric terms on the query batch. Recover
## predvars and contrasts on the fitting rows and verify them against X.
.bamlss_fixed_design <- function(block, train, newdata, env) {
  tt <- block$terms
  if (!inherits(tt, "terms")) stop("BAMLSS prediction needs the fitted model terms.", call. = FALSE)
  labs <- attr(tt, "term.labels")
  sp <- unlist(attr(tt, "specials"), use.names = FALSE)
  if (length(sp))
    labs <- labs[colSums(attr(tt, "factors")[sp, , drop = FALSE] != 0) == 0]
  int <- attr(tt, "intercept") == 1L
  rhs <- if (length(labs)) reformulate(labs, intercept = int) else if (int) ~1 else ~0
  environment(rhs) <- env
  tt <- stats::terms(rhs)
  pv <- attr(tt, "variables")
  mf <- train[, integer(), drop = FALSE]
  for (j in seq_along(pv)[-1L]) {
    e <- pv[[j]]
    nm <- deparse1(e)
    v <- if (nm %in% names(train)) train[[nm]] else eval(e, train, env)
    mf[[nm]] <- v
    pv[[j]] <- stats::makepredictcall(v, e)
  }
  attr(tt, "predvars") <- pv
  attr(mf, "terms") <- tt
  saved <- block$model.matrix
  ctr <- attr(saved, "contrasts")
  X0 <- stats::model.matrix(tt, mf, contrasts.arg = ctr)
  same <- if (!ncol(X0)) is.null(saved) || ncol(saved) == 0L else
    identical(dim(X0), dim(saved)) && identical(colnames(X0), colnames(saved)) &&
      isTRUE(all.equal(as.numeric(X0), as.numeric(saved), tolerance = 1e-12))
  if (!isTRUE(same))
    stop("Cannot reproduce the fitted BAMLSS parametric design; use explicit basis columns.",
         call. = FALSE)
  nd <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                           xlev = stats::.getXlevels(tt, mf))
  X <- stats::model.matrix(tt, nd, contrasts.arg = ctr)
  if (nrow(X) != nrow(newdata) || !identical(colnames(X), colnames(X0)))
    stop("Held-out BAMLSS design does not match the fitted columns.", call. = FALSE)
  X
}

.bamlss_design_draws <- function(fit, newdata, param_names, offsets) {
  train <- fit$model.frame
  if (!is.data.frame(train) || !all(param_names %in% names(fit$x)))
    stop("Held-out BAMLSS predictions require the saved training frame and designs.", call. = FALSE)
  attr(newdata, "terms") <- NULL
  for (nm in intersect(names(train), names(newdata))) {
    v <- train[[nm]]
    if (is.factor(v)) {
      if (any(!is.na(newdata[[nm]]) & !newdata[[nm]] %in% levels(v)))
        stop("New level of BAMLSS factor '", nm, "'.", call. = FALSE)
      newdata[[nm]] <- factor(newdata[[nm]], levels = levels(v), ordered = is.ordered(v))
    }
  }
  nd <- if (is.null(newdata)) train else newdata
  sc <- attr(train, "scale", exact = TRUE)
  if (!is.null(newdata)) for (nm in intersect(names(sc$center), names(nd)))
    nd[[nm]] <- (nd[[nm]] - sc$center[nm])/sc$scale[nm]
  B <- if (!is.null(fit$samples)) as.matrix(bamlss::samples(fit)) else {
    b <- bamlss::parameters(fit, list = FALSE, extract = TRUE)
    if (is.null(dim(b))) matrix(b, 1L, dimnames = list(NULL, names(b))) else as.matrix(b)
  }
  out <- setNames(vector("list", length(param_names)), param_names)
  for (nm in param_names) {
    block <- fit$x[[nm]]
    X <- if (is.null(newdata)) block$model.matrix else
      .bamlss_fixed_design(block, train, nd, environment(fit$formula))
    eta <- matrix(0, nrow(nd), nrow(B))
    upstream <- function(type, ...) {
      v <- predict(fit, newdata = newdata, model = nm, type = type, FUN = identity,
                   drop = FALSE, ...)[[nm]]
      ## A single draw can be a vector; multiple draws are a data frame
      ## with one column per sample, not a vector of scalar predictions.
      if (is.null(dim(v)) && length(v) == length(eta) && 1L %in% dim(eta))
        dim(v) <- dim(eta)
      v <- as.matrix(v)
      if (!identical(dim(v), dim(eta)))
        stop("BAMLSS predictions do not match the held-out rows and coefficient draws.", call. = FALSE)
      v
    }
    if (NCOL(X)) {
      cols <- paste(nm, "p", colnames(X), sep = ".")
      if (!all(cols %in% colnames(B))) cols <- paste(nm, "p.model.matrix", colnames(X), sep = ".")
      if (!all(cols %in% colnames(B)))
        stop("Cannot match BAMLSS design to coefficients for '", nm, "'.", call. = FALSE)
      eta <- X %*% t(B[, cols, drop = FALSE])
    }
    if (!is.null(offsets[[nm]])) eta <- eta + offsets[[nm]]
    smooth <- names(block$smooth.construct)
    if (length(smooth))
      eta <- eta + upstream("link", term = smooth, match.names = FALSE, intercept = FALSE)
    link <- fit$family$links[[nm]]
    if (link %in% c("log", "logit", "identity", "probit", "cauchit", "cloglog", "sqrt", "inverse", "1/mu^2")) {
      inv <- stats::make.link(link)$linkinv
      out[[nm]] <- matrix(inv(as.vector(eta)), nrow(eta), ncol(eta))
    } else {
      ## A bamlss-specific inverse link is upstream's responsibility, but
      ## use it only if upstream used the verified design for every draw.
      lp <- upstream("link")
      if (!isTRUE(all.equal(as.numeric(lp), as.numeric(eta), tolerance = 1e-12)))
        stop("Cannot verify the held-out BAMLSS predictor for link '", link,
             "'; fit with explicit basis columns.", call. = FALSE)
      out[[nm]] <- upstream("parameter")
    }
  }
  out
}
