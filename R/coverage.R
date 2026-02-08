## Predictive coverage for egpd, gamlss, and bamlss fits
## -------------------------------------------------------
## S3 generic + methods for computing simulation-based prediction
## intervals and their empirical coverage.

## ---------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------

#' @keywords internal
.shortest_discrete_set <- function(pmf, level = 0.95) {
  ord <- order(pmf, decreasing = TRUE)
  k <- which(cumsum(pmf[ord]) >= level)[1]
  S <- sort(as.integer(names(pmf))[ord[seq_len(k)]])
  c(min(S), max(S))
}

## Simulate responses from an EGPD-family distribution.
## `par_list` is a named list; each element is either an nsim x n
## matrix (one column per observation, one row per simulation draw)
## or a length-n vector (plug-in case).
## `rfun` is one of regpd / rdiscegpd / rzidiscegpd / rziegpd,
## `dtype` is the integer type code (1/4/5/6), and `nsim` is the
## number of draws per observation.
## Returns an nsim x n matrix.
.sim_responses <- function(par_list, rfun, dtype, nsim) {
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
  ## Replace non-finite simulated values with NA so quantile() can handle them
  ysim[!is.finite(ysim)] <- NA
  alpha <- 1 - level
  if (!use_shortest) {
    if (is_discrete) {
      ## For discrete distributions, quantile-based intervals systematically
      ## overcover because the integer-valued CDF jumps past the target
      ## probability.  Adding uniform jitter U(-0.5, 0.5) creates a
      ## continuous distribution whose quantiles interpolate between the
      ## discrete steps.  Coverage is checked against these continuous
      ## bounds so that boundary values are included/excluded
      ## proportionally, giving well-calibrated coverage on average.
      ysim_j <- ysim + matrix(runif(length(ysim), -0.5, 0.5),
                               nrow = nrow(ysim), ncol = ncol(ysim))
      L <- apply(ysim_j, 2, quantile, probs = alpha / 2, type = 7,
                 names = FALSE, na.rm = TRUE)
      U <- apply(ysim_j, 2, quantile, probs = 1 - alpha / 2, type = 7,
                 names = FALSE, na.rm = TRUE)
    } else {
      L <- apply(ysim, 2, quantile, probs = alpha / 2, type = 7,
                 names = FALSE, na.rm = TRUE)
      U <- apply(ysim, 2, quantile, probs = 1 - alpha / 2, type = 7,
                 names = FALSE, na.rm = TRUE)
    }
  } else {
    L <- U <- numeric(n)
    for (j in seq_len(n)) {
      col_j <- ysim[, j]
      col_j <- col_j[!is.na(col_j)]
      tab <- table(col_j)
      pmf <- as.numeric(tab) / length(col_j)
      names(pmf) <- names(tab)
      lims <- .shortest_discrete_set(pmf, level = level)
      L[j] <- lims[1]; U[j] <- lims[2]
    }
  }
  ## For discrete data, jitter the observed values too (randomised PIT

  ## approach).  This ensures that boundary observations are included
  ## proportionally rather than all-or-nothing, giving smoother and
  ## better-calibrated coverage even when all observations share the
  ## same distribution (intercept-only models).
  if (is_discrete && !use_shortest) {
    y_for_cov <- y + runif(n, -0.5, 0.5)
  } else {
    y_for_cov <- y
  }
  covered <- (y_for_cov >= L) & (y_for_cov <= U)
  ## Round bounds for display (integer intervals for discrete families).
  if (is_discrete) { L <- as.integer(ceiling(L)); U <- as.integer(floor(U)) }
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
#'       \item{\code{method = "parametric"}}{Parameters are perturbed
#'         on the link scale using normal noise (scaled by standard
#'         errors) and then back-transformed.  For \code{egpd} fits
#'         this is a multivariate-normal draw from the joint
#'         coefficient posterior; for \code{gamlss} fits the
#'         perturbation is independent across parameters.}
#'     }
#'   \item \strong{Response simulation.}
#'     \code{nsim} responses are drawn at each observation from the
#'     fitted distribution (using the appropriate \code{r*} function
#'     from the egpd package).
#'   \item \strong{Intervals and coverage.}
#'     Equal-tailed quantile intervals are computed from the simulated
#'     responses.  For discrete families, \code{use_shortest = TRUE}
#'     gives the shortest (highest-density-region) interval instead.
#' }
#'
#' @param fit a fitted model object (class \code{"egpd"},
#'   \code{"gamlss"}, or \code{"bamlss"})
#' @param y numeric vector of observed responses to evaluate coverage
#'   against
#' @param newdata optional data frame of covariate values at which to
#'   compute predictions.  If \code{NULL} (the default), in-sample
#'   predictions are used.
#' @param level nominal coverage level (default 0.95)
#' @param nsim number of parameter/response simulation draws (default
#'   2000)
#' @param method \code{"parametric"} (default) for parameter
#'   uncertainty, or \code{"plug-in"} for fixed-parameter prediction
#' @param use_shortest logical: use shortest (HDR) interval for
#'   discrete families instead of equal-tailed?  Default \code{FALSE}.
#' @param data_fit (gamlss method only) optional data frame passed as
#'   \code{data} to \code{predictAll()}; needed when the original
#'   fitting data cannot be recovered from the fit object.
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
  dtype <- c(1L, 6L, 4L, 5L)[m]
  is_discrete <- family %in% c("degpd", "zidegpd")

  ## choose rng function
  rfun <- switch(family,
    egpd     = regpd,
    degpd    = rdiscegpd,
    zidegpd  = rzidiscegpd,
    stop("predictive_coverage not implemented for egpd family '", family, "'")
  )

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

  if (m == 1) {
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

  ## Clamp parameter draws to valid ranges.
  ## MVN draws on the link scale can produce extreme values after

  ## back-transformation (Inf, 0, negative).
  eps <- 1e-12
  for (nm in c("sigma", "xi", "kappa", "delta")) {
    if (!is.null(par_list[[nm]])) {
      v <- par_list[[nm]]
      v[!is.finite(v) | v <= 0] <- eps
      par_list[[nm]] <- v
    }
  }
  if (!is.null(par_list$prob)) {
    v <- par_list$prob
    v[!is.finite(v)] <- 0.5
    v <- pmin(pmax(v, eps), 1 - eps)
    par_list$prob <- v
  }
  if (!is.null(par_list$pi)) {
    v <- par_list$pi
    v[!is.finite(v)] <- 0.5
    v <- pmin(pmax(v, eps), 1 - eps)
    par_list$pi <- v
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

## Internal helper for gamlss link extraction
.get_link_name_gamlss <- function(fit, p) {
  val <- tryCatch(fit[[paste0(p, ".link")]], error = function(e) NULL)
  if (is.character(val) && length(val) == 1) return(val)
  fam <- fit$family
  if (!is.null(fam) && !is.function(fam)) {
    val <- tryCatch(fam[[paste0(p, ".link")]], error = function(e) NULL)
    if (is.character(val) && length(val) == 1) return(val)
  }
  stop("Could not determine the '", p, "' link for this gamlss fit.")
}

#' @rdname predictive_coverage
#' @export
predictive_coverage.gamlss <- function(fit, y, newdata = NULL, level = 0.95,
                                       nsim = 2000,
                                       method = c("parametric", "plug-in"),
                                       use_shortest = FALSE,
                                       data_fit = NULL, ...) {
  method <- match.arg(method)

  fam_obj  <- fit$family
  fam_name <- tryCatch(fam_obj$family, error = function(e) NULL)
  if (is.null(fam_name)) fam_name <- as.character(fam_obj)[1]

  ## Locate RNG function in the egpd namespace (rDEGPD1, rEGPD1, etc.)
  rng_name <- paste0("r", fam_name)
  egpd_ns <- asNamespace("egpd")
  if (!exists(rng_name, where = egpd_ns, inherits = FALSE)) {
    ## Fallback: try gamlss.dist
    if (requireNamespace("gamlss.dist", quietly = TRUE) &&
        exists(rng_name, where = asNamespace("gamlss.dist"), inherits = FALSE)) {
      rng_fun <- get(rng_name, envir = asNamespace("gamlss.dist"))
    } else {
      stop("No RNG found for family '", fam_name, "' (looked in egpd and gamlss.dist).")
    }
  } else {
    rng_fun <- get(rng_name, envir = egpd_ns)
  }

  ## Determine if discrete/continuous from family name
  is_discrete <- grepl("^(D|ZID)", fam_name, ignore.case = FALSE)

  params <- intersect(c("mu", "sigma", "nu", "tau"), fit$parameters)
  if (!length(params)) stop("Could not detect mu/sigma/nu/tau parameters.")

  n <- if (!is.null(newdata)) nrow(newdata) else length(y)
  if (length(y) != n) stop("length(y) must match the number of observations.")

  link_names <- vapply(params, function(p) .get_link_name_gamlss(fit, p), character(1))
  invlink <- lapply(link_names, function(ln) {
    ml <- tryCatch(make.link(ln), error = function(e) NULL)
    if (!is.null(ml)) return(ml$linkinv)
    ## gamlss-specific links via gamlss package
    if (requireNamespace("gamlss", quietly = TRUE)) {
      mlg_fn <- utils::getFromNamespace("make.link.gamlss", "gamlss")
      ml <- tryCatch(mlg_fn(ln), error = function(e) NULL)
      if (!is.null(ml)) return(ml$linkinv)
    }
    stop("Cannot find inverse link for '", ln, "'")
  })
  names(invlink) <- params

  ## Per-parameter predict helper.
  ## gamlss::predictAll() fails with custom families (e.g. DEGPD1), so we
  ## use per-parameter predict() calls instead.
  ## suppressWarnings: gamlss emits "discrepancy between the original and
  ## the re-fit" warnings for out-of-sample predictions; these are expected.
  .predict_param <- function(p, type = "response", se.fit = FALSE) {
    args <- list(object = fit, what = p, type = type, se.fit = se.fit)
    if (!is.null(newdata)) args$newdata <- newdata
    if (!is.null(data_fit)) args$data <- data_fit
    suppressWarnings(
      tryCatch(do.call(predict, args),
        error = function(e) {
          ## Retry without data argument
          args$data <- NULL
          do.call(predict, args)
        }
      )
    )
  }

  ## Parameter draws on response scale (nsim x n)
  par_draws <- switch(method,
    "plug-in" = {
      out <- list()
      for (p in params) {
        v <- as.numeric(.predict_param(p, type = "response", se.fit = FALSE))
        if (length(v) != n) v <- rep(v, length.out = n)
        out[[p]] <- matrix(rep(v, each = nsim), nrow = nsim, ncol = n)
      }
      out
    },
    "parametric" = {
      draws <- list()
      for (p in params) {
        pred <- .predict_param(p, type = "link", se.fit = TRUE)
        if (is.list(pred) && "fit" %in% names(pred)) {
          m <- as.numeric(pred$fit)
          s <- as.numeric(pred$se.fit)
        } else {
          ## Fallback: no SE available
          m <- as.numeric(pred)
          s <- rep(0, length(m))
        }
        if (length(m) != n) m <- rep(m, length.out = n)
        if (length(s) != n) s <- rep(s, length.out = n)

        m[!is.finite(m)] <- 0
        s[!is.finite(s) | s < 0] <- 0

        eta <- matrix(
          rnorm(nsim * n, mean = rep(m, each = nsim), sd = rep(s, each = nsim)),
          nrow = nsim, ncol = n
        )
        val <- invlink[[p]](eta)
        eps <- 1e-12
        ln <- link_names[[p]]
        if (grepl("^log", ln, ignore.case = TRUE)) {
          val[!is.finite(val) | val <= 0] <- eps
        } else if (grepl("logit|cloglog|probit|cauchit|loglog", ln, ignore.case = TRUE)) {
          val[!is.finite(val)] <- 0.5
          val <- pmin(pmax(val, eps), 1 - eps)
        } else {
          bad <- !is.finite(val)
          if (any(bad)) {
            col_med <- apply(val, 2, function(col) {
              ok <- is.finite(col); if (any(ok)) stats::median(col[ok]) else 0
            })
            val[bad] <- rep(col_med, each = nsim)[bad]
          }
        }
        draws[[p]] <- val
      }
      draws
    }
  )

  ## Simulate responses
  ysim <- matrix(NA_real_, nsim, n)
  for (j in seq_len(n)) {
    args <- list(n = nsim)
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
#' @export
predictive_coverage.bamlss <- function(fit, y, newdata = NULL, level = 0.95,
                                       nsim = 2000,
                                       method = c("parametric", "plug-in"),
                                       use_shortest = FALSE, ...) {
  method <- match.arg(method)

  fam <- fit$family
  fam_name   <- fam$family        # e.g. "egpd1", "degpd1", "zidegpd1"
  param_names <- fam$names        # e.g. c("sigma","xi","kappa")
  links       <- fam$links        # e.g. c(sigma="log", xi="log", kappa="log")

  ## Determine family type, rng function, and discrete flag
  base_fam <- sub("[0-9]+$", "", fam_name)   # "egpd", "degpd", "zidegpd", "ziegpd"
  m_val <- as.integer(sub("^[a-z]+", "", fam_name))
  dtype <- c(1L, 6L, 4L, 5L)[m_val]
  is_discrete <- base_fam %in% c("degpd", "zidegpd")

  rfun <- switch(base_fam,
    egpd     = regpd,
    degpd    = rdiscegpd,
    zidegpd  = rzidiscegpd,
    ziegpd   = rziegpd,
    stop("predictive_coverage not implemented for bamlss family '", fam_name, "'")
  )

  n <- if (!is.null(newdata)) nrow(newdata) else length(y)
  if (length(y) != n) stop("length(y) must match the number of observations.")

  ## Build inverse-link functions from the link names
  make_invlink <- function(ln) {
    switch(ln,
      "log"     = exp,
      "logit"   = function(eta) 1 / (1 + exp(-eta)),
      "probit"  = pnorm,
      "cloglog" = function(eta) 1 - exp(-exp(eta)),
      "identity" = function(eta) eta,
      {
        ml <- tryCatch(make.link(ln), error = function(e) NULL)
        if (!is.null(ml)) ml$linkinv else stop("Unknown link '", ln, "'")
      }
    )
  }
  invlinks <- lapply(links, make_invlink)
  names(invlinks) <- param_names

  ## Build link functions for forward transform
  make_linkfun <- function(ln) {
    switch(ln,
      "log"     = log,
      "logit"   = function(mu) log(mu / (1 - mu)),
      "probit"  = qnorm,
      "cloglog" = function(mu) log(-log(1 - mu)),
      "identity" = function(mu) mu,
      {
        ml <- tryCatch(make.link(ln), error = function(e) NULL)
        if (!is.null(ml)) ml$linkfun else stop("Unknown link '", ln, "'")
      }
    )
  }
  linkfuns <- lapply(links, make_linkfun)
  names(linkfuns) <- param_names

  if (method == "plug-in") {
    ## predict(type="parameter") returns named list of vectors on response scale
    pars <- predict(fit, newdata = newdata, type = "parameter")
    par_mat <- lapply(param_names, function(nm) {
      v <- as.numeric(pars[[nm]])
      matrix(rep(v, each = nsim), nrow = nsim, ncol = n)
    })
    names(par_mat) <- param_names
  } else {
    ## Parametric: get response-scale predictions, transform to link,
    ## estimate SEs from model coefficients, perturb, and back-transform
    pars <- predict(fit, newdata = newdata, type = "parameter")

    par_mat <- list()
    for (nm in param_names) {
      mu_resp <- as.numeric(pars[[nm]])
      mu_link <- linkfuns[[nm]](mu_resp)

      ## Estimate SE on link scale from the model
      ## bamlss stores samples in fit$samples; use their SD as SE proxy
      se_link <- rep(0, n)
      if (!is.null(fit$samples) && !is.null(fit$samples[[nm]])) {
        samp <- fit$samples[[nm]]
        if (is.list(samp)) {
          ## Combine all coefficient samples (intercept + smooth terms)
          all_samps <- do.call(cbind, samp)
          if (!is.null(all_samps) && ncol(all_samps) > 0) {
            ## Get design matrix for this parameter if possible
            ## Use SD of fitted values from MCMC samples as link-scale SE
            coef_sd <- apply(all_samps, 2, sd)
            ## Simple approximation: SE = sqrt(sum(coef_var))
            ## More accurate would use design matrix, but this is reasonable
            se_link <- rep(sqrt(sum(coef_sd^2)), n)
          }
        }
      }

      ## If no samples available, use a small SE based on prediction
      if (all(se_link == 0)) {
        ## Use coefficient vcov if available
        ## Otherwise fall back to a fraction of the absolute link-scale value
        se_link <- rep(abs(mean(mu_link)) * 0.05, n)
        se_link[se_link < 1e-6] <- 1e-6
      }

      ## Perturb on link scale and back-transform
      eta <- matrix(
        rnorm(nsim * n, mean = rep(mu_link, each = nsim),
              sd = rep(se_link, each = nsim)),
        nrow = nsim, ncol = n
      )
      val <- invlinks[[nm]](eta)

      ## Enforce parameter support
      eps <- 1e-12
      ln <- links[[nm]]
      if (ln == "log") {
        val[!is.finite(val) | val <= 0] <- eps
      } else if (ln %in% c("logit", "probit", "cloglog")) {
        val[!is.finite(val)] <- 0.5
        val <- pmin(pmax(val, eps), 1 - eps)
      } else {
        bad <- !is.finite(val)
        if (any(bad)) {
          col_med <- apply(val, 2, function(col) {
            ok <- is.finite(col); if (any(ok)) stats::median(col[ok]) else 0
          })
          val[bad] <- rep(col_med, each = nsim)[bad]
        }
      }
      par_mat[[nm]] <- val
    }
  }

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
