## Neural Bayes inference for bivariate BEGPD
## Requires Julia + JuliaConnectoR + NeuralEstimators (Suggests only)

## Package-level environment to track Julia session state
.begpd_env <- new.env(parent = emptyenv())
.begpd_env$julia_initialized <- FALSE
.begpd_env$npe <- NULL
.begpd_env$nbe <- NULL

.julia_let <- function(...) JuliaConnectoR::juliaLet(...)
.julia_eval <- function(...) JuliaConnectoR::juliaEval(...)
.neural_train <- function(...) NeuralEstimators::train(...)
## Flux state + BSON via the egpd Julia helpers: the same file format as
## NeuralEstimators::savestate()/loadstate() for Flux estimators, without their
## LuxEstimator check (undefined in NeuralEstimators.jl 0.1.x).
.neural_savestate <- function(estimator, filename) {
  .julia_let("egpd_savestate(estimator, filename)", estimator = estimator,
             filename = normalizePath(filename, mustWork = FALSE))
}
.neural_loadstate <- function(estimator, filename) {
  .julia_let("egpd_loadstate(estimator, filename)", estimator = estimator,
             filename = normalizePath(filename, mustWork = TRUE))
}
## Learned log likelihood-to-evidence ratio log r(Z, theta), one value per
## column of theta (network parameter scale).
.neural_logratio <- function(estimator, Z, theta) {
  .julia_let("egpd_logratio(estimator, Z, theta)", estimator = estimator,
             Z = Z, theta = theta)
}
## Julia, NeuralEstimators.jl and Flux.jl versions of the running session
## (NA when Julia is unavailable).
.neural_backend_versions <- function() {
  jl <- tryCatch(as.character(JuliaConnectoR::juliaEval("egpd_backend_versions()")),
                 error = function(e) rep(NA_character_, 3L))
  if (length(jl) != 3L) jl <- rep(NA_character_, 3L)
  jl[jl %in% "unknown"] <- NA_character_
  r_ne <- tryCatch(as.character(utils::packageVersion("NeuralEstimators")),
                   error = function(e) NA_character_)
  c(julia = jl[1L], NeuralEstimators.jl = jl[2L], Flux.jl = jl[3L],
    NeuralEstimators.R = r_ne,
    egpd = tryCatch(as.character(utils::packageVersion("egpd")),
                    error = function(e) NA_character_))
}

## Worker count: 1 on Windows; detectCores() - 1 when unset or invalid; at most
## 2 when _R_CHECK_LIMIT_CORES_ is set (R CMD check --as-cran), using the same
## test as parallel, which otherwise refuses to spawn more processes.
.sanitize_mc_cores <- function(mc.cores) {
  if (.Platform$OS.type == "windows") return(1L)
  if (is.null(mc.cores) || length(mc.cores) != 1L ||
      is.na(mc.cores) || !is.finite(mc.cores) || mc.cores < 1L) {
    dc <- parallel::detectCores()
    if (is.na(dc) || dc < 1L) dc <- 1L
    mc.cores <- as.integer(dc) - 1L
  }
  n <- max(1L, as.integer(mc.cores))
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(chk) && chk != "false") n <- min(n, 2L)
  n
}

.neural_seed <- function(seed) {
  if (!is.null(seed))
    .julia_let("Random.seed!(seed)", seed = as.integer(seed))
  invisible(NULL)
}

.neural_simulate <- function(indices, FUN, mc.cores) {
  mc.cores <- .sanitize_mc_cores(mc.cores)
  seeds <- sample.int(.Machine$integer.max, length(indices))
  parent_seed <- .Random.seed
  on.exit(assign(".Random.seed", parent_seed, envir = .GlobalEnv))
  parallel::mclapply(seq_along(indices), function(i) {
    set.seed(seeds[i])
    FUN(indices[i])
  }, mc.cores = mc.cores, mc.set.seed = FALSE)
}

.neural_model_metadata <- function(family, estimator, data_dim = 2L) {
  cfg <- .neuralbayes_config(family, data_dim)
  bounds <- if (family %in% c("mdgpd", "zimdgpd")) {
    rbind(sigma = c(0.1, 10), xi = c(-0.4, 0.5), rate = c(0.01, 5))
  } else {
    rbind(kappa = c(0.1, 10), sigma = c(0.1, 3), xi = c(0.01, 0.5),
          thL = c(0.1, 20), thU = c(0.1, 20), thw = c(0.01, 0.49))
  }
  if (family %in% c("bzidegpd", "zimdgpd"))
    bounds <- rbind(bounds, pi0 = c(0.01, 0.9))
  param_transform <- if (family %in% c("mdgpd", "zimdgpd")) {
    ifelse(cfg$parnames == "pi0", "logit",
           ifelse(cfg$parnames == "xi", "identity", "log"))
  } else {
    ifelse(cfg$parnames == "pi0", "logit", "log")
  }
  list(contract = 3L, family = family, estimator = estimator,
       data_dim = cfg$data_dim, parnames = cfg$parnames,
       prior_bounds = bounds,
       parameter_transform = param_transform,
       architecture = if (estimator == "npe") "public-ratio-v1" else "point-v1",
       simulator = if (family %in% c("mdgpd", "zimdgpd"))
         "akn-definition-2.1-v3" else "begpd-radial-v1")
}

.write_neural_metadata <- function(path, family, estimator, data_dim = 2L,
                                   training = NULL) {
  metadata <- .neural_model_metadata(family, estimator, data_dim)
  ## Not part of the compatibility contract: provenance used for warnings.
  metadata$training <- training
  metadata$versions <- .neural_backend_versions()
  saveRDS(metadata, paste0(path, ".metadata.rds"))
}

.neural_contract_matches <- function(metadata, family, estimator, data_dim) {
  expected <- .neural_model_metadata(family, estimator, data_dim)
  is.list(metadata) && all(names(expected) %in% names(metadata)) &&
    identical(metadata[names(expected)], expected)
}

.validate_neural_model <- function(path, family, estimator, data_dim) {
  sidecar <- paste0(path, ".metadata.rds")
  metadata <- if (file.exists(sidecar))
    tryCatch(readRDS(sidecar), error = function(e) NULL) else NULL
  if (!.neural_contract_matches(metadata, family, estimator, data_dim))
    stop("Incompatible or unversioned neural asset: ", path,
         ". Retrain with train_begpd(), train_bdegpd(), or train_mdgpd(); ",
         "no bundled .bson models are shipped, and weights without a matching ",
         "'.bson.metadata.rds' sidecar cannot be relabeled as calibrated for this contract.",
         call. = FALSE)
  invisible(metadata)
}

## Warn when the data or the session fall outside what the network saw in
## training: sample size, data support and NeuralEstimators.jl version.
.neural_check_training_range <- function(metadata, n) {
  m <- metadata$training$m
  if (is.numeric(m) && length(m) >= 1L && all(is.finite(m)) &&
      (n < min(m) || n > max(m))) {
    warning(sprintf(paste0(
      "n = %d observations lies outside the training sample-size range ",
      "[%d, %d] of this neural model; amortised estimates are extrapolations."),
      as.integer(n), as.integer(min(m)), as.integer(max(m))), call. = FALSE)
  }
  invisible(NULL)
}

.neural_check_versions <- function(metadata, current = .neural_backend_versions()) {
  rec <- metadata$versions
  if (is.null(rec)) return(invisible(NULL))
  for (nm in c("NeuralEstimators.jl", "Flux.jl")) {
    a <- unname(rec[nm])
    b <- unname(current[nm])
    if (length(a) == 1L && length(b) == 1L && !is.na(a) && !is.na(b) &&
        !identical(a, b)) {
      warning(sprintf(paste0(
        "Neural model was trained with %s %s but %s is loaded; retrain or ",
        "load the training version if estimates look wrong."), nm, a, b),
        call. = FALSE)
    }
  }
  invisible(NULL)
}

.neural_check_support <- function(x, family) {
  bad <- switch(family,
    "begpd" = rowSums(x < 0) > 0,
    "bdegpd" = , "bzidegpd" = rowSums(x < 0 | x != round(x)) > 0,
    "mdgpd" = rowSums(x != round(x)) > 0 | apply(x, 1L, max) < 1,
    "zimdgpd" = rowSums(x != round(x)) > 0 |
      (apply(x, 1L, max) < 1 & rowSums(x != 0) > 0),
    rep(FALSE, nrow(x))
  )
  rule <- switch(family,
    "begpd" = "non-negative values",
    "bdegpd" = , "bzidegpd" = "non-negative integers",
    "mdgpd" = "integer rows with a strictly positive maximum",
    "zimdgpd" = "integer rows that are all zero or have a strictly positive maximum",
    "")
  if (any(bad)) {
    warning(sprintf(paste0(
      "%d of %d rows of 'x' lie outside the support of the %s model the ",
      "network was trained on (%s); neural estimates for such data are ",
      "extrapolations."), sum(bad), nrow(x), family, rule), call. = FALSE)
  }
  invisible(!any(bad))
}

#' Check Julia dependencies are installed
#' @noRd
.check_julia_deps <- function() {
  if (!requireNamespace("JuliaConnectoR", quietly = TRUE))
    stop("Package 'JuliaConnectoR' is required for method='neuralbayes'.\n",
         "Install it with: install.packages('JuliaConnectoR')",
         call. = FALSE)
  if (!requireNamespace("NeuralEstimators", quietly = TRUE))
    stop("Package 'NeuralEstimators' is required for method='neuralbayes'.\n",
         "Install it with: remotes::install_github('msainsburydale/NeuralEstimators')",
         call. = FALSE)
}

#' One-time Julia session initialization
#' @noRd
.init_julia_begpd <- function() {
  if (.begpd_env$julia_initialized) return(invisible(NULL))

  .check_julia_deps()

  ## Load Julia packages: BSON for the egpd_savestate/egpd_loadstate helpers, and
  ## Optimisers because NeuralEstimators::train() runs `using Optimisers: Adam`
  ## and then calls `Optimisers.Adam(...)`, which needs the module bound in Main.
  ## Each must be a direct dependency of the active Julia environment: a package
  ## that is only a dependency of another (BSON of NeuralEstimators.jl 0.1,
  ## Optimisers of Flux) cannot be loaded by `using`. Loading them here makes a
  ## missing one fail before any training data are simulated.
  tryCatch(.julia_eval("using NeuralEstimators, Flux, Random, BSON, Optimisers"),
           error = function(e)
             stop("Loading the Julia packages NeuralEstimators, Flux, BSON and ",
                  "Optimisers failed: ",
                  conditionMessage(e), "\nInstall them in Julia (>= 1.11) with:\n",
                  "  using Pkg; Pkg.add([\"NeuralEstimators\", \"Flux\", \"BSON\", ",
                  "\"Optimisers\"])",
                  call. = FALSE))

  ## Load architecture definition
  arch_file <- system.file("julia", "begpd_architecture.jl", package = "egpd")
  if (arch_file == "")
    stop("Cannot find inst/julia/begpd_architecture.jl. ",
         "Is the egpd package installed correctly?", call. = FALSE)

  arch_code <- paste(readLines(arch_file), collapse = "\n")
  JuliaConnectoR::juliaEval(arch_code)

  .begpd_env$julia_initialized <- TRUE
  invisible(NULL)
}

#' Initialize NPE estimator object in Julia
#' @param d integer; number of parameters (default 6L for BEGPD/BDEGPD, 7L for BZIDEGPD).
#' @param n integer; data dimension (default 2L for bivariate).
#' @noRd
.init_npe_begpd <- function(d = 6L, n = 2L) {
  npe <- .julia_let("initialize_ratio_network(n, d)", d = d, n = n)

  npe
}

#' Initialize NBE estimator object in Julia
#' @param d integer; number of parameters (default 6L for BEGPD/BDEGPD, 7L for BZIDEGPD).
#' @param n integer; data dimension (default 2L for bivariate).
#' @noRd
.init_nbe_begpd <- function(d = 6L, n = 2L) {
  J <- 5L  # ensemble components

  nbe <- .julia_let("
    NeuralEstimators.Ensemble([
      NeuralEstimators.PointEstimator(initialize_point_network(n, d))
      for j in 1:J
    ])
  ", n = n, d = d, J = J)

  nbe
}

#' Variance-stabilizing signed log transform
#' @noRd
.begpd_signed_log <- function(x) {
  sign(x) * log1p(abs(x)) - 1
}

## BEGPD parameter names (in order matching the prior sampler)
.begpd_parnames <- c("kappa", "sigma", "xi", "thL", "thU", "thw")

#' Configuration for each neural-Bayes family
#' @param family character; one of "begpd", "bdegpd", "bzidegpd",
#'   "mdgpd", "zimdgpd".
#' @param data_dim integer; data dimension (default 2L). Only affects MDGPD
#'   families, where a dimension-specific file prefix is used (e.g.,
#'   "MDGPD_" for 2D, "MDGPD3D_" for 3D).
#' @return Named list with d, data_dim, file_prefix, parnames, backtransform.
#' @noRd
.neuralbayes_config <- function(family, data_dim = 2L) {
  ## Dimension-specific file prefix for MDGPD families
  .mdgpd_prefix <- function(base, data_dim) {
    if (data_dim == 2L) base else paste0(base, data_dim, "D_")
  }

  switch(family,
    "begpd" = list(
      d           = 6L,
      data_dim    = 2L,
      file_prefix = "",
      parnames    = c("kappa", "sigma", "xi", "thL", "thU", "thw"),
      backtransform = function(mat) exp(mat)
    ),
    "bdegpd" = list(
      d           = 6L,
      data_dim    = 2L,
      file_prefix = "BDEGPD_",
      parnames    = c("kappa", "sigma", "xi", "thL", "thU", "thw"),
      backtransform = function(mat) exp(mat)
    ),
    "bzidegpd" = list(
      d           = 7L,
      data_dim    = 2L,
      file_prefix = "BZIDEGPD_",
      parnames    = c("kappa", "sigma", "xi", "thL", "thU", "thw", "pi0"),
      backtransform = function(mat) {
        mat[1:6, ] <- exp(mat[1:6, , drop = FALSE])
        mat[7, ]   <- plogis(mat[7, ])
        mat
      }
    ),
    "mdgpd" = list(
      d           = 3L,
      data_dim    = as.integer(data_dim),
      file_prefix = .mdgpd_prefix("MDGPD_", data_dim),
      parnames    = c("sigma", "xi", "rate"),
      backtransform = function(mat) {
        mat[c(1L, 3L), ] <- exp(mat[c(1L, 3L), , drop = FALSE])
        mat
      }
    ),
    "zimdgpd" = list(
      d           = 4L,
      data_dim    = as.integer(data_dim),
      file_prefix = .mdgpd_prefix("ZIMDGPD_", data_dim),
      parnames    = c("sigma", "xi", "rate", "pi0"),
      backtransform = function(mat) {
        mat[c(1L, 3L), ] <- exp(mat[c(1L, 3L), , drop = FALSE])
        mat[4L, ] <- plogis(mat[4L, ])
        mat
      }
    ),
    stop("Unknown neural-Bayes family: ", family, call. = FALSE)
  )
}

#' Draw candidate parameter vectors in the transformed parameter space
#' @noRd
.sample_theta_prior <- function(family, K, data_dim = 2L) {
  K <- as.integer(K)
  out <- switch(family,
    "begpd" = rbind(
      log(runif(K, 0.1, 10)),
      log(runif(K, 0.1, 3)),
      log(runif(K, 0.01, 0.5)),
      log(runif(K, 0.1, 20)),
      log(runif(K, 0.1, 20)),
      log(runif(K, 0.01, 0.49))
    ),
    "bdegpd" = rbind(
      log(runif(K, 0.1, 10)),
      log(runif(K, 0.1, 3)),
      log(runif(K, 0.01, 0.5)),
      log(runif(K, 0.1, 20)),
      log(runif(K, 0.1, 20)),
      log(runif(K, 0.01, 0.49))
    ),
    "bzidegpd" = rbind(
      log(runif(K, 0.1, 10)),
      log(runif(K, 0.1, 3)),
      log(runif(K, 0.01, 0.5)),
      log(runif(K, 0.1, 20)),
      log(runif(K, 0.1, 20)),
      log(runif(K, 0.01, 0.49)),
      qlogis(runif(K, 0.01, 0.9))
    ),
    "mdgpd" = rbind(
      log(runif(K, 0.1, 10)),
      runif(K, -0.4, 0.5),
      log(runif(K, 0.01, 5))
    ),
    "zimdgpd" = rbind(
      log(runif(K, 0.1, 10)),
      runif(K, -0.4, 0.5),
      log(runif(K, 0.01, 5)),
      qlogis(runif(K, 0.01, 0.9))
    ),
    stop("Unknown neural-Bayes family: ", family, call. = FALSE)
  )

  storage.mode(out) <- "double"
  out
}

## Prior support on the network (transformed) parameter scale and the
## natural-scale uniform bounds/transforms it comes from (model contract).
.neural_prior_box <- function(family, data_dim = 2L) {
  md <- .neural_model_metadata(family, "npe", data_dim)
  fwd <- function(x, tr) switch(tr, log = log(x), identity = x, logit = qlogis(x))
  b <- md$prior_bounds
  tr <- unname(md$parameter_transform)
  list(lower = vapply(seq_len(nrow(b)), function(j) fwd(b[j, 1], tr[j]), 0),
       upper = vapply(seq_len(nrow(b)), function(j) fwd(b[j, 2], tr[j]), 0),
       bounds = b, transform = tr)
}

## Log prior density on the network scale of theta_j = g(u_j), u_j ~ U(a_j, b_j):
## log: theta - log(b - a); identity: -log(b - a);
## logit: log plogis(theta) + log plogis(-theta) - log(b - a). -Inf off the box.
.neural_log_prior <- function(theta, box) {
  theta <- as.matrix(theta)
  out <- numeric(ncol(theta))
  for (j in seq_len(nrow(theta))) {
    t <- theta[j, ]
    lw <- log(box$bounds[j, 2] - box$bounds[j, 1])
    lp <- switch(box$transform[j],
      log = t - lw,
      identity = rep(-lw, length(t)),
      logit = plogis(t, log.p = TRUE) + plogis(-t, log.p = TRUE) - lw
    )
    lp[!(t >= box$lower[j] & t <= box$upper[j])] <- -Inf
    out <- out + lp
  }
  unname(out)
}

.log_weights_ess <- function(lw) {
  ok <- is.finite(lw)
  if (!any(ok)) return(0)
  w <- exp(lw[ok] - max(lw[ok]))
  sum(w)^2 / sum(w^2)
}

.normalize_log_weights <- function(lw) {
  w <- numeric(length(lw))
  ok <- is.finite(lw)
  w[ok] <- exp(lw[ok] - max(lw[ok]))
  w / sum(w)
}

.mvt_draw <- function(n, mu, L, df) {
  d <- length(mu)
  z <- matrix(stats::rnorm(n * d), d, n)
  mu + (L %*% z) * rep(sqrt(df / stats::rchisq(n, df)), each = d)
}

## Multivariate t log density with location mu, scale L L^T (L lower triangular).
.mvt_logdens <- function(x, mu, L, df) {
  d <- length(mu)
  y <- forwardsolve(L, x - mu)
  lgamma((df + d) / 2) - lgamma(df / 2) - d / 2 * log(df * pi) -
    sum(log(diag(L))) - (df + d) / 2 * log1p(colSums(y^2) / df)
}

## Heavy-tailed proposal fitted to weighted particles, covariance inflated.
## With fewer than d + 1 effective particles the covariance is not estimable
## and a scale of 5% of the prior box is used instead.
.npe_fit_proposal <- function(theta, w, inflate, box) {
  d <- nrow(theta)
  mu <- as.vector(theta %*% w)
  xc <- theta - mu
  width2 <- (box$upper - box$lower)^2
  if (1 / sum(w^2) < d + 1) {
    Sigma <- diag(0.0025 * width2, d)
  } else {
    Sigma <- (xc * rep(w, each = d)) %*% t(xc) / max(1 - sum(w^2), 1e-8)
    Sigma <- inflate * (Sigma + t(Sigma)) / 2
  }
  diag(Sigma) <- diag(Sigma) + 1e-10 * width2
  L <- tryCatch(t(chol(Sigma)), error = function(e) NULL)
  if (is.null(L)) L <- diag(sqrt(diag(Sigma)), d)
  list(mu = mu, L = L)
}

## K draws from (1 - a) t_df(mu, L L^T) truncated to the prior box + a prior,
## with the exact mixture log density (truncation constant estimated by the
## acceptance rate of the t draws).
.npe_draw_proposal <- function(K, prop, df, defensive, box, family, data_dim) {
  d <- length(prop$mu)
  n_p <- max(1L, as.integer(round(defensive * K)))
  n_t <- K - n_p
  if (n_t < 1L) {
    ## a defensive weight this close to 1 leaves no draw for the fitted t
    theta <- .sample_theta_prior(family, K, data_dim)
    return(list(theta = theta, lq = .neural_log_prior(theta, box)))
  }
  kept <- matrix(0, d, 0L)
  drawn <- 0
  inside <- 0
  while (ncol(kept) < n_t && drawn < 200 * K) {
    z <- .mvt_draw(max(2000L, 2L * (n_t - ncol(kept))), prop$mu, prop$L, df)
    ok <- colSums(z >= box$lower & z <= box$upper) == d
    drawn <- drawn + ncol(z)
    inside <- inside + sum(ok)
    kept <- cbind(kept, z[, ok, drop = FALSE])
  }
  if (ncol(kept) < n_t) {
    ## The fitted t has (almost) no mass on the prior support: use the prior.
    theta <- .sample_theta_prior(family, K, data_dim)
    return(list(theta = theta, lq = .neural_log_prior(theta, box)))
  }
  theta <- cbind(kept[, seq_len(n_t), drop = FALSE],
                 .sample_theta_prior(family, n_p, data_dim))
  lt <- log(n_t / K) + .mvt_logdens(theta, prop$mu, prop$L, df) - log(inside / drawn)
  lp <- log(n_p / K) + .neural_log_prior(theta, box)
  m <- pmax(lt, lp)
  list(theta = theta, lq = m + log(exp(lt - m) + exp(lp - m)))
}

.npe_control <- function(control, nsamples, d) {
  defaults <- list(K = 5000L, K0 = 20000L,
                   ess_target = max(200L, min(as.integer(nsamples), 1000L)),
                   ess_floor = max(200L, 25L * d), max_stages = 30L, df = 5,
                   inflate = 1.5, defensive = 0.1)
  if (!is.list(control) || (length(control) && is.null(names(control))) ||
      !all(names(control) %in% names(defaults)))
    stop("'npe.control' must be a named list with entries in: ",
         paste(names(defaults), collapse = ", "), call. = FALSE)
  out <- utils::modifyList(defaults, control)
  num1 <- function(v) is.numeric(v) && length(v) == 1L && !is.na(v) && is.finite(v)
  whole <- function(v) num1(v) && v == round(v)
  rules <- list(
    K = list(whole(out$K) && out$K >= 1, "a positive integer"),
    K0 = list(whole(out$K0) && out$K0 >= 1, "a positive integer"),
    ess_target = list(num1(out$ess_target) && out$ess_target >= 1, "a number >= 1"),
    ess_floor = list(num1(out$ess_floor) && out$ess_floor >= 1, "a number >= 1"),
    max_stages = list(whole(out$max_stages) && out$max_stages >= 0,
                      "a non-negative integer"),
    df = list(num1(out$df) && out$df > 0, "a positive number"),
    inflate = list(num1(out$inflate) && out$inflate > 0, "a positive number"),
    ## the prior's share of the proposal mixture: at 1 no draw would come from the
    ## fitted proposal and its log density would be NaN
    defensive = list(num1(out$defensive) && out$defensive >= 0 && out$defensive < 1,
                     "a number in [0, 1)"))
  for (nm in names(rules))
    if (!rules[[nm]][[1L]])
      stop(sprintf("'npe.control$%s' must be %s", nm, rules[[nm]][[2L]]),
           call. = FALSE)
  out$K <- as.integer(out$K)
  out$K0 <- as.integer(max(out$K0, out$K))
  if (out$K < 2L * out$ess_floor)
    stop("'npe.control$K' must be at least twice 'ess_floor'", call. = FALSE)
  out
}

## Adaptive importance sampling of p(theta | Z) ~ p(theta) r(Z, theta) on the
## network parameter scale (D3). Stage 0 samples the prior; each later stage
## fits a multivariate-t proposal (inflated covariance, mixed with the prior)
## to the particles reweighted to the tempered target p(theta) r^beta, where
## beta is the largest value keeping ESS >= ess_floor, until the ESS of the
## full-posterior weights reaches ess_target. Returns systematic-resampling
## draws and diagnostics.
.npe_importance_sample <- function(logratio, family, data_dim, nsamples,
                                   control = list()) {
  box <- .neural_prior_box(family, data_dim)
  d <- length(box$lower)
  ctl <- .npe_control(control, nsamples, d)
  K <- ctl$K

  eval_lr <- function(theta) {
    lr <- as.numeric(logratio(theta))
    if (length(lr) != ncol(theta) || anyNA(lr) || any(lr == Inf))
      stop("The neural log-ratio must return one value per parameter vector ",
           "(finite or -Inf).", call. = FALSE)
    lr
  }
  tempered <- function(b, lp, lr, lq) lp + ifelse(is.finite(lr), b * lr, -Inf) - lq

  ## Stage 0 explores the prior with K0 >= K draws.
  theta <- .sample_theta_prior(family, ctl$K0, data_dim)
  lp <- .neural_log_prior(theta, box)
  lq <- lp
  lr <- eval_lr(theta)
  if (!any(is.finite(lr)))
    stop("No prior draw has a finite approximate likelihood; the data are ",
         "outside the range the network can represent.", call. = FALSE)

  beta <- 0
  best <- NULL
  path <- data.frame(stage = integer(), beta = numeric(), ess = numeric())
  for (stage in 0:ctl$max_stages) {
    lw <- tempered(1, lp, lr, lq)
    ess <- .log_weights_ess(lw)
    path[nrow(path) + 1L, ] <- list(stage, beta, ess)
    if (is.null(best) || ess > best$ess)
      best <- list(theta = theta, lw = lw, ess = ess, stage = stage)
    if (ess >= ctl$ess_target || stage == ctl$max_stages) break

    ess_b <- function(b) .log_weights_ess(tempered(b, lp, lr, lq))
    root <- function(lo, hi) stats::uniroot(function(b) ess_b(b) - ctl$ess_floor,
                                            c(lo, hi), tol = 1e-12)$root
    w <- NULL
    if (ess >= ctl$ess_floor) {
      beta <- 1
    } else if (ess_b(beta) >= ctl$ess_floor) {
      beta <- root(beta, 1)
    } else if (ess_b(0) >= ctl$ess_floor) {
      ## the proposal over-concentrated: step back to a wider tempered target
      beta <- root(0, beta)
    } else {
      ## too few particles with a finite target to weight reliably:
      ## cross-entropy step on the best particles
      ok <- which(is.finite(lr))
      top <- ok[order((lp + lr)[ok], decreasing = TRUE)][
        seq_len(min(length(ok), ctl$ess_floor))]
      w <- numeric(ncol(theta))
      w[top] <- 1 / length(top)
    }
    if (is.null(w)) w <- .normalize_log_weights(tempered(beta, lp, lr, lq))
    prop <- .npe_fit_proposal(theta, w, ctl$inflate, box)
    new <- .npe_draw_proposal(K, prop, ctl$df, ctl$defensive, box, family, data_dim)
    theta <- new$theta
    lq <- new$lq
    lp <- .neural_log_prior(theta, box)
    lr <- eval_lr(theta)
  }

  w <- .normalize_log_weights(best$lw)
  u <- (stats::runif(1) + seq_len(nsamples) - 1) / nsamples
  idx <- pmin(findInterval(u, cumsum(w)) + 1L, length(w))
  list(theta = best$theta[, idx, drop = FALSE],
       diagnostics = list(ess = best$ess, ess_target = ctl$ess_target,
                          K = ncol(best$theta), stages = best$stage,
                          unique = length(unique(idx)), nsamples = nsamples,
                          path = path))
}

#' Internal: fit multivariate families via neural Bayes estimation
#'
#' @param x n-by-d matrix or data.frame of observations (d >= 2).
#' @param family character; "begpd", "bdegpd", "bzidegpd", "mdgpd", or "zimdgpd".
#' @param model.path path to .bson model file, or NULL for bundled default
#' @param estimator "npe" or "nbe"
#' @param nsamples number of posterior samples (NPE only)
#' @param call the matched call from fitegpd()
#' @param npe.control optional named list for the NPE adaptive importance
#'   sampler (K, K0, ess_target, ess_floor, max_stages, df, inflate,
#'   defensive).
#' @noRd
.distfit_neuralbayes <- function(x, family = "begpd",
                                  model.path = NULL,
                                  estimator = c("npe", "nbe"),
                                  nsamples = 1000L,
                                  call = NULL,
                                  npe.control = list()) {
  estimator <- match.arg(estimator)
  if (length(nsamples) != 1L || !is.finite(nsamples) ||
      nsamples < 1 || nsamples != floor(nsamples))
    stop("'nsamples' must be a positive integer")

  ## MDGPD families allow d >= 2 columns; others require exactly 2
  mdgpd_families <- c("mdgpd", "zimdgpd")

  ## Validate input
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x))
    stop("For family='", family, "', 'x' must be a numeric matrix or data.frame",
         call. = FALSE)
  data_dim <- ncol(x)
  if (family %in% mdgpd_families) {
    if (data_dim < 2L)
      stop("For family='", family, "', 'x' must have at least 2 columns",
           call. = FALSE)
  } else {
    if (data_dim != 2L)
      stop("For family='", family, "', 'x' must have exactly 2 columns (bivariate data)",
           call. = FALSE)
  }
  if (nrow(x) < 2)
    stop("Need at least 2 observations", call. = FALSE)
  if (any(!is.finite(x))) stop("'x' must contain only finite observations")

  cfg <- .neuralbayes_config(family, data_dim = data_dim)
  n <- nrow(x)

  ## Resolve model path
  if (is.null(model.path)) {
    bson_name <- paste0(cfg$file_prefix, toupper(estimator), ".bson")
    model.path <- system.file("models", bson_name, package = "egpd")
    if (model.path == "") {
      train_fn <- switch(family,
        "begpd" = "train_begpd()", "bdegpd" = , "bzidegpd" = "train_bdegpd()",
        "mdgpd" = , "zimdgpd" = "train_mdgpd()", "train_begpd()")
      stop("No bundled neural model '", bson_name, "' is shipped with egpd. ",
           "Train (or retrain) compatible weights via ", train_fn,
           " (see train_begpd(), train_bdegpd(), train_mdgpd()) and pass the ",
           "saved '.bson' path via 'model.path' (retaining its '.bson.metadata.rds' sidecar).",
           call. = FALSE)
    }
  }
  if (!file.exists(model.path))
    stop("Model file not found: ", model.path, call. = FALSE)
  metadata <- .validate_neural_model(model.path, family, estimator, data_dim)
  .neural_check_support(x, family)
  .neural_check_training_range(metadata, n)

  ## Initialize Julia and architecture
  .check_julia_deps()
  .init_julia_begpd()
  .neural_check_versions(metadata)

  ## Create estimator and load trained state
  if (estimator == "npe") {
    est_obj <- .init_npe_begpd(d = cfg$d, n = data_dim)
  } else {
    est_obj <- .init_nbe_begpd(d = cfg$d, n = data_dim)
  }
  loaded_obj <- .neural_loadstate(est_obj, model.path)
  if (!is.null(loaded_obj)) est_obj <- loaded_obj

  ## Prepare data: transpose to data_dim x n and apply signed_log
  Z <- .begpd_signed_log(t(x))

  ## Inference
  importance <- NULL
  if (estimator == "npe") {
    ## NPE: adaptive importance sampling of prior x learned likelihood ratio
    is_out <- .npe_importance_sample(
      function(theta) .neural_logratio(est_obj, Z, theta),
      family = family, data_dim = data_dim, nsamples = as.integer(nsamples),
      control = npe.control
    )
    importance <- is_out$diagnostics
    if (importance$ess < 200)
      warning(sprintf(paste0(
        "NPE importance sampling reached an effective sample size of only ",
        "%.0f (< 200) after %d adaptive stages; posterior summaries are ",
        "unreliable."), importance$ess, importance$stages), call. = FALSE)
    if (importance$unique < nsamples / 10)
      warning(sprintf(paste0(
        "Only %d distinct parameter vectors among %d posterior draws; ",
        "posterior spread is underestimated."), importance$unique,
        as.integer(nsamples)), call. = FALSE)

    ## Back-transform from unconstrained parameter space
    post_mat <- cfg$backtransform(is_out$theta)
    rownames(post_mat) <- cfg$parnames

    ## Posterior summaries
    estimate <- apply(post_mat, 1, median)
    sd_est   <- apply(post_mat, 1, sd)
    names(estimate) <- cfg$parnames
    names(sd_est)   <- cfg$parnames
    posterior_samples <- post_mat

  } else {
    ## NBE: point estimates
    est_raw <- .julia_let(
      "Float64.(egpd_estimate(estimator, Z; batchsize = batchsize, use_gpu = false))",
      estimator = est_obj,
      Z = Z,
      batchsize = 32L
    )

    ## Result may be a list or matrix; extract and back-transform
    if (is.list(est_raw)) est_raw <- est_raw[[1]]
    raw_vec <- as.numeric(est_raw)
    if (length(raw_vec) != cfg$d || any(!is.finite(raw_vec)))
      stop("NBE returned an incompatible parameter vector")
    ## Back-transform: use a 1-column matrix for cfg$backtransform
    raw_mat <- matrix(raw_vec, ncol = 1)
    estimate <- cfg$backtransform(raw_mat)[, 1]
    names(estimate) <- cfg$parnames
    sd_est <- rep(NA_real_, cfg$d)
    names(sd_est) <- cfg$parnames
    posterior_samples <- NULL
  }

  ## Assemble fitegpd result
  structure(list(
    estimate          = estimate,
    sd                = sd_est,
    vcov              = NULL,
    loglik            = NA_real_,
    aic               = NA_real_,
    bic               = NA_real_,
    n                 = n,
    npar              = cfg$d,
    data              = x,
    type              = NA_integer_,
    family            = family,
    method            = "neuralbayes",
    fix.arg           = NULL,
    convergence       = 0L,
    optim             = NULL,
    call              = call,
    bernstein.m       = NULL,
    bernstein.weights = NULL,
    cpegpd.h          = NULL,
    ## BEGPD-specific fields
    estimator_type    = estimator,
    model_contract    = metadata,
    model.path        = model.path,
    posterior_samples  = posterior_samples,
    nsamples          = if (estimator == "npe") nsamples else NULL,
    importance_sampling = importance
  ), class = if (estimator == "npe") c("distfit_npe", "distfit") else "distfit")
}

.print_npe_diagnostics <- function(diag) {
  if (is.null(diag)) return(invisible(NULL))
  cat(sprintf(paste0(
    "Importance sampling: ESS = %.0f (target %d) from %d proposal draws after ",
    "%d adaptive stage%s; %d distinct of %d posterior draws\n"),
    diag$ess, as.integer(diag$ess_target), as.integer(diag$K),
    as.integer(diag$stages), if (diag$stages == 1L) "" else "s",
    as.integer(diag$unique), as.integer(diag$nsamples)))
  invisible(NULL)
}

#' @export
print.distfit_npe <- function(x, ...) {
  NextMethod()
  .print_npe_diagnostics(x$importance_sampling)
  invisible(x)
}

#' @export
summary.distfit_npe <- function(object, ...) {
  out <- NextMethod()
  out$importance_sampling <- object$importance_sampling
  class(out) <- c("summary.distfit_npe", class(out))
  out
}

#' @export
print.summary.distfit_npe <- function(x, ...) {
  NextMethod()
  .print_npe_diagnostics(x$importance_sampling)
  invisible(x)
}


#' 4-panel diagnostic plot for bivariate BEGPD fits
#' @noRd
.plot_begpd <- function(obj, ...) {
  est <- obj$estimate
  dat <- obj$data

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  ## Panel 1: Observed scatter
  plot(dat[, 1], dat[, 2], pch = ".", cex = 2,
       main = "Observed Data",
       xlab = expression(Y[1]), ylab = expression(Y[2]),
       col = adjustcolor("steelblue", 0.5))

  ## Panel 2: Simulated scatter from fitted model
  n_sim <- min(nrow(dat), 5000L)
  Y_sim <- rbegpd(n_sim,
                   kappa = est["kappa"], sigma = est["sigma"],
                   xi = est["xi"], thL = est["thL"],
                   thU = est["thU"], thw = est["thw"])
  plot(Y_sim[, 1], Y_sim[, 2], pch = ".", cex = 2,
       main = "Simulated from Fit",
       xlab = expression(Y[1]), ylab = expression(Y[2]),
       col = adjustcolor("firebrick", 0.5))

  ## Panel 3: Marginal Q-Q plot for radial component (Y1 + Y2)
  R_obs <- sort(dat[, 1] + dat[, 2])
  R_sim <- sort(Y_sim[, 1] + Y_sim[, 2])
  ## Interpolate to common length
  nn <- min(length(R_obs), length(R_sim))
  q_obs <- quantile(R_obs, probs = seq(0, 1, length.out = nn))
  q_sim <- quantile(R_sim, probs = seq(0, 1, length.out = nn))
  plot(q_sim, q_obs, pch = ".", cex = 2,
       main = "Q-Q Plot (Radial: Y1+Y2)",
       xlab = "Simulated Quantiles", ylab = "Observed Quantiles")
  abline(0, 1, col = "red", lwd = 2)

  ## Panel 4: Posterior densities (NPE) or parameter bar chart (NBE)
  if (!is.null(obj$posterior_samples)) {
    ## NPE: show posterior marginal densities for each parameter
    post <- obj$posterior_samples  # 6 x nsamples
    par_names <- rownames(post)
    cols <- c("steelblue", "firebrick", "forestgreen",
              "darkorange", "purple", "brown")
    ## Pre-compute densities to determine axis ranges
    dens_list <- list()
    max_dens <- 0
    for (i in seq_len(nrow(post))) {
      d <- density(post[i, ], n = 128)
      dens_list[[i]] <- d
      max_dens <- max(max_dens, max(d$y))
    }
    ## Single plot with all densities overlaid
    plot(0, 0, type = "n",
         xlim = range(post),
         ylim = c(0, max_dens * 1.1),
         main = "Posterior Marginals (NPE)",
         xlab = "Parameter value", ylab = "Density")
    for (i in seq_len(nrow(post))) {
      lines(dens_list[[i]], col = cols[i], lwd = 2)
    }
    legend("topright", legend = par_names, col = cols,
           lwd = 2, cex = 0.7, bg = "white")
  } else {
    ## NBE: bar chart of point estimates
    barplot(est, col = "steelblue", main = "Parameter Estimates (NBE)",
            ylab = "Estimate", las = 2)
  }

  invisible(obj)
}


#' 4-panel diagnostic plot for bivariate discrete EGPD fits (BDEGPD/BZIDEGPD)
#' @noRd
.plot_bdegpd <- function(obj, ...) {
  est <- obj$estimate
  dat <- obj$data
  family <- obj$family
  is_zi <- family == "bzidegpd"

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  ## Panel 1: Observed scatter with jitter (integers overlap on grid)
  plot(jitter(dat[, 1]), jitter(dat[, 2]), pch = ".", cex = 2,
       main = "Observed Data (jittered)",
       xlab = expression(Y[1]), ylab = expression(Y[2]),
       col = adjustcolor("steelblue", 0.5))

  ## Panel 2: Simulated scatter with jitter
  n_sim <- min(nrow(dat), 5000L)
  if (is_zi) {
    Y_sim <- rbzidegpd(n_sim,
                        kappa = est["kappa"], sigma = est["sigma"],
                        xi = est["xi"], thL = est["thL"],
                        thU = est["thU"], thw = est["thw"],
                        pi0 = est["pi0"])
  } else {
    Y_sim <- rbdegpd(n_sim,
                      kappa = est["kappa"], sigma = est["sigma"],
                      xi = est["xi"], thL = est["thL"],
                      thU = est["thU"], thw = est["thw"])
  }
  plot(jitter(Y_sim[, 1]), jitter(Y_sim[, 2]), pch = ".", cex = 2,
       main = "Simulated from Fit (jittered)",
       xlab = expression(Y[1]), ylab = expression(Y[2]),
       col = adjustcolor("firebrick", 0.5))

  ## Panel 3: Marginal barplot comparison (observed vs simulated for Y1+Y2)
  R_obs <- dat[, 1] + dat[, 2]
  R_sim <- Y_sim[, 1] + Y_sim[, 2]
  all_vals <- sort(unique(c(R_obs, R_sim)))
  ## Limit to a reasonable number of bars
  if (length(all_vals) > 30) all_vals <- all_vals[seq(1, length(all_vals), length.out = 30)]
  obs_counts <- tabulate(match(R_obs, all_vals), nbins = length(all_vals))
  sim_counts <- tabulate(match(R_sim, all_vals), nbins = length(all_vals))
  bar_mat <- rbind(Observed = obs_counts / length(R_obs),
                   Simulated = sim_counts / length(R_sim))
  colnames(bar_mat) <- all_vals
  barplot(bar_mat, beside = TRUE, col = c("steelblue", "firebrick"),
          main = "Marginal Y1+Y2", xlab = "Y1+Y2", ylab = "Proportion",
          legend.text = TRUE, args.legend = list(cex = 0.7, bg = "white"))

  ## Panel 4: Posterior densities (NPE) or parameter bar chart (NBE)
  if (!is.null(obj$posterior_samples)) {
    post <- obj$posterior_samples  # d x nsamples
    par_names <- rownames(post)
    cols <- c("steelblue", "firebrick", "forestgreen",
              "darkorange", "purple", "brown", "hotpink")[seq_len(nrow(post))]
    dens_list <- list()
    max_dens <- 0
    for (i in seq_len(nrow(post))) {
      d <- density(post[i, ], n = 128)
      dens_list[[i]] <- d
      max_dens <- max(max_dens, max(d$y))
    }
    plot(0, 0, type = "n",
         xlim = range(post),
         ylim = c(0, max_dens * 1.1),
         main = "Posterior Marginals (NPE)",
         xlab = "Parameter value", ylab = "Density")
    for (i in seq_len(nrow(post))) {
      lines(dens_list[[i]], col = cols[i], lwd = 2)
    }
    legend("topright", legend = par_names, col = cols,
           lwd = 2, cex = 0.7, bg = "white")
  } else {
    barplot(est, col = "steelblue", main = "Parameter Estimates (NBE)",
            ylab = "Estimate", las = 2)
  }

  invisible(obj)
}


#' Diagnostic plot for MDGPD/ZIMDGPD fits (d >= 2)
#' @noRd
.plot_mdgpd <- function(obj, ...) {
  est <- obj$estimate
  dat <- obj$data
  fam <- obj$family
  is_zi <- fam == "zimdgpd"
  data_dim <- ncol(dat)
  if (!.neural_contract_matches(obj$model_contract, fam, obj$estimator_type,
                                data_dim) ||
      !identical(names(est), .neuralbayes_config(fam, data_dim)$parnames))
    stop("Legacy MDGPD fit is incompatible with the new simulator; refit after retraining.",
         call. = FALSE)

  ## Simulate from fitted model
  n_sim <- min(nrow(dat), 5000L)
  if (is_zi) {
    Y_sim <- rzimdgpd(n_sim,
                       sigma = est["sigma"], xi = est["xi"],
                       rate = est["rate"],
                       pi0 = est["pi0"], d = data_dim)
  } else {
    Y_sim <- rmdgpd(n_sim,
                     sigma = est["sigma"], xi = est["xi"],
                     rate = est["rate"],
                     d = data_dim)
  }

  if (data_dim == 2L) {
    ## Bivariate: 4-panel layout
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    ## Panel 1: Observed scatter with jitter
    plot(jitter(dat[, 1]), jitter(dat[, 2]), pch = ".", cex = 2,
         main = "Observed Data (jittered)",
         xlab = expression(Y[1]), ylab = expression(Y[2]),
         col = adjustcolor("steelblue", 0.5))

    ## Panel 2: Simulated scatter with jitter
    plot(jitter(Y_sim[, 1]), jitter(Y_sim[, 2]), pch = ".", cex = 2,
         main = "Simulated from Fit (jittered)",
         xlab = expression(Y[1]), ylab = expression(Y[2]),
         col = adjustcolor("firebrick", 0.5))

    ## Panel 3: Marginal barplot comparison (Y1+Y2)
    R_obs <- dat[, 1] + dat[, 2]
    R_sim <- Y_sim[, 1] + Y_sim[, 2]
    all_vals <- sort(unique(c(R_obs, R_sim)))
    if (length(all_vals) > 30) all_vals <- all_vals[seq(1, length(all_vals), length.out = 30)]
    obs_counts <- tabulate(match(R_obs, all_vals), nbins = length(all_vals))
    sim_counts <- tabulate(match(R_sim, all_vals), nbins = length(all_vals))
    bar_mat <- rbind(Observed = obs_counts / length(R_obs),
                     Simulated = sim_counts / length(R_sim))
    colnames(bar_mat) <- all_vals
    barplot(bar_mat, beside = TRUE, col = c("steelblue", "firebrick"),
            main = "Marginal Y1+Y2", xlab = "Y1+Y2", ylab = "Proportion",
            legend.text = TRUE, args.legend = list(cex = 0.7, bg = "white"))

    ## Panel 4: Posterior densities or parameter bar chart
    .plot_mdgpd_posterior(obj, est)

  } else {
    ## d >= 3: pairs-style layout (observed pairwise scatters, simulated
    ## pairwise scatters, marginal sum barplot, and posterior/parameter panel)
    layout_n <- choose(data_dim, 2)
    n_cols <- min(layout_n, 3L)
    n_rows_data <- ceiling(layout_n / n_cols)
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mfrow = c(n_rows_data * 2L + 1L, n_cols),
        mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))

    ## Observed pairwise scatter plots
    for (i in seq_len(data_dim - 1L)) {
      for (j in (i + 1L):data_dim) {
        plot(jitter(dat[, i]), jitter(dat[, j]), pch = ".", cex = 2,
             main = paste0("Obs Y", i, " vs Y", j),
             xlab = paste0("Y", i), ylab = paste0("Y", j),
             col = adjustcolor("steelblue", 0.5))
      }
    }

    ## Simulated pairwise scatter plots
    for (i in seq_len(data_dim - 1L)) {
      for (j in (i + 1L):data_dim) {
        plot(jitter(Y_sim[, i]), jitter(Y_sim[, j]), pch = ".", cex = 2,
             main = paste0("Sim Y", i, " vs Y", j),
             xlab = paste0("Y", i), ylab = paste0("Y", j),
             col = adjustcolor("firebrick", 0.5))
      }
    }

    ## Marginal sum barplot
    R_obs <- rowSums(dat)
    R_sim <- rowSums(Y_sim)
    all_vals <- sort(unique(c(R_obs, R_sim)))
    if (length(all_vals) > 30) all_vals <- all_vals[seq(1, length(all_vals), length.out = 30)]
    obs_counts <- tabulate(match(R_obs, all_vals), nbins = length(all_vals))
    sim_counts <- tabulate(match(R_sim, all_vals), nbins = length(all_vals))
    bar_mat <- rbind(Observed = obs_counts / length(R_obs),
                     Simulated = sim_counts / length(R_sim))
    colnames(bar_mat) <- all_vals
    barplot(bar_mat, beside = TRUE, col = c("steelblue", "firebrick"),
            main = paste0("Marginal sum (", data_dim, "D)"),
            xlab = "Sum", ylab = "Proportion",
            legend.text = TRUE, args.legend = list(cex = 0.7, bg = "white"))

    ## Posterior / parameters
    .plot_mdgpd_posterior(obj, est)

    title(paste0(data_dim, "D MDGPD Diagnostics"), outer = TRUE)
  }

  invisible(obj)
}

#' Posterior density / parameter bar chart panel for MDGPD
#' @noRd
.plot_mdgpd_posterior <- function(obj, est) {
  if (!is.null(obj$posterior_samples)) {
    post <- obj$posterior_samples
    par_names <- rownames(post)
    cols_mdgpd <- c("steelblue", "firebrick", "forestgreen",
                     "darkorange", "purple")[seq_len(nrow(post))]
    dens_list <- list()
    max_dens <- 0
    for (i in seq_len(nrow(post))) {
      d <- density(post[i, ], n = 128)
      dens_list[[i]] <- d
      max_dens <- max(max_dens, max(d$y))
    }
    plot(0, 0, type = "n",
         xlim = range(post),
         ylim = c(0, max_dens * 1.1),
         main = "Posterior Marginals (NPE)",
         xlab = "Parameter value", ylab = "Density")
    for (i in seq_len(nrow(post))) {
      lines(dens_list[[i]], col = cols_mdgpd[i], lwd = 2)
    }
    legend("topright", legend = par_names, col = cols_mdgpd,
           lwd = 2, cex = 0.7, bg = "white")
  } else {
    barplot(est, col = "steelblue",
            main = "Parameter Estimates (NBE)",
            ylab = "Estimate", las = 2)
  }
}
