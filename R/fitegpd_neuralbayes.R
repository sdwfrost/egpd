## Neural Bayes inference for bivariate BEGPD
## Requires Julia + JuliaConnectoR + NeuralEstimators (Suggests only)

## Package-level environment to track Julia session state
.begpd_env <- new.env(parent = emptyenv())
.begpd_env$julia_initialized <- FALSE
.begpd_env$npe <- NULL
.begpd_env$nbe <- NULL

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

  ## Load Julia packages
  JuliaConnectoR::juliaEval("using NeuralEstimators, Flux")

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
#' @noRd
.init_npe_begpd <- function() {
  n <- 2L  # bivariate
  d <- 6L  # number of parameters

  npe <- JuliaConnectoR::juliaLet("
    q = NormalisingFlow(d, 2d)
    network = initializenetwork(n, 2d)
    PosteriorEstimator(q, network)
  ", d = d, n = n)

  npe
}

#' Initialize NBE estimator object in Julia
#' @noRd
.init_nbe_begpd <- function() {
  n <- 2L  # bivariate
  d <- 6L  # number of parameters
  J <- 5L  # ensemble components

  nbe <- JuliaConnectoR::juliaLet("
    Ensemble([PointEstimator(initializenetwork(n, d)) for j in 1:J])
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

#' Internal: fit bivariate BEGPD via neural Bayes estimation
#'
#' @param x n-by-2 matrix or data.frame of bivariate observations
#' @param model.path path to .bson model file, or NULL for bundled default
#' @param estimator "npe" or "nbe"
#' @param nsamples number of posterior samples (NPE only)
#' @param call the matched call from fitegpd()
#' @noRd
.fitegpd_neuralbayes <- function(x, model.path = NULL,
                                  estimator = c("npe", "nbe"),
                                  nsamples = 1000L,
                                  call = NULL) {
  estimator <- match.arg(estimator)

  ## Validate input
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x))
    stop("For family='begpd', 'x' must be an n-by-2 numeric matrix or data.frame",
         call. = FALSE)
  if (ncol(x) != 2)
    stop("For family='begpd', 'x' must have exactly 2 columns (bivariate data)",
         call. = FALSE)
  if (nrow(x) < 2)
    stop("Need at least 2 observations", call. = FALSE)

  n <- nrow(x)

  ## Resolve model path
  if (is.null(model.path)) {
    bson_name <- paste0(toupper(estimator), ".bson")
    model.path <- system.file("models", bson_name, package = "egpd")
    if (model.path == "")
      stop("Bundled model '", bson_name, "' not found. ",
           "Either provide 'model.path' or run train_begpd() first.",
           call. = FALSE)
  }
  if (!file.exists(model.path))
    stop("Model file not found: ", model.path, call. = FALSE)

  ## Initialize Julia and architecture
  .check_julia_deps()
  .init_julia_begpd()

  ## Create estimator and load trained state
  if (estimator == "npe") {
    est_obj <- .init_npe_begpd()
  } else {
    est_obj <- .init_nbe_begpd()
  }
  NeuralEstimators::loadstate(est_obj, model.path)

  ## Prepare data: transpose to 2 x n, apply signed_log, wrap in list
  Z <- t(x)  # 2 x n
  Z <- .begpd_signed_log(Z)
  Z_list <- list(Z)

  ## Inference
  if (estimator == "npe") {
    ## NPE: sample from approximate posterior
    ## sampleposterior(estimator, Z, N) — N is positional, returns list of d x N matrices
    samples_raw <- NeuralEstimators::sampleposterior(est_obj, Z_list,
                                                      as.integer(nsamples))

    ## Result is a list (one element per dataset); extract the matrix
    if (is.list(samples_raw)) {
      post_mat <- samples_raw[[1]]
    } else {
      post_mat <- samples_raw
    }

    ## Exp-transform back from log-parameter space
    post_mat <- exp(as.matrix(post_mat))
    if (nrow(post_mat) != 6) post_mat <- t(post_mat)
    rownames(post_mat) <- .begpd_parnames

    ## Posterior summaries
    estimate <- apply(post_mat, 1, median)
    sd_est   <- apply(post_mat, 1, sd)
    names(estimate) <- .begpd_parnames
    names(sd_est)   <- .begpd_parnames
    posterior_samples <- post_mat

  } else {
    ## NBE: point estimates
    ## estimate() returns a matrix directly
    est_raw <- NeuralEstimators::estimate(est_obj, Z_list)

    ## Result may be a list or matrix; extract and exp-transform
    if (is.list(est_raw)) est_raw <- est_raw[[1]]
    estimate <- exp(as.numeric(est_raw))
    names(estimate) <- .begpd_parnames
    sd_est <- rep(NA_real_, 6)
    names(sd_est) <- .begpd_parnames
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
    npar              = 6L,
    data              = x,
    type              = NA_integer_,
    family            = "begpd",
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
    model.path        = model.path,
    posterior_samples  = posterior_samples,
    nsamples          = if (estimator == "npe") nsamples else NULL
  ), class = "fitegpd")
}


#' 4-panel diagnostic plot for bivariate BEGPD fits
#' @noRd
.plot_begpd <- function(obj, ...) {
  est <- obj$estimate
  dat <- obj$data

  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(op))

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
