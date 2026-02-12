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
#' @param d integer; number of parameters (default 6L for BEGPD/BDEGPD, 7L for BZIDEGPD).
#' @param n integer; data dimension (default 2L for bivariate).
#' @noRd
.init_npe_begpd <- function(d = 6L, n = 2L) {
  npe <- JuliaConnectoR::juliaLet("
    q = NormalisingFlow(d, 2d)
    network = initializenetwork(n, 2d)
    PosteriorEstimator(q, network)
  ", d = d, n = n)

  npe
}

#' Initialize NBE estimator object in Julia
#' @param d integer; number of parameters (default 6L for BEGPD/BDEGPD, 7L for BZIDEGPD).
#' @param n integer; data dimension (default 2L for bivariate).
#' @noRd
.init_nbe_begpd <- function(d = 6L, n = 2L) {
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
      d           = 4L,
      data_dim    = as.integer(data_dim),
      file_prefix = .mdgpd_prefix("MDGPD_", data_dim),
      parnames    = c("sigma", "xi", "lambda", "rho"),
      backtransform = function(mat) {
        ## log for sigma (1), xi (2), lambda (3); logit for rho (4)
        mat[1:3, ] <- exp(mat[1:3, , drop = FALSE])
        mat[4, ]   <- plogis(mat[4, ])
        mat
      }
    ),
    "zimdgpd" = list(
      d           = 5L,
      data_dim    = as.integer(data_dim),
      file_prefix = .mdgpd_prefix("ZIMDGPD_", data_dim),
      parnames    = c("sigma", "xi", "lambda", "rho", "pi0"),
      backtransform = function(mat) {
        ## log for sigma (1), xi (2), lambda (3); logit for rho (4), pi0 (5)
        mat[1:3, ] <- exp(mat[1:3, , drop = FALSE])
        mat[4:5, ] <- plogis(mat[4:5, , drop = FALSE])
        mat
      }
    ),
    stop("Unknown neural-Bayes family: ", family, call. = FALSE)
  )
}

#' Internal: fit multivariate families via neural Bayes estimation
#'
#' @param x n-by-d matrix or data.frame of observations (d >= 2).
#' @param family character; "begpd", "bdegpd", "bzidegpd", "mdgpd", or "zimdgpd".
#' @param model.path path to .bson model file, or NULL for bundled default
#' @param estimator "npe" or "nbe"
#' @param nsamples number of posterior samples (NPE only)
#' @param call the matched call from fitegpd()
#' @noRd
.fitegpd_neuralbayes <- function(x, family = "begpd",
                                  model.path = NULL,
                                  estimator = c("npe", "nbe"),
                                  nsamples = 1000L,
                                  call = NULL) {
  estimator <- match.arg(estimator)

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
      stop("Bundled model '", bson_name, "' not found. ",
           "Either provide 'model.path' or run ", train_fn,
           " first.", call. = FALSE)
    }
  }
  if (!file.exists(model.path))
    stop("Model file not found: ", model.path, call. = FALSE)

  ## Initialize Julia and architecture
  .check_julia_deps()
  .init_julia_begpd()

  ## Create estimator and load trained state
  if (estimator == "npe") {
    est_obj <- .init_npe_begpd(d = cfg$d, n = data_dim)
  } else {
    est_obj <- .init_nbe_begpd(d = cfg$d, n = data_dim)
  }
  NeuralEstimators::loadstate(est_obj, model.path)

  ## Prepare data: transpose to data_dim x n, apply signed_log, wrap in list
  Z <- t(x)  # data_dim x n
  Z <- .begpd_signed_log(Z)
  Z_list <- list(Z)

  ## Inference
  if (estimator == "npe") {
    ## NPE: sample from approximate posterior
    samples_raw <- NeuralEstimators::sampleposterior(est_obj, Z_list,
                                                      as.integer(nsamples))

    ## Result is a list (one element per dataset); extract the matrix
    if (is.list(samples_raw)) {
      post_mat <- samples_raw[[1]]
    } else {
      post_mat <- samples_raw
    }

    ## Back-transform from unconstrained parameter space
    post_mat <- as.matrix(post_mat)
    if (nrow(post_mat) != cfg$d) post_mat <- t(post_mat)
    post_mat <- cfg$backtransform(post_mat)
    rownames(post_mat) <- cfg$parnames

    ## Posterior summaries
    estimate <- apply(post_mat, 1, median)
    sd_est   <- apply(post_mat, 1, sd)
    names(estimate) <- cfg$parnames
    names(sd_est)   <- cfg$parnames
    posterior_samples <- post_mat

  } else {
    ## NBE: point estimates
    est_raw <- NeuralEstimators::estimate(est_obj, Z_list)

    ## Result may be a list or matrix; extract and back-transform
    if (is.list(est_raw)) est_raw <- est_raw[[1]]
    raw_vec <- as.numeric(est_raw)
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


#' 4-panel diagnostic plot for bivariate discrete EGPD fits (BDEGPD/BZIDEGPD)
#' @noRd
.plot_bdegpd <- function(obj, ...) {
  est <- obj$estimate
  dat <- obj$data
  family <- obj$family
  is_zi <- family == "bzidegpd"

  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(op))

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

  ## Simulate from fitted model
  n_sim <- min(nrow(dat), 5000L)
  if (is_zi) {
    Y_sim <- rzimdgpd(n_sim,
                       sigma = est["sigma"], xi = est["xi"],
                       lambda = est["lambda"], rho = est["rho"],
                       pi0 = est["pi0"], d = data_dim)
  } else {
    Y_sim <- rmdgpd(n_sim,
                     sigma = est["sigma"], xi = est["xi"],
                     lambda = est["lambda"], rho = est["rho"],
                     d = data_dim)
  }

  if (data_dim == 2L) {
    ## Bivariate: 4-panel layout
    op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    on.exit(par(op))

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
    ## d >= 3: pairs plot layout
    ## Row 1: observed pairs (lower triangle) + simulated pairs (upper)
    ## Plus one panel for posterior / parameters
    n_panels <- data_dim + 1L
    op <- par(mfrow = c(2, n_panels), mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))
    on.exit(par(op))

    ## Row 1: Observed pairwise scatter plots
    for (i in seq_len(data_dim)) {
      for (j in seq_len(data_dim)) {
        if (i >= j) next
      }
    }

    ## Use pairs-style layout: one row observed, one row simulated
    layout_n <- choose(data_dim, 2)
    n_cols <- min(layout_n, 3L)
    n_rows_data <- ceiling(layout_n / n_cols)
    op <- par(mfrow = c(n_rows_data * 2L + 1L, n_cols),
              mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))
    on.exit(par(op))

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
