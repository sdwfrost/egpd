## Bivariate Discrete EGPD (BDEGPD) and Zero-Inflated BDEGPD (BZIDEGPD)
## Experimental extensions: floor(BEGPD) for discrete bivariate data

#' \[Experimental\] Random generation from the bivariate discrete EGPD
#'
#' Generates random samples from the Bivariate Discrete Extended Generalized
#' Pareto Distribution (BDEGPD) by applying \code{floor()} to continuous
#' bivariate BEGPD samples. This function is pure R and does not require Julia.
#'
#' @param n integer; number of observations to generate.
#' @param kappa positive numeric; EGPD shape parameter (power transform).
#' @param sigma positive numeric; GPD scale parameter.
#' @param xi positive numeric; GPD shape parameter.
#' @param thL positive numeric; lower tail dependence parameter (beta shape).
#' @param thU positive numeric; upper tail dependence parameter (beta shape).
#' @param thw numeric in (0, 0.5); weight mixing parameter.
#'
#' @return An \code{n} by 2 integer matrix with columns \code{Y1} and \code{Y2}.
#'
#' @details
#' The BDEGPD constructs bivariate discrete observations by generating
#' continuous BEGPD samples via \code{\link{rbegpd}} and applying
#' \code{floor()} to obtain non-negative integers.
#'
#' This function is \strong{experimental} and its interface may change.
#'
#' @seealso \code{\link{rbegpd}}, \code{\link{rbzidegpd}}
#'
#' @examples
#' Y <- rbdegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
#' table(Y[,1])
#'
#' @export
rbdegpd <- function(n, kappa, sigma, xi, thL, thU, thw) {
  ## Generate continuous BEGPD and discretize
  Y <- rbegpd(n, kappa = kappa, sigma = sigma, xi = xi,
               thL = thL, thU = thU, thw = thw)
  Y <- floor(Y)
  storage.mode(Y) <- "integer"
  Y
}


#' \[Experimental\] Random generation from the zero-inflated bivariate discrete EGPD
#'
#' Generates random samples from the Zero-Inflated Bivariate Discrete Extended
#' Generalized Pareto Distribution (BZIDEGPD). With probability \code{pi0},
#' a row is set to \code{(0, 0)}; otherwise, it is drawn from \code{\link{rbdegpd}}.
#'
#' @inheritParams rbdegpd
#' @param pi0 numeric in (0, 1); joint zero-inflation probability.
#'
#' @return An \code{n} by 2 integer matrix with columns \code{Y1} and \code{Y2}.
#'
#' @details
#' This function is \strong{experimental} and its interface may change.
#'
#' @seealso \code{\link{rbdegpd}}, \code{\link{rbegpd}}
#'
#' @examples
#' Y <- rbzidegpd(1000, kappa = 2, sigma = 1, xi = 0.1,
#'                thL = 5, thU = 5, thw = 0.2, pi0 = 0.3)
#' mean(Y[,1] == 0 & Y[,2] == 0)  # approximately 0.3
#'
#' @export
rbzidegpd <- function(n, kappa, sigma, xi, thL, thU, thw, pi0) {
  ## Validate pi0
  if (!is.numeric(pi0) || length(pi0) != 1 || pi0 <= 0 || pi0 >= 1)
    stop("'pi0' must be a number in (0, 1)")

  ## Generate BDEGPD
  Y <- rbdegpd(n, kappa = kappa, sigma = sigma, xi = xi,
                thL = thL, thU = thU, thw = thw)

  ## Apply joint zero-inflation
  zi <- runif(n) < pi0
  Y[zi, ] <- 0L
  Y
}


#' \[Experimental\] Train a neural estimator for bivariate discrete EGPD
#'
#' Trains Neural Posterior Estimator (NPE) and/or Neural Bayesian Estimator
#' (NBE) networks for bivariate discrete EGPD inference (BDEGPD or BZIDEGPD).
#' Requires Julia (>= 1.11) with NeuralEstimators.jl and Flux.jl installed,
#' plus the R packages \pkg{JuliaConnectoR} and \pkg{NeuralEstimators}.
#'
#' @param savepath character; directory to save trained model files (.bson).
#'   Created if it does not exist.
#' @param family character; \code{"bdegpd"} (6 parameters) or \code{"bzidegpd"}
#'   (7 parameters, adds pi0).
#' @param estimator character; which estimator(s) to train: \code{"both"}
#'   (default), \code{"npe"}, or \code{"nbe"}.
#' @param K integer; number of training parameter sets to sample from the prior.
#' @param m integer vector; range of sample sizes for simulated datasets.
#' @param epochs integer; maximum number of training epochs.
#' @param stopping_epochs integer; early stopping patience.
#' @param quick logical; if \code{TRUE}, uses reduced training settings.
#' @param mc.cores integer; number of parallel cores for data simulation.
#' @param seed integer or \code{NULL}; random seed for reproducibility.
#' @param verbose logical; if \code{TRUE}, print progress messages.
#'
#' @return A named list of file paths to saved .bson model files (invisible).
#'
#' @details
#' This function is \strong{experimental} and its interface may change.
#'
#' The training procedure mirrors \code{\link{train_begpd}} but uses
#' \code{\link{rbdegpd}} or \code{\link{rbzidegpd}} for data simulation.
#' For BZIDEGPD, the seventh parameter (pi0) uses a logit transform.
#'
#' @seealso \code{\link{train_begpd}}, \code{\link{rbdegpd}}, \code{\link{rbzidegpd}}
#'
#' @examples
#' \dontrun{
#' paths <- train_bdegpd(savepath = tempdir(), family = "bdegpd", quick = TRUE)
#' }
#'
#' @export
train_bdegpd <- function(savepath = "inst/models",
                          family = c("bdegpd", "bzidegpd"),
                          estimator = c("both", "npe", "nbe"),
                          K = NULL, m = NULL,
                          epochs = NULL, stopping_epochs = 10,
                          quick = FALSE,
                          mc.cores = parallel::detectCores() - 1L,
                          seed = 1L, verbose = TRUE) {

  family    <- match.arg(family)
  estimator <- match.arg(estimator)
  is_zi     <- family == "bzidegpd"
  d         <- if (is_zi) 7L else 6L

  ## Windows guard
  if (.Platform$OS.type == "windows") mc.cores <- 1L

  ## Check Julia dependencies
  .check_julia_deps()

  ## Quick mode defaults
  if (quick) {
    K      <- K      %||% 20000L
    m      <- m      %||% 500L:1000L
    epochs <- epochs %||% 10L
  } else {
    K      <- K      %||% 100000L
    m      <- m      %||% 1000L:4000L
    epochs <- epochs %||% 100L
  }

  if (!is.null(seed)) set.seed(seed)

  dir.create(savepath, recursive = TRUE, showWarnings = FALSE)

  ## Initialize Julia session and architecture
  .init_julia_begpd()

  ## --- Prior sampler ---
  sampler <- function(K) {
    kappa <- runif(K, 0.1, 10)
    sigma <- runif(K, 0.1, 3)
    xi    <- runif(K, 0.01, 0.5)
    thL   <- runif(K, 0.1, 20)
    thU   <- runif(K, 0.1, 20)
    thw   <- runif(K, 0.01, 0.49)
    if (is_zi) {
      pi0 <- runif(K, 0.01, 0.9)
      matrix(c(kappa, sigma, xi, thL, thU, thw, pi0), byrow = TRUE, ncol = K)
    } else {
      matrix(c(kappa, sigma, xi, thL, thU, thw), byrow = TRUE, ncol = K)
    }
  }

  ## --- Data simulator ---
  simulator <- function(Theta, m, mc.cores) {
    parallel::mclapply(seq_len(ncol(Theta)), function(k) {
      m_k <- if (length(m) > 1) sample(m, 1) else m
      theta_k <- Theta[, k]
      if (is_zi) {
        z <- t(rbzidegpd(m_k,
                          kappa = theta_k[1], sigma = theta_k[2],
                          xi = theta_k[3], thL = theta_k[4],
                          thU = theta_k[5], thw = theta_k[6],
                          pi0 = theta_k[7]))
      } else {
        z <- t(rbdegpd(m_k,
                        kappa = theta_k[1], sigma = theta_k[2],
                        xi = theta_k[3], thL = theta_k[4],
                        thU = theta_k[5], thw = theta_k[6]))
      }
      z
    }, mc.cores = mc.cores)
  }

  if (verbose) message("Sampling parameters and simulating ", family, " data...")

  ## Generate training and validation sets
  theta_train <- sampler(K)
  theta_val   <- sampler(max(1L, K %/% 10L))
  Z_train <- simulator(theta_train, m, mc.cores)
  Z_val   <- simulator(theta_val, m, mc.cores)

  ## Apply variance-stabilizing transform to data
  signed_log <- .begpd_signed_log
  Z_train <- lapply(Z_train, signed_log)
  Z_val   <- lapply(Z_val, signed_log)

  ## Apply parameter transforms: log for first 6, logit for pi0
  if (is_zi) {
    theta_train[1:6, ] <- log(theta_train[1:6, ])
    theta_train[7, ]   <- qlogis(theta_train[7, ])
    theta_val[1:6, ]   <- log(theta_val[1:6, ])
    theta_val[7, ]     <- qlogis(theta_val[7, ])
  } else {
    theta_train <- log(theta_train)
    theta_val   <- log(theta_val)
  }

  ## File prefix for this family
  cfg <- .neuralbayes_config(family)
  prefix <- cfg$file_prefix

  paths <- list()

  ## --- Train NPE ---
  if (estimator %in% c("both", "npe")) {
    if (verbose) message("Training NPE for ", family, "...")
    npe <- .init_npe_begpd(d = d)
    npe <- NeuralEstimators::train(
      npe,
      theta_train = theta_train,
      theta_val   = theta_val,
      Z_train     = Z_train,
      Z_val       = Z_val,
      epochs      = epochs,
      stopping_epochs = stopping_epochs,
      savepath    = file.path(savepath, paste0(prefix, "NPE"))
    )
    npe_path <- file.path(savepath, paste0(prefix, "NPE.bson"))
    NeuralEstimators::savestate(npe, npe_path)
    paths$npe <- npe_path
    if (verbose) message("NPE saved to ", npe_path)
  }

  ## --- Train NBE ---
  if (estimator %in% c("both", "nbe")) {
    if (verbose) message("Training NBE for ", family, "...")
    nbe <- .init_nbe_begpd(d = d)
    nbe <- NeuralEstimators::train(
      nbe,
      theta_train = theta_train,
      theta_val   = theta_val,
      Z_train     = Z_train,
      Z_val       = Z_val,
      epochs      = epochs,
      stopping_epochs = stopping_epochs,
      savepath    = file.path(savepath, paste0(prefix, "NBE"))
    )
    nbe_path <- file.path(savepath, paste0(prefix, "NBE.bson"))
    NeuralEstimators::savestate(nbe, nbe_path)
    paths$nbe <- nbe_path
    if (verbose) message("NBE saved to ", nbe_path)
  }

  if (verbose) message("Training complete.")
  invisible(paths)
}
