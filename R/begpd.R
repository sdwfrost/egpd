## Multivariate EGPD (BEGPD) — random generation and training
## Bivariate only (d = 2)

#' Random generation from the bivariate BEGPD
#'
#' Generates random samples from the bivariate Multivariate Extended
#' Generalized Pareto Distribution (BEGPD). This function is pure R and
#' does not require Julia.
#'
#' @param n integer; number of observations to generate.
#' @param kappa positive numeric; EGPD shape parameter (power transform).
#' @param sigma positive numeric; GPD scale parameter.
#' @param xi non-negative numeric; GPD shape parameter.
#' @param thL positive numeric; lower tail dependence parameter (beta shape).
#' @param thU positive numeric; upper tail dependence parameter (beta shape).
#' @param thw numeric in (0, 0.5); weight mixing parameter.
#'
#' @return An \code{n} by 2 numeric matrix with columns \code{Y1} and \code{Y2}.
#'
#' @details
#' The BEGPD constructs bivariate observations by:
#' \enumerate{
#'   \item Generating a radial component \eqn{R} from a power-transformed GPD.
#'   \item Generating lower and upper dependence components from symmetric
#'         Beta distributions with parameters \code{thL} and \code{thU}.
#'   \item Mixing the components using a weight function based on \code{thw}.
#' }
#'
#' @references
#' Alotaibi, N., Sainsbury-Dale, M., Naveau, P., Gaetan, C., and Huser, R.
#' (2025). Joint modeling of low and high extremes using a multivariate
#' extended generalized Pareto distribution. \emph{arXiv preprint}
#' arXiv:2509.05982.
#'
#' @examples
#' Y <- rbegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
#' plot(Y[,1], Y[,2], pch = ".", main = "Bivariate BEGPD sample")
#'
#' @export
rbegpd <- function(n, kappa, sigma, xi, thL, thU, thw) {
  ## Input validation
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n < 1 || n != floor(n))
    stop("'n' must be a positive integer")
  if (!is.numeric(kappa) || length(kappa) != 1 || !is.finite(kappa) || kappa <= 0)
    stop("'kappa' must be a positive number")
  if (!is.numeric(sigma) || length(sigma) != 1 || !is.finite(sigma) || sigma <= 0)
    stop("'sigma' must be a positive number")
  if (!is.numeric(xi) || length(xi) != 1 || !is.finite(xi) || xi < 0)
    stop("'xi' must be a non-negative number")
  if (!is.numeric(thL) || length(thL) != 1 || !is.finite(thL) || thL <= 0)
    stop("'thL' must be a positive number")
  if (!is.numeric(thU) || length(thU) != 1 || !is.finite(thU) || thU <= 0)
    stop("'thU' must be a positive number")
  if (!is.numeric(thw) || length(thw) != 1 || !is.finite(thw) || thw <= 0 || thw >= 0.5)
    stop("'thw' must be a number in (0, 0.5)")

  ## Radial component: power-transformed GPD (cancellation-free at large kappa)
  R_unif <- runif(n)
  log_u_kappa <- log(R_unif) / kappa
  one_minus_u <- -expm1(log_u_kappa)
  log_one_minus_u <- ifelse(-log_u_kappa > log(2),
                            log1p(-exp(log_u_kappa)),
                            log(one_minus_u))
  R <- .qgpd_from_logsurv(log_one_minus_u, sigma, xi)

  ## By beta symmetry the smaller share is Q(min(u, 1-u)). Only form
  ## the *larger* share by subtraction: 1 - rbeta(...) could round the
  ## smaller one to zero in the second margin even at shape 0.11.
  angle <- function(shape) {
    u <- runif(n)
    small <- qbeta(pmin(u, 1 - u), shape, shape)
    cbind(Y1 = ifelse(u <= 0.5, small, 1 - small),
          Y2 = ifelse(u <= 0.5, 1 - small, small))
  }
  lo <- angle(thL)
  hi <- angle(thU)
  t <- (R_unif - thw) / (1 - 2 * thw)
  R * (pbeta(t, 3, 3, lower.tail = FALSE) * lo + pbeta(t, 3, 3) * hi)
}


#' Train a neural estimator for bivariate BEGPD
#'
#' Trains Neural Posterior Estimator (NPE) and/or Neural Bayesian Estimator
#' (NBE) networks for bivariate BEGPD inference. Requires Julia (>= 1.11)
#' with NeuralEstimators.jl (>= 0.3 for training), Flux.jl, BSON.jl and
#' Optimisers.jl installed (in Julia,
#' \code{using Pkg; Pkg.add(["NeuralEstimators", "Flux", "BSON", "Optimisers"])}),
#' plus the R packages \pkg{JuliaConnectoR} and \pkg{NeuralEstimators}.
#'
#' @param savepath character; directory to save trained model files (.bson).
#'   Created if it does not exist. The default \code{NULL} uses
#'   \code{file.path(tools::R_user_dir("egpd", "data"), "models")} and reports
#'   it in a message.
#' @param estimator character; which estimator(s) to train: \code{"both"}
#'   (default), \code{"npe"}, or \code{"nbe"}.
#' @param K integer; number of training parameter sets to sample from the prior.
#' @param m integer vector; range of sample sizes for simulated datasets
#'   (e.g., \code{1000:4000}). Recorded in the model sidecar; inference warns
#'   when the number of observations falls outside \code{range(m)}.
#' @param epochs integer; maximum number of training epochs.
#' @param stopping_epochs integer; early stopping patience (epochs without
#'   improvement).
#' @param quick logical; if \code{TRUE}, uses reduced training settings for
#'   quick testing (K=20000, m=500:1000, epochs=10).
#' @param mc.cores integer; number of parallel cores for data simulation.
#'   \code{NULL} (default) uses \code{parallel::detectCores() - 1}. Forced to
#'   1 on Windows and capped at 2 when \code{_R_CHECK_LIMIT_CORES_} is set
#'   (as by \code{R CMD check --as-cran}).
#' @param seed integer or \code{NULL}; seeds R, per-dataset parallel
#'   simulations, and Julia's \code{Random.seed!}. Backend/version-specific
#'   nondeterminism may remain.
#' @param verbose logical; if \code{TRUE}, print progress messages.
#'
#' @return A named list of file paths to saved .bson model files (invisible).
#'
#' @details
#' Newly trained models require their generated \code{.bson.metadata.rds}
#' compatibility sidecars. No bundled .bson models are shipped; train
#' compatible weights and pass \code{model.path}. Compatibility does not
#' establish calibration. The sidecar also records the training sample-size
#' range \code{range(m)} and the Julia, NeuralEstimators.jl, Flux.jl,
#' NeuralEstimators (R) and egpd versions; \code{distfit(method =
#' "neuralbayes")} warns when \code{nrow(x)} is outside the training range,
#' when the data are outside the simulated support, or when the loaded
#' NeuralEstimators.jl or Flux.jl version differs from the training one.
#'
#' The training procedure:
#' \enumerate{
#'   \item Samples parameters from uniform priors.
#'   \item Simulates bivariate BEGPD data using \code{\link{rbegpd}}.
#'   \item Applies variance-stabilizing (signed log) and log parameter transforms.
#'   \item Trains the neural network(s) using \code{NeuralEstimators::train()}.
#'   \item Saves trained states as .bson files.
#' }
#'
#' The NPE is a likelihood-to-evidence ratio network, not the normalising
#' flow of Alotaibi et al. (2025). \code{distfit()} turns it into posterior
#' draws by adaptive importance sampling of prior times learned ratio (see
#' \code{\link{distfit}}), reporting the effective sample size.
#'
#' @section Supported Julia/R versions:
#' The egpd Julia helpers (\code{inst/julia/begpd_architecture.jl}) were run
#' directly with Julia 1.12.7, NeuralEstimators.jl 0.1.0 and Flux 0.14.25
#' (network construction, NBE estimates, log-ratio evaluation, BSON
#' save/load); the log-ratio helper was also run with the
#' \code{RatioEstimator} code of NeuralEstimators.jl 0.1.3 and 0.1.4, which
#' store the network in a \code{network} field. They follow the
#' NeuralEstimators.jl 0.2.1 to 0.3.x APIs (\code{RatioEstimator} constructors
#' and \code{logratio(estimator, Z; grid)}), checked against the upstream
#' sources but not executed. NeuralEstimators.jl 0.2.0 is not supported: its
#' \code{logratio()} has no \code{grid} keyword, so NPE inference fails.
#' Training goes through \code{NeuralEstimators::train()}, and the R package
#' NeuralEstimators 0.2.2 (the GitHub version) needs NeuralEstimators.jl 0.3 or
#' later: its \code{train()} passes \code{shuffle}, which the training methods
#' accept only from 0.3.0 (\code{device} from 0.2.1;
#' \code{freeze_summary_network} is not accepted by 0.1.x). It also runs
#' \code{using Optimisers: Adam}, so Optimisers.jl must be a direct dependency
#' of the Julia environment. Its \code{loadstate()} and \code{savestate()}
#' reference \code{LuxEstimator}, which 0.1.x does not define, so egpd saves
#' and loads weights with its own Flux/BSON helpers (the same file format).
#' Training with these functions therefore requires NeuralEstimators.jl 0.3 or
#' later with R NeuralEstimators 0.2.2.
#'
#' @references
#' Sainsbury-Dale, M., Zammit-Mangion, A., and Huser, R. (2024). Likelihood-free
#' parameter estimation with neural Bayes estimators. \emph{The American
#' Statistician}, 78(1), 1--14.
#'
#' @examples
#' \dontrun{
#' paths <- train_begpd(savepath = tempdir(), quick = TRUE)
#' }
#'
#' @export
train_begpd <- function(savepath = NULL,
                        estimator = c("both", "npe", "nbe"),
                        K = NULL, m = NULL,
                        epochs = NULL, stopping_epochs = 10,
                        quick = FALSE,
                        mc.cores = NULL,
                        seed = 1L, verbose = TRUE) {

  estimator <- match.arg(estimator)

  mc.cores <- .sanitize_mc_cores(mc.cores)

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

  savepath <- .neural_savepath(savepath)

  ## Initialize Julia session and architecture
  .init_julia_begpd()
  .neural_seed(seed)
  if (!is.null(seed)) set.seed(seed)

  ## --- Prior sampler ---
  sampler <- function(K) {
    kappa <- runif(K, 0.1, 10)
    sigma <- runif(K, 0.1, 3)
    xi    <- runif(K, 0.01, 0.5)
    thL   <- runif(K, 0.1, 20)
    thU   <- runif(K, 0.1, 20)
    thw   <- runif(K, 0.01, 0.49)
    matrix(c(kappa, sigma, xi, thL, thU, thw), byrow = TRUE, ncol = K)
  }

  ## --- Data simulator ---
  simulator <- function(Theta, m, mc.cores) {
    .neural_simulate(seq_len(ncol(Theta)), function(k) {
      m_k <- if (length(m) > 1) sample(m, 1) else m
      theta_k <- Theta[, k]
      z <- t(rbegpd(m_k,
                     kappa = theta_k[1], sigma = theta_k[2],
                     xi = theta_k[3], thL = theta_k[4],
                     thU = theta_k[5], thw = theta_k[6]))
      z
    }, mc.cores = mc.cores)
  }

  if (verbose) message("Sampling parameters and simulating data...")

  ## Generate training and validation sets
  theta_train <- sampler(K)
  theta_val   <- sampler(max(1L, K %/% 10L))
  Z_train <- simulator(theta_train, m, mc.cores)
  Z_val   <- simulator(theta_val, m, mc.cores)

  ## Apply variance-stabilizing transform to data
  signed_log <- .begpd_signed_log
  Z_train <- lapply(Z_train, signed_log)
  Z_val   <- lapply(Z_val, signed_log)

  ## Apply log transform to parameters
  theta_train <- log(theta_train)
  theta_val   <- log(theta_val)

  paths <- list()
  training <- .neural_training_record(m, K, epochs, stopping_epochs, seed)

  ## --- Train NPE ---
  if (estimator %in% c("both", "npe")) {
    if (verbose) message("Training NPE...")
    npe <- .init_npe_begpd()
    npe <- .neural_train(
      npe,
      theta_train = theta_train,
      theta_val   = theta_val,
      Z_train     = Z_train,
      Z_val       = Z_val,
      epochs      = epochs,
      stopping_epochs = stopping_epochs,
      savepath    = file.path(savepath, "NPE")
    )
    npe_path <- file.path(savepath, "NPE.bson")
    .neural_savestate(npe, npe_path)
    .write_neural_metadata(npe_path, "begpd", "npe", training = training)
    paths$npe <- npe_path
    if (verbose) message("NPE saved to ", npe_path)
  }

  ## --- Train NBE ---
  if (estimator %in% c("both", "nbe")) {
    if (verbose) message("Training NBE...")
    nbe <- .init_nbe_begpd()
    nbe <- .neural_train(
      nbe,
      theta_train = theta_train,
      theta_val   = theta_val,
      Z_train     = Z_train,
      Z_val       = Z_val,
      epochs      = epochs,
      stopping_epochs = stopping_epochs,
      savepath    = file.path(savepath, "NBE")
    )
    nbe_path <- file.path(savepath, "NBE.bson")
    .neural_savestate(nbe, nbe_path)
    .write_neural_metadata(nbe_path, "begpd", "nbe", training = training)
    paths$nbe <- nbe_path
    if (verbose) message("NBE saved to ", nbe_path)
  }

  if (verbose) message("Training complete.")
  invisible(paths)
}

## Null-coalescing operator (if not already available)
`%||%` <- function(x, y) if (is.null(x)) y else x

## Default directory for trained weights: persistent per-user data directory.
.neural_savepath <- function(savepath) {
  if (is.null(savepath)) {
    savepath <- file.path(tools::R_user_dir("egpd", which = "data"), "models")
    message("Saving trained neural models to ", savepath,
            " (set 'savepath' to choose another directory).")
  }
  if (!is.character(savepath) || length(savepath) != 1L || is.na(savepath) ||
      !nzchar(savepath))
    stop("'savepath' must be a single directory path", call. = FALSE)
  dir.create(savepath, recursive = TRUE, showWarnings = FALSE)
  savepath
}

## Training provenance stored in the model sidecar (not part of the contract).
.neural_training_record <- function(m, K, epochs, stopping_epochs, seed) {
  list(m = range(as.integer(m)), K = as.integer(K), epochs = as.integer(epochs),
       stopping_epochs = as.integer(stopping_epochs),
       seed = if (is.null(seed)) NULL else as.integer(seed))
}
