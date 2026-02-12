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
#' @param xi positive numeric; GPD shape parameter.
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

  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != floor(n))
    stop("'n' must be a positive integer")
  if (!is.numeric(kappa) || length(kappa) != 1 || kappa <= 0)
    stop("'kappa' must be a positive number")
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0)
    stop("'sigma' must be a positive number")
  if (!is.numeric(xi) || length(xi) != 1 || xi <= 0)
    stop("'xi' must be a positive number")
  if (!is.numeric(thL) || length(thL) != 1 || thL <= 0)
    stop("'thL' must be a positive number")
  if (!is.numeric(thU) || length(thU) != 1 || thU <= 0)
    stop("'thU' must be a positive number")
  if (!is.numeric(thw) || length(thw) != 1 || thw <= 0 || thw >= 0.5)
    stop("'thw' must be a number in (0, 0.5)")

  ## Radial component: power-transformed GPD

  R_unif <- runif(n)
  R <- sigma * ((1 - R_unif^(1 / kappa))^(-xi) - 1) / xi

  ## Lower tail dependence (symmetric beta)
  V1 <- rbeta(n, thL, thL)
  L1 <- 1 / V1
  L2 <- 1 / (1 - V1)

  ## Upper tail dependence (symmetric beta)
  U1 <- rbeta(n, thU, thU)
  U2 <- 1 - U1

  ## Mixing weights based on radial quantile
  wR <- pbeta((R_unif - thw) / (1 - 2 * thw), 3, 3)

  ## Bivariate output
  Y1 <- R * ((1 - wR) * L1^(-1) + wR * U1)
  Y2 <- R * ((1 - wR) * L2^(-1) + wR * U2)
  Y <- cbind(Y1 = Y1, Y2 = Y2)

  Y
}


#' Train a neural estimator for bivariate BEGPD
#'
#' Trains Neural Posterior Estimator (NPE) and/or Neural Bayesian Estimator
#' (NBE) networks for bivariate BEGPD inference. Requires Julia (>= 1.11)
#' with NeuralEstimators.jl and Flux.jl installed, plus the R packages
#' \pkg{JuliaConnectoR} and \pkg{NeuralEstimators}.
#'
#' @param savepath character; directory to save trained model files (.bson).
#'   Created if it does not exist.
#' @param estimator character; which estimator(s) to train: \code{"both"}
#'   (default), \code{"npe"}, or \code{"nbe"}.
#' @param K integer; number of training parameter sets to sample from the prior.
#' @param m integer vector; range of sample sizes for simulated datasets
#'   (e.g., \code{1000:4000}).
#' @param epochs integer; maximum number of training epochs.
#' @param stopping_epochs integer; early stopping patience (epochs without
#'   improvement).
#' @param quick logical; if \code{TRUE}, uses reduced training settings for
#'   quick testing (K=20000, m=500:1000, epochs=10).
#' @param mc.cores integer; number of parallel cores for data simulation.
#'   Forced to 1 on Windows.
#' @param seed integer or \code{NULL}; random seed for reproducibility.
#' @param verbose logical; if \code{TRUE}, print progress messages.
#'
#' @return A named list of file paths to saved .bson model files (invisible).
#'
#' @details
#' The training procedure:
#' \enumerate{
#'   \item Samples parameters from uniform priors.
#'   \item Simulates bivariate BEGPD data using \code{\link{rbegpd}}.
#'   \item Applies variance-stabilizing (signed log) and log parameter transforms.
#'   \item Trains the neural network(s) using \code{NeuralEstimators::train()}.
#'   \item Saves trained states as .bson files.
#' }
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
train_begpd <- function(savepath = "inst/models",
                        estimator = c("both", "npe", "nbe"),
                        K = NULL, m = NULL,
                        epochs = NULL, stopping_epochs = 10,
                        quick = FALSE,
                        mc.cores = parallel::detectCores() - 1L,
                        seed = 1L, verbose = TRUE) {

  estimator <- match.arg(estimator)

  ## Windows guard: mclapply doesn't parallelize on Windows

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
    matrix(c(kappa, sigma, xi, thL, thU, thw), byrow = TRUE, ncol = K)
  }

  ## --- Data simulator ---
  simulator <- function(Theta, m, mc.cores) {
    parallel::mclapply(seq_len(ncol(Theta)), function(k) {
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

  ## --- Train NPE ---
  if (estimator %in% c("both", "npe")) {
    if (verbose) message("Training NPE...")
    npe <- .init_npe_begpd()
    npe <- NeuralEstimators::train(
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
    NeuralEstimators::savestate(npe, npe_path)
    paths$npe <- npe_path
    if (verbose) message("NPE saved to ", npe_path)
  }

  ## --- Train NBE ---
  if (estimator %in% c("both", "nbe")) {
    if (verbose) message("Training NBE...")
    nbe <- .init_nbe_begpd()
    nbe <- NeuralEstimators::train(
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
    NeuralEstimators::savestate(nbe, nbe_path)
    paths$nbe <- nbe_path
    if (verbose) message("NBE saved to ", nbe_path)
  }

  if (verbose) message("Training complete.")
  invisible(paths)
}

## Null-coalescing operator (if not already available)
`%||%` <- function(x, y) if (is.null(x)) y else x
