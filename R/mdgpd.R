## Multivariate Discrete Generalized Pareto Distribution (MDGPD)
## Based on Aka, Kratz & Naveau (2025), arXiv:2506.19361
##
## General d-dimensional construction:
##   N_i = G + min(Delta_i, 0)  where  Delta_i = T_i - max_{j != i}(T_j)
##   G ~ Geom(1 - e^{-1}),  T is a d-dimensional equicorrelated Poisson generator.
##
## Non-standard MDGPD (with sigma, xi): apply GPD quantile transform
##   M_i = floor(sigma * (exp(xi * max(N_i, 0)) - 1) / xi)

#' \[Experimental\] Random generation from the multivariate MDGPD
#'
#' Generates random samples from the \code{d}-dimensional Multivariate Discrete
#' Generalized Pareto Distribution (MDGPD) of Aka, Kratz & Naveau (2025).
#' Uses the generator-based construction with an equicorrelated Poisson
#' generator and geometric maximum component. This function is pure R and
#' does not require Julia.
#'
#' @param n integer; number of observations to generate.
#' @param sigma positive numeric; GPD scale parameter (common across dimensions).
#' @param xi non-negative numeric; GPD shape parameter. When \code{xi = 0}
#'   the marginals reduce to scaled geometric.
#' @param lambda positive numeric; Poisson rate for the generator components.
#'   Controls the spread of the dependence structure.
#' @param rho numeric in \[0, 1); equicorrelation of the \code{d}-dimensional
#'   Poisson generator. Higher values give stronger positive dependence:
#'   \code{rho -> 1} gives near-perfect dependence, \code{rho = 0} gives
#'   weakest dependence for the given \code{lambda}.
#' @param d integer >= 2; dimension of the multivariate distribution
#'   (default \code{2L}).
#'
#' @return An \code{n} by \code{d} integer matrix with columns
#'   \code{Y1}, \code{Y2}, \ldots, \code{Yd}.
#'
#' @details
#' The construction follows Aka, Kratz & Naveau (2025), generalised to
#' dimension \code{d}:
#' \enumerate{
#'   \item Generate \code{d}-dimensional equicorrelated Poisson generator
#'     \eqn{T_j = X_j + Z} where \eqn{Z \sim Poisson(\rho\lambda)} (common)
#'     and \eqn{X_j \sim Poisson((1-\rho)\lambda)} (independent), for
#'     \eqn{j = 1, \ldots, d}.
#'   \item For each component \eqn{i}, compute the spectral difference
#'     \eqn{\Delta_i = T_i - \max_{j \neq i} T_j}.
#'   \item Generate \eqn{G \sim Geometric(1 - e^{-1})}, independently.
#'   \item Standard MDGPD: \eqn{N_i = G + \min(\Delta_i, 0)}.
#'   \item Non-standard transform to discrete GPD marginals:
#'     \eqn{M_i = \lfloor \sigma (e^{\xi \max(N_i, 0)} - 1) / \xi \rfloor}.
#' }
#'
#' The parameter \code{rho} controls dependence strength: when \code{rho}
#' is close to 1, all components tend to be equal. The parameter
#' \code{lambda} controls the spread of the spectral differences.
#'
#' For \code{d = 2}, the spectral difference reduces to
#' \eqn{\Delta = T_1 - T_2} and the construction matches the bivariate
#' case in the original paper.
#'
#' This function is \strong{experimental} and its interface may change.
#'
#' @references
#' Aka, S., Kratz, M., and Naveau, P. (2025). Multivariate discrete
#' generalized Pareto distributions: theory, simulation, and applications
#' to dry spells. \emph{arXiv preprint} arXiv:2506.19361.
#'
#' @seealso \code{\link{rzimdgpd}}, \code{\link{rbdegpd}}
#'
#' @examples
#' # Bivariate (default)
#' Y2 <- rmdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5)
#' plot(jitter(Y2[,1]), jitter(Y2[,2]), pch = ".", main = "Bivariate MDGPD")
#'
#' # Trivariate
#' Y3 <- rmdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, d = 3)
#' pairs(Y3 + runif(length(Y3), -0.3, 0.3), pch = ".", main = "Trivariate MDGPD")
#'
#' @export
rmdgpd <- function(n, sigma, xi, lambda, rho, d = 2L) {
  ## Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != floor(n))
    stop("'n' must be a positive integer")
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0)
    stop("'sigma' must be a positive number")
  if (!is.numeric(xi) || length(xi) != 1 || xi < 0)
    stop("'xi' must be a non-negative number")
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
    stop("'lambda' must be a positive number")
  if (!is.numeric(rho) || length(rho) != 1 || rho < 0 || rho >= 1)
    stop("'rho' must be a number in [0, 1)")
  d <- as.integer(d)
  if (d < 2L) stop("'d' must be an integer >= 2")

  ## Step 1: d-dimensional equicorrelated Poisson generator
  ##   T_j = X_j + Z, where Z ~ Pois(rho*lambda), X_j ~ Pois((1-rho)*lambda)
  lambda_z <- rho * lambda
  lambda_indep <- (1 - rho) * lambda
  Z_common <- rpois(n, lambda_z)                 # n-vector (shared)
  X <- matrix(rpois(n * d, lambda_indep), nrow = n, ncol = d)  # n x d
  Tmat <- X + Z_common                           # n x d generator

  ## Step 2: Spectral differences Delta_i = T_i - max_{j != i}(T_j)
  ##   For each row, compute the max of all other columns
  Delta <- matrix(0L, nrow = n, ncol = d)
  for (i in seq_len(d)) {
    max_others <- if (d == 2L) {
      Tmat[, 3L - i]  # the other column
    } else {
      apply(Tmat[, -i, drop = FALSE], 1, max)
    }
    Delta[, i] <- Tmat[, i] - max_others
  }

  ## Step 3: Geometric maximum component
  G <- rgeom(n, prob = 1 - exp(-1))

  ## Step 4: Standard MDGPD (geometric marginals)
  N <- sweep(pmin(Delta, 0L), 1, G, "+")  # N_i = G + min(Delta_i, 0)

  ## Step 5: Non-standard transform to discrete GPD marginals
  N_pos <- pmax(N, 0L)
  if (xi > 1e-8) {
    M <- floor(sigma * (exp(xi * N_pos) - 1) / xi)
  } else {
    M <- floor(sigma * N_pos)
  }

  storage.mode(M) <- "integer"
  colnames(M) <- paste0("Y", seq_len(d))
  M
}


#' \[Experimental\] Random generation from the zero-inflated MDGPD
#'
#' Generates random samples from the Zero-Inflated MDGPD (ZIMDGPD).
#' With probability \code{pi0}, a row is set to all zeros; otherwise,
#' it is drawn from \code{\link{rmdgpd}}.
#'
#' @inheritParams rmdgpd
#' @param pi0 numeric in (0, 1); joint zero-inflation probability.
#'
#' @return An \code{n} by \code{d} integer matrix with columns
#'   \code{Y1}, \code{Y2}, \ldots, \code{Yd}.
#'
#' @details
#' This function is \strong{experimental} and its interface may change.
#'
#' @seealso \code{\link{rmdgpd}}
#'
#' @examples
#' Y <- rzimdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3)
#' mean(rowSums(Y) == 0)  # approximately 0.3 + natural zeros
#'
#' # Trivariate
#' Y3 <- rzimdgpd(1000, sigma = 2, xi = 0.2, lambda = 1, rho = 0.5, pi0 = 0.3, d = 3)
#'
#' @export
rzimdgpd <- function(n, sigma, xi, lambda, rho, pi0, d = 2L) {
  ## Validate pi0
  if (!is.numeric(pi0) || length(pi0) != 1 || pi0 <= 0 || pi0 >= 1)
    stop("'pi0' must be a number in (0, 1)")

  ## Generate MDGPD
  Y <- rmdgpd(n, sigma = sigma, xi = xi, lambda = lambda, rho = rho, d = d)

  ## Apply joint zero-inflation
  zi <- runif(n) < pi0
  Y[zi, ] <- 0L
  Y
}


#' \[Experimental\] Train a neural estimator for MDGPD
#'
#' Trains Neural Posterior Estimator (NPE) and/or Neural Bayesian Estimator
#' (NBE) networks for MDGPD inference. Based on the Aka, Kratz &
#' Naveau (2025) framework. Requires Julia (>= 1.11) with NeuralEstimators.jl
#' and Flux.jl installed, plus the R packages \pkg{JuliaConnectoR} and
#' \pkg{NeuralEstimators}.
#'
#' @param savepath character; directory to save trained model files (.bson).
#'   Created if it does not exist.
#' @param family character; \code{"mdgpd"} (4 parameters) or \code{"zimdgpd"}
#'   (5 parameters, adds pi0).
#' @param data_dim integer >= 2; dimension of the multivariate data
#'   (default \code{2L}). A separate neural network is trained for each
#'   data dimension.
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
#' The training procedure:
#' \enumerate{
#'   \item Samples parameters from uniform priors:
#'     \code{sigma ~ U(0.1, 10)}, \code{xi ~ U(0.01, 0.5)},
#'     \code{lambda ~ U(0.01, 5)}, \code{rho ~ U(0.01, 0.99)},
#'     and for ZIMDGPD: \code{pi0 ~ U(0.01, 0.9)}.
#'   \item Simulates \code{data_dim}-dimensional MDGPD data using
#'     \code{\link{rmdgpd}} or \code{\link{rzimdgpd}}.
#'   \item Applies variance-stabilizing (signed log) and parameter transforms
#'     (\code{log} for sigma/xi/lambda, \code{logit} for rho and pi0).
#'   \item Trains the neural network(s).
#'   \item Saves trained states as .bson files.
#' }
#'
#' Model files are named with a dimension-specific prefix (e.g.,
#' \code{MDGPD_NPE.bson} for \code{data_dim = 2},
#' \code{MDGPD_3D_NPE.bson} for \code{data_dim = 3}).
#'
#' @references
#' Aka, S., Kratz, M., and Naveau, P. (2025). Multivariate discrete
#' generalized Pareto distributions: theory, simulation, and applications
#' to dry spells. \emph{arXiv preprint} arXiv:2506.19361.
#'
#' @seealso \code{\link{rmdgpd}}, \code{\link{rzimdgpd}}, \code{\link{train_begpd}}
#'
#' @examples
#' \dontrun{
#' # 2D (default)
#' paths <- train_mdgpd(savepath = tempdir(), family = "mdgpd", quick = TRUE)
#' # 3D
#' paths3 <- train_mdgpd(savepath = tempdir(), family = "mdgpd",
#'                        data_dim = 3L, quick = TRUE)
#' }
#'
#' @export
train_mdgpd <- function(savepath = "inst/models",
                         family = c("mdgpd", "zimdgpd"),
                         data_dim = 2L,
                         estimator = c("both", "npe", "nbe"),
                         K = NULL, m = NULL,
                         epochs = NULL, stopping_epochs = 10,
                         quick = FALSE,
                         mc.cores = parallel::detectCores() - 1L,
                         seed = 1L, verbose = TRUE) {

  family    <- match.arg(family)
  estimator <- match.arg(estimator)
  is_zi     <- family == "zimdgpd"
  d_params  <- if (is_zi) 5L else 4L
  data_dim  <- as.integer(data_dim)
  if (data_dim < 2L) stop("'data_dim' must be >= 2")

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
    sigma  <- runif(K, 0.1, 10)
    xi     <- runif(K, 0.01, 0.5)
    lambda <- runif(K, 0.01, 5)
    rho    <- runif(K, 0.01, 0.99)
    if (is_zi) {
      pi0 <- runif(K, 0.01, 0.9)
      matrix(c(sigma, xi, lambda, rho, pi0), byrow = TRUE, ncol = K)
    } else {
      matrix(c(sigma, xi, lambda, rho), byrow = TRUE, ncol = K)
    }
  }

  ## --- Data simulator ---
  simulator <- function(Theta, m, mc.cores) {
    parallel::mclapply(seq_len(ncol(Theta)), function(k) {
      m_k <- if (length(m) > 1) sample(m, 1) else m
      theta_k <- Theta[, k]
      if (is_zi) {
        z <- t(rzimdgpd(m_k,
                         sigma = theta_k[1], xi = theta_k[2],
                         lambda = theta_k[3], rho = theta_k[4],
                         pi0 = theta_k[5], d = data_dim))
      } else {
        z <- t(rmdgpd(m_k,
                       sigma = theta_k[1], xi = theta_k[2],
                       lambda = theta_k[3], rho = theta_k[4],
                       d = data_dim))
      }
      z
    }, mc.cores = mc.cores)
  }

  if (verbose) message("Sampling parameters and simulating ", data_dim, "D ",
                       family, " data...")

  ## Generate training and validation sets
  theta_train <- sampler(K)
  theta_val   <- sampler(max(1L, K %/% 10L))
  Z_train <- simulator(theta_train, m, mc.cores)
  Z_val   <- simulator(theta_val, m, mc.cores)

  ## Apply variance-stabilizing transform to data
  signed_log <- .begpd_signed_log
  Z_train <- lapply(Z_train, signed_log)
  Z_val   <- lapply(Z_val, signed_log)

  ## Apply parameter transforms
  theta_train[1:3, ] <- log(theta_train[1:3, ])
  theta_train[4, ]   <- qlogis(theta_train[4, ])
  theta_val[1:3, ]   <- log(theta_val[1:3, ])
  theta_val[4, ]     <- qlogis(theta_val[4, ])
  if (is_zi) {
    theta_train[5, ] <- qlogis(theta_train[5, ])
    theta_val[5, ]   <- qlogis(theta_val[5, ])
  }

  ## File prefix (dimension-specific)
  cfg <- .neuralbayes_config(family, data_dim = data_dim)
  prefix <- cfg$file_prefix

  paths <- list()

  ## --- Train NPE ---
  if (estimator %in% c("both", "npe")) {
    if (verbose) message("Training NPE for ", data_dim, "D ", family, "...")
    npe <- .init_npe_begpd(d = d_params, n = data_dim)
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
    if (verbose) message("Training NBE for ", data_dim, "D ", family, "...")
    nbe <- .init_nbe_begpd(d = d_params, n = data_dim)
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
