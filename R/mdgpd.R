## Multivariate Discrete Generalized Pareto Distribution (MDGPD)
## Based on Aka, Kratz & Naveau (2025), arXiv:2506.19361
##
## Standard: N = T - max(T) + ceiling(E), E ~ Exp(1).
## Nonstandard: M = ceiling(h^{-1}(E + T - max(T))).

#' \[Experimental\] Random generation from the multivariate MDGPD
#'
#' Generates random samples from the \code{d}-dimensional Multivariate Discrete
#' Generalized Pareto Distribution (MDGPD) of Aka, Kratz & Naveau (2025).
#' Uses the generator-based construction with an equicorrelated Poisson
#' generator and geometric maximum component. This function is pure R and
#' does not require Julia.
#'
#' @param n integer; number of observations to generate.
#' @param sigma positive finite numeric scalar or length-\code{d} vector;
#'   GPD scale parameter(s).
#' @param xi finite numeric scalar or length-\code{d} vector; GPD shape
#'   parameter(s) (may be negative, zero, or positive).
#' @param lambda,rho obsolete generator arguments; supplying either is an
#'   error. Migrate explicitly to \code{rate = (1 - rho) * lambda}.
#' @param rate non-negative finite numeric; rate of independent Poisson generator
#'   components (\code{rate = 0} gives the exact comonotonic limit). This is
#'   the identifiable combination of the former \code{lambda} and \code{rho};
#'   smaller rates give more nearly equal components.
#' @param d integer >= 2; dimension of the multivariate distribution
#'   (default \code{2L}).
#'
#' @return An \code{n} by \code{d} integer-valued matrix with columns
#'   \code{Y1}, \code{Y2}, \ldots, \code{Yd}. Storage is double when counts
#'   exceed the R integer range. Components may be negative (threshold excesses);
#'   every row has a strictly positive maximum.
#'
#' @details
#' The construction follows Aka, Kratz & Naveau (2025), generalised to
#' dimension \code{d}:
#' \enumerate{
#'   \item Generate \code{d}-dimensional equicorrelated Poisson generator
#'     \eqn{T_j \sim Poisson(rate)} independently, for \eqn{j=1,\ldots,d}
#'     (with \eqn{T_j \equiv 0} when \code{rate = 0}).
#'   \item Set \eqn{S_i = T_i - \max_j T_j}.
#'   \item Independently generate \eqn{E \sim Exp(1)}. Then
#'     \eqn{G=\lceil E\rceil} has positive geometric support \eqn{1,2,\ldots}.
#'   \item Standard MDGPD: \eqn{N_i=S_i+G} (Theorem 2.1, Eq. (7)).
#'   \item For general scale \eqn{\sigma_i > 0} and shape
#'     \eqn{\xi_i \in \mathbb{R}}, return
#'     \eqn{M_i=\lceil \sigma_i\,\mathrm{expm1}(\xi_i(E+S_i))/\xi_i\rceil},
#'     with limit \eqn{\lceil\sigma_i(E+S_i)\rceil} at zero shape.
#'     When \eqn{\xi_i < 0}, \eqn{M_i \le \lceil -\sigma_i/\xi_i\rceil} has a
#'     finite upper endpoint. When \eqn{\xi_i>0}, the unrounded lower endpoint
#'     \eqn{-\sigma_i/\xi_i} is open: even an integer endpoint is never
#'     attained by the ceiling. If resolving that endpoint would require
#'     distinguishing adjacent integers beyond \eqn{2^{53}}, simulation stops.
#'     This realizes Definition 2.1, not a transform of
#'     the rounded \eqn{N}.
#' }
#'
#' A common Poisson component cancels in \eqn{T-\max(T)}, so its rate
#' cannot be inferred. Only \code{rate} is fitted. Neural estimation
#' (\code{\link{train_mdgpd}}) covers the homogeneous model with scalar
#' \code{sigma} and \code{xi} only; vector values can be simulated but not
#' estimated.
#' Marginal invariance is \emph{conditional}: for standard parameters,
#' \eqn{P(N_i>k\mid N_i>0)=e^{-k}} for nonnegative integers \eqn{k}.
#' In general \eqn{P(M_i>k\mid M_i>0)=(1+\xi_i k/\sigma_i)_+^{-1/\xi_i}}
#' (exponential limit at zero shape). Unconditional marginals are not GPDs.
#'
#' This replaces the former censored, zero-based geometric simulator.
#' Both the support and distribution have changed; old fits and bundled
#' MDGPD neural assets are incompatible and must be retrained, not relabeled.
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
#' Y2 <- rmdgpd(1000, sigma = 2, xi = 0.2, rate = 0.5)
#' plot(jitter(Y2[,1]), jitter(Y2[,2]), pch = ".", main = "Bivariate MDGPD")
#'
#' # Negative shape (bounded upper tail, e.g. dry-spell lengths)
#' Y_neg <- rmdgpd(1000, sigma = c(3, 4), xi = c(-0.18, -0.25), rate = 0.5)
#'
#' # Trivariate
#' Y3 <- rmdgpd(1000, sigma = 2, xi = 0.2, rate = 0.5, d = 3)
#' pairs(Y3 + runif(length(Y3), -0.3, 0.3), pch = ".", main = "Trivariate MDGPD")
#'
#' @export
rmdgpd <- function(n, sigma, xi, lambda = NULL, rho = NULL, d = 2L,
                   rate = NULL) {
  if (!is.null(lambda) || !is.null(rho))
    stop("Legacy 'lambda'/'rho' arguments are not supported. Use rate = (1-rho)*lambda explicitly; the MDGPD simulator has changed and old fits require retraining.",
         call. = FALSE)
  ## Input validation
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n < 1 || n != floor(n))
    stop("'n' must be a positive integer")
  if (!is.numeric(d) || length(d) != 1 || !is.finite(d) ||
      d < 2 || d != floor(d) || d > .Machine$integer.max)
    stop("'d' must be an integer >= 2")
  d <- as.integer(d)
  if (!is.numeric(sigma) || !(length(sigma) %in% c(1L, d)) ||
      any(!is.finite(sigma)) || any(sigma <= 0))
    stop("'sigma' must be a positive finite number or length-'d' vector")
  if (!is.numeric(xi) || !(length(xi) %in% c(1L, d)) || any(!is.finite(xi)))
    stop("'xi' must be a finite number or length-'d' vector")
  if (!is.numeric(rate) || length(rate) != 1 || !is.finite(rate) || rate < 0)
    stop("'rate' must be a non-negative finite number")

  Tmat <- if (rate == 0) {
    matrix(0L, nrow = n, ncol = d)
  } else {
    matrix(rpois(n * d, rate), nrow = n, ncol = d)
  }
  S <- Tmat - apply(Tmat, 1L, max)
  W <- S + rexp(n)

  component <- function(w, s, shape) {
    ## On the negative half-line reflect both the argument and the shape
    ## to use the same range-safe inverse as the univariate GP.
    x <- if (shape == 0) s * w else
      sign(w) * .qgpd_from_logsurv(-abs(w), s, sign(w) * shape)
    m <- ceiling(x)
    m[w > 0 & m == 0] <- 1  # a strictly positive product can underflow
    if (shape < 0) {
      ## A bounded positive draw may round to an integer just below the
      ## true endpoint. Compare to the cell edge before taking the ceiling.
      at <- which(w > 0 & x == m & m == round(-s/shape) & m >= 1 & m < 2^53)
      if (length(at)) {
        h <- -.gpd_logsurv_inside(m[at], rep(s, length(at)), rep(shape, length(at)))
        up <- at[w[at] > h]
        m[up] <- m[up] + 1
      }
    }
    if (shape > 0) {
      ## In the first lattice cell expm1(z) can round to -1. Compare
      ## s*exp(z) to the remainder of s modulo shape instead of subtracting
      ## two nearly equal numbers. A zero remainder is *strictly* exceeded.
      ratio <- s / shape
      at <- which(w < 0 & is.finite(m) & m <= -floor(ratio))
      if (length(at)) {
        if (ratio > 2^53)
          stop("Cannot resolve individual counts at the open lower endpoint beyond ",
               "double integer precision; adjust sigma, xi, or rate.", call. = FALSE)
        ## Split a wide quotient before taking the remainder (%/% and %%
        ## can lose precision at ratios above 1/double.eps).
        unit <- if (ratio >= 2^52) 2^26 else 1
        hi <- if (unit == 1) 0 else s %/% (shape * unit)
        lo <- if (unit == 1) s else s %% (shape * unit)
        k <- hi * unit + lo %/% shape
        r <- lo %% shape
        m[at] <- -k + (r == 0 | shape * w[at] > log(r) - log(s))
      }
    }
    m
  }
  if (length(sigma) == 1L && length(xi) == 1L) {
    M <- component(W, sigma, xi)
  } else {
    sigma <- rep_len(as.numeric(sigma), d)
    xi <- rep_len(as.numeric(xi), d)
    M <- matrix(0, nrow = n, ncol = d)
    for (j in seq_len(d)) {
      M[, j] <- component(W[, j], sigma[j], xi[j])
    }
  }

  if (any(!is.finite(M)))
    stop("Simulated counts exceed finite numeric range; adjust sigma, xi, or rate.")
  if (all(abs(M) <= .Machine$integer.max)) storage.mode(M) <- "integer"
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
#' @return An \code{n} by \code{d} integer-valued matrix with columns
#'   \code{Y1}, \code{Y2}, \ldots, \code{Yd}.
#'
#' @details
#' This function is \strong{experimental} and its interface may change.
#'
#' @seealso \code{\link{rmdgpd}}
#'
#' @examples
#' Y <- rzimdgpd(1000, sigma = 2, xi = 0.2, rate = 0.5, pi0 = 0.3)
#' mean(apply(Y, 1, function(y) all(y == 0)))  # approximately 0.3
#'
#' # Trivariate
#' Y3 <- rzimdgpd(1000, sigma = 2, xi = 0.2, rate = 0.5, pi0 = 0.3, d = 3)
#'
#' @export
rzimdgpd <- function(n, sigma, xi, lambda = NULL, rho = NULL, pi0, d = 2L,
                     rate = NULL) {
  ## Validate pi0
  if (!is.numeric(pi0) || length(pi0) != 1 || !is.finite(pi0) || pi0 <= 0 || pi0 >= 1)
    stop("'pi0' must be a number in (0, 1)")

  ## Generate MDGPD
  Y <- rmdgpd(n, sigma = sigma, xi = xi, lambda = lambda, rho = rho,
              rate = rate, d = d)

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
#' (>= 0.3 for training), Flux.jl, BSON.jl and Optimisers.jl installed (in Julia,
#' \code{using Pkg; Pkg.add(["NeuralEstimators", "Flux", "BSON", "Optimisers"])}),
#' plus the R packages \pkg{JuliaConnectoR} and \pkg{NeuralEstimators} (see
#' "Supported Julia/R versions" in \code{\link{train_begpd}}).
#'
#' @param savepath character; directory to save trained model files (.bson).
#'   Created if it does not exist. The default \code{NULL} uses
#'   \code{file.path(tools::R_user_dir("egpd", "data"), "models")} and reports
#'   it in a message.
#' @param family character; \code{"mdgpd"} (3 parameters: sigma, xi, rate)
#'   or \code{"zimdgpd"} (4 parameters, adds pi0).
#' @param data_dim integer >= 2; dimension of the multivariate data
#'   (default \code{2L}). A separate neural network is trained for each
#'   data dimension.
#' @param estimator character; which estimator(s) to train: \code{"both"}
#'   (default), \code{"npe"}, or \code{"nbe"}.
#' @param K integer; number of training parameter sets to sample from the prior.
#' @param m integer vector; range of sample sizes for simulated datasets.
#'   Recorded in the model sidecar; inference warns when the number of
#'   observations falls outside \code{range(m)}.
#' @param epochs integer; maximum number of training epochs.
#' @param stopping_epochs integer; early stopping patience.
#' @param quick logical; if \code{TRUE}, uses reduced training settings.
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
#' This function is \strong{experimental} and its interface may change.
#'
#' \strong{Homogeneous MDGPD only.} The networks are trained for the
#' exchangeable model with one scale \eqn{\sigma_1 = \cdots = \sigma_d}, one
#' shape \eqn{\xi_1 = \cdots = \xi_d} and a single Poisson generator rate for
#' all components (Aka, Kratz & Naveau 2025, Table 1 types (a)-(b)). The
#' marginal-specific \eqn{(\sigma_j, \xi_j)}, asymmetric generator rates
#' \eqn{\lambda_1 \ne \lambda_2} and shifted generators (Table 1 type (c))
#' fitted in the paper's dry-spell application (Tables 2-3) cannot be
#' represented: \code{\link{rmdgpd}} simulates vector \code{sigma}/\code{xi},
#' but \code{distfit(family = "mdgpd", method = "neuralbayes")} estimates
#' only the common \code{(sigma, xi, rate)}.
#'
#' The training procedure:
#' \enumerate{
#'   \item Samples parameters from uniform priors:
#'     \code{sigma ~ U(0.1, 10)}, \code{xi ~ U(-0.4, 0.5)},
#'     \code{rate ~ U(0.01, 5)},
#'     and for ZIMDGPD: \code{pi0 ~ U(0.01, 0.9)}.
#'   \item Simulates \code{data_dim}-dimensional MDGPD data using
#'     \code{\link{rmdgpd}} or \code{\link{rzimdgpd}}.
#'   \item Applies variance-stabilizing (signed log) and parameter transforms
#'     (\code{log} for sigma/rate, \code{identity} for xi, \code{logit} for pi0).
#'   \item Trains the neural network(s).
#'   \item Saves trained states as .bson files with required
#'     \code{.metadata.rds} compatibility sidecars. No bundled .bson models
#'     are shipped; train compatible weights and pass \code{model.path}.
#' }
#'
#' Model files are named with a dimension-specific prefix (e.g.,
#' \code{MDGPD_NPE.bson} for \code{data_dim = 2},
#' \code{MDGPD_3D_NPE.bson} for \code{data_dim = 3}). The sidecar records
#' the training sample-size range and backend versions; see
#' \code{\link{train_begpd}} for the supported Julia/R version pairings.
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
train_mdgpd <- function(savepath = NULL,
                         family = c("mdgpd", "zimdgpd"),
                         data_dim = 2L,
                         estimator = c("both", "npe", "nbe"),
                         K = NULL, m = NULL,
                         epochs = NULL, stopping_epochs = 10,
                         quick = FALSE,
                         mc.cores = NULL,
                         seed = 1L, verbose = TRUE) {

  family    <- match.arg(family)
  estimator <- match.arg(estimator)
  is_zi     <- family == "zimdgpd"
  d_params  <- if (is_zi) 4L else 3L
  if (length(data_dim) != 1L || !is.finite(data_dim) ||
      data_dim < 2 || data_dim != floor(data_dim))
    stop("'data_dim' must be an integer >= 2")
  data_dim <- as.integer(data_dim)

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
    sigma <- runif(K, 0.1, 10)
    xi    <- runif(K, -0.4, 0.5)
    rate  <- runif(K, 0.01, 5)
    if (is_zi) {
      pi0 <- runif(K, 0.01, 0.9)
      matrix(c(sigma, xi, rate, pi0), byrow = TRUE, ncol = K)
    } else {
      matrix(c(sigma, xi, rate), byrow = TRUE, ncol = K)
    }
  }

  ## --- Data simulator ---
  simulator <- function(Theta, m, mc.cores) {
    .neural_simulate(seq_len(ncol(Theta)), function(k) {
      m_k <- if (length(m) > 1) sample(m, 1) else m
      theta_k <- Theta[, k]
      if (is_zi) {
        z <- t(rzimdgpd(m_k,
                         sigma = theta_k[1], xi = theta_k[2],
                         rate = theta_k[3], pi0 = theta_k[4], d = data_dim))
      } else {
        z <- t(rmdgpd(m_k,
                       sigma = theta_k[1], xi = theta_k[2],
                       rate = theta_k[3],
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

  ## Apply parameter transforms: log for sigma/rate, identity for xi, logit for pi0
  theta_train[c(1L, 3L), ] <- log(theta_train[c(1L, 3L), ])
  theta_val[c(1L, 3L), ]   <- log(theta_val[c(1L, 3L), ])
  if (is_zi) {
    theta_train[4L, ] <- qlogis(theta_train[4L, ])
    theta_val[4L, ]   <- qlogis(theta_val[4L, ])
  }

  ## File prefix (dimension-specific)
  cfg <- .neuralbayes_config(family, data_dim = data_dim)
  prefix <- cfg$file_prefix

  paths <- list()
  training <- .neural_training_record(m, K, epochs, stopping_epochs, seed)

  ## --- Train NPE ---
  if (estimator %in% c("both", "npe")) {
    if (verbose) message("Training NPE for ", data_dim, "D ", family, "...")
    npe <- .init_npe_begpd(d = d_params, n = data_dim)
    npe <- .neural_train(
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
    .neural_savestate(npe, npe_path)
    .write_neural_metadata(npe_path, family, "npe", data_dim,
                           training = training)
    paths$npe <- npe_path
    if (verbose) message("NPE saved to ", npe_path)
  }

  ## --- Train NBE ---
  if (estimator %in% c("both", "nbe")) {
    if (verbose) message("Training NBE for ", data_dim, "D ", family, "...")
    nbe <- .init_nbe_begpd(d = d_params, n = data_dim)
    nbe <- .neural_train(
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
    .neural_savestate(nbe, nbe_path)
    .write_neural_metadata(nbe_path, family, "nbe", data_dim,
                           training = training)
    paths$nbe <- nbe_path
    if (verbose) message("NBE saved to ", nbe_path)
  }

  if (verbose) message("Training complete.")
  invisible(paths)
}
