## Univariate distribution fitting for EGPD families
## MLE and Bernstein polynomial fitting

#' Fit a distribution to data by maximum likelihood
#'
#' Maximum likelihood, Bernstein polynomial, or neural Bayes fitting of the
#' distribution families supported by the package: (discrete, zero-inflated)
#' EGPD, the composite-Pareto model, the mixed-Poisson count families
#' (Poisson-inverse Gaussian and generalised PIG), the Bell distribution, and
#' the bivariate/multivariate EGPD families. Formerly named \code{fitegpd};
#' that name is retained as a deprecated alias.
#'
#' @param x numeric vector of observations (univariate families), or an
#'   n-by-d numeric matrix/data.frame (multivariate families; d=2 for BEGPD,
#'   d >= 2 for MDGPD). Univariate data must lie in the support: positive
#'   values for \code{"egpd"} (the continuous EGPD has no mass at 0, so exact
#'   zeros are an error; use \code{"ziegpd"} or a discrete family for them),
#'   non-negative values for \code{"ziegpd"} and \code{"cpegpd"}, and
#'   non-negative integers for the discrete families; otherwise \code{distfit()}
#'   stops.
#' @param type integer 1-6 specifying the G-transformation type (univariate only).
#' @param family character: "egpd", "degpd", "ziegpd", "zidegpd", "cpegpd",
#'   "cpdegpd", "pig" (Poisson-inverse Gaussian; \code{type} ignored), "zipig"
#'   (zero-inflated PIG), "gpig" (generalised Poisson-inverse Gaussian of
#'   Zhu & Joe 2009, native \code{(a, b, c)} parameterisation; \code{type}
#'   ignored), "zigpig" (zero-inflated GPIG), "bell" (Bell distribution of
#'   Castellares et al. 2018, native \code{theta} parameterisation; \code{type}
#'   ignored), "zibell" (zero-inflated Bell), "gwaring" (generalised Waring;
#'   \code{type} ignored), "begpd" (bivariate EGPD),
#'   "bdegpd" (\[Experimental\] bivariate
#'   discrete EGPD), "bzidegpd" (\[Experimental\] zero-inflated bivariate
#'   discrete EGPD), "mdgpd" (\[Experimental\] multivariate MDGPD via
#'   Aka-Kratz-Naveau construction, d >= 2), or "zimdgpd" (\[Experimental\]
#'   zero-inflated multivariate MDGPD, d >= 2).
#' @param method character: "mle", "bernstein", or "neuralbayes".
#'   \code{"neuralbayes"} is required for bivariate families and requires
#'   Julia dependencies.
#' @param start named list of starting values, or NULL for automatic
#'   (not used for method="neuralbayes"). Every value must carry the name of its
#'   parameter and be a single number inside that parameter's space; otherwise
#'   \code{distfit()} stops.
#' @param fix.arg named list of fixed parameters (not used for
#'   method="neuralbayes"), checked like \code{start}, except that a probability
#'   may also be fixed at 0 or 1 where its distribution allows it: GPIG
#'   \code{a} in (0, 1] and \code{c} in [0, 1], the zero-inflation \code{pi} of
#'   ZIPIG, ZIGPIG and ZIBELL in [0, 1), and the \code{pi} of ZIEGPD and ZIDEGPD
#'   and the type-6 \code{prob} in [0, 1]. In type 2 \code{kappa} can
#'   be fixed at its uniform-carrier limit of zero.
#' @param optim.method optimization method passed to \code{\link{optim}}.
#'   With finite optimizer bounds, use L-BFGS-B unless Brent is requested.
#'   For \code{method = "mle"}, if a feasible step materially improves an
#'   L-BFGS-B solution, the fit retries with smaller steps within the same
#'   bounds rather than accepting the optimizer's apparent convergence.
#' @param hessian logical: compute standard errors via Hessian? MLE and Bernstein
#'   fits warn and leave standard errors unavailable if the Hessian is not
#'   positive definite or the estimate is on an optimizer bound.
#'   Fixed parameters are excluded from estimation and
#'   parameter counts. In exact subfamilies only identifiable coordinates
#'   are fitted: for GPIG with \code{a = 1} fixed and both \code{b} and
#'   \code{c} free, the fitted parameter is \code{mu = b*c} (starts in
#'   either the original coordinates or \code{mu} are accepted). With
#'   equal fixed type-6 powers, the mixture weight is not estimated; a
#'   component of fixed weight zero has no estimated power. For a type-1
#'   ZIDEGPD whose fixed endpoint implies support \eqn{\{0,1\}}, free
#'   \code{kappa} and \code{pi} identify only \code{p1 = P(Y = 1)}.
#'   The fit reports its interior Bernoulli MLE, not separate estimates
#'   or standard errors for the two native parameters. An all-zero or
#'   all-one sample has no interior two-cell estimate without an additional
#'   bound. For \code{cpdegpd}
#'   with two-cell jumps, free intensity and free jump probability, the
#'   aggregate is Poisson and only its mean \code{mu} is fitted. Other two-cell
#'   carriers with multiple unidentified parameters require further pins.
#'   More generally, a bounded count law with \eqn{K} support cells can
#'   identify at most \eqn{K-1} parameter combinations. A fixed endpoint
#'   with too many free parameters requires additional pins. If the
#'   endpoint is estimated this check uses the fitted support, not the
#'   set of observed counts. For \code{cpdegpd} it applies to the
#'   \emph{jump} support, not the aggregate counts: with two jump cells
#'   the aggregate law depends only on \eqn{\lambda P(J=1)}.
#'   On a fitted-support ridge the fit warns and withholds Wald inference;
#'   this is a ridge of the
#'   fitted likelihood, not a reduction of the entire model. See
#'   \code{reduction} in the result. MLE and Bernstein fits also warn when \code{optim} does not
#'   report convergence (e.g. its iteration limit was reached), and stop when no
#'   parameter value visited gives the data a positive likelihood (observations
#'   outside the support of every candidate distribution have log-likelihood
#'   \code{-Inf}).
#'   An all-zero Bell sample has no interior maximum; fitting stops instead
#'   of reporting a Wald SE at an arbitrary approximation to the boundary.
#'   When zero inflation is free it is profiled conditional on the base
#'   parameters. If the observed zero fraction is at or below the base
#'   probability of zero, the conditional MLE is exactly zero unless a
#'   finite logit bound excludes it. A natural or user-imposed boundary
#'   estimate has no ordinary Wald inference; an interior estimate uses
#'   the joint information including the other free parameters.
#'   For a free type-2 \code{kappa} or type-6 mixing weight \code{prob},
#'   fitting also checks the included natural endpoints. A competitive
#'   endpoint is refitted with the remaining parameters free and checked
#'   for an improving direction into the model. At a selected endpoint
#'   the estimate is exact and joint Wald inference is withheld. A power
#'   of a mixture component with fitted weight zero is not identified.
#'   In a continuous type-3 fit, \code{kappa = 0} is excluded. Without a
#'   positive lower bound on \code{kappa}, fitting stops if the likelihood
#'   at the fitted base parameters increases toward that limit.
#'   With two count or jump cells and only the type-6 weight free, the
#'   mean is affine in the weight. The fit uses the exact Bernoulli or
#'   thinned Poisson optimum and its natural-scale information.
#'   For a continuous EGPD with free shape, an interior regular solution
#'   is a local fit, not a guarantee of a global maximum over all shapes.
#'   With fixed or fitted \eqn{\xi < -1}, fitting stops if the continuous
#'   endpoint can approach the largest positive observation from above:
#'   the likelihood is then unbounded. At \eqn{\xi=-1}, if \eqn{u g(u)}
#'   increases for carrier density \eqn{g}, the likelihood increases as
#'   free scale approaches the sample maximum from above. This is checked
#'   for types 1--4 and 6, and for type 5 with \eqn{\delta\le1} or
#'   \eqn{\kappa\ge2}. Its boundary
#'   limit is not an interior MLE; fitting stops instead of inverting a
#'   Hessian across the support boundary. Neither check applies to floored counts.
#'   For floored counts with \eqn{\xi\le-1/2}, the likelihood need not be
#'   twice differentiable when the endpoint crosses the upper edge of the
#'   highest observed cell. If the numerical Hessian spans that edge, the fit
#'   retains its estimates, warns and omits ordinary Wald inference.
#'   For type-1 floored counts without zero inflation, with fixed
#'   \eqn{\xi=-1} and observations in and below the highest cell, scale
#'   is fitted at its conditional maximum for each power, subject to its
#'   bounds. This includes a maximum at the edge of the highest cell.
#'   For a compound fit the same check applies to changes of jump support at
#'   grid edges that contribute to an observed total.
#' @param bernstein.m integer: Bernstein polynomial degree (method="bernstein" only)
#' @param cpegpd.h numeric: gauge resolution \eqn{h} of the cpegpd family. The
#'   cpegpd law is that of \eqn{S_h = \sum_i h \lceil X_i/h \rceil}, a compound
#'   Poisson sum of EGPD amounts each rounded up to a multiple of \eqn{h} (see
#'   \code{\link{dcpegpd}}); for continuous (off-lattice) data use an \eqn{h}
#'   much smaller than typical amounts, otherwise \code{distfit()} warns. The
#'   cpegpd log-likelihood uses the probability of an exact zero for zero
#'   observations and the bin probability divided by \code{cpegpd.h} (a density
#'   approximation) for positive ones, i.e. a likelihood with respect to a point
#'   mass at zero plus Lebesgue measure: comparable across \code{cpegpd.h}, but
#'   not AIC-comparable with continuous families that place a density, not a
#'   mass, at zero.
#' @param model.path character: path to a trained \code{.bson} model file
#'   with its \code{.bson.metadata.rds} sidecar (\code{method="neuralbayes"}).
#' @param estimator character: \code{"npe"} for Neural Posterior Estimation or
#'   \code{"nbe"} for Neural Bayesian Estimation (method="neuralbayes" only).
#' @param nsamples integer: number of posterior samples for NPE
#'   (method="neuralbayes" with estimator="npe" only).
#' @param npe.control optional named list of controls for the adaptive
#'   importance sampler that draws NPE posterior samples
#'   (method="neuralbayes", estimator="npe"): \code{K} (proposal draws per
#'   stage, default 5000), \code{K0} (prior draws in the first stage, default
#'   20000), \code{ess_target} (effective sample size at which sampling stops,
#'   default \code{nsamples} clamped to [200, 1000]), \code{ess_floor}
#'   (minimum ESS kept by the tempering schedule, default
#'   \code{max(200, 25 d)}), \code{max_stages} (default 30), \code{df}
#'   (degrees of freedom of the multivariate-t proposal, default 5),
#'   \code{inflate} (proposal covariance inflation, default 1.5) and
#'   \code{defensive} (weight of the prior in the proposal mixture, default
#'   0.1). \code{K}, \code{K0} and \code{max_stages} must be whole numbers
#'   (\code{max_stages} may be 0), \code{ess_target} and \code{ess_floor} at
#'   least 1, \code{df} and \code{inflate} positive, \code{defensive} in
#'   \eqn{[0, 1)}, and \code{K} at least twice \code{ess_floor}; other values
#'   are an error.
#' @param ... additional arguments passed to \code{\link{optim}}. For the
#'   closed-form Bernoulli and thinned Poisson reductions, scalar \code{lower}
#'   and \code{upper} bound \eqn{\mathrm{logit}(p_1)} or \eqn{\log(\mu)},
#'   respectively. A constrained maximum on a bound has no ordinary Wald SE.
#'
#' @details For \code{method = "neuralbayes"} with \code{estimator = "npe"},
#'   posterior draws are obtained by adaptive importance sampling of the prior
#'   times the learned likelihood-to-evidence ratio: the first stage samples the
#'   prior, later stages fit heavy-tailed proposals to the reweighted particles,
#'   and the final draws are a systematic resample of the importance weights.
#'   The effective sample size (ESS), number of stages and number of distinct
#'   draws are reported, with a warning when the ESS is below 200 or there are
#'   few distinct draws. Neural MDGPD estimation covers the homogeneous model
#'   only (common \eqn{\sigma}, \eqn{\xi} and generator rate).
#'
#'   For \code{method = "bernstein"}, the log-likelihood is evaluated in log
#'   space (no clipping of the GPD probabilities, so tiny observations are exact)
#'   with its exact gradient. With the default Nelder-Mead optimiser, whose
#'   simplex can collapse on the flat ridges of this likelihood, the solution is
#'   refined by alternating BFGS (with the exact gradient) and fresh Nelder-Mead
#'   simplices. When a cycle stalls, the best direction toward a single
#'   mixture component is checked on the weight scale: a near-zero weight
#'   can hide an improvement on the logit scale. \code{convergence} is 0
#'   when neither refinement improves the fit materially and 1 if this
#'   does not happen within 20 cycles. A single free
#'   parameter is fitted by BFGS, as for \code{method = "mle"}. These
#'   unbounded refinements are not used when box constraints are supplied.
#'
#' @return An object of class \code{"distfit"} with components:
#' \describe{
#'   \item{estimate}{named vector of parameter estimates}
#'   \item{sd}{named vector of standard errors (NA if Hessian not computed)}
#'   \item{vcov}{variance-covariance matrix on natural scale}
#'   \item{loglik}{maximized log-likelihood}
#'   \item{aic}{Akaike information criterion}
#'   \item{bic}{Bayesian information criterion}
#'   \item{n}{number of observations}
#'   \item{npar}{number of estimated parameters}
#'   \item{data}{the input data vector}
#'   \item{type}{the G-transformation type}
#'   \item{family}{the distribution family}
#'   \item{method}{the fitting method}
#'   \item{fix.arg}{list of fixed arguments}
#'   \item{reduction}{for a reduced MLE, the reason and representative constants
#'     used only for evaluating the distribution, not parameter estimates}
#'   \item{support.ridge}{whether the fitted count or jump support is too
#'     small to identify all free parameters (ordinary MLE fits)}
#'   \item{support.change}{whether observed-information inference was withheld
#'     because the fitted count support changes across the Hessian stencil}
#'   \item{convergence}{convergence code from optim (0 = success)}
#'   \item{optim}{full optim output}
#'   \item{call}{the matched call}
#'   \item{bernstein.m}{Bernstein degree (NULL for method="mle")}
#'   \item{bernstein.weights}{Bernstein weights (NULL for method="mle")}
#' }
#'
#' For the neural families (\code{method = "neuralbayes"}), additional fields:
#' \describe{
#'   \item{estimator_type}{"npe" or "nbe"}
#'   \item{model.path}{path to the .bson model used}
#'   \item{posterior_samples}{parameter-by-\code{nsamples} matrix of posterior
#'     draws (NPE) or NULL (NBE)}
#'   \item{nsamples}{number of posterior samples (NPE) or NULL (NBE)}
#'   \item{importance_sampling}{NPE only: list with \code{ess},
#'     \code{ess_target}, \code{K}, \code{stages}, \code{unique},
#'     \code{nsamples} and \code{path} of the adaptive importance sampler; NPE
#'     fits have class \code{c("distfit_npe", "distfit")}}
#' }
#'
#' @examples
#' \dontrun{
#' # Univariate fitting
#' x <- regpd(500, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
#' fit <- distfit(x, type = 1)
#' summary(fit)
#' plot(fit)
#'
#' # Bivariate BEGPD (requires Julia and a trained model: no models are bundled;
#' # see train_begpd())
#' Y <- rbegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
#' fit_biv <- distfit(Y, family = "begpd", method = "neuralbayes",
#'                    model.path = "path/to/begpd_model.bson")
#' summary(fit_biv)
#' plot(fit_biv)
#' }
#'
#' @name distfit
#' @aliases fitegpd
#' @export
distfit <- function(x, type = 1,
                    family = c("egpd", "degpd", "ziegpd", "zidegpd",
                               "cpegpd", "cpdegpd", "pig", "zipig",
                               "gpig", "zigpig", "bell", "zibell", "gwaring", "begpd",
                               "bdegpd", "bzidegpd",
                               "mdgpd", "zimdgpd"),
                    method = c("mle", "bernstein", "neuralbayes"),
                    start = NULL, fix.arg = NULL,
                    optim.method = "Nelder-Mead", hessian = TRUE,
                    bernstein.m = 8, cpegpd.h = 0.2,
                    model.path = NULL, estimator = c("npe", "nbe"),
                    nsamples = 1000L, npe.control = list(), ...) {
  cl <- match.call()
  family <- match.arg(family)
  method <- match.arg(method)

  ## Cross-validation: bivariate families <-> neuralbayes
  nb_families <- c("begpd", "bdegpd", "bzidegpd", "mdgpd", "zimdgpd")
  if (family %in% nb_families && method != "neuralbayes")
    stop("family='", family, "' requires method='neuralbayes'", call. = FALSE)
  if (method == "neuralbayes" && !(family %in% nb_families))
    stop("method='neuralbayes' requires a multivariate family (",
         paste(nb_families, collapse = ", "), ")",
         call. = FALSE)

  ## Dispatch to neural Bayes fitting for bivariate families
  if (method == "neuralbayes") {
    estimator <- match.arg(estimator)
    return(.distfit_neuralbayes(x, family = family,
                                 model.path = model.path,
                                 estimator = estimator,
                                 nsamples = as.integer(nsamples),
                                 call = cl,
                                 npe.control = npe.control))
  }

  ## Validate inputs (univariate families only)
  if (!is.numeric(x) || length(x) < 2)
    stop("'x' must be a numeric vector with at least 2 observations")
  if (any(!is.finite(x)))
    stop("'x' must contain only finite observations")
  if (!type %in% 1:6)
    stop("'type' must be an integer from 1 to 6")
  if (method == "bernstein" && family != "egpd")
    stop("Bernstein method is only available for family='egpd'")
  ## Support checks, as in egpd(): an observation outside the support has zero
  ## likelihood under every parameter value, so no fit is meaningful.
  if (family %in% .distfit_discrete_families && any(x != floor(x) | x < 0))
    stop("Discrete family expects non-negative integer data", call.=FALSE)
  if (family %in% c("egpd", "ziegpd", "cpegpd") && any(x < 0))
    stop("family '", family, "' expects non-negative data", call. = FALSE)
  ## The continuous EGPD puts no mass at 0, and its density there is typically 0
  ## or infinite, so one exact zero can pin the fit to the parameters that keep it
  ## finite (type 1: kappa = 1, with a spurious zero SE). Also covers the
  ## Bernstein path.
  if (family == "egpd" && any(x == 0))
    stop("family 'egpd' expects positive data: the continuous EGPD has no mass ",
         "at 0, so exact zeros make its likelihood degenerate. For data with ",
         "exact zeros use family = 'ziegpd' (zero-inflated EGPD), or a discrete ",
         "family ('degpd', 'zidegpd') for counts.", call. = FALSE)
  ## cpegpd is the law of gauge-discretised sums: warn when positive data are off
  ## the cpegpd.h lattice and h is coarse relative to them (D1).
  if (family == "cpegpd") {
    if (!is.numeric(cpegpd.h) || length(cpegpd.h) != 1L || !is.finite(cpegpd.h) ||
        cpegpd.h <= 0)
      stop("'cpegpd.h' must be a single positive finite number", call. = FALSE)
    .aggregate_offgrid_warning(list(positive_samples = list(x[x > 0])), cpegpd.h)
  }
  if (family %in% c("ziegpd", "zidegpd", "zipig", "zigpig", "zibell") && !any(x == 0))
    warning("Zero-inflated family but no zeros observed in data")

  ## Only the outer names identify parameters: coef(pilot)["pi"] is also
  ## a scalar. Its inner name must not turn a fixed `pi` into `pi.pi`.
  if (length(fix.arg) && is.list(fix.arg)) fix.arg <- lapply(fix.arg, unname)
  if (length(start) && is.list(start)) start <- lapply(start, unname)

  ## Dispatch to Bernstein fitting
  if (method == "bernstein") {
    return(.distfit_bernstein(x, type = type, start = start, fix.arg = fix.arg,
                              m = bernstein.m, optim.method = optim.method,
                              hessian = hessian, call = cl, ...))
  }

  ## Parameter specification
  spec <- .distfit_parspec(type, family)

  ## Remove fixed args from spec
  free_spec <- spec
  if (!is.null(fix.arg)) {
    .distfit_check_parlist(fix.arg, "fix.arg", spec, family)
    bad <- setdiff(names(fix.arg), names(spec))
    if (length(bad) > 0)
      stop("Unknown fixed parameters: ", paste(bad, collapse = ", "))
    free_spec <- spec[setdiff(names(spec), names(fix.arg))]
  }
  if (length(free_spec) == 0)
    stop("No free parameters to estimate")

  fixed <- as.list(fix.arg)
  if (family %in% c("bell", "zibell") && !any(x > 0) && "theta" %in% names(free_spec))
    stop("No interior Bell maximum for an all-zero sample; the likelihood is ",
         "maximized only in the boundary limit.", call. = FALSE)
  reduction <- .distfit_reduction(type, family, free_spec, fixed)
  product <- isTRUE(reduction$product)
  binary <- isTRUE(reduction$bernoulli)
  thinned <- isTRUE(reduction$thinned)
  start_vals <- if (binary || thinned) list() else .distfit_start(x, type, family, free_spec, fixed)
  if (!is.null(start)) {
    adds <- if (product || thinned) list(mu = list(transform = "log")) else if (binary)
      list(p1 = list(transform = "logit"))
    .distfit_check_parlist(start, "start", c(spec, adds), family)
    bad <- setdiff(names(start), c(names(free_spec), names(adds)))
    if (length(bad) > 0)
      stop("Start values for fixed or unknown parameters: ", paste(bad, collapse = ", "))
    if (product && "mu" %in% names(start) && any(c("b", "c") %in% names(start)))
      stop("Specify the Poisson mean start `mu` or its factors `b`, `c`, not both.")
    for (nm in names(start)) start_vals[[nm]] <- start[[nm]]
  }
  if (!is.null(reduction)) {
    fixed <- c(fixed, reduction$constants)
    free_spec <- free_spec[setdiff(names(free_spec), names(reduction$constants))]
    if (product) {
      start_vals$mu <- start_vals$mu %||% (start_vals$b * start_vals$c)
      names(free_spec)[names(free_spec) == "b"] <- "mu"
    }
    start_vals <- start_vals[names(free_spec)]
  }
  if (binary) {
    if (any(x > 1)) stop("Observation outside the fixed two-cell support {0, 1}.", call. = FALSE)
    free_spec <- list(p1 = list(transform = "logit"))
    start_vals <- list(p1 = mean(x))
  }
  if (thinned) {
    free_spec <- list(mu = list(transform = "log"))
    start_vals <- list(mu = mean(x))
  }

  ## Transform to unconstrained scale
  start_par <- unlist(start_vals)
  theta0 <- .par_to_theta(start_par, free_spec)
  dots <- list(...)
  continuous <- family %in% c("egpd", "ziegpd") && any(x > 0) &&
    !isTRUE(fixed$pi == 1) && !(isTRUE(fixed$pi == 0) && any(x == 0))
  if (continuous && !is.null(fixed$xi))
    .distfit_check_endpoint(max(x), fixed$xi, fixed$sigma, names(free_spec), dots,
                            monotone = .distfit_monotone_endpoint(type, fixed))

  ## Optimize. Nelder-Mead (the default) is unreliable and noisy for a single
  ## free parameter (e.g. the one-parameter Bell family), so fall back to BFGS
  ## in that case unless the user has explicitly chosen a method.
  opt_method <- .distfit_method(optim.method, length(theta0), dots)
  h_cp <- if (family == "cpegpd") cpegpd.h else NULL
  weight2 <- type == 6L && identical(names(free_spec), "prob") &&
    family %in% c("degpd", "zidegpd", "cpdegpd") && isTRUE(fixed$xi < 0) &&
    .gpd_last_count(fixed$sigma, fixed$xi) == 1
  profile.scale <- type == 1L && (family == "degpd" ||
    (family == "zidegpd" && isTRUE(fixed$pi == 0))) &&
    isTRUE(fixed$xi == -1) && "sigma" %in% names(free_spec) &&
    min(x) < max(x) && max(x) < 2^53 && (dots$control$fnscale %||% 1) > 0
  ## A start at which some observation lies outside the support has zero likelihood.
  ## Move it into the support first, on a surrogate objective that charges such an
  ## observation -log(DBL_MIN) instead of +Inf; the fit itself then uses the exact
  ## likelihood, so it can never settle on parameters that exclude data.
  if (!binary && !thinned && !weight2 && !profile.scale && length(theta0) && family %in% c("egpd", "degpd", "ziegpd", "zidegpd", "cpegpd", "cpdegpd") &&
      .distfit_nll(theta0, x, type, family, fixed, free_spec, h_cp) >= 1e20) {
    pre <- optim(theta0, fn = .distfit_nll, x = x, type = type, family = family,
                 fix.arg = fixed, spec = free_spec, h = h_cp, soft = TRUE,
                 method = opt_method, ...)
    theta0 <- pre$par
  }
  ## Nelder-Mead can stop early on a collapsed simplex (e.g. ZIDEGPD: 0.03 below the
  ## optimum); one restart from its solution, with a fresh simplex, is the standard
  ## remedy, never returns a worse value, and is cheap here.
  nm <- opt_method == "Nelder-Mead"
  profile.pi <- "pi" %in% names(free_spec)
  run <- function(theta, spec = free_spec, pins = fixed, args = dots) {
    method <- .distfit_method(optim.method, length(theta), args)
    if (profile.scale)
      return(.distfit_profile_power_scale(theta, x, family, pins, spec, method, args))
    if ("pi" %in% names(spec))
      return(.distfit_profile_pi(theta, x, type, family, pins, spec, method, args))
    do.call(.distfit_optim, c(list(par = theta, fn = .distfit_nll, x = x, type = type,
                          family = family, fix.arg = pins, spec = spec, h = h_cp,
                          method = method), args))
  }
  opt <- if (binary || thinned) {
    .distfit_reduced_opt(x, theta0, binary, hessian, dots)
  } else if (weight2) {
    lo <- dots$lower %||% -Inf
    hi <- dots$upper %||% Inf
    if (length(lo) != 1L || length(hi) != 1L || anyNA(c(lo, hi)) ||
        lo > hi || lo == Inf || hi == -Inf)
      stop("Invalid optimizer bounds for the two-cell mixture.", call. = FALSE)
    if (family != "cpdegpd" && any(x > 1))
      stop("Observation outside the fixed two-cell support {0, 1}.", call. = FALSE)
    w <- two_cell_weight_cpp(sum(x), length(x), fixed$sigma, fixed$xi, fixed$kappa,
                             fixed$delta, fixed$lambda %||% 1, fixed$pi %||% 0,
                             family == "cpdegpd")
    p <- c(prob = min(hi, max(lo, w[1L])))
    list(par = p, value = .distfit_nll(p, x, type, family, fixed, free_spec),
         convergence = 0L, message = "Closed-form two-cell mixture weight",
         carrier.boundary = if (is.infinite(p)) plogis(p), weight.se = w[2L])
  } else if (length(theta0)) run(theta0) else
    list(par = theta0, value = .distfit_nll(theta0, x, type, family, fixed, free_spec, h_cp),
         convergence = 0L)
  if (!binary && !thinned && !weight2 && length(theta0) && nm && is.finite(opt$value) && opt$value < 1e20)
    opt <- run(opt$par)
  ## A small gradient on the equal-power ridge need not be a maximum.
  ## Refit only if separating the two powers demonstrates a better
  ## *actual* likelihood, not merely an indefinite numerical Hessian.
  if (type == 6 && opt$convergence == 0L && all(c("kappa", "delta") %in% names(opt$par)) &&
      abs(opt$par[["kappa"]] - opt$par[["delta"]]) < .02) {
    p <- fixed$prob %||% plogis(opt$par[["prob"]])
    split <- function(h) {
      th <- opt$par
      th[c("kappa","delta")] <- th[c("kappa","delta")] + h * c(-(1-p), p)
      th
    }
    candidates <- lapply(c(-1, -.1, .1, 1), split)
    nlls <- vapply(candidates, .distfit_nll, numeric(1), x = x, type = type,
                   family = family, fix.arg = fixed, spec = free_spec, h = h_cp)
    i <- which.min(nlls)
    if (nlls[i] + 1e-8 < opt$value) {
      separated <- run(candidates[[i]])
      if (separated$value < opt$value) opt <- separated
    }
  }
  ## Brent calls its objective and returns its solution as unnamed scalars.
  if (is.null(names(opt$par))) names(opt$par) <- names(free_spec)
  if (!binary && !thinned && length(opt$par)) {
    nll <- function(p) .distfit_nll(p, x, type, family, fixed, free_spec, h_cp)
    nm_face <- if (type == 2L) "kappa" else if (type == 6L) "prob"
    if (!weight2 && length(nm_face) && nm_face %in% names(free_spec) && opt$convergence == 0L)
      opt <- .distfit_carrier_face(opt, nm_face, free_spec, fixed, dots, run, nll)
    opt <- .distfit_check_box(opt, if (profile.pi || weight2 || profile.scale) "profile" else opt_method,
                              dots, hessian, nll)
    if (length(opt$carrier.boundary)) {
      opt$bound.active <- TRUE
      if (hessian)
        warning("Carrier MLE (", names(opt$carrier.boundary), ") is on its natural boundary; ",
                "Wald standard errors unavailable.", call. = FALSE)
    }
  }
  if (profile.pi && is.infinite(opt$par[["pi"]])) {
    opt$pi.boundary <- TRUE
    opt$bound.active <- TRUE
    if (hessian)
      warning("Zero-inflation MLE is on its natural boundary; Wald standard errors unavailable.",
              call. = FALSE)
  }
  .distfit_check_optim(opt, opt_method)
  ## A pmf is at most one, so a positive log-likelihood for a discrete family means
  ## that its numerics broke down where the optimiser stopped; as in egpd(), this
  ## is not a fit (the tolerance allows for rounding in the per-observation terms).
  if (family %in% .distfit_discrete_families &&
      -opt$value > sqrt(.Machine$double.eps) * max(1, length(x)))
    stop(sprintf(paste0("the log-likelihood (%.6g) is positive, which is impossible ",
                        "for a discrete family: the likelihood is numerically ",
                        "unreliable at the returned estimates."), -opt$value),
         call. = FALSE)

  ## Back-transform estimates
  estimate <- .theta_to_par(opt$par, free_spec)
  if (profile.scale) estimate["sigma"] <- opt$profile.scale
  par <- c(as.list(estimate), fixed)
  if (continuous)
    .distfit_check_endpoint(max(x), par$xi, par$sigma, names(free_spec), dots,
                            monotone = .distfit_monotone_endpoint(type, par))
  if (continuous && type == 3L && "kappa" %in% names(free_spec) && opt$convergence == 0L) {
    lo <- rep_len(dots$lower %||% -Inf, length(estimate))[match("kappa", names(estimate))]
    if (lo == -Inf) {
      ## The carrier density is proportional to exp((kappa-1)*L).
      ## Its per-row score is L-E[L] and decreases strictly at rate -Var[L].
      s <- .gpd_cdf_surv_std(x[x > 0], par$sigma, par$xi)$surv
      score0 <- mean(log1p(-((15/16)*s)^2)) - .egpd_tb_norm(0)[1L, 2L]
      if (score0 < -1e-12)
        stop("No interior type-3 kappa maximum at the fitted scale and shape: ",
             "the continuous likelihood increases toward the excluded kappa = 0 limit.",
             call. = FALSE)
    }
  }

  ## K cells identify at most K-1 masses (or positive-jump rates).
  ## Count the support of the fitted law, not the distinct observations.
  support.ridge <- FALSE
  if (family %in% c("degpd", "zidegpd", "cpdegpd") && any(c("sigma", "xi") %in% names(estimate))) {
    if (par$xi < 0) {
      cells <- 1 + .gpd_last_count(par$sigma, par$xi)
      support.ridge <- length(estimate) > cells - 1
    }
  }
  if (support.ridge)
    warning(if (family == "cpdegpd") "Fitted jump support" else "Fitted support",
            " has ", cells, " cell(s): the free parameters are not separately ",
            "identifiable. Estimates are non-unique; Wald standard errors are unavailable.",
            call. = FALSE)

  ## At e=max(x)+1 the mass in the highest observed cell starts losing
  ## a term of order (e-max(x)-1)^(-1/xi). For -1/xi <= 2 there is no
  ## common second derivative across the two sides. Do not invert a
  ## numerical Hessian whose stencil crosses that change of support.
  support.change <- FALSE
  if (hessian && !support.ridge && !isTRUE(opt$bound.active) &&
      family %in% c("degpd", "zidegpd", "cpdegpd", "cpegpd") && par$xi <= -.5 &&
      any(c("sigma", "xi") %in% names(estimate))) {
    step <- rep_len(dots$control$ndeps %||% .001, length(estimate)) *
      rep_len(dots$control$parscale %||% 1, length(estimate))
    names(step) <- names(estimate)
    ds <- if ("sigma" %in% names(step)) 2*max(step["sigma"], sqrt(.Machine$double.eps)) else 0
    dx <- if ("xi" %in% names(step)) 2*max(step["xi"], abs(par$xi)*sqrt(.Machine$double.eps)) else 0
    edge <- max(x) + 1
    if (family %in% c("cpdegpd", "cpegpd")) {
      ## A compound sum need not be bounded. Any moving *jump* cell
      ## contributing to an observed total can kink its probability.
      grid <- if (family == "cpegpd") cpegpd.h else 1
      lo <- par$sigma*exp(-ds)/(-par$xi+dx)
      k <- max(1, ceiling(lo/grid))
      limit <- if (family == "cpegpd") ceiling(max(x)/grid) else max(x)+1
      edge <- if (k <= limit) grid*k else NA_real_
    }
    support.change <- is.finite(edge) && par$sigma*exp(-ds) <= (-par$xi+dx)*edge &&
      par$sigma*exp(ds) >= (-par$xi-dx)*edge
    if (support.change)
      warning("The observed-information stencil crosses a change in count support; ",
              "ordinary Wald inference is unavailable.", call. = FALSE)
  }

  ## SE via delta method
  se <- rep(NA_real_, length(estimate))
  names(se) <- names(estimate)
  V_par <- if (hessian && !length(estimate)) matrix(numeric(), 0L, 0L) else NULL

  if (hessian && weight2 && !isTRUE(opt$bound.active)) {
    se[] <- opt$weight.se
    V_par <- matrix(se^2, 1L, 1L, dimnames = list("prob", "prob"))
    opt$hessian <- matrix((plogis(opt$par) * plogis(-opt$par) / se)^2, 1L, 1L,
                          dimnames = dimnames(V_par))
  } else if (hessian && length(estimate) && !support.ridge && !support.change &&
             !isTRUE(opt$bound.active) && !isTRUE(opt$stalled)) {
    if (profile.pi && length(estimate) == 1L) {
      ## At an interior maximum: (1-pi)*Pr_base(Y>0) = n_positive/n.
      opt$hessian <- matrix(sum(x > 0) * estimate[["pi"]]^2 / mean(x == 0),
                            1L, 1L, dimnames = list("pi", "pi"))
    } else if (!binary && !thinned)
      opt$hessian <- stats::optimHess(opt$par, .distfit_nll, gr = dots[["gr"]], x = x,
          type = type, family = family, fix.arg = fixed, spec = free_spec, h = h_cp,
          control = dots$control %||% list())
    R_chol <- tryCatch(chol(opt$hessian), error = function(e) NULL)
    if (!is.null(R_chol)) {
      V_theta <- chol2inv(R_chol)
      jac <- .distfit_jacobian(opt$par, free_spec)
      J <- diag(jac, nrow = length(jac))
      V_par <- J %*% V_theta %*% t(J)
      rownames(V_par) <- colnames(V_par) <- names(estimate)
      se <- sqrt(diag(V_par))
    } else {
      warning("Hessian is not positive definite; standard errors unavailable.",
              call. = FALSE)
    }
  }

  ## Assemble result
  npar <- length(estimate)
  ll <- -opt$value
  n <- length(x)

  structure(list(
    estimate  = estimate,
    sd        = se,
    vcov      = V_par,
    loglik    = ll,
    aic       = -2 * ll + 2 * npar,
    bic       = -2 * ll + log(n) * npar,
    n         = n,
    npar      = npar,
    data      = x,
    type      = type,
    family    = family,
    method    = method,
    fix.arg   = fix.arg,
    reduction = reduction,
    support.ridge = support.ridge,
    support.change = support.change,
    convergence = opt$convergence,
    optim     = opt,
    call      = cl,
    bernstein.m       = NULL,
    bernstein.weights = NULL,
    cpegpd.h  = if (family == "cpegpd") cpegpd.h else NULL
  ), class = "distfit")
}


#' @rdname distfit
#' @export
fitegpd <- function(x, type = 1,
                    family = c("egpd", "degpd", "ziegpd", "zidegpd",
                               "cpegpd", "cpdegpd", "pig", "zipig",
                               "gpig", "zigpig", "bell", "zibell", "gwaring", "begpd",
                               "bdegpd", "bzidegpd",
                               "mdgpd", "zimdgpd"),
                    method = c("mle", "bernstein", "neuralbayes"),
                    start = NULL, fix.arg = NULL,
                    optim.method = "Nelder-Mead", hessian = TRUE,
                    bernstein.m = 8, cpegpd.h = 0.2,
                    model.path = NULL, estimator = c("npe", "nbe"),
                    nsamples = 1000L, npe.control = list(), ...) {
  .Deprecated("distfit", package = "egpd",
              msg = paste("'fitegpd' has been renamed to 'distfit' and is",
                          "deprecated; it will be removed in a future release."))
  cl <- match.call()
  cl[[1L]] <- quote(distfit)
  eval.parent(cl)
}


## ---- Internal helpers ----

.distfit_method <- function(method, npar, dots) {
  method <- match.arg(method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
  if (method != "Brent" && (any(is.finite(dots$lower)) || any(is.finite(dots$upper))))
    return("L-BFGS-B")
  if (npar == 1L && method == "Nelder-Mead") "BFGS" else method
}

## L-BFGS-B's relative-reduction test can report success when a line search
## hits the finite penalty for impossible data and returns to its starting point.
.distfit_optim <- function(par, fn, gr = NULL, ..., method, control = list(),
                           lower = -Inf, upper = Inf) {
  run <- function(p, ctl, objective = fn)
    stats::optim(p, objective, gr, ..., method = method, control = ctl,
                 lower = lower, upper = upper)
  opt <- run(par, control)
  if (method != "L-BFGS-B" || (control$fnscale %||% 1) <= 0) return(opt)
  n <- length(par)
  lo <- rep_len(lower, n)
  hi <- rep_len(upper, n)
  scale <- rep_len(control$parscale %||% 1, n)
  step <- scale * rep_len(control$ndeps %||% .001, n)
  reltol <- max(16, control$factr %||% 1e7) * .Machine$double.eps
  for (attempt in 0:2) {
    if (opt$convergence != 0L || !is.finite(opt$value) || opt$value >= 1e20) return(opt)
    best <- opt$par
    value <- opt$value
    for (j in seq_len(n)) for (sign in c(-1, 1)) {
      p <- opt$par
      p[j] <- min(hi[j], max(lo[j], p[j] + sign * step[j]))
      if (p[j] == opt$par[j]) next
      v <- fn(p, ...)
      if (is.finite(v) && v < value) {
        best <- p
        value <- v
      }
    }
    if (opt$value - value <= reltol * max(control$fnscale %||% 1, abs(opt$value)))
      return(opt)
    opt$par <- best
    opt$value <- value
    if (attempt == 2L) break
    control$parscale <- scale / 4^(attempt + 1L)
    ## An impossible trial must lose to this feasible iterate, but a jump to
    ## 1e20 can collapse the interpolated line-search step to zero again.
    ## Leave every feasible likelihood untouched and verify the returned point.
    penalty <- value + max(1, abs(value))
    search <- function(p, ...) {
      v <- fn(p, ...)
      if (is.finite(v) && v < 1e20) v else penalty
    }
    trial <- run(best, control, search)
    trial$value <- fn(trial$par, ...)
    if (!is.finite(trial$value) || trial$value > value) break
    opt <- trial
  }
  opt$convergence <- 51L
  opt$stalled <- TRUE
  opt$message <- "A feasible step still improved the likelihood; adjust control$parscale."
  opt
}

## Keep optimizer controls in the original coordinate order when profiling.
.distfit_subset_controls <- function(dots, keep, n) {
  for (nm in c("lower", "upper"))
    if (!is.null(dots[[nm]])) dots[[nm]] <- rep_len(dots[[nm]], n)[keep]
  for (nm in c("parscale", "ndeps")) {
    v <- dots$control[[nm]]
    if (!is.null(v)) {
      if (length(v) != n) stop("Wrong length for '", nm, "'.", call. = FALSE)
      dots$control[[nm]] <- v[keep]
    }
  }
  dots
}

.distfit_profile_power_scale <- function(theta, x, family, fixed, spec, method, dots) {
  i <- match("sigma", names(spec))
  n <- length(theta)
  lo <- rep_len(dots$lower %||% -Inf, n)[i]
  hi <- rep_len(dots$upper %||% Inf, n)[i]
  if (anyNA(c(lo, hi)) || lo > hi || lo == Inf || hi == -Inf)
    stop("Invalid optimizer bounds for the scale.", call. = FALSE)
  if (method == "Brent" && n != 1L)
    stop("Brent requires a single free parameter.", call. = FALSE)
  m <- max(x)
  h <- -log1p(-mean(x == m))
  ## For m<s<=m+1, counts below m have total probability (m/s)^k
  ## and their conditional distribution is independent of s. Above m+1
  ## every observed mass decreases with s. Clamp the binomial optimum.
  scale <- function(v) {
    k <- fixed$kappa %||% exp(v[["kappa"]])
    a <- h/k
    s <- if (a >= log1p(1/m)) m+1 else m + m*expm1(a)
    if (s == m) {
      ## The exact optimum is above m; m itself assigns zero probability
      ## to the observed last cell. Use the smallest representable scale above it.
      bin <- 2^floor(log2(m))
      if (bin > m) bin <- bin/2
      s <- m + bin*.Machine$double.eps
    }
    if (log(s) < lo) exp(lo) else if (log(s) > hi) exp(hi) else s
  }
  fn <- function(v) .distfit_nll(v, x, 1L, family, c(fixed, list(sigma = scale(v))), spec[-i])
  full <- function(v) {
    theta[-i] <- v
    theta[i] <- log(scale(v))
    theta
  }
  v <- theta[-i]
  if (length(v)) {
    args <- .distfit_subset_controls(dots, -i, n)
    if (!is.null(gr <- dots$gr))
      args$gr <- function(v) gr(full(v), x = x, type = 1L, family = family,
                                fix.arg = fixed, spec = spec, h = NULL)[-i]
    fit <- do.call(.distfit_optim, c(list(par = v, fn = fn,
                     method = .distfit_method(method, length(v), args)), args))
    v <- fit$par
  } else {
    fit <- list(value = fn(v), convergence = 0L, message = "Closed-form power-carrier scale")
  }
  fit$par <- full(v)
  ## exp(log(m+1)) need not be m+1. Retain the same natural scale used in fn.
  fit$profile.scale <- scale(v)
  fit
}

## Conditional on the base law, the number of positive observations is
## binomial with success probability (1-pi)*Pr_base(Y>0). Maximize this
## factor exactly, including pi=0, while fitting the remaining parameters.
## Return the original working coordinates for the joint Hessian and methods.
.distfit_profile_pi <- function(theta, x, type, family, fixed, spec, method, dots) {
  i <- match("pi", names(spec))
  n <- length(theta)
  lo <- rep_len(dots$lower %||% -Inf, n)
  hi <- rep_len(dots$upper %||% Inf, n)
  if (anyNA(c(lo, hi)) || any(lo > hi) || any(lo == Inf | hi == -Inf))
    stop("Invalid optimizer bounds for the zero-inflated fit.", call. = FALSE)
  if (method == "Brent" && n != 1L)
    stop("Brent requires a single free parameter.", call. = FALSE)
  z <- sum(x == 0)
  zero <- z / length(x)
  obs <- log(mean(x > 0))
  if (z == length(x) && hi[i] == Inf) {
    if (family %in% c("zipig", "zigpig", "zibell"))
      stop("All observations are zero: the likelihood increases to the excluded pi = 1 boundary.",
           call. = FALSE)
    if (n > 1L)
      stop("All observations are zero: at pi = 1 the free base parameters are not identifiable.",
           call. = FALSE)
  }
  positive <- vapply(spec[-i], function(s) s$transform == "log", logical(1))
  full <- function(v) {
    theta[-i] <- v
    base <- .theta_to_par(theta[-i], spec[-i])
    if (any(!is.finite(base)) || any(base[positive] <= 0) || any(base[names(base) == "a"] <= 0)) {
      theta[i] <- NA_real_  # infeasible optimizer proposal
      return(theta)
    }
    if (z == length(x)) mode <- Inf else if (!z) mode <- -Inf else {
      par <- c(as.list(base), fixed)
      ## Compare the smaller observed cell directly. When zeros are rare,
      ## even an exact equality can be lost by a log/complement round trip.
      if (family == "zidegpd") {
        mass <- do.call(pdiscegpd, c(list(q = 0, type = type, lower.tail = zero <= .5,
                                         log.p = zero > .5), par))
      } else {
        lz <- switch(family, ziegpd = -Inf,
          zibell = dbell(0, par$theta, log = TRUE),
          zipig = dpig(0, par$mu, par$sigma, log = TRUE),
          zigpig = dgpig(0, par$a, par$mu %||% par$b, par$c, log = TRUE))
        mass <- if (zero <= .5) exp(lz) else .log1mexp(lz)
      }
      if (zero <= .5) {
        ## pi/(1-pi) = (observed zero fraction - base zero mass)/(1-zero).
        mode <- log(pmax(0, zero - mass)) - log1p(-zero)
      } else {
        l1mpi <- obs - mass
        mode <- if (is.na(l1mpi)) NA_real_ else if (l1mpi >= 0) -Inf else .log1mexp(l1mpi) - l1mpi
      }
      if (family == "zidegpd" && is.finite(mass) &&
          (if (zero <= .5) abs(zero-mass) <= 16*.Machine$double.eps*zero else
             abs(obs-mass) <= 16*.Machine$double.eps*max(1,abs(obs)))) {
        power <- switch(as.integer(type), par$kappa, if (par$kappa == 0) 1,
                        if (par$kappa == 1) 1, if (par$delta == 1) 2,
                        if (par$delta == 1) par$kappa, par$kappa)
        if (!is.null(power)) mode <- profile_zero_odds_cpp(z, length(x), par$sigma, par$xi,
          power, if (type == 6L) par$delta else power, if (type == 6L) par$prob else 1)
      }
    }
    theta[i] <- min(max(mode, lo[i]), hi[i])
    theta
  }
  fn <- function(v) {
    p <- full(v)
    if (is.na(p[i])) return(1e20)
    .distfit_nll(p, x, type, family, fixed, spec)
  }
  if (n == 1L) {
    p <- full(numeric())
    return(list(par = p, value = .distfit_nll(p, x, type, family, fixed, spec),
                convergence = 0L, message = "Closed-form zero-inflation MLE"))
  }
  dots <- .distfit_subset_controls(dots, -i, n)
  if (!is.null(gr <- dots$gr)) {
    dots$gr <- function(v) gr(full(v), x = x, type = type, family = family,
                               fix.arg = fixed, spec = spec, h = NULL)[-i]
  }
  method <- .distfit_method(method, n - 1L, dots)
  opt <- do.call(.distfit_optim, c(list(par = theta[-i], fn = fn, method = method), dots))
  opt$par <- full(opt$par)
  opt
}

.distfit_check_box <- function(opt, method, dots, hessian, nll) {
  n <- length(opt$par)
  lo <- rep_len(dots$lower %||% -Inf, n)
  hi <- rep_len(dots$upper %||% Inf, n)
  if (!any(is.finite(c(lo, hi)))) return(opt)
  if (method == "Brent") {
    ## optimize() searches the open interval. Compare its solution to
    ## the two feasible faces rather than mistaking a nearby optimum for
    ## a boundary solely because it is within the optimizer's tolerance.
    faces <- c(lo, hi)
    val <- vapply(faces, function(x) nll(setNames(x, names(opt$par))), 0)
    tol <- sqrt(.Machine$double.eps)*abs(opt$par) + (dots$control$reltol %||% sqrt(.Machine$double.eps))
    take <- which(is.finite(val) & (val < opt$value |
                    (val == opt$value & abs(faces - opt$par) <= 2*tol)))
    if (length(take)) {
      j <- take[which.min(val[take])]
      opt$par[] <- faces[j]
      opt$value <- val[j]
    }
  }
  on_face <- function(b) is.finite(b) & abs(opt$par - b) <= 8*.Machine$double.eps*abs(b)
  opt$bound.active <- any(on_face(lo) | on_face(hi))
  if (hessian && opt$bound.active)
    warning("MLE is on an optimizer bound; Wald standard errors unavailable.", call. = FALSE)
  opt
}

## Log and logit coordinates never reach these included endpoints at a finite
## parameter value. Compare actual likelihoods, not distances in either scale.
.distfit_carrier_face <- function(opt, nm, spec, fixed, dots, run, nll) {
  theta <- opt$par
  j <- match(nm, names(theta))
  n <- length(theta)
  lo <- rep_len(dots$lower %||% -Inf, n)
  hi <- rep_len(dots$upper %||% Inf, n)
  faces <- c(if (lo[j] == -Inf) 0, if (nm == "prob" && hi[j] == Inf) 1)
  for (face in faces) {
    p <- theta
    p[j] <- if (face == 0) -Inf else Inf
    f0 <- nll(p)
    tol <- 16 * .Machine$double.eps * max(1, abs(opt$value))
    if (!is.finite(f0) || f0 >= 1e20 || f0 > opt$value + tol) next

    ## An inactive power stays only as a representative for evaluating the
    ## fitted mixture. Never ask an optimizer to estimate it on this face.
    idle <- if (nm == "prob") if (face == 0) "kappa" else "delta"
    keep <- !(names(p) %in% c(nm, idle))
    if (any(keep)) {
      pins <- c(fixed, as.list(.theta_to_par(p[!keep], spec[!keep])))
      args <- dots
      if (!is.null(args$lower)) args$lower <- lo[keep]
      if (!is.null(args$upper)) args$upper <- hi[keep]
      for (s in c("parscale", "ndeps")) if (!is.null(args$control[[s]]))
        args$control[[s]] <- rep_len(args$control[[s]], n)[keep]
      ## A supplied gradient is for the original finite working coordinates.
      args$gr <- NULL
      fit <- run(p[keep], spec[keep], pins, args)
      if (fit$convergence != 0L || !is.finite(fit$value) || fit$value > f0 + tol) next
      p[keep] <- fit$par
      f0 <- nll(p)
    }
    if (f0 > opt$value + tol) next

    ## A small *interior* optimum can beat the endpoint even if it beats
    ## the log-scale optimizer too. Check the one-sided derivative in the
    ## natural parameter, and do not call a completely flat slice a face.
    width <- if (nm == "kappa") exp(hi[j]) else if (face == 0)
      plogis(hi[j]) else plogis(lo[j], lower.tail = FALSE)
    h <- min(.Machine$double.eps^0.25, width / 2)
    if (h <= 0) next
    inward <- function(u) {
      p[j] <- if (nm == "kappa") log(u) else qlogis(u, lower.tail = face == 0)
      nll(p) - f0
    }
    d <- vapply(c(h/2, h, 2*h), inward, numeric(1))
    ## Two second-order forward derivatives, in units of h. Their difference
    ## estimates truncation error; tol covers subtraction of total log likelihoods.
    a <- 4*d[1L] - d[2L]
    b <- 2*d[2L] - d[3L]/2
    if (min(d) < -tol || max(d) <= tol || a < -(abs(a-b)/3 + 6*tol)) next
    opt$par <- p
    opt$value <- f0
    opt$carrier.boundary <- setNames(face, nm)
  }
  opt
}

## At xi=-1, f(x;s)=g(x/s)/s increases towards the endpoint whenever
## u*g(u) is increasing. This holds for powers/mixtures and the increasing
## densities of types 2 and 4. In type 3, z=l+(1/2-l)u: even for k<1,
## (z-l)(1-2z)/(z(1-z)) < 1 gives d log(u*g)/d log(u) > 0.
## For delta<=1, G4 is a positive sum of powers, so u*G4'/G4 is
## nondecreasing (a weighted mean of the exponents); so is u*(G4^r)'.
.distfit_monotone_endpoint <- function(type, par) {
  type != 5L || isTRUE(par$delta <= 1) || isTRUE(par$kappa >= 2)
}

## At xi < -1 the continuous density diverges at its upper endpoint.
## Check whether a free coordinate can take that endpoint to xmax,
## respecting any box on the optimizer's (log sigma, identity xi) scales.
.distfit_check_endpoint <- function(xmax, xi, sigma, free, dots, monotone = FALSE) {
  if (!isTRUE(xi <= -1)) return(invisible(NULL))
  reachable <- function(nm, boundary) {
    j <- match(nm, free)
    if (is.na(j) || !is.finite(boundary)) return(FALSE)
    lo <- rep_len(dots$lower %||% -Inf, length(free))[j]
    hi <- rep_len(dots$upper %||% Inf, length(free))[j]
    isTRUE(lo <= boundary && boundary < hi)
  }
  if (xi == -1) {
    if (monotone && reachable("sigma", log(xmax)))
      stop("Nonregular continuous endpoint: at xi = -1 the likelihood ",
           "increases as scale decreases to the largest observation from above. ",
           "There is no interior scale maximum or Wald inference.", call. = FALSE)
    return(invisible(NULL))
  }
  if (reachable("sigma", log(-xi * xmax)) ||
      (isTRUE(sigma > xmax) && reachable("xi", -sigma / xmax)))
    stop("Unbounded continuous likelihood: with xi < -1 the density diverges ",
         "as the upper endpoint approaches the largest observation. No finite ",
         "maximum or Wald inference is available on this branch.", call. = FALSE)
}

## Both reduced likelihoods are concave in their one working coordinate.
.distfit_reduced_opt <- function(x, mode, binary, hessian, dots) {
  lo <- dots[["lower"]] %||% -Inf
  hi <- dots[["upper"]] %||% Inf
  if (!is.numeric(lo) || !is.numeric(hi) || length(lo) != 1L || length(hi) != 1L ||
      is.na(lo) || is.na(hi) || lo >= hi)
    stop("Reduced optimizer 'lower' and 'upper' must be increasing scalar bounds.", call. = FALSE)
  theta <- pmin(pmax(mode, lo), hi)
  if (!is.finite(theta))
    stop(if (binary) "No interior Bernoulli estimate: only one support cell was observed."
         else "No positive Poisson mean estimate: all aggregate counts are zero.", call. = FALSE)
  v <- if (binary) plogis(theta) else exp(theta)
  if (!is.finite(v) || v <= 0 || (binary && v >= 1))
    stop("Reduced optimizer bound has no representable interior parameter value.", call. = FALSE)
  at_bound <- theta == lo || theta == hi
  if (hessian && at_bound)
    warning("Reduced MLE is on an optimizer bound; Wald standard errors unavailable.", call. = FALSE)
  list(par = theta, value = if (binary) -sum(stats::dbinom(x, 1, v, log = TRUE))
                            else -sum(stats::dpois(x, v, log = TRUE)),
       hessian = if (hessian && !at_bound) matrix(length(x) * v * (if (binary) 1 - v else 1)),
       convergence = 0L, bound.active = at_bound,
       message = if (binary) "Closed-form Bernoulli MLE" else "Closed-form thinned Poisson MLE")
}

## Remove only algebraically redundant coordinates, not directions that
## happen to have little observed curvature in a particular sample.
.distfit_reduction <- function(type, family, free, fixed) {
  gp <- family %in% c("gpig", "zigpig")
  mix <- type == 6 && "prob" %in% c(names(free), names(fixed))
  end <- Inf
  if (family %in% c("degpd", "zidegpd", "cpdegpd") && isTRUE(fixed$xi < 0) && !is.null(fixed$sigma))
    end <- .gpd_last_count(fixed$sigma, fixed$xi) + 1
  if ((gp && isTRUE(fixed$c == 0)) ||
      (family %in% c("ziegpd", "zidegpd") && isTRUE(fixed$pi == 1)) || end <= 1) {
    constants <- lapply(free, function(s) switch(s$transform, log = 1, logit = .5, identity = 0))
    return(list(constants = constants, description = "Point mass at zero; no identifiable free parameters."))
  }
  if (gp && isTRUE(fixed$a == 1) && all(c("b","c") %in% names(free)))
    return(list(constants = list(c = 1), product = TRUE,
                description = "Poisson component: only mu = b*c is identifiable."))
  red <- if (mix && "prob" %in% names(free) && isTRUE(fixed$kappa == fixed$delta))
    list(constants = list(prob = .5), description = "Equal fixed powers: mixture weight not identifiable.")
  else if (mix && isTRUE(fixed$prob == 0) && "kappa" %in% names(free))
    list(constants = list(kappa = 1), description = "First component has zero weight; kappa not identifiable.")
  else if (mix && isTRUE(fixed$prob == 1) && "delta" %in% names(free))
    list(constants = list(delta = 1), description = "Second component has zero weight; delta not identifiable.")
  if (end <= 2 && length(setdiff(names(free), names(red$constants))) > 1L) {
    if (family == "cpdegpd" && "lambda" %in% names(free))
      return(list(thinned = TRUE, description = paste("Two-cell jumps: the aggregate is Poisson;",
                  "only mu = lambda * P(jump = 1) is identifiable.")))
    if (type == 1 && family == "zidegpd")
      return(list(bernoulli = TRUE,
                  description = "Two support cells: only p1 = P(Y = 1) is identifiable, not kappa and pi separately."))
    if (family == "cpdegpd")
      stop("Only the jump probability is identifiable with fixed intensity and two-cell jumps; ",
           "fix all but one effective carrier parameter.", call. = FALSE)
    stop("Only P(Y = 1) is identifiable on the fixed two-cell support; ",
         "fix all but one effective carrier or zero-inflation parameter.", call. = FALSE)
  }
  cells <- ceiling(end)
  if (length(setdiff(names(free), names(red$constants))) > cells - 1)
    stop("Fixed ", if (family == "cpdegpd") "jump ", "support has ", cells,
         " cells: at most ", cells - 1, " parameter combinations are identifiable. ",
         "Fix additional parameters.", call. = FALSE)
  red
}

## Univariate families whose observations are counts (a pmf, so log-lik <= 0).
.distfit_discrete_families <- c("degpd", "zidegpd", "cpdegpd", "pig", "zipig",
                                "gpig", "zigpig", "bell", "zibell", "gwaring")

## Check a user-supplied `start` or `fix.arg`: a list (or numeric vector) whose
## values all carry the unique name of their parameter, each a single finite number
## inside that parameter's space -- positive on the log scale (rho > 1 for the
## generalised Waring, whose mean needs it), inside (0, 1) on the logit scale
## (for a fixed value also the ends that its distribution accepts, `closed` in
## the parameter's spec, both by default: a boundary is a model the density may
## accept) and finite on the identity scale. Names outside `spec` are left to the
## caller.
.distfit_check_parlist <- function(lst, what, spec, family) {
  fixed <- what == "fix.arg"
  if (!(is.list(lst) || is.numeric(lst)))
    stop("'", what, "' must be a named list of numbers", call. = FALSE)
  nms <- names(lst)
  if (length(lst) && (is.null(nms) || anyNA(nms) || any(!nzchar(nms))))
    stop("'", what, "' must be a named list: give each value the name of its ",
         "parameter (", paste(names(spec), collapse = ", "), ")", call. = FALSE)
  if (anyDuplicated(nms))
    stop("duplicated names in '", what, "': ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "), call. = FALSE)
  for (nm in intersect(nms, names(spec))) {
    v <- lst[[nm]]
    rho_w <- family == "gwaring" && nm == "rho"
    cl <- spec[[nm]]$closed %||% c(TRUE, TRUE)
    dom <- switch(spec[[nm]]$transform,
                  log = if (rho_w) "a single finite number above 1" else
                    if (fixed && isTRUE(spec[[nm]]$allow_zero)) "a single non-negative finite number" else
                      "a single positive finite number",
                  logit = if (fixed)
                    sprintf("a single number in %s0, 1%s", if (cl[1L]) "[" else "(",
                            if (cl[2L]) "]" else ")") else
                    "a single number strictly between 0 and 1",
                  "a single finite number")
    ok <- is.numeric(v) && length(v) == 1L && !is.na(v) && is.finite(v) &&
      switch(spec[[nm]]$transform,
             log = if (fixed && isTRUE(spec[[nm]]$allow_zero)) v >= 0 else
               v > if (rho_w) 1 else 0,
             logit = if (fixed) (v > 0 || (cl[1L] && v == 0)) &&
                                (v < 1 || (cl[2L] && v == 1)) else v > 0 && v < 1,
             TRUE)
    if (!ok)
      stop(sprintf("'%s' value for %s must be %s", what, nm, dom), call. = FALSE)
  }
  invisible(NULL)
}


#' Parameter specification for each type and family
#' @noRd
.distfit_parspec <- function(type, family) {
  ## PIG / ZIPIG are mixed-Poisson count models (not EGPDs): mu (mean, log),
  ## sigma (dispersion, log), and for ZIPIG pi (zero inflation, logit, in [0, 1)).
  ## The G-transformation `type` does not apply. `closed` gives the ends of a
  ## logit parameter's space at which it may be fixed (both by default).
  if (family %in% c("pig", "zipig")) {
    spec <- list(mu = list(transform = "log"), sigma = list(transform = "log"))
    if (family == "zipig")
      spec <- c(spec, list(pi = list(transform = "logit", closed = c(TRUE, FALSE))))
    return(spec)
  }

  ## GPIG / ZIGPIG (Zhu & Joe 2009), native (a, b, c) parameterisation: a in
  ## (0, 1] tail exponent (logit; a = 1 is Poisson), b > 0 level (log), c in
  ## [0, 1] down-weight (logit; c = 1 is discrete stable); ZIGPIG adds pi in
  ## [0, 1) (logit). The G-transformation `type` does not apply.
  if (family %in% c("gpig", "zigpig")) {
    spec <- list(a = list(transform = "logit", closed = c(FALSE, TRUE)),
                 b = list(transform = "log"),
                 c = list(transform = "logit"))
    if (family == "zigpig")
      spec <- c(spec, list(pi = list(transform = "logit", closed = c(TRUE, FALSE))))
    return(spec)
  }

  ## Generalized Waring, mean parameterisation (as in egpd()): mu = mean (log),
  ## k > 0 (log), rho > 1 (log). The tail is P(Y = y) ~ y^-(rho+1), so xi = 1/rho.
  ## The mean exists only for rho > 1; the likelihood returns +Inf below that, which
  ## keeps the optimiser out of that region without a bespoke transform. The
  ## G-transformation `type` does not apply.
  if (family == "gwaring") {
    return(list(mu  = list(transform = "log"),
                k   = list(transform = "log"),
                rho = list(transform = "log")))
  }

  ## Bell / ZIBell (Castellares et al. 2018), native theta parameterisation:
  ## theta > 0 (log); ZIBell adds pi in [0, 1) (logit). The G-transformation
  ## `type` does not apply.
  if (family %in% c("bell", "zibell")) {
    spec <- list(theta = list(transform = "log"))
    if (family == "zibell")
      spec <- c(spec, list(pi = list(transform = "logit", closed = c(TRUE, FALSE))))
    return(spec)
  }

  ## Base GPD parameters
  spec <- list(
    sigma = list(transform = "log"),
    xi    = list(transform = "identity")
  )

  ## G-transformation parameters by type
  g_params <- switch(as.character(type),
    "1" = list(kappa = list(transform = "log")),
    "2" = list(kappa = list(transform = "log", allow_zero = TRUE)),
    "3" = list(kappa = list(transform = "log")),
    "4" = list(delta = list(transform = "log")),
    "5" = list(delta = list(transform = "log"),
               kappa = list(transform = "log")),
    "6" = list(kappa = list(transform = "log"),
               delta = list(transform = "log"),
               prob  = list(transform = "logit"))
  )

  spec <- c(spec, g_params)

  ## Zero-inflation parameter
  if (family %in% c("ziegpd", "zidegpd")) {
    spec <- c(spec, list(pi = list(transform = "logit")))
  }

  ## Compound Poisson rate parameter
  if (family %in% c("cpegpd", "cpdegpd")) {
    spec <- c(spec, list(lambda = list(transform = "log")))
  }

  spec
}


#' Transform natural-scale parameters to unconstrained scale
#' @noRd
.par_to_theta <- function(par, spec) {
  theta <- numeric(length(par))
  names(theta) <- names(par)
  for (nm in names(par)) {
    theta[nm] <- switch(spec[[nm]]$transform,
      "log"      = log(par[nm]),
      "logit"    = qlogis(par[nm]),
      "identity" = par[nm]
    )
  }
  theta
}


#' Transform unconstrained parameters to natural scale
#' @noRd
.theta_to_par <- function(theta, spec) {
  if (is.null(names(theta))) names(theta) <- names(spec)
  par <- numeric(length(theta))
  names(par) <- names(theta)
  for (nm in names(theta)) {
    par[nm] <- switch(spec[[nm]]$transform,
      "log"      = exp(theta[nm]),
      "logit"    = plogis(theta[nm]),
      "identity" = theta[nm]
    )
  }
  par
}


#' Jacobian diagonal of back-transformation (for delta method)
#' @noRd
.distfit_jacobian <- function(theta, spec) {
  jac <- numeric(length(theta))
  names(jac) <- names(theta)
  for (nm in names(theta)) {
    jac[nm] <- switch(spec[[nm]]$transform,
      "log" = exp(theta[nm]),
      "logit" = {
        p <- plogis(theta[nm])
        p * (1 - p)
      },
      "identity" = 1
    )
  }
  jac
}


#' Automatic starting values
#' @noRd
.distfit_start <- function(x, type, family, spec, fixed = list()) {
  ## Remove non-finite values for moment estimation
  xclean <- x[is.finite(x)]
  if (length(xclean) < 2) xclean <- c(0.1, 0.2)

  ## PIG / ZIPIG: method-of-moments seeds (var = mu + sigma*mu^2).
  if (family %in% c("pig", "zipig")) {
    xm <- mean(xclean); v <- stats::var(xclean)
    if (!is.finite(xm) || xm <= 0) xm <- 1
    if (!is.finite(v)) v <- xm
    start <- list(
      mu    = max(xm, 0.1),
      sigma = min(max((v - xm) / max(xm^2, 1e-8), 0.1), 10)
    )
    if (family == "zipig")
      start$pi <- max(min(mean(x == 0), 0.99), 0.01)
    return(start[names(spec)])
  }

  ## GPIG / ZIGPIG: seed a = 1/2 (PIG case), c from the empirical dispersion
  ## index D via D-1 = (1-a)c/(1-c), and b from the sample mean.
  if (family %in% c("gpig", "zigpig")) {
    xm <- mean(xclean); v <- stats::var(xclean)
    if (!is.finite(xm) || xm <= 0) xm <- 1
    if (!is.finite(v)) v <- xm
    D <- if (xm > 0) max(v / xm, 1.05) else 2
    a0 <- 0.5
    cc <- min(max((D - 1) / ((1 - a0) + (D - 1)), 0.05), 0.95)
    b0 <- max(xm, 0.1) * (1 - cc)^(1 - a0) / (a0 * cc)
    start <- list(a = a0, b = max(b0, 0.1), c = cc)
    if (family == "zigpig")
      start$pi <- max(min(mean(x == 0), 0.99), 0.01)
    return(start[names(spec)])
  }

  ## Generalized Waring: seed the mean from the data, and rho from the index of
  ## dispersion -- a heavier tail shows up as greater overdispersion, so a large
  ## variance-to-mean ratio seeds a smaller rho.
  if (family == "gwaring") {
    xm <- mean(xclean); v <- stats::var(xclean)
    if (!is.finite(xm) || xm <= 0) xm <- 1
    if (!is.finite(v) || v <= xm) v <- 2 * xm
    rho0 <- 1 + 2 / max(log1p(v / xm), 0.2)
    return(list(mu  = max(xm, 0.1),
                k   = 1,
                rho = min(max(rho0, 1.2), 20))[names(spec)])
  }

  ## Bell / ZIBell: seed the native theta at the MLE theta.hat = W0(ybar).
  if (family %in% c("bell", "zibell")) {
    xm <- mean(xclean)
    if (!is.finite(xm) || xm <= 0) xm <- 1
    start <- list(theta = max(bell_W0_cpp(xm), 0.05))
    if (family == "zibell")
      start$pi <- max(min(mean(x == 0), 0.99), 0.01)
    return(start[names(spec)])
  }

  ## Moments of the continuous variable X behind the data. A discrete count is
  ## Y = floor(X), so X ~ Y + U with U uniform on [0, 1): the sample is Y + 1/2 and
  ## its variance gains 1/12 (Sheppard's correction). DEGPD zeros belong to the law
  ## (X < 1) and are kept. The zero-inflated families use their positive
  ## observations only, since part of their zeros is structural: for ZIEGPD these are
  ## draws of X itself; for ZIDEGPD they are X | X >= 1, whose excess X - 1 is GPD
  ## with the same shape and scale sigma + xi (threshold stability), so the excess
  ## Y - 1/2 is used and the scale shifted back below. The compound-Poisson families
  ## use the positive observations, adjusted to a per-event scale below.
  discrete <- family %in% c("degpd", "zidegpd", "cpdegpd")
  xpos <- if (family == "degpd") xclean else xclean[xclean > 0]
  if (discrete) xpos <- xpos + if (family == "zidegpd") -0.5 else 0.5
  if (length(xpos) < 2) xpos <- abs(xclean) + 0.01

  m <- mean(xpos)
  v <- var(xpos) + if (discrete) 1 / 12 else 0
  if (!is.finite(m) || !is.finite(v) || v <= 0) { m <- 1; v <- 1 }

  ## For cpegpd, the data is S = X_1 + ... + X_N where N ~ Poisson(lambda).
  ## E[S|S>0] = lambda*E[X]/(1-exp(-lambda)). Adjust moments to target
  ## the individual severity rather than the aggregate.
  if (family %in% c("cpegpd", "cpdegpd")) {
    p0 <- max(min(mean(xclean == 0), 0.99), 0.01)
    lambda_init <- -log(p0)
    ## Rough per-event mean: E[X] ~ E[S] / lambda = mean(xclean) / lambda
    per_event_mean <- mean(xclean) / lambda_init
    m <- max(per_event_mean, 0.01)
    v <- max(v / lambda_init, 0.01)  # Var(S) ~ lambda * E[X^2], rough scale
  }

  ## GPD moment estimates: with m = sigma / (1 - xi) and
  ## v = sigma^2 / ((1 - xi)^2 (1 - 2 xi)), m^2 / v = 1 - 2 xi, so
  ## xi = (1 - m^2 / v) / 2 and sigma = m (1 - xi). xi is kept in [-0.4, 0.45]
  ## (the variance needs xi < 1/2), which also keeps sigma = m (1 - xi) positive.
  xi_start <- min(max(0.5 * (1 - m^2 / v), -0.4), 0.45)
  sigma_start <- m * (1 - xi_start)
  if (family == "zidegpd") sigma_start <- sigma_start - xi_start   # excess over 1
  sigma_start <- max(sigma_start, 0.01)

  ## Safety: ensure sigma + xi * max(x) > 0
  xmax <- max(xclean)
  if (is.finite(xmax) && xmax > 0 &&
      is.finite(sigma_start + xi_start * xmax) &&
      sigma_start + xi_start * xmax <= 0) {
    xi_start <- max(xi_start, -sigma_start / (2 * xmax))
  }

  ## Grid-search for delta and kappa starting values
  delta_start <- fixed$delta %||% 1
  kappa_start <- fixed$kappa %||% 1
  prob_start <- fixed$prob %||% 0.5

  if ("delta" %in% names(spec) || "kappa" %in% names(spec)) {
    ## Helper to evaluate NLL at candidate shape parameters; the zero-inflated
    ## families are scored with their zero-inflation at its start value
    pi0 <- fixed$pi %||% max(min(mean(x == 0), 0.99), 0.01)
    .eval_nll <- function(sigma, xi, kappa, delta, type, family, x, prob = prob_start) {
      sigma <- fixed$sigma %||% sigma
      xi <- fixed$xi %||% xi
      tryCatch({
        if (family == "zidegpd") {
          -sum(dzidiscegpd(x, pi = pi0, prob = prob, kappa = kappa, delta = delta,
                           sigma = sigma, xi = xi, type = type, log = TRUE))
        } else if (family == "ziegpd") {
          -sum(dziegpd(x, pi = pi0, prob = prob, kappa = kappa, delta = delta,
                       sigma = sigma, xi = xi, type = type, log = TRUE))
        } else if (family %in% c("degpd", "cpdegpd")) {
          -sum(ddiscegpd(x, prob = prob, kappa = kappa, delta = delta,
                         sigma = sigma, xi = xi, type = type, log = TRUE))
        } else {
          ll <- degpd_density(x, prob = prob, kappa = kappa, delta = delta,
                              sigma = sigma, xi = xi, type = type, log = TRUE)
          -sum(ll)
        }
      }, error = function(e) Inf)
    }

    if (type == 6) {
      ## Equal powers have identical scores at weight 1/2: a gradient
      ## search then stays in the one-power submodel even at a saddle.
      seeds <- expand.grid(k = fixed$kappa %||% c(.3, .7, 1.5),
                           d = fixed$delta %||% c(1, 3, 7),
                           p = fixed$prob %||% c(.3, .7))
      seeds <- seeds[seeds$k != seeds$d, , drop = FALSE]
      if (nrow(seeds)) {
        z <- if (family == "cpegpd") x[x > 0] else x
        scores <- vapply(seq_len(nrow(seeds)), function(i)
          .eval_nll(sigma_start, xi_start, seeds$k[i], seeds$d[i], type, family, z, seeds$p[i]),
          numeric(1))
        j <- if (any(is.finite(scores))) which.min(replace(scores, !is.finite(scores), Inf)) else 1L
        kappa_start <- seeds$k[j]; delta_start <- seeds$d[j]; prob_start <- seeds$p[j]
      }
    } else if ("delta" %in% names(spec)) {
      delta_candidates <- c(0.3, 0.5, 1.0, 1.5, 2.5, 4.0)
      kappa_for_grid <- kappa_start
      delta_nlls <- vapply(delta_candidates, function(d) {
        .eval_nll(sigma_start, xi_start, kappa_for_grid, d, type, family, x)
      }, numeric(1))
      best_delta <- which.min(delta_nlls)
      if (length(best_delta) > 0 && is.finite(delta_nlls[best_delta])) {
        delta_start <- delta_candidates[best_delta]
      }
    }

    if (type != 6 && "kappa" %in% names(spec)) {
      kappa_candidates <- c(0.5, 1.0, 1.5, 2.5)
      delta_for_grid <- delta_start
      kappa_nlls <- vapply(kappa_candidates, function(k) {
        .eval_nll(sigma_start, xi_start, k, delta_for_grid, type, family, x)
      }, numeric(1))
      best_kappa <- which.min(kappa_nlls)
      if (length(best_kappa) > 0 && is.finite(kappa_nlls[best_kappa])) {
        kappa_start <- kappa_candidates[best_kappa]
      }
    }
  }

  start <- list()
  for (nm in names(spec)) {
    start[[nm]] <- switch(nm,
      "sigma" = sigma_start,
      "xi"    = xi_start,
      "kappa" = kappa_start,
      "delta" = delta_start,
      "prob"  = prob_start,
      "pi"    = max(min(mean(x == 0), 0.99), 0.01),
      "lambda" = -log(max(min(mean(x == 0), 0.99), 0.01)),
      1 # fallback
    )
  }

  start
}


#' Stop if the optimiser never left the infeasible region, warn if it did not
#' converge. .distfit_nll() returns 1e20 for parameter values with zero likelihood
#' (so that optim() can continue), so an optimum at 1e20 is not a fit.
#' @noRd
.distfit_check_optim <- function(opt, method) {
  if (!is.finite(opt$value) || opt$value >= 1e20)
    stop("no parameter value with a positive likelihood was found: the data lie ",
         "outside the support of every distribution visited from the starting ",
         "values; supply `start` (or check `fix.arg`).", call. = FALSE)
  if (!identical(as.integer(opt$convergence), 0L)) {
    msg <- sprintf(paste0("optim (%s) did not converge (code %d%s); the estimates ",
                          "may not be the maximum likelihood estimates. Increase ",
                          "control$maxit or try another optim.method."),
                   method, as.integer(opt$convergence),
                   if (identical(as.integer(opt$convergence), 1L))
                     ": iteration limit reached" else
                     if (!is.null(opt$message)) paste0(": ", opt$message) else "")
    warning(structure(list(message = msg, call = NULL),
                      class = c("distfit_nonconvergence", "warning", "condition")))
  }
  invisible(opt)
}


#' Negative log-likelihood for fitegpd
#' @noRd
.distfit_nll <- function(theta, x, type, family, fix.arg, spec, h = NULL, soft = FALSE) {
  ## Back-transform to natural scale
  par <- .theta_to_par(theta, spec)

  ## Merge with fixed args
  all_par <- as.list(c(par, unlist(fix.arg)))

  ## Generalized Waring. dgwaring() rejects rho <= 1 outright, so the tryCatch is
  ## what enforces the constraint: an out-of-range proposal scores 1e20 rather than
  ## stopping the optimisation.
  if (family == "gwaring") {
    nll <- tryCatch(
      -sum(dgwaring(x, mu = all_par[["mu"]], k = all_par[["k"]],
                    rho = all_par[["rho"]], log = TRUE)),
      error = function(e) Inf)
    if (!is.finite(nll)) return(1e20)
    return(nll)
  }

  ## PIG / ZIPIG mixed-Poisson count models, by the package's native recursion,
  ## which is exact down to the Poisson limit sigma -> 0 (gamlss.dist::dPIG breaks
  ## down there, giving pmf values far above 1). Parameter values that the native
  ## functions reject (with NaN and a warning) are infeasible here.
  if (family %in% c("pig", "zipig")) {
    mu <- all_par[["mu"]]
    sigma <- all_par[["sigma"]]
    pi_val <- if (family == "zipig") all_par[["pi"]] else 0
    if (!(is.finite(mu) && mu > 0 && is.finite(sigma) && sigma > 0 &&
          is.finite(pi_val) && pi_val >= 0 && pi_val < 1))
      return(1e20)
    nll <- tryCatch({
      if (family == "pig") {
        -sum(dpig(x, mu = mu, sigma = sigma, log = TRUE))
      } else {
        -sum(dzipig(x, mu = mu, sigma = sigma, pi = pi_val, log = TRUE))
      }
    }, error = function(e) Inf)
    if (!is.finite(nll)) return(1e20)
    return(nll)
  }

  ## GPIG / ZIGPIG (native a, b, c) via the package recursion.
  if (family %in% c("gpig", "zigpig")) {
    ## When a=1 with both factors free, estimate the Poisson component mean
    ## directly; c=1 is only its native evaluation representative.
    b <- all_par[["mu"]] %||% all_par[["b"]]
    nll <- tryCatch({
      if (family == "gpig") {
        -sum(dgpig(x, a = all_par[["a"]], b = b, c = all_par[["c"]],
                   log = TRUE))
      } else {
        -sum(dzigpig(x, a = all_par[["a"]], b = b, c = all_par[["c"]],
                     pi = all_par[["pi"]], log = TRUE))
      }
    }, error = function(e) Inf)
    if (!is.finite(nll)) return(1e20)
    return(nll)
  }

  ## Bell / ZIBell (native theta) via the closed-form pmf.
  if (family %in% c("bell", "zibell")) {
    nll <- tryCatch({
      if (family == "bell") {
        -sum(dbell(x, theta = all_par[["theta"]], log = TRUE))
      } else {
        -sum(dzibell(x, theta = all_par[["theta"]], pi = all_par[["pi"]],
                     log = TRUE))
      }
    }, error = function(e) Inf)
    if (!is.finite(nll)) return(1e20)
    return(nll)
  }

  sigma  <- all_par[["sigma"]]
  xi     <- all_par[["xi"]]
  kappa  <- if ("kappa" %in% names(all_par)) all_par[["kappa"]] else NA
  delta  <- if ("delta" %in% names(all_par)) all_par[["delta"]] else NA
  prob   <- if ("prob" %in% names(all_par)) all_par[["prob"]] else NA
  pi_val <- if ("pi" %in% names(all_par)) all_par[["pi"]] else NA
  lambda <- if ("lambda" %in% names(all_par)) all_par[["lambda"]] else NA

  ## Per-observation log-likelihoods. The pmf is not floored: an observation
  ## outside the support (e.g. above the endpoint when xi < 0) has log-likelihood
  ## -Inf, not log(DBL_MIN) = -708.
  ll <- tryCatch({
    switch(family,
      "egpd" = degpd_density(x, prob = prob, kappa = kappa, delta = delta,
                             sigma = sigma, xi = xi, type = type, log = TRUE),
      "degpd" = ddiscegpd(x, prob = prob, kappa = kappa, delta = delta,
                          sigma = sigma, xi = xi, type = type, log = TRUE),
      "ziegpd" = dziegpd(x, pi = pi_val, prob = prob, kappa = kappa, delta = delta,
                         sigma = sigma, xi = xi, type = type, log = TRUE),
      "zidegpd" = dzidiscegpd(x, pi = pi_val, prob = prob, kappa = kappa,
                              delta = delta, sigma = sigma, xi = xi, type = type,
                              log = TRUE),
      "cpegpd" = {
        l <- dcpegpd(x, lambda = lambda, prob = prob, kappa = kappa,
                     delta = delta, sigma = sigma, xi = xi, type = type,
                     h = h, log = TRUE)
        pos <- x > 0
        if (any(pos)) l[pos] <- l[pos] - log(h)
        l
      },
      "cpdegpd" = dcpdegpd(x, lambda = lambda, prob = prob, kappa = kappa,
                           delta = delta, sigma = sigma, xi = xi, type = type,
                           log = TRUE)
    )
  }, error = function(e) NA_real_)

  if (length(ll) != length(x) || !is.numeric(ll)) return(1e20)
  ## `soft`: only to steer an infeasible start into the support (see distfit());
  ## each offending observation costs -log(DBL_MIN) instead of +Inf
  if (soft) ll <- pmax(ifelse(is.na(ll), -Inf, ll), log(.Machine$double.xmin))
  if (any(!is.finite(ll))) return(1e20)
  -sum(ll)
}
