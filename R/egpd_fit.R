## Main egpd fitting function and S3 methods
## Adapted from evgam by Ben Youngman

#' Fit Extended Generalized Pareto Distribution GAMs
#'
#' @param formula a formula or list of formulae, one per distribution parameter. The
#'   first (or the only) formula carries the response; in a named list (for a
#'   family other than \code{"custom"}) any one element may carry it. A list may be named, in any
#'   order, by each parameter's natural name (e.g. \code{sigma}, \code{xi},
#'   \code{kappa}, \code{pi}, \code{lambda}, \code{sdlog}), by its link-prefixed
#'   name (e.g. \code{lsigma}, \code{logscale}, \code{logitpi}), or, for the scale and
#'   shape of the EGPD families, by \code{scale}/\code{psi} and \code{shape}; unnamed
#'   elements fill the remaining slots in order. Unknown names are an error listing
#'   the valid ones. The legacy spellings \code{lpi} (zero-inflation probability) and
#'   \code{lxi} (shape under an identity link) still work but are deprecated. Each
#'   formula is evaluated in its own environment, so smooth arguments may refer to
#'   local variables. As in \pkg{mgcv}, parametric terms, the covariates of smooths
#'   (including \code{by} variables and the factor of an \code{"fs"} smooth) and the
#'   response may be expressions of the columns of \code{data}, such as
#'   \code{log(x)}, \code{factor(g)}, \code{I(x^2)}, \code{poly(x, 2)},
#'   \code{s(log(x))}, \code{s(x, by = factor(g))} or \code{I(y + 1) ~ x}; the
#'   columns they use must be in \code{data} (and in the \code{newdata} of
#'   \code{predict()}). Other symbols in them are looked up in the formula's
#'   environment when fitting, as by \code{mgcv::gam()}; the response is evaluated
#'   once, in the environment of the formula that carries it. As in base \R, an
#'   expression is always evaluated, even if \code{data} has a column named like it
#'   (\code{"log(x)"}); only a name, backquoted if need be, refers to a column.
#'   In a term such as \code{I(x * cfg$gain)}, \code{gain} names a component
#'   of \code{cfg}, not an additional data column.
#'   Censored responses, \code{cens(left, right)} on the left-hand
#'   side (an \pkg{evgam} feature), are not supported by any family and give an
#'   error.
#' @param data a data frame. Incomplete rows, including missing values
#'   introduced by evaluated predictor terms or offsets, are omitted jointly
#'   across parameters before fitting.
#' @param family a character string naming the distribution family, including
#'   \code{"egpd"}, \code{"degpd"}, \code{"zidegpd"}, \code{"pig"}, \code{"zipig"},
#'   \code{"gpig"}, \code{"gpignat"}, \code{"zigpig"}, \code{"zigpignat"},
#'   \code{"bell"}, \code{"bellnat"}, \code{"zibell"}, \code{"zibellnat"},
#'   \code{"cmp"}, \code{"prl"}, \code{"gpois"}, \code{"gwaring"},
#'   \code{"ztdegpd"}, \code{"ztgwaring"}, \code{"plnorm"}, \code{"comppareto"},
#'   or \code{"custom"}.
#'   Continuous zero-inflated EGPD utilities are provided elsewhere in the
#'   package, but are not fitted via \code{egpd()}.
#' @param correctV logical: should variance-covariance matrix account for smoothing parameter uncertainty? Defaults to TRUE
#' @param rho0 initial log smoothing parameters for the smoothing-parameter search:
#'   one per independent parameter not fixed in a term, or one value for all
#'   (with \code{gamma != 1}, of \code{log(sp / gamma)}; see \code{gamma}).
#' @param inits initial parameter values: one value per distribution parameter (the
#'   intercept of its linear predictor), or a list of candidate values per parameter,
#'   of which the best combination is used. Explicit starts are coefficients:
#'   offsets are added to them. Automatic predictor-level starts instead subtract
#'   the mean offset from a free intercept (the median if the mean-centred
#'   start has non-finite likelihood), and never move a fixed predictor.
#'   Starting values at which the likelihood
#'   cannot be evaluated are an error, and so, for the count families whose
#'   likelihoods clamp their parameters (\code{"gpois"}, \code{"gwaring"},
#'   \code{"ztgwaring"}, \code{"plnorm"}, \code{"prl"}, \code{"cmp"}), are values
#'   outside the clamped range, where the fit could not move them: for example
#'   \code{inits = c(0, 0, 0)} for \code{"gwaring"} puts \eqn{\rho} at 1, outside
#'   its parameter space \eqn{\rho > 1}. Without \code{inits}, \code{"ztgwaring"}
#'   starts from the mean, \eqn{k = 1} and \eqn{\rho = 2.5}, unless the intercept-only
#'   fit from there leaves the interior of the parameter space (e.g. runs into
#'   \eqn{\rho = 1}) while the intercept-only maximum likelihood estimate of the
#'   zero-truncated law lies inside it; it then starts from that estimate, as the first
#'   start of \code{\link{fit_hurdle}} does.
#' @param outer outer optimization method: "bfgs" (default), "newton", "fd", or "fixed"
#'   Continuous EGPD models 1 and 3, and DEGPD model 1 with an unbounded log or
#'   identity shape link, supply analytic third derivatives for the REML gradient.
#'   Other models, compact fits (of any family) and a finite \code{xi.max} use
#'   finite differences;
#'   near-zero shape, non-finite analytic derivatives, and DEGPD cells too narrow for
#'   the analytic derivatives (very large counts or scales, where
#'   \eqn{(\sigma + \xi y)/|\xi|} exceeds about \eqn{4.5 \times 10^{10}}) trigger a
#'   warning of class \code{"egpd_d34_fallback"} and a finite-difference retry.
#' @param control a named list with optional \code{outer} and \code{inner} lists.
#'   Each overrides the corresponding optimizer defaults: \code{steptol},
#'   \code{itlim}, \code{fntol}, \code{gradtol}, and \code{stepmax}. The
#'   finite-difference REML derivatives tighten the default inner \code{gradtol}
#'   for their own refits; an inner \code{gradtol} given here is used as is.
#' @param removeData logical: should data be removed from the returned object? Defaults to FALSE
#' @param trace an integer controlling output verbosity
#' @param knots a list of knot values for smooth terms
#' @param maxdata maximum number of data rows
#' @param maxspline maximum number of rows for spline basis construction: with more
#'   rows than this, each smooth basis is set up on a random subsample of
#'   \code{maxspline} rows and then evaluated on all of them. The subsample is
#'   extended by one row of each level of a factor in the formula that it misses, so
#'   factor terms are set up as on the full data.
#' @param compact logical: keep one design row per distinct combination of
#'   evaluated design values and offsets? Terms are evaluated on the fitting
#'   sample before grouping, so compaction preserves their multiplicities
#'   in sample-dependent transformations. Defaults to \code{FALSE}. For custom likelihoods see
#'   the observation map under \code{custom.fns}.
#' @param egpd.args a list of arguments for EGPD family (e.g., m=1)
#' @param degpd.args a list of arguments for the DEGPD family. \code{m} selects the
#'   G-function (model 1-6, default 1). \code{xi.max} (default \code{Inf}) caps the
#'   tail index \eqn{\xi} via a bounded shape link
#'   \eqn{\xi = \texttt{xi.max}\,(1 - \exp(-\exp(\eta)/\texttt{xi.max}))}, which reduces
#'   to the usual unbounded \eqn{\xi = \exp(\eta)} log link as \code{xi.max} tends to
#'   \code{Inf} and saturates at \code{xi.max}; this stabilises fits on short, spiky,
#'   or heavily over-dispersed series. Applies to all DEGPD models 1-6.
#'   Must be a positive scalar; a finite cap cannot be combined with
#'   \code{link = "identity"}. Default and restart shape starts are placed inside
#'   the cap (at most \code{xi.max / 2}), since near the cap, and as \eqn{\xi \to 0},
#'   the likelihood barely depends on \eqn{\eta} and a fit can stop there short of
#'   the optimum (see \code{restarts}).
#' @param zidegpd.args a list of arguments for ZIDEGPD family (e.g., m=1)
#' @param gpig.args a list of arguments controlling how the GPIG/ZIGPIG pmf is
#'   evaluated, passed to \code{\link{gpig_control}}: \code{method} (one of
#'   \code{"hybrid"}, \code{"recursion"}, \code{"trunc"}, \code{"saddlepoint"},
#'   \code{"legacy"}), \code{yswitch}, \code{order}, \code{eps}, \code{sptol},
#'   \code{ymax}, \code{nquad} and \code{normalize}. Applied for the duration of the
#'   fit and then restored. The effective settings (including inherited
#'   settings) are saved as \code{gpig.control} for profile likelihood
#'   intervals. Ignored, with a warning, for other families.
#' @param comppareto.args a list of arguments for the CompPareto GAM family.
#'   Currently this supports \code{spec = "lnorm"}, \code{"gamma"},
#'   \code{"weibull"}, or \code{"exp"} for the body distribution. See
#'   \code{restarts} for splice-point searches and convergence at a splice knot.
#' @param sandwich.args Optional list with \code{id} naming the cluster column
#'   (use a distinct ID per independent observation). Let \eqn{H} be the observed
#'   information at the unpenalised fit and \eqn{J} the cross-product of the
#'   cluster score sums. On the \eqn{p} free coefficients, \code{method="magnitude"}
#'   (default) multiplies the likelihood on the log scale by
#'   \eqn{p / {\rm tr}(H^{-1}J)}. \code{method="curvature"} transforms the
#'   coefficients about the unpenalised estimate so the inverse adjusted curvature
#'   there is \eqn{H^{-1}JH^{-1}}. Smooths are fitted to the adjusted likelihood
#'   with the usual penalty; \code{Vp} inverts the adjusted penalised curvature
#'   at the resulting fit and \code{Vc} includes smoothing-parameter uncertainty.
#'   For BHHH families the unadjusted fit still uses score products, but sandwich
#'   inference differentiates the likelihood scores. Numerically unresolved or
#'   non-positive definite observed information (or a singular \eqn{J} for the
#'   curvature adjustment) is an error, not a fallback to BHHH. For every family
#'   the magnitude adjustment warns if \eqn{J} lacks numerical rank on the free
#'   coefficients, for example when cluster fixed effects make some cluster score
#'   sums zero. Its scalar factor is unchanged and need not give a reliable
#'   standard error for an individual coefficient. The curvature adjustment
#'   requires full-rank positive definite \eqn{H} and \eqn{J} and a numerically
#'   resolved adjusted information matrix at the fitted coefficients. For exact-Hessian
#'   families only, \code{force=TRUE} allows pseudo-inverse adjustments; a forced
#'   singular curvature adjustment does not support inference in the missing
#'   directions.
#' @param custom.fns for \code{family = "custom"}, a list defining the likelihood (as
#'   for \pkg{evgam}'s custom families), with elements
#'   \describe{
#'     \item{\code{d0}}{\code{function(pars, likdata)}: the negative log-likelihood at
#'       the coefficient vector \code{pars}, summed over the rows of
#'       \code{likdata$y} using the linear predictors described below.}
#'     \item{\code{d120}}{\code{function(pars, likdata)}: a matrix with one row per
#'       observation: the first derivatives of its negative log-likelihood with
#'       respect to each linear predictor, then the second derivatives, upper
#'       triangle by row (for two predictors: 1-1, 1-2, 2-2).}
#'     \item{\code{d340}}{optional third and fourth derivatives, used for
#'       smoothing-parameter selection; without them, finite differences are used.}
#'     \item{\code{unlink}}{optional list of inverse-link functions, one per linear
#'       predictor (in formula order), giving the distributional parameter, as used
#'       by \code{fixed}, \code{predict(type = "response")}, \code{simulate()} and
#'       \code{confint()}; each must be monotone (increasing or decreasing). A
#'       missing or \code{NULL} entry uses the link named by the
#'       predictor's prefix (\code{log}, \code{logit}, \code{probit},
#'       \code{cloglog}), and otherwise the identity. An attribute \code{"deriv"},
#'       the derivative of the inverse link, gives the delta-method standard errors
#'       of \code{predict(se.fit = TRUE)}; without it a central difference is used.
#'       For a response-scale pin, an optional function attribute \code{"link"}
#'       maps the fixed value back to a finite linear predictor on the intended
#'       domain. Otherwise a standard forward link or a numerical root bracketed
#'       on a finite grid is used. If a valid inverse cannot be found and checked
#'       against \code{unlink}, fitting stops rather than guessing from the name.}
#'     \item{\code{q}}{optional quantile function, \code{function(p, par1, ...)},
#'       taking the distributional parameters in formula order (at most four); needed for
#'       quantile predictions and \code{simulate(probs = )}.}
#'   }
#'   In each callback, construct predictor \code{i} in the row order of
#'   \code{likdata$y}. Select its coefficients with \code{likdata$idpars},
#'   expand the design product if compact, \emph{then} add the offset:
#'   \preformatted{
#'   eta <- drop(likdata$X[[i]] \%*\% pars[likdata$idpars == i])
#'   if (likdata$duplicate == 1L) eta <- eta[likdata$dupid + 1L]
#'   off <- likdata$offsets[[i]]
#'   if (length(off)) eta <- eta + off
#'   }
#'   The observation-to-design map \code{dupid} is zero-based. The offsets are
#'   already in observation order: do not expand them again. They include formula
#'   offsets and the constant link-scale offsets used by \code{fixed}; for a pin
#'   the design product is zero. Apply the family's inverse link to this complete
#'   predictor when evaluating the likelihood, and return derivative rows in the
#'   same observation order. The design matrices passed to existing callbacks
#'   are not automatically expanded.
#' @param sp optional non-negative smoothing parameters on the reported scale
#'   (see \code{gamma}). Supplying them skips outer optimization. In predictor
#'   order, give one value per independent parameter not already fixed in a
#'   smooth, or a scalar to use for all of them. Alternatively supply one
#'   multiplier per penalty, consistent with all fixed and shared values.
#'   A non-negative term-level \code{sp} fixes that penalty as in \pkg{mgcv}; an \code{id}
#'   shares smoothing parameters only within its predictor (including the
#'   levels of a factor \code{by} smooth). For linked terms specify fixed
#'   values on the first term with that id. A top-level zero removes that
#'   penalty exactly; for a zero in a term, \pkg{mgcv}'s setup supplies a
#'   small positive effective zero, retained in the fit and its REML score.
#'   Custom smooths whose log-penalty map needs a product or power of free
#'   parameters rather than one multiplier per matrix are not supported.
#' @param gamma a positive number; as in \code{mgcv::gam()}, values above 1 make the
#'   REML smoothing-parameter selection favour smoother fits, by selecting the
#'   smoothing parameters with the log-likelihood divided by \code{gamma}. Only the
#'   selection is affected: the coefficients are those that maximise the
#'   log-likelihood minus the penalty at the selected smoothing parameters, and their
#'   covariance matrices \code{Vp} and \code{Vc} (and so \code{summary()},
#'   \code{predict(se.fit = TRUE)}, \code{simulate()} and Wald intervals) are those of
#'   this unscaled penalised log-likelihood, so \code{gamma} does not inflate standard
#'   errors. The reported smoothing parameters \code{sp} are those of the penalty on
#'   the unscaled log-likelihood (\code{gamma} times those of the penalty on the
#'   log-likelihood divided by \code{gamma}), as in \pkg{mgcv}, and a supplied
#'   \code{sp} is read on that scale, so supplying the reported \code{sp} reproduces a
#'   fit whatever \code{gamma} is. Without smooth terms, or with \code{sp} supplied,
#'   \code{gamma} does not change the estimates. The reported log-likelihood is not
#'   divided by \code{gamma} (\code{logLik.adjusted}, the objective used for
#'   estimation, is). Cannot be combined with \code{sandwich.args}.
#' @param fixed a named list of response-scale parameter values to hold fixed
#'   during estimation (for example \code{list(kappa = 1)}). Names are matched as
#'   in a named \code{formula} list, and values are on the response scale
#'   whichever spelling is used (\code{list(lkappa = 1)} also means kappa = 1).
#' @param restarts logical: when \code{inits = NULL}, fit \code{family = "egpd"},
#'   \code{"degpd"} or \code{"zidegpd"} (and \code{"comppareto"}, see below) from several
#'   starting values -- the package
#'   default (for ZIDEGPD a warm start from the nested DEGPD fit) and four more, with
#'   shapes 0.05 and 0.3 combined with two fills of the remaining parameters (0 and
#'   \eqn{\log 2} on their link scales; for the continuous mixture carrier,
#'   \code{egpd.args = list(m = 2)}, weight 0.12 on \eqn{\kappa_1} with
#'   \eqn{\kappa_1 = 1, \kappa_2 = 2}, and weight 0.88 with \eqn{\kappa_1 = 1,
#'   \kappa_2 = 21}), plus the model 2 starts below -- and return the best of the
#'   converged fits? Defaults to \code{TRUE}; \code{FALSE} fits the default start
#'   only. The fits are compared by
#'   the criterion they optimise: without smooth terms, the log-likelihood; with smooth
#'   terms, the REML criterion (the lowest \code{negREML}) when the smoothing
#'   parameters are estimated, so that a wigglier fit with a higher log-likelihood but
#'   a worse REML criterion is not preferred, and the penalised log-likelihood (the
#'   log-likelihood less the penalty) when they are all fixed. Fits within
#'   0.01 of the best by that criterion count as the same optimum, and of those the one
#'   with the highest log-likelihood is returned -- for the continuous EGPD the default
#'   start's fit, so that another start replaces it only if it improves the criterion
#'   by more than 0.01, and a fit whose default start reaches the optimum is unchanged.
#'   The REML criterion is a Laplace approximation, which needs a regular maximum of
#'   the penalised likelihood: if any of the fits has a linear predictor other than the
#'   scale's on a boundary of its link (\eqn{|\eta| > 9}, as below), its shape on a
#'   flat end of the bounded link, or a penalised Hessian that is not positive
#'   definite, the likelihood is flat there and the criterion falls without bound as
#'   that fit drifts along it, so the fits are compared by the log-likelihood instead.
#'   With a sandwich adjustment (\code{sandwich.args}) the fits are always compared by
#'   the log-likelihood: the adjustment is estimated from each start's own unpenalised
#'   fit, so the adjusted criteria of different starts are on different scales.
#'   Starts that differ only for parameters held by
#'   \code{fixed} are fitted once, and every start fits the same data, including the
#'   same \code{maxdata} and \code{maxspline} subsamples.
#'
#'   Model 2 (EGPD, DEGPD and ZIDEGPD) also tries up to three starts that place
#'   the two components at different response quantiles. It screens trial scales,
#'   shapes and mixing weights on at most 512 response ranks; only the best start
#'   from each of three weight ranges is used for a full fit. Screening uses
#'   densities for EGPD and bin probabilities (including zero inflation if fitted)
#'   for counts. This helps escape both nearly equal powers and inferior separated
#'   mixtures; it does not guarantee a global maximum. In continuous model 2
#'   with negative shape, the scale can place the upper support endpoint just
#'   above the largest observation while the larger power grows. The likelihood
#'   can then be unbounded. The starts do not target such spikes.
#'
#'   For EGPD, DEGPD and ZIDEGPD carriers \code{m = 2}, \code{3} and \code{4},
#'   all the starts are tried, even for an interior primary fit.
#'   For the other ZIDEGPD carriers the four extra starts are fitted
#'   only when the default fit fails: it does not converge, its log-likelihood is not
#'   finite, or a fitted linear predictor other than the scale's (and other than
#'   those held by \code{fixed}) reaches a boundary, \eqn{|\eta| > 9} (for
#'   example \eqn{\pi}, \eqn{\kappa} or \eqn{\xi} below about \eqn{10^{-4}}).
#'   Each extra start retained is a full fit. On two intercept-only model 2
#'   examples with 400 observations, eight starts took about 8 to 11 times as
#'   long as a single start. Without the model 2 screen, five starts together
#'   take about twice as long as the default fit alone for intercept-only and
#'   parametric DEGPD models, and about 2.5 to 8 times as long for the continuous
#'   EGPD (about 3 times for a model 1 fit with one smooth). With a
#'   finite \code{degpd.args$xi.max}, if the best of them has its shape on either
#'   flat end of the bounded link (\eqn{\xi} within about \eqn{10^{-4}} relative of the
#'   cap, or below about \eqn{10^{-4}}), a fit with the unbounded log link is mapped
#'   onto the bounded link and tried as one more start (the log-link fit, another
#'   model, only supplies that start and is not one of the fits compared).
#'
#'   For \code{family = "comppareto"} a local maximum can occur at an observation
#'   (a kink in the splice point \eqn{\theta}) or between observations. For constant
#'   predictors a deterministic screen considers up to 128 distinct sample splits and
#'   the adjacent log-midpoints. It profiles a free \eqn{\alpha} analytically and fits
#'   the free body parameters at a coarse set of splits. With covariates, smooths or
#'   varying offsets it instead tries shifts of the fitted log-splice intercept (up
#'   to 32 splits and their midpoints), holding the other fitted predictors intact
#'   for scoring. This searches only a line even for a smooth on \eqn{\theta}; a
#'   pinned splice or a design without a splice intercept skips the screen.
#'   At most two peaks scoring more than 0.01 above the primary log-likelihood
#'   supply intercept starts for full refits. All coefficients and adjustments are
#'   fitted normally from each start. A primary fit with all observations on one
#'   side still tries the four quantile starts (50\%, 75\%, 90\%, 95\%). The
#'   primary fit is kept within 0.01 by the comparison criterion above; a failed
#'   primary is replaced only for a material improvement over the log-likelihood
#'   at its coefficients as well. This bounded search does not guarantee a global
#'   maximum. A selected split with no tail observations and a free unpenalised
#'   tail-shape intercept has no finite tail-shape optimum there.
#'
#'   For an exponential body with one log-splice intercept (possibly offset), the
#'   coefficient optimizer can stop at an observed splice with a nonzero ordinary
#'   gradient. With no smoothing parameters to estimate, the fitter can polish
#'   the smooth free coordinates at the knot and verify positive tangent curvature
#'   and both one-sided slopes of the objective actually optimized. This is a
#'   local convergence check, not a restart or a regularity guarantee for standard
#'   errors. A rate driven to zero is distinguished from a finite kink optimum.
#'   This check does not certify varying splice coefficients, curvature-sandwich
#'   transformations, or the outer REML search.
#'
#' @return An object of class \code{egpd}. Among its components, \code{logLik} is
#'   the log-likelihood of the data at the fitted coefficients (not divided by
#'   \code{gamma} or a sandwich adjustment factor; the adjusted value used for
#'   estimation is \code{logLik.adjusted}), \code{negREML} the REML criterion at the
#'   returned smoothing parameters (or, without smooths, the minimised negative
#'   log-likelihood), and \code{convergence} 0 on success, 1 if an iteration limit
#'   was reached, 2 if the log-likelihood cannot be evaluated at the returned
#'   coefficients (or is positive for a discrete family), and 3 if the coefficient
#'   search stopped with a non-negligible gradient, with the reason in
#'   \code{convergence.report}. A fit that did not converge has \code{logLik}
#'   \code{NA} and \code{AIC = BIC = 1e20}, and \code{egpd()} warns (condition class
#'   \code{"egpd_nonconvergence"}).
#'
#'   For penalised models, \code{sp} lists the independent parameters not
#'   already fixed by the formula (including values supplied at top level).
#'   \code{full.sp} gives the actual multiplier of each penalty matrix,
#'   including term-fixed and repeated values. Each predictor reports its own
#'   slice of both vectors. The correction from \code{Vp} to \code{Vc}
#'   propagates the covariance of the estimated log smoothing parameters
#'   through the first derivatives of the fitted coefficients.
#'
#'   \code{negREML} is the Laplace approximation to the negative log restricted
#'   likelihood (Wood, 2011, \emph{JRSS B} 73, 3--36). With \code{gamma != 1} or
#'   a sandwich adjustment it uses the adjusted log-likelihood. When comparing
#'   \code{negREML} to \code{gcv.ubre} from \code{mgcv::gam(method = "REML")},
#'   match the likelihood, scale and smoothing parameters. Even then the absolute
#'   scores can differ by an additive constant under \code{gamma != 1}, depending
#'   on the \pkg{mgcv} fitter. For example, with \eqn{M_p} unpenalised coefficients,
#'   the fixed-scale Gaussian REML score exceeds \code{negREML} by
#'   \eqn{M_p \log\gamma}, whereas a Gaussian location/scale fit
#'   (\code{mgcv::gaulss()}, \code{gam.fit5} in \pkg{mgcv} 1.9-4) gives
#'   \eqn{(M_p + 1)\log\gamma / 2}. These constants vanish at \code{gamma = 1}
#'   and do not affect smoothing-parameter selection within a model at fixed
#'   \code{gamma}.
#'
#'   \code{df} is the total effective degrees of freedom (also the \code{"df"}
#'   attribute of \code{logLik()}), from which \code{AIC = -2 * logLik + 2 * df} and
#'   \code{BIC = -2 * logLik + log(nobs) * df} are computed, and each parameter's
#'   \code{edf} holds the share of its coefficients. They are the diagonal of
#'   \eqn{(H_0 + S)^{-1} H_0}, where \eqn{H_0} is the Hessian of the negative
#'   log-likelihood at the fitted coefficients (of the adjusted one, with a sandwich
#'   adjustment) and \eqn{S} the penalty (without smooth terms \eqn{S = 0}, and each
#'   identified coefficient counts one). That formula assumes \eqn{H_0} is positive
#'   semi-definite, which it need not be at a penalised optimum -- e.g. with smooth
#'   terms on all the parameters of a zero-truncated DEGPD -- and there it can give a
#'   negative \code{df} and an AIC below \code{-2 * logLik}. When \eqn{H_0} has an
#'   eigenvalue below \eqn{-\sqrt{\epsilon}}{-sqrt(eps)} times its largest
#'   (\eqn{\epsilon}{eps} the machine precision), the edf follow the rule of \pkg{mgcv}
#'   for general families: the eigenvalues of \eqn{H_0} below \eqn{t} times its largest
#'   are raised to that level, with \eqn{t = 0}, then \eqn{10^{-7}}{1e-7},
#'   \eqn{1.1 \times 10^{-6}}{1.1e-6}, \ldots until the clamped \eqn{H_0} plus \eqn{S}
#'   has full rank, and the edf are computed from the clamped \eqn{H_0}. Only the edf,
#'   \code{df}, AIC and BIC (and the smooth-term tests of \code{summary()}, which use the
#'   edf) change: \code{Vp} and \code{Vc} remain the inverse of the exact penalised
#'   Hessian.
#'
#'   \code{plateau} is a named logical vector with one element per linear predictor,
#'   \code{TRUE} where the predictor sits on a plateau of its link, so that the
#'   likelihood no longer depends on it and the parameter is not identified: a
#'   predictor other than the scale's, and not held by \code{fixed}, with a fitted value
#'   beyond \eqn{|\eta| > 9} (e.g. \eqn{\kappa} or \eqn{\xi} below about \eqn{10^{-4}}
#'   or above about \eqn{10^4} under a log link, \eqn{\pi} within about \eqn{10^{-4}}
#'   of 0 or 1 under a logit link), or a DEGPD shape on either flat end of the
#'   bounded link (\code{degpd.args$xi.max}), as for \code{restarts}, which compares
#'   such fits by the log-likelihood instead of \code{negREML}. Where any element is
#'   \code{TRUE}, \code{negREML}, and so the selection of the smoothing parameters, is
#'   ill-defined: the Laplace approximation needs a regular maximum, and its value
#'   depends on how far the coefficient search drifted along the plateau (for a
#'   continuous carrier-5 fit with \eqn{\kappa \to 0}{kappa -> 0}, refitting at the
#'   reported \code{sp} changed \code{negREML} from 447.15 to 463.22 at the same
#'   log-likelihood). \code{egpd()} records this without a warning, and
#'   \code{print()} notes it. \code{plateau} is assessed for the families
#'   \code{"egpd"}, \code{"degpd"}, \code{"zidegpd"}, \code{"ztdegpd"},
#'   \code{"gwaring"} and \code{"ztgwaring"}, and is \code{NA} for the others.
#'
#'   \code{hessian.rank} is \code{c(rank, npar)}: the numerical rank of the Hessian
#'   behind \code{Vp} and \code{Vc} (that of the negative log-likelihood -- the adjusted
#'   one, with a sandwich adjustment -- plus the penalty with smooth terms; \code{Vp} is
#'   its pseudo-inverse) and the number of estimated coefficients it covers (all but the
#'   inert ones of parameters held by \code{fixed}). The rank counts the eigenvalues above
#'   \eqn{\sqrt{\epsilon}}{sqrt(eps)} times the largest eigenvalue of the unpenalised
#'   Hessian, with each coefficient measured in units of the root mean square of its
#'   design column, so that the units of the covariates do not matter. A rank below
#'   \code{npar} means that the Hessian is singular or not positive definite: \code{Vp}
#'   gives the unidentified direction(s) a zero, negative or arbitrarily large variance,
#'   so that standard errors, \code{vcov()}, \code{confint()},
#'   \code{predict(se.fit = TRUE)} and \code{simulate()} do not reflect the uncertainty
#'   along them; \code{print()} and \code{summary()} then add a note. For the families
#'   whose Hessian is the outer product of the per-observation scores (BHHH:
#'   \code{"gwaring"}, \code{"ztgwaring"}, \code{"gpois"}, \code{"plnorm"},
#'   \code{"gpig"}, \code{"gpignat"}, \code{"zigpig"}, \code{"zigpignat"}, \code{"cmp"}
#'   and \code{"prl"}, without a sandwich adjustment) it is singular where the
#'   scores are linearly dependent, which
#'   need not make the likelihood flat: an intercept-only generalized Waring fit whose
#'   estimate has \eqn{a = k}, where \eqn{a = \mu(\rho - 1)/k}, has rank 2 of 3, since
#'   its likelihood is symmetric in \eqn{a} and \eqn{k}.
#'
#'   A converged fit warns (condition class \code{"egpd_singular_hessian"}), naming the
#'   coefficients involved, when the singularity is one that neither its standard errors
#'   nor \code{plateau} show. Among the coefficients of the linear predictors that are not
#'   on a plateau (those are reported by \code{plateau} and \code{print()}, without a
#'   warning), it warns when \code{Vp} gives some direction zero or negative variance --
#'   the pseudo-inverse cut it off, as for collinear terms, or a coefficient carries no
#'   information at all, as for a parameter beyond the range in which its likelihood is
#'   evaluated (e.g. a \code{"gpig"} \eqn{a} on its guard bound), or the Hessian is not
#'   positive definite -- and, for the BHHH families, also when the Hessian is merely
#'   numerically singular, as at \eqn{a = k}, where the outer product of the scores is
#'   singular although the likelihood is not flat. A near-singular exact Hessian, whose
#'   very large variances reflect a flat likelihood, is recorded in \code{hessian.rank} and
#'   noted by \code{print()} and \code{summary()} without a warning. Fixing a parameter
#'   (\code{fixed}), a simpler model, or profile likelihood intervals
#'   (\code{confint(method = "profile")}) do not rely on the singular Hessian. For a fit
#'   that did not converge, the rank is still recorded, but only the nonconvergence
#'   warning is given for the estimates.
#'
#' @export
egpd <- function(formula, data, family="egpd", correctV=TRUE, rho0=0,
inits=NULL, outer="bfgs", control=NULL, removeData=FALSE, trace=0,
knots=NULL, maxdata=1e20, maxspline=1e20, compact=FALSE,
egpd.args=list(), degpd.args=list(), zidegpd.args=list(), comppareto.args=list(),
gpig.args=list(), sandwich.args=list(), custom.fns=list(), sp=NULL, gamma=1,
fixed=list(), restarts=TRUE) {

if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0)
  stop("`gamma` must be a single positive number", call. = FALSE)
xi.max <- if (is.null(degpd.args$xi.max)) Inf else degpd.args$xi.max
if (!is.numeric(xi.max) || length(xi.max) != 1L || is.na(xi.max) || xi.max <= 0)
  stop("`degpd.args$xi.max` must be a positive scalar (or Inf)", call. = FALSE)
if (family %in% c("degpd", "ztdegpd") && identical(degpd.args$link, "identity") && is.finite(xi.max))
  stop("finite `xi.max` cannot be combined with the identity shape link", call. = FALSE)
if (family == "degpd" && !is.null(degpd.args$m) &&
    (!is.numeric(degpd.args$m) || length(degpd.args$m) != 1L ||
     is.na(degpd.args$m) || !(degpd.args$m %in% 1:6)))
  stop("degpd$m must be 1, 2, 3, 4, 5, or 6.", call. = FALSE)

## Set up the family first: it is cheap, and a named formula list has to be put in
## canonical order before anything reads formula[[1]]. The restarts below used to read
## it first, so for list(lxi = ~ 1, lsigma = y ~ x, lkappa = ~ 1) they found no
## response, warned about a non-numeric mean and silently fell back to one start.
family.info <- .setup.family.egpd(family, egpd.args, degpd.args, zidegpd.args, comppareto.args, formula, custom.fns)
## A response written as an expression, e.g. I(y + 1), is computed here, once, in the
## environment of the formula that carries it, and becomes a column of `data` that
## every formula then names (.materialise_response) -- before the response is moved to
## or copied into other formulae, which live in other environments.
.resp <- .materialise_response(formula, data, named = family.info$family != "custom")
formula <- .resp$formula
data <- .resp$data
if (family.info$family != "custom")
  formula <- .reorder_named_formulae(formula, family.info$ptab)
.check_not_all_fixed(fixed, family.info$ptab)

## Multi-start for DEGPD, ZIDEGPD and the continuous EGPD. DEGPD model 1's likelihood
## has a nearly flat ridge along which sigma -> 0 and kappa -> infinity together (they
## enter the upper tail only through kappa * sigma^(1/xi)), and from a single start the
## optimiser can settle on it or on a saddle -- which shows up as a wildly large kappa,
## or as a negative variance for a parameter that is genuinely estimated. ZIDEGPD adds
## the zero-inflation boundary pi -> 0 and inherits from its DEGPD warm start a shape
## driven towards xi -> 0, where progress is slow enough to exhaust the iteration
## limit. The continuous carriers with more than one parameter of their own have
## distinct local optima too: from the default start, model 2 (a mixture of two power
## carriers) stopped 1.65 log-likelihood units below the optimum on 400 EGPD draws,
## and model 4 0.08-0.22 units below it on four of nine simulated data sets. Trying a
## few starts and keeping the best of the converged fits can escape these basins.
## Every candidate is the same model fitted to the same data, so this is plain
## multi-start optimisation, not selection across models, and the fits are compared by
## the criterion each of them optimised (.pick_start): the REML criterion when
## smoothing parameters are estimated, not the log-likelihood. Skipped when the caller
## supplies `inits`, or sets restarts = FALSE.
if (isTRUE(restarts) && is.null(inits) && family %in% c("egpd", "degpd", "zidegpd")) {
  cont <- family == "egpd"
  ## `data` is evaluated here, once, and every start below fits that object: each
  ## start used to re-evaluate the caller's expression for it (match.call()), so an
  ## impure one gave every start different data, and the fit with the highest
  ## log-likelihood was then chosen across data sets.
  yv <- tryCatch({
    fl <- if (is.list(formula)) formula[[1]] else formula
    if (length(fl) == 3L) eval(fl[[2L]], data, environment(fl)) else NULL
  }, error = function(e) NULL)
  ybar <- if (is.numeric(yv)) mean(yv, na.rm = TRUE) else NA_real_
  if (is.finite(ybar) && (!cont || ybar > 0)) {
    zi <- family == "zidegpd"
    ## the continuous EGPD's shape xi has an identity link
    idlink <- cont || isTRUE(family.info$lik.fns$identity.xi)
    ## the shape link the starts are written for: a finite xi.max bounds DEGPD's shape
    ## (ZIDEGPD ignores it), and .shape_eta() keeps those starts inside the cap
    xm <- if (family == "degpd" && !idlink) xi.max else Inf
    nrest <- family.info$npar - 2L - zi
    ## a moderate zero-inflation start: at most half of the observed zeros
    pi0 <- min(max(mean(yv == 0, na.rm = TRUE) / 2, 0.01), 0.5)
    ## A light and a heavier shape, alongside the package default -- and, for the
    ## carriers that have more than one shape parameter of their own, two fills for
    ## them as well. Restarts used to be restricted to model 1, which silently left
    ## every other carrier single-start however `restarts` was set. That flattered
    ## model 1 in any comparison: on a Tycho measles series model 4's default start
    ## reaches logLik -40053.7 with xi-hat 5e-08 and an unbounded interval, while a
    ## start at xi = 0.3 reaches -37672.6 with xi-hat 0.324 [0.295, 0.356] -- the same
    ## optimum the other carriers find, and 2,380 units better.
    fills <- lapply(c(0, log(2)), function(cf) rep(cf, nrest))
    ## The continuous mixture carrier (model 2: kappa1, dkappa = kappa2 - kappa1 and
    ## the weight p of kappa1) puts the minority weight on either component instead:
    ## p = 0.12 with kappa2 = 2, and p = 0.88 with a sharp kappa2 = 21. On 12 simulated
    ## data sets these reached the best known optimum in 7 (the fills above in 3).
    if (cont && "logitp" %in% family.info$nms)
      fills <- list(c(0, 0, -2), c(0, 3, 2))
    ## the continuous default starts its scale at log(mean(y)); counts can average 0
    scale0 <- if (cont) log(ybar) else log(ybar + 1)
    starts <- c(list(NULL), unlist(lapply(c(0.05, 0.3), function(v)
                 lapply(fills, function(cf)
                   .eta_start(c(scale0, .shape_eta(v, idlink, xm), cf,
                                if (zi) stats::qlogis(pi0))))),
                 recursive = FALSE))
    ## Every start also draws the same maxdata and maxspline subsamples: the
    ## random-number state is set back before each one (the subsamples used to be
    ## redrawn for every start, and the best log-likelihood then picked the easiest).
    rng0 <- .rng_freeze(NROW(data) > min(maxdata, maxspline))
    ## Each start is a call of egpd() with the arguments as evaluated here (the
    ## canonically ordered formula list, so that names are matched, and any
    ## deprecation reported, once). Re-evaluating match.call() in the caller's frame
    ## used to re-evaluate the caller's expressions as well (and had to find that
    ## frame from inside lapply(), #6).
    fit_start <- function(s0, dargs = degpd.args)
      egpd(formula, data, family = family, correctV = correctV, rho0 = rho0,
           inits = s0, outer = outer, control = control, removeData = removeData,
           trace = trace, knots = knots, maxdata = maxdata, maxspline = maxspline,
           compact = compact, egpd.args = egpd.args, degpd.args = dargs,
           zidegpd.args = zidegpd.args, comppareto.args = comppareto.args,
           gpig.args = gpig.args, sandwich.args = sandwich.args, custom.fns = custom.fns,
           sp = sp, gamma = gamma, fixed = fixed, restarts = FALSE)
    err <- NULL
    cand_warns <- list()
    run_start <- function(k) {
      .rng_restore(rng0)
      ws <- list()
      val <- tryCatch(
        withCallingHandlers(
          fit_start(starts[[k]]),
          warning = function(w) {
            ws[[length(ws) + 1L]] <<- w
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) { err <<- conditionMessage(e); NULL }
      )
      cand_warns[[k]] <<- ws
      val
    }
    cand <- list(run_start(1L))
    ## An interior ZIDEGPD warm start need not be best for the mixture or the delta
    ## carriers. For the other carriers an interior primary fit is rarely improved.
    zi_m <- if (is.null(zidegpd.args$m)) 1L else as.integer(zidegpd.args$m)
    if (zi && !(zi_m %in% 2:4)) {
      why <- .restart_reason(cand[[1L]], family.info$nms)
      if (is.null(why)) {
        for (w in cand_warns[[1L]]) warning(w)
        return(cand[[1L]])
      }
      if (trace > 0) message("ZIDEGPD restarts: the primary fit ", why, ".")
    }
    if (isTRUE(family.info$lik.fns$m == 2L))
      starts <- c(starts, .m2_restart_starts(yv, cont, zi, idlink, xm, cand[[1L]]))
    ## Starts that differ only for parameters held by `fixed` are the same start (the
    ## pinned parameters' coefficients are inert), so they are fitted once.
    if (length(fixed)) {
      pin <- unique(suppressWarnings(vapply(names(fixed), .match_param_name, integer(1),
                                            ptab = family.info$ptab)))
      starts <- unique(lapply(starts, function(s0) if (is.null(s0)) s0 else replace(s0, pin, 0)))
    }
    for (k in seq_along(starts)[-1L]) cand[k] <- list(run_start(k))
    converged <- function()
      which(vapply(cand, function(z) !is.null(z) && isTRUE(z$convergence == 0L) &&
                     is.finite(as.numeric(z$logLik)), logical(1)))
    ## The best converged fit by the criterion the fits optimised (.pick_start). For the
    ## continuous EGPD the default fit stands unless another start improves on it by
    ## more than 0.01: the other starts mostly reach the same optimum, where a smooth
    ## model's log-likelihood still varies by up to about 4e-3 with the smoothing
    ## parameters the outer search stops at, and a gain below 0.01 is immaterial anyway.
    best_of <- function(keep) .pick_start(cand, keep, prefer_first = cont)
    keep <- converged()
    ## Under the bounded shape link (finite xi.max) a fit can stop on either flat
    ## asymptote of the link, xi -> 0 or xi -> xi.max, where the gradient in eta
    ## vanishes whether or not the likelihood could still increase (.shape_plateau):
    ## the DEGPD Hessian is often indefinite on the way there, and the near-singular
    ## Newton steps that follow can land on a plateau. With the shape starts written
    ## on the log scale, every start did so on 6 of 360 simulated data sets, 0.08-1.4
    ## log-likelihood units below the optimum. The likelihood is the same function of
    ## xi under the log link, which has no plateau at xi.max, so when the best fit is
    ## on a plateau a log-link fit, mapped onto the bounded link, is tried as one more
    ## start: it reaches the interior optimum if there is one, and goes back to the
    ## boundary if the optimum is there. The log-link fit itself is another model (its
    ## smoothing parameters, and so its REML criterion and penalty, belong to another
    ## link) and only supplies that start: it is never compared with the candidates.
    if (is.finite(xm) && length(keep)) {
      best <- best_of(keep)
      why <- .shape_plateau(cand[[best]], xm)
      if (!is.null(why)) {
        if (trace > 0)
          message("DEGPD restarts: the best fit ", why, "; trying the log-link fit as a start.")
        ufit <- tryCatch(suppressWarnings({
          .rng_restore(rng0)
          fit_start(NULL, dargs = utils::modifyList(degpd.args, list(xi.max = Inf)))
        }), error = function(e) NULL)
        if (!is.null(ufit) && isTRUE(ufit$convergence == 0L)) {
          k <- length(starts) + 1L
          starts[[k]] <- .intercept_starts(ufit, shape = function(e)
            .shape_eta(exp(e), xi.max = xm, cap = 1 - 1e-3))
          cand[k] <- list(run_start(k))
          keep <- converged()
        }
      }
    }
    if (length(keep)) {
      best <- best_of(keep)
      for (w in cand_warns[[best]]) warning(w)
      return(cand[[best]])
    }
    ## Never degrade quietly: a silent fall-through to one start is the failure mode
    ## that made #6 invisible, and it can cost thousands of log-likelihood units. Say
    ## whether the starts failed outright: "no restart converged" hid that every start
    ## had stopped with the same error.
    nerr <- sum(vapply(cand, is.null, logical(1)))
    warning("restarts requested but ",
            if (nerr == length(cand)) "every start failed with an error"
            else if (nerr > 0L) sprintf("no restart converged (%d of %d failed with an error)",
                                        nerr, length(cand))
            else "no restart converged",
            "; continuing from a single start", if (!is.null(err)) paste0(" (", err, ")"),
            call. = FALSE)
    .rng_restore(rng0)
  }
}

## Screen on the actual rows of the first fit, after response materialisation,
## complete-case filtering and maxdata. Each full refit repeats the same setup.
if (isTRUE(restarts) && is.null(inits) && family == "comppareto") {
  rng0 <- .rng_freeze(NROW(data) > min(maxdata, maxspline))
  fit_splice_start <- function(s) {
    .rng_restore(rng0)
    ws <- list()
    f <- tryCatch(withCallingHandlers(
      egpd(formula, data, family = family, correctV = correctV, rho0 = rho0,
           inits = s, outer = outer, control = control,
           removeData = removeData, trace = trace, knots = knots, maxdata = maxdata,
           maxspline = maxspline, compact = compact, egpd.args = egpd.args,
           degpd.args = degpd.args, zidegpd.args = zidegpd.args,
           comppareto.args = comppareto.args, gpig.args = gpig.args,
           sandwich.args = sandwich.args, custom.fns = custom.fns, sp = sp,
           gamma = gamma, fixed = fixed, restarts = FALSE),
      warning = function(w) {
        ws[[length(ws) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }), error = identity)
    list(fit = f, warnings = ws)
  }
  primary <- fit_splice_start(NULL)
  f1 <- primary$fit
  if (inherits(f1, "error")) {
    for (w in primary$warnings) warning(w)
    stop(f1)
  }
  starts <- tryCatch(suppressWarnings(.comppareto_split_starts(f1)), error = function(e) {
    if (trace > 0) message("CompPareto split screen: ", conditionMessage(e))
    list()
  })
  why <- .comppareto_restart_reason(f1)
  if (!is.null(why)) {
    if (trace > 0) message("CompPareto restarts: the default fit ", why, ".")
    starts <- c(lapply(.comppareto_restart_starts(f1$likfns$spec, f1$likdata$y[, 1L]),
                       .eta_start), starts)
  }
  pin <- which(vapply(f1$likdata$X, .is_pinned_block, logical(1)))
  starts <- unique(lapply(starts, function(s) replace(s, pin, 0)))
  cand <- c(list(primary), lapply(starts, fit_splice_start))
  fits <- lapply(cand, `[[`, "fit")
  best <- 1L
  if (length(starts)) {
    keep <- which(vapply(fits, function(z) isTRUE(z$convergence == 0L) &&
                           !.egpd_invalid_objective(-z$logLik), logical(1)))
    if (length(keep)) best <- .pick_start(fits, keep, prefer_first = TRUE)
    ## Do not trade a failed primary for a lower, apparently converged, fit.
    if (best != 1L && !isTRUE(f1$convergence == 0L)) {
      nll <- f1$likfns$d0(f1$coefficients, f1$likdata)
      if (!.egpd_invalid_objective(nll) && fits[[best]]$logLik <= -nll + 0.01) best <- 1L
    }
  }
  selected <- .comppareto_check_split(fits[[best]])
  reclassified <- !identical(selected$convergence.report, fits[[best]]$convergence.report)
  ## A split with no finite optimum supersedes the earlier optimizer/rank
  ## diagnosis. Keep unrelated warnings (such as omitted rows).
  for (w in cand[[best]]$warnings)
    if (!reclassified || !inherits(w, c("egpd_nonconvergence", "egpd_singular_hessian"))) warning(w)
  if (reclassified)
    .warn_nonconvergence(selected$convergence, selected$convergence.report, family)
  return(selected)
}

family <- family.info$family
if (!(family %in% c("degpd", "ztdegpd")) || isTRUE(family.info$lik.fns$identity.xi)) xi.max <- Inf

## `fixed`: hold named parameters at known values. This is sugar for the offset-only
## pinning `~ offset(rep(<link>(v), n)) - 1`, which is easy to get wrong by hand because
## the offset must be on the LINK scale (log, logit or identity, according to the
## parameter's link). Values here are given on the RESPONSE scale, so fixed = list(kappa = 1)
## pins kappa itself rather than log kappa.
if (length(fixed)) {
  .fx <- .setup.fixed.egpd(fixed, formula, data, family.info$nms, xi.max = xi.max,
                           ptab = family.info$ptab,
                           unlink = if (family == "custom") family.info$lik.fns$unlink)
  formula <- .fx$formula
  data    <- .fx$data
}

## GPIG/ZIGPIG pmf evaluation route (see gpig_control). Set for the duration of
## this fit only, so a method chosen here cannot leak into later calls.
.gpig.old <- .gpig_apply_control(family, gpig.args)
if (!is.null(.gpig.old)) on.exit(do.call(gpig_control, .gpig.old), add = TRUE)

## setup formulae
formula <- .setup.formulae(formula, family.info$npar, family.info$npar2, data, trace, family.info$nms)
response.name <- attr(formula, "response.name")

## ZIDEGPD nests DEGPD (zero-inflation pi -> 0). Warm-start the richer model from
## the DEGPD MLE so a converged ZIDEGPD can never sit below the simpler model's
## likelihood (the nesting guarantee). The last ZIDEGPD parameter is always
## logitpi, so dropping it yields the matching DEGPD formula (reordering
## kappa/delta for model 4). Seed (sigma, xi, kappa, ...) from DEGPD's
## per-parameter linear-predictor means and pi near 0.
if (family == "zidegpd" && is.null(inits) && length(response.name) == 1) {
  inits <- tryCatch({
    zi_m <- if (is.null(zidegpd.args$m)) 1L else as.integer(zidegpd.args$m)
    dform <- unname(if (zi_m == 4L) formula[c(1L, 2L, 4L, 3L)] else formula[-length(formula)])
    ## The warm start must not move the random-number state: the ZIDEGPD fit then
    ## draws the same maxdata/maxspline subsamples whether or not it is warm-started,
    ## as every restart must (the other starts are not).
    rng <- .rng_freeze(NROW(data) > min(maxdata, maxspline))
    dfit <- tryCatch(suppressWarnings(egpd(
      dform, data, family = "degpd",
      correctV = FALSE, inits = NULL, outer = outer, control = control,
      knots = knots, maxdata = maxdata, maxspline = maxspline, compact = compact,
      degpd.args = list(m = zi_m), trace = 0, gamma = gamma, restarts = FALSE
    )), finally = .rng_restore(rng))
    lp_means <- .intercept_starts(dfit)
    if (zi_m == 4L) lp_means <- lp_means[c(1L, 2L, 4L, 3L)]
    p0 <- min(max(mean(data[[response.name]] == 0, na.rm = TRUE), 1e-3), 0.5)
    .eta_start(c(lp_means, qlogis(p0)), idx = length(lp_means) + 1L)
  }, error = function(e) {
    if (trace > 0) message("ZIDEGPD warm-start from DEGPD failed; using default inits.")
    NULL
  })
}

## setup mgcv objects and data
temp.data <- .setup.data(data, response.name, formula, family, family.info$nms,
  removeData, knots, maxdata, maxspline, compact, sandwich.args,
  tolower(outer), trace, gamma)
data <- temp.data$data

## Only DEGPD's log-shape link supports a finite cap.
temp.data$lik.data$xi.max <- xi.max
temp.data$lik.data$bounded.xi <- is.finite(xi.max)
## The analytic log-determinant contractions combine d340, which has one row per
## observation, with the design rows, of which a compact fit keeps one per distinct
## covariate row -- for every family that has third derivatives (continuous EGPD
## m = 1, 3 as well as DEGPD), not only DEGPD: a continuous compact smooth fit failed
## with "dims do not match", or, when the row counts happened to be commensurate,
## silently used a wrong REML gradient. Compact fits, and DEGPD's bounded shape link,
## therefore use finite-difference REML derivatives.
if (isTRUE(temp.data$lik.data$duplicate == 1) ||
    (family == "degpd" && isTRUE(temp.data$lik.data$bounded.xi) &&
     !isTRUE(family.info$lik.fns$identity.xi)))
  family.info$lik.fns$d340 <- NULL

## Optimizers read lik.data$control, so merge overrides there before initialization.
if (!is.null(control)) {
  if (!is.list(control)) stop("`control' must be a list", call. = FALSE)
  if (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || any(!nzchar(names(control))) ||
      anyDuplicated(names(control))))
    stop("`control' must have unique, nonempty names", call. = FALSE)
  unknown <- setdiff(names(control), c("outer", "inner"))
  if (length(unknown))
    stop("unknown `control' element(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  for (.slot in names(control)) {
    .override <- control[[.slot]]
    if (!is.list(.override) || (length(.override) &&
        (is.null(names(.override)) || anyNA(names(.override)) ||
         any(!nzchar(names(.override))) || anyDuplicated(names(.override)))))
      stop("`control$", .slot, "' must be a named list", call. = FALSE)
    unknown <- setdiff(names(.override), names(temp.data$lik.data$control[[.slot]]))
    if (length(unknown))
      stop("unknown `control$", .slot, "' element(s): ",
           paste(unknown, collapse = ", "), call. = FALSE)
    if (any(!vapply(.override, function(x)
        is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, logical(1))))
      stop("`control$", .slot, "' values must be finite positive numbers", call. = FALSE)
    temp.data$lik.data$control[[.slot]] <-
      utils::modifyList(temp.data$lik.data$control[[.slot]], .override)
  }
  ## an explicit inner gradtol is used as given, including by the finite-difference
  ## REML derivatives, which otherwise tighten the default (see .tighten_inner)
  temp.data$lik.data$gradtol.user <- !is.null(control$inner$gradtol)
}

## initialise inner iteration
beta <- .setup.inner.inits(inits, temp.data$lik.data, family.info$lik.fns, family.info$npar, family,
                           trace = trace)
lik.data <- .sandwich(temp.data$lik.data, beta)
if (trace > 0 & lik.data$adjust > 0) cat(paste("\n Sandwich correct lambda =", signif(lik.data$k, 3), "\n"))

## An fx smooth needs its basis but has no penalty. Conversely, a term-fixed
## penalty still enters the fit even when none of its multipliers are estimated.
has.penalty <- any(vapply(temp.data$gams, function(g) length(g$lsp0) > 0L, logical(1)))

if (has.penalty) {

## initialise outer iteration
S.data <- .joinSmooth(temp.data$gams, gamma)
nsp <- length(attr(S.data, "Sl"))
if (!is.null(sp) || nsp == 0L) {
  rho0 <- if (is.null(sp)) numeric(0) else .fixed_rho(sp, S.data, gamma)
  lik.data$outer <- "fixed"
} else if (is.null(rho0)) {
  rho0 <- vapply(attr(S.data, "Sl"), function(S)
    uniroot(.guess, c(-1e2, 1e2), d=attr(beta, "diagH"), s=diag(S))$root, 0)
} else {
  if (!is.numeric(rho0) || !is.null(dim(rho0)) || !length(rho0) %in% c(1L, nsp) ||
      any(!is.finite(rho0)))
    stop("`rho0` must give a finite start per independent smoothing parameter, or a scalar.",
         call. = FALSE)
  if (length(rho0) == 1) rho0 <- rep(rho0, nsp)
}

lik.data$S <- .makeS(S.data, exp(rho0))

## perform outer iteration
fit.reml <- .outer(rho0, beta, family.info$lik.fns, lik.data, S.data, control, correctV, lik.data$outer, trace)

sp <- exp(fit.reml$par)
lik.data$S <- .makeS(S.data, sp)

} else {

if (length(sp)) stop("`sp` supplied for a model without penalties.", call. = FALSE)
S.data <- NULL
fit.reml <- .outer.nosmooth(beta, family.info$lik.fns, lik.data, control, trace)

}

## covariance matrices
VpVc <- .VpVc(fit.reml, family.info$lik.fns, lik.data, S.data, correctV=correctV, sandwich=temp.data$sandwich, smooths=has.penalty, trace=trace)

## effective degrees of freedom
edf <- .edf(fit.reml$beta, family.info$lik.fns, lik.data, VpVc, temp.data$sandwich)

## update mgcv objects
names(temp.data$gams) <- family.info$nms
gams <- .swap(fit.reml, temp.data$gams, lik.data, VpVc, edf, has.penalty)

## add extra things that make an egpd object
gams <- .finalise(gams, data, family.info$lik.fns, lik.data, S.data, fit.reml, VpVc, family, temp.data$gotsmooth, formula, response.name, removeData, edf, family.info$nms2)
if (family %in% .gpig_families) gams$gpig.control <- .gpig_settings()

## record any `fixed` pins so predict() can re-create their synthetic columns for
## newdata, and so the constraint is visible on the fitted object
if (length(fixed)) {
  gams$fixed.offsets <- .fx$offsets
  gams$fixed         <- .fx$values
}

## Warn on the DEGPD kappa ridge. sigma and kappa enter the upper tail only through
## kappa * sigma^(1/xi), so they are separated by the bulk alone; on a short series
## that separation is weak and the optimum drifts along a flat ridge to very large
## kappa and near-zero sigma. xi usually survives this -- it is not part of the
## confounding -- but sigma and kappa individually do not, so quantiles and any
## interval that leans on them should be treated as provisional.
if (identical(attr(family, "type"), 1) && grepl("degpd", family) && !removeData) {
  kap <- tryCatch(suppressWarnings(max(as.data.frame(
           stats::predict(gams, type = "response"))$kappa, na.rm = TRUE)),
         error = function(e) NA_real_)
  if (is.finite(kap) && kap > 1e4)
    warning(sprintf(paste0("fitted kappa reaches %.3g: the DEGPD scale/carrier pair is ",
                           "effectively unidentified (the sigma-kappa ridge). xi is ",
                           "usually still identified; sigma and kappa are not. Consider ",
                           "fixed = list(kappa = 1) for a plain discretised GPD."), kap),
            call. = FALSE)
}

return(gams)
}

## Why a primary ZIDEGPD fit needs the other restart starts, or NULL if it does not:
## it failed or did not converge, its log-likelihood is not finite, or one of its
## shape-type linear predictors (every parameter but the scale, skipping fixed ones)
## has reached a boundary of its link -- |eta| > `bound`, e.g. pi, kappa or xi below
## 1.2e-4 -- where the fit may be a degenerate limit with a better optimum elsewhere.
.restart_reason <- function(fit, nms, bound = 9) {
  if (is.null(fit)) return("failed")
  if (!isTRUE(fit$convergence == 0L))
    return(sprintf("did not converge (code %s)", format(fit$convergence)))
  if (!is.finite(as.numeric(fit$logLik))) return("has a non-finite log-likelihood")
  fixed <- match(names(fit$fixed), nms)
  for (i in setdiff(seq_len(fit$ngam)[-1L], fixed))
    if (any(abs(fit[[i]]$fitted) > bound))
      return(sprintf("has %s on a boundary (|eta| > %g)", names(fit)[i], bound))
  NULL
}

## Which of the converged fits `fits[keep]` -- one model fitted to the same data from
## different starts -- egpd() (restarts) and fit_hurdle() keep; returns an element of
## `keep`. The fits are compared by the criterion they optimised (.start_criterion), not
## by the log-likelihood: with estimated smoothing parameters a start with a higher
## log-likelihood can be a wigglier fit (smaller sp, more edf) whose REML criterion is
## worse. Fits within `tol` of the best by that criterion reached the same optimum as far
## as it can tell (it still varies with the tolerances at which the searches stop), and
## among them the one with the highest log-likelihood is kept, as it always was, so the
## choice among starts that agree is unchanged -- or, with `prefer_first`, the first fit
## (the default start's) whenever it is among them, so that another start replaces it
## only if it improves the criterion by more than `tol`. Without smooth terms, and with a
## sandwich adjustment, the criterion is the log-likelihood, and the choice is exactly
## what it always was.
.pick_start <- function(fits, keep, prefer_first = FALSE, tol = 0.01) {
  ll <- vapply(fits[keep], function(z) as.numeric(z$logLik), 0)
  crit <- .start_criterion(fits[keep], ll)
  near <- max(crit) - crit <= tol
  if (prefer_first && keep[1L] == 1L && near[1L])
    return(1L)
  keep[near][which.max(ll[near])]
}

## The criterion that converged fits of one model from different starts optimised, one
## value per fit, larger being better (`ll`: their log-likelihoods):
## - without smooth terms, the log-likelihood;
## - with a sandwich adjustment (`sandwich.args`), the log-likelihood as well. The
##   adjustment -- the magnitude k, or the curvature transformation CH -- is estimated
##   from each start's own unpenalised fit (.setup.inner.inits), so it differs between
##   starts, and the adjusted penalised log-likelihood and REML criterion that each fit
##   optimised are on different scales: they preferred the start with the smallest k
##   whatever its fit (a default start with k = 0.47, a log-likelihood 0.61 lower and
##   standard errors 3 times larger, over starts with k = 1.39; or, with estimated sp,
##   a REML criterion 12 units lower of which 10 came from k = 0.98 against 1.01);
## - with the smoothing parameters held fixed (`sp` supplied, or outer = "fixed"), the
##   penalised log-likelihood, which the coefficients maximise (.pen_loglik);
## - with them estimated, minus the REML criterion `negREML`, which the smoothing
##   parameters minimise (with gamma != 1, that of the log-likelihood divided by gamma,
##   as optimised). It is a Laplace approximation, which needs each fit to be a regular
##   maximum of its penalised likelihood (.reml_regular). On a boundary of a link (e.g.
##   DEGPD's xi -> 0 or ZIDEGPD's pi -> 0), or where the penalised Hessian is not
##   positive definite, the likelihood is flat in some direction, and the criterion falls
##   without bound the further a fit drifted along it -- by half a unit per unit of log
##   xi in the DEGPD case -- so it would rank the drift, not the fits. Of five fits of a
##   ZIDEGPD model to data without excess zeros, all ending at pi -> 0, the REML criterion
##   preferred by 6.6 units one whose smooth had collapsed (sp 726, edf 4.0) and whose
##   log-likelihood was 7.6 units lower. If any fit is irregular, the fits are therefore
##   compared by the log-likelihood, as without smooths.
.start_criterion <- function(fits, ll) {
  f1 <- fits[[1L]]
  if (is.null(f1$Sdata) ||
      any(vapply(fits, function(z) isTRUE(z$likdata$adjust > 0), logical(1))))
    return(ll)
  if (identical(f1$likdata$outer, "fixed")) {
    pl <- vapply(fits, .pen_loglik, 0)
    return(if (all(is.finite(pl))) pl else ll)
  }
  if (!all(vapply(fits, .reml_regular, logical(1))))
    return(ll)
  -vapply(fits, function(z) as.numeric(z$negREML), 0)
}

## A fit's penalised log-likelihood at its coefficients and smoothing parameters: the
## log-likelihood less the penalty, on the scale of the log-likelihood. The objective
## minimised (.nllh.pen) divides the log-likelihood by gamma and carries the penalty at
## sp / gamma, so it is multiplied back by gamma: with `sp` supplied, gamma does not
## change the fit, nor so the choice among starts. (Sandwich-adjusted fits, whose
## objective also carries a start-dependent adjustment, are not compared by it.)
.pen_loglik <- function(fit) {
  g <- if (is.null(fit$likdata$gamma)) 1 else fit$likdata$gamma
  v <- .nllh.pen(as.vector(fit$coefficients), fit$likdata, fit$likfns)
  if (.egpd_invalid_objective(v)) NA_real_ else -g * as.numeric(v)
}

## Is a converged fit with estimated smoothing parameters a regular maximum of its
## penalised likelihood, at which its REML criterion can be compared with another fit's?
## Not if the criterion is missing; if a linear predictor other than the scale's (and
## other than one held by `fixed`) has reached a boundary of its link, |eta| > `bound`
## (as for .restart_reason); if a bounded DEGPD shape sits on a plateau of its link
## (.shape_plateau); or if the penalised Hessian at the coefficients is not numerically
## positive definite (the test of .VpVc()).
.reml_regular <- function(fit, bound = 9) {
  if (.egpd_invalid_objective(fit$negREML))
    return(FALSE)
  X <- fit$likdata$X
  for (i in seq_len(fit$ngam)[-1L])
    if (ncol(X[[i]]) && !.is_pinned_block(X[[i]]) && isTRUE(any(abs(fit[[i]]$fitted) > bound)))
      return(FALSE)
  if (!is.null(.shape_plateau(fit, fit$likdata$xi.max)))
    return(FALSE)
  H <- tryCatch(.gH.pen(as.vector(fit$coefficients), fit$likdata, fit$likfns)[[2]],
                error = function(e) NULL)
  if (is.null(H) || anyNA(H))
    return(FALSE)
  ## the inert coefficient of a parameter held by `fixed` has a zero row and column
  inert <- which(rowSums(abs(H)) == 0 & colSums(abs(H)) == 0)
  if (length(inert)) diag(H)[inert] <- 1
  !inherits(try(chol(H), silent = TRUE), "try-error")
}

## The random-number state that each restart starts from, so that every start draws
## the same maxdata and maxspline subsamples, or NULL when no subsample is drawn
## (`need` FALSE). A state is created first if there is none (as simulate() does):
## otherwise every start would seed the generator afresh.
.rng_freeze <- function(need) {
  if (!isTRUE(need))
    return(NULL)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    stats::runif(1)
  get(".Random.seed", envir = .GlobalEnv)
}

.rng_restore <- function(state) {
  if (!is.null(state))
    assign(".Random.seed", state, envir = .GlobalEnv)
  invisible(NULL)
}

## Starting intercepts from a fitted model: each parameter's mean linear predictor,
## without its offsets, since the new fit adds them back (averaging
## predict(type = "link") counted every offset twice). `shape`, if given, first maps
## the shape (second) predictor onto another link. The linear predictors are those the
## fit stores (`$fitted`, with offsets, as predict(type = "link") gives them), so that a
## fit with removeData = TRUE supplies starts too: predict() needs the data, and the
## shape-plateau restart of egpd() crashed on such fits.
.intercept_starts <- function(fit, shape = NULL) {
  eta <- lapply(fit[seq_along(fit$likdata$X)], `[[`, "fitted")
  if (!is.null(shape)) eta[[2L]] <- shape(eta[[2L]])
  .start_coefficients(eta, .start_model(fit$likdata))
}

#' Extract Model Fitted Values
#'
#' @param object a fitted \code{egpd} object
#' @param type \code{"link"} (the default, as in \pkg{evgam}): the fitted linear
#'   predictors, one column per distributional parameter named by its link (e.g.
#'   \code{logscale}, \code{logshape}, \code{logitpi}), including offsets;
#'   \code{"response"}: the fitted distributional parameters on their natural
#'   scale, as \code{predict(object, type = "response")} (columns named by
#'   parameter, e.g. \code{scale}, \code{shape}, plus an \code{xi} alias of
#'   \code{shape} for EGPD families). Neither is a fitted mean of the response.
#' @param ... not used
#'
#' @return A data frame with one row per observation used in the fit.
#'
#' @export
fitted.egpd <- function(object, type = c("link", "response"), ...) {
type <- match.arg(type)
predict(object, type = type)
}

#' Simulations from a fitted \code{egpd} object
#'
#' Draws coefficient vectors from their approximate posterior, a multivariate
#' normal with the fitted mean and covariance \code{Vc} (\code{marginal = TRUE})
#' or \code{Vp}, and returns the implied linear predictors, distributional
#' parameters, or quantiles.
#'
#' @param object a fitted \code{egpd} object
#' @param nsim an integer giving the number of simulations
#' @param seed an integer giving the seed for simulations
#' @param newdata a data frame
#' @param type a character string: \code{"link"}, \code{"response"} or
#'   \code{"quantile"}; supplying \code{probs} implies \code{"quantile"}
#' @param probs a scalar or vector of probabilities: for each coefficient draw, the
#'   quantiles of the response at these probabilities are returned
#' @param threshold a scalar added to the simulated quantiles (e.g. the threshold
#'   of a peaks-over-threshold model); ignored for other types
#' @param marginal logical: should uncertainty integrate out smoothing parameter
#'   uncertainty (use \code{Vc}; the default) rather than condition on it
#'   (\code{Vp})?
#' @param ... additional arguments
#'
#' @return A named list of observation-by-simulation matrices: one per linear
#'   predictor (\code{type = "link"}: draws on the link scale, including model
#'   offsets, named by parameter as for \code{type = "response"}, e.g.
#'   \code{scale} for the \code{logscale} predictor), one
#'   per distributional parameter (\code{type = "response"}, named by parameter;
#'   the fitted inverse links are applied, including a finite DEGPD shape cap and,
#'   for \code{family = "custom"}, \code{custom.fns$unlink}, and
#'   parameters are clamped to the same safe interior as in the likelihood and in
#'   \code{predict()}), or one per probability (\code{type = "quantile"}, named
#'   \code{"q:<prob>"} with the probability to 15 significant digits, e.g.
#'   \code{"q:0.9995"}, with \code{threshold} added).
#'
#' @export
simulate.egpd <- function(object, nsim=1e3, seed=NULL, newdata,
  type="link", probs=NULL, threshold=0, marginal=TRUE, ...) {
if (!is.null(probs))
  type <- "quantile"
type <- match.arg(type, c("link", "response", "quantile"))
if (type == "quantile" && is.null(probs))
  stop("non-NULL `probs' required if `type = \"quantile\"'", call. = FALSE)
if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  runif(1)
if(is.null(seed)) {
  RNGstate <- get(".Random.seed", envir = .GlobalEnv)
} else {
  R.seed <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(seed)
  RNGstate <- structure(seed, kind = as.list(RNGkind()))
  on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
}
family <- object$family
if (marginal) {
  V.type <- "Vc"
} else {
  V.type <- "Vp"
}
B <- .pivchol_rmvn(nsim, object$coefficients, object[[V.type]])
idpars <- object$idpars
X <- predict.egpd(object, newdata, type="lpmatrix")
nms <- names(X)
offsets <- attr(X, "offsets")
B <- lapply(seq_along(X), function(i) B[idpars == i, , drop=FALSE])
X <- lapply(seq_along(X), function(i) {
  eta <- X[[i]] %*% B[[i]]
  if (length(offsets[[i]])) eta <- eta + offsets[[i]]
  eta
})
names(X) <- nms
if (type %in% c("response", "quantile")) {
  for (i in seq_along(nms))
    X[[i]] <- .egpd_invlink(object, nms[i], X[[i]])
}
nms <- .strip_link(nms)
names(X) <- nms
if (type == "link")
  return(X)
X <- .egpd_clamp_response(X, family)
if (type == "response")
  return(X)
## quantiles for every (observation, draw) pair at once
dims <- dim(X[[1]])
q <- .egpd_quantile(object, lapply(X, as.vector), probs)
out <- lapply(seq_along(probs), function(j) matrix(q[, j], dims[1], dims[2],
                                                  dimnames = dimnames(X[[1]])) + threshold)
names(out) <- .egpd_qnames(probs)
out
}

#' Log-likelihood from a fitted \code{egpd} object
#'
#' @param object a fitted \code{egpd} object
#' @param ... not used
#'
#' @return A \code{logLik} object
#'
#' @export
logLik.egpd <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
out <- object$logLik
attr(out, "df") <- attr(object, "df")
attr(out, "nobs") <- nobs(object)
class(out) <- "logLik"
out
}

#' Plot a fitted \code{egpd} object
#'
#' @param x a fitted \code{egpd} object
#' @param onepage logical: should all plots be on one page? Defaults to TRUE
#' @param which a vector of integers identifying which smooths to plot
#' @param main a character string or vector of plot titles, one per plot drawn
#'   (recycled); by default each smooth is titled with its parameter and label
#' @param ask logical: ask to show next plots?
#' @param ... extra arguments to pass to \code{mgcv::plot.gam}
#'
#' @return Plots representing smooth terms. A model without smooth terms has
#'   nothing to plot, which is an error.
#'
#' @export
plot.egpd <- function(x, onepage = TRUE, which = NULL, main, ask = !onepage, ...) {
x <- lapply(x[x$gotsmooth], .gam_prediction_basis)
nsmooth <- sum(unlist(lapply(x, function(x) as.integer(sapply(x$smooth, function(y) y$plot.me)))))
if (nsmooth == 0L)
  stop("No smooth terms to plot: the model has no smooth terms.", call. = FALSE)
if (is.null(which)) {
  nplot <- nsmooth
  which <- seq_len(nplot)
} else {
  nplot <- length(which)
}
if (onepage) {
  omfrow <- par("mfrow")
  on.exit(par(mfrow = omfrow), add = TRUE)
  nmfrow <- rev(n2mfrow(nplot))
  par(mfrow = nmfrow)
}
if (ask) {
  oask <- par("ask")
  if (nplot > prod(par("mfrow")) && dev.interactive()) {
    on.exit(par(ask = oask), add = TRUE)
    par(ask = TRUE)
  }
}
current <- 1
for (i in seq_along(x)) {
  for (j in seq_along(x[[i]]$smooth)) {
    if (current %in% which) {
      if (missing(main)) {
        mgcv::plot.gam(x[[i]], select = j, main = paste(names(x)[i], x[[i]]$smooth[[j]]$label, sep = ": "), ...)
      } else {
        ## `main` is a formal, so it is not in `...`: pass the title of this plot
        mgcv::plot.gam(x[[i]], select = j,
                       main = rep_len(main, length(which))[match(current, which)], ...)
      }
    }
    current <- current + 1
  }
}
invisible(NULL)
}

## A predictor's fitted mgcv object as mgcv's own code (plot.gam) expects it: with the
## coefficients and their covariance matrices in the prediction basis of PredictMat().
## egpd keeps them in the basis the smooths were fitted in, which differs for t2() terms
## (see .predictable.gam); mgcv::gam() maps them in the same way, beta = P beta_fit and
## V = P V_fit P'.
.gam_prediction_basis <- function(g) {
  P <- g[["P"]]
  if (is.null(P)) return(g)
  b <- g$coefficients
  g$coefficients <- stats::setNames(as.vector(P %*% b), names(b))
  for (v in c("Vp", "Vc", "Ve"))
    if (!is.null(g[[v]])) g[[v]] <- P %*% g[[v]] %*% t(P)
  g$P <- NULL
  g
}

#' Summary method for a fitted \code{egpd} object
#'
#' @param object a fitted \code{egpd} object
#' @param ... not used
#'
#' @return A \code{summary.egpd} object: a list whose first element holds, for each
#'   distributional parameter, a table of the parametric (non-smooth) coefficients
#'   with columns \code{Estimate}, \code{Std. Error} (from \code{Vp}),
#'   \code{z value} and \code{Pr(>|z|)}, the two-sided Wald p-value
#'   \eqn{2\Phi(-|z|)} (the likelihood has no free scale parameter, so the
#'   reference distribution is normal, as in \code{mgcv::summary.gam} for
#'   known-scale families); and whose second element holds, for each parameter
#'   with smooth terms, a table with columns \code{edf}, \code{max.df},
#'   \code{Chi.sq} and \code{p-value}. For ordinary smooths these are Wood's
#'   smooth-term tests. For a fully penalised term such as \code{s(g, bs="re")},
#'   zero variance is on the boundary. Its approximate test uses a weighted
#'   sum of chi-squares, conditional on the fitted smoothing parameters, with
#'   nuisance terms from all predictors accounted for (as in \pkg{mgcv}'s
#'   random-effect test). Neither \code{edf} nor \code{max.df} is the reference
#'   chi-square degrees of freedom for this test. If the null law is unavailable
#'   (for example after a cluster adjustment, without usable joint information,
#'   or with no estimable contrast), the statistic and p-value are \code{NA}
#'   and a warning names the term and reason. Its element
#'   \code{hessian.rank} is that of the fit (see \code{\link{egpd}}, Value): when the
#'   Hessian behind the standard errors is singular, printing the summary adds a
#'   note under the tables.
#'
#' @export
summary.egpd <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
out <- list()
out[[1]] <- .parametric.summary.egpd(object)
out[[2]] <- .smooth.summary.egpd(object)
out$hessian.rank <- object$hessian.rank
class(out) <- "summary.egpd"
out
}

#' @param x a \code{summary.egpd} object
#' @rdname summary.egpd
#' @export
print.summary.egpd <- function(x, ...) {
if (!missing(...)) warning("extra arguments discarded")
cat("\n")
cat("** Parametric terms **")
tab <- lapply(x[[1]], .tidyParametricTable)
cat("\n")
for (i in seq_along(tab)) {
cat("\n")
cat(names(tab)[[i]])
cat("\n")
print(tab[[i]])
}
cat("\n")
cat("** Smooth terms **")
tab <- lapply(x[[2]], .tidySmoothTable)
cat("\n")
for (i in seq_along(tab)) {
cat("\n")
cat(names(tab)[[i]])
cat("\n")
print(tab[[i]])
}
note <- .hessian_note(x$hessian.rank)
if (!is.null(note)) cat("\n", note, "\n", sep = "")
invisible(x)
}

#' Print a fitted \code{egpd} object
#'
#' Prints the formulae of the fit, a note if a linear predictor sits on a plateau
#' of its link (component \code{plateau}, see \code{\link{egpd}}), and a note if the
#' Hessian behind the standard errors is singular (component \code{hessian.rank}).
#'
#' @param x a fitted \code{egpd} object
#' @param ... not used
#'
#' @return The call of the object (invisibly)
#'
#' @export
print.egpd <- function(x, ...) {
if (!missing(...)) warning("extra arguments discarded")
print(x$call)
if (isTRUE(any(x$plateau)))
  cat("Note: ", paste(names(x$plateau)[which(x$plateau)], collapse = ", "),
      " on a plateau of the link (`plateau`, see ?egpd): not identified there",
      if (!is.null(x$Sdata)) paste0(", and negREML and the smoothing-parameter ",
                                      "selection are ill-defined"),
      ".\n", sep = "")
note <- .hessian_note(x$hessian.rank)
if (!is.null(note)) cat(note, "\n", sep = "")
invisible(x)
}

#' Bind a list of data frames
#'
#' @param x a list of data frames
#'
#' @return A data frame
#'
#' @export
dfbind <- function(x) {
  if (!is.list(x) || !all(vapply(x, is.data.frame, logical(1))))
    stop("'x' must be a list of data frames.", call. = FALSE)
  if (!length(x)) return(data.frame())
  out <- do.call(rbind, unname(x))
  for (nm in names(out)) {
    if (is.matrix(out[[nm]]) && inherits(x[[1L]][[nm]], "AsIs"))
      out[[nm]] <- I(out[[nm]])
  }
  out
}

#' Moore-Penrose pseudo-inverse of a matrix
#'
#' @param x a matrix
#' @param tol a scalar tolerance. For \code{pinv()}, a negative value (the
#'   default) uses Armadillo's default tolerance and a non-negative value is an
#'   absolute cut-off on the singular values; for \code{ginv.egpd()} it is a
#'   relative tolerance (see Details).
#'
#' @return A matrix
#'
#' @export
pinv <- function(x, tol=-1) {
armapinv(x, tol)
}

#' @rdname pinv
#' @details \code{ginv.egpd()} follows \code{MASS::ginv()}: \code{tol} is a
#'   \emph{relative} tolerance, and singular values \code{d} with
#'   \code{d <= tol * max(d)} are treated as zero, so the result does not
#'   depend on the scale of \code{x}. \code{tol} must be a single non-negative
#'   number. (\code{pinv()} with \code{tol >= 0} uses \code{tol} as an
#'   absolute cut-off.)
#' @export
ginv.egpd <- function(x, tol=sqrt(.Machine$double.eps)) {
if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
  stop("'tol' must be a single non-negative number")
armaginv(as.matrix(x), tol)
}

#' Generate a sequence between a range
#'
#' @param x a 2-vector
#' @param length an integer
#'
#' @return A vector
#'
#' @export
seq_between <- function(x, length=NULL) {
if (is.null(length)) {
    return(seq(x[1], x[2]))
    } else {
    return(seq(x[1], x[2], length.out=length))
    }
}

#' Number of observations
#' @param object an egpd object
#' @param ... not used
#' @return integer
#' @export
nobs.egpd <- function(object, ...) {
object$nobs
}
