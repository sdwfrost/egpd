# egpd 0.1.0 (unreleased)

## Bug fixes

### Model fitting (`egpd()`)

* Native identity-link DEGPD fits retain the probability in a bounded last
  cell even if the continuous endpoint rounds to an integer. Refit affected
  models: the previous likelihood could silently assign this mass to the
  preceding cell.

* Transformed responses and fixed-parameter offsets no longer shadow data
  columns. A list field in a formula such as `I(x * cfg$gain)` is not an
  extra predictor. Refit models affected by these name collisions or by
  omission of rows with missing values in the unrelated column.

* Native fits use complete cases across the evaluated predictor terms
  and offsets. For example, an inline factor transformation that creates
  missing values now gives the same fit as the precomputed column.
* Compact fits evaluate designs and offsets on the fitting sample before
  grouping. A centered exposure is no longer recalculated on just the
  distinct rows. Refit affected compact models.
* `dfbind()` preserves matrix-valued columns and time classes when
  binding observations, without recycling or duplicating response rows.
* `compact = TRUE` now groups exact predictor values, including the
  components of matrix columns made with `I()`. Previously, distinct
  timestamps could be grouped together, changing coefficients and even
  the reported likelihood. Refit affected compact models.
* Multiple offsets on one predictor are summed in the likelihood and in
  prediction, including compact fits and offset-only parameter blocks.
  Previously only the first exposure could enter the model. Refit models
  that used separate offsets on the same predictor.
* Automatic starts now account for offsets on free intercepts. Adding a
  constant exposure offset to the mean of a Poisson-lognormal fit could
  previously trap the dispersion at its lower bound and return a much
  poorer likelihood as converged. A median-centred fallback keeps a single
  outlying offset from spoiling an otherwise finite automatic start.
  Explicit `inits` remain coefficients.
* Sandwich adjustments for score-outer-product (BHHH) families now use
  observed likelihood curvature for inference. A distinct cluster per
  observation no longer forces an identity adjustment. If the required
  curvature cannot be resolved, an error replaces spurious robust SEs.
* An unpenalised `fx = TRUE` smooth keeps its fitted basis for prediction,
  intervals, simulation, summaries and plots even if a different predictor
  has penalties. Previously new-data predictions could revert to the
  intercept. Refit affected models to restore their missing basis metadata.
  `summary()` also handles a smooth with a single coefficient.
* Custom response-scale `fixed` values use the predictor's actual `unlink`
  (or its documented prefix), not a guessed log link. An unverifiable
  inverse gives an error rather than fitting a different value. `?egpd`
  now shows how custom likelihood callbacks expand compact design rows and
  apply formula offsets and pins, and distinguishes mgcv's absolute REML
  score conventions under `gamma != 1`.
* `egpd()` honours smoothing parameters fixed in a term (`s(x, sp = 5)`)
  and parameters shared by terms or factor levels with the same `id`
  within a predictor. Previously a fixed term could be left unpenalised or
  re-estimated, and shared penalties could be estimated separately. `sp`
  reports the independent multipliers not fixed by the formula; `full.sp`
  reports the multiplier of every penalty, including fixed and shared ones.
* `family = "ztgwaring"` (and therefore the default `fit_hurdle()`) now starts
  inside the parameter space and reaches the MLE where it is interior, also
  from `egpd()`'s default start: when that start's intercept-only Newton
  steps would leave the interior (for example run into the lower bound of
  rho) while the intercept-only zero-truncated MLE lies inside it, the fit
  starts from that MLE (one such fit stopped at the bound with convergence 3,
  4.27 log-likelihood units below it). Where the likelihood rises towards
  rho = 1, the fit still stops at the bound with convergence 3. Fits that
  never leave an invalid start, or whose log-likelihood is not finite, are
  reported as failures instead of as converged.
* The inner Newton optimiser no longer stalls when step clipping reverses its
  descent direction. Fits whose Newton iterations stall or run out next to the
  optimum -- common with the outer-product Hessians of the generalized Waring,
  zero-truncated and other count families -- are checked by a short
  gradient-only search before being reported as not converged.
* DEGPD restarts also run when the response formula is not listed first, and
  formulae keep their environments, so smooth arguments such as
  `s(x, k = kk)` may use local variables in every slot.
* Default restarts fit one data set: `data` is evaluated once (each start
  re-evaluated it in the caller's frame), and every start draws the same
  `maxdata` and `maxspline` subsamples (each drew its own, and the highest
  log-likelihood then picked the easiest subsample). Starts that differ only
  in parameters held by `fixed` are fitted once, and a shape held by `fixed`
  under a finite `xi.max` no longer triggers two wasted refits. When the best
  DEGPD fit has its shape on a flat end of the bounded (`xi.max`) link, the
  extra start taken from a log-link fit also works with `removeData = TRUE`
  (it failed).
* Formulae accept expressions of data columns in parametric terms (`log(x)`,
  `factor(g)`, `I(x^2)`, `poly(x, 2)`), in smooth terms (`s(log(x))`,
  `te(log(x), z)`, `s(x, by = factor(g))`, the factor of an `"fs"` smooth)
  and in the response (`I(y + 1) ~ x`), also in `predict(newdata)`; they were
  rejected as variables missing from `data`, or failed with "undefined
  columns selected". As in base R and mgcv, an expression is always
  evaluated, even if `data` has a column named like it: only a name refers to
  a column. A response expression is evaluated once, in the environment of
  the formula that carries it, so a local constant in it is found whatever
  the environments of the other formulae; as in mgcv, a local constant in a
  term (`poly(x, deg)`) is found when fitting but not by `predict(newdata)`.
* `egpd(family = "custom")` accepts a bare formula, for a one-parameter family
  (it set up three parameters and failed with an internal error).
* `egpd(maxspline = )` with fewer rows than the data no longer fails, also
  when its subsample misses a level of a factor in the model, including a
  factor expression in a smooth such as `s(x, by = factor(g))` (the subsample
  is extended by one row of each missing level), and `maxdata`, `maxspline`
  and `compact` work with a data frame holding only the response.
* Compact fits (`compact = TRUE`) of the continuous EGPD (m = 1, 3) with
  smooths no longer fail or use a wrong REML gradient: compact fits of every
  family use finite-difference REML derivatives.
* `t2()` smooths now use mgcv's map from their fit basis to the prediction
  basis (`G$P`), as `mgcv::gam()` does. Predictions, simulations and Wald
  intervals on `newdata`, `compact = TRUE` and `maxspline` fits, and `plot()`
  combined prediction-basis designs with fit-basis coefficients or penalties
  and were off by O(1) (a compact fit's `negREML` was 282.84, against mgcv's
  268.49). `summary()` reports its tests in the prediction basis, as
  `summary.gam()` does (an intercept was 0.076, against mgcv's 0.130). The
  coefficients, their covariance and `predict(type = "lpmatrix")` stay in the
  fit basis.
* `logLik()`, `AIC()` and `BIC()` report the unscaled log-likelihood under a
  sandwich adjustment or `gamma != 1` (the adjusted value is
  `fit$logLik.adjusted`).
* The smoothing-parameter covariance behind `Vc` is computed by
  Richardson-extrapolated finite differences for families without analytic
  third derivatives, judges each curvature against its own error (small but
  well-determined curvatures were dropped, understating some smooth-term
  standard errors by up to 2.6 times), and is never inverted when not
  positive definite.
* The analytic smoothing-parameter searches (`outer = "bfgs"`, the default,
  and `"newton"`) no longer stop silently on the edge of the near-zero-shape
  band where the analytic DEGPD third derivatives are unavailable, reporting
  convergence with smoothing parameters up to 7 times too large: they warn and
  retry with finite differences, as documented.
* Under `sandwich.args = list(method = "curvature")`, the analytic REML
  gradient and the REML Hessian behind `Vc` are correct (the gradient was 8.7%
  off, and the smoothing parameters a few percent).
* The penalty log-determinant `log|S|+` is stable for smooths with several
  penalties (`te()`, `ti()`, `t2()`, `select = TRUE`) however far apart their
  smoothing parameters are (it dropped directions at log-scale spreads of
  about 25, and with three or more penalties failed there), and is taken from
  the right smooth when an unpenalised smooth (`fx = TRUE`) precedes a
  penalised one, whose smoothing parameter was estimated with the wrong
  penalty rank.
* A failed inner refit at a trial smoothing parameter no longer aborts the
  fit; finite-difference REML no longer overrides a user-supplied inner
  `gradtol` and avoids redundant inner refits (about twice as fast at
  n = 2e4).
* `sandwich.args = list(force = TRUE)` works when the Hessian or sandwich
  variance has to be pseudo-inverted (it failed with an internal error).
* `sandwich.args` works with `fixed`: the adjustment uses the free
  coefficients only (it failed, and with `force = TRUE` the magnitude factor
  counted the pinned coefficients and shrank the standard errors).
* A parameter held by `fixed` stays pinned in the intercept-only start fit,
  so its inert coefficient is 0 and `summary()` reports it as `NA` (it showed
  an arbitrary estimate with a zero standard error); a `fixed` value outside a
  count family's evaluable range is an error naming `fixed`, and fixing every
  parameter is a clear error instead of an internal one.
* The ZIDEGPD warm start no longer counts offsets twice.

* `egpd(family = "comppareto")` screens competing splice points, including
  observations and points between adjacent observations. Up to two promising
  splits are refitted; fits stranded with all observations on one side also
  try four quantile starts. With predictors the screen searches shifts of
  the splice intercept, not all possible fits. `restarts = FALSE` and
  explicit `inits` skip the search.
* An exponential-body fit that stops at an observed splice can now converge
  when the smooth parameters are stationary and neither one-sided move
  improves the objective. A zero-rate boundary is not reported as a finite
  optimum. When the selected split has no finite tail-shape maximum, it
  reports that diagnosis once, not warnings about an earlier optimizer stop
  or a singular Hessian followed by a second diagnosis.
  The exponential-body score and curvature now use the same stable tail
  log-ratio as the density when an observation lies just above the splice.

### Analytic derivatives (C++)

* Bell and zero-inflated Bell fits compute curvature directly on the link
  scale. An offset of -400 on one row no longer makes an otherwise estimable
  fit fail with a nonfinite Hessian.
* Continuous EGPD carrier 3 retains its third and fourth derivatives at
  small `delta` by expanding the likelihood in `delta` times the GPD
  cumulative hazard. This fixes the REML gradient even when `delta` is held fixed.
* Fixed the analytic third derivative d^3/(d log(sigma) d xi^2) of the
  continuous EGPD log-likelihood (carriers 1-4; also a fourth derivative of
  carrier 4), which made the default REML gradient of smoothed
  `family = "egpd"` fits inexact.
* Continuous EGPD model 6 (truncated beta): the log-kappa gradient and Hessian
  are exact for large kappa (they were wrong, even in sign, above about 30),
  and a non-finite kappa no longer hangs the normaliser quadrature.
* Continuous EGPD carriers 5 and 6 keep their normalising constant as
  kappa -> 0 (the carrier-6 density integrated to 1.32 at kappa of about
  1e-16, and both carriers returned the 1e20 penalty for smaller kappa), and
  the gradients and Hessians of carriers 3 and 4 are exact for very small
  delta (the log-delta score was 0.36 instead of 0 at log(delta) = -36).
* Continuous EGPD derivatives next to a negative-shape upper endpoint no
  longer pick up the out-of-support penalty (gradients were off by up to
  1e23); observations far below the scale (y/sigma < 1e-4), shapes near zero
  and the far upper tail use exact automatic differentiation.
* The analytic third and fourth derivatives behind the REML gradient are
  accurate in DEGPD (carrier 1) upper-tail cells of very small probability,
  such as the identity-link top attainable bin and far-tail cells (they were
  up to 160% off), and for continuous carrier 3 next to a negative-shape
  upper endpoint and far in the tail, also for small delta (fourth
  derivatives were 45-160% off, and for small delta up to 129 times their
  size, or zero).
* Zero-inflated DEGPD: the y = 0 log-likelihood keeps full precision when the
  zero probability is close to one, is never slightly negative, and is never
  `NaN` when P(0) rounds to one.
* DEGPD and zero-inflated DEGPD likelihoods, gradients and Hessians are
  accurate at very large counts (y >= 65536), where the objective was 5e-3
  off at y = 1e12 and returned the 1e20 penalty once y + 1 == y; this
  includes identity-link gradients and Hessians next to a negative-shape
  upper endpoint.
* DEGPD with a finite `xi.max`: a sign error in the series used for very
  small shapes (below about 1e-4 `xi.max`) made the shape score up to 1.7e-9
  too small (relative).
* DEGPD/ZIDEGPD objective-only evaluations are 2-12 times faster again.

* Continuous EGPD carrier 5, DEGPD/ZIDEGPD carriers 1, 2 and 4 (and carrier
  5 at counts >= 65536) and `family = "ztdegpd"`: likelihoods at log kappa (or
  a log-link log xi) below -708.4, where kappa is subnormal, are computed from
  the coefficient itself and equal their limits. They could exceed the
  maximum (by 21 for an EGPD carrier-5 fit and 55 for a zero-truncated shape)
  or hit the 1e20 penalty. The zero-truncated Hessian is no longer zero for
  log kappa below about -355.

### EGPD distribution functions

* `dcomppareto()` preserves finite tail densities when the sum of the
  observation and splice point overflows, or when the conditional tail
  density overflows before being multiplied by the tail weight.

* Type-4 and type-5 lower tails retain the log of the standardized input
  through the carrier transformation. In particular, type 5 with `delta = 1`
  agrees with the equivalent power carrier in densities, quantiles and draws.

* Bounded count probabilities and quantiles retain a last cell when the
  continuous endpoint rounds to an integer just below its true value.
  Discrete draws no longer acquire mass at an exactly integral open endpoint.
  Power carriers retain representable probabilities, densities and scaled
  draws even when the intermediate base CDF underflows.

* `qbgpd()` inverts narrow blends on the scale of their knots. Blended GP
  and composite-Pareto tail quantiles no longer overflow just because an
  intermediate exponential or doubled scale exceeds the double range.

* GPD CDFs, densities and quantiles retain the exponential limit for
  subnormal shapes, and representable quantiles no longer overflow during
  rescaling. Standalone type-4 and type-5 densities remain finite in the
  interior when `delta` is tiny.
* Discrete and spliced GPD quantiles respect exact survivor jumps, without
  shifting neighboring probabilities by a tolerance. Compound query `NaN`
  values remain `NaN`, distinct from `NA`.

* `pegpd()` now retains the type-3 beta concentration transition at large
  GPD hazards and the type-4/5 lower tail when `delta` is very large.
  Previously these combinations could return probabilities above one.

* Models 3 and 4 retain the lower quantile bracket for large `delta`,
  including count quantiles and inverse draws. Model 6 densities, cell
  probabilities and shape scores now cancel the common beta normalizer
  before evaluation. Discrete GAMLSS scores use the density-weighted
  score within a cell when subtracting its endpoint derivatives loses
  precision.

* For highly concentrated beta carriers (`type = 3`, GAMLSS `EGPD6` and
  `DEGPD6`), quantiles and random draws now use the beta inverse outside
  the narrow range where a local expansion is accurate. With very small
  positive powers, the `type = 6` mixture now retains its upper tail and
  discrete cell masses instead of rounding them to zero during inversion.
* Composite Pareto CDFs retain the body mass above the splice at small
  Pareto shapes. The CDF no longer drops at the splice and the inverse
  no longer returns negative values for these positive-support models.

* Blended GP functions retain the exact nonzero tail shape on both sides
  of zero. The bridge, endpoint, density and quantiles use the scale at
  the lower matching knot instead of switching small shapes to exponential.
* For carriers 4 and 5 (fitting models 3 and 4), small positive `delta` no
  longer collapses an ordinary quantile or inverse draw. The inverse uses
  the smaller beta tail and retains the GPD log survivor. For example, a
  discrete median reported as 14 is now 18, the first count whose CDF reaches 0.5.
* `pegpd()`, `pdiscegpd()`, `ddiscegpd()`, `degpd_density()`, `qegpd()`,
  `qdiscegpd()` and the zero-inflated variants evaluate `log`/`log.p` results
  on the log scale (far-tail values underflowed to `-Inf`/`Inf`), including
  where intermediate values are subnormal, where `q / sigma` overflows, and
  for lower-tail log-probabilities below about -745 (the discrete and
  zero-inflated quantiles returned 0). Zero-inflated log-probabilities never
  exceed 0 (they rounded to +7e-18 where the probability is 1, and the
  `log.p` round trip gave `NaN`).
* `pzidiscegpd()` and `pziegpd()` with `log.p = TRUE` keep full relative
  accuracy close to 0, where the zero-inflated survivor is small, so
  `qzidiscegpd(log.p = TRUE)` at levels within about 1e-13 of 0 is the
  generalised inverse (it was 3 counts too small; see Behaviour changes).
* Lower-tail quantiles of carriers 2, 3 and 6 are accurate for very small
  probabilities; far-tail densities of carriers 4/5 with small `delta` are
  accurate; cancellation was removed from the type-5 upper tail and the
  type-3 lower tail near zero.
* `qegpd()` and `q.G()` are accurate for small kappa: carrier 1 with p > 1/2
  (`qegpd(0.9, kappa = 0.001)` was 0, not 1.7e-46), carrier 5 where
  p^(2/kappa) underflows (quantiles such as 1e-200 were 0), and carrier 3 as
  kappa -> 0, whose normaliser cancelled (`qegpd(pegpd(0.5))` was 673.8 at
  kappa = 1e-15).
* `p.G()` is accurate near u = 0 for carriers 2 and 3 (all digits were lost
  at u = 1e-20) and where its carrier-5 intermediate value underflows, and,
  with `d.G()`, for carrier 3 as kappa -> 0. Carrier 6: `d.G()` and
  `degpd_density()` with `prob` 0 or 1 are no longer `NaN` at zero, and
  `d.G(log = TRUE)` no longer underflows to `-Inf` at tiny u.
* `regpd()`, `rdiscegpd()`, `rzidiscegpd()`, `rziegpd()` and `r.G()` return
  `numeric(0)` for zero-length `n`; `rzidiscegpd(unifsamp = )` is
  reproducible when `pi > 0` and returns exactly `n` draws.
* `regpd()` and `rdiscegpd()` with a huge carrier shape (kappa above about
  1e4) take the GPD quantile from the log survival probability, so draws stay
  accurate and inside the support (they were `Inf` from kappa of about 1e14);
  seeded draws for ordinary parameters are unchanged.

### GAMLSS and BAMLSS families

* Held-out plug-in coverage now predicts formula-based GAMLSS smooths
  such as `lo(~x)` without mistaking their zero placeholders for aliased
  predictors. Real aliases are still checked for estimability.
* Held-out plug-in predictions for backfitted GAMLSS smooths no longer
  confuse a user's `source` column with a training/prediction row marker.
* BAMLSS continuous and floored EGPD means account for concentrated
  mass near a finite endpoint. Zero-inflated means apply the mixture
  weight before an unweighted expectation can overflow.
* New GAMLSS fits evaluate quantile residuals in log tails, including
  randomization within the observed count's probability jump. A finite
  residual is no longer replaced by infinity or a clipped tail value.
  Refit older objects to replace their stored residual expression.
* `GPIG()` and `ZIGPIG()` differentiate the exact count law. Small means
  and tempering parameters near one no longer bias the GAMLSS scores or
  move an already optimal coefficient away from the likelihood maximum.
* BAMLSS quantile residuals randomize over the actual CDF jump, including
  zeros in mixed distributions and counts in held-out samples with no
  observed zero. Response residuals subtract the distributional mean,
  including zero inflation, rather than a transformed scale or component
  mean. Infinite means are flagged rather than used as finite predictions.
* GAMLSS fits can use `family = egpd::BELL()`, `egpd::DEGPD1()` and the
  other egpd families without `library(egpd)`, including when regenerating
  quantile residuals. Separate GAMLSS helpers that look up densities by
  name, including `vcov()`, still require `library(egpd)`.
* `BELL()` and `ZIBELL()` use exact scores instead of finite differences.
  At very small means the old score for a positive count could shift an
  intercept estimate even when GAMLSS reported convergence.
* The GAMLSS EGPD/DEGPD/ZIDEGPD/ZIEGPD families use analytic, numerically
  stable scores. Zero counts with kappa < 1 (and DEGPD4 with kappa < 2) had
  zero mu/xi scores, which collapsed the fitted shape to xi = 0 on
  zero-heavy data.
* `DEGPD5()`/`ZIDEGPD5()` kappa scores no longer underflow to 0 for large
  `nu` or lose accuracy when the scale is much larger than the count, and
  `EGPD4()$dldv` at y = 0 with `tau = 1` is the closed-form limit (it was
  `NaN`, with a false "zero likelihood" warning).
* The numeric score helpers handle parameters at domain boundaries (for
  example ZIBell with `pi -> 0`) instead of returning 0.
* Starting values: `GPIG()` and `ZIGPIG()` start `gamlss()` from the
  intercept-only maximum likelihood estimate, as the other families do (from
  their fixed starts the RS algorithm drifted to a -> 0 on heavy-tailed data
  and stopped 100 or more log-likelihood units below the maximum), and the
  EGPD families start inside (0, 1) for a `"logit"` `nu` or `tau` link (they
  stopped with "Value ... out of range").
* The intercept-only starts of the EGPD and GPIG families are reused only for
  an identical response (a response sharing the summary statistics of the
  previous one reused its start: a `DEGPD1(sigma.link = "identity")` fit then
  failed, and a `GPIG()` fit converged 8.7 log-likelihood units below the
  maximum), and the EGPD start clamps the GPD scale relative to the scale of
  the data and is always checked on the full data (bounded-tail data on a
  scale of about 1e8 with `sigma.link = "identity"` started at zero
  likelihood, and `gamlss()` stopped).
* `bell_bamlss()`, `zibell_bamlss()`, `gpig_bamlss()` and `zigpig_bamlss()`
  pass `lower.tail`, `log.p` and `max.value` on in `$p()` and `$q()` (they
  were silently ignored).

### Count families

* Count quantiles now invert ordinary upper probabilities near one using
  their exact lower-tail complements. Previously a Waring quantile of 13
  could be reported as 0, even when the requested search cap was too small.

* Generalized Poisson upper probabilities and quantiles near unit dispersion
  no longer lose a small survivor to the near-unit mass at zero.

* Extreme lower-tail quantiles now resolve the small upper probability.
  For instance, `qpig(1 - 1e-14, 5, 50)` now gives 11,190 rather than
  10,773, agreeing with inversion on the survivor scale.
* Mean-parameterized GPIG and ZIGPIG distributions now recycle inputs
  together before computing the native level parameter. In vectorized
  `d/p/q/r` calls this keeps observations and draws paired with their
  requested means even when parameter vectors repeat at different periods.
* CMP and Poisson-lognormal upper tails no longer drop to zero at count
  zero just because the zero mass rounded to one. The Poisson-lognormal
  evaluates the small positive tail directly when computing a near-unit
  zero mass and its scores, including at large log-intensity dispersion.
* GPIG lower-tail CDFs stay in `[0,1]`. Near-one GPIG/ZIGPIG quantiles use
  the survivor when it can be resolved and retain the forward log CDF when
  it cannot. A quantile of 794 was incorrectly reported as above the cap
  because a forward sum plateaued below one. At `a = 1/2` the exact PIG
  recurrence avoids quadratic work for long finite-mean ladders.
* GPIG/ZIGPIG: the default "hybrid" likelihood is continuous in the
  parameters (it jumped where the saddlepoint correction was switched off)
  and faster (one cached recursion ladder per parameter row); it is exact
  wherever the saddlepoint has not converged up to `ymax` (default 1000, as
  before), and at every count with `gpig_control(ymax = Inf)`, which removes
  the remaining bias above `ymax` when `c` is close to 1 at a higher cost.
  `gpig_control(normalize = TRUE)` is 40-100 times faster and its score is
  analytic.
* GPIG/ZIGPIG: the zero mass is formed from expm1(a log1p(-c)) where
  (1 - c)^a - 1 cancels, so the pmf is accurate for a tiny `a` (at a = 1e-20,
  b = 1e20 and c = 0.5 it was doubled) and near the Poisson limit of a small
  `c` (P(0) was 2.7e-4 off, relative, at c = 1e-12).
* Generalized Waring: the log-likelihood, its score and `pgwaring()` tails stay
  accurate up to the negative-binomial and Poisson limits (fits on Poisson
  data no longer report positive log-likelihoods), and the fitted likelihood
  equals `dgwaring()` at the reported parameters.
* Generalized Waring: `pgwaring()`/`qgwaring()` upper tails, and the hurdle
  quantiles built on them, are correct when a = mu (rho - 1) / k is tiny (below
  about 1e-16 rho: a tiny mean, rho near 1 or a large k). P(Y > q) used to be
  returned as 1, giving positive log-probabilities, `-Inf` lower-tail
  log-probabilities and `NA` quantiles.
* Generalized Waring: `pgwaring()` and `qgwaring()` are accurate, in both
  tails, for a negative-binomial shape k beyond 1e12, up to the fitting bound
  1e15 and beyond: P(Y > q) was up to 3e-5 off (relative) at k = 1e15 and
  `NaN` from k of about 1e18, and an evaluation took seconds (now
  milliseconds).
* Bell/ZIBell: the log(B_n/n!) table is cached for the session and built in
  O(n) (likelihoods with counts near 1e4 went from seconds to milliseconds).
* CMP: the GAM score is the exact derivative of the likelihood, whose
  interpolation grid is checked against exact solves and replaced by the
  exact sums wherever its error could exceed 1e-5. For strongly
  underdispersed laws (nu above about 10, variance below 1) the grid was
  unnormalised or infinite, giving wrong estimates and log-likelihoods (11
  units too high at nu = 1000); where it validates (everywhere for nu <= 8)
  results are unchanged.
* CMP: `dcmp()`, `pcmp()`, `qcmp()`, `rcmp()` and the likelihood's exact sums
  and grid nodes solve the centred mean equation for log lambda where the
  rounded mean cannot resolve the root, so strongly underdispersed laws with
  an integer mean (e.g. mu = 1, nu = 128) are right: the rounded mean equals
  mu over a whole range of log lambda there, so the root was arbitrary, tail
  log-probabilities were off by 7.6 nats at nu = 128 and 310 at nu = 1000,
  and quantiles by a count. Where the variance underflows (an integer mean
  with nu above about 2000), the distribution functions stop with an error
  instead of returning a wrong value.
* Poisson-lognormal probabilities are accurate to about 1e-13 (they were off by
  up to 2e-3 at y = 0 for large sigma).
* Count distributions: upper tails (`lower.tail = FALSE`) and `log.p` keep full
  relative accuracy (they were computed as 1 - cumsum(pmf)), and quantiles are
  exact generalised inverses, on the survivor scale for upper tails; a
  survivor that cannot be resolved affects only the upper tail, and `pprl()`
  and `qprl()` evaluate the closed-form survivor at the requested counts
  only.
* Count quantiles with `log.p = TRUE` are exact generalised inverses also at
  levels within about 1e-13 of 0: the search shifted every log-scale level by
  a fixed 1.4e-14, so it returned the quantile of a tail probability of that
  size (`qgpois(-1e-20, 3, 0, log.p = TRUE)` was 24, while `qpois()` gives
  30); above log p = -1 its tolerance is now relative to 1 - p. Upper-tail
  log-probabilities keep full relative accuracy where the lower tail is small
  (`pgpois(0:3, 50, 0, lower.tail = FALSE, log.p = TRUE)` was 0). See
  Behaviour changes.
* Zero-inflated count log-probabilities (`dzigpig()`, `dzipig()` and
  `dzibell()` at 0, and the zero-inflated log distribution functions) never
  exceed 0 (they rounded to +6.9e-18, and the `log.p` quantile round trip then
  gave `NaN`), and log P(0) keeps its relative accuracy when it is close to 0.
* `dpig()`, `qpig()` and the ZIPIG functions are computed natively:
  `dpig(-1)` is 0, `NA` inputs propagate, and `qpig()` returns `NA` with a
  warning when `max.value` is not attained.
* PIG/ZIPIG: the pmf recursion is scaled where it under- or overflowed:
  `dpig()` lost accuracy for a mean below about 1e-154 and was `NaN` below
  about 2e-162, where `ppig()` and `qpig()` failed; P(0) was 1, and the other
  probabilities `NaN`, for sigma mu above about 9e307; and the likelihood's
  Hessian was `NaN` for log mu <= -237.
* CompPareto: the likelihood treats y <= 0 as having zero density with a
  consistent score, and `egpd(family = "comppareto")` rejects non-positive
  or non-finite responses.
* CompPareto: `dcomppareto()`, `pcomppareto()`, `qcomppareto()` and the
  fitted likelihood form the splice weight on the log scale where g(theta) or
  G(theta) underflow (they returned `NaN`, and the GAM objective `Inf`).
* CompPareto: below the splice point, `pcomppareto()` takes the upper tail,
  and `log.p` values near 0, from the body distribution's own survivor where
  the survivor is below about 1e-3, and `qcomppareto()` inverts upper-tail
  and `log.p` levels there on it, so both keep full relative accuracy, stay
  monotone and invert each other when the Pareto tail weight is tiny (the
  upper tail was 1 - F, which lost all relative accuracy from a tail weight
  of about 1e-10 down: 118% off at a weight of 3.8e-19, and `-Inf` on the log
  scale at a survivor of e^-787). See Behaviour changes.
* Count and CompPareto random generators return a zero-length vector for a
  zero-length `n`, as base R does (they stopped with an error), and
  `rpoislnorm()` returns doubles instead of `NA` when a draw exceeds the
  integer range (integers otherwise, as `rpois()`).

### Prediction and inference

* Smooth summaries use the prediction design even when a tensor-product and
  random-effect term make the basis mapping singular. Predictor expressions
  and offsets may not expand held-out rows through vector recycling.

* BAMLSS coverage uses the saved fitting offsets, not formula offsets
  ignored by the fitter. Supply varying offsets at `newdata` with `offset`;
  in-sample intervals use the saved values.

* Held-out BAMLSS prediction intervals now use the fitted polynomial basis
  and factor contrasts rather than rebuilding them from the prediction batch.
  Native predictions accept data frames created by `model.frame()` for another
  predictor without treating their terms as the fitted design.

* Shortest empirical intervals count sorted draws and use their exact
  numeric values. They no longer pick an extra observation through
  accumulated rounding or merge distinct endpoints into printed labels.
* Plug-in GAMLSS prediction preserves a predictor named `source` when
  selected by a ridge smoother's `x.vars` argument.

* A missing numeric covariate or offset in one prediction row no longer
  removes the offsets from complete rows in predictions, Wald intervals
  or simulation draws.
* `predictive_coverage(use_shortest = TRUE)` also returns a shortest
  empirical interval for continuous and mixed responses, not just counts.
* GPIG and ZIGPIG profile intervals now use the likelihood evaluator and
  numerical settings that produced the fit. Older fits without these
  settings must be refitted before profiling.
* Response-scale Wald limits use the same count-family bounds as point
  predictions. Delta-method SEs are zero on a flat part of the mean map
  and undefined at a kink with positive link-scale uncertainty. For
  example, a PRL interval can no longer lie below its minimum mean.
* Held-out GAMLSS intervals use the fitted spline knots and factor
  contrasts, and sum all formula offsets. Adding unrelated prediction
  rows no longer changes the fitted spline at the original rows.
* `pinv(x, tol = 0)` now keeps all positive computed singular values,
  even on multiscale matrices. Omitting `tol` still uses the automatic
  tolerance.

* Held-out GAMLSS intervals now calibrate likelihood-curvature steps to
  the model rather than to fixed response units. Both plug-in and
  parametric prediction reject rows that are not estimable from aliased
  training columns. Plug-in prediction does not require a covariance matrix.
* `summary()` now tests fully penalised terms such as random intercepts
  against their boundary null distribution. If the fit cannot support
  that test, it returns `NA` and explains why instead of reporting an
  ordinary smooth-term p-value.
* Clustered magnitude adjustments now warn about rank-deficient scores for
  exact-Hessian families as well as score-product families. This matters
  when the model includes cluster fixed effects: a scalar adjustment can
  substantially understate individual SEs. The default curvature method
  rejects a singular adjustment.

* At `newdata`, parametric GAMLSS prediction intervals use the fitted
  weighted likelihood and its free coefficients, holding fixed parameters
  at their actual values. Changing the link on a fixed shape no longer
  changes the interval. This path does not require `vcov.gamlss()` or
  `library(egpd)`; it still requires identifiable parametric predictors.
* Count-family profile intervals refit free nuisance intercepts within
  their applicable fitting bounds and compare the edge and interior
  branches. A flat score beyond a clamp previously cut GPIG and Waring
  intervals short without a profile warning.
* GAMLSS `predictive_coverage()` supports GAMLSS-specific inverse links
  such as the default degrees-of-freedom link in `TF2()`. Plug-in coverage
  works from fitted response-scale parameters without reconstructing links.
* For grouped binomial GAMLSS fits, both methods of `predictive_coverage()`
  use the fitted number of trials per row instead of simulating a single
  trial. With `newdata`, supply `bd` (a common size or one per row) and use
  success counts for `y`. Ungrouped Bernoulli predictions still default to
  one trial.
* Bounded DEGPD and zero-inflated DEGPD quantiles at probability one now
  return the largest attainable count rather than an earlier count whose
  CDF rounded to one. `rqresid()` uses log tail probabilities so that a
  finite extreme residual isn't lost when the lower CDF rounds to one.
* `confint()` rejects invalid confidence levels for `distfit` models. For
  `egpd` models a 100% profile interval reaches the link boundaries; a
  fixed parameter keeps its point interval, including at `newdata`.
* `predictive_coverage()` requires one observed response per fitted row
  (or per row of `newdata`). The GAMLSS and BAMLSS methods no longer recycle
  in-sample predictions when the response vector has the wrong length.
* `summary()` of `egpd` fits reports two-sided Wald p-values for parametric
  terms (they were one-sided, i.e. halved).
* `confint()` reports slopes and smooth coefficients on the link scale under
  their coefficient names (it passed them through the parameter's inverse
  link), no-intercept predictors no longer silently use their first
  coefficient, and profile intervals that reach the boundary of the parameter
  space report that boundary instead of `NA`.
* `confint(newdata = )` gives `NA` limits where a coefficient has no positive
  variance, as `confint()` does (it gave a zero-width interval).
* `confint()` for `egpd` fits rejects out-of-range confidence levels and
  fractional coefficient indices. A single slope without an intercept no
  longer qualifies for an intercept-only profile interval.
* `confint(method = "profile")` brackets its limits in steps of
  0.1 (|estimate| + 1) when the Hessian behind the standard errors is singular
  or a coefficient's variance is negligible: generalized Waring fits at a = k
  gave k in [0, Inf] (now [1.063, 8.519], matching an independent profile
  likelihood); full-rank fits are unchanged.
* `confint(method = "profile")` no longer cuts limits short where the fit of
  the other coefficients at a profile point stalls or ends in another basin
  (e.g. an EGPD model-3 upper scale limit 1.96 -> 2.385, and fits with xi on
  its boundary): each finite limit is checked by refitting from other starts,
  and the search follows the profile path when the check fails. Model 4 also
  checks the small- and large-delta branches: a 95% lower scale limit was
  0.763 where the likelihood-ratio statistic was only 2.33; the corrected
  limit of 0.482 gives 3.84.
* `predictive_coverage()`: bamlss fits propagate the posterior samples (an
  arbitrary standard error of 5% of the mean linear predictor was used);
  gamlss fits return integer intervals for every discrete family and, with
  `method = "parametric"` and `newdata`, include parameter uncertainty
  through standard errors from `vcov()` and the new-data design (they were
  silently plug-in intervals); and fits made with `removeData = TRUE` explain
  that `newdata` is needed.
* `predictive_coverage()` works for gamlss fits with a fixed parameter
  (`<p>.fix = TRUE`; it stopped with "Could not determine the 'sigma'
  link"): the parameter is held at its fitted value (at `newdata` it must be
  a single value), and the other parameters' standard errors at `newdata`
  are conditional on it.
* `simulate(type = "response")` clamps parameters exactly as `predict()` and
  the likelihood do; `predict()` and `rqresid()` work for `ztgwaring` fits.
* For `family = "custom"`, `simulate()` and `confint()` apply
  `custom.fns$unlink`, as `predict()` does (they returned link-scale values as
  parameters, e.g. negative mean draws and `NaN` quantiles), and
  `predict(se.fit = TRUE)` differentiates an `unlink` without a `"deriv"`
  attribute numerically (it returned the link-scale standard error). A
  decreasing `unlink` (e.g. `1/eta` or `exp(-eta)`) gives non-negative
  standard errors (they were negative, also with a `"deriv"` attribute) and
  `confint()` limits in increasing order, also by profile and with
  `newdata`.
* `plot()` of an `egpd` fit without smooth terms is a clear error (it failed
  in `par(mfrow)`), and graphical parameters are restored on exit.
* `ginv.egpd()` uses a relative singular-value tolerance, like `MASS::ginv()`.

### `distfit()`

* For type-1 discrete fits with `xi = -1`, scale estimation now accounts
  for the change in the highest observed cell's mass at the endpoint.
  A bounded fit could previously report convergence short of this maximum.

* `distfit()` and the discrete POT fitter check for feasible improving
  steps when L-BFGS-B reports success. A line search against the support
  boundary could previously return its starting point as a converged fit.
* Continuous type-3 fits no longer report an interior estimate and Wald
  standard errors when the conditional carrier likelihood increases toward
  the excluded limit `kappa = 0`.

* At the uniform limit of the type-2 carrier or a zero- or unit-weight
  type-6 mixture, a free parameter no longer gets a Wald standard error
  determined by how close the optimizer stopped to the endpoint. For a
  two-cell count or jump law with only the mixing weight free, an exact
  fit also protects small interior weights and their standard errors.

* Bernstein likelihoods and diagnostics preserve representable results
  when the unscaled GPD CDF underflows before the power transformation.

* Compound fits withhold Wald uncertainty when the Hessian crosses a
  non-smooth jump-support boundary. Bernstein fits preserve the exponential
  limit of the likelihood for very small nonzero GPD shapes.

* Near a zero-inflation boundary, power-carrier profiles distinguish an
  exact equality from an adjacent interior optimum. Discrete and POT fits
  withhold Wald inference when the Hessian spans a non-smooth support change.

* Bounded count and jump laws are checked against the number of free
  probabilities on their support, not just for the one- and two-cell cases.
  Equivalent continuous uniform and power carriers no longer report a
  regular scale MLE at an unreachable support endpoint.
* Zero inflation is fitted at its exact conditional maximum, respecting
  logit bounds. A fit at zero no longer reports an interior Wald interval,
  including when the observed zero fraction exactly equals the base mass.

* Unbounded default Bernstein fits check for profitable moves on the
  mixture-weight scale before reporting convergence. Refit models where
  a stalled near-zero weight may have hidden a better fit.

* Bernstein fits with optimizer bounds use a box-constrained method
  throughout, instead of polishing the result outside the feasible box.
  Ordinary MLE and Bernstein fits now withhold Wald inference at active
  optimizer bounds, as for the closed-form reductions.
* When the estimated compound jump support has at most two cells, the
  fit warns if its free parameters cannot be identified separately.
  A two-cell jump law identifies only its thinned Poisson mean, even
  if many different aggregate counts are observed.
* Continuous power-carrier fits at fixed shape `xi = -1` no longer
  report ordinary scale uncertainty when the likelihood increases to the
  excluded support endpoint. For discrete fits on one or two cells, the
  ridge check now covers free carrier and scale parameters as well as
  free scale and shape.
* Closed-form Bernoulli and thinned Poisson fits honor optimizer bounds
  on their reduced coordinates (logit probability or log mean), including
  in the likelihood they report. They withhold Wald SEs on an active bound.
* Continuous fits stop rather than reporting finite estimates and Wald
  uncertainty on an unbounded endpoint branch (`xi < -1`). Fixed endpoints
  and discrete fits are unaffected. Bell fits with no positive counts
  also stop: there is no interior maximum on that sample.
* One-parameter fits with `optim.method = "Brent"` retain the parameter
  name throughout optimization. Previously a Bell fit could report a
  zero log likelihood for positive counts.
* If `cpdegpd` has fixed severity support `{0, 1}` with both intensity
  and jump probability free, the aggregate is exactly Poisson. The fit
  reports its mean, standard error and one degree of freedom rather
  than separate estimates for unidentifiable parameters.

* When the fixed endpoint of a type-1 zero-inflated DEGPD permits only
  counts 0 and 1, jointly free `pi` and `kappa` identify just `p1 = P(Y = 1)`.
  The fit reports that probability, its SE and one degree of freedom.
  Other two-cell carriers that still have multiple unidentified parameters
  require additional fixed values; seeing only 0 and 1 in a *sample* does
  not trigger a reduction.

* Type-6 MLE starts use separated powers and the actual mixture weight.
  A collapsed two-power fit is checked for an improving split instead of
  accepting a symmetric saddle as converged.
* Fixed subfamilies report identifiable parameters and information
  criteria. In GPIG and ZIGPIG with `a` fixed at 1 and both `b,c` free,
  `coef()` reports `mu = b*c`, not two arbitrary factors. Equal fixed
  type-6 powers have no estimated mixture weight. Fix one of `b,c` if
  an estimate for the other is required.
* Named scalar values such as `fix.arg = list(pi = coef(pilot)["pi"])`
  now fit the requested likelihood. Previously the inner name could cause
  the parameter to disappear (for example, fitting Bell instead of ZIBell).
* Observations outside the support have log-likelihood `-Inf` (they cost a
  finite 708 each); fits stop when no feasible parameter value is found and
  warn when the optimiser does not converge.
* `distfit()` and its `plot()` use the native PIG/ZIPIG functions for
  `family = "pig"` and `"zipig"`: near-Poisson samples gave absurd `"pig"`
  fits at the sigma -> 0 boundary (log-likelihood +17718, mu = 1e154)
  reported as converged. A positive maximised log-likelihood for a discrete
  family is an error.
* The shape starts at the GPD method-of-moments estimate (its sign was
  reversed), with moments consistent with discrete and zero-inflated data,
  and Nelder-Mead restarts once from its solution.
* `confint()` for `distfit` fits computes Wald intervals on the log/logit
  scale, so they respect the parameter space; `parm` accepts indices, and
  unknown names or indices are an error (they gave an empty matrix).
* `distfit()` checks `fix.arg` and `start`: unnamed, duplicated, non-scalar
  or out-of-range values are errors (an unnamed `fix.arg` was ignored yet
  reported as fixed, and an unnamed `start` was ignored). A probability may be
  fixed at 0 or 1 only where its family allows it: GPIG `a = 0`, and `pi = 1`
  for ZIPIG, ZIGPIG and ZIBELL, are range errors (they failed later with a
  generic error).
* Bernstein-EGPD quantiles keep their relative accuracy in both tails (they
  were `Inf` above about 1 - 1e-9 and inaccurate below 1e-8).
* Bernstein fits (`distfit(method = "bernstein")`) refine Nelder-Mead
  solutions by alternating BFGS, with the exact gradient, and fresh
  Nelder-Mead runs, so they no longer stop up to 1.3 log-likelihood units
  short while reporting convergence, and are faster; the Bernstein likelihood
  and density are evaluated in log space, without clipping tiny observations,
  and stay accurate for a huge carrier shape kappa (above about 1e15), as does
  the distribution function that `plot()` draws (it rounded to 1 there).
* `plot()` of a discrete `distfit` fit draws the fitted pmf over its bars (it
  was drawn at 1, ..., K, drifting left of them).

### Hurdle, spliced and POT models

* POT fits and threshold scans validate all counts before filtering rows
  or computing thresholds. `fix.arg` supports additional parameter pins;
  reported endpoints and stability tables honor fixed scale and shape.

* Inline terms that supply a finite value for a missing input now retain
  their observations in fitting and prediction, like precomputed columns.

* Both parts of a spliced fit now use the same evaluated response, including
  backquoted columns. Refit models if a column name matched the text of a
  response expression. Missing responses are omitted with the bulk rows;
  tail bases such as `poly()` are trained on the complete exceedances.

* Hurdle quantiles distinguish the zero atom from the first positive count
  even when the log probabilities round together. Default POT threshold
  grids retain exactly represented counts beyond the 32-bit integer limit.
* Compound Poisson functions now report the scalar-parameter requirement
  directly when a vector of rates or severity parameters is supplied.

* Spliced fits check tail rank using factor levels observed in the bulk
  training sample. A factor and equivalent dummy columns now agree when
  a group has no exceedances, instead of returning a fit that cannot
  predict on its own training data.

* Hurdle fits select a joint sample using the evaluated terms and
  offsets of both components. Writing a factor transformation inline
  now gives the same fit as supplying the transformed column.
* Spliced fits omit incomplete bulk rows before selecting exceedances.
  At prediction time a missing bulk covariate or offset now affects only
  its own row. Bulk-only models use exactly `bulk_probs`: `threshold_prob`
  no longer adds a knot when there is no tail. Unknown tail controls and
  starts are rejected; use `start$beta` for log-scale coefficients and
  `fix.arg` to pin tail parameters.
* Spliced bulk quantile regressions now sum multiple offsets. An unseen
  fixed-effect group produces `NA` for that row without changing forecasts
  for known groups in the same batch.
* An unseen fixed-effect factor level gives `NA` for the affected
  prediction rows, without changing estimable rows. This also applies
  when a hurdle model observed a level but no positive count there.
* A discrete POT fit whose estimated support leaves at most two cells
  no longer reports separate Wald uncertainty for free scale and shape.
  The stability table marks the non-unique parameters as unavailable.
* Clarified the existing splice rule: a negative raw threshold puts the
  entire non-negative count distribution in the tail. Only the knots
  used for bulk interpolation are clipped to zero.
* Hurdle quantiles retain the positive gate probability on the log scale.
  Even when the reported zero probability rounds to one, an interior
  upper quantile or the support endpoint is no longer mistaken for zero.
  Distinct requested quantiles near one also have distinct column names.

* `fit_hurdle()` now retains the zero-truncation factor for DEGPD fits
  with very large scales. Previously a converged fit could report the
  untruncated likelihood and estimate a different positive-count shape.

* Spliced DEGPD fits stop if the tail scale design is rank deficient on
  the complete exceedances. Such fits could otherwise give incompatible
  forecasts from equally likely coefficients.
* An existing predictor named `.positive` (hurdle) or `exceedance`
  (splice) is no longer overwritten with the current response during
  fitting. Both fitting and prediction now use the supplied covariate.
* `fit_spliced_degpd()` now pairs each exceedance with its own complete
  design row, reports the number retained, and honors tail-formula offsets
  and omitted intercepts. Predictions reuse fitted polynomial bases and
  also work with just one bulk quantile.
* Spliced forecasts evaluate tail masses directly rather than subtracting
  CDFs rounded to one. Tail quantiles reach the bounded support endpoint
  at probability one and no longer hang where adjacent integers cannot
  be distinguished in double precision.
* `fit_hurdle()` evaluates a transformed response once and takes any
  `maxdata` subsample before splitting. Both components, the joint
  likelihood, and the reported sample size now describe the same rows.
  `degpd_pot_stability()` omits thresholds whose optimizer did not converge
  and names them in its diagnostic rather than including unreliable estimates.
* `fit_hurdle()` no longer returns a failed zero-truncated fit (log-likelihood
  -1e20) as converged: it tries data-driven starts and keeps the best
  converged fit. `family = "degpd"` starts the shape on the scale of the
  chosen link (it put 0.05 on log xi, i.e. xi = 1.05, and the shape usually
  collapsed to 0). Zero-truncated DEGPD standard errors use the exact
  observed information.
* Zero-truncated quantiles (`predict(type = "quantile")` of `fit_hurdle()`
  fits and of `"ztdegpd"`/`"ztgwaring"` fits) invert the base law on the
  log-survivor scale where needed: they were `Inf` when the base law put less
  than 1e-16 on positive counts (e.g. when extrapolating a covariate), and
  inaccurate within about 1e-15 of 1.
* `fit_hurdle()` finds the response on any element of a named formula list,
  as `egpd()` does (it had to come first), also when the names are in
  canonical order (`list(lsigma = ~ x, lxi = y ~ 1, lkappa = ~ 1)`), and the
  default `zero.formula` takes that formula's terms.
* `predict()` for spliced DEGPD fits no longer fails for a single `at` value,
  and for an integer threshold quantile q the tail starts at q + 1 in both
  fitting and prediction (it was fitted from q but given only the mass above
  q).
* `fit_spliced_degpd()` stops when the tail optimiser finds no parameter value
  with a positive likelihood (an infeasible `egpd_control$start` was returned
  as a converged fit with log-likelihood -1e20), and warns (class
  `"egpd_nonconvergence"`, with `tail_fit$loglik` `NA`) when it does not
  converge.

### Multivariate, compound and neural models

* `pcpegpd()` and `pcpdegpd()` now calculate large cumulative probabilities
  from the upper tail, including mass beyond the search grid. Previously the
  implied exceedance probability could be more than twice its true value,
  even when the corresponding quantile was correct.

* Gauge aggregation preserves small positive observations rather than
  snapping them to zero and omitting them from the positive likelihood.
  `rmdgpd()` keeps a strictly positive maximum even at subnormal scales;
  `rmdgpd()` and `rbegpd()` avoid overflow before rescaling finite draws.

* Multiscale fits optimize the shapes directly with box constraints in
  L-BFGS-B so that saturation of a logit does not masquerade as convergence.
  Near a shape bound the moving-block bootstrap compares warm starts with
  automatic starts on the actual resample, rather than reporting spuriously
  narrow intervals from frozen refits.
* Aggregate models honor interior shape starts even when bounds are close
  together or the start is close to a bound. `control = list(maxit = 0)`
  evaluates the requested start instead of moving it a fixed distance.
* Compound upper quantiles, including positive aggregate return levels,
  now invert the survivor directly. Probability beyond the search grid
  is retained, so an extreme but representable quantile inside `K` is
  not rounded to a lower bin or mistaken for grid exhaustion.

* Compound Poisson functions accept the type-2 `kappa = 0` boundary in
  both the gauge and floored severity laws; its carrier is exactly uniform.
* Compound Poisson log probabilities and positive-part likelihoods retain
  jump masses that underflow on the ordinary scale. `pcpdegpd()` floors
  noninteger queries; only gauge-discretised sums snap nearly on-grid
  observations. Blended GP log densities likewise retain finite beta-bulk
  and GP-tail log values below the ordinary probability range.
* `rmdgpd()` and `rzimdgpd()` respect the open lower endpoint at positive
  shape. Legitimate negative excesses are retained.
* Multiscale fits reject polynomial degrees that cannot be identified from
  the distinct observed durations. The monotonicity constraint is on scale
  and intensity: for `cpdegpd`, quantiles conditional on a positive sum can
  cross even when both parameters are non-decreasing.
* `rbegpd()` draws each pair of symmetric beta shares without subtracting
  a rounded near-one draw from one. Small representable angles no longer
  produce spurious zeros in the second continuous margin. Seeded bivariate
  draws change, including the discrete generators that use them.
* `fit_aggregate_egpd()` aggregates finest-scale series exactly and scores
  every window at its own Panjer bin: about 47% of gauge windows were scored
  one bin too high. Automatic starting values are screened (no kappa = 1e9
  or on-bound starts).
* `rcpegpd()` simulates the same gauge-discretised law as `dcpegpd()`,
  `pcpegpd()` and `qcpegpd()` (see Behaviour changes).
* `dcpegpd()` and the cpegpd likelihood take the severity bin masses where
  the CDF exceeds 1/2 from the log-survival, so far-tail lattice
  probabilities keep their precision (log-probabilities were off by up to 2.6
  for probabilities below 1e-17).
* Positive-part quantiles (`predict(type = "quantile")` of
  `fit_aggregate_egpd()` fits) are right for any event rate (they failed or
  were wrong for rates below about 1e-14), and `qcpegpd(1)` and `qcpdegpd(1)`
  are `Inf` for an unbounded support even when `exp(-lambda)` rounds to 1
  (for `qcpdegpd()` judged from the parameters, not from the grid masses,
  which can all underflow).
* Compound Poisson-discrete EGPD: the zero atom exp(-lambda P(Y >= 1)) takes
  P(Y >= 1) from the EGPD survival function, not as 1 - P(Y = 0), once
  P(Y = 0) > 0.999, so `dcpdegpd()`, `pcpdegpd()`, `qcpdegpd()` and the
  positive-part quantiles and log-likelihood of
  `fit_aggregate_egpd(family = "cpdegpd")` stay right when P(Y >= 1) is below
  about 1e-13 (the quantiles failed or were wrong, and the positive-part
  log-likelihood was `-Inf`).
* NPE posterior draws come from adaptive importance sampling on the learned
  log-ratio instead of resampling a fixed 20000-point prior grid, which gave
  an effective sample size of about 1; NPE works with NeuralEstimators.jl
  0.1.3 and 0.1.4 (it failed with an `Inf` log-ratio or a `MethodError`).
* `npe.control` values are range-checked (for example, `defensive` must be
  in [0, 1)): `defensive = 1` left every adaptive stage with an effective
  sample size of 0, and `defensive > 1` failed with a cryptic error. A
  `defensive` weight that leaves the fitted proposal no draws now draws that
  stage from the prior.
* Parallel helpers respect `_R_CHECK_LIMIT_CORES_`.

## Behaviour changes

* `egpd()` warns (class `"egpd_nonconvergence"`) whenever a fit does not
  converge and reports why in `fit$convergence.report`; starting values at
  which the likelihood is not finite are an error.
* Fits that stop where the likelihood is unbounded (a continuous EGPD shape
  below -1 with an observation at the upper end of the support) or on a
  saddle point are no longer reported as converged: they have
  `convergence = 3`, `logLik` `NA` and an `egpd_nonconvergence` warning, and
  `convergence.report` names the cause.
* `egpd()` warns (class `"egpd_singular_hessian"`, naming the coefficients)
  when the Hessian behind a converged fit's standard errors is singular and
  the standard errors cannot show it: `vcov()`, `summary()`, `confint()`,
  `predict(se.fit = TRUE)` and `simulate()` then understate the uncertainty
  along the unidentified direction(s). Examples are intercept-only generalized
  Waring fits at a = k, whose outer-product (BHHH) Hessian is singular because
  the two shape scores coincide (SE(log k) was 1.6e-8), collinear designs,
  collapsed model-2 mixtures and count-family parameters held at a bound.
  Predictors on a plateau of their link (`fit$plateau`) are not warned about,
  as their standard errors are already huge. `fit_hurdle()` reports the
  warning once, from the fit it keeps, and the vague "Final Hessian ... not
  numerically positive definite" message (which also leaked from restart
  candidates) is gone.
* DEGPD smooth fits with cells too narrow for the analytic third derivatives
  (very large counts or scales: (sigma + xi y)/|xi| above about 4.5e10) select
  their smoothing parameters by finite differences, with a warning naming the
  regime. The analytic REML gradient was up to 5e-2 off there (smoothing
  parameters about 5% off), and the REML Hessian behind `Vc` 29-330%. The
  finite-difference fallback warning now has class `"egpd_d34_fallback"`.
* `pgwaring()`, `qgwaring()` and `rgwaring()` warn (class
  `"egpd_gwaring_accuracy"`) when they use the generalized Waring tail
  integral beyond its validated range, `k > 1e15` at counts above 1000, where
  P(Y > q) can be off by up to about 1e-3 (lower-tail probabilities warn only
  where P(Y > q) > 1e-8). Values are unchanged; within the fitting bounds the
  tails are exact to about 1e-11.
* `confint(method = "profile")` warns (class `"egpd_profile_better_fit"`)
  when profiling finds a log-likelihood more than 0.05 above the fitted one,
  so the fit is not the maximum, as `MASS::profile.glm()` reports. The
  intervals are unchanged and relative to that fit; refit, e.g. from other
  `inits` or with `restarts`.
* `restarts` also covers `family = "zidegpd"`: all starts for carriers 2 to 4,
  otherwise only when the default fit fails or reaches a boundary.
* `restarts` also covers the continuous `family = "egpd"` (all carriers,
  m = 1-6). In addition to the five original starts, model 2 screens mixtures
  with separated component medians and tries up to three further full fits.
  This also applies to DEGPD and ZIDEGPD. It improves previously missed fits
  without adding a random search. Finite restarts cannot guarantee a global
  maximum; the continuous mixture with negative shape can even have an
  unbounded likelihood at the sample maximum. `restarts = FALSE` keeps the
  single-start fit.
* `egpd()`'s restarts (DEGPD, ZIDEGPD and the continuous EGPD) and
  `fit_hurdle()`'s starts keep the converged fit that is best by the
  criterion the fits optimise: without smooth terms the log-likelihood, as
  before; with smooth terms the REML criterion (the lowest `negREML`) when the
  smoothing parameters are estimated, and the penalised log-likelihood when
  `sp` is supplied. They were compared by the log-likelihood, which could
  keep a wigglier fit with a worse REML criterion (one continuous EGPD fit:
  0.18 higher log-likelihood, 1.15 worse `negREML`), so some fits with smooths
  now return another start's fit. Fits within 0.01 of the best count as the
  same optimum, and of those the one with the highest log-likelihood is kept
  (for the continuous EGPD, the default start's fit if it is among them). With
  estimated smoothing parameters, if any start has a predictor on a boundary
  or plateau of its link (as `fit$plateau` flags) or a penalised Hessian that
  is not positive definite, where the REML criterion is ill-defined, the fits
  are compared by the log-likelihood; so are sandwich-adjusted fits
  (`sandwich.args`), whose adjustment is estimated per start. When no start of
  `egpd()` converges, its warning says how many failed with an error.
* With `maxdata` or `maxspline` subsampling, a ZIDEGPD fit with default
  `inits` draws the same subsample as its DEGPD warm start, so for a given
  seed it fits a different subsample than before, and default restarts leave
  the random-number stream advanced by one fit's draws, not one per start.
* `gamma` has mgcv semantics: it only steers smoothing-parameter selection.
  `Vp` and `Vc` (and so `summary()`, `predict(se.fit = TRUE)`, `simulate()`
  and Wald intervals) and profile intervals use the unscaled likelihood: `Vp`
  was inflated by a factor `gamma`, and profile intervals solved a likelihood
  ratio of `gamma` times the chi-squared cutoff. The reported `sp` are those
  of the penalty on the unscaled likelihood, `gamma` times the previous
  values, and a supplied `sp` is read on that scale, so `sp = fit$sp`
  reproduces a fit whatever `gamma` is (`rho0` stays on the internal scale,
  log(sp / gamma)). Estimates change only for fits with a supplied `sp`;
  `gamma` must be a single positive number, and `gamma = 1` is unchanged.
* `negREML` is the Laplace approximation to the negative log restricted
  likelihood with mgcv's constants (Wood, 2011): for a likelihood mgcv also
  fits, it equals `gcv.ubre` of `mgcv::gam(method = "REML")` at the same
  smoothing parameters (with `gamma = 1`; see `?egpd`), also for fits with
  supplied `sp` (they reported the penalised negative log-likelihood). Its
  value changes for every fit with smooth terms (it carried a model-dependent
  offset); smoothing parameters, estimates and covariances do not.
* Where the Hessian of the negative log-likelihood is indefinite at a
  penalised optimum (as in zero-truncated DEGPD fits with smooths on all
  parameters, also within `fit_hurdle()`, and some DEGPD and continuous EGPD
  fits), the edf, `df`, AIC and BIC (and so `summary()`'s smooth-term tests
  and `fit_hurdle()`'s `df` and AIC) follow mgcv's eigenvalue clamp for
  general families: the plain formula could give a negative `df` (-51.5 for
  one fit) and an AIC below -2 logLik. `Vp` and `Vc` remain the exact Laplace
  curvature (the inverse of the exact penalised Hessian), and fits whose
  Hessian is positive semi-definite are unchanged.
* `negREML`, and so the selection of the smoothing parameters, is ill-defined
  where a linear predictor sits on a plateau of its link (`fit$plateau`, see
  New features): the Laplace approximation needs a regular maximum, and its
  value depends on how far the search drifted along the plateau (refitting a
  carrier-5 fit with kappa -> 0 at its reported `sp` moved `negREML` from
  447.15 to 463.22 at the same log-likelihood). `?egpd` documents this and
  `print()` notes it, without a warning.
* Named formula lists and `fixed` accept each parameter's natural name (e.g.
  `lambda`, `pi`, `sdlog`) and its link-prefixed forms in any order, and a
  formula list may carry the response on any element; `lpi` and `lxi` are
  deprecated aliases, and unknown names are an error listing the valid ones.
* Censored responses `cens(left, right)`, which no family supported, are an
  error.
* `summary()` labels parametric columns "z value" and "Pr(>|z|)", and smooth
  terms "p-value".
* `vcov()` and `confint()` default to `marginal = TRUE` (the smoothing-parameter
  corrected `Vc`), consistent with `predict()` and `simulate()`.
* `confint()` keeps the generalized Waring `rho` inside `rho > 1`: for
  `egpd()` fits a lower limit below the likelihood's bound 1 + 1e-6 ends there
  (one reached 0.106), and for `distfit()` fits the interval is formed for
  log(rho - 1) (in 13 of 20 ordinary fits the lower limit was below 1).
* `predict(type = "response")` adds its `xi` alias of `shape` only for the EGPD
  families (not CompPareto, whose `shape` is the body shape); select columns
  by name.
* For `family = "custom"` without `custom.fns$unlink`, `predict()` applies
  the inverse link named by a predictor's prefix, as `simulate()` and
  `confint()` do: with `type = "response"` a `logmu` predictor is
  exponentiated (it was returned on the link scale, named `mu`) and its
  `se.fit` is the response-scale delta-method standard error (it was the
  link-scale one), and `type = "quantile"` passes response-scale parameters
  to `custom.fns$q` (it passed the linear predictors, e.g. giving Poisson
  quantiles 1 and 3 instead of 4 and 6).
* Quantile columns of `predict()` and `simulate(probs = )` for `egpd()` fits,
  and of `predict()` for `fit_hurdle()` fits, are named by the probability to
  15 significant digits (`"q:0.9995"`): `round(prob, 3)` named 0.9995 and
  0.9999 both `"q:1"`, so selecting by name returned the wrong quantile.
  Names such as `"q:0.5"` and `"q:0.95"` are unchanged. `predict()` for
  `fit_spliced_degpd()` fits names its quantile and bulk columns the same
  way; it formatted the probabilities together, so 0.5 next to 0.95 was
  `"q:0.50"`.
* `predictive_coverage()` for a gamlss fit with `method = "parametric"` and
  `newdata` is an error when a parameter's formula has additive smoothers
  (`pb()`, `cs()`, `lo()`, `random()`, ...), which have no standard errors at
  new covariate values; use `method = "plug-in"` there.
* `fit_aggregate_egpd()` defaults to `optim.method = "auto"`: L-BFGS-B, which
  imposes the monotonicity constraint as box bounds, and BFGS followed by an
  L-BFGS-B polish, from the same start, keeping the better optimum
  (`fit$optim.runs`). The previous default, Nelder-Mead with 500 evaluations,
  stopped short of the optimum in nearly every test problem (by up to 2.4
  log-likelihood units), sometimes while reporting convergence, and a single
  L-BFGS-B run can stop on the flat ridges of the positive-part likelihood.
  Non-convergence warnings name the run and give the `optim` code and message,
  and `bootstrap_aggregate_egpd()` warns when refits do not converge.
* `rcpegpd()` gains `h` (default 0.2, as in `dcpegpd()`) and returns sums of
  event amounts each rounded up to a multiple of `h`; use `h = 0` for the
  continuous compound sum. `distfit(family = "cpegpd")` warns when positive
  data are off the `cpegpd.h` lattice and `h` is coarse relative to them.
* The EGPD and compound Poisson-EGPD distribution functions stop when a
  parameter required by `type` is missing (some returned `NA`); supplied `NA`
  parameters still give `NA`.
* Count distribution functions follow `stats::dpois`: `NA` gives `NA`, invalid
  parameters and probabilities outside [0, 1] give `NaN` with a warning (an
  `NA` or `NaN` count or probability is still returned silently, as
  `dpois(NA, -1)` is), and non-integer `x` gives 0 with a warning.
* Quantiles and tail probabilities change where they were inaccurate at
  extreme levels (see Bug fixes): count-family quantiles with `log.p = TRUE`
  at levels within about 1e-13 of 0, in either tail
  (`qgpois(-1e-20, 3, 0, log.p = TRUE)` is 30, as from `qpois()`, not 24),
  and `qzidiscegpd(log.p = TRUE)` there (3 counts higher); and, below the
  splice point of a CompPareto law whose Pareto tail weight is below about
  1e-3, the upper tail and `log.p` values near 0 of `pcomppareto()` and the
  `qcomppareto()` quantiles that invert them. Lower-tail quantiles on the
  linear scale, and seeded draws, are unchanged.
* `rpig()`, `rzipig()`, `rbell()`, `rzibell()`, `rBELL()` and `rZIBELL()`
  sample exactly, from the Poisson-inverse Gaussian mixture and from the
  Neyman type A representation of the Bell law, so seeded draws change; they
  inverted the distribution function up to a cap of 1e5, and draws above it
  were `NA` with a warning.
* The generalized Waring GAM bounds on rho and k are 1e15 (were 1e6 and 1e8),
  where the distribution has reached its negative-binomial and Poisson
  limits.
* GAMLSS families: an observation with zero likelihood at the current
  parameters gets a `NaN` score with a warning, so `gamlss()` stops instead of
  silently ignoring it; starting values are maximum-likelihood estimates of
  the intercept-only model.
* `fit_hurdle()` returns the best non-converged fit with a warning, `NA`
  log-likelihood and `converged = FALSE` when no start converges, and stops
  only when no start yields a finite log-likelihood.
* `distfit()` rejects negative data for the continuous families, and exact
  zeros for `family = "egpd"` (also with `method = "bernstein"`), where one
  zero pinned kappa at 1 with a spurious zero standard error: use `"ziegpd"`
  or a discrete family. Its `confint()` intervals are back-transformed from
  the log/logit scale.
* `degpd_pot_stability()` reports skipped thresholds and stops when every
  threshold fails (it returned empty data frames), and
  `fit_degpd_pot()$upper.endpoint` is the largest attainable count.
* `fit_aggregate_egpd()` rejects user starts outside `kappa_bounds` or
  `xi_bounds`.
* Blended GP functions require Beta weight shapes `c1, c2 > 3` (Majumder and
  Richards, p. 12).
* `train_begpd()`, `train_bdegpd()` and `train_mdgpd()` save to
  `tools::R_user_dir("egpd", "data")` by default (not a relative
  `inst/models`), and their `mc.cores` defaults respect the platform and
  `R CMD check` limits.

## New features

* Fits record `hessian.rank` (the rank of the Hessian behind the standard
  errors and the number of coefficients), and `print()`/`summary()` of `egpd`
  and `fit_hurdle()` fits add a one-line note when it is deficient.
* `bootstrap_aggregate_egpd()` gives moving-block bootstrap intervals for
  `fit_aggregate_egpd()` fits (Ailliot et al., Sec. 4.1), and
  `fit_aggregate_egpd(monotone_upto = D)` imposes the paper's monotonicity
  only on durations up to `D`.
* `confint()` gains `newdata` for parameter-level intervals at covariate
  values; `fitted()` gains `type = "response"`; `simulate(probs = )` returns
  simulated quantiles (it stopped with "not yet implemented").
* `predict()`, `summary()` and `nobs()` methods for `fit_hurdle()` objects.
* NPE fits store and print importance-sampling diagnostics (effective sample
  size, distinct draws, stages) and warn when the effective sample size is
  below 200; `distfit()` gains `npe.control`.
* Neural model sidecars record the training sample-size range and backend
  versions; `distfit(method = "neuralbayes")` warns for out-of-range sample
  sizes, out-of-support data and NeuralEstimators.jl/Flux.jl version
  mismatches.
* `rgpig()`, `rzigpig()`, `rGPIG()` and `rZIGPIG()` gain `max.value`, the
  largest count their inversion searches (default 1e5, as before, so default
  draws are unchanged); a draw above it is `NA` with a warning.
* `egpd()` fits gain `plateau`, a named logical vector with one element per
  linear predictor, `TRUE` where the predictor sits on a plateau (a flat end)
  of its link, so that the parameter is not identified: a predictor other
  than the scale's, not held by `fixed`, with a fitted |eta| > 9 (e.g. kappa
  or xi below about 1e-4 under a log link), or a DEGPD shape at either end of
  a bounded (`xi.max`) link. It is assessed for the EGPD, DEGPD, ZIDEGPD,
  ZT-DEGPD and (zero-truncated) generalized Waring families, and is `NA` for
  the others.

## Documentation

* The `doctor-visits`, `gaming-offenses` and `insurance-complaints` vignettes
  report the model comparisons at the fitted mixture likelihoods.
  `?gwaring` documents the tail accuracy limits.
* Help pages are generated from the roxygen comments: 23 hand-written pages
  (the GAMLSS families, `fit_aggregate_egpd()`, `fit_degpd_pot()`,
  `fit_hurdle()`, `fit_spliced_degpd()` and `rqresid()`), whose roxygen
  sources roxygen had been skipping and which had drifted from the code, were
  converted, and the package-level help cites the published versions of the
  methods papers.
* The `distfit-continuous` vignette describes the gauge-discretised
  compound Poisson-EGPD law (Ailliot, Gaetan and Naveau, 2026, Appendix F)
  that `dcpegpd()`/`rcpegpd()` implement, instead of a mean-preserving
  rounding the code no longer uses.
* Neural MDGPD estimation is documented as covering the homogeneous model only
  (common sigma, xi and generator rate), and the supported
  NeuralEstimators.jl/R package version pairings are listed: training through
  the R package NeuralEstimators 0.2.2 needs NeuralEstimators.jl 0.3 or later
  (it passes `shuffle`), and NeuralEstimators.jl 0.2.0 is not supported (its
  `logratio()` has no `grid` keyword). `inst/models/COMPATIBILITY.dcf` agrees
  (it claimed NeuralEstimators.jl 0.1-0.4), as do the README and both neural
  vignettes.
* The Julia set-up in the README, the `multivariate-egpd` and
  `bivariate-discrete-egpd` vignettes, `?train_begpd`, `?train_bdegpd` and
  `?train_mdgpd` installs BSON.jl, which the neural code loads, and
  Optimisers.jl, without which NeuralEstimators' `train()` fails
  (`Pkg.add(["NeuralEstimators", "Flux", "BSON", "Optimisers"])`). egpd loads
  both when it starts Julia, so a missing package is reported before any
  training data are simulated, and a failed Julia package load reports the
  Julia error with that command.

* The `bamlss-comparison`, `gamlss-comparison` and `predictive-coverage`
  vignettes no longer stop with "arguments imply differing number of rows",
  and the `simulation` vignette no longer mislabels its parameter tables.
  `predict(type = "response")` also returns an `xi` column that duplicates
  `shape`, so the vignettes now select its columns by name.
* The `simulation` vignette's Model 2 table reports
  `kappa2 = kappa1 + dkappa` (it showed `dkappa` labelled as `kappa2`) and
  notes that Model 2 is weakly identified.
* The `parameter-coverage` vignette loads precomputed Monte Carlo results
  (`parameter-coverage-results.rds`, generated by
  `vignettes/precompute/parameter-coverage.R` from the vignette's own code)
  instead of running about an hour of simulations during the build. Its
  summary is now computed from the results, and reports how often the
  profile likelihood interval could not be computed. Fits that do not
  converge count as failures.
* `bamlss-comparison` fits `bamlss()` with `sampler = FALSE`, since it only
  compares point estimates; this cuts its build time about tenfold.
* The README lists all 20 vignettes and has corrected citations: Gamet and
  Jalbert (2022) is *Environmetrics* 33(6), e2744; models 5-6 follow Ahmad
  and Arshad (2026, *Journal of Applied Statistics*; arXiv:2409.18719);
  Naveau et al. (2016) is pages 2753-2769; the first author of
  arXiv:2209.04660 is Le Carrer. DESCRIPTION cites the published versions
  by DOI.
* Vignette corrections: in `distfit-continuous`, cpegpd log-likelihoods are
  comparable across `h` (as `?distfit` says); in `gamlss-comparison`,
  `logLik()` of an egpd fit is the unpenalised log-likelihood at the
  penalised estimates; in `predictive-coverage`, `use_shortest = TRUE` gives
  the shortest contiguous interval, which need not be an HDR, even for a
  unimodal pmf, and the out-of-sample gamlss model uses `y ~ x`, as parametric
  coverage at new data requires; `multivariate-egpd` and
  `bivariate-discrete-egpd` cite Naveau et al. (2016) as pages 2753-2769.
* `?egpd` documents `custom.fns` (`d0`, `d120`, the optional `d340`, `unlink`
  and `q`; each `unlink` must be monotone, increasing or decreasing), and
  `?simulate.egpd` says that link-scale draws are named by parameter.
* The `rbzidegpd()` example says the joint-zero proportion is at least `pi0`
  (about 0.73 there), not approximately `pi0`.

## Packaging

* Vignette outputs (`inst/doc`) are no longer committed to git; `R CMD build`
  builds them.
* The long simulation studies moved from `tests/` to `simulations/`
  (excluded from the build), so `R CMD check` no longer runs them for hours.
  They now use the current parameter names (`logitpi`, and `xi` for the
  identity link), which had made their fits fail silently.
* Slow test blocks are skipped on CRAN with `skip_on_cran()`; the full suite
  still runs with `NOT_CRAN=true`.
* The BEGPD and BDEGPD neural tests train a tiny network with
  `train_begpd()`/`train_bdegpd()` when Julia and the Julia packages they need
  are available (not on CRAN), instead of always skipping: they needed bundled
  weights, which are not shipped.
* `.Rbuildignore` now excludes `_pkgdown.yml`, the derivation scripts and
  notes, `simulations/`, `vignettes/precompute/` and `.github/`.
* Added a GitHub Actions `R CMD check` workflow (Ubuntu and macOS, R
  release) and a README badge.
* `Deriv` is no longer a dependency; `numDeriv` and `withr` are suggested
  for the tests.
* DESCRIPTION depends on R >= 4.0.0, for `tools::R_user_dir()` (it said
  3.5.0), and suggests testthat >= 3.1.8, whose mocking of imported functions
  the tests use (it said 3.0.0).
* Added this `NEWS.md`.
