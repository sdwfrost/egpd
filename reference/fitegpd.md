# Fit EGPD distribution to data

Maximum likelihood, Bernstein polynomial, or neural Bayes fitting of
EGPD, discrete EGPD, zero-inflated EGPD, zero-inflated discrete EGPD, or
bivariate multivariate EGPD distributions.

## Usage

``` r
fitegpd(
  x,
  type = 1,
  family = c("egpd", "degpd", "ziegpd", "zidegpd", "cpegpd", "cpdegpd", "begpd",
    "bdegpd", "bzidegpd", "mdgpd", "zimdgpd"),
  method = c("mle", "bernstein", "neuralbayes"),
  start = NULL,
  fix.arg = NULL,
  optim.method = "Nelder-Mead",
  hessian = TRUE,
  bernstein.m = 8,
  cpegpd.h = 0.2,
  model.path = NULL,
  estimator = c("npe", "nbe"),
  nsamples = 1000L,
  ...
)
```

## Arguments

- x:

  numeric vector of observations (univariate families), or an n-by-d
  numeric matrix/data.frame (multivariate families; d=2 for BEGPD, d \>=
  2 for MDGPD).

- type:

  integer 1-6 specifying the G-transformation type (univariate only).

- family:

  character: "egpd", "degpd", "ziegpd", "zidegpd", "cpegpd", "cpdegpd",
  "begpd" (bivariate EGPD), "bdegpd" (\\Experimental\\ bivariate
  discrete EGPD), "bzidegpd" (\\Experimental\\ zero-inflated bivariate
  discrete EGPD), "mdgpd" (\\Experimental\\ multivariate MDGPD via
  Aka-Kratz-Naveau construction, d \>= 2), or "zimdgpd"
  (\\Experimental\\ zero-inflated multivariate MDGPD, d \>= 2).

- method:

  character: "mle", "bernstein", or "neuralbayes". `"neuralbayes"` is
  required for bivariate families and requires Julia dependencies.

- start:

  named list of starting values, or NULL for automatic (not used for
  method="neuralbayes").

- fix.arg:

  named list of fixed parameters (not used for method="neuralbayes").

- optim.method:

  optimization method passed to
  [`optim`](https://rdrr.io/r/stats/optim.html)

- hessian:

  logical: compute standard errors via Hessian?

- bernstein.m:

  integer: Bernstein polynomial degree (method="bernstein" only)

- cpegpd.h:

  numeric: discretization step for cpegpd family.

- model.path:

  character: path to a pre-trained .bson model file
  (method="neuralbayes" only). If NULL, uses bundled model.

- estimator:

  character: `"npe"` for Neural Posterior Estimation or `"nbe"` for
  Neural Bayesian Estimation (method="neuralbayes" only).

- nsamples:

  integer: number of posterior samples for NPE (method="neuralbayes"
  with estimator="npe" only).

- ...:

  additional arguments passed to
  [`optim`](https://rdrr.io/r/stats/optim.html)

## Value

An object of class `"fitegpd"` with components:

- estimate:

  named vector of parameter estimates

- sd:

  named vector of standard errors (NA if Hessian not computed)

- vcov:

  variance-covariance matrix on natural scale

- loglik:

  maximized log-likelihood

- aic:

  Akaike information criterion

- bic:

  Bayesian information criterion

- n:

  number of observations

- npar:

  number of estimated parameters

- data:

  the input data vector

- type:

  the G-transformation type

- family:

  the distribution family

- method:

  the fitting method

- fix.arg:

  list of fixed arguments

- convergence:

  convergence code from optim (0 = success)

- optim:

  full optim output

- call:

  the matched call

- bernstein.m:

  Bernstein degree (NULL for method="mle")

- bernstein.weights:

  Bernstein weights (NULL for method="mle")

For `family="begpd"`, additional fields:

- estimator_type:

  "npe" or "nbe"

- model.path:

  path to the .bson model used

- posterior_samples:

  6 x nsamples matrix of posterior draws (NPE) or NULL (NBE)

- nsamples:

  number of posterior samples (NPE) or NULL (NBE)

## Examples

``` r
if (FALSE) { # \dontrun{
# Univariate fitting
x <- regpd(500, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
fit <- fitegpd(x, type = 1)
summary(fit)
plot(fit)

# Bivariate BEGPD (requires Julia)
Y <- rbegpd(1000, kappa = 2, sigma = 1, xi = 0.1, thL = 5, thU = 5, thw = 0.2)
fit_biv <- fitegpd(Y, family = "begpd", method = "neuralbayes")
summary(fit_biv)
plot(fit_biv)
} # }
```
