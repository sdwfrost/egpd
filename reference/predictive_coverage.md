# Predictive coverage for EGPD models

Compute simulation-based prediction intervals and their empirical
coverage for fitted EGPD-family models. Methods are provided for objects
of class `"egpd"`, `"gamlss"`, and `"bamlss"`.

## Usage

``` r
predictive_coverage(fit, ...)

# S3 method for class 'egpd'
predictive_coverage(
  fit,
  y,
  newdata = NULL,
  level = 0.95,
  nsim = 2000,
  method = c("parametric", "plug-in"),
  use_shortest = FALSE,
  ...
)

# S3 method for class 'gamlss'
predictive_coverage(
  fit,
  y,
  newdata = NULL,
  level = 0.95,
  nsim = 2000,
  method = c("parametric", "plug-in"),
  use_shortest = FALSE,
  data_fit = NULL,
  ...
)

# S3 method for class 'bamlss'
predictive_coverage(
  fit,
  y,
  newdata = NULL,
  level = 0.95,
  nsim = 2000,
  method = c("parametric", "plug-in"),
  use_shortest = FALSE,
  ...
)
```

## Arguments

- fit:

  a fitted model object (class `"egpd"`, `"gamlss"`, or `"bamlss"`)

- ...:

  additional arguments (currently unused)

- y:

  numeric vector of observed responses to evaluate coverage against

- newdata:

  optional data frame of covariate values at which to compute
  predictions. If `NULL` (the default), in-sample predictions are used.

- level:

  nominal coverage level (default 0.95)

- nsim:

  number of parameter/response simulation draws (default 2000)

- method:

  `"parametric"` (default) for parameter uncertainty, or `"plug-in"` for
  fixed-parameter prediction

- use_shortest:

  logical: use shortest (HDR) interval for discrete families instead of
  equal-tailed? Default `FALSE`.

- data_fit:

  (gamlss method only) optional data frame passed as `data` to
  [`predictAll()`](https://rdrr.io/pkg/gamlss/man/predict.gamlss.html);
  needed when the original fitting data cannot be recovered from the fit
  object.

## Value

A list with components:

- `family`:

  character string identifying the distribution family

- `method`:

  the method used (`"parametric"` or `"plug-in"`)

- `level`:

  the nominal coverage level

- `coverage`:

  the empirical coverage (proportion of `y` falling within the
  prediction interval)

- `covered`:

  logical vector indicating which observations are covered

- `L`:

  lower prediction interval bounds

- `U`:

  upper prediction interval bounds

## Details

The algorithm has three steps:

1.  **Parameter draws.**

    - `method = "plug-in"`:

      Fitted parameters are held fixed at their point estimates.

    - `method = "parametric"`:

      Parameters are perturbed on the link scale using normal noise
      (scaled by standard errors) and then back-transformed. For `egpd`
      fits this is a multivariate-normal draw from the joint coefficient
      posterior; for `gamlss` fits the perturbation is independent
      across parameters.

2.  **Response simulation.** `nsim` responses are drawn at each
    observation from the fitted distribution (using the appropriate `r*`
    function from the egpd package).

3.  **Intervals and coverage.** Equal-tailed quantile intervals are
    computed from the simulated responses. For discrete families,
    `use_shortest = TRUE` gives the shortest (highest-density-region)
    interval instead.

## Examples

``` r
if (FALSE) { # \dontrun{
y <- rdiscegpd(500, sigma = 3, xi = 0.15, kappa = 2, type = 1)
df <- data.frame(y = y)

# egpd fit
fit_e <- egpd(list(lsigma = y ~ 1, lxi = ~ 1, lkappa = ~ 1),
              data = df, family = "degpd", degpd.args = list(m = 1))
predictive_coverage(fit_e, y)
} # }
```
