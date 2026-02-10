# Fit EGPD distribution to data

Maximum likelihood or Bernstein polynomial fitting of EGPD, discrete
EGPD, zero-inflated EGPD, zero-inflated discrete EGPD, or Compound
Poisson-EGPD distributions.

## Usage

``` r
fitegpd(
  x,
  type = 1,
  family = c("egpd", "degpd", "ziegpd", "zidegpd", "cpegpd"),
  method = c("mle", "bernstein"),
  start = NULL,
  fix.arg = NULL,
  optim.method = "Nelder-Mead",
  hessian = TRUE,
  bernstein.m = 8,
  cpegpd.h = 0.2,
  ...
)
```

## Arguments

- x:

  numeric vector of observations

- type:

  integer 1-6 specifying the G-transformation type

- family:

  character: `"egpd"`, `"degpd"`, `"ziegpd"`, `"zidegpd"`, or `"cpegpd"`

- method:

  character: `"mle"` or `"bernstein"` (Bernstein only for
  `family="egpd"`)

- start:

  named list of starting values, or `NULL` for automatic

- fix.arg:

  named list of fixed parameters (not optimized)

- optim.method:

  optimization method passed to
  [`optim`](https://rdrr.io/r/stats/optim.html)

- hessian:

  logical: compute standard errors via Hessian?

- bernstein.m:

  integer: Bernstein polynomial degree (`method="bernstein"` only)

- cpegpd.h:

  numeric: discretization step size for the Compound Poisson-EGPD Panjer
  recursion (`family="cpegpd"` only)

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

- cpegpd.h:

  discretization step size (NULL for non-cpegpd families)

## Details

The following families are supported:

- `"egpd"`:

  Continuous Extended Generalized Pareto Distribution

- `"degpd"`:

  Discrete EGPD for non-negative integer data

- `"ziegpd"`:

  Zero-inflated continuous EGPD

- `"zidegpd"`:

  Zero-inflated discrete EGPD

- `"cpegpd"`:

  Compound Poisson-EGPD for aggregated data with point mass at zero,
  fitted via Panjer recursion

For `family="cpegpd"`, the model is \\S = X_1 + \cdots + X_N\\ where \\N
\sim \mathrm{Poisson}(\lambda)\\ and \\X_i \sim \mathrm{EGPD}(\sigma,
\xi, \kappa, \ldots)\\. The `cpegpd.h` argument controls the
discretization grid width for the Panjer recursion; smaller values give
better accuracy at higher computational cost.

The Bernstein method (`method="bernstein"`) replaces the parametric
G-transformation with a flexible Bernstein polynomial density, giving a
semiparametric model.

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuous EGPD
x <- regpd(500, sigma = 2, xi = 0.1, kappa = 1.5, type = 1)
fit <- fitegpd(x, type = 1)
summary(fit)
plot(fit)

# Bernstein EGPD
fit_b <- fitegpd(x, type = 1, method = "bernstein")

# Discrete EGPD
y <- rdiscegpd(500, sigma = 3, xi = 0.1, kappa = 1.5, type = 1)
fit_d <- fitegpd(y, type = 1, family = "degpd")

# Compound Poisson-EGPD
z <- rcpegpd(500, sigma = 2, xi = 0.1, kappa = 1.5, lambda = 2)
fit_cp <- fitegpd(z, type = 1, family = "cpegpd")
} # }
```
