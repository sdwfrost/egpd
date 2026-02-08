# Fit Extended Generalized Pareto Distribution GAMs

Fit Extended Generalized Pareto Distribution GAMs

## Usage

``` r
egpd(
  formula,
  data,
  family = "egpd",
  correctV = TRUE,
  rho0 = 0,
  inits = NULL,
  outer = "bfgs",
  control = NULL,
  removeData = FALSE,
  trace = 0,
  knots = NULL,
  maxdata = 1e+20,
  maxspline = 1e+20,
  compact = FALSE,
  egpd.args = list(),
  degpd.args = list(),
  zidegpd.args = list(),
  sandwich.args = list(),
  custom.fns = list(),
  sp = NULL,
  gamma = 1
)
```

## Arguments

- formula:

  a formula or list of formulae

- data:

  a data frame

- family:

  a character string: "egpd", "degpd", or "zidegpd"

- correctV:

  logical: should variance-covariance matrix account for smoothing
  parameter uncertainty? Defaults to TRUE

- rho0:

  initial log smoothing parameters

- inits:

  initial parameter values

- outer:

  outer optimization method: "bfgs" (default), "newton", "fd", or
  "fixed"

- control:

  a list of control parameters

- removeData:

  logical: should data be removed from the returned object? Defaults to
  FALSE

- trace:

  an integer controlling output verbosity

- knots:

  a list of knot values for smooth terms

- maxdata:

  maximum number of data rows

- maxspline:

  maximum number of rows for spline basis construction

- compact:

  logical: use compact representation? Defaults to FALSE

- egpd.args:

  a list of arguments for EGPD family (e.g., m=1)

- degpd.args:

  a list of arguments for DEGPD family (e.g., m=1)

- zidegpd.args:

  a list of arguments for ZIDEGPD family (e.g., m=1)

- sandwich.args:

  a list of sandwich correction arguments

- custom.fns:

  a list of custom likelihood functions

- sp:

  fixed smoothing parameters (if supplied, outer optimization is
  skipped)

- gamma:

  a gamma multiplier for the likelihood

## Value

An object of class `egpd`
