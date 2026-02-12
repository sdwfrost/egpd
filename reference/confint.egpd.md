# Confidence intervals for egpd model parameters

Confidence intervals for egpd model parameters

## Usage

``` r
# S3 method for class 'egpd'
confint(object, parm = NULL, level = 0.95, method = c("wald", "profile"), ...)
```

## Arguments

- object:

  A fitted `egpd` object.

- parm:

  Which parameters to compute CIs for. `NULL` (default) selects the
  intercept of each distributional parameter. Can also be a character
  vector of stripped parameter names (e.g. `"scale"`, `"shape"`) or an
  integer vector of coefficient indices.

- level:

  Confidence level (default 0.95).

- method:

  Either `"wald"` (default) or `"profile"`. Profile likelihood CIs are
  only available for intercept-only models (no smooth terms or
  covariates beyond an intercept).

- ...:

  Not used.

## Value

A matrix with columns `"lower"` and `"upper"`, on the response (natural)
scale. Rows are named by parameter.
