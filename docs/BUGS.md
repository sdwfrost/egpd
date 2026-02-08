# Bug Fixes in C++ Source Code

This document records bugs found in the C++ source files inherited from
the `degpd-and-zidegpd` prototype (itself a fork of `evgam` v1.0.1) and
fixed in the `egpd` package.

## Fixed Bugs

### 1. `gradHess.cpp:348` – Incorrect Hessian symmetry indices in `.gH6`

**Severity:** High **Affects:** All 6-parameter models (ZIDEGPD models 2
and 4 with smooth terms)

The `.gH6` function computes the gradient and Hessian for 6-parameter
models. After computing lower-triangle blocks `(5,4)` and `(6,4)`, the
code should transpose block `(6,4)` to fill the symmetric block `(4,6)`.
Instead, it wrote to the wrong block `(5,6)` and read wrong column
boundaries.

``` cpp
// BEFORE (buggy):
H.submat(s5, s6, e5, e6) = H.submat(s6, s4, e6, e5).t();

// AFTER (fixed):
H.submat(s4, s6, e4, e6) = H.submat(s6, s4, e6, e4).t();
```

This left the `(4,6)` block of the Hessian uninitialized, producing
incorrect second-derivative information for REML smoothing parameter
estimation in 6-parameter models.

------------------------------------------------------------------------

### 2. `degpd.cpp:372` – Wrong design matrix in `degpd2d12`

**Severity:** High **Affects:** DEGPD model 2 (mixture model) derivative
computation

The `degpd2d12` function (first and second derivatives for DEGPD model
2) used design matrix `X4` (for `lkappa2`) instead of `X5` (for
`logitp`) when computing the mixture proportion linear predictor. The
corresponding `degpd2d0` function correctly used `X5`.

``` cpp
// BEFORE (buggy):
arma::vec logitpvec = X4 * Rcpp::as<arma::vec>(pars[4]);

// AFTER (fixed):
arma::vec logitpvec = X5 * Rcpp::as<arma::vec>(pars[4]);
```

This produced silently wrong gradients and Hessians for DEGPD model 2
whenever `X4` and `X5` differed (e.g., different covariate structures
for `kappa2` and `p`), causing the optimizer to converge to incorrect
parameter estimates.

------------------------------------------------------------------------

### 3. `zi_degpd.cpp:380,458` – Wrong design matrix in `zidegpd2d0` and `zidegpd2d12`

**Severity:** High **Affects:** ZIDEGPD model 2 (mixture model with
zero-inflation) – both likelihood and derivative computation

Both `zidegpd2d0` and `zidegpd2d12` used design matrix `X5` (for
`logitp`) instead of `X6` (for `logitpi`) when computing the
zero-inflation probability. The function signatures correctly accepted
`X6` as a parameter, and `pars[5]` was the correct coefficient vector,
but it was multiplied by the wrong design matrix.

``` cpp
// BEFORE (buggy, in both zidegpd2d0 and zidegpd2d12):
arma::vec logitpivec = X5 * Rcpp::as<arma::vec>(pars[5]);

// AFTER (fixed):
arma::vec logitpivec = X6 * Rcpp::as<arma::vec>(pars[5]);
```

This is the same class of off-by-one design matrix error as Bug \#2,
affecting both the likelihood evaluation and its derivatives. When `X5`
and `X6` happened to be identical (e.g., both intercept-only), the bug
was masked.

------------------------------------------------------------------------

### 4. `zi_degpd.cpp:1664` – `0e5` typo in `zidegpd4d12` (should be `0.5`)

**Severity:** High **Affects:** ZIDEGPD model 4 second derivatives
(lsigma, lsigma) for y \> 0

In the `zidegpd4d12` function, the `(lsigma, lsigma)` second derivative
contained `0e5` instead of `0.5`. The literal `0e5` is valid C++ and
parses as `0.0 * 10^5 = 0.0`, silently zeroing out the entire second
term of the expression.

``` cpp
// BEFORE (buggy):
... * R_pow(ee5,2)) - 0e5 *

// AFTER (fixed):
... * R_pow(ee5,2)) - 0.5 *
```

The corresponding expression in `degpd4d12` (`degpd.cpp:1384`) correctly
uses `0.5`. This bug caused incorrect Hessian values for ZIDEGPD model
4, affecting optimizer convergence and standard error estimates.

------------------------------------------------------------------------

### 5. `families.R` – Incorrect `type` attribute for models 2, 3, and 4

**Severity:** Low (latent) **Affects:** All families (EGPD, DEGPD,
ZIDEGPD) models 2, 3, and 4

The `attr(family, "type")` value set in `.setup.family.egpd()` used the
fitting model number (2, 3, 4) instead of the distribution type number
expected by the `p.G`/`q.G`/`r.G` functions in `distributions.R`. The
distribution type numbering follows the original `degpd_functions.R`
prototype:

| Model | G transformation                       | Was | Should be |
|-------|----------------------------------------|-----|-----------|
| 1     | u^kappa                                | 1   | 1         |
| 2     | p\*u^kappa1 + (1-p)\*u^kappa2          | 2   | 6         |
| 3     | 1 - D_delta((1-u)^delta)               | 3   | 4         |
| 4     | \[1 - D_delta((1-u)^delta)\]^(kappa/2) | 4   | 5         |

The bug was latent because `predict.R` calls the `iG` inverse functions
directly from the likelihood function list rather than dispatching via
`q.G(..., type=...)`. The only runtime use of the `type` attribute was
in `setup_and_fit.R` to set initial values for model 2; this check was
also updated from `type == 2` to `type == 6`.

``` r

# BEFORE (families.R, model 2/3/4 for each family):
attr(family, "type") <- 2  # model 2
attr(family, "type") <- 3  # model 3
attr(family, "type") <- 4  # model 4

# AFTER (fixed):
attr(family, "type") <- 6  # model 2
attr(family, "type") <- 4  # model 3
attr(family, "type") <- 5  # model 4
```

------------------------------------------------------------------------

### 6. Model 2 label switching via reparameterization

**Severity:** High **Affects:** All model 2 variants (EGPD-2, DEGPD-2,
ZIDEGPD-2)

The mixture transformation G(u) = p*u^kappa1 + (1-p)*u^kappa2 is
invariant under swapping (kappa1, p) with (kappa2, 1-p), causing label
switching: the optimizer can converge to either mode, producing bimodal
parameter estimates and poor inference. Simulation tests showed only 34%
usable fits for EGPD-2.

**Fix:** Reparameterized from `(lkappa1, lkappa2)` to
`(lkappa1, ldkappa)` where `ldkappa = log(kappa2 - kappa1)`, enforcing
kappa2 \> kappa1 and breaking the symmetry. The C++ d0 functions compute
`lkappa2 = log(kappa1 + exp(ldkappa))` via numerically stable
log-sum-exp. The d12 derivative functions apply a chain rule
transformation at the end of each observation to convert derivatives
from the old parameterization to the new one.

Files modified: `egpd.cpp`, `degpd.cpp`, `zi_degpd.cpp`, `families.R`,
`setup_and_fit.R`, `predict.R`.

------------------------------------------------------------------------

### 7. `distributions.R` – Wrong type in `q.G` for type 6

**Severity:** Medium **Affects:** `q.G(..., type=6)` inverse
transformation for the mixture model

The root-finding function inside `q.G` for type 6 called
`p.G(..., type=4)` instead of `p.G(..., type=6)`, meaning it inverted
the wrong transformation. This affected `qegpd`, `qdiscegpd`, and
`qzidiscegpd` when called with `type=6`.

``` r

# BEFORE (buggy):
p.G(u = u, prob = prob, kappa = kappa, delta = delta, type = 4)

# AFTER (fixed):
p.G(u = u, prob = prob, kappa = kappa, delta = delta, type = 6)
```

------------------------------------------------------------------------

## Known Limitations (Not Fixed)

### `egpd.cpp` – No handling of xi = 0

All EGPD likelihood functions compute `1/xi` without checking for
`xi = 0`. When `xi = 0`, the GPD reduces to the exponential distribution
and requires a separate limiting form. The unused constant
`const double xieps = 0.0` at line 5 of `egpd.cpp` (and `degpd.cpp`,
`zi_degpd.cpp`) suggests this was planned but never implemented. In
practice, the shape parameter rarely lands exactly on zero during
optimization, but extreme parameter values could trigger `Inf`/`NaN`.

## External Bugs

### `bamlss` – `xbin_fun` declared as `void` causes “converting NULL pointer to R NULL” warnings

**Package:** bamlss (v1.2-5) **File:** `src/bamlss_functions.c:341`

The C function `xbin_fun` is declared as `void` but is registered in
`bamlss_init.c` as returning `SEXP` and called via
[`.Call()`](https://rdrr.io/r/base/CallExternal.html) in
`R/optimizers.R:2086`. The
[`.Call()`](https://rdrr.io/r/base/CallExternal.html) interface expects
all C functions to return an `SEXP` value. Since `xbin_fun` modifies its
`xweights` and `xrres` arguments in-place and has no return statement,
[`.Call()`](https://rdrr.io/r/base/CallExternal.html) receives whatever
happens to be in the return register – typically a null pointer. R \>=
4.4 emits a warning for each such call:

    Warning in .Call("xbin_fun", ...) : converting NULL pointer to R NULL

Since `xbin_fun` is called once per parameter per backfitting iteration,
intercept-only models with 3 parameters and 100 iterations produce ~300
warnings per `bamlss()` call.

**Fix:** Change the signature from `void xbin_fun(...)` to
`SEXP xbin_fun(...)` and add `return R_NilValue;` before the closing
brace. The init file already declares the correct `SEXP` return type.

**Workaround:** Suppress warnings when calling `bamlss()`:

``` r

fit <- suppressWarnings(bamlss(..., verbose = FALSE))
```

------------------------------------------------------------------------

### All files – Missing boundary checks in derivative functions

The `d0` (negative log-likelihood) functions in `degpd.cpp` include
basic boundary checks (e.g., `if (e3 <= 0) { nllh = 1e20; break; }`),
but the corresponding `d12` (derivative) and `d34` (higher derivative)
functions do not. If the optimizer explores regions where the GPD
support condition is violated, the derivative functions may produce
`NaN` values. In practice, the `d0` function returns a large penalty
value that steers the optimizer away before the derivatives are used.

### `egpd.cpp` – Dead code in `egpd2d34` and `egpd4d34`

The third/fourth derivative functions for EGPD models 2 and 4 (~2300
lines) are compiled but never called from R – the likelihood wrapper
sets `d340=NULL` for these models.
